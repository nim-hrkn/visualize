      module m_LMTO
      use m_DATA4GW
      use m_QG
      implicit none
      double precision :: G(3,3),VOL
      double precision,allocatable :: Saug(:,:,:)
      double complex,allocatable :: Spw(:,:,:)

      contains
ccccccccccccccccccccccccccccccccccccccc
      subroutine setup_LMTO()
      implicit none
      double precision :: At(3,3)
      write(6,*) '--- setup_LMTO ---'
      call read_DATA4GW()
      call read_CphiGeig()
      call set_mnla()
      call read_QG()

      call mydet3(plat,VOL)
      VOl=abs(VOL)*alat**3

      call mytranspose(plat,At,3,3)
      call myinv3(At,G)

      write(6,*) 'VOL=',VOL
      write(6,*) 'plat='
      write(6,"(3f10.5)") plat(1:3,1)
      write(6,"(3f10.5)") plat(1:3,2)
      write(6,"(3f10.5)") plat(1:3,3)
      write(6,*) 'G='
      write(6,"(3f10.5)") G(1:3,1)
      write(6,"(3f10.5)") G(1:3,2)
      write(6,"(3f10.5)") G(1:3,3)
c      call set_Saug()
c      call set_Spw()
      end subroutine setup_LMTO
ccccccccccccccccccccccccccccccccccccccc
      subroutine set_Saug()
      implicit none
      double precision :: stmp
      integer :: m1,n1,l1,a1
      integer :: m2,n2,l2,a2
      integer :: ix1,ix2,isp
      integer :: ic

      write(6,*) '--- set_Saug ---'
      allocate(Saug(ldim2,ldim2,nsp))

      Saug(1:ldim2,1:ldim2,1:nsp)=0.0d0
      do isp=1,nsp
        write(6,*) 'isp,nsp = ',isp,nsp
      do ix1=1,ldim2
        m1=mnla(1,ix1)
        n1=mnla(2,ix1)
        l1=mnla(3,ix1)
        a1=mnla(4,ix1)
        do ix2=ix1,ldim2
          m2=mnla(1,ix2)
          n2=mnla(2,ix2)
          l2=mnla(3,ix2)
          a2=mnla(4,ix2)
          if (m1.ne.m2) cycle
          if (l1.ne.l2) cycle
          if (a1.ne.a2) cycle
          ic=iclass(a1)
cgx (nrmx, 0:lmxamx, nphimx, nclass,nsp)
          call GINTxx(gx(1:nr(ic),l1,n1,ic,isp),
     &         gx(1:nr(ic),l1,n2,ic,isp),
     &         aa(ic),bb(ic),nr(ic),stmp)
          Saug(ix1,ix2,isp)=stmp
          Saug(ix2,ix1,isp)=stmp
          write(6,"(4i3,f8.4)") a1,l1,n1,n2,stmp
        enddo
      enddo ! ix1
      enddo ! isp
      end subroutine set_Saug
cccccccccccccccccccccccccccccccccccc
      SUBROUTINE GINTxx(G1,G2,A,B,NR,SUM)
C- Integrate product of two wave functions, Simpson rule
C ----------------------------------------------------------------
c takao \sum_i g(i)*g(i) drdi
      implicit none
      integer :: NR
      double precision :: G1(NR),G2(NR),A,B,SUM
c local
      integer :: IR
      double precision :: EA2,EA4,DRDI,RMPB
      EA2 = DEXP(A + A)
      EA4 = EA2*EA2
      SUM = 0D0
      DRDI = (A*B)*DEXP(A)
      DO  10  IR = 2, NR-1, 2
        SUM = SUM + G1(IR)*G2(IR)*DRDI
        DRDI = DRDI*EA2
  10  CONTINUE
      SUM = SUM + SUM
      DRDI = (A*B)*EA2
      DO  11  IR = 3, NR-2, 2
        SUM = SUM + G1(IR)*G2(IR)*DRDI
        DRDI = DRDI*EA2
  11  CONTINUE
      RMPB = B*DEXP(A*(NR-1))
      SUM = (2*SUM+ G1(1)*G2(1)*(A*B) + G1(NR)*G2(NR)*(A*RMPB))/3d0
      END SUBROUTINE GINTxx
cccccccccccccccccccccccccccccccccccccccc
      subroutine set_Spw()
      implicit none
c local
      integer :: gtmp(3),absgtmp
      integer :: iq,ig1,ig2
      integer :: ibas,ic
      double precision :: rmax,fpi
      double complex :: img
      double precision :: dgtmp(3),G21(3),absG21,G21R

      write(6,*) '--- set_Spw ---'
      allocate(Spw(ngpmx,ngpmx,nqnum))
      Spw(1:ngpmx,1:ngpmx,1:nqnum)=0.0d0
      img=dcmplx(0.0d0,1.0d0)
      fpi=16.0d0*atan(1.0d0)
      do iq=1,nqnum
        do ig2=1,ngp(iq)
          do ig1=ig2,ngp(iq)
            gtmp(1:3)=ngvecp(1:3,ig2,iq)-ngvecp(1:3,ig1,iq)
            dgtmp(1:3)=0.5d0*fpi/alat*dble(gtmp(1:3))
            call mymatvec(G,dgtmp,G21,3,3)
            absG21=sqrt(sum(G21(1:3)**2))

            absgtmp=gtmp(1)**2+gtmp(2)**2+gtmp(3)**2
            if (absgtmp.eq.0) then
              Spw(ig1,ig2,iq)=VOL
            endif
            do ibas=1,nbas
              ic=iclass(ibas)
              rmax = bb(ic)*(exp((nr(ic)-1)*aa(ic))-1d0)
              if (absgtmp.eq.0) then
                Spw(ig1,ig2,iq) = Spw(ig1,ig2,iq) - fpi*rmax**3/3d0
              else
                G21R=absG21*rmax
                Spw(ig1,ig2,iq) = Spw(ig1,ig2,iq)
     &               - exp( img* sum(G21(1:3)*bas(1:3,ibas))*alat )  
     &               * fpi/absG21**3*( -G21R * cos(G21R) + sin(G21R))
              endif
            enddo !ibas

            Spw(ig2,ig1,iq)=conjg(Spw(ig1,ig2,iq))
          enddo !ig2
        enddo !ig1
      enddo !iq
      end subroutine set_Spw
cccccccccccccccccccccccccccccccccccccccc
      end module m_LMTO
