      program calc_wfn
      use m_LMTO
      implicit none
      integer :: nq_wfn
      double precision,allocatable :: q_wfn(:,:)
      integer :: nband_wfn
      integer,allocatable :: bindx_wfn(:)
      double complex,allocatable :: phipw(:,:,:,:,:,:),
     &     phiaug(:,:,:,:,:,:),phitot(:,:,:,:,:,:)
      integer :: npw(3)
      integer :: iq,ib
      logical :: lrho
c Read variables from DATA4GW, CphiGeig, QGphi, QGcou
      call setup_LMTO()

c Read number of bands and k-points for plotting phi
      read(*,*) nq_wfn
      allocate(q_wfn(3,nq_wfn))
      do iq=1,nq_wfn
        read(*,*) q_wfn(1:3,iq)
      enddo
      read(*,*) nband_wfn
      allocate(bindx_wfn(nband_wfn))
      read(*,*) bindx_wfn(1:nband_wfn)
      write(6,"(a,2i5)") '### nq_wfn,nband_wfn=',nq_wfn,nband_wfn

c Setup FFT mesh
      call calc_npw(npw)
c Allocate phi
      allocate(phipw(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp))
      allocate(phiaug(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp))
      allocate(phitot(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp))


      call calc_phiall(nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     &     npw,phipw,phiaug,phitot)


c Dump phi in OpenDX format
      write(*,*) 'dump phi'
      call wfn2dx(alat,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     &     npw,phipw,phiaug,phitot)
      call crystal2dx(alat,plat,nbas,bas,nclass,iclass,zz)


      write(*,*) 'dump rho? [T/F]'
      read(*,*) lrho
      if (lrho) call calc_rho(nq_wfn,nband_wfn,npw,phipw,phiaug,phitot)
      deallocate(phipw)
      deallocate(phiaug)
      deallocate(phitot)
      end program calc_wfn
cccccccccccccccccccccc
      subroutine findq(q,iqin)
      use m_LMTO
      implicit none

      double precision :: q(3)
      integer :: iqin
      double precision :: eps
      parameter (eps=1e-4)

      double precision :: dq(3)
      integer :: iq

      do iq=1,nqtt
        dq(1:3)=qtt(1:3,iq)-q(1:3)
        if (sum(abs(dq(1:3))).le.eps) then
          iqin=iq
          return
        endif
      enddo
      write(6,*) 'Error in findq: cannot find q = ',q(1:3)
      stop 'Error in findq: cannot find q'
      end subroutine findq
cccccccccccccccccccccc
      subroutine calc_phiall(nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     &     npw,phipw,phiaug,phitot)
      use m_LMTO
      use m_FFT3D
      implicit none

c inputs
      integer :: nq_wfn,nband_wfn,bindx_wfn(nband_wfn)
      double precision :: q_wfn(3,nq_wfn)
      integer :: npw(3)
c outputs
      double complex :: 
     &     phipw(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp)

c local
      integer :: isp,iq,iq_wfn,ib,i1,i2,i3
      integer :: augregion(4,npw(1)+1,npw(2)+1,npw(3)+1)

      double complex :: ! temp. array
     &     phipwtmp(npw(1)+1,npw(2)+1,npw(3)+1),
     &     phiaugtmp(npw(1)+1,npw(2)+1,npw(3)+1)

      double complex :: ! e^{ikr} and e^{ikT} (T is lattice vector)
     &     eikr(npw(1)+1,npw(2)+1,npw(3)+1),
     &     eikT(npw(1)+1,npw(2)+1,npw(3)+1)


      write(*,"(a)") '--- calc_phiall ---'
      call fft_init(npw,'B')
      call calc_augregion(npw(1),npw(2),npw(3),augregion)

      do iq_wfn=1,nq_wfn
        call findq(q_wfn(1:3,iq_wfn),iq)   
        call calc_eikreikT(q_wfn(1:3,iq_wfn),npw,augregion,eikr,eikT)
        do isp=1,nsp
        do ib=1,nband_wfn
          write(*,"(a,i2,2i5,3f10.4,i5)")
     &         '# isp,iq_wfn,iq,q,ib=',isp,iq_wfn,iq,qtt(1:3,iq),ib
          call calc_phi(isp,iq,bindx_wfn(ib),
     &         npw,
     &         augregion,phipwtmp,phiaugtmp)

          do i3=1,npw(3)+1
          do i2=1,npw(2)+1
          do i1=1,npw(1)+1
            phipw(i1,i2,i3,ib,iq_wfn,isp)=
     &           eikr(i1,i2,i3)
     &           *phipwtmp(i1,i2,i3)

            phiaug(i1,i2,i3,ib,iq_wfn,isp)=
     &           eikT(i1,i2,i3)
     &           *phiaugtmp(i1,i2,i3)

            phitot(i1,i2,i3,ib,iq_wfn,isp)=
     &           phipw(i1,i2,i3,ib,iq_wfn,isp)+
     &           phiaug(i1,i2,i3,ib,iq_wfn,isp)
          enddo !i3
          enddo !i2
          enddo !i1
        enddo !ib
        enddo !isp
      enddo !iq
      end subroutine calc_phiall
ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_npw(npw)
      use m_LMTO
      implicit none
c output
      integer :: npw(3)
c local
      integer :: iq,ig,id,itmp(3),ntmp(3)
      double precision :: pi,gtmp(3),gcutmax,gcuttmp
      write(*,"(a)") '--- calc_npw ---'
      pi=4.0d0*atan(1.0d0)

      ntmp(1:3)=0

      do iq=1,nqnum
        gcutmax=-1.0d0
        do ig=1,ngp(iq)
          call mymatvec(G,dble(ngvecp(1:3,ig,iq)),gtmp,3,3)
          gtmp(1:3)=gtmp(1:3)+qtt(1:3,iq)
          gtmp(1:3)=gtmp(1:3)*2.0d0*pi/alat
          gcuttmp=sqrt(sum(gtmp(1:3)**2))
          if (gcutmax.lt.gcuttmp) gcutmax=gcuttmp
          do id=1,3
            itmp(id)=abs(ngvecp(id,ig,iq))
            if (ntmp(id).lt.itmp(id)) ntmp(id)=itmp(id)
          enddo
        enddo
        write(*,"(a,i5,a,f10.5)") '# iq=',iq,' gcutmax=',gcutmax
      enddo

      npw(1:3)=2*ntmp(1:3)+2
      write(*,"(a,3i6)") '# npw(1:3)=',npw(1:3)
      end subroutine calc_npw
ccccccccccccccccccccccccccccccccccccc
      subroutine calc_augregion(n1,n2,n3,augregion)
      use m_LMTO
      implicit none
c input
      integer :: n1,n2,n3
c output
      integer :: augregion(4,n1+1,n2+1,n3+1)
c local
      integer :: nshell
      parameter (nshell=4)
      integer :: i1,i2,i3,j1,j2,j3,ibas,ic
      double precision :: rmax,ratom(3),r(3),rtmp(3),dr

      write(*,*) '--- calc_augregion ---'
      augregion(:,:,:,:)=0

      do ibas=1,nbas
        ic=iclass(ibas)
        rmax = bb(ic)*(exp((nr(ic)-1)*aa(ic))-1d0)
        do j1=-nshell,nshell
        do j2=-nshell,nshell
        do j3=-nshell,nshell
          rtmp(1)=j1
          rtmp(2)=j2
          rtmp(3)=j3
          call mymatvec(plat,rtmp,ratom,3,3)
          ratom(1:3)=alat*(ratom(1:3)+bas(1:3,ibas))

          do i3=1,n3+1          
          do i2=1,n2+1
          do i1=1,n1+1

            rtmp(1)=(i1-1)/dble(n1)
            rtmp(2)=(i2-1)/dble(n2)
            rtmp(3)=(i3-1)/dble(n3)
            call mymatvec(plat,rtmp,r,3,3)
            r(1:3)=alat*r(1:3)
            dr=sqrt(sum((r(1:3)-ratom(1:3))**2))
          if (dr.lt.rmax) then
            if (augregion(4,i1,i2,i3).ne.0) then
              stop 'Overlap in augmented region!'
            endif
            augregion(1,i1,i2,i3)=j1
            augregion(2,i1,i2,i3)=j2
            augregion(3,i1,i2,i3)=j3
            augregion(4,i1,i2,i3)=ibas
          endif
        enddo !i1
        enddo !i2
        enddo !i3
        enddo !j3
        enddo !j2
        enddo !j1
      enddo !ibas
      end subroutine calc_augregion
ccccccccccccccccccccccccccccccccccccc
      subroutine calc_eikreikT(kvec,npw,augregion,eikr,eikT)
      use m_LMTO
      implicit none
c input
      double precision :: kvec(3)
      integer :: npw(3),augregion(4,npw(1)+1,npw(2)+1,npw(3)+1)
c output
      double complex :: ! e^{ikr} and e^{ikT} (T is lattice vector)
     &     eikr(npw(1)+1,npw(2)+1,npw(3)+1),
     &     eikT(npw(1)+1,npw(2)+1,npw(3)+1)
c local
      integer :: i1,i2,i3
      double precision :: rtmp(3),r(3),tvec(3)
      double precision :: phase,pi

      pi=4.0d0*atan(1.0d0)

c Calculate e^{ikr}
      do i3=1,npw(3)+1
      do i2=1,npw(2)+1
      do i1=1,npw(1)+1
        rtmp(1)=(i1-1)/dble(npw(1))
        rtmp(2)=(i2-1)/dble(npw(2))
        rtmp(3)=(i3-1)/dble(npw(3))
        call mymatvec(plat,rtmp,r,3,3)
        r(1:3)=alat*r(1:3)
        phase=2.0d0*pi/alat*sum(kvec(1:3)*r(1:3))
        eikr(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
      enddo
      enddo      
      enddo

c Calculate e^{ikT}
      do i3=1,npw(3)+1
      do i2=1,npw(2)+1
      do i1=1,npw(1)+1

        if (augregion(4,i1,i2,i3).ne.0) then
          rtmp(1:3)=augregion(1:3,i1,i2,i3)
          call mymatvec(plat,rtmp,tvec,3,3)          
          tvec(1:3)=alat*tvec(1:3)
          phase=2.0d0*pi/alat*sum(kvec(1:3)*tvec(1:3))
          eikT(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
        else
          eikT(i1,i2,i3)=0.0d0
        endif
      enddo
      enddo      
      enddo

      end subroutine calc_eikreikT
ccccccccccccccccccccccccccccccccccccc
      subroutine calc_phi(isp,iq,iband,npw,
     &     augregion,phipwtmp,phiaugtmp)
      use m_LMTO
      use m_FFT3D
      implicit none
c input
      integer :: isp,iq,iband,npw(3)
      integer :: augregion(4,npw(1)+1,npw(2)+1,npw(3)+1)

c output
      double complex :: 
     &     phipwtmp(npw(1)+1,npw(2)+1,npw(3)+1),
     &     phiaugtmp(npw(1)+1,npw(2)+1,npw(3)+1)

c local
      integer :: itmp(3),ig,id,i1,i2,i3

      double precision :: rtmp(3),r(3) !points to plot
      double precision :: ratom(3) ! atomic points
      double precision :: dr(3)


c      write(*,*) '--- calc_phi ---'

c Plane-wave part
      in_fft(:,:,:)=0.0d0
      do ig=1,ngp(iq)
        do id=1,3
          itmp(id)=ngvecp(id,ig,iq)
          if (itmp(id).lt.0) itmp(id)=itmp(id)+npw(id)
        enddo
        in_fft(itmp(1)+1,itmp(2)+1,itmp(3)+1)=geig(ig,iband,iq,isp)
      enddo
c Perform backward fft
      call myfft()

      do i3=1,npw(3)+1
      do i2=1,npw(2)+1
      do i1=1,npw(1)+1
        phipwtmp(i1,i2,i3)=out_fft(mod(i1-1,npw(1))+1,
     &       mod(i2-1,npw(2))+1,mod(i3-1,npw(3))+1)
      enddo
      enddo
      enddo

c Augmented part
      phiaugtmp(:,:,:)=0.0d0
      do i3=1,npw(3)+1
      do i2=1,npw(2)+1
      do i1=1,npw(1)+1

        if (augregion(4,i1,i2,i3).ne.0) then
c set plane-wave part to zero
          phipwtmp(i1,i2,i3)=0.0d0

          rtmp(1)=(i1-1)/dble(npw(1))
          rtmp(2)=(i2-1)/dble(npw(2))
          rtmp(3)=(i3-1)/dble(npw(3))
          call mymatvec(plat,rtmp,r,3,3)
          r(1:3)=alat*r(1:3)

          rtmp(1:3)=augregion(1:3,i1,i2,i3)
          call mymatvec(plat,rtmp,ratom,3,3)          
          ratom(1:3)=alat*(ratom(1:3)+bas(1:3,augregion(4,i1,i2,i3)))
          dr(1:3)=r(1:3)-ratom(1:3)
          call calc_phiaug(dr,augregion(4,i1,i2,i3),
     &         phiaugtmp(i1,i2,i3),isp,iq,iband)
        endif
      enddo
      enddo
      enddo
      end subroutine calc_phi
cccccccccccccccccccccccccccccccccccc
c Calculate wavefunction phi in the augmented region
c

      subroutine calc_phiaug(dr,ibas,phia,isp,iq,iband)
      use m_LMTO
      implicit none
c input
      double precision :: dr(3)
      integer :: ibas,isp,iq,iband
c output
      double complex :: phia
c
c local
      double precision :: pi
      double precision :: drlength,theta,pphi,sintheta
      integer :: idim,il,mtmp,ntmp,ltmp
c
c parameters for Ylm
      integer :: lmax
      parameter (lmax=6)
      double complex :: Y(2*lmax+1,lmax+1)
      double precision :: Yreal(2*lmax+1,lmax+1)
c function
      double precision :: calc_gxr

c
      pi=4.0d0*atan(1.0d0)
      drlength=sqrt(sum(dr(1:3)**2))
c
c x=r*sin(theta)*cos(pphi)
c y=r*sin(theta)*sin(pphi)
c z=r*cos(theta)
      theta=acos(dr(3)/(drlength+1.0d-15))
      sintheta=sqrt(1.0d0-cos(theta)**2)
      pphi=acos(dr(1)/(drlength*sintheta+1.0d-15))
      if (dr(2).lt.0.0d0) pphi=2*pi-pphi


      do il=0,lmax
        call calc_Ylm(il,theta,pphi,
     &       Y(1:2*il+1,il+1),
     &       Yreal(1:2*il+1,il+1))
      enddo

      phia=0.0d0
      do idim=1,ldim2
        if (mnla(4,idim).eq.ibas) then
          mtmp=mnla(1,idim)
          ntmp=mnla(2,idim)
          ltmp=mnla(3,idim)
          if (ltmp.gt.lmax) then
            stop 'ltmp.gt.lmax!'
          endif
          phia=phia+cphi(idim,iband,iq,isp)
     &         *calc_gxr(drlength,ltmp,ntmp,iclass(ibas),isp)
     &         *Yreal(mtmp+ltmp+1,ltmp+1)
        endif
      enddo
      end 
cccccccccccccccccccccccccccccccc
c Linear interpolation of gx/r
c
      double precision function calc_gxr(r,l,n,ic,isp)
      use m_LMTO      
      implicit none
c input
      double precision :: r
      integer :: l,n,ic,isp
c local
      double precision :: r1,r2
      integer :: ir

      ir=1+int(log(r/bb(ic)+1.0d0)/aa(ic))
      if (ir.lt.1) stop 'ir.lt.1'
      if (ir.gt.nr(ic)-1) stop 'ir.gt.nr(ic)-1'

      r1=bb(ic)*(exp((ir-1)*aa(ic))-1d0)
      r2=bb(ic)*(exp((ir  )*aa(ic))-1d0)
      if (r1.gt.r) stop 'r1.gt.r'
      if (r2.le.r) stop 'r2.le.r'


      calc_gxr=(r-r2)/(r1-r2)*gx(ir,l,n,ic,isp) 
     &     + (r-r1)/(r2-r1)*gx(ir+1,l,n,ic,isp)
      calc_gxr=calc_gxr/(r+1.0d-20)

      end function calc_gxr
cccccccccccccccccccccccccccccc
      subroutine calc_rho(nq_wfn,nband_wfn,npw,phipw,phiaug,phitot)
      use m_LMTO
      implicit none
c input
      integer :: nq_wfn,nband_wfn,npw(3)
      double complex :: 
     &     phipw(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp)

      double complex,allocatable :: rhopw(:,:,:),
     &     rhoaug(:,:,:),rhotot(:,:,:)

c
      integer :: isp,iq,ib,i1,i2,i3
      double precision :: nel

      double precision :: qdum(3)
      integer :: bindxdum
      write(*,*) '--- calc_rho ---'
c Allocate rho
      allocate(rhopw(npw(1)+1,npw(2)+1,npw(3)+1))
      allocate(rhoaug(npw(1)+1,npw(2)+1,npw(3)+1))
      allocate(rhotot(npw(1)+1,npw(2)+1,npw(3)+1))

      rhopw(1:npw(1)+1,1:npw(2)+1,1:npw(3)+1)=0.0d0
      rhoaug(1:npw(1)+1,1:npw(2)+1,1:npw(3)+1)=0.0d0
      rhotot(1:npw(1)+1,1:npw(2)+1,1:npw(3)+1)=0.0d0

      do isp=1,nsp
      do iq=1,nq_wfn
      do ib=1,nband_wfn
        do i3=1,npw(3)+1
        do i2=1,npw(2)+1
        do i1=1,npw(1)+1
          rhopw(i1,i2,i3)=rhopw(i1,i2,i3)+
     &         abs(phipw(i1,i2,i3,ib,iq,isp))**2
          rhoaug(i1,i2,i3)=rhoaug(i1,i2,i3)+
     &         abs(phiaug(i1,i2,i3,ib,iq,isp))**2
          rhotot(i1,i2,i3)=rhotot(i1,i2,i3)+
     &         abs(phitot(i1,i2,i3,ib,iq,isp))**2

        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      rhopw(:,:,:)=
     &     rhopw(:,:,:)/dble(nq_wfn)
      rhoaug(:,:,:)=
     &     rhoaug(:,:,:)/dble(nq_wfn)
      rhotot(:,:,:)=
     &     rhotot(:,:,:)/dble(nq_wfn)


      nel=0.0d0
      do i3=1,npw(3)
      do i2=1,npw(2)
      do i1=1,npw(1)
        nel=nel+rhotot(i1,i2,i3)
      enddo
      enddo
      enddo
      nel=nel*dble(3-nsp)*VOL/dble(npw(1)*npw(2)*npw(3))
      write(*,*) 'nel = ',nel

      qdum(1:3)=0.0d0
      bindxdum=0
      call wfn2dx(alat,plat,1,1,1,qdum,bindxdum,
     &     npw,rhopw,rhoaug,rhotot)

      deallocate(rhopw)
      deallocate(rhoaug)
      deallocate(rhotot)
      end subroutine calc_rho
