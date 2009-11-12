      program calc_wannier4
c 080603, For plotting arbitary region, from calc_wannier
c 071114, Takashi Miyake, from calc_wfn.F
      use m_LMTO
      use m_MLWF
      use m_cubeformat 
      use m_xsfformat 
      implicit none
      integer :: nq_wfn,tlat(3),nfac
      double precision :: qdum(3)
      double precision,allocatable :: q_wfn(:,:)
      integer :: nband_wfn
      integer,allocatable :: bindx_wfn(:)
      double complex,allocatable :: phipw(:,:,:,:,:,:),
     &     phiaug(:,:,:,:,:,:),phitot(:,:,:,:,:,:)
      double complex, allocatable :: wanpw(:,:,:,:,:),
     &     wanaug(:,:,:,:,:),wantot(:,:,:,:,:)
      double precision :: rini(3),rfin(3),ang,alat_ang
      integer :: npw(3),mesh(3)
      integer :: iq,ib
      logical :: lrho

      character(20):: outputformat ='none'

c Read variables from DATA4GW, CphiGeig, QGphi, QGcou
      call setup_LMTO()
c Read variables from MLWU/D
      call setup_MLWF()

      write(6,*) '### ib,bas(1:3,ib)'
      do ib=1,nbas
         write(*,"(i5,3f12.6)")ib,bas(1:3,ib)
      enddo
c      call chk_weight()
c      stop 'chk_weight'

c Read number of bands and k-points for plotting phi
!      read(*,*) nq_wfn
!      allocate(q_wfn(3,nq_wfn))
!      do iq=1,nq_wfn
!        read(*,*) q_wfn(1:3,iq)
!      enddo
      nq_wfn = nqbz
      allocate(q_wfn(3,nq_wfn))
      q_wfn(1:3,1:nqbz) = qtt(1:3,1:nqbz)

      read(*,*) nband_wfn
      allocate(bindx_wfn(nband_wfn))
      read(*,*) bindx_wfn(1:nband_wfn)
!      do ib=1,nband_wfn
!      write(*,*)ib,bindx_wfn(ib)
!      enddo
      write(6,"(a,2i5)") '### nq_wfn,nband_wfn=',nq_wfn,nband_wfn
      if (nband_wfn.gt.nband) stop 'nband_wfn.gt.nband !'
      read(*,*) tlat(1:3)
      write(*,"(a,3i5)")'### tlat',tlat


c Setup FFT mesh
      read(*,*)nfac
      call calc_npw(nfac,npw)
c Allocate phi
      read(*,*)mesh(1),mesh(2),mesh(3)
      read(*,*)rini(1),rini(2),rini(3)
      read(*,*)rfin(1),rfin(2),rfin(3)
      allocate(phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,
     &         nband_wfn,nq_wfn,nsp))
      allocate(phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,
     &         nband_wfn,nq_wfn,nsp))
      allocate(phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,
     &         nband_wfn,nq_wfn,nsp))

      call calc_phiall_2(nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     &     npw,mesh,rini,rfin,phipw,phiaug,phitot)


c from Bloch to Wannier
      allocate(wanpw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp))
      allocate(wanaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp))
      allocate(wantot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp))
      call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phipw,
     &           wanpw)
      call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phiaug,
     &           wanaug)
      call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phitot,
     &           wantot)

c Dump phi in OpenDX format
      outputformat='xsf'

      write(*,*) 'dump phi'
      qdum = 0.0d0
      ang=0.529177d0
      alat_ang=alat*ang
      if (outputformat.eq.'opendx') then
      call wfn2dx_2(alat_ang,plat,nsp,1,nband_wfn,qdum,bindx_wfn,
     &     mesh,rini,rfin,wanpw,wanaug,wantot)
      call crystal2dx_2(alat_ang,plat,rini,rfin,
     &     nbas,bas,nclass,iclass,zz)
      else if (outputformat.eq.'gaussiancube') then

      call wrt_cube(
     i     'psi',
     i     alat,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
c     i     mesh,rini,rfin,phipw,phiaug,phitot  ! for bloch orbital
     i     mesh,rini,rfin,wanpw,wanaug,wantot,  ! for wannier function 
     i     nbas,bas,nclass,iclass,zz )

      call wrt_cube(
     i     'phi',
     i     alat,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
     i     mesh,rini,rfin,phipw,phiaug,phitot,  ! for bloch orbital
c     i     mesh,rini,rfin,wanpw,wanaug,wantot,  ! for wannier function 
     i     nbas,bas,nclass,iclass,zz )

      else i f(outputformat.eq.'xsf') then

      call wrt_xsf(
     i     'psi',
     i     alat_ang,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
c     i     mesh,rini,rfin,phipw,phiaug,phitot  ! for bloch orbital
     i     mesh,rini,rfin,wanpw,wanaug,wantot,  ! for wannier function 
     i     nbas,bas,nclass,iclass,zz )

      call wrt_xsf(
     i     'phi',
     i     alat_ang,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
     i     mesh,rini,rfin,phipw,phiaug,phitot,  ! for bloch orbital
c     i     mesh,rini,rfin,wanpw,wanaug,wantot,  ! for wannier function 
     i     nbas,bas,nclass,iclass,zz )


      endif 

      write(*,*) 'dump rho? [T/F]'
      read(*,*) lrho
!      if (lrho) call calc_rho(nq_wfn,nband_wfn,npw,phipw,phiaug,phitot)
!      if (lrho) call calc_rho(1,nband_wfn,npw,wanpw,wanaug,wantot)
      if (lrho) call calc_rho_2
     &(alat_ang,nq_wfn,nband_wfn,mesh,rini,rfin,wanpw,wanaug,wantot)
      deallocate(phipw)
      deallocate(phiaug)
      deallocate(phitot)
      deallocate(wanpw)
      deallocate(wanaug)
      deallocate(wantot)
      end program calc_wannier4
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
#if 0
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


      write(*,"(a)") '--- calc_phiall (3)---'
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
#endif
ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_npw(nfac,npw)
      use m_LMTO
      implicit none
c input
      integer :: nfac
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

!      npw(1:3)=2*ntmp(1:3)+2
      npw(1:3)=nfac*ntmp(1:3)+2
      write(*,"(a,3i6)") '# npw(1:3)=',npw(1:3)
      end subroutine calc_npw
ccccccccccccccccccccccccccccccccccccc
#if 0
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
#endif
ccccccccccccccccccccccccccccccccccccc
#if 0
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
#endif
ccccccccccccccccccccccccccccccccccccc
#if 0
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
#endif
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
      subroutine b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,npw,
     &                 phi,wan)
      use m_LMTO
      implicit none
c input
      integer :: nq_wfn,nband_wfn,npw(3),bindx_wfn(nband_wfn),tlat(3)
      double precision :: q_wfn(3,nq_wfn),tvec(3),phase,pi,rtmp(3)
      double complex :: 
     &     phi(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp),
     &     wan(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nsp)
     &    ,ephase
c local
      integer :: iq,isp

c debug:
c      wan(:,:,:,:,1) = phi(:,:,:,:,2,1)
c      return

      pi = 4.0d0*atan(1.d0)
      rtmp(:) = dble(tlat(:))
      call mymatvec(plat,rtmp,tvec,3,3)          
      tvec(1:3)=alat*tvec(1:3)

      wan = (0.0d0,0.0d0)
      do isp = 1,nsp
      do iq = 1,nq_wfn
        phase=2.0d0*pi/alat*sum(q_wfn(1:3,iq)*tvec(1:3))
        ephase=dcmplx(cos(phase),-sin(phase))
        wan(:,:,:,:,isp) = wan(:,:,:,:,isp) + phi(:,:,:,:,iq,isp)*ephase
      enddo ! iq 
      enddo ! isp
      wan = wan / dble(nq_wfn)

      end subroutine b2w
cccccccccccccccccccccccccccccccccccccccc
      subroutine chk_weight()
      use m_LMTO
      use m_MLWF
      implicit none

      double complex :: cphi_s(ldim2,nband)
      real(8):: w
      real(8),allocatable:: wbas(:,:)
      integer(4):: ib,iq,ix,ibas,isp,m,n,l

      allocate(wbas(nbas,nband))

      do isp=1,nsp

      cphi_s = 0.0d0
      do iq = 1,nqbz
         cphi_s(:,:) = cphi_s(:,:) + cphi(:,:,iq,isp)
      enddo ! iq
      cphi_s = cphi_s / dble(nqbz)

      wbas = 0.0d0
      do ib=1,nband
         do ix=1,ldim2
            ibas = mnla(4,ix)
            wbas(ibas,ib) = wbas(ibas,ib) + 
     &                 conjg(cphi_s(ix,ib))*cphi_s(ix,ib)
         enddo ! ix
      enddo ! iwf

      write(*,*)'***     ibas,wbas (summed over band)'
      do ibas=1,nbas
         w = 0d0
         do ib=1,nband
            w = w + wbas(ibas,ib)
         enddo
         write(*,*)ibas,w
      enddo

      do ib=1,nband
         write(*,*)'***  band',ib
         do ibas=1,nbas
            write(*,*)ibas,wbas(ibas,ib)
         enddo ! ibas
      enddo ! ib
   
      enddo ! isp
      deallocate(wbas)

      end subroutine chk_weight
cccccccccccccccccccccccccccccccccccccccc
