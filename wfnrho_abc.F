      module m_wfrho_abc 
      contains
      subroutine calc_phiall_abc(nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     &     npw,mesh,phipw,phiaug,phitot)
      use m_LMTO
      use m_FFT3D
      implicit none

c inputs
      integer :: nq_wfn,nband_wfn,bindx_wfn(nband_wfn)
      double precision :: q_wfn(3,nq_wfn)
      integer :: npw(3),mesh(3)
c outputs
      double complex :: 
     &     phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

c local
      integer :: isp,iq,iq_wfn,ib,i1,i2,i3
      integer :: augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)

      double complex :: ! temp. array
     &     phipwtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1),
     &     phiaugtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1)

      double complex :: ! e^{ikr} and e^{ikT} (T is lattice vector)
     &     eikr(mesh(1)+1,mesh(2)+1,mesh(3)+1),
     &     eikT(mesh(1)+1,mesh(2)+1,mesh(3)+1)


      write(*,"(a)") '--- calc_phiall ---'
!      call fft_init(npw,'B')
      call calc_augregion_abc(mesh(1),mesh(2),mesh(3),augregion)

!$omp parallel do private(iq, eikr,eikT, phipwtmp,phiaugtmp )
      do iq_wfn=1,nq_wfn
        call findq(q_wfn(1:3,iq_wfn),iq)   
        call calc_eikreikT_abc
     &      (q_wfn(1:3,iq_wfn),mesh,augregion,eikr,eikT)
        do isp=1,nsp
        do ib=1,nband_wfn
!$omp critical
          write(*,"(a,i2,2i5,3f10.4,i5)")
     &         '# isp,iq_wfn,iq,q,ib=',isp,iq_wfn,iq,qtt(1:3,iq),ib
!$omp end critical
          call calc_phi_abc(isp,iq,bindx_wfn(ib),
     &         npw,mesh,
     &         augregion,phipwtmp,phiaugtmp)

          do i3=1,mesh(3)+1
          do i2=1,mesh(2)+1
          do i1=1,mesh(1)+1
c  bloch function
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
c          write(*,'(6f10.5)') phitot(:,:,:,ib,iq_wfn,isp)
        enddo !ib
        enddo !isp
      enddo !iq
      end subroutine calc_phiall_abc
ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_augregion_abc(n1,n2,n3,augregion)
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
!            call mymatvec(plat,rtmp,r,3,3)
c            r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
            r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
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
      end subroutine calc_augregion_abc
ccccccccccccccccccccccccccccccccccccc
      subroutine calc_eikreikT_abc(kvec,mesh,
     &                           augregion,eikr,eikT)
      use m_LMTO
      implicit none
c input
      double precision :: kvec(3)
      integer :: mesh(3),augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
c output
      double complex :: ! e^{ikr} and e^{ikT} (T is lattice vector)
     &     eikr(mesh(1)+1,mesh(2)+1,mesh(3)+1),
     &     eikT(mesh(1)+1,mesh(2)+1,mesh(3)+1)
c local
      integer :: i1,i2,i3
      double precision :: rtmp(3),r(3),tvec(3)
      double precision :: phase,pi

      pi=4.0d0*atan(1.0d0)

      write(*,*) 'kvec=',kvec (1:3) 
c Calculate e^{ikr}
      do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
      do i1=1,mesh(1)+1
        rtmp(1)=(i1-1)/dble(mesh(1))
        rtmp(2)=(i2-1)/dble(mesh(2))
        rtmp(3)=(i3-1)/dble(mesh(3))
!        call mymatvec(plat,rtmp,r,3,3)
c        r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
        r(:)= plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
        r(1:3)=alat*r(1:3)
        phase=2.0d0*pi/alat*sum(kvec(1:3)*r(1:3))
        eikr(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
      enddo
      enddo      
      enddo

c Calculate e^{ikT}
      do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
      do i1=1,mesh(1)+1

        if (augregion(4,i1,i2,i3).ne.0) then
          rtmp(1:3)=augregion(1:3,i1,i2,i3)
c             tvec(i) =plat(i,j)*rtmp(j)
          call mymatvec(plat,rtmp,tvec,3,3)          
          tvec(1:3)=alat*tvec(1:3)
c  2 pi  k(i)*tvec(i) 
          phase=2.0d0*pi/alat*sum(kvec(1:3)*tvec(1:3))
          eikT(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
        else
          eikT(i1,i2,i3)=0.0d0
        endif
      enddo
      enddo      
      enddo

      end subroutine calc_eikreikT_abc
ccccccccccccccccccccccccccccccccccccc
      subroutine calc_phi_abc(isp,iq,iband,npw,mesh,
     &     augregion,phipwtmp,phiaugtmp)
      use m_LMTO
      use m_FFT3D
      implicit none
c input
      integer :: isp,iq,iband,npw(3),mesh(3)
      integer :: augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
      double precision :: qlat(3,3),tpi
      double complex :: eigr,ci

c output
      double complex :: 
     &     phipwtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1),
     &     phiaugtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1)

c local
      integer :: itmp(3),ig,id,i1,i2,i3,j1,j2,j3,ii

      double precision :: rtmp(3),r(3),r0(3) !points to plot
      double precision :: ratom(3) ! atomic points
      double precision :: dr(3)

      ci = (0.0d0,1.0d0)
      tpi = 8d0*atan(1.0d0)
      call dinv33x(plat,qlat)
      call chkinv33(plat,qlat)

c      write(*,*) '--- calc_phi ---'

c Plane-wave part
!      in_fft(:,:,:)=0.0d0
!      do ig=1,ngp(iq)
!        do id=1,3
!          itmp(id)=ngvecp(id,ig,iq)
!          if (itmp(id).lt.0) itmp(id)=itmp(id)+npw(id)
!        enddo
!        in_fft(itmp(1)+1,itmp(2)+1,itmp(3)+1)=geig(ig,iband,iq,isp)
!      enddo
c Perform backward fft
!      call myfft()
!
!      do i3=1,npw(3)+1
!      do i2=1,npw(2)+1
!      do i1=1,npw(1)+1
!        phipwtmp(i1,i2,i3)=out_fft(mod(i1-1,npw(1))+1,
!     &       mod(i2-1,npw(2))+1,mod(i3-1,npw(3))+1)
!      enddo
!      enddo
!      enddo
!
      phipwtmp = 0d0
      do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
      do i1=1,mesh(1)+1
         rtmp(1)=(i1-1)/dble(mesh(1))
         rtmp(2)=(i2-1)/dble(mesh(2))
         rtmp(3)=(i3-1)/dble(mesh(3))
c         r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
          r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
!         r0(:) = matmul(qlat,r)
         do ii=1,3
           r0(ii) = sum(qlat(:,ii)*r(:))
         enddo ! ii
c   r0(i)=G0(j,i)*r(j)*
c   G(i)= G0(j,i)*nG(i) 
c   exp (i 2 pi  G(i)*r(i) ) 
         do ig=1,ngp(iq)
           eigr=exp(ci*tpi*sum(r0(:)*dble(ngvecp(:,ig,iq))))
           phipwtmp(i1,i2,i3) = phipwtmp(i1,i2,i3)
     &        + eigr*geig(ig,iband,iq,isp)
!        phipwtmp(i1,i2,i3)=out_fft(mod(i1-1,npw(1))+1,
!     &       mod(i2-1,npw(2))+1,mod(i3-1,npw(3))+1)
         enddo ! ig
      enddo ! i1
      enddo ! i2
      enddo ! i3

c Augmented part
      phiaugtmp(:,:,:)=0.0d0
      do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
      do i1=1,mesh(1)+1

        if (augregion(4,i1,i2,i3).ne.0) then
c set plane-wave part to zero
          phipwtmp(i1,i2,i3)=0.0d0
          rtmp(1)=(i1-1)/dble(mesh(1))
          rtmp(2)=(i2-1)/dble(mesh(2))
          rtmp(3)=(i3-1)/dble(mesh(3))
!          call mymatvec(plat,rtmp,r,3,3)
!          r(1:3)=alat*r(1:3)
c          r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
          r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
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
      end subroutine calc_phi_abc
cccccccccccccccccccccccccccccccccccc
      subroutine calc_rho_abc(alat_ang,nq_wfn,nband_wfn,mesh,
     &   phipw,phiaug,phitot)
      use m_LMTO
      implicit none
c input
      integer :: nq_wfn,nband_wfn,mesh(3)
      double precision :: alat_ang
      double complex :: 
     &     phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

      double complex,allocatable :: rhopw(:,:,:),
     &     rhoaug(:,:,:),rhotot(:,:,:)

c
      integer :: isp,iq,ib,i1,i2,i3
      double precision :: nel

      double precision :: qdum(3)
      integer :: bindxdum
      write(*,*) '--- calc_rho ---'
c Allocate rho
      allocate(rhopw(mesh(1)+1,mesh(2)+1,mesh(3)+1))
      allocate(rhoaug(mesh(1)+1,mesh(2)+1,mesh(3)+1))
      allocate(rhotot(mesh(1)+1,mesh(2)+1,mesh(3)+1))

      rhopw(1:mesh(1)+1,1:mesh(2)+1,1:mesh(3)+1)=0.0d0
      rhoaug(1:mesh(1)+1,1:mesh(2)+1,1:mesh(3)+1)=0.0d0
      rhotot(1:mesh(1)+1,1:mesh(2)+1,1:mesh(3)+1)=0.0d0

      do isp=1,nsp
      do iq=1,nq_wfn
      do ib=1,nband_wfn
        do i3=1,mesh(3)+1
        do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1
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
      do i3=1,mesh(3)
      do i2=1,mesh(2)
      do i1=1,mesh(1)
        nel=nel+rhotot(i1,i2,i3)
      enddo
      enddo
      enddo
      nel=nel*dble(3-nsp)*VOL/dble(mesh(1)*mesh(2)*mesh(3))
      write(*,*) 'nel = ',nel

      qdum(1:3)=0.0d0
      bindxdum=0
      call wfn2dx_abc(alat_ang,plat,1,1,1,qdum,bindxdum,
     &     mesh,rhopw,rhoaug,rhotot)

      deallocate(rhopw)
      deallocate(rhoaug)
      deallocate(rhotot)
      end subroutine calc_rho_abc
cccccccccccccccccccccccccccccc
      end module m_wfrho_abc 
