      program calc_wfn2
c for TiFe2
c chane the boundary for the plot
      use m_LMTO
      implicit none
      integer :: nq_wfn
      double precision,allocatable :: q_wfn(:,:)
      integer :: nband_wfn
      integer,allocatable :: bindx_wfn(:)
      double complex,allocatable :: phipw(:,:,:,:,:,:),
     &     phiaug(:,:,:,:,:,:),phitot(:,:,:,:,:,:)
      double precision :: rini(3),rfin(3),ang,alat_ang
      integer :: npw(3),mesh(3)
      integer :: iq,ib,ibi,ibf
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
!      read(*,*) bindx_wfn(1:nband_wfn)
      read(*,*) ibi
      do ib=1,nband_wfn
        bindx_wfn(ib) = ib + ibi - 1
      enddo
      write(6,"(a,2i5)") '### nq_wfn,nband_wfn=',nq_wfn,nband_wfn

c Setup FFT mesh
      call calc_npw(npw)

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


c Dump phi in OpenDX format
      write(*,*) 'dump phi'
      ang=0.529177d0
      alat_ang=alat*ang
      call wfn2dx_2(alat_ang,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     &     mesh,rini,rfin,phipw,phiaug,phitot)
      call crystal2dx_2(alat_ang,plat,rini,rfin,
     &     nbas,bas,nclass,iclass,zz)

      write(*,*) 'dump rho? [T/F]'
      read(*,*) lrho
      if (lrho) call calc_rho_2
     &(alat_ang,nq_wfn,nband_wfn,mesh,rini,rfin,phipw,phiaug,phitot)
      deallocate(phipw)
      deallocate(phiaug)
      deallocate(phitot)
      end program calc_wfn2
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
