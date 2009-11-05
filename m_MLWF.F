      module m_MLWF
      use m_DATA4GW
      use m_QG
      implicit none
      integer :: nqbz,nwf,iko_ix,iko_fx
      double precision :: a
      double complex,allocatable :: dnk(:,:,:,:)

      contains
ccccccccccccccccccccccccccccccccccccccc
      subroutine setup_MLWF()
      implicit none
      integer :: isp,iqbz,ikp,ib,iwf
      double complex,allocatable :: geig2(:,:,:,:),cphi2(:,:,:,:)
      write(6,*) '--- setup_MLWF ---'

      call read_MLWF()

c replace geig and cphi with those for WF
c geig2 = cphi2 = 0 for ik > nqbz
      allocate(geig2(ngpmx,nband,nqtt,nsp))
      allocate(cphi2(ldim2,nband,nqtt,nsp))
      geig2 = geig
      cphi2 = cphi
      deallocate(geig,cphi)
      allocate(geig(ngpmx,nwf,nqtt,nsp))
      allocate(cphi(ldim2,nwf,nqtt,nsp))
      geig = (0.0d0,0.0d0)
      cphi = (0.0d0,0.0d0)
      do ikp = 1,nqbz
      do isp = 1,nsp
         do iwf = 1,nwf
         do ib = iko_ix,iko_fx
            geig(:,iwf,ikp,isp) = geig(:,iwf,ikp,isp) + 
     &           geig2(:,ib,ikp,isp)*dnk(ib,iwf,ikp,isp)
            cphi(:,iwf,ikp,isp) = cphi(:,iwf,ikp,isp) + 
     &           cphi2(:,ib,ikp,isp)*dnk(ib,iwf,ikp,isp)
         enddo ! ib
         enddo ! iwf
      enddo ! isp
      enddo ! ikp

      deallocate(geig2,cphi2,dnk)
      nband = nwf

      end subroutine setup_MLWF
ccccccccccccccccccccccccccccccccccccccc
      subroutine read_MLWF()
      implicit none
      double precision :: q(3),eps
      parameter (eps=1d-4)
      integer :: ifi
      integer :: isp,iqbz
      integer :: iqbz2

      write(6,*) '--- read_MLWF ---'

      do isp=1,nsp
c file open
         ifi = 1000
         if (isp.eq.1) then
           open(ifi,file='MLWU',form='unformatted',status='old',
     &          action='read')
         else
           open(ifi,file='MLWD',form='unformatted',status='old',
     &          action='read')
         endif
c nqbz mesh-points
         read(ifi)nqbz,nwf,iko_ix,iko_fx
         if (isp.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp))
         do iqbz = 1,nqbz
            read(ifi)iqbz2,q(1:3)
            q(:) = q(:) - qtt(:,iqbz)
            if (sum(q(:)**2).gt.eps) stop 'MLWU/D: qbz error'
            if (iqbz2.ne.iqbz) stop 'MLWU/D: iqbz error'
            read(ifi)dnk(iko_ix:iko_fx,1:nwf,iqbz,isp)
         enddo

c fileclose
         close(ifi)
      enddo ! isp

      end subroutine read_MLWF
cccccccccccccccccccccccccccccccccccccccc
      end module m_MLWF
