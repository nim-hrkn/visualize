      module m_FFT3D
      implicit none
#include <fftw3.f>      
      integer :: nfft(3)
      integer*8 :: fftw_plan
      double complex,allocatable :: in_fft(:,:,:),out_fft(:,:,:)
      logical :: linit=.false.
      contains
cccccc
      subroutine fft_init(nin,sign)
      implicit none
      integer :: nin(3)
      character :: sign
      linit=.true.
      nfft(1:3)=nin(1:3)
      write(*,*) '--- fft_init ---'
      allocate(in_fft(nfft(1),nfft(2),nfft(3)),
     &     out_fft(nfft(1),nfft(2),nfft(3)))

      if (sign.eq.'B') then
        call dfftw_plan_dft_3d(fftw_plan,nfft(1),nfft(2),nfft(3),
     &       in_fft,out_fft,FFTW_BACKWARD,FFTW_MEASURE)
      elseif (sign.eq.'F') then
        call dfftw_plan_dft_3d(fftw_plan,nfft(1),nfft(2),nfft(3),
     &       in_fft,out_fft,FFTW_BACKWARD,FFTW_MEASURE)
      else
        stop 'myfft:fft_init unknown sign'
      endif
      end subroutine fft_init
cccccc
      subroutine myfft()
      implicit none
      if (.not.linit) stop 'myfft: not initialized yet'
      call dfftw_execute(fftw_plan)
      end subroutine myfft
cccccc
      end module m_FFT3D
