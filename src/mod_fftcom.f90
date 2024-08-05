      module fftcom

      implicit none

      public


!
!cc   common/fft/ trigs(512),ifax(19)
!t180 common/fft/ trigs(1024),ifax(19)
!t240 common/fft/ trigs(2048),ifax(19)
!t320
      real    trigs(4096)
      integer ifax(19)
      common/fft/trigs,ifax

! reduceg
      real    trigsj(4096,2560)
      integer ifaxj(19,2560)
      common /reducefftm/trigsj,ifaxj
!
      integer, parameter :: ibm_fft = 0 
!
!  working arrays for ibm_fft
!
      logical lessl_fft
!
!  the parameter below is valid only for original rfftmlt used
!
!  length_fft : an option for doing fft with long or short vector
!  length_fft = 0 : use long vector, which favors vector machines like
!                   vpp5000
!               1 : use short vector, which has better performance
!                   on scalar machines
!
      integer, parameter :: length_fft = 1
!xxx 

      save
      end module fftcom
