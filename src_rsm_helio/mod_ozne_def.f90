      module ozne_def
      use machine , only : kind_phys
      implicit none
      save
      integer, parameter :: kozpl=28, kozc=48
      integer latsozp, levozp, timeoz, latsozc, levozc, timeozc
      integer pl_coeff
      real (kind=kind_phys) blatc, dphiozc
      real (kind=kind_phys), allocatable :: PL_LAT(:), PL_Pres(:)
      real (kind=kind_phys), allocatable :: PL_LATC(:),PL_PresC(:)
      real (kind=kind_phys), allocatable :: PL_TIME(:)
      real (kind=kind_phys), allocatable :: ozplin(:,:,:,:)
      real (kind=kind_phys), allocatable :: ozclm(:,:,:)
      end module ozne_def
