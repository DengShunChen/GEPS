      module paramt

      implicit none

      public

!     integer, parameter :: im=1536
!     integer, parameter :: jm=im/2
!     integer, parameter :: lm=60
!     integer, parameter :: jtr=2*((1+(im-1)/3)/2)
!     integer, parameter :: mlm=jtr*(jtr+1)/2

!     integer, parameter :: imax=1536
!     integer, parameter :: jmax=imax/2
!     integer, parameter :: lpx=26
!     integer, parameter :: kvkw=im/2

!     integer, parameter :: ilm=im*lm
!     integer, parameter :: im2=kvkw*(lm+1)

!     integer, parameter :: jmhalf=jm/2
!     integer, parameter :: lmX2=lm*2
!     integer, parameter :: lmX4=lm*4
!     integer, parameter :: lmX10=lm*10
!
      integer  im,jm,lm,jtr,mlm,imax,jmax,lpx,kvkw,ilm,im2, &
               jmhalf,lmX2,lmX4,lmX10

      common   /comparamt/ &
               im,jm,lm,jtr,mlm,imax,jmax,lpx,kvkw,ilm,im2, &
               jmhalf,lmX2,lmX4,lmX10

      end module paramt
