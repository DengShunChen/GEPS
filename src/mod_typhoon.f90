      module mod_typhoon
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use param
      use index,  only: nxp
      use const,  only: RTYPE

      implicit none

      public

      integer  ntyph,ixtyp(5,5),jytyp(5,5),nrec(5)
      common/hurricanI/ntyph,ixtyp,jytyp,nrec

      real clattyp(5),clontyp(5)

!     arrays(cases,fields,typhoons)
      real, dimension(0:168,5,5) ::  tflon,tflat,tensity 
 
!byl      real, dimension(:,:), allocatable, save :: tlon, tydom,  &
!byl                                slp, v850, v700, h850, h500
      real, dimension(:,:), allocatable, save :: tydom
      real(kind=RTYPE), dimension(:,:,:), allocatable,save ::  typtrk

      real, dimension(:),   allocatable, save :: tlon, tlat

      character*15      typname(5)
      common/hurricanC/typname

      logical          nlexist,typhoon,nWlexist,olexist,oWlexist
      common/hurricanL/nlexist,typhoon,nWlexist,olexist,oWlexist

      character typhnam*15,typhpath*120,Wtyphpath*120
      character ntyphfile*24,nWtyphfile*24,otyphfile*24,oWtyphfile*24
      
      common/hurricanC1/typhnam,typhpath,Wtyphpath                  &
                       ,ntyphfile,nWtyphfile,otyphfile,oWtyphfile

     !  tracker
      real :: tau_trk, dt_trk, min_trk_pres=1002.
      logical :: ltrack=.false.

      integer :: write_tau=6, trk_intv=6
      integer :: write_mem=00

      contains 

         subroutine allocate_typhoon_array

           integer  ierr

!byl           allocate (tlon(nx,my),tlat(my),tydom(nx,my),             &
!byl                     slp(nx,my),v850(nx,my),v700(nx,my),h850(nx,my),&
!by                     h500(nx,my), stat=ierr)
           allocate (tlon(nx),tlat(my),tydom(nx,my),             &
                     typtrk(nxp,my_max,5), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_typhoon : allocate fail 1 '
               stop
           end if
!
           return

         end subroutine

         subroutine deallocate_typhoon_array

!byl           deallocate (tlon, tydom, slp, v850, v700, h850, h500, tlat)
           deallocate (tlon, tydom, typtrk)

           return

         end subroutine

      end module mod_typhoon
