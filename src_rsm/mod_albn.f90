      module albn
!------------------------------------------------------------------------------
! for radiation albedo
!------------------------------------------------------------------------------
      use param
      use index
      implicit none
      public
     
      real, dimension (:,:), allocatable, save ::                    &
            alvsf,alvwf,alnsf,alnwf,facsf,facwf

      contains 

         subroutine allocate_alb_array

           integer  ierr

           allocate(alvsf(nxp,my_max),alvwf(nxp,my_max),alnsf(nxp,my_max), &
                    alnwf(nxp,my_max),facsf(nxp,my_max),facwf(nxp,my_max), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_alb : allocate fail 1 '
               stop
           end if

           return

         end subroutine

         subroutine deallocate_alb_array

           deallocate (alvsf,alvwf,alnsf,alnwf,facsf,facwf)                 

           return

         end subroutine

      end module albn

