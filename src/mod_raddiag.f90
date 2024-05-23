      module raddiag
!------------------------------------------------------------------------------
! for radiation diagnose
!------------------------------------------------------------------------------
      use const, only : RTYPE
      use param
      use index
      implicit none
      public
     
      real, dimension(:,:), allocatable, save ::                       &
            asol    , olr    , ss_clr , rs_clr , olr_clr ,             &
            asold   ,                                                  &
            asol_clr, sld_clr, rld_clr 

      real(kind=RTYPE),dimension(:,:,:),allocatable,save ::            &
            fusl   , fdsl   , fuir   , fdir   ,                        &
            fuslr  , fdslr  , fuirr  , fdirr  

      real, dimension(:,:,:), allocatable, save ::                     &
            asl_clr, atl_clr, clds 

      contains 

         subroutine allocate_raddiag_array

           integer  ierr

           allocate (   asol(nxp,my_max),     olr(nxp,my_max), &
                      ss_clr(nxp,my_max),  rs_clr(nxp,my_max), &
                     olr_clr(nxp,my_max),asol_clr(nxp,my_max), &
                       asold(nxp,my_max),                      &
                     sld_clr(nxp,my_max), rld_clr(nxp,my_max), stat=ierr)
        
           if (ierr/= 0) then
               write(6,*) 'mod_raddiag : allocate fail 1 '
               stop
           end if

           allocate ( fusl(nxp,lev+1,my_max)  , fdsl(nxp,lev+1,my_max)  ,     &
                      fuir(nxp,lev+1,my_max)  , fdir(nxp,lev+1,my_max)  ,     &
                      fuslr(nxp,lev+1,my_max) , fdslr(nxp,lev+1,my_max) ,     &
                      fuirr(nxp,lev+1,my_max) , fdirr(nxp,lev+1,my_max) ,     & 
                      asl_clr(nxp,lev,my_max) , atl_clr(nxp,lev,my_max) ,     &
                      clds(nxp,lev,my_max),stat=ierr  )

                     
           if (ierr/= 0) then
               write(6,*) 'mod_raddiag : allocate fail 2 '
               stop
           end if
!      
           asol=0.
           olr=0.
           asold=0.
           ss_clr=0.
           rs_clr=0.
           olr_clr=0.
           asol_clr=0.
           sld_clr=0.
           rld_clr=0.
!      
           fusl=0.
           fdsl=0.
           fuir=0.
           fdir=0.
           fuslr=0.
           fdslr=0.
           fuirr=0.
           fdirr=0.
           asl_clr=0.
           atl_clr=0.
           clds=0.

           return

         end subroutine

         subroutine deallocate_raddiag_array

           deallocate (asol,olr,ss_clr,rs_clr,olr_clr,asol_clr,sld_clr, &
                       rld_clr)
           deallocate (asold)
           deallocate (fusl,fdsl,fuir,fdir,fuslr,fdslr,fuirr,fdirr,     &
                       asl_clr,atl_clr,clds)


           return

         end subroutine

       end module raddiag
