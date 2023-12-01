! retype in
! based on CAM

MODULE shr_kind_mod

   public
   integer,parameter :: SHR_KIND_R8 = selected_real_kind(12)  !8 byte real
!  integer,parameter :: SHR_KIND_R8 = selected_real_kind(6)
   integer,parameter :: SHR_KIND_R4 = selected_real_kind( 6)  !4 byte real 

END MODULE shr_kind_mod
