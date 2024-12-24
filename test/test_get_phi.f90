program test_get_phi
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call get_phi_unit
   call mpe_finalize

end program
subroutine get_phi_unit
   use param, only: lev, my_max, ncld
   use const, only: RTYPE, ptop, cp, rgas, grav
   use index, only: nxjp, nxp, jlistnum, jlist1, nxdef_2d
   use grid, only: sgeo

   implicit none

   real(kind=RTYPE), dimension(nxp, lev, my_max) :: tt, pk, pk2
   real(kind=RTYPE), dimension(nxp, lev*ncld, my_max) :: qt
   real, dimension(nxp, lev + 1, my_max), target :: phii_g, phii_c
   real(kind=RTYPE), dimension(nxp, lev, my_max), target :: phi_g, phi_c
   real, dimension(:, :, :), pointer:: phii
   real(kind=RTYPE), dimension(:, :, :), pointer:: phi

   integer, parameter :: steps = 5
   integer :: i, jj, j, nxj, async_id, seed_size
   integer, allocatable :: seed(:)
   real(kind=RTYPE) abs_diff, rel_diff

   async_id = 1

   call random_seed(size=seed_size)
   allocate (seed(seed_size))
   seed = 123
   call random_seed(put=seed)
   call random_number(tt)
   call random_number(qt)
   call random_number(pk)
   call random_number(pk2)

   phii_c = 0.
   phi_c = 0.
   phii_g = 0.
   phi_g = 0.

   ! << CPU >>
   phii => phii_c
   phi => phi_c
   do i = 1, steps
      do jj = 1, jlistnum
         j = jlist1(jj)
         call get_phi(nxjp(j), nxp, lev, ptop, cp, rgas, grav, sgeo(1, jj), &
                      pk(1, 1, jj), pk2(1, 1, jj), tt(1, 1, jj), qt(1, 1, jj), &
                      phii(1, 1, jj), phi(1, 1, jj))
      end do
   end do

   ! << GPU >>
   phii => phii_g
   phi => phi_g
   !$acc enter data copyin(nxjp, nxp, jlistnum, jlist1, lev, &
   !$acc&                  ptop, cp, rgas, grav, sgeo) async(async_id)
   !$acc enter data copyin(pk, pk2, tt, qt) async(async_id)
   !$acc enter data copyin(phii, phi) async(async_id)
   !$acc wait(async_id)
   do i = 1, steps
      call get_phi_gpu(nxjp, nxp, lev, ptop, cp, rgas, grav, sgeo, &
                       pk, pk2, tt, qt, &
                       phii, phi)
   end do
   !$acc exit data delete(nxjp, nxp, jlistnum, jlist1, lev, &
   !$acc&                 ptop, cp, rgas, grav, sgeo) async(async_id)
   !$acc exit data delete(pk, pk2, tt, qt) async(async_id)
   !$acc exit data copyout(phii, phi) async(async_id)
   !$acc wait(async_id)

#ifdef SP
#else
   call assert_allclose(phii_g, size(phii_g), phii_c, size(phii_c), 1e-8, 1e-8, "Array phii")
   call assert_allclose(phi_g, size(phi_g), phi_c, size(phi_c), 1e-8, 1e-8, "Array phi")
#endif

end subroutine
