module mod_vcopy

   interface vcopy_grid
      module procedure vcopy_grid_r4
      module procedure vcopy_grid_2d_r4
      module procedure vcopy_grid_r8
      module procedure vcopy_grid_2d_r8
   end interface vcopy_grid

   interface vcopy_spec
      module procedure vcopy_spec_r4
      module procedure vcopy_spec_2d_r4
      module procedure vcopy_spec_r8
      module procedure vcopy_spec_2d_r8
   end interface vcopy_spec
contains
   subroutine vcopy_spec_2d_r4(A, B)
      use param, only: jtrun, jtmax
      use index, only: mlist, mlistnum
      real(kind=4), intent(out) :: A(jtrun, jtmax, 2)
      real(kind=4), intent(in) :: B(jtrun, jtmax, 2)
      integer k, i, l, mf, m

      do i = 1, 2
         do m = 1, mlistnum
            mf = mlist(m)
            do l = mf, jtrun
               A(l, m, i) = B(l, m, i)
            end do
         end do
      end do

   end subroutine vcopy_spec_2d_r4

   subroutine vcopy_spec_r4(A, B)
      use param, only: jtrun, jtmax
      use index, only: levp, mlist, mlistnum
      real(kind=4), intent(out) :: A(levp, 2, jtrun, jtmax)
      real(kind=4), intent(in) :: B(levp, 2, jtrun, jtmax)
      integer k, i, l, mf, m

      do m = 1, mlistnum
         mf = mlist(m)
         do l = mf, jtrun
            do i = 1, 2
               do k = 1, levp
                  A(k, i, l, m) = B(k, i, l, m)
               end do
            end do
         end do
      end do

   end subroutine vcopy_spec_r4

   subroutine vcopy_spec_2d_r8(A, B)
      use param, only: jtrun, jtmax
      use index, only: mlist, mlistnum
      real(kind=8), intent(out) :: A(jtrun, jtmax, 2)
      real(kind=8), intent(in) :: B(jtrun, jtmax, 2)
      integer k, i, l, mf, m

      do i = 1, 2
         do m = 1, mlistnum
            mf = mlist(m)
            do l = mf, jtrun
               A(l, m, i) = B(l, m, i)
            end do
         end do
      end do

   end subroutine vcopy_spec_2d_r8

   subroutine vcopy_spec_r8(A, B)
      use param, only: jtrun, jtmax
      use index, only: levp, mlist, mlistnum
      real(kind=8), intent(out) :: A(levp, 2, jtrun, jtmax)
      real(kind=8), intent(in) :: B(levp, 2, jtrun, jtmax)
      integer k, i, l, mf, m

      do m = 1, mlistnum
         mf = mlist(m)
         do l = mf, jtrun
            do i = 1, 2
               do k = 1, levp
                  A(k, i, l, m) = B(k, i, l, m)
               end do
            end do
         end do
      end do

   end subroutine vcopy_spec_r8

   subroutine vcopy_grid_2d_r4(A, B)
      use param, only: lev, my_max
      use index, only: nxp, nxdef_2d, jlist1, jlistnum

      real(kind=4), intent(out) :: A(nxp, my_max)
      real(kind=4), intent(in) ::  B(nxp, my_max)
      integer i, j, k, n, jj, nxj

      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do i = 1, nxj
            A(i, jj) = B(i, jj)
         end do
      end do

   end subroutine vcopy_grid_2d_r4

   subroutine vcopy_grid_r4(A, B, ncld)
      use param, only: lev, my_max
      use index, only: nxp, nxdef_2d, jlist1, jlistnum
      real(kind=4), intent(out) :: A(nxp, lev*ncld, my_max)
      real(kind=4), intent(in) ::  B(nxp, lev*ncld, my_max)
      integer, intent(in):: ncld
      integer i, j, k, n, jj, nxj

      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do k = 1, lev*ncld
            do i = 1, nxj
               A(i, k, jj) = B(i, k, jj)
            end do
         end do
      end do

   end subroutine vcopy_grid_r4

   subroutine vcopy_grid_2d_r8(A, B)
      use param, only: lev, my_max
      use index, only: nxp, nxdef_2d, jlist1, jlistnum

      real(kind=8), intent(out) :: A(nxp, my_max)
      real(kind=8), intent(in) ::  B(nxp, my_max)
      integer i, j, k, n, jj, nxj

      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do i = 1, nxj
            A(i, jj) = B(i, jj)
         end do
      end do

   end subroutine vcopy_grid_2d_r8

   subroutine vcopy_grid_r8(A, B, ncld)
      use param, only: lev, my_max
      use index, only: nxp, nxdef_2d, jlist1, jlistnum

      integer, intent(in):: ncld
      real(kind=8), intent(out) :: A(nxp, lev*ncld, my_max)
      real(kind=8), intent(in) ::  B(nxp, lev*ncld, my_max)
      integer i, j, k, n, jj, nxj

      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do k = 1, lev*ncld
            do i = 1, nxj
               A(i, k, jj) = B(i, k, jj)
            end do
         end do
      end do

   end subroutine vcopy_grid_r8
end module mod_vcopy
