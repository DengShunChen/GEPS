! AMD flang replacements for gfortran libU77 / GNU extensions.
! Only compiled into the ROCm build (see CMakeLists.txt).

subroutine itime(tarray)
  implicit none
  integer, intent(out) :: tarray(3)
  integer :: values(8)
  call date_and_time(values=values)
  tarray(1) = values(5)
  tarray(2) = values(6)
  tarray(3) = values(7)
end subroutine itime

integer function irand(iseed)
  implicit none
  integer, intent(in) :: iseed
  integer, save :: seeded = 0
  integer :: clock, n
  integer, allocatable :: seed_arr(:)
  real :: r
  if (iseed == 0) then
    if (seeded == 0) then
      call system_clock(count=clock)
      call random_seed(size=n)
      allocate(seed_arr(n))
      seed_arr = clock
      call random_seed(put=seed_arr)
      seeded = 1
    end if
  else
    call random_seed(size=n)
    allocate(seed_arr(n))
    seed_arr = iseed
    call random_seed(put=seed_arr)
    seeded = 1
  end if
  call random_number(r)
  irand = int(r * 2147483647.0)
end function irand

subroutine system(command)
  implicit none
  character(len=*), intent(in) :: command
  call execute_command_line(command)
end subroutine system

subroutine flush(unit)
  implicit none
  integer, intent(in) :: unit
  flush (unit)
end subroutine flush
