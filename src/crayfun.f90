      function ismax(n,x,incx)
      dimension x(n)
      mmm = 1
      vmax = x(1)
      do i = 2, n, incx
      if( x(i) .gt. vmax ) then
       vmax = x(i)
       mmm = i
      end if
      end do
      ismax = mmm
      return
      end

      function ismin(n,x,incx)
      dimension x(n)
      mmm = 1
      vmin = x(1)
      do i = 2, n, incx
      if( x(i) .lt. vmin ) then
       vmin = x(i)
       mmm = i
      end if
      end do
      ismin = mmm
      return
      end

      function ilsum(n,x,incx)
      logical x(n)
      num = 0
      do i = 1, n, incx
      if( x(i) ) then
       num = num + 1
      end if
      end do
      ilsum = num
      return
      end

      function isamax(n,x,incx)
      dimension x(n)
      mmm = 1
      vmax = abs(x(1))
      do i = 2, n, incx
      if( abs(x(i)) .gt. vmax ) then
       vmax = abs(x(i))
       mmm = i
      end if
      end do
      isamax = mmm
      return
      end

! Cray rtc() / NCEP timef() — not provided by gfortran
      function rtc()
      real*8 rtc
      integer count, rate
      call system_clock(count, rate)
      if (rate .gt. 0) then
        rtc = dble(count) / dble(rate)
      else
        rtc = 0.d0
      endif
      return
      end

      function timef()
! return wall time in milliseconds (NCEP/W3 convention)
      real timef
      integer count, rate
      call system_clock(count, rate)
      if (rate .gt. 0) then
        timef = real(count) * 1.0e3 / real(rate)
      else
        timef = 0.0
      endif
      return
      end
