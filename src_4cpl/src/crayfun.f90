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
