      subroutine ancrt (m,pl,ancrit)
!
!  ancrt computes the climotological work function as a function of
!  pressure by linear interpolation between values defined at
!  predefined pressures.
!
!  formal parameters:
!
!  ancrit:  output work functions values at pressures pl
!
!  pl:  input pressures
!
!  m: no. of points
!
      implicit none

      integer m,i,ipx
      real    ancrit(m),pl(m),acrit(11),tem,ptem,pxmax

      data ptem/0.01/, pxmax/10.999999/
!
      data acrit/2.0,1.8983,1.2425,.5162,.3252,.1915,.0924,.0577, &
        .0350,.0220,.0150/
!
!  interpolate reference cloud work function to pressure levels
!
      do 110 i=1,m
      ancrit(i)= min(pxmax,1.0+ptem*pl(i))
      ipx= int(ancrit(i))
      tem= ancrit(i)-float(ipx)
      ancrit(i)= acrit(ipx)+tem*(acrit(ipx+1)-acrit(ipx))
  110 continue
!
      return
      end
