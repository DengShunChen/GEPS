      subroutine days (cdtg,julian,hours)

      implicit  none
      integer   iy,im,id,ih,imin,julian,i
      integer   mon(12)
      real      hours
      character cdtg*12
      data mon /0,31,28,31,30,31,30,31,31,30,31,30/
!
      read(cdtg,'(i4,4i2)') iy,im,id,ih,imin
!cc   iy = 1900 + iy
      if (mod(iy,4).eq.0)  mon(3) = 29
      julian = 0
      do 100 i=1,im
      julian = julian + mon(i)
 100  continue
      julian = julian + id
      hours = float(ih)
!xx   hours = float(ih) + 8.0
      return
      end
