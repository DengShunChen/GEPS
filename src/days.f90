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
!
      subroutine leapyear(iy)
      use leapyr
      use rank
      implicit none
      integer iy,iym1
      leap=.false.
      leapm1=.false.
      iym1=iy-1
      if(mod(int(iym1),4).eq.0) then
        leapm1=.true.
        if(mod(int(iym1),100).eq.0) then
          leapm1=.false.
          if(mod(int(iym1),400).eq.0) leapm1=.true.
        endif
      endif
      if(mod(int(iy),4).eq.0) then
        leap=.true.
        if(mod(int(iy),100).eq.0) then
          leap=.false.
          if(mod(int(iy),400).eq.0) leap=.true.
        endif
      endif
      if(leap) then
!        yrjulian=366.d0
         yrd=366
      else
!        yrjulian=365.d0
         yrd=365
      endif
      return
      end
!
      subroutine datecheck(iy,julian,ih)
      use radn,   only: idate
      use leapyr
      use rank
      implicit none
      integer nmon
      parameter (nmon=12)
      integer julian,iday
      integer iy,im,id,ih
      integer nmod(nmon)
      data nmod /31,28,31,30,31,30,31,31,30,31,30,31/
      call leapyear(iy)
      if(leap) then
        nmod(2)=29
      else
        nmod(2)=28
      endif
      id=0
      iday=0
      im=1
  100 id=id+1
      iday=iday+1
      if(id .gt. nmod(im)) then
        id=1
        im=im+1
        if(im .gt. nmon) then
          im=1
          iy=iy+1
          call leapyear(iy)
          if(leap) then
            nmod(2)=29
          else
            nmod(2)=28
          endif
          if(myrank.eq.0) then
            print*, 'leap year or not for',iy,' is',leap
            print*, 'number of days for',iy,' are',yrd
          endif
        endif
      endif
      if(iday.lt.julian) go to 100
!     idate=ih*100+id*10000+im*1000000+iy*100000000
      idate(1)=iy
      idate(2)=im
      idate(3)=id
      idate(5)=ih
      return
      end
