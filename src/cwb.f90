      subroutine sendmsg(str_model,i_fromtau,i_totau,i_status)

      character*3 str_model
      character tauchk*10

      i_status=-1

      tauchk='zzz.tauchk'

      write(tauchk(1:3),'(a3)')str_model

      open(1,file=tauchk,form='formatted',status='unknown', &
           position='append')
      write(1,'(i4.4,x,i4.4)')i_fromtau,i_totau
      close(1)

      i_status=0

      return
      end

!----------------------------------------------------------------------
      subroutine recmsg(model,ifrom,ito,istatus)
      character*3 model
      
      integer tau(512),cnt,ip
      save    tau,cnt,ip

      istatus=0
      if(ip.eq.0)then
        open(1,file=model//'ctl',form='formatted',status='old')
        do i=1,512
           read(1,'(I4)',end=100)tau(i)
           cnt=cnt+1
        enddo
100     close(1)
        ip=1
      endif
      if(ip.ge.cnt) then
        istatus=1
        return
      endif
      ifrom=tau(ip)
      ip=ip+1
      ito=tau(ip)
      return
      end
 
      subroutine dtgfix12(idtg1,idtg2,chg)

      integer*8 idtg1,idtg2
      integer   chg,yyyy,mm,dd,hh,ii,dt,diff,rem
      integer   month(12)
      character*12 cdtg
      data      month/31,28,31,30,31,30,31,31,30,31,30,31/

      write(cdtg,'(i12)')idtg1
      read(cdtg,'(i4,i2,i2,i2,i2)')yyyy,mm,dd,hh,ii

      if(((mod(yyyy,4) .eq.0) .and. (mod(yyyy,100) .ne.0)) .or.  &
         ((mod(yyyy,100) .eq. 0) .and. (mod(yyyy,400) .eq. 0))) then
         month(2)=29
      endif

      dt=hh+chg  
      if(dt .lt. 0) then
	  diff=chg*(-1)/24
          rem=mod(chg*(-1),24)
          hh=hh-rem
	  if (hh .lt. 0)then
             diff=diff+1
	     hh=hh+24
          endif
          dd=dd-diff
          if(dd .lt. 1)then
             mm=mm-1
	     if(mm .gt. 0) then
               dd=dd+month(mm)
             else
               yyyy=yyyy-1
               mm=12
               dd=31
             endif
          endif
       elseif (dt .gt. 23)then
          diff=chg/24
          rem=mod(chg,24)
          hh=hh+rem
	  if(hh .gt. 23)then
	     hh=hh-24
             diff=diff+1
          endif
          dd=dd+diff
          do while (dd.gt.month(mm))
!          if(dd .gt. month(mm))then
            dd=dd-month(mm)
            if(mm .eq. 12)then
	      yyyy=yyyy+1
              mm=1
            else
              mm=mm+1
            endif
            if(((mod(yyyy,4) .eq.0) .and. (mod(yyyy,100) .ne.0)) .or.  &
              ((mod(yyyy,100) .eq. 0) .and. (mod(yyyy,400) .eq. 0))) then
              month(2)=29
            endif
!          endif
          enddo
       else 
         hh=dt 
       endif
!      idtg2=yyyy*100000000+mm*1000000+dd*10000+hh*100+ii
       write(cdtg,'(i4,i2.2,i2.2,i2.2,i2.2)')yyyy,mm,dd,hh,ii
       read(cdtg,'(i12)')idtg2
       return

100    print *,'dtgfix12 error, idtg1=',idtg1
       return
       end
!
!ch 2019-05-10 for prerrtmg potential problem found by Pang-Yen Liu (b827)
      subroutine dtgfix12_new(idtg1,i1,i2,i3,i4,i5,chg)

      implicit none
      integer*8 idtg1
      integer   i1,i2,i3,i4,i5
      integer   chg,yyyy,mm,dd,hh,ii,dt,diff,rem
      integer   month(12)
      character*12 cdtg
      data      month/31,28,31,30,31,30,31,31,30,31,30,31/

      write(cdtg,'(i12)')idtg1
      read(cdtg,'(i4,i2,i2,i2,i2)')yyyy,mm,dd,hh,ii

      if(((mod(yyyy,4) .eq.0) .and. (mod(yyyy,100) .ne.0)) .or.  &
         ((mod(yyyy,100) .eq. 0) .and. (mod(yyyy,400) .eq. 0))) then
         month(2)=29
      endif

      dt=hh+chg
      if(dt .lt. 0) then
          diff=chg*(-1)/24
          rem=mod(chg*(-1),24)
          hh=hh-rem
          if (hh .lt. 0)then
             diff=diff+1
             hh=hh+24
          endif
          dd=dd-diff
          if(dd .lt. 1)then
             mm=mm-1
             if(mm .gt. 0) then
               dd=dd+month(mm)
             else
               yyyy=yyyy-1
               mm=12
               dd=31
             endif
          endif
       elseif (dt .gt. 23)then
          diff=chg/24
          rem=mod(chg,24)
          hh=hh+rem
          if(hh .gt. 23)then
             hh=hh-24
             diff=diff+1
          endif
          dd=dd+diff
          do while (dd.gt.month(mm))
!          if(dd .gt. month(mm))then
            dd=dd-month(mm)
            if(mm .eq. 12)then
              yyyy=yyyy+1
              mm=1
            else
              mm=mm+1
            endif
            if(((mod(yyyy,4) .eq.0) .and. (mod(yyyy,100) .ne.0)) .or.  &
              ((mod(yyyy,100) .eq. 0) .and. (mod(yyyy,400) .eq. 0))) then
              month(2)=29
            endif
!          endif
          enddo
       else
         hh=dt
       endif

       i1=yyyy
       i2=mm
       i3=dd
       i4=hh
       i5=ii

       return

100    print *,'dtgfix12 error, idtg1=',idtg1
       return
       end


      subroutine GETFNAME(pathname,logicfile,truefile,istat)
      character pathname*80,logicfile*80,truefile*80,value*80,blank*80
!     integer getenv

      truefile(1:80)=' '
      value(1:80)=' '
      blank(1:80)=' '
!     i = getenv(pathname,value)
!     if (i .eq. 0) then
      ix=index(pathname,' ')
      call getenv(pathname(1:ix-1),value)
      if (value(1:80) .eq. blank(1:80)) then
 	istat=-1
      else
	j=index(value,' ')
	k=index(logicfile,' ')
	truefile(1:j-1)=value(1:j-1)
	truefile(j:j)='/'
	truefile(j+1:j+k-1)=logicfile(1:k-1)
	istat=0
      endif
      return
      end

! IBM libmassv compatibility

      subroutine vexp(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=exp(x(j))
   10 continue
      return
      end

      subroutine vlog(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=log(x(j))
   10 continue
      return
      end

      subroutine vrec(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=1.d0/x(j)
   10 continue
      return
      end

      subroutine vsqrt(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=sqrt(x(j))
   10 continue
      return
      end
