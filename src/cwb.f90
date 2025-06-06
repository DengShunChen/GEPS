      subroutine sendmsg(str_model,i_fromtau,i_totau,i_status)

      integer i_fromtau,i_totau,i_status
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
      
      integer tau(9999),cnt,ip,istatus,ifrom,ito
      save    tau,cnt,ip

      istatus=0
      if(ip.eq.0)then
        open(1,file=model//'ctl',form='formatted',status='old')
        do i=1,9999
           read(1,*,end=100)tau(i)
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
            else
              month(2)=28
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
            else
              month(2)=28
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
      character pathname*255,logicfile*255,truefile*255,value*255,blank*255
!     integer getenv

      truefile(1:255)=' '
      value(1:255)=' '
      blank(1:255)=' '
!     i = getenv(pathname,value)
!     if (i .eq. 0) then
      ix=index(pathname,' ')
      call getenv(pathname(1:ix-1),value)
      if (value(1:255) .eq. blank(1:255)) then
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
!
      subroutine ptotc(pdrym,lprint)

      use index
      use rank
      use const
      use param
      use grid

      implicit none

      integer i,j,k,n,jj,kk,nxj,nxjf,kn
      logical lprint
      real    sumtot   ,sumwat    ,dsigp                    &
             ,sumtotm  ,sumwatm   ,pdrym                    &
             ,qtot
      real(kind=RTYPE), dimension(:,:), allocatable ::      &
                       sumtotp,sumwatp,sumtottx,sumwattx
      real(kind=RTYPE), dimension(:), allocatable ::        &
                       sumtotpy,sumwatpy,sumtotty,sumwatty


      allocate(sumtotp(nxp,my_max),sumwatp(nxp,my_max),     &
               sumtottx(nx,my_max),sumwattx(nx,my_max))
      ! dry air mass conservation
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
!        nxjf=nxdef(j)
        do i = 1,nxj
          sumtotp(i,jj)=0.
          sumwatp(i,jj)=0.
          do k = 1, lev
            dsigp = dsigma(k,1)*pt(i,jj)+dsigma(k,2)
            qtot=0.
            do n = 1, ncld
               kk=k+(n-1)*lev
               qtot = qtot + qt(i,kk,jj)
            enddo
            sumtotp(i,jj) = sumtotp(i,jj)  + dsigp
            sumwatp(i,jj) = sumwatp(i,jj)  + dsigp * qtot
          enddo
        enddo
      enddo

      call mpe2d_unify_nx(sumtottx,sumtotp)
      call mpe2d_unify_nx(sumwattx,sumwatp)

      deallocate(sumtotp,sumwatp)
      allocate(sumtotpy(my_max),sumwatpy(my_max))

      sumtotpy=0.
      sumwatpy=0.
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxjf=nxdef(j)
        do i = 1, nxjf
          sumtotpy(jj) = sumtotpy(jj) + sumtottx(i,jj)*cosl(j)*nx/nxjf
          sumwatpy(jj) = sumwatpy(jj) + sumwattx(i,jj)*cosl(j)*nx/nxjf
        enddo
      enddo

      deallocate(sumtottx,sumwattx)
      allocate(sumtotty(my),sumwatty(my))

      call mpe2d_unify_my1d(sumtotty,sumtotpy)
      call mpe2d_unify_my1d(sumwatty,sumwatpy)

      deallocate(sumtotpy,sumwatpy)

      sumtot=0.
      sumwat=0.
      do j = 1, my
        sumtot = sumtot + sumtotty(j)
        sumwat = sumwat + sumwatty(j)
      enddo

      deallocate(sumtotty,sumwatty)

      kn      = nx * my
      sumtotm = sumtot / float(kn)
      sumwatm = sumwat / float(kn)
      pdrym   = sumtotm - sumwatm

      if( myrank .eq. 0 .and. lprint ) then
        open(35,file='pdry.txt',form='formatted',status='unknown', &
             position='append')
        write(35,*)pdrym,sumwatm,sumtotm
        close(35)
        print*,'dry air mass = ',pdrym,' hPa'
        print*,'water mass = ',sumwatm,' hPa'
        print*,'total air mass = ',sumtotm,' hPa'
      endif

      return
      end
!
      subroutine adjptq(qt,qp,pnew,pten,nxj,nxp,my_max,lev,ncld,dta)

      use const, only : RTYPE,dsigma

      implicit none

      integer          nxj,nxp,my_max,lev,ncld
      real             dta
      real(kind=RTYPE) qt(nxp,lev*ncld),qp(nxp,lev*ncld), &
                       pnew(nxp),pten(nxp)

      integer i,k,n,jj,kk
      real(kind=RTYPE) pt(nxp)
      real    dsigp    ,qtot      ,qtota                               &
             ,sumtott  ,sumwatt   ,sumwatta                            &
             ,odpondp 
      
      pt(1:nxj) = pnew(1:nxj)
      ! adjustment of surface pressure
        do i = 1,nxj
          sumtott=0.
          sumwatt=0.
          sumwatta=0.
          do k = 1, lev
            dsigp = dsigma(k,1)*pt(i)+dsigma(k,2)
            qtot=0.
            qtota=0.
            do n = 1, ncld
               kk=k+(n-1)*lev
               qtot = qtot  + qp(i,kk)
               qtota= qtota + qt(i,kk)
            enddo
            sumtott = sumtott  + dsigp
            sumwatt = sumwatt  + dsigp * qtot
            sumwatta= sumwatta + dsigp * qtota
          enddo
          pnew(i) = sumtott - sumwatt + sumwatta
        enddo

      ! mass adjustment of all tracers and virtual potential temperature

        do i = 1,nxj
          do k = 1, lev
            odpondp =  (dsigma(k,1)*pt(i)+dsigma(k,2))              &
                     / (dsigma(k,1)*pnew(i)+dsigma(k,2))
            do n = 1, ncld
              kk=k+(n-1)*lev
              qt(i,kk) = qt(i,kk) * odpondp   ! global method
            enddo
          enddo
          pten(i) = ( pnew(i) - pten(i) ) / dta
        enddo

      return
      end

      subroutine adjptqintp(ut,vt,tt,qt,qp,pnew,pten,nxj,nxp,my_max, &
                            lev,ncld,dta)

      use const, only : RTYPE,dsigma,sigma,qmin
      use grid,  only : ndslvvar

      implicit none

      integer          nxj,nxp,my_max,lev,ncld
      real             dta
      real(kind=RTYPE) ut(nxp,lev),vt(nxp,lev),tt(nxp,lev),          &
                       qt(nxp,lev*ncld),qp(nxp,lev*ncld),            &
                       pnew(nxp),pten(nxp),pold(nxp)

      integer          i,k,n,kk,ki
      integer          kuu,kvv,ktt
      real(kind=RTYPE) plnew(lev+1),plold(lev+1)
      real(kind=RTYPE) rqda(lev,ndslvvar),rqnn(lev,ndslvvar)
      real             dsigp,qtot,qtota,odpondp,ptmp

      kuu = 1
      kvv = kuu + 1
      ktt = kvv + 1
      
      pold(1:nxj) = pnew(1:nxj)
      ! new surface pressure
        do i = 1,nxj
          plnew(lev+1) = 0.
          plold(lev+1) = 0.
          ptmp         = 0.
          do k = 1, lev
            ki=lev-k+1
            dsigp = dsigma(k,1)*pold(i)+dsigma(k,2)
            plold(ki) = plold(ki+1)-dsigp
            qtot=0.
            qtota=0.
            do n = 1, ncld
              kk=k+(n-1)*lev
              qtot = qtot  + qp(i,kk)
              qtota= qtota + qt(i,kk)
            enddo
            ptmp = ptmp+dsigp*(1.-qtot+qtota)
          enddo
          pnew(i) = ptmp
          pten(i) = ( pnew(i) - pten(i) ) / dta

          ! mass adjustment of all tracers and virtual potential temperature
          do k = 1, lev
            ki=lev-k+1
            dsigp = dsigma(k,1)*pnew(i)+dsigma(k,2)
            plnew(ki) = plnew(ki+1)-dsigp
            odpondp =  (dsigma(k,1)*pold(i)+dsigma(k,2)) / dsigp 
            rqda(ki,kuu) = ut(i,k) * odpondp
            rqda(ki,kvv) = vt(i,k)
            rqda(ki,ktt) = tt(i,k)
            do n = 1, ncld
              kk=k+(n-1)*lev
              rqda(ki,ktt+n) = max(qt(i,kk) * odpondp,qmin)
            enddo
          enddo

          call fixend_cell_plm_intp(plold,rqda,plnew,rqnn,lev,ndslvvar)
          
          do k = 1, lev
            ki=lev-k+1
            ut(i,k) = rqnn(ki,kuu)
            vt(i,k) = rqnn(ki,kvv)
            tt(i,k) = rqnn(ki,ktt)
            do n = 1, ncld
              kk=k+(n-1)*lev
              qt(i,kk) = rqnn(ki,ktt+n)
            enddo
          enddo

        enddo

      return
      end
