      subroutine rsmout_sigp ( itau,nx,my,my_max,lev,ncld            &
                         , idtg,ptop,rad,grav,cosl                   &
                         , pt,sgeo,snr,gwr,tg,pk                     &
                         , ut,vt,tt,qt,km,smc,stc                    &
                         , ice,land,ocean,xlon,xlat)
#ifdef RSM

#if defined(CWB_MPMD)  || defined(CWBSUM)
#define send_RSM
#else
#define write_RSM
#endif

      use index
      use mpe
      use radn, only : ntoz,ntcw,ntrw,ntiw,ntsw,ntgl
      use noah, only : cice
#if defined(CWB_MPMD) || defined(CWBSUM)
      use rank, only : root_rsm,myrank,itag
#else
      use rank, only : myrank
#endif
      use const, only: rlon1, rlon2, rlat1, rlat2, rgrdsz, rsmsfcmgrhr &
                ,RTYPE

      implicit  none

      logical   ldosfc

      integer   itau,nx,my,my_max,lev,ncld,km

      real      ptop,rad,grav

      real(kind=RTYPE)::ut(nxp,lev,my_max),vt(nxp,lev,my_max) &
          ,tt(nxp,lev,my_max),qt(nxp,lev*ncld,my_max)         &
          ,pt(nxp,my_max),pk(nxp,lev,my_max),sgeo(nxp,my_max)

      real(kind=RTYPE)::wrk1(nxp,my_max),work(nx,my)

      real(kind=RTYPE)      cosl(my)
      real      snr(nxp,my_max),gwr(nxp,my_max),                &
                tg(nxp,my_max),                                 &
                smc(nxp,km,my_max),stc(nxp,km,my_max),          &
                workr8(nx,my)

      logical   land(nxp,my_max),ocean(nxp,my_max),ice(nxp,my_max)

      real      xlon(nx,my_max),xlat(my)
      real      dx,dy

      integer*8 idtg
!
! local array
      integer   i,k,ii,jj,j,nxj,kk,ntrac,nclds
      real      xx
!
      integer   nxmy, nxmyl, nxmys
      integer   nxs, mys, len2
      integer   x1, x2, y1, y2
#ifdef send_RSM
      real,allocatable ::  g3send(:,:,:),g2send(:,:)
#endif
      integer   ierr

      integer   nsig, ndig, recsize
      character cfhour*16,cform*40,cidtg*12

      integer   itag0

! output record      
      integer   kfh,knt,ktt,ksgeo,ksfcp,kuu,kvv
      integer   kqt,kcw,krw,kiw,ksw,kgl,koz
      integer   ktg,ksmc,ksnr,kstc,ksimk,kslmk
      integer   ktrace,kend,krec
      integer   klon,klat
      integer   ntindex(7)

      ldosfc=mod(float(itau)+0.00001, float(rsmsfcmgrhr) ).lt.0.01

! local variable initization
      nxmy=nx*my
      nxmyl=nx*my*lev
      nxmys=nx*my*km
      do i = 1, nxmy
       work(i,1) = 0.
      enddo
!
      if(myrank .eq. 0) print *,' rsmout_sigp : ntoz=',ntoz
      if ( ntoz .gt. 0 ) then
        nclds=ntoz-1
      else
        nclds=ncld
      endif
      if(myrank .eq. 0) print *,' rsmout_sigp : nclds=',nclds

      if(myrank .eq. 0) print *,' rsmout_sigp : ',                     &
                                'ntoz,ntcw,ntrw,ntiw,ntsw,ntgl=',      &
                                 ntoz,ntcw,ntrw,ntiw,ntsw,ntgl
! list tracer index
        ntindex(1)=nclds
        ntindex(2)=ntoz
        ntindex(3)=ntcw
        ntindex(4)=ntrw
        ntindex(5)=ntiw
        ntindex(6)=ntsw
        ntindex(7)=ntgl

#ifdef write_RSM
! output file name      
      nsig=21
#endif
      ndig=max(log10(itau+0.5)+1.,3.)
      write(cform,'("(i",i1,".",i1,")")') ndig,ndig
      write(cfhour,cform) itau
      write(cidtg,'(I12.12)') idtg

! output record sequence 
      ktrace=-999
!lat/lon
      klon=2
      klat=3
      kfh=1
      knt=kfh+1
      ksgeo=knt+1
      ksfcp=ksgeo+1
      ktt=ksfcp+1
      kuu=ktt+lev
      kvv=kuu+lev
      kqt=kvv+lev
      kcw=kqt+lev
      koz=kcw+lev
      kslmk=koz+lev
      ksnr=kslmk+1
      ksimk=ksnr+1
      ktg=ksimk+1
      ksmc=ktg+1
      kstc=ksmc+4
      krw=kstc+4
      kiw=krw+lev
      ksw=kiw+lev
      kgl=ksw+lev
      kend=kgl+lev

      if(myrank.eq.0) then
        print*,'klon,klat='
        print*,klon,klat
        print*,'ktrace,kfh,knt,ktt,ksgeo,ksfcp,kuu,kvv,kqt,kcw,koz='
        print*,ktrace,kfh,knt,ktt,ksgeo,ksfcp,kuu,kvv,kqt,kcw,koz
        print*,'ktg,ksmc,ksnr,kstc,ksimk,kslmk,krw,kiw,ksw,kgl,kend='
        print*,ktg,ksmc,ksnr,kstc,ksimk,kslmk,krw,kiw,ksw,kgl,kend
      endif
! 
! output record, regional domain
      dx=xlon(2,1)-xlon(1,1)
      dy=xlat(2)-xlat(1)
      rgrdsz=min(dx,dy)
      if (myrank.eq.0) print*,'rgrdsz=',rgrdsz
      nxs=nint((rlon2-rlon1)/rgrdsz+1.)
      mys=nint((rlat2-rlat1)/rgrdsz+1.)
      x1=nint(  (rlon1-0.)     /rgrdsz+1.)
      x2=nxs-1+x1
      y1=nint(  (rlat1-(-90.)) /rgrdsz+1.)
      y2=mys-1+y1
      len2=nxs*mys
!
! print/output information of regional domain at itau=0
!
      if(itau.eq.0) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         ii=nxjstart(j)
        do i = 1, nxj
           wrk1(i,jj)=xlon(ii,jj)
           if (xlon(ii,jj).lt. 0) wrk1(i,jj)=(xlon(ii,jj)+360.)
           ii=ii+1
        enddo
        enddo
        call unify_reduceintp(nx,my,my_max,wrk1,work)

        if(myrank.eq.0) then
          open(99,file='rsm_xlon_'//cidtg//'.txt',status='unknown', &
              form='formatted',iostat=ierr)
          if (ierr/=0) stop "rsmout_sigp: open file xlon fail"
          write(99,*) work(x1:x2,1)
#ifdef send_RSM
        workr8=work
        itag=klon
        call mpmd_send(workr8(x1:x2,1),nxs,root_rsm,itag,'R')
#endif
          close(99)
          open(98,file='rsm_xlat_'//cidtg//'.txt',status='unknown', &
              form='formatted',iostat=ierr)
          if (ierr/=0) stop "rsmout_sigp: open file xlat fail"
          write(98,*) xlat(y1:y2)
#ifdef send_RSM
        itag=klat
        call mpmd_send(xlat(y1:y2),mys,root_rsm,itag,'R')
#endif
          close(98)
          print*,'-----output data for rsm-----'
          print*,'--nxs, mys, x1, x2, y1, y2=--'
          print*,'-',nxs, mys, x1, x2, y1, y2,'-'
        endif
      endif  !endif (itau=0)

#ifdef send_RSM
      if(myrank.eq.0) then
! mpmd send tag
      itag0=itag
      endif
#endif


! allocate temporary 
#ifdef send_RSM
      allocate(g3send(nxs,mys,lev),g2send(nxs,mys), &
               stat=ierr)
      if (ierr/=0) stop "rsmout_sigp: allocate fail"
      if (ierr/=0) call mpe_finalize               
#endif


#ifdef write_RSM
      if(myrank.eq.0) then
        inquire(iolength=recsize) workr8(x1:x2,y1:y2)
        open(nsig,file='rsm_data_'//cidtg//'.f'//cfhour,status='unknown', &
            form='unformatted',access='direct',recl=recsize,iostat=ierr)
        if (ierr/=0) stop "rsmout_sigp: open file fail"
      endif
#endif

#ifdef write_RSM
      if(myrank.eq.0) then
        krec=kfh
        write(nsig,rec=krec) float(itau)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        itag=itag0+kfh
        call mpmd_send(float(itau),1,root_rsm,itag,'R')
        print*,'gfs kfh itag=',itag
      endif
#endif

! gfs tracer index
#ifdef write_RSM
      if(myrank.eq.0) then
        krec=knt
        write(nsig,rec=krec) ntindex
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        itag=itag0+knt
        call mpmd_send(ntindex,7,root_rsm,itag,'I')
        print*,'gfs knt itag=',itag
      endif
#endif
!
!  convert terrain geopotential to terrain geopotential hight
!
      do 28 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 28 i=1,nxj
        wrk1(i,jj)=sgeo(i,jj)/grav
 28   continue
      call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=ksgeo
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g2send(:,:)=work(x1:x2,y1:y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        itag=itag0+ksgeo
        call mpmd_send(g2send,len2,root_rsm,itag,'R')
        print*,'gfs ksgeo itag=',itag
      endif
#endif
!
! surface pressure
!
      do 31 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 31 i = 1,nxj
       wrk1(i,jj)=pt(i,jj)+ptop
 31   continue
      call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=ksfcp
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g2send(:,:)=work(x1:x2,y1:y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        itag=itag0+ksfcp
        call mpmd_send(g2send,len2,root_rsm,itag,'R')
        print*,'gfs ksfcp itag=',itag
      endif
#endif
!
!  convert virture potential temperature to temperature
!
      do k=1,lev
        do 20 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
        do 20 i = 1,nxj
          wrk1(i,jj)=tt(i,k,jj)*pk(i,k,jj)/(1.0+0.608*qt(i,k,jj))
 20     continue
        call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=ktt-1+(lev-k+1)
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g3send(:,:,lev-k+1)=work(x1:x2,y1:y2)
      endif
#endif
      enddo
#ifdef send_RSM
      if(myrank.eq.0) then
      do k=1,lev
        itag=itag0+ktt-1+k
        call mpmd_send(g3send(:,:,k),len2,root_rsm,itag,'R')
      enddo
      endif
#endif

!
!  convert gaussain u,v component to normal u,v component
!
!
! u
!
      do k=1,lev
        do 21 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xx = rad/cosl(j)
        do 21 i = 1,nxj
          wrk1(i,jj)=ut(i,k,jj)*xx
 21     continue
        call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=kuu-1+(lev-k+1)
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g3send(:,:,lev-k+1)=work(x1:x2,y1:y2)
      endif
#endif
      enddo
#ifdef send_RSM
      if(myrank.eq.0) then
      do k=1,lev
        itag=itag0+kuu-1+k
        call mpmd_send(g3send(:,:,k),len2,root_rsm,itag,'R')
      enddo
      endif
#endif
!
! v
!
      do k=1,lev
        do 22 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xx = rad/cosl(j)
        do 22 i = 1,nxj
          wrk1(i,jj)=vt(i,k,jj)*xx
 22     continue
        call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=kvv-1+(lev-k+1)
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g3send(:,:,lev-k+1)=work(x1:x2,y1:y2)
      endif
#endif
      enddo
#ifdef send_RSM
      if(myrank.eq.0) then
      do k=1,lev
        itag=itag0+kvv-1+k
        call mpmd_send(g3send(:,:,k),len2,root_rsm,itag,'R')
      enddo
      endif
#endif
!
! specific humidity
!
      do k=1,lev
        do 23 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
        do 23 i = 1,nxj
          wrk1(i,jj)=qt(i,k,jj)
 23     continue
        call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=kqt-1+(lev-k+1)
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g3send(:,:,lev-k+1)=work(x1:x2,y1:y2)
      endif
#endif
      enddo
#ifdef send_RSM
      if(myrank.eq.0) then
      do k=1,lev
        itag=itag0+kqt-1+k
        call mpmd_send(g3send(:,:,k),len2,root_rsm,itag,'R')
      enddo
      endif
#endif
!
! output all hydrometeors and ozone one by one
!
      if( nclds .ge. 2 ) then
        do ntrac=2,nclds
          if(myrank .eq. 0) print *,' rsmout_sigp : ntrac=',ntrac
          do k=1,lev
            kk = (ntrac-1)*lev+k
            do 26 jj = 1, jlistnum
              j=jlist1(jj)
              nxj=nxdef_2d(j)
            do 26 i = 1,nxj
              wrk1(i,jj)=qt(i,kk,jj)
 26         continue
            call unify_reduceintp(nx,my,my_max,wrk1,work)
!
            if(ntrac.eq.ntcw)then
              ktrace=kcw                                ! cloud liquid water content
            else if(ntrac.eq.ntiw)then
              ktrace=kiw                                ! cloud ice content
            else if(ntrac.eq.ntrw)then
              ktrace=krw                                ! rain
            else if(ntrac.eq.ntsw)then
              ktrace=ksw                                ! snow 
            else if(ntrac.eq.ntgl)then
              ktrace=kgl                                ! graupel
            else
              stop "ERROR: rsmout_sigp: wrong tracer index found !"
              call mpe_finalize               
            endif
#ifdef write_RSM
            if(myrank.eq.0) then
              workr8=work
              krec=ktrace-1+(lev-k+1)
              write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
            endif
#endif
#ifdef send_RSM
            if(myrank.eq.0) then
              g3send(:,:,lev-k+1)=work(x1:x2,y1:y2)
            endif
#endif
            enddo   !enddo k=1,lev
#ifdef send_RSM
            if(myrank.eq.0) then
            do k=1,lev
              itag=itag0+ktrace-1+k
              call mpmd_send(g3send(:,:,k),len2,root_rsm,itag,'R')
            enddo
            endif
#endif
!
        enddo   !enddo ntrac=2,nclds
      end if    !endif (nclds.gt.2)
!
! output ozone
!
      if ( ntoz .eq. ncld ) then 
        do k=1,lev
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i = 1,nxj
              wrk1(i,jj)=qt(i,k+(ntoz-1)*lev,jj)
            enddo
          enddo
          call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=koz-1+(lev-k+1)
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g3send(:,:,lev-k+1)=work(x1:x2,y1:y2)
      endif
#endif
      enddo
#ifdef send_RSM
      if(myrank.eq.0) then
      do k=1,lev
        itag=itag0+koz-1+k
        call mpmd_send(g3send(:,:,k),len2,root_rsm,itag,'R')
      enddo
      endif
#endif
!
      endif
!
!----- start to output surface data ------
      if (ldosfc) then

      if(myrank .eq. 0) print *,' rsmout_sigp : output surface file'
! ***land and sea mask***(land=1,sea=0)***
      wrk1(:,:)=0.0
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
      do i=1,nxj
        if(land(i,jj))wrk1(i,jj)=1.0
        if(ocean(i,jj))wrk1(i,jj)=0.0
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=kslmk
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g2send(:,:)=work(x1:x2,y1:y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        itag=itag0+kslmk
        call mpmd_send(g2send,len2,root_rsm,itag,'R')
      endif
#endif
! ***snr***
      wrk1=snr
      call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=ksnr
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g2send(:,:)=work(x1:x2,y1:y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        itag=itag0+ksnr
        call mpmd_send(g2send,len2,root_rsm,itag,'R')
      endif
#endif
!! ***ice***(simk in RSM-csfcfld(:,13))
      wrk1=cice
      call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=ksimk
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g2send(:,:)=work(x1:x2,y1:y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        itag=itag0+ksimk
        call mpmd_send(g2send,len2,root_rsm,itag,'R')
      endif
#endif
! ***tg***
      wrk1=tg
      call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=ktg
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g2send(:,:)=work(x1:x2,y1:y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        itag=itag0+ktg
        call mpmd_send(g2send,len2,root_rsm,itag,'R')
      endif
#endif
!
! ***0-0.1m soil moisture content [fraction]***
! ***0.1-0.4m soil moisture content [fraction]***
! ***0.4-1m soil moisture content [fraction]***
! ***below 1m soil moisture content [fraction]***
      if (km/=4) stop "wrong soil layers, should be 4 layers NOAHLSM"
      do k=1,km
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       wrk1(i,jj)=smc(i,k,jj)
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=ksmc-1+k
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g3send(:,:,k)=work(x1:x2,y1:y2)
      endif
#endif
      enddo
#ifdef send_RSM
      if(myrank.eq.0) then
      do k=1,km
        itag=itag0+ksmc-1+k
        call mpmd_send(g3send(:,:,k),len2,root_rsm,itag,'R')
      enddo
      endif
#endif
!
! ***0-0.1m soil temperature***
! ***0.1-0.4m soil temperature***
! ***0.4-1m soil temperature***
! ***below 1m soil temperature***
      do k=1,km
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       wrk1(i,jj)=stc(i,k,jj)
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk1,work)
!
#ifdef write_RSM
      if(myrank.eq.0) then
        workr8=work
        krec=kstc-1+k
        write(nsig,rec=krec) ((workr8(i,j),i=x1,x2),j=y1,y2)
      endif
#endif
#ifdef send_RSM
      if(myrank.eq.0) then
        g3send(:,:,k)=work(x1:x2,y1:y2)
      endif
#endif
      enddo
#ifdef send_RSM
      if(myrank.eq.0) then
      do k=1,km
        itag=itag0+kstc-1+k
        call mpmd_send(g3send(:,:,k),len2,root_rsm,itag,'R')
      enddo
      endif
#endif
      endif ! (ldosfc)

!--------------
!
#ifdef write_RSM
      if(myrank.eq.0) then
        close(nsig)
      endif
#endif
#ifdef send_RSM
! mpmd send tag
      itag=itag0+kend
#endif
!
#ifdef send_RSM
      deallocate(g3send,g2send,stat=ierr)
      if (ierr/=0) stop "rsmout_sigp: deallocate fail"
#endif

#endif
      return
      end
