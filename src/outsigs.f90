      subroutine outsigs ( itau,nx,my,my_max,lev,ncld                &
                         , idtg,ptop,rad,grav                        &
                         , cp,cosl,pt,sgeo,snr,gwr,tg,pk,pk2         &
                         , ut,vt,tt,qt,phi,km,smc                    &
                         , slc,stc,canopy,zice,ggdef,gmdef)
      use index
      use mpe
      use radn, only : ntoz,ntcw,ntrw,ntiw,ntsw,ntgl,nthl,ntinc,ntrnc
      use const, only : RTYPE,kflag,nmmiph

      implicit  none

      integer   itau,nx,my,my_max,lev,ncld,km,nc,kl

      real      ptop,rad,grav,cp

      real      snr(nxp,my_max),gwr(nxp,my_max),                &
                tg(nxp,my_max),                                 &
                smc(nxp,km,my_max),stc(nxp,km,my_max),          &
                canopy(nxp,my_max),slc(nxp,km,my_max),          &
                zice(nxp,my_max)
      real(kind=RTYPE) ut(nxp,lev,my_max),vt(nxp,lev,my_max),   &
                       tt(nxp,lev,my_max),phi(nxp,lev,my_max),  &
                       qt(nxp,lev*ncld,my_max),                 &
                       pt(nxp,my_max),sgeo(nxp,my_max),         &
                       pk(nxp,lev,my_max),pk2(nxp,lev,my_max),  &
                       mout(nx,my),work(nx,my),cosl(my),        &
                       wrk1(nxp,my_max)
      integer*8 idtg
      character typ*6,mlayer*1
      character*4 ggdef,gmdef
!
      integer   i,lenc,k,jj,j,nxj,istat,kk,iout_b10,ntrac,nclds
      real      xx,capa,pk2top,sfac2,sfac3,sfac4

       work = 0.
       mout = 0.
!
      lenc=nx*my
      if ( ntoz .gt. 0 ) then
        nclds=ntoz-1
      else
        nclds=ncld
      endif

!
!  convert virture potential temperature to temperature
!
      nc= 0
      do k=1,lev
        if ( k .lt. 100 ) then
          write(typ,'("m",i2.2,"100")')k
        else
          write(typ,'("n",i2.2,"100")')mod(k,100)
        endif
        call syslbl_w (typ,idtg,itau,gmdef)
!
        do 20 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
        do 20 i = 1,nxj
          wrk1(i,jj)=tt(i,k,jj)*pk(i,k,jj)/(1.0+0.608*qt(i,k,jj))
 20     continue
        call unify_reduceintp(nx,my,my_max,wrk1,work)
        call split(nx,my,lenc,nc,work,mout)
      enddo
!
      if ( myrank .lt. nc ) call dmswrit_split(nx,my,lenc,kflag,mout,istat)

!
!  convert gaussain u,v component to normal u,v component
!
      nc= 0
      do k=1,lev
        if ( k .lt. 100 ) then
          write(typ,'("m",i2.2,"200")')k
        else
          write(typ,'("n",i2.2,"200")')mod(k,100)
        endif
        call syslbl_w (typ,idtg,itau,gmdef)
!        
        do 21 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xx = rad/cosl(j)
        do 21 i = 1,nxj
          wrk1(i,jj)=ut(i,k,jj)*xx
 21     continue
        call unify_reduceintp(nx,my,my_max,wrk1,work)
        call split(nx,my,lenc,nc,work,mout)
      enddo
!
      if ( myrank .lt. nc ) call dmswrit_split(nx,my,lenc,kflag,mout,istat)
!
      nc= 0
      do k=1,lev
        if ( k .lt. 100 ) then
          write(typ,'("m",i2.2,"210")')k
        else
          write(typ,'("n",i2.2,"210")')mod(k,100)
        endif
!
        call syslbl_w (typ,idtg,itau,gmdef)
        do 22 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xx = rad/cosl(j)
        do 22 i = 1,nxj
          wrk1(i,jj)=vt(i,k,jj)*xx
 22     continue
        call unify_reduceintp(nx,my,my_max,wrk1,work)
        call split(nx,my,lenc,nc,work,mout)
      enddo
!
      if ( myrank .lt. nc ) call dmswrit_split(nx,my,lenc,kflag,mout,istat)
!
      nc= 0
      do k=1,lev
        if ( k .lt. 100 ) then
          write(typ,'("m",i2.2,"500")')k
        else
          write(typ,'("n",i2.2,"500")')mod(k,100)
        endif
!
        call syslbl_w (typ,idtg,itau,gmdef)
        do 23 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
        do 23 i = 1,nxj
          wrk1(i,jj)=qt(i,k,jj)
 23     continue
        call unify_reduceintp(nx,my,my_max,wrk1,work)
        call split(nx,my,lenc,nc,work,mout)
      enddo
!
      if ( myrank .lt. nc ) call dmswrit_split(nx,my,lenc,kflag,mout,istat)
!
! combine cloud water and cloud ice together once there is consideration of 
! more hydrometeors in microphysic scheme.
!
      nc= 0
      do k=1,lev
        if ( k .lt. 100 ) then
          write(typ,'("m",i2.2,"550")')k
        else
          write(typ,'("n",i2.2,"550")')mod(k,100)
        endif
        call syslbl_w (typ,idtg,itau,gmdef)
!
        wrk1=0.
        do ntrac=2,nclds
          if ( .not. (nmmiph.eq.18 .and. ntrac.gt.6) ) then  !exclude ice&rain concentration for 2M Thompson
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i = 1,nxj
              wrk1(i,jj)=wrk1(i,jj)+qt(i,k+(ntrac-1)*lev,jj)
            enddo
          enddo
          endif
        enddo
        call unify_reduceintp(nx,my,my_max,wrk1,work)
        call split(nx,my,lenc,nc,work,mout)
      enddo
!
      if ( myrank .lt. nc ) call dmswrit_split(nx,my,lenc,kflag,mout,istat)
!
! output all hydrometeors and ozone one by one
!
      if( nclds .gt. 2 ) then
        do ntrac=2,nclds
          nc= 0
          do k=1,lev
            if ( k .lt. 100 ) then
              kl     = k
              mlayer = 'm'
            else
              kl     = mod(k,100)
              mlayer = 'n'
            endif
            if(ntrac.eq.ntcw)then
              write(typ,'(A1,i2.2,"551")')mlayer,kl     ! cloud water
            else if(ntrac.eq.ntiw)then
              write(typ,'(A1,i2.2,"552")')mlayer,kl     ! cloud ice
            else if(ntrac.eq.ntrw)then
              write(typ,'(A1,i2.2,"553")')mlayer,kl     ! rain
            else if(ntrac.eq.ntsw)then
              write(typ,'(A1,i2.2,"554")')mlayer,kl     ! snow 
            else if(ntrac.eq.ntgl)then
              write(typ,'(A1,i2.2,"555")')mlayer,kl     ! graupel
            else if(ntrac.eq.nthl .and. nmmiph.eq.16) then
              write(typ,'(A1,i2.2,"556")')mlayer,kl     ! hail
            else if(ntrac.eq.ntinc .and. nmmiph.eq.18) then
              write(typ,'(A1,i2.2,"572")')mlayer,kl     ! ice concentration
            else if(ntrac.eq.ntrnc .and. nmmiph.eq.18) then
              write(typ,'(A1,i2.2,"573")')mlayer,kl     ! rain concentration
            else
              goto 27
            endif
            call syslbl_w (typ,idtg,itau,gmdef)
!
            kk = (ntrac-1)*lev+k
            do 26 jj = 1, jlistnum
              j=jlist1(jj)
              nxj=nxdef_2d(j)
            do 26 i = 1,nxj
              wrk1(i,jj)=qt(i,kk,jj)
 26         continue
            call unify_reduceintp(nx,my,my_max,wrk1,work)
            call split(nx,my,lenc,nc,work,mout)
          enddo
!
          if ( myrank .lt. nc ) call dmswrit_split(nx,my,lenc,kflag,mout,istat)
 27       continue
        enddo
      end if
!
! output ozone
!
      if ( ntoz .eq. ncld ) then
        nc= 0
        do k=1,lev
          if ( k .lt. 100 ) then
            write(typ,'("m",i2.2,"560")')k
          else
            write(typ,'("n",i2.2,"560")')mod(k,100)
          endif
          call syslbl_w (typ,idtg,itau,gmdef)
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i = 1,nxj
              wrk1(i,jj)=qt(i,k+(ntoz-1)*lev,jj)
            enddo
          enddo
          call unify_reduceintp(nx,my,my_max,wrk1,work)
          call split(nx,my,lenc,nc,work,mout)
        enddo
!
        if ( myrank .lt. nc ) call dmswrit_split(nx,my,lenc,kflag,mout,istat)
      endif


!
!  stop outputting rdiv (26/12/2000)
!
!     do 24 jj = 1, jlistnum
!      j=jlist1(jj)
!      nxj=nxdef(j)
!     do 24 i = 1, nxj
!       wrk1(i,jj)=rdiv(i,k,jj)
!24   continue
!
!!byl     call mpe_unify(work,nx,my,2,mpe_double)
!
!     write(typ,'("m",i2.2,"230")')k
!     call syslbl (typ,idtg,itau,gmdef,ihdg)
!      call unify_reduceintp(nx,my,my_max,wrk1,work)
!!byl     if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
!     call dmswrit(nx,my,ihdg,lenc,kflag,work,istat)
!ccc
!
!
! compute geopotential by hydrstatic
!
      capa=1.0/3.5
      pk2top=(ptop/1000.)**capa
!
      do jj =1,jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
!
!  hydrostatic equation
!
      do i=1,nxj
      phi(i,lev,jj)= cp*tt(i,lev,jj)*(pk2(i,lev,jj)-pk(i,lev,jj)) &
                     + sgeo(i,jj)
      enddo
      do k=lev-1,1,-1
      do i=1,nxj
      phi(i,k,jj)= phi(i,k+1,jj)+cp*(tt(i,k,jj)*(pk2(i,k,jj)-pk(i,k,jj)) &
                   + tt(i,k+1,jj)*(pk(i,k+1,jj)-pk2(i,k,jj)))
      enddo
      enddo
!
      enddo
!
!  convert geopotential to geopotential hight
!
      nc= 0
      do 32 k = 1, lev
       if ( k .lt. 100 ) then
         write(typ,'("m",i2.2,"000")')k
       else
         write(typ,'("n",i2.2,"000")')mod(k,100)
       endif
       call syslbl_w (typ,idtg,itau,gmdef)
      do 28 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 28 i=1,nxj
        wrk1(i,jj)=phi(i,k,jj)/grav
 28   continue
!      write(typ,'("m",i2.2,"000")')k
!      call syslbl (typ,idtg,itau,gmdef,ihdg)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
 32   continue
!
      if ( myrank .lt. nc ) call dmswrit_split(nx,my,lenc,kflag,mout,istat)
!
!----- start to output surface data ------
      nc=0
!
      do 31 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 31 i = 1,nxj
       wrk1(i,jj)=pt(i,jj)+ptop
 31   continue
      call syslbl_w ('b00010',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
!  move the output of "b10, "b20" and "b21" to out2d.f (2002/4/29)
!
      iout_b10 = 0
      if( iout_b10 .eq. 0 ) go to 90
!-----------
      do 40 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 40 i = 1,nxj
       wrk1(i,jj) = tt(i,lev,jj)*pk(i,lev,jj)/(1.0+0.608*qt(i,lev,jj))
   40 continue
      call syslbl_w ('b00100',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
      do 41 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       xx = rad/cosl(j)
      do 41 i = 1,nxj
       wrk1(i,jj)=ut(i,lev,jj)*xx
   41 continue
      call syslbl_w ('b00200',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
      do 42 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       xx = rad/cosl(j)
      do 42 i = 1,nxj
       wrk1(i,jj)=vt(i,lev,jj)*xx
   42 continue
      call syslbl_w ('b00210',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
   90 continue
!--------------
!
      wrk1 = snr
      call syslbl_w ('b00650',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)

      wrk1 = gwr
      call syslbl_w ('s005a1',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)


      wrk1 = tg
      call syslbl_w ('s00100',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
! output canopy
!
      wrk1 = canopy
      call syslbl_w ('s005c0',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
      wrk1 = zice
      call syslbl_w ('w00092',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
! s01100 & s015b0
! additional output for gsi (will be remove after gsi modify)
!
      do k=1,1
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       wrk1(i,jj)=smc(i,k,jj)
      enddo
      enddo
      write(typ,'("s0",i1.1,"5b0")')k
      call syslbl_w (typ,idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       wrk1(i,jj)=stc(i,k,jj)
      enddo
      enddo
      write(typ,'("s0",i1.1,"100")')k
      call syslbl_w (typ,idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
      enddo
!
!
! s02100 & s025b0
! additional output for other model's requirement
!
      sfac2=3./19.
      sfac3=6./19.
      sfac4=10./19.
! output smc(2,10-200cm)
      do k=2,2
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       wrk1(i,jj)=smc(i,2,jj)*sfac2+smc(i,3,jj)*sfac3 &
                 +smc(i,4,jj)*sfac4
      enddo
      enddo
      write(typ,'("s0",i1.1,"5b0")')k
      call syslbl_w (typ,idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
! output stc(2,10-200cm)
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       wrk1(i,jj)=stc(i,2,jj)*sfac2+stc(i,3,jj)*sfac3 &
                 +stc(i,4,jj)*sfac4
      enddo
      enddo
      write(typ,'("s0",i1.1,"100")')k
      call syslbl_w (typ,idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
      enddo

      do 200 k = 1, km
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       wrk1(i,jj)=smc(i,k,jj)
      enddo
      enddo
      write(typ,'("l0",i1.1,"5b0")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       wrk1(i,jj)=slc(i,k,jj)
      enddo
      enddo
      write(typ,'("l0",i1.1,"5b1")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       wrk1(i,jj)=stc(i,k,jj)
      enddo
      enddo
      write(typ,'("l0",i1.1,"100")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call unify_reduceintp(nx,my,my_max,wrk1,work)
      call split(nx,my,lenc,nc,work,mout)
!
 200  continue
!
      if ( myrank .lt. nc )             &
         call dmswrit_split(nx,my,lenc,kflag,mout,istat)
!
      return
      end

      subroutine rdpbl(nx,my,my_max,snr,gwr,tg,zice,ifilout,idtg,itau,ggdef)
!
      use index
      use mpe
!
      implicit  none
      integer   nx,my,itau,my_max

      real      snr(nxp,my_max),gwr(nxp,my_max),tg(nxp,my_max),zice(nxp,my_max)
      character ifilout*255
      real      work(nx,my)
      integer*8 idtg
      character typ*6,ggdef*4
!
      integer   lenc,istat,jj,j,ii,nxj,i

      lenc=nx*my
!
      call syslbl_r ('b00650',idtg,itau,ggdef)
      call dmsread(nx,my,lenc,'H',ifilout,work,istat)
      call unify_reducepick(nx,my,my_max,work,snr)

!
!!      call syslbl ('s005a1',idtg,itau,ggdef,ihdg)
!!      call dmsread(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!!!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
!!      do jj=1,jlistnum
!!         j=jlist1(jj)
!!         if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
!!         ii=nxjstart(j)
!!         nxj=nxdef_2d(j)
!!      do i=1,nxj
!!         gwr(i,jj)=work(ii,j)
!!         ii=ii+1
!!      enddo
!!      enddo
!
      call syslbl_r ('s00100',idtg,itau,ggdef)
      call dmsread(nx,my,lenc,'H',ifilout,work,istat)
      call unify_reducepick(nx,my,my_max,work,tg)
!
      call syslbl_r ('w00092',idtg,itau,ggdef)
      call dmsread(nx,my,lenc,'H',ifilout,work,istat)
      call unify_reducepick(nx,my,my_max,work,zice)
!
      return
      end

!-----
!  following subroutine added for new soil model
!-----
      subroutine rdsoil(nx,my,my_max,km,smc,stc,slc,canopy,ifilout &
                       ,idtg,itau,ggdef,gmdef)
!
      use index
      use mpe

      implicit  none

      integer   nx,my,my_max,km,itau

      real        smc(nxp,km,my_max),stc(nxp,km,my_max),canopy(nxp,my_max), &
                  slc(nxp,km,my_max)

      character*4 ggdef,gmdef

      character ifilout*255
      real      work(nx,my),tmp(nxp,my_max)
      integer*8 idtg
      character typ*6
!
      integer   lenc,i,j,k,ii,jj,nxj,istat

      lenc=nx*my
!
      call syslbl_r ('s005c0',idtg,itau,ggdef)
      call dmsread(nx,my,lenc,'H',ifilout,work,istat)
      call unify_reducepick(nx,my,my_max,work,canopy)

!
! for Noah 4-layer land model
! dmskey(1:3) 1st layer :l01
! dmskey(1:3) 2nd layer :l02
! dmskey(1:3) 3rd layer :l03
! dmskey(1:3) 4th layer :l04

      do k = 1,km
!
! read smc  l015b0
!
      write(typ,'("l0",i1.1,"5b0")')k
      call syslbl_r (typ,idtg,itau,gmdef)
      call dmsread(nx,my,lenc,'H',ifilout,work,istat)
      call unify_reducepick(nx,my,my_max,work,tmp)
      smc(:,k,:) = tmp(:,:)
!
! read slc  l015b1
!
      write(typ,'("l0",i1.1,"5b1")')k
      call syslbl_r (typ,idtg,itau,gmdef)
      call dmsread(nx,my,lenc,'H',ifilout,work,istat)
      call unify_reducepick(nx,my,my_max,work,tmp)
      slc(:,k,:) = tmp(:,:)
!
! read stc  l01100
!
      write(typ,'("l0",i1.1,"100")')k
      call syslbl_r (typ,idtg,itau,gmdef)
      call dmsread(nx,my,lenc,'H',ifilout,work,istat)
      call unify_reducepick(nx,my,my_max,work,tmp)
      stc(:,k,:) = tmp(:,:)
!
!
      enddo
!
      return
      end
