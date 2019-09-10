      subroutine outsigs ( itau,nx,my,my_max,lev,ncld                &
                         , idtg,ifilout,ptop,rad,grav                &
                         , cp,cosl,pt,sgeo,snr,gwr,tg,pk,pk2         &
                         , ut,vt,tt,qt,phi,rdiv,km,smc               &
                         , slc,stc,canopy,zice,ggdef,gmdef)
      use index
      use mpe
      use radn, only : ntoz

      implicit  none

      integer   itau,nx,my,my_max,lev,ncld,km

      real      ptop,rad,grav,cp

      real      cosl(my),pt(nx,my_max),sgeo(nxp,my_max),        &
                snr(nxp,my_max),gwr(nxp,my_max),                &
                tg(nxp,my_max),pk(nxp,lev,my_max),              &
                pk2(nxp,lev,my_max),ut(nxp,lev,my_max),         &
                vt(nxp,lev,my_max),tt(nxp,lev,my_max),          &
                qt(nxp,lev*ncld,my_max),phi(nxp,lev,my_max),    &
                rdiv(nxp,lev,my_max),work(nx,my),               &
                smc(nxp,km,my_max),stc(nxp,km,my_max),          &
                canopy(nxp,my_max),slc(nxp,km,my_max),          &
                zice(nxp,my_max)
      integer*8 idtg
      character*80 ifilout
      character typ*6,ihdg*26
      character*4 ggdef,gmdef
!
      integer   i,lenc,k,jj,j,nxj,istat,kk,iout_b10,ntrac
      real      xx,capa,pk2top,sfac2,sfac3,sfac4

      do i = 1, nx*my
       work(i,1) = 0.
      enddo
!
      lenc=nx*my
      do 30 k=1,lev
!
!  convert virture potential temperature to temperature
!
      do 20 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 20 i = 1,nxj
       work(i,j)=tt(i,k,jj)*pk(i,k,jj)/(1.0+0.608*qt(i,k,jj))
 20   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"100")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!  convert gaussain u,v component to normal u,v component
!
      do 21 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       xx = rad/cosl(j)
      do 21 i = 1,nxj
        work(i,j)=ut(i,k,jj)*xx
 21   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"200")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 22 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       xx = rad/cosl(j)
      do 22 i = 1,nxj
        work(i,j)=vt(i,k,jj)*xx
 22   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"210")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 23 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 23 i = 1,nxj
        work(i,j)=qt(i,k,jj)
 23   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"500")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      if( ncld .ge. 2 ) then
      do ntrac=2,ncld
      kk = (ntrac-1)*lev+k
      do 26 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 26 i = 1,nxj
        work(i,j)=qt(i,kk,jj)
 26   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      if(ntrac.eq.2)then
        write(typ,'("m",i2.2,"550")')k     ! cloud liquid water content
      else if(ntrac.eq.ntoz)then
        write(typ,'("m",i2.2,"560")')k     ! ozone
      else
        goto 27
      endif
      call syslbl (typ,idtg,itau,gmdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
 27   continue
      enddo
      end if
!
!  stop outputting rdiv (26/12/2000)
!
!     do 24 jj = 1, jlistnum
!      j=jlist1(jj)
!      nxj=nxdef(j)
!     do 24 i = 1, nxj
!       work(i,j)=rdiv(i,k,jj)
!24   continue
!
!     call mpe_unify(work,nx,my,2,mpe_double)
!
!     write(typ,'("m",i2.2,"230")')k
!     call syslbl (typ,idtg,itau,gmdef,ihdg)
!     if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!ccc
!
 30   continue
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
      do 32 k = 1, lev
      do 28 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 28 i=1,nxj
        work(i,j)=phi(i,k,jj)/grav
 28   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"000")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
 32   continue
!
      do 31 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 31 i = 1,nxj
       work(i,j)=pt(i,jj)+ptop
 31   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      call syslbl ('b00010',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
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
       work(i,j) = tt(i,lev,jj)*pk(i,lev,jj)/(1.0+0.608*qt(i,lev,jj))
   40 continue
      call mpe_unify(work,nx,my,2,mpe_double)
      call syslbl ('b00100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 41 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       xx = rad/cosl(j)
      do 41 i = 1,nxj
       work(i,j)=ut(i,lev,jj)*xx
   41 continue
      call mpe_unify(work,nx,my,2,mpe_double)
      call syslbl ('b00200',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 42 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       xx = rad/cosl(j)
      do 42 i = 1,nxj
       work(i,j)=vt(i,lev,jj)*xx
   42 continue
      call mpe_unify(work,nx,my,2,mpe_double)
      call syslbl ('b00210',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
   90 continue
!--------------
!
      call syslbl ('b00650',idtg,itau,ggdef,ihdg)
!      work = snr
      call mpe2d_unify(work,snr)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call syslbl ('s005a1',idtg,itau,ggdef,ihdg)
!      work = gwr
      call mpe2d_unify(work,gwr)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call syslbl ('s00100',idtg,itau,ggdef,ihdg)
!      work = tg
      call mpe2d_unify(work,tg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
! output canopy
!
      call syslbl ('s005c0',idtg,itau,ggdef,ihdg)
!      work = canopy
      call mpe2d_unify(work,canopy)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      call syslbl ('w00092',idtg,itau,ggdef,ihdg)
!      work = zice
      call mpe2d_unify(work,zice)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
! s01100 & s015b0
! additional output for gsi (will be remove after gsi modify)
!
      do k=1,1
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       work(i,j)=smc(i,k,jj)
      enddo
      enddo
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("s0",i1.1,"5b0")')k
      call syslbl (typ,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       work(i,j)=stc(i,k,jj)
      enddo
      enddo
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("s0",i1.1,"100")')k
      call syslbl (typ,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
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
       work(i,j)=smc(i,2,jj)*sfac2+smc(i,3,jj)*sfac3 &
                 +smc(i,4,jj)*sfac4
      enddo
      enddo
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("s0",i1.1,"5b0")')k
      call syslbl (typ,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
! output stc(2,10-200cm)
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       work(i,j)=stc(i,2,jj)*sfac2+stc(i,3,jj)*sfac3 &
                 +stc(i,4,jj)*sfac4
      enddo
      enddo
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("s0",i1.1,"100")')k
      call syslbl (typ,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      enddo

      do 200 k = 1, km
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       work(i,j)=smc(i,k,jj)
      enddo
      enddo
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("l0",i1.1,"5b0")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       work(i,j)=slc(i,k,jj)
      enddo
      enddo
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("l0",i1.1,"5b1")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       work(i,j)=stc(i,k,jj)
      enddo
      enddo
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("l0",i1.1,"100")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (work(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(work,nx,my,5,mpe_double)
      endif
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
 200  continue
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
      character ifilout*80
      real      work(nx,my)
      integer*8 idtg
      character typ*6,ihdg*26,ggdef*4
!
      integer   lenc,istat,jj,j,ii,nxj,i

      lenc=nx*my
!
      call syslbl ('b00650',idtg,itau,ggdef,ihdg)
      call dmsread(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
      do i=1,nxj
         snr(i,jj)=work(ii,j)
         ii=ii+1
      enddo
      enddo

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
      call syslbl ('s00100',idtg,itau,ggdef,ihdg)
      call dmsread(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
      do i=1,nxj
         tg(i,jj)=work(ii,j)
         ii=ii+1
      enddo
      enddo
!
      call syslbl ('w00092',idtg,itau,ggdef,ihdg)
      call dmsread(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
      do i=1,nxj
         zice(i,jj)=work(ii,j)
         ii=ii+1
      enddo
      enddo
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

      character ifilout*80
      real      work(nx,my)
      integer*8 idtg
      character typ*6,ihdg*26
!
      integer   lenc,i,j,k,ii,jj,nxj,istat

      lenc=nx*my
!
      call syslbl ('s005c0',idtg,itau,ggdef,ihdg)
      call dmsread(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
      do i=1,nxj
         canopy(i,jj)=work(ii,j)
         ii=ii+1
      enddo
      enddo

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
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      call dmsread(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
       ii=nxjstart(j)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       smc(i,k,jj)=work(ii,j)
       ii=ii+1
      enddo
      enddo
!
! read slc  l015b1
!
      write(typ,'("l0",i1.1,"5b1")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      call dmsread(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
       ii=nxjstart(j)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       slc(i,k,jj)=work(ii,j)
       ii=ii+1
      enddo
      enddo
!
! read stc  l01100
!
      write(typ,'("l0",i1.1,"100")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      call dmsread(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
       ii=nxjstart(j)
       nxj=nxdef_2d(j)
      do i = 1,nxj
       stc(i,k,jj)=work(ii,j)
       ii=ii+1
      enddo
      enddo
!
      enddo
!
      return
      end
