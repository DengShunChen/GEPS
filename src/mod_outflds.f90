module mod_outflds
 implicit none
 
 
contains
  subroutine divgout(nx,my,my_max,lpout,lev,itau,ifilout,idtg, &
        plev,num,whtlev,pkout,pk,pklp,rdiv,rdivb,div,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use mod_grb2_param !for write grb2
      use const ,only:out_pres_form

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,ncnt

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max),       &
                rdiv(nxp,lev,my_max),rdivb(nxp,my_max),div(nxp,my_max,lpout),&
                plev(lpout),whtlev(num)
      real      tens(lev+1)
      real      wk1(nx,my),pout(nx,my)

!
      integer*8 idtg
      character*80 ifilout
      character*26 ihdg,ihdg2
      character*6 lrec(lpout)
      character*4 ggdef
!
      integer    k,lpl,lenc,n,num,istat
  
      logical :: lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'230'
      end do
      lrec(lpout) = 'h00230'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,rdiv,rdivb,pkout,div,tens)
!
      lenc= nx*my
      ncnt= -1
!
      do 20 n=1,num
      do 10 k=1,lpout
      if(plev(k).eq.whtlev(n)) then
      ncnt=ncnt+1
      call unify_reduceintp(nx,my,my_max,div(1,1,k),wk1)
      call syslbl(lrec(k),idtg,itau,ggdef,ihdg)
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,wk1,istat)
      if(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,13,6,100,-2,plev(k),wk1)
      endif
      call qmaxn3(wk1,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=wk1
        ihdg2=ihdg
      endif
      go to 20
      endif
   10 continue
   20 continue
      if(lwrite .and. myrank .le. ncnt .and.out_pres_form==1) &
             call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
      return
  end subroutine divgout

  subroutine dragout(nx,my,my_max,lpout,lev,itau,ifilout,idtg, &
        plev,num,whtlev,pkout,pk,pklp,rdrag,rdragb,drag,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use mpe
      use mod_grb2_param !for write grb2
      use const ,only:out_pres_form

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,num,ncnt

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max),               &
                rdrag(nxp,lev,my_max),rdragb(nxp,my_max),drag(nxp,my_max,lpout),&
                plev(lpout),whtlev(num)
      real      tens(lev+1)
      real      wk1(nx,my),pout(nx,my)
!
      integer*8 idtg
      character*80 ifilout
      character*26 ihdg,ihdg2
      character*6 lrec(lpout)
      character*4 ggdef
!
!      real      pk_lev_m1(nx,my)
!

!
      integer   k,lpl,lenc,j,nxj,ij,jj,i,n,istat

      logical lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'290'
      end do
      lrec(lpout) = 'h00290'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,rdrag,rdragb,pkout,drag,tens)
!
      lenc= nx*my
      ncnt= -1
!
!!      do jj =1, jlistnum
!!      j=jlist1(jj)
!!      nxj=nxdef_2d(j)
!!      do i=1,nxj
!!        pk_lev_m1(i,j)=pk(i,lev-1,jj)
!!      enddo
!!      enddo
!
!!      call mpe_unify(pk_lev_m1,nx,my,2,mpe_double)
!
      do 20 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
      ncnt=ncnt+1
      do 45 jj=1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 45 i=1,nxj
!      ij=(j-1)*nx+i
      if(pkout(k).gt.pk(i,lev-1,jj)) drag(i,jj,k)= 0.
   45 continue
!
      call unify_reduceintp(nx,my,my_max,drag(1,1,k),wk1)
      call syslbl(lrec(k),idtg,itau,ggdef,ihdg)
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,wk1,istat)
      if(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,196,6,100,-2,plev(k),wk1)
      endif
      call qmaxn3(wk1,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=wk1
        ihdg2=ihdg
      endif
      go to 20
      endif
   10 continue
   20 continue
      if(lwrite .and. myrank .le. ncnt .and.out_pres_form==1) &
                call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
      return
  end subroutine dragout

  subroutine geopout(nx,my,my_max,lpout,lev,itau,ifilout,idtg,plev, &
        num,whtlev,pkout,pk,pklp,phi,phib,phips,glob,ggdef,phistd,      &
        h850,h500,lwrite)

      use index
      use rank
      use mod_grb2_param !for write grb2
      use const ,only:out_pres_form

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,num,ncnt
!
      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max),            &
                phi(nxp,lev,my_max),phib(nxp,my_max),phips(nxp,my_max,lpout),&
                glob(nx,my),plev(lpout),whtlev(num)

      real      tens(lev+1),phistd(lpout),slp(nx,my),h850(nxp,my_max), &
                h500(nxp,my_max),pout(nx,my),tmp(nxp,my_max)
!

      integer*8    idtg

      character*80 ifilout
      character*26 ihdg,ihdg2
      character*6  lrec(lpout)
      character*4  ggdef
!
      integer      i,k,n,lenc,istat,lpl,jj,j,nxj

      logical lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'000'
      end do
      lrec(lpout) = 'h00000'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,phi,phib,pkout,phips,tens)
!
      lenc= nx*my
      ncnt= -1
!
      do k=1,lpout
        if(plev(k).eq.850.)then
        do jj=1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            h850(i,jj)= phips(i,jj,k)
          enddo
        enddo
!byl          call unify_reduceintp(nx,my,my_max,phips(1,1,k),h850)
!!          do i=1,lenc
!!            h850(i,1)=phips(i,k)
!!          enddo
!
!          if(itau.le.72 .and. lreduce.eq.1)then
!byl            call smth9(nx,my,h850,glob,1)
!byl            h850=glob
!byl            call smth9(nx,my,h850,glob,2)
!byl            h850=glob
!          endif
        else if(plev(k).eq.500.)then
        do jj=1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            h500(i,jj)= phips(i,jj,k)
          enddo
        enddo
!byl          call unify_reduceintp(nx,my,my_max,phips(1,1,k),h500)
!!          do i=1,lenc
!!            h500(i,1)=phips(i,k)
!!          enddo
!          if(itau.le.72 .and. lreduce.eq.1)then
!byl            call smth9(nx,my,h500,glob,1)
!byl            h500=glob
!byl            call smth9(nx,my,h500,glob,2)
!byl            h500=glob
!          endif
        endif
      enddo
!
      do 20 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
      ncnt=ncnt+1
      do 11 jj=1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 11 i=1,nxj
       tmp(i,jj)= phips(i,jj,k)+phistd(k)
   11 continue
      call unify_reduceintp(nx,my,my_max,tmp,glob)
!
      call syslbl(lrec(k),idtg,itau,ggdef,ihdg)
!
!      if(lreduce.eq.1 .and. itau.le.72)then
        call smth9(nx,my,glob,slp,1)
        glob=slp
        call smth9(nx,my,glob,slp,2)
!        glob=slp
!      endif
!
      if(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,3,5,0,100,-2,plev(k),slp)
      endif
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3(slp,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=slp
        ihdg2=ihdg
      endif 
!
      go to 20
      endif
   10 continue
   20 continue
      if(lwrite .and. myrank .le. ncnt.and.out_pres_form==1) &
            call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
      return
  end subroutine geopout

  subroutine shumout(nx,my,my_max,lpout,lev,itau,ifilout,idtg   &
      , plev,num,whtlev,pkout,pk,pklp,dpd,dpdb,dew,glob,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use mod_grb2_param
      use const ,only:out_pres_form

      implicit  none
      integer   nx,my,my_max,lpout,lev,itau,num,ncnt

      real      pkout(lpout),pklp(nxp,my_max)                            &
      , pk(nxp,lev,my_max),dpd(nxp,lev,my_max),dpdb(nxp,my_max)          &
      , dew(nxp,my_max,lpout),glob(nx,my),plev(lpout)                    &
      , whtlev(num),tens(lev+1),pout(nx,my),tmp(nxp,my_max)


      integer   i,k,lpl,n,lenc,istat,jj,j,nxj
!
      integer*8 idtg
      character*80 ifilout
      character*26 ihdg,ihdg2
      character*6 lrec(lpout)
      character*4 ggdef
!
      logical :: lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'510'
      end do
      lrec(lpout) = 'h00510'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,dpd,dpdb,pkout,dew,tens)
!
      lenc= nx*my
      ncnt= -1
!
      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
! relative humidity must be smaller or equal 1.0
!
      ncnt=ncnt+1
      do 20 jj=1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 20 i=1,nxj
       tmp(i,jj)= min(100.,max(dew(i,jj,k)*100.,0.0))
   20 continue
      call unify_reduceintp(nx,my,my_max,tmp,glob)
      if(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,1,2,100,-2,plev(k),glob)
      endif

!!      do 20 i=1,lenc
!!      glob(i,1)= min(100.,max(dew(i,k)*100.,0.0))
!!   20 continue
!
      call syslbl(lrec(k),idtg,itau,ggdef,ihdg)
!
!  reduceintp has been done in voterp (2011/5)
!
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3(glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=glob
        ihdg2=ihdg
      endif
      go to 30
      endif
   10 continue
   30 continue
      if(lwrite .and. myrank .le. ncnt.and.out_pres_form==1) &
          call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
      return
  end subroutine shumout

  subroutine shumout2(nx,my,my_max,lpout,lev,itau,ifilout,idtg   &
      , plev,num,whtlev,pkout,pk,pklp,dpd,dpdb,dew,glob,ggdef,ntrac,lwrite)
!
      use index
      use rank, only : myrank
      use radn, only : ntoz
      use param, only : ncld
      use mod_grb2_param
      use const ,only:out_pres_form

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,num,ntrac,ncnt,ntrchk

      real      pkout(lpout),pklp(nxp,my_max)                        &
      , pk(nxp,lev,my_max),dpd(nxp,lev,my_max),dpdb(nxp,my_max)      &
      , dew(nxp,my_max,lpout),glob(nx,my),plev(lpout)                &
      , whtlev(num),tens(lev+1),pout(nx,my),tmp(nxp,my_max)
!

      integer   i,k,lpl,n,lenc,istat,jj,j,nxj

      integer*8 idtg
      character*80 ifilout
      character*26 ihdg,ihdg2
      character*6 lrec(lpout)
      character*4 ggdef
      character*3 cspec(6)
      logical :: lwrite
!
      cspec=(/'500','551','553','552','554','555'/)
!
      ntrchk=ncld
      if ( ntoz .gt. 0 ) ntrchk=ncld-1
!
      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      if(ntrac.le.ntrchk)then ! all hydrometeors

      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,cspec(ntrac)
      end do
      write( lrec(lpout), '(a3,a3)' ) 'h00',cspec(ntrac)
      Ptp0=0 ;Ptp1=1 ;Ptp2=0 ;Ptp3=6 !grib code
!
      else if(ntrac.eq.ntoz)then
!
        do k = 1, lpout-1
          lpl = int(plev(k)+0.001)
          write( lrec(k), '(i3.3,a3)' ) lpl,'560'   ! ozone
        end do
        lrec(lpout) = 'h00560'
        Ptp0=0 ;Ptp1=14 ;Ptp2=1 ;Ptp3=8 !grib code
!
      else if(ntrac.eq.ncld+1)then
!
        do k = 1, lpout-1
          lpl = int(plev(k)+0.001)
          write( lrec(k), '(i3.3,a3)' ) lpl,'550'   ! combine cloud water and cloud ice together
        end do
        lrec(lpout) = 'h00550'
        Ptp0=0 ;Ptp1=1 ;Ptp2=22 ;Ptp3=8 !grib code 
      else
        goto 40
      endif
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,dpd,dpdb,pkout,dew,tens)
!
      lenc= nx*my
      ncnt= -1
!
      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
!
! relative humidity must be smaller or equal 1.0
!
      ncnt=ncnt+1
      do 20 jj=1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 20 i=1,nxj
       tmp(i,jj)= max(dew(i,jj,k),0.0)
   20 continue
      call unify_reduceintp(nx,my,my_max,tmp,glob)
!!      do 20 i=1,lenc
!!      glob(i,1)= max(dew(i,k),0.0)
!!   20 continue
      if(out_pres_form==2.and.myrank==0)then
          call wrt_grb2(itau,Ptp0,Ptp1,Ptp2,Ptp3,100,-2,plev(k),glob)   !Specit Humility
      endif
!
      call syslbl(lrec(k),idtg,itau,ggdef,ihdg)
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3(glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=glob
        ihdg2=ihdg
      endif
      go to 30
      endif
   10 continue
   30 continue
      if(lwrite .and. myrank .le. ncnt .and. out_pres_form==1) &
           call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
   40 continue

      return
  end subroutine shumout2


  subroutine surfout(nx,my,my_max,ifilout,itau,idtg,taudir,ntau,pdiff  &
       ,pt,ptop,slp,ptend,glob,ggdef,lwrite)
!
      use index
      use mpe
      use rank
      use mod_grb2_param
      use const ,only:out_pres_form
!
      implicit  none
      integer   nx,my,my_max,i,j,jj,kk,n,lev,nxj,itau,ntau,num,lenc,istat

      real      pdiff(nxp,my_max),pt(nxp,my_max),ptend(nxp,my_max),slp(nxp,my_max),glob(nx,my)
      real      tmp(nxp,my_max)
      character*16 taudir(ntau)
      character*4 ggdef
!
      real      ptop,tnshun
!
      integer*8 idtg
      character*80 ifilout
      character*26 lrec
      character*6 label(ntau),labx
   
      logical :: lwrite
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
!byl          slp(i,j)= pt(i,jj)+pdiff(i,jj)
          slp(i,jj)= pt(i,jj)+pdiff(i,jj)
          ptend(i,jj)= ptend(i,jj)*3600.0
        enddo
      enddo
!
!byl      call unify_reduceintp(nx,my,my_max,tmp,slp)
!byl      call mpe_unify(slp,nx,my,2,mpe_double)
!byl      if( lreduce.eq.1 ) call reduceintp (slp,nxdef,nx,my)
!byl      call smth9(nx,my,slp,glob,1)
!byl      slp=glob
!byl      call smth9(nx,my,slp,glob,2)
!byl      slp=glob
!
      num= 0
      do 20 n=1,ntau
        read(taudir(n),'(a6,1x,i4)') labx,lev
        if(lev.eq.0) then
          num= num+1
          label(num)= labx
        endif
   20 continue
      if(num.eq.0) return
!
      tnshun= 1.0
      lenc= nx*my
!
      do 100 kk=1,num
!
!  sea surface level pressure
!

      if(label(kk).eq.'SSL010' .or. label(kk).eq.'ssl010') then
        call unify_reduceintp(nx,my,my_max,slp,glob)
        call syslbl('ssl010',idtg,itau,ggdef,lrec)
        if(out_pres_form==1)then
          if(lwrite) call dmswrit(nx,my,lrec,lenc,'H',ifilout,glob,istat)
        elseif(out_pres_form==2.and.myrank==0)then
          call wrt_grb2(itau,0,3,0,2,101,0,0.,glob)
        endif
        call qmaxn3(glob,lrec(1:14),lrec(15:26),1,1,1,nx,my,1)

!
!  terrain pressure
!
      else if(label(kk).eq.'B00010' .or. label(kk).eq.'b00010') then
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
!byl            glob(i,j) = pt(i,jj) + ptop
            tmp(i,jj) = pt(i,jj) + ptop
          enddo
        enddo
        call unify_reduceintp(nx,my,my_max,tmp,glob)
!byl        call mpe_unify(glob,nx,my,2,mpe_double)
        call syslbl('b00010',idtg,itau,ggdef,lrec)
!byl        if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
        if(lwrite) call dmswrit(nx,my,lrec,lenc,'H',ifilout,glob,istat)
        !write grb2 data 
        if(out_pres_form==2.and.myrank==0)then
          call wrt_grb2(itau,0,3,0,2,1,0,0.,glob)
        endif
        call qmaxn3(glob,lrec(1:14),lrec(15:26),1,1,1,nx,my,1)

!
!  terrain pressure tendency
!
!      else if(label(kk).eq.'B00011' .or. label(kk).eq.'b00011') then
!        do i = 1, lenc
!          glob(i,1) = ptend(i,1)
!        end do
!        call syslbl('b00011',idtg,itau,ggdef,lrec)
!        if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!        if(lwrite) call dmswrit(nx,my,lrec,lenc,'H',ifilout,glob,istat)
!        call qmaxn3(glob,lrec(1:14),lrec(15:26),1,1,1,nx,my,1)
!
      endif
!
  100 continue
      return
  end subroutine surfout

  subroutine tempout(nx,my,my_max,lpout,lev,itau,ifilout,idtg  &
      , plev,num,whtlev,pkout,pk,pklp,tt,ttbot,temp,ggdef,lwrite)
!
      use index
      use rank,   only : myrank
      use mod_grb2_param
      use const ,only:out_pres_form

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,k,i,n,istat
      integer   num,lpl,lenc,ncnt
      real      tnshun

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max)        &
      , tt(nxp,lev,my_max),ttbot(nxp,my_max),temp(nxp,my_max,lpout)     &
      , plev(lpout),whtlev(num)
      real      tens(lev+1),glob(nx,my),slp(nx,my),pout(nx,my)

!
      integer*8 idtg
      character*80 ifilout
      character*26 ihdg
      character*26 ihdg2
      character*6 lrec(lpout)
      character*4 ggdef
!
     logical :: lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev+1)= 0.0
      tens(lev)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'100'
      end do
      lrec(lpout) = 'h00100'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,tt,ttbot,pkout,temp,tens)
!
      lenc= nx*my
!
      tnshun= 1.0
!
      ncnt=-1
!
      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then  
!
      ncnt=ncnt+1
!
      call unify_reduceintp(nx,my,my_max,temp(1,1,k),glob)
!!      do 11 i=1,lenc
!!      glob(i,1)= temp(i,k)
!!   11 continue
!
      call syslbl(lrec(k),idtg,itau,ggdef,ihdg)
!
!  reduceintp has been done in voterp (2011/5)
!
        call smth9(nx,my,glob,slp,1)
        glob=slp
        call smth9(nx,my,glob,slp,2)
!!        glob=slp
!
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      if(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,0,0,2,100,-2,plev(k),slp)
      endif

      call qmaxn3(slp,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=slp
        ihdg2=ihdg
      endif
      go to 30
      endif
   10 continue
   30 continue
      if(lwrite .and. myrank .le. ncnt .and. out_pres_form==1) &
          call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
      return
  end subroutine tempout

  subroutine vortout(nx,my,my_max,lpout,lev,itau,ifilout,idtg     &
      , plev,num,whtlev,pkout,pk,pklp,rvor,rvorb,vor,ggdef,v850,v700,lwrite)
!
      use index
      use rank, only : myrank
      use mod_grb2_param
      use const ,only:out_pres_form

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,num,ncnt

      real      pkout(lpout),pklp(nxp,my_max)                     &
      , pk(nxp,lev,my_max),rvor(nxp,lev,my_max),rvorb(nxp,my_max) &
      , vor(nxp,my_max,lpout),plev(lpout),whtlev(num)
      real      tens(lev+1)
      real      v850(nxp,my_max),v700(nxp,my_max)
      real      wk1(nx,my),pout(nx,my)

!
      integer*8 idtg
      character*80 ifilout
      character*26 ihdg,ihdg2
      character*6 lrec(lpout)
      character*4 ggdef

      integer   k,lpl,lenc,i,n,istat,jj,j,nxj
!
      logical :: lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'240'
      end do
      lrec(lpout) = 'h00240'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,rvor,rvorb,pkout,vor,tens)
!
      lenc= nx*my
      ncnt= -1
!
      do k=1,lpout
        if(plev(k).eq.850.)then
        do jj=1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            v850(i,jj)= vor(i,jj,k)
          enddo
        enddo
!byl        call unify_reduceintp(nx,my,my_max,vor(1,1,k),v850)
!!          do i=1,lenc
!!            v850(i,1)=vor(i,k)
!!          enddo
!
!  reduceintp has been done in voterp (2011/5)
!
        else if(plev(k).eq.700.)then
        do jj=1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            v700(i,jj)= vor(i,jj,k)
          enddo
        enddo
!byl        call unify_reduceintp(nx,my,my_max,vor(1,1,k),v700)
!!          do i=1,lenc
!!            v700(i,1)=vor(i,k)
!!          enddo
!
!  reduceintp has been done in voterp (2011/5)
!
        endif
      enddo
!
      do 20 n=1,num
      do 10 k=1,lpout
      if(plev(k).eq.whtlev(n)) then
      ncnt=ncnt+1
      call unify_reduceintp(nx,my,my_max,vor(1,1,k),wk1)
      call syslbl(lrec(k),idtg,itau,ggdef,ihdg)
!
!  reduceintp has been done in voterp (2011/5)
!
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,wk1,istat)
      if(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,12,6,100,-2,plev(k),wk1)
      endif
      call qmaxn3(wk1,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=wk1
        ihdg2=ihdg
      endif
      go to 20
      endif
   10 continue
   20 continue
      if(lwrite .and. myrank .le. ncnt .and.out_pres_form==1) & 
            call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
      return
  end subroutine vortout

  subroutine windout(nx,my,my_max,lpout,lev,itau,ifilout,idtg      &
      , plev,num,whtlev,cosl,pkout,pk,pklp,ut,vt,sdhat,utb,vtb     &
      , wind,glob,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use mod_grb2_param
      use const ,only:out_pres_form

      implicit none

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max)          &
      , rdiv(nxp,lev,my_max),ut(nxp,lev,my_max),vt(nxp,lev,my_max)        &
      , utb(nxp,my_max),vtb(nxp,my_max),wind(nxp,my_max,lpout),cosl(my)   &
      , glob(nx,my),plev(lpout),whtlev(num),sdhat(nxp,lev,my_max)

      real      tens(lev+1),wtb(nxp,my_max),pout(nx,my),tmp(nxp,my_max)
!
      integer   nx,my,my_max,lpout,lev,itau,jj,nxj,ncnt
      integer   num,k,lenc,lpl,n,i,j,istat

      real      rad,xxx
!lzl +add
!!      real      pklzl(nx,my)
!lzl -end

      integer*8 idtg
      character*80 ifilout
      character*26 ihdg,ihdg2
      character*6 lrec(lpout),krec(lpout),mrec(lpout)
      character*4 ggdef
!
      data rad/6.371e6/
!
     logical :: lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
      lenc= nx*my
      ncnt= -1
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'200'
       write( krec(k), '(i3.3,a3)' ) lpl,'210'
       write( mrec(k), '(i3.3,a3)' ) lpl,'220'
      end do
      lrec(lpout) = 'h00200'
      krec(lpout) = 'h00210'
      mrec(lpout) = 'h00220'
!
! first: do the u components
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,ut,utb,pkout,wind,tens)
!
!!      if( lreduce.eq.1 ) call reduceintp (utb,nxdef,nx,my)
!
!lzl +add======================================================
!!      do i =1,nx
!!      do j =1,my
!!         pklzl(i,j)=pklp(i,j)
!!      end do
!!      end do
!
!!      if( lreduce.eq.1 ) then
!!          call reduceintp(pklzl,nxdef,nx,my)
!!      endif

!lzl -end=====================================================          

      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
!  below ground level extrapolate surface wind downward
!  deweight wind with cos latitude, earth radius
!
        ncnt=ncnt+1
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xxx= rad/cosl(j)
          do i=1,nxj
           if(pkout(k).gt.pklp(i,jj)) wind(i,jj,k)= utb(i,jj)
           tmp(i,jj)=wind(i,jj,k)*xxx
          enddo
        enddo
        call unify_reduceintp(nx,my,my_max,tmp,glob)
!!      do 45 i=1,nx*my
!lzl c      if(pkout(k).gt.pklp(i,1)) wind(i,1,k)= utb(i,1)
!!      if(pkout(k).gt.pklzl(i,1)) wind(i,1,k)= utb(i,1)  !lzl use full grid(pklzl)
!!   45 continue
!
!
!  deweight wind with cos latitude, earth radius
!
!!      do 50 j=1,my
!!      xxx= rad/cosl(j)
!!      do 50 i=1,nx
!!      glob(i,j)= wind(i,j,k)*xxx
!!   50 continue
!
      if(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,2,2,100,-2,plev(k),glob)
      endif

      call syslbl(lrec(k),idtg,itau,ggdef,ihdg)
!
!  reduceintp has been done in voterp (2011/5)
!
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3(glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=glob
        ihdg2=ihdg
      endif
      go to 30
      endif
   10 continue
   30 continue
      if(lwrite .and. myrank .le. ncnt .and.out_pres_form==1) &
           call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
!  now the v components
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,vt,vtb,pkout,wind,tens)
!
!!      if( lreduce.eq.1 ) call reduceintp (vtb,nxdef,nx,my)
!
!lzl +add===============================================================
!!      do i =1,nx
!!      do j =1,my
!!         pklzl(i,j)=pklp(i,j)
!!      end do
!!      end do
!
!!      if( lreduce.eq.1 ) then
!!          call reduceintp(pklzl,nxdef,nx,my)
!!      endif

!lzl -end================================================================
!
      ncnt= -1
      do 40 n=1,num
      do 20 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
!  below ground level extrapolate surface wind downward
!  deweight wind with cos latitude, earth radius
!
        ncnt=ncnt+1
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xxx= rad/cosl(j)
          do i=1,nxj
           if(pkout(k).gt.pklp(i,jj)) wind(i,jj,k)= vtb(i,jj)
           tmp(i,jj)=wind(i,jj,k)*xxx
          enddo
        enddo
        call unify_reduceintp(nx,my,my_max,tmp,glob)
!
!!      do 55 i=1,nx*my
!lzl c      if(pkout(k).gt.pklp(i,1)) wind(i,1,k)= vtb(i,1)
!!      if(pkout(k).gt.pklzl(i,1)) wind(i,1,k)= vtb(i,1) !lzl use full grid (pklzl)
!!   55 continue
!
!  deweight wind with cos latitude, earth radius
!
!!      do 60 j=1,my
!!      xxx= rad/cosl(j)
!!      do 60 i=1,nx
!!      glob(i,j)= wind(i,j,k)*xxx
!!   60 continue
!
      if(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,3,2,100,-2,plev(k),glob)
      endif

      call syslbl(krec(k),idtg,itau,ggdef,ihdg)
!
!  reduceintp has been done in voterp (2011/5)
!
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3(glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=glob
        ihdg2=ihdg
      endif
      go to 40
      endif
   20 continue
   40 continue
      if(lwrite .and. myrank .le. ncnt .and.out_pres_form==1) & 
          call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
!  now the w components
!
      wtb=0.
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,sdhat,wtb,pkout,wind,tens)
!
      ncnt= -1
      do 42 n=1,num
      do 22 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
      ncnt=ncnt+1
      call unify_reduceintp(nx,my,my_max,wind(1,1,k),glob)
!!      do j=1,my
!!      do i=1,nx
!!        glob(i,j)= wind(i,j,k)
!!      enddo
!!      enddo
!
      if(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,8,6,100,-2,plev(k),glob)
      endif

      call syslbl(mrec(k),idtg,itau,ggdef,ihdg)
!
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3(glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if ( myrank .eq. ncnt ) then
        pout=glob*100.
        ihdg2=ihdg
      endif
      go to 42
      endif
   22 continue
   42 continue
      if(lwrite .and. myrank .le. ncnt .and.out_pres_form==1) &
            call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,pout,istat)
!
      return
  end subroutine windout
!
end module mod_outflds
