      subroutine sigful(nx,my,my_max,lev,ncld,lmax,jtrun,jtmax            &
              , cstar,ktrop,idtg,ptop,taux,capa,grav,rgas,rad,cp          &
              , weight,poly,sigma,cosl,phi,tt,ut,vt,sht,o3l,pt,sgeo,pdiff &
!             , weight,poly,sigma,cosl,phi,tt,ut,vt,sht,pt,sgeo,pdiff
              , tsave,t1000,plt,pk,pk2,taup,ggdef,gmdef)
!
!  read and interpolate analysis fields to model's coordinates
!
!  ********************************************************************
!
! **** input ****
!
!  sgeo: gaussian grid terrain geopotential
!
! **** output ****
!
!  phi: 3-d geopotential on gaussian grid and sigma coord.
!  tt: 3-d temperature on gaussian grid and sigma coord.
!  ut: 3-d e-w wind component on gaussian grid and sigma coord.
!  vt: 3-d n-s wind component on gaussian grid and sigma coord.
!  sht: 3-d specific humidity on gaussian grid and sigma coord.
!  pt: 2-d terrain pressure on gaussian grid.
!  pdiff: 2-d difference between initial SLP and terrain pressure.
!  plt: 3-d full level pressure on gaussian grid and sigma coord.
!  pk: 3-d full level exner func on gaussian grid and sigma coord.
!  pk2: 3-d half level exner func on gaussian grid and sigma coord.
!
      use const, only : RTYPE,kflag,qmin,ifilin,keyi,ihdgi
      use mpe
      use rank
      use index
      use radn, only : ntoz,ntcw,ntiw
      use spec, only : trefs

      implicit  none

      integer   nx,my,my_max,lev,ncld,lmax,jtrun,jtmax,KL
      integer   ktrop,nxmy,nxlev,lncrec,lmaxp1,lmaxp2,k,itaux,itaup
      integer   istat,i,ii,jj,j,nxj,kk,lqwset,m,mf,n,llts,ntrac,nclds

      real      taup,cp,rad,rgas,grav,capa,taux,ptop,ppp,fac,ptmp
      real      alaps,rdg,ttt1,ttt2,apha,ttt,sigp,x1,opok,pk800,pk300


      logical cstar
      real      pdiff(nxp,my_max),t1000(nxp,my_max)                     &
              , tsave(nxp,my_max),plt(nxp,lev,my_max)
      real(kind=RTYPE) weight(my),poly(jtrun,jtmax,my/2)                &
              ,        cosl(my),sigma(lev+1,2),pk(nxp,lev,my_max)       &
              ,        pk2(nxp,lev,my_max)
      character*4 ggdef,gmdef
!
!  local work arrays
!
      real      hld1(nx,my),hkd1(nx,lev)
      real     tens(lmax+2),tstd(lmax),utmp(nxp,lev),vtmp(nxp,lev)
      real      puvphi(26),plog(nx,lev),preplt(nx,lmax+2)
!
      real(kind=RTYPE) cc(nx+2,levp,1,my_max)                          &
               ,       ut(nxp,lev,my_max),vt(nxp,lev,my_max)           &
               ,       tt(nxp,lev,my_max),sht(nxp,lev*ncld,my_max)     &
               ,       o3l(nxp,lev,my_max),phi(nxp,lev,my_max)         &
               ,       pt(nxp,my_max),sgeo(nxp,my_max)                 &
               ,       anlslp(nxp,my_max),prett(nx,lmax+2)             &
               ,       ut_tmp(nx,lev)
      real(kind=RTYPE) hld4(nx,levp,ncld,my_max),hld3(nx,levp,my_max)  &
               ,       hld2(nx,my)
      real      work_pr1(lev), work_pr2(lev), work_pr3(lev)
!
      real(kind=RTYPE) plnow(jtrun,jtmax,2),dummy,ww1(nx,my_max)
!
      character*6 typ
      character*3 cspec(6)
      integer      inistat
!dms34
      integer*8 idtg,idtg2
!
!      data puvphi/10.,20.,30.,50.,70., 100.,150.,200.,250.,300.,400.
!     * ,500.,700., 850.,925.,1000./
      data puvphi/10.0,20.0,30.0,50.0,70.0,100.0,150.0,200.0,250.0   &
               ,300.0,350.0,400.0,450.0,500.0,550.0,600.0,650.0      &
               ,700.0,750.0,800.0,850.0,900.0,925.0,950.0,975.0      &
               ,1000.0/
      data cspec/'500','551','553','552','554','555'/
!
      if ( ntoz .gt. 0 ) then
        nclds=ncld-1
      else
        nclds=ncld
      endif

      !new Thompson MP without reading ice/rain number concentration
      if ( ncld .eq. 9 ) then
        nclds = 6
      !Goddard MP without reading hail
      elseif ( ncld .eq. 8 ) then
        nclds = 6
      endif
!
!CWBinit
      plnow=0.
      preplt=0.
      plog=0.
      prett=0.
      anlslp=0.
      sht=0.
      dummy=0.
      cc=0.

      nxmy  = nx*my
      nxlev = nx*lev
      lncrec= nxmy
      lmaxp1 = lmax+1
      lmaxp2 = lmax+2
!
      do 5 k = 1, lmaxp2
       tens(k) = 1.0
    5 continue
      tens(1) = 0.
      tens(lmaxp2) = 0.
!
      itaux = taux + 0.001
      itaup = taup + 0.001
!
      call dtgfix12(idtg,idtg2,-itaup)
!
!  read in pt
!
      call syslbl_r ('b00010',idtg,itaux,ggdef)
      call dmsread (nx,my,lncrec,'H',ifilin,hld1,istat)
!byl      if( lreduce.eq.1 )call reducepick (hld1,nxdef,nx,my)
      do jj = 1, jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
        do i=1,nxj
          pt(i,jj) = hld1(ii,j) - ptop
          ii=ii+1
        enddo
      enddo
      if( lreduce.eq.1 )then
!ch     call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,pt,plnow,nsizey)
        call mpe2d_unify_nx(ww1,pt)
        call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,ww1,plnow,nsizey)
        call transr1(jtrun,jtmax,nx,my,my_max,poly,plnow,pt,nsizey)
      endif
!
!  read in temp at sigma levels
!
      do 71 k = 1, levp
        KL=lev-Llist(k)+1
      if ( KL .lt. 100 ) then
        write (typ, '("m",i2.2,"100")' ) KL
      else
        write (typ, '("n",i2.2,"100")' ) mod(KL,100)
      endif
      call syslbl_r (typ,idtg,itaux,gmdef)
      call dmsread_split(nx,my,lncrec,'H',ifilin,hld1,istat)
      do 71 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef(j)
       if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
       do 71 i = 1, nxj
         hld3(i,k,jj) = hld1(i,j)
  71  continue
      call mpe2d_transpose_ndsl_f2p(hld3,tt, &
            nxp,nx,levf,levp,1,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!
!  read in q at sigma levels
!
      hld4=qmin
      do 73 k = 1, levp
        KL=lev-Llist(k)+1
      if ( KL .lt. 100 ) then
        write (typ, '("m",i2.2,"500")' ) KL
      else
        write (typ, '("n",i2.2,"500")' ) mod(KL,100)
      endif
      call syslbl_r (typ,idtg,itaux,gmdef)
      call dmsread_split (nx,my,lncrec,'H',ifilin,hld1,istat)
      do 73 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef(j)
       if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
       do 73 i = 1, nxj
         if(hld1(i,j).le.1.0e-10) hld1(i,j)=1.0e-10
         hld4(i,k,1,jj) = hld1(i,j)
  73  continue
!
      if( ncld .ge. 2 ) then
!  check initial data of all hydrometeors
!
        if ( nclds .gt. 2 ) then
          inistat=0

          if(col_rank .eq. 0) then
            KL=lev-Llist(1)+1
            do ntrac=2,nclds
              if ( KL .lt. 100 ) then
                write (typ, '("m",i2.2,a3)' ) KL,cspec(ntrac)
              else
                write (typ, '("n",i2.2,a3)' ) mod(KL,100),cspec(ntrac)
              endif
              call syslbl_r (typ,idtg2,itaup,gmdef)
#ifdef I38K
              write(keyi,'(a28,a1,i9.9)') ihdgi,'H',lncrec
#else
              write(keyi,'(a26,a1,i7.7)') ihdgi,'H',lncrec
#endif
              call dmschkr (ifilin,keyi//char(0),istat)
              inistat=inistat+istat
            enddo
          endif

          call mpe_global_sum(inistat,1,mpe_integer)
          if ( myrank .eq. 0 .and. inistat .gt. 0 ) then 
            print *,'========== Warning!!! ==========='
            print *,'not enough initial data for all hydrometeors!!'
          endif
        else
          inistat=1
        endif
!
!  get first guest as initial
!
        if ( inistat .gt. 0) then
          ntrac=2
          do k = 1, levp
            KL=lev-Llist(k)+1
            if ( KL .lt. 100 ) then
              write (typ, '("m",i2.2,"550")' ) KL    ! cloud liquid water content
            else
              write (typ, '("n",i2.2,"550")' ) mod(KL,100)
            endif
            call syslbl_r (typ,idtg2,itaup,gmdef)
            call dmsread_split (nx,my,lncrec,'H',ifilin,hld1,istat)
!
            do jj = 1, jlistnum
              j=jlist1(jj)
              nxj=nxdef(j)
              if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
              do i = 1, nxj
!  no cloud ice data, simple way to split cloud water and ice temporary
                if ( hld3(i,k,jj)-273.15 .le. -15. .and. nclds .ge. 3 ) then
                  ntrac = ntiw
                else
                  ntrac = ntcw
                endif
                hld4(i,k,ntrac,jj) = max(hld1(i,j),qmin)
              end do
            end do
          end do
!
        else
!
          do ntrac=2,nclds
            do k = 1, levp
              KL=lev-Llist(k)+1
            if ( KL .lt. 100 ) then
              write (typ, '("m",i2.2,a3)' ) KL,cspec(ntrac)    ! cloud liquid water content
            else
              write (typ, '("n",i2.2,a3)' ) mod(KL,100),cspec(ntrac)    ! cloud liquid water content
            endif
              call syslbl_r (typ,idtg2,itaup,gmdef)
              call dmsread_split (nx,my,lncrec,'H',ifilin,hld1,istat)
!
              do jj = 1, jlistnum
                j=jlist1(jj)
                nxj=nxdef(j)
                if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
                do i = 1, nxj
                  hld4(i,k,ntrac,jj) = max(hld1(i,j),qmin)
                enddo
              enddo
            enddo
          enddo
!
        endif
!
!  read "observed ozone" at sigma levels for doing ozone forecast
!
      if(ncld.eq.ntoz)then
      ntrac=ntoz
      do k = 1, levp
        KL=lev-Llist(k)+1
        if ( KL .lt. 100 ) then
          write (typ, '("m",i2.2,"560")' ) KL
        else
          write (typ, '("n",i2.2,"560")' ) mod(KL,100)
        endif
        call syslbl_r (typ,idtg,itaux,gmdef)
        call dmsread_split (nx,my,lncrec,'H',ifilin,hld1,istat)
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef(j)
          if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
          do i = 1, nxj
             hld4(i,k,ntrac,jj) = hld1(i,j)
             hld3(i,k,jj) = hld1(i,j)
          enddo
        enddo
      enddo
      endif
!
      endif    ! end of if(ncld.ge.2)
      call mpe2d_transpose_ndsl_f2p(hld4,sht, &
            nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(hld3,o3l, &
            nxp,nx,levf,levp,1,myf,my_max,jlistnum,jlen,nsizex,row_comm)
      
!----
      if( lreduce.eq.1 )then
      call joinrs(cc,tt,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc,trefs,1,nsizey)
      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,trefs,cc,1,nsizey)
      call ujoinsr(cc,tt,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      endif

!------
!
!  ncld > 2 needs to add another cloud micro input
!  "incrini.f" also needs to take care
!
      if( ncld .gt. ntoz ) then
       if(myrank.eq.0)print*,'ncld > 3, incomplete input for cloud micro'
       call mpe_finalize
       call dmsexit(-1)
      end if
!
!  obtain plt, pk, pk2 and phi
!
      do 80 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
!
        call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,pt(1,jj) &
                        ,pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj) )
!
!  derive phi from t and q based on the hydrostatic relation
!
      do 75 i = 1, nxj
       phi(i,lev,jj)= cp*tt(i,lev,jj)*(1.0+0.608*sht(i,lev,jj)) &
                    *(pk2(i,lev,jj)/pk(i,lev,jj)-1.0)+sgeo(i,jj)
   75 continue
!
      do 76 k = lev-1, 1, -1
      do 76 i = 1, nxj
       phi(i,k,jj)= phi(i,k+1,jj) + cp*( tt(i,k,jj)               &
            *(1.0+0.608*sht(i,k,jj))*(pk2(i,k,jj)/pk(i,k,jj)-1.0) &
            +tt(i,k+1,jj)*(1.0+0.608*sht(i,k+1,jj))               &
            *(1.0-pk2(i,k,jj)/pk(i,k+1,jj)) )
   76 continue
!
   80 continue
!
!  compute reference temp and q for horizontal diffusion
!  compute stardard atmospheric temp on pressure and then interpolate it
!  to sigma level and then compute saturated moisture for reference q
!
      call ttstd(lmax, puvphi, tstd)
!
      do 140 k = 1, lmax
      do 140 i = 1, nx
       preplt(i,k+1) = log(puvphi(k))
       prett(i,k+1) = tstd(k)
  140 continue
!
      do 142 i = 1, nx
       preplt(i,1) = log( max(ptop,1.0) )
       prett(i,1) = prett(i,2) + (prett(i,2)-prett(i,3))          &
          * (preplt(i,1)-preplt(i,2))/(preplt(i,2)-preplt(i,3))
  142 continue
!c

!ch> 
      call mpe2d_unify_nx(ww1,pt)
!ch<

      do 170 jj = 1, jlistnum
       j=jlist1(jj)
       ii=nxjstart(j)
       nxj=nxdef_2d(j)
!
      do 145 i = 1, nxdef(j)
!ch    ppp = log( pt(i,jj) + ptop )
       ppp = log( ww1(i,jj) + ptop )
       if( preplt(i,lmaxp1) .ge. ppp ) then
          preplt(i,lmaxp2) = 2.0*preplt(i,lmaxp1) - preplt(i,lmax)
       else
          preplt(i,lmaxp2) = ppp
       end if
  145 continue
!
      do 155 i = 1, nxdef(j)
       prett(i,lmaxp2) =                                            &
               prett(i,lmaxp1) + (prett(i,lmaxp1)-prett(i,lmax)) *  &
               (preplt(i,lmaxp2)-preplt(i,lmaxp1))/                 &
               (preplt(i,lmaxp1)-preplt(i,lmax))
  155 continue
!
      do 157 k = 1, lev
      do 157 i = 1, nxj
       plog(i,k) = log(plt(i,k,jj))
  157 continue
!
      call mpe2d_unify_nx_lev_red(plog,nx,lev,j)
!ch   call vterpj( nx,lmaxp2,lev,preplt,prett,plog,ut(1,1,jj),tens)
      call vterpj( nx,lmaxp2,lev,preplt,prett,plog,ut_tmp,tens)
        do  k = 1, lev
            n=ii
        do  i = 1,nxj
!ch         ut(i,k,jj)=ut_tmp(i,k) 
!cjh        ut(i,k,jj)=ut_tmp(n,k) 
            utmp(i,k)=ut_tmp(n,k)
            ut(i,k,jj) = ut_tmp(n,k)/pk(i,k,jj)
            n=n+1
        enddo
        enddo

!
!      call qsatq_2d( nxjp(j),nxp,lev,utmp,plt(1,1,jj),vtmp)
!
!      do 160 k = 1, lev
!      do 160 i = 1, nxj
!       vt(i,k,jj) = vtmp(i,k)
!  160 continue
!
  170 continue
!
!      do jj =1, jlistnum
!        j=jlist1(jj)
!      if(j .eq. 85) then
!        do k=1,lev
!          work_pr1(k)=plt(58,k,jj)
!          work_pr2(k)=ut (58,k,jj)
!          work_pr3(k)=vt (58,k,jj)
!        enddo
!      call mpe_send_print(work_pr1, lev, 85, mpe_double)
!      call mpe_send_print(work_pr2, lev, 85, mpe_double)
!      call mpe_send_print(work_pr3, lev, 85, mpe_double)
!      endif
!      enddo
!
!      if(myrank .eq. 0) then
!        call mpe_recv_print(work_pr1, lev, 85, mpe_double)
!        call mpe_recv_print(work_pr2, lev, 85, mpe_double)
!        call mpe_recv_print(work_pr3, lev, 85, mpe_double)
!
!        do 168 k=1,lev
!          print 900, k,work_pr1(k),work_pr2(k),work_pr3(k)
!900       format(1x,i3,1x,f6.1,1x,f6.2,1x,e10.3)
!168     continue
!      endif
!
!     do m=1, mlistnum
!       mf=mlist(m)
!     do n=mf, jtrun
!     do k=1, lev*ncld*2
!        qrefs(k,1,n,m)=0.0
!     enddo
!     enddo
!     enddo
!
!!    qrefs=0.0
      call joinrs(cc,ut,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc,trefs,1,nsizey)

!
!!      do m=1, mlistnum
!!        mf=mlist(m)
!!      do n=mf, jtrun
!!        do k=1, levp
!!           qrefs(k,1,n,m)=qtmps(k,1,n,m)
!!           qrefs(k,2,n,m)=qtmps(k,2,n,m)
!!        end do
!!      end do
!!      end do
!
!  compute sea level pressure and write out pdiff
!  The method is based on one used by ecmwf, reseach manual 2 (1988)
!
      alaps = 0.0065
      rdg = rgas/grav
!
!  llts layer's temperature is used to derive an alternative
!  surface skin temperature
!
      llts = lev-5
!
      do 190 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 180 i = 1, nxj
       ttt1 = tt(i,lev,jj) + alaps*rdg*tt(i,lev,jj)*       &
              ((pt(i,jj)+ptop)/plt(i,lev,jj)-1.0)
       ttt2 = tt(i,llts,jj) + alaps*(phi(i,llts,jj)-sgeo(i,jj))/grav
       hld1(i,j) = 0.25*ttt1 + 0.75*ttt2
       hld2(i,j) = hld1(i,j) + alaps*sgeo(i,jj)/grav
  180 continue
      ii=nxjstart(j)
!
      do 185 i = 1, nxj
      if( sgeo(i,jj) .lt. 0.1 ) then
       anlslp(i,jj) = pt(i,jj) + ptop

      elseif( hld1(i,j) .le. 290.5 .and. hld2(i,j) .gt. 290.5 ) then
       apha = rgas*(290.5-hld1(i,j))/sgeo(i,jj)
       ttt = sgeo(i,jj)/(rgas*hld1(i,j))
       anlslp(i,jj) = (pt(i,jj)+ptop)*exp( ttt*(1.0-0.5*apha*ttt+0.333333* &
                     apha*ttt*apha*ttt) )

      elseif( hld1(i,j) .gt. 290.5 .and. hld2(i,j) .gt. 290.5 ) then
       hld1(i,j) = (hld1(i,j)+290.5)*0.5
       anlslp(i,jj) = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,j)) )

      elseif( hld1(i,j) .lt. 255.0 .and. hld2(i,j) .lt. 255.0 ) then
       hld1(i,j) = (hld1(i,j)+255.0)*0.5
       anlslp(i,jj) = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,j)) )

      else
       apha = alaps * rdg
       ttt = sgeo(i,jj)/(rgas*hld1(i,j))
       anlslp(i,jj) = (pt(i,jj)+ptop)*exp( ttt*(1.0-0.5*apha*ttt+0.333333* &
                     apha*ttt*apha*ttt) )
      end if
  185 continue
!
  190 continue
!
      do 195 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 195 i = 1, nxj
       hld1(i,j) = pt(i,jj)
       pdiff(i,jj) = anlslp(i,jj) - pt(i,jj)
!       hld1(i,j) = pdiff(i,jj)
  195 continue
!
      call mpe_unify(hld1,nx,my,2,mpe_double)
      call mpe2d_unify_nx(ww1,anlslp)
      call mpe2d_unify_my(hld2,ww1)
      if(myrank .eq. 0 ) print*,'pt, anlslp at (86,127)= ',hld1(86,127) &
                        ,hld2(86,127)
!
!      call mpe_unify(hld1,nx,my,2,mpe_double)
!      call syslbl ('x00dif',idtg,itaux,ggdef,lrec)
!      if( lreduce.eq.1 ) call reduceintp (hld1,nxdef,nx,my)
!      call dmswrit(nx,my,lrec,lncrec,kflag,ifilout,hld1,istat)
!
!  compute tsave and write out
!
      do 200 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 200 i = 1, nxj
       sigp = sigma(ktrop,1)+sigma(ktrop,2)/pt(i,jj)
       x1= 1.0/(rgas*log(1.0/sigp))
       tsave(i,jj) = x1*(phi(i,ktrop,jj)-sgeo(i,jj))
       hld1(i,j)  = tsave(i,jj)
  200 continue
!
!      call mpe_unify(hld1,nx,my,2,mpe_double)
!      call syslbl ('x00tsv',idtg,itaux,ggdef,lrec)
!      if( lreduce.eq.1 ) call reduceintp (hld1,nxdef,nx,my)
!      call dmswrit (nx,my,lrec,lncrec,kflag,ifilout,hld1,istat)
!
      opok  = 1.0/1000.0**capa
      if( cstar ) then
!
!  construct moisture by assuming a hopefully
!  realistic vertical distribution of moisture profile
!
!
      pk800= opok*exp(.287*log(800.))
      pk300= opok*exp(.287*log(300.))
!
      do 174 jj =1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 173 k = 1, lev
      do 173 i = 1, nxj
      hkd1(i,k) = tt(i,k,jj) - 30.0
      if (pk(i,k,jj).gt.pk300)  hkd1(i,k)= tt(i,k,jj) -10.0
      if (pk(i,k,jj).gt.pk800)  hkd1(i,k)= tt(i,k,jj) -7.0
 173  continue
      call qsatq_2d (nxjp(j),nxp,lev,hkd1(1,1),plt(1,1,jj),vtmp)
      do 175 k = 1, lev
      do 175 i = 1, nxj
       sht(i,k,jj) = vtmp(i,k)
 175  continue
 174  continue
!
      endif
!
!  change real temp to viture potential temp
!

      do 300 jj = 1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 300 k=1,lev
      do 300 i=1,nxj
      tt(i,k,jj)= tt(i,k,jj)*(1.0+0.608*sht(i,k,jj))/pk(i,k,jj)
  300 continue
!
!  read in wind fields at sigma levels
!
      do 320 k = 1, levp
        KL=lev-Llist(k)+1
        if ( KL .lt. 100 ) then
          write (typ, '("m",i2.2,"200")' ) KL
        else
          write (typ, '("n",i2.2,"200")' ) mod(KL,100)
        endif
      call syslbl_r (typ,idtg,itaux,gmdef)
      call dmsread_split (nx,my,lncrec,'H',ifilin,hld1,istat)
      do 320 jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef(j)
        fac = cosl(j)/rad
       if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
       do 320 i = 1,nxj
         hld3(i,k,jj) = hld1(i,j)*fac
  320 continue
      call mpe2d_transpose_ndsl_f2p(hld3,ut, &
            nxp,nx,levf,levp,1,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!
      do 321 k = 1, levp
        KL=lev-Llist(k)+1
        if ( KL .lt. 100 ) then
          write (typ, '("m",i2.2,"210")' ) KL
        else
          write (typ, '("n",i2.2,"210")' ) mod(KL,100)
        endif
      call syslbl_r (typ,idtg,itaux,gmdef)
      call dmsread_split (nx,my,lncrec,'H',ifilin,hld1,istat)
!byl      if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
      do 321 jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef(j)
        fac = cosl(j)/rad
       if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
       do 321 i = 1,nxj
         hld3(i,k,jj) = hld1(i,j)*fac
  321 continue
      call mpe2d_transpose_ndsl_f2p(hld3,vt, &
            nxp,nx,levf,levp,1,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!
!  read in ozone at sigma levels
!  add in june 2010
!
!      do 340 k = 1, lev
!      write (typ, '("m",i2.2,"560")' ) k
!      call syslbl (typ,idtg,itaux,gmdef,lrec)
!      call dmsread (nx,my,lrec,lncrec,'H',ifilin,hld1,istat)
!      if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
!      do 330 jj = 1, jlistnum
!        j=jlist1(jj)
!        nxj=nxdef(j)
!      do 330 i = 1, nxj
!      o3l(i,k,jj) = hld1(i,j)
!  330 continue
!  340 continue
!
      return
      end
