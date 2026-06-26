#define NCCLCHECK(ierr) call nccl_check_helper(ierr, __FILE__, __LINE__)
      subroutine outflds_gpu( itau,nx,my,my_max,lev,ncld             &
             , lmax,numout,idtg                                      &
             , outdir,ktrop,ptop,capa,cp,rgas,grav,sigma,sgeo        &
             , ptend,pt,plt,pk,pk2,phi,ut,vt,vvel                    &
             , tt,qt_org,rdiv,rvor,tg,gwet,z0,hflux,qflux,snr        &
             , raintot,raincu,rainlp,plcl,cumtop,ss,rs,alb,gwclim    &
             , acld,cosl,drag,ugws,vgws,t2,q2,rh2,rh10,u10,v10,gfx,rld,sld &
!byl             , km,smc,slc,stc,canopy,ggdef,slptyp,v850,v700,h850,h500   &
             , km,smc,slc,stc,canopy,ggdef,typtrk                    &
!xb110             , ctot,chig,cmid,clow,hpbl,lwrite,flash,lwritesit)
             , ctot,chig,cmid,clow,hpbl,lwrite,lwritesit)
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
!  output driver subroutine to process sigma level data to standard
!  pressure surfaces and standard grids
!
      use mpe
      use rank
      use index
      use mod_outflds_gpu
      use radn, only : ntcw,ntiw,ntoz
      use const, only : RTYPE,nmmiph,outgrb2,ifilout_grb
      use mod_grb2_param , only :ofdir
      use raddiag, only:clds !cloud fraction on sigma levels
      use mod_qsatq,only:qsatq_gpu,qsatq_3d_gpu
      use openacc
      use cudafor
      implicit  none

      integer   itau,nx,my,my_max,lev,ncld,lmax,numout,ktrop,km
      real      ptop,capa,cp,rgas,grav

      real      pdiff(nxp,my_max)                                                   &
              , t1000(nxp,my_max),plt(nxp,lev,my_max)                               &
!              , qt(nxp,lev*ncld,my_max),rdiv(nxp,lev,my_max)                        &
!              , rvor(nxp,lev,my_max),tg(nxp,my_max),gwet(nxp,my_max)                &
              , tg(nxp,my_max),gwet(nxp,my_max)                                     &
              , z0(nxp,my_max),hflux(nxp,my_max),qflux(nxp,my_max),snr(nxp,my_max)  &
              , raincu(nxp,my_max),rainlp(nxp,my_max),plcl(nxp,my_max),cumtop(nxp,my_max) &
              , ss(nxp,my_max),rs(nxp,my_max),alb(nxp,my_max),gwclim(nxp,my_max)    &
              , acld(lev,my),drag(nxp,lev,my_max)                                   &
              , ugws(nxp,my_max),vgws(nxp,my_max),t2(nxp,my_max)                    &
              , q2(nxp,my_max),rh2(nxp,my_max),rh10(nxp,my_max)                     &
              , u10(nxp,my_max),v10(nxp,my_max),gfx(nxp,my_max),rld(nxp,my_max)     &
              , sld(nxp,my_max),raintot(nxp,my_max)                                 &
!soil
              , smc(nxp,km,my_max),stc(nxp,km,my_max),canopy(nxp,my_max)            &
!noah
              , slc(nxp,km,my_max)                                                  &
! rad-cloud
              , ctot(nxp,my_max),chig(nxp,my_max),cmid(nxp,my_max),clow(nxp,my_max) &
! pbl
              , hpbl(nxp,my_max) 
! river
!byl              , slptyp(nxp,my_max),v850(nx,my),v700(nx,my),h850(nx,my),h500(nx,my)
!byl              , slptyp(nx,my),v850(nx,my),v700(nx,my),h850(nx,my),h500(nx,my)
      real(kind=RTYPE) rdiv(nxp,lev,my_max),rvor(nxp,lev,my_max)    &
                     , ut(nxp,lev,my_max),vt(nxp,lev,my_max)        &
                     , tt(nxp,lev,my_max),qt(nxp,lev*ncld,my_max)   &
                     , qt_org(nxp,lev*ncld,my_max),phi(nxp,lev,my_max) &
                     , sgeo(nxp,my_max),ptend(nxp,my_max)           &
                     , pt(nxp,my_max),vvel(nxp,lev,my_max)          &
                     , sigma(lev+1,2)                               &
                     , pk(nxp,lev,my_max),pk2(nxp,lev,my_max)       &
                     , cosl(my),typtrk(nxp,my_max,5)
!
      character ggdef*4
      integer*8 idtg
!
! local work arrays
!
      real      tmp(nxp,lev,my_max),plog(nxp,lev,my_max),pllp(nxp,my_max)
      real(kind=RTYPE) glob(nx,my)
      real      slp(nxp,my_max)
!
!  pout(16) chnaged into pout(26) to increase p output to 26 levels
!  to respond to the request from regional model
!
      integer,  parameter :: lpout = 31 
      real      wrk1(nxp,lev,my_max),pout(lpout),pkout(lpout),phistd(lpout) &
!              , bt1(nx,my),bt2(nx,my)                               &
              , bt1(nxp,my_max),bt2(nxp,my_max)                      !&
!              , hld1(nxp,my_max),hld2(nxp,my_max) 
!
      real(kind=RTYPE) pres3d(nxp,my_max,lpout)
!
      real(kind=RTYPE) wk_xy(nxp,my_max,12)     ! the last dim is changable
      real(kind=RTYPE) soil_xy(nxp,my_max,12)   ! the last dim is changable
      real      tmpin(nxp,my_max),tmpout(nxp,my_max)
!
!
      real      whtlev(100),whtlevq(100),whtlevz(100)
      character*18 taudir(numout),outdir(numout)
      character*6 labx

      integer   nxmy,nxlev,nxly,ntau,jj,j,nxj,k,i,n,nk,kk,ntrac
      integer,allocatable::   ntt(:),ngq(:)
      integer::ntttmp,ngqtmp
      integer   llts,numz,numq,numt,iqwout,num,nclds
      real      rad,ograv,alaps,rdg,ttb,ttp,ttt,ttt1,ttt2,ttt3,ttt4,anlslp
      real      apha,pl1000,splog,ax,bx,cx,dx,tmid,tsf,tadia,xx,deltap
!
!      integer*8 :: toutsrt, toutend, toutrate
!
      logical :: lwrite,lwritesit

      integer, parameter:: async_id = 1

      integer(kind=cuda_stream_kind) :: stream


!xb110>
!      real      flash(nxp,my_max)         !flash density 
!xb110<


!
!p16  data pout/10.0,20.0,30.0,50.0,70.0,100.0,150.0,200.0,250.0
!p16 1         ,300.0,400.0,500.0,700.0,850.0,925.0,1000.0/
!
!      data pout/10.0,20.0,30.0,50.0,70.0,100.0,150.0,200.0,250.0
      data pout/ 1.0, 2.0, 3.0, 5.0, 7.0                          &
               ,10.0,20.0,30.0,50.0,70.0,100.0,150.0,200.0,250.0  &
               ,300.0,350.0,400.0,450.0,500.0,550.0,600.0,650.0   &
               ,700.0,750.0,800.0,850.0,900.0,925.0,950.0,975.0   &
               ,1000.0/
!
      data rad/6.371e6/

      if(myrank .eq. 0) print*,'   in outflds for tau= ',itau

      call whttau (itau,numout,outdir,ntau,taudir)
      if(ntau.eq.0) return


      allocate( tens(lev+1) )
      do k = 1, lev+1
       tens(k)    = 1.0
      end do
      tens(lev  ) = 0.0
      tens(lev+1) = 0.0
!
!$acc wait(async_id)
!$acc enter data create(tmp,bt1,bt2 ,glob) async(async_id)
!$acc enter data create(pres3d,pkout,plog,pllp,pdiff ) async(async_id)
!$acc enter data create(t1000,wrk1,slp,qt,phistd ) async(async_id)
!$acc enter data create( wk_xy,soil_xy ) async(async_id)
!$acc enter data copyin( tens ) async(async_id)
!$acc enter data copyin(pout) async(async_id)
!$acc wait(async_id)
!

      nxmy = nx*my
      nxlev= nx*lev
      nxly = nx*lev*my
      ograv= 1.0/grav
      alaps = 0.0065
      rdg = rgas/grav
!
      if ( ntoz .gt. 0 ) then
        nclds=ntoz-1
      else
        nclds=ncld
      endif
!  

!$acc parallel loop collapse(2) private( j,nxj ) async(async_id)
      do jj = 1, jlistnum
       do i = 1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)

        do k=1,lpout
        pres3d (i,jj,k) = 0.0
        enddo

        pllp   (i,jj) = 0.0
        bt1    (i,jj) = 0.0
        bt2    (i,jj) = 0.0
       enddo
      enddo

!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
      do k = 1 , 12
       do jj = 1 , my_max
        do i = 1 , nxp
         wk_xy  (i,jj,k) = 0.0
         soil_xy(i,jj,k) = 0.0
        enddo
       enddo
      enddo
!

      if( outgrb2 == 1)then
 134                    format( A  ,A ,I10.10 , i4.4       )
           write(ofdir,134 )trim(ifilout_grb),'/',idtg/100 ,itau
           if(myrank==0) call system("mkdir -p "//trim(ofdir) )
      endif

!
!  copy qt into local qt arrays
!
!$acc parallel loop collapse(3) private( j,nxj ) async(async_id)
      do jj = 1, jlistnum
        do k = 1, lev*ncld
          do i = 1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj) qt(i,k,jj) = qt_org(i,k,jj)  
          enddo
        enddo
      enddo
!
      if(ncld.ge.2)then
        ntrac=2
      else
        ntrac=ncld
      endif
      allocate(ngq(ntrac) , ntt(ntrac)  )
!$acc enter data create(ngq,ntt)  async(async_id)
       ngq(:) = 0
       ntt(:) = 0
!$acc wait(async_id)
!$acc update device(ngq,ntt)  async(async_id)
!$acc wait(async_id)

      
      do n=1,ntrac
        ngqtmp=0
        ntttmp=0
!$acc parallel loop collapse(3) private( j,nxj,kk )  &
!$acc&  reduction(+:ntttmp,ngqtmp )  async(async_id)
        do jj =1,jlistnum
          do k = 1, lev
            do i = 1,nxp
              kk=(n-1)*lev + k
              j=jlist1(jj)
              nxj=nxdef_2d(j)
              if(i<=nxj)then
               if ( qt(i,kk,jj) .lt. 0.0 )  then
                 qt(i,kk,jj) = 0.0
                 !ngq(n) = ngq(n)+1
                 ngqtmp=ngqtmp+1
               endif
               !ntt(n) = ntt(n)+1
               ntttmp=ntttmp+1
             endif
            enddo
          enddo
        enddo
        !$acc wait(async_id)
        ngq(n)=ngqtmp
        ntt(n)=ntttmp
      enddo
   
!!$acc wait(async_id)
!!$acc update self (ngq,ntt) async(async_id)
!$acc wait(async_id)
        call mpe_global_sum(ngq, ntrac ,mpe_integer)
        call mpe_global_sum(ntt, ntrac ,mpe_integer)

!      stream = acc_get_cuda_stream(async_id)
!!$acc host_data use_device(ngq,ntt)
!      NCCLCHECK( ncclAllReduce( ngq , ngq , ntrac, ncclInt32, ncclSum, nccl_comm_gfs, stream ) )
!      NCCLCHECK( ncclAllReduce( ntt , ntt , ntrac, ncclInt32, ncclSum, nccl_comm_gfs, stream ) )
!!$acc end host_data
!!$acc wait(async_id)
!!$acc exit data copyout(ngq,ntt)  async(async_id)
!!$acc wait(async_id)

        if(  myrank.eq.0 .and. ngq(1).ne.0 )  print 901, ngq(1)
        if(  myrank.eq.0 .and. ngq(2).ne.0 )  print 902, ntt(2)-ngq(2)
  901 format ( 1x, " *** warning: in outflds, there are grid points"  &
             , " with q < 0. , total number = ",i10,// )
  902 format ( 1x, " *** checking: in outflds, there are grid points" &
             , " with qc > 0. , total number = ",i10,// )

!$acc wait(async_id)
!$acc exit data delete (ngq,ntt) async(async_id)
!$acc wait(async_id)
      deallocate(ngq,ntt)

!
!                 compute interpolation coeffs
!
!
!  generate structure variables
!
!  first, the log(p) quantities
!
      do k=1,lpout
        pkout(k)= log(pout(k))
      enddo
!$acc update device( pkout ) async(async_id)

!
!  generate std height at p levels
!
      call geostd(lpout,pout,phistd)
!$acc update device( phistd ) async(async_id)
!$acc wait(async_id)


!!
!!  hydrostatic equation
!!
!$acc parallel loop collapse(2) private( j,nxj ) async(async_id)
      do jj =1,jlistnum
        do i=1,nxp
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          if(i<=nxj)then
           phi(i,lev,jj)= cp*tt(i,lev,jj)*(pk2(i,lev,jj)-pk(i,lev,jj)) &
                        + sgeo(i,jj)
           do k=lev-1,1,-1
             phi(i,k,jj)= phi(i,k+1,jj)+cp*(tt(i,k,jj)*(pk2(i,k,jj)-pk(i,k,jj)) &
                        + tt(i,k+1,jj)*(pk(i,k+1,jj)-pk2(i,k,jj)))
           enddo
          endif
        enddo
      enddo
!$acc parallel loop collapse(3) private( j,nxj ) async(async_id)
      do jj =1,jlistnum
        do k=1,lev
          do i=1,nxp
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            if(i<=nxj)then
            plog(i,k,jj)= log(plt(i,k,jj))
            endif
          enddo
        enddo
      enddo
!
!!
!!      if( itau .gt. 0 )then
!!

      call get_prmsl_gpu(nx,my,my_max,lev,ncld &
                    ,sgeo,pt,tt,qt,pk,pk2,plt,phi,pdiff,slp)

!
!      else
!
!  pllp equal to pdiff at the initial time, tau=0
!
!      do jj = 1, jlistnum
!        j=jlist1(jj)
!        nxj=nxdef(j)
!        do i = 1, nxj
!          pllp(i,j) = pdiff(i,jj)
!        enddo
!      enddo
!
!      endif     ! end of ( itau .gt. 0 )
!
! output surface fields
!
!  add terrain pressure output in surfout ( add "ptop" )
!
!$acc wait(async_id)
      if(myrank.eq.0)print*,' outfld : start surfout, lwrite = ',lwrite
      call surfout (nx,my,my_max,itau,idtg,taudir,ntau,pdiff,pt  &
                   ,ptop,typtrk(1,1,1),ptend,glob,ggdef,lwrite)
!
!  obtain the the bottom pressure for the following interpolations
!  bt2 will be used in geoptential interpolations
!
      pl1000= log(1000.1)
!
!$acc parallel loop collapse(2) private( j,nxj ) async(async_id)
      do jj = 1, jlistnum
        do i = 1, nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
          slp(i,jj)= pt(i,jj)+pdiff(i,jj)
          if( slp(i,jj) .le. 1000.1) then
            pllp(i,jj) = 1000.1
          else
            pllp(i,jj) = slp(i,jj)
          endif
         endif
        enddo
      enddo

      call geostd_2d_gpu (nx,my_max,pllp,bt2)
!$acc parallel loop collapse(2) private( j,nxj ) async(async_id)
      do jj = 1, jlistnum
        do i = 1, nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj) pllp(i,jj) = log(pllp(i,jj))
        enddo
      enddo

!
      labx='phi   '
      call whtrec (labx,ntau,taudir,whtlevz,numz)
!
      labx='h2o   '
      call whtrec (labx,ntau,taudir,whtlevq,numq)
!
      labx='tmp   '
      call whtrec (labx,ntau,taudir,whtlev,numt)
!
!$acc wait(async_id)
      if(numz.gt.0.or.numt.gt.0.or.numq.gt.0) then
!
!  temperature output
!
!
!  bottom boundary condition for temperature:
!  limits for 1000 mb boundary condition are (1) standard lapse rate
!  from bottom model level for terrain height and (2) layer mean temp
!  of slp and 1000 mb.  Boundary condition is linear combination
!  of these two cases weighted with subterrainean thickness.
!
      !$acc parallel loop collapse(3) private(j,nxj)  async(async_id)
      do jj =1,jlistnum
        do k = 1, lev
          do i = 1,nxp
            j=jlist1(jj)
            nxj=nxdef_2d(j)
           if(i<=nxj) tmp(i,k,jj) = tt(i,k,jj)*pk(i,k,jj)/(1.0+0.608*qt(i,k,jj))
          enddo
        enddo
      enddo
!
!$acc parallel loop collapse(2) private(i,nxj,ax,bx,cx,dx,tsf, &
!$acc& tadia,splog,tmid) async(async_id)
      do jj =1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
          slp(i,jj) = max( slp(i,jj), pt(i,jj) )
          splog= log(slp(i,jj))
          dx= log(slp(i,jj)/(pt(i,jj)+ptop))
!          bx= tmp(i,lev,jj)*(1.0/sig(lev+1))**0.19023 + sgeo(i,jj)*0.00065
          bx= tmp(i,lev,jj) + sgeo(i,jj)*0.00065
          if ( sgeo(i,jj) .ge. 1000.0 )  then
            tmid = sgeo(i,jj)/(dx*rgas)
            ax = tmp(i,lev,jj)+ (tmid-tmp(i,lev,jj))*(splog-plog(i,lev,jj)) &
                / (log((slp(i,jj)+pt(i,jj)+ptop)*0.5)-plog(i,lev,jj))
          else
            ax = bx
          endif
          cx = (sgeo(i,jj)-1000.0)/(10000.0-1000.0)
          cx = max (0.0, min (1.0, cx))
          dx = cx*cx
          tsf = dx*ax + (1.0-dx)*bx
          tadia= tmp(i,lev,jj)*(slp(i,jj)/plt(i,lev,jj))**0.25
          tsf  = max( tmp(i,lev,jj), min( tsf,tadia ) )
          if ( slp(i,jj) .gt. 1000.1 )  then
            bt1(i,jj)   = tsf
            t1000(i,jj) = tsf
          else
            bt1(i,jj)= tsf + (tsf-tmp(i,lev,jj))*(pl1000-splog)  &
                     / (splog - plog(i,lev,jj))
            t1000(i,jj) = bt1(i,jj)
          endif
         endif
        enddo
      enddo
!
      if(numt.gt.0) then
      if(myrank.eq.0)print*,' outfld : start tempout, lwrite = ',lwrite
        call tempout( nx,my,my_max,lpout,lev,itau,idtg,pout,numt &
               ,whtlev,pkout,plog,pllp,tmp,bt1,pres3d,ggdef,lwrite)
      endif
!
!  moisture output as mixing ratio
!
      if(numq.gt.0) then
!
!  ensure no supersaturated points for output moisture fields
!
      call qsatq_3d_gpu(nx,my_max,lev,tmp,plt,wrk1)
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
      do jj =1,jlistnum
        do k=1, lev
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj) tmp(i,k,jj)= qt(i,k,jj)/wrk1(i,k,jj)
          enddo
        enddo
      enddo
      !$acc parallel loop collapse(2) private(i,nxj) async(async_id)
      do jj =1,jlistnum
        do i=1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj) bt1(i,jj)= tmp(i,lev,jj)
        enddo
      enddo

!!      call mpe_unify(bt1,nx,my,2,mpe_double)
!
      if(myrank.eq.0)print*,' outfld : start shumfout, lwrite = ',lwrite
      call shumout( nx,my,my_max,lpout,lev,itau,idtg,pout,numq &
           ,whtlevq,pkout,plog,pllp,tmp,bt1,pres3d,glob,ggdef,lwrite)
!
!  output clout water content if necessnary
!
      iqwout = 0   ! do output for qw
      if( iqwout .eq. 0 .and. nclds .ge. 2 )then
!  output all hydrometeors
        do ntrac=1,nclds
          if ( ntrac .eq. 1 .or. nclds .gt. 2 ) then
            !$acc parallel loop collapse(3) private(i,nxj,kk) async(async_id)
            do jj = 1, jlistnum
              do k = 1, lev
                do i = 1,nxp
                  j=jlist1(jj)
                  nxj=nxdef_2d(j)
                  if(i<=nxj)then
                    kk = (ntrac-1)*lev+k
                    tmp(i,k,jj)=qt(i,kk,jj)
                  endif
                enddo
              enddo
            enddo
            !$acc parallel loop collapse(2) private(i,nxj) async(async_id)
            do jj = 1, jlistnum
              do i = 1,nxp
                j=jlist1(jj)
                nxj=nxdef_2d(j)
                if(i<=nxj) bt1(i,jj)=tmp(i,lev,jj)
              enddo
            enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
            if(myrank.eq.0)print*,' outfld : start shumout2, lwrite = ',lwrite
            call shumout2( nx,my,my_max,lpout,lev,itau,idtg,pout,numq &
               ,whtlevq,pkout,plog,pllp,tmp,bt1,pres3d,glob,ggdef,ntrac,lwrite )
          endif
        enddo
!
!  output ozone
        if ( ntoz .eq. ncld ) then
         !$acc parallel loop collapse(3) private(i,nxj,kk) async(async_id)
          do jj = 1, jlistnum
            do k = 1, lev
              do i = 1,nxp
              j=jlist1(jj)
              nxj=nxdef_2d(j)
              if(i<=nxj)then
                kk=k+(ntoz-1)*lev
                tmp(i,k,jj)=qt(i,kk,jj)
               endif
              enddo
            enddo
          enddo
         !$acc parallel loop collapse(2) private(i,nxj) async(async_id)
          do jj = 1, jlistnum
            do i = 1,nxp
              j=jlist1(jj)
              nxj=nxdef_2d(j)
              if(i<=nxj) bt1(i,jj)=tmp(i,lev,jj)
            enddo
          enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
          if(myrank.eq.0)print*,' outfld : start shumout2, lwrite =',lwrite
            call shumout2(nx,my,my_max,lpout,lev,itau,idtg,pout,numq &
            ,whtlevq,pkout,plog,pllp,tmp,bt1,pres3d,glob,ggdef,ntoz,lwrite)
        endif
!        
!  output for combination of all condensates
        if ( nmmiph .eq. 18 ) nclds = 6  !do not combine number concentraction for 2M Thompson
        !tmp=0.
!$acc parallel loop collapse(3) private(j,nxj,kk) async(async_id)
        do jj = 1, jlistnum
          do k = 1, lev
            do i = 1,nxp
              j=jlist1(jj)
              nxj=nxdef_2d(j)
              if( i<=nxj)then
                tmp(i,k,jj) = 0.0
                do ntrac=2,nclds
                  kk=k+(ntrac-1)*lev
                  tmp(i,k,jj)=tmp(i,k,jj)+qt(i,kk,jj)
                enddo
              endif
            enddo
          enddo
        enddo
         !$acc parallel loop collapse(2) private(i,nxj) async(async_id)
        do jj = 1, jlistnum
          do i = 1,nxp
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            if(i<=nxj) bt1(i,jj)=tmp(i,lev,jj)
          enddo
        enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
        if(myrank.eq.0)print*,' outfld : start shumout2, lwrite =',lwrite
          call shumout2(nx,my,my_max,lpout,lev,itau,idtg,pout,numq &
          ,whtlevq,pkout,plog,pllp,tmp,bt1,pres3d,glob,ggdef,ncld+1,lwrite)

      endif
!
      endif     ! end of moisture output


!
!  geopotential height output
!
      if(numz.gt.0) then
!$acc parallel loop collapse(2) private( j,nxj ) async(async_id)
        do jj =1,jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj)then
            if(slp(i,jj).le.1000.1) then
              bt1(i,jj)= t1000(i,jj)*rgas*log(slp(i,jj)*0.001)
            else
              bt1(i,jj) = min( sgeo(i,jj), 0.0 )
            endif
           endif
          enddo
        enddo
!
!  compute standard geopotentials, subtract them from sigma
!  level values. vertical interpolation will be done on
!  these deviation from standard values
!
        call geostd_3d_gpu (nx,my_max,lev,plt,wrk1)
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
        do jj =1,jlistnum
          do k=1,lev
            do i=1,nxp
             j=jlist1(jj)
             nxj=nxdef_2d(j)
             if(i<=nxj) tmp(i,k,jj) = ograv*phi(i,k,jj)-wrk1(i,k,jj)
            enddo
          enddo
        enddo
!
!  bottom boundary conditions
!
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
        do jj =1,jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj) bt1(i,jj) = ograv*bt1(i,jj) - bt2(i,jj)
          enddo
        enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
!
        if(myrank.eq.0)print*,' outfld : start geopout, lwrite = ',lwrite
        call geopout (nx,my,my_max,lpout,lev,itau,idtg,pout,numz &
              ,whtlevz,pkout,plog,pllp,tmp,bt1,pres3d,glob,ggdef,phistd  &
              ,typtrk(1,1,4),typtrk(1,1,5),lwrite)
!
      endif   ! end of geopotential height output
!
      endif   ! end of ( itau .gt. 0 )
!
      labx='vor   '
      call whtrec (labx,ntau,taudir,whtlev,num)
      if(num.gt.0) then
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
        do jj =1,jlistnum
          do k=1,lev
            do i=1,nxp
             j=jlist1(jj)
             nxj=nxdef_2d(j)
             if(i<=nxj) tmp(i,k,jj) = rvor(i,k,jj)
            enddo
          enddo
        enddo
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
        do jj =1,jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj) bt1(i,jj)= rvor(i,lev,jj)
          enddo
        enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
        if(myrank.eq.0)print*,' outfld : start vorout, lwrite = ',lwrite
        call vortout (nx,my,my_max,lpout,lev,itau,idtg,pout,num  &
                 ,whtlev,pkout,plog,pllp,tmp,bt1,pres3d,ggdef           &
                 ,typtrk(1,1,2),typtrk(1,1,3),lwrite)
      endif
!
      labx='div   '
      call whtrec (labx,ntau,taudir,whtlev,num)
      if(num.gt.0) then
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
        do jj =1,jlistnum
          do k=1,lev
            do i=1,nxp
             j=jlist1(jj)
             nxj=nxdef_2d(j)
             if(i<=nxj) tmp(i,k,jj) = rdiv(i,k,jj)
            enddo
          enddo
        enddo
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
        do jj =1,jlistnum
          do i=1,nxp
             j=jlist1(jj)
             nxj=nxdef_2d(j)
             if(i<=nxj) bt1(i,jj)= rdiv(i,lev,jj)
          enddo
        enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
        if(myrank.eq.0)print*,' outfld : start divgout, lwrite = ',lwrite
        call divgout (nx,my,my_max,lpout,lev,itau,idtg,pout,num &
                  ,whtlev,pkout,plog,pllp,tmp,bt1,pres3d,ggdef,lwrite)
      endif
!
      labx='wnd   '
      call whtrec (labx,ntau,taudir,whtlev,num)
      if(num.gt.0) then
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
        do jj =1,jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj)then
            bt1(i,jj)= ut(i,lev,jj)
            bt2(i,jj)= vt(i,lev,jj)
           endif
          enddo
        enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
!!        call mpe_unify(bt2,nx,my,2,mpe_double)
        if(myrank.eq.0)print*,' outfld : start windout, lwrite = ',lwrite
        call windout (nx,my,my_max,lpout,lev,itau,idtg,pout,num &
           ,whtlev,cosl,pkout,plog,pllp,ut,vt,vvel,tmp,bt1,bt2,pres3d,glob,ggdef,lwrite)
      endif
!
      labx='dag   '
      call whtrec (labx,ntau,taudir,whtlev,num)
      if(num.gt.0) then
!$acc enter data copyin(drag) async(async_id)
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
        do jj =1,jlistnum
          do k=1,lev
            do i=1,nxp
             j=jlist1(jj)
             nxj=nxdef_2d(j)
             if(i<=nxj) tmp(i,k,jj) = drag(i,k,jj)
            enddo
          enddo
        enddo
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
        do jj =1,jlistnum
          do i=1,nxp
             j=jlist1(jj)
             nxj=nxdef_2d(j)
             if(i<=nxj) bt1(i,jj)= drag(i,lev,jj)
          enddo
        enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
      if(myrank.eq.0)print*,' outfld : start dragout, lwrite = ',lwrite
        call dragout (nx,my,my_max,lpout,lev,itau,idtg,pout,num &
                  ,whtlev,pkout,plog,pllp,tmp,bt1,pres3d,ggdef,lwrite)
      endif
!$acc exit data delete(drag) async(async_id)
!
!!output Fraction of cloud cover on pressure levels
!!clouds output
      labx='cld   '
      call whtrec (labx,ntau,taudir,whtlev,num)
      if(num.gt.0) then
!        tmp=0.
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
        do jj = 1, jlistnum
          do k = 1, lev
            do i = 1,nxp
             j=jlist1(jj)
             nxj=nxdef_2d(j)
             if(i<=nxj) tmp(i,k,jj)=clds(i,k,jj) * 100.0
            enddo
          enddo
        enddo
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
        do jj = 1, jlistnum
          do i = 1,nxp
             j=jlist1(jj)
             nxj=nxdef_2d(j)
             if(i<=nxj) bt1(i,jj)=tmp(i,lev,jj)
          enddo
        enddo
          call cloudout(nx,my,my_max,lpout,lev,itau,idtg,pout,num &
          ,whtlev,pkout,plog,pllp,tmp,bt1,pres3d,glob,ggdef,lwrite)
      endif
!
      if(lwritesit) then
        labx='sit   '
        call whtrec (labx,ntau,taudir,whtlev,num)
        if(num.gt.0) then
          call sitout(nx,my,my_max,itau,idtg,num,whtlev,ggdef)
        endif
      endif

!
!  wk_xy(-,-,1) : temperature at the lowest sigma level
!  wk_xy(-,-,2) : u wind at the lowest sigma level
!  wk_xy(-,-,3) : v wind at the lowest sigma level
!  wk_xy(-,-,4) : precipitable water
!  wk_xy(-,-,5) : relative humidity at the lowest sigma level
!

      if(ncld.ge.2)then
        ntrac=2
      else
        ntrac=ncld
      endif


!$acc parallel loop collapse(2) private(j,nxj,xx) async(async_id)
      do jj = 1 , jlistnum
        do i = 1 , nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj)then
         xx=rad/cosl(j)
         wk_xy(i,jj,1) = tt(i,lev,jj)*pk(i,lev,jj)/(1.0+0.608*qt(i,lev,jj))
         wk_xy(i,jj,2) = ut(i,lev,jj)*xx
         wk_xy(i,jj,3) = vt(i,lev,jj)*xx
        endif
       enddo
      enddo

!$acc parallel loop collapse(2) private(j,nxj,kk,deltap) async(async_id)
      do jj = 1, jlistnum
       do i = 1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj)then
         wk_xy(i,jj,4) = 0.0
         do n = 1, ntrac
          do k = 1, lev
            kk=(n-1)*lev+k
            deltap = ( pt(i,jj)*(sigma(k+1,1)-sigma(k,1))+ &
                     (sigma(k+1,2)-sigma(k,2)) ) * 100./grav
            wk_xy(i,jj,4) = wk_xy(i,jj,4) + qt(i,kk,jj)*deltap
          enddo
         enddo !ntrac
        endif
       end do
      end do


!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
       do i = 1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj)then
         tmpin(i,jj)=wk_xy(i,jj,1)
         bt1  (i,jj)=plt(i,lev,jj)
        endif
       enddo
      enddo
      call qsatq_gpu(nxp,my_max,tmpin ,bt1 ,tmpout)
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
       do i = 1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj)then
          wk_xy(i,jj,5) = tmpout(i,jj)
          wk_xy(i,jj,5) = 100.*(qt(i,lev,jj)/wk_xy(i,jj,5))
          wk_xy(i,jj,5) = min( 100., max( 1., wk_xy(i,jj,5) ) )
        endif
       enddo
      enddo

!
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj =1,jlistnum
       do i=1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj)then
          wk_xy(i,jj,6)  = ctot(i,jj) * 100.0
          wk_xy(i,jj,7)  = chig(i,jj) * 100.0
          wk_xy(i,jj,8)  = cmid(i,jj) * 100.0
          wk_xy(i,jj,9)  = clow(i,jj) * 100.0
          wk_xy(i,jj,10) = hpbl(i,jj)
          wk_xy(i,jj,11) = qt(i,1,jj)
          wk_xy(i,jj,12) = qt(i,lev,jj)
        endif
       enddo
      enddo

!soil variable
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
        do i = 1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj)then
          soil_xy(i,jj,1 ) = smc(i,1,jj)
          soil_xy(i,jj,2 ) = smc(i,2,jj)
          soil_xy(i,jj,3 ) = smc(i,3,jj)
          soil_xy(i,jj,4 ) = smc(i,4,jj)
          soil_xy(i,jj,5 ) = slc(i,1,jj)
          soil_xy(i,jj,6 ) = slc(i,2,jj)
          soil_xy(i,jj,7 ) = slc(i,3,jj)
          soil_xy(i,jj,8 ) = slc(i,4,jj)
          soil_xy(i,jj,9 ) = stc(i,1,jj)
          soil_xy(i,jj,10) = stc(i,2,jj)
          soil_xy(i,jj,11) = stc(i,3,jj)
          soil_xy(i,jj,12) = stc(i,4,jj)
        endif
        enddo
      enddo
!
! output some 2-dimension veriable to dmsfile
!
      if(myrank .eq. 0) print *,'call out2d'
      if ( lwrite )                                                 &
      call out2d_gpu (nx,lev,my,my_max,itau,idtg,taudir,ntau        &
                 ,hflux,qflux,tg,gwet,snr,z0,raintot,raincu,rainlp  &
                 ,plcl,cumtop,ss,rs,alb,gwclim,glob                 &
                 ,acld,ugws,vgws,t2,q2,rh2,rh10,u10,v10,gfx,rld     &
!xb110                 ,sld,wk_xy,soil_xy,canopy,ggdef,lwrite,flash)
                 ,sld,wk_xy,soil_xy,canopy,ggdef)


!
!$acc wait(async_id)
!$acc update self(typtrk) async(async_id)
!$acc exit data delete( wk_xy,soil_xy ) async(async_id)
!$acc exit data delete( tmp,bt1,bt2 ,glob) async(async_id)
!$acc exit data delete( pres3d,pkout,plog,pllp,pdiff) async(async_id)
!$acc exit data delete( t1000,wrk1,slp,qt,phistd ) async(async_id)
!$acc exit data delete( pout) async(async_id)
!$acc exit data delete( tens ) async(async_id)
!$acc wait(async_id)
      deallocate( tens )

      return
      end
