      subroutine rrtmg_gpu &
         ! -------------------------------------------------------------------
         !    -  inputs:
         (sigma, pst, plt, std, tt, qt, o3l, sd, tg, &
          slimsk, cice, xtice, snr, sncover, snoalb, z0, &
          alvsg, alnsg, alvwg, alnwg, facsg, facwg, &
          curate, icsdswg, icsdlwg, &
          sinlj, coslj, xlatj, xlonr, jdat, d2r, xkapa, &
          ptop, dtlw, dtsw, lsswr, lslwr, lssav, &
          nfxr, jjj, &
          nx, nxjdummy, lev, ncld, lprnt, ipt, kdt, &
          uni_cloud, lmfshal, lmfdeep2, &
          deltaq, sup, cnvw, cnvc, &
          ftp, ftp1, fqp, fqp1, nmmiph, &
          !    -  outputs:
          asol, olr, ss, rs, sld, rld, tsflwr, &
          ctot, chig, cmid, clow, &
          cldcov, htrsw, htrlw, &
          fusl, fdsl, fuir, fdir, &
          ! -------------------------------------------------------------------
          fuslr, fdslr, fuirr, fdirr, &
          htrsw0, htrlw0, cosz, &
          asol_clr, olr_clr, ss_clr, rs_clr, &
          sld_clr, rld_clr, sfalb_g, semis_g)
! -------------------------------------------------------------------
! --- for RRTMG scheme :
!
         use physpara
         use module_radiation_driver_gpu, only: grrad_gpu
! -------------------------------------------------------------------
         use mpe
         use rank
         use index
         use radn
         use const, only: RTYPE
         use param, only: my_max, my
         !use nvtx
!      use noah, only:ioutsigr
! -------------------------------------------------------------------
! --- for rrtmg input :
!
         implicit none
         integer i, k, kc, n, myim(my_max), jjj, jj
         integer ntrac, nfxr, nx, nxjdummy, lev, ipt, ncld, kdt, nmmiph, nclds
! --- 3d parameters
!
         real plt(nx, lev, my_max), std(nx, my_max), tg(nx, my_max), &
            tt(nx, lev, my_max), sd(nx, lev, my_max)
         real(kind=RTYPE) qt(nx, lev*ncld, my_max), o3l(nx, lev, my_max), pst(nx, my_max), &
            sigma(lev + 1, 2)
!
! --- 2d parameters
!
         real slimsk(nx, my_max), cice(nx, my_max), xtice(nx, my_max), snr(nx, my_max), sncover(nx, my_max), &
            snoalb(nx, my_max), z0(nx, my_max)
         real alvsg(nx, my_max), alvwg(nx, my_max), alnsg(nx, my_max), alnwg(nx, my_max), facsg(nx, my_max), &
            facwg(nx, my_max), curate(nx, my_max), xlonr(nx, my_max), tsflwr(nx, my_max), cosz(nx, my_max)
         integer icsdlwg(nx, my_max), icsdswg(nx, my_max), jdat(8), j
         real sinlj(my), coslj(my), xlatj(my), ptop, dtlw, dtsw, d2r, xkapa
         logical lsswr, lslwr, lssav, lprnt
         logical uni_cloud, lmfshal, lmfdeep2
         real www, cmax, cmin, imax, imin, tem1, tem2

! --- for grrad input/output (local) :
!
! --- 3d
!
         real prsi(nx, lev + 1, my_max)                      !for even levels
         real prslk(nx, lev, my_max), prsl(nx, lev, my_max)          !for odd levels
         real qgrs(nx, lev, my_max), tgrs(nx, lev, my_max)           !for odd levels
!      real    tracer(nx,lev,ncld),vvl(nx,lev)
         real vvl(nx, lev, my_max)
         real, dimension(:, :, :, :), allocatable ::  tracer

         real slmsk(nx, my_max), xlon(nx, my_max), xlat(nx, my_max), tsfc(nx, my_max), &
            snowd(nx, my_max), sncovr(nx, my_max), zorl(nx, my_max), hprim(nx, my_max), &
            cv(nx,  my_max), cvt(nx, my_max), cvb(nx, my_max), snalb(nx, my_max)

         real alvsf(nx, my_max), alvwf(nx, my_max), alnsf(nx, my_max), alnwf(nx, my_max), &
            facsf(nx, my_max), facwf(nx, my_max), &
            fice(nx, my_max), tisfc(nx, my_max), sinlat(nx, my_max), coslat(nx, my_max)

         real sfalb(nx, my_max), coszen(nx, my_max), coszdg(nx, my_max)
         real tsflw(nx, my_max), semis(nx, my_max)

         integer icsdlw(nx, my_max), icsdsw(nx, my_max)
! --- for pdf cloud
         real sup
         real deltaq(nx, lev, my_max), cnvw(nx, lev, my_max), cnvc(nx, lev, my_max)
! --- for WSM6 & Thompson & GFDL MP
         real ftp(nx, lev, my_max), ftp1(nx, lev, my_max), fqp(nx, lev, my_max), fqp1(nx, lev, my_max)
         real phy3d(nx, lev, 5, my_max)
         real, allocatable, dimension(:, :, :) ::   phy3dnxj

! -------------------------------------------------------------------
! --- for rrtmg output:
         real asol(nx, my_max), olr(nx, my_max), ss(nx, my_max), rs(nx, my_max), sld(nx, my_max), rld(nx, my_max)
         real sfalb_g(nx, my_max), semis_g(nx, my_max)
!
! --- 3d
!
         real htrsw(nx, lev, my_max), htrlw(nx, lev, my_max)
         real(kind=RTYPE) fusl(nx, lev + 1, my_max), fdsl(nx, lev + 1, my_max)
         real(kind=RTYPE) fuir(nx, lev + 1, my_max), fdir(nx, lev + 1, my_max)

         real dummy1(nx, lev, my_max), dummy2(nx, lev, my_max)
         real work1(nx, lev + 1, my_max), work2(nx, lev + 1, my_max)
         real work3(nx, lev + 1, my_max), work4(nx, lev + 1, my_max)
!
! --- 2d
!
!      real    dtrad(nx,lev)
         real ctot(nx, my_max), chig(nx, my_max), cmid(nx, my_max), clow(nx, my_max), csbl(nx)

! for WSM6 & Thompson & GFDL MP
         real cldcov(nx, lev, my_max)   ! input/output layer cloud fraction
         real dummy3(nx, lev, my_max)
         real, allocatable, dimension(:, :) ::   dummy3nxj

!
! --- for clear sky
!
         real asol_clr(nx, my_max), olr_clr(nx, my_max), ss_clr(nx, my_max), rs_clr(nx, my_max)
         real sld_clr(nx, my_max), rld_clr(nx, my_max)
!
! --- 3d
!
         real htrsw0(nx, lev, my_max), htrlw0(nx, lev, my_max)
         real(kind=RTYPE) fuslr(nx, lev + 1, my_max), fdslr(nx, lev + 1, my_max)
         real(kind=RTYPE) fuirr(nx, lev + 1, my_max), fdirr(nx, lev + 1, my_max)
!
         real dummy4(nx, lev, my_max), dummy5(nx, lev, my_max)
         real work5(nx, lev + 1, my_max), work6(nx, lev + 1, my_max)
         real work7(nx, lev + 1, my_max), work8(nx, lev + 1, my_max)
!
! --- 2d
!
         real fluxr(nx, nfxr, my_max)
         character(len=4) :: myrank_str
         integer :: async_id = 1
!
!     if (myrank .eq. 0) print *,'### in rrtmg.f ###'
!     if (myrank .eq. 0) print *,'### j=',j
! -------------------------------------------------------------------
!    to set variables  for grrad input
! -------------------------------------------------------------------
         if (ntoz .eq. 0) then
            ntrac = ncld + 1
         else
            ntrac = ncld
         end if
         allocate (tracer(nx, lev, ntrac, my_max))
         !$acc data create(myim, prsi, prslk, prsl, qgrs, tgrs, vvl, tracer, &
         !$acc&     slmsk, xlon, xlat, tsfc, snowd, sncovr, zorl, hprim, cv, &
         !$acc&     cvt, cvb, snalb, alvsf, alvwf, alnsf, alnwf, facsf, facwf, &
         !$acc&     fice, tisfc, sinlat, coslat, sfalb, coszen, coszdg, tsflw, &
         !$acc&     semis, icsdlw, icsdsw, phy3d, dummy1, dummy2, work1, work2, &
         !$acc&     work3, work4, csbl, fluxr, dummy3, dummy4, dummy5, work5, &
         !$acc&     work6, work7, work8, asol, olr, ss_clr, rs_clr, asol_clr, &
         !$acc&     olr_clr, sld_clr, rld_clr, chig, cmid, clow, ctot, fdirr, &
         !$acc&     fuirr, fdslr, fuslr, fdir, fuir, fdsl, fusl) async(async_id)
         !$acc parallel loop private(j) async(async_id)
         do jj = 1, jlistnum
            j = jlist1(jj)
            myim(jj) = nxjp(j)
         end do

         !$acc parallel loop collapse(2) private(j) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nx
               if (i .le. myim(jj)) then
                  j = jlist1(jj)
                  prsi(i, lev + 1, jj) = (sigma(1, 1)*pst(i, jj) + sigma(1, 2) + ptop)*0.1
               end if
            end do
         end do
         
         !$acc parallel loop gang collapse(2) private(kc) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               kc = lev - k + 1
               !$acc loop vector
               do i = 1, myim(jj)
                  prsi(i, kc, jj) = (sigma(k + 1, 1)*pst(i, jj) + sigma(k + 1, 2) + ptop)*0.1
                  prsl(i, kc, jj) = plt(i, k, jj)*0.1 ! change to cb
                  prslk(i, kc, jj) = (plt(i, k, jj)/1000.)**xkapa
                  tgrs(i, kc, jj) = tt(i, k, jj)
                  qgrs(i, kc, jj) = qt(i, k, jj)
                  vvl(i, kc, jj) = sd(i, k, jj)*0.1                !cb/sec
                  dummy3(i, kc, jj) = cldcov(i, k, jj)
               end do
            end do
         end do

         !$acc parallel loop gang collapse(3) private(kc) async(async_id)
         do jj = 1, jlistnum
            do n = 1, ntrac
               do k = 1, lev
                  kc = lev - k + 1
                  !$acc loop vector
                  do i = 1, myim(jj)
                     tracer(i, kc, n, jj) = qt(i, k + (n - 1)*lev, jj)
                  end do
               end do
            end do
         end do
            !     if (myrank .eq. 0) print *,'### j=',j
            !     if (myrank .eq. 0) print *,'tracer(1,60,1)=',tracer(1,60,1)
            !     if (myrank .eq. 0) print *,'tracer(1,60,2)=',tracer(1,60,2)

            !      do k = 1, lev
            !         kc=lev-k+1
            !      do i = 1, nxj
            !         tracer(i,kc,ntoz) = o3l(i,k)
            !      enddo
            !      enddo
            !
            !      no need the reduction of O3 concentration over model top.
            !      02/23/2024 proposed by Jen-Her Chen
            !
            !      do k=1,8
            !         fac_o3=k*0.1
            !         if(fac_o3 .le. 0.3) fac_o3=0.3
            !         kc=lev-k+1
            !      do i = 1, nxj
   !!       tracer(i,kc,ntoz) = o3l(i,k)*fac_o3
            !       tracer(i,kc,ntoz) = tracer(i,kc,ntoz)*fac_o3
            !      end do
            !      end do
         if (nmmiph .eq. 6 .or. nmmiph .eq. 8 .or. nmmiph .eq. 18) then
            nclds = 3
            ! for MP WSM6 & Thompson effective radius
            !$acc parallel loop gang collapse(2) async(async_id)
            do jj = 1, jlistnum
               do k = 1, lev
                  !$acc loop vector
                  do i = 1, myim(jj)
                     phy3d(i, k, 1, jj) = ftp(i, k, jj)
                     phy3d(i, k, 2, jj) = ftp1(i, k, jj)
                     phy3d(i, k, 3, jj) = fqp(i, k, jj)
                  end do
               end do
            end do
         end if

         if (nmmiph .eq. 11 .or. nmmiph .eq. 12 .or. nmmiph .eq. 13) then
               ! for MP GFDL effective radius
               nclds = 5  ! number of effective cloud condensates used in radiation processes
            !$acc parallel loop gang collapse(2) async(async_id)
            do jj = 1, jlistnum
               do k = 1, lev
                  !$acc loop vector
                  do i = 1, myim(jj)
                     phy3d(i, k, 1, jj) = ftp(i, k, jj)    ! effective radius for liquid water (micron)
                     phy3d(i, k, 2, jj) = ftp1(i, k, jj)   ! effective radius for ice water    (micron)
                     phy3d(i, k, 3, jj) = fqp(i, k, jj)    ! effective radius for snow water   (micron)
                     phy3d(i, k, 4, jj) = fqp1(i, k, jj)   ! effective radius for rain water   (micron)
                  end do
               end do
            end do
         end if

         if (nmmiph .eq. 15 .or. nmmiph .eq. 16) then
               ! for MP Goddard (GCE) effective radius
               nclds = 6
            !$acc parallel loop gang collapse(2) async(async_id)
            do jj = 1, jlistnum
               do k = 1, lev
                  !$acc loop vector
                  do i = 1, myim(jj)
                     phy3d(i, k, 1, jj) = ftp(i, k, jj)    ! effective radius for liquid water (micron)
                     phy3d(i, k, 2, jj) = ftp1(i, k, jj)   ! effective radius for ice water    (micron)
                     phy3d(i, k, 3, jj) = fqp(i, k, jj)    ! effective radius for snow water   (micron)
                     phy3d(i, k, 4, jj) = fqp1(i, k, jj)   ! effective radius for rain water   (micron)
                  end do
               end do
            end do
         end if

         if (nmmiph .eq. 2) then
            nclds = 1
            !$acc parallel loop gang collapse(3) async(async_id)
            do jj = 1, jlistnum
               do kc = 1, 5
                  do k = 1, lev
                     !$acc loop vector
                     do i = 1, myim(jj)
                        phy3d(i, k, kc, jj) = 0.
                     end do
                  end do
               end do
            end do
         end if
            !
            !     if (myrank .eq. 0) print *,'### j=',j
            !     if (myrank .eq. 0) print *,'tracer(1,60,3)=',tracer(1,60,3)

            !--------------------------------------------------------------------
            !    to set xlat(nx), sinlat(nx), coslat(nx) for grrad : input
            !    to set qgrs(nx,lev), tgrs(nx,lev) for grrad : input
            !--------------------------------------------------------------------
         !$acc parallel loop collapse(2) private(j) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nx
               if (i .le. myim(jj)) then
                  j = jlist1(jj)
                  slmsk(i, jj) = slimsk(i, jj)
                  xlon(i, jj) = xlonr(i, jj)
                  xlat(i, jj) = xlatj(j)*d2r
                  !        tsfc(i)=tt(i,lev)          !surface temp in k
                  tsfc(i, jj) = tg(i, jj)
                  snowd(i, jj) = snr(i, jj)
                  sncovr(i, jj) = sncover(i, jj)
                  snalb(i, jj) = snoalb(i, jj)
                  zorl(i, jj) = z0(i, jj)*100.         !surface roughness in cm
                  hprim(i, jj) = std(i, jj)
                  alvsf(i, jj) = alvsg(i, jj)
                  alnsf(i, jj) = alnsg(i, jj)
                  alvwf(i, jj) = alvwg(i, jj)
                  alnwf(i, jj) = alnwg(i, jj)
                  facsf(i, jj) = facsg(i, jj)
                  facwf(i, jj) = facwg(i, jj)
                  fice(i, jj) = cice(i, jj)
                  tisfc(i, jj) = xtice(i, jj)
                  sinlat(i, jj) = sinlj(j)
                  coslat(i, jj) = coslj(j)
                  cvt(i, jj) = 0.
                  cvb(i, jj) = 0.
                  icsdsw(i, jj) = icsdswg(i, jj)
                  icsdlw(i, jj) = icsdlwg(i, jj)
               end if
            end do
         end do

            !     if (myrank .eq. 0) print *,'### j=',j
            !     if (myrank .eq. 0) print *,'xlatj=',xlatj
            !     if (myrank .eq. 0) print *,'xlatj=',xlat(1)

            !     if (myrank .eq. 0) print *,'### j=',j
            !     if (myrank .eq. 0) print *,'sinlj=',sinlj
            !     if (myrank .eq. 0) print *,'coslj=',coslj
            !     if (myrank .eq. 0) print *,'sinlat(1)=',sinlat(1)
            !     if (myrank .eq. 0) print *,'coslat(1)=',coslat(1)

            !     if (myrank .eq. 0) print *,'### j=',j
            !     if (myrank .eq. 0) print *,'qgrs(1,1)=',qgrs(1,1)
            !     if (myrank .eq. 0) print *,'qgrs(1,lev)=',qgrs(1,lev)
            !     if (myrank .eq. 0) print *,'qgrs(ipt,1)=',qgrs(ipt,1)
            !     if (myrank .eq. 0) print *,'qgrs(ipt,lev)=',qgrs(ipt,lev)
            !     if (myrank .eq. 0) print *,'tgrs(1,1)=',tgrs(1,1)
            !     if (myrank .eq. 0) print *,'tgrs(1,lev)=',tgrs(1,lev)
            !     if (myrank .eq. 0) print *,'tgrs(ipt,1)=',tgrs(ipt,1)
            !     if (myrank .eq. 0) print *,'tgrs(ipt,lev)=',tgrs(ipt,lev)
            !------------------------------------------------------------------------
         !!!$acc parallel loop collapse(2) private(cmax, cmin, imin, imax, www) async(async_id)
         !do jj = 1, jlistnum
         !   do i = 1, nx
         !      if (i .le. myim(jj)) then
         !         cmax = -1.0e+25
         !         cmin = 1.0e+25
         !         imin = 1
         !         imax = 1
         !         www = curate(i, jj)*0.0416667*0.1
         !         www = max(www, 5.531e-04)
         !         cv(i, jj) = 0.93 + 0.124*log(www)
         !         if (cv(i, jj) .gt. cmax) then
         !            cmax = cv(i, jj)
         !            imax = i
         !         end if
         !         if (cv(i, jj) .lt. cmin) then
         !            cmin = cv(i, jj)
         !            imin = i
         !         end if
         !      end if
         !   end do
         !end do
         !$acc parallel loop collapse(2) private(www) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nx
               if (i .le. myim(jj)) then
                  www = curate(i, jj)*0.0416667*0.1
                  www = max(www, 5.531e-04)
                  cv(i, jj) = 0.93 + 0.124*log(www)
               end if
            end do
         end do

            !     if (myrank .eq. 0) then
            !     if (myrank .eq. 0) print *,'### j=',j
            !         print *,'### iter=',iter
            !         print *,'### cloud fraction : cv-max =',cmax,' imax=',imax
            !         print *,'### cloud fraction : cv-min =',cmin,' imin=',imin
            !     endif

         tem1 = 0.1
         tem2 = 0.8
         !$acc parallel loop collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nx
               if (i .le. myim(jj)) then
                  if (cv(i, jj) .lt. tem1) cv(i, jj) = 0.0
                  if (cv(i, jj) .gt. tem2) cv(i, jj) = tem2 + 0.01
               end if
            end do
         end do
            !
            !--------------------------------------------------------------------
            !      if (myrank .eq. 0) then
            !     print *,'ptop=', ptop
            !     pstmax=-1.0e10
            !     pstmin=1.0e10
            !     imax=1
            !     imin=1
            !     do i=1,nxj
            !        if (pst(i) .gt. pstmax) then
            !            pstmax=pst(i)
            !            imax=i
            !        endif
            !        if (pst(i) .lt. pstmin ) then
            !            pstmin=pst(i)
            !            imin=i
            !        endif
            !     enddo
            !     if (myrank .eq. 0) print *,'### j=',j
            !     print *,'pst-max=',pstmax,' imax=',imax
            !     print *,'pst-min=',pstmin,' imin=',imin
            !     endif
            !
            !     if (myrank .eq. 0) then
            !     pmax=-1.0e10
            !     pmin=1.0e10
            !     imax=1
            !     imin=1
            !     kmax=1
            !     kmin=1
            !     do k=1,lev+1
            !     do i=1,nxj
            !        if (prsi(i,k) .gt. pmax) then
            !            pmax=prsi(i,k)
            !            imax=i
            !            kmax=k
            !        endif
            !        if (prsi(i,k) .lt. pmin) then
            !            pmin=prsi(i,k)
            !            imin=i
            !            kmin=k
            !        endif
            !     enddo
            !     enddo
            !     if (myrank .eq. 0) print *,'### j=',j
            !     print *,'### prsi : pressure at even level(ptop/pbot)###'
            !     print *,'prsi-max=',pmax,' imax=',imax,' kmax=',kmax
            !     print *,'prsi-min=',pmin,' imin=',imin,' kmin=',kmin
            !     print *,'prsi(ipt,lev+1)=',prsi(ipt,lev+1)
            !     endif
            !--------------------------------------------------------------------
            !     if (myrank .eq. 0) then
            !     pmax=-1.0e10
            !     pmin=1.0e10
            !     imax=1
            !     imin=1
            !     kmax=1
            !     kmin=1
            !     do k=1,lev
            !     do i=1,nxj
            !        if (prsl(i,k) .gt. pmax) then
            !            pmax=prsl(i,k)
            !            imax=i
            !            kmax=k
            !        endif
            !        if (prsl(i,k) .lt. pmin) then
            !            pmin=prsl(i,k)
            !            imin=i
            !            kmin=k
            !        endif
            !     enddo
            !     enddo
            !     if (myrank .eq. 0) print *,'### j=',j
            !     print *,'### prsl : pressure at odd level(mid-level)##'
            !     print *,'prsl-max=',pmax,' imax=',imax,' kmax=',kmax
            !     print *,'prsl-min=',pmin,' imin=',imin,' kmin=',kmin
            !     print *,'prsl(ipt,lev)=',prsl(ipt,lev)
            !     endif
            !--------------------------------------------------------------------
            !     if (myrank .eq. 0) then
            !     pmax=-1.0e10
            !     pmin=1.0e10
            !     imax=1
            !     imin=1
            !     kmax=1
            !     kmin=1
            !     do k=1,lev
            !     do i=1,nxj
            !        if (prslk(i,k) .gt. pmax) then
            !            pmax=prslk(i,k)
            !            imax=i
            !            kmax=k
            !        endif
            !        if (prslk(i,k) .lt. pmin) then
            !            pmin=prslk(i,k)
            !            imin=i
            !            kmin=k
            !        endif
            !     enddo
            !     enddo
            !     if (myrank .eq. 0) print *,'### j=',j
            !     print *,'### prslk : (p_odd/1000.)**xkapa ###'
            !     print *,'prslk-max=',pmax,' imax=',imax,' kmax=',kmax
            !     print *,'prslk-min=',pmin,' imin=',imin,' kmin=',kmin
            !     print *,'prslk(ipt,lev)=',prslk(ipt,lev)
            !     endif
            !--------------------------------------------------------------------
            !     if (myrank .eq. 0) then
            !     pmax=-1.0e10
            !     pmin=1.0e10
            !     imax=1
            !     imin=1
            !     kmax=1
            !     kmin=1
            !     do k=1,lev
            !     do i=1,nxj
            !        if (plt(i,k) .gt. pmax) then
            !            pmax=plt(i,k)
            !            imax=i
            !            kmax=k
            !        endif
            !        if (plt(i,k) .lt. pmin) then
            !            pmin=plt(i,k)
            !            imin=i
            !            kmin=k
            !        endif
            !     enddo
            !     enddo
            !     if (myrank .eq. 0) print *,'### j=',j
            !     print *,'plt-max=',pmax,' imax=',imax,' kmax=',kmax
            !     print *,'plt-min=',pmin,' imin=',imin,' kmin=',kmin
            !     print *,'plt(ipt,lev)=',plt(ipt,lev)
            !     endif
            !--------------------------------------------------------------------
         !write(*,*), 'CCC', nfxr
         !call nvtxStartRange("grrad")
         !write(myrank_str,'(I3)') myrank
         !open(unit=1000, file='grrad_input.'//trim(adjustl(myrank_str)), form='unformatted', &
         !   access='stream', status='replace')
         !write(1000) prsi, prsl, prslk, tgrs, qgrs, tracer, vvl, slmsk, &
         !            xlon, xlat, tsfc, snowd, sncovr, snalb, zorl, hprim, &
         !            alvsf, alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc, &
         !            sinlat, coslat, solhr, jdat, solcon, cv, cvt, cvb, &
         !            icsdsw, icsdlw, ntcw, nclds, ntoz, &
         !            dtlw, dtsw, lsswr, lslwr, lssav, &
         !            nx, myim, lev, lprnt, ipt, kdt, &
         !            ntiw, ntrw, ntsw, ntgl, uni_cloud, lmfshal, lmfdeep2, &
         !            deltaq, sup, cnvw, cnvc, phy3d, dummy3, fluxr
         !close(1000)
         call grrad_gpu(prsi, prsl, prslk, &
                     tgrs, qgrs, tracer, &
                     vvl, slmsk, &
                     xlon, xlat, tsfc, snowd, &
                     sncovr, snalb, zorl, hprim, &
                     alvsf, alnsf, alvwf, alnwf, &
                     facsf, facwf, fice, tisfc, &
                     sinlat, coslat, solhr, jdat, solcon, &
                     cv, cvt, cvb, &
                     icsdsw, icsdlw, ntcw, nclds, ntoz, ntrac, nfxr, &
                     dtlw, dtsw, lsswr, lslwr, lssav, &
                     nx, myim, lev, me, lprnt, ipt, kdt, myrank, &
                     ntiw, ntrw, ntsw, ntgl, uni_cloud, lmfshal, lmfdeep2, &
                     deltaq, sup, cnvw, cnvc, phy3d, -1, &
                     !  ---  outputs:
                     dummy1, sfalb, coszen, coszdg, &
                     dummy2, tsflw, semis, &
                     !  ---  input/output:
                     dummy3, fluxr, &
                     !  ---  optional outputs:
                     dummy4, dummy5, &
                     work1, work2, work3, work4, &
                     work5, work6, work7, work8 &
                     )
         !call nvtxEndRange

         !$acc parallel loop collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nx
               if (i .le. myim(jj)) then
                  !
                  ! total sky
                  !
                  ss(i, jj) = fluxr(i, 4, jj) - fluxr(i, 5, jj)
                  rs(i, jj) = fluxr(i, 7, jj) - fluxr(i, 6, jj)
                  asol(i, jj) = fluxr(i, 1, jj) - fluxr(i, 2, jj)
                  olr(i, jj) = fluxr(i, 3, jj)
                  sld(i, jj) = fluxr(i, 4, jj)
                  rld(i, jj) = fluxr(i, 6, jj)
                  ! clear sky
                  ss_clr(i, jj) = fluxr(i, 24, jj) - fluxr(i, 25, jj)
                  rs_clr(i, jj) = fluxr(i, 27, jj) - fluxr(i, 26, jj)
                  asol_clr(i, jj) = fluxr(i, 1, jj) - fluxr(i, 22, jj)
                  olr_clr(i, jj) = fluxr(i, 23, jj)
                  sld_clr(i, jj) = fluxr(i, 24, jj)
                  rld_clr(i, jj) = fluxr(i, 26, jj)
                  !
                  ! cloud fraction
                  !
                  chig(i, jj) = fluxr(i, 8, jj)   ! high cloud fraction
                  cmid(i, jj) = fluxr(i, 9, jj)   ! middle cloud fraction
                  clow(i, jj) = fluxr(i, 10, jj)  ! low cloud fraction
                  ctot(i, jj) = fluxr(i, 20, jj)  ! total cloud fraction
                  !
                  sfalb_g(i, jj) = sfalb(i, jj)
                  semis_g(i, jj) = semis(i, jj)  ! surface emissivity
               end if
            end do
         end do

         !$acc parallel loop gang collapse(2) private(kc) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               kc = lev - k + 1
               !$acc loop vector
               do i = 1, myim(jj)
                  htrsw(i, kc, jj) = dummy1(i, k, jj)*86400.
                  htrlw(i, kc, jj) = dummy2(i, k, jj)*86400.
                  cldcov(i, kc, jj) = dummy3(i, k, jj)
                  htrsw0(i, kc, jj) = dummy4(i, k, jj)*86400.
                  htrlw0(i, kc, jj) = dummy5(i, k, jj)*86400.
               end do
            end do
         end do

            !       do k = 1, lev
         !$acc parallel loop collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nx
               if (i .le. myim(jj)) then
                  !          dtrad(i,k)=htrsw(i,k)+htrlw(i,k)
                  tsflwr(i, jj) = tsflw(i, jj)
                  cosz(i, jj) = coszen(i, jj)
               end if
            end do
         end do
            !       enddo

            !       if(ioutsigr == 0)then
         !$acc parallel loop gang collapse(2) private(kc) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev + 1
               kc = lev - k + 2
               !$acc loop vector
               do i = 1, myim(jj)
                  fusl(i, kc, jj) = work1(i, k, jj)
                  fdsl(i, kc, jj) = work2(i, k, jj)
                  fuir(i, kc, jj) = work3(i, k, jj)
                  fdir(i, kc, jj) = work4(i, k, jj)
                  fuslr(i, kc, jj) = work5(i, k, jj)
                  fdslr(i, kc, jj) = work6(i, k, jj)
                  fuirr(i, kc, jj) = work7(i, k, jj)
                  fdirr(i, kc, jj) = work8(i, k, jj)
               end do
            end do
         end do
!       endif ! ioutsigr .eq. 0
         !$acc end data
         deallocate (tracer)

         return
      end
