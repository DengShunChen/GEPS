#define ndslfv_monoadvv ndslfv_monoadvv_gpu
#define ndslfv_monoadvv_fgnl ndslfv_monoadvv_fgnl_gpu
#define ndslfv_update  ndslfv_update_gpu
subroutine intgrt_gpu
   !$acc routine(prexp_hybrid_cwb_gpu) vector
   !$acc routine(gridnl_hybrid_ndsl_gpu) vector
   !
   !***********************************************************************
   !  this subroutine is the basic time stepping driver.  it does the
   !  gaussian quadrature integrations to compute the adiabatic spectral
   !  tendencies for the model variables. these tendencies are adjusted
   !  with a semi-implicit method to make the model stable for long
   !  time steps.  time steps are made with a time filter to remove
   !  any computational modes.
   !
   !  modify to f90 by C-H Lee and sort by River Chen in 2015
   !
   !***********************************************************************
   !
   use param
   use mpe
   use rank
   use index
   use const
   use spec
   use grid
   use fftcom
   use phygrid
   use mod_typhoon
   use noah
   use namelist_soilveg
   ! use nvtx
   use cudafor
   use openacc
   use mod_ndslfv_monoadv_gpu, only: ndslfv_monoadvh_fgnl_gpu, &
                                     ndslfv_monoadvh_gpu, &
                                     ndslfv_monoadvv_gpu, ndslfv_monoadvv_fgnl_gpu, &
                                     ndslfv_update_gpu
   !-----------------------------------------------------------------------
   USE mod_sitgrid
   USE mod_sit_vdiff, ONLY: sit_vdiff_end, cal_ratioBlending
   USE mod_sit_control, ONLY: xmissing, lgodas, locaf, lwoa0 &
                              , loutsit24, outsitmean, lsitstart &
                              , ldailysst, ltrigsit
   USE mod_sst, ONLY: deallocate_ocaf_array, deallocate_woa0_array &
                      , deallocate_godas_array, read_dailygodas &
                      , time_weights, nmw1, nmw2, wgt1, wgt2 &
                      , read_dailyFCT, dailyFCTsst, dailyFCTcice &
                      , dailyFCTsndepth, deallocate_dailyFCT_array &
                      , ANAsstT0, dailyClmANAsst, dailyClmFCTsst &
                      , allocate_opgsst_array, read_opgsst &
                      , opgsst, deallocate_opgsst_array &
                      , obswtbnmw1, obswtbnmw2, obswtbwgt1, obswtbwgt2 &
                      , tseaold, tseanow, tseanew, dtseadt &
                      , tseadiffFCT, tseadiffFCT24, outtseadiffFCT24
   USE mo_netcdf, ONLY: lkvl, cleanup_netcdf
   !-----------------------------------------------------------------------
   use raddiag
   use radn
   use albn
   !-----------------------------------------------------------------------
   use mod_stochastic_physics, only: spptout, skebout, &
                                     run_stochastic_physics, &
                                     destroy_stochastic_physics, &
                                     skeb3du, skeb3dv, diss_est, skebfilt, &
                                     keb, kea, skebest
   !-----------------------------------------------------------------------
   implicit none
   !
   !  local work array
   !
   integer nfxr
   !  for Semi-Lagrangian
   !
   real(kind=RTYPE) ndsldtah, dtah, dtahi, dta, &
      diveng(nxp, lev, my_max), &
      qm_sl(nx, levp*ncld, my_max), &
      pten_sl(nx, levp, my_max), &
      uum_sl(nx, levp, my_max), vvm_sl(nx, levp, my_max), &
      ttm_sl(nx, levp, my_max), &
      um(nxp, lev, my_max), vm(nxp, lev, my_max), &
      tm(nxp, lev, my_max), &
      vdmerd(nxp, lev, my_max), vdzonl(nxp, lev, my_max), &
      vdmerdg(nxp, lev, my_max), vdzonlg(nxp, lev, my_max), &
      vdmerdr(nxp, lev, my_max), vdzonlr(nxp, lev, my_max), &
      vdmerdrp(nxp, lev, my_max), vdzonlrp(nxp, lev, my_max), &
      ddtemp(nxp, lev, my_max), &
      pten(nxp, lev, my_max), dummy, &
      ptm(nxp, my_max), &
      deldm(nxp, my_max), sdpbl(nxp, my_max)

   !          integer ierr, ittw, itt, year, yrd
   integer ierr, ittw, itt, year
   !
   real(kind=RTYPE) glob(nx, my)
   real hf24(nxp, my_max), qf24(nxp, my_max), ss24(nxp, my_max), rs24(nxp, my_max), &
      asol24(nxp, my_max), olr24(nxp, my_max), rain24(nxp, my_max), &
      rainlp24(nxp, my_max), totallp(nxp, my_max), &
      drag(nxp, lev, my_max), ugws(nxp, my_max), vgws(nxp, my_max)

   integer kn

   real tmin(nxp, my_max), tmax(nxp, my_max)!,td(nxp,my_max),temp
   !
   real(kind=RTYPE) pltemp(jtrun, jtmax, 2), cc(nx + 2, levp, 1, my_max), &
      ww1(nx, my_max)
   !byl      real      dlgeo(nxp,my_max),dtgeo(nxp,my_max)
   !byl      real      cc3(nx+2,levp,3,my_max),wss3(levp,2,3,jtrun,jtmax)
   !
   character rfile*255, ctau*7, ccore*4
  !!      real      tbar(lev),qbar(lev*ncld),qbrrow(ncld,my)
   integer, parameter :: ktop = 4
   real fac(ktop), wkj(my, 4), wkmf(jtrun), wkmf_local, windmax3
   data windmax3/130./
   !
   logical histim, tchange, flag, forward, fwd
   !
   logical wrestrt
   !
   ! for topographic gravity wave drag
   !
   real hprime_b(nxp, mtnvar, my_max)

   ! for due point temperature
   real tda, tdb
   parameter(tda=17.27, tdb=237.7)

   !
   !   restart  : write(7) work array
   !
   real, dimension(:), allocatable :: work_io
   integer itauezz
   real dtauzz, tautv
   character cmdxx*256
   integer len1, len2
   !-------------------------------------------------------------------
   integer i, j, k, m, n, jj, kk, mf, kw, nxj, nml, lmax, leng, nxmy, &
      jlim, mlst, mlmax2, itaui, itaue, itauo, itaup, &
      ntau, itau, lcwb, lphy, ifromtau, itotau, istat, &
      istst, ii, n_stable, n_unstable, nc_stable, nxjf

   real www, dtx, dtq, thdai, tkei, tpei, dsigp, &
      cosw, tengi, dt24, tg2, dtx_tau, sqhaf, &
      dt1, sptend, wmax, xx, dtaup, hfiltm, &
      sptendmax2, sptendmax1, dt_chg
   integer recn

   ! for io quilting
   character*34 keydoit, keydone
   data keydoit/"DOIT..........................DOIT"/
   data keydone/"DONE..........................DONE"/
   !
   !for sst_restore_tau>0., update sst(W00100), seaice(W00091), snowdepth(B00650)
   integer*8 idtg_sst, idtg1_sst, idtg_temp
   integer icurrenttau, yyyymmdd, hhii
   logical lsstrestore, iceold(nxp, my_max), oceanold(nxp, my_max)
   character*12 cdtg
   real ssttemp, cicetemp, snrtemp
   real sst(nxp, my_max), ssttau, tautemp
   integer yr, mo, dy, hr, mn

   !for opgsst sst
   integer*8 idtg1_opgsst
   integer icurrentyear, inextyear, inexttau
   logical lnewyear
   real tauleft

   !for SIT
   integer icurrentyymmdd, inextyymmdd
   logical lnewday
   integer icurrentyymm, ibeforeyymm
   logical lnewyymm
   integer ic_sit, nc_sit
   logical turn_sit, lrun_sitvdiff
   !pscheckdata
   integer lenc, itautest
   real mout(nx, my)
   integer nc
   integer*8 :: toutsrt, toutend, toutrate      !For CPU Timings
   logical:: lopngrb2 !open grib2 file

   !CWB2021 ndsl single precision test
   real(kind=RTYPE) &
      pdot(nxp, lev + 1, latpart)
   !
   !xb110>
   !byl      real rmr(nxp,lev,my_max),smr(nxp,lev,my_max)
   !for lightning scheme from ECMWF
   real flash(nxp, my_max), flash24(nxp, my_max)
   !xb110<

   integer, parameter:: async_id = 1
   integer(kind=cuda_stream_kind) :: stream
#ifdef TIMING
   ! for timing
   real*8 tm_1, tm_2, tm_use, mpi_wtime
   !CWB2020
   tm_1 = mpi_wtime()
   tm_2 = mpi_wtime()
#endif
   ! ------------------------------------------------------------
   ! << OpenACC ALLOCATE DATA >>
   stream = acc_get_cuda_stream(async_id)
   !$acc enter data copyin(nxjp, nxp, lev, ncld, ndslvvar) async(async_id)
   !$acc enter data copyin(nx, my, my_max, jlistnum, mlistnum, jtrun, jtmax, levp, nsizey, nsizex) async(async_id)
   !$acc enter data copyin(jlist1, nxdef_2d, nxjp_acc, nxptot, &
   !$acc& jlist2_2d, nxjlen_all, nxdef, mlist, Llist, &
   !$acc& gglati, fa1, fa2, fa3, fa4) async(async_id)
  !! << const >>
   !$acc enter data copyin(cp, rad, radsq, sinl, cosl, &
   !$acc& onocos, cor, sigma, dsigma, ptop, hdk1, hdk2) async(async_id)
   !$acc enter data copyin(poly, dpoly, weight, cim, wcfac, wdfac) async(async_id)
   !$acc enter data copyin(ptmeans, spalm, eps4, eigval, evecin ,evectr, &
   !$acc& arrhyd, arsddt) async(async_id)
   !$acc enter data copyin(hfiltx, alphax) async(async_id)
  !! << grid >>
   !$acc enter data copyin(sgeo) async(async_id)
  !! << spec >>
   !$acc enter data copyin(trefs) async(async_id)
   ! ++++++++++++++++++++++++++++++++++++++++++++++++++
  !! << const >>
   !$acc enter data create(pcorr) async(async_id)
  !! << grid >>
   !$acc enter data create(tt, ut, vt, qt, qm) async(async_id)
   !$acc enter data create(dtphi, dlphi, dlpl, dtpl) async(async_id)
   !$acc enter data create(rdivm, rdiv, phi, sd, vvel) async(async_id)
   !$acc enter data create(pt, ptend) async(async_id)
   !$acc enter data create(up, vp, ttp, ptp) async(async_id)
   !$acc enter data create(pk, pk2, plt) async(async_id)
  !! << spec >>
   !$acc enter data create( &
   !$acc& temold, divold, vorold, plold, &
   !$acc& temmid, divmid, vormid, plmid, &
   !$acc& temnow, divnow, vornow, plnow, &
   !$acc& temten, divten, vorten, plten, hldten &
   !$acc& ) async(async_id)
  !! << local >>
   !$acc enter data create(diveng, tm, um, vm, &
   !$acc& vdzonl, vdmerd, vdmerdg, vdzonlg, vdzonlrp, vdmerdrp, &
   !$acc& ddtemp, ptm, deldm, pdot, sdpbl, pltemp) async(async_id)

   !$acc enter data create(cc, ww1, pten, hfiltm) async(async_id)
   !$acc wait(async_id)
   ! ------------------------------------------------------------
   !      fsit=-99.             !fsit>0., turn on sit_vdiff when mod(tau/fsit)<0.001
   !default fsit<=0., turn on sit_vdiff every tau
   ic_sit = -99             !if fsit>0., store now tau is the ic_sit times when sit_vidff is turn on
   nc_sit = 1               !if fsit>0., when mod(tau/fsit)<0.001, turn on sit_vdiff for "nc_sit" timesteps
   turn_sit = .false.       !turn_sit=.true., will run sit_vdiff in some tau
   lrun_sitvdiff = .false.  !lrun_sitvdiff=.true., run sit_vdiff in this tau

   if (dorst) then
      wrestrt = .true.
   else
      wrestrt = .false.
   end if
   !
   ttm_sl = 0.
   pten_sl = 0.
   qm_sl = 0.
   vvm_sl = 0.
   uum_sl = 0.
   !
   lmax = 16
   nfxr = 33
   cc = 0.
   istat = 0
   !
   ! TYW added 20240112
   hprime_b = 0.
   !
   year = idate(1)
   !      yrd  = 365
   !      if ( mod(year,4) .eq. 0 ) yrd = 366
   !
   do jj = 1, jlistnum
      j = jlist1(jj)
      !ch     nxj=nxdef(j)
      call prexp_hybrid_cwb(nxjp(j), nxp, lev, ptop, sigma, pt(1, jj), &
                            pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
   end do
   !      call  outsigs ( 0,nx,my,my_max,lev,ncld                                    &
   !                    , idtg,ifilout,ptop,rad,grav                                 &
   !                    , cp,cosl,pt,sgeo,snr,gwr,tg,pk,pk2                          &
   !                    , ut,vt,tt,qt,phi,rdiv                                       &
   !                    , km_soil,smc,slc,stc,canopy,zice,ggdef,gmdef )

   raintot = 0.
   raincu = 0.
   rainlp = 0.
   raincu6 = 0.
   rainlp6 = 0.
   raincu3 = 0.
   rainlp3 = 0.
   raincu1 = 0.
   rainlp1 = 0.
   gfx = 0.
   if (.not. restrt) then
      rld = 0.
      sld = 0.
   end if
   recn = 1
   rdivm = 0.
   flash = 0.
   pdry = 0.
   !
   !
   ! read mountant variables for topographic gravity wave drag
   !
   if (yesdia .and. dograv .and. nmgwor .eq. 2) then
      call read_mtnvar(nx, my, mtnvar, hprime_b)
      !
      if (myrank .eq. 0) &
         print *, 'read mtnvar=14 hprime_b=', (hprime_b(1, i, 1), i=1, mtnvar)
   end if
   !
   tchange = .false.
   nxmy = nx*my
   nml = nx*my*lev
   mlmax2 = mlmax*2
   leng = mlmax*2*lev
   !
   itaui = taui + 0.1
   itaue = taue + 0.1
   itauo = tauo + 0.1
   itaup = taup + 0.1
   if (.not. restrt) then
      tau = taui
   end if
   dtx = dt
   !
   dta = dtx
   dtah = 0.5*dta
   ndsldtah = dtah/float(itter)
   dtahi = dtah/float(itter)
   !$acc enter data copyin(dta, dtah, ndsldtah, dtahi)
   !
   !jwhwu 201407 add time control
   if (dorst) then
      open (7, file='./timectl', status='old')
      read (7, '(i8)') itauezz
      tautv = 24.           !! history output directory control
      close (7)
      if (myrank .eq. 0) then
         print *, 'the integration will be extended up to ', itauezz, ' hours'
      end if
   end if
   !#endif
   !
   !  compute initial moisture and potential temperature
   !
   qgini = 0.0
   thdai = 0.0
   tkei = 0.0
   tpei = 0.0
   do i = 1, 4
      do j = 1, my
         wkj(j, i) = 0.
      end do
   end do
   !
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do k = 1, lev
         do i = 1, nxj
            dsigp = dsigma(k, 1)*pt(i, jj) + dsigma(k, 2)
            wkj(j, 1) = wkj(j, 1) + qt(i, k, jj)*dsigp*cosl(j)
            wkj(j, 2) = wkj(j, 2) + tt(i, k, jj)*dsigp*cosl(j) &
                        /(1.+0.608*qt(i, k, jj))
            wkj(j, 3) = wkj(j, 3) + (ut(i, k, jj)**2 + vt(i, k, jj)**2)*dsigp &
                        *radsq/(2.*cosl(j))
            wkj(j, 4) = wkj(j, 4) + cp*tt(i, k, jj)*pk(i, k, jj)*dsigp &
                        *cosl(j)
         end do
      end do
   end do

   call mpe_unify(wkj, my, 4, 1, mpe_double)
   !
   cosw = 0.
   do j = 1, my
      nxj = nxdef(j)
      cosw = cosw + cosl(j)*nxj
      qgini = qgini + wkj(j, 1)
      thdai = thdai + wkj(j, 2)
      tkei = tkei + wkj(j, 3)
      tpei = tpei + wkj(j, 4)
   end do
   tengi = tkei + tpei
   qgini = qgini*100./grav/cosw
   cosw = cosw*lev
   thdai = thdai/cosw
   tengi = tengi/cosw
   !
   !
   if (myrank .eq. 0) &
      print *, 'qgini, thdai, tengi= ', qgini, thdai, tengi
   !
   if (myrank .eq. 0) print *, ' beginning integration '
   !
   !  zero out precip arrays
   !
   dt24 = 0.
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do i = 1, nxj
         ss24(i, jj) = 0.
         rs24(i, jj) = 0.
         hf24(i, jj) = 0.
         qf24(i, jj) = 0.
         asol24(i, jj) = 0.
         olr24(i, jj) = 0.
         rain24(i, jj) = 0.
         rainlp24(i, jj) = 0.
         totallp(i, jj) = 0.
         raincu(i, jj) = 0.
         rainlp(i, jj) = 0.
         raincu6(i, jj) = 0.
         rainlp6(i, jj) = 0.
         raincu3(i, jj) = 0.
         rainlp3(i, jj) = 0.
         raincu1(i, jj) = 0.
         rainlp1(i, jj) = 0.
         raintot(i, jj) = 0.
         runoff(i, jj) = 0.  ! soil
         tmax(i, jj) = 0.
         tmin(i, jj) = 0.
         flash24(i, jj) = 0.  !xb110, flash density
      end do
   end do
   !
   if (itaui .eq. 0) then   ! when restart, don't zero out
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do i = 1, nxj
            ss(i, jj) = 0.
            rs(i, jj) = 0.
            sld(i, jj) = 0.
            tg2 = tg(i, jj)*tg(i, jj)     !soil
            rld(i, jj) = stbo*(tg2*tg2)   !soil
         end do
      end do
   end if
   !xb110>
   !byl      rmr = 0.
   !byl      smr = 0.
   !xb110<
   !
   !***********************************************************************
   !     start time integration iterations
   !***********************************************************************
   !$acc update device(dtpl, dlpl) async(async_id)
   !$acc update device(rdiv, ut, vt, tt, qt, pt) async(async_id)
   !$acc update device(vornow, divnow, temnow, plnow) async(async_id)
   !
   ! sppt
   if (.not. restrt) then
      itimestep = 1
   end if
   !
   n_stable = 0
   n_unstable = 0

   if (dta .gt. 720) then
     !!        dt_chg=1800.
      nc_stable = 1
      !        sptendmax2=0.3605
      !        sptendmax1=0.2605
      sptendmax2 = 0.4005
      sptendmax1 = 0.3005
   else if (dta .le. 720 .and. dta .gt. 450) then
     !!        dt_chg=720.
      nc_stable = 2
      sptendmax2 = 0.4005
      sptendmax1 = 0.3005
   else if (dta .le. 450) then
     !!        dt_chg=450.
      nc_stable = 4
      sptendmax2 = 0.4305
      sptendmax1 = 0.3305
   end if

   !
   !      if(typhoon)then
   !        dt_trk=6.
   !      else
   !        dt_trk=tauo
   !      endif
   ! store first idtg to idtg_sst
   icurrenttau = int(tau)
   call dtgfix12(idtg, idtg_sst, icurrenttau)
   if (myrank .eq. 0) print *, 'idtg_sst=', idtg_sst

   ! read opgsst data
   if (lopgsst) then
      call allocate_opgsst_array
      if (myrank .eq. 0) print *, 'myrank=', myrank, 'idtg_sst=', idtg_sst
      call read_opgsst(idtg_sst, ggdef, ocean, ice)
      icurrentyear = idtg_sst/100000000
   end if   !end lopgsst
   !
10 continue

   dtx_tau = dtx/3600.

   if (myrank .eq. 0) then
      print *, 'forcast begin tau=', itaui, ' to tau=', itaue

      !     ! for io quilting
      !      if(io_quilting)then
      !         ntag=ntag+1
      !         call mpe_send_key(keydoit,ntag,istat)
      !      endif

   end if

   !CWB2020 fix the bug for negative timing info
   !#ifdef TIMING
   !      tm_1=mpi_wtime()
   !#endif
100 continue
   !
   tau = tau + dtx/3600.
   ! call nvtxStartRange("Intgrt")
#ifdef TIMING
   tm_use = tm_2 - tm_1
   if (myrank .eq. 0) print 265, tau, tm_use
265 format(' tau= ', f8.3, ',     Timing=', f8.3, ' elapse seconds')
   tm_1 = mpi_wtime()
#else
   if (myrank .eq. 0) print 265, tau
265 format(' tau= ', f8.3)
#endif
   !
   if (do_sit) then
      lsitstart = (itimestep .le. 1)
      if (fsit .le. 0) then
         lrun_sitvdiff = .true.
         ic_sit = -99
         turn_sit = .true.
      else
         if (itimestep .gt. 1) then
            turn_sit = mod(tau + 0.001, fsit) .lt. dtx_tau
            if (turn_sit) ic_sit = 0
            if (turn_sit .or. (ic_sit .ge. 0 .AND. ic_sit .lt. nc_sit)) then
               lrun_sitvdiff = .true.
               ic_sit = ic_sit + 1
            else
               lrun_sitvdiff = .false.
            end if
         end if
      end if
      if (myrank .eq. 0) then
         print *, 'final turn_sit=', turn_sit, ',run_sitvdiff=', lrun_sitvdiff
      end if
   end if
   !
  !!      endif
   !
   !  global mean tempertures (tbar) and specific humid (qbar)
   !

   !$acc host_data use_device(plmid, plnow, divmid, divnow, vormid, vornow, &
   !$acc& temmid, temnow)
   istat = cudaMemcpyAsync(plmid, plnow, size(plmid), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(divmid, divnow, size(divmid), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(vormid, vornow, size(vormid), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(temmid, temnow, size(temmid), cudaMemcpyDeviceToDevice, stream)
   !$acc end host_data

   ! plten = 0.
   ! divten = 0.
   ! vorten = 0.
   ! temten = 0.
   ! hldten = 0.

   !     estimate all field at t+dt/2
   itt = min(itimestep, 2)
   !
   !$acc host_data use_device(up, um, ut, &
   !$acc& vp, vm, vt, ttp, tm, tt, rdivm, rdiv, qm, qt, ptp, ptm, pt)
   istat = cudaMemcpyAsync(up, ut, size(up), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(um, ut, size(um), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(vp, vt, size(vp), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(vm, vt, size(vm), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(ttp, tt, size(ttp), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(tm, tt, size(tm), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(rdivm, rdiv, size(rdivm), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(qm, qt, size(qm), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(ptp, pt, size(ptp), cudaMemcpyDeviceToDevice, stream)
   istat = cudaMemcpyAsync(ptm, pt, size(ptm), cudaMemcpyDeviceToDevice, stream)
   !$acc end host_data

   !
   ! Transfer Spectral to Gridpoint for u,v,t,q,ps at n-1
   !
   !      call transr1(jtrun,jtmax,nx,my,my_max,poly,plmid,ptm,nsizey)
   !      call tranuv(jtrun,jtmax,nx,my,my_max,levp,onocos,wcfac,wdfac    &
   !                 ,poly,dpoly,vormid,divmid,um,vm,nsizey)
   !      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,divmid,cc,1,nsizey)
   !      call ujoinsr(cc,rdivm,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
   !      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,temmid,cc,1,nsizey)
   !      call ujoinsr(cc,tm,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
   forward = .true.
   fwd = .true.

   do itt = 1, itter
      ! call nvtxStartRange("ndsl", mod(itt,2)+1)
      !
      ! for Semi-Lagrangian advection
      !
      !$acc host_data use_device(pdot, vdmerd, vdzonl, &
      !$acc& vdmerdrp, vdzonlrp, ddtemp, pten, deldm)
      istat = cudaMemsetAsync(pdot, 0., size(pdot), stream)
      ! istat = cudaMemsetAsync(vdmerd,   0., size(vdmerd), stream)
      ! istat = cudaMemsetAsync(vdzonl,   0., size(vdzonl), stream)
      ! istat = cudaMemsetAsync(vdmerdr,  0., size(vdmerdr), stream)
      ! istat = cudaMemsetAsync(vdzonlr,  0., size(vdzonlr), stream)
      istat = cudaMemsetAsync(vdmerdrp, 0., size(vdmerdrp), stream)
      istat = cudaMemsetAsync(vdzonlrp, 0., size(vdzonlrp), stream)
      ! istat = cudaMemsetAsync(ddtemp,   0., size(ddtemp), stream)
      istat = cudaMemsetAsync(pten, 0., size(pten), stream)

      istat = cudaMemsetAsync(deldm, 0., size(deldm), stream)
      !$acc end host_data
      !
      !     advet grid non-linear forcing from t-dt/2 to t+dt/2 via NDSL advection
      !
      !$acc wait(async_id)
      call ndslfv_monoadvh_fgnl_gpu(vdzonl, vdmerd, ddtemp, &
                                    ut, vt, tt, um, vm, dtahi, xy, 3, forward)

      ! ************************************************************
      !     calculate vertical velocity at mid-point
      ! ************************************************************
      !$acc parallel loop gang private(j) async(async_id)
      do jj = 1, jlistnum
         j = jlist1(jj)
         ! nxj = nxdef_2d(j)
         ! ************************************************************
         !   new p**capa quantities were computed in previous diabat call
         ! ************************************************************
        !! input: nxjp, nxp, lev, ptop, sigma, ptm
        !! output: pk, pk2, plt
         call prexp_hybrid_cwb_gpu(nxjp(j), nxp, lev, ptop, sigma, ptm(1, jj), &
                                   pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
         !
        !! input: um, vm, rdivm, tm, qt,
        !!        ptm, dtpl, dlpl, pk ,pk2, pten
        !! output: phi, diveng, vdmerdg, vdzonlg, deldm, sdpbl, sd, pdot, vvel
         call gridnl_hybrid_ndsl_gpu(nxjp(j), nxp, lev, ncld &
                                     , cp, radsq, um(1, 1, jj), vm(1, 1, jj), rdivm(1, 1, jj), tm(1, 1, jj) &
                                     , qt(1, 1, jj), phi(1, 1, jj), ptm(1, jj), dtpl(1, jj), dlpl(1, jj), sinl(j) &
                                     , pk(1, 1, jj), pk2(1, 1, jj), dsigma, sigma, onocos(j), cor(j) &
                                     , diveng(1, 1, jj), vdmerdg(1, 1, jj), vdzonlg(1, 1, jj), pten(1, 1, jj) &
                                     !ndy        , deldm(1,jj),sdpbl(1,jj),sd(1,1,jj),pdot(1,1,jj),sgeo(1,jj),2 )
                                     , deldm(1, jj), sdpbl(1, jj), sd(1, 1, jj), pdot(1, 1, jj), vvel(1, 1, jj) &
                                     , sgeo(1, jj))
         !
      end do !jj = 1,jlistnum

      call joinrs_gpu(cc, diveng, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
      call tranrs_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, weight, cc &
                      , hldten, 1, nsizey)
      call trngra3_gpu(jtrun, jtmax, nx, levp, my, my_max, cim, poly, dpoly &
                       , hldten, dlphi, dtphi, nsizey)

      !
      !$acc parallel loop collapse(2) async(async_id) &
      !$acc& private(j,nxj)
      do jj = 1, jlistnum
         do k = 1, lev
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            !$acc loop vector
            do i = 1, nxj
               vdmerdrp(i, k, jj) = vdmerdg(i, k, jj) - dtphi(i, k, jj)/radsq/onocos(j)
               vdzonlrp(i, k, jj) = vdzonlg(i, k, jj) - dlphi(i, k, jj)/radsq
            end do
         end do
      end do !jj = 1,jlistnum

      call ndslfv_update(nxjp, vdzonl, vdmerd, vdzonlrp, vdmerdrp, dtahi, forward)

      call ndslfv_monoadvv_fgnl(vdzonl, vdmerd, ddtemp, pdot, ptm &
                                , nxjp, dtahi, 3, forward)

      !CWB2021 ndsl single precision test

      call mpe2d_unify_nx_gpu(ww1, deldm)
      call tranrs1_gpu(jtrun, jtmax, nx, my, my_max, poly, weight, ww1 &
                       , plten, nsizey)

      if (forward) then
         !$acc parallel loop collapse(3) private(j, nxj) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               do i = 1, nxp
                  j = jlist1(jj)
                  nxj = nxdef_2d(j)
                  if (i .le. nxj) then
                     vdzonl(i, k, jj) = (vdzonl(i, k, jj) - um(i, k, jj))/dtahi
                     vdmerd(i, k, jj) = (vm(i, k, jj) - vdmerd(i, k, jj))/dtahi
                     ddtemp(i, k, jj) = (ddtemp(i, k, jj) - tm(i, k, jj))/dtahi
                  end if
               end do
            end do
         end do
      else
         !$acc parallel loop collapse(3) private(j, nxj) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               do i = 1, nxp
                  j = jlist1(jj)
                  nxj = nxdef_2d(j)
                  if (i .le. nxj) then
                     vdzonl(i, k, jj) = (vdzonl(i, k, jj) - ut(i, k, jj))/dtah
                     vdmerd(i, k, jj) = (vt(i, k, jj) - vdmerd(i, k, jj))/dtah
                     ddtemp(i, k, jj) = (ddtemp(i, k, jj) - tt(i, k, jj))/dtah
                  end if
               end do
            end do
         end do
      end if

      !
     !! input: ddtemp
     !! output: temten, divten, vorten
      call joinrs_gpu(cc, ddtemp, dummy, dummy, dummy, nx, my_max, lev &
                      , jlistnum, 1, 1)
      call tranrs_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, weight, cc &
                      , temten, 1, nsizey)
      call rstrandz_gpu(jtrun, jtmax, nx, my, my_max, levp, vdmerd, vdzonl &
                        , weight, cim, onocos, poly, dpoly, divten, vorten, nsizey)

      if (forward) then
         if (lsimpl) then
           !! input: temmid, divmid, plmid, alphax
           !! output: temten, divten, plten
            call siimpl_gpu(jtrun, jtmax, lev, dtahi, ptmeans, dsigma, spalm, eps4, eigval &
                            , evecin, evectr, arrhyd, arsddt, temmid, divmid, plmid &
                            , temmid, divmid, plmid, temten, divten, plten, alphax)
         end if

      else
         if (lsimpl) then
           !! input: temnow, divnow, plnow, temmid, divmid, plmid, alphax
           !! output: temten, divten, plten
            call siimpl_gpu(jtrun, jtmax, lev, dtah, ptmeans, dsigma, spalm, eps4, eigval &
                            , evecin, evectr, arrhyd, arsddt, temnow, divnow, plnow &
                            , temmid, divmid, plmid, temten, divten, plten, alphax)
         end if
      end if

      mlst = ilist(1)
      if (mlst .ne. 0) then
         !$acc kernels copyin(mlst) async(async_id)
         plten(1, mlst, 1) = 0.0
         plten(1, mlst, 2) = 0.0
         !$acc end kernels
      end if

      !$acc parallel loop collapse(3) private(mf) async(async_id)
      do m = 1, mlistnum
         do i = 1, 2
            do k = 1, levp
               mf = mlist(m)
               if (mf .eq. 1) then
                  divten(k, i, 1, m) = 0.0
                  vorten(k, i, 1, m) = 0.0
               end if
            end do
         end do
      end do

      if (forward) then
         !$acc parallel loop collapse(4) private(mf) async(async_id)
         do m = 1, mlistnum
            do n = 1, jtrun
               do i = 1, 2
                  do k = 1, levp
                     mf = mlist(m)
                     if (n .ge. mf) then
                        vormid(k, i, n, m) = dtahi*vorten(k, i, n, m) + vormid(k, i, n, m)
                        divmid(k, i, n, m) = dtahi*divten(k, i, n, m) + divmid(k, i, n, m)
                        temmid(k, i, n, m) = dtahi*temten(k, i, n, m) + temmid(k, i, n, m)
                     end if
                  end do
               end do
            end do
         end do
         !$acc parallel loop collapse(3) private(mf) async(async_id)
         do i = 1, 2
            do m = 1, mlistnum
               do n = 1, jtrun
                  mf = mlist(m)
                  if (n .ge. mf) then
                     plmid(n, m, i) = dtahi*plten(n, m, i) + plmid(n, m, i)
                  end if
               end do
            end do
         end do
         hfiltm = mwhd*hfiltx
        !! input: hfiltm, rad, trefs, um, vm,
        !! output: vormid, divmid, temmid
         !$acc wait(async_id)
         !$acc update device(hfiltm) async(async_id)
         call whdiffu_gpu(dtahi, my, my_max, nx, jtrun, jtmax, lev, ncld &
                          , hfiltm, rad, cosl, um, vm, vormid, divmid, temmid &
                          , eps4, trefs)
      else
         !$acc parallel loop collapse(4) private(mf) async(async_id)
         do m = 1, mlistnum
            do n = 1, jtrun
               do i = 1, 2
                  do k = 1, levp
                     mf = mlist(m)
                     if (n .ge. mf) then
                        vormid(k, i, n, m) = dtah*vorten(k, i, n, m) + vornow(k, i, n, m)
                        divmid(k, i, n, m) = dtah*divten(k, i, n, m) + divnow(k, i, n, m)
                        temmid(k, i, n, m) = dtah*temten(k, i, n, m) + temnow(k, i, n, m)
                     end if
                  end do
               end do
            end do
         end do
         !$acc parallel loop collapse(3) private(mf) async(async_id)
         do i = 1, 2
            do m = 1, mlistnum
               do n = 1, jtrun
                  mf = mlist(m)
                  if (n .ge. mf) then
                     plmid(n, m, i) = dtah*plten(n, m, i) + plnow(n, m, i)
                  end if
               end do
            end do
         end do

         hfiltm = mwhd*hfiltx
        !! input: hfiltm, rad, trefs, um, vm,
        !! output: vormid, divmid, temmid
         !$acc wait(async_id)
         !$acc update device(hfiltm) async(async_id)
         call whdiffu_gpu(dtah, my, my_max, nx, jtrun, jtmax, lev, ncld &
                          , hfiltm, rad, cosl, um, vm, vormid, divmid, temmid &
                          , eps4, trefs)
      end if
      ! ************************************************************
      !     update all new wind field at mid-point
      ! ************************************************************
     !! input: vormid, divmid
     !! output: um, vm
      call tranuv_gpu(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac &
                      , poly, dpoly, vormid, divmid, um, vm, nsizey)
     !! input: divmid, temmid
     !! output: rdivm, tm
      call transr_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, divmid, cc, 1, nsizey)
      call ujoinsr_gpu(cc, rdivm, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
      call transr_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, temmid, cc, 1, nsizey)
      call ujoinsr_gpu(cc, tm, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
     !! input: plmid
     !! output: ptm, dlpl, dtpl
      call transr1_gpu(jtrun, jtmax, nx, my, my_max, poly, plmid, ptm, nsizey)
      call trngra_gpu(jtrun, jtmax, nx, my, my_max, cim, poly, dpoly, plmid &
                      , dlpl, dtpl, nsizey)

      !$acc wait(async_id)
      ! call nvtxEndRange
      forward = .false.
   end do ! do itt=1,itter

   ! call nvtxStartRange("ndsl-p3", 3)
   !
   !  the gaussian quadrature loop for spectral tendencies.  subroutine
   !  'gridnl' is called for companion mirror image gaussian latitudes
   !  these non-linear contributions are then combined in 'rstran' using
   !  the symmetry properties of the spherical harmonics
   !
   forward = .false.
   !$acc parallel loop gang private(j) async(async_id)
   do jj = 1, jlistnum
      j = jlist1(jj)
      ! nxj = nxdef_2d(j)
      ! ************************************************************
      !   new p**capa quantities were computed in previous diabat call
      ! ************************************************************
     !! input: ptm
     !! output: pk, pk2, plt
      call prexp_hybrid_cwb_gpu(nxjp(j), nxp, lev, ptop, sigma, ptm(1, jj), &
                                pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
      !
      !       Calculate Vertical velocity & Stream Functions
      !
     !!        call gridnl_hybrid_ndsl_2tl (nxjp(j),nxp,lev,ncld              &
     !! input: um, vm, rdivm, tm, qt,
     !!        ptm, dtpl, dlpl, pk ,pk2, pten
     !! output: phi, diveng, vdmerdg, vdzonlg, deldm, sdpbl, sd, pdot, vvel
      call gridnl_hybrid_ndsl_gpu(nxjp(j), nxp, lev, ncld &
                                  , cp, radsq, um(1, 1, jj), vm(1, 1, jj), rdivm(1, 1, jj), tm(1, 1, jj) &
                                  , qm(1, 1, jj), phi(1, 1, jj), ptm(1, jj), dtpl(1, jj), dlpl(1, jj), sinl(j) &
                                  , pk(1, 1, jj), pk2(1, 1, jj), dsigma, sigma, onocos(j), cor(j) &
                                  , diveng(1, 1, jj), vdmerdg(1, 1, jj), vdzonlg(1, 1, jj), pten(1, 1, jj) &
                                  , deldm(1, jj), sdpbl(1, jj), sd(1, 1, jj), pdot(1, 1, jj), vvel(1, 1, jj) &
                                  , sgeo(1, jj))
   end do !jj = 1,jlistnum
   !

   call joinrs_gpu(cc, diveng, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
   call tranrs_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, weight, cc &
                   , hldten, 1, nsizey)
   call trngra3_gpu(jtrun, jtmax, nx, levp, my, my_max, cim, poly, dpoly &
                    , hldten, dlphi, dtphi, nsizey)

   !$acc parallel loop collapse(2) async(async_id) &
   !$acc& private(j,nxj)
   do jj = 1, jlistnum
      do k = 1, lev
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         !$acc loop vector
         do i = 1, nxj
            vdmerdg(i, k, jj) = vdmerdg(i, k, jj) - dtphi(i, k, jj)/radsq/onocos(j)
            vdzonlg(i, k, jj) = vdzonlg(i, k, jj) - dlphi(i, k, jj)/radsq
         end do
      end do
   end do !jj = 1,jlistnum
   !$acc wait(async_id)
   ! ****************************************
   ! Semi-Lagrangian
   !       Horizontal Advection
   ! ****************************************
   call ndslfv_monoadvh_gpu(tt, pten, ut, vt, qt, qm, &
                            um, vm, dtah, xy, forward)
   ! ****************************************
   !       update all horizontal informations
   ! ****************************************
   call ndslfv_update(nxjp, ut, vt, vdzonlg, vdmerdg, dtah, forward)
   ! ****************************************
   !       Vertical Advection
   ! ****************************************
   call ndslfv_monoadvv(tt, qt, ut, vt, pdot, ptm, nxjp, dtah, forward)

   !CWB2021 ndsl single precision test

   call mpe2d_unify_nx_gpu(ww1, deldm)
   call tranrs1_gpu(jtrun, jtmax, nx, my, my_max, poly, weight, ww1 &
                    , plten, nsizey)

   !$acc parallel loop collapse(3) private(mf) async(async_id)
   do i = 1, 2
      do m = 1, mlistnum
         do n = 1, jtrun
            mf = mlist(m)
            if (n .ge. mf) then
               pltemp(n, m, i) = dta*plten(n, m, i) + plnow(n, m, i)
            end if
         end do
      end do
   end do

   call transr1_gpu(jtrun, jtmax, nx, my, my_max, poly, pltemp, pt, nsizey)

   !mass conservation
   !$acc update self(pt, qt) async(async_id)
   !$acc wait(async_id)
  !! input: pt, qt, dpprt
  !! outout: pdry, pcorr
   call ptotc(pdry, dpprt)
   pcorr = (pdryi - pdry)*sqrt(2.)
   !$acc update device(pcorr) async(async_id)
   !
   if (two_loop) then
      !
      !  after phyical parameterization,transform grid point u,v,t,q to
      !  spectrum
      !
      !$acc parallel loop collapse(3) private(j, nxj) async(async_id)
      do jj = 1, jlistnum
         do k = 1, lev
            do i = 1, nxp
               j = jlist1(jj)
               nxj = nxdef_2d(j)
               if (i .le. nxj) then
                  vdzonl(i, k, jj) = (ut(i, k, jj) - up(i, k, jj))/dta
                  vdmerd(i, k, jj) = (vp(i, k, jj) - vt(i, k, jj))/dta
                  ddtemp(i, k, jj) = (tt(i, k, jj) - ttp(i, k, jj))/dta
               end if
            end do
         end do
      end do
      !
      call joinrs_gpu(cc, ddtemp, dummy, dummy, dummy, nx, my_max, lev &
                      , jlistnum, 1, 1)
      call tranrs_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, weight, cc &
                      , temten, 1, nsizey)
      call rstrandz_gpu(jtrun, jtmax, nx, my, my_max, levp, vdmerd, vdzonl &
                        , weight, cim, onocos, poly, dpoly, divten, vorten, nsizey)

      !        if ( mass_dp ) then
      ! sureface pressure global mean correction
      mlst = ilist(1)
      if (mlst .ne. 0) then
         !$acc kernels copyin(mlst, pcorr) async(async_id)
         plten(1, mlst, 1) = plten(1, mlst, 1) + pcorr/dta
         plten(1, mlst, 2) = plten(1, mlst, 2) + pcorr/dta
         !$acc end kernels
      end if
      !        endif
      !
      if (lsimpl) then
         !
         !   apply semi-implicit adjustemts to above explicitly computed
         !   tendencies to stablize integration for long time steps
         !
         !CWB2021
         call siimpl_gpu(jtrun, jtmax, lev, dta, ptmeans, dsigma, spalm, eps4, eigval &
                         , evecin, evectr, arrhyd, arsddt, temnow, divnow, plnow &
                         , temmid, divmid, plmid, temten, divten, plten, alphax)
      end if
      !
      !  zero out global mean tendencies for divergence, vorticity, and
      !  terrain pressure to ensure consistency with gauss's theorem.
      !
      !        if ( .not. mass_dp ) then
      !          mlst=ilist(1)
      !          if(mlst .ne. 0) then
      !            plten(1,mlst,1) = 0.0
      !            plten(1,mlst,2) = 0.0
      !          endif
      !        endif

      !$acc parallel loop collapse(3) private(mf) async(async_id)
      do m = 1, mlistnum
         do i = 1, 2
            do k = 1, levp
               mf = mlist(m)
               if (mf .eq. 1) then
                  divten(k, i, 1, m) = 0.0
                  vorten(k, i, 1, m) = 0.0
               end if
            end do
         end do
      end do
      !
      call transr1_gpu(jtrun, jtmax, nx, my, my_max, poly, plten, ptend, nsizey)
      !$acc enter data create(wkmf) async(async_id)
      !$acc parallel loop gang private(mf, wkmf_local) async(async_id)
      do m = 1, mlistnum
         mf = mlist(m)
         wkmf_local = 0.0
         !$acc loop vector reduction(+:wkmf_local)
         do n = 2, jtrun
            if ((n .ge. mf)) then
               wkmf_local = wkmf_local + plten(n, m, 1)**2 + plten(n, m, 2)**2
            end if
         end do
         wkmf(mf) = wkmf_local
      end do
      !$acc exit data copyout(wkmf) async(async_id)
      !$acc wait(async_id)
      call mpe_unify(wkmf, 1, jtrun, 3, mpe_double)
      sptend = 0.0
      do mf = 1, jtrun
         sptend = sptend + wkmf(mf)
      end do
      !
      sptend = sqrt(0.5*sptend)*3600.0
      if (myrank .eq. 0) &
         print *, '[GPU]  surf pres tend rms =', sptend, ' mb/hrs'
      !
      !  take a time step
      !
      !$acc parallel loop collapse(4) private(mf) async(async_id)
      do m = 1, mlistnum
         do n = 1, jtrun
            do i = 1, 2
               do k = 1, levp
                  mf = mlist(m)
                  if (n .ge. mf) then
                     vorold(k, i, n, m) = vornow(k, i, n, m)
                     divold(k, i, n, m) = divnow(k, i, n, m)
                     temold(k, i, n, m) = temnow(k, i, n, m)
                     vornow(k, i, n, m) = dta*vorten(k, i, n, m) + vorold(k, i, n, m)
                     divnow(k, i, n, m) = dta*divten(k, i, n, m) + divold(k, i, n, m)
                     temnow(k, i, n, m) = dta*temten(k, i, n, m) + temold(k, i, n, m)
                  end if
               end do
            end do
         end do
      end do
      !
      !$acc parallel loop collapse(3) private(mf) async(async_id)
      do i = 1, 2
         do m = 1, mlistnum
            do n = 1, jtrun
               mf = mlist(m)
               if (n .ge. mf) then
                  plold(n, m, i) = plnow(n, m, i)
                  plnow(n, m, i) = dta*plten(n, m, i) + plold(n, m, i)
               end if
            end do
         end do
      end do
      !
      !
      if (hdiff) then
         !$acc update device(hfiltx) async(async_id)
         call hdiffu_gpu(dta, my, my_max, nx, jtrun, jtmax, lev, ncld &
                         , hfiltx, rad, cosl, um, vm, vornow, divnow, temnow &
                         , eps4, trefs)
      end if
      !
      !  for physical parameterization,output spectrum u,v,t,q to grid point
      !
     !! input: temnow, vornow, divnow, plnow
     !! output: tt, ut, vt, pt
      call transr_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, temnow, cc, 1, nsizey)
      call ujoinsr_gpu(cc, tt, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
      call tranuv_gpu(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac &
                      , poly, dpoly, vornow, divnow, ut, vt, nsizey)
      call transr1_gpu(jtrun, jtmax, nx, my, my_max, poly, plnow, pt, nsizey)!
   end if ! two_loop
   !$acc update self(rdivm, ptm, dlpl, dtpl) async(async_id)
   !$acc update self(um, vm, tm) async(async_id)
   !$acc update self(divten, vorten, temten, plten) async(async_id)
   !$acc update self(vormid, divmid, temmid, plmid) async(async_id)
   !$acc update self(vorold, divold, temold, plold) async(async_id)
   !$acc update self(vornow, divnow, temnow, plnow) async(async_id)
   !$acc update self(phi, pdot) async(async_id)
   !$acc update self(vvel, sdpbl, sd) async(async_id)
   !$acc update self(pk ,pk2, plt) async(async_id)

   !$acc update self(tt, ut, vt, qt) async(async_id)
   !$acc update self(pt, pltemp, ptend) async(async_id)
   !$acc update self(ptp, up, vp, ttp, qm) async(async_id)
   !$acc wait(async_id)
   ! call nvtxEndRange
!!!! ---------------------------------------- !!!!
!!!! << END GPU >>                            !!!!
!!!! ---------------------------------------- !!!!

   !  stochastic_physics
   call run_stochastic_physics()
   !
   !  for physical parameterization,output spectrum u,v,t,q to grid point
   !

   if (yesdia) then
      ! call nvtxStartRange("diabat")
      qp(:, :, :) = qt(:, :, :)
      call diabat(fwd, docup, dodry, dolsp, dopbl, dorad, doshl, dograv, tofd &
                  , nx, my, my_max, lev, ncld, nmcup, nmpbl, nmland, nmshl, cgw &
                  !                         , idg, jdg, ldiag, dtx, tau, hours, julian, year, yrd &
                  , idg, jdg, ldiag, dtx, tau, hours, year &
                  , frad, ozon, njump, itypbl, ktcup, ktpbl, ktshl, grav &
                  , rgas, cp, stbo, s0, evaprh, hltm, ptop, sigma, dsigma, il, ib &
                  , cof, xlat, xlon, sgeo, z0, alb, land, ocean, ice &
                  , snr, tg, tgclim, curate, plcl, cumtop, totalp, raintot, raincu &
                  , rainlp, raincu6, rainlp6, raincu3, rainlp3, raincu1, rainlp1 &
                  , hflux, qflux, ustar, tstar, qstar, e &
                  , eps, o3l, dtrad, ss, rs, plt, pk, pk2 &
                  , ptp, up, vp, ttp, qm &
                  , pt, ut, vt, tt, qt &
                  , gwclim, tice, hice, qgini, thdai, tengi &
                  , acld, std, asol, olr, drag, ugws, vgws &
                  , sdpbl, t2, q2, rh2, rh10, u10, v10, gfx &
                  , fm, fh, fm10, fh2, srflag &
                  , rld, km_soil, smc, stc, canopy, runoff &
                  , sigmaf, istyp, ivegtyp, wltsmc, refsmc, maxsmc, dfkt, xktk, dfk &
                  , ftp, fqp, fpsp, ftp1, fqp1, fpsp1, deltaq, cnvwr, cnvcr, pdot &
                  , shdmax, shdmin, snoalb &
                  , slopetyp, sld, slc, zice, cice, xtice, sncover, sndepth &
                  , ctot, chig, cmid, clow, hpbl, asl, atl, cosz &
                  , nmgwor, nmgwcv, hprime_b, mtnvar, docgrav, nmmiph &
                  !--------------------------------------------------------------------------------
                  , fusl, fdsl, fuir, fdir &
                  , fuslr, fdslr, fuirr, fdirr &
                  , asl_clr, atl_clr, clds &
                  , ss_clr, rs_clr, asol_clr, olr_clr, sld_clr, rld_clr &
                  , alvsf, alvwf, alnsf, alnwf, facsf, facwf &
                  , idtg, doo3l, nfxr, sfalb, sfemis, isot, ivegsrc &
                  , itimestep, lrun_sitvdiff, ic_sit &
                  !xb110>
                  !byl                      , rmr,smr,flash)
                  , flash, tsflw, vvel, totallp)
      !xb110<
      !--------------------------------------------------------------------------------
      !
      ! add reynolds stress
      !
      call rayleifr(nx, my, my_max, lev, rad, cosl, dt, ut, vt)

      if (two_loop) then
         ! adjustmen of surface pressure, virtual potential
         ! temperature and all tracers
         if (mass_dp) call adjptq(dta, plnow, pltemp)
         call joinrs(cc, tt, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
         call tranrs(jtrun, jtmax, nx, my, my_max, levp, poly, weight, cc &
                     , temnow, 1, nsizey)
         ! SKEB process
         if (doskeb) call skebest(um, vm)

         call trandv(jtrun, jtmax, nx, my, my_max, lev, ut, vt, weight, cim &
                     , onocos, poly, dpoly, vornow, divnow, nsizey)
      else
         ! adjustmen of surface pressure, virtual potential
         ! temperature and all tracers for one loop
         if (mass_dp) call adjptq(dta, pltemp, plten)

      end if ! two_loop

      ! call nvtxEndRange
   end if    ! end of (yesdia)

   !CWB2021
   itimestep = itimestep + 1
   !        if ( mod(itimestep,2) .eq. 0 ) xy = -1 * xy
   xy = -1*xy
   if (.not. two_loop) then
      !
      !  after phyical parameterization,transform grid point u,v,t,q to
      !  spectrum
      !
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do k = 1, lev
            do i = 1, nxj
               vdzonl(i, k, jj) = (ut(i, k, jj) - up(i, k, jj))/dta
               vdmerd(i, k, jj) = (vp(i, k, jj) - vt(i, k, jj))/dta
               ddtemp(i, k, jj) = (tt(i, k, jj) - ttp(i, k, jj))/dta
            end do
         end do
      end do
      !
      call joinrs(cc, ddtemp, dummy, dummy, dummy, nx, my_max, lev &
                  , jlistnum, 1, 1)
      call tranrs(jtrun, jtmax, nx, my, my_max, levp, poly, weight, cc &
                  , temten, 1, nsizey)
      call rstrandz(jtrun, jtmax, nx, my, my_max, levp, vdmerd, vdzonl &
                    , weight, cim, onocos, poly, dpoly, divten, vorten, nsizey)

      !        if ( mass_dp ) then
      ! sureface pressure global mean correction
      mlst = ilist(1)
      if (mlst .ne. 0) then
         plten(1, mlst, 1) = plten(1, mlst, 1) + pcorr/dta
         plten(1, mlst, 2) = plten(1, mlst, 2) + pcorr/dta
      end if
      !        endif
      !
      if (lsimpl) then
         !
         !   apply semi-implicit adjustemts to above explicitly computed
         !   tendencies to stablize integration for long time steps
         !
         !CWB2021
         call siimpl(jtrun, jtmax, lev, dta, ptmeans, dsigma, spalm, eps4, eigval &
                     , evecin, evectr, arrhyd, arsddt, temnow, divnow, plnow &
                     , temmid, divmid, plmid, temten, divten, plten, alphax)
      end if
      !
      !  zero out global mean tendencies for divergence, vorticity, and
      !  terrain pressure to ensure consistency with gauss's theorem.
      !
      !        if ( .not. mass_dp ) then
      !          mlst=ilist(1)
      !          if(mlst .ne. 0) then
      !            plten(1,mlst,1) = 0.0
      !            plten(1,mlst,2) = 0.0
      !          endif
      !        endif
      do m = 1, mlistnum
         mf = mlist(m)
         if (mf .eq. 1) then
            do i = 1, 2
               do k = 1, levp
                  divten(k, i, 1, m) = 0.0
                  vorten(k, i, 1, m) = 0.0
               end do
            end do
         end if
      end do
      !
      call transr1(jtrun, jtmax, nx, my, my_max, poly, plten, ptend, nsizey)
      do mf = 1, jtrun
         wkmf(mf) = 0.
      end do
      sptend = 0.0
      do m = 1, mlistnum
         mf = mlist(m)
         do n = mf, jtrun
            if (n .ne. 1) then
               wkmf(mf) = wkmf(mf) + plten(n, m, 1)**2 + plten(n, m, 2)**2
            end if
         end do
      end do
      call mpe_unify(wkmf, 1, jtrun, 3, mpe_double)
      do mf = 1, jtrun
         sptend = sptend + wkmf(mf)
      end do
      !
      sptend = sqrt(0.5*sptend)*3600.0
      if (myrank .eq. 0) &
         print *, '  surf pres tend rms =', sptend, ' mb/hrs'
      !
      !  take a time step
      !
      do m = 1, mlistnum
         mf = mlist(m)
         do n = mf, jtrun
            do i = 1, 2
               do k = 1, levp
                  vorold(k, i, n, m) = vornow(k, i, n, m)
                  divold(k, i, n, m) = divnow(k, i, n, m)
                  temold(k, i, n, m) = temnow(k, i, n, m)
                  vornow(k, i, n, m) = dta*vorten(k, i, n, m) + vorold(k, i, n, m)
                  divnow(k, i, n, m) = dta*divten(k, i, n, m) + divold(k, i, n, m)
                  temnow(k, i, n, m) = dta*temten(k, i, n, m) + temold(k, i, n, m)
               end do
            end do
         end do
      end do
      !
      do i = 1, 2
         do m = 1, mlistnum
            mf = mlist(m)
            do n = mf, jtrun
               plold(n, m, i) = plnow(n, m, i)
               plnow(n, m, i) = dta*plten(n, m, i) + plold(n, m, i)
            end do
         end do
      end do
      !
      !
      if (hdiff) call hdiffu(dta, my, my_max, nx, jtrun, jtmax, lev, ncld &
                             , hfiltx, rad, cosl, um, vm, vornow, divnow, temnow &
                             , eps4, trefs)
      ! SKEB process
      if (doskeb) then
         call tranuv(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac &
                     , poly, dpoly, vornow, divnow, ut, vt, nsizey)

         call skebest(um, vm)

         ! compute vorticity and divergence from u and v
         call trandv(jtrun, jtmax, nx, my, my_max, lev, ut, vt, weight, cim &
                     , onocos, poly, dpoly, vornow, divnow, nsizey)
      end if
   end if ! .not. two_loop
   !
   ! update tg, dtsea/dt (W00100)
   !
   if (ldailyFCTsst .OR. ldailyFCTicesndpt .OR. dailyClm_option .ge. 1) then
      if (tau .ge. 24.) then
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do ii = 1, nxj
               i = nxjstart(j) + ii - 1
               if (ocean(ii, jj)) then
                  tseanew(ii, jj) = dta*dtseadt(ii, jj) + tseaold(ii, jj)
                  !                  tseaold(ii,jj)=tseanow(ii,jj) + tfilt*(tseaold(ii,jj)     &
                  !                                 -2.0*tseanow(ii,jj)+tseanew(ii,jj) )
                  tseaold(ii, jj) = tseanew(ii, jj)
                  tseanow(ii, jj) = tseanew(ii, jj)

                  dtaup = mod(tau + 0.001, updatetg)
                  if (dtaup .lt. dtx_tau) then
                     tg(ii, jj) = tseanow(ii, jj)
                  end if
               end if !end if(ocean)
            end do
         end do
      end if
      CALL read_dailyFCT(idtg, tau, dt, tg, cice, sndepth, xlon, xlat, ocean)
   end if
   !
   ! accumulate some flux every time step to output point (24 hour)
   ! 1994 11 11
   !
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do i = 1, nxj
         hf24(i, jj) = hf24(i, jj) + hflux(i, jj)*dtx
         qf24(i, jj) = qf24(i, jj) + qflux(i, jj)*dtx
         ss24(i, jj) = ss24(i, jj) + ss(i, jj)*dtx
         rs24(i, jj) = rs24(i, jj) + rs(i, jj)*dtx
         asol24(i, jj) = asol24(i, jj) + asol(i, jj)*dtx
         olr24(i, jj) = olr24(i, jj) + olr(i, jj)*dtx
         rain24(i, jj) = rain24(i, jj) + totalp(i, jj)
         rainlp24(i, jj) = rainlp24(i, jj) + totallp(i, jj)
         flash24(i, jj) = flash24(i, jj) + flash(i, jj)*dtx
      end do
   end do
   dt24 = dt24 + dtx
   !
   !
   !  from updated spectral variables compute corresponding grid
   !  point fields
   !
  !!      call joinsr(wss3,vornow,divnow,temnow,dummy,jtrun,jtmax,levp &
  !!                 ,mlistnum,3,1)
  !!      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,wss3,cc3,3,nsizey)
  !!      call ujoinsr(cc3,rvor,rdiv,tt,dummy,nx,my_max,lev,jlistnum,3,1)

  !!      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,vornow,cc,1,nsizey)
  !!      call ujoinsr(cc,rvor,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
  !! input: divnow, temnow, plnow, vornow, divnow
  !! output: rdiv, tt, pt, pk, pk2, plt, dlpl, dtpl, ut, vt
   !$acc update device(qt) async(async_id)
   !$acc update device(divnow, plnow, vornow, temnow) async(async_id)
   !$acc update device(sigma, sgeo) async(async_id)
   call transr_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, divnow, cc, 1, nsizey)
   call ujoinsr_gpu(cc, rdiv, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)

   call transr_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, temnow, cc, 1, nsizey)
   call ujoinsr_gpu(cc, tt, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)

   call transr1_gpu(jtrun, jtmax, nx, my, my_max, poly, plnow, pt, nsizey)
   !
   !   computing new p**capa quantities
   !
   !$acc parallel loop gang private(j) async(async_id)
   do jj = 1, jlistnum
      j = jlist1(jj)
      call prexp_hybrid_cwb_gpu(nxjp(j), nxp, lev, ptop, sigma, pt(1, jj), &
                                pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
   end do
   !
   !  zonal and meridional gradients of terrain pressure
   !
   call trngra_gpu(jtrun, jtmax, nx, my, my_max, cim, poly, dpoly, plnow &
                   , dlpl, dtpl, nsizey)
   !
   !  velocity components
   !
   call tranuv_gpu(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac &
                   , poly, dpoly, vornow, divnow, ut, vt, nsizey)
   !$acc update self(rdiv, tt, pt, pk, pk2, plt, &
   !$acc& dlpl, dtpl, ut, vt) async(async_id)
   !$acc wait(async_id)
   !
   !
   !  detact instability occure or not
   !
   if (tau .gt. 12.) then
      !        sptendmax2=0.409
      !        sptendmax1=0.379
      if (sptend .le. sptendmax1) n_stable = n_stable + 1
      if (sptend .gt. sptendmax2) n_unstable = n_unstable + 1
      if (myrank .eq. 0) print *, 'n_stable=', n_stable, ' n_unstable=' &
         , n_unstable
      if (mod(tau + 0.001, 1.) .lt. dtx_tau) then
         if (n_stable .gt. nc_stable) then
            hfiltx = 0.8*hfilt
            alphax = alpha
            if (myrank .eq. 0) print *, '** stable change hfilt=', hfiltx, &
               ' and keep alpha=', alphax
         else if (n_unstable .gt. nc_stable) then
            hfiltx = hfilt
            alphax = alpha + 0.05
            if (myrank .eq. 0) print *, '** unstable change hfilt=', hfiltx, &
               ' and alpha=', alphax
         else
            hfiltx = hfilt
            alphax = alpha
            if (myrank .eq. 0) print *, '** keep hfilt=', hfiltx, &
               ' and alpha=', alphax
         end if
         n_unstable = 0
         n_stable = 0
         !$acc update device(hfiltx, alphax) async(async_id)
      end if
   end if

   ! calculate Tmax Tmin @ 2m from T2
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do i = 1, nxj
         if (tmax(i, jj) .eq. 0.0) tmax(i, jj) = t2(i, jj)
         if (tmin(i, jj) .eq. 0.0) tmin(i, jj) = t2(i, jj)
         if (t2(i, jj) .gt. tmax(i, jj)) then
            tmax(i, jj) = t2(i, jj)
         elseif (t2(i, jj) .lt. tmin(i, jj)) then
            tmin(i, jj) = t2(i, jj)
         end if
      end do
   end do

   ! call nvtxEndRange
   !    ---------------------------------------------------------------
   !     check tau in hourly for output
   dtaup = mod(tau + 0.001, 1.)
   if (dtaup .lt. 0.01) then
      itau = NINT(tau)

      dtaup = mod(tau + 0.001, tauo)
      histim = (dtaup .lt. dtx_tau)
      !       if(myrank.eq.0)print *,'chkhis dtaup,tauo,dtx_tau=',dtaup,tauo,dtx_tau

      !  for tracker
      dt_trk = real(trk_intv)
      dtaup = mod(tau + 0.001, dt_trk)
      ltrack = (dtaup .lt. dtx_tau)
      !       if(myrank.eq.0)print *,'chkltr dtaup,dt_trk,dtx_tau=',dtaup,dt_trk,dtx_tau

      if (myrank == 0) call system_clock(toutsrt)
      if (io_quilting) then
         write (keydoit, '(A6,I4.4,A4,I12.12,A8)') &
            "OPEN..", itau, "....", idtg, "H...DOIT"
         ntag = ntag + 1
         call mpe_send_key(keydoit, ntag, istat)
      end if
      !
      !--- histim (start)
      if (histim .or. ltrack) then
         !itau = tau + 0.1
         itau = NINT(tau)
         if (myrank .eq. 0) print *, ' history file written at tau= ', itau
         !
         ! output gwr or gwet
         ! in new soil model, gwr did not exist
         ! but for output required,
         ! define gwet as near surfce 0.5m soil layer wetness
         ! define gwr as near surfce 0.5m soil layer moist contain
         ! gwr=gwet*gwrcc
         ! gwet(ground wetness) get from smc1*0.2+smc2*0.8
         ! saturate gwrcc set to be 20mm as original land mode setting
         !
         ! xb119 may2022 not need to run every time step ,move to histim inside
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               if (istyp(i, jj) .ne. 0) then
                  www = smc(i, 1, jj)*0.2 + smc(i, 2, jj)*0.8
                  gwet(i, jj) = (www - wltsmc(istyp(i, jj)))/ &
                                (refsmc(istyp(i, jj)) - wltsmc(istyp(i, jj)))
               else
                  gwet(i, jj) = 1.0
               end if
               gwr(i, jj) = max(0., min(1., gwet(i, jj)))*20.0
            end do
         end do
         !      call mpe_unify(gwet,nx,my,2,mpe_double)
         !      call mpe_unify(gwr,nx,my,2,mpe_double)
         !
         !
         if (wrestrt) then
            ! set write out restart at the end of integration
            if (mod(float(itau), float(itauezz)) .lt. 0.01) then
               write (ctau, 800) itau
800            format(i7.7)
               write (ccore, '(i4.4)') myrank
               !
               rfile = trim(cwbout)//'cwbout_'//ctau
               !
               if (myrank .eq. 0) open (unit=7, file=rfile, form='unformatted')
               len1 = 2*levp*jtrun*jtmax*nsize
               len2 = 2*jtrun*jtmax*nsize
               allocate (work_io(len1))
               call mpe_gather_io(work_io, vornow, len1/nsize, nsize)
               if (myrank == 0) write (7) work_io
               call mpe_gather_io(work_io, divnow, len1/nsize, nsize)
               if (myrank == 0) write (7) work_io
               call mpe_gather_io(work_io, temnow, len1/nsize, nsize)
               if (myrank == 0) write (7) work_io
               call mpe_gather_io(work_io, vorold, len1/nsize, nsize)
               if (myrank == 0) write (7) work_io
               call mpe_gather_io(work_io, divold, len1/nsize, nsize)
               if (myrank == 0) write (7) work_io
               call mpe_gather_io(work_io, temold, len1/nsize, nsize)
               if (myrank == 0) write (7) work_io
               call mpe_gather_io(work_io, trefs, len1/nsize, nsize)
               if (myrank == 0) then
                  write (7) work_io
                  call flush (7)
               end if
               deallocate (work_io)
               allocate (work_io(len2))
               call mpe_gather_io(work_io, plnow, len2/nsize, nsize)
               if (myrank == 0) write (7) work_io
               call mpe_gather_io(work_io, plold, len2/nsize, nsize)
               if (myrank == 0) write (7) work_io
               !           call mpe_gather_io(work_io,dsqgeo,len2/nsize,nsize)
               !           call mpe_gather_io(work_io,spgeo,len2/nsize,nsize)
               deallocate (work_io)
               if (myrank == 0) then
                  call flush (7)
                  close (7)
               end if
               cmdxx = 'mkdir -p '//trim(phyout)//'phyout_'//ctau
               call system(trim(cmdxx))
               rfile = trim(phyout)//'phyout_'//ctau//'/'//ccore
               !         if (myrank .eq. 0) &
               !           open (unit=10,file=rfile,form='unformatted')
               i = 200 + myrank
               open (i, file=rfile, form='unformatted')
               !jwhwu 202004 avoid undefined values.
               write (i) land
               write (i) ocean
               write (i) ice
               write (i) alb
               write (i) z0
               write (i) tgclim
               write (i) gwclim
               write (i) sgeo
               write (i) canopy
               write (i) ustar
               write (i) tstar
               write (i) qstar
               write (i) raincu
               write (i) rainlp
               write (i) totalp
               write (i) curate
               write (i) plcl
               write (i) cumtop
               write (i) snr
               write (i) sncover
               write (i) sndepth
               write (i) tg
               write (i) gwr
               write (i) gwet
               write (i) cice
               write (i) xtice
               write (i) zice
               write (i) fpsp
               write (i) fpsp1
               !
               write (i) hflux
               write (i) qflux
               write (i) ss
               write (i) rs
               write (i) asol
               write (i) olr
               write (i) sld
               write (i) rld
               write (i) asold
               !jwhwu 202003
               write (i) sfemis
               write (i) sfalb
               write (i) std
               !jwhwu
               write (i) ctot
               write (i) clds
               write (i) pdiff
               write (i) tsave
               !jwhwu 201702 | add
               write (i) tsflw
               !jwhwu 201702 | add
               write (i) cosz
               write (i) slag, sdec, cdec, solcon, solhr, rsolhr, &
                  ! overcome round off problem in restart
                  tau, hours, itimestep, xy

               write (i) o3l
               write (i) ftp
               write (i) fqp
               write (i) ftp1
               write (i) fqp1
               !jwhwu 201702 ! add
               write (i) asl
               write (i) atl
               write (i) dtrad
               !jwhwu 201910 ! add
               write (i) deltaq
               write (i) cnvwr
               write (i) cnvcr
               !jwhwu 202208 ! add
               write (i) dtcup
               write (i) ducup
               write (i) dvcup
               write (i) dtshl
               write (i) dushl
               write (i) dvshl
               write (i) dtlsp
               write (i) hfiltx
               write (i) alphax
               write (i) pdryi
               !
               write (i) qt
               !         write(i) qp  !  not need
               write (i) smc
               write (i) stc
               write (i) slc
               close (i)
            end if ! end of ( mod(float(itau),float(itauezz)) .lt. 0.01 )
            if (do_sit) then
               if (myrank .eq. 0) print *, 'ready rerun_sitgrid1'
               call rerun_sitgrid1(itau)
               if (myrank .eq. 0) print *, 'ready rerun_sitgrid2'
               call rerun_sitgrid2(itau)
               if (myrank .eq. 0) print *, 'ready rerun_sitgrid3'
               call rerun_sitgrid3(itau)
            end if
         end if   ! end of (wrestrt)
         !
         !        !---output sigma layer data---
         !        dtaup = mod( tauo+0.001, taup )
         !ds        dtaup = mod( tau+0.001, taup )
         !ds       if(tau.le.(taureg+0.001) .and. dtaup.lt.0.01)then
         if (tau .le. (taureg + 0.001) .and. histim) then
            !jh        if( histim ) then
#ifndef NO_OUT
            call outsigs(itau, nx, my, my_max, lev, ncld &
                         , idtg, ptop, rad, grav &
                         , cp, cosl, pt, sgeo, snr, gwr, tg, pk, pk2 &
                         , ut, vt, tt, qt, phi, km_soil, smc &
                         , slc, stc, canopy, zice, ggdef, gmdef)
#endif
         end if
         !-------------------------------------------------------------------------------
         !     !---output sigma layer radiation data---
         !      if (myrank .eq. 0) print *,'ioutsigr =',ioutsigr
         if (ioutsigr .eq. 0 .and. histim) then
            if (myrank .eq. 0) print *, 'outsigr start !!!'
#ifndef NO_OUT
            call outsigr(itau, nx, my, my_max, lev &
                         , idtg &
                         , fusl, fdsl, fuir, fdir &
                         , fuslr, fdslr, fuirr, fdirr &
                         , asl, atl, asl_clr, atl_clr &
                         , dtrad, clds, vvel &
                         , ss, rs, olr, asol, sld, rld &
                         , ss_clr, rs_clr, olr_clr &
                         , asol_clr, sld_clr, rld_clr &
                         , cice, xtice, snr, sncover, snoalb &
                         , ctot, chig, cmid, clow &
                         , ggdef, gmdef)
#endif
            if (myrank .eq. 0) print *, 'outsigr ok !!!'
         end if
         !-------------------------------------------------------------------------------
         !  write operational fields (p-levels and surface)
         !
#ifndef NO_OUT
         !
         !     output cice
         !      call unify_reduceintp(nx,my,my_max,u10,glob)
         !      if ( myrank .eq. 0 ) then
         !         open(35,file='uv10.txt')
         !         write(35,'(i8,1x,i6.6)') idtg/10000,itau
         !         write(35,'(1552f6.2)') (glob(:,j),j=1,my)
         !      endif
         !      call unify_reduceintp(nx,my,my_max,v10,glob)
         !      if ( myrank .eq. 0 ) then
         !         write(35,'(1552f6.2)') (glob(:,j),j=1,my)
         !         close(35)
         !      endif
         !
         !      call unify_reduceintp(nx,my,my_max,cice,glob)
         !      if ( myrank .eq. 0 ) then
         !         open(35,file='cice.txt')
         !         write(35,'(i8,1x,i6.6)') idtg/10000,itau
         !         write(35,'(1552f6.2)') (glob(:,j),j=1,my)
         !         close(35)
         !      endif

         call transr(jtrun, jtmax, nx, my, my_max, levp, poly, vornow, cc, 1, nsizey)
         call ujoinsr(cc, rvor, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
         call outflds(itau, nx, my, my_max, lev, ncld &
                      , lmax, numout, idtg, outdir &
                      , ktrop, ptop, capa, cp, rgas, grav, sigma, sgeo &
                      , ptend, pt, plt, pk, pk2, phi, ut, vt, vvel &
                      , tt, qt, rdiv, rvor, tg, gwr, z0, hflux, qflux, snr &
                      , raintot, raincu, rainlp, asol, olr, ss, rs, alb, gwclim &
                      , acld, cosl, drag, ugws, vgws, t2, q2, rh2, rh10, u10, v10, gfx, rld, sld &
                      !byl                    , km_soil,smc,slc,stc,canopy,ggdef,slp,v850,v700,h850,h500 &
                      , km_soil, smc, slc, stc, canopy, ggdef, typtrk &
                      !xb110                    , ctot,chig,cmid,clow,hpbl,histim,flash,do_sit)
                      , ctot, chig, cmid, clow, hpbl, histim, do_sit)
#endif
         !
         !#ifdef RSM_sigp
#ifdef RSM
#ifdef RSM_sig
         ! for sigma coordinate
         if (outrsm .and. mod(float(itau) + 0.00001, float(rsmoutinv)) .lt. 0.01) then
            if (myrank .eq. 0) print *, ' call rsmout for rsm output at tau=', itau
            call rsmout(idtg, itau, nx, my, my_max, lev, ncld &
                        , ptop, cp, rgas, grav, sgeo, pdiff &
                        , t1000, pt, plt, pk, pk2, phi, ut, vt &
                        , tt, qt, tg, snr, cosl &
                        , km_soil, smc, stc &
                        , ice, land, ocean)
         end if
#else
         ! for sigma-P coordinate
         if (outrsm .and. mod(float(itau) + 0.00001, float(rsmoutinv)) .lt. 0.01) then
            if (myrank .eq. 0) print *, ' call rsmout for rsm output at tau=', itau
            call rsmout_sigp(itau, nx, my, my_max, lev, ncld &
                             , idtg, ptop, rad, grav, cosl &
                             , pt, sgeo, snr, gwr, tg, pk &
                             , ut, vt, tt, qt, km_soil, smc, stc &
                             , ice, land, ocean, xlon, xlat)
         end if
#endif
#endif
         !
         !#ifdef RSM
        !! RSM: output base field ncep-format data for RSM
         !      if(outrsm .and. mod(float(itau)+0.00001, float(rsmoutinv) ) .lt. 0.01)then
         !        if(myrank.eq.0)print*,' call rsmout for rsm output at tau=',itau
         !        call rsmout(idtg,itau,nx,my,my_max,lev,ncld   &
         !                , ptop,cp,rgas,grav,sgeo,pdiff        &
         !                , t1000,pt,plt,pk,pk2,phi,ut,vt       &
         !                , tt,qt,tg,snr,cosl                   &
         !                , km_soil,smc,stc                     &
         !                , ice,land,ocean)
         !      endif
         !#endif
         !
         !       if(typhoon .and. ltrack)then
         if (typhoon .and. ltrack .and. itau .le. 384) then
            if (myrank .eq. 0) print *, ' calling tracking,  tau= ', tau
            !          if(myrank .eq. 0)print *,'dt_trk=',dt_trk
            call tracking(tau, dt_trk, dt, nx, my, &
                          ntyph, typname, ixtyp, jytyp, tlon, tlat, tflon, tflat, idtg, &
                          nrec, typhoon, tensity)
         end if
         !
         !  zero out precip arrays
         !
         !  add 12-hour check to let prep. amount be of 12-hour accumulation for
         !  every 12-hour output, but prep. amount still 24-hour accumulation for
         !  every 24-hour output without 12-hour output points
         !  (in order to be consistent to the rule followed by nfs, 11/30/2000)
         !
         if (mod(float(itau) + 0.00001, 12.) .lt. 0.01) then
            !        if( histim .and. (mod(float(itau)+0.00001, 12.) .lt. 0.01) ) then
            raincu = 0.
            rainlp = 0.
            !         runoff=0.
            !          call zilch (raincu,nxmy)
            !          call zilch (rainlp,nxmy)
            !          call zilch (runoff,nxmy)
         end if
         !
         !
         if (myrank .eq. 0) print *, ' history written at tau=', itau
      end if     ! end of (histim) --- --- ---

      !       output f006 data for FV3
      if (outfv3) then
         if (itau == 6) then
#ifndef NO_OUT
            if (myrank .eq. 0) print *, 'output FV3 data !!!'
            call outflds_fv3(nint(tau), nx, my, my_max, idtg, ggdef &
                             , q2, fm, fh, fm10, fh2, srflag, ustar)
#endif
         end if !(abs(tau+0.00001-6.) .lt. 0.01)
      end if

      !
      !kc             output t2,raintot,u10,v10,ctot at 1 hour interval within 192hr.
      !               if(domfc)then
      !==xb118        change domfc type from logical to real
      !               if(myrank .eq. 0) print*,'domfc at tau,dtaup=',tau,dtaup
      !               if(tau.le.192. .and. dtaup.lt.0.01)then

      if (mod(itau, 1) == 0) then

         if (itau .le. nint(domfc)) then
            if (myrank .eq. 0) print *, 'out1 at tau=', itau
#ifndef NO_OUT
            call out2d_mfc(nx, lev, my, my_max, itau, idtg &
                           , raincu1, rainlp1, raintot, glob, t2, q2, rh2, rh10 &
                           , u10, v10, tmax, tmin, rld, sld, ctot, pt, ggdef)
#endif
         end if !  ( tau <= domfc )
         raincu1 = 0.0
         rainlp1 = 0.0
         tmax = 0.0
         tmin = 0.0
      end if
      !

      !
      !       out green energy plan
      if (out_green .and. mod(itau, nint(otgreen)) == 0) then
#ifndef NO_OUT
         call outflds_green(nint(tau), nx, my, my_max, lev, ncld &
                            , idtg, cp, rgas, grav, t2, u10, v10, ss, pk &
                            , sgeo, pt, plt, ptop, ut, vt, tt, qt, cosl, raincu6, rainlp6)
#endif
         if (mod(itau, 6) == 0) then
            raincu6 = 0.
            rainlp6 = 0.
         end if
      end if
      !
      !       output flux at 24 hour interval
      !
      if (mod(itau, 24) == 0) then
         if (myrank .eq. 0) print *, 'out24 at tau=', itau
         if (myrank .eq. 0) print *, 'julian = ', julian
         !         call mpe_unify(hf24,nx,my,2,mpe_double)
         !         call mpe_unify(qf24,nx,my,2,mpe_double)
         !         call mpe_unify(ss24,nx,my,2,mpe_double)
         !         call mpe_unify(rs24,nx,my,2,mpe_double)
         !         call mpe_unify(asol24,nx,my,2,mpe_double)
         !         call mpe_unify(olr24,nx,my,2,mpe_double)
         !         call mpe_unify(rain24,nx,my,2,mpe_double)
#ifndef NO_OUT
         call out24(nx, my, my_max, hf24, qf24, ss24, rs24, asol24, olr24, rain24, rainlp24, dt24 &
                    , glob, itau, idtg, ggdef, flash24)
#endif
         !
         !         zero set arrays
         !
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               hf24(i, jj) = 0.0
               qf24(i, jj) = 0.0
               ss24(i, jj) = 0.0
               rs24(i, jj) = 0.0
               asol24(i, jj) = 0.0
               olr24(i, jj) = 0.0
               rain24(i, jj) = 0.0
               rainlp24(i, jj) = 0.0
               flash24(i, jj) = 0.0
            end do
         end do
         dt24 = 0.0
      end if ! (  mod( itau , 24 ) == 0  )

      !
      if (mod(itau, 1) == 0) then
         if (dosppt .and. dospptout) then
            call spptout(tau)
         end if
         if (doskeb .and. doskebout) then
            call skebout(tau)
         end if
      end if

      if (do_sit .AND. lgodas .AND. ldailysst) then
         CALL read_dailygodas(idtg, tau, dtx)
      end if

      if (myrank == 0) then
         call system_clock(toutend, toutrate)
         print *, "In output tau: ", itau, " CPU Time: " &
            , dble(toutend - toutsrt)/dble(toutrate)
      end if

      if (histim .or. ltrack) then
         ifromtau = taui
         itotau = taue
         flag = .false.
         if (myrank .eq. 0) then
            flag = .true.
            !CWB2016
            if (.not. io_quilting) then
               !CWB2017             call sendmsg ('gfs',ifromtau,itotau,istat)
               if (itau .eq. itotau) call sendmsg('gfs', ifromtau, itotau, istat)
            else
               istat = 0
            end if
            if (istat .eq. -1) then
               print *, ' SENDMSG ERROR ', 'RANK=', myrank
               call dmsexit(-1)
            end if
         end if
         !ch       call mpe_broadcast(istat,1,flag,mpe_integer)
         !         call mpe_bcast(istat,1,0,mpe_integer)
      end if !histim

   end if !  (mod(tau+0.001, 1.) .lt. 0.01)  hourly for output
   !    ---------------------------------------------------------------
   !
   ! new year, read obs sst
   !
   if (lopgsst) then
      inexttau = int(tau + dtx/3600.+0.001)
      call dtgfix12(idtg, idtg1_sst, inexttau)
      inextyear = idtg1_sst/100000000
      lnewyear = icurrentyear /= inextyear
      if (lnewyear) then
         call read_opgsst(idtg1_sst, ggdef, ocean, ice)
         icurrentyear = idtg1_sst/100000000
      end if
   end if   !end lopgsst
   !
   itau = tau + 0.001
#ifdef TIMING
   tm_2 = mpi_wtime()
#endif
   if (itau .lt. itaue) go to 100
   !
   flag = .false.
   if (myrank .eq. 0) then
      call recmsg('gfs', ifromtau, itotau, istat)
      flag = .true.
   end if
   !ch   call mpe_broadcast(istat,1,flag,mpe_integer)
   call mpe_bcast(istat, 1, 0, mpe_integer)

   if (istat .eq. -1) then
      print *, ' RECMSG ERROR ', 'RANK=', myrank
      call dmsexit(-1)
   else if (istat .eq. 1) then
      if (myrank .eq. 0) then
         print *, ' finished integration '
         ! ------------------------------------------------------------
         ! << openacc deallocate data >>
         !$acc exit data delete(um,vm)
         ! ------------------------------------------------------------

         ! for io quilting
         if (io_quilting) then
            !         ntag=ntag+1
            !         call mpe_send_key(keydoit,ntag,istat)
            ntag = ntag + 1
            call mpe_send_key(keydone, ntag, istat)
         end if

         call dmscls(bckfile, istat)
         call dmscls(ifilin, istst)
         call dmscls(ifilout, istst)
         if (ldailyFCTsst) call dmscls(ifilin_sst, istat)
         if (ldailyFCTicesndpt) call dmscls(ifilin_ncep, istat)
         if (dailyClm_option .ge. 1) then
            call dmscls(ifilin_ClmANA, istat)
            if (dailyClm_option .eq. 2) then
               call dmscls(ifilin_ClmFCT, istat)
            end if
         end if
         close (77)
      end if
      if (ldailyFCTsst .OR. ldailyFCTicesndpt .OR. (dailyClm_option .ge. 1)) then
         call deallocate_dailyFCT_array
      end if
      if (lopgsst) then
         call deallocate_opgsst_array
      end if
      if (do_sit) then
         if (myrank == 0) print *, 'ready sit_vdiff_end'
         if (locaf) call deallocate_ocaf_array
         if (lwoa0) call deallocate_woa0_array
         if (lgodas) call deallocate_godas_array
         call deallocate_sitgrid_array
         call cleanup_netcdf
         call sit_vdiff_end
      end if
      return
   end if
   !
   !ch   call mpe_broadcast(ifromtau,1,flag,mpe_integer)
   !ch   call mpe_broadcast(itotau,1,flag,mpe_integer)
   call mpe_bcast(ifromtau, 1, 0, mpe_integer)
   call mpe_bcast(itotau, 1, 0, mpe_integer)
   !
   taui = float(ifromtau)
   taue = float(itotau)
   tauo = float(itotau)
   itaui = taui + 0.1
   itaue = taue + 0.1
   itauo = tauo + 0.1
   !
   go to 10
   !
   ! finilize stochastic_physics
   call destroy_stochastic_physics()
   close (35)

end subroutine intgrt_gpu
