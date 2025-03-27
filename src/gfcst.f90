#ifndef TIMCOMCPL
      program gfcst
#else
      subroutine gfcst(compid, mpi_comm_mct)
      use gfs_cpl, only:gfs_cpl_init
#endif
!
! main program of CWBGFS
! modify to f90 bt C-H Lee and sort by River Chen in 2015
!
#ifdef USE_CUDA
   use param, only: nx, my, lev, jtrun, ncld, mlmax, &
                    my_max, jtmax
   use index, only: mlistnum, mlist, nlist, &
                    jlistnum, jlist1, jlist2, jlist2_2d, jlist1_sl, &
                    nxdef, mtrundef, &
                    nsizex, nsizey, tcolt_jlist, poly_mlist, &
                    nxp, levp, nxptot, &
                    Llist, nxjp, nxjstart, nxjend, nxjlen, &
                    nxdef_2d, nxjp_acc, nxjlen_all
   use const, only: donnmi, taui, ttl, &
                    cp, sigma, dsigma, ptop, ptmeans, tmeans, spalm, &
                    eigval, evecin, evectr, arrhyd, arsddt, pmcor, tmcor, &
                    ! --------------------
                    poly, dpoly, polyf, dpolyf, weight, cim, &
                    wcfac, wdfac, onocos  ! tranuv, trandv
   use grid, only: rdiv, ut, vt, tt, qt, pt, &
                   dtpl, dlpl, &
                   pk, pk2, plt, &
                   sgeo, phi, dlphi, dtphi, &
                   sd, vvel
   use spec, only: vornow, divnow, temnow, plnow, &
                   vorold, divold, temold, plold, &
                   vorten, divten, temten, &
                   hldten, plten
   use cudafor
#else
   use param
   use const
#endif
!
   implicit none

      integer  no
      integer async_id
#ifdef TIMCOMCPL
      integer, intent(in) :: compid, mpi_comm_mct
#endif
      async_id = 1

!
!  logical io units:
!
!  input namelist file = 'namlsts'
!  date-time-group file = 'crdate'
!  output directives file = 'ocards'
!  model history file = 'cwbout'
!  output model history for diabatic variables='phyout'
!  input file of path/file names='filist'
!
#ifdef TIMCOMCPL
      call mpe_init(mpi_comm_mct)
#else
      call mpe_init
#endif
!
!     get model constants
!
      call cons
#ifdef TIMCOMCPL
      call gfs_cpl_init(compid)
#endif
!
!  read in initial data and prepare for initialization/forecast
!
   call getrdy
!
!  initialization would be done here, taking initial spectral fields
!  out of 'getrdy' and preparing them for 'intgrt'.  for identification
!  purposes only, initialization history file data is assigned
!  "tau"=1.0.
!
   !! << param >>
   !$acc enter data copyin(nx, my, lev, jtrun, ncld, &
   !$acc&                  my_max, jtmax) async(async_id)
   !! << index >>
   !$acc enter data copyin(mlistnum, mlist, nlist, &
   !$acc& jlistnum, jlist1, jlist2, jlist2_2d, jlist1_sl, &
   !$acc& nxdef, mtrundef, &
   !$acc& nsizex, nsizey, tcolt_jlist, poly_mlist, &
   !$acc& nxp, levp, nxptot, &
   !$acc& Llist, nxjp, nxjstart, nxjend, nxjlen, nxdef_2d, nxjp_acc, nxjlen_all &
   !$acc& ) async(async_id)
   !! << const >>
   !$acc enter data copyin(poly, dpoly, polyf, dpolyf, weight, cim, &
   !$acc&                  wcfac, wdfac, onocos) async(async_id)

   !! << spec >>
   !$acc enter data create(temten, vorten, divten, hldten, plten, &
   !$acc&                  temold, vorold, divold, plold, &
   !$acc&                  temnow, vornow, divnow, plnow &
   !$acc& ) async(async_id)
   !! << grid >>
   !$acc enter data create(rdiv, ut, vt, tt, qt, pt, &
   !$acc&                  dtpl, dlpl, &
   !$acc&                  pk, pk2, plt, sgeo, phi, dtphi, dlphi, &
   !$acc&                  sd, vvel &
   !$acc& )async(async_id)

   !$acc update device(pt, ut, vt, tt, qt) async(async_id)
   !$acc update device(rdiv, phi, dtpl, dlpl) async(async_id)
   !$acc update device(vornow, divnow, temnow, plnow) async(async_id)
   !$acc update device(sgeo) async(async_id)
   if (donnmi .and. taui .lt. 1.0) then
      no = 2*((jtrun + 1)/2) + (jtrun/2) + 10
#ifdef USE_CUDA
      call initial_gpu(no, jtrun, jtmax, lev, nx, my, my_max, mlmax)
!      call initial(no, jtrun, jtmax, lev, nx, my, my_max, mlmax)
!!20250307 shian trans to GPU after CPU
!!$acc update device(pt, ut, vt, tt, qt) async(async_id)
!!$acc update device(rdiv, phi, dtpl, dlpl) async(async_id)
!!$acc update device(vornow, divnow, temnow, plnow) async(async_id)
!!$acc update device(sgeo) async(async_id)
      !$acc wait(async_id)
#else
      call initial(no, jtrun, jtmax, lev, nx, my, my_max, mlmax)
#endif
   end if
!
!  recompute eigen value for semi-implicit scheme
!
   call matrix_hybrid_cwb(cp, sigma, dsigma, ptop, ptmeans, tmeans, spalm &
                          , eigval, evecin, evectr, arrhyd, arsddt, pmcor, tmcor)
!
!  time integration
!
#ifdef TIMCOMCPL
   #ifdef USE_CUDA
      call cudaProfilerStart
      call intgrt_gpu(compid)
      call cudaProfilerStop
      call mpe_finalize
   #else
      call intgrt(compid)
      call mpe_finalize
   #endif

   end subroutine gfcst

   
#else
       if ( ttl ) then
  #ifdef USE_CUDA
        call cudaProfilerStart
        call intgrt_gpu
  !      call intgrt
        call cudaProfilerStop
  #else
        call intgrt
  #endif
   else
      call intgrt_3tl
   end if

!
   call mpe_finalize
   call dmsexit(0)
!
      stop
   end program gfcst
#endif
