subroutine initial_gpu(no, jtrun, jtmax, lev, nx, my, my_max, mlmax)

   !
   !  purpose : do nonlinear normal mode initialization
   !  modify to f90 by C-H Lee and sort by River Chen in 2015
   !-----------------------------------------------------------
   !-----------------------------------------------------------------------
   !  **** input *****
   !  - no
   !  - rdiv
   !  - ut
   !  - vt
   !  - tt
   !  - qt
   !  - pt
   !  - dtpl
   !  - dlpl
   !  - phi
   !  - vornow
   !  - divnow
   !  - temnow
   !  - plnow
   !  - rad, omega, eigval
   !  - tmcor, pmcor
   !  - << doincr = ture >
   !    - temold, vorold, divold
   !
   !  **** output *****
   !  - rdiv
   !  - /rvor/
   !  - ut
   !  - vt
   !  - tt
   !  - qt
   !  - pt
   !  - dtpl
   !  - dlpl
   !  - vornow
   !  - divnow
   !  - temnow
   !  - plnow
   !  - vorten
   !  - divten
   !  - vorold, divold, temold, plold
   !  - hldten, dlphi, dtphi, plten
   !  - pk, pk2, plt
   !  - sd, vvel
   !
   !  **** buffer ****
   ! - nnmi_buf, zx_buf, temten1
   ! - spec_cg_buffer
   !-----------------------------------------------------------------------
   !  **** present on device ****
   !  << param >>
   !  - nx, my, lev, jtrun, ncld,
   !  - my_max, jtmax
   ! << index >>
   ! - mlistnum, mlist, nlist,
   ! - jlistnum, jlist1, jlist2, jlist2_2d, jlist1_sl,
   ! - nxdef, mtrundef,
   ! - nsizex, nsizey, tcolt_jlist, poly_mlist,
   ! -  nxp, levp, nxptot,
   ! - Llist, nxjp, nxjstart, nxjend, nxjlen, nxdef_2d, nxjp_acc, nxjlen_all
   ! << const >>
   ! - poly, dpoly, polyf, dpolyf, weight, cim,
   ! - wcfac, wdfac, onocos
   ! << spec >>
   ! - temten, vorten, divten, hldten, plten,
   ! - temold, vorold, divold, plold,
   ! - temnow, vornow, divnow, plnow
   ! << grid >>
   ! - rdiv, ut, vt, tt, qt, pt,
   ! - dtpl, dlpl,
   ! - pk, pk2, plt, sgeo, phi, dtphi, dlphi,
   ! - sd, vvel
   use param, only: ncld
   use mpe
   use rank
   use index
   use const, only: eigval, omega, rad, nnmiit, doincr, evecin, nnmivm, &
                    cutfreq, evectr, pmcor, tmcor, poly, dpoly, cim, &
                    wdfac, wcfac, onocos, coslr, polyf, dpolyf, RTYPE, &
                    ! ----------
                    weight, radsq, cp, sinl, dsigma, sigma, cor, ptop, cosl, &
                    spalm, arrhyd
   use spec, only: temold, vorten, vorold, divten, divold, plnow, temnow, &
                   divnow, vornow, plold, plten, &
                   !byl                        divnow,vornow,qold,qnow,plold
                   ! ----------------------------------------
                   hldten, temten

   use grid, only: pt, tt, rdiv, rvor, dtpl, dlpl, ut, vt, &
                   ! ----------
                   latstr, latlen, lonstr, lonlen, &
                   gglati, fa1, fa2, fa3, fa4, &
                   ! ----------------------------------------
                   qt, plt, pk, pk2, dtpl, dlpl, &
                   sd, vvel, up, vp, ttp, qm, &
                   ut_sl, vt_sl, &
                   phi, dlphi, dtphi

   use fftcom
   use spec_cuda_graph, only: cc_cg, gwk1_cg, wcc_fk_cg, &
                              wc_cg, ws_cg, fj_weight_cg, &
                              transr_fft_cg, transr_lt_cg, &
                              transr1_fft_cg, &
                              trngra_fft_cg, &
                              tranuv_fft_cg, tranuv_lt_cg, &
                              tranrs_fft_cg, tranrs_lt_cg, &
                              tranrs1_fft_cg, &
                              trngra3_fft_cg, trngra3_lt_cg, &
                              trandv_fft_cg, trandv_lt_cg
   use openacc

   implicit none
   integer no, jtrun, jtmax, lev, nx, my, my_max, mlmax
   !
   !  working array
   !
   !     integer, parameter ::  no=2*((jtrun+1)/2)+(jtrun/2)+10
   !byl      real      a(jtrun,jtrun,lev),b(jtrun,jtrun,lev)
   real, dimension(jtrun, jtrun, nnmivm):: a, b, c
   real, dimension(jtrun, jtmax, nnmivm):: h
   real, dimension(no*no, 2*jtmax*nnmivm):: mx, evec, epos
   real, dimension(no, 2*jtmax*nnmivm):: eval
   real, dimension(no, 2, 2*jtmax*nnmivm):: x
   integer nnlist(2*jtmax*nnmivm)

   integer nw(jtrun, jtmax)
   real(kind=RTYPE) phiten(levp, 2, jtrun, jtmax), dummy
   real(kind=RTYPE) temten1(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) zx_buf(lev*6*jtrun*jtmax)
   real nnmi_buf(no, 2, 2*jtmax*nnmivm)
   logical zx_forward_cg_created, zx_backward_cg_created, nnmi_cg_created, &
      tendget_matmul_cg_created
   type(acc_graph_t) zx_forward_cg, zx_backward_cg, nnmi_cg, tendget_matmul_cg

   !byl      real      cc(nx+2,levp,3,my_max),wss(levp,2,3,jtrun,jtmax)
   real bal_tmp(jtrun, nnmivm)
   character lab*10, lrec*16

   integer mlmax2, j, k, l, ic, m, mf, n, ns, na, nbig, kk, kL, jj
   real bal
   logical nnmical
   integer brank   !root rank of row_broadcast
   !
   integer async_id
   async_id = 1

  !! << const >>
   !$acc enter data copyin(cutfreq, eigval, omega, rad, evecin, evectr, pmcor, tmcor, &
   !$acc& coslr, radsq, cosl, spalm, arrhyd &
   !$acc& ) async(async_id)
   !$acc enter data copyin(cp, sinl, dsigma, sigma, cor, ptop) async(async_id)
  !! << index >>
   !$acc enter data copyin(nnmivm, Lstart, Lend) async(async_id)
  !! << grid >>
   !$acc enter data copyin(lonlen, lonstr, latlen, fa1, fa2, fa3, fa4, gglati) async(async_id)

  !! << local >>
   !$acc enter data create(phiten, temten1, &
   !$acc&                  a, b, c, h, x, mx, nnlist, &
   !$acc&                  eval, bal_tmp, &
   !$acc&                  zx_buf, nnmi_buf) async(async_id)

  !! << tendget >>
   !$acc enter data create(rvor, up, vp, ttp, qm, &
   !$acc& ut_sl, vt_sl) async(async_id)

   !
 
   tendget_matmul_cg_created = .false.
   zx_forward_cg_created = .false.
   zx_backward_cg_created = .false.
   nnmi_cg_created = .false.

   mlmax2 = mlmax*2
   if (myrank .eq. 0) print *, 'jrtun=', jtrun, ' mlmax=', mlmax
   lab = 'pt'
   call check(pt, nx, my, my_max, lab)
   !
   !  define constants and comput coefficients for initializatin
   !
   do L = 1, nnmivm
      do m = 1, mlistnum
         mf = mlist(m)
         nbig = jtrun - mf + 1
         !  calculate the size of symmetric and antisymmetric matrix which
         !  include the gravity and rossby wave
         ns = 2*int((nbig + 1)/2) + int(nbig/2)
         na = 2*int(nbig/2) + int((nbig + 1)/2)
         j = 1 + (m - 1)*2 + (L - 1)*2*jtmax
         nnlist(j) = ns
         nnlist(j + 1) = na
      end do
   end do
   call inicons(rad, omega, eigval, lev, jtrun, jtmax, &
                nw, a, b, c, h, nnmivm)
   !$acc update device(h, nnlist) async(async_id)
   !
   !  construct coefficient matrix
   !
   do L = 1, nnmivm
      if ((L .ge. Lstart) .and. (L .le. Lend)) then
         k = L - Lstart + 1
         do m = 1, mlistnum
            mf = mlist(m)

            nbig = jtrun - mf + 1
            ns = 2*int((nbig + 1)/2) + int(nbig/2)
            na = 2*int(nbig/2) + int((nbig + 1)/2)

            j = 1 + (m - 1)*2 + (L - 1)*jtmax*2
            call coftrix_refactor(mf, L, mx(1, j), ns, &
                                  a(1, 1, L), b(1, 1, L), c(1, 1, L), jtrun, lev, +1)
            call coftrix_refactor(mf, L, mx(1, j + 1), na, &
                                  a(1, 1, L), b(1, 1, L), c(1, 1, L), jtrun, lev, -1)
         end do
      end if
   end do
   !$acc update device(mx) async(async_id)
   ! ************************************************************
   do L = 1, nnmivm
      brank = 0
      if ((L .ge. Lstart) .and. (L .le. Lend)) then
         brank = row_rank
      end if
      j = 1 + (L - 1)*jtmax*2
      call mpe2d_row_broadcast_gpu(mx(1, j), no*no*2*jtmax, brank)
   end do
   ! ************************************************************
   !
   !   fine the eigenvector and eigenvalues of the matrix
   !
   call eigen_mx_gpu(mx, eval, nnlist, no)
   ! ************************************************************
   !
   !  begin to iterration, now doing 3 iterrations
   !
  !! ****************************************
   do ic = 1, nnmiit
      !
      if (myrank .eq. 0) print *, 'iteration=', ic
      !$acc parallel loop collapse(2) private(mf) async(async_id)
      do L = 1, nnmivm
         do m = 1, jtrun
            bal_tmp(m, L) = 0.
         end do
      end do
      !
      !  get tendency of vorticity,divergence,geopotential
      !
     !! ****************************************
      call tendget_gpu(phiten, zx_buf(1), zx_buf(1 + 2*lev*jtrun*jtmax), &
                       tendget_matmul_cg_created, &
                       tendget_matmul_cg)
     !! < input >
     !!   pt, ut, vt, tt, qt
     !!   rdiv, phi, dlpl, dtpl
     !! < output >
     !!   phiten, vorten, divten
     !!   hldten, dlphi, dtphi
     !!   plten, pk, pk2, plt
     !!   sd, vvel
     !! ****************************************
      !
      !--------------------------------------------------------
      !  incremental initialization
      !
      if (doincr) then
         print *, "doincr"
         !$acc parallel loop collapse(4) private(mf) async(async_id)
         do m = 1, mlistnum
            do n = 1, jtrun
               do j = 1, 2
                  do k = 1, levp
                     mf = mlist(m)
                     if (n .ge. mf) then
                        phiten(k, j, n, m) = phiten(k, j, n, m) - temold(k, j, n, m)
                        vorten(k, j, n, m) = vorten(k, j, n, m) - vorold(k, j, n, m)
                        divten(k, j, n, m) = divten(k, j, n, m) - divold(k, j, n, m)
                     end if
                  end do
               end do
            end do
         end do
      end if
      !
      !--------------------------------------------------------
      !
      !  do vertical transform
      !
      call zx_gpu(evecin, vorten, divten, phiten, jtrun, jtmax, lev, &
                  zx_buf, zx_forward_cg_created, zx_forward_cg)
      !
      !  create variable vector
      !
      call vartrix_gpu(vorten, divten, phiten, levp, x, no, &
                       rad, omega, h, +2)
      ! ************************************************************
      do L = 1, nnmivm
         brank = 0
         if ((L .ge. Lstart) .and. (L .le. Lend)) then
            brank = row_rank
         end if
         j = 1 + (L - 1)*jtmax*2
         call mpe2d_row_broadcast_gpu(x(1, 1, j), no*4*jtmax, brank)
      end do
      ! ************************************************************
      !
      !   perform nonlinear normal mode initialization
      !
      call nnmi_gpu(x, eval, mx, no, bal_tmp, cutfreq, nnlist, &
                    nnmi_cg_created, nnmi_cg, nnmi_buf)
      ! ****************************************
      !$acc update self(bal_tmp) async(async_id)
      !$acc wait(async_id)
      do L = 1, nnmivm
         if (myrank .eq. 0) print *, 'vertical mode l=', L
         !
         ! set convergence variable
         !
         bal = 0.
         call mpe_unify(bal_tmp(1, L), 1, jtrun, 3, mpe_double)
         do mf = 1, jtrun
            bal = bal + bal_tmp(mf, L)
         end do
         !
         if (myrank .eq. 0) print *, 'bal=', bal
      end do
      ! ****************************************
      !
      !   decompose the variable vector back to 3 individual variable
      !
      call vartrix_gpu(vorten, divten, phiten, levp, x, no, &
                       rad, omega, h, -2)
      !
      !
      !   initial none initialized mode
      !

      !2dMPI>
     !!         kL=1
     !!         kk=Llist(nnmivm+1)
     !!         if(kk .le.levp)then
     !!            kL=nnmivm+1
     !!         endif
      !2dMPI<

      !$acc parallel loop collapse(4) private(mf, kk) async(async_id)
      do m = 1, mlistnum
         do n = 1, jtrun
            do j = 1, 2
!              do k = nnmivm+1,lev
               do k = 1, levp
                  mf = mlist(m)
                  kk = Llist(k)
                  if ((n .ge. mf) .and. (kk .gt. nnmivm)) then
                     vorten(k, j, n, m) = 0.
                     divten(k, j, n, m) = 0.
                     phiten(k, j, n, m) = 0.
                  end if
               end do
            end do
         end do
      end do
     !! ****************************************
      !
      !  conversion  structure of variables (phiten,vorten,divten)
      !
      !   vertical transform back
      !
      call zx_gpu(evectr, vorten, divten, phiten, jtrun, jtmax, lev, &
                  zx_buf, zx_backward_cg_created, zx_backward_cg)
     !! < input/output >
     !!   vorten, divten, phiten
     !! ****************************************

     !! ****************************************
      !
      !   add the correction to variables
      !
      call correct_gpu(vorten, divten, phiten, tmcor, pmcor, vornow, &
                       divnow, temnow, plnow, jtrun, jtmax, lev)
     !! < input >
     !!   vorten, divten, phiten, tmcor, pmcor
     !!   vornow, divnow, temnow, plnow
     !! < output >
     !!   vornow, divnow, temnow, plnow
     !! ****************************************

      !$acc parallel loop collapse(2) private(mf) async(async_id)
      do m = 1, mlistnum
         do k = 1, levp
            mf = mlist(m)
            if (mf .eq. 1) then
               vornow(k, 1, 1, m) = 0.
               vornow(k, 2, 1, m) = 0.
               divnow(k, 1, 1, m) = 0.
               divnow(k, 2, 1, m) = 0.
            end if
         end do
      end do
     !! ****************************************
      !
      !   transform back to phyical space
      !
     !!        call joinsr(wss,vornow,divnow,temnow,dummy,jtrun,jtmax,levp &
     !!                   ,mlistnum,3,1)
     !!        call transr(jtrun,jtmax,nx,my,my_max,levp,poly,wss,cc,3,nsizey)
     !!        call ujoinsr(cc,rvor,rdiv,tt,dummy,nx,my_max,lev,jlistnum,3,1)
      call transr_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, polyf, &
                                 vornow, cc_cg, 1, nsizey, &
                                 gwk1_cg, ws_cg, wc_cg, wcc_fk_cg)
      call ujoinsr_gpu(cc_cg, rvor, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
      call transr_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, polyf, &
                                 divnow, cc_cg, 1, nsizey, &
                                 gwk1_cg, ws_cg, wc_cg, wcc_fk_cg)
      call ujoinsr_gpu(cc_cg, rdiv, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
      call transr_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, polyf, &
                                 temnow, cc_cg, 1, nsizey, &
                                 gwk1_cg, ws_cg, wc_cg, wcc_fk_cg)
      call ujoinsr_gpu(cc_cg, tt, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
      call transr1_gpu(jtrun, jtmax, nx, my, my_max, polyf, plnow, pt, nsizey, &
                       cc_cg, gwk1_cg)
      !
      !  compute zonal and meridional gradients of terrain pressure
      !
      call trngra_gpu(jtrun, jtmax, nx, my, my_max, cim, poly, dpoly, plnow, &
                      dlpl, dtpl, nsizey, cc_cg, gwk1_cg)
      !
      !   transform spectrum vorticity , divergence to physical u , v
      !
      call tranuv_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, &
                                 coslr, wcfac, wdfac, polyf, dpolyf, &
                                 vornow, divnow, ut, vt, nsizey, &
                                 cc_cg, gwk1_cg, ws_cg(1, 1), ws_cg(1, 2), &
                                 wc_cg, wcc_fk_cg, &
                                 fj_weight_cg(1, 1), fj_weight_cg(1, 2))
     !! < input >
     !!   vornow, divnow, temnow, plnow
     !! < output >
     !!   rvor, rdiv, tt, pt, dlpl, dtpl, ut, vt
     !! ****************************************
      !$acc update self(pt) async(async_id)
      !$acc wait(async_id)
      !
      lab = 'pt '
      call check(pt, nx, my, my_max, lab)
   end do
  !! ****************************************
  !! << tendget >>
   !$acc exit data delete(rvor, &
   !$acc& up, vp, ttp, qm, &
   !$acc& ut_sl, vt_sl &
   !$acc& ) async(async_id)
  !! ****************************************

   !$acc parallel loop collapse(3) private(mf) async(async_id)
   do m = 1, mlistnum
      do n = 1, jtrun
         do k = 1, levp
            mf = mlist(m)
            if (n .ge. mf) then
               vorold(k, 1, n, m) = vornow(k, 1, n, m)
               vorold(k, 2, n, m) = vornow(k, 2, n, m)
               divold(k, 1, n, m) = divnow(k, 1, n, m)
               divold(k, 2, n, m) = divnow(k, 2, n, m)
               temold(k, 1, n, m) = temnow(k, 1, n, m)
               temold(k, 2, n, m) = temnow(k, 2, n, m)
               if (k .eq. 1) then
                  plold(n, m, 1) = plnow(n, m, 1)
                  plold(n, m, 2) = plnow(n, m, 2)
               end if
            end if
         end do
      end do
   end do
   !
  !! << const >>
   !$acc exit data delete(eigval, omega, rad, evecin, evectr, pmcor, tmcor, &
   !$acc& coslr, radsq, cosl, spalm, arrhyd &
   !$acc& ) async(async_id)
   !$acc exit data delete(cp, sinl, dsigma, sigma, cor) async(async_id)
  !! << index >>
   !$acc exit data delete(nnmivm, Lstart, Lend) async(async_id)
  !! << grid >>
   !$acc exit data delete(lonlen, lonstr, latlen, fa1, fa2, fa3, fa4, gglati) async(async_id)
  !! << local >>
   !$acc exit data delete(phiten, temten1, &
   !$acc&                  a, b, c, h, x, mx, nnlist, &
   !$acc&                  eval, bal_tmp, &
   !$acc&                  zx_buf, nnmi_buf) async(async_id)
   return
end subroutine initial_gpu
