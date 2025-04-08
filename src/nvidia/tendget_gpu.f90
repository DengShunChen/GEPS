#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)
subroutine tendget_gpu(phiten, temten1, phiten1, cg_created, cg_graph)

   !CWB2021 ndsl single precision test

   !
   !  purpose : compute tendency of vorticity,divergence and geopotential
   !            using subroutine of forecast model
   !-----------------------------------------------------------------------
   !  **** input *****
   !  the input are pass bye those '*.h' files (see below)
   !  - pt
   !  - ut
   !  - vt
   !  - tt
   !  - qt
   !  - rdiv
   !  - phi
   !  - dtpl
   !  - dlpl
   !
   !  **** output *****
   !  phiten   : geopotential tendency array
   !  vorten
   !  divten
   !  other output are pass by 'spec.h'
   !
   !  **** buffer ****
   !  - temten1
   !
   !  **** present on device *****
   !  << output >>
   !  phiten(cg), vorten, divten
   !  << buffer >>
   !  temten1(cg)
   !
   !  << sepc >>
   !  plten(cg), hlten, temten
   !  << grid >>
   !  rdiv, pt, ut, vt, tt, qt,
   !  plt, pk, pk2, dtpl, dlpl,
   !  sd, vvel, up, vp, ttp, qm,
   !  ut_sl, vt_sl,
   !  phi, dlphi, dtphi
   !
   !  << index >>
   !
   !  << const >>
   !  radsq, cp, ptop
   !  poly, dpoly, cim, wdfac, wcfac, onocos, coslr, polyf, dpolyf
   !  weight, sinl, cosl, sigma, dsigma, cor
   !  spalm, arrhyd
   !  << grid >>
   !  sgeo, latstr, latlen, lonstr, lonlen
   !  gglati, fa1, fa2, fa3, fa4
   !  << spec_cuda_graph >>
   !  cc_cg, gwk1_cg, wcc_fk_cg, wc_cg, ws_cg, fj_weight_cg
   !-----------------------------------------------------------------------
   use param
   use mpe
   use rank
   use index
   use const
   use spec
   use grid
   use mod_ndslfv_monoadv_gpu, only: ndslfv_update_gpu, &
                                     ndslfv_monoadvv_gpu

   use spec_cuda_graph, only: cc_cg, gwk1_cg, wcc_fk_cg, &
                              wc_cg, ws_cg, fj_weight_cg
   use openacc
   use cudafor
   use cublas

   implicit none

   integer m, mf, j, jj, k, n, l, mlst, nxj, nk, kk, kl, i, ii, nl, llistnum_fj
   real sqhaf

   real(kind=RTYPE), intent(out) :: phiten(levp, 2, jtrun, jtmax)
   real(kind=RTYPE) temten1(lev, 2, jtrun, jtmax)
   logical, intent(inout) :: cg_created
   type(acc_graph_t), intent(inout):: cg_graph
   !byl      real tbar(lev),qbar(lev*ncld)
   real(kind=RTYPE) sdpbl(nxp, my_max)

   ! for Semi-Lagrangian
   real(kind=RTYPE) deldm(nxp, my_max)

   !CWB2021 ndsl single precision test
   real(kind=RTYPE) uum_sl(nx, levp, my_max), &
      vvm_sl(nx, levp, my_max), &
      ttm_sl(nx, levp, my_max), &
      pten_sl(nx, levp, my_max), &
      qm_sl(nx, levp*ncld, my_max), &
      pdot(nxp, lev + 1, latpart), &
      vdmerdr(nxp, lev, my_max), vdzonlr(nxp, lev, my_max), &
      vdmerd(nxp, lev, my_max), vdzonl(nxp, lev, my_max), &
      ddtemp(nxp, lev, my_max), qvadv(nxp, lev*ncld, my_max), &
      diveng(nxp, lev, my_max), pten(nxp, lev, my_max)

   !byl      real cc(nx+2,levp,1+ncld,my_max)
   real(kind=RTYPE) cc(nx + 2, levp, 1, my_max), dummy
   !byl,wss(levp,2,1+ncld,jtrun,jtmax)
   !
   integer ierr
   real(kind=RTYPE) dta, ONE
   real(kind=RTYPE) ww1(nx, my_max)

   !for 2dMPI
   real(kind=RTYPE) phiten1(lev, 2, jtrun, jtmax)

   logical forward
   integer async_id
   integer(kind=cuda_stream_kind) stream, matmul_stream(jtmax)
   type(cudaEvent) :: spread_event, pack_event
   type(cublashandle) :: handle

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   forward = .false.
   !$acc enter data create(sdpbl, deldm, ww1, &
   !$acc& ttm_sl, pten_sl, uum_sl, vvm_sl, qm_sl, &
   !$acc& pdot, vdmerdr, vdzonlr, vdmerd, vdzonl, &
   !$acc& ddtemp, qvadv, diveng, pten &
   !$acc& ) async(async_id)

   !
   !  global mean tempertures (tbar)
   !
   sqhaf = sqrt(0.5)
   ONE = 1.
   !byl      tbar=0. ; qbar=0.
   !     do 160 m = 1, mlistnum
   !      mf=mlist(m)
   !      if ( mf.eq.1) then
   !       do 161 k=1,lev
   !         tbar(k)= sqhaf*temnow(k,1,1,m)
   ! 161   continue
   !       do 162 k=1,lev*ncld
   !         qbar(k)= sqhaf*qnow(k,1,1,m)
   ! 162   continue
   !      endif
   ! 160 continue
   !

   !2dMPI
   !byl      do 160 m = 1, mlistnum
   !byl       mf=mlist(m)
   !byl       if ( mf.eq.1) then
   !byl        do 161 k=1,levp
   !byl          KK=Llist(k)
   !byl          tbar(KK)= sqhaf*temnow(k,1,1,m)
   !byl  161   continue
   !byl        do 162 n=1,ncld
   !byl           nk=(n-1)*levp
   !byl           nL=(n-1)*lev
   !byl        do 162 k=1,levp
   !byl           KK=nk+k
   !byl           KL=nL+Llist(k)
   !byl          qbar(KL)= sqhaf*qnow(KK,1,1,m)
   !byl  162   continue
   !byl       endif
   !byl  160 continue

   !byl      call mpe_global_sum(tbar,lev,mpe_double)
   !byl      call mpe_global_sum(qbar,lev*ncld,mpe_double)
   !
   !
   !$acc parallel loop collapse(3) private(mf) async(async_id)
   do m = 1, mlistnum
      do n = 1, jtrun
         do k = 1, levp
            mf = mlist(m)
            if (n .ge. mf) then
               if (k .eq. 1) then
                  plten(n, m, 1) = 0.0
                  plten(n, m, 2) = 0.0
               end if
               divten(k, 1, n, m) = 0.0
               divten(k, 2, n, m) = 0.0
               vorten(k, 1, n, m) = 0.0
               vorten(k, 2, n, m) = 0.0
               temten(k, 1, n, m) = 0.0
               temten(k, 2, n, m) = 0.0
               !ncld qten(k,1,n,m) = 0.0
               hldten(k, 1, n, m) = 0.0
               hldten(k, 2, n, m) = 0.0
               phiten(k, 1, n, m) = 0.0
               phiten(k, 2, n, m) = 0.0
            end if
         end do
      end do
   end do
   !
  !!      do 183 m = 1, mlistnum
  !!         mf=mlist(m)
  !!      do 183 n = mf, jtrun
  !!      do 183 k = 1, levp*ncld*2
  !!      qten(k,1,n,m) = 0.0
  !!  183 continue
   !
   !
   dta = 0.5*dt
   !
   !$acc parallel loop collapse(3) private(j, nxj) async(async_id)
   do jj = 1, jlistnum
      do k = 1, lev*ncld
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            if (i .le. nxj) then
               if (k .le. lev) then
                  vdmerd(i, k, jj) = 0.
                  vdzonl(i, k, jj) = 0.
                  ddtemp(i, k, jj) = 0.
                  pten(i, k, jj) = 0.

                  up(i, k, jj) = ut(i, k, jj)
                  vp(i, k, jj) = vt(i, k, jj)
                  ttp(i, k, jj) = tt(i, k, jj)
               end if
               if (k .le. lev + 1) pdot(i, k, jj) = 0.
               qvadv(i, k, jj) = 0.
               qm(i, k, jj) = qt(i, k, jj)
            end if
         end do
      end do
     !!        do i = 1, nxj
     !!          ptp(i,jj)= pt(i,jj)
     !!        enddo
   end do
   !ch>
   ! transpose partial to full: ut -> ut_sl, vt -> vt_sl, ttm -> ttm_sl, qp -> qm_sl

   !#ifdef MULTIPLE
   !      call mpe2d_transpose_ndsl_p2f_multi(ut   ,vt   ,ttp   ,dummy,dummy,qp,    &
   !                                    ut_sl,vt_sl,ttm_sl,dummy,dummy,qm_sl, &
   !                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm,4)
   !#else
   call mpe2d_transpose_ndsl_p2f_gpu(ut, ut_sl, &
                                     nxp, nx, levf, levp, 1, myf, my_max, &
                                     jlistnum, jlen, nsizex, nccl_row_comm)
   call mpe2d_transpose_ndsl_p2f_gpu(vt, vt_sl, &
                                     nxp, nx, levf, levp, 1, myf, my_max, &
                                     jlistnum, jlen, nsizex, nccl_row_comm)
   call mpe2d_transpose_ndsl_p2f_gpu(ttp, ttm_sl, &
                                     nxp, nx, levf, levp, 1, myf, my_max, &
                                     jlistnum, jlen, nsizex, nccl_row_comm)
   call mpe2d_transpose_ndsl_p2f_gpu(qm, qm_sl, &
                                     nxp, nx, levf, levp, ncld, myf, my_max, &
                                     jlistnum, jlen, nsizex, nccl_row_comm)
   !#endif

  !!    call mpe2d_unify_nx(pt_sl,pt)
   !$acc parallel loop collapse(3) private(j, nxj) async(async_id)
   do jj = 1, jlistnum
      do k = 1, levp
         do i = 1, nx
            j = jlist1(jj)
            nxj = nxdef(j)
            if (i .le. nxj) then
               uum_sl(i, k, jj) = ut_sl(i, k, jj)
               vvm_sl(i, k, jj) = vt_sl(i, k, jj)
            end if
         end do
      end do
   end do

  !!    ptp_sl=pt_sl
   !ch<

   !
   !     Semi-Lagrangian
   !       Horizontal Advection
   call ndslfv_monoadvh2_gpu_refactor(ttm_sl, qm_sl, pten_sl, uum_sl, vvm_sl, nxdef, &
                                      dta, levp)

   !ch>
   ! transpose full to partial: ttm_sl -> ddtemp,  pten_sl -> pten, uum_sl -> vdzonl
   !                            vvm_sl -> vdmerd,  qm_sl -> qvadv

   !#ifdef MULTIPLE
   !      call mpe2d_transpose_ndsl_f2p_multi(ttm_sl,pten_sl,uum_sl,vvm_sl,qm_sl, &
   !                                    ddtemp   ,pten   ,vdzonl   ,vdmerd   ,qvadv   , &
   !                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
   !#else
   call mpe2d_transpose_ndsl_f2p_gpu(ttm_sl, ddtemp, &
                                     nxp, nx, levf, levp, 1, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)
   call mpe2d_transpose_ndsl_f2p_gpu(pten_sl, pten, &
                                     nxp, nx, levf, levp, 1, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)
   call mpe2d_transpose_ndsl_f2p_gpu(uum_sl, vdzonl, &
                                     nxp, nx, levf, levp, 1, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)
   call mpe2d_transpose_ndsl_f2p_gpu(vvm_sl, vdmerd, &
                                     nxp, nx, levf, levp, 1, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)
   call mpe2d_transpose_ndsl_f2p_gpu(qm_sl, qvadv, &
                                     nxp, nx, levf, levp, ncld, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)
   !#endif

   !ch<
   !
   !  p**capa quantities
   !
   call prexp_hybrid_cwb_gpu_refactor(nxjp, nxp, lev, ptop, sigma, pt, &
                                      pk, pk2, plt)
   !! input: pt
   !! output: pk, pk2, plt
   !
   ! Calculate Vertical velocity & Stream Functions
   !
   call gridnl_hybrid_ndsl_gpu_refactor(nxjp, nxp, lev, ncld, &
                                        cp, radsq, ut, vt, rdiv, tt, &
                                        qt, phi, pt, dtpl, dlpl, sinl, &
                                        pk, pk2, dsigma, sigma, onocos, cor, &
                                        diveng, vdmerdr, vdzonlr, pten, &
                                        deldm, sdpbl, sd, pdot, vvel, &
                                        sgeo)
   !! input: ut, vt, rdiv, tt, phi, pk, pk2, qt, dtpl, dlpl, pt, pten
   !! output: deldm, sdpbl, diveng, vdmerdr, vdzonlr, sd, vvel, pdot

   ! ------------------------------------------------------------
   !
   !CWB2021 for single precision test
   call joinrs_gpu(cc_cg, diveng, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
   call tranrs_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, polyf, weight, &
                              cc_cg, hldten, 1, nsizey, &
                              gwk1_cg, ws_cg, wc_cg, wcc_fk_cg, fj_weight_cg)

   call trngra3_gpu_cuda_graph(jtrun, jtmax, nx, levp, my, my_max, &
                               cim, polyf, dpolyf, &
                               hldten, dlphi, dtphi, nsizey, &
                               cc_cg, gwk1_cg, ws_cg, wc_cg, wcc_fk_cg)
   ! ------------------------------------------------------------
   !
   !       update all horizontal informations
   !
   !$acc parallel loop collapse(3) private(j, nxj) async(async_id)
   do jj = 1, jlistnum
      do k = 1, lev
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            if (i .le. nxj) then
               vdmerdr(i, k, jj) = vdmerdr(i, k, jj) - dtphi(i, k, jj)/radsq/onocos(j)! &
               !                             -ut(i,k,jj)*cor(j)-(ut(i,k,jj)*ut(i,k,jj)      &
               !                             +vt(i,k,jj)*vt(i,k,jj))*onocos(j)*sinl(j)

               vdzonlr(i, k, jj) = vdzonlr(i, k, jj) - dlphi(i, k, jj)/radsq!           &
               !                             +vt(i,k,jj)*cor(j)
            end if
         end do
      end do
   end do !jj = 1,jlistnum
   call ndslfv_update_gpu(nxjp, vdzonl, vdmerd, vdzonlr, vdmerdr, dta, forward)

   !CWB2021 ndsl single precision test
   !
   !
   !       Vertical Advection
   !
   call ndslfv_monoadvv_gpu(ddtemp, qvadv, vdzonl, vdmerd, pdot, pt, nxjp, dta, forward)

   !CWB2021 ndsl single precision test

   !
   !$acc parallel loop collapse(3) private(j, nxj) async(async_id)
   do jj = 1, jlistnum
      do k = 1, lev
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            if (i .le. nxj) then
               vdzonl(i, k, jj) = (vdzonl(i, k, jj) - up(i, k, jj))/dt
               vdmerd(i, k, jj) = (vdmerd(i, k, jj) - vp(i, k, jj))/dt
               ddtemp(i, k, jj) = (ddtemp(i, k, jj) - ttp(i, k, jj))/dt
            end if
         end do
      end do
   end do
   !
   !  combine non-linear grid point terms via gaussian quadrature
   !
   call joinrs_gpu(cc_cg, ddtemp, dummy, dummy, dummy, &
                   nx, my_max, lev, jlistnum, 1, 1)
   call tranrs_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, polyf, weight, &
                              cc_cg, temten, 1, nsizey, &
                              gwk1_cg, ws_cg, wc_cg, wcc_fk_cg, fj_weight_cg)
  !!      call joinrs(cc,ddtemp,qvadv,dummy,dummy,nx,my_max,lev           &
  !!                 ,jlistnum,2,ncld)
  !!      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc         &
  !!                 ,wss,1+ncld,nsize)
  !!      call ujoinrs(wss,temten,qten,dummy,dummy,jtrun,jtmax,lev        &
  !!                 ,mlistnum,2,ncld)
   call mpe2d_unify_nx_gpu(ww1, deldm)
   call tranrs1_gpu(jtrun, jtmax, nx, my, my_max, polyf, weight, ww1, &
                    plten, nsizey, cc_cg, gwk1_cg)
   !
   !CWB2021 for single precision test
   call trandv_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, vdzonl, vdmerd, &
                              weight, cim, onocos, polyf, dpolyf, vorten, divten, nsizey, &
                              cc_cg, gwk1_cg, ws_cg, wc_cg(1, 1), wc_cg(1, 2), &
                              wcc_fk_cg, fj_weight_cg(1, 1), fj_weight_cg(1, 2))
   !
   !  no tendency of zero mode
   !
   !$acc parallel loop collapse(2) private(mf) async(async_id)
   do m = 1, mlistnum
      do k = 1, levp
         mf = mlist(m)
         if (mf .eq. 1) then
            vorten(k, 1, 1, m) = 0.
            vorten(k, 2, 1, m) = 0.
            divten(k, 1, 1, m) = 0.
            divten(k, 2, 1, m) = 0.
            temten(k, 1, 1, m) = 0.
            temten(k, 2, 1, m) = 0.
         end if
      end do
   end do
   !
   !  compute phiten from temden and plten
   !
   mlst = ilist(1)
   !$acc serial
   if (mlst .ne. 0) then
      plten(1, mlst, 1) = 0.
      plten(1, mlst, 2) = 0.
   end if
   !$acc end serial
   !
   !       do 100 m =1,mlistnum
   !        mf=mlist(m)
   !       do 100 n  = mf,jtrun
   !       do 99 k = 1, lev
   !         phiten(k,1,n,m) = spalm(k)*plten(n,m,1)
   !         phiten(k,2,n,m) = spalm(k)*plten(n,m,2)
   !       do 98 l = 1, lev
   !         phiten(k,1,n,m) = phiten(k,1,n,m)+arrhyd(k,l)*temten(l,1,n,m)
   !         phiten(k,2,n,m) = phiten(k,2,n,m)+arrhyd(k,l)*temten(l,2,n,m)
   !98     continue
   !99     continue
   !100    continue
   !
   !ch>
   call mpe2d_unify_lev_gpu(temten, temten1, lev, levp, jtrun, jtmax, mlistnum, &
                            nsizex, nccl_row_comm)

   if (.not. cg_created) then
      CUDACHECK(cudaEventCreate(spread_event))
      CUDACHECK(cudaEventCreate(pack_event))
      handle = cublasGetHandle()
      do m = 1, mlistnum
         matmul_stream(m) = acc_get_cuda_stream(m + 1)
      end do

      call accx_async_begin_capture(async_id)

      !$acc parallel loop collapse(3) private(mf) async(async_id)
      do m = 1, mlistnum
         do n = 1, jtrun
            do k = 1, lev
               mf = mlist(m)
               if (n .ge. mf) then
                  phiten1(k, 1, n, m) = spalm(k)*plten(n, m, 1)
                  phiten1(k, 2, n, m) = spalm(k)*plten(n, m, 2)
               end if
            end do
         end do
      end do

      CUDACHECK(cudaEventRecord(spread_event, stream))
      do m = 1, mlistnum
         CUDACHECK(cudaStreamWaitEvent(matmul_stream(m), spread_event, 0))
         mf = mlist(m)
         llistnum_fj = 2*(jtrun - mf + 1)
         ierr = cublasSetStream(handle, matmul_stream(m))

         !$acc host_data use_device(arrhyd, temten1, phiten1)
#ifdef SP
         call sgemm('N', 'N', lev, llistnum_fj, lev, &
                    ONE, arrhyd, lev, temten1(1, 1, mf, m), lev, &
                    ONE, phiten1(1, 1, mf, m), lev)
#else
         call dgemm('N', 'N', lev, llistnum_fj, lev, &
                    ONE, arrhyd, lev, temten1(1, 1, mf, m), lev, &
                    ONE, phiten1(1, 1, mf, m), lev)
#endif
         !$acc end host_data
         CUDACHECK(cudaEventRecord(pack_event, matmul_stream(m)))
         CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))
      end do

      !$acc parallel loop collapse(3) private(mf, kk) async(async_id)
      do m = 1, mlistnum
         do n = 1, jtrun
            do k = 1, levp
               mf = mlist(m)
               if (n .ge. mf) then
                  kk = Llist(k)
                  phiten(k, 1, n, m) = phiten1(kk, 1, n, m)
                  phiten(k, 2, n, m) = phiten1(kk, 2, n, m)
               end if
            end do
         end do
      end do
      !ch<

      call accx_async_end_capture(async_id, cg_graph)
      cg_created = .true.
      CUDACHECK(cudaEventDestroy(spread_event))
      CUDACHECK(cudaEventDestroy(pack_event))
   end if
   call accx_graph_launch(cg_graph, async_id)

   !$acc exit data delete(sdpbl, deldm, ww1, &
   !$acc& ttm_sl, pten_sl, uum_sl, vvm_sl, qm_sl, &
   !$acc& pdot, vdmerdr, vdzonlr, vdmerd, vdzonl, &
   !$acc& ddtemp, qvadv, diveng, pten &
   !$acc& ) async(async_id)

   return
end subroutine tendget_gpu
