!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine trngra_gpu(jtrun, jtmax, nx, my, my_max, cim, poly, dpoly, s &
                      , dlpl, dtpl, nsize, cc, gwk1)
!
!  subroutine to transform spectral terrain pressure to grid point
!  fields of zonal and meridional derivatives of terrain pressure
!
! *** input ***
!
!  jtrun: zonal wavenumber resolution limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  cim: zonal wavenumber array
!  poly: legendre polynomials
!  dpoly: d(poly)/d(sin(lat))
!  s: spectral coefficient array
!  nsiz: pe number
!
! **** output ****
!
!  dlpl: d(pt)/d(longitude)
!  dtpl: d(pt)/d(sin(lat))
!
! ****************************************************
!
   use const, only: RTYPE
   use index
   use fftcom
   use spec_cuda_graph, only: fft_cg => trngra_fft_cg
   use openacc
   use cudafor

   implicit none

   integer jtrun, jtmax, nx, my, my_max, nsize
   integer myhalf, m, mf, l, j, jj, i, jtrunj, mm, mp, mlst, nxj

   real(kind=RTYPE) poly(jtrun, my/2, jtmax), dpoly(jtrun, my/2, jtmax) ! Present on device
   real(kind=RTYPE) s(jtrun, jtmax, 2) ! Present on device
   real(kind=RTYPE) cim(jtmax) ! Present on device

   real(kind=RTYPE) dlpl(nxp, my_max), dtpl(nxp, my_max) ! Present on device
   real(kind=RTYPE) cc(nx + 2, 2, my_max)

   real(kind=RTYPE) gwk1(nx + 2, 2, my_max)

   real(kind=RTYPE) twcc_fk(my_max, jtmax*nsize, 2)
   real(kind=RTYPE) twdd_fk(my_max, jtmax*nsize, 2)

   real(kind=RTYPE) wcu_fk(my_max*nsize, jtmax, 2), wcv_fk(my_max*nsize, jtmax, 2)
   real(kind=RTYPE) wcu_t(my, 2), wcv_t(my, 2)
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream
   real(kind=RTYPE) wcu_t1, wcu_t2, wcv_t1, wcv_t2, rt
   integer sign_flip

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   myhalf = my/2

   !$acc data create(wcu_fk, wcv_fk, wcu_t, wcv_t, twcc_fk, twdd_fk) async(async_id)
   !$acc host_data use_device(wcu_fk, wcv_fk)
   CUDACHECK(cudaMemSetAsync(wcu_fk, real(0.0, RTYPE), size(wcu_fk), stream))
   CUDACHECK(cudaMemSetAsync(wcv_fk, real(0.0, RTYPE), size(wcv_fk), stream))
   !$acc end host_data
   !$acc parallel loop gang(dim:2) async(async_id)
   do m = 1, mlistnum
      mf = mlist(m)
      !$acc loop private(wcu_t1, wcu_t2, wcv_t1, wcv_t2)
      do j = 1, my
         wcu_t1 = 0.
         wcu_t2 = 0.
         wcv_t1 = 0.
         wcv_t2 = 0.
         if (j .le. myhalf) then
            if (mf .le. mtrundef(j)) then
               !$acc loop seq
               do l = mf, jtrun
                  wcu_t1 = wcu_t1 + s(l, m, 2)*(cim(m)*poly(l, j, m))
                  wcu_t2 = wcu_t2 - s(l, m, 1)*(cim(m)*poly(l, j, m))
                  wcv_t1 = wcv_t1 - s(l, m, 1)*dpoly(l, j, m)
                  wcv_t2 = wcv_t2 - s(l, m, 2)*dpoly(l, j, m)
               end do
            end if
         else
            jj = my - j + 1
            if (mf .le. mtrundef(jj)) then
               !$acc loop seq
               do l = mf, jtrun
                  sign_flip = (-1)**(l - mf)
                  wcu_t1 = wcu_t1 + sign_flip*s(l, m, 2)*(cim(m)*poly(l, jj, m))
                  wcu_t2 = wcu_t2 - sign_flip*s(l, m, 1)*(cim(m)*poly(l, jj, m))
                  wcv_t1 = wcv_t1 + sign_flip*s(l, m, 1)*dpoly(l, jj, m)
                  wcv_t2 = wcv_t2 + sign_flip*s(l, m, 2)*dpoly(l, jj, m)
               end do
            end if
         end if
         jj = jlist2(j)
         wcu_fk(jj, m, 1) = wcu_t1
         wcu_fk(jj, m, 2) = wcu_t2
         wcv_fk(jj, m, 1) = wcv_t1
         wcv_fk(jj, m, 2) = wcv_t2
      end do
   end do

   call mpe_transpose_rs1_sp_gpu(wcu_fk, twcc_fk, my_max, jtmax, 2, nsize, nccl_col_comm, async_id)
   call mpe_transpose_rs1_sp_gpu(wcv_fk, twdd_fk, my_max, jtmax, 2, nsize, nccl_col_comm, async_id)

   !$acc host_data use_device(gwk1)
   CUDACHECK(cudaMemSetAsync(gwk1, real(0.0, RTYPE), size(gwk1), stream))
   !$acc end host_data

   !$acc parallel loop gang async(async_id)
   do jj = 1, jlistnum
      j = jlist1(jj)
      jtrunj = mtrundef(j)
      !$acc loop vector
      do m = 1, jtrunj
         mm = 2*m - 1
         mp = mm + 1
         mlst = nlist(m)
         gwk1(mm, 1, jj) = twcc_fk(jj, mlst, 1)
         gwk1(mp, 1, jj) = twcc_fk(jj, mlst, 2)
         gwk1(mm, 2, jj) = twdd_fk(jj, mlst, 1)
         gwk1(mp, 2, jj) = twdd_fk(jj, mlst, 2)
      end do
   end do

   if (lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, jlistnum*2, 1)
   else
      if (.not. (fft_cg%created)) then
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, 2, 1, fft_cg%graph)
         CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
         fft_cg%created = .true.
      end if
      CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
   end if

   !$acc parallel loop async(async_id)
   do jj = 1, jlistnum
      j = jlist1(jj)
      dlpl(1:nxjlen(j), jj) = -cc(nxjstart(j):nxjend(j), 1, jj)
      dtpl(1:nxjlen(j), jj) = -cc(nxjstart(j):nxjend(j), 2, jj)
   end do
   !$acc end data
end subroutine trngra_gpu
