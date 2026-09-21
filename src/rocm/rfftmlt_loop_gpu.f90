! ROCm replacement for rfftmlt_loop_identical_cuda_graph (src/nvidia/cufft_wrapper.f90).
! acc2omp splices this over the translated original.
!
! Why: the original creates ONE hipfft PLAN PER LATITUDE (jlistnum = 768) for
! every distinct (jump, batch, isign), so that each latitude's FFT can run on
! its own stream inside a captured graph. On rocFFT a plan of this shape holds
! ~22 MB of internal state (measured 2026-09-16: 768 plans = 16.65 GB, freed
! only by hipfftDestroy), and the model builds ~10 such sets -> ~170 GB of
! VRAM that never comes back, which is the OUT_OF_RESOURCES at integration
! step 2. A plan depends only on (jump, nxj, m, isign), never on the latitude,
! and the graph is not used on ROCm (layer 7 runs the FFT eagerly), so:
!   - one plan per distinct reduced length nxj, via the existing C cache
!     (find_fft_plan / cache_fft_plan are keyed exactly on that tuple);
!   - every latitude executes sequentially on the acc stream, one wait after
!     the loop, then the normalise-and-write-back as a synchronous OpenMP
!     kernel (the layer-7 ordering).
! Concurrency across latitudes is given up; it was already serialised by the
! per-jj cudaStreamSynchronize of the layer-7 eager path.
subroutine rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, jump, m, isign, graph)
   use const, only: RTYPE
   use cudafor
   use cufft
   use openacc
   use iso_c_binding
   use omp_lib
   implicit none

   integer :: jlistnum, jump, m, isign
   integer :: jj, j, nxj
   real(kind=RTYPE), dimension(jump, m, *) :: cc, gwk1 ! Present on device
   real, dimension(4096, *) :: trigsj
   integer, dimension(19, *) :: ifaxj
   integer, dimension(*) :: jlist1, nxdef
   type(cudaGraph) :: graph
   integer(kind=8) :: work_size   ! int_ptr_kind() is rewritten by acc2omp for the original; spliced text is not
   integer :: async_id, istat, i, k
   integer(kind=cuda_stream_kind) :: stream
   integer(4) :: plan_id, plan_jj(jlistnum)
   integer :: ncreated
   ! ---- LRU cache of plan SETS (one set per (jump, m, isign); one plan per
   ! distinct nxj inside a set). rocFFT costs a flat ~20 MB per plan, the
   ! reduced grid has 384 distinct lengths, and the model touches ~9 sets per
   ! time step: unlimited caching is ~70 GB and was the OUT_OF_RESOURCES.
   ! Budget = GEPS_FFT_PLAN_SETS sets (default 8: the model uses 7 distinct
   ! (jump,m,isign) sets per step, ~7 GB each; with 3-4 the LRU thrashed and
   ! rebuilt ~8k plans per step = ~740 s/step, 2026-09-18); the least recently
   ! used set is destroyed when a new one is needed. Re-creating a set costs
   ! ~4 s warm (42 s cold, rocFFT compiles kernels at run time and caches
   ! them on disk), so this trades speed for fitting in 192 GB; a structural
   ! fix (fewer distinct lengths, or one plan covering several) comes later.
   integer, parameter :: MAXSETS = 16, NMAX = 4096
   integer, save :: nsets = 0, clock = 0, budget = -1
   integer, save :: set_key(3, MAXSETS), set_last(MAXSETS), set_count(MAXSETS)
   integer(4), save :: set_plan(0:NMAX, MAXSETS)
   integer :: is, iold, n
   character(len=32) :: envbuf
   ! ---- 2026-09-19: the per-latitude exec loop is host-launch bound (prof8:
   ! ~31k rocFFT kernel launches per step, GPU 30% busy in this phase), so the
   ! latitudes are issued from GEPS_FFT_THREADS host threads (default 4), each
   ! on its own non-blocking stream; every plan has its own work area
   ! (auto-allocation on), so plans never share state. Same kernels, same
   ! data, same order per plan -> bit-identical.
   integer, parameter :: MAXTHR = 16
   integer, save :: nthr = -1
   integer(kind=cuda_stream_kind), save :: thr_stream(MAXTHR)
   integer :: it, dev
   interface
      function geps_hip_set_device(id) bind(C, name="geps_hip_set_device")
         import c_int
         integer(c_int), value :: id
         integer(c_int) :: geps_hip_set_device
      end function
      function geps_hip_get_device() bind(C, name="geps_hip_get_device")
         import c_int
         integer(c_int) :: geps_hip_get_device
      end function
   end interface

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   if (nthr .lt. 0) then
      nthr = 4
      call get_environment_variable("GEPS_FFT_THREADS", envbuf, status=istat)
      if (istat .eq. 0) read (envbuf, *, iostat=istat) nthr
      if (istat .ne. 0 .or. nthr .lt. 1) nthr = 4
      if (nthr .gt. MAXTHR) nthr = MAXTHR
      do it = 1, nthr
         istat = cudaStreamCreatewithFlags(thr_stream(it), cudastreamnonblocking)
      end do
      print *, '[rfftmlt_loop rocm] fft host threads =', nthr
   end if
   if (budget .lt. 0) then
      budget = 8
      call get_environment_variable("GEPS_FFT_PLAN_SETS", envbuf, status=istat)
      if (istat .eq. 0) read (envbuf, *, iostat=istat) budget
      if (istat .ne. 0 .or. budget .lt. 1) budget = 8
      if (budget .gt. MAXSETS) budget = MAXSETS
      print *, '[rfftmlt_loop rocm] plan-set budget =', budget
   end if

   ! find this (jump, m, isign) set, or make room and open a new one
   clock = clock + 1
   is = 0
   do i = 1, nsets
      if (set_key(1, i) .eq. jump .and. set_key(2, i) .eq. m .and. set_key(3, i) .eq. isign) is = i
   end do
   if (is .eq. 0) then
      if (nsets .lt. budget) then
         nsets = nsets + 1
         is = nsets
      else
         iold = 1
         do i = 2, nsets
            if (set_last(i) .lt. set_last(iold)) iold = i
         end do
         do n = 0, NMAX
            if (set_plan(n, iold) .ne. -1) istat = cufftDestroy(set_plan(n, iold))
         end do
         print *, '[rfftmlt_loop rocm] evicted set (jump,m,isign)=', set_key(:, iold), &
            ' plans=', set_count(iold), ' for (', jump, m, isign, ')'
         is = iold
      end if
      set_key(1, is) = jump; set_key(2, is) = m; set_key(3, is) = isign
      set_plan(:, is) = -1
      set_count(is) = 0
   end if
   set_last(is) = clock

   ! one plan per distinct nxj within the set
   ncreated = 0
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      if (nxj .gt. NMAX) then
         print *, '[rfftmlt_loop rocm] nxj exceeds NMAX', nxj
         stop 1
      end if
      plan_id = set_plan(nxj, is)
      if (plan_id .eq. -1) then
         call fft_create_plan(plan_id, 1)
         call fft_make_plan(1, jump, nxj, m, isign, plan_id, work_size)
         set_plan(nxj, is) = plan_id
         set_count(is) = set_count(is) + 1
         ncreated = ncreated + 1
      end if
      plan_jj(jj) = plan_id
   end do
   if (ncreated .gt. 0) print *, '[rfftmlt_loop rocm] new plans:', ncreated, &
      ' (jump,m,isign)=', jump, m, isign, ' sets live=', nsets

   if (mod(jump, 2) .ne. 0) then
      print *, 'fft jump is odd, CWB obsoleted, jump=', jump
      return
   end if

   ! the FFT reads/writes cc/gwk1 that OpenMP kernels produced synchronously;
   ! rocFFT itself is queued on `stream`, so one wait after the loop suffices
   dev = geps_hip_get_device()
   !$omp target data use_device_addr(cc, gwk1)
   !$omp parallel num_threads(nthr) private(jj, it, istat)
   it = omp_get_thread_num() + 1
   ! HIP's current device is per host thread; rank r runs on device r
   if (it .gt. 1) istat = geps_hip_set_device(dev)
   !$omp do schedule(static)
   do jj = 1, jlistnum
      istat = cufftSetStream(plan_jj(jj), thr_stream(it))
      if (isign .eq. 1) then
#ifdef SP
         istat = cufftExecC2R(plan_jj(jj), gwk1(1, 1, jj), cc(1, 1, jj))
#else
         istat = cufftExecZ2D(plan_jj(jj), gwk1(1, 1, jj), cc(1, 1, jj))
#endif
      else
#ifdef SP
         istat = cufftExecR2C(plan_jj(jj), cc(1, 1, jj), gwk1(1, 1, jj))
#else
         istat = cufftExecD2Z(plan_jj(jj), cc(1, 1, jj), gwk1(1, 1, jj))
#endif
      end if
      if (istat .ne. 0) print *, '[rfftmlt_loop rocm] hipfftExec failed', istat, ' jj=', jj
   end do
   !$omp end do
   !$omp end parallel
   !$omp end target data
   do it = 1, nthr
      istat = cudaStreamSynchronize(thr_stream(it))
   end do
   call geps_acc_wait(async_id)

   if (isign .ne. 1) then
      !$omp target teams distribute parallel do collapse(2) private(j, nxj, i)
      do jj = 1, jlistnum
         do k = 1, m
            j = jlist1(jj)
            nxj = nxdef(j)
            do i = 1, nxj + 2
               cc(i, k, jj) = gwk1(i, k, jj)/float(nxj)
            end do
         end do
      end do
   end if
end subroutine rfftmlt_loop_identical_cuda_graph
