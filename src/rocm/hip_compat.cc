/* HIP runtime helpers for GEPS AMD GPU path (CUDA Fortran API subset). */

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <hipfft/hipfft.h>
#include <hipsolver/hipsolver.h>
#include <hipsparse/hipsparse.h>
#include <omp.h>
#include <rccl/rccl.h>

namespace {

std::mutex g_mu;
std::unordered_map<int, hipStream_t> g_streams;
std::vector<hipStream_t> g_all_streams; /* every stream this shim created */
/* streams that received work since the last geps_hip_wait_all(); ~350
   streams exist (prof8), so waiting on all of them cost 147k
   hipStreamSynchronize per step */
std::vector<hipStream_t> g_dirty;
hipStream_t use_stream(int64_t stream) {
  hipStream_t s = (hipStream_t)(uintptr_t)stream;
  std::lock_guard<std::mutex> lock(g_mu);
  for (hipStream_t d : g_dirty)
    if (d == s)
      return s;
  g_dirty.push_back(s);
  return s;
}
std::unordered_map<int, hipfftHandle> g_fft;
int g_fft_next = 1;
hipblasHandle_t g_blas = nullptr;
hipsparseHandle_t g_sparse = nullptr;
hipsolverHandle_t g_solver = nullptr;

hipStream_t stream_for(int async_id) {
  if (async_id <= 0)
    return nullptr; /* default stream */
  std::lock_guard<std::mutex> lock(g_mu);
  auto it = g_streams.find(async_id);
  if (it != g_streams.end())
    return it->second;
  hipStream_t s = nullptr;
  (void)hipStreamCreateWithFlags(&s, hipStreamNonBlocking);
  g_streams[async_id] = s;
  g_all_streams.push_back(s);
  return s;
}

/* Host ptr with an OpenMP mapping → device ptr. Already-device ptr is kept. */
void *device_ptr(void *p) {
  if (!p)
    return p;
  int d = omp_get_default_device();
  if (omp_target_is_present(p, d)) {
    void *m = omp_get_mapped_ptr(p, d);
    if (m)
      return m;
  }
  return p;
}

hipsolverSyevjInfo_t g_syevj = nullptr;

/* #region agent log */
void dbg_ndjson(const char *hyp, const char *loc, const char *msg, int rank,
                long long p0, long long p1, long long p2, int do_sync) {
  /* 2026-09-17: opt-in (GEPS_DBG_NDJSON=1). Each call does hipMemGetInfo +
   * a file append; the src/rocm cyclic_cell wrappers call it ~34x per
   * advection, which is measurable once the model runs at full speed. */
  static int ndjson_on = -1;
  if (ndjson_on < 0)
    ndjson_on = std::getenv("GEPS_DBG_NDJSON") ? 1 : 0;
  if (!ndjson_on)
    return;
  int dev = -1;
  (void)hipGetDevice(&dev);
  hipError_t peek = hipPeekAtLastError();
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  long long ms = (long long)ts.tv_sec * 1000LL + ts.tv_nsec / 1000000LL;
  size_t free_b = 0, tot_b = 0;
  (void)hipMemGetInfo(&free_b, &tot_b);
  FILE *f = fopen("/mlsteam/workspace/data/geps/.cursor/debug-a51a4c.log", "a");
  if (f) {
    fprintf(f,
            "{\"sessionId\":\"a51a4c\",\"runId\":\"pre-fix\",\"hypothesisId\":\"%s\","
            "\"location\":\"%s\",\"message\":\"%s\",\"data\":{\"rank\":%d,\"dev\":%d,"
            "\"peek\":%d,\"sync\":-1,\"p0\":%lld,\"p1\":%lld,\"p2\":%lld,"
            "\"freeB\":%zu,\"usedB\":%zu,\"totB\":%zu},"
            "\"timestamp\":%lld}\n",
            hyp, loc, msg, rank, dev, (int)peek, p0, p1, p2, free_b,
            tot_b > free_b ? tot_b - free_b : 0, tot_b, ms);
    fflush(f);
    fclose(f);
  }
  if (!do_sync)
    return;
  hipError_t syn = hipDeviceSynchronize();
  clock_gettime(CLOCK_REALTIME, &ts);
  ms = (long long)ts.tv_sec * 1000LL + ts.tv_nsec / 1000000LL;
  f = fopen("/mlsteam/workspace/data/geps/.cursor/debug-a51a4c.log", "a");
  if (f) {
    fprintf(f,
            "{\"sessionId\":\"a51a4c\",\"runId\":\"pre-fix\",\"hypothesisId\":\"%s\","
            "\"location\":\"%s\",\"message\":\"%s-sync\",\"data\":{\"rank\":%d,\"dev\":%d,"
            "\"peek\":%d,\"sync\":%d,\"p0\":%lld,\"p1\":%lld,\"p2\":%lld},"
            "\"timestamp\":%lld}\n",
            hyp, loc, msg, rank, dev, (int)peek, (int)syn, p0, p1, p2, ms);
    fflush(f);
    fclose(f);
  }
}
/* #endregion */

} /* namespace */

extern "C" {

/* #region agent log */
void geps_dbg_vram_(int *hyp, int *locid, long long *p0, long long *p1,
                    long long *p2) {
  char loc[64];
  std::snprintf(loc, sizeof(loc), "vram:%d", locid ? *locid : 0);
  const char *h = "A";
  int hv = hyp ? *hyp : 1;
  if (hv == 2)
    h = "B";
  else if (hv == 3)
    h = "C";
  else if (hv == 4)
    h = "D";
  else if (hv >= 5)
    h = "E";
  dbg_ndjson(h, loc, "vram", 0, p0 ? *p0 : 0, p1 ? *p1 : 0, p2 ? *p2 : 0, 0);
}

void geps_dbg_gpu_log_(int *hyp, int *locid, int *rank, long long *p0,
                       long long *p1, long long *p2) {
  const char *h = "E";
  switch (*hyp) {
  case 1:
    h = "A";
    break;
  case 2:
    h = "B";
    break;
  case 3:
    h = "C";
    break;
  case 4:
    h = "D";
    break;
  default:
    h = "E";
    break;
  }
  char loc[64];
  std::snprintf(loc, sizeof(loc), "initial_gpu:%d", *locid);
  dbg_ndjson(h, loc, "ckpt", *rank, *p0, *p1, *p2, 1);
}
/* #endregion */

int64_t geps_acc_get_stream(int async_id) {
  return (int64_t)(uintptr_t)stream_for(async_id);
}

void geps_acc_set_stream(int async_id, int64_t stream) {
  std::lock_guard<std::mutex> lock(g_mu);
  g_streams[async_id] = (hipStream_t)(uintptr_t)stream;
}

void geps_hip_wait(int async_id) {
  (void)hipStreamSynchronize(stream_for(async_id));
}

/* Device-wide wait. acc2omp inserts this after every cudaMemsetAsync /
   cudaMemcpyAsync run: those go to a HIP stream while the translated
   OpenMP kernels run on libomptarget's own queue, so without it a memset
   can land AFTER the kernel that fills the same buffer (2026-09-15: the
   NDSL pack in ndslfv_monoadvh2 lost about half of ddtemp that way). */
void geps_hip_wait_all(void) {
  // 2026-09-18: every HIP async op the model issues (memset/memcpy, RCCL,
  // hipfft, hipblas) goes through a shim function that takes the stream
  // (use_stream marks it dirty), and the translated OpenMP kernels are
  // synchronous. Waiting on the dirty streams is therefore equivalent to a
  // device-wide sync (~16k calls per step). 2026-09-19: waiting on every
  // stream the shim ever created was 350 syncs per call (prof8), hence the
  // dirty list. async_id <= 0 maps to the null stream, always waited.
  std::vector<hipStream_t> ss;
  {
    std::lock_guard<std::mutex> lock(g_mu);
    ss.swap(g_dirty);
  }
  (void)hipStreamSynchronize(nullptr);
  for (hipStream_t s : ss)
    if (s)
      (void)hipStreamSynchronize(s);
}

/* #region agent log: DBGMAP/DBGUNMAP print with the current free VRAM.
   Called with an implicit Fortran interface (character literal + real(8)
   expression), so it needs no declaration in the caller - that matters in
   `!$omp declare target` routines, where a USE after the directive is an
   error. flang passes the character length as a trailing size_t. */
void geps_dbgmap_(const char *tag, const double *gib, size_t taglen) {
  size_t f = 0, t = 0;
  (void)hipMemGetInfo(&f, &t);
  fprintf(stdout, " %.*s GiB= %.4f free= %.2f\n", (int)taglen, tag, *gib,
          (double)f / 1073741824.0);
  fflush(stdout);
}
/* #endregion */

int geps_hip_device_count() {
  int n = 0;
  (void)hipGetDeviceCount(&n);
  if (n <= 1)
    return n < 1 ? 1 : n;
  /* Collapse DPX partitions of the same PCI device (this lab: C5:00.0 / .1). */
  int unique = 0;
  int seen_bus[32], seen_dev[32], seen_dom[32];
  int nseen = 0;
  for (int i = 0; i < n && i < 32; ++i) {
    hipDeviceProp_t p{};
    if (hipGetDeviceProperties(&p, i) != hipSuccess)
      continue;
    int dup = 0;
    for (int k = 0; k < nseen; ++k) {
      if (seen_dom[k] == p.pciDomainID && seen_bus[k] == p.pciBusID &&
          seen_dev[k] == p.pciDeviceID) {
        dup = 1;
        break;
      }
    }
    if (!dup && nseen < 32) {
      seen_dom[nseen] = p.pciDomainID;
      seen_bus[nseen] = p.pciBusID;
      seen_dev[nseen] = p.pciDeviceID;
      nseen++;
      unique++;
    }
  }
  return unique < 1 ? 1 : unique;
}

int geps_hip_set_device(int id) {
  int use = id;
  if (geps_hip_device_count() <= 1)
    use = 0;
  else if (use < 0)
    use = 0;
  int rc = (int)hipSetDevice(use);
  omp_set_default_device(use);
  /* #region agent log */
  {
    int nraw = 0;
    hipDeviceProp_t p{};
    (void)hipGetDeviceCount(&nraw);
    (void)hipGetDeviceProperties(&p, use);
    dbg_ndjson("F", "hip_set_device", "pci", nraw, (long long)p.pciBusID,
               (long long)p.pciDeviceID, (long long)use, 0);
  }
  /* #endregion */
  return rc;
}

int geps_hip_get_device() {
  int id = 0;
  (void)hipGetDevice(&id);
  return id;
}

void geps_hip_error_string(int err, char *buf, int n) {
  const char *s = hipGetErrorString((hipError_t)err);
  if (!s)
    s = "unknown HIP error";
  std::snprintf(buf, n > 0 ? (size_t)n : 0, "%s", s);
}

int geps_hip_memcpy(void *dst, const void *src, int64_t bytes, int /*kind*/) {
  return (int)hipMemcpy(dst, src, (size_t)bytes, hipMemcpyDefault);
}

int geps_hip_memcpy_async(void *dst, const void *src, int64_t bytes, int kind,
                          int64_t stream) {
  // 2026-09-18 (rocprofv3): with hipMemcpyDefault the 14 D2D state copies in
  // intgrt_gpu ran at ~26 MB/s (13 s for a 343 MB array, 87 s for qm) and
  // cost ~190 s per step: the buffers come from libomptarget's allocator,
  // which HIP's pointer classification does not recognise as device memory,
  // so CLR takes the host-staged path. Explicit DtoD skips the lookup.
  // 2026-09-18 (rocprofv3 prof2): hipMemcpyDtoDAsync was no better -- 84
  // calls, 1127 s, 13 s per 343 MB, i.e. the copy still runs at ~26 MB/s.
  // HIP does not know libomptarget's HSA allocations, so it falls back to a
  // CPU memcpy through the large-BAR mapping of device memory. Let
  // libomptarget do the copy (hsa_amd_memory_async_copy on the device):
  // omp_target_memcpy is synchronous, which every caller tolerates (each
  // memcpy run is followed by geps_acc_wait_all anyway).
  if (kind == 3) {
    int dev = omp_get_default_device();
    return omp_target_memcpy(dst, src, (size_t)bytes, 0, 0, dev, dev);
  }
  return (int)hipMemcpyAsync(dst, src, (size_t)bytes, hipMemcpyDefault,
                             use_stream(stream));
}

int geps_hip_memset_async(void *dst, int64_t bytes, int64_t stream) {
  dst = device_ptr(dst);
  return (int)hipMemsetAsync(dst, 0, (size_t)bytes, use_stream(stream));
}

int geps_hip_memset_i32_async(void *dst, int val, int64_t count, int64_t stream) {
  dst = device_ptr(dst);
  return (int)hipMemsetD32Async((hipDeviceptr_t)dst, val, (size_t)count,
                                use_stream(stream));
}

int geps_hip_malloc_async(void **ptr, int64_t bytes, int64_t stream) {
  return (int)hipMallocAsync(ptr, (size_t)bytes, use_stream(stream));
}

int geps_hip_free_async(void *ptr, int64_t stream) {
  return (int)hipFreeAsync(ptr, use_stream(stream));
}

int geps_hip_event_create(void **ev) {
  hipEvent_t e = nullptr;
  int rc = (int)hipEventCreateWithFlags(&e, hipEventDisableTiming);
  *ev = (void *)e;
  return rc;
}

int geps_hip_event_destroy(void *ev) {
  return (int)hipEventDestroy((hipEvent_t)ev);
}

int geps_hip_event_record(void *ev, int64_t stream) {
  return (int)hipEventRecord((hipEvent_t)ev, use_stream(stream));
}

int geps_hip_stream_wait_event(int64_t stream, void *ev, unsigned int flags) {
  return (int)hipStreamWaitEvent(use_stream(stream), (hipEvent_t)ev,
                                 flags);
}

int geps_hip_stream_create_flags(int64_t *stream, unsigned int flags) {
  hipStream_t s = nullptr;
  int rc = (int)hipStreamCreateWithFlags(&s, flags);
  *stream = (int64_t)(uintptr_t)s;
  if (rc == 0) {
    std::lock_guard<std::mutex> lock(g_mu);
    g_all_streams.push_back(s);
  }
  return rc;
}

int geps_hip_stream_sync(int64_t stream) {
  return (int)hipStreamSynchronize(use_stream(stream));
}

int geps_hip_begin_capture(int64_t stream, int mode) {
  return (int)hipStreamBeginCapture(use_stream(stream),
                                    (hipStreamCaptureMode)mode);
}

int geps_hip_end_capture(int64_t stream, void **graph) {
  hipGraph_t g = nullptr;
  int rc = (int)hipStreamEndCapture(use_stream(stream), &g);
  *graph = (void *)g;
  return rc;
}

int geps_hip_graph_instantiate(void **exec, void *graph, int /*unused*/) {
  hipGraphExec_t e = nullptr;
  int rc = (int)hipGraphInstantiateWithFlags(&e, (hipGraph_t)graph, 0);
  *exec = (void *)e;
  return rc;
}

int geps_hip_graph_launch(void *exec, int64_t stream) {
  return (int)hipGraphLaunch((hipGraphExec_t)exec,
                             use_stream(stream));
}

int geps_blas_create(void **handle) {
  hipblasHandle_t h = nullptr;
  int rc = (int)hipblasCreate(&h);
  *handle = (void *)h;
  if (!g_blas)
    g_blas = h;
  return rc;
}

void *geps_blas_default() {
  if (!g_blas)
    hipblasCreate(&g_blas);
  return (void *)g_blas;
}

int geps_blas_set_stream(void *handle, int64_t stream) {
  return (int)hipblasSetStream((hipblasHandle_t)handle,
                               use_stream(stream));
}

/* #region agent log: report the device address OpenMP has mapped a host array
   to. Called from Fortran by reference, so no iso_c_binding is needed:
       call geps_dbg_mapped(nnmi_buf, 1)
   Compare the result with the C pointer printed by [dgemm]. */
extern "C" void geps_dbg_mapped_(double *p, int *tag) {
  int d = omp_get_default_device();
  int present = omp_target_is_present(p, d);
  void *mapped = present ? omp_get_mapped_ptr(p, d) : nullptr;
  fprintf(stderr, "[mapped] tag=%d host=%p is_present=%d device=%p\n",
          tag ? *tag : -1, (void *)p, present, mapped);
  fflush(stderr);
}
/* #endregion */

/* 2026-09-19: acc2omp's Legendre-loop rewrite turns the per-dgemm sync off
   for the duration of a dgemm-only loop and waits all dirty streams after it
   (geps_blas_defer_sync in geps_acc_wait.f90). */
static int g_blas_defer_sync = 0;
void geps_blas_defer_sync_set(int v) { g_blas_defer_sync = v; }

int geps_blas_dgemm(void *handle, int ta, int tb, int m, int n, int k,
                    double alpha, const double *a, int lda, const double *b,
                    int ldb, double beta, double *c, int ldc) {
  /* #region agent log: which buffer does dgemm actually write?
     nnmi_gpu's `bal` implies |wrk| ~ 1 while its x_out implies |wrk| ~ 1e-3.
     If the address dgemm gets from `use_device_addr` differs from the one the
     OpenMP kernel reads (printed by geps_dbg_mapped_), they are two buffers. */
  {
    /* Only nnmi_gpu's shape: dgemm(op, op, nn, 2, nn, ...) - the first few
       dgemm calls in the program come from elsewhere (m=144, n=1..6) and
       their pointers are not comparable. */
    static int seen = 0;
    if (n == 2 && m == k && seen < 6) {
      ++seen;
      fprintf(stderr, "[dgemm] m=%d n=%d k=%d  A=%p B=%p C=%p\n", m, n, k,
              (const void *)a, (const void *)b, (void *)c);
      fflush(stderr);
    }
  }
  /* #endregion */
  int rc = (int)hipblasDgemm((hipblasHandle_t)handle, (hipblasOperation_t)ta,
                             (hipblasOperation_t)tb, m, n, k, &alpha, a, lda, b,
                             ldb, &beta, c, ldc);
  /* 2026-09-16 (layer 15): hipBLAS queues on the handle's current stream and
     returns; translated OpenMP kernels run on libomptarget's own queue with
     no ordering against it. Sites outside the graph blocks (correct_gpu,
     tranuv1_gpu, dgemm_async's 20 callers) read the result in the very next
     kernel with no wait, and with LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0
     that showed up as a 1e-5 shift in sptend and NaN feeding the
     microphysics. Make every dgemm synchronous here, once, instead of
     patching each call site. 2026-09-18: wait on the handle's stream only
     (the dgemm is the last thing queued there); ~2.7k calls per step. */
  if (!g_blas_defer_sync) {
    hipStream_t hs = nullptr;
    if (hipblasGetStream((hipblasHandle_t)handle, &hs) == HIPBLAS_STATUS_SUCCESS)
      (void)hipStreamSynchronize(hs);
    else
      (void)hipDeviceSynchronize();
  }
  return rc;
}

int geps_fft_create(int *plan) {
  hipfftHandle h = nullptr;
  int rc = (int)hipfftCreate(&h);
  std::lock_guard<std::mutex> lock(g_mu);
  int id = g_fft_next++;
  g_fft[id] = h;
  *plan = id;
  return rc;
}

static hipfftHandle fft_of(int plan) {
  std::lock_guard<std::mutex> lock(g_mu);
  auto it = g_fft.find(plan);
  return it == g_fft.end() ? nullptr : it->second;
}

int geps_fft_destroy(int plan) {
  hipfftHandle h = fft_of(plan);
  if (!h)
    return 0;
  int rc = (int)hipfftDestroy(h);
  std::lock_guard<std::mutex> lock(g_mu);
  g_fft.erase(plan);
  return rc;
}

int geps_fft_set_auto_alloc(int plan, int auto_alloc) {
  return (int)hipfftSetAutoAllocation(fft_of(plan), auto_alloc);
}

int geps_fft_make_plan_many(int plan, int rank, int n, int inembed, int istride,
                            int idist, int onembed, int ostride, int odist,
                            int ffttype, int batch, int64_t *work_size) {
  size_t ws = 0;
  int nn = n;
  int ine = inembed;
  int one = onembed;
  /* #region agent log: plans are created per latitude (768 per distinct
     batch size) with auto-allocation ON, and every routine's FIRST call
     costs 15-45 GiB of VRAM that never comes back (2026-09-16 trace). Log
     the reported work size and the real hipMemGetInfo delta per plan set. */
  size_t free0 = 0, tot0 = 0;
  (void)hipMemGetInfo(&free0, &tot0);
  /* #endregion */
  int rc = (int)hipfftMakePlanMany(fft_of(plan), rank, &nn, &ine, istride,
                                   idist, &one, ostride, odist,
                                   (hipfftType)ffttype, batch, &ws);
  if (work_size)
    *work_size = (int64_t)ws;
  /* #region agent log */
  {
    static long long nplans = 0, ws_sum = 0, vram_sum = 0;
    size_t free1 = 0, tot1 = 0;
    (void)hipMemGetInfo(&free1, &tot1);
    ++nplans;
    ws_sum += (long long)ws;
    vram_sum += (long long)free0 - (long long)free1;
    if (nplans <= 3 || nplans % 256 == 0)
      fprintf(stderr,
              "[fftplan] #%lld n=%d batch=%d ws=%.1f MB vram_delta=%.1f MB  "
              "cum ws=%.2f GB cum vram=%.2f GB free=%.1f GB\n",
              nplans, n, batch, ws / 1048576.0,
              ((double)free0 - (double)free1) / 1048576.0,
              ws_sum / 1073741824.0, vram_sum / 1073741824.0,
              free1 / 1073741824.0);
  }
  /* #endregion */
  return rc;
}

int geps_fft_set_stream(int plan, int64_t stream) {
  return (int)hipfftSetStream(fft_of(plan), use_stream(stream));
}

int geps_fft_set_work_area(int plan, void *work) {
  return (int)hipfftSetWorkArea(fft_of(plan), work);
}

int geps_fft_exec_z2d(int plan, void *in, void *out) {
  int rc = (int)hipfftExecZ2D(fft_of(plan), (hipfftDoubleComplex *)in,
                            (double *)out);
  /* 2026-09-18: no per-exec device sync (layer 15 was redundant here: the
   * caller waits on the FFT stream once after its 384-latitude loop; the
   * per-exec sync cost ~3k x 150 us per step). */
  return rc;
}
int geps_fft_exec_d2z(int plan, void *in, void *out) {
  int rc = (int)hipfftExecD2Z(fft_of(plan), (double *)in,
                            (hipfftDoubleComplex *)out);
  /* 2026-09-18: no per-exec device sync (layer 15 was redundant here: the
   * caller waits on the FFT stream once after its 384-latitude loop; the
   * per-exec sync cost ~3k x 150 us per step). */
  return rc;
}
int geps_fft_exec_c2r(int plan, void *in, void *out) {
  int rc = (int)hipfftExecC2R(fft_of(plan), (hipfftComplex *)in, (float *)out);
  /* 2026-09-18: no per-exec device sync (layer 15 was redundant here: the
   * caller waits on the FFT stream once after its 384-latitude loop; the
   * per-exec sync cost ~3k x 150 us per step). */
  return rc;
}
int geps_fft_exec_r2c(int plan, void *in, void *out) {
  int rc = (int)hipfftExecR2C(fft_of(plan), (float *)in, (hipfftComplex *)out);
  /* 2026-09-18: no per-exec device sync (layer 15 was redundant here: the
   * caller waits on the FFT stream once after its 384-latitude loop; the
   * per-exec sync cost ~3k x 150 us per step). */
  return rc;
}

int geps_solver_create(void **handle) {
  hipsolverHandle_t h = nullptr;
  int rc = (int)hipsolverDnCreate(&h);
  *handle = (void *)h;
  if (!g_solver)
    g_solver = h;
  return rc;
}

int geps_solver_destroy(void *handle) {
  return (int)hipsolverDnDestroy((hipsolverHandle_t)handle);
}

int geps_solver_set_stream(void *handle, int64_t stream) {
  return (int)hipsolverDnSetStream((hipsolverHandle_t)handle,
                                   use_stream(stream));
}

int geps_solver_syevj_create(void **info) {
  hipsolverSyevjInfo_t p = nullptr;
  int rc = (int)hipsolverDnCreateSyevjInfo(&p);
  *info = (void *)p;
  return rc;
}

int geps_solver_syevj_set_tol(void *info, double tol) {
  return (int)hipsolverDnXsyevjSetTolerance((hipsolverSyevjInfo_t)info, tol);
}

int geps_solver_syevj_set_sort(void *info, int sort) {
  return (int)hipsolverDnXsyevjSetSortEig((hipsolverSyevjInfo_t)info, sort);
}

int geps_solver_dsyevd_buf(void *handle, int jobz, int uplo, int n,
                           double *A, int lda, double *W, int *lwork) {
  return (int)hipsolverDnDsyevd_bufferSize((hipsolverHandle_t)handle,
                                           (hipsolverEigMode_t)jobz,
                                           (hipblasFillMode_t)uplo, n, A, lda,
                                           W, lwork);
}

int geps_solver_dsyevd(void *handle, int jobz, int uplo, int n, double *A,
                       int lda, double *W, double *work, int lwork,
                       int *devinfo) {
  return (int)hipsolverDnDsyevd((hipsolverHandle_t)handle,
                                (hipsolverEigMode_t)jobz,
                                (hipblasFillMode_t)uplo, n, A, lda, W, work,
                                lwork, devinfo);
}

int geps_solver_dsyevj_buf(void *handle, int jobz, int uplo, int n, double *A,
                           int lda, double *W, int *lwork, void *params) {
  return (int)hipsolverDnDsyevj_bufferSize(
      (hipsolverDnHandle_t)handle, (hipsolverEigMode_t)jobz,
      (hipblasFillMode_t)uplo, n, A, lda, W, lwork,
      (hipsolverSyevjInfo_t)params);
}

int geps_solver_dsyevj(void *handle, int jobz, int uplo, int n, double *A,
                       int lda, double *W, double *work, int lwork,
                       int *devinfo, void *params) {
  /* #region agent log */
  static int syevj_n = 0;
  int seq = syevj_n++;
  if (seq < 2)
    dbg_ndjson("C", "hip_compat:dsyevj", "before", seq, (long long)(uintptr_t)A,
               (long long)n, (long long)lda, 0);
  /* #endregion */
  int rc = (int)hipsolverDnDsyevj(
      (hipsolverDnHandle_t)handle, (hipsolverEigMode_t)jobz,
      (hipblasFillMode_t)uplo, n, A, lda, W, work, lwork, devinfo,
      (hipsolverSyevjInfo_t)params);
  /* #region agent log */
  if (seq < 2)
    dbg_ndjson("C", "hip_compat:dsyevj", "after", seq, (long long)rc,
               (long long)n, (long long)lda, 1);
  /* #endregion */
  return rc;
}

int geps_rocm_dsyevj_dev(double *A, int n, double *W, int *devinfo) {
  if (n <= 0)
    return 0;
  /* #region agent log: is `use_device_addr` actually giving us a device
     address? device_ptr() silently passes a host pointer straight through when
     omp_target_is_present() is false, which would make hipSOLVER write
     somewhere other than the `mx` the model reads back. Print the first few. */
  {
    static int seen = 0;
    if (seen < 4) {
      ++seen;
      int d = omp_get_default_device();
      int present = omp_target_is_present(A, d);
      void *mapped = present ? omp_get_mapped_ptr(A, d) : nullptr;
      void *resolved = device_ptr(A);
      fprintf(stderr,
              "[dsyevj] n=%d A=%p is_present=%d mapped=%p resolved=%p %s\n", n,
              (void *)A, present, mapped, resolved,
              (resolved == (void *)A && !present)
                  ? "<-- passed through unchanged (suspect host address)"
                  : "");
      fflush(stderr);
    }
  }
  /* #endregion */
  A = (double *)device_ptr(A);
  W = (double *)device_ptr(W);
  devinfo = (int *)device_ptr(devinfo);

  hipsolverHandle_t h = nullptr;
  hipsolverSyevjInfo_t params = nullptr;
  {
    std::lock_guard<std::mutex> lock(g_mu);
    if (!g_solver)
      (void)hipsolverDnCreate(&g_solver);
    if (!g_syevj) {
      (void)hipsolverDnCreateSyevjInfo(&g_syevj);
      (void)hipsolverDnXsyevjSetTolerance(g_syevj, 1e-15);
      (void)hipsolverDnXsyevjSetSortEig(g_syevj, 0);
    }
    h = g_solver;
    params = g_syevj;
  }

  int lwork = 0;
  int rc = (int)hipsolverDnDsyevj_bufferSize(
      h, HIPSOLVER_EIG_MODE_VECTOR, HIPBLAS_FILL_MODE_UPPER, n, A, n, W, &lwork,
      params);
  if (rc != 0)
    return rc;
  if (lwork < 1)
    lwork = 1;
  double *work = nullptr;
  if (hipMalloc(&work, sizeof(double) * (size_t)lwork) != hipSuccess)
    return -1;
  rc = (int)hipsolverDnDsyevj(h, HIPSOLVER_EIG_MODE_VECTOR,
                              HIPBLAS_FILL_MODE_UPPER, n, A, n, W, work, lwork,
                              devinfo, params);
  (void)hipDeviceSynchronize();
  (void)hipFree(work);
  return rc;
}

int geps_sparse_create(void **handle) {
  hipsparseHandle_t h = nullptr;
  int rc = (int)hipsparseCreate(&h);
  *handle = (void *)h;
  if (!g_sparse)
    g_sparse = h;
  return rc;
}

int geps_sparse_set_stream(void *handle, int64_t stream) {
  return (int)hipsparseSetStream((hipsparseHandle_t)handle,
                                 use_stream(stream));
}

int geps_sparse_dgtsv_interleaved_buf(void *handle, int algo, int m,
                                      const double *dl, const double *d,
                                      const double *du, const double *x,
                                      int batch, int64_t *bytes) {
  size_t sz = 0;
  int rc = (int)hipsparseDgtsvInterleavedBatch_bufferSizeExt(
      (hipsparseHandle_t)handle, algo, m, dl, d, du, x, batch, &sz);
  if (bytes)
    *bytes = (int64_t)sz;
  return rc;
}

// 2026-09-16: rocSPARSE's Thomas solver (algo 0) overwrites BOTH the main
// diagonal `d` and the upper diagonal `du` (measured with a standalone test:
// d changed 72000/72000, du 71000/72000, dl untouched), whereas the cuSPARSE
// callers in moninedmf_gpu only restore `aug` (du) between solves that share
// alg/adg/aug.  NOTE: those call sites are inside `if (.false.)` today (the
// live path is tridin_gpu/tridi2_gpu), so this is NOT the cause of the
// layer-18 microphysics hang -- it only matters if that path is re-enabled.
// Preserve all three coefficient arrays so both libraries look alike.
static double *g_gtsv_save = nullptr;
static size_t g_gtsv_cap = 0;

int geps_sparse_dgtsv_interleaved(void *handle, int algo, int m, double *dl,
                                  double *d, double *du, double *x, int batch,
                                  void *buf) {
  hipsparseHandle_t h = (hipsparseHandle_t)handle;
  size_t n = (size_t)m * (size_t)batch;
  if (n > g_gtsv_cap) {
    if (g_gtsv_save)
      (void)hipFree(g_gtsv_save);
    g_gtsv_save = nullptr;
    if (hipMalloc(&g_gtsv_save, 3 * n * sizeof(double)) != hipSuccess) {
      g_gtsv_cap = 0;
      return (int)HIPSPARSE_STATUS_ALLOC_FAILED;
    }
    g_gtsv_cap = n;
  }
  hipStream_t s = nullptr;
  (void)hipsparseGetStream(h, &s);
  double *sdl = g_gtsv_save, *sd = g_gtsv_save + n, *sdu = g_gtsv_save + 2 * n;
  (void)hipMemcpyAsync(sdl, dl, n * sizeof(double), hipMemcpyDeviceToDevice, s);
  (void)hipMemcpyAsync(sd, d, n * sizeof(double), hipMemcpyDeviceToDevice, s);
  (void)hipMemcpyAsync(sdu, du, n * sizeof(double), hipMemcpyDeviceToDevice, s);
  int rc = (int)hipsparseDgtsvInterleavedBatch(h, algo, m, dl, d, du, x, batch, buf);
  (void)hipMemcpyAsync(dl, sdl, n * sizeof(double), hipMemcpyDeviceToDevice, s);
  (void)hipMemcpyAsync(d, sd, n * sizeof(double), hipMemcpyDeviceToDevice, s);
  (void)hipMemcpyAsync(du, sdu, n * sizeof(double), hipMemcpyDeviceToDevice, s);
  (void)hipStreamSynchronize(s);
  return rc;
}

int geps_nccl_unique_id(void *id) {
  return (int)ncclGetUniqueId((ncclUniqueId *)id);
}

int geps_nccl_comm_init_rank(void **comm, int nranks, void *id, int rank) {
  ncclComm_t c = nullptr;
  int rc = (int)ncclCommInitRank(&c, nranks, *(ncclUniqueId *)id, rank);
  *comm = (void *)c;
  return rc;
}

int geps_nccl_comm_split(void *comm, int color, int key, void **newcomm) {
  ncclComm_t c = nullptr;
  int rc = (int)ncclCommSplit((ncclComm_t)comm, color, key, &c, nullptr);
  *newcomm = (void *)c;
  return rc;
}

int geps_nccl_comm_destroy(void *comm) {
  return (int)ncclCommDestroy((ncclComm_t)comm);
}

int geps_nccl_group_start() { return (int)ncclGroupStart(); }
int geps_nccl_group_end() { return (int)ncclGroupEnd(); }

int geps_nccl_send(const void *s, int64_t count, int dtype, int peer,
                   void *comm, int64_t stream) {
  s = device_ptr(const_cast<void *>(s));
  return (int)ncclSend(s, (size_t)count, (ncclDataType_t)dtype, peer,
                       (ncclComm_t)comm, use_stream(stream));
}

int geps_nccl_recv(void *r, int64_t count, int dtype, int peer, void *comm,
                   int64_t stream) {
  r = device_ptr(r);
  return (int)ncclRecv(r, (size_t)count, (ncclDataType_t)dtype, peer,
                       (ncclComm_t)comm, use_stream(stream));
}

int geps_nccl_allreduce(const void *s, void *r, int64_t count, int dtype,
                        int op, void *comm, int64_t stream) {
  s = device_ptr(const_cast<void *>(s));
  r = device_ptr(r);
  return (int)ncclAllReduce(s, r, (size_t)count, (ncclDataType_t)dtype,
                            (ncclRedOp_t)op, (ncclComm_t)comm,
                            use_stream(stream));
}

static size_t nccl_elem_size(int dtype) {
  switch ((ncclDataType_t)dtype) {
  case ncclInt8: case ncclUint8: return 1;
  case ncclFloat16: case ncclBfloat16: return 2;
  case ncclInt32: case ncclUint32: case ncclFloat32: return 4;
  default: return 8;
  }
}

int geps_nccl_allgather(const void *s, void *r, int64_t count, int dtype,
                        void *comm, int64_t stream) {
  s = device_ptr(const_cast<void *>(s));
  r = device_ptr(r);
  // 2026-09-18: out-of-place AllGather makes RCCL copy the local chunk with
  // hipMemcpyAsync(hipMemcpyDefault); HIP does not know libomptarget's
  // buffers and takes the CPU/large-BAR path (~26 MB/s: 3-10 s per call,
  // ~30 s per step, seen with an LD_PRELOAD backtrace from
  // mpe2d_unify_lev_gpu). In-place (sendbuff == recvbuff + rank*count)
  // skips that copy; stage the chunk with libomptarget's own D2D copy.
  int rank = 0;
  if (ncclCommUserRank((ncclComm_t)comm, &rank) == ncclSuccess) {
    size_t bytes = (size_t)count * nccl_elem_size(dtype);
    char *slot = (char *)r + (size_t)rank * bytes;
    if ((const void *)slot != s) {
      int dev = omp_get_default_device();
      (void)omp_target_memcpy(slot, s, bytes, 0, 0, dev, dev);
      s = slot;
    }
  }
  return (int)ncclAllGather(s, r, (size_t)count, (ncclDataType_t)dtype,
                            (ncclComm_t)comm, use_stream(stream));
}

int geps_nccl_broadcast(const void *s, void *r, int64_t count, int dtype,
                        int root, void *comm, int64_t stream) {
  /* #region agent log */
  static int bcast_n = 0;
  int seq = bcast_n++;
  int interesting = (seq < 3) || (count >= 100000);
  if (interesting)
    dbg_ndjson("B", "hip_compat:ncclBcast", "before", seq,
               (long long)(uintptr_t)s, (long long)(uintptr_t)r, count, 0);
  /* #endregion */
  s = device_ptr(const_cast<void *>(s));
  r = device_ptr(r);
  // same in-place trick as AllGather: the root's local copy is a slow
  // hipMemcpyAsync otherwise; in-place broadcast (s == r on every rank)
  // needs no copy at all.
  if (s != r) {
    int rank = -1;
    if (ncclCommUserRank((ncclComm_t)comm, &rank) == ncclSuccess) {
      if (rank == root) {
        int dev = omp_get_default_device();
        (void)omp_target_memcpy(r, s, (size_t)count * nccl_elem_size(dtype), 0, 0, dev, dev);
      }
      s = r;
    }
  }
  int rc = (int)ncclBroadcast(s, r, (size_t)count, (ncclDataType_t)dtype, root,
                            (ncclComm_t)comm, use_stream(stream));
  /* #region agent log */
  if (interesting)
    dbg_ndjson("B", "hip_compat:ncclBcast", "after", seq, (long long)rc,
               (long long)root, count, 1);
  /* #endregion */
  return rc;
}

void geps_nccl_error_string(int err, char *buf, int n) {
  const char *s = ncclGetErrorString((ncclResult_t)err);
  if (!s)
    s = "unknown RCCL error";
  std::snprintf(buf, n > 0 ? (size_t)n : 0, "%s", s);
}

void acc_map_data_(void *host, void *device, int *nbytes)
{
  int dev = omp_get_default_device();
  size_t n = nbytes ? (size_t)(*nbytes) : 0;
  (void)omp_target_associate_ptr(host, device, n, 0, dev);
}

void acc_unmap_data_(void *host)
{
  (void)omp_target_disassociate_ptr(host, omp_get_default_device());
}

void acc_map_data(void *host, void *device, int *nbytes)
{
  acc_map_data_(host, device, nbytes);
}

void acc_unmap_data(void *host) { acc_unmap_data_(host); }

} /* extern "C" */
