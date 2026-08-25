/* HIP runtime helpers for GEPS AMD GPU path (CUDA Fortran API subset). */

#include <cstdint>
#include <cstdio>
#include <cstring>
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
  return s;
}

} /* namespace */

extern "C" {

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

int geps_hip_device_count() {
  int n = 0;
  (void)hipGetDeviceCount(&n);
  return n;
}

int geps_hip_set_device(int id) { return (int)hipSetDevice(id); }

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

int geps_hip_memcpy_async(void *dst, const void *src, int64_t bytes, int /*kind*/,
                          int64_t stream) {
  return (int)hipMemcpyAsync(dst, src, (size_t)bytes, hipMemcpyDefault,
                             (hipStream_t)(uintptr_t)stream);
}

int geps_hip_memset_async(void *dst, int64_t bytes, int64_t stream) {
  return (int)hipMemsetAsync(dst, 0, (size_t)bytes, (hipStream_t)(uintptr_t)stream);
}

int geps_hip_memset_i32_async(void *dst, int val, int64_t count, int64_t stream) {
  return (int)hipMemsetD32Async((hipDeviceptr_t)dst, val, (size_t)count,
                                (hipStream_t)(uintptr_t)stream);
}

int geps_hip_malloc_async(void **ptr, int64_t bytes, int64_t stream) {
  return (int)hipMallocAsync(ptr, (size_t)bytes, (hipStream_t)(uintptr_t)stream);
}

int geps_hip_free_async(void *ptr, int64_t stream) {
  return (int)hipFreeAsync(ptr, (hipStream_t)(uintptr_t)stream);
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
  return (int)hipEventRecord((hipEvent_t)ev, (hipStream_t)(uintptr_t)stream);
}

int geps_hip_stream_wait_event(int64_t stream, void *ev, unsigned int flags) {
  return (int)hipStreamWaitEvent((hipStream_t)(uintptr_t)stream, (hipEvent_t)ev,
                                 flags);
}

int geps_hip_stream_create_flags(int64_t *stream, unsigned int flags) {
  hipStream_t s = nullptr;
  int rc = (int)hipStreamCreateWithFlags(&s, flags);
  *stream = (int64_t)(uintptr_t)s;
  return rc;
}

int geps_hip_stream_sync(int64_t stream) {
  return (int)hipStreamSynchronize((hipStream_t)(uintptr_t)stream);
}

int geps_hip_begin_capture(int64_t stream, int mode) {
  return (int)hipStreamBeginCapture((hipStream_t)(uintptr_t)stream,
                                    (hipStreamCaptureMode)mode);
}

int geps_hip_end_capture(int64_t stream, void **graph) {
  hipGraph_t g = nullptr;
  int rc = (int)hipStreamEndCapture((hipStream_t)(uintptr_t)stream, &g);
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
                             (hipStream_t)(uintptr_t)stream);
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
                               (hipStream_t)(uintptr_t)stream);
}

int geps_blas_dgemm(void *handle, int ta, int tb, int m, int n, int k,
                    double alpha, const double *a, int lda, const double *b,
                    int ldb, double beta, double *c, int ldc) {
  return (int)hipblasDgemm((hipblasHandle_t)handle, (hipblasOperation_t)ta,
                           (hipblasOperation_t)tb, m, n, k, &alpha, a, lda, b,
                           ldb, &beta, c, ldc);
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
  int rc = (int)hipfftMakePlanMany(fft_of(plan), rank, &nn, &ine, istride,
                                   idist, &one, ostride, odist,
                                   (hipfftType)ffttype, batch, &ws);
  if (work_size)
    *work_size = (int64_t)ws;
  return rc;
}

int geps_fft_set_stream(int plan, int64_t stream) {
  return (int)hipfftSetStream(fft_of(plan), (hipStream_t)(uintptr_t)stream);
}

int geps_fft_set_work_area(int plan, void *work) {
  return (int)hipfftSetWorkArea(fft_of(plan), work);
}

int geps_fft_exec_z2d(int plan, void *in, void *out) {
  return (int)hipfftExecZ2D(fft_of(plan), (hipfftDoubleComplex *)in,
                            (double *)out);
}
int geps_fft_exec_d2z(int plan, void *in, void *out) {
  return (int)hipfftExecD2Z(fft_of(plan), (double *)in,
                            (hipfftDoubleComplex *)out);
}
int geps_fft_exec_c2r(int plan, void *in, void *out) {
  return (int)hipfftExecC2R(fft_of(plan), (hipfftComplex *)in, (float *)out);
}
int geps_fft_exec_r2c(int plan, void *in, void *out) {
  return (int)hipfftExecR2C(fft_of(plan), (float *)in, (hipfftComplex *)out);
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
                                   (hipStream_t)(uintptr_t)stream);
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
  return (int)hipsolverDnDsyevj(
      (hipsolverDnHandle_t)handle, (hipsolverEigMode_t)jobz,
      (hipblasFillMode_t)uplo, n, A, lda, W, work, lwork, devinfo,
      (hipsolverSyevjInfo_t)params);
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
                                 (hipStream_t)(uintptr_t)stream);
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

int geps_sparse_dgtsv_interleaved(void *handle, int algo, int m, double *dl,
                                  double *d, double *du, double *x, int batch,
                                  void *buf) {
  return (int)hipsparseDgtsvInterleavedBatch((hipsparseHandle_t)handle, algo, m,
                                             dl, d, du, x, batch, buf);
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
  return (int)ncclSend(s, (size_t)count, (ncclDataType_t)dtype, peer,
                       (ncclComm_t)comm, (hipStream_t)(uintptr_t)stream);
}

int geps_nccl_recv(void *r, int64_t count, int dtype, int peer, void *comm,
                   int64_t stream) {
  return (int)ncclRecv(r, (size_t)count, (ncclDataType_t)dtype, peer,
                       (ncclComm_t)comm, (hipStream_t)(uintptr_t)stream);
}

int geps_nccl_allreduce(const void *s, void *r, int64_t count, int dtype,
                        int op, void *comm, int64_t stream) {
  return (int)ncclAllReduce(s, r, (size_t)count, (ncclDataType_t)dtype,
                            (ncclRedOp_t)op, (ncclComm_t)comm,
                            (hipStream_t)(uintptr_t)stream);
}

int geps_nccl_allgather(const void *s, void *r, int64_t count, int dtype,
                        void *comm, int64_t stream) {
  return (int)ncclAllGather(s, r, (size_t)count, (ncclDataType_t)dtype,
                            (ncclComm_t)comm, (hipStream_t)(uintptr_t)stream);
}

int geps_nccl_broadcast(const void *s, void *r, int64_t count, int dtype,
                        int root, void *comm, int64_t stream) {
  return (int)ncclBroadcast(s, r, (size_t)count, (ncclDataType_t)dtype, root,
                            (ncclComm_t)comm, (hipStream_t)(uintptr_t)stream);
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
