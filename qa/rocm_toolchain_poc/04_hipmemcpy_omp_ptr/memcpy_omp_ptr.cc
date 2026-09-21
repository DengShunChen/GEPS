// PoC 04 - hipMemcpy* on device memory that was allocated by libomptarget
// (omp_target_alloc, or any `target enter data`) runs at tens of MB/s.
//
// hipPointerGetAttributes *does* classify the pointer as device memory
// (type=2), yet hipMemcpy / hipMemcpyDtoD / hipMemcpyAsync all fall back to a
// CPU copy through the large-BAR mapping (~50 MB/s) instead of an SDMA or
// blit-kernel copy (~1.4 TB/s).  The bytes are still correct, so nothing
// reports an error - the program is just ~30000x slower than expected.
//
// Secondary observation: omp_target_memcpy between the same buffers reaches
// only ~57 GB/s, 25x below hipMemcpy on hipMalloc'ed buffers.
//
// Interaction: the same thing happens inside any HIP library handed an
// OpenMP-mapped buffer (rocFFT, RCCL out-of-place collectives, ...).
#include <hip/hip_runtime.h>
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define CHK(x) do { hipError_t e = (x); if (e != hipSuccess) { \
   fprintf(stderr, "%s -> %s\n", #x, hipGetErrorString(e)); exit(1); } } while (0)

static double now() {
   return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static void describe(const char *label, const void *p) {
   hipPointerAttribute_t at;
   memset(&at, 0, sizeof at);
   hipError_t e = hipPointerGetAttributes(&at, p);
   if (e == hipSuccess)
      printf("  %-22s %p  hipPointerGetAttributes: type=%d device=%d\n", label, p, (int)at.type, at.device);
   else
      printf("  %-22s %p  hipPointerGetAttributes: %s\n", label, p, hipGetErrorString(e));
}

typedef hipError_t (*copy_fn)(void *, const void *, size_t);

static hipError_t cp_default(void *d, const void *s, size_t n) { return hipMemcpy(d, s, n, hipMemcpyDefault); }
static hipError_t cp_d2d(void *d, const void *s, size_t n)     { return hipMemcpy(d, s, n, hipMemcpyDeviceToDevice); }
static hipError_t cp_dtod(void *d, const void *s, size_t n)    { return hipMemcpyDtoD((hipDeviceptr_t)d, (hipDeviceptr_t)s, n); }
static hipError_t cp_async(void *d, const void *s, size_t n) {
   hipError_t e = hipMemcpyAsync(d, s, n, hipMemcpyDeviceToDevice, 0);
   if (e == hipSuccess) e = hipStreamSynchronize(0);
   return e;
}
static hipError_t cp_omp(void *d, const void *s, size_t n) {
   int dev = omp_get_default_device();
   return omp_target_memcpy(d, s, n, 0, 0, dev, dev) == 0 ? hipSuccess : hipErrorUnknown;
}

static double bench(copy_fn f, void *dst, const void *src, size_t n, int reps) {
   CHK(f(dst, src, n));                       // warm-up
   double t0 = now();
   for (int r = 0; r < reps; r++) CHK(f(dst, src, n));
   CHK(hipDeviceSynchronize());
   double dt = (now() - t0) / reps;
   return n / dt / 1e9;                       // GB/s
}

int main(int argc, char **argv) {
   size_t mb = argc > 1 ? atoi(argv[1]) : 64;
   size_t n = mb << 20;
   int dev = omp_get_default_device();

   void *o_src = omp_target_alloc(n, dev), *o_dst = omp_target_alloc(n, dev);
   void *h_src, *h_dst;
   CHK(hipMalloc(&h_src, n)); CHK(hipMalloc(&h_dst, n));
   if (!o_src || !o_dst) { fprintf(stderr, "omp_target_alloc failed\n"); return 1; }

   printf("buffer %zu MB, device %d\n", mb, dev);
   describe("hipMalloc", h_src);
   describe("omp_target_alloc", o_src);

   // fill via OpenMP so the source content is defined on the device
   unsigned char *s8 = (unsigned char *)o_src;
   #pragma omp target teams distribute parallel for is_device_ptr(s8)
   for (size_t i = 0; i < n; i++) s8[i] = (unsigned char)i;
   CHK(hipMemset(h_src, 0x5a, n));

   struct { const char *name; copy_fn f; } fns[] = {
      { "hipMemcpy(Default)",          cp_default },
      { "hipMemcpy(DeviceToDevice)",   cp_d2d },
      { "hipMemcpyDtoD",               cp_dtod },
      { "hipMemcpyAsync(D2D)+sync",    cp_async },
      { "omp_target_memcpy",           cp_omp },
   };
   printf("\n  %-28s %14s %14s   [GB/s]\n", "copy API", "hipMalloc bufs", "omp bufs");
   int reps_small = 2;
   for (auto &fn : fns) {
      double g_hip = bench(fn.f, h_dst, h_src, n, 10);
      double g_omp = bench(fn.f, o_dst, o_src, n, reps_small);
      printf("  %-28s %14.1f %14.3f%s\n", fn.name, g_hip, g_omp, g_omp < g_hip / 20 ? "   <-- CPU path" : "");
   }

   // correctness check of the slow path: bytes are right, only the speed is wrong
   CHK(hipMemcpy(o_dst, o_src, n, hipMemcpyDeviceToDevice));
   unsigned char *d8 = (unsigned char *)o_dst;
   long bad = 0;
   #pragma omp target teams distribute parallel for is_device_ptr(d8) reduction(+:bad)
   for (size_t i = 0; i < n; i++) bad += d8[i] != (unsigned char)i;
   printf("\n  hipMemcpy on omp buffers: %ld wrong bytes of %zu (data is correct, only slow)\n", bad, n);
   return 0;
}
