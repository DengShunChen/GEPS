// rocFFT/hipFFT exec host-overhead microbenchmark: mirrors rfftmlt's per-latitude
// D2Z plans (batch = m levels, distinct lengths), 1 vs N host threads.
#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>
#include <chrono>
#include <cstdio>
#include <thread>
#include <vector>
#include <cstdlib>
int main(int argc, char **argv) {
  int nlat = argc > 1 ? atoi(argv[1]) : 192;
  int m = argc > 2 ? atoi(argv[2]) : 72;
  int nthr = argc > 3 ? atoi(argv[3]) : 1;
  int reps = argc > 4 ? atoi(argv[4]) : 10;
  // reduced-grid-like lengths: 1024..1536 step 4ish, distinct
  std::vector<int> len(nlat);
  for (int j = 0; j < nlat; ++j) len[j] = 64 + (j * 1500) / nlat;
  int jump = 1536 + 2;
  size_t words = (size_t)jump * m;
  double *cc; hipfftDoubleComplex *gw;
  hipMalloc(&cc, words * nlat * sizeof(double));
  hipMalloc(&gw, words * nlat * sizeof(double));
  hipMemset(cc, 0, words * nlat * sizeof(double));
  std::vector<hipfftHandle> plans(nlat);
  for (int j = 0; j < nlat; ++j) {
    hipfftCreate(&plans[j]);
    int n = len[j], ine = n, one = n;
    size_t ws;
    hipfftMakePlanMany(plans[j], 1, &n, &ine, 1, jump, &one, 1, jump / 2, HIPFFT_D2Z, m, &ws);
  }
  std::vector<hipStream_t> st(nthr);
  for (auto &s : st) hipStreamCreateWithFlags(&s, hipStreamNonBlocking);
  auto run = [&](int t) {
    for (int j = t; j < nlat; j += nthr) {
      hipfftSetStream(plans[j], st[t]);
      hipfftExecD2Z(plans[j], cc + words * j, gw + (words / 2) * j);
    }
  };
  // warm-up
  for (int t = 0; t < nthr; ++t) run(t);
  hipDeviceSynchronize();
  auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    std::vector<std::thread> th;
    for (int t = 0; t < nthr; ++t) th.emplace_back(run, t);
    for (auto &x : th) x.join();
    for (auto &s : st) hipStreamSynchronize(s);
  }
  auto t1 = std::chrono::steady_clock::now();
  double us = std::chrono::duration<double, std::micro>(t1 - t0).count() / reps;
  printf("nlat=%d m=%d threads=%d: %.0f us per loop, %.1f us per exec\n", nlat, m, nthr, us, us / nlat);
  // launch-only cost: same loop but no sync until the end, measure issue time
  t0 = std::chrono::steady_clock::now();
  run(0);
  t1 = std::chrono::steady_clock::now();
  us = std::chrono::duration<double, std::micro>(t1 - t0).count();
  hipDeviceSynchronize();
  printf("single-thread issue only: %.1f us per exec\n", us / ((nlat + nthr - 1) / nthr));
  return 0;
}
