#define _GNU_SOURCE
#include <dlfcn.h>
#include <execinfo.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <hip/hip_runtime.h>
typedef hipError_t (*fn_t)(void*, const void*, size_t, hipMemcpyKind, hipStream_t);
extern "C" hipError_t hipMemcpyAsync(void* dst, const void* src, size_t n, hipMemcpyKind kind, hipStream_t s) {
  static fn_t real = (fn_t)dlsym(RTLD_NEXT, "hipMemcpyAsync");
  static std::set<std::string> seen;
  if (n >= (size_t)16 << 20) {
    void* bt[24]; int d = backtrace(bt, 24);
    char** sy = backtrace_symbols(bt, d);
    std::string key; for (int i = 1; i < d && i < 6; i++) { key += sy[i]; key += "|"; }
    if (!seen.count(key)) { seen.insert(key);
      fprintf(stderr, "[memcpy_bt] hipMemcpyAsync %zu MB kind=%d\n", n >> 20, (int)kind);
      for (int i = 1; i < d; i++) fprintf(stderr, "   %s\n", sy[i]); }
    free(sy);
  }
  return real(dst, src, n, kind, s);
}
