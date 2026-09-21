// Free VRAM as seen by HIP, callable from Fortran.
#include <hip/hip_runtime.h>
extern "C" double freegb(void) { size_t f = 0, t = 0; hipMemGetInfo(&f, &t); return f / 1073741824.0; }
