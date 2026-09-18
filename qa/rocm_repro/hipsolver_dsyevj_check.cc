// Same call sequence as geps_rocm_dsyevj_dev (src/rocm/hip_compat.cc:488),
// inlined so we do not have to link the whole compat layer.
#include <hip/hip_runtime.h>
#include <hipsolver/hipsolver.h>
#include <cstdio>
#include <cmath>

static hipsolverHandle_t g_solver = nullptr;
static hipsolverSyevjInfo_t g_syevj = nullptr;

int dsyevj_dev(double *A, int n, double *W, int *devinfo) {
  if (n <= 0) return 0;
  if (!g_solver) hipsolverDnCreate(&g_solver);
  if (!g_syevj) {
    hipsolverDnCreateSyevjInfo(&g_syevj);
    hipsolverDnXsyevjSetTolerance(g_syevj, 1e-15);
    hipsolverDnXsyevjSetSortEig(g_syevj, 0);
  }
  int lwork = 0;
  int rc = (int)hipsolverDnDsyevj_bufferSize(g_solver, HIPSOLVER_EIG_MODE_VECTOR,
             HIPSOLVER_FILL_MODE_UPPER, n, A, n, W, &lwork, g_syevj);
  printf("  bufferSize rc=%d lwork=%d\n", rc, lwork);
  if (rc != 0) return rc;
  if (lwork < 1) lwork = 1;
  double *work = nullptr;
  if (hipMalloc(&work, sizeof(double)*(size_t)lwork) != hipSuccess) return -1;
  rc = (int)hipsolverDnDsyevj(g_solver, HIPSOLVER_EIG_MODE_VECTOR,
         HIPSOLVER_FILL_MODE_UPPER, n, A, n, W, work, lwork, devinfo, g_syevj);
  hipDeviceSynchronize();
  hipFree(work);
  return rc;
}

int main() {
  const int n = 4;
  double h[16] = { 4,1,0,0,  1,3,1,0,  0,1,2,1,  0,0,1,1 };   // symmetric
  double *dA,*dW; int *dI;
  hipMalloc(&dA,sizeof(h)); hipMalloc(&dW,n*sizeof(double)); hipMalloc(&dI,sizeof(int));
  hipMemcpy(dA,h,sizeof(h),hipMemcpyHostToDevice);
  hipMemset(dI,0,sizeof(int));
  int rc = dsyevj_dev(dA,n,dW,dI);
  double V[16],W[4]; int info=-999;
  hipMemcpy(V,dA,sizeof(V),hipMemcpyDeviceToHost);
  hipMemcpy(W,dW,sizeof(W),hipMemcpyDeviceToHost);
  hipMemcpy(&info,dI,sizeof(int),hipMemcpyDeviceToHost);
  printf("rc=%d  devinfo=%d\n", rc, info);
  printf("eigenvalues: %.6f %.6f %.6f %.6f\n", W[0],W[1],W[2],W[3]);
  double amax=0; for(double v:V) amax=fmax(amax,fabs(v));
  printf("max|eigenvector entry| = %.6g   (must be <= 1)\n", amax);
  for(int c=0;c<n;c++){ double s=0; for(int r=0;r<n;r++) s+=V[c*n+r]*V[c*n+r];
    printf("  column %d norm^2 = %.6f  (must be 1)\n", c, s); }
  return 0;
}
