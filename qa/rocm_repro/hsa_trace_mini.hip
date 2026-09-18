#include <hip/hip_runtime.h>
#include <cstdio>
__global__ void k(float* p, int n){ int i = blockIdx.x*blockDim.x+threadIdx.x; if(i<n) p[i] = p[i]*2.0f + 1.0f; }
int main(){
  const int n = 1<<20; float *d=nullptr, *h=(float*)malloc(n*sizeof(float));
  for(int i=0;i<n;i++) h[i]=1.0f;
  hipMalloc(&d, n*sizeof(float));
  hipMemcpy(d,h,n*sizeof(float),hipMemcpyHostToDevice);
  for(int r=0;r<8;r++) hipLaunchKernelGGL(k, dim3((n+255)/256), dim3(256), 0, 0, d, n);
  hipMemcpy(h,d,n*sizeof(float),hipMemcpyDeviceToHost);
  hipDeviceSynchronize();
  printf("mini ok h[0]=%f\n", h[0]);
  hipFree(d); free(h); return 0;
}
