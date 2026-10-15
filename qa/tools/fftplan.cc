#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>
#include <cstdio>
#include <vector>
static double freeGB(){ size_t f=0,t=0; hipMemGetInfo(&f,&t); return f/1073741824.0; }
int main(int argc,char**argv){
  int nplans = argc>1? atoi(argv[1]) : 768;
  int batch  = argc>2? atoi(argv[2]) : 72;
  int n = 1552, jump = n+2;
  double f0 = freeGB();
  std::vector<hipfftHandle> plans;
  for(int i=0;i<nplans;i++){
    int nxj = 1552 - (i % 64) * 8;          // reduced grid: varying n per latitude
    hipfftHandle p; hipfftCreate(&p);
    int nn = nxj, ine = nxj, one = nxj; size_t ws=0;
    hipfftResult r = hipfftMakePlanMany(p, 1, &nn, &ine, 1, jump, &one, 1, jump/2, HIPFFT_D2Z, batch, &ws);
    if(r!=HIPFFT_SUCCESS){ printf("plan %d failed %d\n", i, r); return 1; }
    plans.push_back(p);
    if(i<3 || (i+1)%128==0) printf("plan %4d n=%d ws=%.1f MB  free=%.2f GB  used=%.2f GB\n", i+1, nxj, ws/1048576.0, freeGB(), f0-freeGB());
  }
  printf("TOTAL %d plans batch=%d: %.2f GB\n", nplans, batch, f0-freeGB());
  for(auto p:plans) hipfftDestroy(p);
  printf("after destroy: used=%.2f GB\n", f0-freeGB());
  return 0;
}
