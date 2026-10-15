#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>
#include <cstdio>
#include <vector>
#include <cstdlib>
static double freeGB(){ size_t f=0,t=0; hipMemGetInfo(&f,&t); return f/1073741824.0; }
double run(const char* tag, int nplans, int batch, int nfixed, int step, int type){
  double f0=freeGB(); std::vector<hipfftHandle> ps;
  for(int i=0;i<nplans;i++){
    int n = nfixed>0 ? nfixed : 1552 - (i%384)*4; int jump=n+2;
    hipfftHandle p; hipfftCreate(&p); size_t ws=0; int nn=n, ine=n, one=n;
    hipfftResult r = type==0 ? hipfftMakePlanMany(p,1,&nn,&ine,1,jump,&one,1,jump/2,HIPFFT_D2Z,batch,&ws)
                             : hipfftMakePlanMany(p,1,&nn,&ine,1,jump/2,&one,1,jump,HIPFFT_Z2D,batch,&ws);
    if(r) { printf("fail %d\n", r); break; } ps.push_back(p);
  }
  double used=f0-freeGB(); printf("%-40s %4d plans batch=%3d : %.2f GB  (%.1f MB/plan)\n", tag, nplans, batch, used, used*1024/nplans);
  for(auto p:ps) hipfftDestroy(p); return used;
}
int main(){
  run("same n=1552, D2Z", 384, 72, 1552, 0, 0);
  run("same n=1024 (pow2), D2Z", 384, 72, 1024, 0, 0);
  run("distinct n (1552-4i), D2Z", 384, 72, 0, 4, 0);
  run("distinct n, batch=1, D2Z", 384, 1, 0, 4, 0);
  run("distinct n, Z2D", 384, 72, 0, 4, 1);
  return 0;
}
