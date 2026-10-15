#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>
#include <cstdio>
#include <vector>
#include <chrono>
int main(){
  for(int rep=0;rep<2;rep++){
    auto t0=std::chrono::steady_clock::now(); std::vector<hipfftHandle> ps;
    for(int i=0;i<384;i++){ int n=1552-i*4, jump=n+2; hipfftHandle p; hipfftCreate(&p); size_t ws; int nn=n,ine=n,one=n;
      hipfftMakePlanMany(p,1,&nn,&ine,1,jump,&one,1,jump/2,HIPFFT_D2Z,72,&ws); ps.push_back(p);}
    auto t1=std::chrono::steady_clock::now();
    for(auto p:ps) hipfftDestroy(p);
    auto t2=std::chrono::steady_clock::now();
    printf("rep %d: create 384 plans %.2f s, destroy %.2f s\n", rep, std::chrono::duration<double>(t1-t0).count(), std::chrono::duration<double>(t2-t1).count());
  }
}
