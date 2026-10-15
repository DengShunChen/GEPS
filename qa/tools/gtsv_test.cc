#include <hip/hip_runtime.h>
#include <hipsparse/hipsparse.h>
#include <cstdio>
#include <vector>
#include <cmath>
int main(){
  const int m=72, batch=1000, n=m*batch;
  std::vector<double> dl(n),d(n),du(n),x(n);
  for(int k=0;k<m;k++) for(int b=0;b<batch;b++){ int i=k*batch+b;
    dl[i]=(k==0)?0.0:-0.3-0.001*b; d[i]=1.0+0.01*k; du[i]=(k==m-1)?0.0:-0.2-0.002*k; x[i]=1.0+k+b*0.5; }
  auto h0=dl, h1=d, h2=du;
  double *ddl,*dd,*ddu,*dx; hipMalloc(&ddl,n*8); hipMalloc(&dd,n*8); hipMalloc(&ddu,n*8); hipMalloc(&dx,n*8);
  hipMemcpy(ddl,dl.data(),n*8,hipMemcpyHostToDevice); hipMemcpy(dd,d.data(),n*8,hipMemcpyHostToDevice);
  hipMemcpy(ddu,du.data(),n*8,hipMemcpyHostToDevice); hipMemcpy(dx,x.data(),n*8,hipMemcpyHostToDevice);
  hipsparseHandle_t h; hipsparseCreate(&h);
  for(int algo=0;algo<3;algo++){
    hipMemcpy(ddl,h0.data(),n*8,hipMemcpyHostToDevice); hipMemcpy(dd,h1.data(),n*8,hipMemcpyHostToDevice);
    hipMemcpy(ddu,h2.data(),n*8,hipMemcpyHostToDevice); hipMemcpy(dx,x.data(),n*8,hipMemcpyHostToDevice);
    size_t sz=0; hipsparseDgtsvInterleavedBatch_bufferSizeExt(h,algo,m,ddl,dd,ddu,dx,batch,&sz);
    void* buf; hipMalloc(&buf,sz);
    int rc=hipsparseDgtsvInterleavedBatch(h,algo,m,ddl,dd,ddu,dx,batch,buf); hipDeviceSynchronize();
    std::vector<double> r0(n),r1(n),r2(n),rx(n);
    hipMemcpy(r0.data(),ddl,n*8,hipMemcpyDeviceToHost); hipMemcpy(r1.data(),dd,n*8,hipMemcpyDeviceToHost);
    hipMemcpy(r2.data(),ddu,n*8,hipMemcpyDeviceToHost); hipMemcpy(rx.data(),dx,n*8,hipMemcpyDeviceToHost);
    int c0=0,c1=0,c2=0; double maxres=0;
    for(int i=0;i<n;i++){ if(r0[i]!=h0[i])c0++; if(r1[i]!=h1[i])c1++; if(r2[i]!=h2[i])c2++; }
    // residual check with original coefficients
    for(int b=0;b<batch;b++) for(int k=0;k<m;k++){ double s=h1[k*batch+b]*rx[k*batch+b];
      if(k>0) s+=h0[k*batch+b]*rx[(k-1)*batch+b]; if(k<m-1) s+=h2[k*batch+b]*rx[(k+1)*batch+b];
      maxres=fmax(maxres,fabs(s-x[k*batch+b])); }
    printf("algo=%d rc=%d bufsize=%zu changed: dl=%d d=%d du=%d (of %d)  max residual=%.3e\n",algo,rc,sz,c0,c1,c2,n,maxres);
    hipFree(buf);
  }
  return 0;
}
