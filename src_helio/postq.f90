      subroutine postq ( nxj,nx,lev,ncld,dsigma,pt,q,qbarrw,cosl )

      implicit  none
      integer   nxj,nx,lev,ncld
      real      cosl,qmin
      real      q(nx,lev*ncld),dsigma(lev,2),pt(nx),qbarrw(ncld)
      data qmin/1.0e-10/

      integer   k,n,kk,i,klev
      real      fac,dsigp
!
      do 10 k = 1,lev-1
      do 15 n=1,ncld
        kk=(n-1)*lev+k
      do 20 i = 1,nxj
        if (q(i,kk).lt.qmin)  then
        fac= (dsigma(k  ,1)*pt(i)+dsigma(k  ,2))/   &
             (dsigma(k+1,1)*pt(i)+dsigma(k+1,2))
          q(i,kk+1) = q(i,kk+1)+fac*q(i,kk)-qmin
          q(i,kk) = qmin
        endif
 20   continue
 15   continue
 10   continue
!
      do n=1,ncld
        klev=(n-1)*lev+lev
      do 30 i = 1,nxj
        if (q(i,klev).lt.qmin)  then
          dsigp = dsigma(lev,1)*pt(i)+dsigma(lev,2)
          qbarrw(n) = qbarrw(n)+dsigp*(qmin-q(i,klev))*cosl
          q(i,klev) = qmin
        endif
 30   continue
      enddo
!
      return
      end
