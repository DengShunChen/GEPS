      subroutine joinsr(wss,s1,s2,s3,s4,jtrun,jtmax,lev                 &
                       ,mlistnum,num,ncld)
!
      implicit  none
      real      wss(*), s1(*), s2(*), s3(*), s4(*)
      integer   jtrun,jtmax,lev,mlistnum,num,ncld
!
      if(num .eq. 2) call join2sr(wss,s1,s2,jtrun,jtmax,lev             &
                                 ,mlistnum,ncld)
      if(num .eq. 3) call join3sr(wss,s1,s2,s3,jtrun,jtmax,lev          &
                                 ,mlistnum,ncld)
      if(num .eq. 4) call join4sr(wss,s1,s2,s3,s4,jtrun,jtmax,lev       &
                                 ,mlistnum,ncld)
!
      return
      end
!
      subroutine join2sr(wss,s1,s2,jtrun,jtmax,lev,mlistnum,ncld)
!
      implicit  none
      real      wss (lev,2,1+ncld,jtrun,jtmax)
      real      s1(lev,2,jtrun,jtmax)
      real      s2(lev*ncld,2,jtrun,jtmax)
      integer   jtrun,jtmax,lev,mlistnum,ncld
      integer   m,l,k,n,nk,kk
!
      wss=0.
      do 10 m=1,mlistnum
      do 10 l=1,jtrun
      do 10 k = 1, lev*2
         wss(k,1,1,l,m) = s1(k,1,l,m)
  10  continue
!
      do 20 m=1,mlistnum
      do 20 l=1,jtrun
      do 20 n=1,ncld
      nk=(n-1)*lev
      do 20 k = 1, lev
         kk=nk+k
         wss(k,1,1+n,l,m) = s2(kk,1,l,m)
         wss(k,2,1+n,l,m) = s2(kk,2,l,m)
  20  continue

!
      return
      end
!
      subroutine join3sr(wss,s1,s2,s3,jtrun,jtmax,lev,mlistnum,ncld)
!
      implicit  none
      real      wss (lev,2,2+ncld,jtrun,jtmax)
      real      s1(lev,2,jtrun,jtmax)
      real      s2(lev,2,jtrun,jtmax)
      real      s3(lev*ncld,2,jtrun,jtmax)
      integer   jtrun,jtmax,lev,mlistnum,ncld
      integer   m,l,k,n,nk,kk
!
      wss=0.
      do 10 m=1,mlistnum
      do 10 l=1,jtrun
      do 10 k = 1, lev*2
         wss(k,1,1,l,m) = s1(k,1,l,m)
         wss(k,1,2,l,m) = s2(k,1,l,m)
  10  continue
!
      do 20 m=1,mlistnum
      do 20 l=1,jtrun
      do 20 n=1,ncld
      nk=(n-1)*lev
      do 20 k = 1, lev
         kk=nk+k
         wss(k,1,2+n,l,m) = s3(kk,1,l,m)
         wss(k,2,2+n,l,m) = s3(kk,2,l,m)
  20  continue
!
      return
      end
!
      subroutine join4sr(wss,s1,s2,s3,s4,jtrun,jtmax,lev,mlistnum,ncld)
!
      implicit  none
      real      wss (lev,2,3+ncld,jtrun,jtmax)
      real      s1(lev,2,jtrun,jtmax)
      real      s2(lev,2,jtrun,jtmax)
      real      s3(lev,2,jtrun,jtmax)
      real      s4(lev*ncld,2,jtrun,jtmax)
      integer   jtrun,jtmax,lev,mlistnum,ncld
      integer   m,l,k,n,nk,kk
!
      wss=0.
      do 10 m=1,mlistnum
      do 10 l=1,jtrun
      do 10 k = 1, lev*2
         wss(k,1,1,l,m) = s1(k,1,l,m)
         wss(k,1,2,l,m) = s2(k,1,l,m)
         wss(k,1,3,l,m) = s3(k,1,l,m)
  10  continue
!
      do 20 m=1,mlistnum
      do 20 l=1,jtrun
      do 20 n=1,ncld
      nk=(n-1)*lev
      do 20 k = 1, lev
         kk=nk+k
         wss(k,1,3+n,l,m) = s4(kk,1,l,m)
         wss(k,2,3+n,l,m) = s4(kk,2,l,m)
  20  continue
!
      return
      end
