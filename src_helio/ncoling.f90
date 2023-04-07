      subroutine ncoling(jtrun,jtmax,nx,my,my_max,lev,radsq      &
      ,                   dta,temnow,eps4,cosl,pk)
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use index
      use mpe
      use const, only: RTYPE

      implicit  none
      integer   jtrun,jtmax,nx,my,my_max,lev

      real      temnow(lev,2,jtrun,jtmax),                       &
                pk(nx,lev,my_max),cosl(my),radsq,dta
      real(kind=RTYPE) eps4(jtrun,jtmax)
!
       integer,  parameter :: levtop=2,lev2=levtop/2
       real      tem(lev,2,jtrun,jtmax),ctime,cool

       integer   k,j,i,jj,nxj,m,mf,n
       real      wt,pk1,pk2,opk,fac,c1
!
!      data ctime/6./
!      data ctime/12./
!      data ctime/4./
      data ctime/24./
!
!
      real gsum(3,my)
!
      cool=1.0/(ctime*3600.*jtrun*(jtrun+1)/radsq)
!
      tem=temnow
!
      do k=1,levtop
!
      wt=0.
      pk1=0.
      pk2=0.
      gsum=0.
!
      do 100 jj =1,jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 100 i=1,nxj
        gsum(1,j)=gsum(1,j)+cosl(j)
        gsum(2,j)=gsum(2,j)+pk(i,k,jj)*cosl(j)
        gsum(3,j)=gsum(3,j)+pk(i,k+1,jj)*cosl(j)
 100  continue
!
      call mpe_unify(gsum,3,my,2,mpe_double)
!
      do j=1,my
        wt =wt +gsum(1,j)
        pk1=pk1+gsum(2,j)
        pk2=pk2+gsum(3,j)
      enddo
!
      pk1=pk1/wt
      pk2=pk2/wt
      opk=pk2/pk1
!
      fac=1.
      c1=dta*cool*fac
!
      do 110 m =1,mlistnum
         mf=mlist(m)
      do 110 n  = mf,jtrun
        temnow(k,1,n,m)=tem(k,1,n,m)+c1*eps4(n,m)*tem(k+1,1,n,m)     &
                       *opk/(1.+c1*eps4(n,m))
        temnow(k,2,n,m)=tem(k,2,n,m)+c1*eps4(n,m)*tem(k+1,2,n,m)     &
                       *opk/(1.+c1*eps4(n,m))
 110  continue
!
      enddo
!
      return
      end
