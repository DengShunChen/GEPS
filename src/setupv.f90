      subroutine setupv (yr,yin,mn,nx,lmaxp,pjy,ipt,tp1)
!
!  setup routine for bicubv parameters.  see bicubv prolog for
!  parameter descriptions
!
      use const, only: RTYPE
!
      implicit  none
      integer   mn,nx,lmaxp,ipt(mn)
      real(kind=RTYPE) pjy(mn,4)
      real      tp1(mn,4),yin(mn),yr(nx,lmaxp)

      integer   n,i,j,k,ii,ix
      real      eps
!
      n= mn/nx
!
!  check out of bounds to ensure interpolation
!
      eps = 1.0e-5
      do 100 i = 1, mn
      ix = i - ((i-1)/nx)*nx
      if ( yin(i).le.yr(ix,1) )      yin(i) = yr(ix,1) + eps
      if ( yin(i).ge.yr(ix,lmaxp) )  yin(i) = yr(ix,lmaxp) - eps
  100 continue
!
      do 10 i=1,nx
      k=1
      do 5 j=1,n
      ii= nx*(j-1)
    6 k= k+1
      if(yin(i+ii).gt.yr(i,k)) go to 6
      ipt(i+ii)= i+nx*(k-1)
      tp1(i+ii,1)= yin(i+ii)-yr(i,k-1)
      tp1(i+ii,2)= yr(i,k)-yr(i,k-1)
      k= k-1
    5 continue
   10 continue
!
      do 90 i=1,mn
      pjy(i,3)= tp1(i,1)/tp1(i,2)
      pjy(i,4)= 1.0-pjy(i,3)
      pjy(i,1)= pjy(i,3)*pjy(i,3)-1.0
      pjy(i,2)= pjy(i,4)*pjy(i,4)-1.0
      pjy(i,1)= pjy(i,1)*tp1(i,1)*tp1(i,2)
      tp1(i,1)= tp1(i,2)-tp1(i,1)
      pjy(i,2)= pjy(i,2)*tp1(i,1)*tp1(i,2)
   90 continue
!
      return
      end
