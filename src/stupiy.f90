      subroutine stupiy ( yr,xin,yin,mn,ix,jy,pix,pjy,ipt,tp1,bitx )
!
      implicit  none

      integer   mn,ix,jy,mn2,i,j,num,ic,ilsum
      real      yin(mn),xin(mn),pix(mn,4),pjy(mn,4),tp1(mn,4)
      real      yr(jy)
      integer   ipt(mn)
      logical   bitx(mn,2)
!
      integer,parameter :: one=1.0
!
!  find x-index of box containing each interpolation point
!  and compute coefficients for spline interpolation
!
      mn2= mn*2
      do 20 i=1,mn
      ipt(i)= xin(i)
      pix(i,3)= ipt(i)
      pix(i,3)= xin(i)-pix(i,3)
      pix(i,4)= one-pix(i,3)
   20 continue
!
      do 30 i=1,mn2
      pix(i,1)= pix(i,3)*pix(i,3)-one
      pix(i,1)= pix(i,1)*pix(i,3)
   30 continue
!
!  find y-index of box containing each interpolation point and compute
!  coefficients for spline interplation
!
      do 45 i=1,mn
      bitx(i,1)= .false.
   45 continue
      num= 0
      do 5 j=2,jy
      if(j.lt.jy) then
      do 55 i=1,mn
      bitx(i,2)= .not.bitx(i,1).and.(yin(i).lt.yr(j))
   55 continue
      else
      do 65 i=1,mn
      bitx(i,2)= .not.bitx(i,1)
   65 continue
      endif
      ic= ilsum(mn,bitx(1,2),1)
      if(ic.eq.0) go to 5
      do 75 i=1,mn
      if(bitx(i,2)) then
      ipt(i)= ipt(i)+ix*(j-1)
      tp1(i,1)= yin(i)-yr(j-1)
      tp1(i,2)= yr(j)-yr(j-1)
      endif
   75 continue
      num= num+ic
      if(num.eq.mn) go to 6
      do 85 i=1,mn
      bitx(i,1)= bitx(i,1).or.bitx(i,2)
   85 continue
    5 continue
    6 continue
!
      do 95 i=1,mn
      pjy(i,3)= tp1(i,1)/tp1(i,2)
      pjy(i,4)= one-pjy(i,3)
   95 continue
!
      do 105 i=1,mn2
      pjy(i,1)= pjy(i,3)*pjy(i,3)-one
  105 continue
!
      do 115 i=1,mn
      pjy(i,1)= pjy(i,1)*tp1(i,1)*tp1(i,2)
      tp1(i,1)= tp1(i,2)-tp1(i,1)
      pjy(i,2)= pjy(i,2)*tp1(i,1)*tp1(i,2)
  115 continue
      return
      end
