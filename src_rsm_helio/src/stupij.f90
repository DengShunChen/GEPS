      subroutine stupij (xin,yin,mn,ix,pix,pjy,ijpt)
!
      implicit  none

      real      xin(mn),yin(mn),pix(mn,4),pjy(mn,4)
      integer   ijpt(mn),mn,ix,i,ipt,jpt
!
      do 20 i=1,mn
      ipt= int(xin(i))
      jpt= int(yin(i))
      pix(i,3)= xin(i)-float(ipt)
      pjy(i,3)= yin(i)-float(jpt)
      pix(i,4)= 1.0-pix(i,3)
      pjy(i,4)= 1.0-pjy(i,3)
      ijpt(i)= ipt+ix*jpt
   20 continue
!
      do 30 i=1,mn*2
      pix(i,1)= pix(i,3)*(pix(i,3)*pix(i,3)-1.0)
      pjy(i,1)= pjy(i,3)*(pjy(i,3)*pjy(i,3)-1.0)
   30 continue
      return
      end
