      subroutine coftrix(m,k,mx,nn,a,b,c,jtrun,lev,isym)
!
!  purpose: to create the symmetric and anti symmetric matrixes for NNMI
!-----------------------------------------------------------------------
!  **** input ****
!  m     : horizontal wave number
!  k     : vertical level
!  nn    : dimension of cofficient matrix
!  a     : coefficient array
!  b     : coefficient array
!  c     : coefficient array
!  jtrun : number of horizontol mode (symmetric or asymmetric)
!  lev   : total vertical levels
!  isym  : = +1  : symmetric case 
!          =-1  : asymmetric case
!  **** output ****
!    mx    : coefficient matrix   (nn by nn) in size
!
! modify to f90 in 2015 by C-H Lee
! sort in 2015 by River Chen
!------------------------------------------------------------------------
      real mx(nn,nn)
!      real a(jtrun,jtrun,lev),b(jtrun,jtrun,lev),c(jtrun,jtrun,lev)
      real a(jtrun,jtrun),b(jtrun,jtrun),c(jtrun,jtrun)
      do 90 i=1,nn
      do 90 j=1,nn
        mx(i,j)=0.
 90   continue
      if(isym.eq.1)then
!
! symmterix case
!
        if(mod(nn,3).eq.0)then
          i=m
          jj=1
          ni=int(nn/3)
          do 110 j=1,ni
            mx(jj,jj)=0.
!byl            mx(jj+1,jj+1)=c(m,i,k)
!byl            mx(jj+2,jj+2)=c(m,i+1,k)
            mx(jj+1,jj+1)=c(m,i)
            mx(jj+2,jj+2)=c(m,i+1)
            i=i+2
            jj=jj+3
 110      continue
          i=m
          jj=1
          ni=int((nn-1)/3)
          do 120 j=1,ni
!byl            mx(jj,jj+1)=b(m,i,k)
!byl            mx(jj+1,jj+2)=a(m,i+1,k)
            mx(jj,jj+1)=b(m,i)
            mx(jj+1,jj+2)=a(m,i+1)
            mx(jj+2,jj+3)=0.
            i=i+2
            jj=jj+3
 120      continue
!byl            mx(nn-2,nn-1)=b(m,i,k)
!byl            mx(nn-1,nn)=a(m,i+1,k)
            mx(nn-2,nn-1)=b(m,i)
            mx(nn-1,nn)=a(m,i+1)
          i=m
          jj=1
          ni=int((nn-2)/3)
          do 130 j=1,ni
            mx(jj,jj+2)=0.
            mx(jj+1,jj+3)=0.
!byl            mx(jj+2,jj+4)=a(m,i+2,k)
            mx(jj+2,jj+4)=a(m,i+2)
            i=i+2
            jj=jj+3
 130      continue
            mx(nn-2,nn)=0.
        else
          i=m
          jj=1
          ni=int(nn/3)
          do 150 j=1,ni
            mx(jj,jj)=0.
!byl            mx(jj+1,jj+1)=c(m,i,k)
!byl            mx(jj+2,jj+2)=c(m,i+1,k)
            mx(jj+1,jj+1)=c(m,i)
            mx(jj+2,jj+2)=c(m,i+1)
            i=i+2
            jj=jj+3
 150      continue
            mx(nn-1,nn-1)=0.
!byl            mx(nn,nn)=c(m,i,k)
            mx(nn,nn)=c(m,i)
          i=m
          jj=1
          ni=int((nn-1)/3)
          do 160 j=1,ni
!byl            mx(jj,jj+1)=b(m,i,k)
!byl            mx(jj+1,jj+2)=a(m,i+1,k)
            mx(jj,jj+1)=b(m,i)
            mx(jj+1,jj+2)=a(m,i+1)
            mx(jj+2,jj+3)=0.
            i=i+2
            jj=jj+3
 160      continue
!byl            mx(nn-1,nn)=b(m,i,k)
            mx(nn-1,nn)=b(m,i)
          i=m
          jj=1
          ni=int((nn-2)/3)
          do 170 j=1,ni
            mx(jj,jj+2)=0.
            mx(jj+1,jj+3)=0.
!byl            mx(jj+2,jj+4)=a(m,i+2,k)
            mx(jj+2,jj+4)=a(m,i+2)
            i=i+2
            jj=jj+3
 170      continue
        endif
!
!  anti-symmetrix case
!
      else if(isym.eq.-1)then
        if(mod(nn,3).eq.0)then
          i=m
          jj=1
          ni=int(nn/3)
          do 190 j=1,ni
!byl            mx(jj,jj)=c(m,i,k)
            mx(jj,jj)=c(m,i)
            mx(jj+1,jj+1)=0.
!byl            mx(jj+2,jj+2)=c(m,i+1,k)
            mx(jj+2,jj+2)=c(m,i+1)
            i=i+2
            jj=jj+3
 190      continue
          i=m
          jj=1
          ni=int((nn-1)/3)
          do 200 j=1,ni
            mx(jj,jj+1)=0.
!byl            mx(jj+1,jj+2)=b(m,i+1,k)
!byl            mx(jj+2,jj+3)=a(m,i+2,k)
            mx(jj+1,jj+2)=b(m,i+1)
            mx(jj+2,jj+3)=a(m,i+2)
            i=i+2
            jj=jj+3
 200      continue
            mx(nn-2,nn-1)=0.
!byl            mx(nn-1,nn)=b(m,i+1,k)
            mx(nn-1,nn)=b(m,i+1)
          i=m
          jj=1
          ni=int((nn-2)/3)
          do 210 j=1,ni
!byl            mx(jj,jj+2)=a(m,i+1,k)
            mx(jj,jj+2)=a(m,i+1)
            mx(jj+1,jj+3)=0.
            mx(jj+2,jj+4)=0.
            i=i+2
            jj=jj+3
 210      continue
!byl            mx(nn-2,nn)=a(m,i+1,k)
            mx(nn-2,nn)=a(m,i+1)
        else
          i=m
          jj=1
          ni=int(nn/3)
          do 230 j=1,ni
!byl            mx(jj,jj)=c(m,i,k)
            mx(jj,jj)=c(m,i)
            mx(jj+1,jj+1)=0.
!byl            mx(jj+2,jj+2)=c(m,i+1,k)
            mx(jj+2,jj+2)=c(m,i+1)
            i=i+2
            jj=jj+3
 230      continue
!byl            mx(nn,nn)=c(m,i,k)
            mx(nn,nn)=c(m,i)
          i=m
          jj=1
          ni=int((nn-1)/3)
          do 240 j=1,ni
            mx(jj,jj+1)=0.
!byl            mx(jj+1,jj+2)=b(m,i+1,k)
!byl            mx(jj+2,jj+3)=a(m,i+2,k)
            mx(jj+1,jj+2)=b(m,i+1)
            mx(jj+2,jj+3)=a(m,i+2)
            i=i+2
            jj=jj+3
 240      continue
          i=m
          jj=1
          ni=int((nn-2)/3)
          do 250 j=1,ni
!byl            mx(jj,jj+2)=a(m,i+1,k)
            mx(jj,jj+2)=a(m,i+1)
            mx(jj+1,jj+3)=0.
            mx(jj+2,jj+4)=0.
            i=i+2
            jj=jj+3
 250      continue
            if (nn.ne.1) then
!byl            mx(nn-3,nn-1)=a(m,i+1,k)
            mx(nn-3,nn-1)=a(m,i+1)
            mx(nn-2,nn)=0.
            endif
        endif
      endif
!
!   copy upper half to lower half
!
      jj=0
      do 270 n=nn-1,1,-1
        jj=jj+1
        do 260 i=1,n
          j=jj+i 
 260      mx(j,i)=mx(i,j)
 270  continue
!      call smetrix(mx,nn)
!
      return
      end
