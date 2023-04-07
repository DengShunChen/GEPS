      subroutine vartrix(vorten,divten,phiten,kk,k,          &
                         x,nn,m,nbig,jtrun,jtmax,isym,iflag)
!
!  purpose: Scatter between vorten, divten, phiten and x vector
!-----------------------------------------------------------------------
!  ****  input ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!  x       : vector matrix( when iflag=-2)
!  kk      : total vertical levels
!  nn      : total number of variables in one dimension
!            i.e. the size of the vector matrixs
!  m       : horizontal wave number
!  nbig    : how many  wave number for initialization
!  isym    : = +1   :  symmetric case
!            = -1   :  asymmetric case
!  iflag   : = +2   :  variable tendency  to x vector (5.43)
!              -2   :  x vector(correction) to (vors,divs,phis)
!
!  **** output ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!  x       : vector matrix(when iglag=+2)
!-----------------------------------------------------------------------
      use index
      use const, only: RTYPE

      implicit none

      integer  iflag,isym,nbig,m,k,kk,jtrun,jtmax,nn,l,j,n,i

      real(kind=RTYPE) vorten(kk,2,jtrun,jtmax),  &
                       divten(kk,2,jtrun,jtmax),  &
                       phiten(kk,2,jtrun,jtmax)
!
        real x(nn,2)
!

        l=mlist(m)

        j=0
        n=nbig/2

        if(iflag.eq.2)then
! symmetric case
          if(isym.eq.1)then
            do 100 i=1,n
              j=j+1
              x(j,1)=phiten(k,1,l,m)
              x(j,2)=phiten(k,2,l,m)
              j=j+1
              x(j,1)=divten(k,1,l,m)
              x(j,2)=divten(k,2,l,m)
              j=j+1
              x(j,1)=vorten(k,1,l+1,m)
              x(j,2)=vorten(k,2,l+1,m)
              l=l+2
 100        continue
            if(mod(nbig,2).ne.0)then
              j=j+1
              x(j,1)=phiten(k,1,l,m)
              x(j,2)=phiten(k,2,l,m)
              j=j+1
              x(j,1)=divten(k,1,l,m)
              x(j,2)=divten(k,2,l,m)
            endif
          else if(isym.eq.-1)then
!  antiaymmetrix case
            do 110 i=1,n
              j=j+1
              x(j,1)=vorten(k,1,l,m)
              x(j,2)=vorten(k,2,l,m)
              j=j+1
              x(j,1)=phiten(k,1,l+1,m)
              x(j,2)=phiten(k,2,l+1,m)
              j=j+1
              x(j,1)=divten(k,1,l+1,m)
              x(j,2)=divten(k,2,l+1,m)
              l=l+2
 110        continue
            if(mod(nbig,2).ne.0)then
              j=j+1
              x(j,1)=vorten(k,1,l,m)
              x(j,2)=vorten(k,2,l,m)
            endif
          endif
        else if(iflag.eq.-2)then
! inverse transform
          if(isym.eq.1)then
! symmmetrix case
            do 120 i=1,n
              j=j+1
              phiten(k,1,l,m)=x(j,1)
              phiten(k,2,l,m)=x(j,2)
              j=j+1
              divten(k,1,l,m)=x(j,1)
              divten(k,2,l,m)=x(j,2)
              j=j+1
              vorten(k,1,l+1,m)=x(j,1)
              vorten(k,2,l+1,m)=x(j,2)
              l=l+2
 120        continue
            if(mod(nbig,2).ne.0)then
              j=j+1
              phiten(k,1,l,m)=x(j,1)
              phiten(k,2,l,m)=x(j,2)
              j=j+1
              divten(k,1,l,m)=x(j,1)
              divten(k,2,l,m)=x(j,2)
            endif
          else if(isym.eq.-1)then
!  antiammetrix case
            do 130 i=1,n
              j=j+1
              vorten(k,1,l,m)=x(j,1)
              vorten(k,2,l,m)=x(j,2)
              j=j+1
              phiten(k,1,l+1,m)=x(j,1)
              phiten(k,2,l+1,m)=x(j,2)
              j=j+1
              divten(k,1,l+1,m)=x(j,1)
              divten(k,2,l+1,m)=x(j,2)
              l=l+2
 130        continue
            if(mod(nbig,2).ne.0)then
              j=j+1
              vorten(k,1,l,m)=x(j,1)
              vorten(k,2,l,m)=x(j,2)
            endif
          endif
        endif
!
        return
        end
