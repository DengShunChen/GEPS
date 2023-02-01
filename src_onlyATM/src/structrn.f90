      subroutine structrn(vorten,divten,phiten,mlmax,lev,   &
                           mlsort,mnsort,jtrun,iflag)
!
!  purpose : change array structure between forecast model and
!  initialization
!            module
!----------------------------------------------------------------------------
!  **** input ****
!  vorten : vorticity tendency(correction) array of spectrum coefficients
!  divten : divergence tendency(correction) array of spectrum coefficients
!  phiten : geopotential tendency(correction)array of spectrum coefficients
!  mlmax  : total number of spectrum coefficient
!  lev    : total vertical levels
!  mlsort : array index of forecast model
!  jtrun  : horizontal wave number truncation limit
!  iflag  : +2  from forecast model to initialization module
!           -2  from initialization module to forecast model
!  **** output ****
!  vorten : vorticity tendency(correction) array of spectrum coefficients
!  divten : divergence tendency(correction) array of spectrum coefficients
!  phiten : geopotential tendency(correction)array of spectrum coefficients
!----------------------------------------------------------------------
      use paramt

      implicit none

      integer  mlmax,lev,iflag,k,m,jtrun,l,ml,mn,mlmax2

      real vorten(mlmax,2,lev),divten(mlmax,2,lev),phiten(mlmax,2,lev)
      integer mlsort(jtrun,jtrun),mnsort(jtrun,jtrun)
      real      vor(mlm,2,lm),div(mlm,2,lm),phe(mlm,2,lm)
!
!  from forecast model to initialization
!
      if(iflag.eq.2)then
!mic$  do all
!mic$1 private (k,m,l,ml,mn)
!mic$2 shared  (lev,jtrun,mlsort,mnsort,vor,div,phe)
!mic$3 shared  (vorten,divten,phiten)
      do 100 k=1,lev
      do 100 m=1,jtrun
      do 100 l=m,jtrun
      ml=mlsort(m,l)
      mn=mnsort(m,l)
      vor(mn,1,k)=vorten(ml,1,k)
      vor(mn,2,k)=vorten(ml,2,k)
      div(mn,1,k)=divten(ml,1,k)
      div(mn,2,k)=divten(ml,2,k)
      phe(mn,1,k)=phiten(ml,1,k)
      phe(mn,2,k)=phiten(ml,2,k)
 100  continue
!
!  form initialization to forecsat model
!
      else
!mic$  do all
!mic$1 private (k,m,l,ml,mn)
!mic$2 shared  (lev,jtrun,mlsort,mnsort,vor,div,phe)
!mic$3 shared  (vorten,divten,phiten)
      do 110 k=1,lev
      do 110 m=1,jtrun
      do 110 l=m,jtrun
      ml=mlsort(m,l)
      mn=mnsort(m,l)
      vor(ml,1,k)=vorten(mn,1,k)
      vor(ml,2,k)=vorten(mn,2,k)
      div(ml,1,k)=divten(mn,1,k)
      div(ml,2,k)=divten(mn,2,k)
      phe(ml,1,k)=phiten(mn,1,k)
      phe(ml,2,k)=phiten(mn,2,k)
 110  continue
      endif
!
      mlmax2=mlmax*2
!mic$  do all
!mic$1 private (k,ml)
!mic$2 shared  (lev,mlmax2,vor,div,phe,vorten,divten,phiten)
      do 120 k=1,lev
      do 120 ml=1,mlmax2
      vorten(ml,1,k)=vor(ml,1,k)
      divten(ml,1,k)=div(ml,1,k)
      phiten(ml,1,k)=phe(ml,1,k)
 120  continue
!
      return
      end
