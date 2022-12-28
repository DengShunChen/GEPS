!fj   subroutine zx(evec,vorten,divten,phiten,mlmax,lev)
      subroutine zx(evec,vorten,divten,phiten,jtrun,jtmax,lev)
!
!  purpose : do vertical transform and inverse transform
!------------------------------------------------------------------------
!  **** input ****
!  evec   : matrix array for transform
!  vorten : vorticity tendency(correction) array of spectrum coefficients
!  divten : divergence tendency(correction) array of spectrum coefficients
!  phiten : geopotential tendency(correction)array of spectrum coefficients
!  lev    : total vertical levels
!  **** output ****
!  vorten : vorticity tendency(correction) array of spectrum coefficients
!  divten : divergence tendency(correction) array of spectrum coefficients
!  phiten : geopotential tendency(correction)array of spectrum coefficients
!---------------------------------------------------------------------------
      use index
      use const, only: RTYPE

      implicit none

      integer  lev,jtrun,jtmax,m,mf,n,j,l,k,KK,KL

      real(kind=RTYPE) evec(lev,lev)
      real(kind=RTYPE) vorten(levp,2,jtrun,jtmax),        &
                       phiten(levp,2,jtrun,jtmax),        &
                       divten(levp,2,jtrun,jtmax)
      real(kind=RTYPE) vor(lev,2,jtrun,jtmax),            &
                       div(lev,2,jtrun,jtmax),            &
                       phe(lev,2,jtrun,jtmax)
!2dMPI
      REAL(kind=RTYPE) wrk(lev,2,3,jtrun,jtmax)
      CALL mpe2d_unify_spec_lev_zx(wrk,vorten,divten,phiten,lev,levp,jtrun,jtmax,mlistnum,nsizex,row_comm)

!
      do 120 m = 1, mlistnum
         mf=mlist(m)
      do 120 n = mf, jtrun
      do 120 j = 1, 2
      do 105 L = 1, lev
!       vor(L,j,n,m)=evec(L,1)*vorten(1,j,n,m)
!       div(L,j,n,m)=evec(L,1)*divten(1,j,n,m)
!       phe(L,j,n,m)=evec(L,1)*phiten(1,j,n,m)
        vor(L,j,n,m)=evec(L,1)*wrk(1,j,1,n,m)
        div(L,j,n,m)=evec(L,1)*wrk(1,j,2,n,m)
        phe(L,j,n,m)=evec(L,1)*wrk(1,j,3,n,m)
      do 100 k = 2, lev
!       vor(L,j,n,m)=vor(L,j,n,m)+evec(L,k)*vorten(k,j,n,m)
!       div(L,j,n,m)=div(L,j,n,m)+evec(L,k)*divten(k,j,n,m)
!       phe(L,j,n,m)=phe(L,j,n,m)+evec(L,k)*phiten(k,j,n,m)
        vor(L,j,n,m)=vor(L,j,n,m)+evec(L,k)*wrk(k,j,1,n,m)
        div(L,j,n,m)=div(L,j,n,m)+evec(L,k)*wrk(k,j,2,n,m)
        phe(L,j,n,m)=phe(L,j,n,m)+evec(L,k)*wrk(k,j,3,n,m)
 100  continue
 105  continue
 120  continue
!
      do 110 m = 1, mlistnum
        mf=mlist(m)
      do 110 n = mf, jtrun
!     do 110 k = 1, lev*2
!       vorten(k,1,n,m)=vor(k,1,n,m)
!       divten(k,1,n,m)=div(k,1,n,m)
!       phiten(k,1,n,m)=phe(k,1,n,m)
      do 110 k = 1, levp
         kk=Llist(k)
         vorten(k,1,n,m)=vor(kk,1,n,m)
         vorten(k,2,n,m)=vor(kk,2,n,m)
         divten(k,1,n,m)=div(kk,1,n,m)
         divten(k,2,n,m)=div(kk,2,n,m)
         phiten(k,1,n,m)=phe(kk,1,n,m)
         phiten(k,2,n,m)=phe(kk,2,n,m)
 110  continue
!
      return
      end
