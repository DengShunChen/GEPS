      subroutine pmatrx_hybrid_cwb ( cp,tmean,ptop,ptmean,spalm,phatk  &
                        , sigma,ai,dcdp )
!
!  pmatrx computes the contributions to the semi-implicit operators
!  that are a result of the differences of the actual terrain pressure
!  and pbl depth from the reference values of these parameters
!
!  ***input***
!
!  cp: specific heat of air
!  tmean: reference temperature profile for linearized atmosphere
!  ptop: pressure at model top
!  ptmean:  reference terrain pressure for linearized atmosphere
!  spalm:  linearized energy conversion term for reference atmosphere
!  phatk: even level exner function
!  sig: sigma coordinates
!  ai: inverse of bi-linear vertical differencing operator array
!
!  ***output***
!
!  dcdp: matrix operator for vertical differencing of energy term
!
! *********************************************************************
!
      use param
      use const, only : RTYPE
!
      implicit  none
  
!      real      sigma(lev+1,2),tmean(lev),a(lev,lev),b(lev,lev),ai(lev,lev),&
      real      dcdp(lev)
      real(kind=RTYPE) sigma(lev+1,2),tmean(lev),a(lev,lev),b(lev,lev),&
                       ai(lev,lev),spalm(lev),phatk(lev)
!
      integer   i,j,k
      real      ptop,cp
      real(kind=RTYPE) ptmean,capa,ps,spamin,tem,tem1,tem2,tem3,phat

      capa= 2.0/7.0
      ps= ptmean+ptop
      call zilch(b,lev*lev)
      spamin= 0.000001
      do 3 k=1,lev-1
      tem2= max(spalm(k+1),spamin)
      tem3= max(spalm(k),spamin)
      phat= sigma(k+1,1)*ptmean+sigma(k+1,2)+ptop
      tem= capa*sigma(k+1,1)/phat
!      phat= sig(k+1)*ptmean+ptop
!      tem= capa*sig(k+1)/phat
      tem= tem*phatk(k)
      tem1= tem*tmean(k+1)/tem2
      tem= tem*tmean(k)/tem3
      b(k,k)= cp*spalm(k)*(1.0-tem)
      b(k,k+1)= cp*spalm(k+1)*(tem1-1.0)
    3 continue
      tem= capa*tmean(lev)/(ps*spalm(lev))
      tem= tem*phatk(lev)
      b(lev,lev)= cp*spalm(lev)*(tem-1.0)
      call mtxmlp(ai,b,a,lev)
      do 5 i=1,lev
      dcdp(i)= 0.0
    5 continue
      do 6 j=1,lev
      do 6 i=1,lev
    6 dcdp(i)= dcdp(i)+a(i,j)
      return
      end
