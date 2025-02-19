      subroutine matrix_hybrid_cwb(cp,sigma,dsigma,ptop,ptmean,tmean,spalm &
                        ,eigval,evecin,evectr,arrhyd,arsddt                &
                        ,pmcor,tmcor )
!
!  ***input***
!
!  cp: specific heat of air
!  sig: sigma coordinate array
!  dsig: sigma coordinate layer thickness array
!  ptop: pressure at model top
!  ptmean: reference terrain pressure for linearization
!  tmean: reference temperature profile for linearization
!  spalm: linearized energy conversion term
!
!  ***output***
!
!  eigval: eigenvalues of gravity wave matrix
!  evectr: eigenvectors of gravity wave matrix
!  evecin: inverse of evectr
!  arrhyd: hydrostatic matrix for reference atmosphere
!  arsddt: matrix operator for mean vertical temperature advection
!  pmcor: vector used in future normal mode initialization
!  tmcor: matrix used in future normal mode initialization
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
! *******************************************************************
!
!  matrix computes the matrix operators used in the semi-implicit
!  time differencing scheme of the model.
!
      use param
      use const, only : RTYPE

      implicit  none
      real      cp,ptop
!
      real(kind=RTYPE) eigval(lev),evecin(lev,lev),evectr(lev,lev),    &
                pmcor(lev),tmcor(lev,lev),arrhyd(lev,lev),             &
                arsddt(lev,lev),sigma(lev+1,2),dsigma(lev,2),          &
                tmean(lev),spalm(lev),ptmean
!
      integer   ipp(lev*2),iwk(lev)
!  in order to keep the stability. two variables below can not be single precision.
      real      dp(lev),dcdp(lev)
!
      real(kind=RTYPE) a(lev,lev),b(lev,lev),asd(lev,lev),ai(lev,lev), &
                csrtn(lev,lev),enorm(lev),tmp(lev),                    &
                phatk(lev+1),pko(lev),thatm(lev),p2(lev+1)
!
!  compute some pressure related variables from reference atmosphere
!
      integer   k,j,ier,jj,i,l
      real      capa,capap1,at,det
#ifdef USE_CUDA
      real evecin_fp64(lev, lev), a_fp64(lev, lev), eigval_fp64(lev), b_fp64(lev, lev), enorm_fp64(lev)
#endif
!
      capa= 1.0/3.5
      capap1= 1.0+capa
      do 135 k=1,lev+1
      p2(k)= ptop+sigma(k,1)*ptmean+sigma(k,2)
!      p2(k)= ptop+sig(k)*ptmean
      phatk(k)= (p2(k)/1000.0)**capa
!jh      sig(k) = (p2(k)-ptop) / ptmean
  135 continue
!
      do 145 k=1,lev
      dp(k)= p2(k+1)-p2(k)
!jh      dsig(k) = sig(k+1) - sig(k)
      pko(k)= (p2(k+1)*phatk(k+1)-p2(k)*phatk(k))/(capap1*dp(k))
  145 continue
!
      do 120 k=1,lev
      tmean(k)= tmean(k)/pko(k)
      spalm(k)= tmean(k)*(sigma(k,1)*(pko(k)-phatk(k))+sigma(k+1,1)     &
       *(phatk(k+1)-pko(k)))/dp(k)
!jh      spalm(k)= tmean(k)*(sig(k)*(pko(k)-phatk(k))+sig(k+1)
!jh     * *(phatk(k+1)-pko(k)))/dp(k)
  120 continue
      do 140 k=1,lev-1
      at= (phatk(k+1)-pko(k))/(pko(k+1)-pko(k))
      thatm(k)= tmean(k)*at+tmean(k+1)*(1.0-at)
  140 continue
!
! computation of hydrostatic matrix operators
!
      call zilch (a,lev*lev)
      call zilch (b,lev*lev)
      do 2 j=1,lev-1
      a(j,j)= -1.0
      a(j,j+1)= 1.0
      b(j,j)= cp*(pko(j)-phatk(j+1))
      b(j,j+1)= cp*(phatk(j+1)-pko(j+1))
    2 continue
      a(lev,lev)= 1.0
      b(lev,lev)= cp*(phatk(lev+1)-pko(lev))
!
      call invmtx (a,lev,ai,lev,lev,det,ipp,ier)
      call mtxmlp (ai,b,arrhyd,lev)
!
      call pmatrx_hybrid_cwb(cp,tmean,ptop,ptmean,spalm,phatk(2),sigma  &
                            ,ai,dcdp)
!
      call zilch (a,lev*lev)
      do 7 k=1,lev-1
      if(k.gt.1) a(k,k)= (tmean(k)-thatm(k-1))/dp(k)
    7 a(k,k+1)= (thatm(k)-tmean(k))/dp(k)
      a(lev,lev)= (tmean(lev)-thatm(lev-1))/dp(lev)
      call zilch (asd,lev*lev)
      do 20 j=2,lev
      jj= j-1
      do 30 i=1,lev
   30   asd(j,i)=  sigma(j,1) *dp(i)
!   30 asd(j,i)= sig(j)*dp(i)
      do 20 i=1,jj
   20   asd(j,i)= (sigma(j,1)-1.0)*dp(i)
!   20 asd(j,i)= (sig(j)-1.0)*dp(i)
      call mtxmlp (a,asd,arsddt,lev)
#ifdef VERBOSE
      tmp=dcdp
      call mtxprt (tmp,lev,1,'dcdp    ','f10.4   ')
      tmp=dp
      call mtxprt (tmp,lev,1,'dp      ','f10.4   ')
#endif
!
      call mtxmlp (arrhyd,arsddt,a,lev)
      call zilch (csrtn,lev*lev)
      do 10 j=1,lev
      spalm(j)= dcdp(j)+cp*spalm(j)
   10 continue
      do 11 i=1,lev
      do 11 j=1,lev
   11 csrtn(j,i)= a(j,i)+spalm(j)*dp(i)
!
      do 15 i=1,lev
      do 15 j=1,lev
      evecin(i,j)= csrtn(i,j)
   15 continue
!
!  find eigenvalues and eigenvectors of gravity wave modes
!
!  evectr in "rg" should be a integer array
!  1/9/2004
!err  call rg (lev,lev,evecin,a,eigval,1,b,evectr,enorm,ier)
#ifdef SP
#ifdef USE_CUDA
      evecin_fp64 = evecin
      a_fp64 = a
      eigval_fp64 = eigval
      b_fp64 = b
      enorm_fp64 = enorm
      call rg (lev,lev,evecin_fp64,a_fp64,eigval_fp64,1,b_fp64,iwk,enorm_fp64,ier)
      evecin = evecin_fp64
      a = a_fp64
      eigval = eigval_fp64
      b = b_fp64
      enorm = enorm_fp64
#else
      !sigle precision version in fujitsu
      call eig1(evecin,lev,lev,0,a,eigval,b,enorm,ier)
#endif
#else
      call rg (lev,lev,evecin,a,eigval,1,b,iwk,enorm,ier)
#endif
!sun  call rg (lev,lev,evecin,eigval,a,1,evectr,b,enorm,ier)
!fuji      call deig1(evecin,lev,lev,0,a,eigval,b,enorm,ier)

!
!  normalize vertical eigenvector matrix
!
!      call orders (2,wrk,a,ipp,lev,1,8,1)
      call indexx (lev,a,ipp)
!
!  normalize vertical eigenvector matrix
!
!sun-no: i.e. cray only  ( make eigval be decending order)
      do 40 j = 1, lev
      jj = lev -j + 1
      eigval(jj)   = a(ipp(j),1)
      do 40 i = 1, lev
      evectr(i,jj) = b(i,ipp(j))
   40 continue
!sun-no: i.e. cray only
      do 45 i=1,lev
      enorm(i)= 0.0
      do 50 j=1,lev
   50 enorm(i)= enorm(i)+(dsigma(j,1)+dsigma(j,2)/ptmean)               &
                         *evectr(j,i)*evectr(j,i)
!jh1   50 enorm(i)= enorm(i)+dsigma(j,1)*evectr(j,i)*evectr(j,i)
!jh   50 enorm(i)= enorm(i)+dsig(j)*evectr(j,i)*evectr(j,i)
      enorm(i)= 1.0/sqrt(enorm(i))
   45 continue
!
      do 55 i=1,lev
      do 55 j=1,lev
   55 evectr(j,i)= evectr(j,i)*enorm(i)
!
      call invmtx (evectr,lev,evecin,lev,lev,det,ipp,ier)
#ifdef VERBOSE
      call mtxprt (arrhyd,lev,lev,'arrhyd  ','20f6.2  ')
      call mtxprt (eigval,lev,1,'eigval  ','f12.4   ')
      call mtxprt (evectr,lev,lev,'evectr  ','20f7.3  ')
      call mtxprt (evecin,lev,lev,'evecin  ','20f7.3  ')
#endif
      call invmtx (csrtn,lev,b,lev,lev,det,ipp,ier)
!
      do 21 k=1,lev
      pmcor(k)= 0.0
      do 21 l=1,lev
   21 pmcor(k)= pmcor(k)+dp(l)*b(l,k)
#ifdef VERBOSE
      call mtxprt (pmcor,lev,1,'pmcor   ','f10.4   ')
#endif
      call mtxmlp (arsddt,b,tmcor,lev)
#ifdef VERBOSE
      call mtxprt (tmcor,lev,lev,'tmcor   ','20f6.3  ')
#endif
!
      return
      end
