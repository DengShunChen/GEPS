      subroutine rcmpir24( nr,nx,lev,bf,ee,ww1,ww2,ww,ttau,ffu,ffd )
!
! **********************************************************************
! In this subroutine, the two- and four- stream combination  scheme  or
! the source function technique (Toon et al. 1989) is used to calculate
! the IR radiative fluxes. The exponential approximation for the Planck
! function in optical depth is used ( Fu, 1991).
! At IR wavelengths, the two-stream results are not exact in the limit
! of no scattering. It also introduces large error in the case of sca-
! ttering. Since the no-scattering limit is of considerable significance
! at IR wavelengths, we have used  the source function technique  that
! would be exact in the limit of the pure absorption and would also en-
! hance the accuracy of the two-stream approach when scattering occurs
! in the IR wavelengths.
! Here, we use nq Gauss points to obtain the fluxes: when nq=2, we use
! double Gaussian quadrature as in Fu and Liou (1993) for  four-stream
! approximation; when nq = 3, we use the regular Gauss quadrature  but
! u1*w1+u2*w2+u3*w3=1.0.
!
! ++ input/output variables:
!    (1) input variables :
!        nr    : number of grid points to calculate ir flux
!        nx    : horizontal dimension of input array
!        lev   : number of vertical layers
!        bf    : blockbody intensity function integrated over
!                the band "ib"
!                bf(1~lev+1) : atmospheric levels
!                bf(lev+2)   : surface level
!        ee    : ground surface emissivity
!        ww1~2 : first two coeff. of legendre polynomials for
!                the scattering phase function
!        ww    : single scattering albedo
!        ttau(lev+1) : total optical depth from top to reference level
!
!    (2) output variables :
!        ffu
!        ffd
!
! ++ auther :
!    original : fu and liou (1992)
!    modified : c. t. fong (1996)
!    modify to f90 by C-H Lee in 2015
!
! **********************************************************************
      use paramt

      implicit  none
!
!fong parameter ( nq = 1 )  ! souce function techn. for two streams
      integer, parameter :: nq = 2 ! .............. for four streams
!fong parameter ( nq = 3 )  ! ..................... for six streams
!
      integer   nr,nx,lev

      real      bf(nx*(lev+2)),ee(nx),ww1(nx*lev),ww2(nx*lev),ww(nx*lev), &
                ttau(nx*(lev+1))
      real      ffu(nx*(lev+1)),ffd(nx*(lev+1))
!
!   local working arrays
!
!     include '../include/paramt.h'
!
      real      w1(ilm),w(ilm),qu0(ilm),f0(ilm),tau(im*(lm+1))
      real      lamdan(ilm)
      real      ug(nq),wg(nq),ugwg(nq)
      real      gamman(ilm),caddn(ilm),                &
                cminn(ilm),caddn0(ilm),cminn0(ilm),    &
                aa(ilm),bb(ilm),expn(ilm)
      real      f(ilm),xn(ilm),yn(ilm),zn(ilm)
      real      etri(im,lm*2,4),alfa(im*(lm+1)),beta(ilm)
      real      fg(ilm),fh(ilm),fj(ilm),fk(ilm)
      real      fiu(im*(lm+1),nq),fid(im*(lm+1),nq),fx(ilm,nq),&
                fz1(ilm,nq),fz2(ilm,nq),fuq1(ilm),fuq2(ilm) 
      real      wk1(im),wk2(im)
!
!6st  data ug / 0.238619, 0.661209, 0.932469 /
!6st  data wg / 0.467914, 0.360762, 0.171324 /
!6st  data ugwg / 0.109475, 0.233886, 0.156639 /
      data ug / 0.2113248, 0.7886752 /
      data wg / 0.5, 0.5 /
      data ugwg / 0.105662, 0.394338 /
!2st  data ug / 0.5 /
!2st  data wg / 1.0 /
!2st  data ugwg / 0.5 /
!
      integer   nrn,nrn1,lev2,i,j,jk,jj,j2n,j2n1,jkes,jke,jje,j2m1,j2m,nst
      real      sqrt3,rsqrt3,pi,ugts1,fw,q1,dtau,q2,x,gamma1,gamma2,z,rsfc
      real      ssfc,wm1,wm2,y1,y,rugts1,xx,xxm1

      nrn  = nr*lev
      nrn1 = nr*(lev+1)
      lev2 = lev*2
      sqrt3= sqrt(3.)
      rsqrt3= 1./sqrt3
      pi   = 4.*atan(1.)
      ugts1 = 0.5
!
!fong01
! incorporate a delta-function adjustment to account for the forward
! diffraction  peak in the context of the two stream approximations.
! w1(n),  w(n), and tau(n+1) are the adjusted parameters.
!
      do 6 i =1, nrn
      f(i) = ww2(i) / 5.0
    6 continue
!
      do 7 i = 1, nr
      tau(i) = 0.
    7 continue
!
      do 10 i = 1, nrn
      fw = 1.0 - f(i) * ww(i)
      w1(i) = ( ww1(i) - 3.0 * f(i) ) / ( 1.0 - f(i) )
      w(i)  = ( 1.0 - f(i) ) * ww(i) / fw
      f0(i) = (ttau(i+nr) - ttau(i)) * fw
   10 continue
!
      do 12 j = 1, lev
      jk = j*nr
      jj = jk-nr
      do 12 i = 1, nr
      tau(i+jk) = tau(i+jj) + f0(i+jj)
   12 continue
!fong01
!
!fong01.1
      do 15 i = 1, nrn
      q1 = log( bf(i+nr)/bf(i) )
      if( abs(q1) .lt. 1.0e-10 ) q1 = 1.0e-10
      dtau = tau(i+nr) - tau(i)
      dtau = max(0.000001, dtau)
!cc   q2 = 1.0 / ( tau(i+nr) - tau(i) )
      q2 = 1.0 / dtau
      f0(i)  = 2.0 * ( 1.0 - w(i) ) * bf(i)
      qu0(i) = - 1.0 / ( q1 * q2 )
      beta(i)= -1.0 / qu0(i)
   15 continue
!fong01.1
!
!fong02
!
      do 20 i = 1, nrn
      w(i) = min ( w(i), 0.999999 )
   20 continue
!
!
! The following is used to calculate the coefficients for hemispheric
! mean two streams.
! The  hemispheric  mean scheme is derived by assuming that the phase
! function is equal to 1 + g in  the forward scattering hemisphere
! and 1 - g  in  the  backward scattering hemisphere where g is the
! asymmetry factor.   The hemispheric mean is only used for infrared
! wavelengths (Toon et al. 1989).
!
      do 40 i = 1, nrn
      x = w(i) * w1(i) / 3.0
      gamma1 = 2.0 - w(i) - x
      gamma2 = w(i) - x
!vpp  lamdan(i) = exp( 0.5*log((gamma1 + gamma2)*(gamma1 - gamma2)) )
      lamdan(i) = sqrt((gamma1 + gamma2)*(gamma1 - gamma2))
      gamman(i) = gamma2 / ( gamma1 + lamdan(i) )
      fw = pi * f0(i)
      x = exp ( beta(i) * ( tau(i+nr) - tau(i) ) )
      z = lamdan(i) * lamdan(i) - beta(i) * beta(i)
      if(abs(z) .lt. 1.0e-6)then
       beta(i) = beta(i)+beta(i)*0.01
       z = lamdan(i) * lamdan(i) - beta(i) * beta(i)
      end if
      caddn0(i) = fw * ( ( gamma1 + beta(i) ) + gamma2 ) / z
      cminn0(i) = fw * ( ( gamma1 - beta(i) ) + gamma2 ) / z
      caddn(i) = caddn0(i) * x
      cminn(i) = cminn0(i) * x
   40 continue
!fong02
!
!fong03
      do 50 i = 1, nrn
      dtau = tau(i+nr) - tau(i)
      expn(i) = exp ( - lamdan(i) * dtau )
      xn(i) = gamman(i) * expn(i)
      yn(i) = ( expn(i) - gamman(i) ) / ( xn(i) - 1.0 )
      zn(i) = ( expn(i) + gamman(i) ) / ( xn(i) + 1.0 )
   50 continue
!fong03
!
!fong04
!
!  etri : tridiagonal matrix for eq(39) in toon et al. (1989)
!  etri(i,k,m) : i=1,nr (horizontal grids), k=1,2*lev (vertical levels)
!                mm=1,4 (1:a, 2:b, 3:d, 4:e, a,b,d,e refer to eq(39))
!
      do 110 i = 1, nr
      etri(i,1,1) = 0.0
      etri(i,1,2) = xn(i) + 1.0
      etri(i,1,3) = xn(i) - 1.0
      etri(i,1,4) = - cminn0(i)
  110 continue
!
      do 150 j = 1, lev-1
      jk = j*nr
      jj = jk-nr
      j2n = j+j
      j2n1= j2n+1
      do 150 i = 1, nr
      etri(i,j2n,1) = 1.0+xn(i+jj)-yn(i+jk)*(gamman(i+jj)+expn(i+jj))
      etri(i,j2n,2) = 1.0-xn(i+jj)-yn(i+jk)*(gamman(i+jj)-expn(i+jj))
      etri(i,j2n,3) = yn(i+jk)*(1.0+xn(i+jk))-expn(i+jk)-gamman(i+jk)
      etri(i,j2n,4) = caddn0(i+jk)-caddn(i+jj)-yn(i+jk) *             &
                ( cminn0(i+jk)-cminn(i+jj) )
      etri(i,j2n1,1) = gamman(i+jj)-expn(i+jj)-zn(i+jj)*(1.0-xn(i+jj))
      etri(i,j2n1,2) = -1.0-xn(i+jk)+zn(i+jj)*(expn(i+jk)+gamman(i+jk))
      etri(i,j2n1,3) = zn(i+jj)*(expn(i+jk)-gamman(i+jk))-xn(i+jk)+1.0
      etri(i,j2n1,4) = cminn0(i+jk)-cminn(i+jj)-zn(i+jj) *            &
                 ( caddn0(i+jk)-caddn(i+jj) )
  150 continue
!fong04
!
!fong05
      jkes= (lev+1)*nr
      jke = lev*nr
      jje = jke-nr
      do 170 i = 1, nr
      rsfc = 1.0 - ee(i)
      ssfc = pi * ( bf(i+jkes) * ee(i) )
      wm1 = 1.0 - rsfc * gamman(i+jje)
      wm2 = xn(i+jje) - rsfc * expn(i+jje)
      etri(i,lev2,1) = wm1 + wm2
      etri(i,lev2,2) = wm1 - wm2
      etri(i,lev2,3) = 0.0
      etri(i,lev2,4) = rsfc * cminn(i+jje) - caddn(i+jje) + ssfc
  170 continue
!fong05
!
!fong06
!  solve a tridiagonal matrix using gaussian elimination
!
      call trigau ( etri, nx, nr, lev2, lev2, 0 )
!fong06
!
!fong07
!
!  aa : Y1 + Y2 as in the table 3 of toon et al. (1989)
!  bb : Y1 - Y2    .......
!
      do 180 j = 1, lev
      jk = j*nr
      jj = jk-nr
      j2m1 = j+j-1
      j2m  = j+j
      do 180 i = 1, nr
      aa(i+jj) = etri(i,j2m1,4) + etri(i,j2m,4)
      bb(i+jj) = etri(i,j2m1,4) - etri(i,j2m,4)
  180 continue
!fong07
!
!fong08
!
!  compute parameters of the source function, eq(53) and (54)
!  in toon et al. (1989)
!  notice : there is a little difference here from original eqs.
!           because of the different approximation to an emitting
!           source of plank function used by fu(1992)
!
!  fg, fh, fj, fk : parameter G, H, J, K of the table 3 respectively
!  fuq1 * (2*pi*bf) : the rest of terms for upward source function
!  fuq2 * (2*pi*bf) : the rest of terms for downward source function
!
      do 210 i =1, nrn
      x = 2.0 * ( 1.0 - w(i) ) * w(i) / ( lamdan(i) *   &
          lamdan(i) - beta(i) * beta(i) )
      y1 = w1(i) / 3.0
      y = 2.0 * ( 1.0 - w(i) * y1 )
      z = -y1 * beta(i)
      fuq1(i) = x * ( y - z ) + 1.0 - w(i)
      fuq2(i) = x * ( y + z ) + 1.0 - w(i)
  210 continue
!
      do 220 i = 1, nrn1
      alfa(i) = 2. * pi * bf(i)
  220 continue
!
      rugts1 = 1.0 / ugts1
      do 230 i = 1, nrn
      y = gamman(i) * ( rugts1 + lamdan(i) )
      z = rugts1 - lamdan(i)
      fg(i) = aa(i) * z
      fh(i) = bb(i) * y
      fj(i) = aa(i) * y
      fk(i) = bb(i) * z
  230 continue
!
      do 250 nst = 1, nq
      do 250 i = 1, nr
      fid(i,nst) = 0.0
  250 continue
!
!  use the source function technique to compute fluxes at each level
!  ,which are eq(55) and (56) of toon et al. (1989)
!
!  form the way down
!
      do 300 nst = 1, nq

      do 280 i = 1, nrn
      fx(i,nst) = exp ( - ( tau(i+nr) - tau(i) ) / ug(nst) )
      xx = lamdan(i) * ug(nst)
      xxm1 = xx - 1.0
      if(abs(xxm1) .lt. 1.0e-6) xx = lamdan(i)*(ug(nst)+ug(nst)*0.01)
      fz1(i,nst) = ( 1.0 - fx(i,nst) * expn(i) ) / ( xx + 1.0 )
      fz2(i,nst) = ( fx(i,nst) - expn(i) ) / ( xx - 1.0 )
      gamman(i) = ug(nst)*beta(i) + 1.0
      if(abs(gamman(i)).lt.1.0e-6)                          &
         gamman(i) = (ug(nst)+ug(nst)*0.01)*beta(i) + 1.0
  280 continue
!
      do 290 j = 1, lev
      jk = j*nr
      jj = jk-nr
      do 290 i = 1, nr
      fid(i+jk,nst) = fid(i+jj,nst)*fx(i+jj,nst) + fj(i+jj)             &
                    * fz1(i+jj,nst) + fk(i+jj)*fz2(i+jj,nst)            &
!org 2              + 1.0/( ug(nst)*beta(i+jj) + 1.0 )
                    + 1.0/gamman(i+jj)                                  &
                    * (alfa(i+jk)-alfa(i+jj)*fx(i+jj,nst)) * fuq2(i+jj)
  290 continue

  300 continue
!
!  fluxes at the surface
!
      do 400 i = 1, nr
      wk1(i) = 0.0
  400 continue
      do 410 nst = 1, nq
      do 410 i = 1, nr
      wk1(i) = wk1(i) + ugwg(nst) * fid(i+jke,nst)
  410 continue
      do 420 i = 1, nr
      wk2(i) = wk1(i)*(1.0-ee(i))*2.0 + 2.*pi*ee(i)*bf(i+jkes)
  420 continue
      do 430 nst = 1, nq
      do 430 i = 1, nr
      fiu(i+jke,nst) = wk2(i)
  430 continue
!
!  from the way up
!
      do 500 nst = 1, nq

      do 480 i = 1, nrn
      gamman(i) = ug(nst)*beta(i) - 1.0
      if(abs(gamman(i)) .lt. 1.0e-6)                   &
         gamman(i)=(ug(nst)+ug(nst)*0.01)*beta(i)-1.0
  480 continue
!
      do 490 j = lev, 1, -1
      jk = j*nr
      jj = jk-nr
      do 490 i = 1, nr
      fiu(i+jj,nst) = fiu(i+jk,nst)*fx(i+jj,nst) + fg(i+jj)             &
                    * fz2(i+jj,nst) + fh(i+jj)*fz1(i+jj,nst)            &
!rog 2              + 1.0/( ug(nst)*beta(i+jj)-1.0 )
                    + 1.0/gamman(i+jj)                                  &
                    * (alfa(i+jk)*fx(i+jj,nst)-alfa(i+jj)) * fuq1(i+jj)
  490 continue

  500 continue
!
      do 600 i = 1, nrn1
      ffu(i) = 0.0
      ffd(i) = 0.0
  600 continue
      do 700 nst = 1, nq
      do 700 i = 1,  nrn1
      ffu(i) = ffu(i) + ugwg(nst) * fiu(i,nst)
      ffd(i) = ffd(i) + ugwg(nst) * fid(i,nst)
  700 continue
!fong08
      return
      end
