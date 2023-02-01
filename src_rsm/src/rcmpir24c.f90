      subroutine rcmpir24c( nr,nx,lev,bf,ee,ww1,ww2,ww,ttau,cld,cldb &
                           ,ffu,ffd )
!
! **********************************************************************
! (1)In this subroutine, the two- and four- stream combination scheme or
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
! (2)The treatment of cloud overlapping for partial clouinesses is
!    included by Fong(1996), which is based on the method proposed
!    by Geleyn and Hollingsworth (1979, contrib. atmos. phys. ).
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
!    c. t. fong and fu (1996, summer)
!    modify to f90 bt C-H Lee in 2015
!
! **********************************************************************
      use paramt
      use fj_pad
      use rank

      implicit  none
!
!fong parameter ( nq = 1 )  ! souce function techn. for two streams
      integer,parameter :: nq = 2 ! ............... for four streams
!fong parameter ( nq = 3 )  ! ..................... for six streams
!
      integer   nr,nx,lev

      real      bf(nx*(lev+2)),ee(nx),ww1(nx*lev,2),ww2(nx*lev,2), &
                ww(nx*lev,2),ttau(nx*(lev+1),2),cld(nx*lev),       &
                cldb(nx*lev,4)
      real      ffu(nx*(lev+1)),ffd(nx*(lev+1))
!
!   local working arrays
!
!     include '../include/paramt.h'
!     include '../include/fj_pad.h'
!     include '../include/rank.h'
!
      real      w1(ilm,2),w(ilm,2),tau(im*(lm+1),2)
      real      ug(nq),wg(nq),ugwg(nq)
      real      lamdan(ilm,2),gamman(ilm,2),caddn(ilm,2),          &
                cminn(ilm,2),caddn0(ilm,2),cminn0(ilm,2),          &
                rl(ilm,2),tl(ilm,2),asir(im),expn(ilm,2)
!fj     3         ,aa(im,4,4)
!fj      dimension amtx(im,lm*4,12),alfa(im*(lm+1)),beta(ilm,2)
      real      amtx(im+npad,lm*4,12),alfa(im*(lm+1)),beta(ilm,2)
      real      fiu(im*(lm+1),2,nq),fid(im*(lm+1),2,nq),           &
                fx(ilm,2,nq),fz1(ilm,2,nq),fz2(ilm,2,nq),          &
                fuq1(ilm,2),fuq2(ilm,2)
      real      wk1(im)
!ibm
      real      tmp1(nr*lev,2),tmp2(nr*lev,2),tmp3(nr*lev,2,nq),   &
                tmp4(nr*lev,2)
!fj
      real      fj_dtau(nr*lev)
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
      integer   nrn,nrn1,lev4,i,m,nst,j,kini,jk,jj,k1,k2
      real      sqrt3,rsqrt3,pi,ugts1
      real      fxx,fw,q1,q2,x,gamma1,gamma2,z,fw1,fw2,xx1,xx2
      real      axx,bxx,cxx,rug,ssfc_f,ssfc_c,rugts1,fp,fq,y1,y
      real      fidm_f,fidm_c,fium_f,fium_c

      nrn  = nr*lev
      nrn1 = nr*(lev+1)
      lev4 = lev*4
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
!  asir : surface albedo for infrared waves
!
      do 8 i = 1, nr
      asir(i) = 1.0 - ee(i)
    8 continue
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
      do 41 m = 1, 2
!ibm
      do i = 1, nr
       tau(i,m) = 0.
      enddo
!
      do i = 1, nrn
         fxx = ww2(i,m) / 5.0
         fw = 1.0 - fxx * ww(i,m)
         w1(i,m) = ( ww1(i,m) - 3.0 * fxx ) / ( 1.0 - fxx )
         w(i,m)  = ( 1.0 - fxx ) * ww(i,m) / fw
         fj_dtau(i) = (ttau(i+nr,m) - ttau(i,m)) * fw
      enddo
!fj
      do i = 1, nrn
         tau(i+nr,m) = tau(i,m) + fj_dtau(i)
      enddo
!--
!ocl simd,swp,norecurrence
      do i = 1, nrn
         tmp2(i,m) = fj_dtau(i)
         fj_dtau(i) = max(0.000001, fj_dtau(i))
         q1 = log( bf(i+nr)/bf(i) )
         if( abs(q1) .lt. 1.0e-10 ) q1 = 1.0e-10
         q2 = 1.0 / fj_dtau(i)
         beta(i,m)= q1 * q2
!--
         if(m.eq.2) w(i,m) = min ( w(i,m), 0.999999 )
         x = w(i,m) * w1(i,m) / 3.0
         gamma1 = 2.0 - w(i,m) - x
         gamma2 = w(i,m) - x
         tmp1(i,m) = (gamma1 + gamma2)*(gamma1 - gamma2)
!
!fj      call vsqrt(lamdan(1,m),tmp1(1,m),nrn)
         lamdan(i,m) = sqrt(tmp1(i,m))
!fj   do i = 1,nrn
!      tmp1(i,m) = -tmp2(i,m)*lamdan(i,m)
         expn(i,m) = -tmp2(i,m)*lamdan(i,m)
         z = lamdan(i,m) * lamdan(i,m) - beta(i,m) * beta(i,m)
         if(abs(z) .lt. 0.000001)then
            beta(i,m) = beta(i,m)+beta(i,m)*0.01
            z = lamdan(i,m) * lamdan(i,m) - beta(i,m) * beta(i,m)
         end if
         tmp1(i,m) = z
!fj      enddo
!     call vexp(expn(1,m),tmp1(1,m),nrn)
!fj      call vexp(expn(1,m),expn(1,m),nrn)
         expn(i,m) = exp(expn(i,m))
!fj      call vrec(tmp1(1,m),tmp1(1,m),nrn)
         tmp1(i,m) = 1.d0/tmp1(i,m)
!ibm
!fj      do 40 i = 1, nrn
      x = w(i,m) * w1(i,m) / 3.0
      gamma1 = 2.0 - w(i,m) - x
      gamma2 = w(i,m) - x
!ibm  dtau = tau(i+nr,m) - tau(i,m)
!ibm  lamdan(i,m) = sqrt((gamma1 + gamma2)*(gamma1 - gamma2))
!vpp  lamdan(i,m) = exp( 0.5*log((gamma1 + gamma2)*(gamma1 - gamma2)) )
      gamman(i,m) = gamma2 / ( gamma1 + lamdan(i,m) )
!ibm  expn(i,m) = exp ( - lamdan(i,m) * dtau )
      fw1 = pi * 2.0 * ( 1.0 - w(i,m) ) * bf(i)
      fw2 = pi * 2.0 * ( 1.0 - w(i,m) ) * bf(i+nr)
      xx1 = ( gamma1 + beta(i,m) ) + gamma2
      xx2 = ( gamma1 - beta(i,m) ) + gamma2
!ibm
!     z = lamdan(i,m) * lamdan(i,m) - beta(i,m) * beta(i,m)
!     if(abs(z) .lt. 0.000001)then
!      beta(i,m) = beta(i,m)+beta(i,m)*0.01
!      z = lamdan(i,m) * lamdan(i,m) - beta(i,m) * beta(i,m)
!     end if
!     caddn0(i,m) = fw1 * xx1 / z
!     cminn0(i,m) = fw1 * xx2 / z
!     caddn(i,m) = fw2 * xx1 / z
!     cminn(i,m) = fw2 * xx2 / z
!ibm
      caddn0(i,m) = fw1 * xx1 * tmp1(i,m)
      cminn0(i,m) = fw1 * xx2 * tmp1(i,m)
      caddn(i,m) = fw2 * xx1 * tmp1(i,m)
      cminn(i,m) = fw2 * xx2 * tmp1(i,m)
      axx = gamman(i,m) * gamman(i,m)
      bxx = expn(i,m) * expn(i,m)
      cxx = axx * bxx - 1.0
      rl(i,m) = gamman(i,m) * ( bxx - 1.0 ) / cxx
      tl(i,m) = expn(i,m) * ( axx - 1.0 ) / cxx
      rl(i,m) = max ( 0., rl(i,m) )
      enddo

!ocl simd,prefetch
      do i = 1, nrn
      z = expn(i,m)*expn(i,m) - gamman(i,m)*gamman(i,m)
      if( abs(z) .lt. 1.0e-10 ) then
         z = 1.0e-10
         tmp4(i,m) = z
      else
         tmp4(i,m) = z
      endif
      enddo
!fj 40   continue
!ibm
!fj      call vrec(tmp4(1,m),tmp4(1,m),nrn)
!ibm
!ocl prefetch
      do i = 1,nrn
         tmp4(i,m) = 1.d0/tmp4(i,m)
      enddo

   41 continue
!fong02
!ibm
!ocl prefetch
      do nst = 1, nq
       rug = 1./ug(nst)
      do i = 1, nrn*2
       tmp3(i,1,nst) = -tmp2(i,1)*rug
       tmp3(i,1,nst) = exp(tmp3(i,1,nst))
      enddo
      enddo
!fj      call vexp(tmp3(1,1,1),tmp3(1,1,1),nrn*2*nq)
!ibm

!
!     do 45 i = 1, nx*lev4*12
!     amtx(i,1,1) = 0.
!  45 continue
!fong03
!fj
!fj j = 1
!fj
      j = 1
      kini = (j-1) * 4
      jk = j*nr
      jj = jk-nr
!ocl prefetch
      do i = 1, nr
         amtx(i,kini+1,3)=0.
         amtx(i,kini+1,4)=0.
         amtx(i,kini+1,5)=0.
         amtx(i,kini+1,7)=0.
         amtx(i,kini+1,8)=0.
         amtx(i,kini+1,9)=0.
         amtx(i,kini+2,3)=0.
         amtx(i,kini+2,4)=0.
         amtx(i,kini+2,5)=0.
         amtx(i,kini+2,7)=0.
         amtx(i,kini+2,8)=0.
         amtx(i,kini+3,3)=0.
         amtx(i,kini+3,4)=0.
         amtx(i,kini+3,5)=0.
         amtx(i,kini+3,7)=0.
         amtx(i,kini+4,3)=0.
         amtx(i,kini+4,4)=0.
         amtx(i,kini+4,5)=0.
         amtx(i,kini+4,9)=0.
         amtx(i,kini+1,6) = 1.0
         amtx(i,kini+1,10) = - tl(i+jj,1) * cldb(i+jj,2)
         amtx(i,kini+1,11) = - tl(i+jj,1) * ( 1.0 - cldb(i+jj,4) )
         amtx(i,kini+2,6) = 1.0
         amtx(i,kini+2,9) =  - tl(i+jj,2) * ( 1.0 - cldb(i+jj,2) )
         amtx(i,kini+2,10) =  - tl(i+jj,2) * cldb(i+jj,4)
         amtx(i,kini+3,6) = 1.0
         amtx(i,kini+3,8) = - rl(i+jj,1) * cldb(i+jj,2)
         amtx(i,kini+3,9) = - rl(i+jj,1) * ( 1.0 - cldb(i+jj,4) )
         amtx(i,kini+4,6) = 1.0
         amtx(i,kini+4,7) =  - rl(i+jj,2) * ( 1.0 - cldb(i+jj,2) )
         amtx(i,kini+4,8) =  - rl(i+jj,2) * cldb(i+jj,4)
         amtx(i,kini+1,12) = ( -cminn0(i+jj,1)*rl(i+jj,1)-caddn(i+jj,1)* &
              tl(i+jj,1)+caddn0(i+jj,1) ) * ( 1.0-cld(i+jj) )
         amtx(i,kini+2,12) = ( -cminn0(i+jj,2)*rl(i+jj,2)-caddn(i+jj,2)* &
              tl(i+jj,2)+caddn0(i+jj,2) ) * cld(i+jj)
         amtx(i,kini+3,12) = ( -cminn0(i+jj,1)*tl(i+jj,1)-caddn(i+jj,1)* &
              rl(i+jj,1)+cminn(i+jj,1) ) * ( 1.0-cld(i+jj) )
         amtx(i,kini+4,12) = ( -cminn0(i+jj,2)*tl(i+jj,2)-caddn(i+jj,2)* &
              rl(i+jj,2)+cminn(i+jj,2) ) * cld(i+jj)
      enddo

!fj
!fj j = 2, lev-1
!fj
      do j = 2,lev-1

         kini = (j-1) * 4
         jk = j*nr
         jj = jk-nr
!
!
         do i = 1,nr
            amtx(i,kini+1,3)=0.
            amtx(i,kini+1,7)=0.
            amtx(i,kini+1,8)=0.
            amtx(i,kini+1,9)=0.
            amtx(i,kini+2,5)=0.
            amtx(i,kini+2,7)=0.
            amtx(i,kini+2,8)=0.
            amtx(i,kini+3,4)=0.
            amtx(i,kini+3,5)=0.
            amtx(i,kini+3,7)=0.
            amtx(i,kini+4,3)=0.
            amtx(i,kini+4,4)=0.
            amtx(i,kini+4,5)=0.
            amtx(i,kini+4,9)=0.
         enddo
!ocl prefetch
         do i = 1,nr
            amtx(i,kini+1,4)  = - rl(i+jj,1) * cldb(i+jj,1)
            amtx(i,kini+1,5)  = - rl(i+jj,1) * ( 1.0 - cldb(i+jj,3) )
            amtx(i,kini+1,10) = - tl(i+jj,1) * cldb(i+jj,2)
            amtx(i,kini+1,11) = - tl(i+jj,1) * ( 1.0 - cldb(i+jj,4) )
            amtx(i,kini+2,3)  = - rl(i+jj,2) * ( 1.0 - cldb(i+jj,1) )
            amtx(i,kini+2,4)  = - rl(i+jj,2) * cldb(i+jj,3)
            amtx(i,kini+2,9)  = - tl(i+jj,2) * ( 1.0 - cldb(i+jj,2) )
            amtx(i,kini+2,10) = - tl(i+jj,2) * cldb(i+jj,4)
            amtx(i,kini+3,2)  = - tl(i+jj,1) * cldb(i+jj,1)
            amtx(i,kini+3,3)  = - tl(i+jj,1) * ( 1.0 - cldb(i+jj,3) )
            amtx(i,kini+3,8)  = - rl(i+jj,1) * cldb(i+jj,2)
            amtx(i,kini+3,9)  = - rl(i+jj,1) * ( 1.0 - cldb(i+jj,4) )
            amtx(i,kini+4,1)  = - tl(i+jj,2) * ( 1.0 - cldb(i+jj,1) )
            amtx(i,kini+4,2)  = - tl(i+jj,2) * cldb(i+jj,3)
            amtx(i,kini+4,7)  = - rl(i+jj,2) * ( 1.0 - cldb(i+jj,2) )
            amtx(i,kini+4,8)  = - rl(i+jj,2) * cldb(i+jj,4)
         enddo
!ocl prefetch
         do i = 1,nr
            amtx(i,kini+1,12) = ( -cminn0(i+jj,1)*rl(i+jj,1)-caddn(i+jj,1)*   &
                 tl(i+jj,1)+caddn0(i+jj,1) ) * ( 1.0-cld(i+jj) )
            amtx(i,kini+2,12) = ( -cminn0(i+jj,2)*rl(i+jj,2)-caddn(i+jj,2)*   &
                 tl(i+jj,2)+caddn0(i+jj,2) ) * cld(i+jj)
            amtx(i,kini+3,12) = ( -cminn0(i+jj,1)*tl(i+jj,1)-caddn(i+jj,1)*   &
                 rl(i+jj,1)+cminn(i+jj,1) ) * ( 1.0-cld(i+jj) )
            amtx(i,kini+4,12) = ( -cminn0(i+jj,2)*tl(i+jj,2)-caddn(i+jj,2)*   &
                rl(i+jj,2)+cminn(i+jj,2) ) * cld(i+jj)
         enddo

         do i = 1,nr
            amtx(i,kini+1,6) = 1.0
            amtx(i,kini+2,6) = 1.0
            amtx(i,kini+3,6) = 1.0
            amtx(i,kini+4,6) = 1.0
         enddo
      enddo

!fj
!fj j = lev
!fj

      j = lev
      kini = (j-1) * 4
      jk = j*nr
      jj = jk-nr
!
!
!ocl prefetch
      do i = 1, nr
         amtx(i,kini+1,3)=0.
         amtx(i,kini+1,7)=0.
         amtx(i,kini+1,8)=0.
         amtx(i,kini+1,9)=0.
         amtx(i,kini+2,5)=0.
         amtx(i,kini+2,7)=0.
         amtx(i,kini+2,8)=0.
         amtx(i,kini+2,9)=0.
         amtx(i,kini+3,4)=0.
         amtx(i,kini+3,5)=0.
         amtx(i,kini+3,7)=0.
         amtx(i,kini+3,8)=0.
         amtx(i,kini+3,9)=0.
         amtx(i,kini+4,3)=0.
         amtx(i,kini+4,4)=0.
         amtx(i,kini+4,5)=0.
         amtx(i,kini+4,7)=0.
         amtx(i,kini+4,8)=0.
         amtx(i,kini+4,9)=0.
         amtx(i,kini+1,6) = 1.0
         amtx(i,kini+1,4) = - rl(i+jj,1) * cldb(i+jj,1)
         amtx(i,kini+1,5) = - rl(i+jj,1) * ( 1.0 - cldb(i+jj,3) )
         amtx(i,kini+1,8) = amtx(i,kini+1,8) +                      &
              asir(i)*(- tl(i+jj,1) * cldb(i+jj,2))
         amtx(i,kini+1,9) = amtx(i,kini+1,9) +                      &
              asir(i)*(- tl(i+jj,1) * ( 1.0 - cldb(i+jj,4) ))
         amtx(i,kini+2,6) = 1.0
         amtx(i,kini+2,3) = - rl(i+jj,2) * ( 1.0 - cldb(i+jj,1) )
         amtx(i,kini+2,4) = - rl(i+jj,2) * cldb(i+jj,3)
         amtx(i,kini+2,7) = amtx(i,kini+2,7) +                      &
              asir(i)*(- tl(i+jj,2) * ( 1.0 - cldb(i+jj,2) ))
         amtx(i,kini+2,8) = amtx(i,kini+2,8) +                      &
              asir(i)*(- tl(i+jj,2) * cldb(i+jj,4))
         amtx(i,kini+3,6) = 1.0
         amtx(i,kini+3,2) = - tl(i+jj,1) * cldb(i+jj,1)
         amtx(i,kini+3,3) = - tl(i+jj,1) * ( 1.0 - cldb(i+jj,3) )
         amtx(i,kini+3,6) = amtx(i,kini+3,6) +                      &
              asir(i)*(- rl(i+jj,1) * cldb(i+jj,2))
         amtx(i,kini+3,7) = amtx(i,kini+3,7) +                      &
              asir(i)*(- rl(i+jj,1) * ( 1.0 - cldb(i+jj,4) ))
         amtx(i,kini+4,6) = 1.0
         amtx(i,kini+4,1) = - tl(i+jj,2) * ( 1.0 - cldb(i+jj,1) )
         amtx(i,kini+4,2) = - tl(i+jj,2) * cldb(i+jj,3)
         amtx(i,kini+4,5) = amtx(i,kini+4,5) +                      &
              asir(i)*(- rl(i+jj,2) * ( 1.0 - cldb(i+jj,2) ))
         amtx(i,kini+4,6) = amtx(i,kini+4,6) +                      &
              asir(i)*(- rl(i+jj,2) * cldb(i+jj,4))
      enddo

!ocl prefetch
      do i = 1, nr
         amtx(i,kini+1,12) = ( -cminn0(i+jj,1)*rl(i+jj,1)-caddn(i+jj,1)* &
              tl(i+jj,1)+caddn0(i+jj,1) ) * ( 1.0-cld(i+jj) )
         amtx(i,kini+2,12) = ( -cminn0(i+jj,2)*rl(i+jj,2)-caddn(i+jj,2)* &
              tl(i+jj,2)+caddn0(i+jj,2) ) * cld(i+jj)
         amtx(i,kini+3,12) = ( -cminn0(i+jj,1)*tl(i+jj,1)-caddn(i+jj,1)* &
              rl(i+jj,1)+cminn(i+jj,1) ) * ( 1.0-cld(i+jj) )
         amtx(i,kini+4,12) = ( -cminn0(i+jj,2)*tl(i+jj,2)-caddn(i+jj,2)* &
              rl(i+jj,2)+cminn(i+jj,2) ) * cld(i+jj)
        enddo

!ocl norecurrence
        do i = 1, nr
           ssfc_f = ee(i) * pi * bf(i+nrn1) * ( 1.0 - cld(i+jj) )
           ssfc_c = ee(i) * pi * bf(i+nrn1) * cld(i+jj) 	
           amtx(i,kini+1,12) = amtx(i,kini+1,12) - ( (- tl(i+jj,1) * cldb(i+jj,2))*ssfc_f +             &
                (- tl(i+jj,1) * ( 1.0 - cldb(i+jj,4) ))*ssfc_c )
           amtx(i,kini+2,12) = amtx(i,kini+2,12) - ( (- tl(i+jj,2) * ( 1.0 - cldb(i+jj,2) ))*ssfc_f +   &
                (- tl(i+jj,2) * cldb(i+jj,4))*ssfc_c )
           amtx(i,kini+3,12) = amtx(i,kini+3,12) - ( (- rl(i+jj,1) * cldb(i+jj,2))*ssfc_f +             &
                (- rl(i+jj,1) * ( 1.0 - cldb(i+jj,4) ))*ssfc_c )
           amtx(i,kini+4,12) = amtx(i,kini+4,12) - ( (- rl(i+jj,2) * ( 1.0 - cldb(i+jj,2) ))*ssfc_f +   &
                (- rl(i+jj,2) * cldb(i+jj,4))*ssfc_c )
        enddo
!fong03
!
!  solve matrix
!
      call mxslvc ( amtx, nx, nr, lev4 )
!
!fong04
!
!fong07
!
!  rl : Y1 + Y2 as in the table 3 of toon et al. (1989)
!  tl : Y1 - Y2    .......
!ibm
!     do 210 i = 1, nrn
!     fz1(i,1,1) = 1.0 - cld(i)
!     fz1(i,2,1) = cld(i)
! 210 continue
!ibm
      rugts1 = 1.0 / ugts1
!
      do 250 m = 1, 2
      do 250 j = 1, lev
      jk = j*nr
      jj = jk-nr
      k1 = (j-1)*4
      k2 = k1 + 2
      do 250 i = 1, nr
      fz1(i+jj,m,1) = (2-m)*(1.0 - cld(i+jj))+(m-1)*cld(i+jj)
      fp = amtx(i,k1+m,12) - caddn0(i+jj,m) * fz1(i+jj,m,1)
      fq = amtx(i,k2+m,12) - cminn(i+jj,m) * fz1(i+jj,m,1)
!ibm  z = expn(i+jj,m)*expn(i+jj,m) - gamman(i+jj,m)*gamman(i+jj,m)
!ibm  if( abs(z) .lt. 1.0e-10 ) z = 1.0e-10
!ibm  rl(i+jj,m) = ( fp * expn(i+jj,m) - fq * gamman(i+jj,m) ) / z
!ibm  tl(i+jj,m) = ( fq * expn(i+jj,m) - fp * gamman(i+jj,m) ) / z
      rl(i+jj,m) = ( fp * expn(i+jj,m) - fq * gamman(i+jj,m) )*tmp4(i+jj,m)
      tl(i+jj,m) = ( fq * expn(i+jj,m) - fp * gamman(i+jj,m) )*tmp4(i+jj,m)
! 250 continue
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
!  caddn0, cminn0, caddn, camin : parameter G, H, J, K of the table 3 i
!  respectively
!  fuq1 * (2*pi*bf) : the rest of terms for upward source function
!  fuq2 * (2*pi*bf) : the rest of terms for downward source function
!
!     do 310 m = 1, 2
!     do 310 i = 1, nrn
!ibm  x = 2.0 * ( 1.0 - w(i+jj,m) ) * w(i+jj,m) / ( lamdan(i+jj,m) *
!ibm 1    lamdan(i+jj,m) - beta(i+jj,m) * beta(i+jj,m) )
      x = 2.0 * ( 1.0 - w(i+jj,m) ) * w(i+jj,m) * tmp1(i+jj,m)
      y1 = w1(i+jj,m) / 3.0
      y = 2.0 * ( 1.0 - w(i+jj,m) * y1 )
      z = -y1 * beta(i+jj,m)
      fuq1(i+jj,m) = ( x * ( y - z ) + 1.0 - w(i+jj,m) ) * fz1(i+jj,m,1)
      fuq2(i+jj,m) = ( x * ( y + z ) + 1.0 - w(i+jj,m) ) * fz1(i+jj,m,1)
! 310 continue
!
!     rugts1 = 1.0 / ugts1
!     do 330 m = 1, 2
!     do 330 i = 1, nrn
      y = gamman(i+jj,m) * ( rugts1 + lamdan(i+jj,m) )
      z = rugts1 - lamdan(i+jj,m)
      caddn0(i+jj,m) = rl(i+jj,m) * z
      cminn0(i+jj,m) = tl(i+jj,m) * y
      caddn(i+jj,m)  = rl(i+jj,m) * y
      cminn(i+jj,m)  = tl(i+jj,m) * z
! 330 continue
  250 continue
!
      do 320 i = 1, nrn1
      alfa(i) = 2. * pi * bf(i)
  320 continue
!
      do 350 nst = 1, nq
      do 350 i = 1, nr
      fid(i,1,nst) = 0.0
      fid(i,2,nst) = 0.0
  350 continue
!
!  use the source function technique to compute fluxes at each level
!  ,which are eq(55) and (56) of toon et al. (1989)
!
!  from the way down
!
!ocl prefetch
      do 400 nst = 1, nq
!ibm
      do i = 1, nrn
       tmp1(i,1) = lamdan(i,1) * ug(nst) + 1.0
       tmp1(i,2) = lamdan(i,2) * ug(nst) + 1.0
       tmp1(i,1) = 1.d0/tmp1(i,1)
       tmp1(i,2) = 1.d0/tmp1(i,2)
      enddo
!fj
!ocl simd,norecurrence
      do i = 1, nrn
       tmp2(i,1) = lamdan(i,1) * ug(nst) - 1.0
       tmp2(i,2) = lamdan(i,2) * ug(nst) - 1.0
       if(abs(tmp2(i,1)).lt.1.0e-6)   &
              tmp2(i,1)=lamdan(i,1)*(ug(nst)+ug(nst)*0.01) - 1.0
       if(abs(tmp2(i,2)).lt.1.0e-6)   &
              tmp2(i,2)=lamdan(i,2)*(ug(nst)+ug(nst)*0.01) - 1.0
       tmp2(i,1) = 1.d0/tmp2(i,1)
       tmp2(i,2) = 1.d0/tmp2(i,2)
      enddo
!
!fj      call vrec(tmp1(1,1),tmp1(1,1),nrn*2)
!fj      call vrec(tmp2(1,1),tmp2(1,1),nrn*2)

!ibm
      do 380 i = 1, nrn
!ibm  fx(i,1,nst) = exp ( - (tau(i+nr,1) - tau(i,1)) / ug(nst) )
!ibm  fx(i,2,nst) = exp ( - (tau(i+nr,2) - tau(i,2)) / ug(nst) )
!ibm  fx ---> tmp3
!fong xx1 = lamdan(i,1) * ug(nst) + 1.0e-10   ! for two streams
!fong xx2 = lamdan(i,2) * ug(nst) + 1.0e-10   ! for two streams
!ibm  xx1 = lamdan(i,1) * ug(nst)
!ibm  xx2 = lamdan(i,2) * ug(nst)
      fz1(i,1,nst) = ( 1.0 - tmp3(i,1,nst) * expn(i,1) )  &
!ibm 1                  / ( xx1 + 1.0 )
                        * tmp1(i,1)
      fz1(i,2,nst) = ( 1.0 - tmp3(i,2,nst) * expn(i,2) )  &
!ibm 1                  / ( xx2 + 1.0 )
                        * tmp1(i,2)
!ibm  xx1m1 = xx1 - 1.
!ibm  xx2m1 = xx2 - 1.
!ibm  if(abs(xx1m1).lt.1.0e-6)xx1=lamdan(i,1)*(ug(nst)+ug(nst)*0.01)
!ibm  if(abs(xx2m1).lt.1.0e-6)xx2=lamdan(i,2)*(ug(nst)+ug(nst)*0.01)
!ibm  fz2(i,1,nst) = ( tmp3(i,1,nst) - expn(i,1) )/( xx1 - 1.0 )
!ibm  fz2(i,2,nst) = ( tmp3(i,2,nst) - expn(i,2) )/( xx2 - 1.0 )
      fz2(i,1,nst) = ( tmp3(i,1,nst) - expn(i,1) )*tmp2(i,1)
      fz2(i,2,nst) = ( tmp3(i,2,nst) - expn(i,2) )*tmp2(i,2)
 380  continue
!fj
!ocl simd,norecurrence
      do i = 1,nrn
      gamman(i,1) = ug(nst)*beta(i,1)+1.0
      gamman(i,2) = ug(nst)*beta(i,2)+1.0
      if(abs(gamman(i,1)).lt.1.0e-6)                        &
         gamman(i,1) = (ug(nst)+ug(nst)*0.01)*beta(i,1)+1.0
      if(abs(gamman(i,2)).lt.1.0e-6)                        &
         gamman(i,2) = (ug(nst)+ug(nst)*0.01)*beta(i,2)+1.0
      tmp1(i,1) = 1.d0/gamman(i,1)
      tmp1(i,2) = 1.d0/gamman(i,2)
      enddo

!fj  380 continue
!ibm
!fj      call vrec(tmp1(1,1),gamman(1,1),nrn)
!fj      call vrec(tmp1(1,2),gamman(1,2),nrn)
!ibm
!
!ocl prefetch
      do 390 j = 1, lev
      jk = j*nr
      jj = jk-nr
      do 390 i = 1, nr
      fidm_f = fid(i+jj,1,nst)*cldb(i+jj,1) + fid(i+jj,2,nst)         &
               *(1.0-cldb(i+jj,3))
      fidm_c = fid(i+jj,1,nst)*(1.0-cldb(i+jj,1)) + fid(i+jj,2,nst)   &
               *cldb(i+jj,3)
      fid(i+jk,1,nst) = fidm_f*tmp3(i+jj,1,nst) + caddn(i+jj,1)       &
                    * fz1(i+jj,1,nst) + cminn(i+jj,1)*fz2(i+jj,1,nst) &
!org 2              + 1.0/( ug(nst)*beta(i+jj,1) + 1.0 )
!ibm 2              + 1.0/gamman(i+jj,1)
                    + tmp1(i+jj,1)                                    &
                    * (alfa(i+jk)-alfa(i+jj)*tmp3(i+jj,1,nst))        &
                    * fuq2(i+jj,1)
      fid(i+jk,2,nst) = fidm_c*tmp3(i+jj,2,nst) + caddn(i+jj,2)       &
                    * fz1(i+jj,2,nst) + cminn(i+jj,2)*fz2(i+jj,2,nst) &
!org 2              + 1.0/( ug(nst)*beta(i+jj,2) + 1.0 )
!ibm 2              + 1.0/gamman(i+jj,2)
                    + tmp1(i+jj,2)                                    &
                    * (alfa(i+jk)-alfa(i+jj)*tmp3(i+jj,2,nst))        &
                    * fuq2(i+jj,2)
  390 continue

  400 continue
!
!  fluxes at the surface
!
      do 405 i = 1, nr
      wk1(i) = 0.0
  405 continue
      do 410 nst = 1, nq
      do 410 i = 1, nr
      wk1(i) = wk1(i) + ugwg(nst)*(fid(i+nrn,1,nst)+fid(i+nrn,2,nst))
  410 continue
      do 420 i = 1, nr
      wk1(i) = wk1(i)*asir(i)*2.0 + 2.*pi*ee(i)*bf(i+nrn1)
  420 continue
      do 430 nst = 1, nq
      do 430 i = 1, nr
      fiu(i+nrn,1,nst) = wk1(i)
      fiu(i+nrn,2,nst) = 0.
  430 continue
!
!  from the way up
!
!ocl prefetch
      do 500 nst = 1, nq

!ocl simd,norecurrence
      do 480 i = 1, nrn
      gamman(i,1) = ug(nst)*beta(i,1)-1.0
      gamman(i,2) = ug(nst)*beta(i,2)-1.0
      if(abs(gamman(i,1)).lt.1.0e-6)                           &
         gamman(i,1) = (ug(nst)+ug(nst)*0.01)*beta(i,1)-1.0
      if(abs(gamman(i,2)).lt.1.0e-6)                           &
         gamman(i,2) = (ug(nst)+ug(nst)*0.01)*beta(i,2)-1.0
      tmp1(i,1) = 1.d0/gamman(i,1)
      tmp1(i,2) = 1.d0/gamman(i,2)
 480  continue

!ibm
!fj      call vrec(tmp1(1,1),gamman(1,1),nrn)
!fj      call vrec(tmp1(1,2),gamman(1,2),nrn)
!ibm
      do 490 j = lev, 1, -1
      jk = j*nr
      jj = jk-nr
      do 490 i = 1, nr
      fium_f = fiu(i+jk,1,nst)*cldb(i+jj,2) + fiu(i+jk,2,nst)*        &
               (1.0-cldb(i+jj,4))
      fium_c = fiu(i+jk,1,nst)*(1.0-cldb(i+jj,2)) + fiu(i+jk,2,nst)*  &
               cldb(i+jj,4)
      fiu(i+jj,1,nst) = fium_f*tmp3(i+jj,1,nst) + caddn0(i+jj,1)      &
                    * fz2(i+jj,1,nst) + cminn0(i+jj,1)*fz1(i+jj,1,nst)&
!org 2              + 1.0/( ug(nst)*beta(i+jj,1)-1.0 )
!ibm 2              + 1.0/gamman(i+jj,1)
                    + tmp1(i+jj,1)                                    &
                    * (alfa(i+jk)*tmp3(i+jj,1,nst)-alfa(i+jj))        &
                    * fuq1(i+jj,1)
      fiu(i+jj,2,nst) = fium_c*tmp3(i+jj,2,nst) + caddn0(i+jj,2)      &
                    * fz2(i+jj,2,nst) + cminn0(i+jj,2)*fz1(i+jj,2,nst)&
!org 2              + 1.0/( ug(nst)*beta(i+jj,2)-1.0 )
!ibm 2              + 1.0/gamman(i+jj,2)
                    + tmp1(i+jj,2)                                    &
                    * (alfa(i+jk)*tmp3(i+jj,2,nst)-alfa(i+jj))        &
                    * fuq1(i+jj,2)
  490 continue

  500 continue
!
      do 600 i = 1, nrn1
      ffu(i) = 0.0
      ffd(i) = 0.0
  600 continue
!ocl prefetch
      do 700 nst = 1, nq
      do 700 i = 1,  nrn1
      ffu(i) = ffu(i) + ugwg(nst) * ( fiu(i,1,nst) + fiu(i,2,nst) )
      ffd(i) = ffd(i) + ugwg(nst) * ( fid(i,1,nst) + fid(i,2,nst) )
  700 continue
!fong08
      return
      end
