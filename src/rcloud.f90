      subroutine rcloud ( nxj,im,lm,ls,ptop,ptopc,pst,pp,tt,rh,cufrac,plcl     &
                        , cumtop,clds,cvclds,rhx,sdpbl,ncld,qtr                &
                        , clwc,ciwc,cre,cde,cldmaxv ) 
!######################################################################c
! a. function:                                                         c
!    (1)compute cloud fraction and cloud optical properties
!    (2)stable cloud fraction is diagnosed by input temperature,
!    specific humidity, and cloud water content, but convective cloud
!    fraction by convective precip, lcl, and cloud heights.
!
! b. common block specification:                                       c
!    none                                                              c
!                                                                      c
! c. parameter specification:                                          c
!    none                                                              c
!                                                                      c
! d. input/output variables:                                           c
!   (1) input variables:                                               c
!     im    : model total number of horizontal grids                   c
!     lm    : model total number of vertical layers                    c
!     ls    : number of layers at model top where no clouds allowed    c
!     ptop  : pressure at top of model                                 c
!     pst    : (terrain pressure -ptop)                                c
!     pp    : model level pressures                                    c
!     tt    : model level temperatures                                 c
!     rh    : model level relative humidity                            c
!     cufrac: convective cloud fractions                               c
!     plcl  : pressure of lifting condensation level (lcl)             c
!     cumtop: pressure at top of deepest convective clouds             c
!     sdpbl : vertical velocity (omega) at the near 800hPa level
!                  to indicate the convergence in the pbl
!     rhx   : work array (reference rh for full grid stable clouds     c
!     qtr   : cloud water content
!                                                                      c
!   (2) output variables:                                              c
!     clds  : stable cloud fractions                                   c
!     cvclds: convective cloud fractions                               c
!     clwc/ciwc : cloud liquid/ice water content for radiation use(g/m3)
!     cre/cde   : cloud optical effective size (um)
!                                                                      c
!   first version : march, 2007
!                                                                      c
!   auther:                                                           c
!    c. t. fong
!    modif to f90 by C-H Lee and sort by River Chen in 2015
!                                                                      c
!######################################################################c
!
      use const, only: RTYPE
!
      implicit  none

      integer   nxj,im,lm,ls,ncld
      real      ptop,ptopc

      real      pst(im),pp(im,lm),rh(im,lm),cufrac(im),                  &
                plcl(im),cumtop(im),clds(im,lm),cvclds(im,lm),rhx(im,lm),&
                qtr(im,lm),                                              &
                clwc(im,lm),ciwc(im,lm),cre(im,lm),cde(im,lm)
      real(kind=RTYPE) sdpbl(im),tt(im,lm)
!
      real      rhc(im,lm)
      real      okk,cre_min,cre_max
      integer   kk
!
      real      qlwc(9),qlwcv(9),qiwc(15),qiwcv(15),deg(8)

      integer i,ls1,ls2,lmm1,k,nn,nwc,nic
      real    cc1,cc2,sq3,cldmaxv,p850,anvil,pps,rrs,cld,fff,tem1,tem2,rgas,t273
      real    t253,qtrc,fcld,ratio,aln10,ttc,cvf,cvf1,xxx,rlwcv,riwcv,rrr,ctt,rlwc,riwc

!
!  create diagnostic formulation of liquid water content as functions
!  of temperature based on the moist thermodynamic relationship proposed
!  by betts and harshvardhan (1987), where assumes lwc proportional to
!  the change of saturation mixing ratio by lifting saturated air at
!  ambient temperature vertically through a small distance(dz) along the
!  moist adiabat.
!     dz = 120(meter) for stratus water clouds
!        = 180(meter) for convective  water clouds
!        = 120(meter) for stratus ice clouds
!        = 180 + 8./70. * (273.-tc)**2 for convective ice clouds
!     tc : cloud's temperature
!  ** the temperature range for the following data set (all with
!     a temperature interval equal to 5 K)
!     qlwc(liquid water content) : 20C to -20C
!     qiwc(ice water content, stratus) : 0C to -70C
!     qiwcv(ice water content, convective) : 0C to -70C
!
      data qlwc /-0.6067,-0.6789,-0.7592,-0.8416,-0.9317,-1.0296,-1.1374       &
               , -1.2552,-1.3836/
      data qlwcv /-0.6067,-0.6789,-0.7592,-0.8416,-0.9317,-1.0296,-1.1374      &
               , -1.2552,-1.3836/
!cc   data qlwcv/-0.4333,-0.5053,-0.5860,-0.6702,-0.7597,-0.8580,-0.9668
!cc  1         , -1.0853,-1.2153/
      data qiwc /-0.9493,-1.0322,-1.1509,-1.2824,-1.4274,-1.5934,-1.7698       &
               , -1.9591,-2.1711,-2.3855,-2.6158,-2.8539,-3.1231,-3.3834       &
               , -3.6670/
      data qiwcv /-0.9493,-1.0322,-1.1509,-1.2824,-1.4274,-1.5934,-1.7698      &
               , -1.9591,-2.1711,-2.3855,-2.6158,-2.8539,-3.1231,-3.3834       &
               , -3.6670/
!cc   data qiwcv/-0.7753,-0.8535,-0.9542,-1.0579,-1.1685,-1.2918,-1.4262
!cc  1         , -1.5714,-1.7344,-1.9080,-2.1012,-2.3038,-2.5388,-2.7737
!cc  2         , -3.0379/
!
!      data deg/ 72.4, 57.6, 82.7, 55.2, 28.8, 27.9, 21.0, 36.4 /
      data deg/ 72.4, 57.6, 82.7, 55.2, 38.8, 32.9, 26.0, 42.4 /   !a4
!      data deg/ 72.4, 57.6, 82.7, 55.2, 38.8, 32.9, 29.7, 49.1 /   !a5
!
!cc   data anvil/0.2/
      data anvil/0.25/
!
      do k = 1, lm
      do i = 1, im
       rhc(i,k)=clds(i,k)
       clds(i,k) = 0.
      enddo
      enddo
!
      cc1= 0.25
      cc2= 1.0 -2.0*cc1
!
!  compute vertically filtered relative humidity and store in array rh
!
!  use rhx and cvclds as work space
!
      ls1= ls+1
      ls2= ls+2
      lmm1 = lm-1
!xx      do 100 k=ls2,lmm1
!xx      do 100 i=1,im
!xx      cvclds(i,k)= cc1*(rh(i,k-1)+rh(i,k+1))
!xx  100 continue
!
!fong
!     do 110 k=ls2,lmm1
!     do 110 i=1,im
!     cvclds(i,k)= cc1*(rh(i,k-1)+rh(i,k+1))
!     rhx(i,k)= cc2*rh(i,k)+cvclds(i,k)
! 110 continue
!     do 120 i=1,im
!     rhx(i,ls1)= rh(i,ls1)
!xx      rhx(i,lm)= rh(i,lm)
!     rhx(i,lm)= (rh(i,lm)+rh(i,lmm1))*0.5
! 120 continue
!fong
      do k = 1, lm
      do i = 1, nxj
       rhx(i,k) = min(0.999, max(0.,rh(i,k)))
      enddo
      enddo
!fong
!
!  stable cloud fractions
!
      sq3= sqrt(3.0)

!
      if( ncld .eq. 1 )then

!      cldmaxv = 0.6   ! max one-layer cloudiness

      do 130 k = ls1, lm
      do 130 i = 1, nxj
      pps = pp(i,k) / (pst(i)+ptop)
      rrs = 1.0 + 2.0*(pps**2 -pps) + sq3*pps*(1.0 -3.0*pps+2.0*pps**2)
      rrs = max( 0.7, min (rrs, 0.9995) )
      cld = (rhx(i,k)-rrs)/(1.0-rrs)
      cld = min (1.0, max(0.0, cld))
      if( pp(i,k) .gt. ptopc ) then
!cc   clds(i,k)= cld*cld
      clds(i,k)= min(cldmaxv,cld*cld)
      end if
  130 continue
!
!  reduce low cloudiness at ocean grids without significant upward
!  motion
!
      do 135 k = lm/2, lmm1
      do 135 i = 1, nxj
!lim
      if(pp(i,k) .gt. 850. .and. sdpbl(i) .gt. -0.000001 ) then
        clds(i,k) = clds(i,k) * 0.10
      end if
!lim
  135 continue

      else     ! (ncld.eq.2)

!      cldmaxv = 0.6   ! max one-layer cloudiness
!      cldmaxv = 0. 85  ! max one-layer cloudiness
      p850 = 850.

      do 131 k = ls1, lm
      do 131 i = 1, nxj
! Xu and Randall (1996)
        fff=-100.*qtr(i,k)/(((1.-rhx(i,k))*rhc(i,k))**0.49)
        cld=rhx(i,k)**0.25*(1.-exp(fff))
! original (fu and Liou)
!         fff = 1000.*qtr(i,k)/(1.-rhx(i,k))
!         cld = rhx(i,k)*(1.-exp(-fff))
!
      if(cld .lt. 0.05) cld=0.    ! avoid too small cloudiness
      if( pp(i,k) .gt. ptopc ) then
!       if(rhx(i,k) .gt. 0.5 .and. qtr(i,k).gt.10e-6)   &
         clds(i,k)= min(cldmaxv,max(0.,cld))
      end if
!  limit cloudiness at lowest two layers
      if( (pp(i,k).gt.p850) .and. (k.ge.lm-1) )then
         clds(i,k) = min (0.2, clds(i,k))
      endif
  131 continue

      endif

!
!  convective cloud fractions, function of convective precip, lifting
!  condensation level, level of vanishing cloud buoynacy, and ice cloud
!  enhancment at low temperatures.  stored in array cvclds.
!
      tem1= 0.1
      tem2= 0.8
      do 140 i = 1, nxj
      if ( cufrac(i) .lt. tem1 )  cufrac(i) = 0.0
      if ( cufrac(i) .gt. tem2 )  cufrac(i) = tem2 + 0.01
 140  continue
!
!     do i = 1, im*lm
!      cvclds(i,1) = 0.
!     enddo
       cvclds = 0.
!
      do i=1,nxj
        kk=0
        do k=lmm1,ls,-1
          if ( (pp(i,k).le.plcl(i)) .and. (pp(i,k).ge.cumtop(i)) ) then
            if(cufrac(i).gt.0.3 .and. pp(i,k).le.400.)then
              kk=kk+1
              okk=1.-0.225*(kk-1)
              cvclds(i,k) = cufrac(i)**okk
            else
              cvclds(i,k) = cufrac(i)
            endif
            cvclds(i,k) = min(tem2,cvclds(i,k))
          endif
        enddo
      enddo
!
!      do 150 i=1,nxj
!      if ( (pp(i,k).le.plcl(i)) .and. (pp(i,k).ge.cumtop(i)) )  &
!            cvclds(i,k) = cufrac(i)
!      if ( (cvclds(i,k).gt.0.0) .and. (tt(i,k).lt.233.0) )      &
!!ccc *      cvclds(i,k) = cvclds(i,k) + anvil
!            cvclds(i,k) = cvclds(i,k) + cvclds(i,k)*anvil
!  150 continue
!
!----
!
!  compute cloud microphysical parameters
!
      rgas= 287.
      t273 = 273.16
      t253 = 253.16
!
!--ncld=2 begin
      if( ncld .ge. 2 ) then
!
!  clwc/ciwc for stable cloud
!  clwc/ciwc are nomalized by cloudiness to take account of real water
!  content in clouds since qtr stands for grid average
!
      do k = 1, lm
      do i = 1, nxj
       if( clds(i,k) .gt. 0.01 )then
       qtrc = qtr(i,k)*1000.*( 100.*pp(i,k)/(rgas*tt(i,k)) )
! test
!      fcld = max(0.1, clds(i,k))
       fcld = 1.
       qtrc = qtrc / fcld
         if( tt(i,k) .ge. t273 )then
          clwc(i,k) = qtrc
          ciwc(i,k) = 0.
         elseif( tt(i,k) .lt. t273 .and. tt(i,k) .ge. t253 )then
          ratio = 1.0 - (t273-tt(i,k))/(t273-t253)
          clwc(i,k) = qtrc*ratio
          ciwc(i,k) = qtrc-clwc(i,k)
         else
          ciwc(i,k) = qtrc
          clwc(i,k) = 0.
         endif
       else
         clwc(i,k) = 0.
         ciwc(i,k) = 0.
       endif
      enddo
      enddo
!
!  for convective clouds
!
       aln10 = log(10.)
!
      do k = 1, lm
      do i = 1, nxj

      ttc = tt(i,k) - t273
      if( cvclds(i,k) .gt. 0. ) then
        cvf = cvclds(i,k) / (clds(i,k)+cvclds(i,k)-clds(i,k)*cvclds(i,k))
!        cvf = cvclds(i,k) / (clds(i,k)+cvclds(i,k))
        cvf1= 1.0 - cvf
!
        if( ttc .ge. 20. ) then
          clwc(i,k) = clwc(i,k)*cvf1 + exp( qlwcv(1)*aln10 )*cvf
!
        else if( ttc .le. -70. ) then
          ciwc(i,k) = ciwc(i,k)*cvf1 + exp( qiwcv(15)*aln10 )*cvf
!
        else if( ttc .ge. 0. ) then
          nwc = int((20.-ttc)/5.+0.0001) + 1
!         if(nwc.gt.9)then
!           nwc=9
          if(nwc.gt.8)then
            nwc=8
          else if(nwc.lt.1)then
            nwc=1
          endif
          xxx = (20.-(nwc-1)*5.) - ttc
          rlwcv = qlwcv(nwc) - (qlwcv(nwc)-qlwcv(nwc+1))*xxx/5.
          clwc(i,k) = clwc(i,k)*cvf1 + exp( rlwcv*aln10 )*cvf
!
        else if( ttc .ge. -20. ) then
          nwc = int((20.-ttc)/5.+0.0001) + 1
!         if(nwc.gt.9)then
!           nwc=9
          if(nwc.gt.8)then
            nwc=8
          else if(nwc.lt.1)then
            nwc=1
          endif
          nic = nwc - 4
          if(nic.lt.1)nic=1
          xxx = (20.-(nwc-1)*5.) - ttc
          rlwcv = qlwcv(nwc) - (qlwcv(nwc)-qlwcv(nwc+1))*xxx/5.
          riwcv = qiwcv(nic) - (qiwcv(nic)-qiwcv(nic+1))*xxx/5.
          rrr = -ttc/20.
          clwc(i,k) = clwc(i,k)*cvf1 + exp( rlwcv*aln10 )*(1.-rrr)*cvf
          ciwc(i,k) = ciwc(i,k)*cvf1 + exp( riwcv*aln10 )*rrr*cvf
!
        else if( ttc .lt. -20. ) then
          nic = int(-ttc/5.+0.0001) + 1
!         if(nic.gt.15)then
!           nic=15
          if(nic.gt.14)then
            nic=14
          else if(nic.lt.1)then
            nic=1
          endif
          xxx = -(nic-1)*5. - ttc
          riwcv = qiwcv(nic) - (qiwcv(nic)-qiwcv(nic+1))*xxx/5.
          ciwc(i,k) = ciwc(i,k)*cvf1 + exp( riwcv*aln10 )*cvf
!
        endif
      endif

      enddo
      enddo
!
      endif
!--ncld=2 end
!
!--ncld=1 begin
      if( ncld .eq. 1)then
!
      do k = 1, lm
      do i = 1, nxj

      ctt = cvclds(i,k)+clds(i,k)
      ttc = tt(i,k)-t273
      if( ctt .gt. 0. ) then
       cvf = cvclds(i,k) / ctt
       cvf1= 1.0 - cvf
!
       if( ttc .ge. 20. ) then
          clwc(i,k) = exp( qlwc(1) *aln10 ) * cvf1 +     &
                      exp( qlwcv(1)*aln10 ) * cvf
!
       else if( ttc .le. -70. ) then
          ciwc(i,k) = exp( qiwc(15) *aln10 ) * cvf1 +    &
                      exp( qiwcv(15)*aln10 ) * cvf
!
       else if( ttc .ge. 0. ) then
          nwc = int((20.-ttc)/5.+0.0001) + 1
!         if(nwc.gt.9)then
!           nwc=9
          if(nwc.gt.8)then
            nwc=8
          else if(nwc.lt.1)then
            nwc=1
          endif
          xxx = (20.-(nwc-1)*5.) - ttc
          rlwc  = qlwc(nwc)  - (qlwc(nwc) -qlwc(nwc+1))*xxx/5.
          rlwcv = qlwcv(nwc) - (qlwcv(nwc)-qlwcv(nwc+1))*xxx/5.
          clwc(i,k) = exp( rlwc *aln10 ) * cvf1 +        &
                      exp( rlwcv*aln10 ) * cvf
!
      else if( ttc .ge. -20. ) then
          nwc = int((20.-ttc)/5.+0.0001) + 1
!         if(nwc.gt.9)then
!           nwc=9
          if(nwc.gt.8)then
            nwc=8
          else if(nwc.lt.1)then
            nwc=1
          endif
          nic = nwc - 4
          if(nic.lt.1)nic=1
          xxx = (20.-(nwc-1)*5.) - ttc
          rlwc = qlwc(nwc) - (qlwc(nwc)-qlwc(nwc+1))*xxx/5.
          riwc = qiwc(nic) - (qiwc(nic)-qiwc(nic+1))*xxx/5.
          rlwcv = qlwcv(nwc) - (qlwcv(nwc)-qlwcv(nwc+1))*xxx/5.
          riwcv = qiwcv(nic) - (qiwcv(nic)-qiwcv(nic+1))*xxx/5.
          rrr = -ttc/20.
          clwc(i,k) = exp( rlwc *aln10 ) * cvf1 +    &
                      exp( rlwcv*aln10 ) * cvf
          ciwc(i,k) = exp( riwc *aln10 ) * cvf1 +    &
                      exp( riwcv*aln10 ) * cvf
          clwc(i,k) = clwc(i,k) * (1.-rrr)
          ciwc(i,k) = ciwc(i,k) * rrr
!
       else if( ttc .lt. -20. ) then
          nic = int(-ttc/5.+0.0001) + 1
!         if(nic.gt.15)then
!           nic=15
          if(nic.gt.14)then
            nic=14
          else if(nic.lt.1)then
            nic=1
          endif
          xxx = -(nic-1)*5. - ttc
          riwc  = qiwc(nic) - (qiwc(nic)-qiwc(nic+1))*xxx/5.
          riwcv = qiwcv(nic) - (qiwcv(nic)-qiwcv(nic+1))*xxx/5.
          ciwc(i,k) = exp( riwc *aln10 ) * cvf1 +    &
                      exp( riwcv*aln10 ) * cvf
!
       endif
      endif

      enddo
      enddo
!
      endif
!--ncld=1 end
!
!  compute the effective size for cloud droplets
!
!  1.the formulation of effective size for water clouds same as the one
!    used by McFarlane et al. (1992) in Canadian climate model
!  2.effective size for ice clouds based on the table2 of Fu (1996), which
!    is organized from data of HP (1984).
!
      do k = 1, lm
      do i = 1, nxj
!      ctt = clds(i,k)+cvclds(i,k)
      ctt = clds(i,k)+cvclds(i,k)-clds(i,k)*cvclds(i,k)
      ttc = tt(i,k)-273.16
      if( ctt .gt. 0. ) then
          cre_min=5.
          cre_max=7.5
          if(ttc.gt.-10.)then
            cre(i,k) = 11.*clwc(i,k) + cre_min
          else if(ttc.gt.-20. .and. ttc.le.-10.)then
            cre(i,k) = 11.*clwc(i,k) + cre_min-cre_min*((ttc+10.)/20.)
          else
            cre(i,k) = 11.*clwc(i,k) + cre_max
          endif
!          cre(i,k) = 11.*clwc(i,k) + 4.
!          cre(i,k) = 11.*clwc(i,k) + 5.
          nn = int( (-ttc-20.)/5. + 0.0001 ) + 1
          nn = min( 8, max(1,nn) )
          cde(i,k) = deg(nn)
      end if
      enddo
      enddo
!
      return
      end
