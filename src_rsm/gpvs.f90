      subroutine gpvs
!
! Abstract: Computes saturation vapor pressure table as a function of
!   temperature for the table lookup function fpvs.
!   Exact saturation vapor pressures are calculated in subprogram fpvsx.
!   The current implementation computes a table with a length
!   of 7501 for temperatures ranging from 180. to 330. Kelvin.
!
      implicit none
      integer jx
      real xmin,xmax,xinc,x,t
      integer nxpvs
      real c1xpvs,c2xpvs,tbpvs,fpvsx
      parameter(nxpvs=7501)
      common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xmin=180.0
      xmax=330.0
      xinc=(xmax-xmin)/(nxpvs-1)
!   c1xpvs=1.-xmin/xinc
      c2xpvs=1./xinc
      c1xpvs=1.-xmin*c2xpvs
      do jx=1,nxpvs
        x=xmin+(jx-1)*xinc
        t=x
        tbpvs(jx)=fpvsx(t)
      enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end subroutine

      function fpvsx(t)
!
! Abstract: Exactly compute saturation vapor pressure from temperature.
!   The saturation vapor pressure over either liquid and ice is computed
!   over liquid for temperatures above the triple point,
!   over ice for temperatures 20 degress below the triple point,
!   and a linear combination of the two for temperatures in between.
!   The water model assumes a perfect gas, constant specific heats
!   for gas, liquid and ice, and neglects the volume of the condensate.
!   The model does account for the variation of the latent heat
!   of condensation and sublimation with temperature.
!   The Clausius-Clapeyron equation is integrated from the triple point
!   to get the formula
!         pvsl=psat*(tr**xa)*exp(xb*(1.-tr))
!   where tr is ttp/t and other values are physical constants.
!   The reference for this computation is Emanuel(1994), pages 116-117.
!   This function should be expanded inline in the calling routine.
!
!   Input argument list:
!       t              Real temperature in Kelvin
!
!   Output argument list:
!       fpvsx        Real saturation vapor pressure in Pascals
!
      implicit none
      real fpvsx
      real t
      real psat,hvap,cvap,cliq,csol,grav,hfus,ttp,rd,rv,cp
      parameter(  psat = 6.1078e+2 )
      parameter(  hvap = 2.5000e+6 )
      parameter(  cvap = 1.8460e+3 )
      parameter(  cliq = 4.1855e+3 )
      parameter(  csol = 2.1060e+3 )
      parameter(  grav = 9.80665e+0 )
      parameter(  hfus = 3.3358e+5 )
      parameter(   ttp = 2.7316e+2 )
      parameter(    rd = 2.8705e+2 )
      parameter(    rv = 4.6150e+2 )
      parameter(    cp = 1.0046e+3 )
      real tliq,tice,dldtl,heatl,xponal,xponbl,dldti,heati,xponai,xponbi
      parameter( tliq=ttp )
      parameter( tice=ttp-20.0 )
      parameter( dldtl=cvap-cliq )
      parameter( heatl=hvap )
      parameter( xponal=-dldtl/rv )
      parameter( xponbl=xponal+heatl/(rv*ttp) )
      parameter( dldti=cvap-csol )
      parameter( heati=hvap+hfus )
      parameter( xponai=-dldti/rv )
      parameter( xponbi=xponai+heati/(rv*ttp) )
      real tr,w,pvl,pvi
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tr=ttp/t
      if(t.ge.tliq) then
        fpvsx=psat*(tr**xponal)*exp(xponbl*(1.-tr))
      elseif(t.lt.tice) then
        fpvsx=psat*(tr**xponai)*exp(xponbi*(1.-tr))
      else
        w=(t-tice)/(tliq-tice)
        pvl=psat*(tr**xponal)*exp(xponbl*(1.-tr))
        pvi=psat*(tr**xponai)*exp(xponbi*(1.-tr))
        fpvsx=w*pvl+(1.-w)*pvi
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end function

