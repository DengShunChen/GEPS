      subroutine cup92 ( nn,nx,lev,ktmin,dt,g,r,cp,hltm,etop,prevap  &
                       , pl,pk,zl,tl,ql,qls,ssl,hhl,hhls,pe,pk2,zl2  &
                       , tl2,ql2,qls2,ssl2,hhl2,hhls2,pmassl,hbase   &
                       , qbase,cubase,e,hc,qc,cu,kcbase,kctop,dhcup  &
                       , dqcup,flxmas,rcup,ncup,nlcl,nnegl           &
                       , nosat,nwork,ntcup,nflx )
!
!***********************************************************************
!
!    1.   function description
!           compute large-scale change due to cumulus convection
!           based upon arakawa-schubert cumulus parameterization scheme
!
!    2.   common block specification
!            none
!
!    3.   parameter specification
!     input :
!       nn        :   number of cumlus points in this latitude ring
!       nx        :   e-w (longitude) dimension of input arrays
!       lev       :   vertical levels
!       ktmin     :   highest sigma level allowed for the top of
!       dt        :   time step for calling cup92
!       g         :   gravity (9.806)
!       r         :   gas constant (287.0)
!       cp        :   heat specific constant of the air (=1004.)
!       hltm      :   latent heat constant of water vapor (2.52e+6)
!       etop      :   coeff for enhanced cloud top entrainment
!       prevap    :   fraction of falling precip evaporated
!       pl(nx,lev):   pressure on model level
!       pk(nx,lev):   p**kapa  on model level
!       pk2(nx,lev):  p**capa on even levels
!       zl(nx,lev) :   height on model level
!       tl(nx,lev) :   temperature on model level
!       ql(nx,lev) :   water vapor mixing ratio on model level
!       qls(nx,lev):   saturation mixing ratio on model level
!       ssl(nx,lev):   static energy on model level
!       hhl(nx,lev):   moist static energy on model level
!       hhls(nx,lev):  saturation moist static energy on model level
!       pe(nx,lev) :   pressure on even level
!       zl2(nx,lev):   height on even level
!       tl2(nx,lev):   temperature on even level
!       ql2(nx,lev):   water vapor mixing ratio on even level
!       qls2(nx,lev):  saturation mixing ratio on even level
!       ssl2(nx,lev):  static energy on even level
!       hhl2(nx,lev):  moist static energy on even level
!       hhls2(nx,lev)  saturation moist static energy on even level
!       pmassl(nx,lev): mass in between model even sigma levels  (kg)
!                      = (pl2(i,l)-pl2(i,l-1))*100.0/g
!       e(nx,lev)  :  normalized entrainment rate
!       hc(nx,lev) :  cloud air moist static energy
!       qc(nx,lev) :  cloud air specific humidity
!       cu(nx,lev) :  cloud air liquid water
!
!     output:
!       kcbase(nx):   cloud base index
!       kctop(nx) :   cloud top index
!       hbase(nx) : cloud base moist static energy
!       qbase(nx) : cloud base specific humidity
!       cubase(nx): cloud base liquid water
!       dhcup(nx,lev): moisture static energy change by cup  (1/call)
!       dqcup(nx,lev): moisture change by cumulus convection (1/call)
!       flxmas(nx,lev): mass flux for each cloud type
!       rcup(nx):   precipitation due to cumulus convection (mm/call)
!       ncup    : no of active convective points in this call
!       nlcl(lev)  : no of clouds at level k with tops above lcl
!       nnegl(lev) : no of clouds at level k with positive lamda
!       nosat(lev) : no of clouds at level k with liquid water at top
!       nwork(lev) : no of clouds at level k with work func above ref
!       ntcup(lev) : no of clouds at level k with non-zero cbase masflx
!       nflx(lev)  : no of clouds exceeding mass flux limit
!
!    local work arrays
!
!       xmb(nx)      : cloud base mass flux
!       pk2(nx,lev)   p**kapa  on even level
!
!    4.  calling modules
!            cupcwb
!
!    5.  usage
!          call cup92 ( nn,nx,lev,ktmin,dt,g,r,cp,hltm,pl,pk,zl,tl
!         1           , ql,qls,ssl,hhl,hhls,pe,zl2,tl2,ql2,qls2
!         2           , ssl2,hhl2,hhls2,pmassl,hbase,qbase,cubase
!         3           , kcbase,kctop,dqcup,rcup,ncup,nlcl,nnegl
!         4           , nosat,nwork,ntcup,nflx )
!
!    6.  modules called
!             cem1, dcbase, evpcup, ancrt
!
!    modify to f90 by C-H Lee and sort by River Chen in 2015
!
!***********************************************************************
      use paramt
!
!     variables passed through argument lists
!
      implicit  none

      integer   nn,nx,lev,ktmin

      real      pl(nx,lev),pk(nx,lev),zl(nx,lev),tl(nx,lev),ql(nx,lev) &
              , qls(nx,lev),ssl(nx,lev),hhl(nx,lev),hhls(nx,lev)       &
              , pe(nx,lev),pk2(nx,lev),zl2(nx,lev),tl2(nx,lev)         &
              , ql2(nx,lev),qls2(nx,lev),ssl2(nx,lev),hhl2(nx,lev)     &
              , hhls2(nx,lev),pmassl(nx,lev)                           &
              , rcup(nx),hbase(nx),qbase(nx),cubase(nx),e(nx,lev)      &
              , hc(nx,lev),qc(nx,lev),cu(nx,lev),flxmas(nx,lev)        &
              , dhcup(nx,lev),dqcup(nx,lev)

      integer   kcbase(nx),kctop(nx),nlcl(lev),nnegl(lev)              &
              , nosat(lev),nwork(lev),ntcup(lev),nflx(lev)

      integer   ncup,negak,nudry,nflxlim,i,k,ik,l
      integer   ltop,lt,ilsum,num1,num2,num3,num4,itype,itop,ibot

      real      prevap,etop,hltm,cp,r,g,dt,rkmin,relaxas,critl,p608
      real      gamfac,pcon,rg,ocp,xkapa,ohltm,const1,rmfmax,ccritl
      real      tem1,tem2,tmp1,tmp2,denom,otest,ftest

!
!     local work arrays
!
! --- cup locally stored large-scale variables

      real      xmb(im), rain(im,lm,lm)                                &
       , dzlc(im,lm),   dzl2c(im,lm),  po(im,lm)                       &
       , to(im,lm),     qo(im,lm),     qswo(im,lm),   so(im,lm)        &
       , ho(im,lm),     hswo(im,lm),   gamwo(im,lm)                    &
       , te(im,lm),     qswe(im,lm),   se(im,lm)                       &
       , he(im,lm),     hswe(im,lm),   gamwe(im,lm)                    &
! --- cup locally calculated large-scale quantities
       , gamqo(im,lm),  gamqe(im,lm)                                   &
       , buof1(im,lm),  buof2(im,lm), ancrit(im)                       &
       , grain(im,lm),  work1(im,lm), ratio(im)                        &
       , dho(im,2), dqo(im,2), dhswo(im,2), dqswo(im,2)                &
! --- cumulus contributions
       , chct(im,lm,lm), cqct(im,lm,lm)                                &
! --- internal variables
       , awp(im,lm), cwork(im,lm), aforce(im,lm), akernl(im,lm)
!
      logical  exist(nx,lev)
!
      real     hoz(im,lm),qoz(im,lm),hswoz(im,lm),qswoz(im,lm)
!
!***********************************************************************
! --- prescribed parameters
!
!lim  data  rkmin/-.005/, relaxas/0.10/
!
!  suggested by Ming-Deng
!
      rkmin = -0.02
      relaxas = dt/7200.
!
! --- constant setting on the first call of cup
!
      chct=0.
      dho=0.
!CWB2015
      gamwe=0.
      ho=0.
      qo=0.
      hswo=0.
      qswo=0.
!
        critl  = 0.
        p608  = 0.608
        gamfac = hltm*5417.9827/cp
        pcon   = 0.002
        rg     = 0.01*g
        ocp    = 1.0/cp
        xkapa= 1.0/3.5
        ohltm  = 1.0/hltm
        const1 = cp*p608*ohltm
!        rmfmax = (500./lev)/(dt/1800.)
!        rmfmax = (500./lev)/(dt/1200.)
!        rmfmax = (150./lev)/(dt/600.)
        rmfmax = (200./lev)/(dt/600.)
!
!cc   linear instability requires: rmfmax*0.001*dt/dp < 1.0
!     constant descriptions
!     cp :     specific heat of dry air at constant pressure  ( j/kg/k )
!     hltm :   latrent heat of condensation ( j/kg )
!     g :   gitational acceleration ( m/s )
!     pcon :   auto - conversion coefficient ( 1/m )
!     gamfacc: hltm/cp
!     critl :  critical ratio for auto - conversion of rainwater cloud
!              liquid water
!     p608 :  = 0.608
!     rg :  0.01*g
!     r :   gas constant ( j/kg/k )
!     rmfmax : maximum value for mass flux (kg/m2/s)
!     ocp                          1./1004.
!     ohltm                        1./2.52*10.**6
!
!***********************************************************************
!
! --- initialization
!
        ccritl= pcon*critl
!
        negak  = 0
        nudry  = 0
        nflxlim= 0
!
      do 5 i=1,nx*lev
        dhcup(i,1)  = 0.
        dqcup(i,1)  = 0.
        flxmas(i,1) = 0.
        exist(i,1)= .true.
        grain(i,1)= 0.0
    5 continue
!
      do 4 k = 1, ktmin-1
      do 4 i = 1, nn
      exist(i,k) = .false.
    4 continue
!
      do  6 i = 1,nn
      xmb(i)     = 0.
      rcup(i)    = 0.
   6  continue
!
      do ik=1,nx*lev*lev
       rain(ik,1,1)=0.
       chct(ik,1,1)=0.
       cqct(ik,1,1)=0.
      enddo
!
!    define the local large-scale environment
!
!      do 7 k = ktmin-1, lev
      do 7 k = 1, lev
      do 7 i = 1,nn
      ql(i,k)   = min (ql(i,k), 0.999*qls(i,k))
      hhl(i,k)  = ssl(i,k) + hltm*ql(i,k)
      po(i,k)   = pl(i,k)
      to(i,k)   = tl(i,k)
      so(i,k)   = ssl(i,k)
      ho(i,k)   = hhl(i,k)
      qo(i,k)   = ql(i,k)
      hswo(i,k) = hhls(i,k)
      qswo(i,k) = qls(i,k)
      qswe(i,k) = qls2(i,k)
      hswe(i,k) = hhls2(i,k)
      se(i,k)   = ssl2(i,k)
      he(i,k)   = hhl2(i,k)
      te(i,k)   = tl2(i,k)
   7  continue
!
!***********************************************************************
!
! cumulus parameterization subroutine
!
!***********************************************************************
!
!  compute a number of cloud dependent parameters for later use
!
      do 12 l = ktmin-1, lev
      do 12 i = 1,nn
        gamwo(i,l)= gamfac*qswo(i,l)/(to(i,l)*to(i,l))
        gamwe(i,l)  = gamfac*qswe(i,l)/(te(i,l)*te(i,l))
        tem1        = const1*te(i,l)
        tem2        = 1.0 + tem1*gamwe(i,l)
        gamqe(i,l)  = ohltm*gamwe(i,l)/(1.0+gamwe(i,l))
!
        dzlc(i,l)   = zl(i,l) - zl2(i,l)
        gamqo(i,l) = ohltm*gamwo(i,l)/(1. + gamwo(i,l))
        dzl2c(i,l)  = zl2(i,l-1) - zl2(i,l)
   12 continue
!
! interpolation for an incomplete cloud top layer:
!
! the values of "ratio" in the following statement must not be too
!   small considering the additional cloud water loading effect.
! if ratio is very close to 1.0, no interpolation is needed.
!
      do 2402 i = 1,nn
       ltop = kctop(i)
       tmp1= hbase(i)-hswe(i,ltop)
       tmp2= hswe(i,ltop-1)-hswe(i,ltop)
       if( tmp1 .gt. 0.9*tmp2 ) then
         ratio(i) = 1.
       else
         if( tmp1 .lt. 0.1*tmp2 ) then
           ratio(i) = 1.
           kctop(i) = ltop+1
         else
           ratio(i)=tmp1/tmp2
         endif
       endif
 2402 continue
!
      do 2401 i = 1,nn
       ltop = kctop(i)
       po(i,ltop)   = pe(i,ltop) - (pe(i,ltop) - po(i,ltop))*ratio(i)
!cc    pk(i,ltop)   = (0.001*po(i,ltop))**xkapa
       pk(i,ltop)   = exp(xkapa*log(0.001*po(i,ltop)))
       dzlc(i,ltop) = dzlc(i,ltop)*ratio(i)
       dzl2c(i,ltop)= dzl2c(i,ltop)*ratio(i)
       so(i,ltop)   = se(i,ltop) + (so(i,ltop) - se(i,ltop))*ratio(i)
       qo(i,ltop)   = (ho(i,ltop) - so(i,ltop))*ohltm
       to(i,ltop)   = te(i,ltop) + (so(i,ltop) - se(i,ltop))*ocp &
                    - g*dzlc(i,ltop)*ocp
       pmassl(i,ltop)= pmassl(i,ltop)*ratio(i)
 2401 continue
!
      do 2407 k = ktmin, lev
      call qsatq (nn,to(1,k),po(1,k),work1(1,k))
 2407 continue
!
      do 2409 i = 1,nn
        ltop = kctop(i)
        qswo(i,ltop) = work1(i,ltop)
        hswo(i,ltop) = so(i,ltop) + hltm*work1(i,ltop)
        gamwo(i,ltop)= gamfac*work1(i,ltop)/(to(i,ltop)*to(i,ltop))
        gamqo(i,ltop)= ohltm*gamwo(i,ltop)/(1. + gamwo(i,ltop))
        qo(i,ltop)   = min (qo(i,ltop), 0.999*qswo(i,ltop))
        ho(i,ltop)   = so(i,ltop) + hltm*qo(i,ltop)
 2409 continue
!
      hoz = ho
      qoz = qo
      hswoz = hswo
      qswoz = qswo
!
!
!  compute weight factors used in quadrature of work function
!  integral ( based on moorthi and suarez)
!
      do 335 k = ktmin+1, lev-1
      do 335 i = 1, nn
      if (k.gt.kcbase(i))  then
        buof1(i,k)= 0.0
        buof2(i,k)= 0.0
        exist(i,k)= .false.
      else
        denom     = 1.0/(pk(i,k)*(1.0+gamwo(i,k)))
        buof1(i,k)= (pk(i,k)-pk2(i,k-1))*denom
        buof2(i,k)= (pk2(i,k)-pk(i,k))*denom
      endif
  335 continue
!
      k = ktmin
      do 337 i = 1, nn
      buof2(i,k)= (pk2(i,k)-pk(i,k))/(pk(i,k)*(1.0+gamwo(i,k)))
  337 continue
!
! ********************************************************
! **  spectral cumulus ensemble model (cloud type: lt)  **
! ********************************************************
!
      do 500 lt = ktmin, lev-1
!
      num1    = ilsum( nn, exist(1,lt),1 )
      nlcl(lt)= num1
      if (num1.eq.0) go to 500
!
!     linear cloud model to calculate alamda  ( moothi and suarez )
!
      call cem1 ( nn,nx,lev,lt,cp,r,g,pcon,ccritl,etop,kcbase,.false. &
                , hswo,hbase, hswo,ho,qo,qswo,hswe,qswe,gamqo,gamqe   &
                , dzlc,dzl2c,e,exist(1,lt),hc,qc,cu,rain(1,1,lt)      &
                , zl,zl2 )
!
      num2     = ilsum( nn,exist(1,lt),1 )
      nnegl(lt)= num2
      if (num2.eq.0) go to 500
!
! --- we required that cloud air must be saturated at cloud top.
!
      do 501 i = 1,nn
      exist(i,lt)= exist(i,lt) .and. (cu(i,lt-1).gt.0.0)
  501 continue
!
      num3     = ilsum( nn,exist(1,lt),1 )
      nosat(lt)= num3
      if (num3.eq.0) go to 500
!
! --- calculate cloud work function (units: j/kg)
!       store the thermal byouancy per unit mass (m/s/s) at even levels
!       effective mass flux at odd levels (normalized by mb)
!       air density at even levels (g/cm/cm/cm)
!       precipitation per unit kg/m**2
!
      do 5020 i=1,nn
        cwork(i,lt)= e(i,lt)*buof2(i,lt)*(hc(i,lt)-hswo(i,lt))
 5020 continue
!
        do 120 l = lt+1, lev-1
        do 122 i = 1, nn
        cwork(i,lt)= cwork(i,lt)+e(i,l-1)*buof1(i,l)*(hc(i,l-1) &
                   - hswo(i,l))+e(i,l)*buof2(i,l)*(hc(i,l)-hswo(i,l))
  122 continue
  120 continue
!
!  reference climotological work function
!
       call ancrt (nn,po(1,lt),ancrit)
!
        do 125 i = 1, nn
        ancrit(i)= ancrit(i)*(pe(i,kcbase(i))-po(i,lt))
        ancrit(i)= max(ancrit(i),0.0)
        aforce(i,lt)= cwork(i,lt)-ancrit(i)
        exist(i,lt) = exist(i,lt) .and. (aforce(i,lt).gt.0.0)
  125 continue
!
      num4     = ilsum( nn,exist(1,lt),1 )
      nwork(lt)= num4
      if (num4.eq.0) go to 500
!
!  donor cell method for environmental budgets, reduces excessive
!  low level drying of normal flux form
!
      do 5030 i = 1,nn
      if (exist(i,lt) ) then
        chct(i,lt,lt) =  e(i,lt-1)*(hc(i,lt-1) - ho(i,lt))
        cqct(i,lt,lt) =  e(i,lt-1)*(qc(i,lt-1) - qo(i,lt))
      endif
 5030   continue
!
      do 140 l = lt+1, lev
      do 140 i=1,nn
      if (exist(i,lt)) then
          chct(i,l,lt) = e(i,l-1)*(ho(i,l-1)-ho(i,l))
          cqct(i,l,lt) = e(i,l-1)*(qo(i,l-1)-qo(i,l))
      endif
  140   continue
!
      do 220 l = lt, lev
      do 220 i = 1, nn
      if ( exist(i,lt) )  then
         chct(i,l,lt) = chct(i,l,lt)/pmassl(i,l)
         cqct(i,l,lt) = cqct(i,l,lt)/pmassl(i,l)
      endif
  220 continue
!
  500 continue
!
!
! *****************************
! **  begin of procedure 5)  **
! *****************************
!
! --- to determine the kernel elements, the following local
!       large-scale variables have to be modified:
!         qo, qswo, ho, hswo,  qswe, hswe, te, se
!       changes of the following parameters are neglected:
!         gamqo, gamqe, buof1, buof2
!         dzl, dzl2
!
      do 510 itype = ktmin, lev-1
!
      num4= ilsum( nn,exist(1,itype),1 )
      if (num4.eq.0) go to 510
        ftest  = 1.0
!
! --- adjust vetical profiles of thermodynamic variables
!       based on the test flux
!
        itop= 1
        ibot= 2
        do 545 i = 1, nn
        dho(i,itop)        = ftest*chct(i,itype,itype)
        dqo(i,itop)        = ftest*cqct(i,itype,itype)
        dhswo(i,itop)= (1. + gamwo(i,itype))*(dho(i,itop) &
                     - hltm*dqo(i,itop))
        dqswo(i,itop)      = dhswo(i,itop)*gamqo(i,itype)
        ho(i,itype) = hoz(i,itype)   + dho(i,itop)
        qo(i,itype) = qoz(i,itype)   + dqo(i,itop)
        hswo(i,itype) = hswoz(i,itype) + dhswo(i,itop)
        qswo(i,itype) = qswoz(i,itype) + dqswo(i,itop)
  545 continue
!
        do 540 l = itype+1, lev
        do 555 i = 1, nn
        dho(i,ibot)      = ftest*chct(i,l,itype)
        dqo(i,ibot)      = ftest*cqct(i,l,itype)
        dhswo(i,ibot)= (1. + gamwo(i,l))*(dho(i,ibot) &
                     - hltm*dqo(i,ibot))
        dqswo(i,ibot)    = dhswo(i,ibot)*gamqo(i,l)
        ho(i,l)     = hhl(i,l)  + dho(i,ibot)
        qo(i,l)     = ql(i,l)  + dqo(i,ibot)
        hswo(i,l)   = hhls(i,l)   + dhswo(i,ibot)
        qswo(i,l)   = qls(i,l)    + dqswo(i,ibot)
        qswe(i,l-1) = qls2(i,l-1) + .5*(dqswo(i,itop) + dqswo(i,ibot))
        hswe(i,l-1) = hhls2(i,l-1)+ .5*(dhswo(i,itop) + dhswo(i,ibot))
  555 continue
      itop = ibot
      ibot = 3-itop
  540 continue
!
! ********************************************************
! **  spectral cumulus ensemble model (cloud type: lt)  **
! ********************************************************
!
!     re-determine the properties of the sub cloud layers
!
!     call dcbase to determine cloud base level and cloud base condition
!
      call dcbase ( nn,nx,lev,hltm,kcbase,zl2,qo,qswe,ho,hswe,gamwe &
                  , e,qc,hc,cu,hbase,qbase,cubase )
!
      call cem1 ( nn,nx,lev,itype,cp,r,g,pcon,ccritl,etop,kcbase,.true.  &
                , hswo,hbase,hswo,ho,qo,qswo,hswe,qswe,gamqo,gamqe       &
                , dzlc,dzl2c,e,exist(1,itype),hc,qc,cu,work1             &
                , zl,zl2 )
!
! --- calculate the cloud work function only (units: j/kg)
!
      lt= itype
!
      do 502 i = 1, nn
      awp(i,lt)= e(i,lt)*buof2(i,lt)*(hc(i,lt)-hswo(i,lt))
 502  continue
!
        do 225 l = lt+1, lev-1
        do 225 i = 1, nn
        awp(i,lt) = awp(i,lt)+e(i,l-1)*buof1(i,l)*(hc(i,l-1) &
                  - hswo(i,l))+e(i,l)*buof2(i,l)*(hc(i,l)-hswo(i,l))
  225   continue
!
      nflxlim = 0
      otest= 1.0/ftest
      do 590 i = 1, nn
      if (exist(i,itype)) then
        akernl(i,itype) = (awp(i,itype) - cwork(i,itype))*otest
!
! --- checking the diagnal component of the mass flux kernel
!       (the criterion is arbitory)
!
      if (akernl(i,itype) .gt. rkmin)  then
        akernl(i,itype) = rkmin
        negak           = negak + 1
      endif
!
      flxmas(i,itype) = -relaxas*aforce(i,itype)/akernl(i,itype)
!
      if (flxmas(i,itype) .gt. rmfmax)  then
        flxmas(i,itype)= rmfmax
        nflxlim        = nflxlim + 1
      endif
      endif
  590 continue
!
      ntcup(itype)= ilsum( nn,exist(1,itype),1 )
      nflx(itype) =  nflxlim
!
  510 continue
!
!  how many convective grid points?
!
      do 520 k = ktmin+1, lev-1
      do 520 i = 1, nn
      exist(i,ktmin)= exist(i,ktmin).or.exist(i,k)
  520 continue
      ncup= ilsum( nn,exist(1,ktmin),1 )
      if (ncup.eq.0)  return
!
!  evaporation of cumulus precipitation
!
      if (prevap.gt.0.0)  then
!
      call evpcup (nn,nx,lev,ktmin,prevap,ql,qls,pmassl,rain,cqct)
!
      endif
!
!  compute rain as equal to integrated column drying
!
      do 525 itype = ktmin, lev-1
      do 535 i = 1, nn
      grain(i,itype)= -cqct(i,lev,itype)*pmassl(i,lev)
  535 continue
!
      do 525 k = itype, lev-1
      do 525 i = 1, nn
      grain(i,itype)= grain(i,itype)-cqct(i,k,itype)*pmassl(i,k)
  525 continue
!
! *****************************
! **  begin of procedure 8)  **
! *****************************
!
! --- total cumulus mass flux at the cloud base level
!     total cumulus precipitation
!
      do 810 itype = ktmin, lev-1
      do 810 i = 1,nn
        xmb(i) = xmb(i)  + flxmas(i,itype)
        rcup(i)= rcup(i) + flxmas(i,itype)* grain(i,itype)
  810 continue
!
! --- total cumulus contributions to temperature and moisture changes
!
      do 840 itype = ktmin, lev-1
      do 840 l = itype, lev
      do 840 i = 1,nn
          dqcup(i,l)  = dqcup(i,l) + flxmas(i,itype)*cqct(i,l,itype)
          dhcup(i,l)  = dhcup(i,l) + flxmas(i,itype)*chct(i,l,itype)
  840   continue
!
!  adjust changes for deepest clouds due to ratio being .lt. 1.0
!
      do 850 i = 1, nn
      dhcup(i,kctop(i))= dhcup(i,kctop(i))*ratio(i)
      dqcup(i,kctop(i))= dqcup(i,kctop(i))*ratio(i)
  850 continue
!
      return
      end
