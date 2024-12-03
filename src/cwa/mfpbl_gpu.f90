!>  \file mfpbl.f
!!  This file contains the subroutine that calculates the updraft properties and mass flux for use in the Hybrid EDMF PBL scheme.

!>  \ingroup PBL
!!  \brief This subroutine is used for calculating the mass flux and updraft properties.
!!
!!  The mfpbl routines works as follows: if the PBL is convective, first, the ascending parcel entrainment rate is calculated as a function of height. Next, a surface parcel is initiated according to surface layer properties and the updraft buoyancy is ca
!!  \param[in] im number of used points
!!  \param[in] ix horizontal dimension
!!  \param[in] km vertical layer dimension
!!  \param[in] ntrac number of tracers
!!  \param[in] delt physics time step
!!  \param[in] cnvflg flag to denote a strongly unstable (convective) PBL
!!  \param[in] zl height of grid centers
!!  \param[in] zm height of grid interfaces
!!  \param[in] thvx virtual potential temperature at grid centers (\f$ K \f$)
!!  \param[in] q1 layer mean tracer concentration (units?)
!!  \param[in] t1 layer mean temperature (\f$ K \f$)
!!  \param[in] u1 u component of layer wind (\f$ m s^{-1} \f$)
!!  \param[in] v1 v component of layer wind (\f$ m s^{-1} \f$)
!!  \param[in,out] hpbl PBL top height (m)
!!  \param[in,out] kpbl PBL top index
!!  \param[in] sflx total surface heat flux (units?)
!!  \param[in] ustar surface friction velocity
!!  \param[in] wstar convective velocity scale
!!  \param[out] xmf updraft mass flux
!!  \param[in,out] tcko updraft temperature (\f$ K \f$)
!!  \param[in,out] qcko updraft tracer concentration (units?)
!!  \param[in,out] ucko updraft u component of horizontal momentum (\f$ m s^{-1} \f$)
!!  \param[in,out] vcko updraft v component of horizontal momentum (\f$ m s^{-1} \f$)
!!
!!  \section general General Algorithm
!!  -# Determine an updraft parcel's entrainment rate, buoyancy, and vertical velocity.
!!  -# Recalculate the PBL height (previously calculated in moninedmf) and the parcel's entrainment rate.
!!  -# Calculate the mass flux profile and updraft properties.
!!  \section detailed Detailed Algorithm
!!  @{
      subroutine mfpbl_gpu(myim,ix,km,ntrac,delt,cnvflg,                      &
         zl,zm,thvx,q1,t1,u1,v1,hpbl,kpbl,                              &
         sflx,ustar,wstar,xmf,tcko,qcko,ucko,vcko)
!
      use machine , only : kind_phys
      use physcons, grav => con_g, cp => con_cp
      use param, only : my_max
      use index, only: jlistnum

!
      implicit none
!
      integer              myim(my_max), ix, km, ntrac
!    &,                    me
      integer              kpbl(ix,my_max)
      logical              cnvflg(ix,my_max)
      real(kind=kind_phys) delt
      real(kind=kind_phys) q1(ix,km,ntrac,my_max), t1(ix,km,my_max),                  &
                           u1(ix,km,my_max),  v1(ix,km,my_max),                       &
                           thvx(ix,km,my_max),                                 &
                           zl(ix,km,my_max),  zm(ix,km+1,my_max),                     &
                           hpbl(ix,my_max),   sflx(ix,my_max),    ustar(ix,my_max),          &
                           wstar(ix,my_max),  xmf(ix,km,my_max),                      &
                           tcko(ix,km,my_max),qcko(ix,km,ntrac,my_max),               &
                           ucko(ix,km,my_max),vcko(ix,km,my_max)
!
!  local variables and arrays
!
      integer   i, j, k, n, kmpbl, jj
!
      real(kind=kind_phys) dt2,     dz,      ce0, &
                           h1,      factor,  gocp, &
                           g,       c1,      d1, &
                           b1,      f1,      bb1,     bb2,   &
                           alp,     a1,      qmin,    zfmin, &
                           xmmx,    rbint,   tau,            &
!    &                     rbint,   tau, &
                           tem,     tem1,    tem2, &
                           ptem,    ptem1,   ptem2, &
                           pgcon
!
      real(kind=kind_phys) sigw1(ix,my_max),   usws3(ix,my_max),  xlamax(ix,my_max), &
                           rbdn(ix,my_max),    rbup(ix,my_max),   delz(ix,my_max)
!
      real(kind=kind_phys) wu2(ix,km,my_max),     xlamue(ix,km,my_max), &
                           thvu(ix,km,my_max),    zi(ix,km,my_max), &
                           buo(ix,km,my_max)
!
      logical totflg, flg(ix,my_max)
!
!  physical parameters
      parameter(g=grav)
      parameter(gocp=g/cp)
!     parameter(ce0=0.37,qmin=1.e-8,alp=1.0,pgcon=0.55)
      parameter(ce0=0.38,qmin=1.e-8,alp=1.0,pgcon=0.55)
      parameter(a1=0.08,b1=0.5,f1=0.15,c1=0.3,d1=2.58,tau=500.)
      parameter(zfmin=1.e-8,h1=0.33333333)
!
!c-----------------------------------------------------------------------
!
!************************************************************************
!

      !jlistnum = 2
      kmpbl = km/2 + 1
      dt2 = delt
!> Since the mfpbl subroutine is called regardless of whether the PBL is convective, a check of the convective PBL flag is performed and the subroutine returns back to moninedmf (with the output variables set to the initialized values) if the PBL is not c
      totflg = .true.
      
      !!$acc parallel loop gang private(jj,i)
      !do jj = 1, jlistnum
      !  !$acc loop vector
      !  do i=1,myim(jj)
      !    !$acc atomic
      !    totflg = totflg .and. (.not. cnvflg(i,jj))
      !  enddo
      !enddo
      !if(totflg) return
  !!
      !$acc enter data create(buo,delz,flg,rbup,rbdn,sigw1,thvu,usws3,wu2, &
      !$acc&                  xlamax,xlamue,zi)
      !$acc parallel loop gang collapse(2) private(jj,k,i)
      do jj = 1, jlistnum
        do k = 1, km
          !$acc loop vector 
          do i=1,myim(jj)
            if (cnvflg(i,jj)) then
              zi(i,k,jj) = zm(i,k+1,jj)
            endif
          enddo
        enddo
      enddo
  !>  ## Determine an updraft parcel's entrainment rate, buoyancy, and vertical velocity.
  !!  Calculate the entrainment rate according to equation 16 in Siebesma et al. (2007) \cite siebesma_et_al_2007 for all levels (xlamue) and a default entrainment rate (xlamax) for use above the PBL top.
      !$acc parallel loop gang private(jj,i,k)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
          if(cnvflg(i,jj)) then
            k = kpbl(i,jj) / 2
            k = max(k, 1)
            delz(i,jj) = zl(i,k+1,jj) - zl(i,k,jj)
            xlamax(i,jj) = ce0 / delz(i,jj)
          endif
        enddo
      enddo
      !$acc parallel loop gang collapse(2) private(jj,k,i,ptem,tem,ptem1)
      do jj = 1, jlistnum
        do k = 1, kmpbl
          !$acc loop vector
          do i=1,myim(jj)
            if(cnvflg(i,jj)) then
              if(k < kpbl(i,jj)) then
                ptem = 1./(zi(i,k,jj)+delz(i,jj))
                tem = max((hpbl(i,jj)-zi(i,k,jj)+delz(i,jj)) ,delz(i,jj))
                ptem1 = 1./tem
                xlamue(i,k,jj) = ce0 * (ptem+ptem1)
              else
                xlamue(i,k,jj) = xlamax(i,jj)
              endif
            endif
          enddo
        enddo
      enddo
  !
  !  compute thermal excess
  !
  !>  Using equations 17 and 7 from Siebesma et al (2007) \cite siebesma_et_al_2007 along with \f$u_*\f$, \f$w_*\f$, and the previously diagnosed PBL height, the initial \f$\theta_v\f$ of the updraft (and its surface buoyancy) is calculated.
      !$acc parallel loop gang private(jj,i,tem,tem1,tem2,ptem,ptem1)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
          if(cnvflg(i,jj)) then
            tem = zl(i,1,jj)/hpbl(i,jj)
            usws3(i,jj) = (ustar(i,jj)/wstar(i,jj))**3.
            tem1 = usws3(i,jj) + 0.6*tem
            tem2 = max((1.-tem), zfmin)
            ptem = (tem1**h1) * sqrt(tem2)
            sigw1(i,jj) = 1.3 * ptem * wstar(i,jj)
            ptem1 = alp * sflx(i,jj) / sigw1(i,jj)
            thvu(i,1,jj) = thvx(i,1,jj) + ptem1
            buo(i,1,jj) = g * (thvu(i,1,jj)/thvx(i,1,jj)-1.)
            !if (jj .eq. 5 .and. i .eq. 3) write(*,*) 'b20',thvu(i,1,jj),thvx(i,1,jj)
            !if (jj .eq. 5 .and. i .eq. 3) write(*,*) 'b21',sflx(i,jj) / sigw1(i,jj)
          endif
        enddo
      enddo
  !
  !  compute potential temperature and buoyancy for updraft air parcel
  !
  !>  From the second level to the middle of the vertical domain, the updraft virtual potential temperature is calculated using the entraining updraft equation as in equation 10 of Siebesma et al (2007) \cite siebesma_et_al_2007, discretized as
  !!  \f[
  !!  \frac{\theta_{v,u}^k - \theta_{v,u}^{k-1}}{\Delta z}=-\epsilon^{k-1}\left[\frac{1}{2}\left(\theta_{v,u}^k + \theta_{v,u}^{k-1}\right)-\frac{1}{2}\left(\overline{\theta_{v}}^k + \overline{\theta_v}^{k-1}\right)\right]
  !!  \f]
  !!  where the superscript \f$k\f$ denotes model level, and subscript \f$u\f$ denotes an updraft property, and the overbar denotes the grid-scale mean value.
      !$acc parallel loop gang private(jj,k,i,dz,tem,ptem,ptem1,tem1)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
          !$acc loop seq
          do k = 2, kmpbl
            if(cnvflg(i,jj)) then
              dz = zl(i,k,jj) - zl(i,k-1,jj)
              tem = xlamue(i,k-1,jj) * dz
              ptem = 2. + tem
              ptem1 = (2. - tem) / ptem
              tem1 = tem  * (thvx(i,k,jj)+thvx(i,k-1,jj)) / ptem
              thvu(i,k,jj) = ptem1 * thvu(i,k-1,jj) + tem1
              buo(i,k,jj) = g * (thvu(i,k,jj)/thvx(i,k,jj)-1.)
            endif
          enddo
        enddo
      enddo
  !
  !  compute updraft velocity square(wu2)
  !
  !>  Rather than use the vertical velocity equation given as equation 15 in Siebesma et al (2007) \cite siebesma_et_al_2007 (which parameterizes the pressure term in terms of the updraft vertical velocity itself), this scheme uses the more widely used form
  !!  \f[
  !!  \frac{w_{u,k}^2 - w_{u,k-1}^2}{\Delta z} = -2b_1\frac{1}{2}\left(\epsilon_k + \epsilon_{k-1}\right)\frac{1}{2}\left(w_{u,k}^2 + w_{u,k-1}^2\right) + 2b_2B
  !!  \f]
  !! The constants used in the scheme are labeled \f$bb1 = 2b_1\f$ and \f$bb2 = 2b_2\f$ and are tuned to be equal to 1.8 and 3.5, respectively, close to the values proposed by Soares et al. (2004) \cite soares_et_al_2004 .
  !     tem = 1.-2.*f1
  !     bb1 = 2. * b1 / tem
  !     bb2 = 2. / tem
  !  from soares et al. (2004,qjrms)
  !     bb1 = 2.
  !     bb2 = 4.
  !
  !  from bretherton et al. (2004, mwr)
  !     bb1 = 4.
  !     bb2 = 2.
  !
  !  from our tuning
         bb1 = 1.8    !(org)
  !      bb1 = 2.8    !(lin)

        bb2 = 3.5
  !
      !$acc parallel loop gang private(jj,i,dz,tem,tem1,ptem1)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(cnvflg(i,jj)) then
  !
  !         tem = zi(i,1)/hpbl(i)
  !         tem1 = usws3(i) + 0.6*tem
  !         tem2 = max((1.-tem), zfmin)
  !         ptem = (tem1**h1) * sqrt(tem2)
  !         ptem1 = 1.3 * ptem * wstar(i)
  !         wu2(i,1) = d1*d1*ptem1*ptem1
  !
            dz   = zi(i,1,jj)
            tem  = 0.5*bb1*xlamue(i,1,jj)*dz
            tem1 = bb2 * buo(i,1,jj) * dz
            ptem1 = 1. + tem
            wu2(i,1,jj) = tem1 / ptem1
            !if (jj .eq. 5 .and. i .eq. 3) write(*,*) 'b1',zi(i,1,jj),xlamue(i,1,jj),buo(i,1,jj)
  !
          endif
        enddo
      enddo
      !$acc parallel loop gang private(jj,k,i,dz,tem,tem1,ptem,ptem1) 
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          !$acc loop seq
          do k = 2, kmpbl
            if(cnvflg(i,jj)) then
              dz    = zi(i,k,jj) - zi(i,k-1,jj)
              tem  = 0.25*bb1*(xlamue(i,k,jj)+xlamue(i,k-1,jj))*dz
              tem1 = bb2 * buo(i,k,jj) * dz
              ptem = (1. - tem) * wu2(i,k-1,jj)
              ptem1 = 1. + tem
              wu2(i,k,jj) = (ptem + tem1) / ptem1
            endif
          enddo
        enddo
      enddo
  !
  !  update pbl height as the height where updraft velocity vanishes
  !
  !>  ## Recalculate the PBL height and the parcel's entrainment rate.
  !!  Find the level where the updraft vertical velocity is less than zero and linearly interpolate to find the height where it would be exactly zero. Set the PBL height to this determined height.
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           flg(i,jj)  = .true.
           if(cnvflg(i,jj)) then
             flg(i,jj)  = .false.
             rbup(i,jj) = wu2(i,1,jj)
           endif
        enddo
      enddo
      !$acc parallel loop gang private(jj,k,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          !$acc loop seq
          do k = 2, kmpbl
            if(.not.flg(i,jj)) then
              rbdn(i,jj) = rbup(i,jj)
              rbup(i,jj) = wu2(i,k,jj)
              kpbl(i,jj) = k
              flg(i,jj)  = rbup(i,jj).le.0.
            endif
          enddo
        enddo
      enddo
      !$acc parallel loop gang private(jj,i,k,rbint)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
          if(cnvflg(i,jj)) then
             k = kpbl(i,jj)
             if(rbdn(i,jj) <= 0.) then
                rbint = 0.
             elseif(rbup(i,jj) >= 0.) then
                rbint = 1.
             else
                rbint = rbdn(i,jj)/(rbdn(i,jj)-rbup(i,jj))
             endif
             hpbl(i,jj) = zi(i,k-1,jj) + rbint*(zi(i,k,jj)-zi(i,k-1,jj))
          endif
        enddo
      enddo
  !
  !>  Recalculate the entrainment rate as before except use the updated value of the PBL height.
      !$acc parallel loop gang private(jj,i,k)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
          if(cnvflg(i,jj)) then
            k = kpbl(i,jj) / 2
            k = max(k, 1)
            delz(i,jj) = zl(i,k+1,jj) - zl(i,k,jj)
            xlamax(i,jj) = ce0 / delz(i,jj)
          endif
        enddo
      enddo
  !
  !  update entrainment rate
  !
  !     do k = 1, kmpbl
  !       do i=1,im
  !         if(cnvflg(i)) then
  !           if(k < kpbl(i)) then
  !             tem = tau * sqrt(wu2(i,k))
  !             tem1 = 1. / tem
  !             ptem = ce0 / zi(i,k)
  !             xlamue(i,k) = max(tem1, ptem)
  !           else
  !             xlamue(i,k) = xlamax(i)
  !           endif
  !         endif
  !       enddo
  !     enddo
  !
      !$acc parallel loop gang collapse(2) private(jj,k,i,ptem,tem,ptem1)
      do jj = 1, jlistnum
        do k = 1, kmpbl
          !$acc loop vector
          do i=1,myim(jj)
            if(cnvflg(i,jj)) then
              if(k < kpbl(i,jj)) then
                ptem = 1./(zi(i,k,jj)+delz(i,jj))
                tem = max((hpbl(i,jj)-zi(i,k,jj)+delz(i,jj)) ,delz(i,jj))
                ptem1 = 1./tem
                xlamue(i,k,jj) = ce0 * (ptem+ptem1)
              else
                xlamue(i,k,jj) = xlamax(i,jj)
              endif
            endif
          enddo
        enddo
      enddo
  !
  !  updraft mass flux as a function of sigmaw
  !   (0.3*sigmaw[square root of vertical turbulence variance])
  !
  !>  ## Calculate the mass flux profile and updraft properties.
  !     do k = 1, kmpbl
  !       do i=1,im
  !         if(cnvflg(i) .and. k < kpbl(i)) then
  !           tem = zi(i,k)/hpbl(i)
  !           tem1 = usws3(i) + 0.6*tem
  !           tem2 = max((1.-tem), zfmin)
  !           ptem = (tem1**h1) * sqrt(tem2)
  !           ptem1 = 1.3 * ptem * wstar(i)
  !           xmf(i,k) = c1 * ptem1
  !         endif
  !       enddo
  !     enddo
  !
  !  updraft mass flux as a function of updraft velocity profile
  !
  !>  Calculate the mass flux:
  !!  \f[
  !!  M = a_uw_u
  !!  \f]
  !!  where \f$a_u\f$ is the tunable parameter that represents the fractional area of updrafts (currently set to 0.08). Limit the computed mass flux to be less than \f$\frac{\Delta z}{\Delta t}\f$. This is different than what is done in Siebesma et al. (200
      !$acc parallel loop gang collapse(2) private(jj,k,i,dz,xmmx)
      do jj = 1, jlistnum
        !write(*,*) 'a4',jj, myim(jj)
        do k = 1, kmpbl
          !$acc loop vector
          do i = 1, myim(jj)
            !if (jj .eq. 2 .and. i .eq. 2 .and. k .eq. 1) write(*,*) 222
            if (cnvflg(i,jj) .and. k < kpbl(i,jj)) then
               xmf(i,k,jj) = a1 * sqrt(wu2(i,k,jj))
               dz   = zl(i,k+1,jj) - zl(i,k,jj)
               xmmx = dz / dt2
               xmf(i,k,jj) = min(xmf(i,k,jj),xmmx)
               !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 1) write(*,*) 'b0',xmf(i,k,jj),wu2(i,k,jj)
            endif
          enddo
        enddo
      enddo
      
  !
  !c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  compute updraft property
  !
  !>  The updraft properties are calculated according to the entraining updraft equation
  !!  \f[
  !!  \frac{\partial \phi}{\partial z}=-\epsilon\left(\phi_u - \overline{\phi}\right)
  !!  \f]
  !!  where \f$\phi\f$ is \f$T\f$ or \f$q\f$. The equation is discretized according to
  !!  \f[
  !!  \frac{\phi_{u,k} - \phi_{u,k-1}}{\Delta z}=-\epsilon_{k-1}\left[\frac{1}{2}\left(\phi_{u,k} + \phi_{u,k-1}\right)-\frac{1}{2}\left(\overline{\phi}_k + \overline{\phi}_{k-1}\right)\right]
  !!  \f]
  !!  The exception is for the horizontal momentum components, which have been modified to account for the updraft-induced pressure gradient force, and use the following equation, following Han and Pan (2006) \cite han_and_pan_2006
  !!  \f[
  !!  \frac{\partial v}{\partial z} = -\epsilon\left(v_u - \overline{v}\right)+d_1\frac{\partial \overline{v}}{\partial z}
  !!  \f]
  !!  where \f$d_1=0.55\f$ is a tunable constant.
      
      !$acc parallel loop gang private(jj,k,i,dz,tem,factor,ptem,ptem1)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          !$acc loop seq
          do k = 2, kmpbl
            if (cnvflg(i,jj) .and. k <= kpbl(i,jj)) then
               dz   = zl(i,k,jj) - zl(i,k-1,jj)
               tem  = 0.5 * xlamue(i,k-1,jj) * dz
               factor = 1. + tem
               ptem = tem + pgcon
               ptem1= tem - pgcon
  !
               tcko(i,k,jj) = ((1.-tem)*tcko(i,k-1,jj)+tem* &
                            (t1(i,k,jj)+t1(i,k-1,jj))-gocp*dz)/factor
               ucko(i,k,jj) = ((1.-tem)*ucko(i,k-1,jj)+ptem*u1(i,k,jj) &
                            +ptem1*u1(i,k-1,jj))/factor
               vcko(i,k,jj) = ((1.-tem)*vcko(i,k-1,jj)+ptem*v1(i,k,jj) &
                            +ptem1*v1(i,k-1,jj))/factor
            endif
          enddo
        enddo
      enddo
      !$acc parallel loop gang collapse(2) private(jj,n,k,i,dz,tem,factor)
      do jj = 1, jlistnum
        do n = 1, ntrac
          !$acc loop vector
          do i = 1, myim(jj)
            !$acc loop seq
            do k = 2, kmpbl
              if (cnvflg(i,jj) .and. k <= kpbl(i,jj)) then
                 dz   = zl(i,k,jj) - zl(i,k-1,jj)
                 tem  = 0.5 * xlamue(i,k-1,jj) * dz
                 factor = 1. + tem

                 qcko(i,k,n,jj) = ((1.-tem)*qcko(i,k-1,n,jj)+tem* &
                              (q1(i,k,n,jj)+q1(i,k-1,n,jj)))/factor
              endif
            enddo
          enddo
        enddo
      enddo
!
      !$acc exit data delete(buo,delz,flg,rbup,rbdn,sigw1,thvu,usws3,wu2, &
      !$acc&                 xlamax,xlamue,zi)

      return
      end


      
