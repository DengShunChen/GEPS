      subroutine sfc_diag(imj,im,ps,u1,v1,t1,q1,                   &
!    &                    tskin,qsurf,f10m,u10m,v10m,t2m,q2m,      &
                          tskin,qsurf,u10m,v10m,t2m,q2m, &
                          prslki,evap,fm,fh,fm10,fh2,fh10,rh2,rh10)
!
      use machine , only : kind_phys
!     use funcphys, only : fpvs
      use physcons, grav => con_g,  cp => con_cp, &
                    eps => con_eps, epsm1 => con_epsm1
      use const,    only : RTYPE
      implicit none
!
      integer              im,imj
      real(kind=RTYPE), dimension(im) :: u1,  v1,  t1,  q1
      real, dimension(im) :: ps,   tskin,  qsurf,                       &
                             f10m, u10m, v10m, t2m, t10m, q2m, q10m,    &
                             prslki,evap,fm,fh,fh10,fm10,fh2,rh2,rh10
!
!     locals
!
      real (kind=kind_phys), parameter :: qmin=1.0e-8
      integer              k,i
!
      real(kind=kind_phys)        fhi, qss, wrk,fpvs
!     real(kind=kind_phys) sig2k, fhi, qss
!
!     real, parameter :: g=grav
!
!     estimate sigma ** k at 2 m
!
!     sig2k = 1. - 4. * g * 2. / (cp * 280.)
!
!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals
!
!!
      do i = 1, imj
        f10m(i) = fm10(i) / fm(i)
!       f10m(i) = min(f10m(i),1.)
        u10m(i) = f10m(i) * u1(i)
        v10m(i) = f10m(i) * v1(i)
        fhi     = min( fh2(i) / fh(i) , 0.99 )
!       t2m(i)  = tskin(i)*(1. - fhi) + t1(i) * prslki(i) * fhi
!       sig2k   = 1. - (grav+grav) / (cp * t2m(i))
!       t2m(i)  = t2m(i) * sig2k
        wrk     = 1.0 - fhi

        t2m(i)  = tskin(i)*wrk + t1(i)*prslki(i)*fhi - (grav+grav)/cp

        if(evap(i) >= 0.) then !  for evaporation>0, use inferred qsurf to deduce q2m
          q2m(i) = qsurf(i)*wrk + max(qmin,q1(i))*fhi
        else                   !  for dew formation, use saturated q at tskin
          qss    = fpvs(tskin(i))
          qss    = eps * qss / (ps(i) + epsm1 * qss)
          q2m(i) = qss*wrk + max(qmin,q1(i))*fhi
        endif
        qss    = fpvs(t2m(i))
        qss    = eps * qss / (ps(i) + epsm1 * qss)
        q2m(i) = min(q2m(i),qss)
!
        rh2(i) = max( q2m(i)/qss , 0.)

        fhi     = min( fh10(i) / fh(i) , 0.99 )
        wrk     = 1.0 - fhi
        t10m(i)  = tskin(i)*wrk + t1(i)*prslki(i)*fhi - (grav+grav)/cp
        if(evap(i) >= 0.) then !  for evaporation>0, use inferred qsurf to deduce q2m
          q10m(i) = qsurf(i)*wrk + max(qmin,q1(i))*fhi
        else                   !  for dew formation, use saturated q at tskin
          qss    = fpvs(tskin(i))
          qss    = eps * qss / (ps(i) + epsm1 * qss)
          q10m(i) = qss*wrk + max(qmin,q1(i))*fhi
        endif
        qss    = fpvs(t10m(i))
        qss    = eps * qss / (ps(i) + epsm1 * qss)
        q10m(i) = min(q10m(i),qss)
!
        rh10(i) = max( q10m(i)/qss , 0.)

      enddo

      return
      end
