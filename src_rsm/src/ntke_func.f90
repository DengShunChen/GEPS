!--------------------------------------------------------
! external functions
!------------------------------------------------------
      real function foealfa(tt)
      use mo_constants, only:rtwat,   &
                             rtice
!     foealfa is calculated to distinguish the three cases:
!
!                       foealfa=1            water phase
!                       foealfa=0            ice phase
!                       0 < foealfa < 1      mixed phase
!
!               input : tt = temperature
!
        implicit none
        real tt
         foealfa = min(1.,((max(rtice,min(rtwat,tt))-rtice) &
     &  /(rtwat-rtice))**2)

        return
      end function foealfa

      real function foelhm(tt)
      use mo_constants, only:als,     &
                             alv
        implicit none
        real tt,foealfa
        foelhm = foealfa(tt)*alv + (1.-foealfa(tt))*als
      return
      end function foelhm

      real function foeewm(tt)
      use mo_constants, only:c4ies,   &
                             c3ies,   &
                             c4les,   &
                             tmelt,   &
                             c3les,   &
                             c2es

        implicit none
        real tt,foealfa
        foeewm  = c2es * &
     &     (foealfa(tt)*exp(c3les*(tt-tmelt)/(tt-c4les))+ &
     &     (1.-foealfa(tt))*exp(c3ies*(tt-tmelt)/(tt-c4ies)))
      return
      end function foeewm

      real function foedem(tt)
      use mo_constants, only:c4ies,   &
                             c4les,   &
                             r5alscp, &
                             r5alvcp

        implicit none
        real tt,foealfa
        foedem  = foealfa(tt)*r5alvcp*(1./(tt-c4les)**2)+ &
     &              (1.-foealfa(tt))*r5alscp*(1./(tt-c4ies)**2)
      return
      end function foedem

      real function foeldcpm(tt)
      use mo_constants, only:ralsdcp, &
                             ralvdcp

        implicit none
        real tt,foealfa
        foeldcpm = foealfa(tt)*ralvdcp+ &
     &        (1.-foealfa(tt))*ralsdcp
      return
      end function  foeldcpm

