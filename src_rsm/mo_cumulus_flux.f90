MODULE mo_cumulus_flux

  use shr_kind_mod, only: r8 => shr_kind_r8
  use param, only: jtrun

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_cumulus_flux* - parameters for cumulus massflux scheme
  !
  ! ----------------------------------------------------------------

  REAL(r8) :: entrpen      !    entrainment rate for penetrative convection
  REAL(r8) :: entrscv      !    entrainment rate for shallow convection
  REAL(r8) :: entrmid      !    entrainment rate for midlevel convection
  REAL(r8) :: entrdd       !    entrainment rate for cumulus downdrafts
  REAL(r8) :: cmfctop      !    relat. cloud massflux at level above nonbuoyanc
  REAL(r8) :: cmfcmax      !    maximum massflux value allowed for
  REAL(r8) :: cmfcmin      !    minimum massflux value (for safety)
  REAL(r8) :: cmfdeps      !    fractional massflux for downdrafts at lfs
  REAL(r8) :: rhcdd        !    relative saturation in downdrafts
  REAL(r8) :: cprcon       !    coefficients for determining conversion
  REAL(r8) :: cprcon_n       !    coefficients for determining conversion
                       !    from cloud water to rain
!> xb110
  REAL(r8) :: zdnoprc  ! deep cloud is thicker than this height (Unit:Pa)
  REAL(r8) :: momtrans ! momentum transport method
!< xb110
  LOGICAL :: lmfpen    !    true if penetrative convection is switched on
  LOGICAL :: lmfscv    !    true if shallow     convection is switched on
  LOGICAL :: lmfmid    !    true if midlevel    convection is switched on
  LOGICAL :: lmfdd     !    true if cumulus downdraft      is switched on
  LOGICAL :: lmfdudv   !    true if cumulus friction       is switched on

CONTAINS

SUBROUTINE cuparam

  ! Description:
  !
  ! Defines disposable parameters for massflux scheme
  !
  ! Method:
  !
  ! This routine is called from *iniphy*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, February 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, Jan 1999, subroutine cuparam -> module mo_cumulus_flux
  ! 
  ! for more details see file AUTHORS
  ! 

! USE mo_control, ONLY: lamip2

  IMPLICIT NONE


  !  Executable Statements 

  lmfpen = .TRUE.
  lmfscv = .TRUE.
  lmfmid = .TRUE.
  lmfdd = .TRUE.
  lmfdudv = .TRUE.

!-- 1. Specify parameters for massflux-scheme

  entrpen = 1.0E-4 ! Average entrainment rate for penetrative convection

!  entrscv = 1.2E-3 ! Average entrainment rate for shallow convection (tie3)
  entrscv = 3.0E-4 ! Average entrainment rate for shallow convection (Tie7 org)
!  entrscv = 6.0E-4 ! Average entrainment rate for shallow convection (Tie8)

  entrmid = 1.0E-4 ! Average entrainment rate for midlevel convection

  entrdd  = 2.0E-4 ! Average entrainment rate for downdrafts

  cmfctop = 0.3 ! Relative cloud massflux at level above nonbuoyancy level

  cmfcmax = 1.0    ! Maximum massflux value allowed for updrafts etc

  cmfcmin = 1.E-10 ! Minimum massflux value (for safety)

  cmfdeps = 0.3    ! Fractional massflux for downdrafts at lfs
  
  zdnoprc = 2.0e4  ! deep cloud is thicker than this height (Unit:Pa)     !xb110

!  cprcon  = 6.E-4  ! Coefficients for determining conversion from cloud water   (org)
!  cprcon  = 1.122E-4  ! Coefficients for determining conversion from cloud water(tie4)
!  cprcon  = 1.0E-3  ! Coefficients for determining conversion from cloud water (tie5)
   cprcon  = 1.0E-4  ! Coefficients for determining conversion from cloud water(tie7)   !xb110, for Tiedtke (1989)
   cprcon_n  = 1.4e-3  !Coefficients for determining conversion from cloud water (MPAS) !xb110, for new Tiedtke
   momtrans = 2  ! momentum transport method !xb110

  if(jtrun.eq.320)then
    cmfctop = 0.35
  endif

  ! Next value is relative saturation in downdrafrs
  ! but is no longer used ( formulation implies saturation)

  rhcdd = 1.

  RETURN
END SUBROUTINE cuparam

END MODULE mo_cumulus_flux
