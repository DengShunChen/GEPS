#define LCWBGFS T
#if defined(LCWBGFS)
#define mpp_pe() myrank
#define mpp_root_pe() 0
#endif

#define DEBUG
#undef CTFREEZE
!!!#undef Prw001 T
!!!#undef Prw01
!!!#define Prw1
#define TMELTS0
#define SALTI0
!!! #define GDCHK2 (lwarning_msg.GE.1).AND.(((mpp_pe().EQ.2).AND.(jrow.EQ.1).AND.(jl.EQ.1)).OR.((mpp_pe().EQ.35).AND.(jrow.EQ.1).AND.(jl.EQ.1)))
!!! #define GDCHK2 (lwarning_msg.GE.1).AND.(((mpp_pe().EQ.2).AND.(jrow.EQ.2).AND.(jl.EQ.6)).OR.((mpp_pe().EQ.35).AND.(jrow.EQ.10).AND.(jl.EQ.3)))
!!!#define GDCHK2 (lwarning_msg.GE.1).AND.(mpp_pe().EQ.97).AND.(jrow.EQ.1).AND.(jl.EQ.15)
!!!#define GDCHK2 (lwarning_msg.GE.1).AND.((mpp_pe().EQ.2).AND.lsitmask(jl)).OR.((mpp_pe().EQ.35).AND.lsitmask(jl))
!!!#define GDCHK2 (lwarning_msg.GE.1).AND.((mpp_pe().EQ.2).OR.(mpp_pe().EQ.35))
! dateline (0N, 180E), if npe=96, pe = npe/2=48. And 48-1 =47. Note pe starts from 0
!!!#define GDCHK2 (lwarning_msg.GE.1).AND.(mpp_pe().EQ.47).AND.(jrow.EQ.43).AND.(jl.EQ.86)
!!!#define GDCHK2 (lwarning_msg.GE.1).AND.(mpp_pe().EQ.47).AND.(jrow.EQ.1)
!!!#define GDCHK2 (lwarning_msg.GE.1).AND.(mpp_pe().EQ.361)
!!!#define GDCHK2 (lwarning_msg.GE.1).AND.(mpp_pe().EQ.47)
!!!#define GDCHK2 (lwarning_msg.GE.1).AND.(mpp_pe().EQ.449)
!!!#define GDCHK2 (lwarning_msg.GE.1).AND.(mpp_pe().EQ.449).AND.(jrow.EQ.1)
!!!#define GDCHK2 (mpp_pe().EQ.216).AND.(jrow.EQ.1).AND.(jl.EQ.1)
!!!#define GDCHK0 (mpp_pe().EQ.216)
!!!#define GDCHK2 (mpp_pe().EQ.216).AND.(jrow.EQ.55)
!!!#define GDCHK2 (mpp_pe().EQ.216)
!!! #define GDCHK2 .FALSE.
!!!#define GDCHK0 (mpp_pe().EQ.399)
!!!#define GDCHK2 (mpp_pe().EQ.399).AND.(jrow.EQ.55)
!!!#define GDCHK0 (mpp_pe().EQ.817)
!!!#define GDCHK2 (mpp_pe().EQ.817).AND.(jrow.EQ.55)

!!!#define GDCHK0 (mpp_pe().EQ.495)
!!!#define GDCHK2 (mpp_pe().EQ.495)

!!!#define GDCHK0 (mpp_pe().EQ.530)
!!!#define GDCHK2 (mpp_pe().EQ.530)

!#define GDCHK0 (mpp_pe().EQ.740)
!#define GDCHK2 (mpp_pe().EQ.740)

!#define GDCHK0 (mpp_pe().EQ.0)
!#define GDCHK2 (mpp_pe().EQ.0)

!!!#define GDCHK0 (mpp_pe().EQ.734)
!!!#define GDCHK2 (mpp_pe().EQ.734)
!!!#define GDCHK2 ( (plon(jl).LE.-175.).AND.(plon(jl).GE.-185.).AND.(plat(jl).LE.-5.).AND.(plat(jl).GE.-13.) )
!!!#define GDCHK0 (mpp_pe().EQ.334)
!!!#define GDCHK1 ( (mpp_pe().EQ.334).AND.(jl.EQ.4) )
!!!#define GDCHK2 ( (mpp_pe().EQ.334).AND.(jl.EQ.4).AND.(jrow.EQ.80) )
!!!#define GDCHK2 ( (mpp_pe().EQ.334).AND.(jl.EQ.4).AND.(jrow.EQ.78) )

!!!Default file number is: 1
!!!X is fixed     Lon = -99.5206  X = 181
!!!Y is fixed     Lat = -30.5  Y = 52
!!!Z is fixed     Lev = 0  Z = 1
!!!T is varying   Time = 01:30Z01JAN1980 to 22:30Z01JAN1980  T = 1 to 8
!!!E is fixed     Ens = 1  E = 1
!!!#define GDCHK2 ( (plon(jl).GT.-100.).AND.(plon(jl).LE.-99.).AND.(plat(jl).LE.-30.).AND.(plat(jl).GE.-31.) )

!!!#define GDCHK0 (mpp_pe().EQ.198)
!!!#define GDCHK1 (mpp_pe().EQ.198).AND.(jrow.EQ.52)
!!!#define GDCHK2  (mpp_pe().EQ.198).AND.(jl.EQ.1).AND.(jrow.EQ.52)
!!!#define GDCHK2 ( (plon(jl).GT.-160.).AND.(plon(jl).LE.-159.).AND.(plat(jl).LE.8.).AND.(plat(jl).GE.7.) )

!!!#define GDCHK0 (mpp_pe().EQ.480)
!!!#define GDCHK1 (mpp_pe().EQ.480).AND.(jrow.EQ.113)
!!!#define GDCHK2 (mpp_pe().EQ.480).AND.(jrow.EQ.113).AND.(jl.EQ.1)

!!!#define GDCHK2 ( (plon(jl).GT.-1.).AND.(plon(jl).LE.-0.).AND.(plat(jl).GE.-46).AND.(plat(jl).LE.-45) )
!#define GDCHK0 (mpp_pe().EQ.171)
!#define GDCHK1 (mpp_pe().EQ.171).AND.(jrow.EQ.37)
!!!#define GDCHK2 (mpp_pe().EQ.171).AND.(jrow.EQ.37).AND.(jl.EQ.10)

!!!#define GDCHK2 ( (plon(jl).GT.-180.5).AND.(plon(jl).LE.-180.).AND.(plat(jl).GE.-5.0).AND.(plat(jl).LE.-4.5) )
!ps
!#define GDCHK0 (mpp_pe().EQ.54)
!#define GDCHK1 (1 .eq. 0).AND.(mpp_pe().EQ.54).AND.(jrow.EQ.5)
!#define GDCHK2 (1 .eq. 0).AND.(mpp_pe().EQ.54).AND.(jrow.EQ.5).AND.(jl.EQ.556)
!#define GDCHK3 (mpp_pe().EQ.54).AND.(jrow.EQ.5).AND.(jl.EQ.556)

#define GDCHK0 (mpp_pe().EQ.2)
#define GDCHK1 (1 .eq. 0).AND.(mpp_pe().EQ.44).AND.(jrow.EQ.3)
!#define GDCHK2 (1 .eq. 0).AND.(mpp_pe().EQ.44).AND.(jrow.EQ.3).AND.(jl.EQ.212)
#define GDCHK2 (1 .eq. 0).AND.(mpp_pe().EQ.254).AND.(jrow.EQ.3).AND.(jl.EQ.1)
!#define GDCHK3 (1 .eq. 1).AND.(mpp_pe().EQ.49).AND.(jrow.EQ.3).AND.(jl.EQ.212)
#define GDCHK3 (1 .eq. 0).AND.(mpp_pe().EQ.209).AND.(jrow.EQ.5).AND.(jl.EQ.12)
!#define GDCHK3 (1 .eq. 1).AND.(mpp_pe().EQ.21).AND.(jrow.EQ.4).AND.(jl.EQ.319)
!#define GDCHK3 (1 .eq. 1).AND.(mpp_pe().EQ.24).AND.(jrow.EQ.3).AND.(jl.EQ.188)
!#define GDCHK3 (mpp_pe().EQ.35).AND.(jrow.EQ.4).AND.(jl.EQ.495) 
!ps#define GDCHK3 ((mpp_pe().EQ.59).OR.(mpp_pe().EQ.60)) .AND. (jrow.EQ.2) .AND. (jl.EQ.151) 
!ps

!#define GDCHK0 (mpp_pe().EQ.741)
!#define GDCHK1 (mpp_pe().EQ.741)
!#define GDCHK2 (mpp_pe().EQ.741)

!#define DBLE REAL   ! EVERYTHING IN REAL
! EVERYTHING IN REAL
#define DBLE(x) REAL(x)

MODULE mod_sit_vdiff

! 1. History
!
!     v0.1: July, 1998,  Ben-Jei Tsuang
!     v0.76g: 18 July 2001, Ben-Jei Tsuang
!     1) correct the bug of molecular heat diffusivity (1.34E-7 m2/s)
!     2) new lwaterlevel logic for fixed water level run
!     3) set heat diffusivity within the surface viscous sublayer to be
!        the molecular heat diffusivity
!     v0.76i: 18 July 2001, Ben-Jei Tsuang
!     1) modify diffusivity of skin layer to be the geometric mean between molecular diff and XK(0)
!     2) modify radiation parameterization of skin layer of thickness "d"
!     v0.76j: 18 July 2001, Ben-Jei Tsuang
!     1) reset the diffusivity of skin layer to be XK(0)
!     v0.76k: 18 July 2001, Ben-Jei Tsuang
!     1) set radiation param. to be dFFN
!     v0.76l: 18 July 2001, Ben-Jei Tsuang
!     1) reset radiation param. to be the thickness of effective thickness
!     v0.76i: 18 July 2001, Ben-Jei Tsuang
!     1) reset radiation param and skin layer diff to be that in v.76i
!     v0.76m: 24 July 2001, Ben-Jei Tsuang
!     1) modify DTDZ to DRHODZ for stabiliy criterion
!     v0.76n: 25 July 2001, Chia-Ying Tu
!     1) ambient diffusion coefficeint is added
!     v0.76o: 25 July 2001, Chia-Ying Tu
!     1) minimum 0.3 m of mixing length is chosen
!     v0.77: 1 Aug 2001,
!     1) calibrate Prandtl number
!     v0.78: emin, EResRatio
!     v0.90: 9 Jan 2003, Chia-Ying Tu
!     1) add mo_thmcln
!     2) xkmmin xkhmin hcoolskin hcoolskin emin are set from arguments
!     3) use Inverse Problem to determine above values
!     v1.00: 29 July 2007, Ben-Jei Tsuang
!     1) port to ECHMA-5.4.00
!     2) modify heice and hesn
!     v1.30: 1 September 2008, Ben-Jei Tsuang
!     1) Modify Water Level Logical for lsoil
!     v3.2: 26 August 2008, Ben-Jei Tsuang
!     1) modify coordinate system according to Tsuang et al. (2009)
!        "A more accurate scheme for calcualting Earth's skin temperature",
!        Climate Dynamics
!     v5.0: June 2009, Ben-Jei Tsuang
!     1) Add below ice surface heat/fresh water fluxes for coupling with embedded ocean model
!     2) Add sufsurface (at about 10-20 m depth) heat/salinity fluxes for coupling with embedded ocean model
!     3) Turn on the ice module with security number for coupling with embedded ocean model
!     v5.1: June 2009, Ben-Jei Tsuang
!     1) A semi-implicity scheme for snow/ice temperature (TSI) coupling with vdiff (atmosphere)
!     v5.3: July 2009, Ben-Jei Tsuang
!     1) A new density function of pressure, temperature and salinity is needed for
!        vertical diffusion calculation
!     2) Bug fixed for density function
!     v7.7: snow/ice temperate set to be the underneach temp if these layers are missing. (bjt,2010/5)
!     v7.9 (bjt, 2010/5/29)
!     1) sit_vdiff.f90: seaice mask: Winner wins!
!     v8.7 (bjt, 2011/8/29)
!      Assunimg no skin layer for momentm since wind shear can enforce on the side. Note that the surface heat flux is from the top. (2011/8/29 bjt)
!     v9.8b          zcor=MAX(ABS(pcoriol(jl)),zepcor)  ! v9.8b: 20130822  (restore the corlios security number)
!     v9.86: bug correction: initialize pwtke
!     v9.84 (Same as v9.83 + xlkmin=0.) (bjt, 2013/9/10)
!     v9.85: salti  (seaice is salty at salti PSU) (bjt, 2013/9/11)
!     v9.86: emin  (truncate tke to 0, for tke < emin, 10-6 ms/s2) (bjt, 2013/9/13)
!     v9.862: rhom(jk)=rho_from_theta(wsm(jk),wtm(jk)-tmelt,0.) (bjt, 2013/9/15), change to potential water density
!     v9.864: CALL nudging_sit_viff_gd(jl,jrow,six_hour,six_hour,six_hour,six_hour,six_hour,.TRUE.) (bjt, 2013/9/16), nudging within 10 m depth as well    
!     v9.865: v9.864+no truncate low tke+no contrain in huge tke, emin=1.0E-6 (bjt, 2013/9/18)
!     v9.866: emin=1.0E-4 !limit min pwtke to 1.0E-4 (v9.866)
!     v9.867: emin=1.0E-5 !limit min pwtke to 1.0E-5 (v9.867)
!     v9.868: emin=1.0E-4 !limit min pwtke to 1.0E-4 (v9.868)
!     v9.869: CALL nudging_sit_viff_gd(jl,jrow,six_hour,six_hour,one_month,six_hour,six_hour,.TRUE.) (bjt, 2013/9/16), nudging > 100 m depth at 1-month time scale
!     v9.871: Add d0 for reaching the bottom d0=0.03  ! zero-displacement (m) (0.02,0.05)
!     v9.873: emin=1.0E-5 !limit min pwtke to 1.0E-5 (v9.873)
!     v9.874: change drho/dz from potentail density at surface to the potential density diff at the level (v9.874)
!             real,PARAMETER::emin=1.0E-4 !limit min pwtke to 1.0E-4 (v9.866),(v9.874)
!     v9.880: 1) bug found for ! bug, v0.9879
! &             -0.5_dp*ce*zdtime*wtkem(level)**(1.5_dp)/pwldisp(jl,level)                   &   ! bug, v0.9879
!             2) emin=1.0E-5 !limit min pwtke to 1.0E-5 (v9.873)
!     v9.884: bug correction for RHOE
!     v9.886: 1)bug correction for pwkm=ck*lk*e^0.5, then pwkh=pwkm/Pr
!             2) introducing steady TKE
!     v9.892: 1) with prho1000 output
!     v9.893: 1) nudging for deep salinity for coupled model during spinup period
!     v9.8999:
!     v10.4:  1) partial open water and partial ice sheet are revised by changing  
!                0.5*hsn to 0.5*hsn/seaice 
!                0.5*hice to 0.5*hice/seaice 
!     v10.601:  1) module sit_vdiff
!                2) steady-TKE  
!     v10.602:  1) friction is added in one-column model (fr=7.E-5/s)
!                2) unsteady-TKE
!
! 11. REFERENCE
!     Gaspar, P., Y. Gregoris, and J.-M. Lefevre,1990: A simple eddy
!         kinetic energy model for simulation of the oceanic vertical
!         mixing: test at station Papa and long-term upper ocean study
!         site. JGR, 95, 16179-16193.
!     Tu, C.-Y.; Tsuang, B.-J., 2005/11: Cool-skin simulation by a one-column
!         ocean model, Geophys. Res. Lett., 32, L22602, doi:10.1029/2005GL024252.
!     Tsuang, B.-J., C.-Y. Tu, J.-L. Tsai, J.A. Dracup, K. Arpe and T. Meyers, 2009:
!         A more accurate scheme for calculating Earth's skin temperature.
!         Climate Dynamics, DOI 10.1007/s00382-008-0479-2, on-line version available.
!
! 12. Bug known
!     Coriols force needed to be corrected for considering the curvature of the Earth
!  ---------------------------------------------------------------------
!

  !
  !  variables for sit coupling only
  !
!!!  USE mo_memory_g3b,     ONLY:                 &
!!!              obsseaice, obstsw, obswsb,       &
!!!              sitmask,                         &
!!!              bathy, ctfreez2, sitwlvl,        &
!!!              ocnmask,obox_mask,               &
!!!              sitwtb, sitwub, sitwvb,          &
!!!              sitwsb, fluxiw, pme2,            &
!!!              subfluxw, wsubsal,               &
!!!              cc,hc,engwac,                    &
!!!              sc,saltwac,                      &
!!!         !              
!!!              obswt,  obsws,                   &            
!!!              obswu,  obswu,                   &            
!!!              sitzsi,  sitsilw, sittsi,        &
!!!              sitwt,   sitwu,   sitwv,         &
!!!              sitww,                           &
!!!              sitws, sitwtke, sitwlmx,         & 
!!!              sitwldisp, sitwkm, sitwkh,       &
!!!              wrho,                            &
!!!              wtfn, wsfn,                      &
!!!              wtfns, wsfns,                    &
!!!              awufl,awvfl,awtfl,               &
!!!              awsfl,awtkefl

!ps  use mpp_mod, only: mpp_init, mpp_pe, mpp_root_pe, &
!ps                   stdlog, mpp_error, NOTE, FATAL, WARNING
      USE mod_eos_ocean,      ONLY: tmaxden,tmelts,rho_from_theta,theta_from_t
      USE mod_eos_ocean,      ONLY: api,argas,avo,ak,stbo,amco2,amch4,amo3 &
                                   ,amn2o,amc11,amc12,amw,amd,cpd,cpv,rd,rv &
                                  ,rcpd,vtmpc1,vtmpc2
      USE mod_eos_ocean,      ONLY: rhoh2o,alv,als,alf,clw,tmelt
      USE mod_eos_ocean,      ONLY: a,omega,g,IDAYLEN
      USE mod_eos_ocean,      ONLY: c1es,c2es,c3les,c3ies,c4les,c4ies,c5les,c5ies,c5alvcp,c5alscp,alvdcp,alsdcp
      USE mo_netcdf
      USE mod_sit_control,    ONLY: lssst,lsit_ice,lsit_salt,lrere,lhd       &
                               ,lsitstart,sit_ice_option,locaf,locaf0,ocaf0_add,lgodas,locn  &
                               ,lsice_nudg,lsit_lw,xmissing                               &
                               ,ssit_restore_time,usit_restore_time,dsit_restore_time     &
                               ,ssits_restore_time,usits_restore_time,dsits_restore_time     &
                               ,ssituv_restore_time,usituv_restore_time,dsituv_restore_time     &
                               ,lwarning_msg,ocn_couple_option                            &
                               ,socn_restore_time,uocn_restore_time,docn_restore_time     &
                               ,obox_restore_time,obox_nudg_flag,ck,ce,Prw,d0,csiced,salti       &
                               ,sit_domain_w,sit_domain_e,sit_domain_s,sit_domain_n,sit_domain_extgrd &
                               ,maskid,lobs_ocn_rerun,lwoa0,ltrigsit,lsteady_TKE,lwave_breaking        &
                               ,lcool_skin,lamip,lmixedlayer,ltimeblending,timebl_option,timebl_start  &
                               ,timebl_allsit,nsit_nudg                          &
                               ,sitbox_nudg_w,sitbox_nudg_e,sitbox_nudg_s,sitbox_nudg_n &
                               ,sitbox_st_restore_time,sitbox_ut_restore_time,sitbox_dt_restore_time &
                               ,sitbox_ss_restore_time,sitbox_us_restore_time,sitbox_ds_restore_time &
                               ,sitbox_suv_restore_time,sitbox_uuv_restore_time,sitbox_duv_restore_time &
                               ,lsftobswt
                              !nn,lsit,
      USE mod_sst,          ONLY: csn,rhosn,xksn,cice,rhoice,xkice,xkw,&
                               omegas,wcri,lou,lov  &
                              ,nodepth,odepths,ot12,os12,ou12,ov12,mixedlayer12   &
                              ,nodepth0,odepth0,ot0,os0,ou0,ov0,mixedlayer0  &
                              ,nwdepth,wdepths,wtfn12, wsfn12,mask1st &
                              ,wgt1,wgt2,nmw1,nmw2,obswtbwgt1,obswtbwgt2,obswtbnmw1,obswtbnmw2 &
                              ,now1,now2,wgto1,wgto2   
!ps							   ,wlvlref,dpthmx
                                   
#if defined (LCWBGFS)
      use rank      ! myrank
#endif

                                   
  IMPLICIT NONE




!--------------------------------------------------------------------------
!  USE convect_tables_mod
!-----------------------

  ! Lookup tables for convective adjustment code
  !
  ! D. Salmond, CRAY (UK), August 1991, original code
  ! L. Kornblueh, MPI, April 2003, cleanup and move of table setup code
  !                                from setphys.f90 in module  
  !

  !!! USE mo_kind,   ONLY: dp

  !!! IMPLICIT NONE
  !!! ! Standard i-o units
  !!! 
  !!! INTEGER, PARAMETER :: nout = 0     ! standard output stream
  !!! INTEGER, PARAMETER :: nerr = 0     ! error output stream
  !!! INTEGER, PARAMETER :: nin  = 5     ! standard input stream
  !!! 
  !!! INTEGER, PARAMETER:: dp=4
  !!! 
  !!! SAVE
  !!! 
  !!! PRIVATE

  ! variables public

  PUBLIC :: jptlucu1          ! lookup table lower bound
  PUBLIC :: jptlucu2          ! lookup table upper bound
  PUBLIC :: tlucua            ! table -- e_s*Rd/Rv
  PUBLIC :: tlucub            ! table -- for derivative calculation: d es/ d t
  PUBLIC :: tlucuc            ! table -- l/cp
  PUBLIC :: tlucuaw           ! table
  PUBLIC :: lookupoverflow    ! lookup table overflow flag

  ! subroutines public

  PUBLIC :: init_convect_tables ! initialize LUTs 
  PUBLIC :: lookuperror         ! error handling routine 

  INTEGER, PARAMETER :: jptlucu1 =  50000  ! lookup table lower bound
  INTEGER, PARAMETER :: jptlucu2 = 400000  ! lookup table upper bound

  LOGICAL :: lookupoverflow = .FALSE.          ! preset with false
  
  real :: tlucua(jptlucu1:jptlucu2)        ! table - e_s*Rd/Rv
  real :: tlucub(jptlucu1:jptlucu2)        ! table - for derivative calculation
  real :: tlucuc(jptlucu1:jptlucu2)        ! table - l/cp
  real :: tlucuaw(jptlucu1:jptlucu2)       ! table

!--------------------------------------------------------------------------
!   USE mo_doctor,         ONLY: nout, nin, nerr
!-----------------------

  ! Standard i-o units

!!!  INTEGER, PARAMETER :: nout = 0     ! standard output stream 
!!  INTEGER, PARAMETER :: nerr = 0     ! error output stream
  INTEGER, PARAMETER :: nin  = 5     ! standard input stream
  INTEGER, PARAMETER :: nout = 6     ! standard output stream 
  INTEGER, PARAMETER :: nerr = 6     ! error output stream


!--------------------------------------------------------------------------
! from mo_exception
!-----------------------

  !!! PUBLIC :: message_text
  !!! PUBLIC :: message, finish
  !!! PUBLIC :: em_none, em_info, em_warn

  INTEGER, PARAMETER :: em_none = 0 
  INTEGER, PARAMETER :: em_info = 1
  INTEGER, PARAMETER :: em_warn = 2

!!!  CHARACTER(512) :: message_text = ''


!--------------------------------------------------------------------------
! from mo_interpro
!-----------------------

!
!    WEIGHTNG FACTORS AND MONTH INDICES FOR
!    INTERPOLATION IN TIME IN *CLSST* AND *RADINT*
!

!ps     real:: WGT1=0.5
!ps     real:: WGT2=0.5
!ps     INTEGER :: NMW1=1
!ps     INTEGER :: NMW2=1
     INTEGER :: NMW1CL=1
     INTEGER :: NMW2CL=1
     real:: WGTD1=0.5
     real:: WGTD2=0.5
     INTEGER :: NDW1=1
     INTEGER :: NDW2=1
!--------------------------------------------------------------------------
! from mo_mpi
!-----------------------
  ! PE identifier

  !!! USE mo_mpi,            ONLY: p_io, p_pe

  !!! PUBLIC :: p_pe, p_io, p_nprocs, p_ocean
  
  ! logical switches

  !!! PUBLIC :: p_parallel
  !!! 
  !!! INTEGER :: p_pe=0                ! this is the PE number of this task
  !!! INTEGER :: p_io=0                ! PE number of PE handling IO
  !!! INTEGER :: p_ocean               ! PE number of PE handling ocean model
  !!! INTEGER :: p_nprocs              ! number of available PEs (processors)
  
  ! public parallel run information

  LOGICAL :: p_parallel = .TRUE.
  
!--------------------------------------------------------------------------
! from mo_parameters
!-----------------------

!----------------------------------------------------------------
!   USE mo_physc2,         ONLY: wicemx, csncri, xicri  
!-----------------------

  real:: wicemx = 0.010               ! maximum water (assume 0.01 m) on the top of a ice layer
  real:: csncri = 5.85036E-3          !  CRITICAL SNOW DEPTH FOR SOIL COMPUTATIONS
  real:: xicri  = 0.005               ! critial depth of ice, causing lost of ice heat flux into water)
  real:: ctfreez= 271.38                        !   temperature at which sea
                                                       !   starts freezing/melting
!----------------------------------------------------------------
! from mo_semi_impl
!-----------------------
  real:: eps  = 0.1 !   time filtering coefficient.

!  USE mo_gaussgrid, ONLY: gl_coriol

!!!  USE mo_geoloc,         ONLY: coriol_2d, philat_2d, philon_2d

!!!  USE ice_model_mod,     ONLY: ice_data_type


  !!! PUBLIC :: sit_vdiff_init,sit_vdiff_ik_init,sit_vdiff,sit_vdiff_end,SICEDFN
  PUBLIC :: sit_vdiff_init,sit_vdiff,sit_vdiff_end,SICEDFN  
!          ,maskid,nwdepth, wdepths, wtfn12, wsfn12, locaf ,lstart  &
!          ,lgodas,sit_nml,nodepth,odepths,ot12,os12,ou12,ov12,lou,lov     &
!          ,lwoa0,nodepth0,odepth0,ot0,os0,ou0,ov0,lwarning_msg    &
!          ,sit_domain_w,sit_domain_e,sit_domain_s,sit_domain_n,sit_domain_extgrd,lamip
            
  PUBLIC :: cal_ratioBlending
 

  !!! real, PARAMETER::tol=1.E-14_dp  ! tol: very small numerical value to prevent numerical error
  real, PARAMETER:: tol=1.E-12  ! tol: very small numerical value to prevent numerical error
                                      ! =1.E-6 (ok), tol=1.E-12 (ok), =1.E-33 (crash), , =1.E-20 (crash)
                                      ! =1.E-16_dp (crash)
                                      ! =1.E-14_dp (ok)
                                      ! =0. (ok)

!!!  real,PARAMETER::emin=0. !limit min pwtke to 0. to avoid too warm below mixing depth (20090710) (v5.4)(v9.913)
  real,PARAMETER:: emin=1.0E-6 !limit min pwtke to 1.0E-6 (v9.865, v9.885) (1.0E-6 in GASPAR ET AL., 1990)
!  real,PARAMETER:: emin=1.0E-4 !limit min pwtke to 1.0E-4 (v9.866),(v9.868)
!  real,PARAMETER:: emin=1.0E-5 !limit min pwtke to 1.0E-5 (v9.867) (v9.873)
  real,PARAMETER:: xkmmin=1.2E-6     ! molecular momentum diffusivity (Paulson and Simpson, 1981; Chia and pwu, 1998; Mellor and Durbin, 1975)
  real,PARAMETER:: xkhmin=1.34E-7    ! molecular heat diffusivity (Paulson and Simpson, 1981; Chia and pwu, 1998; Mellor and Durbin, 1975)
!ps  real,PARAMETER:: e_mixing=1.0E-1      ! default ini TKE within mixing depth (1.0E-3=still some spikes)  
  real,PARAMETER:: e_mixing=1.0E-4      ! default ini TKE within mixing depth (1.0E-3=still some spikes)  
  real,PARAMETER:: wlmx_mixing=50.      ! default ini wlmx within mixing depth (50 m)
  real,PARAMETER:: pwldisp_mixing=50.   ! default ini wldisp within mixing depth (50 m)
  real,PARAMETER:: xkm_mixing=1.0E-1    ! default ini km within mixing depth (1E-3 m2/s=still some spikes)
  real,PARAMETER:: xkh_mixing=1.0E-1    ! default ini kh within mixing depth  (1E-3 m2/s=still some spikes)


! --------------------------------
!
!     2.1 OCEAN PARAMETERS
!
  real,PARAMETER:: rhowcw = rhoh2o*clw
  INTEGER :: wtype = 0
! wtype: water type of the ocean
! wtype = 0, PAULSON AND SIMPSON (1981),  Fairall et al., (1996)
!
 ! WTYPE = 22, two-components, JERLOV'S (1976)
  ! WTYPE = 0, PAULSON AND SIMPSON (1981),  Fairall et al., (1996)
  ! WTYPE = 10, SOLOVIEV AND SCHLUESSEL (1996,type I)
  ! WTYPE = 11, SOLOVIEV AND SCHLUESSEL (1996,type IA)
  ! WTYPE = 12, SOLOVIEV AND SCHLUESSEL (1996,type IB)
  ! WTYPE = 20, SOLOVIEV AND SCHLUESSEL (1996,type II)
  ! WTYPE = 30, SOLOVIEV AND SCHLUESSEL (1996,type III)
  ! WTYPE = 1, SOLOVIEV AND SCHLUESSEL (1996,type 1)
  ! WTYPE = 3, SOLOVIEV AND SCHLUESSEL (1996,type 3)
  ! WTYPE = 5, SOLOVIEV AND SCHLUESSEL (1996,type 5)
  ! WTYPE = 7, SOLOVIEV AND SCHLUESSEL (1996,type 7)
  ! WTYPE = 9, SOLOVIEV AND SCHLUESSEL (1996,type 9)

  real:: mixing_depth=100.          ! mixing depth of ocean, set at 100 m, should be changed as a function of time and location  
!**************************************************	
!    2) SOIL PARAMETERS
!
  real,PARAMETER::xkg = 0.4835E-6
!     heat capacity of soil
!       =(0.5675E-6-0.175E-6*porosity)*soil type factor (de Vries, 1975)
!       =0.4835E-6 ,
!     sandy clay loam at water content , porosity(48%), factor =1
!     CGSOIL=248672.
  real,PARAMETER::rhogcg=1.04E6+0.48*4.19E6
!     area heat capacity of the above soil (Tsuang and Wang, 1994)
!       =rhog*cg*Sqrt(Kg/omega)
!       =(1.04E6+0.48*4.19E6)*0.0815
!     ALBG = 0.3
!     ALBG: ALBEDO OF UNDERNEATH SOIL
!
  real,PARAMETER::hspg=0.         ! HOT SPRING FROM THE BOTTOM (W/M2)                                I
!
!!!  real, PARAMETER:: nudge_depth_10=0.   ! (v8.3 - v9.1)
  real, PARAMETER:: nudge_depth_10=10.     ! (- v8.3 and v9.2 -)
  real, PARAMETER:: nudge_depth_100=100.
!ps
  
!ps
CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE sit_vdiff_init ( kproma, kbdim, jrow, istep,               &
       !  
       ! 1-input only, original ATM/SIT variabels
       !
                  plat,       plon,                                   &
                  psitmask,   pbathy,     pwlvl,          &
                  pocnmask,   obox_mask,                              &
                  psni,       psiced,     ptsi,                       &
                  pobsseaice, pobswtb,    pobswsb,                    &
                  pctfreez2,                                          &
       ! 2-d SIT vars
                  pwtb,       pwub,       pwvb,                       &
                  pwsb,                                               &
                  psubfluxw,  pwsubsal,                               &
                  pcc,        phc,        pengwac,                    &
                  psc,        psaltwac,                               &
                  pwtfns,     pwsfns,                                 &
       ! 3-d SIT vars: snow/ice
                  pzsi,       psilw,      ptsnic,                     &
       ! 3-d SIT vars: water column
                  pobswt,     pobsws,     pobswu,    pobswv,          &
                  pwt,        pwu,        pwv,       pww,             &
                  pws,        pwtke,      pwlmx,                      & 
                  pwldisp,    pwkm,       pwkh,                       &
                  pwrho1000,                                          &
                  pwtfn,      pwsfn,      pwtfn0,    pwsfn0,          &
                  pawufl,     pawvfl,     pawtfl,                     &
                  pawsfl,     pawtfl0,    pawsfl0,                    &
                  pawtkefl,                                           &
       ! 4-output only, original ATM variabels
                  pseaice,                                            &
                  pgrndcapc,  pgrndhflx,  pgrndflux,psftobswt )
!  ---------------------------------------------------------------------
!
!                           SUBROUTINE DESCRIPTION
!
! 1. FUNCTION DESCRIPTION
!
!     THIS SUBROUTINE initialize the sit_vdiff variables including T,S,u,v profiles
!
! 2. CALLING MODULES
!
!          *sit_vdiff_init* is called from *initemp.
!
! 3. PARAMETER SPECIFICATION
!
!      Arguments.
!      ----------
! local dimensions
!  kproma   : IF(jrow==ldc%ngpblks, ldc%npromz,ldc% nproma),
!             number of local longitudes
!  kbdim    : ldc% nproma, gauss grid description
!  jrow     : sequential index for output message
!  istep    : # of time steps
!  pcoriol : coriol factor = 2*omega*sin(latitudes)  (1/s)                         I
!  plat    : latitide (geg)                                                        I
!  plon    : longitude (deg)                                                       I
! 2-d SIT vars
!  pobsseaice  : observed sea ice fraction (fraction)                              I
!  pobswtb     : observed bulk sea surface temperature (K)                         I
!  pobswsb    : observed salinity (PSU, 0/00)                                      I
!  psitmask : mask for sit(1=.TRUE., 0=.FALSE.)                                    I
!  pbathy  : bathymeter (topography or orography) of ocean (m)                     I
!  pctfreez2 : ref water freezing temperature (K)                                  I
!  pwlvl  : current water level (ice/water interface) a water body grid            I/O
!  pocnmask : fractional mask for 3-D ocean grid (0-1)                             I
!  obox_mask : 3-D ocean nudging mask, =0: nudging, = 1 (>0): nudging              I
!  pwtb: bulk water temperature (K)                                                O
!  pwub: bulk water u current (m/s)                                                O
!  pwvb: bulk water v current (m/s)                                                O
!  pwsb: bulk water salinity (PSU)                                                 O
!  psubfluxw: subsurface ocean heat flux (W/m2, + upward)                          O
!  pwsubsal: subsurface ocean salinity flux (m*PSU/s, + upward)                    O
!  pcc: cold content per water fraction (ice sheet+openwater) (J/m2)
!    (energy need to melt snow and ice, i.e., energy below liquid water at tmelt)  I/O
!  phc: heat content per water fraction (ice sheet+openwater) (J/m2)
!    (energy of a water column above tmelt)                                        I/O
!  pengwac: accumulated energy per water fraction (ice sheet+openwater) (+ downward, J/m2)
!    (pfluxw+pfluxi+rain/snow advected energy in respect to liquid water at
!     tmelt)*dt                                                                    I/O
!  psaltwac: accumulated salt into water fraction (+ downward, PSU*m)              I/O
!
!  pgrndcapc: areal heat capacity of the uppermost sit layer (snow/ice/water)
!    (J/m**2/K)                                                                    I/O
!  pgrndhflx: ground heat flux below the surface (W/m**2)
!    (+ upward, into the skin layer)                                               I/O
!  pgrndflx:  acc. ground heat flux below the surface (W/m**2*s)
!    (+ upward, into the skin layer)                                               I/O
! 3-d SIT vars: snow/ice
!  pzsi   :
!        pzsi(jl,0): dry snow water equivalent (m)                                 I/O
!        pzsi(jl,1): dry ice water equivalent (m)                                  I/O
!  psilw  :
!        psilw(jl,0): snow liquid water storage (m)                                I/O
!        psilw(jl,1): ice liquid water storage (m)                                 I/O
!  ptsnic   :
!        ptsnic(jl,0): snow skin temperatrue (jk)                                  I/O
!        ptsnic(jl,1): mean snow temperatrue (jk)                                  I/O
!        ptsnic(jl,2): ice skin temperatrue (jk)                                   I/O
!        ptsnic(jl,3): mean ice temperatrue (jk)                                   I/O
! 3-d SIT vars: water column
!  pobswt: observed potentail water temperature (K)                                I/O
!  pobsws: observed salinity (PSU, 0/00)                                           I/O
!  pobswu: observed u-current (m/s)                                                I/O
!  pobswv: observed v-current (m/s)                                                I/O
!  pwt      : potential water temperature (K)                                      I/O
!  pwu      : u current (m/s)                                                      I/O
!  pwv      : v current (m/s)                                                      I/O
!  pww      : w current (m/s)                                                      I/O
!  pws      : practical salinity (0/00)                                            I/O
!  pwtke    : turbulent kinetic energy (M2/S2)                                     I/O
!  pwlmx    : mixing length (m)                                                    O
!  pwldisp  : dissipation length (m)                                               O
!  pwkm     : eddy diffusivity for momentum (m2/s)                                 O
!  pwkh     : eddy diffusivity for heat (m2/s)                                     O
!  pwrho1000: potential temperature at 1000 m (PSU)                                O
!  pwtfn  : nudging flux into sit-ocean at each level (K, accum. var.)             I/O
!  pwsfn  : nudging salinity flux into sit-ocean at each level (PSU, accum. var.)  I/O
!  pwtfns : nudging flux into sit-ocean for an entire column (J/m**2, accum.var.), i.e.
!     pwtfns=SUM(pwtfn(jl,:))                                                      I/O
!  pwsfn  : nudging salinity flux into sit-ocean for an entire column (PSU*m,accum. var.,)
!     i.e., pwsfns=SUM(pwsfn(jl,:)                                                 I/O
!  pwtfn0  : nudging flux into sit-ocean at each level (K, accum. var.)             I/O
!  pwsfn0  : nudging salinity flux into sit-ocean at each level (PSU, accum. var.)  I/O
!  pawufl : advected u flux at each level (positive into ocean) (m/s/s)            I
!  pawvfl : advected v flux at each level (positive into ocean) (m/s/s)            I
!  pawtfl : advected temperature flux at each level (positive into ocean) (K/s)(12 mons)    I
!  pawsfl : advected salinity flux at each level (positive into ocean) (PSU/s)(12 mons)     I
!  pawtfl0 : advected temperature flux at each level (positive into ocean) (K/s)(1 record)    I
!  pawsfl0 : advected salinity flux at each level (positive into ocean) (PSU/s)(1 record)     I
!  pawtkefl : advected tke at each level (positive into ocean) (m3/s3)             I
! - variables internal to physics
! - water body or ml_ocean variables
!  pseaice   : ice cover (fraction of 1-SLM) (0-1)                                 I/O
!  psni      : snow thickness over ice (m in water equivalent)                     I/O
!  psiced    : ice thickness (m in water equivalent)                               I/O
!  ptsl      : calcuated Earth's skin temperature from vdiff/tsurf at t+dt (K)     I/O
!  ptslm     : calcuated Earth's skin temperature from vdiff/tsurf at t (K)        I/O
!  ptslm1    : calcuated Earth's skin temperature from vdiff/tsurf at t-dt (K)     I/O
!  psftobswt : the shift swt in mixed layer based on difference of obswtb and obswt10m I/O
! 2-d SIT vars

  IMPLICIT NONE

  ! local dimensions
  INTEGER, INTENT(in)                      :: kproma ! number of local longitudes
  INTEGER, INTENT(in)                      :: kbdim  ! number of local longitudes

  ! gauss grid description
  INTEGER, INTENT(in)                      :: jrow   ! sequential index
  INTEGER, INTENT(in)                      :: istep  ! # of time steps
  real, INTENT(in)::                                                      &
       plat(kbdim), plon(kbdim)   
       ! pcoriol(kbdim),   plat(kbdim), plon(kbdim)   

  real, INTENT(in out):: pseaice(kbdim), psni(kbdim), psiced(kbdim),     &
     ptsi(kbdim)
! 2-d SIT vars                                                
  real, INTENT(in out) ::                                                &
       pwtfn(kbdim,0:lkvl+1),  pwsfn(kbdim,0:lkvl+1)
  real, INTENT(in out) ::                                                &
       pwtfn0(kbdim,0:lkvl+1),  pwsfn0(kbdim,0:lkvl+1)
  real, INTENT(in out) ::                                                &
       pwtfns(kbdim),  pwsfns(kbdim)
  real, INTENT(in):: pobsseaice(kbdim)
  real, INTENT(in out):: pobswtb(kbdim),      pobswsb(kbdim)
! 2-d SIT vars                                                
  real, INTENT(in)::   psitmask(kbdim)         ! grid mask for lsit (1 or 0)      
  real, INTENT(in out) ::   pbathy(kbdim)
  real, INTENT(in out) :: pctfreez2(kbdim)
  real, INTENT(in out) :: pwlvl(kbdim)
  real, INTENT(in)::   pocnmask(kbdim)         ! (fractional) grid mask for 3-D ocean (DIECAST)
  real, INTENT(in)::   obox_mask(kbdim)        ! (fractional) ocean iop nudging mask for 3-D ocean (DIECAST)
  real, INTENT(out) :: pwtb(kbdim), pwub(kbdim), pwvb(kbdim), pwsb(kbdim), &
    psubfluxw(kbdim), pwsubsal(kbdim)
  real, INTENT(in out) :: pgrndcapc(kbdim), pgrndhflx(kbdim), pgrndflux(kbdim)
  real, INTENT(in out) :: pcc(kbdim), phc(kbdim), pengwac(kbdim)
  real, INTENT(in out) :: psc(kbdim), psaltwac(kbdim)
! 3-d SIT vars: snow/ice
  real, INTENT(in out) ::                                                &
       pzsi(kbdim,0:1),   psilw(kbdim,0:1), ptsnic(kbdim,0:3)
! 3-d SIT vars: water column
  real, INTENT(in out) ::                                                &
       pobswt(kbdim,0:lkvl+1), pobsws(kbdim,0:lkvl+1),                       &
       pobswu(kbdim,0:lkvl+1), pobswv(kbdim,0:lkvl+1),                       &
       pwt(kbdim,0:lkvl+1), pwu(kbdim,0:lkvl+1), pwv(kbdim,0:lkvl+1),        &
       pww(kbdim,0:lkvl+1),                                                  &
       pws(kbdim,0:lkvl+1), pwtke(kbdim,0:lkvl+1), pwlmx(kbdim,0:lkvl+1),    & 
       pwrho1000(kbdim,0:lkvl+1),                                                & 
       pwldisp(kbdim,0:lkvl+1), pwkm(kbdim,0:lkvl+1), pwkh(kbdim,0:lkvl+1)  

  real, INTENT(in out) :: pawufl(kbdim,0:lkvl+1), pawvfl(kbdim,0:lkvl+1), &
       pawtfl(kbdim,0:lkvl+1), pawsfl(kbdim,0:lkvl+1),              &
       pawtfl0(kbdim,0:lkvl+1), pawsfl0(kbdim,0:lkvl+1),            &
       pawtkefl(kbdim,0:lkvl+1)
  real, INTENT(in out) :: psftobswt(kbdim,0:lkvl+1)

! - variables internal to physics

! - local variables

  !*    1.0 Depth of the coordinates of the water body point. (zdepth = 0 at surface )
!!!  real, PARAMETER:: zdepth(0:lkvl)=(/0.,0.0005_dp,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,&
!!!              20.,30.,50.,75.,100.,125.,150.,200.,250.,300.,400.,500.,&
!!!              600.,700.,800.,900.,1000.,1100.,1200.,1300.,1400.,1500./) 
  !
  !********************
  ! Note that WOA 2005 data are at depths:
  !  depth = 0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600,
  !    700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500 ;
  !********************

  !*    1.1  INFORMATION OF OCEAN GRID (COMLKE)
  !     WATER BODY POINT VERTICAL LEVELS
  INTEGER nls ! first water level -1,
              ! starting water level of the water body point.
              ! z(k)=sitwlvl, IF k is <= nls. 
  INTEGER nle ! nle: last water level (excluding skin level)
              ! nle: maximum water levels of the water body point.
              ! Note nle=-1 IF water level is below point water body bed.
              ! Soil layer always is sitwt(nle+1). If there is no water,
              ! sitwt(nle+1)=sitwt(0,jrow). In addition, z(k)=bathy IF k >= nle.
  !!! INTEGER nlv ! number of water levels (excluding water surface &
  !!!             ! soil layer). Note nlv=-1 IF water level is below
  !!!             ! water bed.  
  real, DIMENSION(0:lkvl+1):: z   ! coordinates of the water body point.
  real, DIMENSION(0:lkvl+1):: zlk ! standard z coordinates of a water body 
  real, DIMENSION(0:lkvl):: hw    ! thickness of each layer
  LOGICAL lshlw  ! .true. = shallow water mode (due to evaportion or freezing)
              ! (one layer water body)
  LOGICAL lsoil  ! .true. = soil grid

  INTEGER:: lstfn     ! index of last fine level
  real:: heice  ! effective skin thickness of ice (m)
  real:: hice   ! ice thickness (m)
  real:: hesn   ! effective skin thickness of snow (m)
  real:: hsn    ! snow thickness (m)
  real:: hew    ! effective skin thickness of water (m)
  real:: xkhskin   ! skin layer heat diffusivity (m2/s)
!***************
! Local variables for backgroud initial ocean profiles
  real:: bg_wt0(kbdim,0:lkvl+1)
  real:: bg_ws0(kbdim,0:lkvl+1)
  real:: bg_wu0(kbdim,0:lkvl+1)
  real:: bg_wv0(kbdim,0:lkvl+1)
!!!!***************
  INTEGER :: jl=-999
  LOGICAL :: lsitmask(kbdim)   ! sit mask
!ps**************
  real,DIMENSION(kbdim):: mixedlayer  !mixed layer depth (m)
  INTEGER:: k10m,kmixdepth
  real:: t10m,t10m1,tmixup,tmixdown,mixlayer,rr

!ps*************


!!!  INTEGER :: istep=1                           ! istep=time step

!!!  RETURN

    lsitmask=psitmask.EQ.1.                          ! convert from REAL to Integer and to Logical
!!!    RETURN
    IF (GDCHK1) THEN
      WRITE(nerr,*) "I am in sit_vdiff_ik_init 1.0:"   
      WRITE(nerr,*) "lsitstart=.TRUE."
!!!!!!      WRITE(nerr,*) "delta_time=",delta_time,"zdtime=",zdtime,"n_sit_step=",n_sit_step
    ENDIF
#ifdef DEBUG
  if (GDCHK1) print *, "I am in sit_vdiff_ik_init 0.0."
  if (GDCHK1) print *, "lat=",plat
  if (GDCHK1) print *, "lon=",plon
!!  if (GDCHK1) print *, "tsw=",ptsw
  if (GDCHK1) print *, "sni=",psni          
  if (GDCHK1) print *, "siced=",psiced          
  if (GDCHK1) print *, "tsi=",ptsi          
  if (GDCHK1) print *, "obsseaice=",pobsseaice          
  if (GDCHK1) print *, "obswtb=",pobswtb          
  if (GDCHK1) print *, "obswsb=",pobswsb
  if (GDCHK1) print *, "ctfreez2=",pctfreez2  
  if (GDCHK1) print *, "seaice=",pseaice
#endif
    DO jl=1,kproma
      IF (GDCHK2) then
        WRITE(nerr,*) ", I am in sit_vdiff_init 1.0: prior to init_sit_viff_gd"
!        IF (lsitmask(jl) ) CALL output2
      ENDIF
      CALL init_sit_viff_gd(jl,jrow)
      IF (GDCHK2) then
        WRITE(nerr,*) ", I am in sit_vdiff_init 2.0: after init_sit_viff_gd"
        CALL output2
      ENDIF
    END DO
!!!    RETURN
!!!    masking RETURN will crash in one day
!!!    run for one time step for self-adjusting before entering to ocean
!!!    RETURN
  CONTAINS

  SUBROUTINE init_sit_viff_gd(jl,jrow)
!
! Initialise sit T,U,V,S on ocean levels (lake, ocean and soil)
!
  !!! USE mo_mpi,           ONLY: p_pe 
   USE mod_eos_ocean,     ONLY: tmelt
   USE mod_sit_control,       ONLY: lgodas,lwoa0
!!!  USE mo_physc2,        ONLY: ctfreez
                          
  IMPLICIT NONE
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  INTEGER  :: jk,kkk
  real:: depth  
  real:: sumxxz,sumxxt,sumxxu,sumxxv,sumxxs
  
  IF (GDCHK2) then
      WRITE(nerr,*) ", I am in init_sit_viff_gd"
  ENDIF

! initialize accumulated variables  
  pwtfns(jl)=0.
  pwsfns(jl)=0.
  pwtfn(jl,:)=0.
  pwsfn(jl,:)=0.
  pwtfn0(jl,:)=0.
  pwsfn0(jl,:)=0.
  pawufl(jl,:)=0.
  pawvfl(jl,:)=0.
  pawtfl(jl,:)=0.
  pawsfl(jl,:)=0.
  pawtfl0(jl,:)=0.
  pawsfl0(jl,:)=0.
  pawtkefl(jl,:)=0.
  psftobswt(jl,:)=0.
!  
  IF (GDCHK2) WRITE(nerr,*) ", I am in init_sit_viff_gd 1.1, lsitmask(",jl,")=",lsitmask(jl)
  IF (.NOT.lsitmask(jl)) THEN
    IF (.FALSE.) THEN    
      ptsnic(jl,:)=tmelt
      pzsi(jl,:)=0.
      
      pobswt(jl,:)=tmelt
      pobsws(jl,:)=0.
      pobswu(jl,:)=0.
      pobswv(jl,:)=0.
      
      pwt(jl,:)=tmelt
      pws(jl,:)=0.
      pwu(jl,:)=0.
      pwv(jl,:)=0.
      pww(jl,:)=0.
          
      pwtke(jl,:)=0.
      pwlmx(jl,:)=0.
      pwldisp(jl,:)=0.
      pwkm(jl,:)=0.
      pwkh(jl,:)=0.
      
      pgrndcapc(jl)=0.
      pgrndhflx(jl)=0.
      pgrndflux(jl)=0.
      pengwac(jl)=0.
      psaltwac(jl)=0.
      pwtb(jl)=tmelt
      pwub(jl)=0.
      pwvb(jl)=0.
      pwsb(jl)=0.
      psubfluxw(jl)=0.
      pwsubsal(jl)=0.
    ELSE
      ptsnic(jl,:)=xmissing
      pzsi(jl,:)=xmissing
      
      pobswt(jl,:)=xmissing
      pobsws(jl,:)=xmissing
      pobswu(jl,:)=xmissing
      pobswv(jl,:)=xmissing
      
      pwt(jl,:)=xmissing
      pws(jl,:)=xmissing
      pwu(jl,:)=xmissing
      pwv(jl,:)=xmissing    
      pww(jl,:)=xmissing
          
      pwtke(jl,:)=xmissing
      pwlmx(jl,:)=xmissing
      pwldisp(jl,:)=xmissing
      pwkm(jl,:)=xmissing
      pwkh(jl,:)=xmissing
      
      pgrndcapc(jl)=xmissing
      pgrndhflx(jl)=xmissing
      pgrndflux(jl)=xmissing
      pengwac(jl)=xmissing
      psaltwac(jl)=xmissing
      pwtb(jl)=xmissing
      pwub(jl)=xmissing
      pwvb(jl)=xmissing
      pwsb(jl)=xmissing
      psubfluxw(jl)=xmissing
      pwsubsal(jl)=xmissing
    ENDIF          
    RETURN
  ENDIF 
!!!  bjt >> 20170703
!     Make all water least 2 layers deep
  pbathy(jl)=MIN(pbathy(jl),pwlvl(jl)-sit_zdepth(2))
  CALL pzcord(pwlvl(jl),pbathy(jl),nls,nle,z,zlk,hw,lsoil,lshlw)
!!!  IF (lsoil) THEN
!!!!!!    pwlvl(jl)=pbathy(jl)+wcri+10.
!!!    pbathy(jl)=pwlvl(jl)-wcri-10.
!!!    !! add 10-m depth water to soil grid.
!!!    !! This 10-m water line should be deleted for more physcial sound,
!!!    CALL pzcord(pwlvl(jl),pbathy(jl),nls,nle,z,zlk,hw,lsoil,lshlw)
!!!  ENDIF
!!
!! Note that the vertical index of g3b starts from 1 although the index in the sit_vdiff 
!!   routine starts from 1,
  IF (ptsi(jl).NE.xmissing) THEN
    ptsnic(jl,:)=ptsi(jl)    ! assume snow temperature = ice temperature
  ELSE
    ptsnic(jl,:)=tmelt       ! assume snow temperature = tmelt
  ENDIF
  IF (GDCHK2) then
      WRITE(nerr,*) ", I am in init_sit_viff_gd 1.8: after setting ptsnic"
      !!! CALL output2
  ENDIF
  
  IF (psni(jl).NE.xmissing) THEN
    pzsi(jl,0)=psni(jl)            ! snow depth
  ELSE
    pzsi(jl,0)=0.
  ENDIF
  IF (psiced(jl).NE.xmissing) THEN
    pzsi(jl,1)=psiced(jl)         ! ice depth
  ELSE
    pzsi(jl,1)=0.
  ENDIF
  IF (GDCHK2) then
    WRITE(nerr,*) "istep=",istep,"pe=",mpp_pe(),"jrow=",jrow,", I am in init_sit_viff_gd 1.9: after setting zsi"
    !!! CALL output2
  ENDIF
!    sitsilw improper chosen might cause larger KM
!  runtoc(jl,jrow)=0.
!  hspg(jl,jrow)=0.
! 
  IF (lobs_ocn_rerun) THEN
    ! modify prog. vars (wt, ws, wu, wv) from rerun file according to obs
    CALL compose_obs_ocean(jl,jrow,.TRUE.,.TRUE.,.TRUE.,.TRUE.)
  ELSE
    psilw(jl,:)=0.                 ! assume no liquid water in snow and ice
    CALL compose_obs_ocean(jl,jrow,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
    IF (GDCHK2) then
      WRITE(nerr,*) ", I am in init_sit_viff_gd 2.1: lobs_ocn_rerun=F"
      CALL output2
    ENDIF
    pwtke(jl,0:nle+1)=emin
    pwlmx(jl,0:nle+1)=0.
    pwldisp(jl,0:nle+1)=0.
    pwkm(jl,0:nle+1)=xkmmin
    pwkh(jl,0:nle+1)=xkhmin
    DO jk=0,nle+1
      depth=(pwlvl(jl)-z(jk))
      IF ( depth.LE.mixing_depth) then
        pwtke(jl,jk)=e_mixing
!ps        pwlmx(jl,jk)=wlmx_mixing
!ps        pwldisp(jl,jk)=pwldisp_mixing
!ps        pwkm(jl,jk)=xkm_mixing
!ps        pwkh(jl,jk)=xkh_mixing
      ELSE
        pwtke(jl,jk)=emin
!ps        pwlmx(jl,jk)=0.
!ps        pwldisp(jl,jk)=0.
!ps        pwkm(jl,jk)=xkmmin
!ps        pwkh(jl,jk)=xkhmin
      ENDIF
    ENDDO
    IF (GDCHK2) then
      WRITE(nerr,*) ", I am in init_sit_viff_gd 2.3: after setting pwkh"
      WRITE(nerr,*) "nle=",nle
!        CALL output2
    ENDIF
  ENDIF  
  !
!
!*   14.2 pgrndcapc, pgrndhflx, current, tsw
!
!!  pseaice(jl)=SEAICEFN(pzsi(jl,1))
!!!  RETURN
  !!! CALL update_snow_ice_property(pzsi(jl,:),hsn,hesn,hice,heice,pseaice(jl))
  CALL update_snow_ice_property(pzsi(jl,0),pzsi(jl,1),hsn,hesn,hice,heice,pseaice(jl))
  !!! hsn=HSNFN(pzsi(jl,0))
  !!! hesn=HEFN(hsn/4.,xksn,omegas)
  !!! hice=HICEFN(pzsi(jl,1))
  !!! heice=HEFN(hice/4.,xkice,omegas)
  xkhskin=SQRT(pwkh(jl,0)*pwkh(jl,nls+1))         ! calc the heat diffusivity for skin layer
  hew=HEFN(hw(0),xkhskin,omegas)

  IF (hesn.GT.csncri) THEN
! snow on top
    pgrndcapc(jl)=rhosn*csn*hesn
  ELSEIF (heice.GT.xicri) THEN
! ice on top
    pgrndcapc(jl)=rhoice*cice*heice
  ELSEIF(.NOT.lsoil) THEN
! water on top
    pgrndcapc(jl)=rhowcw*hew
  ELSE
! soil on top
    pgrndcapc(jl)=rhogcg*SQRT(xkg/omegas)
  ENDIF  
  pgrndhflx(jl)=0.
  pgrndflux(jl)=0.
  pengwac(jl)=0.
  psaltwac(jl)=0.  
  IF (locn) THEN
    sumxxz=0.
    sumxxt=0.
    sumxxu=0.
    sumxxv=0.
    sumxxs=0.
    lstfn=MIN(nls+nfnlvl-1,nle)  ! index of last fine level
    DO jk=nls+1,lstfn
      ! excluidng skin layer
      sumxxz=sumxxz+hw(jk)
      sumxxt=sumxxt+pwt(jl,jk)*hw(jk)
      sumxxu=sumxxu+pwu(jl,jk)*hw(jk)
      sumxxv=sumxxv+pwv(jl,jk)*hw(jk)
      sumxxs=sumxxs+pws(jl,jk)*hw(jk)
      IF ( lwarning_msg.GE.3 ) THEN
        WRITE(nerr,*) 'jk=',jk,'hw=',hw(jk)
      ENDIF
    ENDDO
    pwtb(jl)=sumxxt/sumxxz
    pwub(jl)=sumxxu/sumxxz
    pwvb(jl)=sumxxv/sumxxz
    pwsb(jl)=sumxxs/sumxxz
    psubfluxw(jl)=0.
    pwsubsal(jl)=0.
    IF ( (sumxxz.LE.0.).OR.((pwtb(jl)+pwub(jl)+pwvb(jl)+pwsb(jl)+psubfluxw(jl)+pwsubsal(jl)).LT.-9.E20) ) THEN
      WRITE(nerr,*) ", I am in init_sit_viff_gd 2.6: sumxxz = (<=0)", sumxxz
      WRITE(nerr,*) 'nls+1=',nls+1,'lstfn=',lstfn,'sumxxz=',sumxxz,'sumxxt=',sumxxt,'sumxxu=',sumxxu,'sumxxv=',sumxxv,'sumxxs=',sumxxs
      WRITE(nerr,*) 'pwtb=',pwtb(jl),'pwub=',pwub(jl),'pwvb=',pwvb(jl),'pwsb=',pwsb(jl),'psubfluxw=',psubfluxw(jl),'pwsubsal=',pwsubsal(jl)
      CALL output2      
    ENDIF
  ENDIF
  IF (GDCHK3) then
    WRITE(nerr,*) ", I am in init_sit_viff_gd 2.7: leaving "
    CALL output2
  ENDIF        
  END SUBROUTINE init_sit_viff_gd
! ----------------------------------------------------------------------  
  SUBROUTINE compose_obs_ocean(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
  !
  ! Initialise sit T,U,V,S on ocean levels (lake, ocean and soil)
  ! Change in-situ water temperature to potential water temperature
  !
    !!! USE mo_mpi,           ONLY: mpp_pe() 
    USE mod_eos_ocean,     ONLY: tmelt
    USE mod_sit_control,       ONLY: lgodas,lwoa0
                            
    IMPLICIT NONE
    INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
    LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
    LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
    LOGICAL, INTENT(IN):: l_no_expolation_wu   ! logical for missing data
    LOGICAL, INTENT(IN):: l_no_expolation_wv   ! logical for missing data  
    INTEGER  :: jk,kkk
    real:: depth  
    
  
    IF (GDCHK2) then
      WRITE(nerr,*) ", I am in compose_obs_ocean 1.0: start"
      WRITE(nerr,*) "l_no_expolation_wt=",l_no_expolation_wt
      WRITE(nerr,*) "l_no_expolation_ws=",l_no_expolation_ws
      WRITE(nerr,*) "l_no_expolation_wu=",l_no_expolation_wu
      WRITE(nerr,*) "l_no_expolation_wv=",l_no_expolation_wv
      WRITE(nerr,*) "pobswtb=",pobswtb(jl)
      WRITE(nerr,*) "pobswsb=",pobswsb(jl)
      WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
      !!! WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "obswt,", "obsws", "obswu,", "obswv"
      !!! DO jk = 0, nle+1
      !!!  !!!      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),bg_wt0(jl,jk),bg_ws0(jl,jk),bg_wu0(jl,jk),bg_wv0(jl,jk)
      !!!   WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
      !!! END DO
    ENDIF

    !
    ! 1.0 using observed SST and SSS for the first guess
    !
    IF (.NOT.l_no_expolation_wt) THEN
      DO jk=0,nle+1
        depth=(pwlvl(jl)-z(jk))
        pobswt(jl,jk)=wt_maxden(depth,pobswtb(jl),pobswsb(jl),mixing_depth)
      END DO
    ENDIF
    IF (.NOT.l_no_expolation_ws) pobsws(jl,0:nle+1)=pobswsb(jl)
    IF (.NOT.l_no_expolation_wu) pobswu(jl,0:nle+1)=0.
    IF (.NOT.l_no_expolation_wv) pobswv(jl,0:nle+1)=0.
    IF (GDCHK2) then
      WRITE(nerr,*) ", I am in compose_obs_ocean 2.0: after initial guess"
      WRITE(nerr,*) "pobswtb=",pobswtb(jl)
      WRITE(nerr,*) "pobswsb=",pobswsb(jl)
      WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
      WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "obswt,", "obsws", "obswu,", "obswv"
      DO jk = 0, nle+1
        !!!      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),bg_wt0(jl,jk),bg_ws0(jl,jk),bg_wu0(jl,jk),bg_wv0(jl,jk)
        WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
      END DO
    ENDIF
    !
    ! 2.0 modified again according to monthly OCN profile data such as GODAS
    !    
    IF (lgodas) THEN
      CALL interpolation_godas(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
      ! extrapolation for missing values
      IF (GDCHK2) then
        WRITE(nerr,*) ", I am in compose_obs_ocean 2.0: after monthly OCN profile"
        WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
        WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "obswt,", "obsws", "obswu,", "obswv"
        DO jk = 0, nle+1
          !!!      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),bg_wt0(jl,jk),bg_ws0(jl,jk),bg_wu0(jl,jk),bg_wv0(jl,jk)
          WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
        END DO
      ENDIF
    ENDIF
    !
    ! 3.0 modified again according to initial OCN profile data such as WOA0
    !
    IF (lwoa0) THEN  
      bg_wt0(jl,0:lkvl+1)=xmissing
      bg_ws0(jl,0:lkvl+1)=xmissing
      bg_wu0(jl,0:lkvl+1)=xmissing
      bg_wv0(jl,0:lkvl+1)=xmissing
      ! CALL interpolation_woa0(jl,jrow,.TRUE.,.TRUE.)
      ! using original field in DIECAST for missing values
      ! The above command will crash the model at non-DIECAST grids such as lakes
      CALL interpolation_woa0(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
      pobswt(jl,0:nle+1)=MERGE(bg_wt0(jl,0:nle+1),pobswt(jl,0:nle+1),bg_wt0(jl,0:nle+1).NE.xmissing)
      pobsws(jl,0:nle+1)=MERGE(bg_ws0(jl,0:nle+1),pobsws(jl,0:nle+1),bg_ws0(jl,0:nle+1).NE.xmissing)
      pobswu(jl,0:nle+1)=MERGE(bg_wu0(jl,0:nle+1),pobswu(jl,0:nle+1),bg_wu0(jl,0:nle+1).NE.xmissing)
      pobswv(jl,0:nle+1)=MERGE(bg_wv0(jl,0:nle+1),pobswv(jl,0:nle+1),bg_wv0(jl,0:nle+1).NE.xmissing)
      IF (GDCHK2) then
        WRITE(nerr,*) ", I am in compose_obs_ocean 3.0: after initial OCN profile"
        WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
        WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "obswt,", "obsws", "obswu,", "obswv"
        DO jk = 0, nle+1
          !!!      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),bg_wt0(jl,jk),bg_ws0(jl,jk),bg_wu0(jl,jk),bg_wv0(jl,jk)
          WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
        END DO
      ENDIF
    ENDIF
    !
    ! 4.0 Change from in-situ temperature to potential temperature
    !
    DO jk=0,nle+1
      pobswt(jl,jk)=theta_from_t(pobsws(jl,jk),pobswt(jl,jk)-tmelt,pwlvl(jl)-z(jk),0.)+tmelt
    ENDDO
    !
    IF (GDCHK2) then
      WRITE(nerr,*) ", I am in compose_obs_ocean 4.0: after setting obswt, obsws"
      WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
      WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "obswt,", "obsws", "obswu,", "obswv"
      DO jk = 0, nle+1
        WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
      END DO
    ENDIF
    !
    ! 5.0 Set progonastic variabls pws, pwt, pwu and pwv
    !
    IF (lobs_ocn_rerun) THEN
        IF (GDCHK2) then
          WRITE(nerr,*) ", I am in compose_obs_ocean 5.0: lobs_ocn_rerun=T"
          WRITE(nerr,*) "nle=",nle
        ENDIF
    !! unchanged
    ELSE
    !! coldstart
      !!! 2.0 Set initial pws and pwt according to pobswsb and pobswtb first 
      !!! pws(jl,0:nle+1)=pobswsb(jl)
      pws(jl,0:nle+1)=pobsws(jl,0:nle+1)
      !!! Set some initial values for pwt    
      pwt(jl,0:nle+1)=pobswt(jl,0:nle+1)
      !!! 3.0 Adjust pwt for accounting the ice grid for depth <= 10 m and the limit of pctfreez2
      IF (.NOT.lsoil) pwt(jl,0:nle+1)=MAX(pctfreez2(jl),pwt(jl,0:nle+1))    
      DO jk=0,nle+1
        IF (pzsi(jl,1).GT.0.) THEN
          depth=(pwlvl(jl)-z(jk))
          IF (depth.LE.10.) pwt(jl,jk)=pctfreez2(jl)
        ENDIF
      END DO
      !!! 4.0 Modification according to observations
      IF (lwoa0.OR.lgodas) THEN
        ! reset initial pws and pwt if observation is available  
        pwt(jl,0:nle+1)=MERGE(pobswt(jl,0:nle+1),pwt(jl,0:nle+1),pobswt(jl,0:nle+1).NE.xmissing)
        pws(jl,0:nle+1)=MERGE(pobsws(jl,0:nle+1),pws(jl,0:nle+1),pobsws(jl,0:nle+1).NE.xmissing)
      ENDIF
    ENDIF  
    IF (lobs_ocn_rerun) THEN
      pwu(jl,0:nle+1)=MERGE(pobswu(jl,0:nle+1),pwu(jl,0:nle+1),pobswu(jl,0:nle+1).NE.xmissing)
      pwv(jl,0:nle+1)=MERGE(pobswv(jl,0:nle+1),pwv(jl,0:nle+1),pobswv(jl,0:nle+1).NE.xmissing)
    ELSE
    ! coldstart
      pwu(jl,0:nle+1)=MERGE(pobswu(jl,0:nle+1),0.,pobswu(jl,0:nle+1).NE.xmissing)
      pwv(jl,0:nle+1)=MERGE(pobswv(jl,0:nle+1),0.,pobswv(jl,0:nle+1).NE.xmissing)
      pww(jl,0:nle+1)=0.
    ENDIF
    IF (GDCHK2) then
      WRITE(nerr,*) ", I am in compose_obs_ocean 5.0: leaving "
      CALL output2
    ENDIF
!ps    2300 FORMAT(1X,3A4,     2A9,2A4,1A13,4A9,2A10,10A9)
       2300 FORMAT(1X,A4,A5,A4,2A8,2A4,1A11,4A9,2A10,10A9)
!ps    2301 FORMAT(1X,3(I3,","),2(F8.2,","),2(I3,","),F10.0,9(F10.1,","),2(E10.2,","),&
!ps          1(F8.3,","),2(E9.2,","),10(F8.3,","))  
       2301 FORMAT(1X,(I3,","),(I4,","),(I3,","),2(F7.2,","),2(I3,","),(F10.4,","),9(F8.3,","),2(E10.2,","),&
          1(F8.3,","),2(E9.2,","),10(F8.3,","))  
  END SUBROUTINE compose_obs_ocean
  ! ----------------------------------------------------------------------
  SUBROUTINE interpolation_woa0(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
  !
  !     input
  !      real:: obstsw, ot0, os0, ou0, ov0
  !     output
  !      real:: bg_wt0(0:lkvl+1),bg_ws0(0:lkvl+1)bg_wu0(0:lkvl+1),bg_wv0(0:lkvl+1)
  
    USE mod_sst,       ONLY: nmw1, nmw2, wgt1, wgt2                           
    USE mod_eos_ocean,     ONLY: tmelt
  !!!  USE mo_physc2,        ONLY: ctfreez
    !!! USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, mpp_pe() 
    USE mod_sst,           ONLY: nodepth0, odepth0, ot0, os0, ou0, ov0  
    USE mod_eos_ocean,     ONLY: tmaxden  
    USE mod_sit_control,       ONLY: lwoa0,lsftobswt
    IMPLICIT NONE
    LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
    LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
    LOGICAL, INTENT(IN):: l_no_expolation_wu   ! logical for missing data
    LOGICAL, INTENT(IN):: l_no_expolation_wv   ! logical for missing data
      ! .TRUE. no expolation, missing value returned
      ! .FALSE. expolation enforced, non-missing value returned
    INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
    LOGICAL  :: l_upperdata            ! =TRUE, if data of an upper level is available
    INTEGER  :: jk,kkk
    REAL:: ttt,ttt1,sss,sss1,uuu,uuu1,vvv,vvv1,depth
    INTEGER:: k10m,kmixdepth
    REAL:: t10m,t10m1,pobswt10m,rr,rf_obs,tmixup,mixdepth


    IF (GDCHK2) then
      WRITE(nerr,*) ", I am in interpolation_woa0"
      WRITE(nerr,*) "lwarning_msg=",lwarning_msg
    ENDIF
    IF ((.NOT.lwoa0).AND.((pobswtb(jl).EQ.xmissing).OR.(pobswsb(jl).EQ.xmissing))) THEN
      WRITE(nerr,*) "I am in interpolation_woa0 1.0"
      WRITE(nerr,*) ", FATAL: I am not able to interpolation the grid"
      WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "bg_wt0,", "bg_ws0"
      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(0),pobswtb(jl),pobswsb(jl)
      WRITE(nerr,*) ", pobswtb(jl).or. pobswsb(jl) .EQ.xmissing"
      RETURN
    ENDIF
  !!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
  !!
  !! Note that the vertical index of g3b starts from 1 although the index in the sit_vdiff 
  !!   routine starts from 1,
    IF (GDCHK2) then
      WRITE(nerr,*) "I am in interpolation_woa0 1.0"
      WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
      WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "bg_wt0,", "bg_ws0", "bg_wu0,", "bg_wv0"
      DO jk = 0, nle+1
        WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),bg_wt0(jl,jk),bg_ws0(jl,jk),bg_wu0(jl,jk),bg_wv0(jl,jk)
      END DO
    ENDIF
  
  !!!  IF (GDCHK2) then
  !!!    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "ot0,", "os0"
  !!!  ENDIF   

!find mixing_depth
    IF ( lmixedlayer .AND.              &
       (mixedlayer0(jl,jrow).NE.xmissing)) THEN
       mixedlayer(jl)=max(mixedlayer0(jl,jrow),nudge_depth_10)
       IF(GDCHK3) print *,"woa0:mixedlayer(jl)=",mixedlayer(jl)
    ENDIF
!find mixing_depth abs(t-t10m)<= (0.8) degree  (A.B. Kara et al., 2000)
    k10m=1
    DO WHILE ( (odepths(k10m).LE. 10.) .AND. (k10m.LT.nodepth) )
      k10m=k10m+1
    END DO
    IF (k10m.GT.1) k10m=k10m-1              ! restore back one layer
  
    IF ( (.not.lmixedlayer)   .AND.      &
     (ot0(jl,k10m,jrow).NE.xmissing)    .AND. &
     (ot0(jl,k10m+1,jrow).NE.xmissing)  .AND. &
     (10. .LE.odepths(nodepth)) ) THEN
     t10m=ot0(jl,k10m,jrow)
     t10m1=ot0(jl,k10m+1,jrow)
     pobswt10m=t10m+ &                                    ! pobswt at 10 m
        (10.-odepths(k10m))/(odepths(k10m+1)-odepths(k10m))* (t10m1-t10m)

     kmixdepth=k10m+1
     DO WHILE ( (ot0(jl,kmixdepth,jrow).NE.xmissing) .AND. &
               (kmixdepth.LT.nodepth) )
       tmixdown=ot0(jl,kmixdepth,jrow)
       IF (abs(tmixdown-pobswt10m).LE. 0.8) THEN
         kmixdepth=kmixdepth+1
       ELSE
         exit
       ENDIF
     END DO
     IF(kmixdepth.GT.1) kmixdepth=kmixdepth-1
     tmixup=ot0(jl,kmixdepth,jrow)+ot0(jl,kmixdepth,jrow)
     mixlayer=((pobswt10m-0.8)-tmixup)/(tmixdown-tmixup)*(odepths(kmixdepth+1)-odepths(kmixdepth))&
              +odepths(kmixdepth)
     mixedlayer(jl)=max(mixlayer,nudge_depth_10)
     IF(GDCHK3) print *,"woa0:kmixdepth=",kmixdepth,",mixedlayer=",mixedlayer(jl) &
                       ,",obswt10m=",pobswt10m
    ENDIF


    DO jk=0,nle+1
      depth=(pwlvl(jl)-z(jk))
      IF (lwoa0) THEN
        !! initialized the water profile according to world ocean altas data (woa) data
        kkk=1
        DO WHILE ( (odepth0(kkk).LE.depth) .AND. (kkk.LT.nodepth0) ) 
           kkk=kkk+1
        END DO
        IF (kkk.GT.1) kkk=kkk-1              ! restore back one layer
      ENDIF


  !
  ! 1.0 water salinity
  !
      l_upperdata=.FALSE.
      sss=pobswsb(jl)
      sss1=pobswsb(jl)
      IF (( lwoa0 .AND.                               &
           (os0(jl,kkk,jrow).NE.xmissing)   .AND.     &
           (os0(jl,kkk+1,jrow).NE.xmissing) .AND.     &
           (depth.LE.odepth0(nodepth0)) )) THEN
        ! depth < max. obs. depth
        ! linearly interpolation (no extrapolation)
        l_upperdata=.TRUE.
        sss=os0(jl,kkk,jrow)
        sss1=os0(jl,kkk+1,jrow)
        bg_ws0(jl,jk)=sss+ &
          (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (sss1-sss)
      ELSE
        ! depth > max obs. depth or missing observed data
        ! water salinity
        IF (l_no_expolation_ws) THEN
        ! no extrapolation
          bg_ws0(jl,jk)=xmissing
        ELSE 
        ! assume to be sss1
          bg_ws0(jl,jk)=sss1
  !!!      ! extrapolation
  !!!        bg_ws0(jl,jk)=sss+ &
  !!!           (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (sss1-sss)
        ENDIF
      ENDIF
  !
  ! 1.1 modify pobswsb/freezing temp again according to initial ocean T/S profile.
  !     Note that pobswsb/freezing temp were firstly setup in ioinitial.f90
  !    
      IF (jk.EQ.0) THEN
        pobswsb(jl)=MERGE(bg_ws0(jl,jk),pobswsb(jl),bg_ws0(jl,jk).NE.xmissing)
#if defined (V9897)
          pctfreez2(jl)=tmelt+3.   
#elif defined (TMELTS0)
          pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS1)
          pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1.,tmelt+1.,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS15)
          pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1.5_dp,tmelt+1.5_dp,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS2)
          pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2.,tmelt+2.,pobswsb(jl).NE.xmissing)  ! v9.9003, there is 7600 ice grids, while the observation is 760 ice grid.    
#elif defined (TMELTS25)
          pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2.,tmelt+2.,pobswsb(jl).NE.xmissing)  ! v9.9007.    
#elif defined (TMELTS3)
          pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+3.,tmelt+3.,pobswsb(jl).NE.xmissing)  ! v9.898 (too hot and too salty in S.H.)
#else
          pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)
#endif
      ENDIF
  !
  ! 2.0 water temperature
  !
      l_upperdata=.FALSE.
      ttt=pobswtb(jl)
      ttt1=pobswtb(jl)
      IF (( lwoa0 .AND.                               &
           (ot0(jl,kkk,jrow).NE.xmissing)   .AND.     &
           (ot0(jl,kkk+1,jrow).NE.xmissing) .AND.     &
           (depth.LE.odepth0(nodepth0)) )) THEN
        ! depth < max. obs. depth
        ! linearly interpolation (no extrapolation)
        l_upperdata=.TRUE.
        ttt=ot0(jl,kkk,jrow)
        ttt1=ot0(jl,kkk+1,jrow)
        bg_wt0(jl,jk)=MAX(pctfreez2(jl),ttt+ &
          (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (ttt1-ttt) )
!ps     !shift bg_wt0 to obswtb, and other depth shift same value
        if(jk.eq.0 .AND. pobswtb(jl) .ne. xmissing) then
          psftobswt(jl,0)=pobswtb(jl)-bg_wt0(jl,0)
          if(GDCHK3) then
            WRITE(nerr,*) 'jl=',jl,',bg_wt0(jl,0)=',bg_wt0(jl,0),',pobswtb=',pobswtb(jl) &
                         ,',psftobswt0=',psftobswt(jl,0)
          endif
        endif
        if(psftobswt(jl,0) .ne. 0. )then
          psftobswt(jl,jk)=psftobswt(jl,0)
        endif
        IF(mixedlayer(jl) .GE. nudge_depth_10) then
          IF(depth .LE. mixedlayer(jl) )then
            bg_wt0(jl,jk)=psftobswt(jl,jk)+bg_wt0(jl,jk)
          ELSEIF (depth .GT. mixedlayer(jl) .AND. depth .LT. nudge_depth_100)then            
            rr=(depth-mixedlayer(jl))/(nudge_depth_100-mixedlayer(jl))
            rf_obs=1.-exp(-0.5/rr*exp(0.5/(rr-1.)))
            psftobswt(jl,jk)=rf_obs*psftobswt(jl,0)
            bg_wt0(jl,jk)=psftobswt(jl,jk)+bg_wt0(jl,jk)
          ELSE
            psftobswt(jl,jk)=0.
            bg_wt0(jl,jk)=bg_wt0(jl,jk)
          ENDIF
       
          if(GDCHK3) then
            WRITE(nerr,*) 'jk=',jk,',bg_wt0=',bg_wt0(jl,jk) &
                       ,',psftobswt_jk=',psftobswt(jl,jk)
          endif
        ENDIF
!ps
      ELSE
        ! depth > max obs. depth or missing observed data
        ! water temperature
        IF (l_no_expolation_wt) THEN
        ! no extrapolation
          bg_wt0(jl,jk)=xmissing
        ELSE 
        ! extrapolation
          IF (l_upperdata) THEN
            bg_wt0(jl,jk)=MAX(pctfreez2(jl),ttt+ &
               (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (ttt1-ttt) )
            IF ((tmaxden(pobswsb(jl)).GT.ttt1).AND.(tmaxden(pobswsb(jl)).GT.ttt)) THEN
              bg_wt0(jl,jk)=MIN(tmaxden(bg_ws0(jl,jk)),bg_wt0(jl,jk))
            ELSEIF ((tmaxden(pobswsb(jl)).LT.ttt1).AND.(tmaxden(pobswsb(jl)).LT.ttt)) THEN
              bg_wt0(jl,jk)=MAX(tmaxden(bg_ws0(jl,jk)),bg_wt0(jl,jk))
            ELSE                                      
              bg_wt0(jl,jk)=tmaxden(bg_ws0(jl,jk))
            ENDIF
          ELSE
            ! set initial profile to be expontential decay to tmaxden
            !   =  3.73 C for fresh water s=0
            !   = -4.35 C for ocena water s=36.3 0/00.
            bg_wt0(jl,jk)=wt_maxden(depth,pobswtb(jl),pobswsb(jl),mixing_depth)
          ENDIF
          ! range check
          bg_wt0(jl,jk)=MAX(pctfreez2(jl),bg_wt0(jl,jk))                  
        ENDIF
      ENDIF
  !
  ! 3.0 water u current
  !
      l_upperdata=.FALSE.
      uuu=0.
      uuu1=0.
      IF (( lwoa0 .AND.                               &
           (ou0(jl,kkk,jrow).NE.xmissing)   .AND.     &
           (ou0(jl,kkk+1,jrow).NE.xmissing) .AND.     &
           (depth.LE.odepth0(nodepth0)) )) THEN
        ! depth < max. obs. depth
        ! linearly interpolation (no extrapolation)
        l_upperdata=.TRUE.
        uuu=ou0(jl,kkk,jrow)
        uuu1=ou0(jl,kkk+1,jrow)
        bg_wu0(jl,jk)=uuu+ &
          (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (uuu1-uuu)
      ELSE
        ! depth > max obs. depth or missing observed data
        IF (l_no_expolation_wu) THEN
        ! no extrapolation
          bg_wu0(jl,jk)=xmissing
        ELSE 
        ! extrapolation
          bg_wu0(jl,jk)=0.
  !!!      ! extrapolation
  !!!        bg_wu0(jl,jk)=uuu+ &
  !!!           (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (uuu1-uuu)
        ENDIF
      ENDIF
  !
  ! 4.0 water v current
  !
      l_upperdata=.FALSE.
      vvv=0.
      vvv1=0.
      IF (( lwoa0 .AND.                               &
           (ov0(jl,kkk,jrow).NE.xmissing)   .AND.     &
           (ov0(jl,kkk+1,jrow).NE.xmissing) .AND.     &
           (depth.LE.odepth0(nodepth0)) )) THEN
        ! depth < max. obs. depth
        ! linearly interpolation (no extrapolation)
        l_upperdata=.TRUE.
        vvv=ov0(jl,kkk,jrow)
        vvv1=ov0(jl,kkk+1,jrow)
        bg_wv0(jl,jk)=vvv+ &
          (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (vvv1-vvv)
      ELSE
        ! depth > max obs. depth or missing observed data
        IF (l_no_expolation_wv) THEN
        ! no extrapolation
          bg_wv0(jl,jk)=xmissing
        ELSE 
        ! extrapolation
          bg_wv0(jl,jk)=0.
  !!!      ! extrapolation
  !!!        bg_wv0(jl,jk)=vvv+ &
  !!!           (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (vvv1-vvv)
        ENDIF
      ENDIF
    ENDDO
  
    IF (GDCHK2) then
      WRITE(nerr,*) "I am in interpolation_woa0 2.0"
      WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
      WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "bg_wt0,", "bg_ws0", "bg_wu0,", "bg_wv0"
      DO jk = 0, nle+1
        WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),bg_wt0(jl,jk),bg_ws0(jl,jk),bg_wu0(jl,jk),bg_wv0(jl,jk)
      END DO
    ENDIF
  
  ! Adjust pwt for accounting the ice grid for depth <= 10 m and the limit of pctfreez2
  !!!  IF (.NOT.lsoil) bg_wt0(jl,0:nle+1)=MERGE(MAX(pctfreez2(jl),bg_wt0(jl,0:nle+1)),xmissing,bg_wt0(jl,0:nle+1).NE.xmissing)
    DO jk=0,nle+1
      IF (pobsseaice(jl).GT.0.) THEN
        depth=(pwlvl(jl)-z(jk))
        IF ((depth.LE.10.).AND.(bg_wt0(jl,jk).EQ.xmissing)) bg_wt0(jl,jk)=pctfreez2(jl)
      ENDIF
    END DO
       
    IF (GDCHK3) then
      WRITE(nerr,*) "I am leaving interpolation_woa0"
      WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
      WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "bg_wt0,", "bg_ws0", "bg_wu0,", "bg_wv0,","sftobswt"
      DO jk = 0, nle+1
        WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),bg_wt0(jl,jk),bg_ws0(jl,jk),bg_wu0(jl,jk),bg_wv0(jl,jk),psftobswt(jl,jk)
      END DO
    ENDIF
  
    2300 FORMAT(1X,3A4,2A9,2A4,1A13,4A9,2A10,10A9)
    !ps 2301 FORMAT(1X,3(I3,","),2(F8.2,","),2(I3,","),1(F10.4,","),10(F8.3,","),2(E10.2,","),&
    2301 FORMAT(1X,3(I3,","),2(F8.2,","),2(I3,","),F10.0,9(F10.1,","),2(E10.2,","),&
    &       1(F8.3,","),2(E9.2,","),10(F8.3,","))
  END SUBROUTINE interpolation_woa0
  ! ----------------------------------------------------------------------
  SUBROUTINE interpolation_godas(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
!!!!!!DUPLICATE in sit_vdiff  !!!!!!!!!!!
!
!     input
!      real:: obstsw, ot12, os12, ou12, ov12
!     output
!      real:: obswt(0:lkvl+1),obsws(0:lkvl+1),obswu(0:lkvl+1),obswv(0:lkvl+1)

  USE mod_sst,       ONLY: wgto1, wgto2, now1, now2                 ! GODAS MONTHLY/PENTAD Data                           
  USE mod_sit_control,       ONLY: lgodas,lsftobswt
  USE mod_eos_ocean,     ONLY: tmelt
!!!  USE mo_physc2,        ONLY: ctfreez
  !!! USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, mpp_pe() 
  USE mod_sst,           ONLY: nodepth, odepths, ot12, os12, ou12, ov12
  USE mod_eos_ocean,     ONLY: tmaxden
!ps
  USE mod_sst,        ONLY:obswtbwgt1,obswtbwgt2,obswtbnmw1,obswtbnmw2  
!ps

  IMPLICIT NONE
  LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wu   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wv   ! logical for missing data
    ! .TRUE. no expolation, missing value returned
    ! .FALSE. expolation enforced, non-missing value returned
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  LOGICAL  :: l_upperdata            ! =TRUE, if data of an upper level is available
  INTEGER  :: jk,kkk
  real:: ttt,ttt1,sss,sss1,uuu,uuu1,vvv,vvv1,depth
!ps
!  real:: ttt0,sss0 
  real:: t10m,t10m1,pobswt10m,rr,rf_obs
!ps
  IF (GDCHK2) then
    WRITE(nerr,*) "I am in interpolation_godas1, line 2026"
    WRITE(nerr,*) "lwarning_msg=",lwarning_msg
  ENDIF
  IF ((.NOT.lgodas).AND.(pobswtb(jl).EQ.xmissing)) RETURN
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!!
!! Note that the vertical index of g3b starts from 1 although the index in the sit_vdiff 
!!   routine starts from 1,
  IF (GDCHK2) then
    WRITE(nerr,*) "I am in interpolation_godas1"
    WRITE(nerr,*) "l_no_expolation_wt=",l_no_expolation_wt
    WRITE(nerr,*) "l_no_expolation_ws=",l_no_expolation_ws
    WRITE(nerr,*) "l_no_expolation_wu=",l_no_expolation_wu
    WRITE(nerr,*) "l_no_expolation_wv=",l_no_expolation_wv
    WRITE(nerr,*) "nmw1=",now1,"nmw2=",now2,"wgt1=",wgto1,"wgt2=",wgto2
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "ot1,", "ot2", "os1,", "os2"
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "obswt,", "obsws", "obswu", "obswv"
    DO jk = 0, nle+1
      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
    END DO
  ENDIF


!find mixing_depth
  IF ( lgodas .AND. lmixedlayer .AND.              &
       (mixedlayer12(jl,jrow,now1).NE.xmissing) .AND. &
       (mixedlayer12(jl,jrow,now2).NE.xmissing)) THEN
      mixedlayer(jl)=max(wgto1*mixedlayer12(jl,jrow,now1)+wgto2*mixedlayer12(jl,jrow,now2),nudge_depth_10)
  ENDIF
!find mixing_depth abs(t-t10m)<= (0.8) degree  (A.B. Kara et al., 2000)
  IF (lgodas) then
    k10m=1
    DO WHILE ( (odepths(k10m).LE. 10.) .AND. (k10m.LT.nodepth) )
      k10m=k10m+1
    END DO
    IF (k10m.GT.1) k10m=k10m-1              ! restore back one layer
  ENDIF
  IF ( lgodas .AND. (.not.lmixedlayer)   .AND.      &
     (ot12(jl,k10m,jrow,now1).NE.xmissing)    .AND. &
     (ot12(jl,k10m,jrow,now2).NE.xmissing)    .AND. &
     (ot12(jl,k10m+1,jrow,now1).NE.xmissing)  .AND. &
     (ot12(jl,k10m+1,jrow,now2).NE.xmissing)  .AND. &
     (10. .LE.odepths(nodepth)) ) THEN
    t10m=wgto1*ot12(jl,k10m,jrow,now1)+wgto2*ot12(jl,k10m,jrow,now2)
    t10m1=wgto1*ot12(jl,k10m+1,jrow,now1)+wgto2*ot12(jl,k10m+1,jrow,now2)
    pobswt10m=t10m+ &                                    ! pobswt at 10 m
        (10.-odepths(k10m))/(odepths(k10m+1)-odepths(k10m))* (t10m1-t10m)

    kmixdepth=k10m+1
    DO WHILE ( (ot12(jl,kmixdepth,jrow,now1).NE.xmissing) .AND. &
               (ot12(jl,kmixdepth,jrow,now2).NE.xmissing) .AND. &
               (kmixdepth.LT.nodepth) )
      tmixdown=wgto1*ot12(jl,kmixdepth,jrow,now1)+wgto2*ot12(jl,kmixdepth,jrow,now2)
      IF (abs(tmixdown-pobswt10m).LE. 0.8) THEN
        kmixdepth=kmixdepth+1
      ELSE
        exit
      ENDIF
    END DO
    IF(kmixdepth.GT.1) kmixdepth=kmixdepth-1
    tmixup=wgto1*ot12(jl,kmixdepth,jrow,now1)+wgto2*ot12(jl,kmixdepth,jrow,now2)
    mixlayer=((pobswt10m-0.8)-tmixup)/(tmixdown-tmixup)*(odepths(kmixdepth+1)-odepths(kmixdepth)) &
             +odepths(kmixdepth)
    mixedlayer(jl)=max(mixlayer,nudge_depth_10)
    IF(GDCHK3) print *,"1.kmixdepth=",kmixdepth,",mixedlayer=",mixedlayer(jl) &
                      ,",obswt10m=",pobswt10m
  ENDIF


  DO jk=0,nle+1
    depth=(pwlvl(jl)-z(jk))
    IF (lgodas) THEN
      !! initialized the water profile according to world ocean altas data (woa) data
      kkk=1
      DO WHILE ( (odepths(kkk).LE.depth) .AND. (kkk.LT.nodepth) ) 
         kkk=kkk+1
      END DO
      IF (kkk.GT.1) kkk=kkk-1              ! restore back one layer
    ENDIF

!
! 1.0 water salinity
!
    l_upperdata=.FALSE.
    sss=pobswsb(jl)
    sss1=pobswsb(jl)
    IF ( lgodas .AND.                                    &
         (os12(jl,kkk,jrow,now1).NE.xmissing)      .AND. &
         (os12(jl,kkk,jrow,now2).NE.xmissing)      .AND. &
         (os12(jl,kkk+1,jrow,now1).NE.xmissing)    .AND. &
         (os12(jl,kkk+1,jrow,now2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
      l_upperdata=.TRUE.
!      sss0=obswtbwgt1*os12(jl,1,jrow,obswtbnmw1)+obswtbwgt2*os12(jl,1,jrow,obswtbnmw2)
      sss=wgto1*os12(jl,kkk,jrow,now1)+wgto2*os12(jl,kkk,jrow,now2)
      sss1=wgto1*os12(jl,kkk+1,jrow,now1)+wgto2*os12(jl,kkk+1,jrow,now2)
      pobsws(jl,jk)=sss+ &
        (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (sss1-sss)
    ELSE
      ! depth > max obs. depth or missing observed data
      ! water salinity
      IF (l_no_expolation_ws) THEN
      ! no extrapolation
        pobsws(jl,jk)=xmissing
      ELSE
      ! assume to be sss1
        pobsws(jl,jk)=sss1
!!!   ! extrapolation
!!!      pobsws(jl,jk)=sss+ &
!!!      (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (sss1-sss)
      ENDIF
    ENDIF
!
! 1.1 modify pobswsb/freezing temp again according to initial ocean T/S profile.
!     Note that pobswsb/freezing temp were firstly setup in ioinitial.f90
!    
    IF (lsitstart.AND.jk.EQ.0) THEN
      pobswsb(jl)=MERGE(pobsws(jl,jk),pobswsb(jl),pobsws(jl,jk).NE.xmissing)
#if defined (V9897)
        pctfreez2(jl)=tmelt+3.
#elif defined (TMELTS0)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)           
#elif defined (TMELTS1)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1.,tmelt+1.,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS15)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1.5,tmelt+1.5,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS2)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2.,tmelt+2.,pobswsb(jl).NE.xmissing)  ! v9.9003, there is 7600 ice grids, while the observation is 760 ice grid.    
#elif defined (TMELTS25)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2.,tmelt+2.,pobswsb(jl).NE.xmissing)  ! v9.9007.    
#elif defined (TMELTS3)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+3.,tmelt+3.,pobswsb(jl).NE.xmissing)  ! v9.898 (too hot and too salty in S.H.)
#else
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)
#endif
    ENDIF    
!
! 2.0 water temperature
!
!
    l_upperdata=.FALSE.
    ttt=pobswtb(jl)
    ttt1=pobswtb(jl)    
    IF ( lgodas .AND.                                    &
         (ot12(jl,kkk,jrow,now1).NE.xmissing)      .AND. &
         (ot12(jl,kkk,jrow,now2).NE.xmissing)      .AND. &
         (ot12(jl,kkk+1,jrow,now1).NE.xmissing)    .AND. &
         (ot12(jl,kkk+1,jrow,now2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      l_upperdata=.TRUE.
!      ttt0=obswtbwgt1*ot12(jl,1,jrow,obswtbnmw1)+obswtbwgt2*ot12(jl,1,jrow,obswtbnmw2)
      ttt=wgto1*ot12(jl,kkk,jrow,now1)+wgto2*ot12(jl,kkk,jrow,now2)
      ttt1=wgto1*ot12(jl,kkk+1,jrow,now1)+wgto2*ot12(jl,kkk+1,jrow,now2)

      
      IF(mixedlayer(jl) .ge. nudge_depth_10 .AND. pobswtb(jl) .NE. xmissing   &
         .AND. pobswt10m .NE. xmissing)then

        if(jk .le. 1)then
           pobswt(jl,jk)=MAX(pctfreez2(jl), ttt+ &
            (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
          if(GDCHK3)then
            print*,"int_godas0:jk=",jk,",pobswtb=",pobswtb(jl),",pobswt_jk=",pobswt(jl,jk)
          endif
        endif
        if(lsftobswt) then
          psftobswt(jl,jk)=(pobswtb(jl)-pobswt(jl,0))
        endif

        IF(depth .LE. mixedlayer(jl) ) THEN
          pobswt(jl,jk)=MAX(pctfreez2(jl), psftobswt(jl,jk)+ttt+ &
            (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
        ELSEIF(depth .GT. mixedlayer(jl) .AND. depth .LT. nudge_depth_100) THEN
          rr=(depth-mixedlayer(jl))/(nudge_depth_100-mixedlayer(jl))
          rf_obs=1.-exp(-0.5/rr*exp(0.5/(rr-1.)))
          psftobswt(jl,jk)=rf_obs*psftobswt(jl,jk)
          pobswt(jl,jk)=MAX(pctfreez2(jl), psftobswt(jl,jk)+ttt+ &
            (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
        ELSE
          pobswt(jl,jk)=MAX(pctfreez2(jl), ttt+ &
            (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
        ENDIF
      ELSE
!ps      IF ( depth.LE.10. ) THEN
        IF ( depth .LE. nudge_depth_10 ) THEN
        ! depth < 10 m, note daily SST is avaiable from satellite
          IF (pobswtb(jl).NE.xmissing) THEN
            pobswt(jl,jk)=MAX(pctfreez2(jl),pobswtb(jl))
          ELSE
            pobswt(jl,jk)=xmissing
          ENDIF
        ELSE
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
        pobswt(jl,jk)=MAX(pctfreez2(jl), ttt+ &
          (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
!        pobswt(jl,jk)=MAX(pctfreez2(jl), ttt &
!            +(depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) &
!            +rf_obs*(pobswtb(jl)-ttt0) )
!          IF (GDCHK3) then
!            print *,'now1=',now1,',now2=',now2,'wgto1=',wgto1,',wgto2=',wgto2
!            print *,'kkk=',kkk,',ttt0=',ttt0,',ttt=',ttt,',ttt1=',ttt1
!            print *,'pobsws(',jl,',',jk,')=',pobsws(jl,jk)
!            print *,'pobswt(',jl,',',jk,')=',pobswt(jl,jk)
!            print *,'pobswtb(',jl,')=',pobswtb(jl)
!            print *,'obswtbwgt1=',obswtbwgt1,'obswtbwgt1=',obswtbwgt2
!            print *,'obswtbnmw1=',obswtbnmw1,'obswtbnmw2=',obswtbnmw2
!          ENDIF
         ENDIF
      ENDIF

      IF (GDCHK2) then
        WRITE(nerr,*) "godas1.1, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
      ENDIF
    ELSE
      ! depth > max obs. depth or on observation
      ! water temperature
      IF (l_no_expolation_wt) THEN
      ! no extrapolation
        pobswt(jl,jk)=xmissing
        IF (GDCHK2) then
          WRITE(nerr,*) "godas1.2, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
        ENDIF
      ELSE 
      ! extrapolation
        IF (l_upperdata) THEN
          pobswt(jl,jk)=MAX(pctfreez2(jl),ttt+ &
             (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
          IF (GDCHK2) then
            WRITE(nerr,*) "godas1.3, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
          ENDIF
          IF ((tmaxden(pobswsb(jl)).GT.ttt1).AND.(tmaxden(pobswsb(jl)).GT.ttt)) THEN
            pobswt(jl,jk)=MIN(tmaxden(pobsws(jl,jk)),pobswt(jl,jk))
            IF (GDCHK2) then
              WRITE(nerr,*) "godas1.4, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
            ENDIF
          ELSEIF ((tmaxden(pobswsb(jl)).LT.ttt1).AND.(tmaxden(pobswsb(jl)).LT.ttt)) THEN
            pobswt(jl,jk)=MAX(tmaxden(pobsws(jl,jk)),pobswt(jl,jk))
            IF (GDCHK2) then
              WRITE(nerr,*) "godas1.5, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
            ENDIF
          ELSE                                      
            pobswt(jl,jk)=tmaxden(pobsws(jl,jk))
            IF (GDCHK2) then
              WRITE(nerr,*) "godas1.6, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
            ENDIF
          ENDIF
        ELSE
          ! set initial profile to be expontential decay to tmaxden
          !   =  3.73 C for fresh water s=0
          !   = -4.35 C for ocena water s=36.3 0/00.
          pobswt(jl,jk)=wt_maxden(depth,pobswtb(jl),pobswsb(jl),mixing_depth)                  
          IF (GDCHK2) then
            WRITE(nerr,*) "godas1.7, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
          ENDIF
        ENDIF
        ! range check
        pobswt(jl,jk)=MAX(pctfreez2(jl),pobswt(jl,jk))
        IF (GDCHK2) then
          WRITE(nerr,*) "godas1.8, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
        ENDIF
      ENDIF
    ENDIF
!
! 3.0 water u current
!
    l_upperdata=.FALSE.
    uuu=0.
    uuu1=0.    
    IF ( lgodas .AND.                                    &
         (ou12(jl,kkk,jrow,now1).NE.xmissing)      .AND. &
         (ou12(jl,kkk,jrow,now2).NE.xmissing)      .AND. &
         (ou12(jl,kkk+1,jrow,now1).NE.xmissing)    .AND. &
         (ou12(jl,kkk+1,jrow,now2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      l_upperdata=.TRUE.
      uuu=wgto1*ou12(jl,kkk,jrow,now1)+wgto2*ou12(jl,kkk,jrow,now2)
      uuu1=wgto1*ou12(jl,kkk+1,jrow,now1)+wgto2*ou12(jl,kkk+1,jrow,now2)
      ! linearly interpolation (no extrapolation)
      pobswu(jl,jk)=uuu+ &
        (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (uuu1-uuu)
    ELSE
      IF (l_no_expolation_wu) THEN
        ! no extrapolation
        pobswu(jl,jk)=xmissing 
      ELSE
        pobswu(jl,jk)=0.
!!!     ! extrapolation
!!!       pobswu(jl,jk)=uuu+ &
!!!          (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (uuu1-uuu)
      ENDIF
    ENDIF
!
! 4.0 water v current
!
    l_upperdata=.FALSE.
    vvv=0.
    vvv1=0.    
    IF ( lgodas .AND.                                    &
         (ov12(jl,kkk,jrow,now1).NE.xmissing)      .AND. &
         (ov12(jl,kkk,jrow,now2).NE.xmissing)      .AND. &
         (ov12(jl,kkk+1,jrow,now1).NE.xmissing)    .AND. &
         (ov12(jl,kkk+1,jrow,now2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      l_upperdata=.TRUE.
      vvv=wgto1*ov12(jl,kkk,jrow,now1)+wgto2*ov12(jl,kkk,jrow,now2)
      vvv1=wgto1*ov12(jl,kkk+1,jrow,now1)+wgto2*ov12(jl,kkk+1,jrow,now2)
      ! linearly interpolation (no extrapolation)
      pobswv(jl,jk)=vvv+ &
        (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (vvv1-vvv)
    ELSE
      IF (l_no_expolation_wu) THEN
        ! no extrapolation
        pobswv(jl,jk)=xmissing 
      ELSE
        pobswv(jl,jk)=0.
!!!     ! extrapolation
!!!       pobswv(jl,jk)=vvv+ &
!!!          (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (vvv1-vvv)
      ENDIF
    ENDIF
  ENDDO
!
! 5.0 Adjust pwt for accounting the ice grid for depth <= 10 m and the limit of pctfreez2
!
!  IF (.NOT.lsoil) pobswt(jl,0:nle+1)=MERGE(MAX(pctfreez2(jl),pobswt(jl,0:nle+1)),xmissing,pobswt(jl,0:nle+1).NE.xmissing)
  DO jk=0,nle+1
    IF (pobsseaice(jl).GT.0.) THEN
      depth=(pwlvl(jl)-z(jk))
      IF ((depth.LE.10.).AND.(pobswt(jl,jk).NE.xmissing)) pobswt(jl,jk)=pctfreez2(jl)
    ENDIF
  END DO
  
  IF (GDCHK2) then
    WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "obswt,", "obsws", "obswu", "obswv"
    DO jk = 0, nle+1
      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
    END DO
    WRITE(nerr,*) "I am leaving interpolation_godas"
  ENDIF
 2300 FORMAT(1X,3A4,2A9,2A4,1A11,4A11,2A11,10A9)
 2301 FORMAT(1X,3(I3,","),2(F8.2,","),2(I3,","),(F10.0,","),9(F10.1,","),2(E10.2,","),&
&       1(F8.3,","),2(E9.2,","),10(F8.3,","))
  END SUBROUTINE interpolation_godas  
  !
  !-----------------------------------------------------------------
  !
  SUBROUTINE output2
  !!!! A DUPLICATE VERSION of sit_vdiff_init 
  !
  !-----------------------------------------------------------------
  !
  !*    OUTPUT: WRITE DIAGNOSTIC VARIABLES
  !
  !-----------------------------------------------------------------
    !!! USE mo_mpi,           ONLY: mpp_pe()   
    IMPLICIT NONE
  
    INTEGER jk
    real:: tmp
    
  !
  !!!  CALL pzcord(jl,jrow)   ! bjt 2010/02/21
    WRITE(nerr,*) "kbdim=", kbdim, "kproma=",kproma
    IF (nle.GE.1) THEN
  !     water exists
      WRITE(nerr,*) "Water exist, nle=", nle
    ELSE
  !     soil only
      WRITE(nerr,*) "Soil only, nle=", nle
    ENDIF
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","","","we,","lw,","tsi0,","tsi1,"
    WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,0,0.,pzsi(jl,0),psilw(jl,0),&
  &  ptsnic(jl,0),ptsnic(jl,1)
    WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,1,1.,pzsi(jl,1),psilw(jl,1),&
  &  ptsnic(jl,2),ptsnic(jl,3)
  !  
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,","wt,", "ws,", "wu,", "wv,", &
  &  "wtke,", "wkh,","wldisp,","wlmx,","awtfl,","awtfl0,","awufl,","awvfl,","awsfl,","awsfl0,","awtkefl,","obswt,"
    DO jk = 0, nle+1
      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),&
  &     pwt(jl,jk),pws(jl,jk),pwu(jl,jk),pwv(jl,jk),&
  &     pwtke(jl,jk),pwkh(jl,jk),pwldisp(jl,jk),pwlmx(jl,jk), &
  &     pawtfl(jl,jk),pawtfl0(jl,jk),pawufl(jl,jk),pawvfl(jl,jk),pawsfl(jl,jk),pawsfl0(jl,jk),pawtkefl(jl,jk), &
  &     pobswt(jl,jk)
    END DO
      
  !
    RETURN
  !
   2200 FORMAT(1X,3(A4),2(A9),1(A4),4(A10),2(A9),1(A9), 1(A10))
  !ps 2201 FORMAT(1X,4(I3,","),1(I9,","),2(2(E9.2,","),2(F8.3,",")),&
   2201 FORMAT(1X,3(I3,","),2(F8.2,","),1(I3,","),4(F10.3,","),1(I9,","),2(2(E9.2,","),2(F8.3,",")),&
  &       1(F8.3,","), 1(E10.3,","))
  !
!ps   2300 FORMAT(1X,3A4,      2A9,2A4,1A11,4A9,2A10,10A9)
      2300 FORMAT(1X, A4,A5,A4,2A8,A6,A4,1A11,4A9,2A10,13A9)
!ps   2301 FORMAT(1X,3(I3,","),                 2(F8.2,","),2(I3,","),1(F10.0,","),1(F8.3,","),1(F8.3,","),&
!ps  &       2(F8.2,","),2(E9.2,","),10(F8.3,","))
      2301 FORMAT(1X,(I3,","),(I4,","),(I3,","),2(F7.2,","),(I5,","),(I3,","),1(F10.4,","),1(F8.3,","),1(F8.4,","),&
  &       2(F8.4,","),2(E9.2,","),13(F8.3,","))
  !
  END SUBROUTINE output2
END SUBROUTINE sit_vdiff_init
  
! ****************************************************************************
   SUBROUTINE sit_vdiff ( kproma, kbdim, jrow, istep,  delta_time,    & 
                  plat,       plon, tau, tauhr,                       &
                  pcoriol,    pslm,       plclass,                    &
                  psitmask,   pbathy,     pwlvl,          &
                  pocnmask,   obox_mask,                              &
! - same as lake and ml_ocean
                  pfluxw,     pdfluxs,    psoflw,                     &
                  pfluxi,     psofli,                                 &
! 2-d SIT vars (wind stress)
                  taucx,      taucy,                                  &
! -    water mass variables(rain, snow, evap and runoff):
                  prsf,      pssf,                                    &
                  pevapw,                                             &
                  pdisch,                                             &
! 2-d SIT vars
                  ptemp2,                                             &
                  pwind10w,                                           &
! 2-d SIT vars (sit variables)
                  pobsseaice, pobswtb,    pobswsb,                    &
                  pwtb,       pwub,       pwvb,                       &
                  pwsb,                                               &
                  pfluxiw,    ppme2,                                  &
                  psubfluxw,  pwsubsal,                               &
                  pcc,        phc,        pengwac,                    &
                  psc,        psaltwac,                               &
                  pwtfns,     pwsfns,                                 &
! 3-d SIT vars: snow/ice
                  pzsi,       psilw,      ptsnic,                     &
! 3-d SIT vars: water column
                  pobswt,     pobsws,     pobswu,   pobswv,           &
                  pwt,        pwu,        pwv,      pww,              &
                  pws,        pwtke,      pwlmx,                      & 
                  pwldisp,    pwkm,       pwkh,                       &
                  pwrho1000,                                          &
                  pwtfn,      pwsfn,      pwtfn0,   pwsfn0,           &
                  pawufl,     pawvfl,     pawtfl,                     &
                  pawsfl,     pawtfl0,    pawsfl0,                    &
                  pawtkefl,                                           &
! final output only
                  pseaice,                                            &
                  psni,       psiced,     ptsi,     ptsw,             &    !!
                  ptsl,       ptslm,      ptslm1,                     &
                  pocu,       pocv,       pctfreez2,                  &
! implicit with vdiff
                  pgrndcapc,  pgrndhflx,  pgrndflux,                  &
! store latest 2 steps tsw
                  poldsitwt,  poldsitwu,  poldsitwv,                  &
                  poldsitww,  poldsitws,  poldwtke, pdtswdt,          &
                  psftobswt)
!  ---------------------------------------------------------------------
!
!                           SUBROUTINE DESCRIPTION
!
! 1. FUNCTION DESCRIPTION
!
!     THIS SUBROUTINE CALCULATES WATER BODY THERMOCLINE STRUCTURE
!     AS WELL AS OVER-WATER BODY SNOW/ICE STRUCTURE. THERMOCLINE
!     STRUCTURE IS DETERMINED BY THE TKE METHOD
!     BY COMPUTING pwtke, pwlmx AND pwldisp, AND UPDATE
!     U, T, s DUE TO TURBULENCE MIXING EFFECTS. The design of shallow
!     water mode is to preserve salinity conservation at the conditions
!     water completely evaporated or frozen. The water does not have heat
!     capacity. In addition, it set taucx,taucy to 0 if snow/ice on top
!
! 2. CALLING MODULES
!
!          *sit_vdiff* is called from *physc*.
!
! 3. PARAMETER SPECIFICATION
!     zdtime: time step of the simulation (s), time step of ECHAM     I
!         is used.
!
!      Arguments.
!      ----------
! local dimensions
!  kproma   : IF(jrow==ldc%ngpblks, ldc%npromz,ldc% nproma),
!             number of local longitudes
!  kbdim    : ldc% nproma, gauss grid description
!  jrow     : sequential index for output message
!  istep    : # of time steps
!  pcoriol : coriol factor = 2*omega*sin(latitudes)  (1/s)                         I
!  plat    : latitide (geg)                                                        I
!  plon    : longitude (deg)                                                       I
!  tau     : current forecast time (in hours)                                      I
!  tauhr   : current forecast hr (00-24z)                                          I
! - water body or ml_ocean variables
!  pfluxw    : net surface energy flux over open water per water fraction          I
!              (icesheet+openwater)*(open water fraction) (w/m2) (positive upward)
!              = -(net solar + net longwave - sensible heat - latent heat)*(1-seaice)
!    i.e., fluxw(jl)=-(1.-seaice(jl,jrow))*(zhflw(jl)+zhfsw(jl)+ztrflw(jl)+zsoflw(jl))
!  pdfluxw   : dG/dT (W/m2/K) over water (positive upward)                                    I
!  pfluxi    : net surface energy flux over icesheet per water fraction            I
!              (icesheet+openwater)*(seaice fraction) (w/m2) (positive upward)
!              = -(net solar + net longwave - sensible heat - latent heat)*(seaice)
!    i.e., fluxi(jl)=-seaice(jl,jrow)*(zhfli(jl)+zhfsi(jl)+ztrfli(jl)+zsofli(jl))
!  pdfluxi   : dG/dT (W/m2/K) over ice (positive upward)                                    I
!  psoflw    : net SW flux over open water per water fraction (w/m2) (positive downward)          I
!    i.e., soflw(jl)=(1.-seaice(jl,jrow))*zsoflw(jl)
!  psofli    : net SW flux over over icesheet per water fraction (w/m2) (positive downward)       I
!    i.e., sofli(jl)=seaice(jl,jrow)*zsofli(jl)
!  taucx    : u-stress (Pa) over ice/water at current time step                        I/O
!              (set to tauwx at zero if snow/ice on top)
!  taucy    : v-stress (Pa) over ice/water at current time step                        I
!              (set to tauwy at zero if snow/ice on top)

!  prsf    : large scale + convective rain flux at the surface (kg/m2/s) (positive downward)   I
!  pssf    : large scale + convective snow flux at the surface (kg/m2/s) (positive downward)   I
!  pevapw   : evaporation from water surface (kg/m2/s) (positive downward)         I
!!!!  pruntoc   : surface runoff into ocean (kg/m**2s)                                I
!  pdisch   : surface runoff into ocean (m/s)                                I
!  ptemp2    : 2 m temperature (K) at current time step                            I
!  pwind10w     : 10m windspeed over water (m/s)                                   I
! 2-d SIT vars
!  pobsseaice  : observed sea ice fraction (fraction)                              I
!  pobswtb     : observed bulk sea surface temperature (K)                         I
!  pobswsb    : observed salinity (PSU, 0/00)                                      I
!  psitmask : mask for sit(1=.TRUE., 0=.FALSE.)                                    I
!  pbathy  : bathymeter (topography or orography) of ocean (m)                     I
!  pctfreez2 : ref water freezing temperature (K)                                  I
!  pwlvl  : current water level (ice/water interface) a water body grid            I/O
!  pocnmask : fractional mask for 3-D ocean grid (0-1)                             I
!  obox_mask : 3-D ocean nudging mask, =0: nudging, = 1 (>0): nudging              I
!  pwtb: bulk water temperature (K)                                                O
!  pwub: bulk water u current (m/s)                                                O
!  pwvb: bulk water v current (m/s)                                                O
!  pwsb: bulk water salinity (PSU)                                                 O
!  pfluxiw: over-water net surface heat flux (W/m2, + upward, i.e., from ocean))   O
!  ppme2: net fresh water into ocean (P-E+ice_corr) (m/s, + downward)            O
!  psubfluxw: subsurface ocean heat flux (W/m2, + upward)                          O
!  pwsubsal: subsurface ocean salinity flux (m*PSU/s, + upward)                    O
!  pcc: cold content per water fraction (ice sheet+openwater) (J/m2)
!    (energy need to melt snow and ice, i.e., energy below liquid water at tmelt)  I/O
!  phc: heat content per water fraction (ice sheet+openwater) (J/m2)
!    (energy of a water column above tmelt)                                        I/O
!  pengwac: accumulated energy per water fraction (ice sheet+openwater) (+ downward, J/m2)
!    (pfluxw+pfluxi+rain/snow advected energy in respect to liquid water at
!     tmelt)*dt                                                                    I/O
!!!!  pengw: mean net surface heat flux per water fraction (ice sheet+openwater)
!!!!    (+ downward, W/m2)                                                            I/O
!!!!  pengw2: mean snow corrected net surface heat flux per water fraction
!!!!    (ice sheet+openwater) (+ downward, W/m2)
!!!!    (pfluxw+pfluxi+rain/snow advected heat flux in respect to liquid water at
!!!!     tmelt)                                                                       I/O
!  psc: salinity content per water fraction (ice sheet+openwater) (PSU*m)          I/O
!  psaltwac: accumulated salt into water fraction (+ downward, PSU*m)              I/O
!
! 3-d SIT vars: snow/ice
!  pzsi   :
!        pzsi(jl,0): dry snow water equivalent (m)                                 I/O
!        pzsi(jl,1): dry ice water equivalent (m)                                  I/O
!  psilw  :
!        psilw(jl,0): snow liquid water storage (m)                                I/O
!        psilw(jl,1): ice liquid water storage (m)                                 I/O
!  ptsnic   :
!        ptsnic(jl,0): snow skin temperatrue (jk)                                  I/O
!        ptsnic(jl,1): mean snow temperatrue (jk)                                  I/O
!        ptsnic(jl,2): ice skin temperatrue (jk)                                   I/O
!        ptsnic(jl,3): mean ice temperatrue (jk)                                   I/O
! 3-d SIT vars: water column
!  pobswt: observed potentail water temperature (K)                                I/O
!  pobsws: observed salinity (PSU, 0/00)                                           I/O
!  pobswu: observed u-current (m/s)                                                I/O
!  pobswv: observed v-current (m/s)                                                I/O
!  pwt      : potential water temperature (K)                                      I/O
!  pwu      : u current (m/s)                                                      I/O
!  pwv      : v current (m/s)                                                      I/O
!  pww      : w current (m/s)                                                      I/O
!  pws      : practical salinity (0/00)                                            I/O
!  pwtke    : turbulent kinetic energy (M2/S2)                                     I/O
!  pwlmx    : mixing length (m)                                                    O
!  pwldisp  : dissipation length (m)                                               O
!  pwkm     : eddy diffusivity for momentum (m2/s)                                 O
!  pwkh     : eddy diffusivity for heat (m2/s)                                     O
!  pwrho1000: potential temperature at 1000 m (PSU)                                O
!  pwtfn  : nudging flux into sit-ocean at each level (W/m**2)                     I/O
!  pwsfn  : nudging salinity flux into sit-ocean at each level (PSU*m/s)           I/O
!  pwtfn0  : nudging flux into sit-ocean at each level (W/m**2)                     I/O
!  pwsfn0  : nudging salinity flux into sit-ocean at each level (PSU*m/s)           I/O
!  pwtfns : nudging flux into sit-ocean for an entire column (W/m**2), i.e.
!     pwtfns=SUM(pwtfn(jl,:))                                                      I/O
!  pwsfn  : nudging salinity flux into sit-ocean for an entire column (PSU*m/s),
!     i.e., pwsfns=SUM(pwsfn(jl,:)                                                 I/O
!  pawufl : advected u flux at each level (positive into ocean) (m/s/s)            I
!  pawvfl : advected v flux at each level (positive into ocean) (m/s/s)            I
!  pawtfl : advected temperature flux at each level (positive into ocean) (K/s)(12 mons)    I
!  pawsfl : advected salinity flux at each level (positive into ocean) (PSU/s)(12 mons)     I
!  pawtfl0 : advected temperature flux at each level (positive into ocean) (K/s)(1 recored)    I
!  pawsfl0 : advected salinity flux at each level (positive into ocean) (PSU/s)(1 record)     I
!  pawtkefl : advected tke at each level (positive into ocean) (m3/s3)             I
! - variables internal to physics
!  pslm: land fraction [0-1]                                                       I
!  pseaice   : ice cover (fraction of 1-SLM) (0-1)                                 I/O
!  plclass    :
!  psni      : snow thickness over ice (m in water equivalent)                     I/O
!  psiced    : ice thickness (m in water equivalent)                               I/O
!  ptsw      : skin temperatrue (K) over water                                     O
!  ptsl      : calcuated Earth's skin temperature from vdiff/tsurf at t+dt (K)     I/O
!  ptslm     : calcuated Earth's skin temperature from vdiff/tsurf at t (K)        I/O
!  ptslm1    : calcuated Earth's skin temperature from vdiff/tsurf at t-dt (K)     I/O
! 2-d SIT vars
!  pocu      : ocean eastw. velocity (m/s)                                         O
!  pocv      : ocean northw. velocity (m/s)                                        O
!
!  pgrndcapc: areal heat capacity of the uppermost sit layer (snow/ice/water)
!    (J/m**2/K)                                                                    I/O
!  pgrndhflx: ground heat flux below the surface (W/m**2)
!    (+ upward, into the skin layer)                                               I/O
!  pgrndflx:  acc. ground heat flux below the surface (W/m**2*s)
!    (+ upward, into the skin layer)                                               I/O
!  poldtsw    : store latest 1 steps tsw                                           I
!  poldtsitwt : store latest 2 steps sitwt                                         I/O
!  poldtsitwu : store latest 2 steps sitwu                                         I/O
!  poldtsitwv : store latest 2 steps sitwv                                         I/O
!  poldtsitww : store latest 2 steps sitww                                         I/O
!  poldtsitws : store latest 2 steps sitws                                         I/O
!  poldtwtke  : store latest 2 steps wtke                                          I/O
!  pdtswdt   : tsw tendency                                                       O
!  psftobswt : the shift swt in mixed layer based on difference of obswtb and obswt10m I/O
!
! 4. LOCAL VARIABLE  (STAGGERED GRID)
!
!     tsim: ptsnic temperatures at pervious time step, respectively (jk)
!     wum: VELOCITY (kk) AT PAST T LVL (M/s)
!     wtm: TEMPERATURE (kk) AT PAST T LVL (jk)
!     wsm: salinity (kk) at past T LVL (KG/KG)
!     wtkem: TURBULENCE KINETIC ENERGY (kk) AT PAST LVL (M2/S2)
!     nle: VERTICAL DIMENSION
!     AA: INPUT TRI-DIAGONAL MATRIX COEF.        (MN,kk,4)
!     RHS: RIGHT HAND SIDE OF THE MATRIX (FORCING)
!     SOL: SOLUTION OF THE MATRIX
!     z: ELEVATION AT THE CENTER OF A GIRD                  (MN,kk+1)
!     z(0): ELEVATION OF THE SURFACE OF THE WATER BODY, ZERO IS DEFAULT (M)
!     z(nle+1): ELEVATION OF THE BOTTOM OF THE WATER BODY (M)
!     beta : PARAMETER TO CONTROL TIME SCHEME ;
!         1.0 -> BACKWARD, 1/2. -> CRANK-NICHOLSON, 0. -> FORWARD.
!     beta2: PARAMETER TO CONTROL TIME SCHEME OF TKE;
!         1.0 -> BACKWARD, 1/2. -> CRANK-NICHOLSON, 0. -> FORWARD.
!
!
!
! 6. USAGE
!
!      CALL thermocline(LKID,jl,JL,istep,zdtime,SRFL,pfluxw,pdfluxw,
!     &                prsf,pssf,
!     &                TSN,SN,istep,lsit_debug,
!     &                pwkm,pwlmx,pwldisp)
!
! 7. FUNCTIONS CALLED
!
!     FFN
!     rhofn
!     LU,LU2
!     pzcord(LKID,jl,nls,nle,z,zlk,hw)
!     LKDIFC(mas,nle,zdtime,z,pwkm,hw,X,Y)
!     lkerr
!
! 8. LIMITATION
!
!     1) Sometimes, the model is unable to judge the ice is completely
!        melted or still some ice remains (recursive more than 4 times).
!        It is because the solar radiation might penetrate into certain
!        depth. That causes the surface energy budget of ice on top
!        slightly differ from water on top. While ice on top, all the
!        solar radiation will be used to melt ice. If water on top, there is
!        some solar radiation penetrate into deeper layer causing the skin
!        colder
!        than tmelt, and water starts to refreeze.
!     2) If the depths of snow, ice & water below their critial depths
!        (csncri,xicri, wcri), their temperature are assumed to be their
!        underneath temperatures in order to prevent "divide by zero error".
!        However, this introduces artifical energies. More rigous tests are
!        needed to decide their existence.
!
!


  IMPLICIT NONE
  !
  ! Arguments
  !
  ! local dimensions
  INTEGER, INTENT(in)                      :: kproma ! number of local longitudes
  INTEGER, INTENT(in)                      :: kbdim  ! number of local longitudes

  ! gauss grid description
  INTEGER, INTENT(in)                      :: jrow   ! sequential index
  INTEGER, INTENT(in)                      :: istep  ! # of time steps
  real, INTENT(in):: delta_time                       ! time step of the parent routine (s)

  real, INTENT(in):: plat(kbdim), plon(kbdim)   
  real, INTENT(in):: tau,  tauhr             ! current forecast time (tau:in hours, tauhr: 00-24z)
  real, INTENT(in):: pcoriol(kbdim)
  real, INTENT(in):: pslm(kbdim), plclass(kbdim)
  real, INTENT(in out):: psitmask(kbdim)         ! grid mask for lsit (1 or 0)      
  real, INTENT(in out):: pbathy(kbdim)
  real, INTENT(in out):: pwlvl(kbdim)
  real, INTENT(in):: pocnmask(kbdim)         ! (fractional) grid mask for 3-D ocean (DIECAST)
  real, INTENT(in):: obox_mask(kbdim)        ! (fractional) ocean iop nudging mask for 3-D ocean (DIECAST)

! - same as lake and ml_ocean
  real, INTENT(in):: pfluxw(kbdim), pdfluxs(kbdim), psoflw(kbdim)
  real, INTENT(in):: pfluxi(kbdim), psofli(kbdim)

! 2-d SIT vars (wind stress)
  real, INTENT(in out):: taucx(kbdim),    taucy(kbdim)

! -    water mass variables(rain, snow, evap and runoff):
  real, INTENT(in):: prsf(kbdim), pssf(kbdim)
  real, INTENT(in):: pevapw(kbdim)
  real, INTENT(in):: pdisch(kbdim)

! 2-d SIT vars
  real, INTENT(in):: ptemp2(kbdim), pwind10w(kbdim)

! 2-d SIT vars (sit variables)
  real, INTENT(in):: pobsseaice(kbdim)
  real, INTENT(in out):: pobswtb(kbdim),      pobswsb(kbdim)
  real, INTENT(out):: pwtb(kbdim), pwub(kbdim), pwvb(kbdim), pwsb(kbdim)
  real, INTENT(out):: pfluxiw(kbdim),ppme2(kbdim)
  real, INTENT(out):: psubfluxw(kbdim), pwsubsal(kbdim)

  real, INTENT(in out) :: pcc(kbdim), phc(kbdim), pengwac(kbdim)
  real, INTENT(in out) :: psc(kbdim), psaltwac(kbdim)
  real, INTENT(in out) :: pwtfns(kbdim),  pwsfns(kbdim)
! 3-d SIT vars: snow/ice
  real, INTENT(in out)::                                                &
       pzsi(kbdim,0:1),   psilw(kbdim,0:1), ptsnic(kbdim,0:3)
! 3-d SIT vars: water column
  real, INTENT(in out)::                                                &
       pobswt(kbdim,0:lkvl+1), pobsws(kbdim,0:lkvl+1),                       &
       pobswu(kbdim,0:lkvl+1), pobswv(kbdim,0:lkvl+1),                       &
       pwt(kbdim,0:lkvl+1), pwu(kbdim,0:lkvl+1), pwv(kbdim,0:lkvl+1),        &
       pww(kbdim,0:lkvl+1),                                                  &
       pws(kbdim,0:lkvl+1), pwtke(kbdim,0:lkvl+1), pwlmx(kbdim,0:lkvl+1),    & 
       pwldisp(kbdim,0:lkvl+1), pwkm(kbdim,0:lkvl+1), pwkh(kbdim,0:lkvl+1),  &
       pwrho1000(kbdim,0:lkvl+1)
  real, INTENT(in out):: pwtfn(kbdim,0:lkvl+1),  pwsfn(kbdim,0:lkvl+1)
  real, INTENT(in out):: pwtfn0(kbdim,0:lkvl+1), pwsfn0(kbdim,0:lkvl+1)
  real, INTENT(in out):: pawufl(kbdim,0:lkvl+1), pawvfl(kbdim,0:lkvl+1), pawtfl(kbdim,0:lkvl+1)
  real, INTENT(in out):: pawtfl0(kbdim,0:lkvl+1),pawsfl0(kbdim,0:lkvl+1)
  real, INTENT(in out):: pawsfl(kbdim,0:lkvl+1), pawtkefl(kbdim,0:lkvl+1)


! final output only
  real, INTENT(in out):: pseaice(kbdim)
  real, INTENT(in out):: psni(kbdim), psiced(kbdim), ptsi(kbdim),ptsw(kbdim)
  real, INTENT(in out):: ptsl(kbdim), ptslm(kbdim), ptslm1(kbdim)
  real, INTENT(in out):: pocu(kbdim), pocv(kbdim), pctfreez2(kbdim)

! implicit with vdiff
  real, INTENT(in out) :: pgrndcapc(kbdim), pgrndhflx(kbdim), pgrndflux(kbdim)

! store latest 2 step tsw
  real, INTENT(in out) :: poldsitwt(kbdim,0:lkvl+1,0:1)
  real, INTENT(in out) :: poldsitwu(kbdim,0:lkvl+1,0:1)
  real, INTENT(in out) :: poldsitwv(kbdim,0:lkvl+1,0:1)
  real, INTENT(in out) :: poldsitww(kbdim,0:lkvl+1,0:1)
  real, INTENT(in out) :: poldsitws(kbdim,0:lkvl+1,0:1)
  real, INTENT(in out) :: poldwtke(kbdim,0:lkvl+1,0:1)
  real, INTENT(in out) :: pdtswdt(kbdim)
  real, INTENT(in out) :: psftobswt(kbdim,0:lkvl+1)

! - local variables


  !*    1.0 Depth of the coordinates of the water body point. (zdepth = 0 at surface )
!!!  real, PARAMETER:: zdepth(0:lkvl)=(/0.,0.0005_dp,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,&
!!!              20.,30.,50.,75.,100.,125.,150.,200.,250.,300.,400.,500.,&
!!!              600.,700.,800.,900.,1000.,1100.,1200.,1300.,1400.,1500./) 
  !
  !********************
  ! Note that WOA 2005 data are at depths:
  !  depth = 0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600,
  !    700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500 ;
  !********************

  !*    1.1  INFORMATION OF OCEAN GRID (COMLKE)
  !     WATER BODY POINT VERTICAL LEVELS
  INTEGER nls ! first water level -1,
              ! starting water level of the water body point.
              ! z(k)=sitwlvl, IF k is <= nls. 
  INTEGER nle ! nle: last water level (excluding skin level)
              ! nle: maximum water levels of the water body point.
              ! Note nle=-1 IF water level is below point water body bed.
              ! Soil layer always is sitwt(nle+1). If there is no water,
              ! sitwt(nle+1)=sitwt(0,jrow). In addition, z(k)=bathy IF k >= nle.
  !!! INTEGER nlv ! number of water levels (excluding water surface &
  !!!             ! soil layer). Note nlv=-1 IF water level is below
  !!!             ! water bed.  
  real, DIMENSION(0:lkvl+1):: z   ! coordinates of the water body point.
  real, DIMENSION(0:lkvl+1):: zlk ! standard z coordinates of a water body 
  real, DIMENSION(0:lkvl):: hw    ! thickness of each layer
  LOGICAL lshlw  ! .true. = shallow water mode (due to evaportion or freezing)
              ! (one layer water body)
  LOGICAL lsoil  ! .true. = soil grid
!
!*      1.1 CONSTANTS (LOCAL PARAMETERS)
!    1) eddy DIFFUSION PARAMETERS (GASPAR ET AL., 1990)
!ps  real,PARAMETER::ck=0.1    ! (0.1 in GASPAR ET AL., 1990) (v9.887)
!!!#if defined (Prw001)
!!!  real,PARAMETER::ck=0.1    ! (0.1 in GASPAR ET AL., 1990) (v9.887)
!!!  real,PARAMETER::Prw=0.01_dp ! 1., Prt_molecular=8.96, (40.,50.)  (v8.6)  0.0006 m2/s ocnkvm 0.016
!!!#elif defined (Prw01)
!!!  real,PARAMETER::ck=0.1    ! (0.1 in GASPAR ET AL., 1990) (v9.887)
!!!  real,PARAMETER::Prw=0.1_dp ! 1., Prt_molecular=8.96, (40.,50.)  (v8.6)  0.0006 m2/s ocnkvm 0.016
!!!#elif defined (Prw1)
!!!  real,PARAMETER::ck=0.1    ! (0.1 in GASPAR ET AL., 1990) (v9.887)
!!!  real,PARAMETER::Prw=1. ! 1., Prt_molecular=8.96, (40.,50.)  (v8.6)  0.0006 m2/s ocnkvm 0.016
!!!#else
!!!  real,PARAMETER::ck=1.0        ! (v9.885 - v9.886 for sensitivity test)  (OK)
!!!  real,PARAMETER::Prw=1. ! 1., Prt_molecular=8.96, (40.,50.) ! (Cold tongue is vanished) ! v9.4
!!!#endif
!ps  real,PARAMETER::ce=0.7    ! (0.7 in Bougeault and Lacarrere, 1989)
!  real,PARAMETER::Prw=10. ! 1., Prt_molecular=8.96, (40.,50.) ! V9.2 v9.3
!  real,PARAMETER::Prw=1. ! 1., Prt_molecular=8.96, (40.,50.) ! (Cold tongue is vanished) ! v9.4
!  real,PARAMETER::Prw=0.04_dp ! 1., Prt_molecular=8.96, (40.,50.)  (v8.6)  0.0006 m2/s ocnkvm 0.016
!  real,PARAMETER::Prw=0.01_dp ! 1., Prt_molecular=8.96, (40.,50.)  (v8.6)  0.0006 m2/s ocnkvm 0.016
!!!  real,PARAMETER::emin=0. !limit min pwtke to 0. to avoid too warm below mixing depth (20090710) (v5.4)(v9.913)
  real,PARAMETER::emin=1.0E-6 !limit min pwtke to 1.0E-6 (v9.865, v9.885) (1.0E-6 in GASPAR ET AL., 1990)
!  real,PARAMETER::emin=1.0E-4 !limit min pwtke to 1.0E-4 (v9.866),(v9.868)
!  real,PARAMETER::emin=1.0E-5 !limit min pwtke to 1.0E-5 (v9.867) (v9.873)
  real,PARAMETER::xkmmin=1.2E-6     ! molecular momentum diffusivity (Paulson and Simpson, 1981; Chia and pwu, 1998; Mellor and Durbin, 1975)
  real,PARAMETER::xkhmin=1.34E-7    ! molecular heat diffusivity (Paulson and Simpson, 1981; Chia and pwu, 1998; Mellor and Durbin, 1975)
!  real,PARAMETER::hcoolskin=4.E-4 ! thickness of conductive sublayer (m), where only molecular transfer exists (m) (Khundzuha et al., 1977; Paulson and Simpson, 1981)
!ps  real,PARAMETER::d0=0.03  ! zero-displacement (m) (0.02,0.05) (v0.9871)
!  real,PARAMETER::d0=0.    ! The wt(nle+1) at 30-50 m is too high, 
                                   ! due to the very small vertical diffusivity but still strong solar radiaion  (v0.9871)
!  real,PARAMETER::xlkmin=0.3        ! minimum pwlmx (0.5,0.8) produced unreasoable high SST near Equator bjt 2007/8, v09.83
  real,PARAMETER::xlkmin=0.        ! minimum pwlmx (0.5,0.8) produced unreasoable high SST near Equator bjt 2013/9, v0.984
  real,PARAMETER::xldispmin=0.3 ! minimum pwldisp
!  real,PARAMETER::xldispmin=0. ! 1., minimum pwldisp, (0.5,0.6)
! optimal (Lotus: Case 2: WT10, emin=3.E-5, xldispmin=2., xlkmin=0., stderr=0.282), no warm-layer
! optimal (Lotus: Case 3: WT10, emin=4.E-5, xldispmin=0.4, xlkmin=0., stderr=0.397)
! optimal (Lotus: Case 4: WT0, emin=3.E-5, xldispmin=0.6, xlkmin=0., stderr=0.260)
! optimal (Lotus: Case 5: WT0, emin=2.E-5, xldispmin=0.03, xlkmin=0.05, stderr=0.529)
! optimal (Lotus: Case 6: WT0, emin=3.E-5, xldispmin=2., xlkmin=0.,d0=0.001, stderr=0.307), no warm-layer
! optimal (Lotus: Case 7: WT0, emin=1.E-5, xldispmin=0.03, xlkmin=0.,d0=0.001, stderr=1.45), has warm-layer
! optimal (TOGA: Case 1: WT0, emin=1.E-6, xldispmin=0., xlkmin=0.)
! optimal (Lotus: WT0, emin=3.E-5, Prw=40., stderr=0.368)
!       v.78
!v90  real,PARAMETER::EResRatio=0.01 ! emin=EResRatio*E(0), 1% of TKE(0) (0.5%,1%) 
!v90  real,PARAMETER::EMPOW=.8               !limit min pwtke to 1.0E-6 (0.5,1)
!v90  real,PARAMETER::EMINMAX=3.E-5 ! 1., minimum pwldisp, (1.E-5,3.E-5)
!v90  real,PARAMETER::EMINMIN=3.E-5 ! 1., minimum pwldisp, (0.5,0.6)
! optimal (Lotus: Case 8: WT0, EResRatio=3%,EMINMAX=1.E-4,EMINMIN=1.E-7,EMPOW=1, xldispmin=1., xlkmin=0.,d0=0.02, stderr=0.54), has warm-layer
! optimal (Lotus: Case 9: WT0, EResRatio=3%,EMINMAX=1.E-4,EMINMIN=2.E-5,EMPOW=1., xldispmin=1., xlkmin=0.,d0=0.02,stderr=0.461), no warm-layer
! optimal (Lotus: Case 10: WT0, EResRatio=2%,EMINMAX=1.E-4,EMINMIN=1.E-6,EMPOW=0.8, xldispmin=1., xlkmin=0.,d0=0.02,stderr=0.537), no warm-layer

! optimal (Lotus: Case 1: WT0, EResRatio=1%, d0=0.03, stderr=1.90)
! optimal (TOGA: Case 1: WT0, EResRatio=1%, d0<0.05)
!  REAL, PARAMETER:: zepcor=5.e-05           ! minimum corilol force, 5.e-05 at 20 deg. There are still huge current within 20S-20N
!  REAL, PARAMETER:: zepcor=2*omegas*sin(30./180.*api) ! minimum corilol force  30 deg: 7.29212E-05
!  REAL, PARAMETER:: zepcor=7.3e-05           ! minimum corilol force, 7.3 e-05 at 30 deg
  real, PARAMETER:: zepcor=0.                 ! minimum corilol force, 7.3 e-05 at 30 deg
  real:: fr=7.E-5                             ! = friction factor for current (1/s). It is due to that one-column model neglect horizontal diffusion.
                                              ! It should be 0. for 3-D ocean. Otherwise oscilliation and current > 3 m/s may occurs.
  
!pwldisp
! Security number
  real, PARAMETER:: SALT_MAX=41. ! maximun salinity
#if defined(SALTI0)
  real, PARAMETER:: salti=0. ! salinity of ice (PSU)  
!!!#else  
!!!  real, PARAMETER:: salti=20. ! salinity of ice (PSU)  
#endif
!
!*    4) PARAMETER TO CONTROL TIME SCHEME; 
!
  real, PARAMETER:: beta=1.
  real, PARAMETER:: beta2=0.5
  real, PARAMETER:: zero_hour=0.
  real, PARAMETER:: one_hour=3600.
  real, PARAMETER:: six_hour=6.*3600.
  real, PARAMETER:: one_day=86400.  
  real, PARAMETER:: one_month=30.*86400.  
!
!     beta: 1.-> BACKWARD,  1/2.-> CRANK-NICOLSON, 0.-> FORWARD.
!
!*    5) Logics to control the run:
!
  LOGICAL, PARAMETER::  lv81=.TRUE.
  LOGICAL, PARAMETER::  lwaterlevel=.FALSE.
! .TRUE. for handeling water level change
! .FALSE. for not fixed water level
  LOGICAL:: lpenetrative_convection=.FALSE.  ! .true. = set penetrative convection
  INTEGER, PARAMETER:: debug_level=1
!!#ifdef ARGCHECK
!!  LOGICAL, PARAMETER:: lssst=.False.
!!#else
!!  LOGICAL, PARAMETER:: lssst=.TRUE.
!!#endif
! .TRUE. for with thermocline skin layer
! .FALSE. for without thermocline skin layer
!!  LOGICAL, PARAMETER:: lsit_ice=.FALSE.
!!  LOGICAL, PARAMETER:: lsit_salt=.FALSE.
  !! .FALSE.=turn off the salinity module
  LOGICAL, PARAMETER:: ldiag=.FALSE.
! LOGICAL, PARAMETER:: lwarning_msg=.TRUE.
!!!  INTEGER, PARAMETER:: sit_ice_option=1
  LOGICAL, PARAMETER:: ldeep_water_nudg=.FALSE.
!! deep water column nudging where obsevation data are not available
  real, PARAMETER:: QFLTI=0.2
  real, PARAMETER:: QFLTA=0.2
    ! QFLTI+QFLTA should be a value within 0-1., a security number
    ! for preventing ocillation. It should be modified with
    ! implicit coupling with the atmospehre for getting 1st order accuracy.
    ! QFLTI: weighting of current time step
    ! QFLTA: weighting of previous time step
    ! (1-QFLTI-QFLTA): weighting of previous time step   
!
!     6) LOCAL VARIABLES, ARRAY
!
  INTEGER:: lstfn     ! index of last fine level
  real, DIMENSION(0:3)::      tsim
  real, DIMENSION(0:lkvl+1):: wtm,wum,wvm,wsm
  real, DIMENSION(0:lkvl+1):: wtkem
  ! potential water density at the surface for T and S at old time step
  real, DIMENSION(0:lkvl+1):: rhom                           
  ! potential water density at one level higher (denoted as "h")
  real, DIMENSION(0:lkvl+1):: pwrhoh,rhomh,pwrho
  real:: pgrndhflx_int(kbdim)
  real:: pfluxiw_int(kbdim),ppme2_int(kbdim),pwsubflux_int(kbdim), pwsubsal_int(kbdim)
  real:: zsf     ! salinity flux (PSU*m) (positive upward)  
  real:: wlvlm  ! old water level (m in elevation)
  real:: heice  ! effective skin thickness of ice (m)
  real:: hice   ! ice thickness (m)
  real:: hesn   ! effective skin thickness of snow (m)
  real:: hsn    ! snow thickness (m)
  real:: hew    ! effective skin thickness of water (m)
  real:: xkhskin   ! skin layer heat diffusivity (m2/s)
  real:: pfluxwm  ! original surface energy flux over open water per water fraction (icesheet+openwater) (W/m2)
  real:: pfluxim  ! original surface energy flux over icesheet per water fraction (icesheet+openwater)(W/m2)
  real:: pfluxw2  ! same as pfluxw but (adv energy included) (w/m2) (positive upward) 
  real:: pfluxi2  ! same as pfluxi but (adv energy included) (w/m2) (positive upward)
  real:: ccm      ! cold content at previous time step (J/m2) 
  real:: utauw         ! water-side friction velocity (m/s)     
!  real:: tauwx,    tauwy   ! wind stress over water, They will 0 while iced.
!!!!***************
!!!! Local variables for backgroud initial ocean profiles
!!!  real:: bg_wt0(kbdim,0:lkvl+1)
!!!  real:: bg_ws0(kbdim,0:lkvl+1)
!!!  real:: bg_wu0(kbdim,0:lkvl+1)
!!!  real:: bg_wv0(kbdim,0:lkvl+1)
!***************
  INTEGER :: jl
! bjt
  INTEGER, DIMENSION(1) :: imax, imin
  LOGICAL :: lsitmask(kbdim)   ! sit mask
!!!  INTEGER :: istep                          ! istep=time step
  INTEGER :: i_sit_step, n_sit_step         ! istep=time step
! bjt
  real:: zdtime                         ! sit time step (s)
  real:: acc_time                       ! accum. time (s)
  real:: tmaxb,tminb,tmax,tmin  

!ps  
  real, DIMENSION(kbdim):: mixedlayer  
  real, DIMENSION(kbdim):: poldtsw    !store latest 1 steps tsw  
  real:: st_restore_jl,ut_restore_jl,dt_restore_jl 
  real:: ss_restore_jl,us_restore_jl,ds_restore_jl 
  real:: suv_restore_jl,uuv_restore_jl,duv_restore_jl 
  integer:: isit_nudg
  real:: mixlayer, tmixup, tmixdown
! ----------------------------------------------------------------------

! store latest tsw (pskuo)
      poldtsw=ptsw
      if(GDCHK3) then
        print *,"oldtsw=",poldtsw &
               ,",pwt=",pwt 
      endif
      poldsitwt(:,:,0)=MERGE(pwt(:,:),poldsitwt(:,:,1),poldsitwt(:,:,1).EQ.xmissing)
      poldsitwu(:,:,0)=MERGE(pwu(:,:),poldsitwu(:,:,1),poldsitwu(:,:,1).EQ.xmissing)
      poldsitwv(:,:,0)=MERGE(pwv(:,:),poldsitwv(:,:,1),poldsitwv(:,:,1).EQ.xmissing)
      poldsitww(:,:,0)=MERGE(pww(:,:),poldsitww(:,:,1),poldsitww(:,:,1).EQ.xmissing)
      poldsitws(:,:,0)=MERGE(pws(:,:),poldsitws(:,:,1),poldsitws(:,:,1).EQ.xmissing)
      poldwtke(:,:,0)=MERGE(pwtke(:,:),poldwtke(:,:,1),poldwtke(:,:,1).EQ.xmissing)

      poldsitwt(:,:,1)=pwt(:,:)
      poldsitwu(:,:,1)=pwu(:,:)
      poldsitwv(:,:,1)=pwv(:,:)
      poldsitww(:,:,1)=pwv(:,:)
      poldsitws(:,:,1)=pws(:,:)
      poldwtke(:,:,1)=pwtke(:,:)

! -----------------------------------------------
!!!  istep=999
!  istep=get_time_step()
  lsitmask=psitmask.EQ.1.                          ! convert from REAL to Integer and to Logical
!
! 1.0 Determine time step
!     
  n_sit_step = CEILING(delta_time/900.)
  zdtime=delta_time/DBLE( n_sit_step)

!
! 2.0 Initialization
!     
  IF (GDCHK1) then
     WRITE(nerr,*) ", I am in sit_vdiff" 
  ENDIF

  IF (GDCHK1) THEN
    WRITE(nerr,*) "lsitstart=.FALSE."
    WRITE(nerr,*) "sit_vdiff: jrow=",jrow   
    WRITE(nerr,*) "delta_time=",delta_time,"zdtime=",zdtime,"n_sit_step=",n_sit_step
    WRITE(nerr,*) "lsitmask=",lsitmask    

     tmax=maxval(ptsw(1:kproma),mask=lsitmask(1:kproma))
     tmin=minval(ptsw(1:kproma),mask=lsitmask(1:kproma))
     imax=maxloc(ptsw(1:kproma),mask=lsitmask(1:kproma))
     imin=minloc(ptsw(1:kproma),mask=lsitmask(1:kproma))
!!!     IF ((tmax.GE.tmelt+60.).or.(tmin.LE.tmelt-100.)) THEN
        !!! WRITE(nerr,*) "istep=",istep,"jrow=,",jrow,"istep=",istep
        WRITE(nerr,*) "tmax=",tmax,"at",imax,                          &
           "lat=",plat(imax),"lon=",plon(imax),    &
           "slm=",pslm(imax),"lclass=",plclass(imax),    &
           "sit=",lsitmask(imax)
        WRITE(nerr,*) "tmin=",tmin,"at",imin,                          &
           "lat=",plat(imin),"lon=",plon(imin),    &
           "slm=",pslm(imin),"lclass=",plclass(imin),    &
           "sit=",lsitmask(imin)
        tmaxb=maxval(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        tminb=minval(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        imax=maxloc(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        imin=minloc(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        WRITE(nerr,*) "tmaxb=",tmaxb,"at",imax,                          &
           "lat=",plat(imax),"lon=",plon(imax),    &
           "slm=",pslm(imax),"lclass=",plclass(imax),    &
           "sit=",lsitmask(imax)
        WRITE(nerr,*) "tminb=",tminb,"at",imin,                          &
           "lat=",plat(imin),"lon=",plon(imin),    &
           "slm=",pslm(imin),"lclass=",plclass(imin),    &
           "sit=",lsitmask(imin)
        WRITE(nerr,*) "sitmask=",lsitmask(1:kproma)
        WRITE(nerr,*) "lclass=",plclass(1:kproma)
        WRITE(nerr,*) "slm=",pslm(1:kproma)
        WRITE(nerr,*) "seaice=",pseaice(1:kproma)
        WRITE(nerr,*) "tsw=",ptsw(1:kproma)
        WRITE(nerr,*) "obstsw=",pobswtb(1:kproma)
        WRITE(nerr,*) "obswsb=",pobswsb(1:kproma)
        WRITE(nerr,*) "fluxi=",pfluxi(1:kproma)
        WRITE(nerr,*) "dfluxs=",pdfluxs(1:kproma)
        WRITE(nerr,*) "sofli=",psofli(1:kproma)
        WRITE(nerr,*) "fluxw=",pfluxw(1:kproma)
        WRITE(nerr,*) "soflw=",psoflw(1:kproma)
        WRITE(nerr,*) "wtfn(10)=",pwtfn(1:kproma,10)
        WRITE(nerr,*) "wsfn(10)=",pwsfn(1:kproma,10)
        WRITE(nerr,*) "wtfn0(10)=",pwtfn0(1:kproma,10)
        WRITE(nerr,*) "wsfn0(10)=",pwsfn0(1:kproma,10)
        WRITE(nerr,*) "disch=",pdisch(1:kproma)
        WRITE(nerr,*) "temp2=",ptemp2(1:kproma)
        WRITE(nerr,*) "wind10w=",pwind10w(1:kproma)
        WRITE(nerr,*) "wlvl=",pwlvl(1:kproma)
        WRITE(nerr,*) "evapw=",pevapw(1:kproma)
        WRITE(nerr,*) "rsf=",prsf(1:kproma)
        WRITE(nerr,*) "ssf=",pssf(1:kproma)
!!!        CALL mp_barrier()
!!!        call mpp_error ('sit_vdiff', 'unreasonable ptsw', FATAL)
! 
!!!     ENDIF
  ENDIF


!
! 3.0 Start to work
!     

  DO jl=1,kproma
  IF (locaf0) CALL interpolation_ocaf0(jl,jrow,.TRUE.,.TRUE.)
    IF (lsitmask(jl)) THEN     ! sitmask true
      IF(GDCHK3) THEN
 3300 FORMAT(1X,27(A11,E13.5))

        WRITE(nerr,3300)           &
           "1sitvdifftau,",tau,",ptsw,",ptsw(jl)                       &
          ,",pfluxw,",pfluxw(jl),",pdfluxs,",pdfluxs(jl)     &
          ,",psoflw,",psoflw(jl),",pfluxi,",pfluxi(jl),",psofli,",psofli(jl) &
          ,",taucx,",taucx(jl),",taucy,",taucy(jl),",prsf,",prsf(jl)   &
          ,",pssf,",pssf(jl),",pevapw,",pevapw(jl),"ptemp2",ptemp2(jl) &
          ,",pwind10w,",pwind10w(jl),",pobswtb,",pobswtb(jl)           &    
          ,",pobswsb,",pobswsb(jl),",ptsw,",ptsw(jl),",pdtswdt,",pdtswdt(jl)&
          ,",pwt0,",pwt(jl,0),",pws0,",pws(jl,0)  &
          ,",plat,",plat(jl),",plon,",plon(jl)  &
          ,",sbox_w(1),",sitbox_nudg_w(1),",sbox_e(1),",sitbox_nudg_e(1) &
          ,",sbox_s(1),",sitbox_nudg_s(1),",sbox_n(1),",sitbox_nudg_n(1)

        WRITE(nerr,*) "nsit_nudg=", nsit_nudg

      ENDIF
      CALL pzcord(pwlvl(jl),pbathy(jl),nls,nle,z,zlk,hw,lsoil,lshlw)  ! bjt 2010/2/21
      IF (GDCHK2) THEN
        WRITE(nerr,*) ", I am in sit_vdiff 3.1: entering"
        print*,"ps: 3.0 psilw(jl,:)=",psilw(jl,:)
        CALL output2
      ENDIF
      pgrndhflx(jl)=0.
      pfluxiw(jl)=0.
      ppme2(jl)=0.
      IF (locn) THEN
        psubfluxw(jl)=0.
        pwsubsal(jl)=0.
      ENDIF
      IF (locaf) CALL interpolation_ocaf(jl,jrow,.TRUE.,.TRUE.)
!      IF (locaf0) CALL interpolation_ocaf0(jl,jrow,.TRUE.,.TRUE.)
      IF (GDCHK2) THEN
        WRITE(nerr,*) ", I am in sit_vdiff 3.2: after interpolation_ocaf"
        WRITE(nerr,*) "pawtfl=",pawtfl(jl,:)
        WRITE(nerr,*) "pawtfl0=",pawtfl0(jl,:)
        CALL output2
      ENDIF
      acc_time=0.
      DO i_sit_step=1, n_sit_step
!        IF (GDCHK3) THEN
!          WRITE(nerr,*) "i_sit_step=",i_sit_step,",n_sit_step=",n_sit_step
!        ENDIF
        acc_time=acc_time+zdtime
        CALL thermocline(jl,jrow)
        CALL acc_flux
        IF (GDCHK2) THEN
          WRITE(nerr,*) ", I am in sit_vdiff 3.3: after thermocline"
          CALL output2
        ENDIF
        ! nudging the sit_vdiff grid accoridng to water temperature
!sitbox
        st_restore_jl =ssit_restore_time
        ut_restore_jl =usit_restore_time
        dt_restore_jl =dsit_restore_time   
        ss_restore_jl =ssits_restore_time
        us_restore_jl =usits_restore_time
        ds_restore_jl =dsits_restore_time 
        suv_restore_jl=ssituv_restore_time
        uuv_restore_jl=usituv_restore_time
        duv_restore_jl=dsituv_restore_time
        IF (nsit_nudg .GE. 1) THEN
          DO isit_nudg= 1, nsit_nudg
            IF( (plon(jl) .GE. sitbox_nudg_w(isit_nudg)) .AND. &
                (plon(jl) .LT. sitbox_nudg_e(isit_nudg)) .AND. &
                (plat(jl) .GE. sitbox_nudg_s(isit_nudg)) .AND. &
                (plat(jl) .LT. sitbox_nudg_n(isit_nudg)) ) THEN
              st_restore_jl =sitbox_st_restore_time(isit_nudg)
              ut_restore_jl =sitbox_ut_restore_time(isit_nudg)
              dt_restore_jl =sitbox_dt_restore_time(isit_nudg)
              ss_restore_jl =sitbox_ss_restore_time(isit_nudg)
              us_restore_jl =sitbox_us_restore_time(isit_nudg)
              ds_restore_jl =sitbox_ds_restore_time(isit_nudg)
              suv_restore_jl=sitbox_suv_restore_time(isit_nudg)
              uuv_restore_jl=sitbox_uuv_restore_time(isit_nudg)
              duv_restore_jl=sitbox_duv_restore_time(isit_nudg)
              IF (GDCHK3) THEN
                WRITE(nerr,*) "in sitbox:",  &
                  ",st_restore_jl=",st_restore_jl,",ut_jl=",ut_restore_jl, &
                  ",dt_jl=",dt_restore_jl,",ss_jl=",ss_restore_jl, &
                  ",us_jl=",us_restore_jl,",ds_jl=",ds_restore_jl
              ENDIF
            ENDIF 
          END DO
        ENDIF
        CALL nudging_sit_viff_gd(jl,jrow,obox_restore_time,       &
          socn_restore_time,uocn_restore_time,docn_restore_time,   &
          ssit_restore_time,ssits_restore_time,ssituv_restore_time, &
          st_restore_jl,ut_restore_jl,dt_restore_jl,   &
          ss_restore_jl,us_restore_jl,ds_restore_jl,  &
          suv_restore_jl,uuv_restore_jl,duv_restore_jl,tau, tauhr)
        IF (GDCHK3) THEN
          WRITE(nerr,*) ", I am in sit_vdiff 3.4: after nudging_sit_viff_gd"
          WRITE(nerr,*) "lsice_nudg=",lsice_nudg,",delta_time=",delta_time
          WRITE(nerr,*) "lsice_nudg=",lsice_nudg,",n_sit_step=",n_sit_step
          WRITE(nerr,*) "zdtime=",zdtime 
          CALL output2
        ENDIF
      END DO
      CALL final
    ENDIF
!
! 4.0 Warning message
!         
!    IF (GDCHK2) then
!       WRITE(nerr,*) ", I am leaving sit_vdiff"
!!!       CALL output2
!    ENDIF
  END DO

  IF (GDCHK1) THEN
     WRITE(nerr,*) ", I am leaving sit_vdiff II:"
     tmax=maxval(ptsw(1:kproma),mask=lsitmask(1:kproma))
     tmin=minval(ptsw(1:kproma),mask=lsitmask(1:kproma))
     imax=maxloc(ptsw(1:kproma),mask=lsitmask(1:kproma))
     imin=minloc(ptsw(1:kproma),mask=lsitmask(1:kproma))
     IF ((tmax.GE.tmelt+60.).or.(tmin.LE.tmelt-100.)) THEN
        !!! WRITE(nerr,*) "istep=",istep,"jrow=,",jrow,"istep=",istep
        WRITE(nerr,*) "tmax=",tmax,"at",imax,                          &
           "lat=",plat(imax),"lon=",plon(imax),    &
           "slm=",pslm(imax),"lclass=",plclass(imax),    &
           "sit=",lsitmask(imax)
        WRITE(nerr,*) "tmin=",tmin,"at",imin,                          &
           "lat=",plat(imin),"lon=",plon(imin),    &
           "slm=",pslm(imin),"lclass=",plclass(imin),    &
           "sit=",lsitmask(imin)
        tmaxb=maxval(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        tminb=minval(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        imax=maxloc(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        imin=minloc(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        WRITE(nerr,*) "tmaxb=",tmaxb,"at",imax,                          &
           "lat=",plat(imax),"lon=",plon(imax),    &
           "slm=",pslm(imax),"lclass=",plclass(imax),    &
           "sit=",lsitmask(imax)
        WRITE(nerr,*) "tminb=",tminb,"at",imin,                          &
           "lat=",plat(imin),"lon=",plon(imin),    &
           "slm=",pslm(imin),"lclass=",plclass(imin),    &
           "sit=",lsitmask(imin)
        WRITE(nerr,*) "sitmask=",lsitmask(1:kproma)
        WRITE(nerr,*) "lclass=",plclass(1:kproma)
        WRITE(nerr,*) "slm=",pslm(1:kproma)
        WRITE(nerr,*) "seaice=",pseaice(1:kproma)
        WRITE(nerr,*) "tsw=",ptsw(1:kproma)
        WRITE(nerr,*) "obstsw=",pobswtb(1:kproma)
        WRITE(nerr,*) "obswsb=",pobswsb(1:kproma)
        WRITE(nerr,*) "fluxi=",pfluxi(1:kproma)
        WRITE(nerr,*) "sofli=",psofli(1:kproma)
        WRITE(nerr,*) "fluxw=",pfluxw(1:kproma)
        WRITE(nerr,*) "dfluxs=",pdfluxs(1:kproma)
        WRITE(nerr,*) "soflw=",psoflw(1:kproma)
        WRITE(nerr,*) "wtfn(10)=",pwtfn(1:kproma,10)
        WRITE(nerr,*) "wsfn(10)=",pwsfn(1:kproma,10)
        WRITE(nerr,*) "wtfn0(10)=",pwtfn0(1:kproma,10)
        WRITE(nerr,*) "wsfn0(10)=",pwsfn0(1:kproma,10)
        WRITE(nerr,*) "disch=",pdisch(1:kproma)
        WRITE(nerr,*) "temp2=",ptemp2(1:kproma)
        WRITE(nerr,*) "wind10w=",pwind10w(1:kproma)
        WRITE(nerr,*) "wlvl=",pwlvl(1:kproma)
        WRITE(nerr,*) "evapw=",pevapw(1:kproma)
        WRITE(nerr,*) "rsf=",prsf(1:kproma)
        WRITE(nerr,*) "ssf=",pssf(1:kproma)
!!!        CALL mp_barrier()
#if defined (LCWBGFS)
        WRITE(nerr,*), "sit_vdiff, unreasonable ptsw, myrank=", myrank 
        call mpe_finalize
#else
        call mpp_error ('sit_vdiff', 'unreasonable ptsw', FATAL)
#endif        
! 
     ENDIF
  ENDIF
  RETURN
! **********************************************************************
CONTAINS
  !----------------------------------------------------------  
  !*    5.0  SUBROUTINES
  ! **********************************************************************
  SUBROUTINE final
  !
    !!! USE mo_semi_impl,        ONLY: eps
  !!!   USE convect_tables_mod,   ONLY: jptlucu1, jptlucu2
  !   INTEGER, PARAMETER:: jptlucu1 =  50000  ! lookup table lower bound (50K)
  !   INTEGER, PARAMETER:: jptlucu2 = 400000  ! lookup table upper bound (400K) 
    IMPLICIT NONE  
    real:: sumxxz,sumxxt,sumxxu,sumxxv,sumxxs
    real, PARAMETER:: tmin_table=jptlucu1/1000.
    real, PARAMETER:: tmax_table=jptlucu2/1000.
  !!!  real, PARAMETER:: eps=0.001_dp
    INTEGER:: jk
! location blending
    real:: rrn,rrs,rre,rrw,rr_sn,rr_we,rrtmp,wtr
!    real, PARAMETER:: rf=10._dp
    real, PARAMETER:: rf=5.
!ps    

#if defined(LCWBGFS)
    Do jk=nls,nle+1
       pwt(jl,jk)=0.2*pwt(jl,jk)+0.3*poldsitwt(jl,jk,1)+0.5*poldsitwt(jl,jk,0)
       pwu(jl,jk)=0.2*pwu(jl,jk)+0.3*poldsitwu(jl,jk,1)+0.5*poldsitwu(jl,jk,0)
       pwv(jl,jk)=0.2*pwv(jl,jk)+0.3*poldsitwv(jl,jk,1)+0.5*poldsitwv(jl,jk,0)
       pww(jl,jk)=0.2*pww(jl,jk)+0.3*poldsitww(jl,jk,1)+0.5*poldsitww(jl,jk,0)
       pws(jl,jk)=0.2*pws(jl,jk)+0.3*poldsitws(jl,jk,1)+0.5*poldsitws(jl,jk,0)
       pwtke(jl,jk)=0.2*pwtke(jl,jk)+0.3*poldwtke(jl,jk,1)+0.5*poldwtke(jl,jk,0)
!       pwt(jl,jk)=1.*pwt(jl,jk)+0.*poldsitwt(jl,jk,1)+0.*poldsitwt(jl,jk,0)
!       pwu(jl,jk)=1.*pwu(jl,jk)+0.*poldsitwu(jl,jk,1)+0.*poldsitwu(jl,jk,0)
!       pwv(jl,jk)=1.*pwv(jl,jk)+0.*poldsitwv(jl,jk,1)+0.*poldsitwv(jl,jk,0)
!       pww(jl,jk)=1.*pww(jl,jk)+0.*poldsitww(jl,jk,1)+0.*poldsitww(jl,jk,0)
!       pws(jl,jk)=1.*pws(jl,jk)+0.*poldsitws(jl,jk,1)+0.*poldsitws(jl,jk,0)
!       pwtke(jl,jk)=1.*pwtke(jl,jk)+0.*poldwtke(jl,jk,1)+0.*poldwtke(jl,jk,0)
    ENDDO
      if(GDCHK3) print *,"final,jk=0,oldwt(0)=",poldsitwt(jl,0,0),",oldwt(1)=",poldsitwt(jl,0,1)
      if(GDCHK3) print *,"final,jk=1,oldwt(0)=",poldsitwt(jl,1,0),",oldwt(1)=",poldsitwt(jl,1,1)
#endif


    
  IF (ltrigsit) THEN
  !!! couple with SIT: return SIT sst and ice to atmosphere model
  !   
  !*     5.1 pgrndcapc, pgrndhflx, current, tsw
  !   
      IF (hesn.GT.csncri) THEN
  !   snow on top
        pgrndcapc(jl)=rhosn*csn*hesn
        ptsl(jl)=MAX(MIN(ptsl(jl),tmelt),tmin_table)
        ptslm(jl)=MAX(MIN(ptslm(jl),tmelt),tmin_table)
        ptslm1(jl)=MAX(MIN(ptslm1(jl),tmelt),tmin_table)
      ELSEIF (heice.GT.xicri) THEN
  !   ice on top
        pgrndcapc(jl)=rhoice*cice*heice
        ptsl(jl)=MAX(MIN(ptsl(jl),tmelt),tmin_table)
        ptslm(jl)=MAX(MIN(ptslm(jl),tmelt),tmin_table)
        ptslm1(jl)=MAX(MIN(ptslm1(jl),tmelt),tmin_table)
      ELSEIF(.NOT.lsoil) THEN
  !   water on top
        pgrndcapc(jl)=rhowcw*hew
        ptsl(jl)=MIN(MAX(ptsl(jl),pctfreez2(jl)),tmelt+100.)
        ptslm(jl)=MIN(MAX(ptslm(jl),pctfreez2(jl)),tmelt+100.)
        ptslm1(jl)=MIN(MAX(ptslm1(jl),pctfreez2(jl)),tmelt+100.)
      ELSE
  !   soil on top
        pgrndcapc(jl)=rhogcg*SQRT(xkg/omegas)
      ENDIF
      pocu(jl)=pwu(jl,0)
      pocv(jl)=pwv(jl,0)
  !*       5.3     Time filter for surface temperature
      IF (.NOT.lsitstart) THEN
        ptslm1(jl)=ptslm(jl)+eps*(ptslm1(jl)-2.*ptslm(jl)+ptsl(jl))
        ptslm(jl)=ptsl(jl)        
      ELSE
        ptslm1(jl)=ptslm(jl)
      ENDIF
  !   
  !*    5.2 snow/ice properties
  !   
      IF(lsit_ice) THEN
        psni(jl)=pzsi(jl,0)
        !! an extra varible for snow is needed for partial water/partial land,
        !! and modification is need for subroutine albedo for distinquish snow on ice
        !! or snow on land
        psiced(jl)=pzsi(jl,1)
  !!  !
  !!  
  !!  ! change to lognormal distribution:
  !!  !
        pseaice(jl)=SEAICEFN(psiced(jl))    
        !
        IF (hesn.GT.csncri) THEN
  !     snow on top
          IF(sit_ice_option.EQ.0) THEN
            ptsnic(jl,0)=ptslm1(jl)
          ELSE IF (sit_ice_option.EQ.1) THEN
            ptsnic(jl,0)=MAX(ptsnic(jl,0),tmelt-10.)
            ! tmelt-10.: security number
            ! for preventing ocillation
            ! It should be modified with implicit coupling
            ! with the atmospehre.
          ELSE IF (sit_ice_option.EQ.2) THEN
            ptsnic(jl,0)=ptsnic(jl,0)
            ! This also crash after few time steps.
          ELSE IF (sit_ice_option.EQ.3) THEN
            ptsnic(jl,0)=ptsnic(jl,1)
            ! This also crash after few time steps.
          ELSE IF (sit_ice_option.EQ.4) THEN
            ptsnic(jl,0)=MAX(QFLTI*ptsnic(jl,0)+QFLTA*ptemp2(jl)+(1.-QFLTI-QFLTA)*tsim(0),tmelt-50.)
            ! This also crash after few time steps.
            ! QFLTI should be a value within 0-1., a security number
            ! for preventing ocillation. It should be modified with
            ! implicit coupling with the atmospehre for getting 1st order accuracy.
          ENDIF
          ptsi(jl)=ptsnic(jl,0)
          ptsw(jl)=pwt(jl,0)
        ELSEIF (heice.GT.xicri) THEN
  !     ice on top
          IF(sit_ice_option.EQ.0) THEN
            ptsnic(jl,2)=ptslm1(jl)
          ELSE IF (sit_ice_option.EQ.1) THEN
            ptsnic(jl,2)=MAX(ptsnic(jl,2),tmelt-10.)
             ! tmelt-10.: security number
             ! for preventing ocillation
             ! It should be modified with implicit coupling
             ! with the atmospehre.
          ELSE IF (sit_ice_option.EQ.2) THEN
            ptsnic(jl,2)=ptsnic(jl,2)
              ! This will creash in few time steps
          ELSE IF (sit_ice_option.EQ.3) THEN
            ptsnic(jl,2)=ptsnic(jl,3)
              ! This will creash in few time steps
          ELSE IF (sit_ice_option.EQ.4) THEN
            ptsnic(jl,2)=MAX(QFLTI*ptsnic(jl,2)+QFLTA*ptemp2(jl)+(1.-QFLTI-QFLTA)*tsim(2),tmelt-50.)
              ! QFLTI should be a value within 0-1., a security number
              ! for preventing ocillation. It should be modified with
              ! implicit coupling with the atmospehre for getting 1st order accuracy.
          ENDIF
          ptsi(jl)=ptsnic(jl,2)
          ptsw(jl)=pwt(jl,0)
        ELSEIF(.NOT.lsoil) THEN
  !     water on top
          IF(sit_ice_option.EQ.0) THEN
            !!! ptsw(jl)=ptslm1(jl)
            !!! rather using pwt than ptslm1 for ptsw
            ptsw(jl)=pwt(jl,0)
          ELSE
            ptsw(jl)=pwt(jl,0)
          ENDIF
          ptsi(jl)=tmelt
        ELSE
  !     soil on top
          ptsi(jl)=tmelt
          ptsw(jl)=pwt(jl,nle+1)
        ENDIF
      ELSE
        IF(.NOT.lsoil) THEN
          ptsw(jl)=pwt(jl,0) 
        ELSE
  !     soil on top
          ptsw(jl)=pwt(jl,nle+1)
        ENDIF
      ENDIF
  !!!????    pobswsb(jl)=pws(jl,0)
  !   
  !!!   IF ( (ptsw(jl).GE.400_dp).OR.(ptsw(jl).LT.50.) ) THEN
  !!!     CALL output(hesn,hew,heice,fcew,pfluxwm,wtm,wum,wvm,wsm,wtkem,tsim,mas,mae)
  !!!   ENDIF
!ps
#if defined(LCWBGFS)
!!! location blending
    IF(ptsw(jl) .lt. pctfreez2(jl) .or. ptsw(jl) .gt. 350.) then
      ptsw(jl)=pobswtb(jl)
!      pdtswdt(jl)=(ptsw(jl)-poldtsw(jl))/(n_sit_step*zdtime)
!      pdtswdt(jl)=(pwt(jl,0)-poldsitwt(jl,0,0))/(2.*n_sit_step*zdtime)
      pdtswdt(jl)=(pwt(jl,0)-poldsitwt(jl,0,1))/(n_sit_step*zdtime)
      psitmask(jl)=0.
    ELSE 
      IF( plon(jl).GT.sit_domain_e) THEN
        rre=min(sit_domain_e+sit_domain_extgrd,360.)
        rrw=sit_domain_e
      ELSEIF ( plon(jl).LT.sit_domain_w) THEN
        rre=max(sit_domain_w-sit_domain_extgrd,0.)
        rrw=sit_domain_w
      ELSE
        rre=1.
        rrw=1.
      ENDIF

      IF( plat(jl).GT.sit_domain_n) THEN
        rrn=min(sit_domain_n+sit_domain_extgrd,90.)
        rrs=sit_domain_n
      ELSEIF( plat(jl).LT.sit_domain_s) THEN
        rrn=max(sit_domain_s-sit_domain_extgrd,-90.)
        rrs=sit_domain_s
      ELSE
        rrn=1.
        rrs=1.
      ENDIF
    
      IF(rrw .EQ. 1.) THEN
        rr_we=1.
      ELSE
        rr_we=1.-ABS((plon(jl)-rrw)/(rre-rrw))
      ENDIF
      IF(rrs .EQ. 1.) THEN
        rr_sn=1.
      ELSE
        rr_sn=1.-ABS((plat(jl)-rrs)/(rrn-rrs))
      ENDIF

!      rrtmp=sqrt(rr_sn**2.0+rr_we**2.0)
      rrtmp=0.5*(rr_sn+rr_we)

      IF(rrtmp .EQ. 1.) THEN
        wtr=1.
      ELSE if(rrtmp .eq. 0.) then
        wtr=0.
      ELSE
        wtr=exp(-rf/rrtmp*exp(1./(rrtmp-1.))) 
      ENDIF

      ptsw(jl)= wtr*ptsw(jl)+(1.-wtr)*pobswtb(jl)
!      pwt(jl,0)=ptsw(jl)
!      pdtswdt(jl)=(ptsw(jl)-poldtsw(jl))/(n_sit_step*zdtime)
!      pdtswdt(jl)=(pwt(jl,0)-poldsitwt(jl,0,0))/(2.*n_sit_step*zdtime)
      pdtswdt(jl)=(pwt(jl,0)-poldsitwt(jl,0,1))/(n_sit_step*zdtime)
    ENDIF
!!!end location blending

#endif
!ps  
  ELSE
    ptsw(jl)=pobswtb(jl)
!    pwt(jl,0)=ptsw(jl)
!    pdtswdt(jl)=(ptsw(jl)-poldtsw(jl))/(n_sit_step*zdtime)
    pdtswdt(jl)=(pwt(jl,0)-poldsitwt(jl,0,1))/(n_sit_step*zdtime)
  ENDIF   !ENDIF ltrigsit
    IF(GDCHK3) then
      print *,"sitvdiff:ltrigsit=",ltrigsit,",poldtsw(jl)="       &
             ,poldtsw(jl),",ptsw(jl)=",ptsw(jl) &
             ,",pdtswdt=",pdtswdt(jl),",ltrigsit=",ltrigsit       &
             ,",rrs=",rrs,",rrn=",rrn,",rrw=",rrw,",rre=",rre     &
             ,",rr_sn=",rr_sn,",rr_we=",rr_we,",wtr=",wtr

 3301 FORMAT(1X,9(A11,E13.5))

      WRITE(nerr,3301)           &
        "2sitvdifftau,",tau,",ptsw,",ptsw(jl)                     &
       ,",pobswtb,",pobswtb(jl),",pobswsb,",pobswsb(jl)           &
       ,",ptsw,",ptsw(jl),",pdtswdt,",pdtswdt(jl)                 &
       ,",pwt0,",pwt(jl,0),",pws0,",pws(jl,0),",zsf,",zsf
      ENDIF

  !
  !*   14.4 Calc the uppermost bulk layer properties for coupling with 3-D ocean
  !
  !
    IF (locn) THEN
      ! Not valid for lwaterlevel=.TRUE.
      sumxxz=0.
      sumxxt=0.
      sumxxu=0.
      sumxxv=0.
      sumxxs=0.
      DO jk=nls+1,lstfn
        ! excluidng skin layer
        sumxxz=sumxxz+hw(jk)
        sumxxt=sumxxt+pwt(jl,jk)*hw(jk)
        sumxxu=sumxxu+pwu(jl,jk)*hw(jk)
        sumxxv=sumxxv+pwv(jl,jk)*hw(jk)
        sumxxs=sumxxs+pws(jl,jk)*hw(jk)
        IF ( lwarning_msg.GE.3 ) THEN
          WRITE(nerr,*) 'jk=',jk,'hw=',hw(jk)
        ENDIF
      ENDDO
      pwtb(jl)=sumxxt/sumxxz
      pwub(jl)=sumxxu/sumxxz
      pwvb(jl)=sumxxv/sumxxz
      pwsb(jl)=sumxxs/sumxxz
  
      IF ( (sumxxz.LE.0.).OR.((pwtb(jl)+pwub(jl)+pwvb(jl)+pwsb(jl)+psubfluxw(jl)+pwsubsal(jl)).LT.-9.E20) ) THEN
        WRITE(nerr,*) ", I am in thermocline: sumxxz = (<=0)", sumxxz
        WRITE(nerr,*) 'nls+1=',nls+1,'lstfn=',lstfn,'sumxxz=',sumxxz,'sumxxt=',sumxxt,'sumxxu=',sumxxu,'sumxxv=',sumxxv,'sumxxs=',sumxxs
        WRITE(nerr,*) 'pwtb=',pwtb(jl),'pwub=',pwub(jl),'pwvb=',pwvb(jl),'pwsb=',pwsb(jl),'psubfluxw=',psubfluxw(jl),'pwsubsal=',pwsubsal(jl)
        WRITE(nerr,*) 'pobswtb=',pobswtb(jl)
        CALL output2
      ENDIF
    ENDIF
  !
  !*  14.5 Calc ground heat flux for coupling with vdiff/surftemp,
  !        net surface heat/fresh water flux into ocean, and 
  !        subsurface heat/salinity fluxes
  !  
    pgrndhflx(jl)=pgrndhflx(jl)/acc_time
    pfluxiw(jl)=pfluxiw(jl)/acc_time
    ppme2(jl)=ppme2(jl)/acc_time
    IF (locn) THEN
      psubfluxw(jl)=psubfluxw(jl)/acc_time
      pwsubsal(jl)=pwsubsal(jl)/acc_time
    ENDIF
  !
  !*  14.6 Calc heat content of a water column above tmelt
  ! 
    IF (nle.GE.1)THEN
  !     water exists
      phc(jl)=rhowcw*DOT_PRODUCT((pwt(jl,nls+1:nle)-tmelt),MAX(hw(nls+1:nle),0.))
      psc(jl)=DOT_PRODUCT(pws(jl,nls+1:nle),MAX(hw(nls+1:nle),0.))
    ELSE
  !     soil only
      phc(jl)=0.
      psc(jl)=0.
    ENDIF
    IF (lwaterlevel) THEN
      pengwac(jl)=pengwac(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl)        &
        +prsf(jl)*clw*(tmelt-ptemp2(jl))                    &
        +pssf(jl)*(alf+csn*(tmelt-ptemp2(jl)))  )
  !!!  ! Assuming that rain temp to be wtm(nls+1), that of snowfall to be tsim(1) 
  !!!    pengwac(jl)=pengwac(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl)        &
  !!!      +prsf(jl)*clw*(tmelt-wtm(nls+1))                    &
  !!!      +pssf(jl)*(alf+csn*(tmelt-tsim(1)))  )
    ELSE
    ! Neglected adveced heat flux 
      pengwac(jl)=pengwac(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl)        &
        +pssf(jl)*(alf)  )
    ENDIF
  !!!  pengw(jl)=pengw(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl) )
  !!!  pengw2(jl)=pengw2(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl)+pssf(jl)*(alf)  )
      
    IF (nle.GE.1)THEN
  !     water exists
      pwtfns(jl)=rhowcw*DOT_PRODUCT(pwtfn(jl,nls+1:nle),hw(nls+1:nle))
      pwsfns(jl)=DOT_PRODUCT(pwsfn(jl,nls+1:nle),hw(nls+1:nle))
!ps
      IF (GDCHK3) then
        WRITE(nerr,*) ", water exist"
        WRITE(nerr,*) "wtfn(",jl,")=",pwtfn(jl,nls:nle)
        WRITE(nerr,*) "hw=",hw(nls+1:nle)
        WRITE(nerr,*) "wtfns(",jl,")=",pwtfns(jl),"wsfns(",jl,")=",pwsfns(jl)
      ENDIF
!ps
    ELSE
  !     soil only
      pwtfns(jl)=0.
      pwsfns(jl)=DOT_PRODUCT(pwsfn(jl,nls+1:nle),hw(nls+1:nle))
    ENDIF
  END SUBROUTINE final
!----------------------------------------------------------  
  SUBROUTINE acc_flux
! ----------------------------------------------------------------------
!
!*   accumlate fluxes
!
! ----------------------------------------------------------------------
!     pfluxiw: net surface heat flux into ocean (W/m2, + upward)
!     ppme2: net fresh water into ocean (m/s, + downward)
!
!
!*   14.2 pgrndcapc, pgrndhflx, current, tsw
!
   pgrndhflx(jl)=pgrndhflx(jl)+pgrndhflx_int(jl)*zdtime
   pgrndflux(jl)=pgrndflux(jl)+pgrndhflx_int(jl)*(1.-pslm(jl))*zdtime
   pfluxiw(jl)=pfluxiw(jl)+pfluxiw_int(jl)*zdtime
   ppme2(jl)=ppme2(jl)+ppme2_int(jl)*zdtime
   IF (locn) THEN
     psubfluxw(jl)=psubfluxw(jl)+pwsubflux_int(jl)*zdtime
     pwsubsal(jl)=pwsubsal(jl)+pwsubsal(jl)*zdtime
   ENDIF
  END SUBROUTINE acc_flux
!----------------------------------------------------------  
SUBROUTINE thermocline(jl,jrow)
!!    input from SURF or calling routine
!v77  INTEGER LKID ! Ocean ID (Caspian Sea, Great Lake have their own ID)
!v77  INTEGER jl  ! Ocean grid (point) ID (In land mask, each ocean grid
               !   has its own ID) 
!v77  INTEGER JL   ! the corresponding longitude ID for jl
!v77  INTEGER istep ! the corresponding latitude ID for jl
!v77  real:: zdtime   ! time step in sec
!v77  INTEGER istep ! (YYMMDDHH) or (MMDDHHMM) an I8 integer to be used in 
!                Subroutine OUTPUT to add the time stamp of the resut.
!!
!! INPUT & OUTPUT VARIABLE
!!
!! Note that, there are two extra variables from SURF
!! pfluxw & pdfluxs. It is located in module mo_pgrads in this
!! routine. While coupling with ECHAM4, pfluxw and pdfluxs should
!! be included in the Calling variables.
!
!-----------------------------------------------------------------
! 
  USE mod_eos_ocean,      ONLY: tmelt, rhoh2o, alf, clw, g, alv, cpd, stbo
  IMPLICIT NONE
!
! 0.0 Calling Variables

  INTEGER, INTENT(IN):: jl, jrow    ! lonitude and latitude index
!-----------------------------------------------------------------
! 1.1 Local Integer
      
  INTEGER :: jk
  INTEGER mas ! mas: AA matrix staring level
  INTEGER mae ! mae: AA/AAC matrix last level
  INTEGER iderr ! iderr: error id. See Subroutine lkerr for details.
  INTEGER irsv  ! irsv: recursive number 
  INTEGER :: levelm,level,levelp
!
  LOGICAL lwf ! lwf: water freeze/ice melt logic at the interface  
  LOGICAL lim ! lim: ice melt logic at the surface
  LOGICAL lsm ! lsm: snow melt logic at the surface

!
  real, DIMENSION(0:lkvl+5):: X,Y,RHS,SOL
  real, DIMENSION(0:lkvl+5,3):: AA
  COMPLEX, DIMENSION(0:lkvl+1)::wumc,pawuflc
  COMPLEX, DIMENSION(0:lkvl+5)::RHSC,SOLCMPLX
  COMPLEX, DIMENSION(0:lkvl+5,3)::AAC
!
! 1.3 LOCAL VARIABLES
!
  real:: zsoflw ! net solar radiation flux over openwater per water fraction (w/m2) (positive downward)
  real:: zsofli ! net solar radiation flux over ice sheet per water fraction (w/m2) (positive downward)
!v90  real:: emin  !limit min pwtke to 3.0E-5              v78
  real:: fcew   ! calculate FCE of ice/water interface (ice melt positive)
  real:: fcei,fces ! potential phase change energy of ice and snow
  real:: xife   ! ice freezing energy due to liquid water on ice freezed (j/m2)
  real:: sfe    ! snow freezing energy due to liquid water in snow freezed (j/m2)
  real:: fcewm,tmp ! temporary working variables. ??M means previous variable.
  real:: fcewf  ! final phase change energy for water.
  
  real:: sfm    ! temporary zsf
!  real:: depth  
  REAL::hw1m
  
  COMPLEX :: tauc ! tauc: shear stress at the surface (N/m2)
  COMPLEX :: fr_corc  ! friction and corilois factor (1/s)
!
  real:: totice ! total ice above water (snow+snowfall+ice-sublimation) (m)
!
  !
  real:: pzsim       !
  real:: s      ! effective water content
  real:: hair   ! porosity of snow
  real:: hi     ! hight of irreducible water content
!
  real:: gzero  ! test of non-solar energy in Equation (12)
!!!  real:: g0          ! test of non-solar energy in Equation (12)
!
!  real:: alphaffn  !Ratio of soalr absorbtion within skin layer, v92
!
  real:: hcoolskin
! thickness of conductive sublayer (m), where only molecular transfer exists (m) 
!   (Khundzuha et al., 1977; Paulson and Simpson, 1981)
!
!
  real:: zcor          ! coriolis force   
  real:: dz            ! distance between two density levels/thickness for freeze (m)


!*    2. Initialization
!
!
!
!*    2.2 VERTICAL LEVELS (lkvl=11)
!
!v77      CALL pzcord(LKID,jl,nls,nle,z,zlk,hw,lsoil,lshlw)
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!
!     2.3 Initialization
!

!      pfluxw(jl)=-(pahflw(jl)+pahfsw(jl)+ptrflw(jl)+psoflw(jl))   ! zfluw: positive upward, which is oppsitive to lake subroutine
      zsf=0.
!     
      lwf=.FALSE.
      lim=.FALSE.
      lsm=.FALSE.
!     default all phase change logics are false
      irsv = 0   ! recursive number, this is beginning
!
!     2.4 Store Old Value (Note Soil always at nle+1 level)
!
      pfluxwm=pfluxw(jl)
      if(locaf0) then
!ps        pfluxw2=pfluxw(jl)
        pfluxw2=pfluxw(jl)-(pawtfl0(jl,0)+ocaf0_add)
      else
        pfluxw2=pfluxw(jl)
      endif
      pfluxim=pfluxi(jl)
      pfluxi2=pfluxi(jl)
      wtm(0:lkvl+1)=pwt(jl,0:lkvl+1)
      wum(0:lkvl+1)=pwu(jl,0:lkvl+1)
      wvm(0:lkvl+1)=pwv(jl,0:lkvl+1)
      wsm(0:lkvl+1)=pws(jl,0:lkvl+1)
      wtkem(0:lkvl+1)=pwtke(jl,0:lkvl+1)
      
      IF(lsice_nudg.AND.(i_sit_step.EQ.1)) THEN
        IF(ptsi(jl).NE.xmissing) ptsnic(jl,0:3)=ptsi(jl)
        IF(psni(jl).NE.xmissing) pzsi(jl,0)=psni(jl)
        IF(psiced(jl).NE.xmissing) pzsi(jl,1)=psiced(jl)
      ENDIF
!!!      IF(lsit_ice.AND.(i_sit_step.EQ.1)) THEN
!!!        IF(ptsi(jl).NE.xmissing) ptsnic(jl,0:3)=ptsi(jl)
!!!        IF(psni(jl).NE.xmissing) pzsi(jl,0)=psni(jl)
!!!        IF(psiced(jl).NE.xmissing) pzsi(jl,1)=psiced(jl)
!!!      ENDIF
      tsim(0:3)=ptsnic(jl,0:3)
 
      wlvlm=pwlvl(jl)
      ppme2_int(jl)=0.
      IF (.NOT.lsoil) THEN
        hw1m=hw(nls+1)
        IF (lhd) THEN
!!!          hw(nls+1)=hw(nls+1)+pruntoc(jl)*zdtime/rhoh2o     
          hw(nls+1)=hw(nls+1)+pdisch(jl)*zdtime
        ENDIF
      ENDIF
!
!     calc cold content ccm (J/m2): the energy to melt snow + ice
!
      IF (lwaterlevel) THEN
      !! take advected heat flux from snow into account
        ccm=pcc(jl)+ zdtime*(                                           &
&         ( prsf(jl) )*clw*                                  &
&         ( tmelt-ptemp2(jl) ) +                                        &
&         ( pssf(jl) )*                                      &
          ( alf+csn*(tsim(1)-ptemp2(jl)) )                              &
         )
!       Since the advection energy of snowfall has been included in
!       pfluxi2, the calculated Tsn will equal to temp2 (temp of snowfall).
      ELSE
      !! assuming no advected heat flux from snowfall and rainfall
        ccm=( pssf(jl) )*zdtime*( alf ) &
          +pcc(jl)
      ENDIF
!-----------------------------------------------------------------
!
!*    3. Precipitation & Evaporation (Sublimation) Events
!
!-----------------------------------------------------------------
  300 CONTINUE
!      PRINT *, ", I am in Precipitation & Evaporation Events."
!     3.1 Calc total ice above water (soil)
      totice=pzsi(jl,0)+pzsi(jl,1)&
&       +( pssf(jl)+pevapw(jl))*zdtime/rhoh2o
!
      IF ((totice-pzsi(jl,1).GT.0.).AND.&
&             ( (pzsi(jl,1).GT.0.).OR.lsoil)&
&                )  THEN
!
!*    3.2 snow over ice/ snow over land 
!
        mas=0
        IF (lwaterlevel) THEN
        !! take advected heat flux from snow into account
          pfluxi2=pfluxi2+&
&           ( prsf(jl) )*clw*&
&           (tmelt-ptemp2(jl) ) +&
&           ( pssf(jl) )*csn*&
&           (tsim(1)-ptemp2(jl))
!         Since the advection energy of snowfall has been included in
!         pfluxi2, the calculated Tsn will equal to temp2 (temp of snowfall).
        ENDIF
        pzsi(jl,0)=totice-pzsi(jl,1)
        psilw(jl,0)=psilw(jl,0)+prsf(jl)*zdtime/rhoh2o
      ELSEIF ((totice .GT. 0.).AND.(pzsi(jl,1).GT.0.)) THEN
!
!*    3.3 Ice on top or snow sublimates in this time step completely
!
        mas=2
        IF (lwaterlevel) THEN
        !! take advected heat flux from snow into account
          pfluxi2=pfluxi2                                         &
&          +prsf(jl)*clw*(tmelt-ptemp2(jl))          &
&          +pssf(jl)*                                &
&           ( csn*(tmelt-ptemp2(jl))                              &
&               +cice*(tsim(3)-tmelt) )                           &
&          +pzsi(jl,0)*rhoh2o/zdtime*                             &
&           ( csn*(tmelt-tsim(1))+cice*(tsim(3)-tmelt) )
        ELSE
          pfluxi2=pfluxi2                                         &
&          +pzsi(jl,0)*rhoh2o/zdtime*                             &
&           ( csn*(tmelt-tsim(1))+cice*(tsim(3)-tmelt) )
        ENDIF
        pzsi(jl,1)=totice
!       merge snow, snowfall into ice & assume sublimation only
        psilw(jl,1)=psilw(jl,1)+psilw(jl,0)&
&         +prsf(jl)*zdtime/rhoh2o
!       If there was snow, reset snow progonastic variables & ptsw.
        IF (pzsi(jl,0).GT.0.) THEN
          pzsi(jl,0)=0.
          psilw(jl,0)=0.
!!!          ptsnic(jl,0)=tmelt
!!!          ptsnic(jl,1)=tmelt
!!! set to be undeneath temp instead (v7.7)
        ENDIF
      ELSE
!
!*    3.4 Water/Soil on top, snow & ice sublimate completely in this
!         time step, or snow on water  
!         starting level
        mas=4
        IF (.NOT.lsoil) THEN
!       water on top (Normal & Shallow Water)
!         ref temp of water (1st layer water temp)
          tmp=wtm(nls+1)
        ELSE
!       soil on top
!         ref temp of soil (soil skin temp)
          tmp=wtm(nle+1)
        ENDIF 
        IF (lwaterlevel) THEN
        !! take advected heat flux from snow into account
#ifdef DEBUG
        if(GDCHK2) print*,"ps: 3.4.0 pfluxw2=",pfluxw2 
#endif
          pfluxw2=pfluxi2+pfluxw2                                      &
&           +( prsf(jl) )*clw*                              &
&              (tmp-ptemp2(jl))                                        &
&           +( pssf(jl) )*                                  &
&                ( alf+csn*(tmelt-ptemp2(jl))                          &
&                     +clw*(tmp-tmelt) )                               &
&           +pzsi(jl,0)*rhoh2o/zdtime*                                 &
&                ( alf+csn*(tmelt-tsim(1))+clw*(tmp-tmelt) )           &
&           +psilw(jl,0)*rhoh2o/zdtime* clw*(tmp-tmelt)                &
&           +pzsi(jl,1)*rhoh2o/zdtime*                                 &
&                ( alf+cice*(tmelt-tsim(3))+clw*(tmp-tmelt) )          &
&           +psilw(jl,1)*rhoh2o/zdtime*clw*(tmp-tmelt)
#ifdef DEBUG
       if(GDCHK2) print*,"ps: 3.4.1 pfluxi2=",pfluxi2,",prsf(jl)=",prsf(jl),",pssf(jl)=",pssf(jl)
       if(GDCHK2) print*,"ps: 3.4.1 tmp=",tmp,",ptemp2(jl)=",ptemp2(jl),",clw=",clw
       if(GDCHK2) print*,"ps: 3.4.1 alf=",alf,",csn=",csn,",tmelt=",tmelt,",rhoh2o/zdtime=",rhoh2o/zdtime
       if(GDCHK2) print*,"ps: 3.4.1 psilw(jl,0)=",psilw(jl,0),",psilw(jl,1)=",psilw(jl,1)
       if(GDCHK2) print*,"ps: 3.4.1 pzsi(jl,0)=",pzsi(jl,0),",pzsi(jl,1)=",pzsi(jl,1),",cice=",cice
       if(GDCHK2) print*,"ps: 3.4.1 tsim(1)=",tsim(1),",tsim(3)=",tsim(3)
       if(GDCHK2) print*,"ps: 3.4.1 pfluxw2=",pfluxw2 
#endif
        ELSE
#ifdef DEBUG
        if(GDCHK2) print*,"ps: 3.4.2 pfluxw2=",pfluxw2 
#endif
          pfluxw2=pfluxi2+pfluxw2                                      &
&           +( pssf(jl) )* alf                              &
&           +pzsi(jl,0)*rhoh2o/zdtime*                                 &
&                ( alf+csn*(tmelt-tsim(1))+clw*(tmp-tmelt) )           &
&           +psilw(jl,0)*rhoh2o/zdtime* clw*(tmp-tmelt)                &
&           +pzsi(jl,1)*rhoh2o/zdtime*                                 &
&                ( alf+cice*(tmelt-tsim(3))+clw*(tmp-tmelt) )          &
&           +psilw(jl,1)*rhoh2o/zdtime*clw*(tmp-tmelt)
#ifdef DEBUG
       if(GDCHK2) print*,"ps: 3.4.3 pfluxi2=",pfluxi2,",prsf(jl)=",prsf(jl),",pssf(jl)=",pssf(jl)
       if(GDCHK2) print*,"ps: 3.4.3 tmp=",tmp,",ptemp2(jl)=",ptemp2(jl),",clw=",clw
       if(GDCHK2) print*,"ps: 3.4.3 alf=",alf,",csn=",csn,",tmelt=",tmelt,",rhoh2o/zdtime=",rhoh2o/zdtime
       if(GDCHK2) print*,"ps: 3.4.3 psilw(jl,0)=",psilw(jl,0),",psilw(jl,1)=",psilw(jl,1)
       if(GDCHK2) print*,"ps: 3.4.3 pzsi(jl,0)=",pzsi(jl,0),",pzsi(jl,1)=",pzsi(jl,1),",cice=",cice
       if(GDCHK2) print*,"ps: 3.4.3 tsim(1)=",tsim(1),",tsim(3)=",tsim(3)
       if(GDCHK2) print*,"ps: 3.4.3 pfluxw2=",pfluxw2 
#endif
        ENDIF
        pfluxi2=0.
!         includes advection flux (positive upward) 
!         assume precipitation temperature is temp2 &
!         snowfall melted released latent heat
!         assume wtm(nls+1) is the temperature of outflow
!         merge snow/ice layers into water
!
        IF (.NOT.lsoil) THEN
          hw(nls+1)=hw(nls+1)+&
&           (prsf(jl)+pssf(jl)+pevapw(jl))*zdtime/rhoh2o&
&           +pzsi(jl,0)+psilw(jl,0)+pzsi(jl,1)+psilw(jl,1)
          IF (.FALSE.) THEN
            zsf=zsf+(MAX(hw(nls+1),0.)-hw1m)*wsm(nls+1)
          ELSE
            zsf=zsf+&
&             ((prsf(jl)+pssf(jl)+pevapw(jl))*zdtime/rhoh2o+pzsi(jl,0)+psilw(jl,0))*wsm(nls+1)&
&             +(pzsi(jl,1)+psilw(jl,1))*salti*wsm(nls+1)
          ENDIF
!
          DO WHILE (hw(nls+1).LE.0.AND.nle.GE.nls+2)
            hw1m=hw(nls+2)
            IF ((nle-nls).EQ.2) THEN
              hw(nls+2)=MAX(hw1m+hw(nls+1),wcri)
!             preserve wcri for shallow water mode
!             Although it is improper for water conservation, 
!             it's important for salinity conservation.
            ELSE
              hw(nls+2)=hw1m+hw(nls+1)
            ENDIF
            zsf=zsf+(MAX(hw(nls+2),0.)-hw1m)*wsm(nls+2)
            hw(nls+1)=0.
            nls=nls+1
          ENDDO
          IF ((nle-nls).EQ.1) THEN
            lshlw = .TRUE.
!           shallow water body
          ENDIF
!       the minmum thickness of wcri is reserved to perserve 
!       salinity conservation.
        ELSE
          IF (lwaterlevel) THEN
            pwlvl(jl)=pwlvl(jl)+&
&             (prsf(jl)+pssf(jl)+pevapw(jl))*zdtime/rhoh2o&
&             +pzsi(jl,0)+psilw(jl,0)+pzsi(jl,1)+psilw(jl,1)
          ENDIF
        ENDIF
!
!       If there was snow or ice, reset snow/ice progonastic
!       variables & ptsw.
        IF (SUM(pzsi(jl,0:1))+SUM(psilw(jl,0:1)).GT.0.) THEN
          pzsi(jl,:)=0.
          psilw(jl,:)=0.
!!!          ptsnic(jl,:)=tmelt
!!!          ptsnic(jl,:)=tmelt
!!! set to be undeneath temp instead (v7.7)
        ENDIF
      ENDIF

!*    3.42 Calculate Effective thickness of snow and ice. Note the the effective thickness can 
!         be zero although its physcial thickness is not due to numerical error.

      CALL update_snow_ice_property(pzsi(jl,0),pzsi(jl,1),hsn,hesn,hice,heice,pseaice(jl))


      !!! hsn=HSNFN(pzsi(jl,0))
      !!! hesn=HEFN(hsn/4.,xksn,omegas)
      !!! hice=HICEFN(pzsi(jl,1))
      !!! heice=HEFN(hice/4.,xkice,omegas)
      !!! pseaice(jl)=SEAICEFN(pzsi(jl,1))          
!
!*    3.5 Determine Solar radaition on Water. This value is fixed at
!         current stage to prevent recursive. The model might not be
!         able to determine icemelt or water freeze IF it change
!         accordingly. 
!
!     surface net solar radiation flux over water
!
      
      IF ( (heice.LE.xicri).AND.(hesn.LE.csncri) ) THEN
!     no/or only thin snow and ice on top
        zsoflw=psofli(jl)+psoflw(jl)
        zsofli=0.
      ELSE
        zsoflw=psoflw(jl)
        zsofli=psofli(jl)
      ENDIF
 
      IF (mas.EQ.0) THEN
! ----------------------------------------------------------------------
!
!*    4. Snow & Ice Melt Runoff
!        Calc Liquid Water Balance in Snow
!        One-time step is assumed to allow liquid water becoming runoff
!        before it refreezes.
!
! ----------------------------------------------------------------------
  400 CONTINUE
!      PRINT *, ", I am in PAHSE Change Energy."
!
!*    4.1 Chk consistency
!
        IF (hesn.LE.csncri) THEN
!       Small value of snow thickness. Assume all liquid water becoming 
!       runoff to prevent numerical error.
          psilw(jl,1)=psilw(jl,1)+MAX(psilw(jl,0),0.)
!         reset snow progonastic variables
          psilw(jl,0)=0.
        ELSE
!
!     4.2 Calc porosity (hair), irreducible water content (hi),
!         effective water content (s)
!
!         Calc porosity
          hair=pzsi(jl,0)*(rhoh2o/rhosn-rhoh2o/rhoice)
!         irreducible water content
!         csi=0.03-0.07. A value of 0.05 is used
!         hi = hight of irreducible water content
          hi = 0.05*hair
!     4.3 Chk water content > irreducible water content?
          IF (psilw(jl,0).GE.hi) THEN
            IF( psilw(jl,0).LT.hair) THEN
              s=(psilw(jl,0)-hi)/(hair-hi)
            ELSE
!             more liquid water than porosity
!             assume it becomes runoff & store in ice layer
              psilw(jl,1)=psilw(jl,1)+(psilw(jl,0)-hair)
              psilw(jl,0)=hair
              s=1.
            ENDIF
!!
!!    4.4 Calc Runoff (tmp) according to Darcy's Law (Shimizu's Formula)
!!
            tmp=MIN(5.47E6*(0.077*0.9E-3**2)*&
&               Exp(-7.8*rhosn/rhoh2o)*s**3*rhoh2o*zdtime,&
&               psilw(jl,0)-hi)
!           where grain size of 0.9E-3 M is assumed.
!           DTY*G/MU=ALPHA=5.47E6
            psilw(jl,0)=psilw(jl,0)-tmp
            psilw(jl,1)=psilw(jl,1)+tmp
          ELSE
!           water content less than irreducible water content. No runoff.
          ENDIF
        ENDIF
      ENDIF
!
      IF(mas.LE.2) THEN
! ----------------------------------------------------------------------
!
!*    5. Calc Liquid Water Balance in Ice
!
! ----------------------------------------------------------------------
!     5.1 Chk LAGER THAN MAX WATER ON THE TOP OF ICE LAYER
!
        IF (heice.GT.xicri) THEN
          IF (psilw(jl,1).GT.wicemx) THEN
!           ice melt runoff occurs
            IF(.NOT.lsoil) THEN
              hw(nls+1)=hw(nls+1)+psilw(jl,1)-wicemx
              zsf=zsf+(psilw(jl,1)-wicemx)*salti*wsm(nls+1)
            ELSE
              IF(lwaterlevel) THEN
                pwlvl(jl)=pwlvl(jl)+psilw(jl,1)-wicemx
              ENDIF
            ENDIF
            psilw(jl,1)=wicemx
          ENDIF
        ENDIF
      ENDIF
!
! ----------------------------------------------------------------------
!
!*    6. Determine phase change energy of snow and ice & modify snow & ice
!        thickness due to refreezing of liquid water. Assume the refreezing
!        occurs at the center (not surface) of snow, and at the surface of
!        ice. The refreezing will change temperature profiles, which will
!        be evaluated in section 8. Note that, putting this routine after
!        section 6 implies allowing one time step for snow and ice melt
!        becoming runoff before it might refreeze in this section.
!
! ----------------------------------------------------------------------
!
  600 CONTINUE
!      PRINT *, "I am update SN & ICE due to refreeze."
!     6.1 Snow

      IF ( (tsim(1).LT.tmelt).AND.(psilw(jl,0).GT.0.) ) THEN
        sfe=MIN( (tmelt-tsim(1))*rhoh2o*pzsi(jl,0)*csn,&
&         psilw(jl,0)*rhoh2o*alf)
!       refreezing amount
        pzsi(jl,0)=pzsi(jl,0)+sfe/rhoh2o/alf
        !!! hsn=HSNFN(pzsi(jl,0))
        !!! hesn=HEFN(hsn/4.,xksn,omegas)
        psilw(jl,0)=MAX(psilw(jl,0)-sfe/rhoh2o/alf,0.)
      ELSE
        sfe=0.
      ENDIF
!
!     6.2 Ice
!
 620  CONTINUE
      IF ((psilw(jl,1).GT.0.).AND.(heice.GT.xicri)) THEN
        IF (pzsi(jl,0).GT. csncri) THEN
          xife=&
&           +zdtime*rhosn*csn*xksn*(-(tsim(1)-tmelt) )&
&             /(0.5*hsn/pseaice(jl))&
&           +zdtime*rhoh2o*cice*xkice*(tmelt-tsim(3))&
&             /(pzsi(jl,1)*rhoh2o/rhoice/2.)
        ELSE
          xife=&
&           +zdtime*rhoh2o*cice*xkice*(tmelt-tsim(3))&
&             /(pzsi(jl,1)*rhoh2o/rhoice/2.)
        ENDIF
        xife=MIN( xife, psilw(jl,1)*rhoh2o*alf)
!       refreezing amount
        pzsi(jl,1)=pzsi(jl,1)+xife/rhoh2o/alf
        !!! hice=HICEFN(pzsi(jl,1))
        !!! heice=HEFN(hice/4.,xkice,omegas)     
        psilw(jl,1)=MAX(psilw(jl,1)-xife/rhoh2o/alf,0.)
      ELSE
        xife=0.
      ENDIF
      CALL update_snow_ice_property(pzsi(jl,0),pzsi(jl,1),hsn,hesn,hice,heice,pseaice(jl))

      
!
!-----------------------------------------------------------------
!
!*    7.  COMPUTE eddy mixing coeff. pwkm, pwlmx, pwldisp
!
!-----------------------------------------------------------------
  700 CONTINUE
!      PRINT *, ", I am in computing eddy diffusivity."
!
      CALL eddy(hcoolskin,utauw,wtkem,wtm,wsm)
!
!-----------------------------------------------------------------
!
!*    8. PREPARE T MATRIX
!
!-----------------------------------------------------------------
  800 CONTINUE
!      PRINT *, ", I am in Prepare T matrix."
!
!       first level matrix updated by this routine
!
!*    8.1 SNOW LAYER
!
!!! (v7.7)
!!!
  810 CONTINUE
!!!
!!!  810 IF ( mas.EQ.0) THEN
!     snow exists
      IF (hesn.GT.csncri) THEN
!     thick snow
        Y(0) = zdtime/hesn*xksn/(0.5*hsn/pseaice(jl))
        RHS(0)= tsim(0)-(1.-beta)*Y(0)*(tsim(0)-tsim(1))&
&         +zdtime*( -pfluxi2+pdfluxs(jl)*tsim(0) )/rhosn/csn/hesn
        AA(0,1)=0.
        AA(0,3)=-beta*Y(0)
        AA(0,2)= 1.+zdtime*pdfluxs(jl)/rhosn/csn/hesn                        &
&          - AA(0,1) - AA(0,3)
!       First Layer
        Y(1) = zdtime/hsn*xksn/(0.5*hsn/pseaice(jl))
        X(1)= Y(1)
        RHS(1)= tsim(1)-hesn/hsn*tsim(0) &
&         +sfe/rhosn/csn/hsn+(1.-beta)*&
&         ( X(1)*(tsim(0)-tsim(1))-Y(1)*(tsim(1)-tsim(2)) )
        AA(1,1)=-beta*X(1)-hesn/hsn
        AA(1,3)=-beta*Y(1)
        AA(1,2)= 1. +beta*X(1)+beta*Y(1)
!
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 810 thick, zdtime=",zdtime,", hesn=",hesn,", hsn=",hsn
        if(GDCHK2) print *,"ps: 810 thick, pseaice(jl)=", pseaice(jl),",tsim(0)=", tsim(0)
        if(GDCHK2) print *,"ps: 810 thick, beta=",beta,", tsim(1)=", tsim(1)
        if(GDCHK2) print *,"ps: 810 thick, pfluxi2=", pfluxi2, ",pdfluxi(jl)=",pfluxi(jl)
        if(GDCHK2) print *,"ps: 810 thick, rhosn=",rhosn,", csn=",csn,", xksn=",xksn,", sfe=",sfe 
#endif
      ELSE
!     thin snow
!     assume snow temperature equals to the temperature underneath
        RHS(0)=0.
        AA(0,1)=0.
        AA(0,3)=-1.
        AA(0,2)=1.
!
        RHS(1)=0.
        AA(1,1)=0.
        AA(1,3)=-1.
        AA(1,2)=1.
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 810 thin"
#endif
      ENDIF
!!! (v7.7)
!!!      ENDIF
!
!*    8.2 ICE LAYER
!
  820 CONTINUE
!!!  
!!!  820 IF (mas.LE.2) THEN
      IF (heice.GT.xicri) THEN
!     thick ice
        IF (hesn.GT.csncri) THEN
!       thick snow on top
!          hsn=HSNFN(pzsi(jl,0))
          X(2) = zdtime/(rhoice*cice*heice)*&
&           (rhosn*csn*xksn)/(0.5*hsn/pseaice(jl))
          Y(2) = zdtime/heice*xkice/(0.5*hice/pseaice(jl))
          RHS(2)= tsim(2)+(1.-beta)*&
&           ( X(2)*(tsim(1)-tsim(2))-Y(2)*(tsim(2)-tsim(3)) )&
&           +xife/rhoice/cice/heice
          AA(2,1)=-beta*X(2)
          AA(2,3)=-beta*Y(2)
          AA(2,2)= 1. - AA(2,1) - AA(2,3)
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 820 thick, zdtime=",zdtime,", rhoice=",rhoice,", cice=",cice,", heice=",heice
        if(GDCHK2) print *,"ps: 820 thick, pseaice(jl)=", pseaice(jl),", beta=",beta
        if(GDCHK2) print *,"ps: 820 thick, tsim(1)=", tsim(1),", tsim(2)=", tsim(2),", tsim(3)=", tsim(3)
        if(GDCHK2) print *,"ps: 820 thick, hsn=",hsn
        if(GDCHK2) print *,"ps: 820 thick, rhosn=",rhosn,", csn=",csn,", xksn=",xksn,", xife=",xife 
#endif
        ELSE
!       thin snow or no snow on top
          Y(2) = zdtime/heice*xkice/(0.5*hice/pseaice(jl))
          RHS(2)= tsim(2)-(1.-beta)*Y(2)*(tsim(2)-tsim(3))               &
            +(  rhosn*csn*hsn*tsim(1)                                      &
                +zdtime*(-pfluxi2+pdfluxs(jl)*tsim(2))+xife  )              &
             /rhoice/cice/heice                                             
          AA(2,1)=0.                                                     
          AA(2,3)=-beta*Y(2)                                                
          AA(2,2)=1.-AA(2,1)-AA(2,3)                                     &
            +(rhosn*csn*hsn+zdtime*pdfluxs(jl))/rhoice/cice/heice
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 820 thin, zdtime=",zdtime,", heice=",heice,", xkice=",xkice,", hice=",hice
        if(GDCHK2) print *,"ps: 820 thin, pseaice(jl)=", pseaice(jl),", beta=",beta
        if(GDCHK2) print *,"ps: 820 thin, tsim(1)=", tsim(1),", tsim(2)=", tsim(2),", tsim(3)=", tsim(3)
        if(GDCHK2) print *,"ps: 820 thin, hsn=",hsn,", pfluxi2=", pfluxi2, ",pdfluxi(jl)=",pfluxi(jl)
        if(GDCHK2) print *,"ps: 820 thin, rhosn=",rhosn,", csn=",csn,", xksn=",xksn,", xife=",xife 
#endif
!
        ENDIF
!       First Layer
        Y(3) = zdtime/hice*xkice/(0.5*hice/pseaice(jl))
        X(3)= Y(3)
        RHS(3)= tsim(3)-heice/hice*tsim(2)+(1.-beta)*&
&         ( X(3)*(tsim(2)-tsim(3))-Y(3)*(tsim(3)-wtm(0)) )
        AA(3,1)=-beta*X(3)-heice/hice
        AA(3,3)=-beta*Y(3)
        AA(3,2)= 1. +beta*X(3)+beta*Y(3)
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 820 first, wtm(0)=", wtm(0)
#endif
      ELSE
!     thin ice
!     assume ice temperature equals to the temperature underneath
        RHS(2)=0.
        AA(2,1)=0.
        AA(2,3)=-1.
        AA(2,2)=1.
!
        RHS(3)=0.
        AA(3,1)=0.
        AA(3,3)=-1.
        AA(3,2)=1.
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 820 else"
#endif
      ENDIF
!!!      ENDIF
!
!*    8.3 WATER & SOIL LAYER 
!
  830 IF (.NOT.lsoil) THEN
  !
  !*   Total Fresh Water Inflow = (pwlvl(jl)-wlvlm)/zdtime
  !
        ppme2_int(jl)=ppme2_int(jl)+(pbathy(jl)+SUM(hw(nls+1:nle))-wlvlm)/zdtime
        IF (lwaterlevel) THEN
          pwlvl(jl)=pbathy(jl)+SUM(hw(nls+1:nle))
        ELSE
        ! restore back to original water level, and associtaed hw
          pwlvl(jl)=wlvlm
        ENDIF
        CALL pzcord(pwlvl(jl),pbathy(jl),nls,nle,z,zlk,hw,lsoil,lshlw)  ! bjt 2010/2/21
!
        CALL LKDIFKH(hew,X,Y)
!       Skin layer
        IF (lssst) THEN
          IF ( (hesn.LE.csncri).AND.(heice.LE.xicri)) THEN
!           water on top or only thin snow/ice layer exists
            RHS(4)= wtm(0)-(1.-beta)*Y(4)*(wtm(0)-wtm(nls+1))              &
                    +(  rhosn*csn*hsn*tsim(1)+rhoice*cice*hice*tsim(3)        &
                    +zdtime*( zsoflw*(FFN(0.)-FFN(zlk(nls+1)-pwlvl(jl)))   &
                      -(pfluxw2+zsoflw)+pdfluxs(jl)*wtm(0) )                  &
                     )/rhoh2o/clw/hew+zdtime*pawtfl(jl,0)
            AA(4,1)=0.
            AA(4,3)=-beta*Y(4)
            AA(4,2)= 1.-AA(4,1)-AA(4,3)                                    &
                    +(  rhosn*csn*hsn+rhoice*cice*hice+zdtime*pdfluxs(jl)  )  &
                    /rhoh2o/clw/hew
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 830 .1, Y(4)=",Y(4)
        if(GDCHK2) print *,"ps: 830 .1, nls=",nls,", wtm(0)=",wtm(0),", wtm(nls+1)=",wtm(nls+1)
        if(GDCHK2) print *,"ps: 830 .1, zdtime=",zdtime,", hice=",hice,", zsoflw=",zsoflw
        if(GDCHK2) print *,"ps: 830 .1, zlk(nls+1)=",zlk(nls+1),", pwlvl(jl)=", pwlvl(jl)
        if(GDCHK2) print *,"ps: 830 .1, beta=", beta,", tsim(3)=", tsim(3),", rhoice=",rhoice,", cice=",cice
        if(GDCHK2) print *,"ps: 830 .1, hsn=",hsn,", pfluxw2=", pfluxw2, ",pdfluxw(jl)=",pfluxw(jl)
        if(GDCHK2) print *,"ps: 830 .1, rhosn=",rhosn,", csn=",csn,", rhoh2o=",rhoh2o,", clw=",clw,", hew=",hew 
#endif
!!!           IF (lsit_lw) THEN
!!!             AA(4,3)=AA(4,3)-4.*stbo*zdtime*wtm(nls+1)**3/rhoh2o/clw/hew
!!!             AA(4,2)=AA(4,2)+4.*stbo*zdtime*wtm(0)**3/rhoh2o/clw/hew
!!!           ENDIF
          ELSE
!         thick ice or thick snow on top
!           set pwt(0) =pctfreez2(jl)
            RHS(4)= pctfreez2(jl)
            AA(4,1)=0.
            AA(4,3)=0.
            AA(4,2)=1. - AA(4,1) - AA(4,3)
#ifdef DEBUG
            if(GDCHK2) print *,"ps: 830 .2"
#endif		
          ENDIF

        ELSE
!         Without skin layer, skin temperature = first layer temperature
!         so AA(4,3)=-1 & AA(4,2)= 1.  
          RHS(4)=0.
          AA(4,1)=0.
          AA(4,3)=-1.
          AA(4,2)= 1.
        ENDIF
!       First Layer
        jk = 1
        IF (lssst) THEN
          RHS(4+jk)= wtm(nls+jk)-hew/hw(nls+jk)*wtm(0)+(1.-beta)*                &
              ( X(4+jk)*(wtm(0)-wtm(nls+jk))-Y(4+jk)*(wtm(nls+jk)-wtm(nls+jk+1)) )  &
            +zdtime*(                                                               &
              zsoflw*( FFN(zlk(nls+jk)-pwlvl(jl))-FFN(zlk(nls+jk+1)-pwlvl(jl)) )    &
                /rhoh2o/clw/hw(nls+jk)                                              &
              +pawtfl(jl,nls+jk)                     &
              )

          AA(4+jk,1)=-beta*X(4+jk)-hew/hw(nls+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1. +beta*X(4+jk)+beta*Y(4+jk)
          IF (lsit_lw) THEN
            AA(4+jk,3)=AA(4+jk,3)-4.*stbo*zdtime*wtm(nls+jk+1)**3/rhoh2o/clw/hw(nls+jk)
            AA(4+jk,2)=AA(4+jk,2)+4.*stbo*zdtime*wtm(nls+jk)**3/rhoh2o/clw/hw(nls+jk)
          ENDIF
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 830 .3, jk=",jk
        if(GDCHK2) print *,"ps: 830 .3, X(4+jk)=",X(4+jk),", Y(4+jk)=",Y(4+jk)
        if(GDCHK2) print *,"ps: 830 .3, nls=",nls,", wtm(0)=",wtm(0),", wtm(nls+jk)=",wtm(nls+jk),", wtm(nls+jk+1)=",wtm(nls+jk+1)
        if(GDCHK2) print *,"ps: 830 .3, zdtime=",zdtime,", zsoflw=",zsoflw,", hew=",hew
        if(GDCHK2) print *,"ps: 830 .3, FFN(zlk(nls+jk)-pwlvl(jl))=",FFN(zlk(nls+jk)-pwlvl(jl)),", FFN(zlk(nls+jk+1)-pwlvl(jl))=",FFN(zlk(nls+jk+1)-pwlvl(jl))
        if(GDCHK2) print *,"ps: 830 .3, beta=", beta,", hw(nls+jk)=", hw(nls+jk)
        if(GDCHK2) print *,"ps: 830 .3, pawtfl(jl,nls+jk)=", pawtfl(jl,nls+jk)
        if(GDCHK2) print *,"ps: 830 .3, pawtfl0(jl,nls+jk)=", pawtfl0(jl,nls+jk)
        if(GDCHK2) print *,"ps: 830 .3, rhosn=",rhosn,", csn=",csn,", rhoh2o=",rhoh2o,", clw=",clw
        if(GDCHK2) print *,"ps: 830 .3, RHS(4+jk)=",RHS(4+jk) 
#endif
        ELSE
          IF ( (hesn.LE.csncri).AND.(heice.LE.xicri)) THEN
            RHS(4+jk)= wtm(nls+jk)-(1.-beta)*Y(4+jk)*(wtm(nls+jk)-wtm(nls+jk+1))  &
                    +(  rhosn*csn*hsn*tsim(1)+rhoice*cice*hice*tsim(3)               &
                       +zdtime*( zsoflw*(FFN(0.)-FFN(zlk(nls+jk+1)-pwlvl(jl)))    &
                         -(pfluxw2+zsoflw)+pdfluxs(jl)*wtm(0)                        &
                        )                                                            &
                      )/rhoh2o/clw/hw(nls+jk)                                        &
                       +zdtime*pawtfl(jl,nls+jk)

            AA(4+jk,1)=0.
            AA(4+jk,3)=-beta*Y(4+jk)
            AA(4+jk,2)= 1. - AA(4+jk,1) - AA(4+jk,3)                              &
                    +(  rhosn*csn*hsn+rhoice*cice*hice+zdtime*pdfluxs(jl)  )         &
                     /rhoh2o/clw/hw(nls+jk)
            IF (lsit_lw) THEN
              AA(4+jk,3)=AA(4+jk,3)-4.*stbo*zdtime*wtm(nls+jk+1)**3/rhoh2o/clw/hw(nls+jk)
              AA(4+jk,2)=AA(4+jk,2)+4.*stbo*zdtime*wtm(nls+jk)**3/rhoh2o/clw/hw(nls+jk)
            ENDIF
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 830 .4, X(4+jk)=",X(4+jk),", Y(4+jk)=",Y(4+jk)
        if(GDCHK2) print *,"ps: 830 .4, nls=",nls,", wtm(0)=",wtm(0),", wtm(nls+jk)=",wtm(nls+jk),", wtm(nls+jk+1)=",wtm(nls+jk+1)
        IF(GDCHK2) print *,"ps: 830 .4, tsim(3)=", tsim(3),", rhoice=",rhoice,", cice=", cice,", hice=",hice
        if(GDCHK2) print *,"ps: 830 .4, zdtime=",zdtime,", zsoflw=",zsoflw,", hew=",hew,", hw=",hw
        if(GDCHK2) print *,"ps: 830 .4, FFN(0.)=",FFN(0.),", FFN(zlk(nls+jk+1)-pwlvl(jl))=",FFN(zlk(nls+jk+1)-pwlvl(jl))
        if(GDCHK2) print *,"ps: 830 .4, beta=", beta,", hw(nls+jk)=", hw(nls+jk)
        if(GDCHK2) print *,"ps: 830 .4, pfluxw2=",pfluxw2,", pdfluxs(jl)=",pdfluxs(jl)
        if(GDCHK2) print *,"ps: 830 .4, rhosn=",rhosn,", csn=",csn,", rhoh2o=",rhoh2o,", clw=",clw,", hsn=",hsn
        if(GDCHK2) print *,"ps: 830 .4, RHS(4+jk)=",RHS(4+jk) 
#endif
          ELSE
!         thick ice or thick snow on top
!           set pwt(jk) =pctfreez2(jl)
            RHS(4+jk)= pctfreez2(jl)
            AA(4+jk,1)=0.
            AA(4+jk,3)=0.
            AA(4+jk,2)= 1. - AA(4+jk,1) - AA(4+jk,3)
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 830 .5"
#endif
          ENDIF
        ENDIF
!       Middle & Bottom Layers
        DO jk =2,nle-nls
          RHS(4+jk)= wtm(nls+jk)+(1.-beta)*&
&             ( X(4+jk)*(wtm(nls+jk-1)-wtm(nls+jk))-&
&               Y(4+jk)*(wtm(nls+jk)-wtm(nls+jk+1)) )+&
&           zdtime*( zsoflw* &
&             ( FFN(zlk(nls+jk)-pwlvl(jl))-FFN(zlk(nls+jk+1)-pwlvl(jl)) )&
&               /rhoh2o/clw/hw(nls+jk)                                   &
&             +pawtfl(jl,nls+jk) )


          AA(4+jk,1)=-beta*X(4+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1. - AA(4+jk,1) - AA(4+jk,3)
          IF (lsit_lw) THEN
            AA(4+jk,1)=AA(4+jk,1)-4.*stbo*zdtime*wtm(nls+jk-1)**3/rhoh2o/clw/hw(nls+jk)
            AA(4+jk,3)=AA(4+jk,3)-4.*stbo*zdtime*wtm(nls+jk+1)**3/rhoh2o/clw/hw(nls+jk)
            AA(4+jk,2)=AA(4+jk,2)+8.*stbo*zdtime*wtm(nls+jk)**3/rhoh2o/clw/hw(nls+jk)
          ENDIF
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 830 .6, jk=",jk
        if(GDCHK2) print *,"ps: 830 .6, X(4+jk)=",X(4+jk),", Y(4+jk)=",Y(4+jk)
        if(GDCHK2) print *,"ps: 830 .6, nls=",nls,", nle=",nle,", wtm(nls+jk-1)=",wtm(nls+jk-1),", wtm(nls+jk)=",wtm(nls+jk),", wtm(nls+jk+1)=",wtm(nls+jk+1)
        if(GDCHK2) print *,"ps: 830 .6, tsim(3)=", tsim(3),", rhoice=",rhoice,", cice=", cice,", hice=",hice
        if(GDCHK2) print *,"ps: 830 .6, zdtime=",zdtime,", zsoflw=",zsoflw,", hew=",hew,", hw=",hw
        if(GDCHK2) print *,"ps: 830 .6, FFN(zlk(nls+jk)-pwlvl(jl))=",FFN(zlk(nls+jk)-pwlvl(jl)),", FFN(zlk(nls+jk+1)-pwlvl(jl))=",FFN(zlk(nls+jk+1)-pwlvl(jl))
        if(GDCHK2) print *,"ps: 830 .6, beta=", beta,", hw(nls+jk)=", hw(nls+jk)
        if(GDCHK2) print *,"ps: 830 .6, pfluxw2=",pfluxw2,", pdfluxs(jl)=",pdfluxs(jl),", pawtfl(jl,nls+jk)=",pawtfl(jl,nls+jk),", pawtfl0(jl,nls+jk)=",pawtfl0(jl,nls+jk)
        if(GDCHK2) print *,"ps: 830 .6, rhosn=",rhosn,", csn=",csn,", rhoh2o=",rhoh2o,", clw=",clw,", hsn=",hsn,", stbo=",stbo
        if(GDCHK2) print *,"ps: 830 .6, RHS(4+jk)=",RHS(4+jk) 
#endif
        ENDDO
!       soil layer
        jk=nle+1-nls
!       
        RHS(4+jk)= wtm(nls+jk)&
&         +(1.-beta)*X(4+jk)*(wtm(nls+jk-1)-wtm(nls+jk))&
&         +zdtime*(zsoflw*FFN(zlk(nls+jk)-pwlvl(jl))+hspg)&
&          /(rhogcg*SQRT(xkg/omegas))  
        AA(4+jk,1)=-beta*X(4+jk)
        AA(4+jk,3)=0.
        AA(4+jk,2)=1. - AA(4+jk,1) - AA(4+jk,3)
        IF (lsit_lw) THEN
          AA(4+jk,1)=AA(4+jk,1)-4.*stbo*zdtime*wtm(nls+jk-1)**3/(rhogcg*SQRT(xkg/omegas))
          AA(4+jk,2)=AA(4+jk,2)+4.*stbo*zdtime*wtm(nls+jk)**3/(rhogcg*SQRT(xkg/omegas))
        ENDIF
!       last level matrix updated by this routine
        mae=5+nle-nls
      ELSE
!     Soil Mode
!
!       water completely frozen or evaporated
!       soil layer (note nle=-1)
        IF  ((hesn.LE.csncri).AND.(heice.LE.xicri)) THEN
!       soil on top
          RHS(4)= wtm(nle+1)&
&           +zdtime*(-pfluxw2+hspg+pdfluxs(jl)*wtm(nle+1))&
&            /(rhogcg*SQRT(xkg/omegas))  
          AA(4,1)=0.
          AA(4,3)=0.
          AA(4,2)=1.+zdtime*pdfluxs(jl)/(rhogcg*SQRT(xkg/omegas))&
&           -AA(4,1)-AA(4,3)
        ELSEIF (heice.GT.xicri) THEN
!       water completely frozen, with ice above
          X(4) = zdtime/(rhogcg*SQRT(xkg/omegas))&
&           *rhoice*cice*xkice/(pzsi(jl,1)*rhoh2o/rhoice/2.)
          RHS(4)= wtm(nle+1)&
&           +(1.-beta)*X(4)*(tsim(3)-wtm(nle+1))&
&           +zdtime*hspg/(rhogcg*SQRT(xkg/omegas))  
          AA(4,1)=-beta*X(4)
          AA(4,3)=0.
          AA(4,2)=1. - AA(4,1) - AA(4,3)
        ELSEIF (pzsi(jl,0).GT.csncri) THEN
!       snow on top
          X(4) = zdtime/(rhogcg*SQRT(xkg/omegas))&
&           *rhosn*csn*xksn/(0.5*hsn/pseaice(jl))
          RHS(4)= wtm(nle+1)&
&           +(1.-beta)*X(4)*(tsim(1)-wtm(nle+1))&
&           +zdtime*hspg/(rhogcg*SQRT(xkg/omegas))  
          AA(4,1)=-beta*X(4)
          AA(4,3)=0.
          AA(4,2)=1. - AA(4,1) - AA(4,3)
        ENDIF
        mae=4
      ENDIF
!
!-----------------------------------------------------------------
!
!*    9. Update T by Solving a Tri-Diagonal Matrix.
!         Determine Melting Rate in Snow & Ice, Freezing Rate in Water
!         and adjust their water storage (not water equivalent)
!
!-----------------------------------------------------------------
  900 CONTINUE
!      PRINT *, "I am update Temperature."
!
!*    9.1  SOLVE THE TRIDIAGONAL MXTRIX
!
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 900 9.1 before LU, mae=",mae
        if(GDCHK2) print *,"ps: 900 9.1 before LU, AA(:,1)=",AA(:,1)
        if(GDCHK2) print *,"ps: 900 9.1 before LU, AA(:,2)=",AA(:,2)
        if(GDCHK2) print *,"ps: 900 9.1 before LU, AA(:,3)=",AA(:,3)
        if(GDCHK2) print *,"ps: 900 9.1 before LU, RHS=",RHS(:)
        if(GDCHK2) print *,"ps: 900 9.1 before LU, SOL=",SOL(:)
#endif

      CALL LU(0,mae,AA,RHS,SOL)
!!!      CALL LU(mas,mae,AA,RHS,SOL)
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 900 9.1 LU, mae=",mae
        if(GDCHK2) print *,"ps: 900 9.1 LU, AA(:,1)=",AA(:,1)
        if(GDCHK2) print *,"ps: 900 9.1 LU, AA(:,2)=",AA(:,2)
        if(GDCHK2) print *,"ps: 900 9.1 LU, AA(:,3)=",AA(:,3)
        if(GDCHK2) print *,"ps: 900 9.1 LU, RHS=",RHS(:)
        if(GDCHK2) print *,"ps: 900 9.1 LU, SOL=",SOL(:)
#endif
!
      IF (mas.GE.4.AND..NOT.lsoil) THEN
!     9.2 Water on top: CHK Water Skin Temperature pwt(jl,0) = SOL(4)
        IF ( SOL(4).LT.pctfreez2(jl)-tol) THEN
!       freeze
          lwf=.TRUE.
!         set pwt(0) =tmelt
          RHS(4)= pctfreez2(jl)
          AA(4,1)=0.
          AA(4,3)=0.
          AA(4,2)=1. - AA(4,1) - AA(4,3)
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 900 9.2 before LU, mae=",mae
        if(GDCHK2) print *,"ps: 900 9.2 before LU, AA(:,1)=",AA(:,1)
        if(GDCHK2) print *,"ps: 900 9.2 before LU, AA(:,2)=",AA(:,2)
        if(GDCHK2) print *,"ps: 900 9.2 before LU, AA(:,3)=",AA(:,3)
        if(GDCHK2) print *,"ps: 900 9.2 before LU, RHS=",RHS(:)
        if(GDCHK2) print *,"ps: 900 9.2 before LU, SOL=",SOL(:)
#endif

          CALL LU(0,mae,AA,RHS,SOL)
!!!          CALL LU(mas,mae,AA,RHS,SOL)
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 900 9.2 LU, mae=",mae
        if(GDCHK2) print *,"ps: 900 9.2 LU, AA(:,1)=",AA(:,1)
        if(GDCHK2) print *,"ps: 900 9.2 LU, AA(:,2)=",AA(:,2)
        if(GDCHK2) print *,"ps: 900 9.2 LU, AA(:,3)=",AA(:,3)
        if(GDCHK2) print *,"ps: 900 9.2 LU, RHS=",RHS(:)
        if(GDCHK2) print *,"ps: 900 9.2 LU, SOL=",SOL(:)
#endif

!         calculate FCE of ice/water interface (ice melt positive)
!v.76i
          fcew=rhowcw*hew*( -pctfreez2(jl)+wtm(0)&
&             -beta*Y(4)*(pctfreez2(jl)-SOL(5+nls))&
&             -(1.-beta)*Y(4)*(wtm(0)-wtm(nls+1))  )&
&           +zdtime*(  zsoflw* ( FFN(0.)-FFN(-hw(0)) )&
&             -(pfluxw2+zsoflw) )
!
          iderr=11
          CALL lkerr("Water starts to freeze.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
        ELSE
          lwf=.FALSE.
          fcew=0.
        ENDIF
!
      ELSEIF ( ( (mas.GE.1) .OR.&
&        ((mas.EQ.0).AND.(hesn.LE.csncri))  )&
&      .AND.(heice.GT.xicri) )THEN
!
!     9.3 Ice on top: CHK Thick Ice Skin Temperature store in SOL(1)
!
        IF (SOL(2) .GT. (tmelt+tol)) THEN
          lim=.TRUE.
!         set ptsw(jl) =tmelt
          RHS(2)= tmelt
          AA(2,1)=0.
          AA(2,3)=0.
          AA(2,2)=1.
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 900 9.3 before LU, mae=",mae
        if(GDCHK2) print *,"ps: 900 9.3 before LU, AA(:,1)=",AA(:,1)
        if(GDCHK2) print *,"ps: 900 9.3 before LU, AA(:,2)=",AA(:,2)
        if(GDCHK2) print *,"ps: 900 9.3 before LU, AA(:,3)=",AA(:,3)
        if(GDCHK2) print *,"ps: 900 9.3 before LU, RHS=",RHS(:)
        if(GDCHK2) print *,"ps: 900 9.3 before LU, SOL=",SOL(:)
#endif

          CALL LU(0,mae,AA,RHS,SOL)
!!!          CALL LU(mas,mae,AA,RHS,SOL)
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 900 9.3 LU, mae=",mae
        if(GDCHK2) print *,"ps: 900 9.3 LU, AA(:,1)=",AA(:,1)
        if(GDCHK2) print *,"ps: 900 9.3 LU, AA(:,2)=",AA(:,2)
        if(GDCHK2) print *,"ps: 900 9.3 LU, AA(:,3)=",AA(:,3)
        if(GDCHK2) print *,"ps: 900 9.3 LU, RHS=",RHS(:)
        if(GDCHK2) print *,"ps: 900 9.3 LU, SOL=",SOL(:)
#endif

!         calculate FCE of ice surface (ice melt positive)
          fcei=rhoice*cice*heice*&
&           ( -tmelt+tsim(2)-beta*Y(2)*(tmelt-SOL(3))&
&                       -(1.-beta)*Y(2)*(tsim(2)-tsim(3)) )&
&           +zdtime*( -pfluxi2 )+xife
        ELSE
          lim=.FALSE.
          fcei=0.
        ENDIF
      ELSEIF (pzsi(jl,0).GT.csncri) THEN
!
!     9.4 Snow on top: CHK Thick Snow Skin Temperature store in SOL(0)
!
        IF (SOL(0) .GT. (tmelt+tol) ) THEN
          lsm=.TRUE.
!         set ptsw(jl) =tmelt
          RHS(0)= tmelt
          AA(0,1)=0.
          AA(0,3)=0.
          AA(0,2)= 1.
!         recalculate temp profile
          CALL LU(0,mae,AA,RHS,SOL)
!!!          CALL LU(mas,mae,AA,RHS,SOL)
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 900 9.4 LU, mae=",mae
        if(GDCHK2) print *,"ps: 900 9.4 LU, AA(:,1)=",AA(:,1)
        if(GDCHK2) print *,"ps: 900 9.4 LU, AA(:,2)=",AA(:,2)
        if(GDCHK2) print *,"ps: 900 9.4 LU, AA(:,3)=",AA(:,3)
        if(GDCHK2) print *,"ps: 900 9.4 LU, RHS=",RHS(:)
        if(GDCHK2) print *,"ps: 900 9.4 LU, SOL=",SOL(:)
#endif

!         calculate FCE of snow surface (snow melt positive)
          fces=rhosn*csn*hesn*&
&           ( -tmelt+tsim(0)-beta*Y(0)*(tmelt-SOL(1))&
&                       -(1.-beta)*Y(0)*(tsim(0)-tsim(1)) )&
&           +zdtime*( -pfluxi2 )
        ELSE
          lsm=.FALSE.
          fces=0.
        ENDIF
      ENDIF
!
!     9.5 Determine Phase Change Energy in the ice/water interface
!
      IF ((mas.LT.4).AND.(.NOT.lsoil)) THEN
        lwf=.TRUE.
        IF ( ((mas.EQ.2).AND.(heice.LE.xicri)).OR.&
&            ((mas.EQ.0).AND.(heice.LE.xicri).AND.&
&             (hesn.LE.csncri)) ) THEN
!       thin ice layer only or thin ice & thin snow layer
          fcew=rhowcw*hew*( -SOL(4)+wtm(0)&
&             -beta*Y(4)*(SOL(4)-SOL(5+nls))&
&             -(1.-beta)*Y(4)*(wtm(0)-wtm(nls+1))  )&
&           +zdtime*(  zsoflw* ( FFN(0.)-FFN(-hw(0)) )&
&             -(pfluxi2+pfluxw2+zsoflw) )

!!!          fcew=rhowcw*hew*( -tmelt+wtm(0)&
!!!&             -beta*Y(4)*(tmelt-SOL(5+nls))&
!!!&             -(1.-beta)*Y(4)*(wtm(0)-wtm(nls+1))  )&
!!!&           +zdtime*(  zsoflw* ( FFN(0.)-FFN(-hw(0)) )&
!!!&             -(pfluxw2+zsoflw) )
        ELSEIF (heice.GT.xicri) THEN
!       thick ice layer
          fcew=+(                                                            &
&           (1.-beta)*( rhoice*cice*xkice*(tsim(3)-wtm(0))/(0.5*hice/pseaice(jl))     &
&            -rhowcw*pwkh(jl,nls+1)*(wtm(0)-wtm(nls+1))/(z(0)-z(nls+1))      &
&                     )                                                      &
&           +beta*( rhoice*cice*xkice*(SOL(3)-SOL(4))/(0.5*hice/pseaice(jl))             &
&             -rhowcw*pwkh(jl,nls+1)*(SOL(4)-SOL(5+nls))/(z(0)-z(nls+1))     &
&                 )                                                          &
&           )*zdtime                                                         &
&           +zdtime*(  zsoflw* ( FFN(0.)-FFN(-hw(0)) )                    &
&             -(pfluxw2+zsoflw) )
!         note pwt(jl,0)=wtm(0)=pctfreez2(jl).
        ELSE
!       thin ice layer & thick snow layer
          fcew=+(                                                            &
&           (1.-beta)*( rhosn*csn*xksn*(tsim(1)-tsim(2))/(0.5*hsn/pseaice(jl))        &
&            -rhowcw*pwkh(jl,nls+1)*(wtm(0)-wtm(nls+1))/(z(0)-z(nls+1)) )    &
&           +beta*( rhosn*csn*xksn*(SOL(1)-SOL(2))/(0.5*hsn/pseaice(jl))                 &
&            -rhowcw*pwkh(jl,nls+1)*(SOL(4)-SOL(5+nls))/(z(0)-z(nls+1)) )    &
&           )*zdtime                                                         &
&           +zdtime*(  zsoflw* ( FFN(0.)-FFN(-hw(0)) )                    &
&             -(pfluxw2+zsoflw) )
        ENDIF
        IF (lsit_lw) THEN
          fcew=fcew+zdtime*4.*stbo*(-wtm(nls+1)**3*SOL(4)+wtm(nls+2)**3*SOL(5+nls))
        ENDIF
      ENDIF
!
!*    9.6 Restore Temperatures
!
!*    9.6.1 WATER & SOIL LAYER 
!
      IF (.NOT.lsoil) THEN
        pwt(jl,0) = SOL(4)
        pwt(jl,nls+1:nle+1) = SOL(5:nle+5-nls)
!       Assume the missing layer temperatures to be the first
!         available water temperature, i.e., =pwt(nls+1).
        DO jk = 1, nls, 1
          pwt(jl,jk) = pwt(jl,nls+1)
        END DO
#ifdef DEBUG
          IF (GDCHK2) print *, "ps:9.6.1, pwt(",jl,",0)=", pwt(jl,0)
          IF (GDCHK2) print *, "ps:9.6.1, pwt(",jl,",",nls+1,":",nle+1,")=", pwt(jl,nls+1:nle+1)
#endif
      ELSE
        pwt(jl,nle+1)=SOL(4)
#ifdef DEBUG
          IF (GDCHK2) print *, "ps:9.6.1 else, pwt(",jl,",",nle+1,")=", pwt(jl,nle+1)
#endif
!       SOL(4) is for soil while bare soil or water body completely frozen/
!       evaporated
        DO jk = 0, nle,1
          pwt(jl,jk) = pwt(jl,nle+1)
        END DO
!       Set temp of the missing water to be the temp underneath, soil
!       temp at this case.
      ENDIF
!
!     9.6.2 Ice, Snow & Skin Temp
!!!      IF (mas.LE.2) THEN
!!!        ptsnic(jl,2) = SOL(2)
!!!        ptsnic(jl,3) = SOL(3)
!!!        IF (mas.EQ.0) THEN
!!!          ptsnic(jl,0) = SOL(0)
!!!          ptsnic(jl,1) = SOL(1)
!!!        ENDIF
!!!      ENDIF
      ptsnic(jl,0:3) = SOL(0:3)
!
!v.76i
       gzero=rhoh2o*clw*hew*(pwt(jl,0)-wtm(0))/zdtime-&
&                      zsoflw*( FFN(0.)-FFN(-hw(0)) )+&
&                      rhoh2o*clw*pwkh(jl,nls+1)*(wtm(0)-wtm(1))/(z(0)-z(nls+1))
!v1
!       gzero=-zsoflw*( FFN(0.)-FFN(-hw(0)) )+&
!&                      rhoh2o*clw*pwkh(jl,nls+1)*(wtm(0)-wtm(1))/(z(0)-z(nls+1))
!!!        g0=-(pfluxw2+zsoflw)
!
!-----------------------------------------------------------------
!
!*    10. Determine Melting Rate in Snow & Ice, Freezing Rate in Water
!         and adjust their water storage (not water equivalent)
!         (Note melted water won't refreeze untill next time step)
!
!-----------------------------------------------------------------
 1000 CONTINUE
!      PRINT *, ", I am in 10). PAHSE Change Energy ."
!
!*    10.1 Calc Net Freezing Rate of Water
!
      IF (lwf)THEN
!     water freeze/ice melt in the water-ice interface
        lwf=.FALSE.
!
        IF (fcew.LT.0.) THEN
!       freeze
!-----------------------------------------------------------------------
          fcewf=fcew
          DO WHILE ((fcewf.LE.0.).AND.(nle.GE.nls+1))
            IF (nle.EQ.nls+1) THEN
              dz=hw(nle)-wcri
!             remain minimum water layer of thickness wcri
!             to perserve salinity properity
            ELSE
              dz=hw(nls+1)
            ENDIF
            hw1m=hw(nls+1)
            fcewm=fcewf
            sfm=zsf
            pzsim=pzsi(jl,1)
            tmp=ptsnic(jl,3)
!
            hw(nls+1)=hw1m-dz
            fcewf=fcewf+dz*rhoh2o*&
&             (alf+clw*(pwt(jl,nls+1)-pctfreez2(jl)))
            zsf=zsf-dz*salti*wsm(nls+1)
            pzsi(jl,1)=pzsi(jl,1)+dz
!           modify ice mean temp due to adding in new freezed
!           ice with temp = pctfreez2(jl)
            ptsnic(jl,3)=ptsnic(jl,3)+dz*(pctfreez2(jl)-ptsnic(jl,3))/pzsi(jl,1)
            nls=nls+1
          ENDDO
!         roll back one layer to point correct index
          nls=nls-1
          
          IF(fcewf.GT.0.) THEN
            dz=-fcewm/rhoh2o/(alf+clw*(pwt(jl,nls+1)-pctfreez2(jl)))
            hw(nls+1)=hw1m-dz
            zsf=sfm-dz*salti*wsm(nls+1)
            pzsi(jl,1)=pzsim+dz
            ptsnic(jl,3)=tmp+dz*(pctfreez2(jl)-tmp)/pzsi(jl,1)
          ENDIF
          hice=HICEFN(pzsi(jl,1))
          heice=HEFN(hice/4.,xkice,omegas)     
          IF(fcewf.LT.0.) THEN
!         water body freezes completely
            lshlw=.TRUE.
            iderr=4
            CALL lkerr("Water freezes completely.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
            pfluxi2=pfluxi2-(fcewf-fcew)/zdtime
            pfluxw2=0.
            lim=.FALSE.
            lsm=.FALSE.
            irsv=irsv+1
!!!            CONTINUE
            IF (irsv.LT.4) THEN
!           recursive less than 4 times
              GOTO 800
!             prepare lvl3 matrix & recalculate temperature profile
            ELSE
              iderr=5
              CALL lkerr("Recursive more than 4 times.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
!             The model is unable to judge melt or freeze, just arbitraily
!             pick one.
              CONTINUE
            ENDIF
          ENDIF
        ELSE
!       ice melt
          tmp=fcew/rhoh2o/(alf+cice*(tmelt-ptsnic(jl,3)))
!         tmp is the melting amount of pzsi(jl,1)
          IF (tmp.GE.pzsi(jl,1)) THEN
!         ice melts completely
!         causing snow falling into the water
#ifdef DEBUG
       if(GDCHK2) print*,"ps: 10.1.0 pfluxw2=",pfluxw2
#endif
            pfluxw2=pfluxw2+pfluxi2+(  pzsi(jl,0)*&
&                ( alf+csn*(tmelt-tsim(1))+clw*(wtm(nls+1)-tmelt) )&
&             +pzsi(jl,1)*&
&                ( alf+cice*(tmelt-tsim(3))+clw*(wtm(nls+1)-tmelt) )&
&             +psilw(jl,0)*clw*(wtm(nls+1)-tmelt)&
&             +psilw(jl,1)*clw*(wtm(nls+1)-tmelt)&
&             )*rhoh2o/zdtime
            pfluxi2=0.
#ifdef DEBUG
       if(GDCHK2) print*,"ps: 10.1.1 pfluxi2=",pfluxi2,",prsf(jl)=",prsf(jl),",pssf(jl)=",pssf(jl)
       if(GDCHK2) print*,"ps: 10.1.1 tmp=",tmp,",ptemp2(jl)=",ptemp2(jl),",clw=",clw
       if(GDCHK2) print*,"ps: 10.1.1 alf=",alf,",csn=",csn,",tmelt=",tmelt,",rhoh2o/zdtime=",rhoh2o/zdtime
       if(GDCHK2) print*,"ps: 10.1.1 psilw(jl,0)=",psilw(jl,0),",psilw(jl,1)=",psilw(jl,1)
       if(GDCHK2) print*,"ps: 10.1.1 pzsi(jl,0)=",pzsi(jl,0),",pzsi(jl,1)=",pzsi(jl,1),",cice=",cice
       if(GDCHK2) print*,"ps: 10.1.1 tsim(1)=",tsim(1),",tsim(3)=",tsim(3)
       if(GDCHK2) print*,"ps: 10.1.1 pfluxw2=",pfluxw2
#endif
            zsoflw=zsoflw+zsofli
            zsofli=0.
            hw(nls+1)=hw(nls+1)&
&             +pzsi(jl,0)+pzsi(jl,1)+psilw(jl,0)+psilw(jl,1)
            zsf=zsf+(pzsi(jl,0)+psilw(jl,0))*wsm(nls+1)&
&             +(pzsi(jl,1)+psilw(jl,1))*salti*wsm(nls+1)
!           reset snow/ice progonastic variables
!!!            DO jk=0,1
!!!              pzsi(jl,jk)=0.
!!!              psilw(jl,jk)=0.
!!!              ptsnic(jl,2*jk)=tmelt
!!!              ptsnic(jl,2*jk+1)=tmelt
!!!            ENDDO
            pzsi(jl,:)=0.
            psilw(jl,:)=0.
!!! set to be underneath temp instead
            ptsnic(jl,:)=pwt(jl,0)
!!!
            hice=HICEFN(pzsi(jl,1))
            heice=HEFN(hice/4.,xkice,omegas)
            iderr=7
            CALL lkerr("Ice melts completely due to level 4 forcing.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)

            lim=.FALSE.
            lsm=.FALSE.
            irsv=irsv+1
!!!            CONTINUE
            IF (irsv.LT.4) THEN
!           recursive less than 4 times
              mas=4
              GOTO 830
!             prepare lvl3 matrix & recalculate temperature profile
            ELSE
              iderr=5
              CALL lkerr("Recursive more than 4 times.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
!             The model is unable to judge melt or freeze, just arbitraily
!             pick one.
              CONTINUE
            ENDIF
          ELSE
            pzsi(jl,1)=pzsi(jl,1)-tmp
            hw(nls+1)=hw(nls+1)+tmp
            zsf=zsf+tmp*salti*wsm(nls+1)
            hice=HICEFN(pzsi(jl,1))
            heice=HEFN(hice/4.,xkice,omegas)     
          ENDIF
        ENDIF
      ENDIF
!
!*    10.2 Calc Melting Rate of ice
!
      IF (lim) THEN
!     ice melt on the surface
        lim=.FALSE.
!       melts only
        IF (fcei.LT.0.) THEN
!         numerical error
          iderr=2
          CALL lkerr("Ice melts but phase change energy < 0.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
          fcei=0.
        ENDIF
        tmp=fcei/rhoh2o/(alf+clw*(tmelt-ptsnic(jl,3)))
!       surface ice melting amount
        IF (tmp.GE.pzsi(jl,1)) THEN
!         ice melts completely
!         causing snow falling into the water
          iderr=8
          CALL lkerr("Ice melts completely due to level 3 forcing.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
          pfluxw2=pfluxw2+pfluxi2+(  pzsi(jl,0)*&
&              ( alf+csn*(tmelt-tsim(1))+clw*(wtm(nls+1)-tmelt) )&
&           +pzsi(jl,1)*&
&              ( alf+cice*(tmelt-tsim(3))+clw*(wtm(nls+1)-tmelt) )&
&           +psilw(jl,0)*clw*(wtm(nls+1)-tmelt)&
&           +psilw(jl,1)*clw*(wtm(nls+1)-tmelt)&
&           )*rhoh2o/zdtime
          pfluxi2=0.
          zsoflw=zsoflw+zsofli
          zsofli=0.
          IF (.NOT.lsoil) THEN
            hw(nls+1)=hw(nls+1)+pzsi(jl,0)+pzsi(jl,1)+&
&             psilw(jl,0)+psilw(jl,1)
            zsf=zsf+(pzsi(jl,0)+psilw(jl,0))*wsm(nls+1)&
&             +(pzsi(jl,1)+psilw(jl,1))*salti*wsm(nls+1)
          ELSE
            IF (lwaterlevel) THEN
              pwlvl(jl)=pwlvl(jl)+pzsi(jl,0)+pzsi(jl,1)+&
&               psilw(jl,0)+psilw(jl,1)
            ENDIF
          ENDIF
!         reset snow/ice progonastic variables
!!!          DO jk=0,1
!!!            pzsi(jl,jk)=0.
!!!            psilw(jl,jk)=0.
!!!            ptsnic(jl,2*jk)=tmelt
!!!            ptsnic(jl,2*jk+1)=tmelt
!!!          ENDDO
          pzsi(jl,:)=0.
          psilw(jl,:)=0.
!!! set to be underneath temp instead
          ptsnic(jl,:)=pwt(jl,0)
!!!
          CALL update_snow_ice_property(pzsi(jl,0),pzsi(jl,1),hsn,hesn,hice,heice,pseaice(jl))


          !!! hsn=HSNFN(pzsi(jl,0))
          !!! hesn=HEFN(hsn/4.,xksn,omegas)
          !!! hice=HICEFN(pzsi(jl,1))
          !!! heice=HEFN(hice/4.,xkice,omegas)             
!
          lsm=.FALSE.
          irsv=irsv+1
!!!          CONTINUE
          IF (irsv.LT.4) THEN
!         recursive less than 4 times
            mas=4
            GOTO 830
!           prepare lvl3 matrix & recalculate temperature profile
          ELSE
            iderr=5
            CALL lkerr("Recursive more than 4 times.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
!           The model is unable to judge melt or freeze, just arbitraily
!           pick one.
            CONTINUE
          ENDIF
        ELSE
          pzsi(jl,1)=pzsi(jl,1)-tmp
          psilw(jl,1)=psilw(jl,1)+tmp
          CALL update_snow_ice_property(pzsi(jl,0),pzsi(jl,1),hsn,hesn,hice,heice,pseaice(jl))

          !!! hice=HICEFN(pzsi(jl,1))
          !!! heice=HEFN(hice/4.,xkice,omegas)
        ENDIF
      ENDIF
!
!*    10.3 Calc Melting Rate of snow
!
      IF (lsm) THEN
!     snow melt on the surface
        lsm=.FALSE.
!       melts only
        IF (fces.LT.0.) THEN
!         numerical error
          iderr=3
          CALL lkerr("Snow melts but phase change energy <0.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
          fces=0.
        ENDIF
        tmp= fces/rhoh2o/(alf+clw*(tmelt-ptsnic(jl,1)))
!       tmp is melted amount of SIWE
!       similar error chk routine as 11.1
        IF (tmp.GE.pzsi(jl,0)) THEN
!         snow melts completely
          iderr=12
!          CALL lkerr("Snow melts completely.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
          pfluxi2=pfluxi2+(  pzsi(jl,0)*( alf+csn*(tmelt-tsim(1)) )&
&             )*rhoh2o/zdtime
          psilw(jl,1)=psilw(jl,1)+pzsi(jl,0)+psilw(jl,0)
!         reset snow progonastic variables
          DO jk=0,0
            pzsi(jl,jk)=0.
            psilw(jl,jk)=0.
!!            ptsnic(jl,2*jk)=tmelt
!!            ptsnic(jl,2*jk+1)=tmelt
          ENDDO
!!! set to be underneath temp instead
          ptsnic(jl,0:1)=ptsnic(jl,2)
!
          CALL update_snow_ice_property(pzsi(jl,0),pzsi(jl,1),hsn,hesn,hice,heice,pseaice(jl))

          !!! hsn=HSNFN(pzsi(jl,0))
          !!! hesn=HEFN(hsn/4.,xksn,omegas)
!
          irsv=irsv+1
          IF (irsv.LT.4) THEN
!           recursive less than 4 times
!!!            CONTINUE
            IF (pzsi(jl,1) .GT.xicri) THEN
              mas=2
              GOTO 820
            ELSE
              mas=4
              GOTO 830
            ENDIF
!           prepare lvl2 matrix & recalculate temperature profile
          ELSE
            iderr=5
            CALL lkerr("Recursive more than 4 times.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,pctfreez2(jl),mas,mae,SOL)
!           The model is unable to judge melt or freeze, just arbitraily
!           pick one.
            CONTINUE
          ENDIF
        ELSE
          pzsi(jl,0)=pzsi(jl,0)-tmp
          psilw(jl,0)=psilw(jl,0)+tmp
          CALL update_snow_ice_property(pzsi(jl,0),pzsi(jl,1),hsn,hesn,hice,heice,pseaice(jl))

          !!! hsn=HSNFN(pzsi(jl,0))
          !!! hesn=HEFN(hsn/4.,xksn,omegas)
        ENDIF
      ENDIF
      
! ----------------------------------------------------------------------
      IF ((.NOT.lsoil).AND.(lsit_salt)) THEN
! ----------------------------------------------------------------------
!
!*   11. UPDATE SALINITY DUE TO VERTICAL MIXING
!        BY SOLVING A TRI-DIAGONAL MATRIX:
!
! ----------------------------------------------------------------------
 1100 CONTINUE
!      PRINT *, ", I am in 11) UPDATE Salinity."
        CALL LKDIFKH(hew,X,Y)  
!
!
!* 11.0 Calc. the salinity flux
!
!      
!*   pwlvl(jl)=pbathy(jl)+SUM(hw(nls+1:nle))
!*   Total Fresh Water Inflow = (pwlvl(jl)-wlvlm)/zdtime
!
!!!     ppme2_int(jl)=(pbathy(jl)+SUM(hw(nls+1:nle))-wlvlm)/zdtime
!!! bjt 2010/2/22
!!!     IF(.NOT.lwaterlevel) THEN
!!!       zsf=-ppme2_int(jl)*wsm(nls+1)
!!!     ENDIF
!
     psaltwac(jl)=psaltwac(jl)-zsf
!
!* 11.1 SKIN LAYER
!
        jk = 0
        IF (.FALSE.) THEN
        !!! v7.5 for including runoff salt flux. (2010/5/8) (bjt)
        !!! Since runoff from river usually enters a water column in the top few meters,
        !!! not only in the skin layer
!!!        IF (lssst) THEN
          RHS(4)= wsm(0)-(1.-beta)*Y(4)*(wsm(0)-wsm(nls+1))-zsf/hew
!         Total Fresh Water Inflow = (pwlvl(jl)-wlvlm)/zdtime
          AA(4,1)=0.
          AA(4,3)=-beta*Y(4)
          AA(4,2)= 1. - AA(4,1) - AA(4,3)
        ELSE
          RHS(4)=0.
          AA(4,1)=0.
          AA(4,3)=-1.
          AA(4,2)= 1.
        ENDIF
!
!* 11.2 FIRST LAYER
!
        jk=1
        IF (.FALSE.) THEN
        !!! v7.5 for including runoff salt flux. (2010/5/8) (bjt)
        !!! Since runoff from river usually enters a water column in the top few meters,
        !!! not only in the skin layer
!!!        IF (lssst) THEN
          RHS(4+jk)= wsm(nls+jk)-hew/hw(nls+jk)*wsm(0)+ &
&          (1.-beta)*&
&          ( X(4+jk)*(wsm(0)-wsm(nls+jk))-&
&            Y(4+jk)*(wsm(nls+jk)-wsm(nls+jk+1)) )+ &
&          zdtime*pawsfl(jl,nls+jk)+zdtime*pawsfl0(jl,nls+jk)
          AA(4+jk,1)=-beta*X(4+jk)-hew/hw(nls+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1. +beta*X(4+jk)+beta*Y(4+jk)
        ELSE
          RHS(4+jk)= wsm(nls+jk)-(1.-beta)*Y(4+jk)*(wsm(0)-wsm(nls+jk)) &
&          -zsf/hw(nls+jk)+zdtime*pawsfl(jl,nls+jk)+zdtime*pawsfl0(jl,nls+jk)
          AA(4+jk,1)=0.
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1. - AA(4+jk,1) - AA(4+jk,3)
        ENDIF
!
!       Note AA Matrix from jk=5 to nle is the same as pwt matrix
!
!* 11.3 MIDDLE LAYERS (jk=2,nle)
!
        DO jk=2,nle-nls
!         If nle less than nls+2, do-loop won't do anything.
          RHS(jk+4)= wsm(nls+jk)+(1.-beta)*(X(jk+4)*(wsm(nls+jk-1)-wsm(nls+jk))-&
&           Y(jk+4)*(wsm(nls+jk)-wsm(nls+jk+1)))+ &
&           zdtime*pawsfl(jl,nls+jk)+zdtime*pawsfl0(jl,nls+jk)
          AA(4+jk,1)=-beta*X(4+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1. - AA(4+jk,1) - AA(4+jk,3)
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 11.3, nle=",nle,",nls=",nls 
        if(GDCHK2) print *,"ps: 11.3, wsm(nls+jk)=",wsm(nls+jk),",X(jk+4)=",X(jk+4) 
        if(GDCHK2) print *,"ps: 11.3, wsm(nls+jk-1)=",wsm(nls+jk-1),", Y(jk+4)=",Y(jk+4)
        if(GDCHK2) print *,"ps: 11.3, zdtime=",zdtime,", pawsfl(jl,nls+jk)=",pawsfl(jl,nls+jk)
        if(GDCHK2) print *,"ps: 11.3, zdtime=",zdtime,", pawsfl0(jl,nls+jk)=",pawsfl0(jl,nls+jk)
#endif
        END DO
!
!* 11.4 SOIL LAYER (jk=nle+1=G)
!
!       assume equal to the salinity of previous layer, LVL nle
        RHS(5+nle-nls)=0.
        AA(5+nle-nls,1)=-1.
        AA(5+nle-nls,3)=0.
        AA(5+nle-nls,2)=1.
!
!* 11.5 SOLVE THE TRIDIAGONAL MXTRIX
!
#ifdef DEBUG
        if(GDCHK2) print *,"ps: 11.5, RHS=", RHS(:)
        if(GDCHK2) print *,"ps: 11.5, AA(:,1)=", AA(:,1)
        if(GDCHK2) print *,"ps: 11.5, AA(:,2)=", AA(:,2)
        if(GDCHK2) print *,"ps: 11.5, AA(:,3)=", AA(:,3)
#endif
        CALL LU(4,5+nle-nls,AA,RHS,SOL)
!
!* 11.6 Restore salinity profile
!
#if defined DEBUG
        if(GDCHK2) print *,"ps: 11.6, SOL=", SOL(:)
#endif

        pws(jl,0) = MIN(MAX(SOL(4),0.),SALT_MAX)
        DO jk = 1, nle+1-nls
          pws(jl,nls+jk) = MIN(MAX(SOL(jk+4),0.),SALT_MAX)
        END DO
!       Assume the missing layer salinity to be the first
!       available salinity, i.e., =pws(nls+1).
        DO jk = 1, nls, 1
          pws(jl,jk) = pws(jl,nls+1)
        END DO
!
      ENDIF
!
!* 11.7 Update seawater density for new salinity and temperature
!
     DO jk=0,nle+1
       !! pwrho(jk)=rho_from_theta(pws(jl,jk),pwt(jl,jk)-tmelt,0.)     ! density at surface (0 m depth)
       pwrho(jk)=rho_from_theta(pws(jl,jk),pwt(jl,jk)-tmelt,pwlvl(jl)-z(jk))     ! in situ density
       pwrho1000(jl,jk)=rho_from_theta(pws(jl,jk),pwt(jl,jk)-tmelt,1000.)     ! potential temperature at 1000 m depth
       IF (jk.GE.1) pwrhoh(jk)=rho_from_theta(pws(jl,jk-1),pwt(jl,jk-1)-tmelt,pwlvl(jl)-z(jk))     ! in situ density
     ENDDO

!
!*   12. UPDATE U, V DUE TO VERTICAL MIXING
! ----------------------------------------------------------------------
!      IF (.NOT.(lsoil.OR.(locn.AND.(pocnmask(jl).GT.0.).AND.(ocn_couple_option.NE.11)))) THEN
      IF ( .NOT.lsoil.AND.(.NOT.locn.OR.(pocnmask(jl).LE.0.).OR.(ocn_couple_option.EQ.11)) ) THEN
! ----------------------------------------------------------------------
!
!*   12. UPDATE U, V DUE TO VERTICAL MIXING
!          BY SOLVING A COMPLEX(dp) :: TRI-DIAGONAL MATRIX:
!
! ----------------------------------------------------------------------
!      PRINT *, ", I am in 12) UPDATE U,V."
!
!*   12.1 Update Outflow Velocity & Water Level of Ocean
!
!      OUTFLW(jl)=OUTFLM(jl)+MAX(pwlvl(jl)-z(0),0.)
!      pwlvl(jl)=MIN(pwlvl(jl),z(0))
!
        CALL LKDIFKM(hew,X,Y)
!
!*    this line is moved to Ocean Sub
!
!!!        IF (lsit_ice.AND.mas.LT.4) THEN
        IF (mas.LT.4) THEN
!         set taucx and taucy to zero if snow/ice on top
          taucx(jl)=0.
          taucy(jl)=0.
        ENDIF
        tauc=taucx(jl)*(1.,0.)+taucy(jl)*(0.,1.)
!!
!
!*   12.2 PREPARE CURRENT VECTOR
!
 1220   CONTINUE 
        IF(locn.AND.(ocn_couple_option.EQ.11).AND.(pocnmask(jl).GT.0.)) THEN
        ! ocean grids: taking care in mo_ocean
          fr=0.
          zcor=0.
        ELSE
!          zcor=pcoriol(jl)
!          zcor=MAX(ABS(pcoriol(jl)),zepcor)  ! v9.8b: 20130822
          ! determining coriol, same as vdiff, but strange.
          ! It should have different in different hemisphere
          zcor=SIGN(MAX(ABS(pcoriol(jl)),zepcor),pcoriol(jl))  ! v9.8b: 20130822
          fr=7.E-5                                     ! = friction factor for current (1/s).
                                                       ! It is due to that one-column model neglect horizontal diffusion.
        ENDIF
!        PRINT *, "istep=",istep,"pe=",mpp_pe(),"jrow=",jrow,", I am in UPDATE U,V., 12.2."
        DO jk=0,nle+1
          wumc(jk)=wum(jk)*(1.,0.)+wvm(jk)*(0.,1.)
          pawuflc(jk)=pawufl(jl,jk)*(1.,0.)+pawvfl(jl,jk)*(0.,1.)
        ENDDO
        fr_corc=fr*(1.,0.)+zcor*(0.,1.)
!
!*   12.3 SKIN LAYER
!
!      PRINT *, "istep=",istep,"pe=",mpp_pe(),"jrow=",jrow,", I am in UPDATE U,V., 12.3."
        jk=0
        IF (lssst) THEN
          IF (lv81) THEN
            RHSC(4)= wumc(0)&
&                   -(1.-beta2)*zdtime*fr_corc*wumc(0)&
&                   +zdtime*tauc/rhoh2o/hew&
&                   -(1.-beta)*Y(4)*(wumc(0)-wumc(nls+1))
            AAC(4,1)=(0.,0.)
            AAC(4,3)=-beta*Y(4)*(1.,0.)
            AAC(4,2)= (1.,0.) - AAC(4,1) - AAC(4,3)&
&                   +beta2*zdtime*fr_corc
          ELSE
            RHSC(4)= wumc(0)&
&                   +zdtime*tauc/rhoh2o/hew&
&                   -(1.-beta)*Y(4)*(wumc(0)-wumc(nls+1))
            AAC(4,1)=(0.,0.)
            AAC(4,3)=-beta*Y(4)*(1.,0.)
            AAC(4,2)= (1.,0.) - AAC(4,1) - AAC(4,3)
          ENDIF
        ELSE
          RHSC(4)=(0.,0.)
          AAC(4,1)=(0.,0.)
          AAC(4,3)=(-1.,0.)
          AAC(4,2)= (1.,0.)
        ENDIF
!
!     Note: In the equator, gl_coriol is too small to balance the wind
!     shear. Extra friction might be needed.
!
!*   12.5 FIRST LAYER
!
!      PRINT *, "istep=",istep,"pe=",mpp_pe(),"jrow=",jrow,", I am in UPDATE U,V., 12.4."
        jk=1
        IF (lssst) THEN
          IF(lv81) THEN
            RHSC(4+jk)= wumc(nls+jk)-hew/hw(nls+jk)*wumc(0)&
&             -(1.-beta2)*zdtime*fr_corc*wumc(nls+jk)*0.75&
&             +(1.-beta)*(X(4+jk)*(wumc(0)-wumc(nls+jk))&
&             -Y(4+jk)*(wumc(nls+jk)-wumc(nls+jk+1)))&
&             +zdtime*pawuflc(nls+jk)
            AAC(4+jk,1)=(-beta*X(4+jk)-hew/hw(nls+jk))*(1.,0.)
            AAC(4+jk,3)=(-beta*Y(4+jk))*(1.,0.)
            AAC(4+jk,2)= (1.,0.) + (beta*X(4+jk)+beta*Y(4+jk))*(1.,0.)&
&             +beta2*zdtime*fr_corc*0.75
          ELSE
            RHSC(4+jk)= wumc(nls+jk)-hew/hw(nls+jk)*wumc(0)&
&             -(1.-beta2)*zdtime*fr_corc*wumc(nls+jk) &
&             +(1.-beta)*(X(4+jk)*(wumc(0)-wumc(nls+jk))&
&             -Y(4+jk)*(wumc(nls+jk)-wumc(nls+jk+1)))&
&             +zdtime*pawuflc(nls+jk)
            AAC(4+jk,1)=(-beta*X(4+jk)-hew/hw(nls+jk))*(1.,0.)
            AAC(4+jk,3)=(-beta*Y(4+jk))*(1.,0.)
            AAC(4+jk,2)= (1.,0.) + (beta*X(4+jk)+beta*Y(4+jk))*(1.,0.)&
&             +beta2*zdtime*fr_corc
          ENDIF
        ELSE
          IF(lv81) THEN
            RHSC(4+jk)= wumc(nls+jk)-(1.-beta2)*zdtime*fr_corc*wumc(nls+jk)&
&             -(1.-beta)*Y(4+jk)*(wumc(0)-wumc(nls+jk))&
&             +zdtime*pawuflc(nls+jk)
            AAC(4+jk,1)=(0.,0.)
            AAC(4+jk,3)=(-beta*Y(4+jk))*(1.,0.)
            AAC(4+jk,2)= (1.,0.) - AAC(4+jk,1) - AAC(4+jk,3)&
&             +beta2*zdtime*fr_corc
          ELSE
            RHSC(4+jk)= wumc(nls+jk)-(1.-beta2)*zdtime*fr_corc*wumc(nls+jk)&
&             +zdtime*tauc/rhoh2o/hw(nls+jk)&
&             -(1.-beta)*Y(4+jk)*(wumc(0)-wumc(nls+jk))&
&             +zdtime*pawuflc(nls+jk)
            AAC(4+jk,1)=(0.,0.)
            AAC(4+jk,3)=(-beta*Y(4+jk))*(1.,0.)
            AAC(4+jk,2)= (1.,0.) - AAC(4+jk,1) - AAC(4+jk,3)&
&             +beta2*zdtime*fr_corc
          ENDIF
        ENDIF
!
!       or
! ----------------------------------------------------------------------
!         RHSC(5)= wumc(1)
!     1    -(1.-beta2)*zdtime*fr_corc*wumc(1)
!     2    +zdtime*tauc/rhoh2o/hw(1)
!     3    -(1.-beta)*Y(5)*(wumc(1)-wumc(2))
!
!        AAC(5,1)=(0.,0.)
!        AAC(5,3)=-beta*Y(5)*(1.,0.)
!        AAC(5,2)= (1.,0.) - AAC(5,1) - AAC(5,3)
!     1    +beta2*zdtime*fr_corc
! ----------------------------------------------------------------------
!
!       Note both of the above two equations are identical. The one we
!       choose is still a tridigonal matrix while coupled with ATM PBL
!       directly. i.e. tauc can be implicitly determined by a coupled ATM
!       and ocean momentum matrix by solving a tri-diagonal matrix only. A
!       value of 0.75 to adjust coriolis force due to the physical thickness
!       of skin layer is 0.25 hw(1). Therefore only 75% remained to be
!       counted here.
!
!*   12.6 MIDDLE LAYERS (jk=2,nle)
!
!      PRINT *, "istep=",istep,"pe=",mpp_pe(),"jrow=",jrow,", I am in UPDATE U,V., 12.5."
        DO jk=2,nle-nls
          RHSC(4+jk)= wumc(nls+jk)&
&           -(1.-beta2)*zdtime*fr_corc*wumc(nls+jk)&
&           +(1.-beta)*(X(4+jk)*(wumc(nls+jk-1)-wumc(nls+jk))&
&           -Y(4+jk)*(wumc(nls+jk)-wumc(nls+jk+1)))&
&           +zdtime*pawuflc(nls+jk)
          AAC(4+jk,1)=-beta*X(4+jk)*(1.,0.)
          AAC(4+jk,3)=-beta*Y(4+jk)*(1.,0.)
          AAC(4+jk,2)= (1.,0.) -AAC(4+jk,1) - AAC(4+jk,3)&
&           +beta2*zdtime*fr_corc
!         X,Y is the same as temperature's
        END DO
!
!*   12.7 SOIL LAYER (jk=nle+1=G)
!     
!     zero velocity is assumed.
        jk=nle+1-nls
        RHSC(4+jk)= (0.,0.)
        AAC(4+jk,1)=(0.,0.)
        AAC(4+jk,3)=(0.,0.)
        AAC(4+jk,2)= (1.,0.) - AAC(4+jk,1) - AAC(4+jk,3)
!
!*   12.8 SOLVE THE TRIDIAGONAL MATRIX
!
!      PRINT *, "I am going to LU2."
        CALL LU2(4,5+nle-nls,AAC,RHSC,SOLCMPLX)
!      PRINT *, "I left LU2."
!
!*   12.9 Restore velocity profile
!
!        IF(GDCHK3) THEN
!          WRITE(nerr,*) "I am in 12.9, 1:"
!          WRITE(nerr,*) "pwu=",pwu(jl,:)
!          WRITE(nerr,*) "pwv=",pwv(jl,:)
!        ENDIF 
        pwu(jl,0) = REAL(SOLCMPLX(4))
        pwv(jl,0) = AIMAG(SOLCMPLX(4))
        DO jk = 1, nle+1-nls
          pwu(jl,nls+jk) = REAL(SOLCMPLX(jk+4))
          pwv(jl,nls+jk) = AIMAG(SOLCMPLX(jk+4))
        END DO
!        IF(GDCHK3) THEN
!          WRITE(nerr,*) "I am in 12.9, 2:"
!          WRITE(nerr,*) "pwu=",pwu(jl,:)
!          WRITE(nerr,*) "pwv=",pwv(jl,:)
!        ENDIF 
!         Assume the missing layer current to be the first
!         available current, i.e., =pwu(nls+1), pwv(nls+1).
        DO jk = 1, nls, 1
          pwu(jl,jk) = pwu(jl,nls+1)
          pwv(jl,jk) = pwv(jl,nls+1)
        END DO
! ----------------------------------------------------------------------
      ENDIF
! ----------------------------------------------------------------------
!
!*   13. UPDATE TKE DUE TO VERTICAL MIXING
!            BY SOLVING A TRI-DIAGONAL MATRIX:
!
! ----------------------------------------------------------------------
      IF (.NOT.lsoil) THEN
! ----------------------------------------------------------------------
 1300   CONTINUE
!      PRINT *, "istep=",istep,"pe=",mpp_pe(),"jrow=",jrow,", I am in UPDATE TKE"
        IF (lsteady_TKE) THEN
          CALL calc_steady_tke
        ELSE
          CALL calc_unsteady_tke
        ENDIF
! ----------------------------------------------------------------------
      ENDIF
! ----------------------------------------------------------------------
!*   14. This is needed unless resolution is better than about 10 km; otherwise
!        rotation does not allow proper and full convective adjustment. Should
!        be applied stochastically, as the strong events leading to
! ----------------------------------------------------------------------
      IF (lpenetrative_convection) THEN !!v9.83 (MINOR=3), bjt, 20130824
        CALL penetrative_convection(jl)
      ENDIF
! ----------------------------------------------------------------------
!
!*   15. FINAL: PRINT DIAGNOSTIC VARIABLES
!
! ----------------------------------------------------------------------
!      
!*   15.1 Fluxes
!
      CALL calc_flux
!      
!*   15.2 Water Level
!
!!!      IF(.NOT.lsoil) THEN
!!!        IF (lwaterlevel) THEN      
!!!          pwlvl(jl)=pbathy(jl)+SUM(hw(nls+1:nle))
!!!        ENDIF
!!!      ENDIF
!!!
!!!  do not change the progonastic variable hw, lsoil during the iteration 
!!!      IF (lwaterlevel) THEN
!!!      ! for both water and soil grid      
!!!          CALL pzcord(jl,jrow)                    ! bjt 2010/2/21
!!!      ENDIF
  RETURN
 2100 FORMAT(1X,2(I4,1X),1(I8,1X),(I4,1X),5F8.3,2E9.2)
!
END SUBROUTINE thermocline
! **********************************************************************
! **********************************************************************  
  real FUNCTION calc_tkewave(utauw,depth)
  !  utauw   : friction velocity of the water side (m/s)              I
  !  depth: depth (m), 0 at the surface and plus downward (+, always)
  !  calc_tkewave: tke due to wave (m2/s2)
  !  Mellor and Blumberg (2004)
  !
    IMPLICIT NONE
    real, INTENT(in):: utauw
    real, INTENT(in):: depth
    real,PARAMETER:: ccb=100.      ! constant for wave breaking proposed by Craig and Banner (1994) in Mellor and Blumberg (2004)
                                          ! Mellor and Blumberg (2004), 150 in Stacey (1999)
    real:: lambda
    IF (.NOT.lwave_breaking) THEN
      calc_tkewave=0.
    ELSE
      lambda=1.15925e-5*G/utauw**2.
      calc_tkewave=0.5*(15.8*ccb)**(2./3.)*utauw**2.*EXP(-lambda*depth)
    ENDIF
    IF(GDCHK2) THEN
      print *,"ccb=",ccb,",utauw=",utauw
      print *,"lambda=",lambda,",depth=",depth
    ENDIF
      
  END FUNCTION calc_tkewave
! **********************************************************************
  SUBROUTINE calc_steady_tke
  ! ----------------------------------------------------------------------
  !*   calc TKE at steady state
  ! ----------------------------------------------------------------------          
    IMPLICIT NONE
    INTEGER:: jk,levelm,level
    DO jk=1,nle+1-nls
      levelm=nls+jk-1
      IF (levelm.EQ.nls) levelm=0
      level=nls+jk
      pwtke(jl,level) = ck*pwlmx(jl,level)*pwldisp(jl,level)/ce*                            &
     &        (  ( (pwu(jl,levelm)-pwu(jl,level))/(z(levelm)-z(level)) )**2.+               &
     &           ( (pwv(jl,levelm)-pwv(jl,level))/(z(levelm)-z(level)) )**2.+               &
     &           G/rhoh2o/Prw*( pwrhoh(level)-pwrho(level) )/(z(levelm)-z(level))  )
#ifdef DEBUG     
      IF (GDCHK2) THEN
        print *, "ck=",ck,"pwlmx=",pwlmx(jl,level),"pwldisp=",pwldisp(jl,level)
        print *, "ce=",ce,"pwu(jl,levelm)=",pwu(jl,levelm),"pwu(jl,level)=",pwu(jl,level)
        print *, "z(levelm)=",z(levelm),"z(level)=",z(level)
        print *, "G=",G,"rhoh2o=",rhoh2o,"Prw=",Prw
        print *, "pwrhoh(level)=",pwrhoh(level),"pwrho(level)=",pwrho(level)
        print *, "pwtke(",jl,",",level,")=", pwtke(jl,level)
        print *, "tkewave=",calc_tkewave(utauw,pwlvl(jl)-zlk(level))
        print *, "utauw=",utauw
        print *, "depth=",pwlvl(jl)-zlk(level)
      ENDIF
#endif
      pwtke(jl,level) = MAX(  pwtke(jl,level), emin )+ calc_tkewave(utauw,pwlvl(jl)-zlk(level))
    END DO
    !     Assume the skin layer and the missing layer TKE to be the first
    !     available TKE, i.e., =pwtke(nls+1).
    DO jk = 0, nls, 1
      pwtke(jl,jk) = pwtke(jl,nls+1)
      IF(GDCHK2) THEN
        print *, "pwtke(",jl,",",level,")=", pwtke(jl,level)
      ENDIF
    END DO
  END SUBROUTINE calc_steady_tke
!***********************************************************************
  SUBROUTINE calc_unsteady_tke
  ! ----------------------------------------------------------------------
  !*   calc TKE at unsteady state
  ! ----------------------------------------------------------------------          
    IMPLICIT NONE
    INTEGER:: jk,levelm,level,levelp
    real, DIMENSION(0:lkvl+5):: X,Y,RHS,SOL
    real, DIMENSION(0:lkvl+5,3):: AA
!
!*   13.1 Skin layer (jk = 0)
!
        jk = 0
        IF (lwave_breaking) THEN
        ! set pwtke(0) to be those in Mellor and Blumberg (2004)
          RHS(4)= calc_tkewave(utauw,0.)
          AA(4,1)=0.
          AA(4,3)=0.
          AA(4,2)=1.              
        ELSE
        ! assume to equal to the value of the first layer
          RHS(4)=0.
          AA(4,1)=0.
          AA(4,3)=-1.
          AA(4,2)= 1.
        ENDIF
#ifdef DEBUG
          IF (GDCHK2) print *, "calc_tkewave=",RHS(4)
#endif        
!
!*   13.2 For first to bottom layers (jk = 1 - nle+1)  
!

        DO jk=1,nle+1-nls
          levelm=nls+jk-1
          IF (levelm.EQ.nls) levelm=0
          level=nls+jk
          levelp=nls+jk+1
!!!          IF (levelm.EQ.nls) THEN
!!!          !! first layer
!!!            X(4+jk)= 0.
!!!            Y(4+jk)= zdtime*pwkm(jl,level)/&
!!!&             ( (z(levelm)-z(level))*(zlk(level)-zlk(levelp)) )
!!!            RHS(4+jk)= wtkem(level)&
!!!&             +(1.-beta)*(                                       &
!!!&                         -Y(4+jk)*(wtkem(level)-wtkem(levelp)) )&
!!!&             +beta2* (&
!!!&             +pwkm(jl,level)*zdtime*((wvm(levelm)-wvm(level))/(z(levelm)-z(level)))**2.   &
!!!&             +pwkh(jl,level)*zdtime*G/rhoh2o*( rhomh(level)-            &
!!!&               rhom(level) )/(z(levelm)-z(level)) )                      &
!!!&             +(1-beta2)*(                                                                    &
!!!&             +pwkm(jl,level)*zdtime*((pwu(jl,levelm)-pwu(jl,level))                          &
!!!&                            /(z(levelm)-z(level)))**2.                                    &
!!!&             +pwkm(jl,level)*zdtime*((pwv(jl,levelm)-pwv(jl,level))                          &
!!!&                            /(z(levelm)-z(level)))**2.                                    &
!!!&             +pwkh(jl,level)*zdtime*G/rhoh2o*( pwrhoh(level)-      &
!!!&               pwrho(level) )/(z(levelm)-z(level))    )             &
!!!&             -0.5_dp*ce*zdtime*wtkem(level)**(1.5_dp)/pwldisp(jl,level)                   &   ! bug, v0.9879
!!!&             +zdtime*pawtkefl(jl,level)
!!!
!!!          ELSE IF (levelp.EQ.(nle+2)) THEN
          IF (levelp.EQ.(nle+2)) THEN
          !! bottom layer
            X(4+jk)= zdtime*pwkm(jl,level)/&
&             ( (z(levelm)-z(level))*(zlk(levelm)-zlk(level)) )
            Y(4+jk)= 0.
            RHS(4+jk)= wtkem(level)&
&             +(1.-beta)*( X(4+jk)*(wtkem(levelm)-wtkem(level))&
&                                                               )&
&             +beta2* (&
&             +pwkm(jl,level)*zdtime*((wum(levelm)-wum(level))/(z(levelm)-z(level)))**2.   &
&             +pwkm(jl,level)*zdtime*((wvm(levelm)-wvm(level))/(z(levelm)-z(level)))**2.   &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( rhomh(level)-            &
&               rhom(level) )/(z(levelm)-z(level)) )                      &
&             +(1-beta2)*(                                                                    &
&             +pwkm(jl,level)*zdtime*((pwu(jl,levelm)-pwu(jl,level))                          &
&                            /(z(levelm)-z(level)))**2.                                    &
&             +pwkm(jl,level)*zdtime*((pwv(jl,levelm)-pwv(jl,level))                          &
&                            /(z(levelm)-z(level)))**2.                                    &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( pwrhoh(level)-      &
&               pwrho(level) )/(z(levelm)-z(level))    )             &
&             -0.5*ce*zdtime*wtkem(level)**(1.5)/pwldisp(jl,level)                   &             ! bug, v0.9879
&             +zdtime*pawtkefl(jl,level)
          ELSE
          !! middle layers
            X(4+jk)= zdtime*pwkm(jl,level)/                                                   &
&             ( (z(levelm)-z(level))*(zlk(levelm)-zlk(level)) )
            Y(4+jk)= zdtime*pwkm(jl,level)/&
&             ( (z(levelm)-z(level))*(zlk(level)-zlk(levelp)) )
            RHS(4+jk)= wtkem(level)&
&             +(1.-beta)*( X(4+jk)*(wtkem(levelm)-wtkem(level))                            &
&                         -Y(4+jk)*(wtkem(level)-wtkem(levelp)) )                             &
&             +beta2* (&
&             +pwkm(jl,level)*zdtime*((wum(levelm)-wum(level))/(z(levelm)-z(level)))**2.   &
&             +pwkm(jl,level)*zdtime*((wvm(levelm)-wvm(level))/(z(levelm)-z(level)))**2.   &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( rhomh(level)-            &
&               rhom(level) )/(z(levelm)-z(level)) )                      &
&             +(1-beta2)*(                                                                    &
&             +pwkm(jl,level)*zdtime*((pwu(jl,levelm)-pwu(jl,level))                          &
&                            /(z(levelm)-z(level)))**2.                                    &
&             +pwkm(jl,level)*zdtime*((pwv(jl,levelm)-pwv(jl,level))                          &
&                            /(z(levelm)-z(level)))**2.                                    &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( pwrhoh(level)-      &
&               pwrho(level) )/(z(levelm)-z(level))    )             &
&             -0.5*ce*zdtime*wtkem(level)**(1.5)/pwldisp(jl,level)                   &               ! bug, v0.9879
&             +zdtime*pawtkefl(jl,level)

!            IF(GDCHK3)  then
!              print *,'jk=',jk,',level=',level,',levelm=',levelm
!              print *,'wtkem(levelm)=',wtkem(levelm),',wtkem(level)=',wtkem(level)  
!              print *,'wtkem(levelp)=',wtkem(levelp),',pwkm(jl,level)=',pwkm(jl,level)    
!              print *,'wum(levelm)=',wum(levelm),',wum(level)=',wum(level)              
!              print *,'wvm(levelm)=',wvm(levelm),',wvm(level)=',wvm(level)              
!              print *,'pwkh(jl,level)=',pwkh(jl,level),',rhomh(level)=',rhomh(level)    
!              print *,'rhom(level)=',rhom(level),',z(levelm)=',z(levelm),',z(level)=',z(level)
!              print *,'pwu(jl,levelm)=',pwu(jl,levelm),',pwu(jl,level)=',pwu(jl,level)
!              print *,'pwv(jl,levelm)=',pwv(jl,levelm),',pwv(jl,level)=',pwv(jl,level)
!              print *,'pwrhoh(level)=',pwrhoh(level),',pwrho(level)=',pwrho(level) 
!              print *,'pawtkefl(jl,level)=',pawtkefl(jl,level),',pwldisp(jl,level)=',pwldisp(jl,level)
!              print *,'zdtime=',zdtime,',ce=',ce,',beta=',beta,',beta2=',beta2
!              print *,'pstmp1=',(1.-beta)*( X(4+jk)*(wtkem(levelm)-wtkem(level))                            &
!&                         -Y(4+jk)*(wtkem(level)-wtkem(levelp)) )
!              print *,'pstmp2=',beta2* (&
!&                         +pwkm(jl,level)*zdtime*((wum(levelm)-wum(level))/(z(levelm)-z(level)))**2.   &
!&                         +pwkm(jl,level)*zdtime*((wvm(levelm)-wvm(level))/(z(levelm)-z(level)))**2.   &
!&                         +pwkh(jl,level)*zdtime*G/rhoh2o*( rhomh(level)-            &
!&                           rhom(level) )/(z(levelm)-z(level)) )
!             print *,'pstmp3=',(1-beta2)*(&
!&                         +pwkm(jl,level)*zdtime*((pwu(jl,levelm)-pwu(jl,level))&
!&                            /(z(levelm)-z(level)))**2.&
!&                         +pwkm(jl,level)*zdtime*((pwv(jl,levelm)-pwv(jl,level))&
!&                            /(z(levelm)-z(level)))**2.&
!&                         +pwkh(jl,level)*zdtime*G/rhoh2o*( pwrhoh(level)-      &
!&                           pwrho(level) )/(z(levelm)-z(level))    )
!              print *,'pstmp4=',0.5*ce*zdtime*wtkem(level)**(1.5)/pwldisp(jl,level)
!              print *,'RHS(4+jk)=',RHS(4+jk)
!            ENDIF
          ENDIF
!
!
!     -- DO NOT ALLOW   zdtime*(TOTAL DAMPING) > pwtke   --
!
          RHS(4+jk)=MAX(RHS(4+jk),0.)
!
          AA(4+jk,1)=-beta*X(4+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1.+1.5*ce*zdtime*SQRT(wtkem(level))/pwldisp(jl,level)&
&          -AA(4+jk,1)- AA(4+jk,3)
        ENDDO
!
!*   13.4 SOLVE THE TRIDIAGONAL MATRIX
!
!        IF(GDCHK3) then 
!          print *,'before LU, RHS(:)=',RHS(:)
!          print *,'AA(:,1)=',AA(:,1)
!          print *,'AA(:,2)=',AA(:,2)
!          print *,'AA(:,3)=',AA(:,3)
!          print *,'SOL(:)=',SOL(:)
!        ENDIF
        CALL LU(4,5+nle-nls,AA,RHS,SOL)
!        IF(GDCHK3) then
!          print *,'after LU, RHS(:)=',RHS(:)
!          print *,'AA(:,1)=',AA(:,1)
!          print *,'AA(:,2)=',AA(:,2)
!          print *,'AA(:,3)=',AA(:,3)
!          print *,'SOL(:)=',SOL(:)
!        ENDIF
!
!*   13.5 Restore TKE profile
!     limit tom emin (=1.0E-6) (GASPAR ET AL., 1990)
!
        pwtke(jl,0) = MAX(SOL(4)+emin, tol)
        DO jk = 1, nle+1-nls
          pwtke(jl,nls+jk) = SOL(jk+4)
#ifdef DEBUG
          IF (GDCHK2) print *, "pwtke(",jl,",",nls+jk,")=", pwtke(jl,nls+jk)
#endif
          pwtke(jl,nls+jk) = MAX(pwtke(jl,nls+jk),emin)
!          IF (GDCHK3) print *, "pwtke(",jl,",",nls+jk,")=", pwtke(jl,nls+jk)
        END DO
        
!        IF(GDCHK3) print *,'nls=',nls,',nle=',nle
!         Assume the skin layer and the missing layer TKE to be the first
!         available TKE, i.e., =pwtke(nls+1).
        DO jk = 0, nls, 1
          pwtke(jl,jk) = pwtke(jl,nls+1)
        END DO
  END SUBROUTINE calc_unsteady_tke        
!***********************************************************************
  SUBROUTINE calc_flux
  ! ----------------------------------------------------------------------
  !*   calc fluxes
  ! ----------------------------------------------------------------------
!*   14.2 pgrndcapc, pgrndhflx_int, current, tsw
!
   IF (hesn.GT.csncri) THEN
!  snow on top
     pgrndhflx_int(jl)=-rhosn*csn*xksn*                               &  
       (  beta*(ptsnic(jl,0)-ptsnic(jl,1))                            &
         +(1.-beta)*(tsim(0)-tsim(1))                              &
         )/(0.5*hsn/pseaice(jl))
   ELSEIF (heice.GT.xicri) THEN
!  ice on top
     pgrndhflx_int(jl)=-rhoice*cice*xkice*                            &  
       (  beta*(ptsnic(jl,2)-ptsnic(jl,3))                            &
         +(1.-beta)*(tsim(2)-tsim(3))                              &
         )/(0.5*hice/pseaice(jl))
   ELSEIF(.NOT.lsoil) THEN
!  water on top
     pgrndhflx_int(jl)=-rhowcw*pwkh(jl,nls+1)*                        &  
       (  beta*(pwt(jl,0)-pwt(jl,nls+1))                              &
         +(1.-beta)*(wtm(0)-wtm(nls+1))                            &
         )/(z(0)-z(nls+1))-psoflw(jl)*FFN(-hw(0))       
   ELSE
!  soil on top
     pgrndhflx_int(jl)=0.
     ! the flux below a soil layer of the infinity thickness is zero
   ENDIF
!
!    14.4 Calc cold content ccm (J/m2): the energy to melt snow + ice
!         energy below liquid water at tmelt
!
   pcc(jl)=pzsi(jl,0)*rhoh2o*( alf+csn*(tmelt-tsim(1)) )          &
     +pzsi(jl,1)*rhoh2o*( alf+cice*(tmelt-tsim(3)) )
!
!*   14.5 Calc net surface heat/fresh water flux into ocean
!
!!!     pfluxiw_int: net surface heat flux into ocean (W/m2, + downward)
!!!     pfluxiw_int(jl)=-pfluxi2-pfluxw2+(cc-ccm)/zdtime
!
!    pfluxiw_int: net surface heat flux from ocean (W/m2, + upward)     
!    ppme2_int: net fresh water into ocean (m/s, + downward)
     pfluxiw_int(jl)=pfluxi2+pfluxw2-(pcc(jl)-ccm)/zdtime

!!!     IF (pzsi(jl,1) .GT. 0.) THEN
!!!       pfluxiw_int(jl)=-pfluxi2-pfluxw2+(cc-ccm)/zdtime
!!!!!!       pfluxiw_int(jl)= rhowcw*hw(nls+1)*(pwt(jl,nls+1)-wtm(nls+1))/zdtime           &
!!!!!!         +rhowcw*pwkh(jl,nls+2)*( beta*(pwt(jl,nls+1)-pwt(jl,nls+2))                  &
!!!!!!            +(1.-beta)*(wtm(nls+1)-wtm(nls+2)) )/(z(nls+1)-z(nls+2))
!!!!!!
!!!!!!    Although the below Eq. is correcnt, it is numerically unstable. The above Eq is also correct.
!!!!!!
!!!!!!    pfluxiw_int(jl)= rhowcw*hew*(pwt(jl,0)-wtm(0))/zdtime          &
!!!!!!      rhowcw*pwkh(jl,nls+1)*( beta*(pctfreez2(jl)-pwt(jl,nls+1))    &
!!!!!!         +(1.-beta)*(wtm(0)-wtm(nls+1)) )/(z(0)-z(nls+1))
!!!!!!
!!!!!!    After second thought, the above Eq. is correct
!!!!!!
!!!!!!    the above Eq. is used for freeze/melt seaice
!!!!!!     Since wt(jl,0) is fixed at pctfreez2, pfluxiw_int(jl)=0.
!!!!!!
!!!!!!       pfluxiw_int(jl)=0.
!!!!!!
!!!     ELSE
!!!       pfluxiw_int(jl)=-pfluxi2-pfluxw2
!!!     ENDIF
!!!     IF (abs(pfluxiw_int(jl)).GT.2000.)  CALL output("fluxw > 2000 W/m2. ",hesn,hew,heice,pfluxwm,wtm,wsm,tsim)            
!!!     ppme2_int(jl)=-zsf/wsm(nls+1)/zdtime
!      
!*   pwlvl(jl)=pbathy(jl)+SUM(hw(nls+1:nle))
!*   Total Fresh Water Inflow = (pwlvl(jl)-wlvlm)/zdtime
!
!!!     ppme2_int(jl)=(pbathy(jl)+SUM(hw(nls+1:nle))-wlvlm)/zdtime

!
!*   14.6 Calc the uppermost bulk layer properties for coupling with 3-D ocean
!
!
   IF (locn) THEN
     ! Not valid for lwaterlevel=.TRUE.
     lstfn=MIN(nls+nfnlvl-1,nle)  ! index of last fine level
!
!
!    14.7 Calc net subsurface flux
!
!     pwsubflux_int: subsurface heat flux (W/m2, + upward)
!     pwsubsal_int: subsurface salinity flux (m*PSU/s, + upward)
!
     pwsubflux_int(jl)=-rhowcw*pwkh(jl,lstfn+1)*           &  
         (  beta*(pwt(jl,lstfn)-pwt(jl,lstfn+1))      &
           +(1.-beta)*(wtm(lstfn)-wtm(lstfn+1))    &
           )/(z(lstfn)-z(lstfn+1))                    &
         -psoflw(jl)*FFN(zlk(nls+lstfn+1)-pwlvl(jl))
     pwsubsal_int(jl)=-pwkh(jl,lstfn+1)*                  &    ! unit: 
         ( beta*(pws(jl,lstfn)-pws(jl,lstfn+1))       &
           +(1.-beta)*(wsm(lstfn)-wsm(lstfn+1))    &
           )/(z(lstfn)-z(lstfn+1)) 
   ENDIF
  END SUBROUTINE calc_flux
  !----------------------------------------------------------  
  SUBROUTINE eddy(hcoolskin,utauw,wtkem,wtm,wsm)
  ! ----------------------------------------------------------------------
   real, DIMENSION(0:lkvl+1), INTENT(in):: wtm,wsm
   real, DIMENSION(0:lkvl+1), INTENT(in):: wtkem
   real, INTENT(out):: utauw,hcoolskin
   real:: XLD    ! downward mixing length   (m)
   real:: XLU    ! upward mixing length (m)
   real:: TKESTR ! =TKE*rho/G in (rho*z) (m2/s2*kg/m3/(m/s2))=(kg/m3*m)
   real:: POT    ! rho*z (potential energy difference) (kg/m3 *m)
   real:: wse,wte,RHOE   ! salinity, potential temp, density at TKE level (kg/m3)
   real:: dz     ! distance between two density levels/thickness for freeze (m)
   INTEGER  :: kk, jk  

   IF (.NOT.lsoil) THEN
   !
   !   1.0 calc friction velocity of water side and thickness of cool skin (hcoolskin)
   !           = lamda * nu / uf
     IF(GDCHK2) THEN
       print *,"before frictionvelocity, taucx=",taucx(jl) &
              ,",taucx=",taucx(jl)
     ENDIF
     utauw=frictionvelocity(taucx(jl),taucy(jl))
     IF (lcool_skin) THEN
       hcoolskin=cool_skin(pwind10w(jl),utauw)
     ELSE
       hcoolskin=0. 
     ENDIF
     IF(GDCHK2) THEN
       print *,"after frictionvelocity, utauw=",utauw
     ENDIF
   !
   !   2.0 For the skin layer and the missing layers
   !
     pwldisp(jl,0:nls)=xldispmin
     pwlmx(jl,0:nls)=xlkmin
     pwkm(jl,0:nls)=xkmmin
     pwkh(jl,0:nls)=xkhmin

     DO jk=0,nle+1
!!!      rhom(jk)=rhofn2(wtm(jk),wsm(jk),pwlvl(jl)-z(jk))
!!!       rhom(jk)=rho_from_theta(wsm(jk),wtm(jk)-tmelt,0.)     
       rhom(jk)=rho_from_theta(wsm(jk),wtm(jk)-tmelt,pwlvl(jl)-z(jk))     
       ! potential water density of one level higher (denoted as "h") at level jk
       IF (jk.GE.1) rhomh(jk)=rho_from_theta(wsm(jk-1),wtm(jk-1)-tmelt,pwlvl(jl)-z(jk))     
       IF(GDCHK2) THEN
         print *,"jk=",jk,",rhom=",rhom(jk),",rhomh=",rhomh(jk)
       ENDIF     
     ENDDO
   !
   !   3.0 For layere nls+1 to nle+1
   !
     DO jk=nls+1,nle+1,1
       TKESTR         = wtkem(jk)*rhoh2o/G
!      density at TKE level
       IF (jk.EQ.(nls+1)) THEN
         wse     = (wsm(0)+wsm(jk))*.5
         wte     = (wtm(0)+wtm(jk))*.5
       ELSE
         wse     = (wsm(jk-1)+wsm(jk))*.5
         wte     = (wtm(jk-1)+wtm(jk))*.5
       ENDIF
       IF (GDCHK2) THEN
         print *,"jk=",jk,",TKESTR=",TKESTR,",wtkem=",wtkem(jk)
         print *,"rhoh2o=",rhoh2o,",G=",G,",wse=",wse,",wte=",wte
       ENDIF
!       Calc LD
       XLD        = 0.
       POT        = 0.
       kk         = jk
       DO WHILE ((POT.LT.TKESTR).AND.(kk.LE.nle+1))
        IF ((kk.LE.nls).AND.(kk.NE.0)) THEN
!         Do nothing & skip the layer!
          kk=kk+1
        ELSE
          IF (kk.EQ.nle+1) THEN
            dz=zlk(kk)-z(nle+1)     ! +, always
          ELSE
            dz=zlk(kk)-zlk(kk+1)    ! +, always
          ENDIF
          RHOE=rho_from_theta(wse,wte-tmelt,pwlvl(jl)-z(kk))     
          XLD=XLD+dz
          POT=POT+dz*(rhom(kk)-RHOE)  ! increase along the depth
          kk=kk+1
        ENDIF
        IF(GDCHK2) THEN
          print *,"kk-1=",kk-1,",RHOE=",RHOE,",XLD=",XLD,",POT=",POT
        ENDIF
       ENDDO 
       IF (POT.GT.TKESTR) THEN
       !! restore back slightlty
         IF(kk.EQ.0) THEN
           XLD  = MAX(XLD-(POT-TKESTR)/&
&            (rhom(0)-RHOE),tol)
         ELSE
           XLD  = MAX(XLD-(POT-TKESTR)/&
&            (rhom(kk-1)-RHOE),tol)
           IF(GDCHK2) print *,"1: XLD=",XLD
         ENDIF
       ELSE
         XLD  = MAX(XLD,tol)
         IF(GDCHK2) print *,"2: XLD=",XLD
       ENDIF
       IF(GDCHK2) print *,"3: XLD=",XLD
!       
       ! add zero-dosplacement for reaching the bottom         
       XLD  = MAX(XLD,d0)
!      
!      Calc LU
       XLU        = 0.
       POT        = 0.
       kk         = jk
       DO WHILE ((POT.LT.TKESTR).AND.(kk.GE.nls+1))
         IF ((kk.LE.nls).AND.(kk.NE.0)) THEN
!          Do nothing & skip the layer!
           kk=kk-1
         ELSE
           RHOE=rho_from_theta(wse,wte-tmelt,pwlvl(jl)-z(kk-1))     
           dz=zlk(kk-1)-zlk(kk)
!           dz=zlk(kk)-zlk(kk+1)
           XLU=XLU+dz
           POT=POT+dz*(RHOE-rhom(kk-1))
           kk=kk-1
         ENDIF
        IF(GDCHK2) THEN
          print *,"kk+1=",kk+1,",RHOE=",RHOE,",XLU=",XLU,",POT=",POT
        ENDIF
       ENDDO
       IF (POT.GT.TKESTR) THEN
       !! restore back slightlty
         XLU  = MAX(XLU-(POT-TKESTR)/&
&            (RHOE-rhom(kk+1)),tol)
         IF(GDCHK2) print *,"1: XLU=",XLU
       ELSE
         XLU  = MAX(XLU,tol)
         IF(GDCHK2) print *,"2: XLU=",XLU
       ENDIF
       IF(GDCHK2) print *,"3: XLU=",XLU
       ! add zero-dosplacement for reaching the top         
       XLU  = MAX(XLU,d0)    ! add zero-dosplacement
       pwldisp(jl,jk)  = MAX(SQRT(XLD*XLU),xldispmin)        ! v.76, v76p
       pwlmx(jl,jk)    = MIN(XLD,XLU)
!      pwlmx(jl,jk)    = 2./(1./XLD+1./XLU) 
      pwlmx(jl,jk)    = MAX(pwlmx(jl,jk),xlkmin)            !v.76p
!      
!      4.0 Calc Momentum Diffusivity
!
!      Assunimg no skin layer for momentm since wind shear can enforce on the side. Note that the surface heat flux is from the top. (2011/8/29 bjt)
!!!       pwkm(jl,jk) = MIN(ck*pwlmx(jl,jk)*SQRT(wtkem(jk)),0.1_dp)*Prw+xkmmin                 !v9.7 or v.77
!!!       pwkm(jl,jk) = ck*pwlmx(jl,jk)*SQRT(wtkem(jk))*Prw+xkmmin                             !v9.8: bjt, 2013/5/2, crash
       pwkm(jl,jk) = MIN(ck*pwlmx(jl,jk)*SQRT(wtkem(jk)),100.)             !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s        
!!!       IF ( (z(0)-zlk(jk)).LE. (2.*hcoolskin) ) THEN                       ! v.91, 2004.10.18, (pwu, 1971)
!!!         ! momentum within viscous layer(1mm)
!!!         pwkm(jl,jk)  = xkmmin
!!!       ELSE
!!!         pwkm(jl,jk) = MIN(ck*pwlmx(jl,jk)*SQRT(wtkem(jk)),0.1_dp)*Prw+xkmmin                 ! v.77
!!!       ENDIF
!!       PRINT *,"pwlmx=",pwlmx(jl,jk),"TKE=",wtkem(jk),"pwkm=",pwkm(jl,jk)
!      A maximum value of 10 m2/s is imposed on KM to prevent
!      numerical error.
!      1.4e-6 is the molecular momentum diffusivity of water
!       (Pond & Pickard, 1983).
!      1.34e-6 is the molecular momentum diffusivity of water (Chia and pwu, 1998; Mellor and Durbin, 1975)
!      1.34e-7 is the molecular heat diffusivity of water (Chia and pwu, 1998; Mellor and Durbin, 1975)
!
!      5.0 Calc Heat Diffusivity
!
       IF ( (z(0)-zlk(jk)).LT.hcoolskin) THEN                              !v.77
         ! heat within conduction sublayer (0.4mm)
         pwkh(jl,jk)   = xkhmin
       ELSE
         IF (.TRUE.) THEN
         !! KE due to molecular diffusion has been considered in TKESTR
!ps           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/Prw,100.)                             !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s
           IF( (z(0)-zlk(jk)).LE. mixing_depth) THEN 
             pwkh(jl,jk)   = MIN(MAX(pwkm(jl,jk)/Prw,xkhmin),100.)                             !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s
           ELSE
             pwkh(jl,jk)   = MIN(pwkm(jl,jk)/Prw,100.)
           ENDIF
         ELSEIF ( (z(0)-zlk(jk)).LE.10.) THEN                              !v.77
         ! heat within conduction sublayer (0.4mm)
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/Prw,100.)+xkhmin                      ! otherwise it will crash
         ELSEIF (.TRUE.) THEN
         !! KE due to molecular diffusion has been considered in TKESTR
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/Prw,100.)                             !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s
         ELSEIF (.TRUE.) THEN
           IF (jk.GE.nle) THEN
             pwkh(jl,jk)   = pwkm(jl,jk)/Prw                                        !v9.8994: last 2 levels remain the same
           ELSE
             pwkh(jl,jk)   = SQRT(pwkm(jl,jk)*pwkm(jl,jk+1))/Prw                    !v9.8994: Artifically using diffusivity one-level below (usually very small)
                                                                                        !         (or geometric mean) to preventing heat cross thermocline
           ENDIF
         ELSEIF (.TRUE.) THEN                                                           !v9.83 (MINOR=3), bjt, 20130824
           IF (jk.EQ.nle+1) THEN
             pwkh(jl,jk)   = pwkm(jl,jk)/Prw+xkhmin                                 !v9.8994: last level remain the same
           ELSE 
             pwkh(jl,jk)   = SQRT(pwkm(jl,jk)*pwkm(jl,jk+1))/Prw+xkhmin             !v9.8994: Artifically using diffusivity one-level below (usually very small)
                                                                                        !         (or geometric mean) to preventing heat cross thermocline
           ENDIF
         ELSEIF (.TRUE.) THEN                                                               !v9.83, v9.98998 (MINOR=3), bjt, 20130824
           pwkh(jl,jk)   = pwkm(jl,jk)/Prw+xkhmin                                   !v9.8: bjt, 2013/5/2, crash
         ELSEIF (.TRUE.) THEN                                                           !v9.83 (MINOR=3), bjt, 20130824
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/Prw,100.)+xkhmin                       !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s
         ELSEIF (lv81) THEN
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/Prw,0.1)+xkhmin                       !v9.7 or v.77
         ELSE
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/Prw,100.)+xkhmin                       !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s
         ENDIF
       ENDIF
     END DO
   ENDIF
  END SUBROUTINE eddy
  ! **********************************************************************
  SUBROUTINE REGRID(  pbathy,     pctfreez2,  pwlvl,                      &
                      pzsi,       psilw,      ptsnic,                     &
                      pwt,        pwu,        pwv,                        &
                      pws,     pwtke,      pwlmx,                      & 
                      pwldisp,    pwkm,       pwkh)
  ! ----------------------------------------------------------------------
!!!    USE mo_sst,            ONLY: csn,rhosn,xksn,cice,rhoice,xkice,xkw,                   &
!!!                                 omegas,wcri,tol,wlvlref,dpthmx
    IMPLICIT NONE
    real, INTENT(in):: pbathy(kbdim),   pctfreez2(kbdim)
  ! 3-d SIT vars: water column
    real, INTENT(in out) :: pwlvl(kbdim),                                   &
         pzsi(kbdim,0:1),   psilw(kbdim,0:1), ptsnic(kbdim,0:3),                &
         pwt(kbdim,0:lkvl+1), pwu(kbdim,0:lkvl+1), pwv(kbdim,0:lkvl+1),         &
         pws(kbdim,0:lkvl+1), pwtke(kbdim,0:lkvl+1), pwlmx(kbdim,0:lkvl+1),  & 
         pwldisp(kbdim,0:lkvl+1), pwkm(kbdim,0:lkvl+1), pwkh(kbdim,0:lkvl+1)  
  
  END SUBROUTINE REGRID
  ! ----------------------------------------------------------------------
  real FUNCTION frictionvelocity(taucx,taucy)
    !  taucx    : u-stress (Pa) over water                              I  
    !  taucy    : v-stress (Pa) over water                              I
    !  frictionvelocity: friction velocity of water side (m/s)          O 
    !   Tau/rho=N/m^2/(kg/m3)=(kg*m/s^2/m^2)/(kg/m3)=m^2/s^2
    !   frictionvelocity=sqrt(abs(taucx/rhoh2o)+abs(taucy/rhoh2o))
    !   or, frictionvelocity=sqrt((ustrw**2 + vstrw**2))
    !
    IMPLICIT NONE
    real, INTENT(in):: taucx,taucy
    frictionvelocity=sqrt(abs(taucx/rhoh2o)+abs(taucy/rhoh2o))
    IF(GDCHK2) THEN
     print *,"taucx=",taucx,",taucy=",taucy,",rhoh2o=",rhoh2o
    ENDIF
  END FUNCTION frictionvelocity
  ! ---------------------------------------------------------------------- 
  real FUNCTION cool_skin(wind10w,utauw)
    ! Description:
    !
    ! Calculates the thickness of cool skin by Saunders (1967) and Artale (2002)
    !
    !  wind10w  : 10m wind speed (m/s) over water                I
    !  utauw    : friction velocity of water side (m/s)          I    
    !  cool_skin: thickness of cool skin (m)                     O
    !
    ! Method:
    !
    !   A. Gamma
    !      = 0.2u + 0.5  ; u <= 7.5m/s
    !      = 1.6u - 10   ; 7.5m/s < u < 10m/s
    !      = 6           ; 10 <= u
    !   Where u is the 10m wind speed
    !
    !   B. lamda
    !       = ( uf * k * C )/(gamma * rhoh2o * cw * h * nu )
    !
    !  Where ( for sea water of 20 degree Celsius at salinity 35 [g/kg] )
    !     uf = frictionl velocity of water (or stress)
    !     k = thermal conductivity = 0.596 [W/m K]
    !     C = 86400 s in a day
    !     rhoh2o = density of water = 1024.75 [kg/m**3]
    !     cw = specific heat of water = 3993  [kj/kg] 
    !     h = reference depth = 10m
    !     nu = kinematic viscosity = 1.05*10**-6 [m**2/s]
    !
    !   C. hcoolskin  (Cool skin thickness)
    !       = lamda * nu / uf
    !
    !   D. Temperature difference across cool skin
    !       = Qn * hcoolskin / k
    !   Where 
    !       Qn = net surface heat flux [W/m**2]
    !
    ! *skinsst* is called from *physc*.
    !
    ! Authors:
    !
    ! N. Keenlyside & Chia-Ying Tu, IFM-GEOMAR, June 2006, original source
  
    IMPLICIT NONE
    real, INTENT(in):: wind10w,utauw
    
    ! Local Physical Parameters
    real, PARAMETER:: k = 0.596      ! thermal conductivity  [W/m K]
    INTEGER, PARAMETER  :: C = 86400      ! s in a day
  !!!  real, PARAMETER:: rhoh2o = 1024.75 ! density of water [kg/m**3]
  !!!  real, PARAMETER:: cw = 3993      ! specific heat of water [kj/kg] 
    real, PARAMETER:: refh = 10       
    ! reference depth [m] in Artale's formula of lamda_d (Artale, 2002)
    real, PARAMETER:: nu = 1.05E-6   ! kinematic viscosity of water [m**2/s]
  !  real:: Qb                       !CYTu, v91
  !  real:: gamma            !CYTu, v91
  !  real:: lamda_d          ! nodimensional constant in determing cool skin depth (Sauners, 1967)
  !                                                       ! lamda_d = 5.8 (Jin sitwu, 1971)
  !  real:: lamda_h          
  !  real:: nuw !/1.14E-6/           ! kinematic viscosity of seawater (m2 s-1) (Pauson and Simpson, 1981)
  !
    INTEGER, PARAMETER::  lsaunders=1     ! Artale et al.,2002 scheme
  
    ! Local variables
    real:: gamma, lamda, frictionvelocity, Qnet, Tdiff
    !
    !   A. Gamma
    !      = 0.2u + 0.5  ; u <= 7.5m/s
    !      = 1.6u - 10   ; 7.5m/s < u < 10m/s
    !      = 6           ; 10 <= u
    !   Where u is the 10m wind speed
    IF (lsaunders .EQ. 1) THEN !   Saunders Constant (Artale et al.,2002)
       IF ( wind10w <= 7.5 ) THEN
          gamma = 0.2 * wind10w  + 0.5
       ELSEIF (( 7.5 < wind10w).and.(wind10w < 10 )) THEN
          gamma = 1.6 * wind10w - 10
       ELSE
          gamma=6
       ENDIF
  !   B. lamda
  !           =( uf * k * C )/(gamma * rhoh2o * cw * h * nu )
  !                tauc/rhoh2o
       lamda=(utauw*0.596*86400)/(gamma*rhoh2o*clw*refh*nu)
  !
  !!!!      wsstrw=wsstra*SQRT(1.26/1020.)        ! density ratio of air/seawater (Paulson and Simpson, 1981)
    ELSE IF (lsaunders .EQ. 2)THEN !      Saunders Constant (pwu,1971)
       lamda=5.8
    ELSE IF (lsaunders .EQ. 3)THEN !      Saunders Constant (Paulson and Simpson,1981)
       lamda=6.5
    ELSE IF (lsaunders .EQ. 4)THEN !      Saunders Constant (pwu,1985)
       IF (wind10w.LE.7.)THEN
          lamda=2.+5.*wind10w/7.
       ELSE
          lamda=7.
       ENDIF
  !  ELSE IF (lsaunders .EQ. 5)THEN !     Saunders Constant (Fairall,1996)
  !     Qb=(-Rld+LE+H+Rlu)+((0.026*4.19E-3)/(alv*3.E-4))*H 
  !     lamda=6.*(1.+((Qb*16.*g*3.E-4*rhowcw*(nuw**3.))/((wsstrw**4.)*(0.59**2.)))**(3./4.))**(-1./3.)
    ENDIF
  !
  ! Thickness of Thermal Sublayer (Saunders, 1967)
  !
  !   C. hcoolskin  (Cool skin thickness)
  !           = lamda * nu / uf
    cool_skin=lamda*nu/utauw
  !  !   D. Temperature difference across cool skin
  !  !           = Qn * hcoolskin / k
  !  Tdiff=-pfluxw2*hcoolskin/0.596
  !  !
  !  ptsw(jl)=pobswtb(jl)+Tdiff
    RETURN
  END FUNCTION cool_skin
  ! **********************************************************************
  SUBROUTINE LKDIFKH(hew,X,Y)
  !
  ! Calculate X, Y coefficient matrix for tmeperature and salinity with kh
  !
  ! Output: X,Y,hew
  !
  !*    0. Locate Space
  !
    IMPLICIT NONE
  !
    real, INTENT(out):: hew    ! effective skin thickness of water (m)
    real, DIMENSION(0:lkvl+5), INTENT(out):: X,Y
    
    INTEGER :: jk,levelm,level,levelp
    
  ! For the skin layer 
  
    xkhskin=SQRT(pwkh(jl,0)*pwkh(jl,nls+1))         ! calc the heat diffusivity for skin layer
    hew=HEFN(hw(0),xkhskin,omegas)
  !
  ! For skin layer, middle layers and bottom soil layers
  !
    DO jk=0,nle+1-nls
      levelm=nls+jk-1
      IF (levelm.EQ.nls) levelm=0
      level=nls+jk
      IF (level.EQ.nls) level=0
      levelp=nls+jk+1
      
      IF (level.EQ.0) THEN
      !! Skin layer
        X(4+jk)=0.
        Y(4) = zdtime/hew*pwkh(jl,levelp)/(z(level)-z(levelp))  
      ELSE IF (level.EQ.nle+1) THEN
      !! Soil layer
        X(4+jk)= zdtime/(rhogcg*SQRT(xkg/omegas))*  &
          rhoh2o*clw*pwkh(jl,level)/(z(levelm)-z(level))
        Y(4+jk)=0.
      ELSE
      !! Middle water layers
        X(4+jk)= zdtime/hw(level)*pwkh(jl,level)/(z(levelm)-z(level))
        Y(4+jk)= zdtime/hw(level)*pwkh(jl,levelp)/(z(level)-z(levelp))
      ENDIF
    ENDDO
  END SUBROUTINE LKDIFKH
  ! ----------------------------------------------------------------------
  SUBROUTINE LKDIFKM(hew,X,Y)
  !
  ! Calculate X, Y coefficient (dimensionless) matrix for TKE with Km
  !
  ! Output: X,Y,hew
  !
  !*    0. Locate Space
  !
    IMPLICIT NONE
    real, INTENT(out):: hew    ! effective skin thickness of water (m)
    real, DIMENSION(0:lkvl+5), INTENT(out):: X,Y
  !
    INTEGER :: jk,levelm,level,levelp
    real:: xkmskin   ! skin layer momentum diffusivity (m2/s)
    
  ! For the skin layer 
  
    xkmskin=SQRT(xkmmin*pwkm(jl,0))         ! calc the momentum diffusivity for skin layer
    hew=HEFN(hw(0),xkmskin,omegas)
  !
  ! For skin layer, middle layers and bottom soil layers
  !
    DO jk=0,nle-nls
      levelm=nls+jk-1
      IF (levelm.EQ.nls) levelm=0
      level=nls+jk
      IF (level.EQ.nls) level=0
      levelp=nls+jk+1
      
      IF (level.EQ.0) THEN
      !! Skin layer
        X(4+jk)=0.
        Y(4) = zdtime/hew*pwkm(jl,levelp)/(z(level)-z(levelp))  
      ELSE
      !! Middle water layers
        X(4+jk)= zdtime/hw(level)*pwkm(jl,level)/(z(levelm)-z(level))
        Y(4+jk)= zdtime/hw(level)*pwkm(jl,levelp)/(z(level)-z(levelp))
      ENDIF
    ENDDO
  END SUBROUTINE LKDIFKM
  ! ----------------------------------------------------------------------
  SUBROUTINE LU(mas,mae,AA,RHS,SOL)
  ! ----------------------------------------------------------------------
  !  AA : MATRIX OF COEFFICIENTS A  (MN,kk,3)           I
  !       1ST COLUMN = ELEMENTS LEFT TO THE DIAGONAL,
  !       2ND COLUMN = DIAGONAL ELEMENTS,
  !       3RD COLUMN = ELEMENTS RIGHT TO THE DIAGONAL,
  !  RHS: RIGHT HAND SIDE OF THE MATRIX (FORCING)
  !  SOL: SOLUTION OF THE MATRIX
  !  WKK : WORKING ARRAY               (MN,kk,4)
  ! ----------------------------------------------------------------------
  !
!!!    USE mo_kind
    IMPLICIT NONE
    real:: RHS(0:lkvl+5),AA(0:lkvl+5,3),SOL(0:lkvl+5)
    real:: WKK(0:lkvl+5,4)
    INTEGER jk,mas,mae
    DO jk=mas,mae
      WKK(jk,1)=DBLE(AA(jk,1))
      WKK(jk,2)=DBLE(AA(jk,2))
      WKK(jk,3)=DBLE(AA(jk,3))
      WKK(jk,4)=DBLE(RHS(jk))
    ENDDO
  !
  ! ----------------------------------------------------------------------
  !
  !*   SOLVE THE TRIDIAGONAL MXTRIX
  !       
  !   [ 1  ]       [DR   ] [23  ]
  !   [ L1 ] * [ DR  ]=[123 ]
  !   [ L1 ]       [      DR ] [ 123]
  !-----------------------------------------------------------------
  !
    DO jk = mas+1, mae
      WKK(jk,1) = WKK(jk,1) / WKK(jk-1,2)
      WKK(jk,2) = WKK(jk,2) - WKK(jk,1)*WKK(jk-1,3)
      WKK(jk,4) = WKK(jk,4) - WKK(jk,1)*WKK(jk-1,4)
    ENDDO
  !
  !     BACK SUBSTITUTION
  !
    WKK(mae,4) = WKK(mae,4) /WKK(mae,2)
  !
    DO jk = mae-1, mas,-1
      WKK(jk,4) = (WKK(jk,4)-WKK(jk,3)*WKK(jk+1,4)) /WKK(jk,2)
    ENDDO
    DO jk=mas,mae
      SOL(jk)=REAL(WKK(jk,4))
    ENDDO
  END     SUBROUTINE LU
  ! ----------------------------------------------------------------------
  ! ----------------------------------------------------------------------
  SUBROUTINE LU2(mas,mae,AA,RHS,SOL)
  ! ----------------------------------------------------------------------
  !     same as LU, except for complex matrix
  ! ----------------------------------------------------------------------
  !
    IMPLICIT NONE
    INTEGER mas,mae
    COMPLEX :: RHS(0:lkvl+5),AA(0:lkvl+5,3),SOL(0:lkvl+5),WKK(0:lkvl+5,4)
    INTEGER jk
    DO jk=mas,mae
          WKK(jk,1)=AA(jk,1)
      WKK(jk,2)=AA(jk,2)
      WKK(jk,3)=AA(jk,3)
      WKK(jk,4)=RHS(jk)
    ENDDO
  !
  !-----------------------------------------------------------------
  !*   SOLVE THE TRIDIAGONAL MXTRIX
  !       
  !       [ 1      ]       [DR   ] [23  ]
  !    [ L1        ] * [ DR  ]=[123 ]
  !    [  L1 ]     [      DR ] [ 123]
  !-----------------------------------------------------------------
  !
    DO jk = mas+1, mae
      WKK(jk,1) = WKK(jk,1) / WKK(jk-1,2)
      WKK(jk,2) = WKK(jk,2) - WKK(jk,1)*WKK(jk-1,3)
      WKK(jk,4) = WKK(jk,4) - WKK(jk,1)*WKK(jk-1,4)
    ENDDO
  !
  !     BACK SUBSTITUTION
  !
      WKK(mae,4) = WKK(mae,4) /WKK(mae,2)
  !
    DO jk = mae-1, mas,-1
      WKK(jk,4) = (WKK(jk,4)-WKK(jk,3)*WKK(jk+1,4)) /WKK(jk,2)
    ENDDO
    DO jk=mas,mae
      SOL(jk)=WKK(jk,4)
    ENDDO
  END SUBROUTINE LU2
  ! ----------------------------------------------------------------------
  real FUNCTION levitus_t(tb,z)
  !  z=0. m at the surface (+ upward)
  !  tb: bulk sea temperature (K)
  !-----------------------------------------------------------------------
    IMPLICIT NONE
    real, INTENT(in)::  z, tb
    levitus_t   = tmelt+4+(tb-tmelt-4)*EXP(z/100.)  ! set initial profile to be expontential decay to tmelt+4 K
  END FUNCTION levitus_t
  ! ----------------------------------------------------------------------
  real FUNCTION levitus_s(obswsb,z)
  !  z=0. m at the surface (+ upward)
  !  obswsb: bulk sea salinity (PSU)
  !-----------------------------------------------------------------------
    IMPLICIT NONE
    real, INTENT(in)::  z, obswsb
    levitus_s   = obswsb
  END FUNCTION levitus_s
  ! ----------------------------------------------------------------------
  real FUNCTION FFN(z)
  !*** 2 components ***
  !*****************************************************************************
  ! CALC PENETRATION COEFFICIENT OF SOLAR RADIATION PAULSON AND SIMPSON 1977
  ! J. OF PHYSICAL OCEANOGRAPHY, VOL. 7, PP. 952- 956 
  !*****************************************************************************
  ! CONSTANT (JERLOV's (1976) OPTICAL WATER TYPE I)
  !  z=0. M at the surface (+ upward)
  !-----------------------------------------------------------------------
  !  IMPLICIT NONE
  !  real:: z,R,D1,D2
  !  R    = 0.58
  !  D1   = 0.35
  !  D2   = 23.
  !  FFN   = R*EXP(+z/D1)+(1-R)*EXP(+z/D2)
  !END FUNCTION FFN
  !
  !*** 9 components ***
  !*****************************************************************************
  ! THE TEMPERATURE DIFFERENCE ACROSS COOL SKIN : PAULSON AND SIMPSON 1981
  ! J. OF GEOPHYSICAL RESEARCH, VOL. 86, PP. 11,044- 11,054 
  !*****************************************************************************
  ! EVOLUTION OF COOL SKIN AND DIRECT SEA-AIR GAS TRANSFER COEFFICEINT DURING 
  ! DAYTIME : SOLOVIEV AND SCHLUESSEL 1996
  ! BOUNDARY LAYER METEOROLOGY, VOL. 77, PP. 45- 68
  !*****************************************************************************
  ! SOLOVIEV AND SCHLUESSEL (1996) MODIFIED THE PAULSON AND SIMPSON (1981) EQUATION
  ! OF CLEAR WATER FRO VARIOUS TYPES OF WATER DEFINED BY JERLOV (1976)
  !*****************************************************************************
    IMPLICIT NONE
    real:: z ! depth (m), plus upward
    real:: F(9)      ! spectral distribution (ratio)
    real:: D(9)      ! absorption length (m)
    D(1:9)=(/34.795,2.27,3.15E-2,5.48E-3,8.32E-4,1.26E-4,3.13E-4,7.82E-5,1.44E-5/)
    F(1:9)=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/)
  
    IF(wtype .EQ. 10)THEN
          D(1)    = 15.152                ! SOLOVIEV AND SCHLUESSEL (1996,type I)
    ELSE IF(wtype .EQ. 11)THEN
          D(1)    = 13.158                ! SOLOVIEV AND SCHLUESSEL (1996,type IA)
    ELSE IF(wtype .EQ. 12)THEN
          D(1)    = 11.364                ! SOLOVIEV AND SCHLUESSEL (1996,type IB)
    ELSE IF(wtype .EQ. 20)THEN
          D(1)    = 7.576                 ! SOLOVIEV AND SCHLUESSEL (1996,type II)
    ELSE IF(wtype .EQ. 30)THEN
          D(1)    = 2.618                 ! SOLOVIEV AND SCHLUESSEL (1996,type III)
    ELSE IF(wtype .EQ. 1)THEN
          D(1)    = 2.041                 ! SOLOVIEV AND SCHLUESSEL (1996,type 1)
    ELSE IF(wtype .EQ. 3)THEN
          D(1)    = 1.429                 ! SOLOVIEV AND SCHLUESSEL (1996,type 3)
    ELSE IF(wtype .EQ. 5)THEN
          D(1)    = 1.                   ! SOLOVIEV AND SCHLUESSEL (1996,type 5)
    ELSE IF(wtype .EQ. 7)THEN
          D(1)    = 0.917                 ! SOLOVIEV AND SCHLUESSEL (1996,type 7)
    ELSE IF(wtype .EQ. 9)THEN
          D(1)    = 0.625                 ! SOLOVIEV AND SCHLUESSEL (1996,type 9)
    ELSE
          D(1)    = 34.795                ! PAULSON AND SIMPSON (1981)
    ENDIF
    FFN   = SUM(F*EXP(z/D)) 
  !  dFFN   = (F1/D1)+(F2/D2)+(F3/D3)+(F4/D4)+(F5/D5)+(F6/D6) &
  !                       +(F7/D7)+(F8/D8)+(F9/D9)
  
  END FUNCTION FFN
  ! ----------------------------------------------------------------------
  ! ----------------------------------------------------------------------
  real FUNCTION DFFN(z)
  !*** 2 components ***
  !*****************************************************************************
  ! CALC Derivative of PENETRATION COEFFICIENT OF SOLAR RADIATION PAULSON AND SIMPSON 1977
  ! J. OF PHYSICAL OCEANOGRAPHY, VOL. 7, PP. 952- 956 
  !*****************************************************************************
  ! CONSTANT (JERLOV's (1976) OPTICAL WATER TYPE I)
  !  z=0. M at the surface (+ upward)
  !-----------------------------------------------------------------------
  !  IMPLICIT NONE
  !  real:: z,R,D1,D2
  !  R    = 0.58
  !  D1   = 0.35
  !  D2   = 23.
  !  FFN   = R*EXP(+z/D1)+(1-R)*EXP(+z/D2)
  !END FUNCTION FFN
  !
  !*** 9 components ***
  !*****************************************************************************
  ! THE TEMPERATURE DIFFERENCE ACROSS COOL SKIN : PAULSON AND SIMPSON 1981
  ! J. OF GEOPHYSICAL RESEARCH, VOL. 86, PP. 11,044- 11,054 
  !*****************************************************************************
  ! EVOLUTION OF COOL SKIN AND DIRECT SEA-AIR GAS TRANSFER COEFFICEINT DURING 
  ! DAYTIME : SOLOVIEV AND SCHLUESSEL 1996
  ! BOUNDARY LAYER METEOROLOGY, VOL. 77, PP. 45- 68
  !*****************************************************************************
  ! SOLOVIEV AND SCHLUESSEL (1996) MODIFIED THE PAULSON AND SIMPSON (1981) EQUATION
  ! OF CLEAR WATER FRO VARIOUS TYPES OF WATER DEFINED BY JERLOV (1976)
  !*****************************************************************************
    IMPLICIT NONE
    real:: z ! depth (m), plus upward
    real:: F(9)      ! spectral distribution (ratio)
    real:: D(9)      ! absorption length (m)
    D(1:9)=(/34.795,2.27,3.15E-2,5.48E-3,8.32E-4,1.26E-4,3.13E-4,7.82E-5,1.44E-5/)
    F(1:9)=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/)
  
    IF(wtype .EQ. 10)THEN
          D(1)    = 15.152                ! SOLOVIEV AND SCHLUESSEL (1996,type I)
    ELSE IF(wtype .EQ. 11)THEN
          D(1)    = 13.158                ! SOLOVIEV AND SCHLUESSEL (1996,type IA)
    ELSE IF(wtype .EQ. 12)THEN
          D(1)    = 11.364                ! SOLOVIEV AND SCHLUESSEL (1996,type IB)
    ELSE IF(wtype .EQ. 20)THEN
          D(1)    = 7.576                 ! SOLOVIEV AND SCHLUESSEL (1996,type II)
    ELSE IF(wtype .EQ. 30)THEN
          D(1)    = 2.618                 ! SOLOVIEV AND SCHLUESSEL (1996,type III)
    ELSE IF(wtype .EQ. 1)THEN
          D(1)    = 2.041                 ! SOLOVIEV AND SCHLUESSEL (1996,type 1)
    ELSE IF(wtype .EQ. 3)THEN
          D(1)    = 1.429                 ! SOLOVIEV AND SCHLUESSEL (1996,type 3)
    ELSE IF(wtype .EQ. 5)THEN
          D(1)    = 1.                   ! SOLOVIEV AND SCHLUESSEL (1996,type 5)
    ELSE IF(wtype .EQ. 7)THEN
          D(1)    = 0.917                 ! SOLOVIEV AND SCHLUESSEL (1996,type 7)
    ELSE IF(wtype .EQ. 9)THEN
          D(1)    = 0.625                 ! SOLOVIEV AND SCHLUESSEL (1996,type 9)
    ELSE
          D(1)    = 34.795                ! PAULSON AND SIMPSON (1981)
    ENDIF
    DFFN=SUM(F/D*EXP(z/D))
  !  DFFN   = (F1/D1)+(F2/D2)+(F3/D3)+(F4/D4)+(F5/D5)+(F6/D6) &
  !                       +(F7/D7)+(F8/D8)+(F9/D9)
  
  END FUNCTION DFFN
  !!! bjt >> too cold
  !!!  real, PARAMETER:: csiced=0.8     ! threshold sea ice depth, (= 0.8 m)
  !!!                                           ! seaice[siced]=1-Exp[-siced/csiced]
  !!!                                           ! seaice[2.]=0.91795
  !!! 1)
  !!!    ! Assuming sea ice fraction increases with siced
  !!!    ! sea ice fraction = 100% if siced > csiced (=1 m)
  !!!    ! Otherwise increase linearly with siced/csiced
  !!!    IF (pzsi(jl,1) .GT. csiced) THEN
  !!!      pseaice(jl)=1.
  !!!    ELSE
  !!!      pseaice(jl)=psiced(jl)/csiced
  !!!    ENDIF
  !!! 2)
  !!!    pseaice(jl)=1.-EXP(-psiced(jl)/csiced)
  !!!    seaice2[siced_, m_, c_] := N[Erf[c*(siced - m)]/2 - Erf[-c*m]/2]
  !!!    c=2
  !!!    m=1
  !!! cold bias is found for the above Eq.
  ! ----------------------------------------------------------------------
  SUBROUTINE lkerr(err_message,hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,ctfreeze,mas,mae,SOL)
  !   USE mo_memory_g3b,    ONLY: obswt,obsws
    IMPLICIT NONE
    real, INTENT(in):: hesn,heice,hew,fcew,pfluxwm,ctfreeze
    real, DIMENSION(0:3), INTENT(in):: tsim
    real, DIMENSION(0:lkvl+1), INTENT(in):: wtm,wsm
    real, DIMENSION(0:lkvl+5), INTENT(in):: SOL
    INTEGER, INTENT(in)  :: mas,mae
    CHARACTER(len = *) :: err_message
  !  CHARACTER(len = 80) :: err_message(12)={"Water freezes but phase change energy > 0.",         &  !  1
  !                                       "Ice melts but phase change energy < 0.",             &  !  2
  !                                       "Snow melts but phase change energy <0.",             &  !  3
  !                                       "Water freezes completely.",                           &  !  4
  !                                       "Recursive more than 4 times.",                       &  !  5
  !                                       "Depth is less than 0.2 M. Soil only is assumed.",    &  !  6
  !                                       "Ice melts completely due to level 4 forcing.",       &  !  7
  !                                       "Ice melts completely due to level 3 forcing.",       &  !  8
  !                                       "Recursive goto 830 skin water layer setup.",         &  !  9
  !                                       "Recursive goto 840 ice layer setup.",                &  ! 10 
  !                                       "Water starts to freeze.",                             &  ! 11
  !                                       "Snow melts completely." }                               ! 12
#ifdef DEBUG
    IF (GDCHK2) THEN
#else
    IF (lwarning_msg.GE.3) THEN
#endif
       WRITE(nerr,*) "sit_vdiff: ","jl=",jl,"jrow=",jrow,"istep=",istep,                        &
           "sit_step=",i_sit_step,                                                                &
           "lat=",plat(jl),"lon=",plon(jl), err_message,                      &
           "hew=",hew,"hesn=",hesn,                          &
           "rhoh2o=",rhoh2o,"rhosn=",rhosn,"xksn=",xksn,"omegas=",omegas,                         &
           "heice=",heice,"fcew=",fcew,"ctfreeze=",ctfreeze,                     &
           "pfluxw=",pfluxw2,"sn=",pzsi(jl,0),"ice=",pzsi(jl,1),"tsw=",pwt(jl,0),              &
           "pfluxwm=",pfluxwm,"pdfluxs=",pdfluxs(jl),"soflw=",psoflw(jl),                         &
           "evapw=",pevapw(jl),"rsf=",prsf(jl),"ssf=",pssf(jl),             &
           "silw(0)=",psilw(jl,0),"silw(1)=",psilw(jl,1),                       &
           "rainfall heat flux=",( prsf(jl) )*clw*(wtm(nls+1)-ptemp2(jl)),             &
           "snowfall heat flux=",( pssf(jl) )*( alf+csn*(tmelt-ptemp2(jl))+clw*(wtm(nls+1)-tmelt) ), &
           "snow heat flux=",pzsi(jl,0)*rhoh2o/zdtime*( alf+csn*(tmelt-tsim(1))+clw*(wtm(nls+1)-tmelt) ) &
             +psilw(jl,0)*rhoh2o/zdtime*clw*(wtm(nls+1)-tmelt),                                   &
           "ice flux=", pzsi(jl,1)*rhoh2o/zdtime*( alf+cice*(tmelt-tsim(3))+clw*(wtm(nls+1)-tmelt) )     &
             +psilw(jl,1)*rhoh2o/zdtime*clw*(wtm(nls+1)-tmelt),                                   & 
           "temp2=",ptemp2(jl),"nls=",nls,"nle=",nle,                                             &
           "tsim(:)=",tsim(:),"tsi=",ptsi(jl),                                                    &
           "sittsi(:)=",ptsnic(jl,:),                                                             &
           "wtm(:)=",wtm(:),"sitwt(:)=",pwt(jl,:),                                                &
           "mas=",mas,"mae=",mae,"SOL=",SOL
  !!!         "mas=",mas,"mae=",mae,"AA(0)=",AA(0,:),"RHS=",RHS,"SQL=",SOL
  
    ELSE
       RETURN
    ENDIF
        
  !  CALL write_date(current_date,' sit_vdiff Current date: ')
  !  WRITE(nerr,*) "sit_vdiff Pt. (jl,jrow,istep)=(",jl,jrow, istep,")"
  !  IF (iderr.EQ.1) THEN
  !    WRITE(nerr,*) "Water freezes but phase change energy > 0."
  !  ELSEIF (iderr.EQ.2) THEN
  !    WRITE(nerr,*) "Ice melts but phase change energy < 0."
  !  ELSEIF (iderr.EQ.3) THEN
  !    WRITE(nerr,*) "Snow melts but phase change energy <0."
  !  ELSEIF (iderr.EQ.4) THEN
  !    WRITE(nerr,*) "Water freezes completely."
  !  ELSEIF (iderr.EQ.5) THEN
  !    WRITE(nerr,*) "Recursive more than 4 times."
  !  ELSEIF (iderr.EQ.6) THEN
  !    WRITE(nerr,*) "Depth is less than 0.2 M. Soil only is assumed."
  !  ELSEIF (iderr.EQ.7) THEN
  !    WRITE(nerr,*) "Ice melts completely due to level 4 forcing."
  !  ELSEIF (iderr.EQ.8) THEN
  !    WRITE(nerr,*) "Ice melts completely due to level 3 forcing."
  !  ELSEIF (iderr.EQ.9) THEN
  !    WRITE(nerr,*) "Recursive goto 830 skin water layer setup."
  !  ELSEIF (iderr.EQ.10) THEN
  !    WRITE(nerr,*) "Recursive goto 840 ice layer setup."
  !  ELSEIF (iderr.EQ.11) THEN
  !    WRITE(nerr,*) "Water starts to freeze."
  !  ELSEIF (iderr.EQ.12) THEN
  !    WRITE(nerr,*) "Snow melts completely."
  !  ENDIF
  END SUBROUTINE lkerr
  ! ----------------------------------------------------------------------
  !!!SUBROUTINE OUTPUT(hesn,hew,heice,fcew,pfluxwm,wtm,wum,wvm,wsm,wtkem,tsim,mas,mae)
  SUBROUTINE output(err_message,hesn,hew,heice,pfluxwm,wtm,wsm,tsim)
  !
  !-----------------------------------------------------------------
  !
  !*    OUTPUT: WRITE DIAGNOSTIC VARIABLES
  !
  !-----------------------------------------------------------------
    !!! USE mo_mpi,           ONLY: mpp_pe()   
    IMPLICIT NONE
    real, INTENT(in):: hesn,heice,hew,pfluxwm
    real, DIMENSION(0:3), INTENT(in):: tsim
    real, DIMENSION(0:lkvl+1), INTENT(in):: wtm,wsm
    CHARACTER(len = *) :: err_message  
  !!!  real, INTENT(in):: hesn,heice,hew,fcew,pfluxwm
  !!!  real, DIMENSION(0:3), INTENT(in):: tsim
  !!!  real, DIMENSION(0:lkvl+1), INTENT(in):: wtm,wum,wvm,wsm,wtkem
  !!!  INTEGER, INTENT(in)  :: mas,mae
    INTEGER jk
    real:: tmp
    
  !!!  WRITE(nerr,*) "sit_vdiff: ","pe=",mpp_pe(),"jl=,",jl,"jrow=,",jrow,"istep=",istep,   &
  !!!     "lat=",plat(jl),"lon=",plon(jl), "SST overflow"
!!!    CALL write_date(current_date,' sit_vdiff Current date: ')
  !  WRITE(nerr,*) "sit_vdiff Current date: ",current_date
    WRITE(nerr,*) "sit_vdiff: ","it=",istep,"pe=",mpp_pe(),"jl=",jl,"jrow=",jrow,                  &
        "sit_step=",i_sit_step,                                                                &
        "lat=",plat(jl),"lon=",plon(jl), err_message,                      &
        "fluxiw=",pfluxiw_int(jl),                                                              &
        "pfluxw=",pfluxw2,"pfluxwm=",pfluxwm,"pdfluxs=",pdfluxs(jl),"soflw=",psoflw(jl),       &
        "temp2=",ptemp2(jl),"tsw=",pwt(jl,0),"tsi=",ptsi(jl),                                  &
        "sn=",pzsi(jl,0),"ice=",pzsi(jl,1),"hew=",hew,"hesn=",hesn,                            &
        "rhoh2o=",rhoh2o,"rhosn=",rhosn,"xksn=",xksn,"omegas=",omegas,                         &
        "heice=",heice,                                                                        &
        "evapw=",pevapw(jl),"rsf=",prsf(jl),"ssf=",pssf(jl),                                   &
        "silw(0)=",psilw(jl,0),"silw(1)=",psilw(jl,1),                                         &
        "rainfall heat flux=",( prsf(jl) )*clw*(wtm(nls+1)-ptemp2(jl)),             &
        "snowfall heat flux=",( pssf(jl) )*( alf+csn*(tmelt-ptemp2(jl))+clw*(wtm(nls+1)-tmelt) ), &
        "snow heat flux=",pzsi(jl,0)*rhoh2o/zdtime*( alf+csn*(tmelt-tsim(1))+clw*(wtm(nls+1)-tmelt) ) &
          +psilw(jl,0)*rhoh2o/zdtime*clw*(wtm(nls+1)-tmelt),                                   &
        "ice flux=", pzsi(jl,1)*rhoh2o/zdtime*( alf+cice*(tmelt-tsim(3))+clw*(wtm(nls+1)-tmelt) )     &
          +psilw(jl,1)*rhoh2o/zdtime*clw*(wtm(nls+1)-tmelt)
  
    WRITE(nerr,*) "sitmask=",psitmask(jl).EQ.1.
    WRITE(nerr,*) "lclass=",plclass(jl)
    WRITE(nerr,*) "slm=",pslm(jl)
    WRITE(nerr,*) "seaice=",pseaice(jl)
    WRITE(nerr,*) "siced=",psiced(jl)
    WRITE(nerr,*) "lsoil=",lsoil
    WRITE(nerr,*) "temp2=",ptemp2(jl)
    WRITE(nerr,*) "tsi=",ptsi(jl)
    WRITE(nerr,*) "tsw=",ptsw(jl)
    WRITE(nerr,*) "obstsw=",pobswtb(jl)
    WRITE(nerr,*) "tmelt=",pctfreez2(jl)
    WRITE(nerr,*) "fluxw=",pfluxw(jl)
    WRITE(nerr,*) "fluxw2=",pfluxw2
    WRITE(nerr,*) "fluxi=",pfluxi(jl)
    WRITE(nerr,*) "fluxi=",pfluxi2
    WRITE(nerr,*) "dfluxs=",pdfluxs(jl)
    WRITE(nerr,*) "soflw=",psoflw(jl)
    WRITE(nerr,*) "wtfn(10)=",pwtfn(jl,10)
    WRITE(nerr,*) "wsfn(10)=",pwsfn(jl,10)
    WRITE(nerr,*) "wtfn0(10)=",pwtfn0(jl,10)
    WRITE(nerr,*) "wsfn0(10)=",pwsfn0(jl,10)
    WRITE(nerr,*) "runtoc=",pdisch(jl)
    WRITE(nerr,*) "wind10w=",pwind10w(jl)
    WRITE(nerr,*) "wlvl=",pwlvl(jl)
    WRITE(nerr,*) "evapw=",pevapw(jl)
    WRITE(nerr,*) "rsf=",prsf(jl)
    WRITE(nerr,*) "ssf=",pssf(jl)
    WRITE(nerr,*) "nls=",nls,"nle=",nle,                                                       &
        "tsim(:)=",tsim(:),                                                                    &
        "sittsi(:)=",ptsnic(jl,:),                                                             &
        "wtm(:)=",wtm(:),"sitwt(:)=",pwt(jl,:)
  !
    CALL output2
  !    
    WRITE(nerr,*) "Data in previous time step"
    DO jk = 0, nle+1
      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),&
  &     wtm(jk),wum(jk),wvm(jk),wsm(jk),&
  &     wtkem(jk)
    END DO
  !
  !     calc. column mean temp 
    IF (nle.GE.1)THEN
  !     water exists
      tmp=DOT_PRODUCT(pwt(jl,nls+1:nle),MAX(hw(nls+1:nle),0.))/SUM(MAX(hw(nls+1:nle),0.))
    ELSE
  !     soil only
      tmp=pwt(jl,0)
    ENDIF
    WRITE(nerr,*) "depth=",(z(0)-z(nle+1))
    WRITE(nerr,*) "mean column temp=",tmp
  !
  !ps 2301 FORMAT(1X,4(I3,","),1(I9,","),1(I3,","),1(F10.4,","),1(F8.3,","),2(F8.4,","),&
   2301 FORMAT(1X,3(I3,","),2(F8.2,","),2(I3,","),1(F10.0,","),1(F8.3,","),2(F8.4,","),&
  &       1(F8.3,","),2(E9.2,","),10(F8.3,","))
  
  !
    RETURN
  !
  !
  END SUBROUTINE output
  !
  ! ----------------------------------------------------------------------
  !
  SUBROUTINE penetrative_convection(jl)
  !
  !
  !-----------------------------------------------------------------
  ! Parameterize penetrative convection by mixing unstable surface
  ! layer regions directly with deep water layer where density starts
  ! increasing with depth.
  ! This is needed unless resolution is better than about 10 km; otherwise
  ! rotation does not allow proper and full convective adjustment. Should
  ! be applied stochastically, as the strong events leading to
  ! vigorous convection are short lived and infrequent.
  ! This is especially designed for 18 deg water formation in Sargasso Sea.
  !
  !-----------------------------------------------------------------
  ! Corresponding variables:
  !  TIMCOM: SIT
  !  KB(I,J): nle
  !  K=1: K=nls+1
  ! original code:
  !  u40bjt00@alps6:/work/j07tyh00/PRODUCTION_RUNS/GLOBAL_1DEG_Z61_NOCN_HNUDGa/global.F
  !  line 759-774
  !
    IMPLICIT NONE
  ! 0.0 Calling Variables
    INTEGER, INTENT(IN):: jl    ! lonitude and latitude index
  ! 0.1 Local Variables
    real:: eps               ! exchange ratio for each zdtime time_step
    real:: sdif,tdif
    INTEGER  :: jk
  !!!      hw(nls+1)=Z(3)-Z(1)    ! calculated somewhere else
  ! one-day by half exchange time scale may be too fast for NUMERICAL stability
  ! Set one-day to desired vertical exchange time scale (such as six_hour, one_hour, ...)
    eps=(1.-0.5**(zdtime/one_day))
  ! find unstable top layer points
    IF (pwrho(nls+1).GT.pwrho(nls+2)) THEN
      jk=nls+1
      jk=jk+1
  ! find first layer below top layer having density more than top layer
      DO WHILE (pwrho(nls+1).GT.pwrho(jk).AND.jk.LT.nle) 
        jk=jk+1
      END DO
  ! Vertically exchange heat and salt.
  ! This will usually destabilize next layer down, but that is more easily
  ! convectively adjusted by model due to deep layers being thicker than
  ! those near surface.
      sdif=pws(jl,jk)-pws(jl,nls+1)
      pws(jl,nls+1)=pws(jl,nls+1)+eps*sdif
      pws(jl,jk)=pws(jl,jk)-eps*sdif*hw(nls+1)/hw(nls+jk)
      tdif=pwt(jl,jk)-pwt(jl,nls+1)
      pwt(jl,nls+1)=pwt(jl,nls+1)+eps*tdif
      pwt(jl,jk)=pwt(jl,jk)-eps*tdif*hw(nls+1)/hw(nls+jk)
    ENDIF
  END SUBROUTINE penetrative_convection
  ! ----------------------------------------------------------------------

  SUBROUTINE interpolation_godas(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
!
!     input
!      real:: obstsw, ot12, os12, ou12, ov12
!     output
!      real:: obswt(0:lkvl+1),obsws(0:lkvl+1),obswu(0:lkvl+1),obswv(0:lkvl+1)

  !!! USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2, ndw1, ndw2, wgtd1, wgtd2                           
!!!  USE mo_control,       ONLY: lgodas
!!!  USE sit_constants_mod,     ONLY: tmelt
!!!  USE mo_physc2,        ONLY: ctfreez
  !!! USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, mpp_pe() 
  USE mod_sst,           ONLY: nodepth, odepths, ot12, os12, ou12, ov12
  USE mod_eos_ocean,     ONLY: tmaxden
!ps
  USE mod_sst,        ONLY:now1,now2,wgto1,wgto2,obswtbnmw1,obswtbnmw2,obswtbwgt1,obswtbwgt2  
!ps

  IMPLICIT NONE
  LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wu   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wv   ! logical for missing data
    ! .TRUE. no expolation, missing value returned
    ! .FALSE. expolation enforced, non-missing value returned
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  LOGICAL  :: l_upperdata            ! =TRUE, if data of an upper level is available
  INTEGER  :: jk,kkk
  real:: ttt,ttt1,sss,sss1,uuu,uuu1,vvv,vvv1,depth
!ps
!  real:: ttt0,sss0
  real:: t10m,t10m1,pobswt10m,rr,rf_obs
  INTEGER:: k10m,kmixdepth
!ps
  IF (GDCHK2) then
    WRITE(nerr,*) "I am in interpolation_godas2"
    WRITE(nerr,*) "lwarning_msg=",lwarning_msg
  ENDIF
  IF ((.NOT.lgodas).AND.(pobswtb(jl).EQ.xmissing)) RETURN
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!!
!! Note that the vertical index of g3b starts from 1 although the index in the sit_vdiff 
!!   routine starts from 1,
  IF (GDCHK2) then
    WRITE(nerr,*) "I am in interpolation_godas2, line 6441"
    WRITE(nerr,*) "l_no_expolation_wt=",l_no_expolation_wt
    WRITE(nerr,*) "l_no_expolation_ws=",l_no_expolation_ws
    WRITE(nerr,*) "l_no_expolation_wu=",l_no_expolation_wu
    WRITE(nerr,*) "l_no_expolation_wv=",l_no_expolation_wv
    WRITE(nerr,*) "godas2: nmw1=",nmw1,",nmw2=",nmw2,",wgt1=",wgt1,",wgt2=",wgt2 &
                          ,",now1=",now1,",now2=",now2,",wgto1=",wgto1,",wgto2=",wgto2
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "ot1,", "ot2", "os1,", "os2"
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "obswt,", "obsws", "obswu", "obswv"
    DO jk = 0, nle+1
      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
    END DO
  ENDIF

!find mixing_depth
  IF ( lgodas .AND. lmixedlayer .AND.              &
     (mixedlayer12(jl,jrow,now1).NE.xmissing) .AND. &
     (mixedlayer12(jl,jrow,now2).NE.xmissing)) THEN
    mixedlayer(jl)=max(wgto1*mixedlayer12(jl,jrow,now1)+wgto2*mixedlayer12(jl,jrow,now2),nudge_depth_10)
  ENDIF
!find mixing_depth abs(t-t10m)<= (0.8) degree  (A.B. Kara et al., 2000)
  IF (lgodas) then
    k10m=1
    DO WHILE ( (odepths(k10m).LE. 10.) .AND. (k10m.LT.nodepth) )
      k10m=k10m+1
    END DO
    IF (k10m.GT.1) k10m=k10m-1              ! restore back one layer
  ENDIF
  IF ( lgodas .AND. (.not.lmixedlayer) .AND.              &
       (ot12(jl,k10m,jrow,now1).NE.xmissing)      .AND. &
       (ot12(jl,k10m,jrow,now2).NE.xmissing)      .AND. &
       (ot12(jl,k10m+1,jrow,now1).NE.xmissing)    .AND. &
       (ot12(jl,k10m+1,jrow,now2).NE.xmissing)    .AND. &
       (10. .LE.odepths(nodepth)) ) THEN
    t10m=wgto1*ot12(jl,k10m,jrow,now1)+wgto2*ot12(jl,k10m,jrow,now2)
    t10m1=wgto1*ot12(jl,k10m+1,jrow,now1)+wgto2*ot12(jl,k10m+1,jrow,now2)
    pobswt10m=t10m+ &                                    ! pobswt at 10 m
          (10.-odepths(k10m))/(odepths(k10m+1)-odepths(k10m))* (t10m1-t10m)

    kmixdepth=k10m+1
    DO WHILE ( (ot12(jl,kmixdepth,jrow,now1).NE.xmissing) .AND. &
               (ot12(jl,kmixdepth,jrow,now2).NE.xmissing) .AND. &
               (kmixdepth.LT.nodepth) )
      tmixdown=wgto1*ot12(jl,kmixdepth,jrow,now1)+wgto2*ot12(jl,kmixdepth,jrow,now2)
      IF (abs(tmixdown-pobswt10m).LE. 0.8) THEN
        kmixdepth=kmixdepth+1
      ELSE
        exit
      ENDIF
    END DO
    IF(kmixdepth.GT.1) kmixdepth=kmixdepth-1
      tmixup=wgto1*ot12(jl,kmixdepth,jrow,now1)+wgto2*ot12(jl,kmixdepth,jrow,now2)
      mixlayer=((pobswt10m-0.8)-tmixup)/(tmixdown-tmixup)*(odepths(kmixdepth+1)-odepths(kmixdepth))&
                +odepths(kmixdepth)
      mixedlayer(jl)=max(mixlayer,nudge_depth_10)
    IF(GDCHK3) print *,"1. kmixdepth=",kmixdepth,",mixedlayer=",mixedlayer(jl) &
                      ,",obswt10m=",pobswt10m
  ENDIF


  DO jk=0,nle+1
    depth=(pwlvl(jl)-z(jk))
    IF (lgodas) THEN
      !! initialized the water profile according to world ocean altas data (woa) data
      kkk=1
      DO WHILE ( (odepths(kkk).LE.depth) .AND. (kkk.LT.nodepth) ) 
         kkk=kkk+1
      END DO
      IF (kkk.GT.1) kkk=kkk-1              ! restore back one layer
    ENDIF
!
! 1.0 water salinity
!
    l_upperdata=.FALSE.
    sss=pobswsb(jl)
    sss1=pobswsb(jl)
    IF ( lgodas .AND.                                    &
         (os12(jl,kkk,jrow,now1).NE.xmissing)      .AND. &
         (os12(jl,kkk,jrow,now2).NE.xmissing)      .AND. &
         (os12(jl,kkk+1,jrow,now1).NE.xmissing)    .AND. &
         (os12(jl,kkk+1,jrow,now2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
      l_upperdata=.TRUE.
!      sss0=obswtbwgt1*os12(jl,1,jrow,obswtbnmw1)+obswtbwgt2*os12(jl,1,jrow,obswtbnmw2)
      sss=wgto1*os12(jl,kkk,jrow,now1)+wgto2*os12(jl,kkk,jrow,now2)
      sss1=wgto1*os12(jl,kkk+1,jrow,now1)+wgto2*os12(jl,kkk+1,jrow,now2)
      pobsws(jl,jk)=sss+ &
        (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (sss1-sss)
    ELSE
      ! depth > max obs. depth or missing observed data
      ! water salinity
      IF (l_no_expolation_ws) THEN
      ! no extrapolation
        pobsws(jl,jk)=xmissing
      ELSE 
      ! assume to be sss1
        pobsws(jl,jk)=sss1
!!!      ! extrapolation
!!!        pobsws(jl,jk)=sss+ &
!!!           (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (sss1-sss)
      ENDIF
    ENDIF
!
! 1.1 modify pobswsb/freezing temp again according to initial ocean T/S profile.
!     Note that pobswsb/freezing temp were firstly setup in ioinitial.f90
!    
    IF (lsitstart.AND.jk.EQ.0) THEN
      pobswsb(jl)=MERGE(pobsws(jl,jk),pobswsb(jl),pobsws(jl,jk).NE.xmissing)
#if defined (V9897)
        pctfreez2(jl)=tmelt+3.
#elif defined (TMELTS0)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)           
#elif defined (TMELTS1)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1.,tmelt+1.,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS15)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1.5_dp,tmelt+1.5_dp,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS2)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2.,tmelt+2.,pobswsb(jl).NE.xmissing)  ! v9.9003, there is 7600 ice grids, while the observation is 760 ice grid.    
#elif defined (TMELTS25)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2.,tmelt+2.,pobswsb(jl).NE.xmissing)  ! v9.9007.    
#elif defined (TMELTS3)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+3.,tmelt+3.,pobswsb(jl).NE.xmissing)  ! v9.898 (too hot and too salty in S.H.)
#else
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)
#endif
    ENDIF    
!
! 2.0 water temperature
!
!
!
    l_upperdata=.FALSE.
    ttt=pobswtb(jl)
    ttt1=pobswtb(jl)    
    IF(GDCHK2) then
      print *,"lgodas=",lgodas,",ttt=",ttt,",kkk=",kkk          &
             ,",ot12(jl,kkk,jrow,now1)=",ot12(jl,kkk,jrow,now1) &
             ,",ot12(jl,kkk,jrow,now2)=",ot12(jl,kkk,jrow,now2) &
             ,",ot12(jl,kkk+1,jrow,now1)=",ot12(jl,kkk+1,jrow,now1) &
             ,",ot12(jl,kkk+1,jrow,now2)=",ot12(jl,kkk+1,jrow,now2) &
             ,",depth=",depth,",nodepth=",nodepth                   &
             ,",odepths(nodepth)=",odepths(nodepth)

    ENDIF
    IF ( lgodas .AND.                                    &
         (ot12(jl,kkk,jrow,now1).NE.xmissing)      .AND. &
         (ot12(jl,kkk,jrow,now2).NE.xmissing)      .AND. &
         (ot12(jl,kkk+1,jrow,now1).NE.xmissing)    .AND. &
         (ot12(jl,kkk+1,jrow,now2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      l_upperdata=.TRUE.
!      ttt0=obswtbwgt1*ot12(jl,1,jrow,obswtbnmw1)+obswtbwgt2*ot12(jl,1,jrow,obswtbnmw2)
      ttt=wgto1*ot12(jl,kkk,jrow,now1)+wgto2*ot12(jl,kkk,jrow,now2)
      ttt1=wgto1*ot12(jl,kkk+1,jrow,now1)+wgto2*ot12(jl,kkk+1,jrow,now2)

      IF(mixedlayer(jl) .ge. nudge_depth_10 .AND. pobswtb(jl) .NE. xmissing   &
         .AND. pobswt10m .NE. xmissing)then
        IF(depth .LE. nudge_depth_10)THEN
          rf_obs=1. 
        ELSEIF(depth .GT. mixedlayer(jl) .AND. depth .LT. nudge_depth_100) THEN
          rr=(depth-mixedlayer(jl))/(nudge_depth_100-mixedlayer(jl))
          rf_obs=1.-exp(-0.5/rr*exp(0.5/(rr-1.)))
        ELSE
          rf_obs=0.
        ENDIF

        if(jk .le. 0)then
           pobswt(jl,jk)=MAX(pctfreez2(jl), ttt+ &
            (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
          if(lsftobswt) then
            psftobswt(jl,0)=(pobswtb(jl)-pobswt(jl,0))
          endif
          if(GDCHK3)then
            print*,"int_godas: jk=",jk,",pobswtb=",pobswtb(jl),",pobswt_jk=",pobswt(jl,jk)
          endif
        endif
        if(lsftobswt) then
          psftobswt(jl,jk)=psftobswt(jl,0)
        endif

        IF(depth .LE. mixedlayer(jl) ) THEN
          pobswt(jl,jk)=MAX(pctfreez2(jl), psftobswt(jl,jk)+ttt+ &
            (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
        ELSEIF(depth .GT. mixedlayer(jl) .AND. depth .LT. nudge_depth_100) THEN
          rr=(depth-mixedlayer(jl))/(nudge_depth_100-mixedlayer(jl))
          rf_obs=1.-exp(-0.5/rr*exp(0.5/(rr-1.)))
          psftobswt(jl,jk)=rf_obs*psftobswt(jl,jk)
          pobswt(jl,jk)=MAX(pctfreez2(jl), psftobswt(jl,jk)+ttt+ &
            (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
        ELSE
          psftobswt(jl,jk)=0.
          pobswt(jl,jk)=MAX(pctfreez2(jl), ttt+ &
            (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
        ENDIF
      ELSE
        IF ( depth.LE. nudge_depth_10 ) THEN
      ! depth < 10 m, note daily SST is avaiable from satellite
          IF (pobswtb(jl).NE.xmissing) THEN
            pobswt(jl,jk)=MAX(pctfreez2(jl),pobswtb(jl))
          ELSE
           pobswt(jl,jk)=xmissing
         ENDIF
        ELSE
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
        pobswt(jl,jk)=MAX(pctfreez2(jl), ttt+ &
          (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
!ps        pobswt(jl,jk)=MAX(pctfreez2(jl), ttt &
!ps            +(depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) &
!ps            +rf_obs*(pobswtb(jl)-ttt0) )

!          IF (GDCHK3) then
!            print *,'now1=',now1,',now2=',now2,'wgto1=',wgto1,',wgto2=',wgto2
!            print *,'kkk=',kkk,',ttt0=',ttt0,',ttt=',ttt,',ttt1=',ttt1
!            print *,'pobsws(',jl,',',jk,')=',pobsws(jl,jk)
!            print *,'pobswt(',jl,',',jk,')=',pobswt(jl,jk)
!            print *,'pobswtb(',jl,')=',pobswtb(jl)
!            print *,'obswtbwgt1=',obswtbwgt1,'obswtbwgt2=',obswtbwgt2
!            print *,'obswtbnmw1=',obswtbnmw1,'obswtbnmw2=',obswtbnmw2
!          ENDIF
        ENDIF
      ENDIF
      IF (GDCHK2) then
        WRITE(nerr,*) "godas2.1, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
      ENDIF
    ELSE
      ! depth > max obs. depth or on observation
      ! water temperature
      IF (l_no_expolation_wt) THEN
      ! no extrapolation
        pobswt(jl,jk)=xmissing
        IF (GDCHK2) then
          WRITE(nerr,*) "godas2.2, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
        ENDIF
      ELSE 
      ! extrapolation
        IF (l_upperdata) THEN
          pobswt(jl,jk)=MAX(pctfreez2(jl),ttt+ &
             (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
          IF (GDCHK2) then
            WRITE(nerr,*) "godas2.3, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
          ENDIF
          IF ((tmaxden(pobswsb(jl)).GT.ttt1).AND.(tmaxden(pobswsb(jl)).GT.ttt)) THEN
            pobswt(jl,jk)=MIN(tmaxden(pobsws(jl,jk)),pobswt(jl,jk))
            IF (GDCHK2) then
              WRITE(nerr,*) "godas2.4, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
            ENDIF
          ELSEIF ((tmaxden(pobswsb(jl)).LT.ttt1).AND.(tmaxden(pobswsb(jl)).LT.ttt)) THEN
            pobswt(jl,jk)=MAX(tmaxden(pobsws(jl,jk)),pobswt(jl,jk))
            IF (GDCHK2) then
              WRITE(nerr,*) "godas2.5, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
            ENDIF
          ELSE                                      
            pobswt(jl,jk)=tmaxden(pobsws(jl,jk))
            IF (GDCHK2) then
              WRITE(nerr,*) "godas2.6, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
            ENDIF
          ENDIF
        ELSE          
          ! set initial profile to be expontential decay to tmaxden
          pobswt(jl,jk)=wt_maxden(depth,pobswtb(jl),pobswsb(jl),mixing_depth)
          IF (GDCHK2) then
            WRITE(nerr,*) "godas2.7, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
          ENDIF
        ENDIF
        ! range check
        pobswt(jl,jk)=MAX(pctfreez2(jl),pobswt(jl,jk))
        IF (GDCHK2) then
          WRITE(nerr,*) "godas2.8, pobswt(",jl,",",jk,")=",pobswt(jl,jk)
        ENDIF
      ENDIF
    ENDIF
!=================================================
!
! 3.0 water u current
!
    l_upperdata=.FALSE.
    uuu=0.
    uuu1=0.    
    IF ( lgodas .AND.                                    &
         (ou12(jl,kkk,jrow,now1).NE.xmissing)      .AND. &
         (ou12(jl,kkk,jrow,now2).NE.xmissing)      .AND. &
         (ou12(jl,kkk+1,jrow,now1).NE.xmissing)    .AND. &
         (ou12(jl,kkk+1,jrow,now2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      l_upperdata=.TRUE.
      uuu=wgto1*ou12(jl,kkk,jrow,now1)+wgto2*ou12(jl,kkk,jrow,now2)
      uuu1=wgto1*ou12(jl,kkk+1,jrow,now1)+wgto2*ou12(jl,kkk+1,jrow,now2)
      ! linearly interpolation (no extrapolation)
      pobswu(jl,jk)=uuu+ &
        (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (uuu1-uuu)
    ELSE
      IF (l_no_expolation_wu) THEN
        ! no extrapolation
        pobswu(jl,jk)=xmissing 
      ELSE
        pobswu(jl,jk)=0.
!!!     ! extrapolation
!!!       pobswu(jl,jk)=uuu+ &
!!!          (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (uuu1-uuu)
      ENDIF
    ENDIF
!
! 4.0 water v current
!
    l_upperdata=.FALSE.
    vvv=0.
    vvv1=0.    
    IF ( lgodas .AND.                                    &
         (ov12(jl,kkk,jrow,now1).NE.xmissing)      .AND. &
         (ov12(jl,kkk,jrow,now2).NE.xmissing)      .AND. &
         (ov12(jl,kkk+1,jrow,now1).NE.xmissing)    .AND. &
         (ov12(jl,kkk+1,jrow,now2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      l_upperdata=.TRUE.
      vvv=wgto1*ov12(jl,kkk,jrow,now1)+wgto2*ov12(jl,kkk,jrow,now2)
      vvv1=wgto1*ov12(jl,kkk+1,jrow,now1)+wgto2*ov12(jl,kkk+1,jrow,now2)
      ! linearly interpolation (no extrapolation)
      pobswv(jl,jk)=vvv+ &
        (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (vvv1-vvv)
    ELSE
      IF (l_no_expolation_wu) THEN
        ! no extrapolation
        pobswv(jl,jk)=xmissing 
      ELSE
        pobswv(jl,jk)=0.
!!!     ! extrapolation
!!!       pobswv(jl,jk)=vvv+ &
!!!          (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (vvv1-vvv)
      ENDIF
    ENDIF
  ENDDO
!
! 5.0 Adjust pwt for accounting the ice grid for depth <= 10 m and the limit of pctfreez2
!
!  IF (.NOT.lsoil) pobswt(jl,0:nle+1)=MERGE(MAX(pctfreez2(jl),pobswt(jl,0:nle+1)),xmissing,pobswt(jl,0:nle+1).NE.xmissing)
  DO jk=0,nle+1
    IF (pobsseaice(jl).GT.0.) THEN
      depth=(pwlvl(jl)-z(jk))
      IF ((depth.LE.10.).AND.(pobswt(jl,jk).NE.xmissing)) pobswt(jl,jk)=pctfreez2(jl)
    ENDIF
  END DO
  
  IF (GDCHK2) then
    WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "obswt,", "obsws", "obswu", "obswv"
    DO jk = 0, nle+1
      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
    END DO
    WRITE(nerr,*) "I am leaving interpolation_godas"
  ENDIF
 2300 FORMAT(1X,3A4,2A9,2A4,1A11,4A9,2A10,10A9)
 2301 FORMAT(1X,3(I3,","),2(F8.2,","),2(I3,","),1(F10.0,","),10(F8.3,","),2(E10.2,","),&
&       1(F8.3,","),2(E9.2,","),10(F8.3,","))
  END SUBROUTINE interpolation_godas  
! **********************************************************************
  SUBROUTINE interpolation_ocaf(jl,jrow,l_no_expolation_wt,l_no_expolation_ws)
!
!     input
!      real:: wtfn12, wsfn12
!     output
!      real:: awtfl(0:lkvl+1),awsfl(0:lkvl+1)

  USE mod_sst,       ONLY: nmw1, nmw2, wgt1, wgt2                           
  USE mod_sit_control,       ONLY: locaf
!!!  USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe 
  USE mod_sst,           ONLY: nwdepth, wdepths, wtfn12, wsfn12

  
  IMPLICIT NONE
  LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
    ! .TRUE. no expolation, missing value returned
    ! .FALSE. expolation enforced, non-missing value returned
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  INTEGER  :: jk,kkk
  real:: ttt,ttt1,sss,sss1,depth
   
  IF(lwarning_msg.GE.4) then
    WRITE(nerr,*) ", I am in interpolation_ocaf"
  ENDIF
  IF (.NOT.locaf) RETURN
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!!
!! Note that the vertical index of g3b starts from 1 although the index in the sit_vdiff 
!!   routine starts from 1,
  DO jk=0,nle+1
    depth=(pwlvl(jl)-z(jk))
    !! initialized the water profile according to world ocean altas data (woa) data
    kkk=1
    DO WHILE ( (wdepths(kkk).LE.depth) .AND. (kkk.LT.nwdepth) ) 
       kkk=kkk+1
    END DO
    IF (kkk.GT.1) kkk=kkk-1  ! restore back one layer
    IF (  &
         (wtfn12(jl,kkk,jrow,nmw1).NE.xmissing) .AND. &
         (wtfn12(jl,kkk,jrow,nmw2).NE.xmissing) .AND. &
         (wtfn12(jl,kkk+1,jrow,nmw1).NE.xmissing) .AND. &
         (wtfn12(jl,kkk+1,jrow,nmw2).NE.xmissing) .AND. &
         (wsfn12(jl,kkk,jrow,nmw1).NE.xmissing) .AND. &
         (wsfn12(jl,kkk,jrow,nmw2).NE.xmissing) .AND. &
         (wsfn12(jl,kkk+1,jrow,nmw1).NE.xmissing) .AND. &
         (wsfn12(jl,kkk+1,jrow,nmw2).NE.xmissing) ) THEN
      ttt=wgt1*wtfn12(jl,kkk,jrow,nmw1)+wgt2*wtfn12(jl,kkk,jrow,nmw2)
      ttt1=wgt1*wtfn12(jl,kkk+1,jrow,nmw1)+wgt2*wtfn12(jl,kkk+1,jrow,nmw2)
      sss=wgt1*wsfn12(jl,kkk,jrow,nmw1)+wgt2*wsfn12(jl,kkk,jrow,nmw2)
      sss1=wgt1*wsfn12(jl,kkk+1,jrow,nmw1)+wgt2*wsfn12(jl,kkk+1,jrow,nmw2)
      IF ( (depth.LE.wdepths(nwdepth)) ) THEN
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
        pawtfl(jl,jk)=ttt+ &
          (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (ttt1-ttt)
        pawsfl(jl,jk)=sss+ &
          (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (sss1-sss)
      ELSE
      ! depth > max obs. depth
      ! water temperature flux
        IF (l_no_expolation_wt) THEN
        ! no extrapolation
          pawtfl(jl,jk)=0.
        ELSE 
        ! extrapolation
          pawtfl(jl,jk)=ttt+ &
             (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (ttt1-ttt)
        ENDIF
      ! water salinity flux
        IF (l_no_expolation_ws) THEN
        ! no extrapolation
          pawsfl(jl,jk)=0.
        ELSE 
        ! extrapolation
          pawsfl(jl,jk)=sss+ &
             (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (sss1-sss)
        ENDIF
      ENDIF
    ELSE
      IF (l_no_expolation_wt) THEN
        ! no extrapolation
        pawtfl(jl,jk)=0.
      ELSE
        pawtfl(jl,jk)=0.
      ENDIF
      IF (l_no_expolation_ws) THEN
        ! no extrapolation
        pawsfl(jl,jk)=0. 
      ELSE
        pawsfl(jl,jk)=0.
      ENDIF
    ENDIF
  ENDDO
  IF (GDCHK2) then
    DO jk=1,nle
      WRITE(nerr,*) 'pe=',mpp_pe(),',jl=',jl,',z=',z(jk),',pawtfl=',pawtfl(jl,jk),',pawsfl=',pawsfl(jl,jk)
    ENDDO
    WRITE(nerr,*) "I am leaving interpolation_ocaf"
  ENDIF
  END SUBROUTINE interpolation_ocaf  

! **********************************************************************
  SUBROUTINE interpolation_ocaf0(jl,jrow,l_no_expolation_wt,l_no_expolation_ws)
!
!     input
!      real:: wtfn12, wsfn12
!     output
!      real:: awtfl0(0:lkvl+1),awsfl(0:lkvl+1)

  USE mod_sst,       ONLY: nmw1, nmw2, wgt1, wgt2
  USE mod_sit_control,       ONLY: locaf0
!!!  USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe
  USE mod_sst,           ONLY: nwdepth, wdepths, wtfn1st, mask1st


  IMPLICIT NONE
  LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
    ! .TRUE. no expolation, missing value returned
    ! .FALSE. expolation enforced, non-missing value returned
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  INTEGER  :: jk,kkk
  real:: ttt,ttt1,sss,sss1,depth

  IF(lwarning_msg.GE.4) then
    WRITE(nerr,*) ", I am in interpolation_ocaf"
  ENDIF
  IF (.NOT.locaf0) RETURN
!!   routine starts from 1,
!ps  DO jk=0,nle+1
!ps    depth=(pwlvl(jl)-z(jk))
!ps    !! initialized the water profile according to world ocean altas data (woa)
!ps    !data
!ps    kkk=1
!ps    DO WHILE ( (wdepths(kkk).LE.depth) .AND. (kkk.LT.nwdepth) )
!ps       kkk=kkk+1
!ps    END DO
!ps    IF (kkk.GT.1) kkk=kkk-1  ! restore back one layer
!ps    IF (  &
!ps         (wtfn0(jl,kkk,jrow).NE.xmissing) .AND. &
!ps         (wtfn0(jl,kkk+1,jrow).NE.xmissing) .AND. &
!ps         (wsfn0(jl,kkk,jrow).NE.xmissing) .AND. &
!ps         (wsfn0(jl,kkk+1,jrow).NE.xmissing) ) THEN
!ps      ttt=wtfn0(jl,kkk,jrow)
!ps      ttt1=wtfn0(jl,kkk+1,jrow)
!ps      sss=wsfn0(jl,kkk,jrow)
!ps      sss1=wsfn0(jl,kkk+1,jrow)
!ps      IF ( (depth.LE.wdepths(nwdepth)) ) THEN
!ps      ! depth < max. obs. depth
!ps      ! linearly interpolation (no extrapolation)
!ps        pawtfl0(jl,jk)=ttt+ &
!ps          (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (ttt1-ttt)
!ps        pawsfl0(jl,jk)=sss+ &
!ps          (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (sss1-sss)
!ps      ELSE
!ps      ! depth > max obs. depth
!ps      ! water temperature flux
!ps        IF (l_no_expolation_wt) THEN
!ps        ! no extrapolation
!ps          pawtfl0(jl,jk)=0.
!ps        ELSE
!ps        ! extrapolation
!ps          pawtfl0(jl,jk)=ttt+ &
!ps             (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))*(ttt1-ttt)
!ps        ENDIF
!ps      ! water salinity flux
!ps        IF (l_no_expolation_ws) THEN
!ps        ! no extrapolation
!ps          pawsfl0(jl,jk)=0.
!ps        ELSE
!ps        ! extrapolation
!ps          pawsfl0(jl,jk)=sss+ &
!ps             (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (sss1-sss)
!ps        ENDIF
!ps      ENDIF
!ps    ELSE
!ps      IF (l_no_expolation_wt) THEN
!ps        ! no extrapolation
!ps        pawtfl0(jl,jk)=0.
!ps      ELSE
!ps        pawtfl0(jl,jk)=0.
!ps      ENDIF
!ps      IF (l_no_expolation_ws) THEN
!ps        ! no extrapolation
!ps        pawsfl0(jl,jk)=0.
!ps      ELSE
!ps        pawsfl0(jl,jk)=0.
!ps      ENDIF
!ps    ENDIF
!ps  ENDDO
  pawtfl0(jl,0)=wtfn1st(jl,1,jrow)
  IF (GDCHK3) then
    DO jk=0,0
      WRITE(nerr,*)'pe=',mpp_pe(),',jl=',jl,',z=',z(jk),',pawtfl0=',pawtfl0(jl,jk),',pawsfl0=',pawsfl0(jl,jk)
    ENDDO
    WRITE(nerr,*) "I am leaving interpolation_ocaf0"
  ENDIF

  END SUBROUTINE interpolation_ocaf0


! **********************************************************************                                                                        
  SUBROUTINE nudging_sit_viff_gd_sfc(jl,jk,st_restore_time,ss_restore_time,restore_temp,restore_salt &
                                     ,suv_restore_time,restore_u,restore_v,fratio    &
                                     ,st_restore_time_all,ss_restore_time_all,suv_restore_time_all)
    IMPLICIT NONE
    INTEGER, INTENT(IN):: jl,jk    ! lonitude index, and level
    real, INTENT(IN):: st_restore_time, ss_restore_time, suv_restore_time          
    real, INTENT(IN):: st_restore_time_all, ss_restore_time_all, suv_restore_time_all
    real, INTENT(OUT) :: restore_temp
    real, INTENT(OUT) :: restore_salt    
    real, INTENT(OUT) :: restore_u 
    real, INTENT(OUT) :: restore_v
    real, INTENT(IN):: fratio
    real :: restore_temp2,restore_salt2,restore_u2,restore_v2 
    logical :: lrestore_mix

    IF(GDCHK2) then
          print *,"in nudging_sfc:jl=",jl,",jk=",jk,",plon=",plon(jl),",plat=",plat(jl) &
             ,",pwt=",pwt(jl,jk),",fratio=",fratio
    ENDIF 

    IF ((pobswtb(jl).NE.xmissing).AND.(st_restore_time.GT.0.)) THEN
      restore_temp=(pobswtb(jl)-pwt(jl,jk))*(1.-0.5**(zdtime/st_restore_time))
      IF(lgodas .AND. ltimeblending .AND. timebl_option .EQ. 1 ) THEN
!        restore_temp=(pobswt(jl,jk)-pwt(jl,jk))*(1.-fratio)
        restore_temp=(pobswt(jl,jk)-pwt(jl,jk))*(1.-0.5**(zdtime/st_restore_time))
        restore_temp=(pobswt(jl,jk)-pwt(jl,jk))*(1.-fratio)+restore_temp*fratio
      ELSE
        IF(ltimeblending .AND. timebl_option .EQ. 1 ) THEN
!        restore_temp=(pobswtb(jl)-pwt(jl,jk))*(1.-fratio)
          restore_temp=(pobswtb(jl)-pwt(jl,jk))*(1.-fratio)+restore_temp*fratio
        ENDIF
      ENDIF
    ELSEIF ((pobswtb(jl).NE.xmissing).AND.(st_restore_time.EQ.0.)) THEN
      restore_temp=(pobswtb(jl)-pwt(jl,jk))
    ELSEIF ((st_restore_time.LT.0.)) THEN
      restore_temp=0.
    ELSE
      ! no pobswtb data
      IF (.NOT.ldeep_water_nudg) THEN
      ! no nudging
        restore_temp=0.
      ELSE
      ! assuming restore_temp to that of previous level
      ENDIF
    ENDIF


    IF ((pobswsb(jl).NE.xmissing).AND.(ss_restore_time.GT.0.)) THEN
      restore_salt=(pobswsb(jl)-pws(jl,jk))*(1.-0.5**(zdtime/ss_restore_time))
      IF(ltimeblending .AND. timebl_option .EQ. 1) THEN
!        restore_salt=(pobswsb(jl)-pws(jl,jk))*(1.-fratio)
        restore_salt=(pobswsb(jl)-pws(jl,jk))*(1.-fratio)+restore_salt*fratio
      ENDIF
    ELSEIF ((pobswsb(jl).NE.xmissing).AND.(ss_restore_time.EQ.0.)) THEN
      restore_salt=(pobswsb(jl)-pws(jl,jk))
    ELSEIF ((ss_restore_time.LT.0.)) THEN
      restore_salt=0.
    ELSE
    ! no pobswsb data
      IF (.NOT.ldeep_water_nudg) THEN
      ! no nudging
        restore_salt=0.
      ELSE
      ! assuming restore_salt to that of previous level
      ENDIF
    ENDIF


    IF ((pobswu(jl,jk).NE.xmissing).AND.(suv_restore_time.GT.0.)) THEN
       restore_u=(pobswu(jl,jk)-pwu(jl,jk))*(1.-0.5**(zdtime/suv_restore_time))
       IF(ltimeblending .AND. timebl_option .EQ. 1) THEN
!         restore_u=(pobswu(jl,jk)-pwu(jl,jk))*(1.-fratio)
         restore_u=(pobswu(jl,jk)-pwu(jl,jk))*(1.-fratio)+restore_u*fratio
       ENDIF
    ELSEIF ((pobswu(jl,jk).NE.xmissing).AND.(suv_restore_time.EQ.0.)) THEN
       restore_u=(pobswu(jl,jk)-pwu(jl,jk))
    ELSEIF ((suv_restore_time.LT.0.)) THEN
       restore_u=0.
    ELSE
      ! no pobswtb data
      IF (.NOT.ldeep_water_nudg) THEN
      ! no nudging
        restore_u=0.
      ELSE
      ! assuming restore_temp to that of previous level
      ENDIF
    ENDIF


    IF ((pobswv(jl,jk).NE.xmissing).AND.(suv_restore_time.GT.0.)) THEN
       restore_v=(pobswv(jl,jk)-pwv(jl,jk))*(1.-0.5**(zdtime/suv_restore_time))
       IF(ltimeblending .AND. timebl_option .EQ. 1) THEN
!         restore_v=(pobswv(jl,jk)-pwv(jl,jk))*(1.-fratio)
         restore_v=(pobswv(jl,jk)-pwv(jl,jk))*(1.-fratio)+restore_v*fratio
       ENDIF
    ELSEIF ((pobswv(jl,jk).NE.xmissing).AND.(suv_restore_time.EQ.0.)) THEN
       restore_v=(pobswv(jl,jk)-pwv(jl,jk))
    ELSEIF ((suv_restore_time.LT.0.)) THEN
       restore_v=0.
    ELSE
    ! no pobswtb data
      IF (.NOT.ldeep_water_nudg) THEN
      ! no nudging
       restore_v=0.
      ELSE
      ! assuming restore_temp to that of previous level
      ENDIF
    ENDIF

    lrestore_mix=.false. 
    IF (nsit_nudg .GE. 1) THEN
      DO isit_nudg= 1, nsit_nudg
        IF( (plon(jl) .GE. sitbox_nudg_w(isit_nudg)) .AND. &
                (plon(jl) .LT. sitbox_nudg_e(isit_nudg)) .AND. &
                (plat(jl) .GE. sitbox_nudg_s(isit_nudg)) .AND. &
                (plat(jl) .LT. sitbox_nudg_n(isit_nudg)) ) THEN

          IF( ( (plon(jl) .GE. sitbox_nudg_w(isit_nudg)) .AND. &
                (plon(jl) .LT. sitbox_nudg_w(isit_nudg)+1.) ) .OR. &
              ( (plon(jl) .GE. sitbox_nudg_e(isit_nudg)-1.) .AND. &
                (plon(jl) .LT. sitbox_nudg_e(isit_nudg)) ) ) THEN
            lrestore_mix=.true.
          ENDIF
          IF( ( (plat(jl) .GE. sitbox_nudg_s(isit_nudg)) .AND. &
                    (plat(jl) .LT. sitbox_nudg_s(isit_nudg)+1.) ) .OR. &
                  ( (plat(jl) .GE. sitbox_nudg_n(isit_nudg)-1.) .AND. &
                    (plat(jl) .LT. sitbox_nudg_n(isit_nudg)) ) ) THEN
            lrestore_mix=.true.
          ENDIF
        ENDIF
      ENDDO 
    ENDIF

    IF(lrestore_mix) then
      IF ((pobswtb(jl).NE.xmissing).AND.(st_restore_time_all.GT.0.)) THEN
         restore_temp2=(pobswtb(jl)-pwt(jl,jk))*(1.-0.5**(zdtime/st_restore_time_all))
         IF(ltimeblending .AND. timebl_option .EQ. 1) THEN
!           restore_temp2=(pobswtb(jl)-pwt(jl,jk))*(1.-fratio)
           restore_temp2=(pobswtb(jl)-pwt(jl,jk))*(1.-fratio)+restore_temp2*fratio
         ENDIF
      ELSEIF ((pobswtb(jl).NE.xmissing).AND.(st_restore_time_all.EQ.0.)) THEN
         restore_temp2=(pobswtb(jl)-pwt(jl,jk))
      ELSEIF ((st_restore_time_all.LT.0.)) THEN
         restore_temp2=0.
      ELSE
        ! no pobswtb data
        IF (.NOT.ldeep_water_nudg) THEN
        ! no nudging
          restore_temp2=0.
        ELSE
        ! assuming restore_temp to that of previous level
        ENDIF
      ENDIF

      IF ((pobswsb(jl).NE.xmissing).AND.(ss_restore_time_all.GT.0.)) THEN
        restore_salt2=(pobswsb(jl)-pws(jl,jk))*(1.-0.5**(zdtime/ss_restore_time_all))
        IF(ltimeblending .AND. timebl_option .EQ. 1) THEN
!          restore_salt2=(pobswsb(jl)-pws(jl,jk))*(1.-fratio)
          restore_salt2=(pobswsb(jl)-pws(jl,jk))*(1.-fratio)+restore_salt2*fratio
        ENDIF
      ELSEIF ((pobswsb(jl).NE.xmissing).AND.(ss_restore_time_all.EQ.0.)) THEN
        restore_salt2=(pobswsb(jl)-pws(jl,jk))
      ELSEIF ((ss_restore_time_all.LT.0.)) THEN
        restore_salt2=0.
      ELSE
        ! no pobswsb data
        IF (.NOT.ldeep_water_nudg) THEN
        ! no nudging
          restore_salt2=0.
        ELSE
        ! assuming restore_salt to that of previous level
        ENDIF
      ENDIF

      IF ((pobswu(jl,jk).NE.xmissing).AND.(suv_restore_time_all.GT.0.)) THEN
         restore_u2=(pobswu(jl,jk)-pwu(jl,jk))*(1.-0.5**(zdtime/suv_restore_time_all))
         IF(ltimeblending .AND. timebl_option .EQ. 1) THEN
!           restore_u2=(pobswu(jl,jk)-pwu(jl,jk))*(1.-fratio)
           restore_u2=(pobswu(jl,jk)-pwu(jl,jk))*(1.-fratio)+restore_u2*fratio
         ENDIF
      ELSEIF ((pobswu(jl,jk).NE.xmissing).AND.(suv_restore_time_all.EQ.0.)) THEN
         restore_u2=(pobswu(jl,jk)-pwu(jl,jk))
      ELSEIF ((suv_restore_time_all.LT.0.)) THEN
         restore_u2=0.
      ELSE
        ! no pobswtb data
        IF (.NOT.ldeep_water_nudg) THEN
        ! no nudging
         restore_u2=0.
        ELSE
        ! assuming restore_temp to that of previous level
        ENDIF
      ENDIF

      IF ((pobswv(jl,jk).NE.xmissing).AND.(suv_restore_time_all.GT.0.)) THEN
         restore_v2=(pobswv(jl,jk)-pwv(jl,jk))*(1.-0.5**(zdtime/suv_restore_time_all))
         IF(ltimeblending .AND. timebl_option .EQ. 1) THEN
!           restore_v2=(pobswv(jl,jk)-pwv(jl,jk))*(1.-fratio)
           restore_v2=(pobswv(jl,jk)-pwv(jl,jk))*(1.-fratio)+restore_v2*fratio
         ENDIF
      ELSEIF ((pobswv(jl,jk).NE.xmissing).AND.(suv_restore_time_all.EQ.0.)) THEN
         restore_v2=(pobswv(jl,jk)-pwv(jl,jk))
      ELSEIF ((suv_restore_time_all.LT.0.)) THEN
         restore_v2=0.
      ELSE
        ! no pobswtb data
        IF (.NOT.ldeep_water_nudg) THEN
        ! no nudging
          restore_v2=0.
        ELSE
        ! assuming restore_temp to that of previous level
        ENDIF
      ENDIF

      restore_temp=0.5*(restore_temp+restore_temp2) 
      restore_salt=0.5*(restore_salt+restore_salt2) 
      restore_u=0.5*(restore_u+restore_u2) 
      restore_v=0.5*(restore_v+restore_v2) 
    ENDIF

!    IF(GDCHK3) then
!      print *,"nudging_sfc: jl=",jl,",jk=",jk,",plon=",plon(jl),",plat=",plat(jl) &
!             ,",pwt=",pwt(jl,jk) &
!             ,",fratio=",fratio,",lrestore_mix=",lrestore_mix &
!             ,",st_restore_time=",st_restore_time  &
!             ,",st_restore_time_all=",st_restore_time_all &
!             ,",ss_restore_time=",ss_restore_time  &
!             ,",ss_restore_time_all=",ss_restore_time_all &
!             ,",suv_restore_time=",suv_restore_time  &
!             ,",suv_restore_time_all=",suv_restore_time_all &
!             ,",restore_temp=",restore_temp,",restore_temp2=",restore_temp2 &
!             ,",restore_salt=",restore_salt,",restore_salt2=",restore_salt2 &
!             ,",restore_u=",restore_u,",restore_u2=",restore_u2 &
!             ,",restore_v=",restore_v,",restore_v2=",restore_v2
!    ENDIF



  END SUBROUTINE nudging_sit_viff_gd_sfc
! **********************************************************************                                                                        
  SUBROUTINE nudging_sit_viff_gd(jl,jrow,obox_restore_time, &
    socn_restore_time,uocn_restore_time,docn_restore_time,   &
    ssit_restore_time,ssits_restore_time,ssituv_restore_time, &       
    pssit_restore_time,pusit_restore_time,pdsit_restore_time,   &
    pssits_restore_time,pusits_restore_time,pdsits_restore_time, &
    pssituv_restore_time,pusituv_restore_time,pdsituv_restore_time, tau, tauhr) 

    IMPLICIT NONE
    INTEGER, INTENT(IN):: jl,jrow    ! lonitude index, and latitude index
    real, INTENT(IN):: obox_restore_time  ! nudging restore time for obox_mask(jl).GT.0. grids
    real, INTENT(IN):: socn_restore_time ! surface [0m,10m) ocean restore time for pocnmask(jl).GT.0. grids
    real, INTENT(IN):: uocn_restore_time ! upper [10m,100m) ocean restore time for pocnmask(jl).GT.0. grids 
    real, INTENT(IN):: docn_restore_time ! deep (>100 m) ocean restore time for pocnmask(jl).GT.0. grids
    real, INTENT(IN):: ssit_restore_time   ! surface [0m,10m) ocean temp/salt restore time for the rest SIT grids
    real, INTENT(IN):: ssits_restore_time   ! surface [0m,10m) ocean temp/salt restore time for the rest SIT grids for salinity
    real, INTENT(IN):: ssituv_restore_time   ! surface [0m,10m) ocean temp/salt restore time for the rest SIT grids for u, v
    real, INTENT(IN):: pssit_restore_time   ! surface [0m,10m) ocean temp/salt restore time for the rest SIT grids
    real, INTENT(IN):: pusit_restore_time   ! upper [10m,100m) ocean temp/salt restore time for the rest SIT grids
    real, INTENT(IN):: pdsit_restore_time   ! deep (>100 m) ocean temp/salt restore time for the rest SIT grids
    real, INTENT(IN):: pssits_restore_time   ! surface [0m,10m) ocean temp/salt restore time for the rest SIT grids for salinity
    real, INTENT(IN):: pusits_restore_time   ! upper [10m,100m) ocean temp/salt restore time for the rest SIT grids for salinity
    real, INTENT(IN):: pdsits_restore_time   ! deep (>100 m) ocean temp/salt restore time for the rest SIT grids for salinity
    real, INTENT(IN):: pssituv_restore_time   ! surface [0m,10m) ocean temp/salt restore time for the rest SIT grids for u, v
    real, INTENT(IN):: pusituv_restore_time   ! upper [10m,100m) ocean temp/salt restore time for the rest SIT grids for u, v
    real, INTENT(IN):: pdsituv_restore_time   ! deep (>100 m) ocean temp/salt restore time for the rest SIT grids for u, v
    INTEGER :: jk,jk10,jk100,jkmixing
    real:: restore_temp 
    real:: restore_salt
    real:: restore_u 
    real:: restore_v
    real:: st_restore_time, ss_restore_time, suv_restore_time
    real:: st_restore_time_all, ss_restore_time_all, suv_restore_time_all
!ps
! time blending
    real:: tau,tauhr,tmpt0,tmpt1,tmpt2,xtau,gm,fratio
    real:: urand
    logical:: lbefore_bl

!   nudging from surface to nudge_depth=100 m depth (default)
    IF (GDCHK2) then
       WRITE(nerr,*) ", I am in nudging_sit_viff_gd"

       WRITE(nerr,*) "ssit_restore_time=",pssit_restore_time
       WRITE(nerr,*) "usit_restore_time=",pusit_restore_time
       WRITE(nerr,*) "dsit_restore_time=",pdsit_restore_time
       
       WRITE(nerr,*) "socn_restore_time=",socn_restore_time
       WRITE(nerr,*) "uocn_restore_time=",uocn_restore_time
       WRITE(nerr,*) "docn_restore_time=",docn_restore_time
       WRITE(nerr,*) "pocnmask(",jl,",)=",pocnmask(jl)
    ENDIF
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!!!    IF (pobswtb(jl).NE.xmissing) pobswtb(jl)=MAX(pobswtb(jl),pctfreez2(jl))
    IF (lsoil) RETURN
    IF (lgodas) THEN
      CALL interpolation_godas(jl,jrow,.TRUE.,.TRUE.,.TRUE.,.TRUE.)
    ELSE
      CALL interpolation_godas(jl,jrow,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
    ENDIF
      IF (pobsws(jl,1).NE.xmissing) pobswsb(jl)=pobsws(jl,1)     

    IF ((obox_restore_time.LT.0.).AND.(socn_restore_time.LT.0.).AND.    &
        ( uocn_restore_time.LT.0.).AND.( docn_restore_time.LT.0.).AND. &
        ( pssit_restore_time.LT.0.).AND.( pusit_restore_time.LT.0.).AND. &
        ( pdsit_restore_time.LT.0.) ) RETURN                               ! no nudging    
!
!   1.0 Find jk at 10 m depth and nudge_depth (=100 m)
    jk=nls+1
    DO WHILE ( ((pwlvl(jl)-z(jk)).LT.nudge_depth_10).AND.(jk.LT.nle) ) 
!      nudging at depth just >= 10 m or at level nle (last water level)
       jk=jk+1
    END DO
    jk10=jk
    DO WHILE ( ((pwlvl(jl)-z(jk)).LT.mixedlayer(jl)).AND.(jk.LT.nle) )
!      nudging at depth just>= mixedlayer (55 m)
       jk=jk+1
    END DO
    jkmixing=jk
    IF (GDCHK3) then
      WRITE(nerr,*) ",mixedlayer=",mixedlayer(jl),",jkmixing=",jkmixing
    ENDIF
    DO WHILE ( ((pwlvl(jl)-z(jk)).LT.nudge_depth_100).AND.(jk.LT.nle) ) 
!      nudging at depth >= nudge_depth (100 m)
       jk=jk+1
    END DO
    jk100=jk


!!! time blending
    tmpt0=0.
    fratio=1.
    IF(ltimeblending)THEN
! unit day
      tmpt1=timebl_start*24. !(convert from day to hour)

      IF(tau .le. tmpt1)then
        fratio=0.
        gm=-1.
        lbefore_bl=.true.
      ELSE
        lbefore_bl=.false.
        IF(timebl_option .EQ. 0) THEN
          tmpt0=(timebl_start+timebl_allsit)/2.*24.
          tmpt2=tmpt0-tmpt1
          xtau=tau-tmpt1

          gm=abs(0.5*xtau/tmpt2-1.)

          if(xtau .ge. timebl_allsit*24.)then
            fratio=1.                           !! use st/ss/suv_restore_time
          else if(gm.gt.0. .and. gm.lt.1.)then
            fratio=1.-exp(-1./gm*exp(0.5/(gm-1.)))
          else
            fratio=0.
          endif
        ENDIF

        IF(timebl_option .EQ. 1)then
          fratio=(sin((tauhr-6.)/12.*api)+1.)*0.5
        ENDIF
      ENDIF
      IF(fratio .eq. 0.) then
        ltrigsit=.false.
      else
        ltrigsit=.true.
      endif
      IF(GDCHK3) then
       WRITE(nerr,*) "ltimeblending=",ltimeblending  &
                  ,",timebl_start=",timebl_start,",timebl_allsit=",timebl_allsit  &
                  ,",tmpt1=",tmpt1,",tmpt0=",tmpt0,",nls=",nls,",fratio=",fratio
      ENDIF
    ENDIF

!!! end time blending

!    IF(GDCHK3) print *,'in nudge_gd: nudge_depth_10=',nudge_depth_10,',jk10=',jk10        
!   3.0 Start nudging from 10 m (at nle, 10 m and >=100 m) 
    restore_temp=0.
    restore_salt=0.    
    restore_u=0.    
    restore_v=0.    
!ps    jk=nls+1
    jk=nls
!ps    DO WHILE ( jk.LE.nle )
    DO WHILE ( jk.LE.nle+1 )
      ! [0m, 10m (or waterbed) )
      ! .?: (first digit): T,S nudging. 0: no nudging, 1: nudging inside square boxes, 2: nudging outside the square boxes
      ! ?.: (decimal digit): high_damping coefficient for u,v currents in Sulu Sea, no nudging, 1: nudging inside square boxes, 2: nudging outside the square boxes
      ! obox_nudg_flag
      ! 0             |    0         |    0: no nudging
      ! 1:  0- 10 m   |    1: T,S    |    1: nudging inside
      ! 2: 10-100 m   |    2: U,V    |    2: nudging outside
      ! 4    >100 m   |              |    
      IF(lmixedlayer .and. jk .LT. jkmixing) THEN
      ! [0m, mixing (or waterbed) )
        IF ((obox_mask(jl).GT.0.).AND.(MOD(INT(obox_nudg_flag/100),2).EQ.1)) THEN
        ! obox grids
        ! ie., obox_mask(jl).GT.0. .AND. obox_nudg_flag=1xx
        !
          st_restore_time=obox_restore_time
          ss_restore_time=obox_restore_time
          suv_restore_time=obox_restore_time
        ELSEIF (pocnmask(jl).GT.0.) THEN
        ! ocn grids
          st_restore_time=socn_restore_time
          ss_restore_time=socn_restore_time
          suv_restore_time=socn_restore_time
        ELSE
        ! other SIT grids
!ps          ss_restore_time=ssit_restore_time
          st_restore_time=pssit_restore_time
          ss_restore_time=pssits_restore_time
          suv_restore_time=pssituv_restore_time
          st_restore_time_all=ssit_restore_time
          ss_restore_time_all=ssits_restore_time
          suv_restore_time_all=ssituv_restore_time
        ENDIF
      ELSEIF ( jk.LT.jk10) THEN
      ! [0m, 10m (or waterbed) )
        IF ((obox_mask(jl).GT.0.).AND.(MOD(INT(obox_nudg_flag/100),2).EQ.1)) THEN
        ! obox grids
        ! ie., obox_mask(jl).GT.0. .AND. obox_nudg_flag=1xx
        !        
          st_restore_time=obox_restore_time
          ss_restore_time=obox_restore_time
          suv_restore_time=obox_restore_time
        ELSEIF (pocnmask(jl).GT.0.) THEN
        ! ocn grids
          st_restore_time=socn_restore_time
          ss_restore_time=socn_restore_time
          suv_restore_time=socn_restore_time
        ELSE
        ! other SIT grids
!ps          ss_restore_time=ssit_restore_time
          st_restore_time=pssit_restore_time
          ss_restore_time=pssits_restore_time
          suv_restore_time=pssituv_restore_time
          st_restore_time_all=ssit_restore_time
          ss_restore_time_all=ssits_restore_time
          suv_restore_time_all=ssituv_restore_time
        ENDIF
      ELSEIF ( jk.LT.jk100) THEN
      ! [10m, 100m(waterbed) )
        IF ((obox_mask(jl).GT.0.).AND.(MOD(INT(INT(obox_nudg_flag/100)/2),2).EQ.1)) THEN
        ! ie., obox_mask(jl).GT.0. .AND. obox_nudg_flag=2xx
        !
          st_restore_time=obox_restore_time
          ss_restore_time=obox_restore_time
        ELSEIF (pocnmask(jl).GT.0.) THEN
        ! no nudging for depth < 100 m (default)
          st_restore_time=uocn_restore_time
          ss_restore_time=uocn_restore_time
        ELSE  
!ps          ss_restore_time=usit_restore_time
          st_restore_time=pusit_restore_time
          ss_restore_time=pusits_restore_time
          suv_restore_time=pusituv_restore_time
        ENDIF
      ELSE
      ! >= 100 m or at nle
        IF ((obox_mask(jl).GT.0.).AND.(MOD(INT(INT(INT(obox_nudg_flag/100)/2)/2),2).EQ.1)) THEN
        ! ie., obox_mask(jl).GT.0. .AND. obox_nudg_flag=4xx
        !
          st_restore_time=obox_restore_time
          ss_restore_time=obox_restore_time
        ELSEIF (pocnmask(jl).GT.0.) THEN
        ! no nudging
          st_restore_time=docn_restore_time
          ss_restore_time=docn_restore_time
        ELSE
!ps          ss_restore_time=dsit_restore_time
          st_restore_time=pdsit_restore_time
          ss_restore_time=pdsits_restore_time
          suv_restore_time=pdsituv_restore_time
        ENDIF
      ENDIF

      IF(ltimeblending) THEN
        IF(st_restore_time .GE. 0.) THEN
          st_restore_time=fratio*st_restore_time
        ELSE
          IF(timebl_option .EQ. 0 .AND. fratio .LT. 1) THEN
            st_restore_time=fratio*(2.*30.*24.*3600.) !2m
          ENDIF
          IF(timebl_option .EQ. 1 .AND. lbefore_bl ) THEN
            st_restore_time=0.
          ENDIF
        ENDIF
        IF(st_restore_time_all .GE. 0.) THEN
          st_restore_time_all=fratio*st_restore_time_all
        ELSE
          IF(timebl_option .EQ. 0 .AND. fratio .LT. 1) THEN
            st_restore_time_all=fratio*(2.*30.*24.*3600.) !2m
          ENDIF
          IF(timebl_option .EQ. 1 .AND. lbefore_bl ) THEN
            st_restore_time_all=0.
          ENDIF
        ENDIF


        IF(ss_restore_time .GE. 0.) THEN
          ss_restore_time=fratio*ss_restore_time
        ELSE
          IF(timebl_option .EQ. 0 .AND. fratio .LT. 1.) THEN
            ss_restore_time=fratio*(2.*30.*24.*3600.) !2m
          ENDIF
          IF(timebl_option .EQ. 1 .AND. lbefore_bl ) THEN
            ss_restore_time=0.
          ENDIF
        ENDIF
        IF(ss_restore_time_all .GE. 0.) THEN
          ss_restore_time_all=fratio*ss_restore_time_all
        ELSE
          IF(timebl_option .EQ. 0 .AND. fratio .LT. 1.) THEN
            ss_restore_time_all=fratio*(2.*30.*24.*3600.) !2m
          ENDIF
          IF(timebl_option .EQ. 1 .AND. lbefore_bl ) THEN
            ss_restore_time_all=0.
          ENDIF
        ENDIF


        IF(suv_restore_time .GE. 0.) THEN
          suv_restore_time=fratio*suv_restore_time
        ELSE
          IF(timebl_option .EQ. 0 .AND. fratio .LT. 1.) THEN
            suv_restore_time=fratio*(2.*30.*24.*3600.) !2m
          ENDIF
          IF(timebl_option .EQ. 1 .AND. lbefore_bl ) THEN
            suv_restore_time=0.
          ENDIF
        ENDIF
        IF(suv_restore_time_all .GE. 0.) THEN
          suv_restore_time_all=fratio*suv_restore_time_all
        ELSE
          IF(timebl_option .EQ. 0 .AND. fratio .LT. 1.) THEN
            suv_restore_time_all=fratio*(2.*30.*24.*3600.) !2m
          ENDIF
          IF(timebl_option .EQ. 1 .AND. lbefore_bl ) THEN
            suv_restore_time_all=0.
          ENDIF
        ENDIF


!        IF(GDCHK3 .AND. jk.EQ.1) print *,"ltimeblending=",ltimeblending, &
!                       ",fratio=",fratio,",jk=",jk, &
!                       ",st_restore_time=",st_restore_time, &
!                       ",ss_restore_time=",ss_restore_time, &
!                       ",suv_restore_time=",suv_restore_time, &
!                       ",st_restore_time_all=",st_restore_time_all, &
!                       ",ss_restore_time_all=",ss_restore_time_all, &
!                       ",suv_restore_time_all=",suv_restore_time_all
      ENDIF

!      
!     3.1 Determining restore temperature at each level
      IF(lmixedlayer .and. jk .LT. jkmixing) THEN
      ! nudging to obstsw at [0m, 10m] or at waterbed(=nle)
        CALL nudging_sit_viff_gd_sfc(jl,jk,st_restore_time,ss_restore_time,restore_temp,restore_salt,&
                                     suv_restore_time,restore_u,restore_v,fratio,                    &
                                     st_restore_time_all,ss_restore_time_all,suv_restore_time_all)
!ps      ELSEIF ( jk.LE.jk10) THEN
      ELSEIF ( jk.LT.jk10) THEN
      ! nudging to obstsw at [0m, 10m] or at waterbed(=nle)
        CALL nudging_sit_viff_gd_sfc(jl,jk,st_restore_time,ss_restore_time,restore_temp,restore_salt,&
                                     suv_restore_time,restore_u,restore_v,fratio,                    &
                                     st_restore_time_all,ss_restore_time_all,suv_restore_time_all)
      ELSE
      ! for depth > 10 m
        IF (lgodas .OR. lwoa0) THEN
          IF ((pobswt(jl,jk).NE.xmissing).AND.(st_restore_time.GT.0.)) THEN
            restore_temp=(pobswt(jl,jk)-pwt(jl,jk))*(1.-0.5**(zdtime/st_restore_time))
            ! obsws  vertical index start at 1
          ELSEIF ((pobswt(jl,jk).NE.xmissing).AND.(st_restore_time.EQ.0.)) THEN
            restore_temp=(pobswt(jl,jk)-pwt(jl,jk))
          ELSEIF (st_restore_time.LT.0.) THEN
          !! no nudging
            restore_temp=0.
          ELSE
            IF (.NOT.ldeep_water_nudg) THEN
              restore_temp=0.
            ELSE
            !! no pobswt data, assume restore_temp to be that of the previous level
            ENDIF
          ENDIF
          IF ((pobsws(jl,jk).NE.xmissing).AND.(ss_restore_time.GT.0.)) THEN
            restore_salt=(pobsws(jl,jk)-pws(jl,jk))*(1.-0.5**(zdtime/ss_restore_time))
            ! obsws  vertical index start at 1
          ELSEIF ((pobsws(jl,jk).NE.xmissing).AND.(ss_restore_time.EQ.0.)) THEN
            restore_salt=(pobsws(jl,jk)-pws(jl,jk))
          ELSEIF (ss_restore_time.LT.0.) THEN
            restore_salt=0.
          ELSE
            IF (.NOT.ldeep_water_nudg) THEN
              restore_salt=0.
            ELSE
            !! no pobsws data, assume restore_temp to be that of the previous level
            ENDIF
          ENDIF
          IF ((pobswu(jl,jk).NE.xmissing).AND.(suv_restore_time.GT.0.)) THEN
            restore_u=(pobswu(jl,jk)-pwu(jl,jk))*(1.-0.5**(zdtime/suv_restore_time))
            ! obsws  vertical index start at 1
          ELSEIF ((pobswu(jl,jk).NE.xmissing).AND.(suv_restore_time.EQ.0.)) THEN
            restore_u=(pobswu(jl,jk)-pwu(jl,jk))
          ELSEIF (suv_restore_time.LT.0.) THEN
            restore_u=0.
          ELSE
            IF (.NOT.ldeep_water_nudg) THEN
              restore_u=0.
            ELSE
            !! no pobsws data, assume restore_temp to be that of the previous
            !level
            ENDIF
          ENDIF
          IF ((pobswv(jl,jk).NE.xmissing).AND.(suv_restore_time.GT.0.)) THEN
            restore_v=(pobswv(jl,jk)-pwv(jl,jk))*(1.-0.5**(zdtime/suv_restore_time))
            ! obsws  vertical index start at 1
          ELSEIF ((pobswv(jl,jk).NE.xmissing).AND.(suv_restore_time.EQ.0.)) THEN
            restore_v=(pobswv(jl,jk)-pwv(jl,jk))
          ELSEIF (suv_restore_time.LT.0.) THEN
            restore_v=0.
          ELSE
            IF (.NOT.ldeep_water_nudg) THEN
              restore_v=0.
            ELSE
            !! no pobsws data, assume restore_temp to be that of the previous
            !level
            ENDIF
          ENDIF
        ELSE
         ! no woa data
           IF (st_restore_time.LT.0.) THEN
           !! no nudging
             restore_temp=0.
           ELSE
             IF (.NOT.ldeep_water_nudg) THEN
             !! no nudging
               restore_temp=0.
             ELSE
             !! no pobswt data, assume restore_temp to be that of the previous level
             ENDIF
           ENDIF
           IF (ss_restore_time.LT.0.) THEN
             restore_salt=0.
           ELSE
             IF (.NOT.ldeep_water_nudg) THEN
             !! no nudging
               restore_salt=0.
             ELSE
             !! no pobsws data, assume restore_temp to be that of the previous level
             ENDIF
           ENDIF
         ENDIF
      ENDIF
    
!     3.2 Nudging at each level     
!      IF(GDCHK3) print *,"before restore,jk=",jk,",pwt(jl,jk)=",pwt(jl,jk),",restore_temp=",restore_temp
      pwt(jl,jk)= pwt(jl,jk)+restore_temp
!      IF(GDCHK3) print *,"after restore,jk=",jk,",pwt(jl,jk)=",pwt(jl,jk),",restore_temp=",restore_temp
      pws(jl,jk)= pws(jl,jk)+restore_salt

      pwu(jl,jk)= pwu(jl,jk)+restore_u
      pwv(jl,jk)= pwv(jl,jk)+restore_v

!ps      IF(GDCHK3) WRITE(nerr,*) "before random:pwu(",jl,",",jk,")=",pwu(jl,jk),",pwv=",pwv(jl,jk)
!ps      call random_number(urand)
!ps      urand=0.01*(urand*2.-1.)      !-0.05<=urand<0.05
!ps      pwu(jl,jk)=pwu(jl,jk)*(1.+urand)
!ps      pwv(jl,jk)=pwv(jl,jk)*(1.+urand)
!ps      IF(GDCHK3) then
!ps        WRITE(nerr,*) "after random: urand=",urand,",pwu(",jl,",",jk,")=",pwu(jl,jk),",pwv=",pwv(jl,jk)
!ps      ENDIF

!     calc acc. nudging flux (positive into the water column)
!     Note that energy in skin layer (jk=0) has been accounted in the first water layer (jk=1)
!ps      IF ((jk.NE.0).AND.(hw(jk).NE.xmissing)) THEN
!ps
!      IF(GDCHK3) then
!        WRITE(nerr,*) "before,myrank=",myrank,",pwtfn(",jl,",",jk,")=",pwtfn(jl,jk)&
!                ,",restore_temp=",restore_temp,"pawtfl(jl,jk)=",pawtfl(jl,jk)&
!                ,",zdtime=",zdtime,",pobswt(jl,jk)=",pobswt(jl,jk)&
!                ,"pwt(jl,jk)=",pwt(jl,jk),"st_restore_time=",st_restore_time
!      ENDIF
!ps
        pwtfn(jl,jk)=pwtfn(jl,jk)+restore_temp+pawtfl(jl,jk)*zdtime           ! v9.8 including added pawtfl
        pwsfn(jl,jk)=pwsfn(jl,jk)+restore_salt+pawsfl(jl,jk)*zdtime           ! v9.8 including added pawsfl
        pwtfn0(jl,jk)=pwtfn0(jl,jk)+restore_temp                              ! v9.8 including added pawtfl
        pwsfn0(jl,jk)=pwsfn0(jl,jk)+restore_salt                              ! v9.8 including added pawsfl
!ps      ENDIF
!ps
      IF(GDCHK2) then
        WRITE(nerr,*) "after,myrank=",myrank,",pwtfn(",jl,",",jk,")=",pwtfn(jl,jk) &
                     ,",pawtfl=",pawtfl(jl,jk)
        WRITE(nerr,*) "after,myrank=",myrank,",pwtfn0(",jl,",",jk,")=",pwtfn0(jl,jk)
      ENDIF
!ps
      jk=jk+1
    END DO

    IF (GDCHK2) then    
       WRITE(nerr,*) "I am leaving nudging_sit_viff_gd"
    ENDIF
  END SUBROUTINE nudging_sit_viff_gd
  ! ----------------------------------------------------------------------
  SUBROUTINE output2
  !!!! A DUPLICATE VERSION of sit_vdiff_init 
  !
  !-----------------------------------------------------------------
  !
  !*    OUTPUT: WRITE DIAGNOSTIC VARIABLES
  !
  !-----------------------------------------------------------------
    !!! USE mo_mpi,           ONLY: mpp_pe()   
    IMPLICIT NONE
  
    INTEGER jk
    real:: tmp
    
  !
  !!!  CALL pzcord(jl,jrow)   ! bjt 2010/02/21
    WRITE(nerr,*) "kbdim=", kbdim, "kproma=",kproma
    IF (nle.GE.1) THEN
  !     water exists
      WRITE(nerr,*) "Water exist, nle=", nle
    ELSE
  !     soil only
      WRITE(nerr,*) "Soil only, nle=", nle
    ENDIF
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","","","we,","lw,","tsi0,","tsi1,","","","obswtb,","ptsw,"
    WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,0,0.,pzsi(jl,0),psilw(jl,0),&
  &  ptsnic(jl,0),ptsnic(jl,1),0.,0.,pobswtb(jl),ptsw(jl)

    WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,1,1.,pzsi(jl,1),psilw(jl,1),&
  &  ptsnic(jl,2),ptsnic(jl,3),0.,0.,pobswtb(jl),ptsw(jl)
  !  
    WRITE(nerr,2300) "pe,","jl,","row,","lat,","lon,","step,","k,","z,", "wt,", "ws,", "wu,", "wv,",&
  &  "wtke,", "wkh,","wldisp,","wlmx,","awtfl,","awtfl0,","awufl,","awvfl,","awsfl,","awtkefl,","obswt,",",rhom,","sftobswt,"
    DO jk = 0, nle+1
      WRITE(nerr,2301) mpp_pe(),jl,jrow,plat(jl),plon(jl),istep,jk,z(jk),&
  &     pwt(jl,jk),pws(jl,jk),pwu(jl,jk),pwv(jl,jk),&
  &     pwtke(jl,jk),pwkh(jl,jk),pwldisp(jl,jk),pwlmx(jl,jk), &
  &     pawtfl(jl,jk)*1000.,pawtfl0(jl,jk),pawufl(jl,jk)*1000.,pawvfl(jl,jk)*1000.,pawsfl(jl,jk)*1000.,pawtkefl(jl,jk), &
  &     pobswt(jl,jk),rhom(jk),psftobswt(jl,jk)
    END DO
      
  !
    RETURN
  !
   2200 FORMAT(1X,3(A4),2(A9),1(A4),2(2(A10),2(A9)),1(A9), 1(A10))
  !ps 2201 FORMAT(1X,4(I3,","),1(I9,","),2(2(E9.2,","),2(F8.3,",")),&
   2201 FORMAT(1X,3(I3,","),2(F8.2,","),1(I3,","),4(F10.3,","),1(I9,","),2(2(E9.2,","),2(F8.3,",")),&
  &       1(F8.3,","), 1(E10.3,","))
  !
!ps   2300 FORMAT(1X,4A4,     1A10,1A4,1A11,4A9,2A10,10A9)
      2300 FORMAT(1X, A4,A5,A4,2A8,A6,A4,1A11,4A9,2A10,11A9)
  !ps 2301 FORMAT(1X,4(I3,","),1(I9,","),1(I3,","),1(F10.4,","),1(F8.3,","),2(F8.4,","),&
  !!! 2301 FORMAT(1X,3(I3,","),2(F8.2,","),2(I3,","),1(F10.0,","),1(F8.3,","),1(F8.4,","),&
  !!!&       2(F8.3,","),2(E9.2,","),10(F8.3,","))
!ps   2301 FORMAT(1X,3(I3,","),2(F8.2,","),2(I3,","),1(F10.0,","),1(F8.3,","),1(F8.3,","),&
!ps  &       2(F8.2,","),2(E9.2,","),10(F8.3,","))  
      2301 FORMAT(1X,(I3,","),(I4,","),(I3,","),2(F7.2,","),(I5,","),(I3,","),1(F10.4,","),1(F8.3,","),1(F8.4,","),&
  &       2(F8.4,","),2(E9.2,","),11(F8.3,","))
  !
  END SUBROUTINE output2
END SUBROUTINE sit_vdiff
! **********************************************************************

!!!  SUBROUTINE sit_vdiff_end (Sit)
  SUBROUTINE sit_vdiff_end ()
    IMPLICIT NONE  
    !!! type (sit_vdiff_data_type), intent(inout) :: Sit
    integer           :: k
    integer           :: unit
    character(len=22) :: restart='RESTART/sit_vdiff.res'

    !!!if (conservation_check) then
    !!!   do k=1,4
    !!!      CALL mpp_sum(h2o(k))
    !!!      CALL mpp_sum(heat(k))
    !!!   end do
    !!!   if (mpp_pe()==mpp_root_pe()) then
    !!!      print *
    !!!      print '(a10,5a13)',   'ICE MODEL ','   AT START  ', &
    !!!           ' TOP FLUX DN.', &
    !!!           ' BOT FLUX DN.', &
    !!!           '   AT END    ', &
    !!!           '   ERROR     '
    !!!      print '(a10,5es13.5)','WATER     ', h2o , h2o (4)-(h2o (1)+h2o (2)-h2o (3))
    !!!      print '(a10,5es13.5)','HEAT      ', heat, heat(4)-(heat(1)+heat(2)-heat(3))
    !!!      print *
    !!!   end if
    !!!end if

    CALL sit_vdiff_restart()

    !--- release memory ------------------------------------------------  

    
    !!!! Interpolation array     
    !!!!nxj,                    !  nxj      : number of local longitudes (echam, =nproma=nbdim)
    !!!!!!  gauss grid description
    !!!!krow,                   !  krow     : sequential index (=1)
    !!!!!! - water body or ml_ocean variables
    !!!!echam, check GFS    DEALLOCATE (Sit%slm)             !  pslm      : sea land fraction (1. = land, 0. = sea/lakes)
    !!!!echam, check GFS    DEALLOCATE (Sit%seaice)          !  pseaice   : ice cover (fraction of 1-SLM) (0-1)                                                I/O
    !!!!echam, check GFS    DEALLOCATE (Sit%alake)           !  palake    : lake fraction of grid box                                                          
    !!!!echam, check GFS    DEALLOCATE (Sit%sni)             !  psni      : snow thickness over ice (m in water equivalent)                                    I/O
    !!!!echam, check GFS    DEALLOCATE (Sit%siced)           !  psiced    : ice thickness (m in water equivalent)                                              I/O
    !!!!echam, check GFS    DEALLOCATE (Sit%tsi)             !  ptsi      : surface temperature of ice (K)                                                     
    !!!!echam, check GFS    DEALLOCATE (Sit%tsw)             !  ptsw      : skin temperatrue (K) over water                                                    O     
    !!!!echam, check GFS    DEALLOCATE (Sit%ssw)             !  pssw      : skin salinity (PSU) over water                                                    O     
    !!!!echam, check GFS    DEALLOCATE (Sit%tsl)             !  ptsl      : calcuated Earth's skin temperature from vdiff/tsurf at t+dt (K)                    I/O
    !!!!echam, check GFS    DEALLOCATE (Sit%tslm)            !  ptslm     : calcuated Earth's skin temperature from vdiff/tsurf at t (K)                       I/O
    !!!!echam, check GFS    DEALLOCATE (Sit%tslm1)           !  ptslm1    : calcuated Earth's skin temperature from vdiff/tsurf at t-dt (K)                    I/O
    !!!!echam, check GFS    DEALLOCATE (Sit%fluxw)           !  pfluxw    : net surface energy flux over open water per water fraction                         I
    !!!!echam, check GFS                                 !              (icesheet+openwater) (w/m2) (positive upward)                                      
    !!!!echam, check GFS                                 !              = -(net solar + net longwave - sensible heat - latent heat)                        
    !!!!echam, check GFS    DEALLOCATE (Sit%dfluxs)          !  pdfluxs   : dG/dT (W/m2/K) (positive upward)                                                   I
    !!!!echam, check GFS    DEALLOCATE (Sit%soflw)           !  psoflw    : net SW flux over open water per water fraction (w/m2) (positive downward)          I         
    !!!!echam, check GFS    DEALLOCATE (Sit%fluxi)           !  pfluxi    : net surface energy flux over icesheet per water fraction                           I
    !!!!echam, check GFS                                 !              (icesheet+openwater) (w/m2) (positive upward)
    !!!!echam, check GFS                                 !              = -(net solar + net longwave - sensible heat - latent heat)
    !!!!echam, check GFS    DEALLOCATE (Sit%sofli)           !  psofli    : net SW flux over over icesheet per water fraction (w/m2) (positive downward)       I
    !!!!echam, check GFS      
    !!!!echam, check GFS!!! 2-d SIT vars 
    !!!!echam, check GFS    DEALLOCATE (Sit%disch)           !  pdisch   : surface runoff into ocean (m/s)                                                     I
    !!!!echam, check GFS    DEALLOCATE (Sit%ustrw)           !  taucx    : u-stress (Pa) over water at current time step                                       I/O
    !!!!echam, check GFS                                 !              (set to zero if snow/ice on top)
    !!!!echam, check GFS    DEALLOCATE (Sit%vstrw)           !  taucy    : v-stress (Pa) over water at current time step                                       I
    !!!!echam, check GFS                                 !              (set to zero if snow/ice on top)
    !!!!echam, check GFS    DEALLOCATE (Sit%temp2)           !  ptemp2    : 2 m temperature (K) at current time step                                           I
    !!!!echam, check GFS    DEALLOCATE (Sit%wind10w)         !  pwind10w     : 10m windspeed over water (m/s)                                                  I
    !!!!echam, check GFS    DEALLOCATE (Sit%ocu)             !  pocu      : ocean eastw. velocity (m/s)                                                        O
    !!!!echam, check GFS    DEALLOCATE (Sit%ocv)             !  pocv      : ocean northw. velocity (m/s)                                                       O
    !!!                         
    !!!!!! 2-d SIT vars (sit variable, allocate as 2D) 
    !!!DEALLOCATE (Sit%obsseaice)       !  pobsseaice  : observed sea ice fraction (fraction)                                             I
    !!!DEALLOCATE (Sit%obstsw)          !  pobswtb     : observed bulk sea surface temperature (K)                                        I
    !!!DEALLOCATE (Sit%obswsb)          !  pobswsb    : observed salinity (PSU, 0/00)                                                     I
    !!!DEALLOCATE (Sit%sitmask)         !  psitmask : mask for sit(1=.TRUE., 0=.FALSE.                                                    I
    !!!DEALLOCATE (Sit%bathy)           !  pbathy  : bathymeter (topography or orography) of ocean (m)                                    I  
    !!!DEALLOCATE (Sit%ctfreez2)        !  pctfreez2 : ref water freezing temperature (K)                                                 I
    !!!DEALLOCATE (Sit%wlvl)         !  pwlvl  : current water level (ice/water interface) a water body grid                           I/O
    !!!DEALLOCATE (Sit%ocnmask)         !  pocnmask : fractional mask for 3-D ocean grid (0-1)                                            I
    !!!DEALLOCATE (Sit%obox_mask)       !  obox_mask : 3-D ocean nudging mask, =0: nudging, = 1 (>0): nudging                             I
    !!!DEALLOCATE (Sit%wtb)          !  pwtb: bulk water temperature (K)                                                               O
    !!!DEALLOCATE (Sit%wub)          !  pwub: bulk water u current (m/s)                                                               O
    !!!DEALLOCATE (Sit%wvb)          !  pwvb: bulk water v current (m/s)                                                               O
    !!!DEALLOCATE (Sit%wsb)          !  pwsb: bulk water salinity (PSU)                                                                O
    !!!DEALLOCATE (Sit%fluxiw)          !  pfluxiw: over-water net surface heat flux (W/m2, + upward, i.e., from ocean))                  O
    !!!DEALLOCATE (Sit%pme2)            !  ppme2: net fresh water into ocean (P-E+ice_corr) (m/s, + downward)                             O
    !!!DEALLOCATE (Sit%subfluxw)        !  psubfluxw: subsurface ocean heat flux (W/m2, + upward)                                         O
    !!!DEALLOCATE (Sit%wsubsal)         !  pwsubsal: subsurface ocean salinity flux (m*PSU/s, + upward)                                   O
    !!!DEALLOCATE (Sit%cc)              !  pcc: cold content per water fraction (ice sheet+openwater) (J/m2)                              
    !!!                                !    (energy need to melt snow and ice, i.e., energy below liquid water at tmelt)                 I/O
    !!!DEALLOCATE (Sit%hc)              !  phc: heat content per water fraction (ice sheet+openwater) (J/m2)                              
    !!!                                !    (energy of a water column above tmelt)                                                       I/O
    !!!DEALLOCATE (Sit%engwac)          !  pengwac: accumulated energy per water fraction (ice sheet+openwater) (+ downward, J/m2)
    !!!                             !           (pfluxw+pfluxi+rain/snow advected energy in respect to liquid water at
    !!!                             !           tmelt)*dt                                                                             I/O
    !!!                             !!!  pengw: mean net surface heat flux per water fraction (ice sheet+openwater)
    !!!                             !!!    (+ downward, W/m2)                                                                         I/O
    !!!                             !!!  pengw2: mean snow corrected net surface heat flux per water fraction
    !!!                             !!!    (ice sheet+openwater) (+ downward, W/m2)
    !!!                             !!!    (pfluxw+pfluxi+rain/snow advected heat flux in respect to liquid water at
    !!!                             !!!     tmelt)                                                                                    I/O
    !!!!!!              engw(:,:)engw2(:,:)                              &
    !!!DEALLOCATE (Sit%sc)              !  psc: salinity content per water fraction (ice sheet+openwater) (PSU*m)                         I/O
    !!!DEALLOCATE (Sit%saltwac)         !  psaltwac: accumulated salt into water fraction (+ downward, PSU*m)                             I/O
    !!!
    !!!! implicit with vdiff
    !!!!echam, check GFS    DEALLOCATE (Sit%slm)             !  pslm: land fraction [0-1]                                                                      I
    !!!!echam, check GFS    DEALLOCATE (Sit%grndcapc)        !  pgrndcapc: areal heat capacity of the uppermost sit layer (snow/ice/water)
    !!!!echam, check GFS                                    !       (J/m**2/K)                                                                                I/O
    !!!!echam, check GFS    DEALLOCATE (Sit%grndhflx)        !  pgrndhflx: ground heat flux below the surface (W/m**2)
    !!!!echam, check GFS                                    !          (+ upward, into the skin layer)                                                        I/O
    !!!!echam, check GFS    DEALLOCATE (Sit%grndflux)        !  pgrndflx:  acc. ground heat flux below the surface (W/m**2*s)                                  
    !!!                                    !          (+ upward, into the skin layer)                                                        I/O
    !!!               
    !!!! - 2D from mo_memory_g3b (sit variables, allocate as 3D)
    !!!DEALLOCATE (Sit%obswt)         !  pobswt: observed potentail water temperature (K)                                          I/O
    !!!DEALLOCATE (Sit%obsws)         !  pobsws: observed salinity (PSU, 0/00)                                                     I/O
    !!!DEALLOCATE (Sit%obswu)         !  pobswu: observed u-current (m/s)                                                          I/O
    !!!DEALLOCATE (Sit%obswv)         !  pobswv: observed v-current (m/s)                                                          I/O           
    !!!DEALLOCATE (Sit%zsi)        !  pzsi  :                                                                                   
    !!!                           !        pzsi(jl,0): dry snow water equivalent (m)                                           I/O
    !!!                           !        pzsi(jl,1): dry ice water equivalent (m)                                            I/O
    !!!DEALLOCATE (Sit%silw)       !  psilw  :                                                                                  
    !!!                           !        psilw(jl,0): snow liquid water storage (m)                                          I/O
    !!!                           !        psilw(jl,1): ice liquid water storage (m)                                           I/O
    !!!DEALLOCATE (Sit%tsi)        !  ptsnic   :                                                                      
    !!!                           !        ptsnic(jl,0): snow skin temperatrue (jk)                                            I/O
    !!!                           !        ptsnic(jl,1): mean snow temperatrue (jk)                                            I/O
    !!!                           !        ptsnic(jl,2): ice skin temperatrue (jk)                                             I/O
    !!!                           !        ptsnic(jl,3): mean ice temperatrue (jk)                                             I/O
    !!!                                                                                                
    !!!DEALLOCATE (Sit%wt)         !  pwt      : potential water temperature (K)                                                I/O   
    !!!DEALLOCATE (Sit%wu)         !  pwu      : u current (m/s)                                                                I/O
    !!!DEALLOCATE (Sit%wv)         !  pwv      : v current (m/s)                                                                I/O
    !!!DEALLOCATE (Sit%ww)         !  pww      : w current (m/s)                                                                I/O
    !!!DEALLOCATE (Sit%ws)         !  pws      : practical salinity (0/00)                                                      I/O
    !!!DEALLOCATE (Sit%wtke)       !  pwtke    : turbulent kinetic energy (M2/S2)                                               I/O
    !!!DEALLOCATE (Sit%wlmx)       !  pwlmx    : mixing length (m)                                                              O
    !!!DEALLOCATE (Sit%wldisp)     !  pwldisp  : dissipation length (m)                                                         O
    !!!DEALLOCATE (Sit%wkm)        !  pwkm     : eddy diffusivity for momentum (m2/s)                                           O
    !!!DEALLOCATE (Sit%wkh)        !  pwkh     : eddy diffusivity for heat (m2/s)                                               O
    !!!DEALLOCATE (Sit%wrho)          !  pwrho1000: potential temperature at 1000 m (PSU)                                          O
    !!!DEALLOCATE (Sit%wtfn)          !  pwtfn  : nudging flux into sit-ocean at each level (W/m**2)                               I/O
    !!!DEALLOCATE (Sit%wsfn)          !  pwsfn  : nudging salinity flux into sit-ocean at each level (PSU*m/s)                     I/O
    !!!DEALLOCATE (Sit%wtfns)         !  pwtfns : nudging flux into sit-ocean for an entire column (W/m**2), i.e.
    !!!                           !        pwtfns=SUM(pwtfn(jl,:))                                                             I/O
    !!!DEALLOCATE (Sit%wsfns)         !  pwsfn  : nudging salinity flux into sit-ocean for an entire column (PSU*m/s),             
    !!!                           !        i.e., pwsfns=SUM(pwsfn(jl,:)                                                        I/O
    !!!DEALLOCATE (Sit%awufl)         !  pawufl : advected u flux at each level (positive into ocean) (m/s*m/s)                    I
    !!!DEALLOCATE (Sit%awvfl)         !  pawvfl : advected v flux at each level (positive into ocean) (m/s*m/s)                    I
    !!!DEALLOCATE (Sit%awtfl)         !  pawtfl : advected temperature flux at each level (positive into ocean) (W/m2)             I
    !!!DEALLOCATE (Sit%awsfl)         !  pawsfl : advected salinity flux at each level (positive into ocean) (PSU*m/s)             I
    !!!DEALLOCATE (Sit%awtkefl)       !  pawtkefl : advected tke at each level (positive into ocean) (m3/s3)                       I
    !!!                                                                                                                                  
    !!!!!! - variables internal to physics                                                                                                     
    !!!!!!echam, check GFS    DEALLOCATE (Sit%zevapw)        !  pevapw   : evaporation from water surface (kg/m2/s) (positive downward)                   I
    !!!!!!echam, check GFS    DEALLOCATE (Sit%rsfl)          !  prsfl    : large scale rain flux at the surface (kg/m2/s) (positive downward)             I
    !!!!!!echam, check GFS    DEALLOCATE (Sit%rsfc)          !  prsfc    : convective rain flux at the surface (kg/m2/s) (positive downward)              I
    !!!!!!echam, check GFS    DEALLOCATE (Sit%ssfl)          !  pssfl    : large scale snow flux at the surface (kg/m2/s) (positive downward)             I
    !!!!!!echam, check GFS    DEALLOCATE (Sit%ssfc)          !  pssfc    : convective snow flux at the surface (kg/m2/s) (positive downward)              I
    !!!!end sit_vdiff variable                                             
    !#######################################################################
    CONTAINS   
    !#######################################################################
    ! <SUBROUTINE NAME="sit_vdiff_restart">
    ! <DESCRIPTION>
    !  Write out restart files registered through register_restart_file
    ! </DESCRIPTION>
    !SUBROUTINE sit_vdiff_restart(Sit, time_stamp)
    SUBROUTINE sit_vdiff_restart( time_stamp)
      IMPLICIT NONE  
      !!! type (sit_vdiff_data_type), intent(inout), optional :: Sit
      character(len=*),         intent(in), optional :: time_stamp
    
      !!! CALL save_restart(Sit_restart, time_stamp)
     !CALL icebergs_save_restart(Ice%icebergs)
     ! This should go here but since "Ice" is not available we have to
     ! rely on the restart written via sit_vdiff_end() -AJA
    
    END SUBROUTINE sit_vdiff_restart
    ! </SUBROUTINE>
    !#######################################################################
  END SUBROUTINE sit_vdiff_end
  !#######################################################################

  SUBROUTINE init_convect_tables

    !!! USE sit_constants_mod, ONLY: alv, als, cpd, rd, rv, tmelt, &
    !!!                         c3les, c3ies, c4les, c4ies, c5les, c5ies

    real, PARAMETER :: zavl1 = -6096.9385
    real, PARAMETER :: zavl2 =    21.2409642
    real, PARAMETER :: zavl3 =    -2.711193
    real, PARAMETER :: zavl4 =     1.673952
    real, PARAMETER :: zavl5 =     2.433502 

    real, PARAMETER :: zavi1 = -6024.5282
    real, PARAMETER :: zavi2 =    29.32707
    real, PARAMETER :: zavi3 =     1.0613868
    real, PARAMETER :: zavi4 =    -1.3198825
    real, PARAMETER :: zavi5 =    -0.49382577        

    real :: z5alvcp, z5alscp, zalvdcp, zalsdcp
    real :: ztt, zldcp
    real :: zcvm3, zcvm4, zcvm5
    real :: zavm1, zavm2, zavm3, zavm4, zavm5

    INTEGER :: it

    z5alvcp = c5les*alv/cpd
    z5alscp = c5ies*als/cpd

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    DO it = jptlucu1, jptlucu2
      ztt = 0.001*it
      IF ((ztt-tmelt) > 0.0) THEN
        zcvm3 = c3les
        zcvm4 = c4les
        zcvm5 = z5alvcp
        zldcp = zalvdcp
        zavm1 = zavl1
        zavm2 = zavl2
        zavm3 = zavl3
        zavm4 = zavl4
        zavm5 = zavl5
      ELSE
        zcvm3 = c3ies
        zcvm4 = c4ies
        zcvm5 = z5alscp
        zldcp = zalsdcp
        zavm1 = zavi1
        zavm2 = zavi2
        zavm3 = zavi3
        zavm4 = zavi4
        zavm5 = zavi5
      END IF
      tlucuc(it)  = zldcp
      tlucua(it)  = EXP((zavm1/ztt+zavm2+zavm3*0.01*ztt+zavm4*ztt*ztt*1.e-5+zavm5*LOG(ztt)))*rd/rv
      tlucub(it)  = zcvm5*(1.0/(ztt-zcvm4))**2
      tlucuaw(it) = EXP((zavl1/ztt+zavl2+zavl3*0.01*ztt+zavl4*ztt*ztt*1.e-5+zavl5*LOG(ztt)))*rd/rv
    END DO
    
  END SUBROUTINE init_convect_tables

  SUBROUTINE lookuperror (name)

    !!! USE mo_exception,  ONLY: message, finish

    CHARACTER(len=*), INTENT(in) :: name

    ! normal run informs only
    ! CALL message (name, ' lookup table overflow')
    ! debug run, so stop at problem
    !!! CALL finish (name, ' lookup table overflow')
    WRITE (nerr,*) name, ' lookup table overflow'
    ! reset lookupoverflow for next test  

    lookupoverflow = .FALSE.

  END SUBROUTINE lookuperror

  !#######################################################################

  
SUBROUTINE update_snow_ice_property(swe,iwe,hsn,hesn,hice,heice,seaice)
  IMPLICIT NONE
  real, INTENT(IN):: swe    ! snow water equivalent of the grid (m)
  real, INTENT(IN):: iwe    ! iwe water equivalent of the grid (m)
  real, INTENT(OUT):: hsn,    & ! physical snow thickness of the grid (m)
                          hesn,   & ! effictive snow thickness of the grid (m)
                          hice,   & ! physical ice thickness of the grid (m)
                          heice,  & ! effective ice thickness of the grid (m)
                          seaice    ! sea ice fraction of the grid [0-1]
  hsn=HSNFN(swe)
  hesn=HEFN(hsn/4.,xksn,omegas)
  hice=HICEFN(iwe)
  heice=HEFN(hice/4.,xkice,omegas)
  seaice=SEAICEFN(iwe)
END SUBROUTINE update_snow_ice_property

! ----------------------------------------------------------------------
real FUNCTION HASTRFN(hstr)
!*********************
! CALCUL Effective thickness
! h0: physical thickness of the skin layer (m)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! HEFN: effective thickness of the skin layer (m)
!*********************
  IMPLICIT NONE
  real:: hstr
  hastrfn=SQRT( 1.-2.*COS(hstr)*EXP(-hstr)+EXP(-hstr)**2. )/SQRT(2.)
END FUNCTION HASTRFN
! ----------------------------------------------------------------------
real FUNCTION HEFN(h0,kh,omegas)
!*********************
! CALCUL Effective thickness
! h0: physical thickness of the skin layer (m)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! HEFN: effective thickness of the skin layer (m)
!*********************
  IMPLICIT NONE
  real:: kh,omegas,h0
  real:: h_ref,hstr,hastr
  h_ref=SQRT(2.*kh/omegas)
  hstr=h0/h_ref
  hastr=HASTRFN(hstr)
  HEFN=h_ref*hastr
END FUNCTION HEFN
! ----------------------------------------------------------------------
real FUNCTION TASTRFN(hstr)
!*********************
! CALCUL Effective thickness
! h0: physical thickness of the skin layer (m)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! HEFN: effective thickness of the skin layer (m)
!*********************
  USE mod_eos_ocean,    ONLY: api
  IMPLICIT NONE
  real:: hstr
!!!  TASTRFN=api/4.-ARCTAN( SIN(hstr)*EXP(-hstr)/(1.-Cos(hstr)*EXP(-hstr)) )
  TASTRFN=api/4.-ATAN( SIN(hstr)/(EXP(hstr)-Cos(hstr)) )
END FUNCTION TASTRFN
! ----------------------------------------------------------------------
real FUNCTION SSTRFN(zkm1str,zkstr,zkp1str)
!*********************
! CALCUL dimensionless elasticity of numerical layer k
! zkm1,zk,zkp1: vertical coordinates of k-1, k, k+1 (m) (positive upward)
! s: elasticity (dGk/dTk-dGk+1/dTk)
! sstr: dimensionless elasticity
! sstr=s/rhogcg/SQRT(kh*omegas)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
!*********************
  IMPLICIT NONE
  real, INTENT(in):: zkm1str,zkstr,zkp1str
  SSTRFN=(1./(zkm1str-zkstr)+1./(zkstr-zkp1str))/SQRT(2.)
END FUNCTION SSTRFN
! ----------------------------------------------------------------------
real FUNCTION dGdTFN(em,ra,rc,T0,ps,rhoa)
!*********************
! CALCUL dG/dT of land surface (W/m2/K)
! em: emissivity of land surface
! ra: aerodynamic resistsnce (s/m)
! rc: canopy resistance (s/m)
! T0: land skin temperature (K)
! ps: surface pressure (Pa)
! rhoa: air density (kg/m3)
! dGdT: 
!*********************
  USE mod_eos_ocean,      ONLY: stbo,cpd,alv
!!!  USE convect_tables_mod,   ONLY: jptlucu1,jptlucu2,tlucub  
!!!  USE convect_tables_mod, ONLY : lookuperror, lookupoverflow, jptlucu1    &
!!!                              , jptlucu2, tlucua, tlucub, tlucuaw
  IMPLICIT NONE
  real, INTENT(in):: em,ra,rc,T0,ps,rhoa
  INTEGER  :: it
  real:: dqsatdT
  it = MAX(MIN(NINT(T0*1000.),jptlucu2),jptlucu1)
  dqsatdT=0.622*tlucub(it)/ps
  dGdTFN=4.*em*stbo*T0**3+rhoa*cpd/ra+rhoa*alv/(ra+rc)*dqsatdT
END FUNCTION dGdTFN
! ----------------------------------------------------------------------
real FUNCTION SSTR0FN(rhogcg,kh,omegas,dGdT,z0str,z1str)
!*********************
! CALCUL dimensionless elasticity of numerical layer k
! rhogcg: volume heat capacity (kg/m3*J/kg/K)=(J/m3/K)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! dGdT: 
! z0str,z1str: dimensionless vertical coordinates of k=0, k=1 (positive upward)
! s: elasticity (dGk/dTk-dGk+1/dTk)
! sstr: dimensionless elasticity
! sstr=s/rhogcg/SQRT(kh*omegas)
!*********************
  IMPLICIT NONE
  real, INTENT(in):: rhogcg,kh,omegas,dGdT,z0str,z1str
  SSTR0FN=dGdT/rhogcg/SQRT(kh*omegas)+1./(z0str-z1str)/SQRT(2.)
END FUNCTION SSTR0FN
! ----------------------------------------------------------------------
real FUNCTION HEPSTRFN(ha,ht,s,ta)
  IMPLICIT NONE
  real:: ha,ht,s,ta
  HEPSTRFN=( 2.*ha**2-EXP(-ht)**2*s**2                                                  &
      +SQRT(4.*ha**4+4.*Cos(2.*(ta-ht))*EXP(-ht)**2*s**2+EXP(-ht)**4*s**4)        &
    )/(4.*Cos(ta-ht)*EXP(-ht)*ha) 
END FUNCTION HEPSTRFN
! ----------------------------------------------------------------------
real FUNCTION HEPFN(h,ht,s,kh,omegas,rhogcg)
!*********************
! CALCUL Effective thickness
! h: physical thickness of the numerial layer (m)
! ht: center of temperature below upper interface (m)
! s: elasticity (dGk/dTk-dGk+1/dTk)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! rhogcg: volume heat capacity (kg/m3*J/kg/K)=(J/m3/K)
! HEPFN: effective thickness of the numerical layer (m)
!*********************
  IMPLICIT NONE
  real, INTENT(in):: h,ht,s,kh,omegas,rhogcg
  real:: h_ref,hstr,htstr,hastr,tastr,sstr,hepstr
  h_ref=SQRT(2.*kh/omegas)
  hstr=h/h_ref
  htstr=ht/h_ref
  hastr=HASTRFN(hstr)
  tastr=TASTRFN(hstr)
  sstr=s/rhogcg/SQRT(kh*omegas)
  hepstr=HEPSTRFN(hastr,htstr,sstr,tastr)
  HEPFN=h_ref*hepstr
END FUNCTION HEPFN
! ----------------------------------------------------------------------
real FUNCTION HE0PFN(h,s,kh,omegas,rhogcg)
!*********************
! CALCUL Effective thickness of the surface numerical layer
! h: physical thickness of the numerial layer (m)
! s: elasticity (dGk/dTk-dGk+1/dTk)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! rhogcg: volume heat capacity (kg/m3*J/kg/K)=(J/m3/K)
! HE0PFN: optimal effective thickness of the skin layer (m)
!*********************
  IMPLICIT NONE
  real, INTENT(in):: h,s,kh,omegas,rhogcg
  HE0PFN=HEPFN(h,0.,s,kh,omegas,rhogcg)
END FUNCTION HE0PFN
! ----------------------------------------------------------------------
real FUNCTION SEAICEFN(siced)
!*********************
! CALCUL sea ice cover fraction of a grid basen on mean sea ice depth
! siced: mean sea ice depth of a grid (m, swe)
! SEAICEFN: sea ice cover fraction [0-1] (dimensionless)
! assume to be lognormal distribution
! seaice[0.]=0.
! seaice[csiced]=0.5
! seaice[2*csiced]=0.84
!*********************
  IMPLICIT NONE
! threshold sea ice depth, (= 2 m)
! this value should be resolution dependent
! 0 to be infinity fine resolution
! seaice[csiced]=50%
!  real:: csiced      ! seaice[csiced]=0.5

  real, INTENT(in):: siced
  IF (csiced.GT.TOL) THEN
#ifdef __ibm__
    SEAICEFN=0.5*(ERF(LOG(siced/csiced))+1.)
#else            
    SEAICEFN=0.5*(ERF(LOG(siced/csiced))+1.)
    !!! SEAICEFN=0.5_dp*(DERF(LOG(siced/csiced))+1.)
#endif
  ELSEIF (.FALSE.) THEN
    ! csiced=3.      !    
    ! csiced=2.      ! seaice fraction is estimated to be 2.8%, while the observation 3.7% (T31)
    ! csiced=0.5_dp     ! seaice fraction is estimated to be 6.8%, while the observation 3.7% m (v9.9003, T31)                      
    ! csiced=1.0_dp     ! seaice fraction is estimated to be still as high as 6.8%, while the observation 3.7% m (v9.9004, T31, T63)
      csiced=2.0
    !!!   IF (GDCHK2) then
    !!!     WRITE(nerr,*) "csiced: Truncation is not tested in T',nn,' runs."
    !!!   ENDIF
#ifdef __ibm__
    SEAICEFN=0.5*(ERF(LOG(siced/csiced))+1.)
#else            
    !!! SEAICEFN=0.5_dp*(DERF(LOG(siced/csiced))+1.)
    SEAICEFN=0.5*(ERF(LOG(siced/csiced))+1.)
#endif
  ELSE
    IF (siced.GT.0.) THEN  
    ! seaice mask: Winner wins!
    ! seaice fraction is estimated to be 5.6%, while the observation 3.7% (T31)(cob10dnnnn, cob10d10dnn, cob10d10d10d)
    ! there is 660 ice grids, while the observation is 760 ice grid.    
      SEAICEFN=1.
    ELSE
      SEAICEFN=0.
    ENDIF
  ENDIF
END FUNCTION SEAICEFN
! ----------------------------------------------------------------------
real FUNCTION SICEDFN(seaice)
!*********************
! CALCUL sea ice cover fraction of a grid basen on mean sea ice depth
! siced: mean sea ice depth of a grid (m, swe)
! SEAICEFN: sea ice cover fraction [0-1] (dimensionless)
! assume to be lognormal distribution
! seaice[0.]=0.
! seaice[csiced]=0.5
! seaice[2*csiced]=0.84
!*********************
  IMPLICIT NONE
! threshold sea ice depth, (= 2 m)
! this value should be resolution dependent
! 0 to be infinity fine resolution
! seaice[csiced]=50%
!!!  real:: csiced      ! seaice[csiced]=0.5

  real, INTENT(in):: seaice
!!!  IF (csiced.GT.TOL) THEN
!!!#ifdef __ibm__
!!!    !!! SICEDFN=csiced*EXP(ERFINV(2.*seaice-1.))
!!!    SICEDFN=csiced*EXP(ERFINV(2.*seaice-1.))
!!!#else            
!!!    !!! SICEDFN=csiced*EXP(DERFINV(2.*seaice-1.))
!!!    SICEDFN=csiced*EXP(ERFINV(2.*seaice-1.))
!!!#endif
!!!  ELSE
    IF (seaice.GT.0.5) THEN  
    !!! IF (seaice.GT.TOL) THEN  
    ! seaice mask: Winner wins!
    ! seaice fraction is estimated to be 5.6%, while the observation 3.7% (T31)(cob10dnnnn, cob10d10dnn, cob10d10d10d)
    ! there is 660 ice grids, while the observation is 760 ice grid.    
      SICEDFN=1.
    ELSE
      SICEDFN=0.
    ENDIF
!!!  ENDIF
END FUNCTION SICEDFN
! ----------------------------------------------------------------------
  real FUNCTION HSNFN(sni)
  !*********************
  ! CALCU snow thickness of the grid (m), based on grid mean snow depth in swe (sni)
  !*********************
    IMPLICIT NONE
    real, INTENT(in):: sni  
    HSNFN=sni*rhoh2o/rhosn
  END FUNCTION HSNFN
  !---------------------
  real FUNCTION HICEFN(siced)
  !*********************
  ! CALCUL mean ice thickness of the grid (m), based on grid mean sea ice depth in swe (siced)
  !*********************
    IMPLICIT NONE
    real, INTENT(in):: siced  
    HICEFN=siced*rhoh2o/rhoice
  END FUNCTION HICEFN
! ----------------------------------------------------------------------
  real FUNCTION wt_maxden(depth,obswtb,obswsb,mixing_depth)
  !*********************
  ! CALCUL water temperature based on max density water temperature approach
  ! set initial profile to be expontential decay to tmaxden
  !   =  3.73 C for fresh water s=0
  !   = -4.35 C for ocena water s=36.3 0/00.
  !*********************
    IMPLICIT NONE
    real, INTENT(in):: depth           ! depth (m, positive downward)
    real, INTENT(in):: obswtb          ! observed bulk water temperature (K)
    real, INTENT(in):: obswsb          ! observed bulk salinity (PSU)
    real, INTENT(in):: mixing_depth    ! mixing depth
!ps    IF ( depth.LE.mixing_depth) then
!ps      wt_maxden=obswtb
!ps    ELSE
!ps      wt_maxden=tmaxden(obswsb)+(obswtb-tmaxden(obswsb))*EXP(-(depth-mixing_depth)/100.)
!ps      wt_maxden=tmaxden(obswsb)+(obswtb-tmaxden(obswsb))*EXP(-(depth-mixing_depth)/mixing_depth)
      wt_maxden=tmaxden(obswsb)+(obswtb-tmaxden(obswsb))*EXP(-depth/mixing_depth)
!ps    ENDIF  
  END FUNCTION wt_maxden
! ----------------------------------------------------------------------        
  SUBROUTINE pzcord(zwlvl,bathy,nls,nle,z,zlk,hw,lsoil,lshlw)
!
!     WATER BODY POINT VERTICAL LEVELS
!     zlk: coordinates of a water body flux (m ASL).
!     z: z coordinates of the water body temperature (m ASL).
!     nls: starting water level of the water body point.
!       z(k)=zwlvl, IF k is <= nls. 
!     nle: maximum water levels of the water body point.
!        Note nle=-1 IF water level is below point water body bed.
!        Soil layer always is sitwt(nle+1). If there is no water,
!        sitwt(nle+1)=sitwt(0,jrow). In addition, z(k)=bathy IF k >= nle.
!     hw: thickness of each layer
!     lshlw = shallow water mode (due to evaportion or freezing)
!       (one layer water body)
!     lsoil = soil grid
!     zwlvl  : current water level (ice/water interface) a water body grid            I
!----------------------------------------------------------
!*    0 Locate Space
      IMPLICIT NONE
!     0.1 Calling Variables
      real, INTENT(IN):: zwlvl,bathy
      
!     input
!v77      INTEGER oceanid,jl
!     output
      INTEGER, INTENT(OUT):: nls,nle
      real, INTENT(OUT):: z(0:lkvl+1),zlk(0:lkvl+1),hw(0:lkvl)
      LOGICAL, INTENT(OUT):: lsoil, lshlw
!     0.2 local variables
      INTEGER k
      REAL :: tmp
      real:: zdepth(0:lkvl)
      zdepth(0:lkvl)=sit_zdepth(0:lkvl)
!     0.3 Initial	(default)
      nls=-999
      nle=-999
      z=xmissing
      zlk=xmissing
      hw=xmissing
      lshlw = .FALSE.
      lsoil = .FALSE.
!
!*    1.0 Determine Ocean Coordinate
!
!
      IF ( (zwlvl.GE.bathy+wcri) )  THEN
!
!     2. Water 
!
        nls=0
        nle=lkvl
        z(lkvl+1)=bathy       
        DO k = 0,lkvl
          tmp=zwlvl-zdepth(k)
          IF ( tmp.GT.(bathy+wcri) ) THEN
            z(k)=tmp
            nle=k
          ELSE
            z(k)=bathy       
          ENDIF
        ENDDO
!
!*    3. Determine Thickness (hw) and zlk of Each Layer
!
        zlk(0)=zwlvl
        DO k = 1,nls
          !! vanisih top layers when waterlevel is lower than the reference level
          zlk(k)=zwlvl
        ENDDO
        DO k = nls+1,lkvl+1
!!!          zlk(k)=(z(k-1)+z(k))/2.
          zlk(k)=zwlvl-sit_fluxdepth(k)
        ENDDO
                
        hw(0)=zlk(0)-zlk(nls+1)
        DO k = 1,nls
          hw(k)=0.
        ENDDO
        hw(nls+1)=zlk(0)-zlk(nls+2)
        DO k = nls+2,nle-1
          hw(k)=zlk(k)-zlk(k+1)
        ENDDO
        hw(nle)=zlk(nle)-z(nle+1)

!
!*    4. Determine Current Existence of Water
!
        IF ((nle-nls).EQ.1) THEN
          lshlw = .TRUE.
!         shallow water body
        ENDIF
!
      ELSE
!
!     5.0 Not a Ocean Grid (Soil)
!
        nle=-1
        z(nle+1)=bathy 
        lsoil = .TRUE.
      ENDIF
   IF(lwarning_msg.GE.4) then
      WRITE(nerr,*) "lshlw=",lshlw,"lsoil=",lsoil,"nls=",nls,"nle=",nle
   ENDIF
      
  END SUBROUTINE pzcord
! **********************************************************************

  SUBROUTINE cal_ratioBlending(pmyrank,pii,pjj,plat,inlon,outratio)

      use mod_sit_control,  only: sit_domain_w,sit_domain_e &
                               ,sit_domain_s,sit_domain_n &
                               ,sit_domain_extgrd
! location blending
      real:: plat,inlon
      real:: plon
      real:: outratio
      real:: rrn,rrs,rre,rrw,rr_sn,rr_we,rrtmp
!    real, PARAMETER:: rf=10._dp
      real, PARAMETER:: rf=5.
      integer:: pmyrank,pii,pjj


      IF(inlon .LT. 0.) THEN
        plon=inlon + 360.
      else
        plon=inlon
      ENDIF

      IF((plon .GE. sit_domain_w) .AND. &
         (plon .LE. sit_domain_e) ) then      
        rre=1.
        rrw=1.
      ELSEIF( (plon .GT. sit_domain_e) .AND. &
              (plon .LE. sit_domain_e+sit_domain_extgrd) ) THEN
        rre=min(sit_domain_e+sit_domain_extgrd,360.)
        rrw=sit_domain_e
      ELSEIF( (plon .GT. sit_domain_w-sit_domain_extgrd) .AND. &
              (plon .LE. sit_domain_w) ) THEN
        rre=max(sit_domain_w-sit_domain_extgrd,0.)
        rrw=sit_domain_w
      ELSE
        rre=0.
        rrw=0.
      ENDIF

      IF((plat .GE. sit_domain_s) .AND. &
         (plat .LE. sit_domain_n) ) then      
        rrn=1.
        rrs=1.
      ELSEIF( (plat .GT. sit_domain_n) .AND. &
              (plat .LE. sit_domain_n+sit_domain_extgrd)  ) THEN
        rrn=min(sit_domain_n+sit_domain_extgrd,90.)
        rrs=sit_domain_n
      ELSEIF( (plat .GE. sit_domain_s-sit_domain_extgrd) .AND. &
              (plat .LT. sit_domain_s) ) THEN
        rrn=max(sit_domain_s-sit_domain_extgrd,-90.)
        rrs=sit_domain_s
      ELSE
        rrn=0.
        rrs=0.
      ENDIF

      IF(rrw .EQ. 1.) THEN
        rr_we=1.
      ELSEIF(rrw .EQ. 0.) THEN
        rr_we=0.
      ELSE
        rr_we=1.-ABS((plon-rrw)/(rre-rrw))
      ENDIF

      IF(rrs .EQ. 1.) THEN
        rr_sn=1.
      ELSEIF(rrs .EQ. 0.) THEN
        rr_sn=0.
      ELSE
        rr_sn=1.-ABS((plat-rrs)/(rrn-rrs))
      ENDIF

      if(rr_sn .eq. 1.) then
        rrtmp=rr_we
      elseif(rr_we .eq. 1.)then
        rrtmp=rr_sn
      else
        rrtmp=0.5*(rr_sn+rr_we)
      endif

      IF(rrtmp .EQ. 1.) THEN
        outratio=1.
      ELSE if(rrtmp .eq. 0.) then
        outratio=0.
      ELSE
        outratio=exp(-rf/rrtmp*exp(1./(rrtmp-1.)))
      ENDIF

      if(pmyrank.eq.6 .AND. pii.eq.1 .AND. pjj.eq.12) then
        print *,'in ratio: myrank=',pmyrank,',ii=',pii &
               ,',jj=',pjj,',rre=',rre,',rrw=',rrw    &
               ,',rr_we=',rr_we,',rrn=',rrn,',rrs=',rrs &
               ,',rr_sn=',rr_sn,',rrtmp=',rrtmp      &
               ,',outratio=',outratio
      endif


  END SUBROUTINE cal_ratioBlending 

END MODULE mod_sit_vdiff


