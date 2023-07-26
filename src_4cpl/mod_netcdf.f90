#define LCWBGFS T
#if defined(LCWBGFS)
#define mpp_pe() myrank
#define p_pe myrank
#define p_io 0
#define p_parallel_io (myrank .eq. 0)
#define p_parallel_ocean (myrank .eq. 0)
#endif

      MODULE mo_netcdf
	  
        USE mpe
        USE const,           ONLY:ifilin_nc,do_sit
        USE mod_sit_control,  ONLY:nocn_sv,zocn_option,ocn_tlz,ocn_k1 &
                                  ,lssst,lwarning_msg

        IMPLICIT NONE
	  
        PRIVATE                      ! used in
        PUBLIC :: io_get_varindx     ! mo_memory_base
        PUBLIC :: add_dim            !
        PUBLIC :: add_unknown_dim    !
        PUBLIC :: max_dim_name       ! mo_memory_gl
        PUBLIC :: FILE_INFO          ! mo_io
        PUBLIC :: nf_noerr           ! mo_io
        PUBLIC :: nf_enotnc          ! mo_io
        PUBLIC :: nf_nowrite         ! mo_io
        PUBLIC :: nf_write           ! mo_io
        PUBLIC :: nf_nofill          ! mo_io
        PUBLIC :: io_ndim_ids        ! mo_io
        PUBLIC :: io_dim_ids         ! mo_io
        PUBLIC :: io_dim             ! io_dim_ids
        PUBLIC :: nf_max_name        ! mo_io
        PUBLIC :: nf_inq_varid       ! mo_io
        PUBLIC :: nf_inq_dimid       ! mo_io
        PUBLIC :: nf_global          ! mo_io
        PUBLIC :: nf_double          ! mo_io
        PUBLIC :: nf_create          ! mo_io
        PUBLIC :: nf__create         ! mo_io, version allows setting I/O buffer size
        PUBLIC :: nf_set_fill        ! mo_io
        PUBLIC :: nf_strerror        ! mo_io
        PUBLIC :: nf_open            ! mo_io
        PUBLIC :: nf__open           ! mo_io, version allows setting I/O buffer size
        PUBLIC :: nf_clobber         ! mo_io
        PUBLIC :: nf_close           ! mo_io
#if defined (HAVE_LIBNETCDF64)
        PUBLIC :: nf_64bit_offset    ! mo_io
#endif
        PUBLIC :: io_inq_dimid       ! mo_o3clim
        PUBLIC :: io_inq_dimlen      ! mo_o3clim
        PUBLIC :: io_inq_varid       ! mo_o3clim
        PUBLIC :: io_get_var_double  ! mo_o3clim
        PUBLIC :: io_get_var_double1
        PUBLIC :: io_get_vara_double ! mo_o3clim
        PUBLIC :: io_info_print      ! mo_sst
        PUBLIC :: message_text       ! mo_sst
!        PUBLIC :: io_init_dims       ! inictl
        PUBLIC :: io_get_att_int     ! mo_io
        PUBLIC :: io_get_att_double  ! mo_io
        PUBLIC :: io_get_att_text    ! mo_io
        PUBLIC :: io_def_var         ! mo_io
        PUBLIC :: io_put_var_double  ! mo_io
        PUBLIC :: io_put_var_double1
        PUBLIC :: io_put_att_text    ! mo_io
        PUBLIC :: io_put_att_int     ! mo_io
        PUBLIC :: io_put_att_double  ! mo_io
        PUBLIC :: io_enddef          ! mo_io
        PUBLIC :: io_def_dim         ! mo_io
        PUBLIC :: io_info_construct  ! mo_io
        PUBLIC :: t_att_text         ! text attributes data type
        PUBLIC :: put_att            ! store attribute in attributes data type
        PUBLIC :: global_att         ! global text attributes
        
        PUBLIC :: initialsize        ! initial size of netCDF file
        PUBLIC :: chunksize          ! preferred netCDF I/O buffer size 
        
        PUBLIC :: cleanup_netcdf     ! deallocate module variables

        INCLUDE 'netcdf.inc'

        PUBLIC :: nfnlvl,lkvl
        PUBLIC :: diecast_zdepth,diecast_fluxdepth,sit_zdepth,sit_fluxdepth,ocn_z
        PUBLIC :: set_ocndepth
!ps
        PUBLIC :: io_open
        PUBLIC :: io_open_unit
        PUBLIC :: io_close
        
        !
        INTEGER, PARAMETER :: IO_READ = 1, IO_WRITE = 2
        CHARACTER(512) :: message_text = ''
        INTEGER, PARAMETER :: nerr = 6     ! error output stream
        LOGICAL, PARAMETER :: ldebugio =.FALSE.
        INTEGER, PARAMETER :: NETCDF64=4
        INTEGER, PARAMETER :: NETCDF=2
        INTEGER:: rerun_filetype= NETCDF
        PUBLIC:: NETCDF,NETCDF64
!ps		
#if defined (__SX__) || defined (ES) || defined (__uxp__) || defined (__PGI)
  INTEGER :: initialsize = 33554432      ! that's 32 MByte   
  INTEGER :: chunksize   = 33554432      ! too
#else
  INTEGER :: initialsize =    32768      ! that's 32 kByte   
  INTEGER :: chunksize   =    32768      ! too
#endif

!------------------------------------------------------------------------------
  TYPE FILE_INFO

    LOGICAL :: opened                              ! open = .true. or closed = .FALSE.

    INTEGER :: file_id                             ! netCDF file id 
    INTEGER :: access_mode                         ! access mode for that file
    INTEGER :: ncdims(NF_MAX_VAR_DIMS) 
    INTEGER :: format                              ! file format NETCDF

    CHARACTER (NF_MAX_NAME) :: creation_program    ! name of this program
    CHARACTER (NF_MAX_NAME) :: creation_user       ! who has run this program
    CHARACTER (NF_MAX_NAME) :: creation_date       ! date of creation of netCDF file
    CHARACTER (NF_MAX_NAME) :: binary_source       ! binary data type of src (CRAY/IEEE)
    CHARACTER (NF_MAX_NAME) :: file_type           ! initital or restart file ...
    CHARACTER (NF_MAX_NAME) :: file_name           ! nc file name
    CHARACTER (NF_MAX_NAME) :: title
  END TYPE FILE_INFO
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: max_dim_name = 32
!------------------------------------------------------------------------------
  TYPE IO_dim
    INTEGER                      :: dim_id   =   0      ! temporary NetCDF id
    INTEGER                      :: var_id   =   0      ! temporary NetCDF id
    INTEGER                      :: dim_len  =  -1
    CHARACTER (len=max_dim_name) :: dim_name =  ''
    CHARACTER (len=64)           :: longname =  ''
    CHARACTER (len=32)           :: units    =  ''
    INTEGER                      :: levtyp   =   0      ! GRIB level type
    LOGICAL                      :: single   =  .FALSE. ! single level
    REAL  ,POINTER           :: value(:) => NULL()
  END TYPE IO_dim
!------------------------------------------------------------------------------
  INTEGER       ,PARAMETER    :: max_dim_ids = 50
  INTEGER                     :: IO_ndim_ids
  TYPE (IO_dim) ,TARGET ,SAVE :: IO_dim_ids (max_dim_ids)

  INTEGER :: nfnlvl        ! =12, number of fine layers in sit for the uppermost ocean layer
  INTEGER :: lkvl          ! =49, =ocn_k1+nfnlvl-3, i.e., lkvl=49, ocn_k1-2=lkvl+1-nfnlvl
                           ! number of vertical levels of sit ocean model 
!!!  INTEGER :: ocn_k1        != 40, number of ocean ocean levels
  REAL, POINTER  :: ocn_z(:)
  REAL, POINTER  :: diecast_zdepth(:),diecast_fluxdepth(:)     ! unit: m
  REAL, POINTER  :: sit_zdepth(:),sit_fluxdepth(:)             ! unit: m
  REAL, POINTER  :: sit_gribzdepth(:),sit_gribfluxdepth(:)     ! sit_gribzdepth (for grib format)! bjt

!------------------------------------------------------------------------------
! data type to hold file attributes
!----------------------------------
  TYPE t_att_text
    CHARACTER (len= 32) :: name = ''
    CHARACTER (len=128) :: text = ''
  END TYPE t_att_text

  TYPE (t_att_text) ,SAVE :: global_att (20)
!------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE put_att (att, name, value)
  TYPE (t_att_text) ,INTENT(inout) :: att(:)
  CHARACTER (len=*) ,INTENT(in)    :: name
  CHARACTER (len=*) ,INTENT(in)    :: value
  !----------------------------------------------
  ! store an attribute in the attribute data type
  !----------------------------------------------
    INTEGER :: i
    DO i=1,SIZE(att)
      IF (att(i)% name == '' .OR. att(i)% name == name) THEN
        att(i)% name = name
        att(i)% text = value
        EXIT
      ENDIF
    END DO
  END SUBROUTINE put_att
!------------------------------------------------------------------------------
  SUBROUTINE IO_info_construct(fileinfo)
    TYPE(FILE_INFO), INTENT(inout)  :: fileinfo

    fileinfo%opened           = .FALSE.
    fileinfo%file_id          = 0
    fileinfo%access_mode      = 0
    fileinfo%ncdims(:)        = 0
    fileinfo%creation_program = '?'
    fileinfo%creation_user    = '?'
    fileinfo%creation_date    = '?'
    fileinfo%binary_source    = '?'
    fileinfo%file_type        = '?'
    fileinfo%file_name        = '?'
    fileinfo%title            = '?'

  END SUBROUTINE IO_info_construct
!------------------------------------------------------------------------------
  FUNCTION IO_get_varindx (dimname) RESULT(indx)

    CHARACTER (*), INTENT(in) :: dimname
    INTEGER :: indx

    DO indx = 1, IO_ndim_ids
      IF (IO_dim_ids(indx)%dim_name == dimname) RETURN
    END DO

    indx = -1

!    DO indx = 1, IO_ndim_ids
!      IF (IO_dim_ids(indx)%dim_name(1:1) /= ' ') THEN
!        WRITE(nerr,*) 'IO_get_varindx', indx, IO_dim_ids(indx)%dim_name
!      END IF
!    END DO
!    WRITE (nerr,*) 'element ',dimname,' not available ...'
!    CALL finish('IO_get_varindx','IO_dim_ids error')

  END FUNCTION IO_get_varindx
!------------------------------------------------------------------------------
  SUBROUTINE add_dim (name, len, longname, units, levtyp, single, value, indx)
  !--------------------------------------------------
  ! define a new dimension for Netcdf and GRIB output
  !--------------------------------------------------
  CHARACTER (len=*) ,INTENT(in)            :: name      ! mnemonic
  INTEGER           ,INTENT(in)            :: len       ! size of dimension
  CHARACTER (len=*) ,INTENT(in)  ,OPTIONAL :: longname  ! long name
  CHARACTER (len=*) ,INTENT(in)  ,OPTIONAL :: units     ! units
  INTEGER           ,INTENT(in)  ,OPTIONAL :: levtyp    ! GRIB level type
  LOGICAL           ,INTENT(in)  ,OPTIONAL :: single    ! single layer flag
  REAL          ,INTENT(in)  ,OPTIONAL :: value (:) ! coordinates
  INTEGER           ,INTENT(out) ,OPTIONAL :: indx      ! index in IO_dim_ids 

    INTEGER                :: idx, i
    TYPE (IO_dim) ,POINTER :: dim
    !--------------------
    ! check for zero size
    !--------------------
    if (len == 0) then
      WRITE(nerr,*) 'add_dim',TRIM(name)//' len = 0'
      stop
    endif
    !-----------------------------------------------------------
    ! if dimension is already defined it must have the same size
    !-----------------------------------------------------------
    idx = 0
    do i=1, IO_ndim_ids
      if (IO_dim_ids(i)% dim_name == name) then
        idx = i
        if (IO_dim_ids(i)% dim_len /= len) &
          WRITE(nerr,*) 'add_dim',TRIM(name)//' redifined with different len'
          stop
        exit
      endif
    end do
    !--------------------------------------
    ! increase number of dimensions defined
    !--------------------------------------
    if (idx == 0) then
      IO_ndim_ids = IO_ndim_ids + 1
      IF (IO_ndim_ids > SIZE(IO_dim_ids)) &
        WRITE(nerr,*) 'add_dim','increase size of IO_dim_ids'
        stop
      idx = IO_ndim_ids
      dim => IO_dim_ids(idx)
    endif
    !---------------------------------------------------------
    ! define attributes, chose defaults or optional parameters
    !---------------------------------------------------------
    dim% dim_id   =   0   
    dim% var_id   =   0
    dim% dim_len  =  len
    dim% dim_name =  name
    dim% longname =  ''     ; IF (PRESENT (longname)) dim% longname = longname
    dim% units    =  ''     ; IF (PRESENT (units   )) dim% units    = units
    dim% levtyp   =   0     ; IF (PRESENT (levtyp  )) dim% levtyp   = levtyp
    dim% single   =  .FALSE.; IF (PRESENT (single  )) dim% single   = single
    IF (ASSOCIATED (dim% value)) DEALLOCATE (dim% value)
    IF (PRESENT (value)) THEN
      IF (SIZE (value) /= dim% dim_len) &
        WRITE(nerr,*) 'add_dim',TRIM(name)//' len /= size (value)'
        stop
      ALLOCATE (dim% value (SIZE (value)))
      dim% value = value
    ENDIF
    IF (PRESENT(indx)) indx = idx
  END SUBROUTINE add_dim
!------------------------------------------------------------------------------
  SUBROUTINE add_unknown_dim (len, indx)
  !-----------------------------------------
  ! add a dummy entry in the dimension table
  !-----------------------------------------

  INTEGER   ,INTENT(in)  :: len  ! size of dimension
  INTEGER   ,INTENT(out) :: indx ! index of entry in table

    CHARACTER(len=4) :: nam
    CHARACTER(len=4) :: form
    INTEGER          :: i
    SELECT CASE (len)
    CASE (1:9)
      form ='(i1)'
    CASE (10:99)
      form ='(i2)'
    CASE (100:999)
      form ='(i3)'
    CASE default
      WRITE(nerr,*) 'add_unknown_dim','len < 1 or len > 999'
      stop
  
    END SELECT
    nam = 'n' 
    WRITE (nam(2:4),form) len
    DO i=1, IO_ndim_ids
      IF (IO_dim_ids(i)% dim_name == nam) THEN
        indx = i
        RETURN
      ENDIF
    END DO
    CALL add_dim (nam, len, levtyp=109, indx=indx)    
  END SUBROUTINE add_unknown_dim
!------------------------------------------------------------------------------
  SUBROUTINE IO_DEF_DIM (ncid, name, len, dimid)

    INTEGER :: ncid, len, dimid, status, indx
    CHARACTER*(*) name


    status = NF_DEF_DIM (ncid, name, len, dimid)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_DEF_DIM :', ncid, name, len
      WRITE(nerr,*) 'IO_DEF_DIM', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_DEF_DIM', 'Run terminated.'
      stop  
    END IF

    DO indx = 1, IO_ndim_ids
      IF (name == IO_dim_ids(indx)%dim_name) THEN
        IO_dim_ids(indx)%dim_id = dimid
        RETURN
      END IF
    END DO

    WRITE (nerr,*) 'element ',name,' not available ...'
    WRITE (nerr,*) 'IO_DEF_DIM','IO_dim_ids error'
    stop  


  END SUBROUTINE IO_DEF_DIM
!------------------------------------------------------------------------------
  SUBROUTINE IO_INQ_DIMID (ncid, name, dimid)

    INTEGER :: ncid, dimid, status
    CHARACTER*(*) name


    status = NF_INQ_DIMID (ncid, name, dimid)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_INQ_DIMID :', name
      WRITE(nerr,*) 'IO_INQ_DIMID', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_INQ_DIMID', 'Run terminated.'
      stop  
    END IF

  END SUBROUTINE IO_INQ_DIMID
!------------------------------------------------------------------------------
  SUBROUTINE IO_INQ_DIMLEN (ncid, dimid, len)

    INTEGER :: ncid, dimid, len, status


    status = NF_INQ_DIMLEN (ncid, dimid, len)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_INQ_DIMLEN :', ncid, dimid
      WRITE(nerr,*) 'IO_INQ_DIMLEN', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_INQ_DIMLEN', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_INQ_DIMLEN
!------------------------------------------------------------------------------
  SUBROUTINE IO_INQ_VARID (ncid, name, varid)

    INTEGER :: ncid, varid, status
    CHARACTER*(*) name


    status = NF_INQ_VARID (ncid, name, varid)

    IF (ldebugio) THEN
      WRITE(nerr,*) 'IO_INQ_VARID : ','     Id=',ncid,' varid=',varid,'   ',name
    END IF

    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_INQ_VARID :', ncid, name
      WRITE(nerr,*) 'IO_INQ_VARID', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_INQ_VARID', 'Run terminated.'
      stop
    END IF

  END SUBROUTINE IO_INQ_VARID
!------------------------------------------------------------------------------
  SUBROUTINE IO_DEF_VAR (ncid, name, xtype, nvdims, vdims, varid)

    INTEGER :: ncid, varid, xtype, nvdims
    CHARACTER*(*) name
    INTEGER :: vdims(*)
    INTEGER :: status


    status = NF_DEF_VAR (ncid, name, xtype, nvdims, vdims, varid)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_DEF_VAR :', ncid, name, xtype, nvdims,vdims(1:nvdims)
      WRITE(nerr,*) 'IO_DEF_VAR', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_DEF_VAR', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_DEF_VAR
!------------------------------------------------------------------------------
  SUBROUTINE IO_GET_ATT_TEXT (ncid, varid, name, text)

    INTEGER :: ncid, varid
    CHARACTER*(*) name, text
    INTEGER :: status


    status = NF_GET_ATT_TEXT (ncid, varid, name, text)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_GET_ATT_TEXT :', ncid, varid, name
      WRITE(nerr,*) 'IO_GET_ATT_TEXT', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_GET_ATT_TEXT', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_GET_ATT_TEXT
!------------------------------------------------------------------------------
  SUBROUTINE IO_PUT_ATT_TEXT (ncid, varid, name, text)

    INTEGER :: ncid, varid
    CHARACTER*(*) name, text
    INTEGER :: status
    INTEGER :: lentext


    lentext = len(TRIM(text))

    status = NF_PUT_ATT_TEXT (ncid, varid, name, lentext, text)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_PUT_ATT_TEXT :', ncid, varid, name, lentext, text
      WRITE(nerr,*) 'IO_PUT_ATT_TEXT', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_PUT_ATT_TEXT', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_PUT_ATT_TEXT
!------------------------------------------------------------------------------
  SUBROUTINE IO_GET_ATT_INT (ncid, varid, name, ival)

    INTEGER :: ncid, varid
    CHARACTER*(*) name
    INTEGER :: ival
    INTEGER :: status


    status = NF_GET_ATT_INT (ncid, varid, name, ival)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_GET_ATT_INT :', ncid, varid, name
      WRITE(nerr,*) 'IO_GET_ATT_INT', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_GET_ATT_INT', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_GET_ATT_INT
!------------------------------------------------------------------------------
  SUBROUTINE IO_PUT_ATT_INT (ncid, varid, name, ival)

    INTEGER :: ncid, varid
    CHARACTER*(*) name
    INTEGER :: ival
    INTEGER :: xtype, len
    INTEGER :: status


    len = 1
    xtype = NF_INT
    status = NF_PUT_ATT_INT (ncid, varid, name, xtype, len, ival)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_PUT_ATT_INT :', ncid, varid, name
      WRITE(nerr,*) 'IO_PUT_ATT_INT', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_PUT_ATT_INT', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_PUT_ATT_INT
!------------------------------------------------------------------------------
  SUBROUTINE IO_GET_ATT_DOUBLE (ncid, varid, name, rval)

    INTEGER :: ncid, varid
    CHARACTER*(*) name
    REAL :: rval
    INTEGER :: status


    status = NF_GET_ATT_DOUBLE (ncid, varid, name, rval)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_GET_ATT_DOUBLE :', ncid, varid,  name
      WRITE(nerr,*) 'IO_GET_ATT_DOUBLE', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_GET_ATT_DOUBLE', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_GET_ATT_DOUBLE
!------------------------------------------------------------------------------
  SUBROUTINE IO_PUT_ATT_DOUBLE (ncid, varid, name, rval)

    INTEGER :: ncid, varid
    CHARACTER*(*) name
    REAL :: rval
    INTEGER :: xtype, len
    INTEGER :: status


    len = 1
    xtype = NF_DOUBLE

    status = NF_PUT_ATT_DOUBLE (ncid, varid, name, xtype, len, rval)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_PUT_ATT_DOUBLE :', ncid, varid, name
      WRITE(nerr,*) 'IO_PUT_ATT_DOUBLE', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_PUT_ATT_DOUBLE', 'Run terminated.'
      stop  
    END IF

  END SUBROUTINE IO_PUT_ATT_DOUBLE
!------------------------------------------------------------------------------
  SUBROUTINE IO_GET_VAR_DOUBLE (ncid, varid, dvals)

    INTEGER :: ncid, varid
    REAL :: dvals(*)
    INTEGER :: status


    status = NF_GET_VAR_DOUBLE (ncid, varid, dvals)

    IF (ldebugio) THEN
      WRITE(nerr,*) 'IO_GET_VAR_DOUBLE :',' Id=',ncid,' varid=',varid
    END IF

    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_GET_VAR_DOUBLE :', ncid, varid
      WRITE(nerr,*) 'IO_GET_VAR_DOUBLE', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_GET_VAR_DOUBLE', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_GET_VAR_DOUBLE
!------------------------------------------------------------------------------
  SUBROUTINE IO_GET_VAR_DOUBLE1 (ncid, varid, dval)

    INTEGER :: ncid, varid
    REAL :: dval
    INTEGER :: status


    status = NF_GET_VAR_DOUBLE (ncid, varid, dval)

    IF (ldebugio) THEN
      WRITE(nerr,*) 'IO_GET_VAR_DOUBLE :',' Id=',ncid,' varid=',varid
    END IF

    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_GET_VAR_DOUBLE :', ncid, varid
      WRITE(nerr,*) 'IO_GET_VAR_DOUBLE', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_GET_VAR_DOUBLE', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_GET_VAR_DOUBLE1
!------------------------------------------------------------------------------
  SUBROUTINE IO_GET_VARA_DOUBLE (ncid, varid, start, count, dvals)

    INTEGER :: ncid, varid
    INTEGER :: start(*), COUNT(*)
    REAL :: dvals(*)
    INTEGER :: status

    status = NF_GET_VARA_DOUBLE (ncid, varid, start, count, dvals)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_GET_VARA_DOUBLE :', ncid, varid, start(1:4), COUNT(1:4)
      WRITE(nerr,*) 'IO_GET_VARA_DOUBLE', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_GET_VARA_DOUBLE', 'Run terminated.'
      stop
    END IF

  END SUBROUTINE IO_GET_VARA_DOUBLE
!------------------------------------------------------------------------------
  SUBROUTINE IO_PUT_VARA_DOUBLE (ncid, varid, start, count, dvals)

    INTEGER :: ncid, varid
    INTEGER :: start(*), COUNT(*)
    REAL :: dvals(*)
    INTEGER :: status
    CHARACTER (80) :: varname


    status = NF_PUT_VARA_DOUBLE (ncid, varid, start, count, dvals)
    IF (status /= NF_NOERR) THEN
      status = NF_INQ_VARNAME (ncid, varid, varname)
      WRITE(nerr,*) 'IO_PUT_VARA_DOUBLE :', ncid, varid, varname
      WRITE(nerr,*) 'IO_PUT_VARA_DOUBLE', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_PUT_VARA_DOUBLE', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_PUT_VARA_DOUBLE
!------------------------------------------------------------------------------
  SUBROUTINE IO_PUT_VAR_DOUBLE (ncid, varid, dvals)

    INTEGER :: ncid, varid
    REAL :: dvals(*)
    INTEGER :: status
    CHARACTER (80) :: varname


    IF (ldebugio) THEN
      status = NF_INQ_VARNAME (ncid, varid, varname)
      WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid, varname
      IF (status /= NF_NOERR) THEN
        WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid
        WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE', NF_STRERROR(status)
        WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE', 'Run terminated.'
        stop 
      END IF
    END IF

    status = NF_PUT_VAR_DOUBLE (ncid, varid, dvals)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid
      WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE', 'Run terminated.'
      stop
    END IF

  END SUBROUTINE IO_PUT_VAR_DOUBLE
!------------------------------------------------------------------------------
  SUBROUTINE IO_PUT_VAR_DOUBLE1 (ncid, varid, dval)

    INTEGER :: ncid, varid
    REAL :: dval
    INTEGER :: status
    CHARACTER (80) :: varname


    IF (ldebugio) THEN
      status = NF_INQ_VARNAME (ncid, varid, varname)
      WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid, varname
      IF (status /= NF_NOERR) THEN
        WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid
        WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE', NF_STRERROR(status)
        WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE', 'Run terminated.'
        stop 
      END IF
    END IF

    status = NF_PUT_VAR_DOUBLE (ncid, varid, dval)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid
      WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE', 'Run terminated.'
      stop
    END IF

  END SUBROUTINE IO_PUT_VAR_DOUBLE1
!------------------------------------------------------------------------------
  SUBROUTINE IO_ENDDEF (ncid)

    INTEGER :: ncid, status


    status = NF_ENDDEF (ncid)
    IF (status /= NF_NOERR) THEN
      WRITE(nerr,*) 'IO_ENDDEF', NF_STRERROR(status)
      WRITE(nerr,*) 'IO_ENDDEF', 'Run terminated.'
      stop 
    END IF

  END SUBROUTINE IO_ENDDEF
!------------------------------------------------------------------------------
  SUBROUTINE IO_info_print(fileinfo)
!    USE mo_filename,    ONLY: NETCDF
    TYPE(FILE_INFO) :: fileinfo
    INTEGER :: i

!    CALL message('',' ')
    WRITE(message_text,*) &
         'File name ',TRIM(fileinfo%file_name), &
         ' type ',     TRIM(fileinfo%file_type), &
         ' title ',    TRIM(fileinfo%title)
    WRITE(nerr,*) '',message_text
    SELECT CASE (fileinfo%format)
    CASE (NETCDF)
      WRITE(nerr,*)'',' File type   : NETCDF'
    CASE default
      WRITE(nerr,*)'',' File type   : UNKNOWN'
    END SELECT
    IF (fileinfo%opened) THEN
       WRITE(nerr,*)'',' File open.'
    ELSE
       WRITE(nerr,*)'',' File closed.'
    END IF
    WRITE(message_text,*) &
         'File id :',   fileinfo%file_id,&
         ' access mode :',fileinfo%access_mode
    WRITE(nerr,*)'',message_text
    DO i=1,NF_MAX_VAR_DIMS
       IF (fileinfo%ncdims(i) /= 0) THEN
          WRITE(message_text,*) 'Variable id :',i,' field dimension :',fileinfo%ncdims(i)
          WRITE(nerr,*)'',message_text
       END IF
    END DO
    WRITE(message_text,*) &
         'Created by ',TRIM(fileinfo%creation_user), &
         ' with ',     TRIM(fileinfo%creation_program), &
         ' at ',       TRIM(fileinfo%creation_date), &
         ' source ',   TRIM(fileinfo%binary_source)
    WRITE(nerr,*)'',message_text
    WRITE(nerr,*)'',' '

  END SUBROUTINE IO_info_print


!------------------------------------------------------------------------------
  SUBROUTINE IO_open(filename, fileinfo, mode, iostat)

!    USE mo_filename,    ONLY: rerun_filetype, NETCDF64

    INTEGER,            INTENT(IN)     :: mode
    CHARACTER (*),      INTENT(IN)     :: filename
    TYPE (FILE_INFO), INTENT(INOUT)    :: fileinfo
    INTEGER, OPTIONAL,  INTENT(OUT)    :: iostat

    INTEGER :: status, fill_mode, write_mode
!ps

    IF (PRESENT (iostat)) iostat = 0
    IF (p_pe == p_io) THEN
      IF (fileinfo%opened) THEN
        WRITE(nerr,*) 'IO_open : file ',fileinfo%file_name,' already open'
        WRITE(nerr,*) 'IO_open', 'Run terminated.'
        stop
      END IF

      IF (mode /= IO_READ .AND. mode /= IO_WRITE) THEN
        WRITE (nerr,*) 'IO_open', 'unexpected mode'
        WRITE (nerr,*) 'IO_open', 'Run terminated.'
        stop
      END IF

      ! clean up all elements of netCDF structure

      CALL IO_info_construct(fileinfo)

      fileinfo%file_name = filename

        IF (mode == IO_READ) THEN
          status = NF__OPEN (filename, NF_NOWRITE, chunksize, fileinfo%file_id)
        ELSE
          write_mode = NF_CLOBBER

          IF ( rerun_filetype == NETCDF64 ) THEN
#if defined (HAVE_LIBNETCDF64)
            write_mode = write_mode + NF_64BIT_OFFSET
#else
            WRITE(nerr,*) 'IO_open', 'NETCDF64 not available!'
            stop
#endif /*HAVE_LIBNETCDF64*/
          END IF

          status = NF__CREATE (filename, write_mode, &
                               initialsize, chunksize, fileinfo%file_id)
          status = NF_SET_FILL (fileinfo%file_id, NF_NOFILL, fill_mode)
        END IF

      IF (ldebugio) THEN
        WRITE(nerr,*) 'IO_open  : ', filename, ' mode=', mode, ' Id=', fileinfo%file_id
      END IF

      IF (status == NF_NOERR) THEN
        fileinfo%opened = .TRUE.
      ELSE
        WRITE(nerr,*) 'IO_open failed on', TRIM(filename)
        IF ( fileinfo%format == NETCDF ) THEN
          WRITE(nerr,*) 'IO_open', NF_STRERROR(status)
        ELSE
          WRITE(nerr,*) 'IO_open  : error status = ', status
        END IF
        IF (PRESENT (iostat)) THEN
          iostat = status
        ELSE
          WRITE(nerr,*) 'IO_open', 'Run terminated.'
          stop
        ENDIF
      END IF

    END IF

!    CALL IO_info_print(fileinfo)

  END SUBROUTINE IO_open

  
!------------------------------------------------------------------------------
  SUBROUTINE IO_close(fileinfo)

    TYPE (FILE_INFO), INTENT(INOUT)  :: fileinfo
    INTEGER :: status


    IF (p_pe == p_io) THEN

      IF (ldebugio) THEN
        WRITE(nerr,*) 'IO_close : ',fileinfo%file_name(1:30), ' Id=',fileinfo%file_id
      END IF

      IF (fileinfo%opened) THEN
          status = NF_CLOSE(fileinfo%file_id)

          IF (status == NF_NOERR) THEN
            fileinfo%opened = .FALSE.
          ELSE
            IF ( fileinfo%format == NETCDF ) THEN
              WRITE(nerr,*) 'IO_close', NF_STRERROR(status)
            END IF
              WRITE(nerr,*) 'IO_close', 'Run terminated.'
              stop
          END IF
      ELSE
        WRITE(nerr,*) 'IO_close','file is not open: '//TRIM(fileinfo%file_name)
      ENDIF

      IO_dim_ids(:)%dim_id = -1

    END IF

  END SUBROUTINE IO_close

  
!------------------------------------------------------------------------------
  SUBROUTINE IO_open_unit(unit, fileinfo, mode)

    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN) :: mode
    TYPE (FILE_INFO), INTENT(INOUT)  :: fileinfo
!ps    CHARACTER (13) :: filename
    CHARACTER (100) :: filename
    INTEGER :: lnc


    IF (p_pe == p_io) THEN
      IF (fileinfo%opened) THEN
        WRITE(nerr,*) 'IO_open_unit: unit ',unit,' allready assigned to ', &
             fileinfo%file_name
        WRITE(nerr,*) 'IO_open_unit', 'Run terminated.'
        stop
      END IF

      IF (unit < 10 .OR. unit > 99) THEN
        WRITE(nerr,*) 'IO_open_unit: unit ',unit,' out of range'
        WRITE(nerr,*) 'IO_open_unit', 'Run terminated.'
        stop
      END IF
    END IF

!ps      WRITE (filename,'(A5,I2.2)') 'unit.',unit 
       call chlen (ifilin_nc,80,lnc)
       WRITE (filename, '(A,A1,A5,I2.2)') ifilin_nc(1:lnc),'/','unit.',unit

    CALL IO_open(trim(filename), fileinfo, mode)

  END SUBROUTINE IO_open_unit

  
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_netcdf
    !
    ! deallocate module variables
    !
    INTEGER :: i
    DO i=1, size(IO_dim_ids)
      IF (ASSOCIATED (IO_dim_ids(i)% value)) DEALLOCATE (IO_dim_ids(i)% value)
    END DO
!ps    IF (lsit) THEN
    IF (do_sit) THEN
      IF (ASSOCIATED(diecast_zdepth)) DEALLOCATE (diecast_zdepth)
      IF (ASSOCIATED(diecast_fluxdepth)) DEALLOCATE (diecast_fluxdepth)
      IF (ASSOCIATED(sit_zdepth)) DEALLOCATE (sit_zdepth)
      IF (ASSOCIATED(sit_fluxdepth)) DEALLOCATE (sit_fluxdepth)
      IF (ASSOCIATED(sit_gribzdepth)) DEALLOCATE (sit_gribzdepth)
      IF (ASSOCIATED(sit_gribfluxdepth)) DEALLOCATE (sit_gribfluxdepth)
      IF (ASSOCIATED(ocn_z)) DEALLOCATE (ocn_z)
    ENDIF

  END SUBROUTINE cleanup_netcdf

  
!----------------------------------------------------------------------------------------------------
      SUBROUTINE set_ocndepth()        
      IMPLICIT NONE
      INTEGER :: ocn_k0                    ! vertical dimension parameter (number of layer interfaces, equal
                                         ! the number of layers or pressure levels plus 1, does not include
                                         ! ghost zone)  
      real:: tmp,dzz,b,c,d,aa,zz
      INTEGER:: jk,k
      IF (zocn_option.EQ.1) THEN
     ! RESET ocn_k1
        ocn_k1=50
        WRITE (nerr,*) 'zocn_option=',zocn_option
        WRITE (nerr,*) 'RESET ocn_k1=',ocn_k1
      ELSEIF (zocn_option.EQ.2) THEN
     ! RESET ocn_k1
       ocn_k1=40
       WRITE (nerr,*) 'zocn_option=',zocn_option
       WRITE (nerr,*) 'RESET ocn_k1=',ocn_k1
      ENDIF
      ocn_k0=ocn_k1+1   ! = 41  number of ocean ocean levels plus 1
      ALLOCATE (ocn_z(ocn_k0+ocn_k1))
      ALLOCATE (diecast_zdepth(1:ocn_k1))
      ALLOCATE (diecast_fluxdepth(1:ocn_k1))
      IF (zocn_option.EQ.1) THEN
      ! fine in thermocline, 10 m apart within 100-225 m, but coarse (~200 m) in deep-water formation depths (>500 m).
        diecast_zdepth=(/                                                                              &
          1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,    &
         12.,  15.,  25.,  35.,  45.,  55.,  65.,  75.,  85.,  95.,    &
        105., 115., 125., 135., 145., 155., 165., 175., 185., 195.,    &
        205., 215., 225., 238., 262., 303., 366., 459., 584., 747.,    &
        949.,1193.,1479.,1807.,2174.,2579.,3016.,3483.,3972.,4478./)
        diecast_fluxdepth(1)=(0.+diecast_zdepth(1))/2.
        DO jk = 2, ocn_k1
          diecast_fluxdepth(jk)=(diecast_zdepth(jk-1)+diecast_zdepth(jk))/2.
        ENDDO
        DO jk = 1, ocn_k1
          ocn_z(2*jk-1)=diecast_fluxdepth(jk)    
          ocn_z(2*jk)=diecast_zdepth(jk)
        ENDDO
        ocn_z(ocn_k0+ocn_k1)=ocn_tlz
       ELSEIF (zocn_option.EQ.2) THEN
       ! fine in thermocline, 10 m apart within 100-225 m, but coarse (~200 m) in deep-water formation depths (>500 m).
         diecast_zdepth=(/                                                                              &
          5.,  15.,  25.,  35.,  45.,  55.,  65.,  75.,  85.,  95.,    &
        105., 115., 125., 135., 145., 155., 165., 175., 185., 195.,    &
        205., 215., 225., 238., 262., 303., 366., 459., 584., 747.,    &
        949.,1193.,1479.,1807.,2174.,2579.,3016.,3483.,3972.,4478./)
        diecast_fluxdepth(1)=(0.+diecast_zdepth(1))/2.
        DO jk = 2, ocn_k1
          diecast_fluxdepth(jk)=(diecast_zdepth(jk-1)+diecast_zdepth(jk))/2.
        ENDDO
        DO jk = 1, ocn_k1
          ocn_z(2*jk-1)=diecast_fluxdepth(jk)
          ocn_z(2*jk)=diecast_zdepth(jk)
        ENDDO
        ocn_z(ocn_k0+ocn_k1)=ocn_tlz
      ELSE 
        CALL DIECAST_ZCORD(ocn_k0,ocn_k1,ocn_z)
        DO jk = 1, ocn_k1
          diecast_zdepth(jk)=ocn_z(2*jk)
          diecast_fluxdepth(jk)=ocn_z(2*jk-1)
       ENDDO
      ENDIF
      IF (lssst) THEN
  !   nfnlvl=12     ! number of fine layers in sit for the uppermost ocean layer
        nfnlvl=MIN(INT(diecast_zdepth(2)-1.),10)+2     ! number of fine layers in sit for the uppermost ocean layer
                                                        ! 0,0.0005,1,2,...,MIN(INT(diecast_zdepth(2)),10.)
      ELSE
        nfnlvl=2      ! number of fine layers in sit for the uppermost ocean layer (0m, 0.0005m)
      ENDIF
      lkvl=ocn_k1+nfnlvl-3    ! number of vertical levels of sit ocean model, =49, lssst=T; =39 , lssst=F 
#ifdef DEBUG
      if(myrank.eq.0) print*,"ps: lkvl=",lkvl,",ocn_k1=",ocn_k1,",nfnlvl=",nfnlvl
#endif
      ALLOCATE (sit_zdepth(0:lkvl+1))
      ALLOCATE (sit_fluxdepth(0:lkvl+1))
      ALLOCATE (sit_gribzdepth(1:lkvl+2))
      ALLOCATE (sit_gribfluxdepth(1:lkvl+2))
  ! 
  ! for capture the diurnal cycle of SST, we need the following z-cord for upper ocean
  !
      IF (lssst) THEN
        sit_zdepth(0:1)=(/0.,0.0005/)         !=zdepth(0)*10=-0 m in dm
       DO k=1,nfnlvl-2
        !sit_zdepth(k+1)=REAL(k,dp)
        sit_zdepth(k+1)=REAL(k)
       ENDDO
      ELSE
       sit_zdepth(0:nfnlvl-1)=(/0.,diecast_zdepth(1)/)
      ENDIF
  !
  ! for combining with 3-D ocean model (DIECAST), we need the following z-cord for ocean
  !
    sit_zdepth(nfnlvl:lkvl+1)=diecast_zdepth(2:ocn_k1)                 ! convert unit from m to m
                                                                       ! for the index, it can be seen that lkvl+1-nfnlvl=ocn_k1-2
                                                                       ! i.e., lkvl=ocn_k1+nfnlvl-3
      IF (.TRUE.) THEN
       sit_fluxdepth(0)=sit_zdepth(0)
       DO jk = 1,nfnlvl !i.e., jk = 1,12
        sit_fluxdepth(jk)=(sit_zdepth(jk-1)+sit_zdepth(jk))/2.
       ENDDO
       sit_fluxdepth(nfnlvl+1:lkvl+1)=diecast_fluxdepth(3:ocn_k1)           ! convert unit from m to m
      ! modify diecast_fluxdepth(2) and ocn_z(3) based on lssst option 
       diecast_fluxdepth(2)=sit_fluxdepth(nfnlvl)
       ocn_z(3)=diecast_fluxdepth(2)
      ELSE
       sit_fluxdepth(0)=sit_zdepth(0)
       DO jk = 1,nfnlvl-1 !i.e., jk = 1,11
        sit_fluxdepth(jk)=(sit_zdepth(jk-1)+sit_zdepth(jk))/2.
       ENDDO
       sit_fluxdepth(nfnlvl:lkvl+1)=diecast_fluxdepth(2:ocn_k1)           ! convert unit from m to m  
      ENDIF
      IF (p_parallel_ocean) THEN
       OPEN (UNIT=nocn_sv,FILE="diecast_zdepth.txt",ACTION="WRITE",STATUS="REPLACE")
       WRITE (nocn_sv,*) 'zaxistype = depth_below_sea'
       WRITE (nocn_sv,*) 'size =',  size(diecast_zdepth)
       WRITE (nocn_sv,*) 'levels =',diecast_zdepth
       CLOSE (nocn_sv)
       OPEN (UNIT=nocn_sv,FILE="diecast_fluxdepth.txt",ACTION="WRITE",STATUS="REPLACE")
       WRITE (nocn_sv,*) 'zaxistype = depth_below_sea'
       WRITE (nocn_sv,*) 'size =',  size(diecast_fluxdepth)
       WRITE (nocn_sv,*) 'levels =',diecast_fluxdepth
       CLOSE (nocn_sv)
       OPEN (UNIT=nocn_sv,FILE="zdepth.txt",ACTION="WRITE",STATUS="REPLACE")
       WRITE (nocn_sv,*) 'zaxistype = depth_below_sea'
       WRITE (nocn_sv,*) 'size =',  size(sit_zdepth)
       WRITE (nocn_sv,*) 'levels =',sit_zdepth
       CLOSE (nocn_sv)
       OPEN (UNIT=nocn_sv,FILE="fluxdepth.txt",ACTION="WRITE",STATUS="REPLACE")
       WRITE (nocn_sv,*) 'zaxistype = depth_below_sea'
       WRITE (nocn_sv,*) 'size =',  size(sit_fluxdepth)
       WRITE (nocn_sv,*) 'levels =',sit_fluxdepth
       CLOSE (nocn_sv)
       OPEN (UNIT=nocn_sv,FILE="ocnz_depth.txt",ACTION="WRITE",STATUS="REPLACE")
       WRITE (nocn_sv,*) 'zaxistype = depth_below_sea'
       WRITE (nocn_sv,*) 'size =',  size(ocn_z)
       WRITE (nocn_sv,*) 'levels =',ocn_z
       CLOSE (nocn_sv)
      ENDIF
      IF (lssst) THEN    
       sit_gribzdepth(1:2)=(/0.,                                   &    !=zdepth(0)*10=-0 m in dm
        1./)                                                          !-0.0005 msit_zdepth(0)
       sit_gribzdepth(3:lkvl+2)=sit_zdepth(2:lkvl+1)*10.                ! convert unit from m to dm (i.e., 0.1 m)
!     
       sit_gribfluxdepth(1:2)=(/0.,                               &    !=-0 m                               
        1./)                                                          !=(zdepth(0)+zdepth(1))/2 =-0.00025 m
       sit_gribfluxdepth(3:lkvl+2)=sit_fluxdepth(2:lkvl+1)*10.         !=(zdepth(1)+zdepth(2))/2 =-0.50025 m to dm (i.e., 0.1 m)
      ELSE
       sit_gribzdepth(1:lkvl+2)=sit_zdepth(0:lkvl+1)*10.                ! convert unit from m to dm (i.e., 0.1 m)
       sit_gribfluxdepth(1:lkvl+2)=sit_fluxdepth(0:lkvl+1)*10.         !=(zdepth(1)+zdepth(2))/2 =-0.50025 m to dm (i.e., 0.1 m)
      ENDIF
      IF ( (lwarning_msg.GE.2).AND. p_parallel_io ) THEN
       WRITE(nerr,*) "pe=",mpp_pe(),"I am in mo_netcdf 0.0: IO_init_dims"    
       WRITE (nerr,*) 'diecast_zdepth(m)= ',diecast_zdepth
       WRITE (nerr,*) 'sit_zdepth(m)= ',sit_zdepth
!!!      WRITE (nerr,*) 'sit_gribzdepth(dm)= ',sit_gribzdepth(lkvl-1:lkvl+2)
       WRITE (nerr,*) 'sit_gribzdepth(dm)= ',sit_gribzdepth
       WRITE (nerr,*) 'sit_fluxdepth(m)= ',sit_fluxdepth
       WRITE (nerr,*) 'sit_gribfluxdepth(dm)= ',sit_gribfluxdepth
      ENDIF   

      CONTAINS	
      SUBROUTINE DIECAST_ZCORD(k0,k1,ocn_z)
!   generate combined linear-exponential expanding coordinate Z(zz)
!   Z=aa+b*exp(c*zz)+d*zz, 0 < zz < tlz
!   where zz is the linear (unstretched) space coordinate
!   NOTE: for expanding coordinate, 0 < d < 1
!   at zz=0, Z=aa+b*exp(0)=0, so aa+b=0.
!   at zz=tlz: Z=aa+b*exp(c*tlz)+d*tlz=tlz
!   so aa*(1-exp(c*tlz))=tlz-d*tlz, aa=tlz*(1-d)/(1-exp(c*tlz))
!  convert tlz to meters
      INTEGER, INTENT(IN):: k0,k1
      REAL, INTENT(out):: ocn_z(k0+k1)
      INTEGER:: k
      REAL:: tmp,dzz,b,c,d,aa,zz
      tmp=ocn_tlz
      dzz=tmp/(2.*k1)
!   exponential part
      c=.001
!   linear part
      d=0.03
      aa=tmp*(1.-d)/(1.-exp(c*tmp))
      b=-aa
      zz=0.
      DO k=1,k0+k1
!   z(k) in m
       ocn_z(k)=(aa+b*exp(c*zz)+d*zz)
       zz=zz+dzz
      ENDDO
      END SUBROUTINE DIECAST_ZCORD

      END SUBROUTINE set_ocndepth
 
	  
      END MODULE mo_netcdf
