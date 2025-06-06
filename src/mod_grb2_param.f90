module mod_grb2_param
!************************************************************************
!* ABSTRACT
!*
!* This module is for CWB/GFS T511 input/output data condition.
!*
!* All the parameter and the corresponding subroutine used to create
!*
!* GRIB2 data are listing bellow.
!*
!* The mark '**' before the variables means the variables have to
!setted*
!* in the main program. The mark '##' before the variables means the
!*
!* output variabels.
!*
!*
!* EDITED BY
!*
!* Jia-ying Wu 2017,09
!*
!************************************************************************
implicit none
public
integer*4 :: ierr
!GRIB1 SECTION 0 & 1
integer*4 :: listsec0(2),listsec1(13)
!GRIB2 SECTION 2 SKIP
!GRIB2 SECTION 3
integer*4 :: igds(5),idefnum,ideflist
integer*4,parameter :: igdstmplen=19 !Max dimension of igdstmpl()
integer*4 :: igdstmpl(igdstmplen)
!GRIB2 SECTION 4
integer*4 :: ipdsnum ,ipdsnum_acc
integer*4,parameter :: ipdstmplen=18 !Max dimension of ipdstmpl()
integer*4 :: ipdstmpl(ipdstmplen)
integer*4,parameter :: ipdstmplen8=32 !Max dimension of ipdstmpl()
integer*4 :: ipdstmpl8(ipdstmplen8)
integer*4,parameter :: numcoord=1 !number of values in array
                                      !coordlist.
      real :: pcoord
      real*4 :: coordlist(numcoord)    !Array containg floating point
                                     !values
                                     ! intended to document
      !integer,parameter :: numcoord2=60, plst26=26, plst16=16, plst31=31
      !real :: coordlist2(numcoord2),                               &
      !        pcoord26(plst26),pcoord16(plst16),pcoord31(plst31)
      !integer,parameter :: llst5=5
      !real :: lcoord5(llst5)
      !GRIB2 SECTION 5
      integer*4 :: idrsnum40,idrsnum0 ,idrsnum3
                 !Max dimension of idrstmpl()
      integer*4,parameter :: idrstmplen40=7, idrstmplen0=5 ,idrstmplen3=18
      integer*4 :: idrstmpl40(idrstmplen40),idrstmpl0(idrstmplen0),idrstmpl3(idrstmplen3)
      integer*4 :: ngrdpts
      !GRIB2 SECTION 6
      integer*4 :: ibmap

      character,allocatable,save :: cgrib(:)*1
      integer*4:: lcgrib=2*1e8,lengrib
      logical*1,allocatable,save :: bmap(:)
      !add default 
      integer*4  :: grbnxmy
      !grib2 file name
      character::grbfile*255 ,ofdir*255
      integer*4::grbid=23,grbparid=24
      integer*8::grb_idtg
      !for ensemble
      integer::grbmem=-1 ,grbnumm
      !=======================================================================
      !  call baopenw(g2num,g2name,ierr)                                     !
      !=======================================================================
      !**   g2num=      !Integer num. of created Grib2 file
      !**   g2name=     !Charater of created Grib2 file name
      !=======================================================================
      !  call gribcreate(cgrib,lcgrib,listsec0,listsec1,ierr)                !
      !=======================================================================
      !##     cgrib     Character array to contain the GRIB2 message
      !**     lcgrib    Maximum length (bytes) of array cgrib.

      ! Grib2 section 0 (Must be dimensioned >= 2)
       data listsec0/0,2/
      ! listsec0(1)=0  !Discipline-GRIB Master Table Number (Code Table 0.0)
      ! listsec0(2)=2  !GRIB Edition Number (currently 2)

      ! Grib2 section 1 (Must be dimensioned >= 13)
       data listsec1/138, 10,  2,  1,  1, -1, -1, -1, -1,  0,  0,  0, -1/
      ! listsec1(1)=137    !Id of orginating centre (Code Table C-1)
      ! listsec1(2)=0      !Id of orginating sub-centre 
      ! listsec1(3)=9      !GRIB Master Tables Version Number (Code Table 1.0)
      ! listsec1(4)=0      !GRIB Local Tables Version Number  (Code Table 1.1)
      ! listsec1(5)=1      !Significance of Reference Time    (Code Table 1.2)
      ! listsec1(6)=yy     !Reference Time - Year (4 digits)
      ! listsec1(7)=mon    !Reference Time - Month
      ! listsec1(8)=dd     !Reference Time - Day
      ! listsec1(9)=hh     !Reference Time - Hour
      ! listsec1(10)=0     !Reference Time - Minute
      ! listsec1(11)=0     !Reference Time - Second
      ! listsec1(12)=0     !Production status of data        (Code Table 1.3)
      ! listsec1(13)=tproc !Type of processed data           (Code Table 1.4)

      !=======================================================================
      !  call addlocal(cgrib,lcgrib,csec2,lcsec2,ierr)                       !
      !=======================================================================
      ! Grib2 section 2 skip!!!

      !=======================================================================
      ! call addgrid(cgrib,lcgrib,igds,igdstmpl,igdstmplen,ideflist,idefnum,ierr)
      !=======================================================================
      !*#     cgrib
      !**     lcgrib

      ! Grib2 section 3 (Must be dimensioned >= 5)
       data igds/0,-1,0,0,-1/
      !   igds(1)=0     !Source of grid definition (see Code Table 3.0)
      !** igds(2)=nx*my !Number of grid points in the defined grid.
      !   igds(3)=0     !Number of octets needed for each additional grid points definition.
      !                 !Used to define number of points in each row ( or column ) 
      !                 !for non-regular grids.  = 0, if using regular grid.
      !   igds(4)=0     !Interpretation of list for optional points definition. (Code Table 3.11)
      !   igds(5)=40    !Grid Definition Template Number (Code Table 3.1)
                        !igds(5)=40 :Gaussian latitude/longitude
                        !igds(5)=0  :Equal latitude/longitude

      ! Contains the data values for the specified Grid Definition Template
      ! ( NN=igds(5) ).  Each element of this integer array contains an entry
      ! (in the order specified) of Grid Defintion Template 3.NN
      !      data igdstmpl//
       data igdstmpl/ 6,        0,  0 , 0,  0,  &
                      0,        0, -1, -1,  0,  &
                      0,-90000000,  0, 48, -1,  &
                     -1,       -1, -1, 64/
      !                     0,-90000000,  0,00000000, -1,  &
      !                    -1,       -1, -1,01000000/
      ! igdstmpl(1)=6      !Shpape of the Earth( See Code talbe 3.2)
      ! igdstmpl(2)=0      !Scale factor of radius of spherical Earth
      ! igdstmpl(3)=0      !Scale value of radius of spherical Earth
      ! igdstmpl(4)=0      !Scale factor of major axis of oblate spheroid Earth
      ! igdstmpl(5)=0      !Scale value of major axis of oblate spheroid Earth
      ! igdstmpl(6)=0      !Scale factor of minor axis of oblate spheroid Earth
      ! igdstmpl(7)=0      !Scale value of minor axis of oblate spheroid Earth
      ! igdstmpl(8)=1536   !Ni - number of points of a parallel
      ! igdstmpl(9)=768    !Nj - number of points of a meridian
      ! igdstmpl(10)=0     !Basic angle of the inital production domain
      ! igdstmpl(11)=0     !Subdivisions of basic angle used to define extreme lon. and lat., and direction increments
      ! igdstmpl(12)=-90000000  !La1 - latitude of first grid point
      ! igdstmpl(13)=0          !Lo1 - longitude of first grid point
      ! igdstmpl(14)=00000000   !Resolution and component flags ( See Flage Table 3.3 )
      ! igdstmpl(15)=89765625   !La2 - latitude of last grid point
      ! igdstmpl(16)=359765625  !Lo2 - longitude of last grid point
      ! igdstmpl(17)=234375     !Di - i direction increment
      ! igdstmpl(18)=768        !N - number of parallels between a pole and the equator
      ! igdstmpl(19)=01000000  &!Scanning mode( See Flage Table 3.4 )
        data ideflist/0/
      ! ideflist=0     !(Used if igds(3) .ne. 0)  This array contains the
      !                     number of grid points contained in each row ( or column )
        data idefnum/1/
      ! idefnum=1      !(Used if igds(3) .ne. 0)  The number of entries in array
      !                 ideflist.  i.e. number of rows ( or columns ) for which
      !                 optional grid points are defined.

      !=======================================================================
      !  call addfield(cgrib,lcgrib,ipdsnum,ipdstmpl,ipdstmplen,coordlist,   !
      !                numcoord,idrsnum,idrstmpl,idrstmplen,fld,ngrdpts,     !
      !                ibmap,bmap,ierr)                                      !
      !=======================================================================
      !*#     cgrib
      !**     lcgrib

      ! Grib2 section 4
       data ipdsnum/0/
       data ipdsnum_acc/8/
      !  ipdsnum=0    !Product Definition Template Number ( see Code Table 4.0)
      !  ipdsnum=0  :Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
      ! Contains the data values for the specified Product Definition Template
      ! ( N=ipdsnum ).  Each element of this integer array  contains an entry
      ! (in the order specified) of Product Defintion Template 4.N
       data ipdstmpl/ -1, -1, -1,  0, 53,  0,  0,  1, -1, -1, -1, -1,255,  0,  0, 0, 0, 0/
      ! ipdstmpl(1)=    !Parameter category ( See Code Table 4.1 )
      ! ipdstmpl(2)=    !Parameter number ( See Code Table 4.2 )
      ! ipdstmpl(3)=    !Type of generating process ( See Code Table 4.3 )
      ! ipdstmpl(4)=0   !Background generating process identifier
      ! ipdstmpl(5)=    !Analysis or forecast generating process identified ( See Code ON388 Table A )
      ! ipdstmpl(5)=81 : Analysis from GFS(Global Forecast System)
      ! ipdstmpl(6)=0   !Hours of observational data cutoff after reference time
      ! ipdstmpl(7)=0   !Minutes of observational data cutoff after reference time
      ! ipdstmpl(8)=1   !Indicator of unit of time range ( See Code Table 4.4 )
      ! ipdstmpl(9)=6   !Forecast time in units defined by ipdstmpl(8)
      ! ipdstmpl(10)=   !Type of first fixed surface( See Code Table 4.5 )
      ! ipdstmpl(11)=0  !Scale factor of first fixed surface
      ! ipdstmpl(12)=0  !Scaled value of first fixed surface
      ! ipdstmpl(13)=255 !Type of second fixed surface( See Code Table 4.5 )
      ! ipdstmpl(14)=0  !Scale factor of second fixed surface
      ! ipdstmpl(15)=0  !Scaled value of second fixed surface

       data ipdstmpl8/ -1, -1, -1,  0, 53,  0,  0,  1, -1, -1,  &
                       -1, -1,255,  0,  0, -1, -1, -1,          &
                       -1, -1, -1, -1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
      !  ipdstmpl(16)= 0 ! Year   | Time of end of overall time interval
      !  ipdstmpl(17)= 0 ! Month  | Time of end of overall time interval
      !  ipdstmpl(18)= 0 ! Day    | Time of end of overall time interval
      !  ipdstmpl(19)= 0 ! Hour   | Time of end of overall time interval
      !  ipdstmpl(20)= 0 ! Minute | Time of end of overall time interval
      !  ipdstmpl(21)= 0 ! Second | Time of end of overall time interval
      !  ipdstmpl(22)= 1 ! n - Number of time range specifications describing
                         !     the time intervals used to calculate the
                         !     statistically processed field
      !  ipdstmpl(23)= 0 ! Total number of data values missing in statistical process.
      !  ipdstmpl(24)= 0 ! Statistical process used to calculate the processed
                         ! field from the field at each time increment during
                         ! the time range (Code Table 4.10)= 0
                         ! 0 --> Average  1 --> Accum
      !  ipdstmpl(25)= 2 ! Type of time increment between successive fields
                         ! used in the statistical processing (Code Table 4.11)
                         ! 2 --> Sucessive times processed have same start
                         ! time of forecast, forecast time is incremented
      !  ipdstmpl(26)= 1 ! Indicator of unit of time for time range over which
                         ! statistical processing is done (Code Table 4.4)
                         ! 1 --> hour
      !  ipdstmpl(27)= 1 ! Length of the time range over which statistical
                         ! processing is done, in units defined by the
                         ! previous octet
      !  ipdstmpl(28)= 1 ! Indicator of unit of time for the increment between
                         ! the successive fields used (Code Table 4.4)
                         ! 1 -- Hour
      !  ipdstmpl(29)= 6 ! Time increment between successive fields, in units
                         ! defined by the previous octet (see Note 3 and 4)

       data coordlist/1/
      !coordlist(1)=0 !Array containg floating point values intended to
                      !document the vertical discretisation associated to
                      !model
                      !data on hybrid coordinate vertical levels.
      ! data coordlist2/          0.00085400,0.00182570,0.00293106, &
      !     0.00418811,0.00561725,0.00724145,0.00908662,0.01118186, &
      !     0.01355986,0.01625721,0.01931478,0.02277809,0.02669766, &
      !     0.03112936,0.03613468,0.04178099,0.04814158,0.05529569, &
      !     0.06332831,0.07232963,0.08239436,0.09362052,0.10610791, &
      !     0.11995598,0.13526124,0.15211401,0.17059470,0.19076949, &
      !     0.21268569,0.23636683,0.26180792,0.28897095,0.31778133, &
      !     0.34812545,0.37984984,0.41276224,0.44663487,0.48120970, &
      !     0.51620575,0.55132793,0.58627691,0.62075919,0.65449686, &
      !     0.68723622,0.71875478,0.74886625,0.77742342,0.80431894, &
      !     0.82948416,0.85288658,0.87452598,0.89442998,0.91264919, &
      !     0.92925232,0.94432163,0.95794864,0.97023054,0.98126698, &
      !     0.99058760,1.00000000/
      !
      ! data pcoord16/ 10.0, 20.0, 30.0, 50.0, 70.0,    &
      !               100.0,150.0,200.0,250.0,300.0,    &
      !               400.0,500.0,700.0,850.0,925.0,1000.0/
      !
      ! data pcoord26/ 10.0, 20.0, 30.0, 50.0, 70.0,100.0,150.0,200.0, &
      !               250.0,300.0,350.0,400.0,450.0,500.0,550.0,600.0, &
      !               650.0,700.0,750.0,800.0,850.0,900.0,925.0,950.0, &
      !               975.0,1000.0/
      !
      ! data pcoord31/  1.0,  2.0,  3.0,  5.0,  7.0, 10.0, 20.0, 30.0, &
      !                      50.0, 70.0,100.0,150.0,200.0,250.0,300.0,350.0, &
      !                     400.0,450.0,500.0,550.0,600.0,650.0,700.0,750.0, &
      !                     800.0,850.0,900.0,925.0,950.0,975.0,1000.0/
      !
      ! data lcoord5/0.0,10.0,40.0,100.0,200.0/
      !=======================================================================
      ! Grib2 section 5
       data idrsnum40/40/  !jpeg compress
       data idrsnum0  /0/  !normal
       data idrsnum3  /3/  !complex 3 compress 
      !idrsnum=0    !Data Representation Template Number ( see Code Table 5.0 )
      !idrsnum=40   !JPEG 2000 Code Stream Format

      !Contains the data values for the specified Data Representation Template
      ! ( N=idrsnum ). Each element of this integer array contains an entry
      ! (in the order specified) of Data Representation Template 5.N.
      ! Note that some values in this template (eg. reference values, number
      ! of bits,
      ! etc...) may be changed by the data packing algorithms. Use this to
      ! specify
      ! scaling factors and order of spatial differencing, if desired.
        data idrstmpl40/0,0,2,0,0,0,255/
      ! idrstmpl(1)=0  !Reference Value (R)(IEEE 32-bit folating-point value)
      ! idrstmpl(2)=0  !Binary scale factor (E)
      ! idrstmpl(3)=0  !Decimal scale factor (D)
      ! idrstmpl(4)=8  !Number of bits required to hold the resulting scaled and
      !          referenced data values.(i.e. The depth of the grayscale image.)
      ! idrstmpl(5)=    !Type of original field values( see Code Table 5.1 )
                        ! 0:folating points, 1:integer
      ! idrstmpl(6)=0    !Type of Compression used.( see Code Table 5.40 )
      ! idrstmpl(7)=255  !Target compression ratio, M:1.
                     !with respect to the bit-depth specified in
                     !idrstmpl(4), when idrstmpl(6) indicates Lossy
                     !Compression. Otherwise, set to missing.
        data idrstmpl0/0,0,2,0,0/
      ! idrstmpl(1)=0  !Reference Value (R)(IEEE 32-bit folating-point value)
      ! idrstmpl(2)=0  !Binary scale factor (E)
      ! idrstmpl(3)=0  !Decimal scale factor (D)
      ! idrstmpl(4)=8  !Number of bits required to hold the resulting scaled and
      !            referenced data values.(i.e. The depth of the grayscale image.)
      ! idrstmpl(5)=    !Type of original field values( see Code Table 5.1 )
                        ! 0 -> folating points, 1 -> integer
      !  ngrdpts=nx*my  !Number of data points in grid. i.e. size of fld and bmap.
       data idrstmpl3/0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0/
      !=======================================================================
      ! grib2 section 6
        data ibmap/255/
      !      ibmap=255  !Bitmap indicator ( see Code Table 6.0 )
      !                 0 = bitmap applies and is included in Section 6.
      !                 1-253 = Predefined bitmap applies
      !                 254 = Previously defined bitmap applies to this field
      !                 255 = Bit map does not apply to this product.
      ! bmap(nx*my)=    Logical*1 array containing bitmap to be added. ( if ibmap=0 or ibmap=254)
      !=======================================================================
      ! grib2 section 7
      ! fld(nx*my)          !Array of data points to pack.
      !=======================================================================
      ! grib2 section 8
      !
      ! call gribend(cgrib,lcgrib,lengrib,ierr)                             
      !*#     cgrib
      !**     lcgrib
      !##     lengrib   !Length of the final GRIB2 message in octets (bytes)
      !=======================================================================
      !end module cwbgfs_param
      !=======================================================================
      contains
      !=======================================================================
            subroutine seclist01(idtg,itau)
            !use cwbgfs_param

            integer*4 yy,mm,dd,hh
            integer*8::idtg
            integer   itau
            character::cdtg*12
            grb_idtg=idtg
            write(cdtg,'(I12.12)')idtg
            read(cdtg(1:4),'(I4)')yy
            read(cdtg(5:6),'(I2)')mm
            read(cdtg(7:8),'(I2)')dd
            read(cdtg(9:10),'(I2)')hh
      ! Create GRIB2 title (section 0&1)
             listsec0(1)=0 ! sec00

             listsec1(6)=yy     !Reference Time - Year (4 digits)
             listsec1(7)=mm     !Reference Time - Month
             listsec1(8)=dd     !Reference Time - Day
             listsec1(9)=hh     !Reference Time - Hour
            if( itau .eq. 0 )then
             listsec1(13)=0 !Analysis data (Code Table 1.4)
            else
             listsec1(13)=1 !Forecast data
            endif

            return
            end subroutine
      !=======================================================================
            subroutine wrt_grb2(itau,t0,t1,t2,p3,t10,t11,t12,fld)
            use param, only : io_quilting
            use rank,only:ntag
            use const, only: RTYPE
            integer::  t0,t1,t2,t10,t11,p3
            real::  t12
            integer::ptp0(9)
            integer::itau,ist
            real(kind=RTYPE)::fld(grbnxmy)
            real*4::r4out(grbnxmy)
            character:: keydoit*34
            data keydoit/"DOIT..........................DOIT"/
            if(io_quilting)then
              ptp0=(/t0,t1,t2,p3,t10,t11,nint(t12),-999,-999/)
              ntag=ntag+1
              call mpe_send_key( keydoit  ,ntag,ist)
              ntag=ntag+1
              call mpe_send_int( ptp0 , 9 ,ntag,ist)
              ntag=ntag+1
              call mpe_send_data(fld ,grbnxmy ,ntag,ist)
            else
              r4out(:)=fld(:)
              call  wrt_grb2_io(itau,t0,t1,t2,p3,t10,t11,nint(t12), r4out )
            endif
            end
      !=======================================================================
            subroutine wrt_grb2_v2(itau,t0,t1,t2,p3,t10,t11,t12,fld)
            use param, only : nx,my
            use const, only: idtg,RTYPE,ihdgo2
            integer::  t0,t1,t2,t10,t11,p3
            integer::  t12
            integer::ptp0(9)
            integer::itau,ist,istat,i,ia
            real(kind=RTYPE)::fld(nx*my)
            real*4::r4out(nx*my)
!            character:: ihdg2*26,clen*7

             do i = 1, 26
              ia=ichar(ihdgo2(i:i))
              if((ia.ge.97).and.(ia.le.122))then
                ia=ia-32
                ihdgo2(i:i)=char(ia)
              endif
             enddo

 133                      format( A  ,A1 ,A6,A4 )
#ifdef O38K
              write(grbfile,133 )trim(ofdir),'/',ihdgo2(1:6),ihdgo2(13:16)
#else
              write(grbfile,133 )trim(ofdir),'/',ihdgo2(1:6),ihdgo2(11:14)
#endif
              r4out(:)=fld(:)
              call opn_grb2(grbparid,nx,my,idtg,itau,istat)              
              call wrt_grb2_io(itau,t0,t1,t2,p3,t10,t11,t12, r4out )
              call cls_grb2(grbparid,istat)
            end
      !=======================================================================
            subroutine wrt_grb2_accu(itau,t0,t1,t2,p3,t10,t11,t12,t24,t27,fld)
            use param, only : io_quilting
            use rank,only:ntag
            use const, only: RTYPE
            integer::  t0,t1,t2,t10,t11,p3,t24,t27
            real::  t12
            integer::ptp0(9)
            integer::itau,ist
            real(kind=RTYPE)::fld(grbnxmy)
            real*4::r4out(grbnxmy)
            character:: keydoit*34
            data keydoit/"DOIT..........................DOIT"/
            if(io_quilting)then
              ptp0=(/t0,t1,t2,p3,t10,t11,nint(t12), t24 ,t27/)
              ntag=ntag+1
              call mpe_send_key( keydoit  ,ntag,ist)
              ntag=ntag+1
              call mpe_send_int( ptp0 , 9 ,ntag,ist)
              ntag=ntag+1
              call mpe_send_data(fld ,grbnxmy ,ntag,ist)
            else
              r4out(:)=fld(:)
              call wrt_grb2_accu_io(itau,t0,t1,t2,p3,t10,t11,nint(t12),t24,t27,r4out)
            endif
            end
      !=======================================================================
            subroutine wrt_grb2_accu_v2(itau,t0,t1,t2,p3,t10,t11,t12,t24,t27,fld)
            use param, only : nx,my
            use const, only: idtg,RTYPE,ihdgo2
            integer::  t0,t1,t2,t10,t11,p3,t24,t27
            integer::  t12
            integer::ptp0(9)
            integer::itau,ist,istat,i,ia
            real(kind=RTYPE)::fld(nx*my)
            real*4::r4out(nx*my)
!            character:: ihdg2*26,clen*7
             do i = 1, 26
              ia=ichar(ihdgo2(i:i))
              if((ia.ge.97).and.(ia.le.122))then
                ia=ia-32
                ihdgo2(i:i)=char(ia)
              endif
             enddo
 133                      format( A  ,A1 ,A6,A4 )
#ifdef O38K
              write(grbfile,133 )trim(ofdir),'/',ihdgo2(1:6),ihdgo2(13:16)
#else
              write(grbfile,133 )trim(ofdir),'/',ihdgo2(1:6),ihdgo2(11:14)
#endif
              r4out(:)=fld(:)
              call opn_grb2(grbparid,nx,my,idtg,itau,istat)              
              call wrt_grb2_accu_io(itau,t0,t1,t2,p3,t10,t11,t12,t24,t27,r4out)
              call cls_grb2(grbparid,istat)
            end
      !=======================================================================
            !subroutine seclist45(itau,t1,t2,t10,t11,t12,t13,t14,t15,p3,p5)
            !subroutine seclist45(itau,t0,t1,t2,p3,t10,t11,t12,r4out)
            subroutine wrt_grb2_io(itau,t0,t1,t2,p3,t10,t11,t12,fld)
            !itau :tau
      !===varitabel set ===
            !t0   :Product Discipline ( Code Table 0.0 )
            !t1   :Parameter category ( See Code Table 4.1 )
            !t2   :Parameter number ( See Code Table 4.2 )
      !======compress 
            !p3   :Decimal scale factor (D) [[  set varibale to int( var*10**p3 )  ]]
      !======layer set
            !t10  :Type of surface (See Code Table 4.5) 
            !     : 1:GroundSurface 100:prssure leyer(Pa) 101:MeanSeaLevel 103:Height Above Ground(m)
            !t11  :Scale factor         final value is result=t12*10**t11
            !t12  :Value   (ex: T2M  t12=2    850hPa t12=85000 )
            use grib_mod
      !     use mod_typhoon, only: write_mem    !for read ensemble member num
            implicit none
            integer*4::  t0,t1,t2,t10,t11,t12,t13,t14,p3,p5
            real::  t15
            integer::itau
            real*4::fld(grbnxmy)

            listsec0(1)=t0   !Product Discipline ( Code Table 0.0 )
      ! Add data info. (section 4)
      !     ipdsnum=0  ! Product Definition Template Number (Code Table 4.0)
            ipdstmpl(1)=t1   !Parameter category ( See Code Table 4.1 )
            ipdstmpl(2)=t2   !Parameter number ( See Code Table 4.2 )
      !
            if(itau.eq.0)then
            ipdstmpl(3)=0    !Type of generating process (See Code Table 4.3)
            else
            ipdstmpl(3)=2
            endif
            ipdstmpl(9)=itau !Forecast time in units defined by ipdstmpl(8)
            ipdstmpl(10)=t10 !Type of first fixed surface(See Code Table 4.5)
            ipdstmpl(11)=t11 !Scale factor of first fixed surface
            ipdstmpl(12)=t12 !Scaled value of first fixed surface
            ipdstmpl(13)=255!t13
            ipdstmpl(14)=0  !t14
            ipdstmpl(15)=0  !t15
            if(ipdsnum==1)then
             !Type of ensemble forecast (see Code table 4.6)
             !0:hi-res-ctl 1:low-res-ctl 3:member 
             if(grbmem==0)then
               ipdstmpl(16)=1 !control run 
             else 
               ipdstmpl(16)=3 
             endif
             ipdstmpl(17)=grbmem  !Perturbation number
             ipdstmpl(18)=grbnumm !Number of forecasts in ensemble
            endif
      ! Add packing info. (section 5)
           idrstmpl3=(/0,0,p3,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0/)

            call gribcreate(cgrib,lcgrib,listsec0,listsec1,ierr)
            call addgrid(cgrib,lcgrib,igds,igdstmpl,igdstmplen,ideflist,idefnum,ierr)
            call addfield(cgrib,lcgrib,ipdsnum,ipdstmpl,ipdstmplen,    &
                 coordlist,numcoord,idrsnum3,idrstmpl3,idrstmplen3, &
                 fld,grbnxmy,ibmap,bmap,ierr)
            call gribend(cgrib,lcgrib,lengrib,ierr)
            call wryte(grbparid,lengrib,cgrib)
            return
            end subroutine
      !=======================================================================

      !      subroutine seclist4_85(ita,t1,t2,t10,t11,t12,t13,t14,t15,t16,    &
      !                             t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,  &
      !                             t27,t28,t29,p3,p5)
            subroutine wrt_grb2_accu_io(itau,t0,t1,t2,p3,t10,t11,t12,t24,t27,fld)
            !itau :tau
      !===varitabel set ===
            !t0   :Product Discipline ( Code Table 0.0 )
            !t1   :Parameter category ( See Code Table 4.1 )
            !t2   :Parameter number ( See Code Table 4.2 )
      !======compress 
            !p3   :Decimal scale factor (D) [[  set varibale to int( var*10**p3 )  ]]
      !======layer set
            !t10  :Type of surface (See Code Table 4.5) 
            !     : 1:GroundSurface 100:prssure leyer(Pa) 101:MeanSeaLevel 103:Height Above Ground(m)
            !t11  :Scale factor         final value is result=t12*10**t11
            !t12  :Value   (ex: T2M  t12=2    850hPa t12=85000 )
            !t24  ;0 ave 1 accu
            !t27  :average of time range 
            use grib_mod
            integer itau,t0,t1,t2,t10,t11,t12,t13,t14,t16,t17,t18,t19,t20,t21,   &
                        t22,t23,t24,t25,t26,t27,t28,t29,p3,p5
            real t15
            real*4::fld(grbnxmy)
            integer*8::idtg2
            character:: cdtg*12

            if( (itau-t27)  <  0 )return
            listsec0(1)=t0   !Product Discipline ( Code Table 0.0 )
      ! Add data info. (section 4)
      !     ipdsnum_acc=8  ! Product Definition Template Number (Code Table 4.0)
            ipdstmpl8(1)=t1   !Parameter category ( See Code Table 4.1 )
            ipdstmpl8(2)=t2   !Parameter number ( See Code Table 4.2 )
            if(itau.eq.0)then
            ipdstmpl8(3)=0    !Type of generating process (See Code Table 4.3)
            else
            ipdstmpl8(3)=2
            endif
            ipdstmpl8(8 )=1   !Forecast time in units defined by ipdstmpl(8)
            ipdstmpl8(9 )= itau-t27  !Forecast time in units defined by ipdstmpl(8)
            ipdstmpl8(10)=t10 !Type of first fixed surface(See Code Table 4.5)
            ipdstmpl8(11)=t11 !Scale factor of first fixed surface
            ipdstmpl8(12)=t12 !Scaled value of first fixed surface
            ipdstmpl8(13)=255!t13
            ipdstmpl8(14)=0  !t14
            ipdstmpl8(15)=0  !t15
            call dtgfix12(grb_idtg,idtg2 ,itau)
            write(cdtg,'(i12.12)')idtg2
            if(ipdsnum_acc==8)then
              read(cdtg,'(i4,i2,i2,i2)')ipdstmpl8(16:19)!yyyy mm dd hh
              ipdstmpl8(20)=0 !t20 ! minute
              ipdstmpl8(21)=0 !t21 ! Second | Time of end of overall time interval
              ipdstmpl8(22)=1 !t22
              ipdstmpl8(23)=0 !t23
              ipdstmpl8(24)=t24 !t24  0 ave 1 accu
              ipdstmpl8(25)=2 !t25 code 4.11
              ipdstmpl8(26)=1 !t26 ! unit of time for time range. 1:hour, 2:day
              ipdstmpl8(27)=t27  ! time interval value
              ipdstmpl8(28)=0  !t28 fcst,dt unit 0:minut 1:hour 2:day
              ipdstmpl8(29)=0  !t29
            elseif(ipdsnum_acc==11)then
             if(grbmem==0)then
               ipdstmpl8(16)=1 !control run 
             else 
               ipdstmpl8(16)=3 !Type of ensemble forecast (see Code table 4.6)
             endif
             ipdstmpl8(17)=grbmem  !Perturbation number
             ipdstmpl8(18)=grbnumm !Number of forecasts in ensemble
              read(cdtg,'(i4,i2,i2,i2)')ipdstmpl8(19:22)!yyyy mm dd hh
              ipdstmpl8(23)=0 !t20 ! minute
              ipdstmpl8(24)=0 !t21 ! Second | Time of end of overall time interval
              ipdstmpl8(25)=1 !t22
              ipdstmpl8(26)=0 !t23
              ipdstmpl8(27)=t24 !t24  0 ave 1 accu
              ipdstmpl8(28)=2 !t25 code 4.11
              ipdstmpl8(29)=1 !t26 ! unit of time for time range. 1:hour, 2:day
              ipdstmpl8(30)=t27  ! time interval value
              ipdstmpl8(31)=0  !t28 fcst,dt unit 0:minut 1:hour 2:day
              ipdstmpl8(32)=0  !t29
            endif
      ! Add packing info. (section 5) 
           idrstmpl3=(/0,0,p3,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0/)

            call gribcreate(cgrib,lcgrib,listsec0,listsec1,ierr)
            call addgrid(cgrib,lcgrib,igds,igdstmpl,igdstmplen,ideflist,idefnum,ierr)
            call addfield(cgrib,lcgrib,ipdsnum_acc,ipdstmpl8,ipdstmplen8,    &
                 coordlist,numcoord,idrsnum3,idrstmpl3,idrstmplen3, &
                 fld,grbnxmy,ibmap,bmap,ierr)
            call gribend(cgrib,lcgrib,lengrib,ierr)
            call wryte(grbparid,lengrib,cgrib)
            return
            !end
            end subroutine
      !=======================================================================
            !subroutine latlong(im,jm,la2,lo2,din) !gaussian grid
            subroutine latlong(im,jm)!,la2,lo2,din) !gaussian grid
            !for golbal grid (start at 0E,90S)

            integer*4 im,jm,xxyy,la2,lo2,din

            din=360000000/im
      !     la2=90000000-din
            la2=90000000
            lo2=360000000-din
      !     lo2=360000000
      !
             xxyy=im*jm
             grbnxmy=xxyy
      ! Add GRIB2 dimention info. (section 3)
             igds(2)=xxyy  !Number of grid points in the defined grid.
             igds(5)=40    !Gaussian latitude/longitude
             igdstmpl(8)=im    !Ni  - number of points of a parallel
             igdstmpl(9)=jm    !Nj  - number of points of a meridian
             igdstmpl(15)=la2  !La2 - latitude of last grid point
             igdstmpl(16)=lo2  !Lo2 - longitude of last grid point
             igdstmpl(17)=din  !Di  - i direction increment
      !      igdstmpl(18)=jm   !N   - number of parallels between a pole and
             igdstmpl(18)=jm/2 !N   - number of parallels between a pole and the equator

           !CWB define processID 
           if( im < 1560)then  !GEPS is 1552 GI0G
            ipdstmpl (5)=53
            ipdstmpl8(5)=53
           else !GFS is GJ0G
            ipdstmpl (5)=34
            ipdstmpl8(5)=34
           endif           
          
            return
            end subroutine
      !=======================================================================
            !subroutine latlonc(im,jm,la2,lo2,din)
            subroutine latlonc(im,jm)!,la2,lo2,din) !lat-lon grid
            !for golbal grid (start at 0E,90S)

            integer im,jm,xxyy,la2,lo2,din

            din=360000000/im
      !     la2=90000000-din
            la2=90000000
            lo2=360000000-din
      !     lo2=360000000

             xxyy=im*jm
             grbnxmy=xxyy
      ! Add GRIB2 dimention info. (section 3)
             igds(2)=xxyy  !Number of grid points in the defined grid.
             igds(5)=0     !Equal latitude/longitude
             igdstmpl(8)=im    !Ni  - number of points of a parallel
             igdstmpl(9)=jm    !Nj  - number of points of a meridian
             igdstmpl(15)=la2  !La2 - latitude of last grid point
             igdstmpl(16)=lo2  !Lo2 - longitude of last grid point
             igdstmpl(17)=din  !Di  - i direction increment
             igdstmpl(18)=din  !Di  - j direction increment
      !       igdstmpl(18)=jm   !N   - number of parallels between a pole and the equator

            return
            end subroutine
      !=======================================================================
            subroutine opn_grb2(grbidin,nx,my,idtg,itau ,ierr)
            use grib_mod
                 integer::nx,my,itau,ierr,grbidin
                 integer*8::idtg
                 if( .not. allocated( cgrib ) )allocate( cgrib(lcgrib) , bmap(nx*my) )
                 lengrib=0
                 cgrib=''
                 call baopenw(grbidin,trim(grbfile),ierr)
                 call latlong(nx,my)
                 call seclist01(idtg,itau)
                 if(grbmem>=0)then
                   ipdsnum=1
                   ipdsnum_acc=11
                 endif
            end subroutine 
            subroutine cls_grb2 (grbidin,ierr)
            use grib_mod
              integer::ierr,grbidin
              call baclose(grbidin,ierr)
      !       deallocate(cgrib,bmap)
            end subroutine 

      end module mod_grb2_param

