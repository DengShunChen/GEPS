subroutine tracking(tau,dt_trk,dt,nx,my,slp,v850,v700,h850,h500, &
                    ntyph,typname,ixtyp,jytyp,tlon,tlat,tflon,tflat,idtg, &
                    nrec,typhoon,tensity)
!---------------------------------------------------------------------------!
!  This subroutine is doing to find posiotn of vortex center
!
!  Input :
!       tau     : forcast time 
!       dt_trk  : tracking frequency  
!       dt      : delta t in seconds  
!       nx      : x-asixs of model points  
!       ny      : y-asixs of model points 
!       slp     : Sea level pressure filed 
!       v850    : 850 hPa vorticity
!       v700    : 700 hPa vorticity
!       h850    : 850 hPa geopotential height
!       h500    : 500 hPa geopotential height
!       ntyph   : numbers of typhoon 
!       typname : typhoon name 
!       idtg    : date time 
!       tlon    : lon valuse at model grid point      
!       tlat    : lat valuse at model grid point
!
!  Output :
!       tflon   : typhoon center longitude        
!       tflat   : typhoon center latitude 
!
!---------------------------------------------------------------------------!
!  Thanks to TWRF team, CWB for Technical supporting. 
!
!                                       Created by Chen, Jen-Her
!                                       Modified by Chen, Deng-Shun
!                                       Date : 7 Aug, 2014
!       
!  log:
!   2014-08-27  Chen, Deng-Shun   Be able to output every 6-hour 
!   2014-09-02  Chen, Deng-Shun   Changed track data format 
!   2014-10-13  Chen, Deng-Shun   Changed search domain size from 4 grid points to 6 grid points,
!                                 in order to consider the high moving speed within mid-lattitude area.
!   2015-02-06  Chen, Deng-Shun   Bug fixed. This bugs will cause run-time error
!   2015-11-10  Chen, Deng-Shun   Changed from Fortran 77 code to 90
!
!---------------------------------------------------------------------------!
!
  use mpe
  use rank
  use index
  use mod_typhoon,only:write_mem
  use const,only:ifilout
!  use mod_outflds,only:ifilout
!  use param
  implicit none

! parameters     
  integer, parameter   :: ntau=168, nvar=5, write_tau=6
  real, parameter      :: undef=-99.999
  real, parameter      :: nodata=99999.
  integer,parameter    :: datalength=25000
!--- 
  integer :: nx,my, ntyph, ndt, write_mem2
  character*2 :: mem
  real :: tau,dt_trk,dt 

  integer :: ixtyp(nvar,ntyph),jytyp(nvar,ntyph),nrec(ntyph)
  real :: slp(nx,my),v850(nx,my),v700(nx,my),h850(nx,my),h500(nx,my)
  real :: field(nx,my,nvar)
  real :: tlon(nx,my),tlat(my)

  real :: tflon(0:ntau,nvar,ntyph),tflat(0:ntau,nvar,ntyph)
  real :: tclat(ntyph),tclon(ntyph)
  real :: p6lat(ntyph),p6lon(ntyph),tcslp(ntyph)
  real :: smxv(ntyph),smr30(ntyph),smr50(ntyph)
  real :: slat(0:ntau),slon(0:ntau)
  real :: tensity(0:ntau,nvar,ntyph) 
  real :: rixtyp(nvar),rjytyp(nvar)
  
  logical :: typhoon
  character(15)  :: typname(ntyph)
  integer(8) idtg,idtg8
  character*80 trkpath
  character*150 trkfilename
!---
  logical :: lfound(nvar,ntyph)
  logical :: WriteTrack=.false. ,DoFindTrack=.false.
  character(15) :: trackfile,tensfile
  character*100 line, headerline, tauline
  integer ityp, ist
  integer il,nty,n_h_lat,n_h_lon
  logical l_loop
  integer istat,a
!
  integer itauo(0:ntau,10)  ! The array size should be specified use SAVE is applied 
  save itauo
   
  real  dtx_tau,dtaup
  integer n,i,j,nc,ip, nt, ntime
!DMS variable
      character dmsdb*10
      character dfile*60
      character cdtg*10,idtgc*12
      character epsno*4
      character epstype*1
      character work(datalength)*16,tytrack*9
      character dmshead*3
      character dmstail*8
      integer nstm ! the number of forecasted typhoon
      character domain1*16
      character dmskeytrack*34
!    
      data tytrack /'TYTRACKGT'/
      data domain1 /'CWB GFS  T511L60'/
      data dmstail/'X0025000'/
      data dmsdb/'test'/
      data epsno/'00'/

!---------------------------------------------------------------------------c
! Using five fields to the TC center
  field(:,:,1) =  slp(:,:)      ! sea level pressure
  field(:,:,2) = v850(:,:)      ! 850hPa vorticity 
  field(:,:,3) = v700(:,:)      ! 700hPa vorticity
  field(:,:,4) = h850(:,:)      ! 850hPa height
  field(:,:,5) = h500(:,:)      ! 500hPa height
  
! initialize
  lfound=.false.
  dtx_tau=dt/3600.
  ndt=int((6./dt_trk)+0.001)

! checking write track time
!  dtaup= mod(tau+0.001, real(write_tau))
  dtaup= mod(tau+0.001, 6.)   ! every 6 hrs will output
  WriteTrack=(dtaup .lt. dtx_tau)
!  print*,'WriteTrack',dtaup,dtx_tau,WriteTrack

! do tracking
  dtaup= mod(tau+0.001, tau)
  DoFindTrack=(dtaup .lt. dtx_tau)
!  print*,'DoFindTrack',dtaup,dtx_tau,DoFindTrack
  if(myrank.eq.0)print *,' in tracking tau=',tau

  if(DoFindTrack)then
  !for multi-typhoon 
    do n=1,ntyph

      itauo(0,n)=0
      nrec(n)=nrec(n)+1
      nc=nrec(n)
      itauo(nc,n)=int(tau+0.01)
!    
      do ip=1,nvar
        call findtrk(field(:,:,ip),nx,my,ixtyp(ip,n),jytyp(ip,n),rixtyp(ip),rjytyp(ip),tlon,tlat,ip,lfound(ip,n))
      enddo  
!
      do ip=1,nvar
        if(lfound(ip,n))then
          call xy2ll(rixtyp(ip),rjytyp(ip),tflon(nc,ip,n),tflat(nc,ip,n),tlon,tlat,nx,my)
          i=ixtyp(ip,n) ; j=jytyp(ip,n) 
          if(ip .eq. 4) then 
            tensity(nc,ip,n)=field(i,j,ip)+1457.0 
          elseif(ip .eq. 5) then
            tensity(nc,ip,n)=field(i,j,ip)+5574.0
          else
            tensity(nc,ip,n)=field(i,j,ip)
          endif
        else
          tflon(nc,ip,n)=undef ; tflat(nc,ip,n)=undef
          tensity(nc,ip,n)=undef
        endif
      enddo
!
      if(myrank.eq.0)then
        print *,' tau=',tau,' typhoon=',typname(n),' nrec=',nrec(n)
        print *,' lfound=',lfound(1,n),' fcst slp =',tensity(nc,1,n)
        print *,' at ',tflon(nc,1,n),tflat(nc,1,n)
        print *,' itauo=',itauo(nc,n)
      endif
    enddo ! end of do n typhoons 
!
    DoFindTrack=.false.
  endif
!
!
     dfile=ifilout

     write(idtgc,'(i12)') idtg
     write(mem,'(i2.2)') write_mem
     cdtg=idtgc(1:10)
 if(WriteTrack)then
    if(myrank.eq.0)then
     idtg8=(idtg-200000000000)/100
    call dmsmsg('ERR',ist)
      print *,'dmsdb= ',dfile,'  ist= ',ist,'cdtg=',cdtg,'mem=',mem
!
     do il=1,5
      if (il.NE.4) then

      if (il.eq.1) dmshead='SSL'
      if (il.eq.2) dmshead='850'
      if (il.eq.3) dmshead='700'
      if (il.eq.5) dmshead='500'

      do i=1,datalength
        write(work(i),'(i16.0)') nint(nodata*10000)
      enddo

      write(work(1),'(6x,a10)') cdtg
      nstm=ntyph
      write(work(2),'(i16.0)') nstm
      write(work(3),'(a16)') domain1
!
      do nty=1,nstm ! output all the number of forecasted typhoon

        if(nty.eq.1) n=100
        if(nty.eq.2) n=5100
        if(nty.eq.3) n=10100
        if(nty.eq.4) n=15100
!
        p6lat(nty)=nodata
        p6lon(nty)=nodata
        tcslp(nty)=nodata
        smxv(nty)=nodata
        smr30(nty)=nodata
        smr50(nty)=nodata
        tclat(nty)=tflat(0,1,nty)
        if (tclat(nty) .eq. undef) tclat(nty)=99999.
        tclon(nty)=tflon(0,1,nty)
        if (tclon(nty) .eq. undef) tclon(nty)=99999.
        write(work(n+1),'(1x,a15)') typname(nty)
        write(work(n+2),'(i16.0)') nint(p6lat(nty)*10000)
        write(work(n+3),'(i16.0)') nint(p6lon(nty)*10000)
        write(work(n+4),'(i16.0)') nint(tclat(nty)*10000)
        write(work(n+5),'(i16.0)') nint(tclon(nty)*10000)
        write(work(n+6),'(i16.0)') nint(tcslp(nty)*10000)
        write(work(n+7),'(i16.0)') nint(smxv(nty)*10000)
        write(work(n+8),'(i16.0)') nint(smr30(nty)*10000)
        write(work(n+9),'(i16.0)') nint(smr50(nty)*10000)
!
! ... latitude
        if(nty.eq.1) n_h_lat=110
        if(nty.eq.2) n_h_lat=5110
        if(nty.eq.3) n_h_lat=10110
        if(nty.eq.4) n_h_lat=15110
!
        do i=0,nrec(nty),ndt
          slat(i)=tflat(i,il,nty)
          if (slat(i) .eq. undef) slat(i)=99999.
          write(work(n_h_lat),'(i16.0)') nint(slat(i)*10000)
          n_h_lat=n_h_lat+int(dt_trk) ! every 6 hrs output
        enddo
!
! ... longitude
        if(nty.eq.1) n_h_lon=831
        if(nty.eq.2) n_h_lon=5831
        if(nty.eq.3) n_h_lon=10831
        if(nty.eq.4) n_h_lon=15831
!
        do i=0,nrec(nty),ndt
          slon(i)=tflon(i,il,nty)
          if (slon(i) .eq. undef) slon(i)=99999.
            write(work(n_h_lon),'(i16.0)') nint(slon(i)*10000)
            n_h_lon=n_h_lon+int(dt_trk) ! every 6 hrs output
        enddo
!
!      if (trim(epsno) .eq. 'mean') then
!      dmskeytrack=dmshead//tytrack//'MN'//cdtg//'00'//dmstail
!      else
      dmskeytrack=dmshead//tytrack//trim(mem)//cdtg//'00'//dmstail
!      endif

      call dmsput (dfile,dmskeytrack//char(0),work,istat)
      if (istat. ne. 0) then
        print *,'dmskey : ',dmskeytrack,' put error !!'
        print *,'abort the tytrack program !!'
      call dmsexit(1)
      else
!        print *,'dmskey : ',dmskeytrack,' put OK !!'
      endif
!
      enddo

      endif ! end of il=4 
     enddo ! end of il loop
!
!
      if(myrank.eq.0)print *,'idtg,idtg8=',idtg,idtg8

      ! for track file
      write(trackfile,'(a3,i8.8,a4)')'trk',idtg8,'.dat'
      open(15,file=trackfile,form='formatted',status='unknown')
      if(myrank.eq.0)print *,'open typhoon track file=',trackfile
      write(15,'(i8.8,a24,i2)')idtg8,' number of typhoons= ',ntyph

      ! for intensity file
      write(tensfile,'(a3,i8.8,a4)')'ten',idtg8,'.dat'
      open(16,file=tensfile,form='formatted',status='unknown')
      if(myrank.eq.0)print *,'open typhoon intensity file=',tensfile
      write(16,'(i8.8,a24,i2)')idtg8,' number of typhoons= ',ntyph


      do n=1,ntyph
        write(15,"('number = ',i2,' typh-name= ',a15)")n,typname(n)
        write(16,"('number = ',i2,' typh-name= ',a15)")n,typname(n)
        do i=2,5
          tflat(0,i,n)=tflat(0,1,n)
          tflon(0,i,n)=tflon(0,1,n)
        enddo
        do nt=0,nrec(n),ndt
          ntime=itauo(nt,n)
          write(15,1000)ntime,tflat(nt,2,n),tflon(nt,2,n), &
                              tflat(nt,1,n),tflon(nt,1,n), &
                              tflat(nt,4,n),tflon(nt,4,n), &
                              tflat(nt,3,n),tflon(nt,3,n), &
                              tflat(nt,5,n),tflon(nt,5,n)
 1000  format(i3.3,10f8.3)
!
          write(16,1002)ntime,tensity(nt,2,n), &
                              tensity(nt,1,n), &
                              tensity(nt,4,n), &
                              tensity(nt,3,n), &
                              tensity(nt,5,n)
 1002  format(i3.3,1x,e10.4,1x,f7.2,1x,f7.2,1x,e10.4,1x,f8.2)
!
        enddo
      enddo
    endif
!    typhoon=.false.  ! close for every 6 hrs output
    close(15)
    close(16)
  endif
!
return
end
!
subroutine findtrk(fld,nx,my,ix,iy,rx,ry,tlon,tlat,index,lfound)
!-------------------------------------------------------------------------
!       Finding Tropical Cyclone Center Position 
!
!       fld     : fields
!       nx      : total grid points at x-axis 
!       my      : total grid points at y-axis
!       ix      : center position at x-axis with integer 
!       iy      : center position at y-axis with integer
!       rx      : center position at x-axis with float number 
!       ry      : center position at y-axis with float number
!       index   : different field 
!       lfound  : logical flag of finding typhoon
!
!       hx      : Newton's Metod ajustment at x-axis 
!       hy      : Newton's Metod ajustment at y-axis 
! 
!-------------------------------------------------------------------------
!
  use rank
  use mod_typhoon, only : min_trk_pres

  implicit none

  integer :: i,j
  integer :: nx,my,ix,iy,index
  integer :: ib,ie,jb,je
  real    :: fld(nx,my),fldavg
  real    :: tlon(nx,my),tlat(my)
  logical :: lfound

  integer :: ixyrange
  real    :: f0,f1,f2,f3,f4
  real    :: xxx,yyy,hx,hy,rx,ry
  real    :: max_value, min_value
!-------------------------------------------------------------------------
!     max. Typhoon moving speed is around 50 km/h, which in t320(about
!     50 km) resolution, the serch area should be bigger than
!     50*(tracking interval)  
!-------------------------------------------------------------------------
! data ixyrange/4/    ! for 3 hours tracking interval 
! data ixyrange/6/    ! for 3 hours tracking interval, grid points. T320
!  data ixyrange/12/    ! for 3 hours tracking interval, grid points. T512 
! data ixyrange/16/
! data ixyrange/8/    ! river's original setup
!
 if (tlat(iy) .le. 30. ) then
!   ixyrange=8
   ixyrange=12
 else
!   ixyrange=14
   ixyrange=20
 endif
 
 if(myrank.eq.0) print*,'tlat = ',tlat(iy),' ixyrange = ',ixyrange 

!-- set search domain 
  ib=ix-ixyrange
  ie=ix+ixyrange
  jb=iy-ixyrange
  je=iy+ixyrange
!
  if(index.eq.1 .or. index.eq.4 .or. index.eq.5)then
    min_value=99999.
    do j=jb,je
    do i=ib,ie
      fldavg=fld(i,j)
      if(fldavg.lt.min_value)then
        ix=i
        iy=j
        min_value=fldavg
      endif
    enddo
    enddo

    if(myrank.eq.0) then
      write(6,*)'ib/ie/jb/je/ix/iy/min_value/index = ',ib,ie,jb,je,ix,iy,min_value,index
    endif

!  if the vortex is too weak 
    if(index.eq.1) then 
!      if(min_value.le.999. )then  
      if(min_value .le. min_trk_pres)then  
        lfound=.true.
      else
        lfound=.false.       
        return
      endif
!    if(myrank.eq.0) then
!    print*,'min_trk_pres=',min_trk_pres
!    endif
    endif

!--- check center position 
    f0=fld(ix,iy)
    f1=fld(ix+1,iy)
    f2=fld(ix,iy+1)
    f3=fld(ix-1,iy)
    f4=fld(ix,iy-1)

    if ( f0.gt.f1 .or. f0.gt.f2 .or. f0.gt.f3 .or. f0.gt.f4 ) then
      if(myrank.eq.0) then
        write(6,*)'findtrk : f0 greater than surrounding values ,index = ',index
        write(6,*)'ix/iy/f0/f1/f2/f3/f4/index = ',ix,iy,f0,f1,f2,f3,f4,index
      endif 
      lfound=.false.
    else
      xxx=f1-2.*f0+f3
      yyy=f2-2.*f0+f4

      ! x-direction 
      if(xxx.ne.0.) then
        hx = 0.5*(f1-f3)/(f1-2.*f0+f3)
      else
        hx=0.
      endif
   
      ! y-direction 
      if(yyy.ne.0.) then
        hy = 0.5*(f2-f4)/(f2-2.*f0+f4)
      else
        hy=0.
      endif
      rx=real(ix)-hx
      ry=real(iy)-hy
      !ix=int(rx+0.5)
      !iy=int(ry+0.5)
      lfound=.true.
    endif
!------------------------------------------------------------
  else if(index.eq.2 .or. index.eq.3)then
    max_value=-99999.
    do  j=jb,je
    do  i=ib,ie
      fldavg=fld(i,j)
      if(fldavg.gt.max_value)then
        ix=i
        iy=j
        max_value=fldavg
      endif
    enddo
    enddo
    if(myrank.eq.0) then
      write(6,*)'ib/ie/jb/je/ix/iy/max_value/index = ',ib,ie,jb,je,ix,iy,max_value,index
    endif
    !--- check center position
    f0=fld(ix,iy)
    f1=fld(ix+1,iy)
    f2=fld(ix,iy+1)
    f3=fld(ix-1,iy)
    f4=fld(ix,iy-1)
    if ( f0 .lt. f1 .or. f0 .lt. f2 .or. f0 .lt. f3 .or. f0 .lt. f4 ) then
      if(myrank.eq.0) then
        write(6,*)'findtrk : f0 lower than surrounding ,index=',index
        write(6,*)'ix/iy/f0/f1/f2/f3/f4/index = ',ix,iy,f0,f1,f2,f3,f4,index
      endif 
      lfound=.false.
    else
      xxx=f1-2.*f0+f3
      yyy=f2-2.*f0+f4
      ! x-direction 
      if(xxx.ne.0.) then
        hx = 0.5*(f1-f3)/(f1-2.*f0+f3)
      else
        hx=0.
      endif
      ! y-direction 
      if(yyy.ne.0.) then
        hy = 0.5*(f2-f4)/(f2-2.*f0+f4)
      else
        hy=0.
      endif
      rx=real(ix)-hx
      ry=real(iy)-hy
      ix=int(rx+0.5)
      iy=int(ry+0.5)
      lfound=.true.
    endif
  endif
!
return
end subroutine findtrk

subroutine xy2ll (rx,ry,lon,lat,tlon,tlat,nx,my)
!------------------------------------------------------------------c
!     convert model x-y to lat-lon position
!------------------------------------------------------------------c
  use rank   
  implicit none

  integer nx,my,ix,iy
  real rx,ry, dx, dy
  real lon,lat, dlon, dlat
  real tlon(nx,my), tlat(my)
!
  ix=int(rx)
  iy=int(ry)        

  dx=rx-real(ix)
  dy=ry-real(iy)
!
  if (ix+1 .gt. nx) then
    if(myrank.eq.0) print*,'Warning !! ix+1 greater than nx !!'
    dlon=tlon(1,iy)-tlon(ix,iy)
  else
    dlon=tlon(ix+1,iy)-tlon(ix,iy)
  endif

  if (iy+1 .gt. my) then
    if(myrank.eq.0) print*,'Warning !! iy+1 greater than my !!'
    dlat=0.
  else
    dlat=tlat(iy+1)-tlat(iy)
  endif
!
  lon=tlon(ix,iy)+(dx*dlon)
  lat=tlat(iy)+(dy*dlat)
  if(myrank.eq.0) print*,'xy2ll : x/y/lon/lat =',rx,ry,lon,lat 

return 
end subroutine xy2ll
!
