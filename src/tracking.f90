subroutine tracking(tau,dt_trk,dt,nx,my,                                  &
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
!       typtrk  : (:,:,1) Sea level pressure filed 
!               : (:,:,2) 850 hPa vorticity
!               : (:,:,3) 700 hPa vorticity
!               : (:,:,4) 850 hPa geopotential height
!               : (:,:,5) 500 hPa geopotential height
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
!   2019-10-22  Liu, Pang-Yen     Reducing memory
!
!---------------------------------------------------------------------------!
!
  use mpe
  use rank
  use index
  use mod_typhoon,only:write_mem,write_tau,typtrk
  use const,only:ifilout,RTYPE,KLENO
  use param,only:my_max
!  use mod_outflds,only:ifilout
!  use param
  implicit none

! parameters     
  integer, parameter   :: ntau=168, nvar=5
  real, parameter      :: undef=-99.999
  real, parameter      :: nodata=99999.
  integer,parameter    :: datalength=25000
!--- 
  integer :: nx,my, ntyph, ndt, write_mem2
  character*2 :: mem
  real :: tau,dt_trk,dt

  integer :: ixtyp(nvar,ntyph),jytyp(nvar,ntyph),nrec(ntyph)
!byl  real :: slp(nx,my),v850(nx,my),v700(nx,my),h850(nx,my),h500(nx,my)
!byl  real :: field(nx,my,nvar)
  real(kind=RTYPE) :: field(nx,my)
  real :: tlon(nx),tlat(my)

  real :: tflon(0:ntau,nvar,ntyph),tflat(0:ntau,nvar,ntyph)
  real :: tclat(ntyph),tclon(ntyph)
  real :: p6lat(ntyph),p6lon(ntyph),tcslp(ntyph)
  real :: smxv(ntyph),smr30(ntyph),smr50(ntyph)
  real :: slat(0:ntau),slon(0:ntau)
  real :: sten(0:ntau),speed1(0:ntau)
  real :: tensity(0:ntau,nvar,ntyph) 
  real :: rixtyp(nvar),rjytyp(nvar)
  real :: dist,speed,slatpre,slonpre
  logical :: typhoon
  character(15)  :: typname(ntyph)
  integer(8) idtg,idtg8
  character*255 trkpath
  character*150 trkfilename
!---
  logical :: lfound(nvar,ntyph)
  logical :: WriteTrack=.false. ,DoFindTrack=.false.
  character(15) :: trackfile,tensfile
  character*100 line, headerline, tauline
  integer ityp, ist
  integer il,nty,n_h_lat,n_h_lon,n_h_ten,n_h_spd
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
#ifdef O38K
      character work(datalength)*16,tytrack*11
      character dmstail*10
#else
      character work(datalength)*16,tytrack*9 
      character dmstail*8
#endif
      character dmshead*3
      integer nstm ! the number of forecasted typhoon
      character domain1*16
      character(len=KLENO) dmskeytrack
!   
#ifdef O38K 
      data tytrack /'TYPHTRACKGT'/
      data dmstail/'X002500000'/
#else
      data tytrack /'TYTRACKGT'/
      data dmstail/'X0025000'/
#endif
      data domain1 /'CWB GFS  T511L60'/
      data dmsdb/'test'/
      data epsno/'00'/

! for judging the undef value
  real ::  min_trk_pres,distpre
  real, save    :: rixtyp_savep(60,5,10), rjytyp_savep(60,5,10)
  logical, save :: lfound_save(60,5,10)
  integer :: range_le30, range_gt30

!---------------------------------------------------------------------------c
! get resolution dependent tracking range
  !range_le30 = int(50*dt_trk*nx/36000)+2
!20240812 wei
  range_le30 = int(46.0*dt_trk*nx/36000)+2 
  range_gt30 = int(87.5*dt_trk*nx/36000)+2

! Using five fields to the TC center
!byl  field(:,:,1) =  slp(:,:)      ! sea level pressure
!byl  field(:,:,2) = v850(:,:)      ! 850hPa vorticity 
!byl  field(:,:,3) = v700(:,:)      ! 700hPa vorticity
!byl  field(:,:,4) = h850(:,:)      ! 850hPa height
!byl  field(:,:,5) = h500(:,:)      ! 500hPa height
  
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
  !if ( tau < 0.001 )    dtaup=0.0
  
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
      if (nc .le. 72/dt_trk+1) then                                     ! 3-day fcst
        min_trk_pres=1007.
      else if (nc .gt. 72/dt_trk+1 .and. nc .le. 120/dt_trk+1) then     ! 3~5-day fcst
        min_trk_pres=1003.
      else                                                              ! 5~7-day fcst
        min_trk_pres=1000.
      endif

      do ip=1,nvar
        call unify_reduceintp(nx,my,my_max,typtrk(1,1,ip),field)
!byl        call findtrk(field(:,:,ip),nx,my,ixtyp(ip,n),jytyp(ip,n),rixtyp(ip),rjytyp(ip),tlon,tlat,ip,lfound(ip,n),min_trk_pres,range_le30,range_gt30)
        call findtrk(field,nx,my,ixtyp(ip,n),jytyp(ip,n),rixtyp(ip),rjytyp(ip),tlon,tlat,ip,lfound(ip,n),min_trk_pres,range_le30,range_gt30)
!      enddo

! for pressure -999*3 check
      if (nc .ge. 3) then

!      do ip=1,nvar

      if (lfound(ip,n)) then
        call xy2ll(rixtyp(ip),rjytyp(ip),tflon(nc,ip,n),tflat(nc,ip,n),tlon,tlat,nx,my)

        if (lfound_save(nc-1,ip,n)) then

        call xy2ll(rixtyp_savep(nc-1,ip,n),rjytyp_savep(nc-1,ip,n),tflon(nc-1,ip,n),tflat(nc-1,ip,n),tlon,tlat,nx,my)
        call greatcir(tflon(nc,ip,n),tflat(nc,ip,n),tflon(nc-1,ip,n),tflat(nc-1,ip,n),distpre)

!        if(myrank .eq. 0) print*,'dist pre-1 = ',distpre
        if (distpre .gt. dt_trk*85) then
           lfound(ip,n)=.false.
        endif

        else ! (lfound_save(nc-1,1,n) == .false.) then

          if (lfound_save(nc-2,ip,n)) then

            call xy2ll(rixtyp_savep(nc-2,ip,n),rjytyp_savep(nc-2,ip,n),tflon(nc-2,ip,n),tflat(nc-2,ip,n),tlon,tlat,nx,my)
            call greatcir(tflon(nc,ip,n),tflat(nc,ip,n),tflon(nc-2,ip,n),tflat(nc-2,ip,n),distpre)

!            if(myrank .eq. 0)print*,'dist pre-2 = ',distpre
            if (distpre .gt. dt_trk*85*2) then
              lfound(ip,n)=.false.
            endif

          else ! lfound_save(nc-2,1,n) == .false.
            lfound(ip,n)=.false.
          endif
        endif

      endif

!      enddo
      endif ! nc > 3

!      do ip=1,nvar
         rixtyp_savep(nc,ip,n)=rixtyp(ip)
         rjytyp_savep(nc,ip,n)=rjytyp(ip)
         lfound_save(nc,ip,n)=lfound(ip,n)
!      enddo
! for pressure -999*3 check   

!      do ip=1,nvar
        if(lfound(ip,n))then
          call xy2ll(rixtyp(ip),rjytyp(ip),tflon(nc,ip,n),tflat(nc,ip,n),tlon,tlat,nx,my)
          i=ixtyp(ip,n) ; j=jytyp(ip,n) 
          if(ip .eq. 4) then 
            tensity(nc,ip,n)=field(i,j)+1457.0 
          elseif(ip .eq. 5) then
            tensity(nc,ip,n)=field(i,j)+5574.0
          else
            tensity(nc,ip,n)=field(i,j)
          endif
        else
          tflon(nc,ip,n)=undef ; tflat(nc,ip,n)=undef
          tensity(nc,ip,n)=undef
        endif
      enddo ! ip=1,nvar
!
!cjh
      if   (tensity(nc,1,n) .eq. undef) then
        if (tensity(nc,4,n) .ne. undef) then
          if(tensity(nc-1,1,n) .ne. undef )then
           tflon(nc,1,n)=0.5*(tflon(nc-1,1,n)+tflon(nc,4,n))
           tflat(nc,1,n)=0.5*(tflat(nc-1,1,n)+tflat(nc,4,n))
           tensity(nc,1,n)=tensity(nc-1,1,n)
          endif
        elseif (tensity(nc,4,n) .eq. undef) then
          if( tensity(nc,2,n) .ne. undef)  then
            if( abs(tflon(nc,2,n)-tflon(nc-1,1,n)) .lt. 4. .and. &
                abs(tflat(nc,2,n)-tflat(nc-1,1,n)) .lt. 4.) then
             if(tensity(nc-1,4,n) .ne. undef )then
              tflon(nc,1,n)=0.5*(tflon(nc-1,1,n)+tflon(nc,2,n))
              tflat(nc,1,n)=0.5*(tflat(nc-1,1,n)+tflat(nc,2,n))
              tensity(nc,1,n)=tensity(nc-1,1,n)
              tflon(nc,4,n)=0.5*(tflon(nc-1,4,n)+tflon(nc,2,n))
              tflat(nc,4,n)=0.5*(tflat(nc-1,4,n)+tflat(nc,2,n))
              tensity(nc,4,n)=tensity(nc-1,4,n)
             endif
            endif
          endif
        endif
      endif
      if( (tensity(nc,4,n) .eq. undef) .and. (tensity(nc,1,n).ne.undef) )then
       if(tensity(nc-1,4,n) .ne. undef )then
        tflon(nc,4,n)=0.5*(tflon(nc-1,4,n)+tflon(nc,1,n))
        tflat(nc,4,n)=0.5*(tflat(nc-1,4,n)+tflat(nc,1,n))
        tensity(nc,4,n)=tensity(nc-1,4,n)
       endif
      endif
      if( (tensity(nc,2,n) .eq. undef) .and. (tensity(nc,3,n).ne.undef) )then
       if(tensity(nc-1,2,n) .ne. undef )then
        tflon(nc,2,n)=0.5*(tflon(nc-1,2,n)+tflon(nc,3,n))
        tflat(nc,2,n)=0.5*(tflat(nc-1,2,n)+tflat(nc,3,n))
        tensity(nc,2,n)=tensity(nc-1,2,n)
       endif
      endif
      if( (tensity(nc,3,n) .eq. undef) .and. (tensity(nc,2,n).ne.undef) )then
       if(tensity(nc-1,3,n) .ne. undef )then
        tflon(nc,3,n)=0.5*(tflon(nc-1,3,n)+tflon(nc,2,n))
        tflat(nc,3,n)=0.5*(tflat(nc-1,3,n)+tflat(nc,2,n))
        tensity(nc,3,n)=tensity(nc-1,3,n)
       endif
      endif

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
     if(idtg.ge.200000000000_8)then
     idtg8=(idtg-200000000000_8)/100
     else
     idtg8=(idtg-190000000000_8)/100
     endif
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
          n_h_lat=n_h_lat+int(write_tau) ! every 6 hrs output
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
            n_h_lon=n_h_lon+int(write_tau) ! every 6 hrs output
        enddo
!
!.. 6 hours average wind speed
        if(nty.eq.1) n_h_spd=2273
        if(nty.eq.2) n_h_spd=7273
        if(nty.eq.3) n_h_spd=12273
        if(nty.eq.4) n_h_spd=17273

        do i=0,nrec(nty),ndt
          if ((i-ndt) .lt. 0) then
            speed=undef
            dist=undef
          else
            slatpre=tflat(i-ndt,il,nty)
            slonpre=tflon(i-ndt,il,nty)

            if (slon(i) .eq. undef .or. slat(i) .eq. undef .or. &
                slonpre .eq.undef .or. slatpre .eq. undef) then
              speed=undef
              dist=undef
            else
              call greatcir(slonpre,slatpre,slon(i),slat(i),dist)
              speed=dist/int(write_tau) ! every 6 hrs average
            endif
          endif

          if (speed .eq. undef .or. speed .gt. 100.) then
             speed1(i)=99999.
          else 
             speed1(i)=speed
          endif
!          if(myrank .eq. 0) print*,"nty = ",nty,"dist(",i,") = ",dist,"speed(",i,") = ",speed1(i)
            write(work(n_h_spd),'(i16.0)') nint(speed1(i)*10000)
            n_h_spd=n_h_spd+int(write_tau) ! every 6 hrs output
        enddo
!
! ... center min-pressure
        if(nty.eq.1) n_h_ten=1552
        if(nty.eq.2) n_h_ten=6552
        if(nty.eq.3) n_h_ten=11552
        if(nty.eq.4) n_h_ten=16552

        do i=0,nrec(nty),ndt
          sten(i)=tensity(i,1,nty)
          if (sten(i) .eq. undef) sten(i)=99999.
!            print*,"tensity=",sten(i)
            write(work(n_h_ten),'(i16.0)') nint(sten(i)*10000)
            n_h_ten=n_h_ten+int(write_tau) ! every 6 hrs output
        enddo 
               
!
!      if (trim(epsno) .eq. 'mean') then
!      dmskeytrack=dmshead//tytrack//'MN'//cdtg//'00'//dmstail
!      else
      dmskeytrack=dmshead//tytrack//trim(mem)//cdtg//'00'//dmstail
!      endif

      call dmsput (dfile,dmskeytrack//char(0),work,istat)
      if (istat .ne. 0) then
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
subroutine findtrk(fld,nx,my,ix,iy,rx,ry,tlon,tlat,index,lfound,min_trk_pres,range_le30,range_gt30)
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
  use const, only: RTYPE
!  use mod_typhoon, only : min_trk_pres

  implicit none

  integer :: i,j
  integer :: nx,my,ix,iy,index
  integer :: ib,ie,jb,je
  real(kind=RTYPE) :: fld(nx,my)
  real    :: fldavg
  real    :: tlon(nx),tlat(my)
  logical :: lfound

  integer :: ixyrange
  real    :: f0,f1,f2,f3,f4
  real    :: xxx,yyy,hx,hy,rx,ry
  real    :: max_value, min_value, min_trk_pres
  integer :: range_le30, range_gt30

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
 !if (tlat(iy) .le. 30. ) then
 !20240812 wei
 if (tlat(iy) .le. 40. ) then
   ixyrange=range_le30
!   ixyrange=8
 else
   ixyrange=range_gt30
!   ixyrange=14
 endif

! print*, 'range_le30 = ',range_le30
! print*, 'range_gt30 = ',range_gt30
 
 if(myrank.eq.0) print*,'tlat = ',tlat(iy),' ixyrange = ',ixyrange 

!-- set search domain 
  ib=max(ix-ixyrange,1)
  ie=min(ix+ixyrange,nx)
  jb=max(iy-ixyrange,1)
  je=min(iy+ixyrange,my)
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
    f1=fld( min(ix+1,nx) , iy )
    f2=fld( ix , min(iy+1,my) )
    f3=fld( max(ix-1, 1) , iy )
    f4=fld( ix , max(iy-1, 1) )

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
    f1=fld( min(ix+1,nx) , iy )
    f2=fld( ix , min(iy+1,my) )
    f3=fld( max(ix-1, 1) , iy )
    f4=fld( ix , max(iy-1, 1) )

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
  real tlon(nx), tlat(my)
!
  ix=int(rx)
  iy=int(ry)        

  dx=rx-real(ix)
  dy=ry-real(iy)
!
  if (ix+1 .gt. nx) then
    if(myrank.eq.0) print*,'Warning !! ix+1 greater than nx !!'
!    dlon=tlon(1)-tlon(ix)
    if ( ix .eq. nx ) dlon=tlon(ix+1-nx)-tlon(ix)
    if ( ix .gt. nx ) dlon=tlon(ix+1-nx)-tlon(ix-nx)
  else if ( ix .lt. 1 ) then
    if(myrank.eq.0) print*,'Warning !! ix smaller than 1 !!'
    if ( ix+1 .eq. 1  ) dlon=tlon(ix+1)-tlon(ix+nx)
    if ( ix+1 .lt. 1  ) dlon=tlon(ix+1+nx)-tlon(ix+nx)
  else
    dlon=tlon(ix+1)-tlon(ix)
  endif

  if (iy+1 .gt. my) then
    if(myrank.eq.0) print*,'Warning !! iy+1 greater than my !!'
    dlat=0.
  else
    dlat=tlat(iy+1)-tlat(iy)
  endif
!
  lon=tlon(ix)+(dx*dlon)
  lat=tlat(iy)+(dy*dlat)
  if(myrank.eq.0) print*,'xy2ll : x/y/lon/lat =',rx,ry,lon,lat 

return 
end subroutine xy2ll

subroutine greatcir(lon1,lat1,lon2,lat2,dist)
!--------------------------------------------------------------c
!    calculate distance of lon1,lon2 lat1,lat2
!==============================================================c
      integer np
      real dist
      real lon1,lat1,lon2,lat2
      real dlon,dlat
!
      pi=4.0*atan(1.0)
      rad=6370. ! earth radius in km
      torad=pi/180.
      dlon=(lon2-lon1)*torad
      dlat=(lat2-lat1)*torad
      a=sin(dlat/2.)*sin(dlat/2.)+cos(lat1*torad)*cos(lat2*torad)*sin(dlon/2.)*sin(dlon/2.)
      c=2.*atan2(sqrt(a), sqrt(1-a))
      dist=rad*c
!      print *,'The distance between (lon1,lat1) and (lon2,lat2) = ',dist
!      print *,'lon1,lat1,lon2,lat2= ',lon1,lat1,lon2,lat2
!
return
end subroutine greatcir
!
