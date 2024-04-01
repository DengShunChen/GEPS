      subroutine ioserver(nxmy)

!CWB2016 for gfs io_quilting

      use const, only : ifilout,KLENO,keyo

      implicit none

      integer, parameter :: NMAX=2048

      integer  nxmy,i,j,ist,ist1,ist2,iL,Len,ncnt,ntag
      real*8   z(nxmy,NMAX)
      character(len=KLENO) keys(NMAX)

!CWB2016
      integer  ifromtau,itotau,istat
      ifromtau=0

      call dmsmsg("ALL",ist)
      call dmsopn(ifilout,"w",ist)
      if(ist.ne.0)then
        print*,'ioserver: dmsfile open error, dmsfile and status=',trim(ifilout),ist
        call dmsexit (1)
        stop
      endif

      ntag=0
      ncnt=0

      do while (.true.)
        ntag=ntag+1
        call mpe_recv_key(keyo,ntag,ist1)
        if(keyo(1:4).eq."DONE")goto 100
        if(keyo(1:4).eq."DOIT")then

!CWB2018 bug fixed
        if(ncnt.gt.0)then

          do j=1,ncnt
             call dmsput(ifilout,keys(j)//char(0),z(1,j),ist)
             if(ist.ne.0)then
               print *,' '
               print *,'** Bad I/O: ',keys(j),' dmsput error:',trim(ifilout)
               call dmscls (ifilout,ist)
               stop
             endif
          enddo

!CWB2016
#ifdef O38K
          read(keys(ncnt)(7:12),'(i6)')itotau
#else
          read(keys(ncnt)(7:10),'(i4)')itotau
#endif

!CWB20160927 for NWP control
          if(itotau == 9) call sleep(20)

          if((itotau /= 0).and.(itotau /= ifromtau))then
          call sendmsg ('gfs',ifromtau,itotau,istat)
          if (istat.eq.-1) then
              print *,'ioserver:  SENDMSG ERROR '
              call dmsexit(-1)
          endif
          ifromtau=itotau
          endif

          ncnt=0

!CWB2018 bug fixed
        endif

        else
          ncnt=ncnt+1
          if(ncnt.gt.NMAX)then
            print *,'<<< ioserver : flushing dms buffer to prevent overflow >>>'
            do j=1,NMAX
               call dmsput(ifilout,keys(j)//char(0),z(1,j),ist)
               if(ist.ne.0)then
                 print *,' '
                 print *,'** Bad I/O: ',keys(j),' dmsput error:',trim(ifilout)
                 call dmscls (ifilout,ist)
                 stop
               endif
            enddo
            ncnt=1
          endif

          keys(ncnt)=keyo
          ntag=ntag+1
          call mpe_recv_data(z(1,ncnt),nxmy,ntag,ist2)
        endif
      enddo

100   continue
      call dmscls (ifilout,ist)

      call mpe_finalize
      call dmsexit(0)
      stop

      return
      end subroutine ioserver

!-----------------------------

subroutine ioserver_grb2(nx,my)
use const, only : ifilout_grb,RTYPE,keyo,KLENO2,cleno
use mod_grb2_param
implicit none
integer :: NMAX=500
integer ::nx,my,nxmy,i,j,ist,ncnt,ntag,itau,istat
integer*8::idtg
!integer::ptp0(9,NMAX)
!real*8   z(nx*my,NMAX)
real(kind=RTYPE)::fld(nx*my)
integer,allocatable::ptp0(:,:)
real*4,allocatable ::z(:,:)
integer::t12
integer  ifromtau,itotau,istat
ifromtau=0
nxmy=nx*my

!        memery GB                          core       real-4 
!NMAX=  (   26.    * 1024. * 1024. * 1024. / 48. ) / ( 4.  *  nxmy )
!print*,'in ioserver NMAX=',NMAX
allocate( ptp0(9,NMAX) , z(nx*my,NMAX) )

grbid=233  
ntag=0
 133                      format( A  ,A ,I10.10 ,A , i4.4 ,A      )
do while (.true.)
  ntag=ntag+1
  call mpe_recv_key(keyo,ntag,ist)
  if(keyo(1:4).eq."DONE")exit !goto 100
  if(keyo(1:4).eq."OPEN")then
    read(keyo(cleno:KLENO2),'(I12)')idtg
#ifdef O38K
    read(keyo(7:12),'(I6)')itau
#else
    read(keyo(7:10),'(I4)')itau
#endif
    write(grbfile,133 )trim(ifilout_grb),'/GFS_',idtg/100 ,'_',itau,'.grb2'
    print*,'OutFileName= ',trim(grbfile)
    call opn_grb2(grbid,nx,my,idtg,itau,istat)

    ncnt=0
    ntag=ntag+1
    call mpe_recv_key(keyo,ntag,ist)
    do while ( keyo(1:4)=='DOIT' ) 

      do while ( ncnt < NMAX .and. keyo(1:4)=='DOIT' )
        ncnt=ncnt+1
        ntag=ntag+1
        call mpe_recv_int( ptp0(:,ncnt) , 9 ,ntag,ist)
        ntag=ntag+1
        call mpe_recv_data(fld(:),nx*my,ntag,ist)
        z(:,ncnt)=fld(:)

        ntag=ntag+1
        call mpe_recv_key(keyo,ntag,ist)
        !if(key(1:4)=="CLSE")goto 101
      enddo
      101   continue
      do j=1,ncnt
        t12=ptp0(7,j)
        if( ptp0(8,j)==-999 .and. ptp0(9,j)==-999 )then
          call wrt_grb2_io(itau,ptp0(1,j),ptp0(2,j),ptp0(3,j) &
                  ,ptp0(4,j),ptp0(5,j) ,ptp0(6,j),t12, z(:,j) )
        else
          call wrt_grb2_accu_io(itau,ptp0(1,j),ptp0(2,j),ptp0(3,j),ptp0(4,j), &
                   ptp0(5,j),ptp0(6,j),t12,ptp0(8,j),ptp0(9,j), z(:,j) )
        endif
      enddo
      ncnt=0
    end do ! while ( key(1:4)='DOIT' )

    call cls_grb2(grbid,istat)

    !CWB20160927 for NWP control
#ifdef O38K
    read(keyo(7:12),'(i6)')itotau
#else
    read(keyo(7:10),'(i4)')itotau
#endif
    !if(itotau == 9) call sleep(20)
    if((itotau /= 0).and.(itotau /= ifromtau))then
    call sendmsg ('gfs',ifromtau,itotau,istat)
    if (istat.eq.-1) then
        print *,'ioserver:  SENDMSG ERROR '
        call dmsexit(-1)
    endif
    ifromtau=itotau
    endif

  endif ! (key(1:4).eq."OPEN")
enddo

100   continue
deallocate( z , ptp0 ) 
call mpe_finalize
call dmsexit(0)
stop
return
end subroutine ioserver_grb2




