      subroutine ioserver(nxmy)

!CWB2016 for gfs io_quilting

      use const, only : ifilout

      implicit none

      integer, parameter :: NMAX=2048

      integer  nxmy,i,j,ist,ist1,ist2,iL,Len,ncnt,ntag
      real*8   z(nxmy,NMAX)
      character*34 key,keys(NMAX)

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
        call mpe_recv_key(key,ntag,ist1)
        if(key(1:4).eq."DONE")goto 100
        if(key(1:4).eq."DOIT")then

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
          read(keys(ncnt)(7:10),'(i4)')itotau

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

          keys(ncnt)=key
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
      end
