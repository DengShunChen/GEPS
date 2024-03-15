      subroutine dmswrit(nx,my,lenc,kflag,z,istat)

!CWB2021 single precision test, writing dms output in 32 bits float format

!
!  subroutine to read data in pressure level fields
!
! **** input ****
!
!  lrec: field identifaction label
!  lenc: no. of data values in a 2-d lrec field
!  ifile: path/file name containing initial fields
!
! **** output ****
!
!  z: array containing gaussian grid 2-d field
!  istat: status code. .ne. zero means bad read
!
      use param, only : io_quilting
      use mpe
      use rank
      use const, only : RTYPE,ifilout,keyo,ihdgo
!     use index

      implicit  none

      integer   nx,my,lenc,istat
      logical   t_flg
!CWB2021
      real(kind=RTYPE) z(nx,my)
      character kflag*1
!
! working array
!
!
      write(keyo,1000)ihdgo,kflag,lenc
#ifdef O38K
 1000 format(a28,a1,i9.9)
#else
 1000 format(a26,a1,i7.7)
#endif
!
      t_flg=.false.

!< remove io_quilting
!
!      if(io_quilting)then
!
!        if(myrank .eq. 0) then
!          ntag=ntag+1
!          call mpe_send_key(key,ntag,istat)
!          ntag=ntag+1
!          call mpe_send_data(z,nx*my,ntag,istat)
!#ifdef VERBOSE
!          print *,'dmsput key=',key,' ok'
!#endif
!        endif
!
!      else
!>

       if(myrank .eq. 0) then
!CWB2021
!       if(key(27:27).eq.'R')then
!          z4=z
!          call dmsput(ifile,key//char(0),z4,istat)
!       endif
!       if(key(27:27).eq.'H')then
          call dmsput(ifilout,keyo//char(0),z,istat)
!       endif
       t_flg=.true.
       endif
 
       call mpe_bcast(istat,1,0,mpe_integer)
!
       if(istat.ne.0)then
         if(myrank .eq. 0)print *,'dmsput key=',keyo,' error'
         call mpe_finalize
         call dmsexit(-1)
       else
#ifdef VERBOSE
         if(myrank .eq. 0) print *,'dmsput key=',keyo,' ok'
#endif
       endif

!       endif  !remove io_quilting
!
      return
      end

!---------------------------------------------------------
!CWB2016
! write dmsdata by rank 0, bypassing io_quilting server

      subroutine dmswrit_mfc(nx,my,lenc,kflag,z,istat)
!
!  subroutine to read data in pressure level fields
!
! **** input ****
!
!  lrec: field identifaction label
!  lenc: no. of data values in a 2-d lrec field
!  ifile: path/file name containing initial fields
!
! **** output ****
!
!  z: array containing gaussian grid 2-d field
!  istat: status code. .ne. zero means bad read
!
      use param, only : io_quilting
      use mpe
      use rank
      use const, only : RTYPE,ifilout,keyo,ihdgo
!     use index

      implicit  none

      integer   nx,my,lenc,istat
      logical   t_flg
      real(kind=RTYPE) z(nx,my)
!CWB2021
      character kflag*1
!
! working array
!
!
      write(keyo,1000)ihdgo,kflag,lenc
#ifdef O38K
 1000 format(a28,a1,i9.9)
#else
 1000 format(a26,a1,i7.7)
#endif
!
      t_flg=.false.

       if(myrank .eq. 0) then
!CWB2021
!       if(key(27:27).eq.'R')then
!          z4=z
!          call dmsput(ifile,key//char(0),z4,istat)
!       endif
!       if(key(27:27).eq.'H')then
          call dmsput(ifilout,keyo//char(0),z,istat)
!       endif
       t_flg=.true.
       endif
 
       call mpe_bcast(istat,1,0,mpe_integer)
!
       if(istat.ne.0)then
         if(myrank .eq. 0)print *,'dmsput key=',keyo,' error'
         call mpe_finalize
         call dmsexit(-1)
       else

#ifdef VERBOSE
         if(myrank .eq. 0) print *,'dmsput key=',keyo,' ok'
#endif
       endif

!
      return
      end
!---------------------------------------------------------
!CWB2016
! write dmsdata by rank 0, bypassing io_quilting server

      subroutine dmswrit_split(nx,my,lenc,kflag,z,istat)
!
!  subroutine to read data in pressure level fields
!
! **** input ****
!
!  lrec: field identifaction label
!  lenc: no. of data values in a 2-d lrec field
!  ifile: path/file name containing initial fields
!
! **** output ****
!
!  z: array containing gaussian grid 2-d field
!  istat: status code. .ne. zero means bad read
!
      use param, only : io_quilting
      use mpe
      use rank
      use index, only : col_rank
      use const, only : RTYPE,ifilout,keyo,ihdgo2

      implicit  none

      integer   nx,my,lenc,istat
      logical   t_flg
      real(kind=RTYPE) z(nx,my)
!CWB2021
      character kflag*1
!
! working array
!
!
      write(keyo,1000)ihdgo2,kflag,lenc
#ifdef O38K
 1000 format(a28,a1,i9.9)
#else
 1000 format(a26,a1,i7.7)
#endif
!
      t_flg=.false.

!       if(myrank .eq. iroot) then
!CWB2021
!       if(key(27:27).eq.'R')then
!          z4=z
!          call dmsput(ifile,key//char(0),z4,istat)
!       endif
!       if(key(27:27).eq.'H')then
          call dmsput(ifilout,keyo//char(0),z,istat)
!       endif
       t_flg=.true.
!       endif
 
!       call mpe_bcast_col(istat,1,0,mpe_integer)
!
       if(istat.ne.0)then
!!         if(col_rank .eq. 0)print *,'dmsput key=',key,' error'
         print *,'dmsput key=',keyo,' error'
!!         call mpe_finalize
         call dmsexit(-1)
       else
!!         if(col_rank .eq. 0) print *,'dmsput key=',key,' ok'
#ifdef VERBOSE
         print *,'dmsput key=',keyo,' ok'
#endif
       endif

!
      return
      end

