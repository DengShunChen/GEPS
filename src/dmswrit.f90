      subroutine dmswrit(nx,my,lrec,lenc,kflag,ifile,z,istat)
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
!     use index
      use const, only : OutR4key 

      implicit  none

      integer   nx,my,lenc,istat
      logical   t_flg
      real      z(nx,my)
      character lrec*26,ifile*80,kflag*1
      real*4,allocatable:: r4out(:,:)
      character::rflag*1
!
! working array
!
      character key*34
!   
      rflag=kflag
      if (OutR4key)rflag='R' 
      write(key,1000)lrec,rflag,lenc
 1000 format(a26,a1,i7.7)
!
      t_flg=.false.

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

       if(myrank .eq. 0) then
         if(OutR4key)then
           allocate(r4out(nx,my))
           r4out(:,:)=z(:,:)
           call dmsput(ifile,key//char(0),r4out,istat)
           deallocate(r4out)
         else
           call dmsput(ifile,key//char(0),z,istat)
        endif
       t_flg=.true.
       endif
 
!ch    call mpe_broadcast(istat,1,t_flg,mpe_integer)
       call mpe_bcast(istat,1,0,mpe_integer)
!
       if(istat.ne.0)then
         if(myrank .eq. 0)print *,'dmsput key=',key,' error'
         call mpe_finalize
         call dmsexit(-1)
       else
#ifdef VERBOSE
         if(myrank .eq. 0) print *,'dmsput key=',key,' ok'
#endif
       endif

!       endif
!
      return
      end

!---------------------------------------------------------
!CWB2016
! write dmsdata by rank 0, bypassing io_quilting server

      subroutine dmswrit_mfc(nx,my,lrec,lenc,kflag,ifile,z,istat)
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
!     use index

      implicit  none

      integer   nx,my,lenc,istat
      logical   t_flg
      real      z(nx,my)
      character lrec*26,ifile*80,kflag*1
!
! working array
!
      character key*34
!
      write(key,1000)lrec,kflag,lenc
 1000 format(a26,a1,i7.7)
!
      t_flg=.false.

       if(myrank .eq. 0) then
        call dmsput(ifile,key//char(0),z,istat)
        t_flg=.true.
       endif
 
!ch    call mpe_broadcast(istat,1,t_flg,mpe_integer)
       call mpe_bcast(istat,1,0,mpe_integer)
!
       if(istat.ne.0)then
         if(myrank .eq. 0)print *,'dmsput key=',key,' error'
         call mpe_finalize
         call dmsexit(-1)
       else

#ifdef VERBOSE
         if(myrank .eq. 0) print *,'dmsput key=',key,' ok'
#endif
       endif

!
      return
      end
!---------------------------------------------------------
!CWB2016
! write dmsdata by rank 0, bypassing io_quilting server

      subroutine dmswrit_split(nx,my,lrec,lenc,kflag,ifile,z,istat)
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
      use const, only : OutR4key 

      implicit  none

      integer   nx,my,lenc,istat
      logical   t_flg
      real      z(nx,my)
      character lrec*26,ifile*80,kflag*1
      real*4,allocatable:: r4out(:,:)
      character::rflag*1
!
! working array
!
      character key*34
!
!
      rflag=kflag
      if (OutR4key)rflag='R' 
      write(key,1000)lrec,rflag,lenc
 1000 format(a26,a1,i7.7)
!
      t_flg=.false.

!       if(myrank .eq. iroot) then
      if(OutR4key)then
        allocate(r4out(nx,my))
        r4out(:,:)=z(:,:)
        call dmsput(ifile,key//char(0),r4out,istat)
        deallocate(r4out)
      else
        call dmsput(ifile,key//char(0),z,istat)
      endif
       t_flg=.true.
!       endif
 
!ch    call mpe_broadcast(istat,1,t_flg,mpe_integer)
!       call mpe_bcast_col(istat,1,0,mpe_integer)
!
       if(istat.ne.0)then
!!         if(col_rank .eq. 0)print *,'dmsput key=',key,' error'
         print *,'dmsput key=',key,' error'
!!         call mpe_finalize
         call dmsexit(-1)
       else
!!         if(col_rank .eq. 0) print *,'dmsput key=',key,' ok'
#ifdef VERBOSE
         print *,'dmsput key=',key,' ok'
#endif
       endif

!
      return
      end

