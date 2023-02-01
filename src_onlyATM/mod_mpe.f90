      module mpe
 
      use rank

      implicit none

      public

      integer, parameter :: mpe_integer=0
      integer, parameter :: mpe_double=1
      integer, parameter :: mpe_logical=2
!CWB 2019-06-25
      integer, parameter :: mpe_single=3

!
      interface mpe_broadcast
         module procedure mpe_broadcast_i_scalar
         module procedure mpe_broadcast_i
         module procedure mpe_broadcast_i_2d
         module procedure mpe_broadcast_r4_scalar
         module procedure mpe_broadcast_r4
         module procedure mpe_broadcast_r4_2d
         module procedure mpe_broadcast_r8_scalar
         module procedure mpe_broadcast_r8
         module procedure mpe_broadcast_r8_2d
      end interface

      interface mpe_bcast
         module procedure mpe_bcast_i_scalar
         module procedure mpe_bcast_i
         module procedure mpe_bcast_i_2d
         module procedure mpe_bcast_r4_scalar
         module procedure mpe_bcast_r4
         module procedure mpe_bcast_r4_2d
         module procedure mpe_bcast_r8_scalar
         module procedure mpe_bcast_r8
         module procedure mpe_bcast_r8_2d
      end interface

      interface mpe_bcast_col
         module procedure mpe_bcast_col_i_scalar
         module procedure mpe_bcast_col_i
         module procedure mpe_bcast_col_i_2d
         module procedure mpe_bcast_col_r4_scalar
         module procedure mpe_bcast_col_r4
         module procedure mpe_bcast_col_r4_2d
         module procedure mpe_bcast_col_r8_scalar
         module procedure mpe_bcast_col_r8
         module procedure mpe_bcast_col_r8_2d
      end interface

      interface mpe_global_sum
         module procedure mpe_global_sum_i_scalar
         module procedure mpe_global_sum_i
         module procedure mpe_global_sum_r4_scalar
         module procedure mpe_global_sum_r4
         module procedure mpe_global_sum_r8_scalar
         module procedure mpe_global_sum_r8
      end interface

      interface mpe_global_max
         module procedure mpe_global_max_i_scalar
         module procedure mpe_global_max_i
         module procedure mpe_global_max_r4_scalar
         module procedure mpe_global_max_r4
         module procedure mpe_global_max_r8_scalar
         module procedure mpe_global_max_r8
      end interface

      interface mpe_global_maxloc
         module procedure mpe_global_maxloc_i
         module procedure mpe_global_maxloc_r4
         module procedure mpe_global_maxloc_r8
      end interface

      interface mpe_global_min
         module procedure mpe_global_min_i_scalar
         module procedure mpe_global_min_i
         module procedure mpe_global_min_r4_scalar
         module procedure mpe_global_min_r4
         module procedure mpe_global_min_r8_scalar
         module procedure mpe_global_min_r8
      end interface

      interface mpe_global_minloc
         module procedure mpe_global_minloc_i
         module procedure mpe_global_minloc_r4
         module procedure mpe_global_minloc_r8
      end interface

      interface mpe_send_print
         module procedure mpe_send_print_i
         module procedure mpe_send_print_r4
         module procedure mpe_send_print_r8
      end interface

      interface mpe_recv_print
         module procedure mpe_recv_print_i
         module procedure mpe_recv_print_r4
         module procedure mpe_recv_print_r8
      end interface
 
      contains

!--------------------------------------------------------------
      subroutine mpe_broadcast_i_scalar(buf,n,flag,type)
 
      use rank, only : myrank,MPI_COMM_gfs
      use mpi

      integer  buf,n,iroot_in,iroot,IERR
      logical  flag
      integer, optional :: type
 
      iroot_in=0
      if (flag) iroot_in = myrank+1
 
      call MPI_ALLREDUCE( iroot_in, iroot, 1, MPI_INTEGER,    &
                          MPI_SUM, MPI_COMM_gfs, IERR )

      call MPI_BCAST( BUF, N, MPI_INTEGER, IROOT-1,           &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_broadcast_i_scalar
!--------------------------------------------------------------
      subroutine mpe_broadcast_i(buf,n,flag,type)
 
      use rank, only : myrank,MPI_COMM_gfs
      use mpi

      integer  buf(n),n,iroot_in,iroot,IERR
      logical  flag
      integer, optional :: type
 
      iroot_in=0
      if (flag) iroot_in = myrank+1
 
      call MPI_ALLREDUCE( iroot_in, iroot, 1, MPI_INTEGER,    &
                          MPI_SUM, MPI_COMM_gfs, IERR )

      call MPI_BCAST( BUF, N, MPI_INTEGER, IROOT-1,           &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_broadcast_i
!--------------------------------------------------------------
      subroutine mpe_broadcast_i_2d(buf,n,flag,type)
 
      use rank, only : myrank,MPI_COMM_gfs
      use mpi

      integer  buf(n,1),n,iroot_in,iroot,IERR
      logical  flag
      integer, optional :: type
 
      iroot_in=0
      if (flag) iroot_in = myrank+1
 
      call MPI_ALLREDUCE( iroot_in, iroot, 1, MPI_INTEGER,    &
                          MPI_SUM, MPI_COMM_gfs, IERR )

      call MPI_BCAST( BUF, N, MPI_INTEGER, IROOT-1,           &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_broadcast_i_2d
!--------------------------------------------------------------
      subroutine mpe_broadcast_r4_scalar(buf,n,flag,type)
 
      use rank, only : myrank,MPI_COMM_gfs
      use mpi

      real*4   buf
      integer  n,iroot_in,iroot,IERR
      logical  flag
      integer, optional :: type
 
      iroot_in=0
      if (flag) iroot_in = myrank+1
 
      call MPI_ALLREDUCE( iroot_in, iroot, 1, MPI_INTEGER,    &
                          MPI_SUM, MPI_COMM_gfs, IERR )

      call MPI_BCAST( BUF, N, MPI_REAL4           , IROOT-1,  &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_broadcast_r4_scalar
!--------------------------------------------------------------
      subroutine mpe_broadcast_r8_scalar(buf,n,flag,type)
 
      use rank, only : myrank,MPI_COMM_gfs
      use mpi

      real*8   buf
      integer  n,iroot_in,iroot,IERR
      logical  flag
      integer, optional :: type
 
      iroot_in=0
      if (flag) iroot_in = myrank+1
 
      call MPI_ALLREDUCE( iroot_in, iroot, 1, MPI_INTEGER,    &
                          MPI_SUM, MPI_COMM_gfs, IERR )

      call MPI_BCAST( BUF, N, MPI_REAL8           , IROOT-1,  &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_broadcast_r8_scalar
!--------------------------------------------------------------
      subroutine mpe_broadcast_r4(buf,n,flag,type)
 
      use rank, only : myrank,MPI_COMM_gfs
      use mpi

      real*4   buf(n)
      integer  n,iroot_in,iroot,IERR
      logical  flag
      integer, optional :: type
 
      iroot_in=0
      if (flag) iroot_in = myrank+1
 
      call MPI_ALLREDUCE( iroot_in, iroot, 1, MPI_INTEGER,    &
                          MPI_SUM, MPI_COMM_gfs, IERR )

      call MPI_BCAST( BUF, N, MPI_REAL4           , IROOT-1,  &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_broadcast_r4
!--------------------------------------------------------------
      subroutine mpe_broadcast_r8(buf,n,flag,type)
 
      use rank, only : myrank,MPI_COMM_gfs
      use mpi

      real*8   buf(n)
      integer  n,iroot_in,iroot,IERR
      logical  flag
      integer, optional :: type
 
      iroot_in=0
      if (flag) iroot_in = myrank+1
 
      call MPI_ALLREDUCE( iroot_in, iroot, 1, MPI_INTEGER,    &
                          MPI_SUM, MPI_COMM_gfs, IERR )

      call MPI_BCAST( BUF, N, MPI_REAL8           , IROOT-1,  &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_broadcast_r8
!--------------------------------------------------------------
      subroutine mpe_broadcast_r4_2d(buf,n,flag,type)
 
      use rank, only : myrank,MPI_COMM_gfs
      use mpi

      real*4   buf(n,1)
      integer  n,iroot_in,iroot,IERR
      logical  flag
      integer, optional :: type
 
      iroot_in=0
      if (flag) iroot_in = myrank+1
 
      call MPI_ALLREDUCE( iroot_in, iroot, 1, MPI_INTEGER,    &
                          MPI_SUM, MPI_COMM_gfs, IERR )

      call MPI_BCAST( BUF, N, MPI_REAL4           , IROOT-1,  &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_broadcast_r4_2d
!--------------------------------------------------------------
      subroutine mpe_broadcast_r8_2d(buf,n,flag,type)
 
      use rank, only : myrank,MPI_COMM_gfs
      use mpi

      real*8   buf(n,1)
      integer  n,iroot_in,iroot,IERR
      logical  flag
      integer, optional :: type
 
      iroot_in=0
      if (flag) iroot_in = myrank+1
 
      call MPI_ALLREDUCE( iroot_in, iroot, 1, MPI_INTEGER,    &
                          MPI_SUM, MPI_COMM_gfs, IERR )

      call MPI_BCAST( BUF, N, MPI_REAL8           , IROOT-1,  &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_broadcast_r8_2d
!--------------------------------------------------------------
      subroutine mpe_bcast_i_scalar(buf,n,iroot,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      integer  buf,n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_INTEGER, IROOT,             &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_bcast_i_scalar
!--------------------------------------------------------------
      subroutine mpe_bcast_i(buf,n,iroot,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      integer  buf(n),n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_INTEGER, IROOT,             &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_bcast_i
!--------------------------------------------------------------
      subroutine mpe_bcast_i_2d(buf,n,iroot,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      integer  buf(n,1),n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_INTEGER, IROOT,             &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_bcast_i_2d
!--------------------------------------------------------------
      subroutine mpe_bcast_r4_scalar(buf,n,iroot,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      real*4   buf
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL4, IROOT,    &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_bcast_r4_scalar
!--------------------------------------------------------------
      subroutine mpe_bcast_r8_scalar(buf,n,iroot,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      real*8   buf
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL8, IROOT,    &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_bcast_r8_scalar
!--------------------------------------------------------------
      subroutine mpe_bcast_r4(buf,n,iroot,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      real*4   buf(n)
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL4, IROOT,    &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_bcast_r4
!--------------------------------------------------------------
      subroutine mpe_bcast_r8(buf,n,iroot,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      real*8   buf(n)
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL8, IROOT,    &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_bcast_r8
!--------------------------------------------------------------
      subroutine mpe_bcast_r4_2d(buf,n,iroot,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      real*4   buf(n,1)
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL4, IROOT,    &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_bcast_r4_2d
!--------------------------------------------------------------
      subroutine mpe_bcast_r8_2d(buf,n,iroot,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      real*8   buf(n,1)
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL8, IROOT,    &
                      MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_bcast_r8_2d
!--------------------------------------------------------------
      subroutine mpe_bcast_col_i_scalar(buf,n,iroot,type)
 
      use index, only : col_comm
      use mpi

      integer  buf,n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_INTEGER, IROOT,             &
                      col_comm, IERR )
 
      return

      end subroutine mpe_bcast_col_i_scalar
!--------------------------------------------------------------
      subroutine mpe_bcast_col_i(buf,n,iroot,type)
 
      use index, only : col_comm
      use mpi

      integer  buf(n),n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_INTEGER, IROOT,             &
                      col_comm, IERR )
 
      return

      end subroutine mpe_bcast_col_i
!--------------------------------------------------------------
      subroutine mpe_bcast_col_i_2d(buf,n,iroot,type)
 
      use index, only : col_comm
      use mpi

      integer  buf(n,1),n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_INTEGER, IROOT,             &
                      col_comm, IERR )
 
      return

      end subroutine mpe_bcast_col_i_2d
!--------------------------------------------------------------
      subroutine mpe_bcast_col_r4_scalar(buf,n,iroot,type)
 
      use index, only : col_comm
      use mpi

      real*4   buf
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL4, IROOT,    &
                      col_comm, IERR )
 
      return

      end subroutine mpe_bcast_col_r4_scalar
!--------------------------------------------------------------
      subroutine mpe_bcast_col_r8_scalar(buf,n,iroot,type)
 
      use index, only : col_comm
      use mpi

      real*8   buf
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL8, IROOT,    &
                      col_comm, IERR )
 
      return

      end subroutine mpe_bcast_col_r8_scalar
!--------------------------------------------------------------
      subroutine mpe_bcast_col_r4(buf,n,iroot,type)
 
      use index, only : col_comm
      use mpi

      real*4   buf(n)
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL4, IROOT,    &
                      col_comm, IERR )
 
      return

      end subroutine mpe_bcast_col_r4
!--------------------------------------------------------------
      subroutine mpe_bcast_col_r8(buf,n,iroot,type)
 
      use index, only : col_comm
      use mpi

      real*8   buf(n)
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL8, IROOT,    &
                      col_comm, IERR )
 
      return

      end subroutine mpe_bcast_col_r8
!--------------------------------------------------------------
      subroutine mpe_bcast_col_r4_2d(buf,n,iroot,type)
 
      use index, only : col_comm
      use mpi

      real*4   buf(n,1)
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL4, IROOT,    &
                      col_comm, IERR )
 
      return

      end subroutine mpe_bcast_col_r4_2d
!--------------------------------------------------------------
      subroutine mpe_bcast_col_r8_2d(buf,n,iroot,type)
 
      use index, only : col_comm
      use mpi

      real*8   buf(n,1)
      integer  n,iroot,IERR
      integer, optional :: type
 
      call MPI_BCAST( BUF, N, MPI_REAL8, IROOT,    &
                      col_comm, IERR )
 
      return

      end subroutine mpe_bcast_col_r8_2d
!--------------------------------------------------------------
      subroutine mpe_global_sum_i_scalar(swork,n,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      integer  swork,rwork,i,n,IERR
      integer, optional   :: type
 
      call MPI_ALLREDUCE( SWORK, RWORK, N, MPI_INTEGER,           &
                          MPI_SUM, MPI_COMM_gfs, IERR )
 
      swork=rwork
 
      return

      end subroutine mpe_global_sum_i_scalar
!--------------------------------------------------------------
      subroutine mpe_global_sum_i(swork,n,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      integer  swork(n),rwork(n),i,n,IERR
      integer, optional   :: type
 
      call MPI_ALLREDUCE( SWORK, RWORK, N, MPI_INTEGER,           &
                          MPI_SUM, MPI_COMM_gfs, IERR )
 
      do i=1,n
         swork(i)=rwork(i)
      enddo
 
      return

      end subroutine mpe_global_sum_i
!--------------------------------------------------------------
      subroutine mpe_global_sum_r4_scalar(swork,n,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      real*4 swork, rwork
      integer  i,n,IERR
      integer, optional   :: type
 
      call MPI_ALLREDUCE( SWORK, RWORK, N, MPI_REAL4,  &
                          MPI_SUM, MPI_COMM_gfs, IERR )
 
      swork=rwork
 
      return

      end subroutine mpe_global_sum_r4_scalar
!--------------------------------------------------------------
      subroutine mpe_global_sum_r8_scalar(swork,n,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi

      real*8 swork, rwork
      integer  i,n,IERR
      integer, optional   :: type
 
      call MPI_ALLREDUCE( SWORK, RWORK, N, MPI_REAL8,  &
                          MPI_SUM, MPI_COMM_gfs, IERR )
 
      swork=rwork
 
      return

      end subroutine mpe_global_sum_r8_scalar
!--------------------------------------------------------------
      subroutine mpe_global_sum_r4(swork,n,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      real*4 swork(n), rwork(n)
      integer  i,n,IERR
      integer, optional   :: type
 
      rwork=0.

      call MPI_ALLREDUCE( SWORK, RWORK, N, MPI_REAL4,  &
                          MPI_SUM, MPI_COMM_gfs, IERR )
 
      do i=1,n
         swork(i)=rwork(i)
      enddo
 
      return

      end subroutine mpe_global_sum_r4
!--------------------------------------------------------------
      subroutine mpe_global_sum_r8(swork,n,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      real*8 swork(n), rwork(n)
      integer  i,n,IERR
      integer, optional   :: type
 
      rwork=0.

      call MPI_ALLREDUCE( SWORK, RWORK, N, MPI_REAL8,  &
                          MPI_SUM, MPI_COMM_gfs, IERR )
 
      do i=1,n
         swork(i)=rwork(i)
      enddo
 
      return

      end subroutine mpe_global_sum_r8
!--------------------------------------------------------------
      subroutine mpe_global_max_i_scalar(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      integer  sbuf,rbuf,n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_INTEGER,      &
                       MPI_MAX, MPI_COMM_gfs, IERR )
 
      sbuf=rbuf
 
      return

      end subroutine mpe_global_max_i_scalar
!--------------------------------------------------------------
      subroutine mpe_global_max_i(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      integer  sbuf(n),rbuf(n),i,n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_INTEGER,      &
                       MPI_MAX, MPI_COMM_gfs, IERR )
 
      do i=1,n
        sbuf(i)=rbuf(i)
      enddo
 
      return

      end subroutine mpe_global_max_i
!--------------------------------------------------------------
      subroutine mpe_global_max_r4_scalar(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      real*4   sbuf, rbuf
      integer  n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_REAL4,  &
                       MPI_MAX, MPI_COMM_gfs, IERR )
 
      sbuf=rbuf
 
      return
      end subroutine mpe_global_max_r4_scalar
!--------------------------------------------------------------
      subroutine mpe_global_max_r8_scalar(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      real*8   sbuf, rbuf
      integer  n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_REAL8,  &
                       MPI_MAX, MPI_COMM_gfs, IERR )
 
      sbuf=rbuf
 
      return
      end subroutine mpe_global_max_r8_scalar
!--------------------------------------------------------------
      subroutine mpe_global_max_r4(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      real*4   sbuf(n), rbuf(n)
      integer  i,n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_REAL4,  &
                       MPI_MAX, MPI_COMM_gfs, IERR )
 
      do i=1,n
        sbuf(i)=rbuf(i)
      enddo
 
      return

      end subroutine mpe_global_max_r4
!--------------------------------------------------------------
      subroutine mpe_global_max_r8(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      real*8   sbuf(n), rbuf(n)
      integer  i,n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_REAL8,  &
                       MPI_MAX, MPI_COMM_gfs, IERR )
 
      do i=1,n
        sbuf(i)=rbuf(i)
      enddo
 
      return

      end subroutine mpe_global_max_r8
!--------------------------------------------------------------
      subroutine mpe_global_maxloc_i(a,n,i1,i2,i3,i4,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      integer swork(2), rwork(2)
      integer a,n,i1,i2,i3,i4,i,iroot,IERR
      integer ibuf(n)
      integer, optional  :: type

      swork(1)=a
      swork(2)=myrank
      if(n.ge.1) ibuf(1)=i1
      if(n.ge.2) ibuf(2)=i2
      if(n.ge.3) ibuf(3)=i3
      if(n.ge.4) ibuf(4)=i4
 
      call MPI_ALLREDUCE( SWORK, RWORK, 1, MPI_2INTEGER,    &
                       MPI_MAXLOC, MPI_COMM_gfs, IERR )
 
      iroot=rwork(2)
 
      call MPI_BCAST( IBUF, N, MPI_INTEGER, IROOT,          &
                      MPI_COMM_gfs, IERR )
 
      a =rwork(1)
      if(n.ge.1) i1=ibuf(1)
      if(n.ge.2) i2=ibuf(2)
      if(n.ge.3) i3=ibuf(3)
      if(n.ge.4) i4=ibuf(4)
 
      return 

      end subroutine mpe_global_maxloc_i
!--------------------------------------------------------------
      subroutine mpe_global_maxloc_r4(a,n,i1,i2,i3,i4,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      real*4 swork(2), rwork(2)
      real*4 a
      integer ibuf(n)
      integer n,i1,i2,i3,i4,i,iroot,IERR
      integer, optional  :: type

      swork(1)=a
      swork(2)=myrank
      if(n.ge.1) ibuf(1)=i1
      if(n.ge.2) ibuf(2)=i2
      if(n.ge.3) ibuf(3)=i3
      if(n.ge.4) ibuf(4)=i4
 
      call MPI_ALLREDUCE( SWORK, RWORK, 1, MPI_2REAL,  &
                       MPI_MAXLOC, MPI_COMM_gfs, IERR )
 
      iroot=rwork(2)
 
      call MPI_BCAST( IBUF, N, MPI_INTEGER, IROOT,                 &
                      MPI_COMM_gfs, IERR )
 
      a =rwork(1)
      if(n.ge.1) i1=ibuf(1)
      if(n.ge.2) i2=ibuf(2)
      if(n.ge.3) i3=ibuf(3)
      if(n.ge.4) i4=ibuf(4)
 
      return

      end subroutine mpe_global_maxloc_r4
!--------------------------------------------------------------
      subroutine mpe_global_maxloc_r8(a,n,i1,i2,i3,i4,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      real*8 swork(2), rwork(2)
      real*8 a
      integer ibuf(n)
      integer n,i1,i2,i3,i4,i,iroot,IERR
      integer, optional  :: type

      swork(1)=a
      swork(2)=myrank
      if(n.ge.1) ibuf(1)=i1
      if(n.ge.2) ibuf(2)=i2
      if(n.ge.3) ibuf(3)=i3
      if(n.ge.4) ibuf(4)=i4
 
      call MPI_ALLREDUCE( SWORK, RWORK, 1, MPI_2DOUBLE_PRECISION,  &
                       MPI_MAXLOC, MPI_COMM_gfs, IERR )
 
      iroot=rwork(2)
 
      call MPI_BCAST( IBUF, N, MPI_INTEGER, IROOT,                 &
                      MPI_COMM_gfs, IERR )
 
      a =rwork(1)
      if(n.ge.1) i1=ibuf(1)
      if(n.ge.2) i2=ibuf(2)
      if(n.ge.3) i3=ibuf(3)
      if(n.ge.4) i4=ibuf(4)
 
      return

      end subroutine mpe_global_maxloc_r8
!--------------------------------------------------------------
      subroutine mpe_global_min_i_scalar(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      integer  sbuf,rbuf,n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_INTEGER,      &
                       MPI_MIN, MPI_COMM_gfs, IERR )
 
      sbuf=rbuf
 
      return

      end subroutine mpe_global_min_i_scalar
!--------------------------------------------------------------
      subroutine mpe_global_min_i(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      integer  sbuf(n),rbuf(n),i,n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_INTEGER,      &
                       MPI_MIN, MPI_COMM_gfs, IERR )
 
      do i=1,n
        sbuf(i)=rbuf(i)
      enddo
 
      return

      end subroutine mpe_global_min_i
!--------------------------------------------------------------
      subroutine mpe_global_min_r4_scalar(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      real*4   sbuf, rbuf
      integer  n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_REAL4,  &
                       MPI_MIN, MPI_COMM_gfs, IERR )
 
      sbuf=rbuf
 
      return

      end subroutine mpe_global_min_r4_scalar
!--------------------------------------------------------------
      subroutine mpe_global_min_r8_scalar(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      real*8   sbuf, rbuf
      integer  n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_REAL8,  &
                       MPI_MIN, MPI_COMM_gfs, IERR )
 
      sbuf=rbuf
 
      return

      end subroutine mpe_global_min_r8_scalar
!--------------------------------------------------------------
      subroutine mpe_global_min_r4(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      real*4   sbuf(n), rbuf(n)
      integer  i,n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_REAL4,  &
                       MPI_MIN, MPI_COMM_gfs, IERR )
 
      do i=1,n
        sbuf(i)=rbuf(i)
      enddo
 
      return

      end subroutine mpe_global_min_r4
!--------------------------------------------------------------
      subroutine mpe_global_min_r8(sbuf,n,type)

      use rank, only : MPI_COMM_gfs
      use mpi
      real*8   sbuf(n), rbuf(n)
      integer  i,n,IERR
      integer, optional  :: type
 
      call MPI_ALLREDUCE( SBUF, RBUF, N, MPI_REAL8,  &
                       MPI_MIN, MPI_COMM_gfs, IERR )
 
      do i=1,n
        sbuf(i)=rbuf(i)
      enddo
 
      return

      end subroutine mpe_global_min_r8

      subroutine mpe_global_minloc_i(a,n,i1,i2,i3,i4,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      integer swork(2), rwork(2)
      integer a,n,i1,i2,i3,i4,i,iroot,IERR
      integer ibuf(n)
      integer, optional  :: type

      swork(1)=a
      swork(2)=myrank
      if(n.ge.1) ibuf(1)=i1
      if(n.ge.2) ibuf(2)=i2
      if(n.ge.3) ibuf(3)=i3
      if(n.ge.4) ibuf(4)=i4
 
      call MPI_ALLREDUCE( SWORK, RWORK, 1, MPI_2INTEGER,    &
                       MPI_MINLOC, MPI_COMM_gfs, IERR )
 
      iroot=rwork(2)
 
      call MPI_BCAST( IBUF, N, MPI_INTEGER, IROOT,          &
                      MPI_COMM_gfs, IERR )
 
      a =rwork(1)
      if(n.ge.1) i1=ibuf(1)
      if(n.ge.2) i2=ibuf(2)
      if(n.ge.3) i3=ibuf(3)
      if(n.ge.4) i4=ibuf(4)
 
      return 

      end subroutine mpe_global_minloc_i
!--------------------------------------------------------------
      subroutine mpe_global_minloc_r4(a,n,i1,i2,i3,i4,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      real*4 swork(2), rwork(2)
      real*4 a
      integer ibuf(n)
      integer n,i1,i2,i3,i4,i,iroot,IERR
      integer, optional  :: type

      swork(1)=a
      swork(2)=myrank
      if(n.ge.1) ibuf(1)=i1
      if(n.ge.2) ibuf(2)=i2
      if(n.ge.3) ibuf(3)=i3
      if(n.ge.4) ibuf(4)=i4
 
      call MPI_ALLREDUCE( SWORK, RWORK, 1, MPI_2REAL,  &
                       MPI_MINLOC, MPI_COMM_gfs, IERR )
 
      iroot=rwork(2)
 
      call MPI_BCAST( IBUF, N, MPI_INTEGER, IROOT,                 &
                      MPI_COMM_gfs, IERR )
 
      a =rwork(1)
      if(n.ge.1) i1=ibuf(1)
      if(n.ge.2) i2=ibuf(2)
      if(n.ge.3) i3=ibuf(3)
      if(n.ge.4) i4=ibuf(4)
 
      return

      end subroutine mpe_global_minloc_r4
!--------------------------------------------------------------
      subroutine mpe_global_minloc_r8(a,n,i1,i2,i3,i4,type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      real*8 swork(2), rwork(2)
      real*8 a
      integer ibuf(n)
      integer n,i1,i2,i3,i4,i,iroot,IERR
      integer, optional  :: type

      swork(1)=a
      swork(2)=myrank
      if(n.ge.1) ibuf(1)=i1
      if(n.ge.2) ibuf(2)=i2
      if(n.ge.3) ibuf(3)=i3
      if(n.ge.4) ibuf(4)=i4
 
      call MPI_ALLREDUCE( SWORK, RWORK, 1, MPI_2DOUBLE_PRECISION,  &
                       MPI_MINLOC, MPI_COMM_gfs, IERR )
 
      iroot=rwork(2)
 
      call MPI_BCAST( IBUF, N, MPI_INTEGER, IROOT,                 &
                      MPI_COMM_gfs, IERR )
 
      a =rwork(1)
      if(n.ge.1) i1=ibuf(1)
      if(n.ge.2) i2=ibuf(2)
      if(n.ge.3) i3=ibuf(3)
      if(n.ge.4) i4=ibuf(4)
 
      return

      end subroutine mpe_global_minloc_r8
!--------------------------------------------------------------
      subroutine mpe_send_print_i(sbuf, n, j, type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      integer  sbuf(n)
      integer  n,j,idst,IERR
      integer, optional   :: type
 
      idst=0
      call MPI_SEND( SBUF, N, MPI_INTEGER, IDST, J, MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_send_print_i
!--------------------------------------------------------------
      subroutine mpe_send_print_r4(sbuf, n, j, type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      real*4   sbuf(n,1)
      integer  n,j,idst,IERR
      integer, optional   :: type
 
      idst=0
      call MPI_SEND( SBUF, N, MPI_REAL4, IDST, J, MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_send_print_r4
!--------------------------------------------------------------
      subroutine mpe_send_print_r8(sbuf, n, j, type)
 
      use rank, only : MPI_COMM_gfs
      use mpi
      real*8   sbuf(n,1)
      integer  n,j,idst,IERR
      integer, optional   :: type
 
      idst=0
      call MPI_SEND( SBUF, N, MPI_REAL8, IDST, J, MPI_COMM_gfs, IERR )
 
      return

      end subroutine mpe_send_print_r8
!--------------------------------------------------------------
      subroutine mpe_recv_print_i(rbuf, n, j, type)
      use rank, only : MPI_COMM_gfs
      include 'mpif.h'
!     use mpi
 
      integer  rbuf(n)
      integer  n,j,ISTATUS,IERR
      integer, optional   :: type
 
      call MPI_RECV( RBUF, N, MPI_INTEGER, MPI_ANY_SOURCE,   &
                     J, MPI_COMM_gfs, ISTATUS, IERR )
 
      return
      end subroutine mpe_recv_print_i
!--------------------------------------------------------------
      subroutine mpe_recv_print_r4(rbuf, n, j, type)
      use rank, only : MPI_COMM_gfs
      include 'mpif.h'
!     use mpi
 
      real*4   rbuf(n,1)
      integer  n,j,ISTATUS,IERR
      integer, optional   :: type
 
      call MPI_RECV( RBUF, N, MPI_REAL4, MPI_ANY_SOURCE,  &
                     J, MPI_COMM_gfs, ISTATUS, IERR )
 
      return
      end subroutine mpe_recv_print_r4
!--------------------------------------------------------------
      subroutine mpe_recv_print_r8(rbuf, n, j, type)
      use rank, only : MPI_COMM_gfs
      include 'mpif.h'
!     use mpi
 
      real*8   rbuf(n,1)
      integer  n,j,ISTATUS,IERR
      integer, optional   :: type
 
      call MPI_RECV( RBUF, N, MPI_REAL8, MPI_ANY_SOURCE,  &
                     J, MPI_COMM_gfs, ISTATUS, IERR )
 
      return
      end subroutine mpe_recv_print_r8
!--------------------------------------------------------------

      end module mpe
