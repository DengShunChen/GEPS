      subroutine mpe_recv_data(rbuf, n, tag, ierr)
!
#if defined(RSM) && defined(CWB_MPMD)
      use rank, only : root_gfs,MPI_COMM_gfs_all
#else
      use rank, only : root_gfs,MPI_COMM_atm
#endif
      use mpi
      integer n,tag,isrc,ierr,ISTATUS(MPI_STATUS_SIZE)
      real*8  rbuf(n)

#if defined(RSM) && defined(CWB_MPMD)
      call MPI_RECV( RBUF, n, MPI_DOUBLE_PRECISION, root_gfs, &
                     tag, MPI_COMM_gfs_all, ISTATUS,  IERR )
#else
      call MPI_RECV( RBUF, n, MPI_DOUBLE_PRECISION, root_gfs, &
                     tag, MPI_COMM_atm, ISTATUS,  IERR )
#endif
      return
      end
