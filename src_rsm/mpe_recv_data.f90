      subroutine mpe_recv_data(rbuf, n, tag, ierr)
!
#if defined(RSM) && defined(CWB_MPMD)
      use rank, only : root_gfs,MPI_COMM_gfs_all
#else
      use rank, only : root_gfs,MPI_COMM_atm
#endif
      use mpi
      integer n,tag,isrc,ierr,ISTATUS(MPI_STATUS_SIZE)
#ifdef SP
      real*4  rbuf(n)
#else
      real*8  rbuf(n)
#endif

#if defined(RSM) && defined(CWB_MPMD)
#ifdef SP
      call MPI_RECV( RBUF, n, MPI_REAL, root_gfs, &
                     tag, MPI_COMM_gfs_all, ISTATUS,  IERR )
#else
      call MPI_RECV( RBUF, n, MPI_DOUBLE_PRECISION, root_gfs, &
                     tag, MPI_COMM_gfs_all, ISTATUS,  IERR )
#endif
#else
#ifdef SP
      call MPI_RECV( RBUF, n, MPI_REAL, root_gfs, &
                     tag, MPI_COMM_WORLD, ISTATUS,  IERR )
#else
      call MPI_RECV( RBUF, n, MPI_DOUBLE_PRECISION, root_gfs, &
                     tag, MPI_COMM_atm, ISTATUS,  IERR )
#endif
#endif
      return
      end
