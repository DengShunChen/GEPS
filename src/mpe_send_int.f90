      subroutine mpe_send_int( itau ,n , itag , ist)
#if defined(RSM) && defined(CWB_MPMD)
      use rank, only : root_io,MPI_COMM_gfs_all
#else
      use rank, only : root_io
#endif
      use mpi
      integer::n
      integer :: itau(n)  ,itag

#if defined(RSM) && defined(CWB_MPMD)
      call MPI_SEND( itau, n, MPI_integer , root_io, itag , &
                     MPI_COMM_gfs_all, ist )
#else
      call MPI_SEND( itau, n, MPI_integer , root_io, itag , &
                     MPI_COMM_WORLD, ist )
#endif
      return
      end
!=================================
      subroutine mpe_recv_int( itau ,n , itag , ist)
#if defined(RSM) && defined(CWB_MPMD)
      use rank, only : root_gfs,MPI_COMM_gfs_all
#else
      use rank, only : root_gfs
#endif
      use mpi
      integer::n
      integer :: itau(n)  ,itag
      integer:: ierr ,ISTATUS(MPI_STATUS_SIZE)

#if defined(RSM) && defined(CWB_MPMD)
      call MPI_RECV( itau , n, MPI_INTEGER , root_gfs, &
                     itag, MPI_COMM_gfs_all, ISTATUS, ierr )
#else
      call MPI_RECV( itau , n, MPI_INTEGER , root_gfs, &
                     itag, MPI_COMM_WORLD, ISTATUS, ierr )
#endif
      return
      end
