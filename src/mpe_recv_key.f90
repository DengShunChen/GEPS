      subroutine mpe_recv_key(key, tag, ierr)
!
#if defined(RSM) && defined(CWB_MPMD)
      use rank, only : root_gfs,MPI_COMM_gfs_all
#else
      use rank, only : root_gfs,MPI_COMM_atm
#endif
      use const, only : keyo,KLENO
      use mpi
      integer tag,isrc,ierr,ISTATUS(MPI_STATUS_SIZE)

#if defined(RSM) && defined(CWB_MPMD)
      call MPI_RECV( keyo, KLENO, MPI_CHARACTER, root_gfs, &
                     tag, MPI_COMM_gfs_all, ISTATUS, ierr )
#else
      call MPI_RECV( keyo, KLENO, MPI_CHARACTER, root_gfs, &
                     tag, MPI_COMM_WORLD, ISTATUS, ierr )
#endif
      return
      end
