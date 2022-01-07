      subroutine mpe_recv_key(key, tag, ierr)
!
#if defined(RSM) && defined(CWB_MPMD)
      use rank, only : root_gfs,MPI_COMM_gfs_all
#else
      use rank, only : root_gfs
#endif
      use mpi
      character*38 key
      integer tag,isrc,ierr,ISTATUS(MPI_STATUS_SIZE)

#if defined(RSM) && defined(CWB_MPMD)
      call MPI_RECV( key, 38, MPI_CHARACTER, root_gfs, &
                     tag, MPI_COMM_gfs_all, ISTATUS, ierr )
#else
      call MPI_RECV( key, 38, MPI_CHARACTER, root_gfs, &
                     tag, MPI_COMM_WORLD, ISTATUS, ierr )
#endif
      return
      end
