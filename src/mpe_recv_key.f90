      subroutine mpe_recv_key(key, tag, ierr)
!
      use rank, only : root_gfs
      use mpi
      character*34 key
      integer tag,isrc,ierr,ISTATUS(MPI_STATUS_SIZE)

      call MPI_RECV( key, 34, MPI_CHARACTER, root_gfs, &
                     tag, MPI_COMM_WORLD, ISTATUS, ierr )
      return
      end
