      subroutine mpe_send_key(key, j, ist)
!
      use rank, only : root_io
      use mpi
      character*34 key

      call MPI_SEND( key, 34, MPI_CHARACTER, root_io, J, &
                     MPI_COMM_WORLD, ist )
      return
      end
