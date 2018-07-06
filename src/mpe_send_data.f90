      subroutine mpe_send_data(sbuf, n, j, ist)
!
      use rank, only : root_io
      use mpi
      real*8  sbuf(n)

      call MPI_SEND( SBUF, n, MPI_DOUBLE_PRECISION, root_io, J, &
                     MPI_COMM_WORLD, ist )
 
      return
      end
