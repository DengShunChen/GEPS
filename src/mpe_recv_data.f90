      subroutine mpe_recv_data(rbuf, n, tag, ierr)
!
      use rank, only : root_gfs
      use mpi
      integer n,tag,isrc,ierr,ISTATUS(MPI_STATUS_SIZE)
      real*8  rbuf(n)

      call MPI_RECV( RBUF, n, MPI_DOUBLE_PRECISION, root_gfs, &
                     tag, MPI_COMM_WORLD, ISTATUS,  IERR )
      return
      end
