      subroutine mpe_barrier
!
!     include 'mpif.h'
      use rank
      use mpi

      implicit none

      integer  ierr
!
      call MPI_BARRIER(MPI_COMM_gfs, IERR )
!
      return
      end
