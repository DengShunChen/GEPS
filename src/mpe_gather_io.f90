      subroutine mpe_gather_io(a,b,len,nsize)
!     include 'mpif.h'
      use rank, only : MPI_COMM_gfs
      use mpi
!
      implicit  none
      real      a(len*nsize)
      real      b(len)
      integer   len,nsize,iroot,ierr
!
      iroot=0
!
      call MPI_BARRIER(MPI_COMM_gfs, IERR)
!
      call MPI_GATHER( B,LEN,       MPI_DOUBLE_PRECISION,     &
                       A,LEN,       MPI_DOUBLE_PRECISION,     &
                       IROOT,       MPI_COMM_gfs, IERR )
!
      return
      end
