      subroutine mpe_scatter_io(a,b,len,nsize)
!     include 'mpif.h'
      use rank, only :MPI_COMM_gfs
      use mpi
!
      implicit  none
      real      a(len*nsize)
      real      b(len)
      integer   len,nsize,iroot,ierr
!
      iroot=0
      call MPI_SCATTER(A,LEN,       MPI_DOUBLE_PRECISION,     &
                       B,LEN,       MPI_DOUBLE_PRECISION,     &
                       IROOT,       MPI_COMM_gfs,IERR )
!
      return
      end
