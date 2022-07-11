      subroutine mpe_gather_io(a,b,len,nsize)
!     include 'mpif.h'
      use rank, only : MPI_COMM_gfs
      use mpi
      use const, only: RTYPE
!
      implicit  none
      real(kind=RTYPE) a(len*nsize)
      real(kind=RTYPE) b(len)
      integer   len,nsize,iroot,ierr
!
      iroot=0
!
      call MPI_BARRIER(MPI_COMM_gfs, IERR)
!
#ifdef SP
      call MPI_GATHER( B,LEN,       MPI_REAL4,                &
                       A,LEN,       MPI_REAL4,                &
                       IROOT,       MPI_COMM_gfs, IERR )
#else               
      call MPI_GATHER( B,LEN,       MPI_DOUBLE_PRECISION,     &
                       A,LEN,       MPI_DOUBLE_PRECISION,     &
                       IROOT,       MPI_COMM_gfs, IERR )
#endif
!
      return
      end
