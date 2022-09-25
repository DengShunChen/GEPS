      subroutine mpe_gather_io(a,b,len,nsize)
!     include 'mpif.h'
      use rank, only : MPI_COMM_gfs
      use mpi
      use const, only: RTYPE,MPI_RTYPE
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
      call MPI_GATHER( B,LEN,       MPI_RTYPE,                &
                       A,LEN,       MPI_RTYPE,                &
                       IROOT,       MPI_COMM_gfs, IERR )
!
      return
      end
