      subroutine mpe_scatter_io(a,b,len,nsize)
!     include 'mpif.h'
      use rank, only :MPI_COMM_gfs
      use mpi
      use const, only: RTYPE
!
      implicit  none
      real(kind=RTYPE) a(len*nsize)
      real(kind=RTYPE) b(len)
      integer   len,nsize,iroot,ierr
!
      iroot=0
#ifdef SP
      call MPI_SCATTER(A,LEN,       MPI_REAL4,                &
                       B,LEN,       MPI_REAL4,                &
                       IROOT,       MPI_COMM_gfs,IERR )
#else      
      call MPI_SCATTER(A,LEN,       MPI_DOUBLE_PRECISION,     &
                       B,LEN,       MPI_DOUBLE_PRECISION,     &
                       IROOT,       MPI_COMM_gfs,IERR )
#endif
!
      return
      end
!
      subroutine mpe_scatter_sppt(a,b,len,nsize)
!     include 'mpif.h'
      use index, only :col_comm
      use mpi
      use const, only :RTYPE,MPI_RTYPE
!
      implicit  none
      real      a(len*nsize)
      real(kind=RTYPE) b(len)
      integer   len,nsize,iroot,ierr
!
      iroot=0
      call MPI_SCATTER(A,LEN,       MPI_DOUBLE_PRECISION,     &
                       B,LEN,       MPI_RTYPE,                &
                       IROOT,       col_comm,IERR )
!
      return
      end
