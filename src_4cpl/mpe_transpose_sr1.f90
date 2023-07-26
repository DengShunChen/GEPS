      subroutine mpe_transpose_sr1(sbuf,rbuf,n,m,lev,nsize,comm)
!
      use mpi

      implicit none
      integer  n,m,lev,nsize,i,j,k,ii,len_tr,ierr,comm

      real sbuf(n,m*nsize,lev),rbuf(n,nsize,m,lev)
#ifdef MPISP
#   ifdef UNIFY_ONLY
      real*8 swork(lev,n,m*nsize),rwork(lev,n,m,nsize)
#   else
      real*4 swork(lev,n,m*nsize),rwork(lev,n,m,nsize)
#   endif
#else
      real*8 swork(lev,n,m*nsize),rwork(lev,n,m,nsize)
#endif

      len_tr=m*n
!
      do j=1,m*nsize
      do i=1,n
      do k=1,lev
       swork(k,i,j)=sbuf(i,j,k)
      enddo
      enddo
      enddo
!
#ifdef MPISP
#   ifdef UNIFY_ONLY
      call MPI_ALLTOALL( SWORK, LEN_TR*LEV, MPI_REAL8,   &
                         RWORK, LEN_TR*LEV, MPI_REAL8,   &
                         comm,           IERR )
#   else
      call MPI_ALLTOALL( SWORK, LEN_TR*LEV, MPI_REAL4,   &
                         RWORK, LEN_TR*LEV, MPI_REAL4,   &
                         comm,           IERR )
#   endif
#else
      call MPI_ALLTOALL( SWORK, LEN_TR*LEV, MPI_REAL8,   &
                         RWORK, LEN_TR*LEV, MPI_REAL8,   &
                         comm,           IERR )
#endif
!
      do k=1,lev
      do j=1,m
      do ii=1,nsize
      do i=1,n
        rbuf(i,ii,j,k)=rwork(k,i,j,ii)
      enddo
      enddo
      enddo
      enddo
!
      return
      end
