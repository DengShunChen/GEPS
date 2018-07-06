      subroutine mpe_transpose_rs(sbuf,rbuf,lev,n,m,nsize,comm)
!
      use mpi

      implicit none
      integer  n,m,lev,nsize,i,j,k,ii,len_tr,ierr,comm

      real sbuf (lev,n,nsize,m)
      real rbuf (lev,n,m*nsize)
#ifdef MPISP
#   ifdef UNIFY_ONLY
      real*8 swork(lev,n,m,nsize)
#   else
      real*4 swork(lev,n,m,nsize), rwork (lev,n,m*nsize)
#   endif
#else
      real*8 swork(lev,n,m,nsize)
#endif
!
      len_tr=m*n
!
      do j=1,m
      do ii=1,nsize
      do i=1,n
      do k=1,lev
        swork(k,i,j,ii)=sbuf(k,i,ii,j)
      enddo
      enddo
      enddo
      enddo
!
#ifdef MPISP
#   ifdef UNIFY_ONLY
         call MPI_ALLTOALL( SWORK, LEN_TR*LEV, MPI_REAL8,   &
                            RBUF , LEN_TR*LEV, MPI_REAL8,   &
                            comm,              IERR )
#   else
         call MPI_ALLTOALL( SWORK, LEN_TR*LEV, MPI_REAL4,   &
                            rwork, LEN_TR*LEV, MPI_REAL4,   &
                            comm,              IERR )
         rbuf=rwork
#   endif
#else
         call MPI_ALLTOALL( SWORK, LEN_TR*LEV, MPI_REAL8,   &
                            RBUF , LEN_TR*LEV, MPI_REAL8,   &
                            comm,              IERR )
#endif
!
      return
      end
