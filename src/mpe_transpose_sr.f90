      subroutine mpe_transpose_sr(sbuf,rbuf,lev,n,m,nsize,comm)
!
      use mpi

      implicit none
      integer  n,m,lev,nsize,i,j,k,ii,len_tr,ierr,comm

      real sbuf (lev,n,m*nsize)
      real rbuf (lev,n,nsize,m)
#ifdef MPISP
#   ifdef UNIFY_ONLY
      real*8 rwork(lev,n,m,nsize)
#   else
      real*4 rwork(lev,n,m,nsize),sbuf1(lev,n,m*nsize)
#   endif
#else
      real*8 rwork(lev,n,m,nsize)
#endif
!
!
      len_tr=m*n

#ifdef MPISP
#   ifdef UNIFY_ONLY
         call MPI_ALLTOALL( SBUF,  LEN_TR*LEV, MPI_REAL8,   &
                            RWORK, LEN_TR*LEV, MPI_REAL8,   &
                            comm,              IERR )
#   else
         sbuf1=sbuf
         call MPI_ALLTOALL( sbuf1, LEN_TR*LEV, MPI_REAL4,   &
                            RWORK, LEN_TR*LEV, MPI_REAL4,   &
                            comm,              IERR )
#   endif
#else
         call MPI_ALLTOALL( SBUF,  LEN_TR*LEV, MPI_REAL8,   &
                            RWORK, LEN_TR*LEV, MPI_REAL8,   &
                            comm,              IERR )
#endif
!
      do j=1,m
      do ii=1,nsize
      do i=1,n
      do k=1,lev
        rbuf(k,i,ii,j)=rwork(k,i,j,ii)
      enddo
      enddo
      enddo
      enddo
!
      return
      end
!------------------------------------------------------------
!CWB2021 for single precision test
      subroutine mpe_transpose_sr_sp(sbuf,rbuf,lev,n,m,nsize,comm)
!
      use const, only : RTYPE,MPI_RTYPE
      use mpi

      implicit none
      integer  n,m,lev,nsize,i,j,k,ii,len_tr,ierr,comm

      real(kind=RTYPE) sbuf(lev,n,m*nsize),rbuf(lev,n,nsize,m)
      real(kind=RTYPE) rwork(lev,n,m,nsize)
 
      len_tr=m*n

      call MPI_ALLTOALL( SBUF,  LEN_TR*LEV, MPI_RTYPE,   &
                         RWORK, LEN_TR*LEV, MPI_RTYPE,   &
                         comm,              IERR )
 
      do j=1,m
      do ii=1,nsize
      do i=1,n
      do k=1,lev
        rbuf(k,i,ii,j)=rwork(k,i,j,ii)
      enddo
      enddo
      enddo
      enddo
!
      return
      end
