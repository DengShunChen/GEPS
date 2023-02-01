      subroutine mpe_unify(a,n,m,idcmp,type)

!2dMPI version

      use param
      use mpe
      use rank
      use index
      use mpi

      integer n,m,idcmp,type
      dimension a(*)
!
      if(idcmp .eq. 1) then
          call mpe_unify1_r(a,n,m,my_max,nsize)
      else if(idcmp .eq. 2) then
        if(type .eq. mpe_integer) then
          call mpe_unify2_i(a,n,m,my_max,nsize)
!CWB2021
        else if(type .eq. mpe_single) then
          call mpe_unify2_r_sp(a,n,m,my_max,nsize)
        else if(type .eq. mpe_double) then
          call mpe_unify2_r(a,n,m,my_max,nsize)
        else if(type .eq. mpe_logical) then
          call mpe_unify2_l(a,n,m,my_max,nsize)
        else
          write(6,*) 'mpe_unify: Argument(type) Error  RANK=',myrank
        endif
      else if(idcmp .eq. 3) then
        call mpe_unify3_r(a,n,m,jtmax,nsizey)
      else if(idcmp .eq. 4) then
        call mpe_unify4_r(a,n,m,my_max,nsize)
      else if(idcmp .eq. 5) then
        if(type .eq. mpe_integer) then
          call mpe_unify5_i(a,n,m,my_max,nsize)
        else if(type .eq. mpe_double) then
          call mpe_unify5_r(a,n,m,my_max,nsize)
        else if(type .eq. mpe_logical) then
          call mpe_unify5_l(a,n,m,my_max,nsize)
        else
          write(6,*) 'mpe_unify: Argument(type) Error  RANK=',myrank
        endif
      else
        write(6,*) 'mpe_unify: Argument(idcmp) Error  RANK=',myrank
      endif
!
      return
      end
!-------------------------------------------------------------------------
      subroutine mpe_unify1_r(a,m,n,mx,nsize)

      use rank, only : MPI_COMM_gfs
      use index
      use mpi

      real a(m,n)
#ifdef MPISP
      real*4 b1(n,mx),b2(n,mx*nsize)
#else
      real*8 b1(n,mx),b2(n,mx*nsize)
#endif

      do jj=1,jlistnum
        j1=jlist1(jj)
      do i=1,n
        b1(i,jj)=a(j1,i)
      enddo
      enddo

#ifdef MPISP
      call MPI_ALLGATHER( B1,n*mx,   MPI_REAL4,     &
                          B2,n*mx,   MPI_REAL4,     &
                          MPI_COMM_gfs,  IERR )
#else
      call MPI_ALLGATHER( B1,n*mx,   MPI_REAL8,     &
                          B2,n*mx,   MPI_REAL8,     &
                          MPI_COMM_gfs,  IERR )
#endif
!
      a=0.
      do j=1,m
      do ii=1,nsizex
         jf=jlist2_2d(ii,j)
      do i=1,n
        a(j,i)=a(j,i)+b2(i,jf)
      enddo
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe_unify2_i(a,n,m,mx,nsize)
 
      use rank, only : MPI_COMM_gfs
      use index
      use mpi

      integer a(n,m)
      integer b1(nxp,mx)
      integer b2(nxp,mx*nsize)

      do jj=1,jlistnum
         j=jlist1(jj)
      do i=1,nxjlen(j)
         b1(i,jj)=a(i,j)
      enddo
      enddo

      call MPI_ALLGATHER( B1,nxp*mx, MPI_INTEGER, &
                          B2,nxp*mx, MPI_INTEGER, &
                          MPI_COMM_gfs,   IERR )

      do j=1,m
         ii=1
      do i=1,nsizex
         jf=jlist2_2d(i,j)
         nn=nxjlen_all(i,j)
         a(ii:ii+nn-1,j)=b2(1:nn,jf)
         ii=ii+nn
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe_unify2_r(a,n,m,mx,nsize)
 
      use rank, only : MPI_COMM_gfs
      use index
      use mpi

      real a(n,m)
#ifdef MPISP
      real*4 b1(nxp,mx),b2(nxp,mx*nsize)
#else
      real*8 b1(nxp,mx),b2(nxp,mx*nsize)
#endif

      do jj=1,jlistnum
         j=jlist1(jj)
      do i=1,nxjlen(j)
         b1(i,jj)=a(i,j)
      enddo
      enddo

#ifdef MPISP
      call MPI_ALLGATHER( B1,nxp*mx,   MPI_REAL4,    &
                          B2,nxp*mx,   MPI_REAL4,    &
                          MPI_COMM_gfs,  IERR )
#else
      call MPI_ALLGATHER( B1,nxp*mx,   MPI_REAL8,    &
                          B2,nxp*mx,   MPI_REAL8,    &
                          MPI_COMM_gfs,  IERR )
#endif
!
      do j=1,m
         ii=1
      do i=1,nsizex
         jf=jlist2_2d(i,j)
         nn=nxjlen_all(i,j)
         a(ii:ii+nn-1,j)=b2(1:nn,jf)
         ii=ii+nn
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
!CWB2021 for single precision test
      subroutine mpe_unify2_r_sp(a,n,m,mx,nsize)
 
      use rank, only : MPI_COMM_gfs
      use index
      use mpi
      use const, only : RTYPE,MPI_RTYPE

      real(kind=RTYPE) a(n,m),b1(nxp,mx),b2(nxp,mx*nsize)

      do jj=1,jlistnum
         j=jlist1(jj)
      do i=1,nxjlen(j)
         b1(i,jj)=a(i,j)
      enddo
      enddo

      call MPI_ALLGATHER( B1,nxp*mx,   MPI_RTYPE,    &
                          B2,nxp*mx,   MPI_RTYPE,    &
                          MPI_COMM_gfs,  IERR )
 
      do j=1,m
         ii=1
      do i=1,nsizex
         jf=jlist2_2d(i,j)
         nn=nxjlen_all(i,j)
         a(ii:ii+nn-1,j)=b2(1:nn,jf)
         ii=ii+nn
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      subroutine mpe_unify2_l(a,n,m,mx,nsize)
 
      use rank, only : MPI_COMM_gfs
      use index
      use mpi

      logical a(n,m)
      logical b1(nxp,mx)
      logical b2(nxp,mx*nsize)

      do jj=1,jlistnum
         j=jlist1(jj)
      do i=1,nxjlen(j)
         b1(i,jj)=a(i,j)
      enddo
      enddo

      call MPI_ALLGATHER( B1,nxp*mx, MPI_LOGICAL, &
                          B2,nxp*mx, MPI_LOGICAL, &
                          MPI_COMM_gfs,   IERR )

      do j=1,m
         ii=1
      do i=1,nsizex
         jf=jlist2_2d(i,j)
         nn=nxjlen_all(i,j)
         a(ii:ii+nn-1,j)=b2(1:nn,jf)
         ii=ii+nn
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe_unify3_r(a,n,m,jx,nsize)

      use rank, only : MPI_COMM_gfs
      use index
      use mpi
 
      real a(n,m)
      real b1(n,jx)
      real b2(n,jx*nsize)

      do jj=1,mlistnum
        j1=mlist(jj)
      do i=1,n
        b1(i,jj)=a(i,j1)
      enddo
      enddo
 
      call MPI_ALLGATHER( B1,N*jx,   MPI_DOUBLE_PRECISION, &
                          B2,N*jx,   MPI_DOUBLE_PRECISION, &
                          col_comm, IERR )
 
      do j=1,m
        jf=nlist(j)
      do i=1,n
        a(i,j)=b2(i,jf)
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe_unify4_r(a,n,m,mx,nsize)
 
      use rank, only : MPI_COMM_gfs
      use index
      use mpi

      real a(n,m)
      real b1(n,mx)
      real b2(n,mx*nsize)

      do jj=1,jlistnum
         j=jlist1(jj)
         b1(1:n,jj)=a(1:n,j)
      enddo

      call MPI_ALLGATHER( B1,n*mx, MPI_DOUBLE_PRECISION, &
                          B2,n*mx, MPI_DOUBLE_PRECISION, &
                          MPI_COMM_gfs,   IERR )

      a=0.
      do j=1,m
      do i=1,nsizex
         jf=jlist2_2d(i,j)
         a(:,j)=a(:,j)+b2(:,jf)
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe_unify5_i(a,n,m,mx,nsize)
 
      use rank, only : MPI_COMM_gfs
      use index
      use mpi

      integer a(n,m)
      integer b1(n,mx)
      integer b2(n,mx*nsize)

      do jj=1,jlistnum
         j=jlist1(jj)
         b1(:,jj)=a(:,j)
      enddo

      call MPI_ALLGATHER( B1,n*mx, MPI_INTEGER, &
                          B2,n*mx, MPI_INTEGER, &
                          MPI_COMM_gfs,   IERR )

      do j=1,m
      do i=1,nsizex
         jf=jlist2_2d(i,j)
         a(:,j)=b2(:,jf)
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe_unify5_r(a,n,m,mx,nsize)
 
      use rank, only : MPI_COMM_gfs
      use index
      use mpi

      real a(n,m)
#ifdef MPISP
      real*4 b1(n,mx),b2(n,mx*nsize)
#else
      real*8 b1(n,mx),b2(n,mx*nsize)
#endif

      do jj=1,jlistnum
         j=jlist1(jj)
         b1(:,jj)=a(:,j)
      enddo

#ifdef MPISP
      call MPI_ALLGATHER( B1,n*mx,   MPI_REAL4,    &
                          B2,n*mx,   MPI_REAL4,    &
                          MPI_COMM_gfs,  IERR )
#else
      call MPI_ALLGATHER( B1,n*mx,   MPI_REAL8,    &
                          B2,n*mx,   MPI_REAL8,    &
                          MPI_COMM_gfs,  IERR )
#endif
!
      do j=1,m
      do i=1,nsizex
         jf=jlist2_2d(i,j)
         a(:,j)=b2(:,jf)
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe_unify5_l(a,n,m,mx,nsize)
 
      use rank, only : MPI_COMM_gfs
      use index
      use mpi

      logical a(n,m)
      logical b1(n,mx)
      logical b2(n,mx*nsize)

      do jj=1,jlistnum
         j=jlist1(jj)
         b1(:,jj)=a(:,j)
      enddo

      call MPI_ALLGATHER( B1,n*mx, MPI_LOGICAL, &
                          B2,n*mx, MPI_LOGICAL, &
                          MPI_COMM_gfs,   IERR )

      do j=1,m
      do i=1,nsizex
         jf=jlist2_2d(i,j)
         a(:,j)=b2(:,jf)
      enddo
      enddo

      return
      end
