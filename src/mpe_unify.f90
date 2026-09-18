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

      use rank, only : MPI_COMM_gfs, myrank
      use index
      use mpi

      real a(n,m)
      real, allocatable :: b1(:), b2(:)
      integer nn,mm,jxx,ncount,i,j,jj,j1,jf,ierr
      integer(kind=8) :: dbg0, dbg1, dbg2
      interface
        subroutine geps_dbg_log(hyp, locid, irank, p0, p1, p2)
          integer hyp, locid, irank
          integer(kind=8) p0, p1, p2
        end subroutine
        subroutine geps_mpi_allgather_dbl(send, recv, count, fcomm, ierr)
          real send, recv
          integer count, fcomm, ierr
        end subroutine
      end interface

      nn = n
      mm = m
      jxx = jx
      ncount = nn*jxx
      allocate(b1(ncount), b2(ncount*nsize))
      b1 = 0.

      do jj=1,mlistnum
        j1=mlist(jj)
      do i=1,nn
        b1(i+(jj-1)*nn)=a(i,j1)
      enddo
      enddo

      ! #region agent log
      if (nn .eq. 1) then
        dbg0 = loc(b1(1))
        dbg1 = loc(b2(1))
        dbg2 = ncount
        call geps_dbg_log(18, 261, myrank, dbg0, dbg1, dbg2)
      endif
      ! #endregion
      call geps_mpi_allgather_dbl(b1(1), b2(1), ncount, col_comm, ierr)
      ! #region agent log
      if (nn .eq. 1) then
        dbg0 = ierr
        dbg1 = loc(b2(1))
        dbg2 = loc(a)
        call geps_dbg_log(18, 262, myrank, dbg0, dbg1, dbg2)
      endif
      ! #endregion

      do j=1,mm
        jf=nlist(j)
      do i=1,nn
        a(i,j)=b2(i+(jf-1)*nn)
      enddo
      enddo
      deallocate(b1, b2)

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
