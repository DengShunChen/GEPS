!----------------------------------------------------------------------------------------
!new routines created for 2Dmpi
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_nxp_lev(ain,aout,nxp,nx,lev,levp,num,my,my_max,jlistnum,jlen,nsizex,comm)

! transpose (nx partial,lev full) to (nx full,lev partial), num variable packed

      use index, only : jlist1,nxjlen_all,nxjlen
      use const, only : RTYPE,MPI_RTYPE

      implicit none

      include 'mpif.h'
      integer  nx,nxp,lev,levp,my,my_max,jlen,nsizex,comm
      real(kind=RTYPE) ain(nxp,lev,num,my_max),aout(nx,levp,num,my_max)
      real(kind=RTYPE) b1(nxp,jlen,num,lev),b2(nxp,jlen,num,levp,nsizex)
      integer  nlen,j,i,k,ierr,jlistnum,num,n,i1,i2,j1

      b1=0.
      b2=0.
      aout=0.
      i2=0

      do k=1,lev
      do n=1,num
      do j=1,jlistnum
         j1=jlist1(j)
      do i=1,nxjlen(j1)
         b1(i,j,n,k)=ain(i,k,n,j)
      enddo
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen*num
      call MPI_ALLTOALL( b1 ,nlen, MPI_RTYPE, &
                         b2, nlen, MPI_RTYPE, &
                         comm, IERR )

      do n=1,num
      do k=1,levp
      do j=1,jlistnum
         j1=jlist1(j)
         i1=1
      do i=1,nsizex
         i2=nxjlen_all(i,j1)
         aout(i1:i1+i2-1,k,n,j)=b2(1:i2,j,n,k,i)
         i1=i1+i2
      enddo
      enddo
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------------------------------
!CWB2021 for single precision test

      subroutine mpe2d_transpose_nxp_lev_sp(ain,aout,nxp,nx,lev,levp,num,my,my_max,jlistnum,jlen,nsizex,comm)

! transpose (nx partial,lev full) to (nx full,lev partial), num variable packed

      use index, only : jlist1,nxjlen_all,nxjlen

      implicit none

      include 'mpif.h'
      integer  nx,nxp,lev,levp,my,my_max,jlen,nsizex,comm
      real*4   ain(nxp,lev,num,my_max),aout(nx,levp,num,my_max)
      real*4   b1(nxp,jlen,num,lev),b2(nxp,jlen,num,levp,nsizex)
      integer  nlen,j,i,k,ierr,jlistnum,num,n,i1,i2,j1

      b1=0.
      b2=0.
      aout=0.
      i2=0

      do k=1,lev
      do n=1,num
      do j=1,jlistnum
         j1=jlist1(j)
      do i=1,nxjlen(j1)
         b1(i,j,n,k)=ain(i,k,n,j)
      enddo
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen*num
      call MPI_ALLTOALL( b1 ,nlen, MPI_REAL4, &
                         b2, nlen, MPI_REAL4, &
                         comm, IERR )

      do n=1,num
      do k=1,levp
      do j=1,jlistnum
         j1=jlist1(j)
         i1=1
      do i=1,nsizex
         i2=nxjlen_all(i,j1)
         aout(i1:i1+i2-1,k,n,j)=b2(1:i2,j,n,k,i)
         i1=i1+i2
      enddo
      enddo
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_nxp_lev_ncld(ain,aout,nxp,nx,lev,levp,num,my,my_max,jlistnum,jlen,nsizex,comm,ncld)

! transpose (nx partial,lev full) to (nx full,lev partial), num variable packed,ncld version

      use index, only : jlist1,nxjlen_all,nxjlen

      implicit none

      include 'mpif.h'
      integer  nx,nxp,lev,levp,my,my_max,jlen,nsizex,comm
      real*8   ain(nxp,lev,ncld,num,my_max),aout(nx,levp,ncld,num,my_max)
      real*8   b1(nxp,jlen,num,ncld,lev),b2(nxp,jlen,num,ncld,levp,nsizex)
      integer  nlen,j,i,k,ierr,jlistnum,num,n,nn,ncld,i1,i2,j1

      b1=0.
      b2=0.
      aout=0.

      do j=1,jlistnum
         j1=jlist1(j)
      do nn=1,ncld
      do k=1,lev
      do i=1,nxjlen(j1)
         b1(i,j,1,nn,k)=ain(i,k,nn,1,j)
         b1(i,j,2,nn,k)=ain(i,k,nn,2,j)
      enddo
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen*num*ncld
      call MPI_ALLTOALL( b1 ,nlen, MPI_DOUBLE_PRECISION, &
                         b2, nlen, MPI_DOUBLE_PRECISION, &
                         comm, IERR )

      do k=1,levp
      do nn=1,ncld
      do j=1,jlistnum
         j1=jlist1(j)
         i1=1
      do i=1,nsizex
         i2=nxjlen_all(i,j1)
         aout(i1:i1+i2-1,k,nn,1,j)=b2(1:i2,j,1,nn,k,i)
         aout(i1:i1+i2-1,k,nn,2,j)=b2(1:i2,j,2,nn,k,i)
         i1=i1+i2
      enddo
      enddo
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------------------------------
      subroutine mpe2d_transpose1_nxp_lev(ain,aout,nx,nxp,lev,levp,my,my_max,jlistnum,jlen,nsizex,comm)

! transpose (nx partial,lev full) to (nx full,lev partial), only 1 variable 

      use index, only : lreduce,nxjlen_all,jlist1,nxjlen

      implicit none

      include 'mpif.h'
      integer  nx,nxp,lev,levp,my,my_max,jlen,nsizex,comm
      real*8   ain(nxp,lev,my_max),aout(nx,levp,my_max)
      real*8   b1(nxp,jlen,lev),b2(nxp,jlen,levp,nsizex)
      integer  nlen,j,i,k,ierr,jlistnum,nn,jj,i1,i2,j1

      b1=0.
      b2=0.
      aout=0.

      do k=1,lev
      do j=1,jlistnum
         j1=jlist1(j)
      do i=1,nxjlen(j1)
         b1(i,j,k)=ain(i,k,j)
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen
      call MPI_ALLTOALL( b1 ,nlen, MPI_DOUBLE_PRECISION, &
                         b2, nlen, MPI_DOUBLE_PRECISION, &
                         comm, IERR )

      do k=1,levp
      do j=1,jlistnum
         j1=jlist1(j)
         i1=1
      do i=1,nsizex
         i2=nxjlen_all(i,j1)
         aout(i1:i1+i2-1,k,j)=b2(1:i2,j,k,i)
         i1=i1+i2
      enddo
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_nx_levp(ain,aout,nxp,nx,lev,levp,num,my,my_max,jlistnum,jlen,nsizex,comm)

! transpose (nx full,lev partial) to (nx partial,lev full), num variables packed

      use index, only : jlist1,nxjlen_all
      use const, only : RTYPE,MPI_RTYPE

      implicit none

      include 'mpif.h'
      integer  nxp,nx,lev,levp,my,my_max,jlen,nsizex,comm
      real(kind=RTYPE)   ain(nx,levp,num,my_max),aout(nxp,lev,num,my_max)
      real(kind=RTYPE)   b1(levp,num,jlen,nxp,nsizex),b2(levp,num,jlen,nxp,nsizex)
      integer  nlen,j,jj,i,k,KL,ierr,jlistnum,num,n,i1,i2,j1

      b1=0.
      b2=0.
      aout=0.

      do j=1,jlistnum
         j1=jlist1(j)
      do n=1,num
      do k=1,levp
         i1=1
      do i=1,nsizex
         i2=nxjlen_all(i,j1)
         b1(k,n,j,1:i2,i)=ain(i1:i1+i2-1,k,n,j)
         i1=i1+i2
      enddo
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen*num
      call MPI_ALLTOALL( b1 ,nlen, MPI_RTYPE, &
                         b2, nlen, MPI_RTYPE, &
                         comm, IERR )
      do jj=1,jlistnum
      do n=1,num
      do i=1,nxp
         k=1
      do j=1,nsizex
         aout(i,k:k+levp-1,n,jj)=b2(1:levp,n,jj,i,j)
         k=k+levp
      enddo
      enddo
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_nx_levp_sp(ain,aout,nxp,nx,lev,levp,num,my,my_max,jlistnum,jlen,nsizex,comm)

! transpose (nx full,lev partial) to (nx partial,lev full), num variables packed

      use index, only : jlist1,nxjlen_all

      implicit none

      include 'mpif.h'
      integer  nxp,nx,lev,levp,my,my_max,jlen,nsizex,comm
      real*4   ain(nx,levp,num,my_max),aout(nxp,lev,num,my_max)
      real*4   b1(levp,num,jlen,nxp,nsizex),b2(levp,num,jlen,nxp,nsizex)
      integer  nlen,j,jj,i,k,KL,ierr,jlistnum,num,n,i1,i2,j1

      b1=0.
      b2=0.
      aout=0.

      do j=1,jlistnum
         j1=jlist1(j)
      do n=1,num
      do k=1,levp
         i1=1
      do i=1,nsizex
         i2=nxjlen_all(i,j1)
         b1(k,n,j,1:i2,i)=ain(i1:i1+i2-1,k,n,j)
         i1=i1+i2
      enddo
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen*num

      call MPI_ALLTOALL( b1 ,nlen, MPI_REAL4, &
                         b2, nlen, MPI_REAL4, &
                         comm, IERR )

      do jj=1,jlistnum
      do n=1,num
      do i=1,nxp
         k=1
      do j=1,nsizex
         aout(i,k:k+levp-1,n,jj)=b2(1:levp,n,jj,i,j)
         k=k+levp
      enddo
      enddo
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------------------------------
      subroutine mpe2d_reshape_eps4(eps4in,eps4out)

! reshape eps4, used in siimpl

      use mpi
      use param
      use index
      use const, only: RTYPE

      implicit none

      integer i,j,m,mf,nl

      real(kind=RTYPE) eps4in(jtrun,jtmax)
      real(kind=RTYPE) eps4out(jtp)
      real(kind=RTYPE) b1(jtf)

      b1=0.
      eps4out=0.

      i=1

      do m=1,mlistnum
         mf=mlist(m)
         NL=jtrun-mf+1
         b1(i:i+nl-1)=eps4in(mf:jtrun,m)
         i=i+NL
      enddo

      i=0
      do j=jtstart,jtend
         i=i+1
         eps4out(i)=b1(j)
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_get_allrow(n,a,a_all)

! gather all data of a row

      use index,only : row_comm
      use mpi

      integer a(*),a_all(*)

      call MPI_ALLGATHER( a,     n, MPI_INTEGER, &
                          a_all, n, MPI_INTEGER, &
                          row_comm,  IERR )
      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_row_broadcast(buf,n,brank)

! broadcast for row

      use mpi
      use index

      real*8 buf(n)
      integer brank

!     call MPI_BCAST( BUF, N, MPI_DOUBLE_PRECISION, 0, row_comm, IERR )
      call MPI_BCAST( BUF, N, MPI_DOUBLE_PRECISION, brank, row_comm, IERR )

      return
      end
!---------------------------------------------------------------------------------------
      subroutine mpe2d_unify(work,a,opt)

! unify a(nx_partial,my_partial) into work(nx_full,my_full)

      use param
      use rank
      use index
      use mpi
      use const, only: RTYPE,MPI_RTYPE

      real(kind=RTYPE) work(nx,my)
      real(kind=RTYPE) a(nxp,my_max)
      real(kind=RTYPE) b(nxp,my_max*nsize)
      logical, optional :: opt

      b=0.
      work=0.

      call MPI_ALLGATHER( a,nxp*my_max, MPI_RTYPE, &
                          b,nxp*my_max, MPI_RTYPE, &
                          MPI_COMM_gfs, IERR )

      do j=1,my
         ii=1
      do i=1,nsizex
         jf=jlist2_2d(i,j)
         nn=nxjlen_all(i,j)
         if(present(opt)) nn=nxp
         work(ii:ii+nn-1,j)=b(1:nn,jf)
         ii=ii+nn
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe2d_unify_nx(work,a)

! unify a(nx_partial,my_partial) to work(nx_full,my_partial)

      use param
      use index
      use mpi
      use const, only: RTYPE,MPI_RTYPE

      real(kind=RTYPE) work(nx,my_max)
      real(kind=RTYPE) a(nxp,my_max)
      real(kind=RTYPE) b(nxp,my_max,nsizex)

      work=0.
      b=0.

      call MPI_ALLGATHER( a,nxp*my_max, MPI_RTYPE, &
                          b,nxp*my_max, MPI_RTYPE, &
                          row_comm, IERR )

      do jj=1,jlistnum
         j=jlist1(jj)
         ii=1
      do i=1,nsizex
         nn=nxjlen_all(i,j)
         work(ii:ii+nn-1,jj)=b(1:nn,jj,i)
         ii=ii+nn
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe2d_unify_my(work,a)

! unify a(nx_full,my_partial) to work(nx_full,my_full)

      use param
      use index
      use mpi
      use const, only: RTYPE,MPI_RTYPE

      real(kind=RTYPE) work(nx,my)
      real(kind=RTYPE) a(nx,my_max)
      real(kind=RTYPE) b(nx,my_max*nsizey)

      work=0.
      b=0.

      call MPI_ALLGATHER( a,nx*my_max, MPI_RTYPE, &
                          b,nx*my_max, MPI_RTYPE, &
                          col_comm, IERR )


      do j=1,my
        jj=jlist2(j)
        work(1:nx,j)=b(1:nx,jj)
      enddo

      return
      end
!-------------------------------------------------------------------------
      subroutine mpe2d_unify_my1d(work,a)

! unify a(nx_full,my_partial) to work(nx_full,my_full)

      use param
      use index
      use mpi
      use const, only: RTYPE,MPI_RTYPE

      real(kind=RTYPE) work(my)
      real(kind=RTYPE) a(my_max)
      real(kind=RTYPE) b(my_max*nsizey)

      work=0.
      b=0.

      call MPI_ALLGATHER( a,my_max, MPI_RTYPE, &
                          b,my_max, MPI_RTYPE, &
                          col_comm, IERR )


      do j=1,my
        jj=jlist2(j)
        work(j)=b(jj)
      enddo

      return
      end
!----------------------------------------------------------------------------------------------
      subroutine mpe2d_unify_spec_lev(ain,aout,lev,levp,jtrun,jtmax,mlistnum,proc,comm)

! unify lev of spectrum
 
      use mpi
      use const, only : RTYPE,MPI_RTYPE
 
      implicit none

      integer lev,levp,jtrun,jtmax,mlistnum,proc,comm,ii,j,jj,k,kk,ierr

      real(kind=RTYPE) ain(levp,2,jtrun,jtmax)
      real(kind=RTYPE) aout(lev,2,jtrun,jtmax)
      real(kind=RTYPE) b1(levp,2,jtrun,jtmax)
      real(kind=RTYPE) b2(levp,2,jtrun,jtmax,proc)

      aout = 0.
      b2   = 0.

      call MPI_ALLGATHER( ain, levp*2*jtrun*jtmax, MPI_RTYPE, &
                           b2, levp*2*jtrun*jtmax, MPI_RTYPE, &
                          comm,IERR )

      do jj=1,mlistnum
      do j=1,jtrun
         kk=1
      do ii=1,proc
         aout(kk:kk+levp-1,1,j,jj)=b2(1:levp,1,j,jj,ii)
         aout(kk:kk+levp-1,2,j,jj)=b2(1:levp,2,j,jj,ii)
         kk=kk+levp
      enddo
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------------------------------
      subroutine mpe2d_unify_spec_lev_sp(ain,aout,lev,levp,jtrun,jtmax,mlistnum,proc,comm)

! unify lev of spectrum
 
      use mpi
 
      implicit none

      integer lev,levp,jtrun,jtmax,mlistnum,proc,comm,ii,j,jj,k,kk,ierr

      real*4 ain(levp,2,jtrun,jtmax)
      real*4 aout(lev,2,jtrun,jtmax)
      real*4 b1(levp,2,jtrun,jtmax)
      real*4 b2(levp,2,jtrun,jtmax,proc)

      aout = 0.
      b2   = 0.

      call MPI_ALLGATHER( ain, levp*2*jtrun*jtmax, MPI_REAL4, &
                           b2, levp*2*jtrun*jtmax, MPI_REAL4, &
                          comm,IERR )

      do jj=1,mlistnum
      do j=1,jtrun
         kk=1
      do ii=1,proc
         aout(kk:kk+levp-1,1,j,jj)=b2(1:levp,1,j,jj,ii)
         aout(kk:kk+levp-1,2,j,jj)=b2(1:levp,2,j,jj,ii)
         kk=kk+levp
      enddo
      enddo
      enddo

      return
      end
!-------------------------------------------------------------------------------------------------
      subroutine mpe2d_unify_spec_lev_zx(aout,a1,a2,a3,lev,levp,jtrun,jtmax,mlistnum,proc,comm)

! unify lev of spectrum for zx
 
      use mpi
      use const,only : RTYPE,MPI_RTYPE
 
      implicit none

      integer lev,levp,jtrun,jtmax,mlistnum,proc,comm,ii,j,jj,k,kk,n,ierr

      real(kind=RTYPE) a1(levp,2,jtrun,jtmax),a2(levp,2,jtrun,jtmax),  &
                       a3(levp,2,jtrun,jtmax)
      real(kind=RTYPE) ain(levp,2,3,jtrun,jtmax)
      real(kind=RTYPE) aout(lev,2,3,jtrun,jtmax)
      real(kind=RTYPE) b1(levp,2,jtrun,jtmax)
      real(kind=RTYPE) b2(levp,2,3,jtrun,jtmax,proc)

      ain  = 0.
      aout = 0.
      b2   = 0.

      ain(:,:,1,:,:)=a1(:,:,:,:)
      ain(:,:,2,:,:)=a2(:,:,:,:)
      ain(:,:,3,:,:)=a3(:,:,:,:)

      call MPI_ALLGATHER( ain, levp*2*3*jtrun*jtmax, MPI_RTYPE, &
                           b2, levp*2*3*jtrun*jtmax, MPI_RTYPE, &
                          comm,IERR )

      do jj=1,mlistnum
      do j=1,jtrun
         kk=1
      do ii=1,proc
         do n=1,3
         aout(kk:kk+levp-1,1,n,j,jj)=b2(1:levp,1,n,j,jj,ii)
         aout(kk:kk+levp-1,2,n,j,jj)=b2(1:levp,2,n,j,jj,ii)
         enddo
         kk=kk+levp
      enddo
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_unify_uzmean(a,my,lev)

! unify uzmean for 2dMPI

      use param, only : my_max
      use rank,  only : nsize, MPI_COMM_gfs
      use index
      use mpi

      implicit none

      real a(my,lev)
      real b1(levp,my_max)
      real b2(levp,my_max*nsize)
      integer i,j,ii,jj,j1,jf,kk,my,lev,ierr
 
      a=0.

      do jj=1,jlistnum
      do i=1,levp
         b1(i,jj)=a(jj,i)
      enddo
      enddo

      call MPI_ALLGATHER( b1,levp*my_max, MPI_DOUBLE_PRECISION, &
                          b2,levp*my_max, MPI_DOUBLE_PRECISION, &
                          MPI_COMM_gfs, IERR )

      do j=1,my
         kk=1
      do ii=1,nsizex
         jf=jlist2_2d(ii,j)
      do i=1,levp
         a(j,kk+i-1)=       b2(i,jf)
      enddo
         kk=kk+levp
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_reshape_pl(plin,plout)

! reshape spec coef 

      use mpi
      use param
      use index
      use const,only : RTYPE

      implicit none

      integer i,j,m,mf,nl

      real(kind=RTYPE) plin(jtrun,jtmax,2)
      real(kind=RTYPE) plout(jtp,2)
      real(kind=RTYPE) b1(jtf,2)

      plout=0.
      b1=0.
      i=1

      do m=1,mlistnum
         mf=mlist(m)
         NL=jtrun-mf+1
         b1(i:i+nl-1,1)=plin(mf:jtrun,m,1)
         b1(i:i+nl-1,2)=plin(mf:jtrun,m,2)
         i=i+NL
      enddo

      i=0
      do j=jtstart,jtend
         i=i+1
         plout(i,1)=b1(j,1)
         plout(i,2)=b1(j,2)
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_reshape_pl_multi(plin1,plin2,plin3,plout1,plout2,plout3)

! reshape spec coef (3 arrays)

      use mpi
      use param
      use index
      use const,only : RTYPE

      implicit none

      integer i,j,m,mf,nl

      real(kind=RTYPE) plin1(jtrun,jtmax,2),plin2(jtrun,jtmax,2),plin3(jtrun,jtmax,2)
      real(kind=RTYPE) plout1(jtp,2),plout2(jtp,2),plout3(jtp,2)
      real(kind=RTYPE) b1(jtf,6)

      b1=0.
      plout1=0.
      plout2=0.
      plout3=0.

      i=1

      do m=1,mlistnum
         mf=mlist(m)
         NL=jtrun-mf+1
         b1(i:i+nl-1,1)=plin1(mf:jtrun,m,1)
         b1(i:i+nl-1,2)=plin1(mf:jtrun,m,2)
         b1(i:i+nl-1,3)=plin2(mf:jtrun,m,1)
         b1(i:i+nl-1,4)=plin2(mf:jtrun,m,2)
         b1(i:i+nl-1,5)=plin3(mf:jtrun,m,1)
         b1(i:i+nl-1,6)=plin3(mf:jtrun,m,2)
         i=i+NL
      enddo

      i=0
      do j=jtstart,jtend
         i=i+1
         plout1(i,1)=b1(j,1)
         plout1(i,2)=b1(j,2)
         plout2(i,1)=b1(j,3)
         plout2(i,2)=b1(j,4)
         plout3(i,1)=b1(j,5)
         plout3(i,2)=b1(j,6)
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_reshape_pl_back(plin,plout)

! reshape spec coef back

      use mpi
      use param
      use index
      use const,only : RTYPE,MPI_RTYPE

      implicit none

      integer n,i,j,m,mf,ierr

      real(kind=RTYPE) plin(jtp,2)
      real(kind=RTYPE) plout(jtrun,jtmax,2)
      real(kind=RTYPE) b2(jtp,2,nsizex)

      plout=0.
      b2=0.

      call MPI_ALLGATHER( plin, jtp*2, MPI_RTYPE, &
                          b2,   jtp*2, MPI_RTYPE, &
                          row_comm, IERR )

      i=1
      j=1

      do m=1,mlistnum
         mf=mlist(m)
      do n=mf,jtrun
         plout(n,m,1)=b2(i,1,j)
         plout(n,m,2)=b2(i,2,j)
         i=i+1
         if(i.gt.jtlen_all(j))then
            i=1
            j=j+1
         endif
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_unify_nx_lev(a,nx,lev)

! unify nx of lev
 
      use index
      use mpi
 
      implicit none

      real a(nx,lev)
      real b1(nxp,lev)
      real b2(nxp,lev,nsizex)
      integer nx,lev,i,ii,k,IERR

      b1=0.
      b2=0.

      do k=1,lev
         b1(1:nxp,k)=a(1:nxp,k)
      enddo

      call MPI_ALLGATHER( b1, nxp*lev, MPI_DOUBLE_PRECISION, &
                          b2, nxp*lev, MPI_DOUBLE_PRECISION, &
                          row_comm, IERR )

      do k=1,lev
         ii=1
      do i=1,nsizex
         a(ii:ii+nxp-1,k)=b2(1:nxp,k,i)
         ii=ii+nxp
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_unify_nx_lev_red(a,nx,lev,j)

! unify nx of lev for reduce mode
 
      use index
      use mpi
      use const, only: RTYPE
 
      implicit none

      real a(nx,lev)
      real b1(nxp,lev)
      real b2(nxp,lev,nsizex)
      integer j,nx,lev,i,i2,ii,k,IERR

      b1=0.
      b2=0.

      do k=1,lev
         b1(1:nxp,k)=a(1:nxp,k)
      enddo

      call MPI_ALLGATHER( b1, nxp*lev, MPI_DOUBLE_PRECISION, &
                          b2, nxp*lev, MPI_DOUBLE_PRECISION, &
                          row_comm, IERR )

      do k=1,lev
         ii=1
      do i=1,nsizex
         i2=nxjlen_all(i,j)
         a(ii:ii+i2-1,k)=b2(1:i2,k,i)
         ii=ii+i2
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_siimpl_multi(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6, &
                 levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)

! transpose siimpl spec (6 arrays)

      use mpi
      use const, only: RTYPE,MPI_RTYPE

      implicit none

      integer levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist(jtrun),nsizex,row_comm
      integer m,mf,nl,i,levp2,nlen,ierr,j,k

      real(kind=RTYPE) a1(levp,2,jtrun,jtmax)
      real(kind=RTYPE) a2(levp,2,jtrun,jtmax)
      real(kind=RTYPE) a3(levp,2,jtrun,jtmax)
      real(kind=RTYPE) a4(levp,2,jtrun,jtmax)
      real(kind=RTYPE) a5(levp,2,jtrun,jtmax)
      real(kind=RTYPE) a6(levp,2,jtrun,jtmax)
      real(kind=RTYPE) b1(lev,2,jtp)
      real(kind=RTYPE) b2(lev,2,jtp)
      real(kind=RTYPE) b3(lev,2,jtp)
      real(kind=RTYPE) b4(lev,2,jtp)
      real(kind=RTYPE) b5(lev,2,jtp)
      real(kind=RTYPE) b6(lev,2,jtp)
      real(kind=RTYPE) c1(levp,2,6,jtf)
      real(kind=RTYPE) c2(levp,2,6,jtp,nsizex)

      b1=0.
      b2=0.
      b3=0.
      b4=0.
      b5=0.
      b6=0.
      c1=0.
      c2=0.

      levp2=levp*2
      i=1

      do m=1,mlistnum
         mf=mlist(m)
         NL=jtrun-mf+1
         c1(1:levp,1,1,i:i+nl-1)=a1(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,1,i:i+nl-1)=a1(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,2,i:i+nl-1)=a2(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,2,i:i+nl-1)=a2(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,3,i:i+nl-1)=a3(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,3,i:i+nl-1)=a3(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,4,i:i+nl-1)=a4(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,4,i:i+nl-1)=a4(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,5,i:i+nl-1)=a5(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,5,i:i+nl-1)=a5(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,6,i:i+nl-1)=a6(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,6,i:i+nl-1)=a6(1:levp,2,mf:jtrun,m)
         i=i+NL
      enddo

      nlen=levp2*6*jtp
      call MPI_ALLTOALL( c1 ,nlen, MPI_RTYPE, &
                         c2, nlen, MPI_RTYPE, &
                         row_comm, IERR )

      do j=1,jtp
         k=1
      do i=1,nsizex
         b1(k:k+levp-1,1,j)=c2(1:levp,1,1,j,i)
         b1(k:k+levp-1,2,j)=c2(1:levp,2,1,j,i)
         b2(k:k+levp-1,1,j)=c2(1:levp,1,2,j,i)
         b2(k:k+levp-1,2,j)=c2(1:levp,2,2,j,i)
         b3(k:k+levp-1,1,j)=c2(1:levp,1,3,j,i)
         b3(k:k+levp-1,2,j)=c2(1:levp,2,3,j,i)
         b4(k:k+levp-1,1,j)=c2(1:levp,1,4,j,i)
         b4(k:k+levp-1,2,j)=c2(1:levp,2,4,j,i)
         b5(k:k+levp-1,1,j)=c2(1:levp,1,5,j,i)
         b5(k:k+levp-1,2,j)=c2(1:levp,2,5,j,i)
         b6(k:k+levp-1,1,j)=c2(1:levp,1,6,j,i)
         b6(k:k+levp-1,2,j)=c2(1:levp,2,6,j,i)
         k=k+levp
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_siimpl(ain,aout, &
                 levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)

! transpose siimpl spec

      use mpi
      use const, only: RTYPE,MPI_RTYPE

      implicit none

      integer levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist(jtrun),nsizex,row_comm
      integer m,mf,nl,i,levp2,nlen,ierr,j,k

      real(kind=RTYPE) ain(levp,2,jtrun,jtmax)
      real(kind=RTYPE) aout(lev,2,jtp)
      real(kind=RTYPE) c1(levp,2,jtf)
      real(kind=RTYPE) c2(levp,2,jtp,nsizex)

      aout=0.
      c1=0.
      c2=0.

      levp2=levp*2
      i=1

      do m=1,mlistnum
         mf=mlist(m)
         NL=jtrun-mf+1
         c1(1:levp,1,i:i+nl-1)=ain(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,i:i+nl-1)=ain(1:levp,2,mf:jtrun,m)
         i=i+NL
      enddo

      nlen=levp2*jtp
      call MPI_ALLTOALL( c1 ,nlen, MPI_RTYPE, &
                         c2, nlen, MPI_RTYPE, &
                         row_comm, IERR )

      do j=1,jtp
         k=1
      do i=1,nsizex
         aout(k:k+levp-1,1,j)=c2(1:levp,1,j,i)
         aout(k:k+levp-1,2,j)=c2(1:levp,2,j,i)
         k=k+levp
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_siimpl_back(ain,aout,levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)

! transpose siimpl spec back

      use mpi
      use index,only : jtlen_all
      use const,only : RTYPE,MPI_RTYPE

      implicit none

      integer levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist(jtrun),nsizex,row_comm
      integer m,n,mf,nl,i,nlen,ierr,j,k

      real(kind=RTYPE) ain(lev,2,jtp)
      real(kind=RTYPE) aout(levp,2,jtrun,jtmax)

      real(kind=RTYPE) b1(jtp,2,lev)
      real(kind=RTYPE) b2(jtp,2,levp,nsizex)

      b1=0.
      b2=0.
      aout=0.

      do j=1,jtp
      do k=1,lev
         b1(j,1,k)=ain(k,1,j)
         b1(j,2,k)=ain(k,2,j)
      enddo
      enddo

      nlen=jtp*2*levp
      call MPI_ALLTOALL( b1 ,nlen, MPI_RTYPE, &
                         b2, nlen, MPI_RTYPE, &
                         row_comm, IERR )

      i=1
      j=1

      do m=1,mlistnum
         mf=mlist(m)
      do n=mf,jtrun
      do k=1,levp
         aout(k,1,n,m)=b2(i,1,k,j)
         aout(k,2,n,m)=b2(i,2,k,j)
      enddo
         i=i+1
         if(i.gt.jtlen_all(j))then
            i=1
            j=j+1
         endif
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_siimpl_back_multi(ain1,ain2,aout1,aout2,levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)

! transpose siimpl spec back (2 arrays)

      use mpi
      use index,only : jtlen_all
      use const,only : RTYPE,MPI_RTYPE

      implicit none

      integer levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist(jtrun),nsizex,row_comm
      integer m,n,mf,nl,i,nlen,ierr,j,k

      real(kind=RTYPE) ain1(lev,2,jtp),ain2(lev,2,jtp)
      real(kind=RTYPE) aout1(levp,2,jtrun,jtmax),aout2(levp,2,jtrun,jtmax)

      real(kind=RTYPE) b1(jtp,4,lev)
      real(kind=RTYPE) b2(jtp,4,levp,nsizex)
   
      b1=0.
      b2=0.
      aout1=0.
      aout2=0. 

      do j=1,jtp
      do k=1,lev
         b1(j,1,k)=ain1(k,1,j)
         b1(j,2,k)=ain1(k,2,j)
         b1(j,3,k)=ain2(k,1,j)
         b1(j,4,k)=ain2(k,2,j)
      enddo
      enddo

      nlen=jtp*4*levp
      call MPI_ALLTOALL( b1 ,nlen, MPI_RTYPE, &
                         b2, nlen, MPI_RTYPE, &
                         row_comm, IERR )

      i=1
      j=1

      do m=1,mlistnum
         mf=mlist(m)
      do n=mf,jtrun
      do k=1,levp
         aout1(k,1,n,m)=b2(i,1,k,j)
         aout1(k,2,n,m)=b2(i,2,k,j)
         aout2(k,1,n,m)=b2(i,3,k,j)
         aout2(k,2,n,m)=b2(i,4,k,j)
      enddo
         i=i+1
         if(i.gt.jtlen_all(j))then
            i=1
            j=j+1
         endif
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_unify_lev(ain,aout,lev,levp,jtrun,jtmax,mlistnum,proc,comm)

! unify lev of spectrum var

      use mpi
      use const,only : RTYPE,MPI_RTYPE

      implicit none

      integer lev,levp,jtrun,jtmax,mlistnum,proc,comm,ii,j,jj,k,kk,ierr

      real(kind=RTYPE) ain(levp,2,jtrun,jtmax)
      real(kind=RTYPE) aout(lev,2,jtrun,jtmax)
      real(kind=RTYPE) b1(levp,2,jtrun,jtmax)
      real(kind=RTYPE) b2(levp,2,jtrun,jtmax,proc)

      b2=0.
      aout=0.

      call MPI_ALLGATHER( ain, levp*2*jtrun*jtmax, MPI_RTYPE, &
                           b2, levp*2*jtrun*jtmax, MPI_RTYPE, &
                          comm,IERR )

      do jj=1,mlistnum
      do j=1,jtrun
         kk=1
      do ii=1,proc
         aout(kk:kk+levp-1,1,j,jj)=b2(1:levp,1,j,jj,ii)
         aout(kk:kk+levp-1,2,j,jj)=b2(1:levp,2,j,jj,ii)
         kk=kk+levp
      enddo
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_ndsl_p2f_multi(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,nxp,nx,lev, &
                                          levp,ncld,my,my_max,jlistnum,jlen,nsizex,comm,num)

! transpose (nx partial,lev full) to (nx full,lev partial) for NDSL

      use index, only : lreduce,nxjlen_all,jlist1,nxjlen

      implicit none

      include 'mpif.h'

      real*8   a1(nxp,lev,my_max),a2(nxp,lev,my_max),a3(nxp,lev,my_max),     &
               a4(nxp,lev,my_max),a5(nxp,lev,my_max),a6(nxp,lev,ncld,my_max)

      real*8   b1(nx,levp,my_max),b2(nx,levp,my_max),b3(nx,levp,my_max),     &
               b4(nx,levp,my_max),b5(nx,levp,my_max),b6(nx,levp,ncld,my_max)

      real*8   c1(nxp,jlen,5+ncld,lev)
      real*8   c2(nxp,jlen,5+ncld,levp,nsizex)
      real*8   d1(nxp,jlen,3+ncld,lev)
      real*8   d2(nxp,jlen,3+ncld,levp,nsizex)

      integer  nxp,nx,lev,levp,ncld,my,my_max,jlen,nsizex,comm,num
      integer  nlen,ii,j,i,k,kk,ierr,jlistnum,nn,jj,n

      b1=0.
      b2=0.
      b3=0.
      b4=0.
      b5=0.
      b6=0.
      c1=0.
      c2=0.
      d1=0.
      d2=0.

      if(num.eq.6)then  ! 6 vars in intgrt

      do k=1,lev
         kk=lev-k+1
      do j=1,jlistnum
         jj=jlist1(j)
      do i=1,nxjlen(jj)
         c1(i,j,1,kk)=a1(i,k,j)
         c1(i,j,2,kk)=a2(i,k,j)
         c1(i,j,3,kk)=a3(i,k,j)
         c1(i,j,4,kk)=a4(i,k,j)
         c1(i,j,5,kk)=a5(i,k,j)
         do n=1,ncld
            c1(i,j,5+n,kk)=a6(i,k,n,j)
         enddo
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen*(5+ncld)

      call MPI_ALLTOALL( c1 ,nlen, MPI_DOUBLE_PRECISION, &
                         c2, nlen, MPI_DOUBLE_PRECISION, &
                         comm, IERR )

      do k=1,levp
      do j=1,jlistnum
         jj=jlist1(j)
         ii=1
      do i=1,nsizex
         nn=nxjlen_all(i,jj)
         b1(ii:ii+nn-1,k,j)=c2(1:nn,j,1,k,i)
         b2(ii:ii+nn-1,k,j)=c2(1:nn,j,2,k,i)
         b3(ii:ii+nn-1,k,j)=c2(1:nn,j,3,k,i)
         b4(ii:ii+nn-1,k,j)=c2(1:nn,j,4,k,i)
         b5(ii:ii+nn-1,k,j)=c2(1:nn,j,5,k,i)
         do n=1,ncld
            b6(ii:ii+nn-1,k,n,j)=c2(1:nn,j,5+n,k,i)
         enddo
         ii=ii+nn
      enddo
      enddo
      enddo

      else              ! 4 vars in tendget

      do k=1,lev
         kk=lev-k+1
      do j=1,jlistnum
         jj=jlist1(j)
      do i=1,nxjlen(jj)
         d1(i,j,1,kk)=a1(i,k,j)
         d1(i,j,2,kk)=a2(i,k,j)
         d1(i,j,3,kk)=a3(i,k,j)
         do n=1,ncld
            d1(i,j,3+n,kk)=a6(i,k,n,j)
         enddo
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen*(3+ncld)

      call MPI_ALLTOALL( d1 ,nlen, MPI_DOUBLE_PRECISION, &
                         d2, nlen, MPI_DOUBLE_PRECISION, &
                         comm, IERR )

      do k=1,levp
      do j=1,jlistnum
         jj=jlist1(j)
         ii=1
      do i=1,nsizex
         nn=nxjlen_all(i,jj)
         b1(ii:ii+nn-1,k,j)=d2(1:nn,j,1,k,i)
         b2(ii:ii+nn-1,k,j)=d2(1:nn,j,2,k,i)
         b3(ii:ii+nn-1,k,j)=d2(1:nn,j,3,k,i)
         do n=1,ncld
            b6(ii:ii+nn-1,k,n,j)=d2(1:nn,j,3+n,k,i)
         enddo
         ii=ii+nn
      enddo
      enddo
      enddo

      endif

      return
      end

!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_ndsl_f2p_multi(a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,nxp,nx,lev,  &
                                          levp,ncld,my,my_max,jlistnum,jlen,proc,comm)

! transpose (nx full,lev partial) to (nx partial,lev full) for NDSL

      use index, only : jlist1,nsizex,lreduce,nxjstart_all,nxjend_all,nxjlen_all,nxjp

      implicit none

      include 'mpif.h'

      integer  nxp,nx,lev,levp,ncld,my,my_max,jlen,proc,comm
      integer  nlen,ii,j,jj,i,k,KL,ierr,jlistnum,num,n,i1,i2,i3,i4

      real*8   a1(nx,levp,my_max),a2(nx,levp,my_max),a3(nx,levp,my_max),     &
               a4(nx,levp,my_max),a5(nx,levp,ncld,my_max)

      real*8   b1(nxp,lev,my_max),b2(nxp,lev,my_max),b3(nxp,lev,my_max),     &
               b4(nxp,lev,my_max),b5(nxp,lev,ncld,my_max)

      real*8   c1(levp,4+ncld,jlen,nxp,proc)
      real*8   c2(levp,4+ncld,jlen,nxp,proc)

      b1=0.
      b2=0.
      b3=0.
      b4=0.
      b5=0.
      c1=0.
      c2=0.

      do j=1,jlistnum
         jj=jlist1(j)
      do k=1,levp
         i1=1
      do i=1,nsizex
         i2=nxjlen_all(i,jj)
         c1(k,1,j,1:i2,i)=a1(i1:i1+i2-1,k,j)
         c1(k,2,j,1:i2,i)=a2(i1:i1+i2-1,k,j)
         c1(k,3,j,1:i2,i)=a3(i1:i1+i2-1,k,j)
         c1(k,4,j,1:i2,i)=a4(i1:i1+i2-1,k,j)
         do n=1,ncld
            c1(k,4+n,j,1:i2,i)=a5(i1:i1+i2-1,k,n,j)
         enddo
         i1=i1+i2
      enddo
      enddo
      enddo


      nlen=nxp*levp*jlen*(4+ncld)

      call MPI_ALLTOALL( c1 ,nlen, MPI_DOUBLE_PRECISION, &
                         c2, nlen, MPI_DOUBLE_PRECISION, &
                         comm, IERR )

      do jj=1,jlistnum

      ii=jlist1(jj)
      do i=1,nxjp(ii)

         KL=lev
      do j=1,proc
      do k=1,levp
         b1(i,KL,jj)=c2(k,1,jj,i,j)
         b2(i,KL,jj)=c2(k,2,jj,i,j)
         b3(i,KL,jj)=c2(k,3,jj,i,j)
         b4(i,KL,jj)=c2(k,4,jj,i,j)
         do n=1,ncld
            b5(i,KL,n,jj)=c2(k,4+n,jj,i,j)
         enddo
         KL=KL-1
      enddo
      enddo
      enddo
      enddo

      return
      end
!-----------------------------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_ndsl_p2f(ain,aout,nxp,nx,lev,levp,ncld,my,my_max,jlistnum,jlen,nsizex,comm)

! transpose (nx partial,lev full) to (nx full,lev partial) for NDSL

      use index, only : lreduce,nxjlen_all,jlist1,nxjlen
      use const, only : RTYPE,MPI_RTYPE

      implicit none

      include 'mpif.h'

      real(kind=RTYPE) ain(nxp,lev,ncld,my_max),        &
                       aout(nx,levp,ncld,my_max),       &
                       c1(nxp,jlen,ncld,lev),           &
                       c2(nxp,jlen,ncld,levp,nsizex)

      integer  nxp,nx,lev,levp,ncld,my,my_max,jlen,nsizex,comm
      integer  nlen,ii,j,i,k,kk,ierr,jlistnum,nn,jj,n

      c1=0.
      c2=0.
      aout=0.

      do k=1,lev
         kk=lev-k+1
      do j=1,jlistnum
         jj=jlist1(j)
      do i=1,nxjlen(jj)
         do n=1,ncld
            c1(i,j,n,kk)=ain(i,k,n,j)
         enddo
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen*ncld

      call MPI_ALLTOALL( c1 ,nlen, MPI_RTYPE, &
                         c2, nlen, MPI_RTYPE, &
                         comm, IERR )

      do k=1,levp
      do j=1,jlistnum
         jj=jlist1(j)
         ii=1
      do i=1,nsizex
         nn=nxjlen_all(i,jj)
         do n=1,ncld
            aout(ii:ii+nn-1,k,n,j)=c2(1:nn,j,n,k,i)
         enddo
         ii=ii+nn
      enddo
      enddo
      enddo


      return
      end
!---------------------------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_ndsl_f2p(ain,aout,nxp,nx,lev,levp,ncld,my,my_max,jlistnum,jlen,proc,comm)

! transpose (nx full,lev partial) to (nx partial,lev full) for NDSL

      use index, only : jlist1,nsizex,lreduce,nxjstart_all,nxjend_all,nxjlen_all,nxjp
      use const, only : RTYPE,MPI_RTYPE

      implicit none

      include 'mpif.h'

      integer  nxp,nx,lev,levp,ncld,my,my_max,jlen,proc,comm
      real(kind=RTYPE) ain(nx,levp,ncld,my_max),      &
                       aout(nxp,lev,ncld,my_max),     &
                       c1(levp,ncld,jlen,nxp,proc),   &
                       c2(levp,ncld,jlen,nxp,proc)

      integer  nlen,ii,j,jj,i,k,KL,ierr,jlistnum,n,i1,i2,i3,i4

      c1=0.
      c2=0.
      aout=0.

      do j=1,jlistnum
         jj=jlist1(j)
      do k=1,levp
         i1=1
      do i=1,nsizex
         i2=nxjlen_all(i,jj)
         do n=1,ncld
            c1(k,n,j,1:i2,i)=ain(i1:i1+i2-1,k,n,j)
         enddo
         i1=i1+i2
      enddo
      enddo
      enddo

      nlen=nxp*levp*jlen*ncld

      call MPI_ALLTOALL( c1 ,nlen, MPI_RTYPE, &
                         c2, nlen, MPI_RTYPE, &
                         comm, IERR )

      do jj=1,jlistnum

      ii=jlist1(jj)
      do i=1,nxjp(ii)

         KL=lev
      do j=1,proc
      do k=1,levp
         do n=1,ncld
            aout(i,KL,n,jj)=c2(k,n,jj,i,j)
         enddo
         KL=KL-1
      enddo
      enddo
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
!CWB2021
      subroutine mpe2d_reshape_pl_sp(plin,plout)

! reshape spec coef 

      use mpi
      use param
      use index

      implicit none

      integer i,j,m,mf,nl

      real*4 plin(jtrun,jtmax,2)
      real*4 plout(jtp,2)
      real*4 b1(jtf,2)

      b1=0.
      plout=0.

      i=1

      do m=1,mlistnum
         mf=mlist(m)
         NL=jtrun-mf+1
         b1(i:i+nl-1,1)=plin(mf:jtrun,m,1)
         b1(i:i+nl-1,2)=plin(mf:jtrun,m,2)
         i=i+NL
      enddo

      i=0
      do j=jtstart,jtend
         i=i+1
         plout(i,1)=b1(j,1)
         plout(i,2)=b1(j,2)
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_reshape_pl_multi_sp(plin1,plin2,plin3,plout1,plout2,plout3)

! reshape spec coef (3 arrays)

      use mpi
      use param
      use index

      implicit none

      integer i,j,m,mf,nl

      real*4 plin1(jtrun,jtmax,2),plin2(jtrun,jtmax,2),plin3(jtrun,jtmax,2)
      real*4 plout1(jtp,2),plout2(jtp,2),plout3(jtp,2)
      real*4 b1(jtf,6)

      b1=0.
      plout1=0.
      plout2=0.
      plout3=0.

      i=1

      do m=1,mlistnum
         mf=mlist(m)
         NL=jtrun-mf+1
         b1(i:i+nl-1,1)=plin1(mf:jtrun,m,1)
         b1(i:i+nl-1,2)=plin1(mf:jtrun,m,2)
         b1(i:i+nl-1,3)=plin2(mf:jtrun,m,1)
         b1(i:i+nl-1,4)=plin2(mf:jtrun,m,2)
         b1(i:i+nl-1,5)=plin3(mf:jtrun,m,1)
         b1(i:i+nl-1,6)=plin3(mf:jtrun,m,2)
         i=i+NL
      enddo

      i=0
      do j=jtstart,jtend
         i=i+1
         plout1(i,1)=b1(j,1)
         plout1(i,2)=b1(j,2)
         plout2(i,1)=b1(j,3)
         plout2(i,2)=b1(j,4)
         plout3(i,1)=b1(j,5)
         plout3(i,2)=b1(j,6)
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_reshape_pl_back_sp(plin,plout)

! reshape spec coef back

      use mpi
      use param
      use index

      implicit none

      integer n,i,j,m,mf,ierr

      real*4 plin(jtp,2)
      real*4 plout(jtrun,jtmax,2)
      real*4 b2(jtp,2,nsizex)

      b2=0.
      plout=0.

      call MPI_ALLGATHER( plin, jtp*2, MPI_REAL4, &
                          b2,   jtp*2, MPI_REAL4, &
                          row_comm, IERR )

      i=1
      j=1

      do m=1,mlistnum
         mf=mlist(m)
      do n=mf,jtrun
         plout(n,m,1)=b2(i,1,j)
         plout(n,m,2)=b2(i,2,j)
         i=i+1
         if(i.gt.jtlen_all(j))then
            i=1
            j=j+1
         endif
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_siimpl_multi_sp(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6, &
                 levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)

! transpose siimpl spec (6 arrays)

      use mpi

      implicit none

      integer levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist(jtrun),nsizex,row_comm
      integer m,mf,nl,i,levp2,nlen,ierr,j,k

      real*4 a1(levp,2,jtrun,jtmax)
      real*4 a2(levp,2,jtrun,jtmax)
      real*4 a3(levp,2,jtrun,jtmax)
      real*4 a4(levp,2,jtrun,jtmax)
      real*4 a5(levp,2,jtrun,jtmax)
      real*4 a6(levp,2,jtrun,jtmax)
      real*4 b1(lev,2,jtp)
      real*4 b2(lev,2,jtp)
      real*4 b3(lev,2,jtp)
      real*4 b4(lev,2,jtp)
      real*4 b5(lev,2,jtp)
      real*4 b6(lev,2,jtp)
      real*4 c1(levp,2,6,jtf)
      real*4 c2(levp,2,6,jtp,nsizex)

      b1=0.
      b2=0.
      b3=0.
      b4=0.
      b5=0.
      b6=0.
      c1=0.
      c2=0.

      levp2=levp*2
      i=1

      do m=1,mlistnum
         mf=mlist(m)
         NL=jtrun-mf+1
         c1(1:levp,1,1,i:i+nl-1)=a1(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,1,i:i+nl-1)=a1(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,2,i:i+nl-1)=a2(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,2,i:i+nl-1)=a2(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,3,i:i+nl-1)=a3(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,3,i:i+nl-1)=a3(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,4,i:i+nl-1)=a4(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,4,i:i+nl-1)=a4(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,5,i:i+nl-1)=a5(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,5,i:i+nl-1)=a5(1:levp,2,mf:jtrun,m)
         c1(1:levp,1,6,i:i+nl-1)=a6(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,6,i:i+nl-1)=a6(1:levp,2,mf:jtrun,m)
         i=i+NL
      enddo

      nlen=levp2*6*jtp

      call MPI_ALLTOALL( c1 ,nlen, MPI_REAL4, &
                         c2, nlen, MPI_REAL4, &
                         row_comm, IERR )

      do j=1,jtp
         k=1
      do i=1,nsizex
         b1(k:k+levp-1,1,j)=c2(1:levp,1,1,j,i)
         b1(k:k+levp-1,2,j)=c2(1:levp,2,1,j,i)
         b2(k:k+levp-1,1,j)=c2(1:levp,1,2,j,i)
         b2(k:k+levp-1,2,j)=c2(1:levp,2,2,j,i)
         b3(k:k+levp-1,1,j)=c2(1:levp,1,3,j,i)
         b3(k:k+levp-1,2,j)=c2(1:levp,2,3,j,i)
         b4(k:k+levp-1,1,j)=c2(1:levp,1,4,j,i)
         b4(k:k+levp-1,2,j)=c2(1:levp,2,4,j,i)
         b5(k:k+levp-1,1,j)=c2(1:levp,1,5,j,i)
         b5(k:k+levp-1,2,j)=c2(1:levp,2,5,j,i)
         b6(k:k+levp-1,1,j)=c2(1:levp,1,6,j,i)
         b6(k:k+levp-1,2,j)=c2(1:levp,2,6,j,i)
         k=k+levp
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_siimpl_sp(ain,aout, &
                 levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)

! transpose siimpl spec

      use mpi

      implicit none

      integer levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist(jtrun),nsizex,row_comm
      integer m,mf,nl,i,levp2,nlen,ierr,j,k

      real*4 ain(levp,2,jtrun,jtmax)
      real*4 aout(lev,2,jtp)
      real*4 c1(levp,2,jtf)
      real*4 c2(levp,2,jtp,nsizex)

      c1=0.
      c2=0.
      aout=0.

      levp2=levp*2
      i=1

      do m=1,mlistnum
         mf=mlist(m)
         NL=jtrun-mf+1
         c1(1:levp,1,i:i+nl-1)=ain(1:levp,1,mf:jtrun,m)
         c1(1:levp,2,i:i+nl-1)=ain(1:levp,2,mf:jtrun,m)
         i=i+NL
      enddo

      nlen=levp2*jtp

      call MPI_ALLTOALL( c1 ,nlen, MPI_REAL4, &
                         c2, nlen, MPI_REAL4, &
                         row_comm, IERR )

      do j=1,jtp
         k=1
      do i=1,nsizex
         aout(k:k+levp-1,1,j)=c2(1:levp,1,j,i)
         aout(k:k+levp-1,2,j)=c2(1:levp,2,j,i)
         k=k+levp
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_siimpl_back_sp(ain,aout,levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)

! transpose siimpl spec back

      use mpi
      use index,only : jtlen_all

      implicit none

      integer levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist(jtrun),nsizex,row_comm
      integer m,n,mf,nl,i,nlen,ierr,j,k

      real*4 ain(lev,2,jtp)
      real*4 aout(levp,2,jtrun,jtmax)
      real*4 b1(jtp,2,lev)
      real*4 b2(jtp,2,levp,nsizex)

      b1=0.
      b2=0.
      aout=0.

      do j=1,jtp
      do k=1,lev
         b1(j,1,k)=ain(k,1,j)
         b1(j,2,k)=ain(k,2,j)
      enddo
      enddo

      nlen=jtp*2*levp

      call MPI_ALLTOALL( b1 ,nlen, MPI_REAL4, &
                         b2, nlen, MPI_REAL4, &
                         row_comm, IERR )

      i=1
      j=1

      do m=1,mlistnum
         mf=mlist(m)
      do n=mf,jtrun
      do k=1,levp
         aout(k,1,n,m)=b2(i,1,k,j)
         aout(k,2,n,m)=b2(i,2,k,j)
      enddo
         i=i+1
         if(i.gt.jtlen_all(j))then
            i=1
            j=j+1
         endif
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
      subroutine mpe2d_transpose_siimpl_back_multi_sp(ain1,ain2,aout1,aout2,levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)

! transpose siimpl spec back (2 arrays)

      use mpi
      use index,only : jtlen_all

      implicit none

      integer levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist(jtrun),nsizex,row_comm
      integer m,n,mf,nl,i,nlen,ierr,j,k

      real*4 ain1(lev,2,jtp),ain2(lev,2,jtp)
      real*4 aout1(levp,2,jtrun,jtmax),aout2(levp,2,jtrun,jtmax)

      real*4 b1(jtp,4,lev)
      real*4 b2(jtp,4,levp,nsizex)

      b1=0.
      b2=0.
      aout1=0.
      aout2=0.
      nlen=0

      do j=1,jtp
      do k=1,lev
         b1(j,1,k)=ain1(k,1,j)
         b1(j,2,k)=ain1(k,2,j)
         b1(j,3,k)=ain2(k,1,j)
         b1(j,4,k)=ain2(k,2,j)
      enddo
      enddo

      nlen=jtp*4*levp

      call MPI_ALLTOALL( b1 ,nlen, MPI_REAL4, &
                         b2, nlen, MPI_REAL4, &
                         row_comm, IERR )

      i=1
      j=1

      do m=1,mlistnum
         mf=mlist(m)
      do n=mf,jtrun
      do k=1,levp
         aout1(k,1,n,m)=b2(i,1,k,j)
         aout1(k,2,n,m)=b2(i,2,k,j)
         aout2(k,1,n,m)=b2(i,3,k,j)
         aout2(k,2,n,m)=b2(i,4,k,j)
      enddo
         i=i+1
         if(i.gt.jtlen_all(j))then
            i=1
            j=j+1
         endif
      enddo
      enddo

      return
      end
!----------------------------------------------------------------------------------------
