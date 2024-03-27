      subroutine mpe2d_transpose_nx_levp_gpu(ain,aout,nxp,nx,lev,levp,num,my,my_max,jlistnum,jlen,nsizex,comm)

! transpose (nx full,lev partial) to (nx partial,lev full), num variables packed

      use index, only : jlist1,nxjlen_all
      use const, only : RTYPE,MPI_RTYPE

      implicit none

      include 'mpif.h'
      integer  nxp,nx,lev,levp,my,my_max,jlen,nsizex,comm
      real(kind=RTYPE)   ain(nx,levp,num,my_max),aout(nxp,lev,num,my_max)
      real(kind=RTYPE)   b1(levp,num,jlen,nxp,nsizex),b2(levp,num,jlen,nxp,nsizex)
      integer  nlen,j,jj,i,k,KL,ierr,jlistnum,num,n,i1,i2,j1

      !$acc data create(b1, b2) copyin(jlist1, nxjlen_all)
      !$acc kernels present(b1, b2)
      b1=0.
      b2=0.
      aout=0.
      !$acc end kernels

      !$acc parallel loop gang present(jlist1, nxjlen_all, b1, ain)
      do j=1,jlistnum
         j1=jlist1(j)
      !$acc loop vector collapse(2)
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
      !$acc host_data use_device(b1, b2)
      call MPI_ALLTOALL( b1 ,nlen, MPI_RTYPE, &
                         b2, nlen, MPI_RTYPE, &
                         comm, IERR )
      !$acc end host_data
      !$acc parallel loop collapse(2) present(aout, b2)
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
      !$acc end data

      return
      end
