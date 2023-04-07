      subroutine para_we2ns(a,b,levs,latg)
!
! mpi transport from full dimension of west-east to full dimension of
! north-south with latitude shuffl. created by hann-ming henry juang
!
! program log
! 2011 02 20 : henry juang, created for ndsl advection
!
!
      use const   , only : RTYPE  ,MPI_RTYPE
      use grid    , only : lonfull,lonhalf,lonpart,lonlenmax,mylonlen, &
                           latfull,lathalf,latpart,latlenmax,mylatlen, &
                           latstr ,latlen ,lonstr ,lonlen
      use rank    
      use index
      use mpi
!      use gfs_dyn_mpi_def
!      use gfs_dyn_layout1
      implicit none
!
      integer levs,latg
      real(kind=RTYPE) a(lonfull,levs,latpart)
      real(kind=RTYPE) b(latfull,levs,lonpart)
!      integer global_lats_a(latg)
!
!     real(kind=RTYPE) (kind=kind_mpi_r) works(2*levs*mylatlen*lonhalf)
!     real(kind=RTYPE) (kind=kind_mpi_r) workr(2*levs*mylonlen*lathalf)
      real(kind=RTYPE) works(2,levs,lonlenmax*latlenmax,nsizey)
      real(kind=RTYPE) workr(2,levs,lonlenmax*latlenmax,nsizey)
      integer lensend(nsizey),lenrecv(nsizey)
      integer locsend(nsizey),locrecv(nsizey)
!      integer i,j,k,n,mn,jj,lat1,lat2,ierr
      integer i,j,k,n,mn,lat1,lat2,ierr
!
!     print *,' enter ndslfv_we2ns '

!$omp parallel do private(n,mn,i,j,k) &
!$omp schedule(dynamic)
      do n=1,nsizey
        mn=0
        do j=1,mylatlen
          do i=lonstr(n),lonstr(n)+lonlen(n)-1
            mn = mn + 1
            do k=1,levs
              works(1,k,mn,n)=a(i        ,k,j)
              works(2,k,mn,n)=a(i+lonhalf,k,j)
            enddo
          enddo
        enddo
        lensend(n) = mn * 2 * levs
        locsend(n) = (n-1)*lonlenmax*latlenmax*2*levs
        lenrecv(n)=latlen(n)*mylonlen*2*levs
        locrecv(n)=locsend(n)
      enddo
!$omp end parallel do
!
      call mpi_barrier (col_comm,ierr)
!
!     call mpi_alltoallv(works,lensend,locsend,MPI_REAL8_r,
!    &                   workr,lenrecv,locrecv,MPI_REAL8_r,
!    &                   MPI_COMM_WORLD,ierr)

!     call mpi_alltoallv(works,lensend,locsend,MPI_REAL8, &
!                        workr,lenrecv,locrecv,MPI_REAL8, &
!                        col_comm,    ierr)

      call mpi_alltoallv(works,lensend,locsend,MPI_RTYPE, &
                         workr,lenrecv,locrecv,MPI_RTYPE, &
                         col_comm,    ierr)
!
!$omp parallel do private(n,mn,i,j,lat1,lat2,k) &
!$omp schedule(dynamic)
      do n=1,nsizey
        mn=0
        do j=1,latlen(n)
!          jj = latstr(n) + j - 1
!          lat1=global_lats_a(jj)
          lat1=jlist1_sl(j,n)
!          lat1=latg-jlist1_sl(j,n)+1
          lat2=latfull+1-lat1
          do i=1,mylonlen
            mn = mn + 1
            do k=1,levs
              b(lat1,k,i) = workr(1,k,mn,n)
              b(lat2,k,i) = workr(2,k,mn,n)
            enddo
          enddo
        enddo
      enddo
!$omp end parallel do
!
!     print *,' end of ndslfv_we2ns '

      return
      end

! ======================================================================
      subroutine para_ns2we(a,b,levs,latg)
!
! mpi transport from full dimension of west-east to full dimension of
! north-south with latitude shuffl.
!
      use const   , only : RTYPE  ,MPI_RTYPE
      use grid    , only : lonfull,lonhalf,lonpart,lonlenmax,mylonlen, &
                           latfull,lathalf,latpart,latlenmax,mylatlen, &
                           latstr ,latlen ,lonstr ,lonlen
      use rank
      use index
      use mpi
!      use gfs_dyn_mpi_def
!      use gfs_dyn_layout1
      implicit none
!
      integer levs,latg
      real(kind=RTYPE) a(latfull,levs,lonpart)
      real(kind=RTYPE) b(lonfull,levs,latpart)
!      integer global_lats_a(latg)
!
!     real(kind=RTYPE) (kind=kind_mpi_r) works(2*levs*mylonlen*lathalf)
!     real(kind=RTYPE) (kind=kind_mpi_r) workr(2*levs*mylatlen*lonhalf)
      real(kind=RTYPE) works(2,levs,lonlenmax*latlenmax,nsizey)
      real(kind=RTYPE) workr(2,levs,lonlenmax*latlenmax,nsizey)
      integer lensend(nsizey),lenrecv(nsizey)
      integer locsend(nsizey),locrecv(nsizey)
!      integer i,j,k,n,mn,jj,lat1,lat2,ierr
      integer i,j,k,n,mn,lat1,lat2,ierr
!
!     print *,' enter ndslfv_ns2we '

!$omp parallel do private(n,mn,i,j,lat1,lat2,k) &
!$omp schedule(dynamic)
      do n=1,nsizey
        mn=0
        do j=1,latlen(n)
!          jj=latstr(n) + j -1
!          lat1=global_lats_a(jj)
!          lat1=latg-jlist1_sl(j,n)+1
          lat1=jlist1_sl(j,n)
          lat2=latfull+1-lat1
          do i=1,mylonlen
            mn = mn + 1
            do k=1,levs
              works(1,k,mn,n) = a(lat1,k,i)
              works(2,k,mn,n) = a(lat2,k,i)
            enddo
          enddo
        enddo
        lensend(n)=mn*2*levs
        locsend(n)=(n-1)*lonlenmax*latlenmax*2*levs
        lenrecv(n)=mylatlen*lonlen(n)*2*levs
        locrecv(n)=locsend(n)
      enddo
!$omp end parallel do
!
      call mpi_barrier (col_comm,ierr)
!     call mpi_alltoallv(works,lensend,locsend,MPI_REAL8_r,
!    &                   workr,lenrecv,locrecv,MPI_REAL8_r,
!    &                   MPI_COMM_WORLD,ierr)

!     call mpi_alltoallv(works,lensend,locsend,MPI_REAL8, &
!                        workr,lenrecv,locrecv,MPI_REAL8, &
!                        col_comm,    ierr)

      call mpi_alltoallv(works,lensend,locsend,MPI_RTYPE, &
                         workr,lenrecv,locrecv,MPI_RTYPE, &
                         col_comm,    ierr)
!
!$omp parallel do private(n,mn,i,j,k) &
!$omp schedule(dynamic)
      do n=1,nsizey
        mn=0
        do j=1,mylatlen
          do i=lonstr(n),lonstr(n)+lonlen(n)-1
            mn = mn + 1
            do k=1,levs
              b(i        ,k,j) = workr(1,k,mn,n)
              b(i+lonhalf,k,j) = workr(2,k,mn,n)
            enddo
          enddo
        enddo
      enddo
!$omp end parallel do
!
!     print *,' end of ndslfv_ns2we '
!
      return
      end
