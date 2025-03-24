      subroutine mpe_finalize
!
!     include 'mpif.h'
      use mpi

      use param
      use index
      use rank
      use const
      use grid
      use spec
      use phygrid
      use noah
      use mod_typhoon
      use albn
      use radn
      use ozne_def
      use raddiag
! #ifdef USE_CUDA
!       use mod_ndslfv_monoadv_gpu, only: deallocate_ndslfv_array_gpu
! #endif

      implicit none
      integer  ierr
#ifdef W3TAG
      if (myrank_all==0) call w3tage('CWBGFS')
#endif
      if(io_quilting)then
! method 1
      if(myrank_all .le. Ngfs-1)then

! method 2
!     if(myrank_all .ne. 1)then

! deallocate dynamic arrays
#ifdef USE_CUDA
      ! call deallocate_ndslfv_array_gpu
      call fft_destroy_plan
      call nccl_destroy
#endif
      call deallocate_grid_array
      call deallocate_phygrid_array
      call deallocate_noah_array
      call deallocate_const_array
      call deallocate_spec_array
      call deallocate_typhoon_array
      call deallocate_index_array

!rad
      call deallocate_alb_array
      call deallocate_raddiag_array
#ifdef Readaeroclx
      call deallocate_aerogrid_array
#endif

      endif

      else ! normal mode, i.e. no io_quilting

! deallocate dynamic arrays
#ifdef USE_CUDA
      ! call deallocate_ndslfv_array_gpu
#endif
      call deallocate_grid_array
      call deallocate_phygrid_array
      call deallocate_noah_array
      call deallocate_const_array
      call deallocate_spec_array
      call deallocate_typhoon_array
      call deallocate_index_array

!rad
      call deallocate_alb_array
      call deallocate_raddiag_array

      endif

!     call MPI_BARRIER(MPI_COMM_WORLD, IERR )
#ifndef TIMCOMCPL
      call MPI_FINALIZE(IERR)
#endif
      return
      end
