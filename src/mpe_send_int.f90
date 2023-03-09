      subroutine mpe_send_int( itau ,n , itag , ist)
      use rank, only : root_io
      use mpi
      integer::n
      integer :: itau(n)  ,itag

      call MPI_SEND( itau, n, MPI_integer , root_io, itag , &
                     MPI_COMM_WORLD, ist )
      return
      end

      subroutine mpe_recv_int( itau ,n , itag , ist)
      use rank, only : root_gfs
      use mpi
      integer::n
      integer :: itau(n)  ,itag
      integer:: ierr ,ISTATUS(MPI_STATUS_SIZE)

      call MPI_RECV( itau , n, MPI_INTEGER , root_gfs, &
                     itag, MPI_COMM_WORLD, ISTATUS, ierr )
      return
      end
