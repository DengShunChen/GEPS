#ifdef TIMCOMCPL
  subroutine mpe_init(mpi_comm_mct)
#else
  subroutine mpe_init
#endif
! CWB2016 io_quilting version

! CWB2017 2dMPI  V1   version
! 1. grid dimension is (isd:ied,lev,my_max)
! 2. nx should be devided by npex, i.e. mod(nx,npex) = 0
! 3. i direction loop example : do i=isd, nxj

! CWB2018/06/12 2dMPI  V2   version
! 1. grid dimension is (nxp,lev,my_max), 'nxp' means nx partial
! 2. nx is no longer needed to be devided by npex, i.e. mod(nx,npex) >= 0
! 3. i direction loop example : do i=1, nxj

     use mpi

     use param
     use index
     use rank
     use const, only: allocate_const_array
     use grid, only: allocate_grid_array
     use spec, only: allocate_spec_array
     use phygrid, only: allocate_phygrid_array
     use noah, only: allocate_noah_array
     use mod_typhoon, only: allocate_typhoon_array
#ifdef USE_CUDA
     use spec_cuda_graph, only: allocate_spec_cg_buffer
#endif
     !>>>rad------------------------------------------------------------------
     use radn
     use ozne_def
     use albn, only: allocate_alb_array
     use raddiag, only: allocate_raddiag_array
     !<<<rad------------------------------------------------------------------

     implicit none

#ifdef TIMCOMCPL
    integer, intent(in), optional :: mpi_comm_mct
#endif
    integer i,ierr,istat,iworld,igfs,iio,mini,m,n

     integer, dimension(:), allocatable :: ranks_gfs, ranks_io

#if defined(RSM) && defined(CWB_MPMD)
     !for MPMD mode
     call mpmd_init(nsize_all, myrank_all, MPI_COMM_gfs_all, root_rsm, istat)
     if (istat .ne. 0) stop'mpmd_init fail !'
#else
    ! the whole group, (gfs + io)
      #ifdef TIMCOMCPL
      MPI_COMM_atm = mpi_comm_mct
      #else
      MPI_COMM_atm = MPI_COMM_WORLD
      call MPI_INIT( ierr )
      #endif
    call MPI_COMM_RANK( MPI_COMM_atm, myrank_all, ierr )
    call MPI_COMM_SIZE( MPI_COMM_atm, nsize_all,  ierr )
    root_rsm = nsize_all
#endif
#ifdef USE_CUDA
     call device_init(myrank_all, nsize_all)
#endif
#ifdef W3TAG
     if (myrank_all == 0) call w3tagb('TCoGFS', 2021, 1721, 067, 'GFS')
#endif
     call get_model_param

     if (io_quilting) then
        Ngfs = nsize_all - 1
        Nio = 1

        allocate (ranks_gfs(Ngfs), ranks_io(Nio), stat=ierr)
        if (ierr /= 0) then
           write (6, *) 'mpe_init : allocate fail 1 '
           stop
        end if

        ! gfs     : 0,1,2 .....nsize_all-1
        ! io      : nsize_all
        do i = 1, Ngfs
           ranks_gfs(i) = i - 1
        end do
        do i = 1, Nio
           ranks_io(i) = (i - 1) + Ngfs
        end do

#if defined(RSM) && defined(CWB_MPMD)
        call MPI_COMM_GROUP(MPI_COMM_gfs_all, iworld, ierr)
#else
        call MPI_COMM_GROUP(MPI_COMM_atm, iworld, ierr)
#endif
        call MPI_GROUP_excl(iworld, Nio, ranks_io, igfs, ierr)
        call MPI_GROUP_excl(iworld, Ngfs, ranks_gfs, iio, ierr)

#if defined(RSM) && defined(CWB_MPMD)
        ! create the sub_group(gfs)
        call MPI_COMM_create(MPI_COMM_gfs_all, igfs, MPI_COMM_gfs, ierr)
        ! create the sub_group(io)
        call MPI_COMM_create(MPI_COMM_gfs_all, iio, MPI_COMM_io, ierr)
#else
        ! create the sub_group(gfs)
        call MPI_COMM_create(MPI_COMM_atm, igfs, MPI_COMM_gfs, ierr)
        ! create the sub_group(io)
        call MPI_COMM_create(MPI_COMM_atm, iio, MPI_COMM_io, ierr)
#endif

        if (myrank_all .le. Ngfs - 1) then
           ! gfs group goes here
           call MPI_COMM_RANK(MPI_COMM_gfs, myrank_gfs, ierr)
           call MPI_COMM_SIZE(MPI_COMM_gfs, nsize_gfs, ierr)
           myrank = myrank_gfs
           nsize = nsize_gfs
        else
           ! io  group goes here
           call MPI_COMM_SIZE(MPI_COMM_io, nsize_io, ierr)
           call MPI_COMM_RANK(MPI_COMM_io, myrank_io, ierr)
           myrank = myrank_io
           nsize = nsize_io
        end if

        root_gfs = 0
        root_io = Ngfs

        call MPI_GROUP_free(iworld, IERR)
        call MPI_GROUP_free(igfs, IERR)
        call MPI_GROUP_free(iio, IERR)

        deallocate (ranks_gfs, ranks_io)

        ! from now on, gfs and io groups can
        !              communicate with whole group through MPI_COMM_WORLD, or
        !              communicate with gfs  group through MPI_COMM_gfs,  or
        !              communicate with io    group through MPI_COMM_io
        npe = nsize
        ntag = 0

        if (myrank_all .le. Ngfs - 1) then
           if (npex .eq. -1 .or. npey .eq. -1) then
              ! get number of procs in i and j directions
              mini = 2*nsize
              nsizex = 1
              nsizey = nsize
              do m = 1, nsize
                 if (mod(nsize, m) == 0) then
                    n = nsize/m
                    if (abs(m - n) < mini) then
                       mini = abs(m - n)
                       nsizex = m
                       nsizey = n
                    end if
                 end if
              end do
           else
              ! using namlsts setting
              nsizex = npex
              nsizey = npey
           end if
           if (npe .lt. lev) then
              if (myrank == 0) print *, 'fatal error : npe  .lt. lev !'
!        call MPI_FINALIZE(IERR)
              stop ! force abort
           end if

           if ((nsizex*nsizey) /= npe) then
              if (myrank == 0) print *, 'fatal error : npex*npey  .ne. npe !'
!       call MPI_FINALIZE(IERR)
              stop ! force abort
           end if
           mrow = myrank/nsizex
           ncol = mod(myrank, nsizex)
           call MPI_Comm_split(MPI_COMM_gfs, mrow, ncol, row_comm, ierr)
           call MPI_Comm_split(MPI_COMM_gfs, ncol, mrow, col_comm, ierr)
           call MPI_Comm_rank(row_comm, row_rank, ierr)
           call MPI_Comm_rank(col_comm, col_rank, ierr)

           nxp = nx/nsizex       !nx partial
           n = mod(nx, nsizex)
           if (n .ne. 0) nxp = nxp + 1

           my_max = my/nsizey + 1
           jtmax = jtrun/nsizey + 1

           nxf = nx       !nx full
           myf = my       !my full
           jlen = my_max

           if (mod(lev, nsizex) /= 0) then
              if (myrank == 0) print *, 'fatal error : lev not devided by nsizex !'
!          call MPI_FINALIZE(IERR)
              stop
           else
              levf = lev  !lev full
              Llen = lev/nsizex
              levp = Llen
              Lstart = row_rank*Llen + 1
              Lend = Lstart + Llen - 1

              ! for lev*ncld array
              Lstart_ncld = row_rank*((lev*ncld)/nsizex) + 1
              Lend_ncld = Lstart_ncld + ((lev*ncld)/nsizex) - 1
           end if

           ! allocate dynamic arrays
           call allocate_grid_array
           call allocate_phygrid_array
           call allocate_noah_array
           call allocate_const_array
           call allocate_spec_array
           call allocate_typhoon_array
           call allocate_index_array
           !rad
           call allocate_alb_array
           call allocate_raddiag_array
           ! initial block data
           call init_block
        else
!       call ioserver(nx*my)
           call ioserver_grb2(nx, my)
        end if
     else ! non io_quilting
        nsize = nsize_all
        myrank = myrank_all
#if defined(RSM) && defined(CWB_MPMD)
        MPI_COMM_gfs = MPI_COMM_gfs_all
#else
        MPI_COMM_gfs = MPI_COMM_atm
#endif

        npe = nsize
        if (npex .eq. -1 .or. npey .eq. -1) then
           ! get number of procs in i and j directions
           mini = 2*nsize
           nsizex = 1
           nsizey = nsize
           do m = 1, nsize
              if (mod(nsize, m) == 0) then
                 n = nsize/m
                 if (abs(m - n) < mini) then
                    mini = abs(m - n)
                    nsizex = m
                    nsizey = n
                 end if
              end if
           end do
        else
           ! using namlsts setting
           nsizex = npex
           nsizey = npey
        end if

        if ((nsizex*nsizey) /= npe) then
           if (myrank == 0) print *, 'fatal error : npex*npey  .ne. npe !'
!        call MPI_FINALIZE(IERR)
           stop ! force abort
        end if

        mrow = myrank/nsizex
        ncol = mod(myrank, nsizex)
        call MPI_Comm_split(MPI_COMM_gfs, mrow, ncol, row_comm, ierr)
        call MPI_Comm_split(MPI_COMM_gfs, ncol, mrow, col_comm, ierr)

        call MPI_Comm_rank(row_comm, row_rank, ierr)
        call MPI_Comm_rank(col_comm, col_rank, ierr)

        nxp = nx/nsizex       !nx partial
        n = mod(nx, nsizex)
        if (n .ne. 0) nxp = nxp + 1

        my_max = my/nsizey + 1
        jtmax = jtrun/nsizey + 1

        nxf = nx       !nx full
        myf = my       !my full
        jlen = my_max

        if (mod(lev, nsizex) /= 0) then
           if (myrank == 0) print *, 'fatal error : lev not devided by nsizex !'
!        call MPI_FINALIZE(IERR)
           stop
        else
           levf = lev  !lev full
           Llen = lev/nsizex
           levp = Llen
           Lstart = row_rank*Llen + 1
           Lend = Lstart + Llen - 1

           ! for lev*ncld array
           Lstart_ncld = row_rank*((lev*ncld)/nsizex) + 1
           Lend_ncld = Lstart_ncld + ((lev*ncld)/nsizex) - 1
        end if

        ! allocate dynamic arrays
        call allocate_grid_array
        call allocate_phygrid_array
        call allocate_noah_array
        call allocate_const_array
        call allocate_spec_array
        call allocate_typhoon_array
        call allocate_index_array
        !rad
        call allocate_alb_array
        call allocate_raddiag_array
        ! initial block data
        call init_block
#ifdef USE_CUDA
        call nccl_init
        call allocate_spec_cg_buffer
#endif
     end if

  end subroutine mpe_init
