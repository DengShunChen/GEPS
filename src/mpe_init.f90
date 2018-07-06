      subroutine mpe_init

! CWB2016 io_quilting version

! CWB2017 2dMPI  V1   version
! 1. grid dimension is (isd:ied,lev,my_max)
! 2. nx should be devided by npex, i.e. mod(nx,npex) = 0
! 3. i direction loop example : do i=isd, nxj

! CWB2018/06/12 2dMPI  V2   version
! 1. grid dimension is (nxp,lev,my_max), 'nxp' means nx partial 
! 2. nx is no longer needed to be devided by npex, i.e. mod(nx,npex) >= 0
! 3. i direction loop example : do i=1, nxj

!
      use mpi
!
      use param
      use index
      use rank
      use const
      use grid
      use spec
      use phygrid
      use noah
      use mod_typhoon
!rad------------------------------------------------------------------
      use radn
      use ozne_def
      use albn
      use raddiag
!rad------------------------------------------------------------------

      implicit none

      integer i,ierr,iworld,igfs,iio,mini,m,n

      integer,dimension(:),allocatable :: ranks_gfs,ranks_io

! the whole group, (gfs + io)
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank_all, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nsize_all,  ierr )

      call get_model_param

      if(io_quilting)then

      Ngfs=nsize_all-1
      Nio =1

      allocate (ranks_gfs(Ngfs) ,ranks_io(Nio) ,stat=ierr)
      if (ierr/= 0) then
         write(6,*) 'mpe_init : allocate fail 1 '
         stop
      end if

! gfs     : 0,1,2 .....nsize_all-1
! io      : nsize_all
      do i=1,Ngfs
         ranks_gfs(i)=i-1
      enddo
      do i=1,Nio
         ranks_io(i)=(i-1)+Ngfs
      enddo

      call MPI_COMM_GROUP( MPI_COMM_WORLD, iworld, ierr )
      call MPI_GROUP_excl( iworld, Nio,  ranks_io,  igfs, ierr )
      call MPI_GROUP_excl( iworld, Ngfs ,ranks_gfs, iio, ierr )

! create the sub_group(gfs)
      call MPI_COMM_create( MPI_COMM_WORLD,igfs,MPI_COMM_gfs,ierr )
! create the sub_group(io)
      call MPI_COMM_create( MPI_COMM_WORLD,iio, MPI_COMM_io,ierr )

      if(myrank_all .le. Ngfs-1)then

! gfs group goes here
         call MPI_COMM_RANK( MPI_COMM_gfs, myrank_gfs, ierr )
         call MPI_COMM_SIZE( MPI_COMM_gfs, nsize_gfs,  ierr )
         myrank=myrank_gfs
         nsize=nsize_gfs
!        print *, 'in mpe init gfs group : nsize_all myrank_all nsize myrank= ',nsize_all, myrank_all,nsize, myrank

      else

! io  group goes here
         call MPI_COMM_SIZE( MPI_COMM_io, nsize_io, ierr )
         call MPI_COMM_RANK( MPI_COMM_io, myrank_io,ierr )
         myrank=myrank_io
         nsize=nsize_io
!        print *, 'in mpe init io  group : nsize_all myrank_all nsize myrank= ',nsize_all, myrank_all,nsize, myrank

      endif

      root_gfs=0

      root_io=Ngfs

      call MPI_GROUP_free( iworld, IERR )
      call MPI_GROUP_free( igfs,  IERR )
      call MPI_GROUP_free( iio  ,  IERR )

      deallocate (ranks_gfs,ranks_io)

! from now on, gfs and io groups can
!              communicate with whole group through MPI_COMM_WORLD, or
!              communicate with gfs  group through MPI_COMM_gfs,  or
!              communicate with io    group through MPI_COMM_io       


      npe=nsize
!1d   jtmax=jtrun/npe+1
!1d   my_max=my/npe+1

      ntag=0

      if(myrank_all .le. Ngfs-1)then

!ch>
      if(npex .eq. -1 .or. npey .eq. -1 ) then

! get number of procs in i and j directions
      mini = 2*nsize
      nsizex = 1
      nsizey = nsize
      do m = 1, nsize
        if ( mod( nsize, m ) == 0 ) then
          n = nsize / m
          if ( abs(m-n) < mini  ) then
            mini = abs(m-n)
            nsizex = m
            nsizey = n
          end if
        end if
      end do

      else

! using namlsts setting
      nsizex=npex
      nsizey=npey

      endif

      if ( (nsizex*nsizey) /= npe ) then
         if(myrank == 0)print *,'fatal error : npex*npey  .ne. npe !'
!        call MPI_FINALIZE(IERR)
         stop ! force abort
      endif

         mrow = myrank/nsizex
         ncol = mod(myrank, nsizex)
         call MPI_Comm_split(MPI_COMM_gfs, mrow, ncol, row_comm, ierr)
         call MPI_Comm_split(MPI_COMM_gfs, ncol, mrow, col_comm, ierr)

         call MPI_Comm_rank(row_comm, row_rank, ierr)
         call MPI_Comm_rank(col_comm, col_rank, ierr)

         nxp=nx/nsizex       !nx partial
         n=mod(nx,nsizex)
         if(n.ne.0)nxp=nxp+1

!        if(myrank .eq. 0) then
!         print *,'nsizex nsizey =',nsizex,nsizey
!        endif

!        do i=0,nsize-1
!        if(myrank.eq.i)then
!          print 10,myrank,mrow,ncol,row_rank,col_rank
!        endif
!10      format(5i8)
!        call MPI_Barrier(MPI_COMM_gfs, ierr)
!        enddo

!2d
      my_max=my/nsizey+1
      jtmax=jtrun/nsizey+1

      nxf=nx       !nx full
      myf=my       !my full
      jlen=my_max

      if ( mod( lev, nsizex ) /= 0 ) then
         if(myrank == 0)print *,'fatal error : lev not devided by nsizex !'
!        call MPI_FINALIZE(IERR)
         stop
      else
         levf=lev  !lev full
         Llen=lev/nsizex
         levp=Llen
         Lstart=row_rank*Llen+1
         Lend=Lstart+Llen-1

! for lev*ncld array
         Lstart_ncld=row_rank*((lev*ncld)/nsizex)+1
         Lend_ncld=Lstart_ncld+((lev*ncld)/nsizex)-1
      endif
!ch<

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

      call ioserver(nx*my)

      endif

! non io_quilting
      else

      nsize=nsize_all
      myrank=myrank_all
      MPI_COMM_gfs=MPI_COMM_WORLD

      npe=nsize
!1d   jtmax=jtrun/npe+1
!1d   my_max=my/npe+1

!ch>
      if(npex .eq. -1 .or. npey .eq. -1 ) then

! get number of procs in i and j directions
      mini = 2*nsize
      nsizex = 1
      nsizey = nsize
      do m = 1, nsize
        if ( mod( nsize, m ) == 0 ) then
          n = nsize / m
          if ( abs(m-n) < mini  ) then
            mini = abs(m-n)
            nsizex = m
            nsizey = n
          end if
        end if
      end do

      else

! using namlsts setting
      nsizex=npex
      nsizey=npey

      endif

      if ( (nsizex*nsizey) /= npe ) then
         if(myrank == 0)print *,'fatal error : npex*npey  .ne. npe !'
!        call MPI_FINALIZE(IERR)
         stop ! force abort
      endif

         mrow = myrank/nsizex
         ncol = mod(myrank, nsizex)
         call MPI_Comm_split(MPI_COMM_WORLD, mrow, ncol, row_comm, ierr)
         call MPI_Comm_split(MPI_COMM_WORLD, ncol, mrow, col_comm, ierr)

         call MPI_Comm_rank(row_comm, row_rank, ierr)
         call MPI_Comm_rank(col_comm, col_rank, ierr)

         nxp=nx/nsizex       !nx partial
         n=mod(nx,nsizex)
         if(n.ne.0)nxp=nxp+1

!        if(myrank .eq. 0) then
!         print *,'nsizex nsizey =',nsizex,nsizey
!        endif

         do i=0,nsize-1
!        if(myrank.eq.i)then
!          print 10,myrank,mrow,ncol,row_rank,col_rank
!        endif
!        call MPI_Barrier(MPI_COMM_WORLD, ierr)
         enddo

!2d
      my_max=my/nsizey+1
      jtmax=jtrun/nsizey+1

      nxf=nx       !nx full
      myf=my       !my full
      jlen=my_max

      if ( mod( lev, nsizex ) /= 0 ) then
         if(myrank == 0)print *,'fatal error : lev not devided by nsizex !'
!        call MPI_FINALIZE(IERR)
         stop
      else
         levf=lev  !lev full
         Llen=lev/nsizex
         levp=Llen
         Lstart=row_rank*Llen+1
         Lend=Lstart+Llen-1

! for lev*ncld array
         Lstart_ncld=row_rank*((lev*ncld)/nsizex)+1
         Lend_ncld=Lstart_ncld+((lev*ncld)/nsizex)-1
      endif
!ch<

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

      endif

      return
      end
