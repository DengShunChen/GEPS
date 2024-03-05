      subroutine dmsread(nx,my,lenc,kflag,ifile,z,istat)
!
!  subroutine to read data in pressure level fields
!
! **** input ****
!
!  lrec: field identifaction label
!  lenc: no. of data values in a 2-d lrec field
!  ifile: path/file name containing initial fields
!
! **** output ****
!
!  z: array containing gaussian grid 2-d field
!  istat: status code. .ne. zero means bad read
!
      use mpe
      use rank
      use index
      use const, only:keyi,ihdgi

      implicit  none
      integer   nx,my,lenc,istat

      real      z(nx,my)
!
      logical t_flg
!
      character ifile*255,kflag*1
!
      character crmk*88
!
      write(keyi,1000)ihdgi,kflag,lenc
#ifdef I38K
 1000 format(a28,a1,i9.9)
#else
 1000 format(a26,a1,i7.7)
#endif
!
      t_flg=.false.
!
      if(myrank .eq. 0) then
      call dmsget(ifile,keyi//char(0),z,istat)
      t_flg=.true.
      endif
!
!ch   call mpe_broadcast(z,nx*my,t_flg,mpe_double)
!ch   call mpe_broadcast(istat,1,t_flg,mpe_integer)
      call mpe_bcast(z,nx*my,0,mpe_double)
      call mpe_bcast(istat,1,0,mpe_integer)
!
      if(istat.ne.0) then
!
      write(crmk,100) keyi
  100 format('#######  record ',a38,' missing  ######')
!
      if(myrank .eq. 0) then
      print*, crmk
      print *,'istat=',istat
!CWB2016 force mpi code abort
      stop'dmsread: fatal error !'
      endif
      call mpe_finalize
      call dmsexit(-1)
!
      else
!
#ifdef VERBOSE
      if(myrank .eq. 0) print *,'dms key=',keyi,' found'
#endif
!
      endif
!
      return
      end


      subroutine dmsreadi(nx,my,lenc,kflag,ifile,z,istat)
!
!  subroutine to read data in pressure level fields
!
! **** input ****
!
!  lrec: field identifaction label
!  lenc: no. of data values in a 2-d lrec field
!  ifile: path/file name containing initial fields
!
! **** output ****
!
!  z: array containing gaussian grid 2-d field
!  istat: status code. .ne. zero means bad read
!
      use mpe
      use rank
      use index
      use const, only:keyi,ihdgi

      implicit  none
      integer   nx,my,lenc,istat

      integer z(nx,my)
!
      logical t_flg
!
      character ifile*255,kflag*1
!
      character crmk*88
!
      write(keyi,1000)ihdgi,kflag,lenc
#ifdef I38K
 1000 format(a28,a1,i9.9)
#else
 1000 format(a26,a1,i7.7)
#endif
!
      t_flg=.false.
!
      if(myrank .eq. 0) then
      call dmsget(ifile,keyi//char(0),z,istat)
      t_flg=.true.
      endif
!
!ch   call mpe_broadcast(z,nx*my,t_flg,mpe_integer)
!ch   call mpe_broadcast(istat,1,t_flg,mpe_integer)
      call mpe_bcast(z,nx*my,0,mpe_integer)
      call mpe_bcast(istat,1,0,mpe_integer)
!
      if(istat.ne.0) then
!
      write(crmk,100) keyi
  100 format('#######  record ',a38,' missing  ######')
!
      if(myrank .eq. 0) then
      print*, crmk
      print *,'istat=',istat
!CWB2016 force mpi code abort
      stop'dmsreadi: fatal error !'
      endif
      call mpe_finalize
      call dmsexit(-1)
!
      else
!
#ifdef VERBOSE
      if(myrank .eq. 0) print *,'dms key=',keyi,' found'
#endif
!
      endif
!
      return
      end
!-------------------------------------------------
      subroutine dmsread_split(nx,my,lenc,kflag,ifile,z,istat)
!
!  subroutine to read data in pressure level fields
!
! **** input ****
!
!  lrec: field identifaction label
!  lenc: no. of data values in a 2-d lrec field
!  ifile: path/file name containing initial fields
!
! **** output ****
!
!  z: array containing gaussian grid 2-d field
!  istat: status code. .ne. zero means bad read
!
      use mpe
      use rank
      use index
      use const, only:keyi,ihdgi

      implicit  none
      integer   nx,my,lenc,istat

      real      z(nx,my)
!
      logical t_flg
!
      character ifile*255,kflag*1
!
      character crmk*88
!
      write(keyi,1000)ihdgi,kflag,lenc
#ifdef I38K
 1000 format(a28,a1,i9.9)
#else
 1000 format(a26,a1,i7.7)
#endif
!
      t_flg=.false.
!
      if(col_rank .eq. 0) then
      call dmsget(ifile,keyi//char(0),z,istat)
      t_flg=.true.
      endif
!
!ch   call mpe_broadcast(z,nx*my,t_flg,mpe_double)
!ch   call mpe_broadcast(istat,1,t_flg,mpe_integer)
      call mpe_bcast_col(z,nx*my,0,mpe_double)
      call mpe_bcast_col(istat,1,0,mpe_integer)
!
      if(istat.ne.0) then
!
      write(crmk,100) keyi
  100 format('#######  record ',a38,' missing  ######')
!
      if(col_rank .eq. 0) then
      print*, crmk
      print *,'istat=',istat
!CWB2016 force mpi code abort
      stop'dmsread: fatal error !'
      endif
      call mpe_finalize
      call dmsexit(-1)
!
      else
!
#ifdef VERBOSE
      if(col_rank .eq. 0) print *,'dms key=',keyi,' found'
#endif
!
      endif
!
      return
      end

