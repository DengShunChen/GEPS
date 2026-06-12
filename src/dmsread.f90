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
      use const, only:keyi,ihdgi,KLENI

      implicit  none
      integer   nx,my,lenc,istat

      real      z(nx,my)
!
      logical t_flg
!
      character kflag*1, ifn*512
!
      character crmk*88
      character(len=*)ifile
      integer :: i,ia
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
#ifdef NDMS
      do i = 1, KLENI
      ia = ichar(keyi(i:i))
       if((ia.ge.97).and.(ia.le.122))then
         ia=ia-32
         keyi(i:i)=char(ia)
       endif
      enddo
  #ifdef O38K
      ifn=trim(ifile)//'/'//keyi(17:28)//keyi(7:12)//'/' &
                  //keyi(1:6)//keyi(13:16)//'H'//keyi(30:38)
  #else
      ifn=trim(ifile)//'/'//keyi(15:26)//keyi(7:10)//'/' &
                  //keyi(1:6)//keyi(11:14)//'H'//keyi(28:34)
  #endif
      open(12,file=trim(ifn),access='stream',form='unformatted' &
             ,action='read',convert='little_endian',iostat=istat)
      read(12)z(:,:)
      close(12)
#else
      call dmsget(ifile,keyi//char(0),z,istat)
#endif
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
      use const, only:keyi,ihdgi,KLENI

      implicit  none
      integer   nx,my,lenc,istat

      integer z(nx,my)
!
      logical t_flg
!
      character kflag*1, ifn*512
!
      character crmk*88

      character(len=*)ifile
      integer :: i,ia
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
#ifdef NDMS
      do i = 1, KLENI
      ia = ichar(keyi(i:i))
       if((ia.ge.97).and.(ia.le.122))then
         ia=ia-32
         keyi(i:i)=char(ia)
       endif
      enddo
  #ifdef O38K
      ifn=trim(ifile)//'/'//keyi(17:28)//keyi(7:12)//'/' &
                  //keyi(1:6)//keyi(13:16)//'I'//keyi(30:38)
  #else
      ifn=trim(ifile)//'/'//keyi(15:26)//keyi(7:10)//'/' &
                  //keyi(1:6)//keyi(11:14)//'I'//keyi(28:34)
  #endif
      open(12,file=trim(ifn),access='stream',form='unformatted' &
             ,action='read',convert='little_endian',iostat=istat)
      read(12)z(:,:)
      close(12)
#else
      call dmsget(ifile,keyi//char(0),z,istat)
#endif
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
      use const, only:keyi,ihdgi,KLENI

      implicit  none
      integer   nx,my,lenc,istat

      real      z(nx,my)
!
      logical t_flg
!
      character kflag*1, ifn*512
!
      character crmk*88
      character(len=*)ifile
      integer :: i,ia
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
#ifdef NDMS
      do i = 1, KLENI
      ia = ichar(keyi(i:i))
       if((ia.ge.97).and.(ia.le.122))then
         ia=ia-32
         keyi(i:i)=char(ia)
       endif
      enddo
  #ifdef O38K
      ifn=trim(ifile)//'/'//keyi(17:28)//keyi(7:12)//'/' &
                  //keyi(1:6)//keyi(13:16)//'H'//keyi(30:38)
  #else
      ifn=trim(ifile)//'/'//keyi(15:26)//keyi(7:10)//'/' &
                  //keyi(1:6)//keyi(11:14)//'H'//keyi(28:34)
  #endif
      open(12,file=trim(ifn),access='stream',form='unformatted' &
             ,action='read',convert='little_endian',iostat=istat)
      read(12)z(:,:)
      close(12)
#else
      call dmsget(ifile,keyi//char(0),z,istat)
#endif
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

