#ifdef NDMS

  subroutine dmsexit(istat)
    integer istat
    if (istat == 0 )then
      !call mpe_finalize
      stop 0
    else
      stop -1
    endif
  end subroutine
  !------------------------------
  subroutine dmscls(bckfile,istat)
    character(len=*) bckfile
    integer :: istat
    istat = 0
  end subroutine
  !------------------------------
  subroutine dmsmsg(bckfile,istat)
    character(len=*) bckfile
    integer :: istat
    istat = 0
  end subroutine
  !------------------------------
  subroutine  dmscfg(type_w,argument,istat_w)
    character(len=*) type_w,argument
    integer::istat_w
    istat_w = 0
  end subroutine
  !------------------------------
  subroutine dmschkr (ifilin,keyi,istat)
    use const, only:KLENI
    character(len=*) ifilin,keyi
    character::ifn*1024
    integer :: i,ia,istat
    LOGICAL       :: lex0

    do i = 1, KLENI
    ia = ichar(keyi(i:i))
     if((ia.ge.97).and.(ia.le.122))then
       ia=ia-32
       keyi(i:i)=char(ia)
     endif
    enddo

#ifdef O38K
        ifn=trim(ifilin)//'/'//keyi(17:28)//keyi(7:12)//'/' &
                    //keyi(1:6)//keyi(13:16)//'H'//keyi(30:38)
#else
        ifn=trim(ifilin)//'/'//keyi(15:26)//keyi(7:10)//'/' &
                    //keyi(1:6)//keyi(11:14)//'H'//keyi(28:34)
#endif
    INQUIRE (file=trim(ifn), exist=lex0)
    if ( lex0 ) THEN
      istat = 0
    else
      istat = 1
    endif
  end subroutine
  !------------------------------
  subroutine dmsopn(bckfile,rwstat,istat)
    character(len=*) bckfile
    character::ifn*1024
    character::ipath*1024,idata*1024
    integer :: i,ierr,istat
    logical :: rungetpath ,isdir

    istat = 0
    !origin namelist path
    !inquire(file=trim(bckfile)//'/.', exist=isdir)
    isdir=.false.
    if(bckfile(1:1) =='/'.or.bckfile(1:1)=='.')isdir=.true.
    if(isdir)then 
      istat=0
      return
    endif

    rungetpath = .false.
    if ( index( bckfile ,'@' ) > 0 )then
       rungetpath = .true.
    endif

    !get env variable change into path
    ! exp: bckfile = 'BCKGEPS'
    if( .not. rungetpath ) then
      rungetpath = .false.
      call getenv(bckfile,ipath)
      bckfile=trim(ipath)
      ! ipath format have '@' exp: 'MASOPS@NWPDB'
      if ( index( ipath , '@' ) > 0 )then
        rungetpath = .true.
      endif
    endif

    if (.not. rungetpath )then
      !inquire(file=trim(bckfile)//'/.', exist=isdir)
      isdir=.false.
      if(bckfile(1:1) =='/'.or.bckfile(1:1)=='.')isdir=.true.
      if(isdir)then
       istat = 0
      else
       istat = 1
      endif
      return
    endif

    ! bckfile format have '@' exp: 'MASOPS@NWPDB'
    !get DMS DB path
    call getenv('HOME',ipath)
    ifn=trim(ipath)//'/.dmsrc'
    open(11,file=trim(ifn),action='read',status='old',iostat=ierr )
    if( ierr .ne. 0 )then
      istat = 1
    endif
    
    do while( istat == 0 )
      read(11,'(A512)',end=99)idata
      if( index( idata , '[DMSPATH]' ) > 0 )then
        do i = 1 , 5
         read(11,'(A512)',end=99)idata
         if( index(idata,'/') == 1 ) exit
         if( index(idata,'[') > 0 )then
           istat=1
           return
         endif
        enddo
        exit
      endif
    end do
 
    close(11)
    if( index(idata,'/') == 1 ) then
      ipath = trim(idata)
      i = index( bckfile , '@')
      i = i - 1
      ifn=trim(bckfile(1:i))
      idata=trim( bckfile(i+2:) )
      if (  index( idata , '/' ) > 0 )then
        bckfile =trim(idata)//'.ufs/'//trim(ifn)
      else
        bckfile =trim(ipath)//'/'//trim(idata)//'.ufs/'//trim(ifn)
      endif
      !check path exist
      !exp: /nwpr/gfs/username/.DMSPATH/NWPDB.ufs/MASOPS_eps000
      isdir=.false.
      if(bckfile(1:1) =='/'.or.bckfile(1:1)=='.')isdir=.true.
      !inquire(file=trim(bckfile)//'/.', exist=isdir)
      if(.not.isdir)istat=1
    else
      bckfile=''
      istat=1
    endif

    return

    !dms open error 
    99 continue
    istat = 1
    close(11)
    return

  end subroutine

#endif
