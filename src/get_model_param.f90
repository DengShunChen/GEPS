
      subroutine get_model_param
!
      use param

      use const, only : ifilin,cwbout,bckfile,namlsts,   &
                        ifilout,crdate,ocards,phyout,cntrl, &
                        ifilin_ncep,ifilin_sst,ifilin_nc,&
                        ifilin_ClmANA,ifilin_ClmFCT,ifilout_grb, &
                        ifilin_aero

      use paramt

      use rank
      use mod_grb2_param, only:grbmem,grbnumm
!
      implicit  none

!
      namelist /model_param/nco,lev,ncld,octahedral,nout &
                           ,io_quilting,npex,npey !2dMPI
!
      integer istat
      character*80 filist
      character*255 pathname,logicname,truefile
      namelist /filst/ ifilin,cwbout,bckfile,namlsts &
                     , ifilout,crdate,ocards,phyout,cntrl &
                     , ifilin_ncep,ifilin_sst,ifilin_nc &
                     , ifilin_ClmANA,ifilin_ClmFCT,ifilout_grb &
                     , ifilin_aero

      namelist /grb_conf/ grbmem,grbnumm

      data pathname/'NWPETCGLB'/
      data logicname/'filist'/
!
!  read namlist of path/file name(operation)
!
      call getfname(pathname,logicname,truefile,istat)
      if(istat.ne.0)then
        print *,'getfname : error','RANK=',myrank
        call mpe_finalize
        call dmsexit(-1)
!      else
!        if(myrank .eq. 0) print *, 'check truefile=',truefile
      endif
!
      open (unit=12,file=trim(truefile),form='formatted')
!
      read (12,filst,end=110)
!
  110 continue
      close(12)
!
!  read namelist for model parameters
!
      open (unit=1,file=trim(namlsts),form='formatted')
!
      read (1,model_param,end=120)
!
  120 continue
      read (1,grb_conf,end=121)
  121 continue
      close(1)


!     if(myrank .eq. 0) print model_param
      if(myrank_all .eq. 0) print model_param

      if ( octahedral ) then
        my=2*nco
        nx=16+4*nco
        jtrun= 2*((1+(4*nco-1)/4)/2)

        im=16+4*nco
        jm=2*nco
        jtr= 2*((1+(4*nco-1)/4)/2)
      else
        nx=4*nco
        my=nx/2
        jtrun= 2*((1+(nx-1)/2)/2)
        
        im=nx
        jm=im/2
        jtr= 2*((1+(im-1)/2)/2)
      endif
      mlmax= jtrun*(jtrun+1)/2

      lm=lev
      mlm=jtr*(jtr+1)/2

      imax=nx
      jmax=imax/2
      lpx=26
      kvkw=im/2

      ilm=im*lm
      im2=kvkw*(lm+1)

      jmhalf=jm/2
      lmX2=lm*2
      lmX4=lm*4
      lmX10=lm*10
 
      return
      end
