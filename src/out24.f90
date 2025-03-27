      subroutine out24 (nx,my,my_max,hf24,qf24,ss24,rs24,asol24,olr24  &
                      ,rain24,rainlp24,dt24,glob,itau,idtg,ggdef,flash24)
!
      use index
      use rank
      use mpe
      use mod_grb2_param , only :ofdir,wrt_grb2_accu_v2
      use const ,only:outdms ,outgrb2 ,ifilout_grb, RTYPE,kflag &
                     ,do_sit,ldailyFCTsst,dailyClm_option,ldailyFCTicesndpt &
                     ,ihdgo,ihdgo2
      USE mo_netcdf,           ONLY:lkvl
      use mod_sitgrid, only:ratioSIT,dtsit24
      USE mod_sit_control,only:loutsit24
      use mod_sst ,only:outtseadiffFCT24

      implicit  none

      integer   nx,my,my_max,itau
      real      dt24

      real      hf24(nxp,my_max),qf24(nxp,my_max),ss24(nxp,my_max),rs24(nxp,my_max), &
                asol24(nxp,my_max),olr24(nxp,my_max),rain24(nxp,my_max)              &
               ,flash24(nxp,my_max),rainlp24(nxp,my_max)

      real(kind=RTYPE) glob(nx,my),wrk(nxp,my_max)
!
      integer*8 idtg
      character*4  ggdef

      integer   imax,jmax,lenc,j,nxj,i,istat,jj
!
      imax=nx
      jmax=my
      lenc= imax*jmax
!
      if( outgrb2 == 1)then
 134                    format( A  ,A ,I10.10 , i4.4       )
           write(ofdir,134 )trim(ifilout_grb),'/',idtg/100 ,itau
           if(myrank==0) call system("mkdir -p "//trim(ofdir) )
      endif
!===
!  Latent heat flux at the surface (W/m**2)
      do jj = 1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
         wrk(i,jj)=qf24(i,jj)/dt24
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('s0043f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0)then
         ihdgo2 = ihdgo
         call wrt_grb2_accu_v2(itau,0,0,10,2,1,0,0,0,24,glob)
      endif
! Sensible heat flux at the surface (W/m**2)
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=hf24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('s0042f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0)then
         ihdgo2 = ihdgo
         call wrt_grb2_accu_v2(itau,0,0,11,2,1,0,0,0,24,glob)
      endif
! Net shortwave (solar) flux at the surface (W/m**2) (positive : downward flux)
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=ss24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('s0031f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0)then
        ihdgo2 = ihdgo
        call wrt_grb2_accu_v2(itau,0,4,9,2,1,0,0,0,24,glob)
      endif
! net surface longwave radiation
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=rs24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('s0032f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0)then
        ihdgo2 = ihdgo
        call wrt_grb2_accu_v2(itau,0,5,5,2,1,0,0,0,24,glob)
      endif
!
!  Total precipitation  24-hours
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=rain24(i,jj)
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('b00626',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0)then
        ihdgo2 = ihdgo
        call wrt_grb2_accu_v2(itau,0,1,8,2,103,0,0,1,24,glob)
      endif
!
!  Total precipitation  24-hours
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=rainlp24(i,jj)
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('b00646',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)

!  The average of latent heat flux release for total precipitation within 24-hours
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=rain24(i,jj)/dt24*2.5e+6
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('b0062f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)

! model top of net solor shortwave radiation
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=asol24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('x0033f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0)then
        ihdgo2 = ihdgo
        call wrt_grb2_accu_v2(itau,0,4,1,2,8,0,0,0,24,glob)
      endif

! model top of Outgoing longwave radiation (OLR)
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=olr24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('x0034f',idtg,itau,ggdef)
      if(outdms.gt.0)  call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0)then
        ihdgo2 = ihdgo
        call wrt_grb2_accu_v2(itau,0,5,5,2,8,0,0,0,24,glob)
      endif
!
!xb110>>
! 24hr average flash density (km-2day-1)
      do jj = 1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
         wrk(i,jj)=flash24(i,jj)/dt24 
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('x00999',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0)then
        ihdgo2 = ihdgo
        call wrt_grb2_accu_v2(itau,0,17,4,6,10,0,0,0,24,glob)
      endif
!xb110<<

!Ocean SIT daily mean output
      if(ldailyFCTsst .OR. ldailyFCTicesndpt .OR. (dailyClm_option.ge.1)) then
        call outtseadiffFCT24(nx,my,my_max,dt24,itau,idtg,ggdef)
      endif

      if(do_sit)then
        if(ldailyFCTsst .OR. ldailyFCTicesndpt .OR. (dailyClm_option.ge.1)) then
          call outtseadiffSIT24(nx,my,my_max,ratioSIT,dtsit24,itau,idtg,ggdef)
        endif
        if(loutsit24)then
          call outsit24(nx,my,my_max,lkvl,itau,idtg,ggdef)
        endif
      endif !do_sit


      return
      end
