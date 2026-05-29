      subroutine out24_gpu (nx,my,my_max,hf24,qf24,ss24,rs24,asol24,olr24  &
                      ,rain24,rainlp24,dt24,glob,itau,idtg,ggdef,flash24)
!
      use index
      use rank
      use mpe
      use mod_grb2_param , only :ofdir,wrtgrb2_accu_v2_gpu
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

      real(kind=RTYPE) wrk(nxp,my_max), glob(nx,my), mout(nx,my)
!
      integer*8 idtg
      character*4  ggdef

      integer   imax,jmax,lenc,j,nxj,i,istat,jj
      integer:: ptp0(9),ptp1(9) ,nc
      integer, parameter:: async_id = 1

!$acc wait(async_id)
!$acc enter data create( wrk,glob,mout ) async(async_id)
!
      imax=nx
      jmax=my
      lenc= imax*jmax

      nc = 0
!
      if( outgrb2 == 1)then
 134                    format( A  ,A ,I10.10 , i4.4       )
           write(ofdir,134 )trim(ifilout_grb),'/',idtg/100 ,itau
           if(myrank==0) call system("mkdir -p "//trim(ofdir) )
      endif
!===
!  Latent heat flux at the surface (W/m**2)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
          wrk(i,jj)=qf24(i,jj)/dt24
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('s0043f',idtg,itau,ggdef)
      ptp0=(/0,0,10,2,1,0,0,0,24/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
! Sensible heat flux at the surface (W/m**2)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj=1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
           wrk(i,jj)=hf24(i,jj)/dt24
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('s0042f',idtg,itau,ggdef)
      ptp0=(/0,0,11,2,1,0,0,0,24/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
! Net shortwave (solar) flux at the surface (W/m**2) (positive : downward flux)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj=1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
           wrk(i,jj)=ss24(i,jj)/dt24
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('s0031f',idtg,itau,ggdef)
      ptp0=(/0,4,9,2,1,0,0,0,24/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
! net surface longwave radiation
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj=1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
           wrk(i,jj)=rs24(i,jj)/dt24
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('s0032f',idtg,itau,ggdef)
      ptp0=(/0,5,5,2,1,0,0,0,24/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!
!  Total precipitation  24-hours
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj=1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
           wrk(i,jj)=rain24(i,jj)
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('b00626',idtg,itau,ggdef)
      ptp0=(/0,1,8,2,103,0,0,1,24/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!
!  Total precipitation  24-hours
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj=1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
           wrk(i,jj)=rainlp24(i,jj)
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('b00646',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)

!  The average of latent heat flux release for total precipitation within 24-hours
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj=1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
           wrk(i,jj)=rain24(i,jj)/dt24*2.5e+6
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('b0062f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)

! model top of net solor shortwave radiation
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj=1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
          wrk(i,jj)=asol24(i,jj)/dt24
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('x0033f',idtg,itau,ggdef)
      ptp0=(/0,4,1,2,8,0,0,0,24/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

! model top of Outgoing longwave radiation (OLR)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj=1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
          wrk(i,jj)=olr24(i,jj)/dt24
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('x0034f',idtg,itau,ggdef)
      ptp0=(/0,5,5,2,8,0,0,0,24/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!
!xb110>>
! 24hr average flash density (km-2day-1)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1,jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
          wrk(i,jj)=flash24(i,jj)/dt24 
         endif
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call syslbl_w ('x00999',idtg,itau,ggdef)
      ptp0=(/0,17,4,6,10,0,0,0,24/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!xb110<<

      if(outdms.gt.0)then
      if ( myrank .lt. nc ) then
         !$acc update self(mout) async(async_id)
         !$acc wait(async_id)
         call dmswrit_split(nx,my,lenc,kflag,mout,istat)
      endif
      endif

      if(outgrb2 == 1 )then
       if ( myrank .lt. nc ) then
        call wrtgrb2_accu_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),ptp1(8),ptp1(9),mout)
       endif
      endif
!=====================

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

!$acc wait(async_id)
!$acc exit data delete( wrk,glob,mout ) async(async_id)
!$acc wait(async_id)

      return
      end
