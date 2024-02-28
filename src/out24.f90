      subroutine out24 (nx,my,my_max,hf24,qf24,ss24,rs24,asol24,olr24  &
                      ,rain24,rainlp24,dt24,glob,itau,idtg,ggdef,flash24)
!
      use index
      use mpe
      use mod_grb2_param  !for write grib2 data
      use const ,only:outdms ,outgrb2 ,ifilout_grb, RTYPE,kflag

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
      call syslbl ('s0043f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) call wrt_grb2_accu(itau,0,0,10,2,1,0,0.,0,24,glob)
! Sensible heat flux at the surface (W/m**2)
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=hf24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('s0042f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) call wrt_grb2_accu(itau,0,0,11,2,1,0,0.,0,24,glob)
! Net shortwave (solar) flux at the surface (W/m**2) (positive : downward flux)
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=ss24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('s0031f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) call wrt_grb2_accu(itau,0,4,9,2,1,0,0.,0,24,glob)
! net surface longwave radiation
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=rs24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('s0032f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) call wrt_grb2_accu(itau,0,5,5,2,1,0,0.,0,24,glob)
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
      call syslbl ('b00626',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) call wrt_grb2_accu(itau,0,1,8,2,103,0,0.,1,24,glob)
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
      call syslbl ('b00646',idtg,itau,ggdef)
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
      call syslbl ('b0062f',idtg,itau,ggdef)
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
      call syslbl ('x0033f',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) call wrt_grb2_accu(itau,0,4,1,2,8,0,0.,0,24,glob)

! model top of Outgoing longwave radiation (OLR)
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=olr24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('x0034f',idtg,itau,ggdef)
      if(outdms.gt.0)  call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0)call wrt_grb2_accu(itau,0,5,5,2,8,0,0.,0,24,glob)
!
!xb110>>
      do jj = 1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
         wrk(i,jj)=flash24(i,jj)/dt24 !xb110, 24hr average flash density (km-2day-1)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('x00999',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) call wrt_grb2_accu(itau,0,17,4,6,7,0,0.,0,24,glob)
!xb110<<


      return
      end
