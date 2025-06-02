      subroutine out2d (nx,lev,my,my_max,itau,idtg,taudir,ntau         &
       ,hflux,qflux,tg,gwet,snr,z0,raintot,raincu,rainlp               &
       ,plcl,cumtop,ss,rs,alb,gwclim,glob,acld                         &
       ,ugws,vgws,t2,q2,rh2,rh10,u10,v10,gfx,rld,sld,wk_xy             &
!       ,soil_xy,canopy,ggdef,lwrite,flash)
       ,soil_xy,canopy,ggdef)
!
      use rank
      use mpe
      use index
      use mod_outflds
      use mod_grb2_param , only :ofdir,wrt_grb2_v2,wrt_grb2_accu_v2
      use const ,only:outgrb2 ,outdms ,domfc ,RTYPE,kflag
!
      implicit  none

      integer   nx,lev,my,my_max,itau,ntau 

      real      hflux(nxp,my_max),qflux(nxp,my_max),                       &
                tg(nxp,my_max),gwet(nxp,my_max),                           &
                snr(nxp,my_max),z0(nxp,my_max),                            &
                raintot(nxp,my_max),raincu(nxp,my_max),rainlp(nxp,my_max), &
                plcl(nxp,my_max),cumtop(nxp,my_max),                       &
                ss(nxp,my_max),rs(nxp,my_max),                             &
                alb(nxp,my_max),gwclim(nxp,my_max),acld(lev,my),           &
                ugws(nxp,my_max),vgws(nxp,my_max),t2(nxp,my_max),          &
                q2(nxp,my_max),rh2(nxp,my_max),rh10(nxp,my_max),           &
                u10(nxp,my_max),v10(nxp,my_max),gfx(nxp,my_max),           &
                rld(nxp,my_max),sld(nxp,my_max)

      character*16 taudir(ntau)
      character*4 ggdef
      character::varkey*6
      integer*8 idtg
!
      real(kind=RTYPE) glob(nx,my),mout(nx,my)
!
!  local array
!
      real(kind=RTYPE) globp(nxp,my_max)
!
      real(kind=RTYPE) wk_xy(nxp,my_max,12),soil_xy(nxp,my_max,12)
!soil
      real      canopy(nxp,my_max)
!
      character*6  label(ntau),labx
      character*1  lflg
!
      integer   num,n,levz,lenc,lenc2,i,ia,kk,j,jj,nxj,istat,nc,ilev
      real      tnshun
      integer   kfdb,jfdb
!kc >
      real      sfac2,sfac3,sfac4
      integer::praint
      integer:: ptp0(9),ptp1(9)
!kc <
!      logical lwrite
!xb110>
!      real flash(nxp,my_max)  !flash density (in flashes km^-2 day^-1)
!xb110<

      num= 0
      nc = 0
      if(myrank .eq. 0) print *,' out2d ntau=',ntau
      do 20 n=1,ntau
      read(taudir(n),'(a6,1x,i4)') labx,levz
      if(levz.eq.0) then
      num= num+1
      label(num)= labx
      if(myrank .eq. 0) print *,'num= ',num,' out2d labx=',labx
      endif
   20 continue

      if(num.eq.0) return

!
      tnshun= 1.0
      lenc= nx*my
      lenc2= lev*my
!
!  change the letter from uppercase to lowercase
!
      do n = 1, num
      labx = label(n)  
      do i = 1, 6
       ia=ichar(labx(i:i))
       if((ia.ge.65).and.(ia.le.90))then
         ia=ia+32
         labx(i:i)=char(ia)
       endif
      enddo
      label(n) = labx
      enddo
!
      do 30 kk=1,num
!
!surface albedo  0 - 1.0
      if(label(kk).eq.'s00030') then
      globp=alb
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00030',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,19,1,3,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!surface roughness(m)
      if(label(kk).eq.'s00040') then
      globp=z0
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00040',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,1,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif

!   ---------- Precipitation ----------
!
      if(label(kk).eq.'b00620')then
      !precipitaion accum. interval info for grib2
      praint=mod(itau,12)
      if(praint==0)praint=12
!
!cumulus parameterization precipitation
      globp=raincu
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b00630',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,10,2,103,0,0,1,praint/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
 
      globp=rainlp
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b00640',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,47,2,103,0,0,1,praint/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!total precipitation  every 12 hour reset to zero (mm)
      call syslbl_w ('b00620',idtg,itau,ggdef)
      do 98 jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 98 i=1,nxj
        globp(i,jj)=raincu(i,jj)+rainlp(i,jj)
 98   continue
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,8,2,103,0,0,1,praint/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!accu. total precipitation from tau 0
      if( itau==0 .or.  itau .gt. nint(domfc) )then
!
      globp=raintot
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b0062t',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,49,0,103,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      endif !domfc
      go to 30
      endif
!snow depth of water state at the surface (mm)
      if(label(kk).eq.'b00650') then
      globp=snr
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b00650',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,60,2,103,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!atmosphere column precipitable water (mm)
      if(label(kk).eq.'x00590') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,4),glob)
      call syslbl_w ('x00590',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,3,2,10,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!   ---------- short wave Radiation ----------
!net shortwave solar flux at the surface
      if(label(kk).eq.'s00310') then
      globp=ss
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00310',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,4,9,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!net radiation flux at the surface
      if(label(kk).eq.'s00300') then
      do jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
        globp(i,jj)=ss(i,jj)-rs(i,jj)
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00300',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,4,0,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!
!donwward shortwave radiation at the surface
      if(label(kk).eq.'s003u0') then
      globp=sld
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s003u0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,4,7,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!solar absorption at the top atmosphere
      if(label(kk).eq.'x00330') then
      globp=plcl
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('x00330',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,4,1,2,8,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!   ---------- long wave Radiation ----------
!net longwave infrared flux at the surface
      if(label(kk).eq.'s00320') then
      globp=rs
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00320',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,5,5,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!downward longwave radiation at the surface
      if(label(kk).eq.'s003x0') then
      globp=rld
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s003x0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,5,3,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!outgoing longwave radiation (OLR) 
!net longwave flux at the top atmosphere
      if(label(kk).eq.'x00340') then
      globp=cumtop
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('x00340',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,5,5,2,8,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!   ---------- Heat flux ----------
!sensible heat flux at the surface
      if(label(kk).eq.'s00420') then
      globp=hflux
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00420',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,0,11,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!lentent heat flux at the surface
      if(label(kk).eq.'s00430') then
      globp=qflux
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00430',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,0,10,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!ground heat flux
      if(label(kk).eq.'s00440') then
      globp=gfx
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00440',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,10,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!u conponent of surface grivity wave stress
      if(label(kk).eq.'s00450') then
      globp=ugws
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00450',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,3,16,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!v conponent of surface grivity wave stress
      if(label(kk).eq.'s00460') then
      globp=vgws
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('s00460',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,3,17,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!   ---------- temperature ----------
!land suface tempaerature or sea surface temperature
      if(label(kk).eq.'s00100') then
      globp=tg
!      call unify_reduceintp_idw(nx,my,my_max,globp,glob)
      call unify_reduceintp_idw(nx,my,my_max,globp,glob)
      call syslbl_w ('s00100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,0,17,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!land skin air temperature (the model lowest)
      if(label(kk).eq.'b00100') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,1),glob)
      call syslbl_w ('b00100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,0,0,2,103,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 2m temerature
      if(label(kk).eq.'b02100') then
      if( itau==0 .or. itau .gt. nint(domfc) )then
      globp=t2
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b02100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,0,0,2,103,0,2,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      endif !domfc
      go to 30
      endif
!   ---------- wind component ----------
!skin u component ( model lowest)
      if(label(kk).eq.'b00200') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,2),glob)
      call syslbl_w ('b00200',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,2,2,2,103,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!skinvu component ( model lowest)
      if(label(kk).eq.'b00210') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,3),glob)
      call syslbl_w ('b00210',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,2,3,2,103,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!10m u component
      if(label(kk).eq.'b10200') then
      if( itau==0 .or. itau .gt. nint(domfc) )then
      globp=u10
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b10200',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,2,2,2,103,0,10,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      endif !domfc
      go to 30
      endif
!10m v component
      if(label(kk).eq.'b10210') then
      if( itau==0 .or.  itau .gt. nint(domfc) )then
      globp=v10
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b10210',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,2,3,2,103,0,10,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      endif !domfc
      go to 30
      endif
!   ---------- humidity ----------
!skin relative humidity (model lowest)
      if(label(kk).eq.'b00510') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,5),glob)
      call syslbl_w ('b00510',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,1,2,103,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!2m specific humidity
      if(label(kk).eq.'b02500') then
      if( itau==0 .or. itau .gt. nint(domfc) )then
      globp=q2
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b02500',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,0,6,103,0,2,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      endif !domfc
      go to 30
      endif
!2m relative humidity
      if(label(kk).eq.'b02510') then
      if( itau==0 .or. itau .gt. nint(domfc) )then
      do 37 jj=1,jlistnum
        j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 37 i=1,nxj
        globp(i,jj)=rh2(i,jj) * 100.0
 37   continue
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b02510',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,1,2,103,0,2,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      endif !domfc
      go to 30
      endif
!10m relative humidity
      if(label(kk).eq.'b10510') then
      do 38 jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 38 i=1,nxj
        globp(i,jj)=rh10(i,jj) * 100.0
 38   continue
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('b10510',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,1,2,103,0,10,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!
!   ---------- Soil and Land model ----------
!ground wetness 0 - 1.0
      if(label(kk).eq.'s005a0') then
      do 36 jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 36 i=1,nxj
        globp(i,jj)=gwet(i,jj)/20.
 36   continue
      call unify_reduceintp_idw(nx,my,my_max,globp,glob)
      call syslbl_w ('s005a0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,3,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!soil moisture content (mm)
      if(label(kk).eq.'s005a1') then
      globp=gwet
      call unify_reduceintp_idw(nx,my,my_max,globp,glob)
      call syslbl_w ('s005a1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,3,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!canopy moisture content(mm)
      if(label(kk).eq.'s005c0') then
      globp=canopy
      call unify_reduceintp_idw(nx,my,my_max,globp,glob)
      call syslbl_w ('s005c0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,19,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif

! s**5b0
!
!   s01??? = sa1??? = L01???   0-10cm
!   s02??? =                  10-200cm
!   s03??? = sa2??? = L02???  10-40cm
!   s04??? = sa3??? = L03???  40-100cm
!   s05??? = sa4??? = L04???  100-200cm
!
!
! 0-10cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'sa15b0') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,1),glob)
      call syslbl_w ('sa15b0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,2,151,0,1,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif

! 10-40cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'sa25b0') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,2),glob)
      call syslbl_w ('sa25b0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,2,151,0,2,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 40-100cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'sa35b0') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,3),glob)
      call syslbl_w ('sa35b0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,2,151,0,3,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 100-200cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'sa45b0') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,4),glob)
      call syslbl_w ('sa45b0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,2,151,0,4,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!kc >
! 0-10cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s015b0') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,1),glob)
      call syslbl_w ('s015b0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,2,151,0,1,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif

      sfac2=3./19.
      sfac3=6./19.
      sfac4=10./19.
!  10-200cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s025b0') then
      do jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
       globp(i,jj)=soil_xy(i,jj,2)*sfac2+soil_xy(i,jj,3)*sfac3 &
                  +soil_xy(i,jj,4)*sfac4
      end do
      end do
      call unify_reduceintp_idw(nx,my,my_max,globp,glob)
      call syslbl_w ('s025b0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,2,151,0,2,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!xb13 <
! 10-40cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s035b0') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,2),glob)
      call syslbl_w ('s035b0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,2,151,0,3,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif

! 40-100cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s045b0') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,3),glob)
      call syslbl_w ('s045b0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,2,151,0,4,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 100-200cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s055b0') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,4),glob)
      call syslbl_w ('s055b0',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,9,2,151,0,5,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!xb13 >
! s**5b1
! 0-10cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'sa15b1') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,5),glob)
      call syslbl_w ('sa15b1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,10,2,151,0,1,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 10-40cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'sa25b1') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,6),glob)
      call syslbl_w ('sa25b1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,10,2,151,0,2,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 40-100cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'sa35b1') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,7),glob)
      call syslbl_w ('sa35b1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,10,2,151,0,3,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 100-200cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'sa45b1') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,8),glob)
      call syslbl_w ('sa45b1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,10,2,151,0,4,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!xb13
! 0-10cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s015b1') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,5),glob)
      call syslbl_w ('s015b1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,10,2,151,0,1,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif

! 10-200cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s025b1') then
      do jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
       globp(i,jj)=soil_xy(i,jj,6)*sfac2+soil_xy(i,jj,7)*sfac3 &
                  +soil_xy(i,jj,8)*sfac4
      end do
      end do
      call unify_reduceintp_idw(nx,my,my_max,globp,glob)
      call syslbl_w ('s025b1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,10,2,151,0,2,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 10-40cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s035b1') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,6),glob)
      call syslbl_w ('s035b1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,10,2,151,0,3,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 40-100cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s045b1') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,7),glob)
      call syslbl_w ('s045b1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,10,2,151,0,4,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 100-200cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s055b1') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,8),glob)
      call syslbl_w ('s055b1',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,10,2,151,0,5,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!<xb13
! s**100
! 0-10cm Volumetric soil temperature(K)
      if(label(kk).eq.'sa1100') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,9),glob)
      call syslbl_w ('sa1100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,18,2,151,0,1,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 10-40cm Volumetric soil temperature(K)
      if(label(kk).eq.'sa2100') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,10),glob)
      call syslbl_w ('sa2100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,18,2,151,0,2,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 40-100cm Volumetric soil temperature(K)
      if(label(kk).eq.'sa3100') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,11),glob)
      call syslbl_w ('sa3100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,18,2,151,0,3,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 100-200cm Volumetric soil temperature(K)
      if(label(kk).eq.'sa4100') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,12),glob)
      call syslbl_w ('sa4100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,18,2,151,0,4,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!xb13>
! 0-10cm Volumetric soil temperature(K)
      if(label(kk).eq.'s01100') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,9),glob)
      call syslbl_w ('s01100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,18,2,151,0,1,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 10-200cm Volumetric soil temperature(K)
      if(label(kk).eq.'s02100') then

      do jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
       globp(i,jj)=soil_xy(i,jj,10)*sfac2+soil_xy(i,jj,11)*sfac3 &
                  +soil_xy(i,jj,12)*sfac4
      end do
      end do
      call unify_reduceintp_idw(nx,my,my_max,globp,glob)
      call syslbl_w ('s02100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,18,2,151,0,2,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 10-40cm Volumetric soil temperature(K)
      if(label(kk).eq.'s03100') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,10),glob)
      call syslbl_w ('s03100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,18,2,151,0,3,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 40-100cm Volumetric soil temperature(K)
      if(label(kk).eq.'s04100') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,11),glob)
      call syslbl_w ('s04100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,18,2,151,0,4,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! 100-200cm Volumetric soil temperature(K)
      if(label(kk).eq.'s05100') then
      call unify_reduceintp_idw(nx,my,my_max,soil_xy(1,1,12),glob)
      call syslbl_w ('s05100',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,3,18,2,151,0,5,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!   ---------- cloud cover ----------
!< xb13
! ctot_total cloud fraction
      if(label(kk).eq.'x00770') then
      if( itau==0 .or. itau .gt. nint(domfc) )then
      do  jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do  i=1,nxj
        globp(i,jj)=wk_xy(i,jj,6) * 100.0
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('x00770',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,6,1,2,10,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      endif !domfc
      go to 30
      endif
! chig_high cloud fraction
      if(label(kk).eq.'x00760') then
      do  jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do  i=1,nxj
        globp(i,jj)=wk_xy(i,jj,7) * 100.0
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('x00760',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,6,5,2,10,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! cmid_middle cloud fraction
      if(label(kk).eq.'x00750') then
      do  jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do  i=1,nxj
        globp(i,jj)=wk_xy(i,jj,8) * 100.0
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('x00750',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,6,4,2,10,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
! clow_low cloud fraction
      if(label(kk).eq.'x00740') then
      do  jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do  i=1,nxj
        globp(i,jj)=wk_xy(i,jj,9) * 100.0
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl_w ('x00740',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,6,3,2,10,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!zonal mean cloudiness of Y-Z cross section 0-1 
      if(label(kk).eq.'x00730') then
      call mpe_unify(acld,lev,my,2,mpe_double)
      !use 'mout' could effect another output
      do jfdb=1,my
      do kfdb=1,lev
       glob(kfdb,jfdb)=acld(kfdb,jfdb)*100.
      end do
      end do
      call syslbl_w ('x00730',idtg,itau,ggdef)
      if(outdms.gt.0) call dmswrit(lev,my,lenc2,kflag,glob,istat)
      !if(outgrb2==1.and.myrank==0) call wrt_grb2(itau,0,6,22,2,10,0,0.,glob)
      go to 30
      endif

!   ---------- other ----------
! pbl height
      if(label(kk).eq.'pbl000') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,10),glob)
      call syslbl_w ('pbl000',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,3,18,1,10,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!
! specific humidity 
      if(label(kk).eq.'m01500') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,11),glob)
      call syslbl_w ('m01500',idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,0,9,104,0,1,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif
!
      if(label(kk).eq.'m60500') then
      if ( lev .lt. 100 ) then
        ilev=lev
        lflg='m'
      else
        ilev=mod(lev,100)
        lflg='n'
      endif
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,12),glob)
      write(varkey,'(A1,I2.2,A3)')lflg,ilev,'500'
      call syslbl_w (varkey,idtg,itau,ggdef)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,0,9,104,0,lev,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      go to 30
      endif

!move to  out24.f90
!xb110> flash density
!      if(label(kk).eq.'fshden') then
!      call unify_reduceintp(nx,my,my_max,flash,glob)
!      call syslbl ('x00999',idtg,itau,ggdef,ihdg)
!      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
!      call split(nx,my,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
!      go to 30
!      endif
!xb110<

   30 continue  !============end do ( kk=1,num )

!
      if(outdms.gt.0)then
       if ( myrank .lt. nc )             &
          call dmswrit_split(nx,my,lenc,kflag,mout,istat)
      endif
      if(outgrb2 == 1 )then
       if ( myrank .lt. nc ) then
        if(ptp1(8)==-999)then
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),mout)
        else
        call wrt_grb2_accu_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),ptp1(8),ptp1(9),mout)
        endif
       endif
      endif
!
      return
      end
!---------------------------------------------------------------
      subroutine split(nx,my,lenc,nc,glob,mout)
!
      use rank
      use const, only: RTYPE,kflag,ihdgo,ihdgo2
!
      implicit none
!
      integer   nx,my,nc,istat,lenc,lev
      real(kind=RTYPE) glob(nx,my),mout(nx,my)

      if (myrank .eq. nc) then
        mout  = glob
        ihdgo2 = ihdgo
      endif
        nc    = nc + 1
!
      if ( nc .eq. nsize .and. myrank .lt. nc ) then
        call dmswrit_split(nx,my,lenc,kflag,mout,istat)
        nc    = 0
      endif
!
      return
      end
!---------------------------------------------------------------
      subroutine split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!
      use rank
      use const, only: RTYPE,kflag,ihdgo,ihdgo2,outgrb2 ,outdms
      use mod_grb2_param , only :wrt_grb2_v2 , wrt_grb2_accu_v2
!
      implicit none
!
      integer   nx,my,nc,istat,lenc,lev
      real(kind=RTYPE) glob(nx,my),mout(nx,my)
!      character*26 ihdg,ihdg2
!      character*80 ifilout
      integer::ptp0(9),ptp1(9)
      integer::itau_spl

      if (myrank .eq. nc) then
        mout   = glob
        ihdgo2 = ihdgo
        ptp1   = ptp0
      endif
        nc    = nc + 1
!
      if ( nc .eq. nsize .and. myrank .lt. nc ) then
        if(outdms.gt.0)then
          call dmswrit_split(nx,my,lenc,kflag,mout,istat)
        endif
        if(outgrb2 == 1 )then
#ifdef O38K
         read(ihdgo2(7:12),'(I6)')itau_spl
#else
         read(ihdgo2(7:10),'(I4)')itau_spl
#endif
         if(ptp1(8)==-999)then
          call wrt_grb2_v2(itau_spl,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
              ,ptp1(6),ptp1(7),mout)
         else
          call wrt_grb2_accu_v2(itau_spl,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
              ,ptp1(6),ptp1(7),ptp1(8),ptp1(9),mout)
         endif
        endif
        nc    = 0
      endif
!
      return
      end subroutine split2
