      subroutine out2d (nx,lev,my,my_max,ifilout,itau,idtg,taudir,ntau &
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
!hcwei  write grb2
      use mod_grb2_param
      use const ,only:out_pres_form,domfc

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
      integer*8 idtg
!
      real      glob(nx,my),mout(nx,my)
!
!  local array
!
!byl      real      glob1(nx,my),glob2(nx,my),glob3(nx,my)
      real      globp(nxp,my_max)
!
      real      wk_xy(nxp,my_max,12)     
!soil
      real      soil_xy(nxp,my_max,12),canopy(nxp,my_max)
!
      character*80 ifilout
      character*26 ihdg,ihdg2
      character*6 label(ntau),labx
!
      integer   num,n,levz,lenc,lenc2,i,ia,kk,j,jj,nxj,istat,nc
      real      tnshun
!kc >
      real      sfac2,sfac3,sfac4
      integer::praint
!kc <
!      logical lwrite
!xb110>
!      real flash(nxp,my_max)  !flash density (in flashes km^-2 day^-1)
!xb110<
!xb119

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
      do 30 kk=1,num   ! ---  start  ---


!surface albedo  0 - 1.0
      if(label(kk).eq.'s00030') then
      call unify_reduceintp(nx,my,my_max,alb,glob)
!byl      call mpe2d_unify(glob,alb)
      call syslbl ('s00030',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,19,1,3,1,0,0.,glob)
      endif
      go to 30
      endif
!surface roughness(m)
      if(label(kk).eq.'s00040') then
      call unify_reduceintp(nx,my,my_max,z0,glob)
!byl      call mpe2d_unify(glob,z0)
      call syslbl ('s00040',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,1,2,1,0,0.,glob)
      endif
      go to 30
      endif

!   ---------- Precipitation ----------
!
      if(label(kk).eq.'b00620')then
!
!byl      call mpe2d_unify(glob1,raincu)
!byl      call mpe2d_unify(glob2,rainlp)

!byl      glob=glob1
!cumulus parameterization precipitation
      call unify_reduceintp(nx,my,my_max,raincu,glob)
      call syslbl ('b00630',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,10,2,1,0,0.,glob)
      endif
 
!byl      glob=glob2
!large scale precipitation
      call unify_reduceintp(nx,my,my_max,rainlp,glob)
      call syslbl ('b00640',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,9,2,1,0,0.,glob)
      endif
!total precipitation  every 12 hour reset to zero (mm)
      call syslbl ('b00620',idtg,itau,ggdef,ihdg)
!byl      do 98 j=1,my
      do 98 jj=1,jlistnum
!byl       nxj=nxdef(j)
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 98 i=1,nxj
!byl        glob(i,j)=glob1(i,j)+glob2(i,j)
        globp(i,jj)=raincu(i,jj)+rainlp(i,jj)
 98   continue
      call unify_reduceintp(nx,my,my_max,globp,glob)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        !call wrt_grb2(itau,0,1,7,2,1,0,0.,glob)
        praint=mod(itau,12)
        if(praint==0)praint=12
        call wrt_grb2_accu(itau,0,1,7,2,1,0,0.,1,praint,glob)
      endif
!accu. total precipitation from tau 0
      if( itau==0 .or.  itau .ge. nint(domfc) )then
      call unify_reduceintp(nx,my,my_max,raintot,glob)
!byl      call mpe2d_unify(glob,raintot)
      call syslbl ('b0062t',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,8,0,1,0,0.,glob)
      endif
      endif !domfc
      go to 30
      endif
!snow depth of water state at the surface (mm)
      if(label(kk).eq.'b00650') then
      call unify_reduceintp(nx,my,my_max,snr,glob)
!byl      call mpe2d_unify(glob,snr)
      call syslbl ('b00650',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,60,2,103,0,0.,glob)
      endif
      go to 30
      endif
!atmosphere column precipitable water (mm)
      if(label(kk).eq.'x00590') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,4),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,4))
      call syslbl ('x00590',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,3,2,7,0,0.,glob)
      endif
      go to 30
      endif
!   ---------- short wave Radiation ----------
!net shortwave solar flux at the surface
      if(label(kk).eq.'s00310') then
      call unify_reduceintp(nx,my,my_max,ss,glob)
!byl      call mpe2d_unify(glob,ss)
      call syslbl ('s00310',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,4,0,2,1,0,0.,glob)
      endif
      go to 30
      endif
!net radiation flux at the surface
      if(label(kk).eq.'s00300') then
!byl      call mpe2d_unify(glob1,ss)
!byl      call mpe2d_unify(glob2,rs)
!byl      do j=1,my
      do jj=1,jlistnum
!byl       nxj=nxdef(j)
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
!byl        glob(i,j)=glob1(i,j)-glob2(i,j)
        globp(i,jj)=ss(i,jj)-rs(i,jj)
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00300',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,4,9,2,1,0,0.,glob)
      endif
      go to 30
      endif
!donwward shortwave radiation at the surface
      if(label(kk).eq.'s003u0') then
      call unify_reduceintp(nx,my,my_max,sld,glob)
!byl      call mpe2d_unify(glob,sld)
      call syslbl ('s003u0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,4,7,2,1,0,0.,glob)
      endif
      go to 30
      endif
!solar absorption at the top atmosphere
      if(label(kk).eq.'x00330') then
      call unify_reduceintp(nx,my,my_max,plcl,glob)
!byl      call mpe2d_unify(glob,plcl)
      call syslbl ('x00330',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,4,1,2,1,0,0.,glob)
      endif
      go to 30
      endif
!   ---------- long wave Radiation ----------
!net longwave infrared flux at the surface
      if(label(kk).eq.'s00320') then
      call unify_reduceintp(nx,my,my_max,rs,glob)
!byl      call mpe2d_unify(glob,rs)
      call syslbl ('s00320',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,5,0,2,1,0,0.,glob)
      endif
      go to 30
      endif
!downward longwave radiation at the surface
      if(label(kk).eq.'s003x0') then
      call unify_reduceintp(nx,my,my_max,rld,glob)
!byl      call mpe2d_unify(glob,rld)
      call syslbl ('s003x0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,5,3,2,1,0,0.,glob)
      endif
      go to 30
      endif
!outgoing longwave radiation (OLR) 
!net longwave flux at the top atmosphere
      if(label(kk).eq.'x00340') then
      call unify_reduceintp(nx,my,my_max,cumtop,glob)
!byl      call mpe2d_unify(glob,cumtop)
      call syslbl ('x00340',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,5,1,2,1,0,0.,glob)
      endif
      go to 30
      endif
!   ---------- Heat flux ----------
!sensible heat flux at the surface
      if(label(kk).eq.'s00420') then
      call unify_reduceintp(nx,my,my_max,hflux,glob)
!byl      call mpe2d_unify(glob,hflux)
      call syslbl ('s00420',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,0,11,2,1,0,0.,glob)
      endif
      go to 30
      endif
!lentent heat flux at the surface
      if(label(kk).eq.'s00430') then
      call unify_reduceintp(nx,my,my_max,qflux,glob)
!byl      call mpe2d_unify(glob,qflux)
      call syslbl ('s00430',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,0,10,2,1,0,0.,glob)
      endif
      go to 30
      endif
!ground heat flux
      if(label(kk).eq.'s00440') then
      call unify_reduceintp(nx,my,my_max,gfx,glob)
!byl      call mpe2d_unify(glob,gfx)
      call syslbl ('s00440',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,10,2,1,0,0.,glob)
      endif
      go to 30
      endif
!u conponent of surface grivity wave stress
      if(label(kk).eq.'s00450') then
      call unify_reduceintp(nx,my,my_max,ugws,glob)
!byl      call mpe2d_unify(glob,ugws)
      call syslbl ('s00450',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,3,16,2,1,0,0.,glob)
      endif
      go to 30
      endif
!v conponent of surface grivity wave stress
      if(label(kk).eq.'s00460') then
      call unify_reduceintp(nx,my,my_max,vgws,glob)
!byl      call mpe2d_unify(glob,vgws)
      call syslbl ('s00460',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,3,17,2,1,0,0.,glob)
      endif
      go to 30
      endif
!   ---------- temperature ----------
!land suface tempaerature or sea surface temperature
      if(label(kk).eq.'s00100') then
      call unify_reduceintp(nx,my,my_max,tg,glob)
!byl      call mpe2d_unify(glob,tg)
      call syslbl ('s00100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,0,0,2,1,0,0.,glob)
      endif
      go to 30
      endif
!land skin air temperature (the model lowest)
      if(label(kk).eq.'b00100') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,1),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,1))
     call syslbl ('b00100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,0,0,2,103,0,0.,glob)
      endif
      go to 30
      endif
! 2m temerature
      if(label(kk).eq.'b02100') then
      if( itau==0 .or. itau .ge. nint(domfc) )then
      call unify_reduceintp(nx,my,my_max,t2,glob)
!byl      call mpe2d_unify(glob,t2)
      call syslbl ('b02100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,0,0,2,103,0,2.,glob)
      endif
      go to 30
      endif !domfc
      endif
!   ---------- wind component ----------
!skin u component ( model lowest)
      if(label(kk).eq.'b00200') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,2),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,2))
      call syslbl ('b00200',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,2,2,103,0,0.,glob)
      endif
      go to 30
      endif
!skinvu component ( model lowest)
      if(label(kk).eq.'b00210') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,3),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,3))
      call syslbl ('b00210',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,3,2,103,0,0.,glob)
      endif
      go to 30
      endif
!10m u component
      if(label(kk).eq.'b10200') then
      if( itau==0 .or. itau .ge. nint(domfc) )then
      call unify_reduceintp(nx,my,my_max,u10,glob)
!byl      call mpe2d_unify(glob,u10)
      call syslbl ('b10200',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,2,2,103,0,10.,glob)
      endif
      go to 30
      endif !domfc
      endif
!10m v component
      if(label(kk).eq.'b10210') then
      if( itau==0 .or.  itau .ge. nint(domfc) )then
      call unify_reduceintp(nx,my,my_max,v10,glob)
!byl      call mpe2d_unify(glob,v10)
      call syslbl ('b10210',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,2,3,2,103,0,10.,glob)
      endif
      go to 30
      endif !domfc
      endif
!   ---------- humidity ----------
!skin relative humidity (model lowest)
      if(label(kk).eq.'b00510') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,5),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,5))
      call syslbl ('b00510',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,1,2,103,0,0.,glob)
      endif
      go to 30
      endif
!2m specific humidity
      if(label(kk).eq.'b02500') then
      if( itau==0 .or. itau .ge. nint(domfc) )then
      call unify_reduceintp(nx,my,my_max,q2,glob)
!byl      call mpe2d_unify(glob,q2)
      call syslbl ('b02500',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,0,6,103,0,2.,glob)
      endif
      go to 30
      endif !domfc
      endif
!2m relative humidity
      if(label(kk).eq.'b02510') then
      if( itau==0 .or. itau .ge. nint(domfc) )then
      do 37 jj=1,jlistnum
        j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 37 i=1,nxj
        globp(i,jj)=rh2(i,jj) * 100.0
 37   continue
      call unify_reduceintp(nx,my,my_max,globp,glob)
!hcw      call unify_reduceintp(nx,my,my_max,rh2,glob)
!byl      call mpe2d_unify(glob,rh2)
      call syslbl ('b02510',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,1,2,103,0,2.,glob)
      endif
      go to 30
      endif !domfc
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
!      call mpe2d_unify(glob,rh10)
      call syslbl ('b10510',idtg,itau,ggdef,ihdg)
!      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,1,2,103,0,10.,glob)
      endif
      go to 30
      endif
!
!      if(label(kk).eq.'b10510') then
!      call mpe2d_unify(glob,wk_xy(1,1,5))
!      call syslbl ('b10510',idtg,itau,ggdef,ihdg)
!      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
!      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
!      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
!      go to 30
!      endif

!   ---------- Soil and Land model ----------
!ground wetness 0 - 1.0
      if(label(kk).eq.'s005a0') then
!byl      call mpe2d_unify(glob,gwet)
!byl      do 36 j=1,my
      do 36 jj=1,jlistnum
!byl       nxj=nxdef(j)
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 36 i=1,nxj
!byl        glob(i,j)=glob(i,j)/20.
        globp(i,jj)=gwet(i,jj)/20.
 36   continue
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s005a0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,3,1,0,0.,glob)
      endif
      go to 30
      endif
!soil moisture content (mm)
      if(label(kk).eq.'s005a1') then
      call unify_reduceintp(nx,my,my_max,gwet,glob)
!byl      call mpe2d_unify(glob,gwet)
      call syslbl ('s005a1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,25,2,1,0,0.,glob)
      endif
      go to 30
      endif
!canopy moisture content(mm)
      if(label(kk).eq.'s005c0') then
      call unify_reduceintp(nx,my,my_max,canopy,glob)
!byl      call mpe2d_unify(glob,canopy)
      call syslbl ('s005c0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,19,2,1,0,0.,glob)
      endif
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
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,1),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,1))
      call syslbl ('sa15b0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,2,151,0,1.,glob)
      endif
      go to 30
      endif

! 10-40cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'sa25b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,2),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,2))
      call syslbl ('sa25b0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,2,151,0,2.,glob)
      endif
      go to 30
      endif
! 40-100cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'sa35b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,3),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,3))
      call syslbl ('sa35b0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,2,151,0,3.,glob)
      endif
      go to 30
      endif
! 100-200cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'sa45b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,4),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,4))
      call syslbl ('sa45b0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,2,151,0,4.,glob)
      endif
      go to 30
      endif
!kc >
! 0-10cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s015b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,1),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,1))
      call syslbl ('s015b0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,2,151,0,1.,glob)
      endif
      go to 30
      endif

      sfac2=3./19.
      sfac3=6./19.
      sfac4=10./19.
!  10-200cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s025b0') then
!byl      call mpe2d_unify(glob1,soil_xy(1,1,2))
!byl      call mpe2d_unify(glob2,soil_xy(1,1,3))
!byl      call mpe2d_unify(glob3,soil_xy(1,1,4))
!byl      do j=1,my
      do jj=1,jlistnum
!byl       nxj=nxdef(j)
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
!byl       glob(i,j)=glob1(i,j)*sfac2+glob2(i,j)*sfac3 &
!byl                  +glob3(i,j)*sfac4
       globp(i,jj)=soil_xy(i,jj,2)*sfac2+soil_xy(i,jj,3)*sfac3 &
                  +soil_xy(i,jj,4)*sfac4
      end do
      end do
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s025b0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,2,151,0,2.,glob)
      endif
      go to 30
      endif
!xb13 <
! 10-40cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s035b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,2),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,2))
      call syslbl ('s035b0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,2,151,0,3.,glob)
      endif
      go to 30
      endif

! 40-100cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s045b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,3),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,3))
      call syslbl ('s045b0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,2,151,0,4.,glob)
      endif
      go to 30
      endif
! 100-200cm Volumetric soil moisture fraction (0-1.0)
      if(label(kk).eq.'s055b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,4),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,4))
      call syslbl ('s055b0',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,0,9,2,151,0,5.,glob)
      endif
      go to 30
      endif
!xb13 >
! s**5b1
! 0-10cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'sa15b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,5),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,5))
      call syslbl ('sa15b1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,5,2,151,0,1.,glob)
      endif
      go to 30
      endif
! 10-40cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'sa25b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,6),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,6))
      call syslbl ('sa25b1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,5,2,151,0,2.,glob)
      endif
      go to 30
      endif
! 40-100cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'sa35b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,7),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,7))
      call syslbl ('sa35b1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,5,2,151,0,3.,glob)
      endif
      go to 30
      endif
! 100-200cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'sa45b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,8),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,8))
      call syslbl ('sa45b1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,5,2,151,0,4.,glob)
      endif
      go to 30
      endif
!xb13
! 0-10cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s015b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,5),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,5))
      call syslbl ('s015b1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,5,2,151,0,1.,glob)
      endif
      go to 30
      endif

! 10-200cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s025b1') then
!byl      call mpe2d_unify(glob1,soil_xy(1,1,6))
!byl      call mpe2d_unify(glob2,soil_xy(1,1,7))
!byl      call mpe2d_unify(glob3,soil_xy(1,1,8))
!byl      do j=1,my
      do jj=1,jlistnum
!byl       nxj=nxdef(j)
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
!byl       glob(i,j)=glob1(i,j)*sfac2+glob2(i,j)*sfac3 &
!byl                  +glob3(i,j)*sfac4
       globp(i,jj)=soil_xy(i,jj,6)*sfac2+soil_xy(i,jj,7)*sfac3 &
                  +soil_xy(i,jj,8)*sfac4
      end do
      end do
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s025b1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,5,2,151,0,2.,glob)
      endif
      go to 30
      endif
! 10-40cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s035b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,6),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,6))
      call syslbl ('s035b1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,5,2,151,0,3.,glob)
      endif
      go to 30
      endif
! 40-100cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s045b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,7),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,7))
      call syslbl ('s045b1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,5,2,151,0,4.,glob)
      endif
      go to 30
      endif
! 100-200cm Unfrozen(liquid) soil moisture content(volumetric fraction)
      if(label(kk).eq.'s055b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,8),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,8))
      call syslbl ('s055b1',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,5,2,151,0,5.,glob)
      endif
      go to 30
      endif
!<xb13
! s**100
! 0-10cm Volumetric soil temperature(K)
      if(label(kk).eq.'sa1100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,9),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,9))
      call syslbl ('sa1100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,18,2,151,0,1.,glob)
      endif
      go to 30
      endif
! 10-40cm Volumetric soil temperature(K)
      if(label(kk).eq.'sa2100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,10),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,10))
      call syslbl ('sa2100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,18,2,151,0,2.,glob)
      endif
      go to 30
      endif
! 40-100cm Volumetric soil temperature(K)
      if(label(kk).eq.'sa3100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,11),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,11))
      call syslbl ('sa3100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,18,2,151,0,3.,glob)
      endif
      go to 30
      endif
! 100-200cm Volumetric soil temperature(K)
      if(label(kk).eq.'sa4100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,12),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,12))
      call syslbl ('sa4100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,18,2,151,0,4.,glob)
      endif
      go to 30
      endif
!xb13>
! 0-10cm Volumetric soil temperature(K)
      if(label(kk).eq.'s01100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,9),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,9))
      call syslbl ('s01100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,18,2,151,0,1.,glob)
      endif
      go to 30
      endif
! 10-200cm Volumetric soil temperature(K)
      if(label(kk).eq.'s02100') then

!byl      call mpe2d_unify(glob1,soil_xy(1,1,10))
!byl      call mpe2d_unify(glob2,soil_xy(1,1,11))
!byl      call mpe2d_unify(glob3,soil_xy(1,1,12))
!byl      do j=1,my
      do jj=1,jlistnum
!byl       nxj=nxdef(j)
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
!byl       glob(i,j)=glob1(i,j)*sfac2+glob2(i,j)*sfac3 &
!byl                  +glob3(i,j)*sfac4
       globp(i,jj)=soil_xy(i,jj,10)*sfac2+soil_xy(i,jj,11)*sfac3 &
                  +soil_xy(i,jj,12)*sfac4
      end do
      end do
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s02100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,18,2,151,0,2.,glob)
      endif
      go to 30
      endif
! 10-40cm Volumetric soil temperature(K)
      if(label(kk).eq.'s03100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,10),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,10))
      call syslbl ('s03100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,18,2,151,0,3.,glob)
      endif
      go to 30
      endif
! 40-100cm Volumetric soil temperature(K)
      if(label(kk).eq.'s04100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,11),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,11))
      call syslbl ('s04100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,18,2,151,0,4.,glob)
      endif
      go to 30
      endif
! 100-200cm Volumetric soil temperature(K)
      if(label(kk).eq.'s05100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,12),glob)
!byl      call mpe2d_unify(glob,soil_xy(1,1,12))
      call syslbl ('s05100',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,2,3,18,2,151,0,5.,glob)
      endif
      go to 30
      endif
!   ---------- cloud cover ----------
!< xb13
! ctot_total cloud fraction
      if(label(kk).eq.'x00770') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,6),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,6))
      call syslbl ('x00770',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,6,1,3,7,0,0.,glob)
      endif
      go to 30
      endif
! chig_high cloud fraction
      if(label(kk).eq.'x00760') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,7),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,7))
      call syslbl ('x00760',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,6,5,3,7,0,0.,glob)
      endif
      go to 30
      endif
! cmid_middle cloud fraction
      if(label(kk).eq.'x00750') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,8),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,8))
      call syslbl ('x00750',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,6,4,3,7,0,0.,glob)
      endif
      go to 30
      endif
! clow_low cloud fraction
      if(label(kk).eq.'x00740') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,9),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,9))
      call syslbl ('x00740',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,6,3,3,7,0,0.,glob)
      endif
      go to 30
      endif
!zonal mean cloudiness of Y-Z cross section 0-1 
      if(label(kk).eq.'x00730') then
      call mpe_unify(acld,lev,my,2,mpe_double)
      do i=1,lev*my
       acld(i,1)=acld(i,1)*100.
      end do
      call syslbl ('x00730',idtg,itau,ggdef,ihdg)
      if(out_pres_form==1)then
        call dmswrit(lev,my,ihdg,lenc2,'H',ifilout,acld,istat)
      !elseif(out_pres_form==2.and.myrank==0)then
      !  call wrt_grb2(itau,0,6,22,2,7,0,0.,glob)
      endif
      go to 30
      endif

!   ---------- other ----------
! pbl height
      if(label(kk).eq.'pbl000') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,10),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,10))
      call syslbl ('pbl000',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,3,18,1,7,0,0.,glob)
      endif
      go to 30
      endif
!
! specific humidity 
      if(label(kk).eq.'m01500') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,11),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,11))
      call syslbl ('m01500',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,0,9,104,0,1.,glob)
      endif
      go to 30
      endif
!
      if(label(kk).eq.'m60500') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,12),glob)
!byl      call mpe2d_unify(glob,wk_xy(1,1,12))
      call syslbl ('m60500',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      if(out_pres_form==1)then
        call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      elseif(out_pres_form==2.and.myrank==0)then
        call wrt_grb2(itau,0,1,0,9,104,0,72.,glob)
      endif
      go to 30
      endif

!move to  out24.f90
!xb110> flash density
!      if(label(kk).eq.'fshden') then
!      call unify_reduceintp(nx,my,my_max,flash,glob)
!!byl      call mpe2d_unify(glob,flash)
!      call syslbl ('x00999',idtg,itau,ggdef,ihdg)
!!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!!byl      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
!      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
!      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
!      go to 30
!      endif
!xb110<

   30 continue  !============end do ( kk=1,num )

!
      if(out_pres_form==1)then
      if ( myrank .lt. nc )             &
         call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,mout,istat)
      endif
!
      return
      end
!---------------------------------------------------------------
      subroutine split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
!
      use rank
!
      implicit none
!
      integer   nx,my,nc,istat,lenc,lev
      real      glob(nx,my),mout(nx,my)
      character*26 ihdg,ihdg2
      character*80 ifilout

      if (myrank .eq. nc) then
        mout  = glob
        ihdg2 = ihdg
      endif
        nc    = nc + 1
!
      if ( nc .eq. lev .and. myrank .lt. nc ) then
        call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,mout,istat)
        nc    = 0
      endif
!
      return
      end
