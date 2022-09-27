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
      use const, only: RTYPE,kflag
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
      character*80 ifilout
      character*26 ihdg,ihdg2
      character*6 label(ntau),labx
!
      integer   num,n,levz,lenc,lenc2,i,ia,kk,j,jj,nxj,istat,nc
      real      tnshun
!kc >
      real      sfac2,sfac3,sfac4
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
      if(label(kk).eq.'s00430') then
      globp=qflux
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00430',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00420') then
      globp=hflux
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00420',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00100') then
      globp=tg
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00030') then
      globp=alb
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00030',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s005a0') then
      do 36 jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 36 i=1,nxj
        globp(i,jj)=gwet(i,jj)/20.
 36   continue
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s005a0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s005a1') then
      globp=gwet
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s005a1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b00650') then
      globp=snr
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b00650',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00040') then
      globp=z0
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00040',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b00620')then
      globp=raincu
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b00630',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
 

      globp=rainlp
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b00640',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
!
      call syslbl ('b00620',idtg,itau,ggdef,ihdg)
      do 98 jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 98 i=1,nxj
        globp(i,jj)=raincu(i,jj)+rainlp(i,jj)
 98   continue
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
!
      globp=rainlp
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b0062t',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00310') then
      globp=ss
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00310',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00320') then
      globp=rs
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00320',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00300') then
      do jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
        globp(i,jj)=ss(i,jj)-rs(i,jj)
      enddo
      enddo
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00300',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s003x0') then
      globp=rld
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s003x0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s003u0') then
      globp=sld
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s003u0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'x00330') then
      globp=plcl
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('x00330',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'x00340') then
      globp=cumtop
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('x00340',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00440') then
      globp=gfx
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00440',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00450') then
      globp=ugws
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00450',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s00460') then
      globp=vgws
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s00460',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'x00730') then
      call mpe_unify(acld,lev,my,2,mpe_double)
      do i=1,lev*my
       acld(i,1)=acld(i,1)*100.
      end do
      call syslbl ('x00730',idtg,itau,ggdef,ihdg)
      call dmswrit(lev,my,ihdg,lenc2,kflag,ifilout,acld,istat)
      go to 30
      endif
!
      if(label(kk).eq.'b00100') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,1),glob)
      call syslbl ('b00100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b00200') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,2),glob)
      call syslbl ('b00200',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b00210') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,3),glob)
      call syslbl ('b00210',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'x00590') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,4),glob)
      call syslbl ('x00590',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b00510') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,5),glob)
      call syslbl ('b00510',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
!      if(label(kk).eq.'b10510') then
!      call mpe2d_unify(glob,wk_xy(1,1,5))
!      call syslbl ('b10510',idtg,itau,ggdef,ihdg)
!      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
!      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
!      go to 30
!      endif
!
      if(label(kk).eq.'b02100') then
      globp=t2
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b02100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b02500') then
      globp=q2
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b02500',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b02510') then
      globp=rh2
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b02510',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b10200') then
      globp=u10
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b10200',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b10210') then
      globp=v10
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b10210',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'b10510') then
      globp=rh10
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('b10510',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'s005c0') then
      globp=canopy
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s005c0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
! s**5b0
!
!   s01??? = sa1??? = L01???   0-10cm
!   s02??? =                  10-200cm
!   s03??? = sa2??? = L02???  10-40cm
!   s04??? = sa3??? = L03???  40-100cm
!   s05??? = sa4??? = L04???  100-200cm
!
!
! 0-10cm
      if(label(kk).eq.'sa15b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,1),glob)
      call syslbl ('sa15b0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif

! 10-40cm
      if(label(kk).eq.'sa25b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,2),glob)
      call syslbl ('sa25b0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 40-100cm
      if(label(kk).eq.'sa35b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,3),glob)
      call syslbl ('sa35b0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 100-200cm
      if(label(kk).eq.'sa45b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,4),glob)
      call syslbl ('sa45b0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!kc >
! 0-10cm
      if(label(kk).eq.'s015b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,1),glob)
      call syslbl ('s015b0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif

      sfac2=3./19.
      sfac3=6./19.
      sfac4=10./19.
! output 10-200cm
      if(label(kk).eq.'s025b0') then
      do jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
       globp(i,jj)=soil_xy(i,jj,2)*sfac2+soil_xy(i,jj,3)*sfac3 &
                  +soil_xy(i,jj,4)*sfac4
      end do
      end do
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s025b0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!xb13 <
! 10-40cm
      if(label(kk).eq.'s035b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,2),glob)
      call syslbl ('s035b0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif

! 40-100cm
      if(label(kk).eq.'s045b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,3),glob)
      call syslbl ('s045b0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 100-200cm
      if(label(kk).eq.'s055b0') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,4),glob)
      call syslbl ('s055b0',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!xb13 >
!
! s**5b1
! 0-10cm
      if(label(kk).eq.'sa15b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,5),glob)
      call syslbl ('sa15b1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 10-40cm
      if(label(kk).eq.'sa25b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,6),glob)
      call syslbl ('sa25b1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 40-100cm
      if(label(kk).eq.'sa35b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,7),glob)
      call syslbl ('sa35b1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 100-200cm
      if(label(kk).eq.'sa45b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,8),glob)
      call syslbl ('sa45b1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!xb13
! 0-10cm
      if(label(kk).eq.'s015b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,5),glob)
      call syslbl ('s015b1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif

! 10-200cm
      if(label(kk).eq.'s025b1') then
      do jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
       globp(i,jj)=soil_xy(i,jj,6)*sfac2+soil_xy(i,jj,7)*sfac3 &
                  +soil_xy(i,jj,8)*sfac4
      end do
      end do
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s025b1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 10-40cm
      if(label(kk).eq.'s035b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,6),glob)
      call syslbl ('s035b1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 40-100cm
      if(label(kk).eq.'s045b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,7),glob)
      call syslbl ('s045b1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 100-200cm
      if(label(kk).eq.'s055b1') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,8),glob)
      call syslbl ('s055b1',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!<xb13
!
! s**100
! 0-10cm
      if(label(kk).eq.'sa1100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,9),glob)
      call syslbl ('sa1100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!10-40cm
      if(label(kk).eq.'sa2100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,10),glob)
      call syslbl ('sa2100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 40-100cm
      if(label(kk).eq.'sa3100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,11),glob)
      call syslbl ('sa3100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 100-200cm
      if(label(kk).eq.'sa4100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,12),glob)
      call syslbl ('sa4100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!xb13>
! 0-10cm
      if(label(kk).eq.'s01100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,9),glob)
      call syslbl ('s01100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 10-200cm
      if(label(kk).eq.'s02100') then

      do jj=1,jlistnum
         j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i=1,nxj
       globp(i,jj)=soil_xy(i,jj,10)*sfac2+soil_xy(i,jj,11)*sfac3 &
                  +soil_xy(i,jj,12)*sfac4
      end do
      end do
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call syslbl ('s02100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 10-40cm

      if(label(kk).eq.'s03100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,10),glob)
      call syslbl ('s03100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 40-100cm
      if(label(kk).eq.'s04100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,11),glob)
      call syslbl ('s04100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! 100-200cm
      if(label(kk).eq.'s05100') then
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,12),glob)
      call syslbl ('s05100',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!< xb13
!
! ctot_total cloud fraction
      if(label(kk).eq.'x00770') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,6),glob)
      call syslbl ('x00770',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! chig_high cloud fraction
      if(label(kk).eq.'x00760') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,7),glob)
      call syslbl ('x00760',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! cmid_middle cloud fraction
      if(label(kk).eq.'x00750') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,8),glob)
      call syslbl ('x00750',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
! clow_low cloud fraction
      if(label(kk).eq.'x00740') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,9),glob)
      call syslbl ('x00740',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
! pbl hight
      if(label(kk).eq.'pbl000') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,10),glob)
      call syslbl ('pbl000',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
! specific humidity 
      if(label(kk).eq.'m01500') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,11),glob)
      call syslbl ('m01500',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!
      if(label(kk).eq.'m60500') then
      call unify_reduceintp(nx,my,my_max,wk_xy(1,1,12),glob)
      call syslbl ('m60500',idtg,itau,ggdef,ihdg)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
      go to 30
      endif
!xb110> flash density
!      if(label(kk).eq.'fshden') then
!      call unify_reduceintp(nx,my,my_max,flash,glob)
!      call syslbl ('x00999',idtg,itau,ggdef,ihdg)
!      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
!      call split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
!      go to 30
!      endif
!xb110<

   30 continue
!
      if ( myrank .lt. nc )             &
         call dmswrit_split(nx,my,ihdg2,lenc,kflag,ifilout,mout,istat)
!
      return
      end
!---------------------------------------------------------------
      subroutine split(nx,my,lev,lenc,ifilout,nc,glob,mout,ihdg,ihdg2)
!
      use rank
      use const, only: RTYPE,kflag
!
      implicit none
!
      integer   nx,my,nc,istat,lenc,lev
      real(kind=RTYPE) glob(nx,my),mout(nx,my)
      character*26 ihdg,ihdg2
      character*80 ifilout

      if (myrank .eq. nc) then
        mout  = glob
        ihdg2 = ihdg
      endif
        nc    = nc + 1
!
      if ( nc .eq. lev .and. myrank .lt. nc ) then
        call dmswrit_split(nx,my,ihdg2,lenc,kflag,ifilout,mout,istat)
        nc    = 0
      endif
!
      return
      end
