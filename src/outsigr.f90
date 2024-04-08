      subroutine outsigr ( itau,nx,my,my_max,lev                         &
                         , idtg                                          & 
                         , fusl,fdsl,fuir,fdir                           & 
                         , fuslr,fdslr,fuirr,fdirr                       &
                         , asl,atl,asl_clr,atl_clr                       & 
                         , dtrad,cldcov,sd                               &
                         , ss,rs,olr,asol,sld,rld                        &
                         , ss_clr,rs_clr,olr_clr                         &
                         , asol_clr,sld_clr,rld_clr                      &
                         , cice,xtice,snr,sncover,snoalb                 &
                         , ctot,chig,cmid,clow                           &
                         , ggdef,gmdef)
      use index
      use mpe
      use const, only: RTYPE,kflag
!
      real(kind=RTYPE)  fusl(nxp,lev+1,my_max) ,fdsl(nxp,lev+1,my_max)   &
                      , fuir(nxp,lev+1,my_max) ,fdir(nxp,lev+1,my_max)   &
                      , fuslr(nxp,lev+1,my_max),fdslr(nxp,lev+1,my_max)  &
                      , fuirr(nxp,lev+1,my_max),fdirr(nxp,lev+1,my_max)  
      dimension  asl(nxp,lev,my_max),atl(nxp,lev,my_max)         &
              , asl_clr(nxp,lev,my_max),atl_clr(nxp,lev,my_max)  &
              , dtrad(nxp,lev,my_max),dtrad0(nxp,lev,my_max)     &
              , cldcov(nxp,lev,my_max)                           &
              , ss(nxp,my_max),rs(nxp,my_max)                    &
              , olr(nxp,my_max),asol(nxp,my_max)                 &
              , sld(nxp,my_max),rld(nxp,my_max)                  &
              , ss_clr(nxp,my_max),rs_clr(nxp,my_max)            &
              , olr_clr(nxp,my_max),asol_clr(nxp,my_max)         &
              , sld_clr(nxp,my_max),rld_clr(nxp,my_max)          &
              , cice(nxp,my_max),xtice(nxp,my_max)               &
              , snr(nxp,my_max),sncover(nxp,my_max)              &
              , snoalb(nxp,my_max),ctot(nxp,my_max)              &
              , chig(nxp,my_max),cmid(nxp,my_max)                &
              , clow(nxp,my_max) 
      real(kind=RTYPE) sd(nxp,lev,my_max),work(nx,my),work1(nxp,my_max)

      integer*8 idtg
      character*6 typ
      character*4 ggdef,gmdef
!
      do i = 1, nx*my
         work(i,1) = 0.
      enddo
!
      lenc=nx*my
!
!------------------------------------------------------------
! write 3D radiation flux
!------------------------------------------------------------
      do 30 k=1,lev+1
!
!     output upward solar flux for total sky
!
      do 20 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 20 i = 1,nxj
       work1(i,jj)=fusl(i,k,jj)
 20   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3A0")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
! 
!     output downward solar flux for total sky
!
      do 21 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 21 i = 1,nxj
       work1(i,jj)=fdsl(i,k,jj)
 21   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3B0")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output upward ir flux for total sky
! 
      do 22 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 22 i = 1,nxj
       work1(i,jj)=fuir(i,k,jj)
 22   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3C0")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output downward ir flux for total sky
! 
      do 23 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 23 i = 1,nxj
       work1(i,jj)=fdir(i,k,jj)
 23   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3D0")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output upward solar flux for clear sky
! 
      do 24 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 24 i = 1,nxj
       work1(i,jj)=fuslr(i,k,jj)
 24   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3A1")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output downward solar flux for clear sky
! 
      do 25 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 25 i = 1,nxj
       work1(i,jj)=fdslr(i,k,jj)
 25   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3B1")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output upward ir flux for clear sky
! 
      do 26 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 26 i = 1,nxj
       work1(i,jj)=fuirr(i,k,jj)
 26   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3C1")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output downward ir flux for clear sky
! 
      do 27 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 27 i = 1,nxj
       work1(i,jj)=fdirr(i,k,jj)
 27   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3D1")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output downward net solar flux for total sky
! 
      do 281 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 281 i = 1,nxj
       work1(i,jj)=fdsl(i,k,jj)-fusl(i,k,jj)
 281   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"310")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output downward net solar flux for clear sky
! 
      do 282 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 282 i = 1,nxj
       work1(i,jj)=fdslr(i,k,jj)-fuslr(i,k,jj)
 282   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"311")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output upward net ir flux for total sky
! 
      do 283 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 283 i = 1,nxj
       work1(i,jj)=fuir(i,k,jj)-fdir(i,k,jj)
 283   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"320")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output upward net ir flux for clear sky
! 
      do 284 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 284 i = 1,nxj
       work1(i,jj)=fuirr(i,k,jj)-fdirr(i,k,jj)
 284   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"321")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output downward net radiation flux for total sky
! 
      do 285 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 285 i = 1,nxj
       work1(i,jj)=fdsl(i,k,jj)-fusl(i,k,jj)+fdir(i,k,jj)-fuir(i,k,jj)
 285   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"300")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
!     output downward net radiation flux for clear sky
! 
      do 286 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 286 i = 1,nxj
       work1(i,jj)=fdslr(i,k,jj)-fuslr(i,k,jj)+fdirr(i,k,jj)-fuirr(i,k,jj)
 286   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"301")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)


 30   continue

!
      do 40 k  = 1, lev
      do 41 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 41 i = 1,nxj
       work1(i,jj)=asl(i,k,jj)
 41   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3E0")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
      do 42 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 42 i = 1,nxj
       work1(i,jj)=atl(i,k,jj)
 42   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3F0")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
      do 43 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 43 i = 1,nxj
       work1(i,jj)=asl_clr(i,k,jj)
 43   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3E1")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
      do 44 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 44 i = 1,nxj
       work1(i,jj)=atl_clr(i,k,jj)
 44   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3F1")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
      do 45 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 45 i = 1,nxj
       work1(i,jj)=dtrad(i,k,jj)
 45   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3G0")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
      do 46 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 46 i = 1,nxj
       work1(i,jj)=asl_clr(i,k,jj)+atl_clr(i,k,jj)
 46   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3G1")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
      do 47 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 47 i = 1,nxj
       work1(i,jj)=cldcov(i,k,jj)
 47   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"3H0")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
      do 48 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 48 i = 1,nxj
       work1(i,jj)=sd(i,k,jj)
 48   continue

      call unify_reduceintp(nx,my,my_max,work1,work)
      write(typ,'("m",i2.2,"220")')k
      call syslbl_w (typ,idtg,itau,gmdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
 40   continue
!------------------------------------------------------------
!  write 2D field
!------------------------------------------------------------
      work1=ss
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('s00310',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=rs
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('s00320',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=asol
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('X00330',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=olr
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('X00340',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=sld
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('S003U0',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=rld
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('S003X0',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=ss_clr
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('s0031C',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=rs_clr
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('s0032C',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=asol_clr
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('X0033C',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=olr_clr
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('X0034C',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=sld_clr
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('S003UC',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=rld_clr
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('S003XC',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=ctot
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('X00770',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=chig
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('X00760',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=cmid
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('X00750',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=clow
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('X00740',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=cice
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('W00093',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=xtice
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('W00094',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)

      work1=snr
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('b00650',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
 
      work1=sncover
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('B00651',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
 
      work1=snoalb
      call unify_reduceintp(nx,my,my_max,work1,work)
      call syslbl_w ('S0003X',idtg,itau,ggdef)
      call dmswrit(nx,my,lenc,kflag,work,istat)
!
      return
      end
