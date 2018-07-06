      subroutine outsigr ( itau,nx,my,my_max,lev                         &
                         , idtg,ifilout                                  & 
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
!
      dimension fusl(nxp,lev+1,my_max),fdsl(nxp,lev+1,my_max)    &
              , fuir(nxp,lev+1,my_max),fdir(nxp,lev+1,my_max)    &
              , fuslr(nxp,lev+1,my_max),fdslr(nxp,lev+1,my_max)  &
              , fuirr(nxp,lev+1,my_max),fdirr(nxp,lev+1,my_max)  &
              , asl(nxp,lev,my_max),atl(nxp,lev,my_max)          &
              , asl_clr(nxp,lev,my_max),atl_clr(nxp,lev,my_max)  &
              , dtrad(nxp,lev,my_max),dtrad0(nxp,lev,my_max)     &
              , cldcov(nxp,lev,my_max),sd(nxp,lev,my_max)        &
              , ss(nxp,my_max),rs(nxp,my_max)                    &
              , olr(nxp,my_max),asol(nxp,my_max)                 &
              , sld(nxp,my_max),rld(nxp,my_max)                  &
              , ss_clr(nxp,my_max),rs_clr(nxp,my_max)            &
              , olr_clr(nxp,my_max),asol_clr(nxp,my_max)         &
              , sld_clr(nxp,my_max),rld_clr(nxp,my_max)          &
              , cice(nxp,my_max),xtice(nxp,my_max),snr(nxp,my_max)  &
              , sncover(nxp,my_max),snoalb(nxp,my_max)           &
              , ctot(nxp,my_max),chig(nxp,my_max),cmid(nxp,my_max),clow(nxp,my_max) &
              , work(nx,my)

      integer*8 idtg
      character*80 ifilout
      character typ*6,ihdg*26
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
       work(i,j)=fusl(i,k,jj)
 20   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3A0")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
! 
!     output downward solar flux for total sky
!
      do 21 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 21 i = 1,nxj
       work(i,j)=fdsl(i,k,jj)
 21   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3B0")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output upward ir flux for total sky
! 
      do 22 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 22 i = 1,nxj
       work(i,j)=fuir(i,k,jj)
 22   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3C0")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output downward ir flux for total sky
! 
      do 23 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 23 i = 1,nxj
       work(i,j)=fdir(i,k,jj)
 23   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3D0")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output upward solar flux for clear sky
! 
      do 24 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 24 i = 1,nxj
       work(i,j)=fuslr(i,k,jj)
 24   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3A1")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output downward solar flux for clear sky
! 
      do 25 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 25 i = 1,nxj
       work(i,j)=fdslr(i,k,jj)
 25   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3B1")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output upward ir flux for clear sky
! 
      do 26 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 26 i = 1,nxj
       work(i,j)=fuirr(i,k,jj)
 26   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3C1")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output downward ir flux for clear sky
! 
      do 27 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 27 i = 1,nxj
       work(i,j)=fdirr(i,k,jj)
 27   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3D1")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output downward net solar flux for total sky
! 
      do 281 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 281 i = 1,nxj
       work(i,j)=fdsl(i,k,jj)-fusl(i,k,jj)
 281   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"310")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output downward net solar flux for clear sky
! 
      do 282 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 282 i = 1,nxj
       work(i,j)=fdslr(i,k,jj)-fuslr(i,k,jj)
 282   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"311")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output upward net ir flux for total sky
! 
      do 283 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 283 i = 1,nxj
       work(i,j)=fuir(i,k,jj)-fdir(i,k,jj)
 283   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"320")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output upward net ir flux for clear sky
! 
      do 284 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 284 i = 1,nxj
       work(i,j)=fuirr(i,k,jj)-fdirr(i,k,jj)
 284   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"321")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output downward net radiation flux for total sky
! 
      do 285 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 285 i = 1,nxj
       work(i,j)=fdsl(i,k,jj)-fusl(i,k,jj)+fdir(i,k,jj)-fuir(i,k,jj)
 285   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"300")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
!     output downward net radiation flux for clear sky
! 
      do 286 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 286 i = 1,nxj
       work(i,j)=fdslr(i,k,jj)-fuslr(i,k,jj)+fdirr(i,k,jj)-fuirr(i,k,jj)
 286   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"301")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)


 30   continue

!
      do 40 k  = 1, lev
      do 41 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
!
      do 41 i = 1,nxj
       work(i,j)=asl(i,k,jj)
 41   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3E0")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 42 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 42 i = 1,nxj
       work(i,j)=atl(i,k,jj)
 42   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3F0")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 43 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 43 i = 1,nxj
       work(i,j)=asl_clr(i,k,jj)
 43   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3E1")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 44 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 44 i = 1,nxj
       work(i,j)=atl_clr(i,k,jj)
 44   continue
      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3F1")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 45 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 45 i = 1,nxj
       work(i,j)=dtrad(i,k,jj)
 45   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3G0")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 46 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 46 i = 1,nxj
       work(i,j)=asl_clr(i,k,jj)+atl_clr(i,k,jj)
 46   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3G1")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 47 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 47 i = 1,nxj
       work(i,j)=cldcov(i,k,jj)
 47   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"3H0")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      do 48 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 48 i = 1,nxj
       work(i,j)=sd(i,k,jj)
 48   continue

      call mpe_unify(work,nx,my,2,mpe_double)
      write(typ,'("m",i2.2,"220")')k
      call syslbl (typ,idtg,itau,gmdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
 40   continue
!------------------------------------------------------------
!  write 2D field
!------------------------------------------------------------
      call mpe2d_unify(work,ss)
      call syslbl ('s00310',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,rs)
      call syslbl ('s00320',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,asol)
      call syslbl ('X00330',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,olr)
      call syslbl ('X00340',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,sld)
      call syslbl ('S003U0',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,rld)
      call syslbl ('S003X0',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,ss_clr)
      call syslbl ('s0031C',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,rs_clr)
      call syslbl ('s0032C',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,asol_clr)
      call syslbl ('X0033C',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,olr_clr)
      call syslbl ('X0034C',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,sld_clr)
      call syslbl ('S003UC',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,rld_clr)
      call syslbl ('S003XC',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,ctot)
      call syslbl ('X00770',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,chig)
      call syslbl ('X00760',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,cmid)
      call syslbl ('X00750',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,clow)
      call syslbl ('X00740',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,cice)
      call syslbl ('W00093',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,xtice)
      call syslbl ('W00094',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)

      call mpe2d_unify(work,snr)
      call syslbl ('b00650',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
 
      call mpe2d_unify(work,sncover)
      call syslbl ('B00651',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
 
      call mpe2d_unify(work,snoalb)
      call syslbl ('S0003X',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (work,nxdef,nx,my)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,work,istat)
!
      return
      end
