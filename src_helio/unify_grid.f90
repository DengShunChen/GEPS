      subroutine unify_grid                                             &
          ( snr,gwr,tg,tm1,tm2,ss,rs,tm3,tm4,ustar,tstar,qstar          &
          , hflux,qflux,raincu,rainlp,totalp,curate,plcl,cumtop         &
          , tgclim,gwet,z0,alb,land,ice,ocean,gwclim,acld               &
          , tmc1,tmc2,tmc3,tmc4,tmc5,tmc6,fpsp,ftp,fqp,fpsp1,ftp1,fqp1  &
          , e,eps,o3l,dtrad,pt,ptend,ugws,vgws,nx,my,my_max,lev         &
!soil
          , smc,slc,stc,canopy,sigmaf,istyp,ivegtyp,km_soil,temp1,temp2 &
!    6    , temp3,rld)
          , temp3,rld,zice,asl,atl,sfalb,sfemis)

!
!     include "../include/index.h"
!     include "../include/mpe.h"
      use index
      use mpe
!
      dimension e(nx,lev,my_max),eps(nx,lev,my_max)                     &
      ,         o3l(nx,lev,my_max),dtrad(nx,lev,my_max)                 &
      ,         ftp(nx,lev,my_max)                                      &
      ,         ftp1(nx,lev,my_max)                                     &
      ,         fqp(nx,lev,my_max)                                      &
      ,         fqp1(nx,lev,my_max)                                     &
!
! add for new_sas_pbl
      ,         asl(nx,lev,my_max)                                      &
      ,         atl(nx,lev,my_max)
!soil
      dimension smc(nx,km_soil,my_max),stc(nx,km_soil,my_max)           &
               ,temp1(nx,km_soil,my),temp2(nx,km_soil,my)               &
               ,canopy(nx,my_max),sigmaf(nx,my_max),istyp(nx,my)        &
               ,ivegtyp(nx,my)                                          &
! noah
               ,slc(nx,km_soil,my_max),temp3(nx,km_soil,my)             &
               ,zice(nx,my_max),sfalb(nx,my_max),sfemis(nx,my_max)
!
      dimension snr(nx,my_max),gwr(nx,my_max),tg(nx,my_max)             &
      , ts(nx,my_max),ss(nx,my_max),rs(nx,my_max)                       &
      , ustar(nx,my_max),tstar(nx,my_max),qstar(nx,my_max)              &
      , hflux(nx,my_max),qflux(nx,my_max),raincu(nx,my_max),rainlp(nx,my_max)   &
      , totalp(nx,my_max),curate(nx,my_max),plcl(nx,my_max),cumtop(nx,my_max)   &
      , tgclim(nx,my_max),gwet(nx,my_max),z0(nx,my_max),alb(nx,my_max),land(nx,my)  &
      , ice(nx,my),ocean(nx,my),gwclim(nx,my_max)                       &
      , acld(lev,my),rld(nx,my_max)                                     &
      , fpsp(nx,my_max),fpsp1(nx,my_max)
!
      dimension tm1(nx,lev,my), tm2(nx,lev,my), tm3(nx,lev,my)          &
              , tm4(nx,lev,my),tm2d(nx,my,40)
! hmhj for ncld=2
      dimension tmc1(nx,lev,my), tmc2(nx,lev,my)
      dimension tmc3(nx,lev,my), tmc4(nx,lev,my)
      dimension tmc5(nx,lev,my), tmc6(nx,lev,my)
!
      dimension pt(nx,my),ptend(nx,my)
      dimension ugws(nx,my),vgws(nx,my)
!
      logical land, ice, ocean
!
      call mpe_unify_1(tm2d(1,1,1),snr,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,2),gwr,nx,my,2,mpe_double)

      call mpe_unify_1(tm2d(1,1,3),tg,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,4),ts,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,5),ss,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,6),rs,nx,my,2,mpe_double)
!soil
      call mpe_unify_1(tm2d(1,1,7),rld,nx,my,2,mpe_double)
!soil
      call mpe_unify_1(tm2d(1,1,8),ustar,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,9),tstar,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,10),qstar,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,11),hflux,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,12),qflux,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,13),raincu,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,14),rainlp,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,15),totalp,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,16),curate,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,17),cumtop,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,18),plcl,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,19),tgclim,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,20),gwet,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,21),z0,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,22),alb,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,23),fpsp,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,24),fpsp1,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,25),zice,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,26),gwclim,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,27),ptend,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,28),pt,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,29),ugws,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,30),vgws,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,31),canopy,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,32),sigmaf,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,33),sfalb,nx,my,2,mpe_double)
!
      call mpe_unify_1(tm2d(1,1,34),sfemis,nx,my,2,mpe_double)
!
      call mpe_unify(land,nx,my,2,mpe_logical)
!
      call mpe_unify(ice,nx,my,2,mpe_logical)
!
      call mpe_unify(ocean,nx,my,2,mpe_logical)
!
      call mpe_unify(acld,lev,my,2,mpe_double)
!
      do jj =1, jlistnum
      j=jlist1(jj)
      do k = 1, lev
      do i = 1, nx
        tm1(i,k,j)  = e(i,k,jj)
        tm2(i,k,j)  = eps(i,k,jj)
        tm3(i,k,j)  = o3l(i,k,jj)
        tm4(i,k,j)  = dtrad(i,k,jj)
        tmc1(i,k,j)  = ftp(i,k,jj)
        tmc2(i,k,j)  = fqp(i,k,jj)
        tmc3(i,k,j)  = ftp1(i,k,jj)
        tmc4(i,k,j)  = fqp1(i,k,jj)
        tmc5(i,k,j)  = asl(i,k,jj)
        tmc6(i,k,j)  = atl(i,k,jj)

      enddo
      enddo
      enddo
!
      call mpe_unify(tm1,nx*lev,my,2,mpe_double)
      call mpe_unify(tm2,nx*lev,my,2,mpe_double)
      call mpe_unify(tm3,nx*lev,my,2,mpe_double)
      call mpe_unify(tm4,nx*lev,my,2,mpe_double)
      call mpe_unify(tmc1,nx*lev,my,2,mpe_double)
      call mpe_unify(tmc2,nx*lev,my,2,mpe_double)
      call mpe_unify(tmc3,nx*lev,my,2,mpe_double)
      call mpe_unify(tmc4,nx*lev,my,2,mpe_double)
      call mpe_unify(tmc5,nx*lev,my,2,mpe_double)
      call mpe_unify(tmc6,nx*lev,my,2,mpe_double)
!
!soil
      call mpe_unify(ivegtyp,nx,my,2,mpe_integer)
      call mpe_unify(istyp,nx,my,2,mpe_integer)
!
      do jj =1, jlistnum
      j=jlist1(jj)
      do k = 1,km_soil
      do i = 1, nx
        temp1(i,k,j)  = smc(i,k,jj)
        temp2(i,k,j)  = stc(i,k,jj)
        temp3(i,k,j)  = slc(i,k,jj)
      enddo
      enddo
      enddo
      call mpe_unify(temp1,nx*km_soil,my,2,mpe_double)
      call mpe_unify(temp2,nx*km_soil,my,2,mpe_double)
      call mpe_unify(temp3,nx*km_soil,my,2,mpe_double)
!soil
      return
      end
