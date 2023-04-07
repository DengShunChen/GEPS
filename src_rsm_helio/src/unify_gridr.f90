      subroutine unify_gridr                                            &
      ( tmr1,tmr2,tmr3,tmr4,tmr5,tmr6,tmr7,tmr8,tms1,tms2,tms3,tms4     &
      , fusl,fdsl,fuir,fdir,fuslr,fdslr,fuirr,fdirr,cldcov,sd           &
      , asl_clr,atl_clr,asol,olr,sld                                    &
      , asol_clr,olr_clr,rld_clr,sld_clr,ss_clr,rs_clr                  &
      , nx,my,my_max,lev) 
!------------------------------------------------------------------------
      use index
      use mpe
!
      dimension asl_clr(nx,lev,my_max),atl_clr(nx,lev,my_max),          &
                cldcov(nx,lev,my_max),sd(nx,lev,my_max)
      
      dimension fusl(nx,lev+1,my_max),fdsl(nx,lev+1,my_max),            &
                fuir(nx,lev+1,my_max),fdir(nx,lev+1,my_max),            &
                fuslr(nx,lev+1,my_max),fdslr(nx,lev+1,my_max),          &
                fuirr(nx,lev+1,my_max),fdirr(nx,lev+1,my_max)

      dimension ss_clr(nx,my_max),rs_clr(nx,my_max),                    &
                asol_clr(nx,my_max),olr_clr(nx,my_max),                 &
                rld_clr(nx,my_max),sld_clr(nx,my_max),                  &
                asol(nx,my_max),olr(nx,my_max),sld(nx,my_max)

      dimension tmr1(nx,lev+1,my), tmr2(nx,lev+1,my)
      dimension tmr3(nx,lev+1,my), tmr4(nx,lev+1,my)
      dimension tmr5(nx,lev+1,my), tmr6(nx,lev+1,my)
      dimension tmr7(nx,lev+1,my), tmr8(nx,lev+1,my)
      dimension tmr9(nx,my,9)

      dimension tms1(nx,lev,my), tms2(nx,lev,my)
      dimension tms3(nx,lev,my), tms4(nx,lev,my)
! 
      call mpe_unify_1(tmr9(1,1,1),asol,nx,my,2,mpe_double)
      call mpe_unify_1(tmr9(1,1,2),olr,nx,my,2,mpe_double)
      call mpe_unify_1(tmr9(1,1,3),sld,nx,my,2,mpe_double)
      call mpe_unify_1(tmr9(1,1,4),ss_clr,nx,my,2,mpe_double)
      call mpe_unify_1(tmr9(1,1,5),rs_clr,nx,my,2,mpe_double)
      call mpe_unify_1(tmr9(1,1,6),asol_clr,nx,my,2,mpe_double)
      call mpe_unify_1(tmr9(1,1,7),olr_clr,nx,my,2,mpe_double)
      call mpe_unify_1(tmr9(1,1,8),rld_clr,nx,my,2,mpe_double)
      call mpe_unify_1(tmr9(1,1,9),sld_clr,nx,my,2,mpe_double)

      do jj =1, jlistnum
      j=jlist1(jj)
      do k = 1, lev+1
      do i = 1, nx
        tmr1(i,k,j)  = fusl(i,k,jj)
        tmr2(i,k,j)  = fdsl(i,k,jj)
        tmr3(i,k,j)  = fuir(i,k,jj)
        tmr4(i,k,j)  = fdir(i,k,jj)
        tmr5(i,k,j)  = fuslr(i,k,jj)
        tmr6(i,k,j)  = fdslr(i,k,jj)
        tmr7(i,k,j)  = fuirr(i,k,jj)
        tmr8(i,k,j)  = fdirr(i,k,jj)
      enddo
      enddo
      enddo
!
      do jj =1, jlistnum
      j=jlist1(jj)
      do k = 1, lev
      do i = 1, nx
        tms1(i,k,j)  = cldcov(i,k,jj)
        tms2(i,k,j)  = asl_clr(i,k,jj)
        tms3(i,k,j)  = atl_clr(i,k,jj)
        tms4(i,k,j)  = sd(i,k,jj)
      enddo
      enddo
      enddo
!
      call mpe_unify(tmr1,nx*lev+1,my,2,mpe_double)
      call mpe_unify(tmr2,nx*lev+1,my,2,mpe_double)
      call mpe_unify(tmr3,nx*lev+1,my,2,mpe_double)
      call mpe_unify(tmr4,nx*lev+1,my,2,mpe_double)
      call mpe_unify(tmr5,nx*lev+1,my,2,mpe_double)
      call mpe_unify(tmr6,nx*lev+1,my,2,mpe_double)
      call mpe_unify(tmr7,nx*lev+1,my,2,mpe_double)
      call mpe_unify(tmr8,nx*lev+1,my,2,mpe_double)
!
      call mpe_unify(tms1,nx*lev,my,2,mpe_double)
      call mpe_unify(tms2,nx*lev,my,2,mpe_double)
      call mpe_unify(tms3,nx*lev,my,2,mpe_double)
      call mpe_unify(tms4,nx*lev,my,2,mpe_double)

      return
      end
