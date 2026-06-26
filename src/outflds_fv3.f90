      subroutine outflds_fv3(itau,nx,my,my_max,idtg,ggdef       &
                            ,q2,fm,fh,fm10,fh2,srflag,ustar)
!
!  output additional variables for FV3
!
      use mpe
      use rank
      use index
      use const, only: RTYPE,kflag

      implicit  none

      integer   itau,nx,my,my_max,istat,lenc

      real      q2(nxp,my_max),ustar(nxp,my_max)        &
              , fm(nxp,my_max),fh(nxp,my_max)           &
              , fm10(nxp,my_max),fh2(nxp,my_max)        &
              , srflag(nxp,my_max)
         
      character ggdef*4
      integer*8 idtg
!
! local work arrays
!
      real(kind=RTYPE) glob(nx,my),globp(nxp,my_max)
!      integer   int_glob(nx,my)
      
!
      integer   jj,j,nxj,k,i
!
      character wtemp*6

!$acc enter data create(globp,glob)
      lenc =nx*my
!=======================================================================
      if(myrank .eq. 0) print*,'   in outflds_fv3 for tau= ',itau
!----------------------------------------------------------------------

!output q2
      write(wtemp,'(a6)')'B02500'
      call syslbl_w(wtemp,idtg,itau,ggdef)
      globp=q2
      !$acc update device(globp)
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call dmswrit(nx,my,lenc,kflag,glob,istat)

!output fm
      write(wtemp,'(a6)'),"S004F1"
      call syslbl_w(wtemp,idtg,itau,ggdef)
      globp=fm
      !$acc update device(globp)
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call dmswrit(nx,my,lenc,kflag,glob,istat)

!output fm10
      write(wtemp,'(a6)'),"S004F2"
      call syslbl_w(wtemp,idtg,itau,ggdef)
      globp=fm10
      !$acc update device(globp)
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call dmswrit(nx,my,lenc,kflag,glob,istat)

!output fh
      write(wtemp,'(a6)'),"S004F3"
      call syslbl_w(wtemp,idtg,itau,ggdef)
      globp=fh
      !$acc update device(globp)
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call dmswrit(nx,my,lenc,kflag,glob,istat)

!output fh2
      write(wtemp,'(a6)'),"S004F4"
      call syslbl_w(wtemp,idtg,itau,ggdef)
      globp=fh2
      !$acc update device(globp)
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call dmswrit(nx,my,lenc,kflag,glob,istat)

!output ustar
      write(wtemp,'(a6)'),"S004F5"
      call syslbl_w(wtemp,idtg,itau,ggdef)
      globp=ustar
      !$acc update device(globp)
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call dmswrit(nx,my,lenc,kflag,glob,istat)

!output srflag
      write(wtemp,'(a6)'),"S001A0"
      call syslbl_w(wtemp,idtg,itau,ggdef)
      globp=srflag
      !$acc update device(globp)
      call unify_reduceintp(nx,my,my_max,globp,glob)
      call dmswrit(nx,my,lenc,kflag,glob,istat)

!=======================================================================
!$acc exit data delete(globp,glob)
      return
      end
