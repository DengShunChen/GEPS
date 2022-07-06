      subroutine outflds_fv3(itau,nx,my,my_max,idtg,ggdef,ifilout       &
                            ,q2,fm,fh,fm10,fh2,srflag,ustar)
!
!  output additional variables for FV3
!
      use mpe
      use rank
      use index

      implicit  none

      integer   itau,nx,my,my_max,istat,lenc

      real      q2(nxp,my_max),ustar(nxp,my_max)        &
              , fm(nxp,my_max),fh(nxp,my_max)           &
              , fm10(nxp,my_max),fh2(nxp,my_max)        &
              , srflag(nxp,my_max)
         
      character ifilout*60, ggdef*4, ihdg*28
      integer*8 idtg
!
! local work arrays
!
      real      glob(nx,my)
!      integer   int_glob(nx,my)
      
!
      integer   jj,j,nxj,k,i
!
      character wtemp*6

      lenc =nx*my
!=======================================================================
      if(myrank .eq. 0) print*,'   in outflds_fv3 for tau= ',itau
!----------------------------------------------------------------------

!output q2
!byl      call mpe2d_unify(glob,q2)
      write(wtemp,'(a6)')'B02500'
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call unify_reduceintp(nx,my,my_max,q2,glob)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

!output fm
!byl      call mpe2d_unify(glob,fm)
      write(wtemp,'(a6)'),"S004F1"
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call unify_reduceintp(nx,my,my_max,fm,glob)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

!output fm10
!byl      call mpe2d_unify(glob,fm10)
      write(wtemp,'(a6)'),"S004F2"
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call unify_reduceintp(nx,my,my_max,fm10,glob)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

!output fh
!byl      call mpe2d_unify(glob,fh)
      write(wtemp,'(a6)'),"S004F3"
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call unify_reduceintp(nx,my,my_max,fh,glob)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

!output fh2
!byl      call mpe2d_unify(glob,fh2)
      write(wtemp,'(a6)'),"S004F4"
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call unify_reduceintp(nx,my,my_max,fh2,glob)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

!output ustar
!byl      call mpe2d_unify(glob,ustar)
      write(wtemp,'(a6)'),"S004F5"
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call unify_reduceintp(nx,my,my_max,ustar,glob)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

!output srflag
!byl      call mpe2d_unify(glob,srflag)
      write(wtemp,'(a6)'),"S001A0"
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call unify_reduceintp(nx,my,my_max,srflag,glob)
      call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

!=======================================================================
      return
      end
