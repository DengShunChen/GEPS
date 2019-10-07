      subroutine out2d_mfc (nx,lev,my,my_max,ifilout,itau,idtg,ntau     &
                            ,rain1,raintot,glob,t2,q2,rh2,rh10,u10,v10  &
                            ,tmax,tmin,td,rld,sld,ctot,slpty,ggdef )
!
      use rank
      use mpe
      use index
!
      implicit  none

      integer   nx,lev,my,my_max,itau,ntau

      real      raintot(nxp,my_max),t2(nxp,my_max),u10(nxp,my_max),   &
                v10(nxp,my_max),ctot(nxp,my_max),slpty(nxp,my_max)

      real rain1(nxp,my_max),q2(nxp,my_max),rh2(nxp,my_max),          &
           rh10(nxp,my_max),tmax(nxp,my_max),tmin(nxp,my_max),        &
           td(nxp,my_max),rld(nxp,my_max),sld(nxp,my_max)

      character*4 ggdef
      integer*8 idtg
!
      real      glob(nx,my)
!
      character*80 ifilout
      character*26 ihdg
!
      integer   num,n,levz,lenc,lenc2,i,ia,kk,j,nxj,istat,jj
      real      tnshun

!
      tnshun= 1.0
      lenc= nx*my
      lenc2= lev*my
!
! Total Precp.
      call mpe2d_unify(glob,raintot)
      call syslbl ('b0062t',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! 1hr Precp.
      call mpe2d_unify(glob,rain1)
      call syslbl ('b00621',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! T2
      call mpe2d_unify(glob,t2)
      call syslbl ('b02100',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! T2max
      call mpe2d_unify(glob,tmax)
      call syslbl ('b02171',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! T2min
      call mpe2d_unify(glob,tmin)
      call syslbl ('b02181',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! T2d
      call mpe2d_unify(glob,td)
      call syslbl ('b02150',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! q2
      call mpe2d_unify(glob,q2)
      call syslbl ('b02500',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! rh2
      call mpe2d_unify(glob,rh2)
      call syslbl ('b02510',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

!! rh10
!      call mpe2d_unify(glob,rh10)
!      call syslbl ('b10510',idtg,ntau,ggdef,ihdg)
!!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!      if( lreduce.eq.1 ) then
!        do jj = 1, jlistnum
!         j=jlist1(jj)
!         call reduceintp (glob(1,j),nxdef(j),nx,1)
!        enddo
!        call mpe_unify(glob,nx,my,5,mpe_double)
!      endif
!!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
!      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! U wind
      call mpe2d_unify(glob,u10)
      call syslbl ('b10200',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! V wind
      call mpe2d_unify(glob,v10)
      call syslbl ('b10210',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! ctot_total cloud fraction
      call mpe2d_unify(glob,ctot)
      call syslbl ('x00770',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
!
! Sea level pressure(hPa)
      call mpe2d_unify(glob,slpty)
      call syslbl ('ssl010',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! downward SW
      call mpe2d_unify(glob,sld)
      call syslbl ('s003u0',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! downwrd LW
      call mpe2d_unify(glob,rld)
      call syslbl ('s003x0',idtg,ntau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if( lreduce.eq.1 ) then
        do jj = 1, jlistnum
         j=jlist1(jj)
         call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
      endif
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

      return
      end

