      subroutine out2d_mfc (nx,lev,my,my_max,ifilout,itau,idtg,ntau &
       ,raintot,glob,t2,u10,v10,ctot,slpty,ggdef) 
!
      use rank
      use mpe
      use index
!
      implicit  none

      integer   nx,lev,my,my_max,itau,ntau 

      real      raintot(nxp,my_max),t2(nxp,my_max),u10(nxp,my_max),&
                v10(nxp,my_max),ctot(nxp,my_max),slpty(nxp,my_max)

      character*4 ggdef
      integer*8 idtg
!
      real      glob(nx,my)
!
      character*80 ifilout
      character*26 ihdg
!
      integer   num,n,levz,lenc,lenc2,i,ia,kk,j,nxj,istat
      real      tnshun

!
      tnshun= 1.0
      lenc= nx*my
      lenc2= lev*my
!

! Total Precp.
      call mpe2d_unify(glob,raintot)
      call syslbl ('b0062t',idtg,ntau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! T2
      call mpe2d_unify(glob,t2)
      call syslbl ('b02100',idtg,ntau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! U wind
      call mpe2d_unify(glob,u10)
      call syslbl ('b10200',idtg,ntau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! V wind
      call mpe2d_unify(glob,v10)
      call syslbl ('b10210',idtg,ntau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)

! ctot_total cloud fraction
      call mpe2d_unify(glob,ctot)
      call syslbl ('x00770',idtg,ntau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
!
! Sea level pressure(hPa)
      call mpe2d_unify(glob,slpty)
      call syslbl ('ssl010',idtg,ntau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      return
      end
