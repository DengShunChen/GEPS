      subroutine out24 (nx,my,my_max,hf24,qf24,ss24,rs24,asol24,olr24  &
                      ,rain24,dt24,ifilout,glob,itau,idtg,ggdef,flash24)
!
      use index
      use mpe
      use const, only: RTYPE,kflag

      implicit  none

      integer   nx,my,my_max,itau
      real      dt24

      real      hf24(nxp,my_max),qf24(nxp,my_max),ss24(nxp,my_max),rs24(nxp,my_max), &
                asol24(nxp,my_max),olr24(nxp,my_max),rain24(nxp,my_max)              &
               ,flash24(nxp,my_max)

      real(kind=RTYPE) glob(nx,my),wrk(nxp,my_max)
!
      integer*8 idtg
      character*80 ifilout
      character*26 ihdg
      character*4  ggdef

      integer   imax,jmax,lenc,j,nxj,i,istat,jj
!
      imax=nx
      jmax=my
      lenc= imax*jmax
!
      do jj = 1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
         wrk(i,jj)=qf24(i,jj)/dt24
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('s0043f',idtg,itau,ggdef,ihdg)
      call dmswrit(imax,jmax,ihdg,lenc,kflag,ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=hf24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('s0042f',idtg,itau,ggdef,ihdg)
      call dmswrit(imax,jmax,ihdg,lenc,kflag,ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=ss24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('s0031f',idtg,itau,ggdef,ihdg)
      call dmswrit(imax,jmax,ihdg,lenc,kflag,ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=rs24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('s0032f',idtg,itau,ggdef,ihdg)
      call dmswrit(imax,jmax,ihdg,lenc,kflag,ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=rain24(i,jj)/dt24*2.5e+6
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('b0062f',idtg,itau,ggdef,ihdg)
      call dmswrit(imax,jmax,ihdg,lenc,kflag,ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=asol24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('x0033f',idtg,itau,ggdef,ihdg)
      call dmswrit(imax,jmax,ihdg,lenc,kflag,ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          wrk(i,jj)=olr24(i,jj)/dt24
         enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('x0034f',idtg,itau,ggdef,ihdg)
      call dmswrit(imax,jmax,ihdg,lenc,kflag,ifilout,glob,istat)
!
!xb110>>
      do jj = 1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
         wrk(i,jj)=flash24(i,jj)/dt24 !xb110, 24hr average flash density (km-2day-1)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl ('x00999',idtg,itau,ggdef,ihdg)
      call dmswrit(imax,jmax,ihdg,lenc,kflag,ifilout,glob,istat)
!xb110<<
      return
      end
