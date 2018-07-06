      subroutine out24 (nx,my,my_max,hf24,qf24,ss24,rs24,asol24,olr24  &
                      ,rain24,dt24,ifilout,glob,itau,idtg,ggdef)
!
      use index
      use mpe

      implicit  none

      integer   nx,my,my_max,itau
      real      dt24

      real      hf24(nxp,my_max),qf24(nxp,my_max),ss24(nxp,my_max),rs24(nxp,my_max), &
                asol24(nxp,my_max),olr24(nxp,my_max),rain24(nxp,my_max)

      real      glob(nx,my)
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
         glob(i,j)=qf24(i,jj)/dt24
        enddo
      enddo
      call mpe_unify(glob,nx,my,2,mpe_double)
      call syslbl ('s0043f',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call dmswrit(imax,jmax,ihdg,lenc,'H',ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          glob(i,j)=hf24(i,jj)/dt24
         enddo
      enddo
      call mpe_unify(glob,nx,my,2,mpe_double)
      call syslbl ('s0042f',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call dmswrit(imax,jmax,ihdg,lenc,'H',ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          glob(i,j)=ss24(i,jj)/dt24
         enddo
      enddo
      call mpe_unify(glob,nx,my,2,mpe_double)
      call syslbl ('s0031f',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call dmswrit(imax,jmax,ihdg,lenc,'H',ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          glob(i,j)=rs24(i,jj)/dt24
         enddo
      enddo
      call mpe_unify(glob,nx,my,2,mpe_double)
      call syslbl ('s0032f',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call dmswrit(imax,jmax,ihdg,lenc,'H',ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          glob(i,j)=rain24(i,jj)/dt24*2.5e+6
         enddo
      enddo
      call mpe_unify(glob,nx,my,2,mpe_double)
      call syslbl ('b0062f',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call dmswrit(imax,jmax,ihdg,lenc,'H',ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          glob(i,j)=asol24(i,jj)/dt24
         enddo
      enddo
      call mpe_unify(glob,nx,my,2,mpe_double)
      call syslbl ('x0033f',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call dmswrit(imax,jmax,ihdg,lenc,'H',ifilout,glob,istat)
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
          glob(i,j)=olr24(i,jj)/dt24
         enddo
      enddo
      call mpe_unify(glob,nx,my,2,mpe_double)
      call syslbl ('x0034f',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call dmswrit(imax,jmax,ihdg,lenc,'H',ifilout,glob,istat)
!
      return
      end
