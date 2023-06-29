      subroutine rayleifr(nx,my,my_max,lev,rad,cosl,dt,ut,vt)

      use mpe
      use rank
      use index
      use const, only: RTYPE

      implicit  none
      integer   nx,my,my_max,lev
      real      rad

      real(kind=RTYPE) cosl(my)
      real(kind=RTYPE) ut(nxp,lev,my_max),vt(nxp,lev,my_max)
      real      wmax(lev),wmaxtmp
!
      integer, parameter :: levtop=6,lev2=3

      real     ckdx,day,windmax
!     data ckdx/20./, day/1./
      data ckdx/40./, day/1./

      data windmax/140./
!     data windmax/130./

      logical dofric

      integer jj,j,nxj,k,i   
      real    xx,dt,ckdy,frictime,fac,fricd,cdx,cdy
!
      wmax(1:lev)=0.
      wmaxtmp=windmax
      do 10 jj =1,jlistnum
        j=jlist1(jj)
!ch     nxj=nxdef(j)
        nxj=nxdef_2d(j)
        xx=rad/cosl(j)
      do 10 k=1,lev
      do 10 i=1,nxj
        wmax(k)= max(wmax(k),xx*sqrt(ut(i,k,jj)**2+vt(i,k,jj)**2))
 10   continue
      call mpe_global_max(wmax,lev,mpe_double)
      dofric=.false.
      do k=1,levtop
        wmaxtmp=max(wmax(k),wmaxtmp)
        if( wmax(k) .gt. windmax )then
!          dofric=.true.
          if(myrank.eq.0)print *,' rayleifr k=',k,' wmax=',wmax(k)
!           ' dofric=',dofric
        endif
      enddo
      if(wmaxtmp.gt.windmax)  dofric=.true.
      if(myrank.eq.0)print *,' rayleifr  dofric=',dofric
!
      ckdy=ckdx
      frictime=1./(day*86400.)
      if(dofric)then
        do k=1,levtop
          fac=1.+exp(-1.0*(k-lev2)/lev2)
          if(fac.gt.3.)fac=3.
          fricd=frictime*fac
          cdx=ckdx*fricd
          cdy=ckdy*fricd
          do jj = 1, jlistnum
            j=jlist1(jj)
!ch         nxj=nxdef(j)
            nxj=nxdef_2d(j)
          do i = 1, nxj
            ut(i,k,jj)=ut(i,k,jj)-cdx*ut(i,k,jj)
            vt(i,k,jj)=vt(i,k,jj)-cdy*vt(i,k,jj)
          enddo
          enddo
        enddo
      endif
!
      return
      end
