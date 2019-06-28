      subroutine incrini
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use param
      use rank
      use index
      use const
      use grid
      use spec
      use fftcom
!
      implicit  none

      real      hld1(nx,my),pt1(nx,my_max)
      character typ*6,lrec*26
      real      cc(nx+2,levp,1,my_max)
!!      real      cc(nx+2,levp,3+ncld,my_max),wss(levp,2,3+ncld,jtrun,jtmax)

      integer   k,ii,jj,j,nxj,lmax,itaup,lncrec,istat,i,kk,m,n,mf,ntrac
      real      fac,dummy
!
      lmax=16
!
      itaup = taup + 0.001
      lncrec = nx*my
!
      do k = 1, lev
        write (typ, '("m",i2.2,"200")' ) k
        call syslbl (typ,idtg2,itaup,gmdef,lrec)
        call dmsread (nx,my,lrec,lncrec,'H',ifilin,hld1,istat)
        if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
        do jj = 1, jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
          fac = cosl(j)/rad
          do i = 1,nxj
            ut(i,k,jj) = hld1(ii,j)*fac
            ii=ii+1
          enddo
        enddo
      enddo
!
      do k = 1, lev
        write (typ, '("m",i2.2,"210")' ) k
        call syslbl (typ,idtg2,itaup,gmdef,lrec)
        call dmsread (nx,my,lrec,lncrec,'H',ifilin,hld1,istat)
        if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
        do jj = 1, jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
          fac = cosl(j)/rad
          do i = 1,nxj
            vt(i,k,jj) = hld1(ii,j)*fac
            ii=ii+1
          enddo
        enddo
      enddo
!
      call syslbl ('B00010',idtg2,itaup,ggdef,lrec)
      call dmsread (nx,my,lrec,lncrec,'H',ifilin,hld1,istat)
      if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
      do jj = 1, jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        do i = 1, nxj
          pt(i,jj) = hld1(ii,j) - ptop
          ii=ii+1
        enddo
!ch
        pt1(:,jj) = hld1(:,j) - ptop
      enddo
!
      do k = 1, lev
        write (typ, '("m",i2.2,"100")' ) k
        call syslbl (typ,idtg2,itaup,gmdef,lrec)
        call dmsread (nx,my,lrec,lncrec,'H',ifilin,hld1,istat)
        if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
        do jj = 1, jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
          do i = 1,nxj
            tt(i,k,jj) = hld1(ii,j)
            ii=ii+1
          enddo
        enddo
      enddo
!
      do k = 1, lev
        write (typ, '("m",i2.2,"500")' ) k
        call syslbl (typ,idtg2,itaup,gmdef,lrec)
        call dmsread (nx,my,lrec,lncrec,'H',ifilin,hld1,istat)
        if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
        do jj = 1, jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
          do i = 1,nxj
            qt(i,k,jj) = hld1(ii,j)
            ii=ii+1
          enddo
        enddo
      enddo
!
      if( ncld .ge. 2 ) then
        ntrac=2
        do k = 1, lev
          kk = (ntrac-1)*lev+k
          write (typ, '("m",i2.2,"550")' ) k     ! cloud liquid water content
          call syslbl (typ,idtg2,itaup,gmdef,lrec)
          call dmsread (nx,my,lrec,lncrec,'H',ifilin,hld1,istat)
          if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
          do jj = 1, jlistnum
            j=jlist1(jj)
            ii=nxjstart(j)
            nxj=nxdef_2d(j)
            do i = 1,nxj
              qt(i,kk,jj) = hld1(ii,j)
              ii=ii+1
            enddo
          enddo
        enddo
        if(ncld.ge.3)then
        ntrac=3
        do k = 1, lev
          kk = (ntrac-1)*lev+k
          write (typ, '("m",i2.2,"560")' ) k
          call syslbl (typ,idtg2,itaup,gmdef,lrec)
          call dmsread (nx,my,lrec,lncrec,'H',ifilin,hld1,istat)
          if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
          do jj = 1, jlistnum
            j=jlist1(jj)
            ii=nxjstart(j)
            nxj=nxdef_2d(j)
            do i = 1,nxj
              qt(i,kk,jj) = hld1(ii,j)
              ii=ii+1
            enddo
          enddo
        enddo
        endif
      endif
!
      do k = 1, lev
        write (typ, '("m",i2.2,"000")' ) k
        call syslbl (typ,idtg2,itaup,gmdef,lrec)
        call dmsread (nx,my,lrec,lncrec,'H',ifilin,hld1,istat)
        if( lreduce.eq.1 ) call reducepick (hld1,nxdef,nx,my)
        do jj = 1, jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
          do i = 1,nxj
            phi(i,k,jj) = hld1(ii,j)*grav
            ii=ii+1
          enddo
        enddo
      enddo
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,pt(1,jj) &
                        ,pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj) )
      enddo
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do k=1,lev
          do i=1,nxj
            tt(i,k,jj)= tt(i,k,jj)*(1.0+0.608*qt(i,k,jj))/pk(i,k,jj)
          enddo
        enddo
      enddo

      call joinrs(cc,tt,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc   &
                   ,temnow,1,nsizey)
!ch   call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,pt      &
      call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,pt1     &
                 ,plnow,nsizey)
      call trandv(jtrun,jtmax,nx,my,my_max,lev,ut,vt,weight,cim &
                 ,onocos,poly,dpoly,vornow,divnow,nsizey)
!
      do m=1,mlistnum
        mf=mlist(m)
        if (mf.eq.1) then
          do k = 1, levp
            divnow(k,1,1,m)= 0.0
            divnow(k,2,1,m)= 0.0
            vornow(k,1,1,m)= 0.0
            vornow(k,2,1,m)= 0.0
          enddo
        endif
      enddo
!
!!      call joinsr(wss,vornow,divnow,temnow,dummy,jtrun,jtmax,levp  &
!!                 ,mlistnum,3,1)
!!      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,wss,cc,3,nsizey)
!!      call ujoinsr(cc,rvor,rdiv,tt,dummy,nx,my_max,lev,jlistnum,3,1)
      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,vornow,cc,1,nsizey)
      call ujoinsr(cc,rvor,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,divnow,cc,1,nsizey)
      call ujoinsr(cc,rdiv,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,temnow,cc,1,nsizey)
      call ujoinsr(cc,tt,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call transr1(jtrun,jtmax,nx,my,my_max,poly,plnow,pt,nsizey)
      call tranuv(jtrun,jtmax,nx,my,my_max,levp,onocos,wcfac,wdfac &
                 ,poly,dpoly,vornow,divnow,ut,vt,nsizey)
      call trngra(jtrun,jtmax,nx,my,my_max,cim,poly,dpoly,plnow   &
                 ,dlpl,dtpl,nsizey)
!
      call tendget (temold)
!
      do m=1,mlistnum
        mf=mlist(m)
        do n=mf,jtrun
          do k = 1, levp
            divold(k,1,n,m) = divten(k,1,n,m) 
            divold(k,2,n,m) = divten(k,2,n,m) 
            vorold(k,1,n,m) = vorten(k,1,n,m) 
            vorold(k,2,n,m) = vorten(k,2,n,m) 
          enddo
        enddo
      enddo
!
      return
      end
