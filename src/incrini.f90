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
      use const, only: RTYPE
!
      implicit  none

      real hld1(nx,my)
      real(kind=RTYPE) pt1(nx,my_max)
      character typ*6
      real(kind=RTYPE) cc(nx+2,levp,1,my_max),dummy,                  &
                       hld3(nx,levp,my_max),hld4(nx,levp,ncld,my_max)
!!      real      cc(nx+2,levp,3+ncld,my_max),wss(levp,2,3+ncld,jtrun,jtmax)

      integer   k,ii,jj,j,nxj,lmax,itaup,lncrec,istat,i,KL,m,n,mf,ntrac
      real      fac
!
      lmax=16
!
      itaup = taup + 0.001
      lncrec = nx*my
!
      do k = 1, levp
        KL=lev-Llist(k)+1
        write (typ, '("m",i2.2,"200")' ) KL
        call syslbl_r (typ,idtg2,itaup,gmdef)
        call dmsread_split(nx,my,lncrec,'H',ifilin,hld1,istat)
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef(j)
          fac = cosl(j)/rad
          if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
          do i = 1, nxj
            hld3(i,k,jj) = hld1(i,j)*fac
          enddo
        enddo
      enddo
      call mpe2d_transpose_ndsl_f2p(hld3,ut, &
            nxp,nx,levf,levp,1,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!
      do k = 1, levp
        KL=lev-Llist(k)+1
        write (typ, '("m",i2.2,"210")' ) KL
        call syslbl_r (typ,idtg2,itaup,gmdef)
        call dmsread_split(nx,my,lncrec,'H',ifilin,hld1,istat)
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef(j)
          fac = cosl(j)/rad
          if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
          do i = 1,nxj
            hld3(i,k,jj) = hld1(i,j)*fac
          enddo
        enddo
      enddo
      call mpe2d_transpose_ndsl_f2p(hld3,vt, &
            nxp,nx,levf,levp,1,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!
      call syslbl_r ('B00010',idtg2,itaup,ggdef)
      call dmsread (nx,my,lncrec,'H',ifilin,hld1,istat)
      do jj = 1, jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
        do i = 1, nxj
          pt(i,jj) = hld1(ii,j) - ptop
          ii=ii+1
        enddo
!ch
        pt1(:,jj) = hld1(:,j) - ptop
      enddo
!
      do k = 1, levp
        KL=lev-Llist(k)+1
        write (typ, '("m",i2.2,"100")' ) KL
        call syslbl_r (typ,idtg2,itaup,gmdef)
        call dmsread_split(nx,my,lncrec,'H',ifilin,hld1,istat)
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef(j)
          if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
          do i = 1,nxj
            hld3(i,k,jj) = hld1(i,j)
          enddo
        enddo
      enddo
      call mpe2d_transpose_ndsl_f2p(hld3,tt, &
            nxp,nx,levf,levp,1,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!
      do k = 1, levp
        KL=lev-Llist(k)+1
        write (typ, '("m",i2.2,"500")' ) KL
        call syslbl_r (typ,idtg2,itaup,gmdef)
        call dmsread_split(nx,my,lncrec,'H',ifilin,hld1,istat)
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef(j)
          if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
          do i = 1,nxj
            hld4(i,k,1,jj) = hld1(i,j)
          enddo
        enddo
      enddo
!
      if( ncld .ge. 2 ) then
        ntrac=2
        do k = 1, levp
          KL=lev-Llist(k)+1
          write (typ, '("m",i2.2,"550")' ) KL     ! cloud liquid water content
          call syslbl_r (typ,idtg2,itaup,gmdef)
          call dmsread_split(nx,my,lncrec,'H',ifilin,hld1,istat)
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef(j)
            if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
            do i = 1,nxj
              hld4(i,k,ntrac,jj) = hld1(i,j)
            enddo
          enddo
        enddo
        if(ncld.ge.3)then
        ntrac=3
        do k = 1, levp
          KL=lev-Llist(k)+1
          write (typ, '("m",i2.2,"560")' ) KL
          call syslbl_r (typ,idtg2,itaup,gmdef)
          call dmsread_split(nx,my,lncrec,'H',ifilin,hld1,istat)
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef(j)
            if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
            do i = 1,nxj
              hld4(i,k,ntrac,jj) = hld1(i,j)
            enddo
          enddo
        enddo
        endif
      endif
      call mpe2d_transpose_ndsl_f2p(hld4,qt, &
            nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!
      do k = 1, levp
        KL=lev-Llist(k)+1
        write (typ, '("m",i2.2,"000")' ) KL
        call syslbl_r (typ,idtg2,itaup,gmdef)
        call dmsread_split(nx,my,lncrec,'H',ifilin,hld1,istat)
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef(j)
          if( lreduce.eq.1 )call reducepick (hld1(1,j),nxdef(j),nx,1)
          do i = 1,nxj
            hld3(i,k,jj) = hld1(i,j)*grav
          enddo
        enddo
      enddo
      call mpe2d_transpose_ndsl_f2p(hld3,phi, &
            nxp,nx,levf,levp,1,myf,my_max,jlistnum,jlen,nsizex,row_comm)
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
