      subroutine rstrantq ( jtrun,jtmax,nx,my,my_max,lev,tadvv,tadvu  &
                  ,w,cim,onocos,poly,dpoly,temten,nsize )
!
      use index
      use paramt
      use fftcom
!
      implicit  none

      integer   jtrun,jtmax,nx,my,my_max,lev,nsize
      integer   myhalf,lev2,jj,j,nxj,k,i,m,mm,mp,mlst,mf,j2,j1,l
      integer   lchk,lle,jlistnum_fj,j_fj,l_fj,llistnum_fj

      real      poly(jtrun,my/2,jtmax),dpoly(jtrun,my/2,jtmax)   &
               ,cim(jtmax),onocos(my),w(my)
      real      tadvv(nxp ,levf,my_max),tadvu(nxp ,levf,my_max),dummy
      real      temten(lev,2,jtrun,jtmax)
!
      real      gwk1(nx+2,lev,2,my_max)
      real      wss(lev,2,1,jtrun)
      real      wcc_fk(lev,2,2,jtmax,my_max*nsize)
      real      twcc_fk(lev,2,2,jtmax*nsize,my_max)
      real      cc(nx+2,lev,2,my_max)
!
      real      wcc2(lev,2,1,my/2)
      real      wcc3(lev,2,1,my/2)
      real      wcc4(lev,2,1,my/2)
      real      wcc5(lev,2,1,my/2)

!
      real      wp(my/2,jtrun),wd(my/2,jtrun)

      real      fj_wcc2(lev*2*1,my/2)
      real      fj_wcc3(lev*2*1,my/2)
      real      fj_wcc4(lev*2*1,my/2)
      real      fj_wcc5(lev*2*1,my/2)
      real      fj_wd2(my/2,jtrun),fj_wp3(my/2,jtrun)
      real      fj_wd4(my/2,jtrun),fj_wp5(my/2,jtrun)
      real      fj_wss23(lev*2*1,jtrun)
      real      fj_wss45(lev*2*1,jtrun)
      integer   jlist_fj(my/2)
!CWBinit
      twcc_fk=0.
      cc=0.

!CWB2014
      gwk1=0.

      myhalf=my/2
      lev2=lev*2
!
!     do 23 jj =1, jlistnum
!       j= jlist1(jj)
!       nxj=nxdef(j)
!     do 23 k=1,lev
!!ocl novrec
!      do 23 i=1,nxj+3
!      cc(i,k,1,jj)= tadvv(i,k,jj)
!      cc(i,k,2,jj)= tadvu(i,k,jj)
!   23 continue

! replace loop 23 for 2dMPI >>
      call joinrs(cc,tadvv,tadvu,dummy,dummy,nx,my_max,levf,jlistnum,2,1)

      if( length_fft .eq. 0 .and. lreduce.eq.0 ) then
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*2,-1)
      else
!$omp  parallel do default(none)                            &
!$omp  private(jj,j,nxj,gwk1)                               &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx,lev) &
!$omp  schedule(dynamic)
      do jj=1,jlistnum
        j= jlist1(jj)
        nxj=nxdef(j)
        call rfftmlt(cc(1,1,1,jj),gwk1,trigsj(1,j),ifaxj(1,j), &
                     1,nx+2,nxj,lev*2,-1)
      enddo
!$omp end parallel do
      end if

      do j =1, jlistnum
      do m=1,jtrun
      do k=1,lev
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
         twcc_fk(k,1,1,mlst,j)=cc(mm,k,1,j)
         twcc_fk(k,2,1,mlst,j)=cc(mp,k,1,j)
         twcc_fk(k,1,2,mlst,j)=cc(mm,k,2,j)
         twcc_fk(k,2,2,mlst,j)=cc(mp,k,2,j)
      enddo
      enddo
      enddo
!
      call mpe_transpose_rs(twcc_fk,wcc_fk,lev*2*2,jtmax,my_max,nsize,col_comm)
!
      do m=1,mlistnum
         mf=mlist(m)

!ocl scalar
      do j=1,myhalf
      j1=jlist2(j); j2=jlist2(my-j+1)
      if( mf.le.mtrundef(j) ) then
!ocl vector
      do k=1,lev

      wcc2(k,1,1,j) = +(wcc_fk(k,1,1,m,j1)-wcc_fk(k,1,1,m,j2))
      wcc2(k,2,1,j) = +(wcc_fk(k,2,1,m,j1)-wcc_fk(k,2,1,m,j2))
      wcc3(k,1,1,j) = +(wcc_fk(k,2,2,m,j1)+wcc_fk(k,2,2,m,j2))
      wcc3(k,2,1,j) = -(wcc_fk(k,1,2,m,j1)+wcc_fk(k,1,2,m,j2))
      wcc4(k,1,1,j) = +(wcc_fk(k,1,1,m,j1)+wcc_fk(k,1,1,m,j2))
      wcc4(k,2,1,j) = +(wcc_fk(k,2,1,m,j1)+wcc_fk(k,2,1,m,j2))
      wcc5(k,1,1,j) = +(wcc_fk(k,2,2,m,j1)-wcc_fk(k,2,2,m,j2))
      wcc5(k,2,1,j) = -(wcc_fk(k,1,2,m,j1)-wcc_fk(k,1,2,m,j2))

      enddo
      endif
      enddo

      do l=mf,jtrun
      do j=1,myhalf
        wd(j,l)=w(j)*dpoly(l,j,m)
        wp(j,l)=w(j)*onocos(j)*cim(m)*poly(l,j,m)
      enddo
      enddo
!
      do L=mf,jtrun
      do k=1,lev*2
        wss(k,1,1,L) = 0
      enddo
      enddo

      fj_wcc2 = 0.0
      fj_wcc3 = 0.0
      fj_wcc4 = 0.0
      fj_wcc5 = 0.0
      fj_wd2  = 0.0
      fj_wp3  = 0.0
      fj_wd4  = 0.0
      fj_wp5  = 0.0
      fj_wss23 = 0.0
      fj_wss45 = 0.0

      lchk = iand(jtrun - mf + 1,1)
      lle = jtrun - lchk

      if(jtrun .gt. mf) then
         jlistnum_fj = 0
         do j = 1,myhalf
            if(mf.le.mtrundef(j)) then
               jlistnum_fj = jlistnum_fj + 1
               jlist_fj(jlistnum_fj) = j
            endif
         enddo

         do l = mf,lle,2
            l_fj = (l -mf + 2)/2
            do j_fj = 1,jlistnum_fj
               j = jlist_fj(j_fj)
               fj_wd2(j_fj,l_fj) = wd(j,l)
               fj_wp3(j_fj,l_fj) = wp(j,l)
               fj_wd4(j_fj,l_fj) = wd(j,l+1)
               fj_wp5(j_fj,l_fj) = wp(j,l+1)
            enddo
         enddo

         llistnum_fj = l_fj

         do j_fj = 1,jlistnum_fj
            j = jlist_fj(j_fj)
            do k = 1,lev*2
               fj_wcc2(k,j_fj) = wcc2(k,1,1,j)
               fj_wcc3(k,j_fj) = wcc3(k,1,1,j)
               fj_wcc4(k,j_fj) = wcc4(k,1,1,j)
               fj_wcc5(k,j_fj) = wcc5(k,1,1,j)
            enddo
         enddo

         call dgemm('n','n',lev*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc2,lev*2,fj_wd2,myhalf, &
              1.0d+0,fj_wss23,lev*2)

         call dgemm('n','n',lev*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc3,lev*2,fj_wp3,myhalf, &
              1.0d+0,fj_wss23,lev*2)

         do l = mf,lle,2
            l_fj = (l - mf + 2)/2
            do k = 1,lev*2
               wss(k,1,1,l) = wss(k,1,1,l) + fj_wss23(k,l_fj)
            enddo
         enddo

         call dgemm('n','n',lev*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc4,lev*2,fj_wd4,myhalf, &
              1.0d+0,fj_wss45,lev*2)

         call dgemm('n','n',lev*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc5,lev*2,fj_wp5,myhalf, &
              1.0d+0,fj_wss45,lev*2)

         do l = mf,lle,2
            l_fj = (l - mf + 2)/2
            do k = 1,lev*2
               wss(k,1,1,l+1) = wss(k,1,1,l+1) + fj_wss45(k,l_fj)
            enddo
         enddo

      endif

      if(lchk .eq. 1) then
         L = jtrun
         do j = 1,myhalf
            if(mf.le.mtrundef(j)) then
               do k = 1,lev*2
                  wss(k,1,1,L  )=wss(k,1,1,L  )+wcc2(k,1,1,j)*wd(j,L) &
                       +wcc3(k,1,1,j)*wp(j,L)
               enddo
            endif
         enddo
      endif
!

      do L=mf,jtrun
      do k=1,lev*2
        temten(k,1,L,m) = temten(k,1,L,m)+wss(k,1,1,L)
      enddo
      enddo

      enddo
!
      return
      end
!
! for 2dMPI -----------------------------------
!
      subroutine rstrantq_ncld( jtrun,jtmax,nx,my,my_max,lev,qadvv,qadvu  &
                  ,w,cim,onocos,poly,dpoly,temten,nsize )
!
      use param, only : ncld
      use index
      use paramt
      use fftcom
!
      implicit  none

      integer   jtrun,jtmax,nx,my,my_max,lev,nsize
      integer   myhalf,lev2,jj,j,nxj,k,i,m,mm,mp,mlst,mf,j2,j1,l
      integer   lchk,lle,jlistnum_fj,j_fj,l_fj,llistnum_fj

      real      poly(jtrun,my/2,jtmax),dpoly(jtrun,my/2,jtmax)   &
               ,cim(jtmax),onocos(my),w(my)
      real      qadvv(nxp ,levf*ncld,my_max),qadvu(nxp ,levf*ncld,my_max),dummy
      real      temten(lev,2,jtrun,jtmax)
!
      real      gwk1(nx+2,lev,2,my_max)
      real      wss(lev,2,1,jtrun)
      real      wcc_fk(lev,2,2,jtmax,my_max*nsize)
      real      twcc_fk(lev,2,2,jtmax*nsize,my_max)
      real      cc(nx+2,lev,2,my_max)
!
      real      wcc2(lev,2,1,my/2)
      real      wcc3(lev,2,1,my/2)
      real      wcc4(lev,2,1,my/2)
      real      wcc5(lev,2,1,my/2)

!
      real      wp(my/2,jtrun),wd(my/2,jtrun)

      real      fj_wcc2(lev*2*1,my/2)
      real      fj_wcc3(lev*2*1,my/2)
      real      fj_wcc4(lev*2*1,my/2)
      real      fj_wcc5(lev*2*1,my/2)
      real      fj_wd2(my/2,jtrun),fj_wp3(my/2,jtrun)
      real      fj_wd4(my/2,jtrun),fj_wp5(my/2,jtrun)
      real      fj_wss23(lev*2*1,jtrun)
      real      fj_wss45(lev*2*1,jtrun)
      integer   jlist_fj(my/2)

      real      bufA(nxp,levf*ncld, 2,my_max)
      real      bufB(nx ,levp*ncld, 2,my_max)

!CWBinit
      twcc_fk=0.
      cc=0.

!CWB2014
      gwk1=0.

      myhalf=my/2
      lev2=lev*2

      do jj =1, jlistnum
      do k=1,levf*ncld
      do i=1,nxp
         bufA(i,k,1,jj)=qadvv(i,k,jj)
         bufA(i,k,2,jj)=qadvu(i,k,jj)
      enddo
      enddo
      enddo

      call mpe2d_transpose_nxp_lev_ncld(bufA,bufB,nxp,nx,levf,levp,2,myf,my_max,jlistnum,jlen,nsizex,row_comm,ncld)

      do jj =1, jlistnum
      do k=1,lev  ! lev here is levp*ncld
      do i=1,nx
         cc(i,k,1,jj)= bufB(i,k,1,jj)
         cc(i,k,2,jj)= bufB(i,k,2,jj)
      enddo
      enddo
      enddo

      if( length_fft .eq. 0 .and. lreduce.eq.0 ) then
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*2,-1)
      else
!$omp  parallel do default(none)                            &
!$omp  private(jj,j,nxj,gwk1)                               &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx,lev) &
!$omp  schedule(dynamic)
      do jj=1,jlistnum
        j= jlist1(jj)
        nxj=nxdef(j)
        call rfftmlt(cc(1,1,1,jj),gwk1,trigsj(1,j),ifaxj(1,j), &
                     1,nx+2,nxj,lev*2,-1)
      enddo
!$omp end parallel do
      end if

      do j =1, jlistnum
      do m=1,jtrun
      do k=1,lev
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
         twcc_fk(k,1,1,mlst,j)=cc(mm,k,1,j)
         twcc_fk(k,2,1,mlst,j)=cc(mp,k,1,j)
         twcc_fk(k,1,2,mlst,j)=cc(mm,k,2,j)
         twcc_fk(k,2,2,mlst,j)=cc(mp,k,2,j)
      enddo
      enddo
      enddo
!
      call mpe_transpose_rs(twcc_fk,wcc_fk,lev*2*2,jtmax,my_max,nsize,col_comm)
!
      do m=1,mlistnum
         mf=mlist(m)

!ocl scalar
      do j=1,myhalf
      j1=jlist2(j); j2=jlist2(my-j+1)
      if( mf.le.mtrundef(j) ) then
!ocl vector
      do k=1,lev

      wcc2(k,1,1,j) = +(wcc_fk(k,1,1,m,j1)-wcc_fk(k,1,1,m,j2))
      wcc2(k,2,1,j) = +(wcc_fk(k,2,1,m,j1)-wcc_fk(k,2,1,m,j2))
      wcc3(k,1,1,j) = +(wcc_fk(k,2,2,m,j1)+wcc_fk(k,2,2,m,j2))
      wcc3(k,2,1,j) = -(wcc_fk(k,1,2,m,j1)+wcc_fk(k,1,2,m,j2))
      wcc4(k,1,1,j) = +(wcc_fk(k,1,1,m,j1)+wcc_fk(k,1,1,m,j2))
      wcc4(k,2,1,j) = +(wcc_fk(k,2,1,m,j1)+wcc_fk(k,2,1,m,j2))
      wcc5(k,1,1,j) = +(wcc_fk(k,2,2,m,j1)-wcc_fk(k,2,2,m,j2))
      wcc5(k,2,1,j) = -(wcc_fk(k,1,2,m,j1)-wcc_fk(k,1,2,m,j2))

      enddo
      endif
      enddo

      do l=mf,jtrun
      do j=1,myhalf
        wd(j,l)=w(j)*dpoly(l,j,m)
        wp(j,l)=w(j)*onocos(j)*cim(m)*poly(l,j,m)
      enddo
      enddo
!
      do L=mf,jtrun
      do k=1,lev*2
        wss(k,1,1,L) = 0
      enddo
      enddo

      fj_wcc2 = 0.0
      fj_wcc3 = 0.0
      fj_wcc4 = 0.0
      fj_wcc5 = 0.0
      fj_wd2  = 0.0
      fj_wp3  = 0.0
      fj_wd4  = 0.0
      fj_wp5  = 0.0
      fj_wss23 = 0.0
      fj_wss45 = 0.0

      lchk = iand(jtrun - mf + 1,1)
      lle = jtrun - lchk

      if(jtrun .gt. mf) then
         jlistnum_fj = 0
         do j = 1,myhalf
            if(mf.le.mtrundef(j)) then
               jlistnum_fj = jlistnum_fj + 1
               jlist_fj(jlistnum_fj) = j
            endif
         enddo

         do l = mf,lle,2
            l_fj = (l -mf + 2)/2
            do j_fj = 1,jlistnum_fj
               j = jlist_fj(j_fj)
               fj_wd2(j_fj,l_fj) = wd(j,l)
               fj_wp3(j_fj,l_fj) = wp(j,l)
               fj_wd4(j_fj,l_fj) = wd(j,l+1)
               fj_wp5(j_fj,l_fj) = wp(j,l+1)
            enddo
         enddo

         llistnum_fj = l_fj

         do j_fj = 1,jlistnum_fj
            j = jlist_fj(j_fj)
            do k = 1,lev*2
               fj_wcc2(k,j_fj) = wcc2(k,1,1,j)
               fj_wcc3(k,j_fj) = wcc3(k,1,1,j)
               fj_wcc4(k,j_fj) = wcc4(k,1,1,j)
               fj_wcc5(k,j_fj) = wcc5(k,1,1,j)
            enddo
         enddo

         call dgemm('n','n',lev*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc2,lev*2,fj_wd2,myhalf, &
              1.0d+0,fj_wss23,lev*2)

         call dgemm('n','n',lev*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc3,lev*2,fj_wp3,myhalf, &
              1.0d+0,fj_wss23,lev*2)

         do l = mf,lle,2
            l_fj = (l - mf + 2)/2
            do k = 1,lev*2
               wss(k,1,1,l) = wss(k,1,1,l) + fj_wss23(k,l_fj)
            enddo
         enddo

         call dgemm('n','n',lev*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc4,lev*2,fj_wd4,myhalf, &
              1.0d+0,fj_wss45,lev*2)

         call dgemm('n','n',lev*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc5,lev*2,fj_wp5,myhalf, &
              1.0d+0,fj_wss45,lev*2)

         do l = mf,lle,2
            l_fj = (l - mf + 2)/2
            do k = 1,lev*2
               wss(k,1,1,l+1) = wss(k,1,1,l+1) + fj_wss45(k,l_fj)
            enddo
         enddo

      endif

      if(lchk .eq. 1) then
         L = jtrun
         do j = 1,myhalf
            if(mf.le.mtrundef(j)) then
               do k = 1,lev*2
                  wss(k,1,1,L  )=wss(k,1,1,L  )+wcc2(k,1,1,j)*wd(j,L) &
                       +wcc3(k,1,1,j)*wp(j,L)
               enddo
            endif
         enddo
      endif
!

      do L=mf,jtrun
      do k=1,lev*2
        temten(k,1,L,m) = temten(k,1,L,m)+wss(k,1,1,L)
      enddo
      enddo

      enddo
!
      return
      end
