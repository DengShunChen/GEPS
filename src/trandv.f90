      subroutine trandv (jtrun,jtmax,nx,my,my_max,lev,ut,vt,w,cim  &
                        ,onocos,poly,dpoly,vor,div,nsize)
!
!  subroutine to do grid point velocities to spectral vorticity
!  and divergence
!
! **** input ****
!
!  jtrun: zonal wavenumber truncation limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-w dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  lev: number of vertical levels to be transformed
!  ut: e-w velocity component
!  vt: n-s velocity component
!  w: gaussian quadrature weights
!  cim: zonal wavenumber array
!  onocos: 1.o/(cos(lat)**2)
!  poly: legendre polynomials
!  dpoly: d(poly)/d(sin(lat))
!
! *** output ***
!
!  vor: spectral vorticity
!  div: spectral divergence
!
! ********************************************
!
      use const, only : RTYPE
      use param, only : ncld
      use index
!     use paramt
      use fftcom

      implicit none
!
      integer   jtrun,jtmax,nx,my,my_max,lev,nsize
      integer   mlx,myhalf,levp2,jj,i,j,k,nxj,m,mm,mp,mlst,mf,j2,j1
      integer   l,lchk,lle,j_fj,i_fj,jlistnum_fj,l_fj,llistnum_fj

      real(kind=RTYPE)      poly(jtrun,my/2,jtmax),                  &
                            dpoly(jtrun,my/2,jtmax),cim(jtmax),      &
                            onocos(my),w(my)

      real(kind=RTYPE)      ut(nxp,lev,my_max),vt(nxp,lev,my_max) 
      real(kind=RTYPE)      vor(levp,2,jtrun,jtmax),div(levp,2,jtrun,jtmax)
!
      real(kind=RTYPE)      gwk1(nx+2,levp,2,my_max)
      real(kind=RTYPE)      wss (levp,2,2,jtrun)
      real(kind=RTYPE)      wcc_fk (levp,2,2,jtmax,my_max*nsize)
      real(kind=RTYPE)      twcc_fk(levp,2,2,jtmax*nsize,my_max)
      real(kind=RTYPE)      cc(nx+2,levp,2,my_max),dummy
!
      real(kind=RTYPE)      wcc2(levp,2,2,my/2)
      real(kind=RTYPE)      wcc3(levp,2,2,my/2)
      real(kind=RTYPE)      wcc4(levp,2,2,my/2)
      real(kind=RTYPE)      wcc5(levp,2,2,my/2)
!
      real(kind=RTYPE)      wp(my/2,jtrun),wd(my/2,jtrun)

      real      fj_wcc2(levp*2*2,my/2)
      real      fj_wcc3(levp*2*2,my/2)
      real      fj_wcc4(levp*2*2,my/2)
      real      fj_wcc5(levp*2*2,my/2)
      real      fj_wss23(levp*2*2,jtrun)
      real      fj_wss45(levp*2*2,jtrun)
      real      fj_wd2(my/2,jtrun),fj_wp3(my/2,jtrun)
      real      fj_wd4(my/2,jtrun),fj_wp5(my/2,jtrun)
      integer   jlist_fj(my/2)
!
!CWBinit
      twcc_fk=0.

!CWB2014
      gwk1=0.

      mlx= (jtrun/2)*((jtrun+1)/2)
      myhalf=my/2
      levp2=levp*2
!
!     do 23 jj =1, jlistnum
!       j= jlist1(jj)
!       nxj=nxdef(j)
!     do 23 k=1,lev
!     do 23 i=1,nxj
!     cc(i,k,1,jj)= ut(i,k,jj)
!     cc(i,k,2,jj)= vt(i,k,jj)
!  23 continue

!2dMPI
      call joinrs(cc,ut,vt,dummy,dummy,nx,my_max,lev,jlistnum,2,1)

      if( length_fft .eq. 0 .and. lreduce.eq.0 )then
#ifdef SP
      call rfftmlt_sp(cc,gwk1,trigs,ifax,1,nx+2,nx,levp*jlistnum*2,-1)
#else
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,levp*jlistnum*2,-1)
#endif
      else
!$omp  parallel do default(none)  &
!$omp  private(jj,j,nxj,gwk1)     &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx,levp) &
!$omp  schedule(dynamic)
      do jj=1,jlistnum
        j= jlist1(jj)
        nxj=nxdef(j)
#ifdef SP
        call rfftmlt_sp(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j),  &
#else
        call rfftmlt(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j),  &
#endif
                     1,nx+2,nxj,levp*2,-1)
      enddo
!$omp end parallel do
      end if
!
      do jj =1, jlistnum
      do m=1,jtrun
      do k=1,levp
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
         twcc_fk(k,1,1,mlst,jj)=cc(mm,k,1,jj)
         twcc_fk(k,2,1,mlst,jj)=cc(mp,k,1,jj)
         twcc_fk(k,1,2,mlst,jj)=cc(mm,k,2,jj)
         twcc_fk(k,2,2,mlst,jj)=cc(mp,k,2,jj)
      enddo
      enddo
      enddo

      call mpe_transpose_rs_sp(twcc_fk,wcc_fk,levp*2*2,jtmax,my_max,nsize,col_comm)

      do m=1,mlistnum
         mf=mlist(m)

      do j=1,myhalf
      j1=jlist2(j); j2=jlist2(my-j+1)

      if( mf.le.mtrundef(j) ) then
      do k=1,levp
      wcc2(k,1,1,j) = +(wcc_fk(k,1,1,m,j1)-wcc_fk(k,1,1,m,j2))
      wcc2(k,2,1,j) = +(wcc_fk(k,2,1,m,j1)-wcc_fk(k,2,1,m,j2))
      wcc2(k,1,2,j) = -(wcc_fk(k,1,2,m,j1)-wcc_fk(k,1,2,m,j2))
      wcc2(k,2,2,j) = -(wcc_fk(k,2,2,m,j1)-wcc_fk(k,2,2,m,j2))
      wcc3(k,1,1,j) = -(wcc_fk(k,2,2,m,j1)+wcc_fk(k,2,2,m,j2))
      wcc3(k,2,1,j) = +(wcc_fk(k,1,2,m,j1)+wcc_fk(k,1,2,m,j2))
      wcc3(k,1,2,j) = -(wcc_fk(k,2,1,m,j1)+wcc_fk(k,2,1,m,j2))
      wcc3(k,2,2,j) = +(wcc_fk(k,1,1,m,j1)+wcc_fk(k,1,1,m,j2))
      wcc4(k,1,1,j) = +(wcc_fk(k,1,1,m,j1)+wcc_fk(k,1,1,m,j2))
      wcc4(k,2,1,j) = +(wcc_fk(k,2,1,m,j1)+wcc_fk(k,2,1,m,j2))
      wcc4(k,1,2,j) = -(wcc_fk(k,1,2,m,j1)+wcc_fk(k,1,2,m,j2))
      wcc4(k,2,2,j) = -(wcc_fk(k,2,2,m,j1)+wcc_fk(k,2,2,m,j2))
      wcc5(k,1,1,j) = -(wcc_fk(k,2,2,m,j1)-wcc_fk(k,2,2,m,j2))
      wcc5(k,2,1,j) = +(wcc_fk(k,1,2,m,j1)-wcc_fk(k,1,2,m,j2))
      wcc5(k,1,2,j) = -(wcc_fk(k,2,1,m,j1)-wcc_fk(k,2,1,m,j2))
      wcc5(k,2,2,j) = +(wcc_fk(k,1,1,m,j1)-wcc_fk(k,1,1,m,j2))
      enddo
      endif

      enddo

!----
      do j=1,myhalf
      do l=mf,jtrun
        wd(j,l)=w(j)*dpoly(l,j,m)
        wp(j,l)=w(j)*onocos(j)*cim(m)*poly(l,j,m)
      enddo
      enddo

      do l=mf,jtrun
      do k=1,levp2*2
        wss(k,1,1,l) = 0.
      enddo
      enddo

      lchk = iand(jtrun - mf + 1,1)
      lle = jtrun - lchk

      if(jtrun .gt. mf) then

         fj_wss23 = 0.0
         fj_wss45 = 0.0
         fj_wcc2  = 0.0
         fj_wcc3  = 0.0
         fj_wcc4  = 0.0
         fj_wcc5  = 0.0
         fj_wd2   = 0.0
         fj_wp3   = 0.0
         fj_wd4   = 0.0
         fj_wp5   = 0.0

         jlistnum_fj = 0
         do j = 1,myhalf
            if(mf .le. mtrundef(j)) then
               jlistnum_fj = jlistnum_fj + 1
               jlist_fj(jlistnum_fj) = j
            endif
         enddo

         do j_fj = 1,jlistnum_fj
            j = jlist_fj(j_fj)
            do k=1,levp*2*2
               fj_wcc2(k,j_fj) = wcc2(k,1,1,j)
               fj_wcc3(k,j_fj) = wcc3(k,1,1,j)
               fj_wcc4(k,j_fj) = wcc4(k,1,1,j)
               fj_wcc5(k,j_fj) = wcc5(k,1,1,j)
            enddo
         enddo

         do l = mf,lle,2
            l_fj = (l - mf + 2)/2
            do j_fj = 1,jlistnum_fj
               j = jlist_fj(j_fj)
               fj_wd2(j_fj,l_fj) = wd(j,l)
               fj_wp3(j_fj,l_fj) = wp(j,l)
               fj_wd4(j_fj,l_fj) = wd(j,l+1)
               fj_wp5(j_fj,l_fj) = wp(j,l+1)
            enddo
         enddo

         llistnum_fj = l_fj

         call dgemm('n','n',levp*2*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc2,levp*2*2,&
              fj_wd2,myhalf,1.0d+0,fj_wss23,levp*2*2)

         call dgemm('n','n',levp*2*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc3,levp*2*2,&
              fj_wp3,myhalf,1.0d+0,fj_wss23,levp*2*2)

         do l = mf,lle,2
            l_fj = (l - mf + 2)/2
            do k = 1,levp*2*2
               wss(k,1,1,l) = wss(k,1,1,l) + fj_wss23(k,l_fj)
            enddo
         enddo

         call dgemm('n','n',levp*2*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc4,levp*2*2,&
              fj_wd4,myhalf,1.0d+0,fj_wss45,levp*2*2)

         call dgemm('n','n',levp*2*2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wcc5,levp*2*2,&
              fj_wp5,myhalf,1.0d+0,fj_wss45,levp*2*2)


         do l = mf,lle,2
            l_fj = (l - mf + 2)/2
            do k = 1,levp*2*2
               wss(k,1,1,l+1) = wss(k,1,1,l+1) + fj_wss45(k,l_fj)
            enddo
         enddo

      endif

      if(lchk .eq. 1) then
         l = jtrun
         do j=1,myhalf
            if(mf .le. mtrundef(j)) then
               do k=1,levp*2*2
                  wss(k,1,1,l  )= wss(k,1,1,l  )+wcc2(k,1,1,j)*wd(j,l) &
                                                +wcc3(k,1,1,j)*wp(j,l)
               enddo
            endif
         enddo
      endif

      do l=mf,jtrun
      do k=1,levp*2
      vor(k,1,l,m) = wss(k,1,1,l)
      div(k,1,l,m) = wss(k,1,2,l)
      enddo
      enddo

      enddo
!
      return
      end
