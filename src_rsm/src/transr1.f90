      subroutine transr1(jtrun,jtmax,nx,my,my_max,poly,s,r,nsize)
!
!  subroutine to transform a spectral coefficient field to
!  grid point form
!
! *** input ***
!
!  jtrun: zonal wavenumber resolution limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  poly: legendre polynomials
!  s: spectral coefficient array to transform
!
! *** output ***
!
!  r: 3-d output grid point fields
!
!  **************************************
!
      use const, only : RTYPE
      use index
!     use paramt
      use fftcom
!
      implicit  none

      integer   jtrun,jtmax,nx,my,my_max,nsize

      integer   mlx,myhalf,m,mf,l,jlistnum_fj,j,l_fj,j_fj
      integer   llistnum_fj,jj,i,jtrunj,mm,mp,mlst,nxj

      real(kind=RTYPE)      r(nxp,my_max)
!
      real(kind=RTYPE)      s(jtrun,jtmax,2)
      real(kind=RTYPE)      gwk1(nx+2,my_max)
!
      real(kind=RTYPE)      wcc_fk(my_max*nsize,jtmax,2), twcc_fk(my_max,jtmax*nsize,2)
      real                  wss(jtrun,2)

      real(kind=RTYPE)      poly(jtrun,my/2,jtmax)
      real(kind=RTYPE)      cc(nx+2,my_max)
      real                  tcc(my,2)
!
      real                  ws2(jtrun,2)
!
      integer   jlist_fj(my/2)
      real                  fj_poly(jtrun,my/2)
      real                  fj_tcc(2,my)
      real                  fj_wss(2,jtrun)
      real                  fj_ws2(2,jtrun)

      mlx= (jtrun/2)*((jtrun+1)/2)
      myhalf=my/2

!CWB2014
      wcc_fk= 0.0

!
      do m=1,mlistnum
      mf=mlist(m)

      do l=1,jtrun
        wss(l,1) = s(l,m,1)
        wss(l,2) = s(l,m,2)
      enddo

!     do j=1,my_max*nsize
!     wcc_fk(j,m,1)= 0.0
!     wcc_fk(j,m,2)= 0.0
!     enddo

      do l=mf,jtrun,2
      ws2(l,1) = wss(l,1)
      ws2(l,2) = wss(l,2)
      enddo
      do l=mf+1,jtrun,2
      ws2(l,1) = -wss(l,1)
      ws2(l,2) = -wss(l,2)
      enddo

      tcc = 0.0
      fj_tcc = 0.0
      fj_poly = 0.0
      fj_wss = 0.0
      fj_ws2 = 0.0

      jlistnum_fj = 0
      do j=1,myhalf
         if( mf.le.mtrundef(j) ) then
            jlistnum_fj = jlistnum_fj + 1
            jlist_fj(jlistnum_fj) = j
         endif
      enddo

      do l = mf,jtrun
         l_fj = l - mf + 1
         do j_fj = 1,jlistnum_fj
            j = jlist_fj(j_fj)
            fj_poly(l_fj,j_fj) = poly(l,j,m)
         enddo
      enddo

      do l = mf,jtrun
         l_fj = l - mf + 1
         fj_wss(1,l_fj) = wss(l,1)
         fj_wss(2,l_fj) = wss(l,2)
      enddo

      do l = mf,jtrun
         l_fj = l - mf + 1
         fj_ws2(1,l_fj) = ws2(l,1)
         fj_ws2(2,l_fj) = ws2(l,2)
      enddo

      llistnum_fj = jtrun - mf + 1
      call dgemm('n','n',2,jlistnum_fj,llistnum_fj,1.0d+0,fj_wss,2 &
           ,fj_poly,jtrun,1.0d+0,fj_tcc,2)

      do j_fj = 1,jlistnum_fj
         j = jlist_fj(j_fj)
         tcc(j,1)=tcc(j,1)+fj_tcc(1,j_fj)
         tcc(j,2)=tcc(j,2)+fj_tcc(2,j_fj)
      enddo

      fj_tcc = 0.0

      llistnum_fj = jtrun - mf + 1
      call dgemm('n','n',2,jlistnum_fj,llistnum_fj,1.0d+0,fj_ws2,2 &
           ,fj_poly,jtrun,1.0d+0,fj_tcc,2)

      do j_fj = 1,jlistnum_fj
         j = jlist_fj(j_fj)
         jj=my-j+1
         tcc(jj,1)=tcc(jj,1)+fj_tcc(1,j_fj)
         tcc(jj,2)=tcc(jj,2)+fj_tcc(2,j_fj)
      enddo

      do j=1, my
       wcc_fk(jlist2(j),m,1) = tcc(j,1)
       wcc_fk(jlist2(j),m,2) = tcc(j,2)
      enddo

      enddo

!*** r1 start ***

      call mpe_transpose_rs1_sp(wcc_fk,twcc_fk,my_max,jtmax,2,nsize,col_comm)

      do jj =1,jlistnum
      do i=1,nx+2
      cc(i,jj)=0.
      enddo
      enddo

      do jj =1,jlistnum
        j= jlist1(jj)
        jtrunj = mtrundef(j)
!!ocl repeat(jtr)
      do m=1,jtrunj
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
         cc(mm,jj)=twcc_fk(jj,mlst,1)
         cc(mp,jj)=twcc_fk(jj,mlst,2)
      enddo
      enddo

!
      if( lreduce.eq.0 ) then
#ifdef SP
      call rfftmlt_sp(cc,gwk1,trigs,ifax,1,nx+2,nx,jlistnum,1)
#else
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,jlistnum,1)
#endif
      else
!$omp  parallel do default(none)                            &
!$omp  private(jj,j,nxj,gwk1)                               &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx    ) &
!$omp  schedule(dynamic)
      do jj=1,jlistnum
        j= jlist1(jj)
        nxj=nxdef(j)
#ifdef SP
        call rfftmlt_sp(cc(1,jj),gwk1(1,jj),trigsj(1,j),ifaxj(1,j),  &
#else
        call rfftmlt(cc(1,jj),gwk1(1,jj),trigsj(1,j),ifaxj(1,j),  &
#endif
                     1,nx+2,nxj,1,1)
      enddo
!$omp end parallel do
      endif
!
      do 22 jj =1,jlistnum
      j=jlist1(jj)
!      nxj=nxdef(j)
!      do 22 i=1,nxj
!      r(i,jj)= cc(i,jj)
      r(1:nxjlen(j),jj)= cc(nxjstart(j):nxjend(j),jj)
   22 continue

!*** r1  end  ***

      return
      end
