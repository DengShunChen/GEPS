      subroutine tranrs1(jtrun,jtmax,nx,my,my_max,poly,w,r,s,nsize)
!
!  subroutine to transform a scalar grid point field to spectral
!  coefficients
!
! *** input ***
!
!  jtrun: zonal wavenumber truncation limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  poly: legendre polynomials
!  w: gaussian quadrature weights
!  r: 3-dim input grid pt. field to be transformed
!
! *** output ***
!
!  s: spectral coefficient fields
!
!  **********************************
!
      use const, only : RTYPE
      use index
      use paramt
      use fftcom
!
      implicit none

      integer jtrun,jtmax,nx,my,my_max,nsize
      integer mlx,myhalf,jj,j,nxj,i,jtrunj,m,mm,mp,mlst,mf
      integer l,i1,i2,i3,j1,j2

      real(kind=RTYPE)    poly(jtrun,my/2,jtmax),w(my)
      real(kind=RTYPE)    r(nx,my_max)
      real(kind=RTYPE)    s(jtrun,jtmax,2)
!
!      real(kind=RTYPE)    gwk1(nx+2,1,6,my_max)
      real(kind=RTYPE)    gwk1(nx+2,my_max)
!
      real(kind=RTYPE)    wcc_fk(jtmax,my_max*nsize,2)
      real(kind=RTYPE)    twcc_fk(jtmax*nsize,my_max,2)
      real(kind=RTYPE)    wss(jtrun,2)
      real(kind=RTYPE)    cc(nx+2,my_max)
!
      real(kind=RTYPE)    wccSUM(my,2)
      real(kind=RTYPE)    wccDIF(my,2)
      real(kind=RTYPE)    fj_polyw(my/2,jtrun)
      logical wfirst
      data wfirst/.true./
      save wfirst

!CWBinit
      twcc_fk=0.

!CWB2014
      gwk1=0.


!     mlx= (jtrun/2)*((jtrun+1)/2)
      myhalf=my/2
!
!  put grid point fields into two dimensional horizontal array
!
!*** r1 start ***
      do 23 jj =1, jlistnum
      j=jlist1(jj)
      nxj=nxdef(j)
!ocl repeat(im)
      do 23 i=1,nxj
      cc(i,jj)= r(i,jj)
   23 continue
!
!  fft for each guassian latitude of 2-d field
!
      if( lreduce.eq.0 ) then
#ifdef SP
      call rfftmlt_sp(cc,gwk1,trigs,ifax,1,nx+2,nx,jlistnum,-1)
#else
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,jlistnum,-1)
#endif
      else
!$omp  parallel do default(none)                         &
!$omp  private(jj,j,nxj,gwk1)                            &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx)  &
!$omp  schedule(dynamic)
      do jj=1,jlistnum
        j= jlist1(jj)
        nxj=nxdef(j)
#ifdef SP
        call rfftmlt_sp(cc(1,jj),gwk1(1,jj),trigsj(1,j),ifaxj(1,j), &
#else
        call rfftmlt(cc(1,jj),gwk1(1,jj),trigsj(1,j),ifaxj(1,j), &
#endif
                     1,nx+2,nxj,1,-1)
      enddo
!$omp end parallel do
      endif
!
      do j =1, jlistnum
         jj = jlist1(j)
         jtrunj = mtrundef(jj)
!ocl repeat(jtr)
      do m=1,jtrunj
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
         twcc_fk(mlst,j,1)=cc(mm,j)
         twcc_fk(mlst,j,2)=cc(mp,j)
      enddo
      enddo

!ch    call mpe_transpose_rs1(twcc_fk,wcc_fk,jtmax,my_max,2,nsize,col_comm)
       call mpe_transpose_rs1_sp(twcc_fk,wcc_fk,jtmax,my_max,2,nsize,col_comm)


!*** r1  end  ***
!ibm---beg

      do m=1,mlistnum
         mf=mlist(m)

      do L=mf,jtrun
        wss(L,1) = 0.
        wss(L,2) = 0.
      enddo

      do j = 1, myhalf
         do L = mf, jtrun
            fj_polyw(j,L) = poly(L,j,m)*w(j)
         end do
      end do

      i1=(jtrun-mf+1)/4
      i2=(jtrun-mf+1-i1*4)/2
      i3= jtrun-mf+1-i1*4-i2*2

      do j=1,myhalf
      j1=jlist2(j)
      j2=jlist2(my-j+1)
      if( mf.le.mtrundef(j) ) then
       wccSUM(j,1)=wcc_fk(m,j1,1)+wcc_fk(m,j2,1)
       wccDIF(j,1)=wcc_fk(m,j1,1)-wcc_fk(m,j2,1)
       wccSUM(j,2)=wcc_fk(m,j1,2)+wcc_fk(m,j2,2)
       wccDIF(j,2)=wcc_fk(m,j1,2)-wcc_fk(m,j2,2)
      endif
      enddo

      do i=1,i1
      L=mf+(i-1)*4
      do j=1,myhalf
      jj=my-j+1
       if( mf.le.mtrundef(j) ) then

        wss(L  ,1)=wss(L  ,1)+fj_polyw(j,L)*wccSUM(j,1)
        wss(L+1,1)=wss(L+1,1)+fj_polyw(j,L+1)*wccDIF(j,1)
        wss(L+2,1)=wss(L+2,1)+fj_polyw(j,L+2)*wccSUM(j,1)
        wss(L+3,1)=wss(L+3,1)+fj_polyw(j,L+3)*wccDIF(j,1)
        wss(L  ,2)=wss(L  ,2)+fj_polyw(j,L)*wccSUM(j,2)
        wss(L+1,2)=wss(L+1,2)+fj_polyw(j,L+1)*wccDIF(j,2)
        wss(L+2,2)=wss(L+2,2)+fj_polyw(j,L+2)*wccSUM(j,2)
        wss(L+3,2)=wss(L+3,2)+fj_polyw(j,L+3)*wccDIF(j,2)
        endif
      enddo
      enddo

      do i=1,i2
      L=mf+i1*4+(i-1)*2
      do j=1,myhalf
      jj=my-j+1
       if( mf.le.mtrundef(j) ) then
        wss(L  ,1)=wss(L  ,1)+fj_polyw(j,L)*wccSUM(j,1)
        wss(L+1,1)=wss(L+1,1)+fj_polyw(j,L+1)*wccDIF(j,1)
        wss(L  ,2)=wss(L  ,2)+fj_polyw(j,L)*wccSUM(j,2)
        wss(L+1,2)=wss(L+1,2)+fj_polyw(j,L+1)*wccDIF(j,2)
        endif
      enddo
      enddo

      do i=1,i3
      L=jtrun
      do j=1,myhalf
      jj=my-j+1
       if( mf.le.mtrundef(j) ) then
        wss(L  ,1)=wss(L  ,1)+fj_polyw(j,L)*wccSUM(j,1)
        wss(L  ,2)=wss(L  ,2)+fj_polyw(j,L)*wccSUM(j,2)
        endif
      enddo
      enddo

      do L=mf,jtrun
        s(L,m,1)=wss(L,1)
        s(L,m,2)=wss(L,2)
      enddo

      enddo

      return
      end
