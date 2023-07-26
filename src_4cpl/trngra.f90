      subroutine trngra (jtrun,jtmax,nx,my,my_max,cim,poly,dpoly,s &
                        ,dlpl,dtpl,nsize)
!
!  subroutine to transform spectral terrain pressure to grid point
!  fields of zonal and meridional derivatives of terrain pressure
!
! *** input ***
!
!  jtrun: zonal wavenumber resolution limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  cim: zonal wavenumber array
!  poly: legendre polynomials
!  dpoly: d(poly)/d(sin(lat))
!  s: spectral coefficient array
!  nsiz: pe number
!
! **** output ****
!
!  dlpl: d(pt)/d(longitude)
!  dtpl: d(pt)/d(sin(lat))
!
! ****************************************************
!
      use const, only : RTYPE
      use index
!     use paramt
      use fftcom

      implicit  none
!
      integer   jtrun,jtmax,nx,my,my_max,nsize
      integer   myhalf,m,mf,l,j,jj,i,jtrunj,mm,mp,mlst,nxj

      real(kind=RTYPE) poly(jtrun,my/2,jtmax),dpoly(jtrun,my/2,jtmax)
      real(kind=RTYPE) s(jtrun,jtmax,2)
      real(kind=RTYPE) cim(jtmax)
!
      real(kind=RTYPE) dlpl(nxp,my_max),dtpl(nxp,my_max)
      real(kind=RTYPE) cc(nx+2,2,my_max)
!
      real(kind=RTYPE) gwk1(nx+2,2,my_max)
!
      real(kind=RTYPE) twcc_fk(my_max,jtmax*nsize,2)
      real(kind=RTYPE) twdd_fk(my_max,jtmax*nsize,2)

      real(kind=RTYPE) wcu_fk(my_max*nsize,jtmax,2),wcv_fk(my_max*nsize,jtmax,2)
      real(kind=RTYPE) wcu_t(my,2),wcv_t(my,2)

      real      ws3(jtrun,2,2)
      real      ws4(jtrun,2,2)

!
!CWBinit
      wcu_fk=0.
      wcv_fk=0.

      myhalf=my/2

      do m=1,mlistnum
         mf=mlist(m)

!     do 55 j=1,my
!     wcu_fk(j,m,1)= 0.0
!     wcu_fk(j,m,2)= 0.0
!     wcv_fk(j,m,1)= 0.0
!     wcv_fk(j,m,2)= 0.0
!  55 continue
!
      do l=mf,jtrun
      ws3(l,1,1) = +s(l,m,2)
      ws3(l,2,1) = -s(l,m,1)
      ws4(l,1,2) = -s(l,m,1)
      ws4(l,2,2) = -s(l,m,2)
      enddo

      do j=1,my*2
      wcu_t(j,1)=0.
      wcv_t(j,1)=0.
      enddo

      do l=mf,jtrun
!!ocl loop,repeat(jmhalf)
      do j=1,myhalf
      if ( mf.le.mtrundef(j) ) then
      wcu_t(j,1) = wcu_t(j,1) + ws3(l,1,1)*(cim(m)* poly(l,j,m))
      wcu_t(j,2) = wcu_t(j,2) + ws3(l,2,1)*(cim(m)* poly(l,j,m))
      wcv_t(j,1) = wcv_t(j,1) + ws4(l,1,2)*dpoly(l,j,m)
      wcv_t(j,2) = wcv_t(j,2) + ws4(l,2,2)*dpoly(l,j,m)
      endif
      enddo
      enddo

      do l=mf,jtrun,2
      ws3(l,1,1) = +s(l,m,2)
      ws3(l,2,1) = -s(l,m,1)
      ws4(l,1,2) = +s(l,m,1)
      ws4(l,2,2) = +s(l,m,2)
      enddo

      do l=mf+1,jtrun,2
      ws3(l,1,1) = -s(l,m,2)
      ws3(l,2,1) = +s(l,m,1)
      ws4(l,1,2) = -s(l,m,1)
      ws4(l,2,2) = -s(l,m,2)
      enddo

      do l=mf,jtrun
!!ocl loop,repeat(jmhalf)
      do jj=myhalf+1,my
      j=my-jj+1
      if ( mf.le.mtrundef(j) ) then
      wcu_t(jj,1) = wcu_t(jj,1) + ws3(l,1,1)*(cim(m)* poly(l,j,m))
      wcu_t(jj,2) = wcu_t(jj,2) + ws3(l,2,1)*(cim(m)* poly(l,j,m))
      wcv_t(jj,1) = wcv_t(jj,1) + ws4(l,1,2)*dpoly(l,j,m)
      wcv_t(jj,2) = wcv_t(jj,2) + ws4(l,2,2)*dpoly(l,j,m)
      endif
      enddo
      enddo

      do j=1,my
      jj=jlist2(j)
      wcu_fk(jj,m,1)=wcu_t(j,1)
      wcu_fk(jj,m,2)=wcu_t(j,2)
      wcv_fk(jj,m,1)=wcv_t(j,1)
      wcv_fk(jj,m,2)=wcv_t(j,2)
      enddo

      enddo

!      call mpe_transpose_rs1(wcu_fk,twcc_fk,my_max,jtmax,2,nsize)
!      call mpe_transpose_rs1(wcv_fk,twdd_fk,my_max,jtmax,2,nsize)
       call mpe_transpose_rs1_sp(wcu_fk,twcc_fk,my_max,jtmax,2,nsize,col_comm)
       call mpe_transpose_rs1_sp(wcv_fk,twdd_fk,my_max,jtmax,2,nsize,col_comm)

      do jj=1,jlistnum
      do i=1,nx+2
      cc(i,1,jj)= 0.
      cc(i,2,jj)= 0.
      enddo
      enddo

      do jj=1,jlistnum
        j      = jlist1(jj)
        jtrunj = mtrundef(j)
      do m=1,jtrunj
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
         cc(mm,1,jj)=twcc_fk(jj,mlst,1)
         cc(mp,1,jj)=twcc_fk(jj,mlst,2)
         cc(mm,2,jj)=twdd_fk(jj,mlst,1)
         cc(mp,2,jj)=twdd_fk(jj,mlst,2)
      enddo
      enddo

      if( lreduce.eq.0 ) then
#ifdef SP
      call rfftmlt_sp(cc,gwk1,trigs,ifax,1,nx+2,nx,jlistnum*2,1)
#else
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,jlistnum*2,1)
#endif
      else
!$omp  parallel do default(none)                          &
!$omp  private(jj,j,nxj,gwk1)                             &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx)
      do jj=1,jlistnum
        j= jlist1(jj)
        nxj=nxdef(j)
#ifdef SP
        call rfftmlt_sp(cc(1,1,jj),gwk1(1,1,jj),trigsj(1,j),ifaxj(1,j),1,nx+2,nxj,2,1)
#else
        call rfftmlt(cc(1,1,jj),gwk1(1,1,jj),trigsj(1,j),ifaxj(1,j),1,nx+2,nxj,2,1)
#endif
      enddo
!$omp end parallel do
      endif

      do 22 jj=1,jlistnum
      j=jlist1(jj)
!      nxj=nxdef(j)
!!ocl novrec
!      do 22 i=1,nxj
!      dlpl(i,jj)= -cc(i,1,jj)
!      dtpl(i,jj)= -cc(i,2,jj)
      dlpl(1:nxjlen(j),jj)= -cc(nxjstart(j):nxjend(j),1,jj)
      dtpl(1:nxjlen(j),jj)= -cc(nxjstart(j):nxjend(j),2,jj)
   22 continue
!
      return
      end
