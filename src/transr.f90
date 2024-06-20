      subroutine transr (jtrun,jtmax,nx,my,my_max,lev,poly,wss  &
                        ,cc,num,nsize)
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
!  lev: number of levels to transform
!  poly: legendre polynomials
!  wss: spectral coefficient array to transform
!  num: number of variables grouped together
!
! *** output ***
!
!  cc_r8: 3-d output grid point fields
!
!  **************************************
!
      use const, only : RTYPE
      use index
!     use paramt
      use fftcom
!
      implicit  none

      integer   jtrun,jtmax,nx,my,my_max,lev,num,nsize
      integer   mlx,myhalf,lev2,m,mf,lmax,lchk,j,k,l
      integer   nb,jchk,jje,jlistnum_fj,j_fj,l_fj,llistnum_fj
      integer   kk,ll,jj,jx,j2,ii,i,jtrunj,mchk,mm,mp,mlst
      integer   mm1,mp1,mlst1,mm2,mp2,mlst2,mm3,mp3,mlst3,nxj

      real(kind=RTYPE)      sa00,sa10,sb00,sb10

      real(kind=RTYPE)      poly(jtrun,my/2,jtmax)
      real(kind=RTYPE)      cc(nx+2,lev,num,my_max),wss(lev,2,num,jtrun,jtmax)
!
      real(kind=RTYPE)      gwk1(nx+2,lev,num,my_max)
!
      real(kind=RTYPE)      wcc_fk (lev,2,num,jtmax,my_max*nsize)
      real(kind=RTYPE)      twcc_fk(lev,2,num,jtmax*nsize,my_max)
      real(kind=RTYPE)      tcc(lev,2,num,my/2),tc2(lev,2,num,my/2)
      real      ws2(lev,2,num,jtrun)

      real      fj_poly(jtrun,my/2)
      real      fj_tcc(lev*2*num,my/2)
      real      fj_tc2(lev*2*num,my/2)
      real      fj_ws2(lev*2*num,jtrun)
      real      fj_wss(lev*2*num,jtrun)
      integer   jlist_fj(my/2)

!CWBinit
      cc=0.
      wcc_fk=0.

      mlx= (jtrun/2)*((jtrun+1)/2)
      myhalf=my/2
      lev2=lev*2
!
!-- do start
!
      do m=1,mlistnum
        mf=mlist(m)
        lmax=jtrun-mf+1
        lchk=iand(lmax,1)

      do l=mf,jtrun-1,2
      do k=1,lev2*num
      ws2(k,1,1,l) = wss(k,1,1,l,m)
      ws2(k,1,1,l+1) = -wss(k,1,1,l+1,m)
      enddo
      enddo
      if (lchk.eq.1) then
        l=jtrun
        do k=1,lev2*num
          ws2(k,1,1,l) = wss(k,1,1,l,m)
        enddo
      endif
!
        do k=1,lev2*num*myhalf
         tcc(k,1,1,1)=0.
         tc2(k,1,1,1)=0.
        enddo
!
!  nb : multiple of 2
!
        nb=32
        jchk= iand(myhalf, 1)
        jje = myhalf-jchk
!
        fj_tcc = 0.0
        fj_tc2 = 0.0
        fj_wss = 0.0
        fj_ws2 = 0.0
        fj_poly = 0.0

        jlistnum_fj = 0
        do j = 1,jje
           if( mf.le.mtrundef(j) ) then
              jlistnum_fj = jlistnum_fj + 1
              jlist_fj(jlistnum_fj) = j
           endif
        enddo

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do l = mf,jtrun
              l_fj = l - mf + 1
              fj_poly(l_fj,j_fj) = poly(l,j,m)
           enddo
        enddo

        do l = mf,jtrun
           l_fj = l - mf + 1
           do k = 1,lev2*num
              fj_wss(k,l_fj) = wss(k,1,1,l,m)
           enddo
        enddo

        do l = mf,jtrun
           l_fj = l - mf + 1
           do k = 1,lev2*num
              fj_ws2(k,l_fj) = ws2(k,1,1,l)
           enddo
        enddo

        llistnum_fj = jtrun - mf + 1

        call dgemm('n','n',lev2*num,jlistnum_fj,llistnum_fj   &
             ,1.0d+0,fj_wss,lev2*num,fj_poly,jtrun,1.0d+0,fj_tcc,lev2*num)

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do k = 1,lev2*num
              tcc(k  ,1,1,j) = tcc(k  ,1,1,j)+ fj_tcc(k,j_fj)
           enddo
        enddo

        call dgemm('n','n',lev2*num,jlistnum_fj,llistnum_fj   &
             ,1.0d+0,fj_ws2,lev2*num,fj_poly,jtrun,1.0d+0,fj_tc2,lev2*num)

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do k = 1,lev2*num
              tc2(k  ,1,1,j) = tc2(k  ,1,1,j)+ fj_tc2(k,j_fj)
           enddo
        enddo

!
! odd number
!
        if( jchk.eq.1 )then
          j = myhalf
        if( mf.le.mtrundef(j) ) then
        do kk=1,lev2*num,nb
        do ll=mf,jtrun,nb
          do k = kk,min(kk+nb-1,lev2*num),2
            sa00=tcc(k,1,1,j)
            sa10=tcc(k+1,1,1,j)
            sb00=tc2(k,1,1,j)
            sb10=tc2(k+1,1,1,j)
            do l = ll,min(ll+nb-1,jtrun)
              sa00=sa00+poly(l,j,  m)*wss(k,1,1,l,m)
              sa10=sa10+poly(l,j,  m)*wss(k+1,1,1,l,m)
              sb00=sb00+poly(l,j,  m)*ws2(k,1,1,l)
              sb10=sb10+poly(l,j,  m)*ws2(k+1,1,1,l)
            enddo
            tcc(k,1,1,j  )=sa00
            tcc(k+1,1,1,j  )=sa10
            tc2(k,1,1,j)=sb00
            tc2(k+1,1,1,j)=sb10
          enddo
        enddo
        enddo
        endif
        endif

      do j=1, myhalf
      jj=jlist2(j)
      jx=my-j+1
      j2=jlist2(jx)
      do k=1,lev2*num
       wcc_fk(k,1,1,m,jj)=tcc(k,1,1,j)
       wcc_fk(k,1,1,m,j2)=tc2(k,1,1,j)
      enddo
      enddo

      enddo

      call mpe_transpose_sr_sp(wcc_fk,twcc_fk,lev*2*num,jtmax,my_max,nsize,col_comm)
!      call mpe_transpose_sr(wcc_fk,twcc_fk,lev*2*num,jtmax,my_max,nsize,col_comm)

      do jj =1,jlistnum

      do ii=1,num
      do k=1,lev
      do i=1,nx+2
       cc(i,k,ii,jj)=0.
      enddo
      enddo
      enddo

      j = jlist1(jj)
      jtrunj = mtrundef(j)
      mchk=iand(jtrunj,3)

      do m=1,mchk
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
      do ii=1,num
      do k=1,lev
         cc(mm,k,ii,jj)=twcc_fk(k,1,ii,mlst,jj)
         cc(mp,k,ii,jj)=twcc_fk(k,2,ii,mlst,jj)
      enddo
      enddo
      enddo

      do m=mchk+1,jtrunj,4
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
         mm1= 2*(m+1)-1
         mp1= mm1+1
         mlst1=nlist(m+1)
         mm2= 2*(m+2)-1
         mp2= mm2+1
         mlst2=nlist(m+2)
         mm3= 2*(m+3)-1
         mp3= mm3+1
         mlst3=nlist(m+3)
      do ii=1,num
      do k=1,lev
         cc(mm,k,ii,jj)=twcc_fk(k,1,ii,mlst,jj)
         cc(mp,k,ii,jj)=twcc_fk(k,2,ii,mlst,jj)
         cc(mm1,k,ii,jj)=twcc_fk(k,1,ii,mlst1,jj)
         cc(mp1,k,ii,jj)=twcc_fk(k,2,ii,mlst1,jj)
         cc(mm2,k,ii,jj)=twcc_fk(k,1,ii,mlst2,jj)
         cc(mp2,k,ii,jj)=twcc_fk(k,2,ii,mlst2,jj)
         cc(mm3,k,ii,jj)=twcc_fk(k,1,ii,mlst3,jj)
         cc(mp3,k,ii,jj)=twcc_fk(k,2,ii,mlst3,jj)
      enddo
      enddo
      enddo

      enddo

!
      if( length_fft .eq. 0 .and. lreduce.eq.0 )then
#ifdef SP
      call rfftmlt_sp(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*num,1)
#else
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*num,1)
#endif
      else
!$omp  parallel do default(none)                                      &
!$omp  private(jj,j,nxj,gwk1)                                         &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx,lev,num)       &
!$omp  schedule(dynamic)
         do jj = 1, jlistnum
            j= jlist1(jj)
            nxj=nxdef(j)
#ifdef SP
            call rfftmlt_sp(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j), &
#else
            call rfftmlt(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j), &
#endif
                 1,nx+2,nxj,lev*num,1)
         end do
!$omp end parallel do
      end if
!
   20 continue

!CWB2021

      return
      end
