      subroutine tranuv1(jtrun,jtmax,nx,my,my_max,lev,onocos,wcfac &
                        ,wdfac,poly,dpoly,vor,div,ut,vt,nsize)
!
!  subroutine to transform vorticity and divergence to velocity
!  components
!
! *** input ***
!
!  jtrun: zonal wavenumber truncation limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  lev: number of veritical levels to transform
!  onocos: 1.0/(cos(lat)**2)
!  wcfac: constants defined in cons
!  wdfac: constants defined in cons
!  poly: legendre polynomials
!  dpoly: d(poly)/d(sin(lat))
!  vor: spectral vorticity
!  div: spectral divergence
!
! *** output ***
!
!  ut: e-w velocity component
!  vt: n-s velocity component
!
!  ****************************************
!
      use const, only : RTYPE
      use index
!     use paramt
      use fftcom

      implicit  none

      integer   jtrun,jtmax,nx,my,my_max,lev,nsize
      integer   myhalf,lev2,mlx,j,m,mf,l,k,nb,jchk
      integer   jje,jlistnum_fj,l_fj,j_fj,llistnum_fj
      integer   kk,ll,jj,jx,j2,i,jtrunj,mchk,mm,mp,mlst
      integer   mm1,mp1,mlst1,mm2,mp2,mlst2,mm3,mp3,mlst3,nxj,ierr

      real      sa00,sa10,sa20,sa30,dummy
!
      real(kind=RTYPE) onocos(my),wcfac(jtrun,jtmax),wdfac(jtrun,jtmax)   &
               ,poly(jtrun,my/2,jtmax),dpoly(jtrun,my/2,jtmax)    &
               ,vor(jtrun,jtmax,2),div(jtrun,jtmax,2)             &
               ,ut(nxp,my_max),vt(nxp,my_max)
!
      real(kind=RTYPE)      gwk1(nx+2,lev,2,my_max)
      real(kind=RTYPE)      wcc_fk (lev,2,2,jtmax,my_max*nsize)
      real(kind=RTYPE)     twcc_fk(lev,2,2,jtmax*nsize,my_max)
      real(kind=RTYPE)     cc(nx+2,lev,2,my_max)
      real(kind=RTYPE)      tcc(lev,2,2,my)
      real(kind=RTYPE)      ws3(lev,2,2,jtrun)
      real(kind=RTYPE)      ws4(lev,2,2,jtrun)
!
      real(kind=RTYPE)      tc2(lev,2,2,my)
      real(kind=RTYPE)      wc(jtrun,my/2),wd(jtrun,my/2)

!CWB2015
!     real      coslr(jm)
!     save coslr
      real(kind=RTYPE), dimension(:), allocatable, save ::  coslr

!
      logical lfirst
      data lfirst/.true./
      save lfirst
!
      integer   jlist_fj(my/2)
      real      fj_ws3(lev*2*2,jtrun)
      real      fj_ws4(lev*2*2,jtrun)
      real      fj_tcc(lev*2*2,my)
      real      fj_wc(jtrun,my/2),fj_wd(jtrun,my/2)
      real      fj_tc2(lev*2*2,my)
!

!CWBinit
      wcc_fk=0.

      myhalf=my/2
      lev2=lev*2
      mlx= (jtrun/2)*((jtrun+1)/2)
!
      if (lfirst) then
!CWB2015 >>>
!       allocate(coslr(jm),stat=ierr)
        allocate(coslr(my),stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'tranuv : allocate fail '
               stop
           end if
!CWB2015 <<<
        do j=1,my
          coslr(j)=1./onocos(j)
        enddo
        lfirst=.false.
      endif

      do m=1,mlistnum
         mf=mlist(m)
!
        do j=1,myhalf
        if( mf.le.mtrundef(j) ) then
        do l=mf,jtrun
          wc(l,j)=wcfac(l,m)* poly(l,j,m)
          wd(l,j)=wdfac(l,m)*dpoly(l,j,m)*coslr(j)
        enddo
        endif
        enddo
!

        do l=mf,jtrun
        do k=1,lev
          ws3(k,1,1,l) = +div(l,m,2)
          ws3(k,2,1,l) = -div(l,m,1)
          ws3(k,1,2,l) = +vor(l,m,2)
          ws3(k,2,2,l) = -vor(l,m,1)
          ws4(k,1,1,l) = +vor(l,m,1)
          ws4(k,2,1,l) = +vor(l,m,2)
          ws4(k,1,2,l) = -div(l,m,1)
          ws4(k,2,2,l) = -div(l,m,2)
        enddo
        enddo

       do k=1,lev*2*2*myhalf
         tcc(k,1,1,1)=0.
       enddo

        nb=32
        jchk= iand(myhalf, 1)
        jje = myhalf-jchk

        fj_ws3 = 0.0
        fj_ws4 = 0.0
        fj_tcc = 0.0
        fj_wc  = 0.0
        fj_wd  = 0.0

        jlistnum_fj = 0
        do j = 1,jje
           if(mf .le. mtrundef(j)) then
              jlistnum_fj = jlistnum_fj + 1
              jlist_fj(jlistnum_fj) = j
           endif
        enddo

        do l = mf,jtrun
           l_fj = l - mf + 1
           do k = 1,lev2*2
              fj_ws3(k,l_fj) = ws3(k,1,1,l)
           enddo
        enddo

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do l = mf,jtrun
              l_fj = l - mf + 1
              fj_wc(l_fj,j_fj) = wc(l,j)
           enddo
        enddo

        llistnum_fj = jtrun - mf + 1
        call dgemm('n','n',lev*2*2,jlistnum_fj,llistnum_fj,1.0d+0,fj_ws3,lev*2*2,fj_wc, &
             jtrun,1.0d+0,fj_tcc,lev*2*2)


        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do k = 1,lev2*2
              tcc(k,1,1,j) = tcc(k,1,1,j) + fj_tcc(k,j_fj)
           enddo
        enddo


        do l = mf,jtrun
           l_fj = l - mf + 1
           do k = 1,lev2*2
              fj_ws4(k,l_fj) = ws4(k,1,1,l)
           enddo
        enddo

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do l = mf,jtrun
              l_fj = l - mf + 1
              fj_wd(l_fj,j_fj) = wd(l,j)
           enddo
        enddo

        fj_tcc = 0.0

        llistnum_fj = jtrun - mf + 1
        call dgemm('n','n',lev*2*2,jlistnum_fj,llistnum_fj,1.0d+0,fj_ws4,lev*2*2,fj_wd, &
             jtrun,1.0d+0,fj_tcc,lev*2*2)


        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do k = 1,lev2*2
              tcc(k,1,1,j) = tcc(k,1,1,j) + fj_tcc(k,j_fj)
           enddo
        enddo


!
! odd number
!
        if( jchk.eq.1 )then
          j = myhalf
        if( mf.le.mtrundef(j) ) then
        do kk=1,lev2*2,nb
        do ll=mf,jtrun,nb
          do k = kk,min(kk+nb-1,lev2*2),4
            sa00 = tcc(k  ,1,1,j  )
            sa10 = tcc(k+1,1,1,j  )
            sa20 = tcc(k+2,1,1,j  )
            sa30 = tcc(k+3,1,1,j  )
          do l = ll,min(ll+nb-1,jtrun)
            sa00 = sa00+ws3(k,1,1,l)*wc(l,j)+ws4(k,1,1,l)*wd(l,j)
            sa10 = sa10+ws3(k+1,1,1,l)*wc(l,j)+ws4(k+1,1,1,l)*wd(l,j)
            sa20 = sa20+ws3(k+2,1,1,l)*wc(l,j)+ws4(k+2,1,1,l)*wd(l,j)
            sa30 = sa30+ws3(k+3,1,1,l)*wc(l,j)+ws4(k+3,1,1,l)*wd(l,j)
          enddo
            tcc(k  ,1,1,j  ) = sa00
            tcc(k+1,1,1,j  ) = sa10
            tcc(k+2,1,1,j  ) = sa20
            tcc(k+3,1,1,j  ) = sa30
          enddo
        enddo
        enddo
        endif
        endif

        do l=mf,jtrun,2
        do k=1,lev
          ws3(k,1,1,l) = +div(l,m,2)
          ws3(k,2,1,l) = -div(l,m,1)
          ws3(k,1,2,l) = +vor(l,m,2)
          ws3(k,2,2,l) = -vor(l,m,1)
          ws4(k,1,1,l) = -vor(l,m,1)
          ws4(k,2,1,l) = -vor(l,m,2)
          ws4(k,1,2,l) = +div(l,m,1)
          ws4(k,2,2,l) = +div(l,m,2)
        enddo
        enddo
        do l=mf+1,jtrun,2
        do k=1,lev
          ws3(k,1,1,l) = -div(l,m,2)
          ws3(k,2,1,l) = +div(l,m,1)
          ws3(k,1,2,l) = -vor(l,m,2)
          ws3(k,2,2,l) = +vor(l,m,1)
          ws4(k,1,1,l) = +vor(l,m,1)
          ws4(k,2,1,l) = +vor(l,m,2)
          ws4(k,1,2,l) = -div(l,m,1)
          ws4(k,2,2,l) = -div(l,m,2)
        enddo
        enddo
!

        do k=1,lev*2*2*myhalf
         tc2(k,1,1,1)=0.
        enddo

!
        nb=32
        jchk= iand(myhalf, 1)

        fj_ws3 = 0.0
        fj_ws4 = 0.0
        fj_tc2 = 0.0
        fj_wc  = 0.0
        fj_wd  = 0.0

        do l = mf,jtrun
           l_fj = l - mf + 1
           do k = 1,lev2*2
              fj_ws3(k,l_fj) = ws3(k,1,1,l)
           enddo
        enddo

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do l = mf,jtrun
              l_fj = l - mf + 1
              fj_wc(l_fj,j_fj) = wc(l,j)
           enddo
        enddo

        llistnum_fj = jtrun - mf + 1
        call dgemm('n','n',lev*2*2,jlistnum_fj,llistnum_fj,1.0d+0,fj_ws3,lev*2*2,fj_wc, &
             jtrun,1.0d+0,fj_tc2,lev*2*2)

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do k = 1,lev2*2
              tc2(k,1,1,j) = tc2(k,1,1,j) + fj_tc2(k,j_fj)
           enddo
        enddo

        do l = mf,jtrun
           l_fj = l - mf + 1
           do k = 1,lev2*2
              fj_ws4(k,l_fj) = ws4(k,1,1,l)
           enddo
        enddo

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do l = mf,jtrun
              l_fj = l - mf + 1
              fj_wd(l_fj,j_fj) = wd(l,j)
           enddo
        enddo

        fj_tc2 = 0.0

        llistnum_fj = jtrun - mf + 1
        call dgemm('n','n',lev*2*2,jlistnum_fj,llistnum_fj,1.0d+0,fj_ws4,lev*2*2,fj_wd, &
             jtrun,1.0d+0,fj_tc2,lev*2*2)

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do k = 1,lev2*2
              tc2(k,1,1,j) = tc2(k,1,1,j) + fj_tc2(k,j_fj)
           enddo
        enddo

!
! odd number
!
        if( jchk.eq.1 )then
          j = myhalf
        if( mf.le.mtrundef(j) ) then
        do kk=1,lev2*2,nb
        do ll=mf,jtrun,nb
          do k = kk,min(kk+nb-1,lev2*2),4
            sa00 = tc2(k  ,1,1,j  )
            sa10 = tc2(k+1,1,1,j  )
            sa20 = tc2(k+2,1,1,j  )
            sa30 = tc2(k+3,1,1,j  )
          do l = ll,min(ll+nb-1,jtrun)
            sa00 = sa00+ws3(k,1,1,l)*wc(l,j)+ws4(k,1,1,l)*wd(l,j)
            sa10 = sa10+ws3(k+1,1,1,l)*wc(l,j)+ws4(k+1,1,1,l)*wd(l,j)
            sa20 = sa20+ws3(k+2,1,1,l)*wc(l,j)+ws4(k+2,1,1,l)*wd(l,j)
            sa30 = sa30+ws3(k+3,1,1,l)*wc(l,j)+ws4(k+3,1,1,l)*wd(l,j)
          enddo
            tc2(k  ,1,1,j  ) = sa00
            tc2(k+1,1,1,j  ) = sa10
            tc2(k+2,1,1,j  ) = sa20
            tc2(k+3,1,1,j  ) = sa30
          enddo
        enddo
        enddo
        endif
        endif

        do j=1,myhalf
          jj=jlist2(j)
          jx=my-j+1
          j2=jlist2(jx)
          do k=1,lev*2*2
            wcc_fk(k,1,1,m,jj)=tcc(k,1,1,j)
            wcc_fk(k,1,1,m,j2)=tc2(k,1,1,j)
          enddo
        enddo
!
      enddo   ! end of big m loop


!ch   call mpe_transpose_sr(wcc_fk,twcc_fk,lev*2*2,jtmax,my_max,nsize)
!     call mpe_transpose_sr(wcc_fk,twcc_fk,lev*2*2,jtmax,my_max,nsize,col_comm)
      call mpe_transpose_sr_sp(wcc_fk,twcc_fk,lev*2*2,jtmax,my_max,nsize,col_comm)

      do jj=1,jlistnum

      do i=1,(nx+2)*lev*2
       cc(i,1,1,jj)= 0.
      enddo
!
      j      = jlist1(jj)
      jtrunj = mtrundef(j)
      mchk=iand(jtrunj,3)

      do m=1,mchk
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
      do k=1,lev
         cc(mm,k,1,jj)=twcc_fk(k,1,1,mlst,jj)
         cc(mp,k,1,jj)=twcc_fk(k,2,1,mlst,jj)
         cc(mm,k,2,jj)=twcc_fk(k,1,2,mlst,jj)
         cc(mp,k,2,jj)=twcc_fk(k,2,2,mlst,jj)
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
      do k=1,lev
         cc(mm,k,1,jj)=twcc_fk(k,1,1,mlst,jj)
         cc(mp,k,1,jj)=twcc_fk(k,2,1,mlst,jj)
         cc(mm,k,2,jj)=twcc_fk(k,1,2,mlst,jj)
         cc(mp,k,2,jj)=twcc_fk(k,2,2,mlst,jj)
         cc(mm1,k,1,jj)=twcc_fk(k,1,1,mlst1,jj)
         cc(mp1,k,1,jj)=twcc_fk(k,2,1,mlst1,jj)
         cc(mm1,k,2,jj)=twcc_fk(k,1,2,mlst1,jj)
         cc(mp1,k,2,jj)=twcc_fk(k,2,2,mlst1,jj)
         cc(mm2,k,1,jj)=twcc_fk(k,1,1,mlst2,jj)
         cc(mp2,k,1,jj)=twcc_fk(k,2,1,mlst2,jj)
         cc(mm2,k,2,jj)=twcc_fk(k,1,2,mlst2,jj)
         cc(mp2,k,2,jj)=twcc_fk(k,2,2,mlst2,jj)
         cc(mm3,k,1,jj)=twcc_fk(k,1,1,mlst3,jj)
         cc(mp3,k,1,jj)=twcc_fk(k,2,1,mlst3,jj)
         cc(mm3,k,2,jj)=twcc_fk(k,1,2,mlst3,jj)
         cc(mp3,k,2,jj)=twcc_fk(k,2,2,mlst3,jj)
      enddo
      enddo

      enddo
!
      if( length_fft .eq. 0 .and. lreduce.eq.0 )then
#ifdef SP
      call rfftmlt_sp(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*2,1) ! CWB2015
#else
!     call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*num,1)
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*2,1) ! CWB2015
#endif
      else
!$omp  parallel do default(none)                            &
!$omp  private(jj,j,nxj,gwk1)                               &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx,lev) &
!$omp  schedule(dynamic)
      do jj = 1, jlistnum
        j= jlist1(jj)
        nxj=nxdef(j)
#ifdef SP
        call rfftmlt_sp(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j), &
                     1,nx+2,nxj,lev*2,1)
#else
        call rfftmlt(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j), &
                     1,nx+2,nxj,lev*2,1)
#endif
      end do
!$omp end parallel do
      end if
!
      do 22 jj =1,jlistnum
      j=jlist1(jj)
      ut(1:nxjlen(j),jj)= cc(nxjstart(j):nxjend(j),1,1,jj)
      vt(1:nxjlen(j),jj)= cc(nxjstart(j):nxjend(j),1,2,jj)
   22 continue
!
!     do 22 jj=1,jlistnum
!       j= jlist1(jj)
!       nxj=nxdef(j)
!     do 22 k=1,lev
!     do 22 i=1,nxj
!     ut(i,k,jj)= cc(i,k,1,jj)
!     vt(i,k,jj)= cc(i,k,2,jj)
!  22 continue

!2dMPI
!      call ujoinsr(cc,ut,vt,dummy,dummy,nx,my_max,levF,jlistnum,2,1)

      return
      end
