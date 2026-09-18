      subroutine tranrs (jtrun,jtmax,nx,my,my_max,lev,poly,w,cc  &
                        ,wss,num,nsize)
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
!  lev: number of vertical levels to transform
!  poly: legendre polynomials
!  w: gaussian quadrature weights
!  cc_r8: 3-dim input grid pt. field to be transformed
!  num: number of variables grouped together
!
! *** output ***
!
!  wss: spectral coefficient fields
!
!  **********************************
!
      use const, only : RTYPE
      use index
!     use paramt
      use fftcom
      use fj_pad
!
      ! #region agent log
      use rank, only : myrank, MPI_COMM_gfs
      use mpi
      ! #endregion
      implicit none
      ! #region agent log
      ! Clean comparison points against the GPU: twcc_fk and wss are both
      ! explicitly zeroed here (L57/L58) and on the GPU side twcc_fk is
      ! cudaMemsetAsync'd while wss is fully overwritten by dgemm beta=0.0,
      ! so neither carries uninitialised padding. `cc` deliberately NOT
      ! measured: its whole-array sum is NaN on CPU (handoff 7.1 rule 2).
      real(kind=8) :: dbgtw, dbgwf, dbgws
      ! cc is cc(nx+2, lev, num, my_max) but only i<=nxj / jj<=jlistnum are
      ! ever written, so a whole-array sum reads uninitialised padding and
      ! comes back NaN. Sum the USED region only and count non-finite values
      ! separately, so a stray NaN shows up instead of poisoning the total.
      real(kind=8) :: dbgcb, dbgca, dbgv
      integer :: dbi, dbk, dbii, dbjj, dbj, dbnxj, dbnb, dbna
      ! Post-transpose quantities are split by WAVENUMBER across ranks, and
      ! spectral energy is concentrated at low m, so "rank 0 x nsize" is not
      ! a valid global estimate for them (handoff, 2026-09-15). Instead:
      !   - per-m sums (wcc_fk over all latitudes, wss over l>=mf only, the
      !     region the GPU dgemm actually writes) gathered to rank 0 and
      !     printed by mf, so the GPU's single rank can be matched 1:1;
      !   - ALLREDUCE'd global totals.
      ! dbg_call counts collective calls so the analyser can pick the record
      ! belonging to a DBGMARK (tranrs is shared by getrdy/intgrt/incrini).
      ! MPI_REAL8, NOT MPI_DOUBLE_PRECISION: this OpenMPI was configured with
      ! -fdefault-real-8 alone, so its DOUBLE PRECISION is 16 bytes
      ! (OMPI_SIZEOF_FORTRAN_DOUBLE_PRECISION 16) while the model's is 8.
      ! With MPI_DOUBLE_PRECISION every call moved 2x the bytes: half the
      ! gathered rows were garbage and the receive overrun corrupted the heap
      ! ("double free or corruption" at the deallocate). Measured 2026-09-15.
      real(kind=8), allocatable :: dbg_loc(:,:), dbg_all(:,:)
      real(kind=8) :: dbg_g(3), dbg_gl(3)
      integer :: dbg_m, dbg_ierr, dbg_np
      integer, save :: dbg_call = 0
      ! #endregion

      integer jtrun,jtmax,nx,my,my_max,lev,num,nsize
      integer mlx,myhalf,lev2,nxj,j,jj,mchk,jtrunj,m,mm,mp
      integer mlst,ii,k,mm1,mp1,mlst1,mm2,mp2,mlst2,mm3,mp3,mlst3,mf
      integer lchk,lle,jlistnum_fj,j_fj,j1,j2,l,l_fj,llistnum_fj

      real(kind=RTYPE)    poly(jtrun,my/2,jtmax),w(my)
      real(kind=RTYPE)    wss(lev,2,num,jtrun,jtmax)
!
      real(kind=RTYPE)    gwk1(nx+2,lev,num,my_max)
!
      real(kind=RTYPE)    wcc_fk (lev,2,num,jtmax,my_max*nsize)
      real(kind=RTYPE)    twcc_fk(lev,2,num,jtmax*nsize,my_max)
      real(kind=RTYPE)    cc(nx+2,lev,num,my_max)
!
      real    fj_wss_sum(lev*2*num,jtrun)
      real    fj_wss_dif(lev*2*num,jtrun)
      real    fj_polyw_sum(my/2+npad,jtrun)
      real    fj_polyw_dif(my/2+npad,jtrun)
      real    fj_wccSUM(lev*2*num,my/2)
      real    fj_wccDIF(lev*2*num,my/2)
      integer jlist_fj(my/2)
!
!CWBinit
      twcc_fk=0.
      wss=0.
      ! #region agent log
      dbg_call = dbg_call + 1
      if (myrank .eq. 0) print *,'DBGTRC call=',dbg_call
      ! #endregion

!CWB2014
      gwk1=0.

!CWB2021 for single precision test

      mlx= (jtrun/2)*((jtrun+1)/2)
      myhalf=my/2
      lev2=lev*2
!
!  fft for each guassian latitude of 2-d field
!
      ! #region agent log
      dbgcb = 0.0d0
      dbnb = 0
      do dbjj = 1, jlistnum
         dbj = jlist1(dbjj)
         dbnxj = nxdef(dbj)
         do dbii = 1, num
            do dbk = 1, lev
               do dbi = 1, dbnxj
                  dbgv = real(cc(dbi, dbk, dbii, dbjj), kind=8)
                  if (dbgv /= dbgv) then
                     dbnb = dbnb + 1
                  else
                     dbgcb = dbgcb + dbgv*dbgv
                  end if
               end do
            end do
         end do
      end do
      if (myrank .eq. 0) print *,'DBGCB cc_before=',dbgcb,' nonfinite=',dbnb
      ! #endregion
      if( length_fft .eq. 0 .and. lreduce.eq.0 )then
#ifdef SP
      call rfftmlt_sp(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*num,-1)
#else
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*num,-1)
#endif
      else
!$omp  parallel do default(none)                                &
!$omp  private(jj,j,nxj,gwk1)                                   &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx,lev,num) &
!$omp  schedule(dynamic)
      do jj=1,jlistnum
        j= jlist1(jj)
        nxj=nxdef(j)
#ifdef SP
        call rfftmlt_sp(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j), &
#else
        call rfftmlt(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j),  &
#endif
                     1,nx+2,nxj,lev*num,-1)
      end do
!$omp end parallel do
      end if
      ! #region agent log
      dbgca = 0.0d0
      dbna = 0
      do dbjj = 1, jlistnum
         dbj = jlist1(dbjj)
         dbnxj = nxdef(dbj)
         do dbii = 1, num
            do dbk = 1, lev
               do dbi = 1, dbnxj + 2
                  dbgv = real(cc(dbi, dbk, dbii, dbjj), kind=8)
                  if (dbgv /= dbgv) then
                     dbna = dbna + 1
                  else
                     dbgca = dbgca + dbgv*dbgv
                  end if
               end do
            end do
         end do
      end do
      if (myrank .eq. 0) print *,'DBGCB cc_after =',dbgca,' nonfinite=',dbna
      ! #endregion
!
      mchk=iand(jtrun,3)

      do j =1, jlistnum
        jj = jlist1(j)
        jtrunj = mtrundef(jj)
        mchk=iand(jtrunj,3)

        do m=1,mchk
          mm= 2*m-1
          mp= mm+1
          mlst=nlist(m)
          do ii=1,num
            do k=1,lev
              twcc_fk(k,1,ii,mlst,j)=cc(mm,k,ii,j)
              twcc_fk(k,2,ii,mlst,j)=cc(mp,k,ii,j)
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
              twcc_fk(k,1,ii,mlst,j)=cc(mm,k,ii,j)
              twcc_fk(k,2,ii,mlst,j)=cc(mp,k,ii,j)
              twcc_fk(k,1,ii,mlst1,j)=cc(mm1,k,ii,j)
              twcc_fk(k,2,ii,mlst1,j)=cc(mp1,k,ii,j)
              twcc_fk(k,1,ii,mlst2,j)=cc(mm2,k,ii,j)
              twcc_fk(k,2,ii,mlst2,j)=cc(mp2,k,ii,j)
              twcc_fk(k,1,ii,mlst3,j)=cc(mm3,k,ii,j)
              twcc_fk(k,2,ii,mlst3,j)=cc(mp3,k,ii,j)
            enddo
          enddo
        enddo

      enddo

      ! #region agent log
      dbgtw = sum(real(twcc_fk, kind=8)**2)
      if (myrank .eq. 0) print *,'DBGTR stwcc_fk=',dbgtw
      ! #endregion
      call mpe_transpose_rs_sp(twcc_fk,wcc_fk,lev*2*num,jtmax,my_max,nsize,col_comm)
      ! #region agent log
      dbgwf = sum(real(wcc_fk, kind=8)**2)
      if (myrank .eq. 0) print *,'DBGTR swcc_fk=',dbgwf
      ! #endregion
!      call mpe_transpose_rs(twcc_fk,wcc_fk,lev*2*num,jtmax,my_max,nsize,col_comm)

      do m=1,mlistnum
         mf=mlist(m)

      fj_wss_sum = 0.0
      fj_wss_dif = 0.0
      fj_polyw_sum = 0.0
      fj_polyw_dif = 0.0
      fj_wccSUM = 0.0
      fj_wccDIF = 0.0

      lchk = iand(jtrun-mf+1, 1)
      lle = jtrun -lchk

      jlistnum_fj = 0
      do j = 1,myhalf
         if(mf .le. mtrundef(j)) then
            jlistnum_fj = jlistnum_fj + 1
            jlist_fj(jlistnum_fj) = j
         endif
      enddo

      do j_fj = 1,jlistnum_fj
         j = jlist_fj(j_fj)
         j1=jlist2(j)
         j2=jlist2(my-j+1)
         do k = 1,lev2*num
            fj_wccSUM(k,j_fj) = (wcc_fk(k,1,1,m,j1)+wcc_fk(k,1,1,m,j2))
            fj_wccDIF(k,j_fj) = (wcc_fk(k,1,1,m,j1)-wcc_fk(k,1,1,m,j2))
         enddo
      enddo

      if(jtrun .gt. mf) then

         do l = mf,lle,2
            l_fj = (l - mf + 2) / 2
            do j_fj = 1,jlistnum_fj
               j = jlist_fj(j_fj)
               fj_polyw_sum(j_fj,l_fj) = poly(l,j,m)*w(j)
               fj_polyw_dif(j_fj,l_fj) = poly(l+1,j,m)*w(j)
            enddo
         enddo

         llistnum_fj = l_fj

         call dgemm('n','n',lev2*num,llistnum_fj,jlistnum_fj,1.0d+0,fj_wccSUM,lev*2*num, &
              fj_polyw_sum,my/2+npad,1.0d+0,fj_wss_sum,lev*2*num)

         do l = mf,lle,2
            l_fj = (l - mf + 2) / 2
            do k = 1,lev2*num
               wss(k,1,1,l,m) = wss(k,1,1,l,m) + fj_wss_sum(k,l_fj)
            enddo
         enddo


         call dgemm('n','n',lev2*num,llistnum_fj,jlistnum_fj,1.0d+0,fj_wccDIF,lev*2*num, &
              fj_polyw_dif,my/2+npad,1.0d+0,fj_wss_dif,lev*2*num)

         do l = mf,lle,2
            l_fj = (l - mf + 2) / 2
            do k = 1,lev2*num
               wss(k,1,1,l+1,m) = wss(k,1,1,l+1,m) + fj_wss_dif(k,l_fj)
            enddo
         enddo

      endif

      if( lchk .eq. 1) then
         l=jtrun
         do j_fj = 1,jlistnum_fj
            j = jlist_fj(j_fj)
            fj_polyw_sum(j_fj,l) = poly(l,j,m)*w(j)
         enddo

         do j_fj = 1,jlistnum_fj
            do k=1,lev2*num
               wss(k,1,1,l  ,m)=wss(k,1,1,l,m)+fj_polyw_sum(j_fj,l)*fj_wccSUM(k,j_fj)
            enddo
         enddo

      endif
!
      enddo

      ! #region agent log
      dbgws = sum(real(wss, kind=8)**2)
      if (myrank .eq. 0) print *,'DBGTR swss=',dbgws
      allocate(dbg_loc(3, jtmax))
      dbg_loc = 0.0d0
      do m = 1, mlistnum
         mf = mlist(m)
         dbg_loc(1, m) = real(mf, kind=8)
         dbg_loc(2, m) = sum(real(wcc_fk(:, :, :, m, :), kind=8)**2)
         dbg_loc(3, m) = sum(real(wss(:, :, :, mf:jtrun, m), kind=8)**2)
      end do
      dbg_gl(1) = sum(dbg_loc(2, :))
      dbg_gl(2) = sum(dbg_loc(3, :))
      dbg_gl(3) = dbgws
      call MPI_COMM_SIZE(MPI_COMM_gfs, dbg_np, dbg_ierr)
      allocate(dbg_all(3, jtmax*dbg_np))
      dbg_all = 0.0d0
      if (dbg_np .eq. 1) then
         ! single rank (the GPU binary's getrdy path): nothing to gather.
         ! 2026-09-15: the 1-rank GPU run died with a glibc top-chunk
         ! assertion right after this block; bypassing MPI here is the
         ! discriminating experiment (standalone gather test was clean).
         dbg_g = dbg_gl
         dbg_all(:, 1:jtmax) = dbg_loc
      else
         call MPI_ALLREDUCE(dbg_gl, dbg_g, 3, MPI_REAL8, MPI_SUM, &
                            MPI_COMM_gfs, dbg_ierr)
         call MPI_GATHER(dbg_loc, 3*jtmax, MPI_REAL8, &
                         dbg_all, 3*jtmax, MPI_REAL8, 0, &
                         MPI_COMM_gfs, dbg_ierr)
      end if
      if (myrank .eq. 0) then
         print *,'DBGTRG call=',dbg_call,' gwcc_fk=',dbg_g(1), &
                 ' gwss_lgem=',dbg_g(2),' gwss_all=',dbg_g(3)
         do dbg_m = 1, jtmax*dbg_np
            if (dbg_all(1, dbg_m) .gt. 0.5d0) &
               print *,'DBGTRM',dbg_call,nint(dbg_all(1, dbg_m)), &
                       dbg_all(2, dbg_m),dbg_all(3, dbg_m)
         end do
      end if
      deallocate(dbg_loc, dbg_all)
      ! #endregion
      return
      end
!      
      subroutine tranrs_dp (jtrun,jtmax,nx,my,my_max,lev,poly,w,cc  &
                        ,wss,num,nsize)
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
!  lev: number of vertical levels to transform
!  poly: legendre polynomials
!  w: gaussian quadrature weights
!  cc_r8: 3-dim input grid pt. field to be transformed
!  num: number of variables grouped together
!
! *** output ***
!
!  wss: spectral coefficient fields
!
!  **********************************
!
      use const, only : RTYPE
      use index
!     use paramt
      use fftcom
      use fj_pad
!
      implicit none

      integer jtrun,jtmax,nx,my,my_max,lev,num,nsize
      integer mlx,myhalf,lev2,nxj,j,jj,mchk,jtrunj,m,mm,mp
      integer mlst,ii,k,mm1,mp1,mlst1,mm2,mp2,mlst2,mm3,mp3,mlst3,mf
      integer lchk,lle,jlistnum_fj,j_fj,j1,j2,l,l_fj,llistnum_fj

      real(kind=RTYPE)    poly(jtrun,my/2,jtmax),w(my)
      real    wss(lev,2,num,jtrun,jtmax)
!
      real    gwk1(nx+2,lev,num,my_max)
!
      real    wcc_fk (lev,2,num,jtmax,my_max*nsize)
      real    twcc_fk(lev,2,num,jtmax*nsize,my_max)
      real(kind=RTYPE)    cc(nx+2,lev,num,my_max)
      real    cc8(nx+2,lev,num,my_max)
!
      real    fj_wss_sum(lev*2*num,jtrun)
      real    fj_wss_dif(lev*2*num,jtrun)
      real    fj_polyw_sum(my/2+npad,jtrun)
      real    fj_polyw_dif(my/2+npad,jtrun)
      real    fj_wccSUM(lev*2*num,my/2)
      real    fj_wccDIF(lev*2*num,my/2)
      integer jlist_fj(my/2)
!
!CWBinit
      twcc_fk=0.
      wss=0.

!CWB2014
      gwk1=0.

!CWB2021 for single precision test

      mlx= (jtrun/2)*((jtrun+1)/2)
      myhalf=my/2
      lev2=lev*2
      cc8=cc
!
!  fft for each guassian latitude of 2-d field
!
      if( length_fft .eq. 0 .and. lreduce.eq.0 )then
!#ifdef SP
!      call rfftmlt_sp(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*num,-1)
!#else
      call rfftmlt(cc8,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*num,-1)
!#endif
      else
!$omp  parallel do default(none)                                &
!$omp  private(jj,j,nxj,gwk1)                                   &
!$omp  shared(jlistnum,jlist1,nxdef,cc8,trigsj,ifaxj,nx,lev,num) &
!$omp  schedule(dynamic)
      do jj=1,jlistnum
        j= jlist1(jj)
        nxj=nxdef(j)
!#ifdef SP
!        call rfftmlt_sp(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j), &
!#else
        call rfftmlt(cc8(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j),  &
!#endif
                     1,nx+2,nxj,lev*num,-1)
      end do
!$omp end parallel do
      end if
!
      mchk=iand(jtrun,3)

      do j =1, jlistnum
        jj = jlist1(j)
        jtrunj = mtrundef(jj)
        mchk=iand(jtrunj,3)

        do m=1,mchk
          mm= 2*m-1
          mp= mm+1
          mlst=nlist(m)
          do ii=1,num
            do k=1,lev
              twcc_fk(k,1,ii,mlst,j)=cc8(mm,k,ii,j)
              twcc_fk(k,2,ii,mlst,j)=cc8(mp,k,ii,j)
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
              twcc_fk(k,1,ii,mlst,j)=cc8(mm,k,ii,j)
              twcc_fk(k,2,ii,mlst,j)=cc8(mp,k,ii,j)
              twcc_fk(k,1,ii,mlst1,j)=cc8(mm1,k,ii,j)
              twcc_fk(k,2,ii,mlst1,j)=cc8(mp1,k,ii,j)
              twcc_fk(k,1,ii,mlst2,j)=cc8(mm2,k,ii,j)
              twcc_fk(k,2,ii,mlst2,j)=cc8(mp2,k,ii,j)
              twcc_fk(k,1,ii,mlst3,j)=cc8(mm3,k,ii,j)
              twcc_fk(k,2,ii,mlst3,j)=cc8(mp3,k,ii,j)
            enddo
          enddo
        enddo

      enddo

!      call mpe_transpose_rs_sp(twcc_fk,wcc_fk,lev*2*num,jtmax,my_max,nsize,col_comm)
      call mpe_transpose_rs(twcc_fk,wcc_fk,lev*2*num,jtmax,my_max,nsize,col_comm)

      do m=1,mlistnum
         mf=mlist(m)

      fj_wss_sum = 0.0
      fj_wss_dif = 0.0
      fj_polyw_sum = 0.0
      fj_polyw_dif = 0.0
      fj_wccSUM = 0.0
      fj_wccDIF = 0.0

      lchk = iand(jtrun-mf+1, 1)
      lle = jtrun -lchk

      jlistnum_fj = 0
      do j = 1,myhalf
         if(mf .le. mtrundef(j)) then
            jlistnum_fj = jlistnum_fj + 1
            jlist_fj(jlistnum_fj) = j
         endif
      enddo

      do j_fj = 1,jlistnum_fj
         j = jlist_fj(j_fj)
         j1=jlist2(j)
         j2=jlist2(my-j+1)
         do k = 1,lev2*num
            fj_wccSUM(k,j_fj) = (wcc_fk(k,1,1,m,j1)+wcc_fk(k,1,1,m,j2))
            fj_wccDIF(k,j_fj) = (wcc_fk(k,1,1,m,j1)-wcc_fk(k,1,1,m,j2))
         enddo
      enddo

      if(jtrun .gt. mf) then

         do l = mf,lle,2
            l_fj = (l - mf + 2) / 2
            do j_fj = 1,jlistnum_fj
               j = jlist_fj(j_fj)
               fj_polyw_sum(j_fj,l_fj) = poly(l,j,m)*w(j)
               fj_polyw_dif(j_fj,l_fj) = poly(l+1,j,m)*w(j)
            enddo
         enddo

         llistnum_fj = l_fj

         call dgemm('n','n',lev2*num,llistnum_fj,jlistnum_fj,1.0d+0,fj_wccSUM,lev*2*num, &
              fj_polyw_sum,my/2+npad,1.0d+0,fj_wss_sum,lev*2*num)

         do l = mf,lle,2
            l_fj = (l - mf + 2) / 2
            do k = 1,lev2*num
               wss(k,1,1,l,m) = wss(k,1,1,l,m) + fj_wss_sum(k,l_fj)
            enddo
         enddo


         call dgemm('n','n',lev2*num,llistnum_fj,jlistnum_fj,1.0d+0,fj_wccDIF,lev*2*num, &
              fj_polyw_dif,my/2+npad,1.0d+0,fj_wss_dif,lev*2*num)

         do l = mf,lle,2
            l_fj = (l - mf + 2) / 2
            do k = 1,lev2*num
               wss(k,1,1,l+1,m) = wss(k,1,1,l+1,m) + fj_wss_dif(k,l_fj)
            enddo
         enddo

      endif

      if( lchk .eq. 1) then
         l=jtrun
         do j_fj = 1,jlistnum_fj
            j = jlist_fj(j_fj)
            fj_polyw_sum(j_fj,l) = poly(l,j,m)*w(j)
         enddo

         do j_fj = 1,jlistnum_fj
            do k=1,lev2*num
               wss(k,1,1,l  ,m)=wss(k,1,1,l,m)+fj_polyw_sum(j_fj,l)*fj_wccSUM(k,j_fj)
            enddo
         enddo

      endif
!
      enddo

      return
      end
