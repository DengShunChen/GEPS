      subroutine tendget (phiten)

!CWB2021 ndsl single precision test

!
!  purpose : compute tendency of vorticity,divergence and geopotential
!            using subroutine of forecast model
!-----------------------------------------------------------------------
!  **** input *****
!  the input are pass bye those '*.h' files (see below)
!
!  **** output *****
!  phiten   : geopotential tendency array
!  other output are pass by 'spec.h'
!-----------------------------------------------------------------------
      use param
      use mpe
      use rank
      use index
      use const
      use spec
      use grid

      implicit none

      integer  m,mf,j,jj,k,n,l,mlst,nxj,nk,kk,kl,i,ii,nl
      real     sqhaf

      real(kind=RTYPE) phiten(levp,2,jtrun,jtmax)
!byl      real tbar(lev),qbar(lev*ncld)
      real(kind=RTYPE) sdpbl(nxp,my_max)

! for Semi-Lagrangian
      real(kind=RTYPE) deldm(nxp,my_max)

!CWB2021 ndsl single precision test
      real(kind=RTYPE) uum_sl(nx,levp,my_max)              &
      ,vvm_sl(nx,levp,my_max)                              &
      ,ttm_sl(nx,levp,my_max)                              &
      ,pten_sl(nx,levp,my_max)                             &
      ,qm_sl(nx,levp*ncld,my_max)                          &
      ,pdot(nxp,lev+1,latpart)                             &
      ,vdmerdr(nxp,lev,my_max),vdzonlr(nxp,lev,my_max)     &
      ,vdmerd(nxp,lev,my_max),vdzonl(nxp,lev,my_max)       &
      ,ddtemp(nxp,lev,my_max),qvadv(nxp,lev*ncld,my_max)   &
      ,diveng(nxp,lev,my_max),pten(nxp,lev,my_max)

!byl      real cc(nx+2,levp,1+ncld,my_max)
      real(kind=RTYPE) cc(nx+2,levp,1,my_max),dummy
!byl,wss(levp,2,1+ncld,jtrun,jtmax)
!
      integer   ierr
      real(kind=RTYPE) dta
      real(kind=RTYPE) ww1(nx,my_max)

!for 2dMPI
      real(kind=RTYPE) temten1(lev,2,jtrun,jtmax)
      real(kind=RTYPE) phiten1(lev,2,jtrun,jtmax)
      ! #region agent log
      real(kind=8) :: dbgtpl, dbgtsp, dbgtar, dbgtt1, dbgtp1
      real(kind=8) :: dbgcc0, dbgcc1, dbgtm0, dbgtm1
      real(kind=8) :: dbgdv, dbgdc0, dbgdc1
      real(kind=8) :: dbgdvu, dbgdvv, dbgdvn
      integer :: dvi, dvk, dvjj, dvnb, dvj, dvnxj
      real(kind=8) :: dbgccj, dbgccn, dbgccv
      integer :: ccnb
      ! #endregion
      logical forward
      forward = .false.
      phiten1=0.

!
!  global mean tempertures (tbar)
!
      sqhaf= sqrt(0.5)
!byl      tbar=0. ; qbar=0.
!     do 160 m = 1, mlistnum
!      mf=mlist(m)
!      if ( mf.eq.1) then
!       do 161 k=1,lev
!         tbar(k)= sqhaf*temnow(k,1,1,m)
! 161   continue
!       do 162 k=1,lev*ncld
!         qbar(k)= sqhaf*qnow(k,1,1,m)
! 162   continue
!      endif
! 160 continue
!

!2dMPI
!byl      do 160 m = 1, mlistnum
!byl       mf=mlist(m)
!byl       if ( mf.eq.1) then
!byl        do 161 k=1,levp
!byl          KK=Llist(k)
!byl          tbar(KK)= sqhaf*temnow(k,1,1,m)
!byl  161   continue
!byl        do 162 n=1,ncld
!byl           nk=(n-1)*levp
!byl           nL=(n-1)*lev 
!byl        do 162 k=1,levp
!byl           KK=nk+k
!byl           KL=nL+Llist(k)
!byl          qbar(KL)= sqhaf*qnow(KK,1,1,m)
!byl  162   continue
!byl       endif
!byl  160 continue

!byl      call mpe_global_sum(tbar,lev,mpe_double)
!byl      call mpe_global_sum(qbar,lev*ncld,mpe_double)
!
      do 170 m = 1, mlistnum
         mf=mlist(m)
      do 170 n = mf, jtrun
      plten(n,m,1) = 0.0
      plten(n,m,2) = 0.0
  170 continue
!
      do 180 m = 1, mlistnum
         mf=mlist(m)
      do 180 n = mf, jtrun
      do 180 k = 1, levp
      divten(k,1,n,m) = 0.0
      divten(k,2,n,m) = 0.0
      vorten(k,1,n,m) = 0.0
      vorten(k,2,n,m) = 0.0
      temten(k,1,n,m) = 0.0
      temten(k,2,n,m) = 0.0
!ncld qten(k,1,n,m) = 0.0
      hldten(k,1,n,m) = 0.0
      hldten(k,2,n,m) = 0.0
      phiten(k,1,n,m) = 0.0
      phiten(k,2,n,m) = 0.0
  180 continue
!
!!      do 183 m = 1, mlistnum
!!         mf=mlist(m)
!!      do 183 n = mf, jtrun
!!      do 183 k = 1, levp*ncld*2
!!      qten(k,1,n,m) = 0.0
!!  183 continue
!
       pdot=0.
       vdmerd=0.
       vdzonl=0.
       qvadv=0.
       ddtemp=0.
       pten=0.
!
       dta=0.5*dt
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do k = 1, lev
          do i = 1, nxj
            up(i,k,jj) = ut(i,k,jj)
            vp(i,k,jj) = vt(i,k,jj)
            ttp(i,k,jj)= tt(i,k,jj)
          enddo
        enddo
        do k = 1, lev*ncld
          do i = 1, nxj
            qm(i,k,jj) = qt(i,k,jj)
          enddo
        enddo
!!        do i = 1, nxj
!!          ptp(i,jj)= pt(i,jj)
!!        enddo
      enddo
!ch>
! transpose partial to full: ut -> ut_sl, vt -> vt_sl, ttm -> ttm_sl, qp -> qm_sl

!#ifdef MULTIPLE
!      call mpe2d_transpose_ndsl_p2f_multi(ut   ,vt   ,ttp   ,dummy,dummy,qp,    &
!                                    ut_sl,vt_sl,ttm_sl,dummy,dummy,qm_sl, &
!                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm,4)
!#else
      call mpe2d_transpose_ndsl_p2f(ut,ut_sl,    &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(vt,vt_sl,    &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(ttp,ttm_sl,  &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(qm,qm_sl,    &
                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!#endif

!!    call mpe2d_unify_nx(pt_sl,pt)
      uum_sl=ut_sl
      vvm_sl=vt_sl
!!    ptp_sl=pt_sl
!ch<

!
!     Semi-Lagrangian
!       Horizontal Advection

          ! #region agent log
          ! ddtemp/vdzonl/vdmerd are already wrong when they reach the
          ! vertical advection (GPU/CPU 0.48 / 0.54 / 1.38, 2026-09-15), and
          ! they are the f2p transposes of ttm_sl/uum_sl/vvm_sl -- so bracket
          ! the horizontal semi-Lagrangian advection that produces those.
          call geps_dbg_ssq_full(ttm_sl,  nx, levp, my_max, 'ttm_sl_pre_hadv')
          call geps_dbg_ssq_full(uum_sl,  nx, levp, my_max, 'uum_sl_pre_hadv')
          call geps_dbg_ssq_full(vvm_sl,  nx, levp, my_max, 'vvm_sl_pre_hadv')
          call geps_dbg_ssq_full(pten_sl, nx, levp, my_max, 'pten_sl_pre_hadv')
          ! #endregion
        call ndslfv_monoadvh2(ttm_sl,qm_sl,pten_sl,uum_sl,vvm_sl,nxdef &
                             ,dta,levp)
          ! #region agent log
          call geps_dbg_ssq_full(ttm_sl,  nx, levp, my_max, 'ttm_sl_post_hadv')
          call geps_dbg_ssq_full(uum_sl,  nx, levp, my_max, 'uum_sl_post_hadv')
          call geps_dbg_ssq_full(vvm_sl,  nx, levp, my_max, 'vvm_sl_post_hadv')
          call geps_dbg_ssq_full(pten_sl, nx, levp, my_max, 'pten_sl_post_hadv')
          ! #endregion

!ch>
! transpose full to partial: ttm_sl -> ddtemp,  pten_sl -> pten, uum_sl -> vdzonl
!                            vvm_sl -> vdmerd,  qm_sl -> qvadv

!#ifdef MULTIPLE
!      call mpe2d_transpose_ndsl_f2p_multi(ttm_sl,pten_sl,uum_sl,vvm_sl,qm_sl, &
!                                    ddtemp   ,pten   ,vdzonl   ,vdmerd   ,qvadv   , &
!                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!#else
      call mpe2d_transpose_ndsl_f2p(ttm_sl,ddtemp, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
          ! #region agent log
          call geps_dbg_ssq_grid(ddtemp, nxp, lev, my_max, 'ddtemp_post_f2p')
          call geps_dbg_ssq_prof(ddtemp, nxp, lev, my_max, 'ddtemp_post_f2p')
          ! #endregion
      call mpe2d_transpose_ndsl_f2p(pten_sl,pten, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(uum_sl,vdzonl, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(vvm_sl,vdmerd, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(qm_sl,qvadv,   &
                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!#endif

!ch<

      do jj =1, jlistnum
      j=jlist1(jj)
!     nxj=nxdef(j)
!
!  p**capa quantities
!
       call prexp_hybrid_cwb (nxjp(j),nxp,lev,ptop,sigma,pt(1,jj),   &
                          pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj))
!
! Calculate Vertical velocity & Stream Functions
!
        call gridnl_hybrid_ndsl (nxjp(j),nxp,lev,ncld                  &
        , cp,radsq,ut(1,1,jj),vt(1,1,jj),rdiv(1,1,jj),tt(1,1,jj)       &
        , qt(1,1,jj),phi(1,1,jj),pt(1,jj),dtpl(1,jj),dlpl(1,jj),sinl(j)&
        , pk(1,1,jj),pk2(1,1,jj),dsigma,sigma,onocos(j),cor(j)         &
        , diveng(1,1,jj),vdmerdr(1,1,jj),vdzonlr(1,1,jj),pten(1,1,jj)  &
        , deldm(1,jj),sdpbl(1,jj),sd(1,1,jj),pdot(1,1,jj),vvel(1,1,jj) &
        , sgeo(1,jj) )

      enddo !jj = 1,jlistnum
!
!CWB2021 for single precision test
          ! #region agent log
          ! Baseline for GPU probes 435/436. The FFT fix stopped tranrs
          ! call 1 being zeroed, revealing its whole chain at 1e15;
          ! diveng is what feeds it. Expected GPU/CPU ratio 32.
          ! NOTE: sum(cc**2) is NaN on CPU (uninitialised padding), so
          ! only diveng is comparable here -- cc is printed to confirm
          ! the NaN, not to compare.
          ! RETRACTED the whole-array version: diveng(nxp,lev,my_max) is
          ! only written for jj<=jlistnum, so sum() reads uninitialised
          ! padding -- its first print was literally NaN, and the 4.876e16
          ! second print (which gave the bogus 'ratio 31.99, diveng is
          ! fine') was very likely padding-dominated too.
          ! join1rs reads r1 over i=1..nxp, k=1..lev, jj=1..jlistnum,
          ! so that is the region that matters.
          dbgdv  = sum(real(diveng, kind=8)**2)   ! kept for comparison
          ! TWO sums over different extents, because the two sides of the
          ! comparison were not using the same one:
          !   dbgdvu : i <= nxp   (what join1rs's read loop covers)
          !   dbgdvn : i <= nxj   (the physically meaningful reduced-grid
          !                        count, and what the cc probe in tranrs uses)
          ! If diveng carries leftover junk in nxj < i <= nxp then dbgdvu is
          ! inflated by it and every ratio computed from it was meaningless.
          dbgdvu = 0.0d0
          dbgdvn = 0.0d0
          dvnb = 0
          do dvjj = 1, jlistnum
             dvj = jlist1(dvjj)
             dvnxj = nxdef(dvj)
             do dvk = 1, lev
                do dvi = 1, nxp
                   dbgdvv = real(diveng(dvi, dvk, dvjj), kind=8)
                   if (dbgdvv /= dbgdvv) then
                      dvnb = dvnb + 1
                   else
                      dbgdvu = dbgdvu + dbgdvv*dbgdvv
                      if (dvi .le. dvnxj) dbgdvn = dbgdvn + dbgdvv*dbgdvv
                   end if
                end do
             end do
          end do
          dbgdc0 = sum(real(cc,     kind=8)**2)
        call joinrs(cc,diveng,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
          ! cc measured with join1rs's OWN write extents:
          !   cc(nx+2, levp, ncld, my_max), written for i<=nx (the
          !   transpose fills only i<=nxjlen and aout=0. zeroes the rest),
          !   k<=levp, n<=ncld(=1), jj<=jlistnum.
          ! The probe inside tranrs summed i<=nxj instead, so if the two
          ! disagree the earlier cc_before number was over the wrong region.
          ! Both extents are accumulated in one pass so no run-to-run
          ! difference can creep in between them.
          dbgdc1 = sum(real(cc,     kind=8)**2)   ! whole array (expect NaN)
          dbgccn = 0.0d0
          dbgccj = 0.0d0
          ccnb = 0
          do dvjj = 1, jlistnum
             dvj = jlist1(dvjj)
             dvnxj = nxdef(dvj)
             do dvk = 1, levp
                do dvi = 1, nx
                   dbgccv = real(cc(dvi, dvk, 1, dvjj), kind=8)
                   if (dbgccv /= dbgccv) then
                      ccnb = ccnb + 1
                   else
                      dbgccn = dbgccn + dbgccv*dbgccv
                      if (dvi .le. dvnxj) dbgccj = dbgccj + dbgccv*dbgccv
                   end if
                end do
             end do
          end do
          if (myrank .eq. 0) print *,'DBGCJ cc_i_le_nx=',dbgccn,' nonfinite=',ccnb
          if (myrank .eq. 0) print *,'DBGCJ cc_i_le_nxj=',dbgccj
          if (myrank .eq. 0) print *,'DBGDV sdiveng_all=',dbgdv
          if (myrank .eq. 0) print *,'DBGDV sdiveng_nxp=',dbgdvu,' nonfinite=',dvnb
          if (myrank .eq. 0) print *,'DBGDV sdiveng_nxj=',dbgdvn
          if (myrank .eq. 0) print *,'DBGDV scc_before=',dbgdc0,' scc_after=',dbgdc1
          ! #endregion
          ! #region agent log
          ! tranrs is a SHARED routine: getrdy.f90:779, intgrt.f90 (seven
          ! sites), incrini.f90:186 and tendget all call it, and getrdy runs
          ! first. 39 DBGCB records appear in one run, so "the first match"
          ! picked getrdy's call, not this one -- that invalidated several
          ! rounds of cc_before comparisons. This marker lets the analyser
          ! take the record that belongs to THIS call site.
          if (myrank .eq. 0) print *,'DBGMARK tendget_tranrs_1'
          ! #endregion
        call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc       &
                   ,hldten,1,nsizey)

          ! #region agent log
          ! hldten = the FIRST tranrs pair's output, same tranrs code as
          ! the one producing temten. Measured before trngra3 consumes it.
          dbgcc0 = sum(real(hldten, kind=8)**2)
          if (myrank .eq. 0) print *,'DBGHL shldten=',dbgcc0
          ! #endregion
        call trngra3(jtrun,jtmax,nx,levp,my,my_max,cim,poly,dpoly      &
                   ,hldten,dlphi,dtphi,nsizey)
!
!       update all horizontal informations
!
          ! #region agent log
          ! With the layer-10 memset race fixed, ddtemp is exact but
          ! vdzonl/vdmerd are still 1.9x (2026-09-15) although uum_sl/vvm_sl
          ! (their f2p sources) are exact: bracket the tendency update.
          call geps_dbg_ssq_grid(vdzonl,  nxp, lev, my_max, 'vdzonl_post_f2p')
          call geps_dbg_ssq_grid(vdmerd,  nxp, lev, my_max, 'vdmerd_post_f2p')
          call geps_dbg_ssq_grid(vdzonlr, nxp, lev, my_max, 'vdzonlr_pre_upd')
          call geps_dbg_ssq_grid(vdmerdr, nxp, lev, my_max, 'vdmerdr_pre_upd')
          call geps_dbg_ssq_grid(dlphi,   nxp, lev, my_max, 'dlphi_pre_upd')
          call geps_dbg_ssq_grid(dtphi,   nxp, lev, my_max, 'dtphi_pre_upd')
          ! #endregion
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            do k=1,lev
              vdmerdr(i,k,jj)=vdmerdr(i,k,jj)-dtphi(i,k,jj)/radsq/onocos(j)! &
!                             -ut(i,k,jj)*cor(j)-(ut(i,k,jj)*ut(i,k,jj)      &
!                             +vt(i,k,jj)*vt(i,k,jj))*onocos(j)*sinl(j) 

              vdzonlr(i,k,jj)=vdzonlr(i,k,jj)-dlphi(i,k,jj)/radsq!           &
!                             +vt(i,k,jj)*cor(j)
            enddo
          enddo
        enddo !jj = 1,jlistnum
          ! #region agent log
          call geps_dbg_ssq_grid(vdzonlr, nxp, lev, my_max, 'vdzonlr_post_loop')
          call geps_dbg_ssq_grid(vdmerdr, nxp, lev, my_max, 'vdmerdr_post_loop')
          ! #endregion
      call ndslfv_update(nxjp,vdzonl,vdmerd,vdzonlr,vdmerdr,dta,forward)
          ! #region agent log
          call geps_dbg_ssq_grid(vdzonl,  nxp, lev, my_max, 'vdzonl_post_upd')
          call geps_dbg_ssq_grid(vdmerd,  nxp, lev, my_max, 'vdmerd_post_upd')
          ! #endregion

!CWB2021 ndsl single precision test
!
!
!       Vertical Advection
!
          ! #region agent log
          ! tendget's 2nd tranrs input is 1.3e5x too large on the GPU while
          ! the 1st is exact (2026-09-15 per-mf comparison), so the stage
          ! that makes ddtemp is bracketed here: before / after the vertical
          ! advection, and after the (ddtemp-ttp)/dt step, plus its inputs.
          call geps_dbg_ssq_grid(ddtemp, nxp, lev, my_max, 'ddtemp_pre_vadv')
          call geps_dbg_ssq_grid(vdzonl, nxp, lev, my_max, 'vdzonl_pre_vadv')
          call geps_dbg_ssq_grid(vdmerd, nxp, lev, my_max, 'vdmerd_pre_vadv')
          call geps_dbg_ssq_grid(pdot,   nxp, lev+1, latpart, 'pdot_pre_vadv')
          call geps_dbg_ssq_grid(pt,     nxp, 1, my_max, 'pt_pre_vadv')
          call geps_dbg_ssq_grid(ttp,    nxp, lev, my_max, 'ttp_pre_vadv')
          ! #endregion
      call ndslfv_monoadvv(ddtemp,qvadv,vdzonl,vdmerd,pdot,pt,nxjp,dta,forward)
          ! #region agent log
          call geps_dbg_ssq_grid(ddtemp, nxp, lev, my_max, 'ddtemp_post_vadv')
          call geps_dbg_ssq_grid(vdzonl, nxp, lev, my_max, 'vdzonl_post_vadv')
          call geps_dbg_ssq_grid(vdmerd, nxp, lev, my_max, 'vdmerd_post_vadv')
          ! #endregion

!CWB2021 ndsl single precision test

!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do k = 1, lev
          do i = 1, nxj
            vdzonl(i,k,jj) = (vdzonl(i,k,jj)- up(i,k,jj))/dt
            vdmerd(i,k,jj) = (vp(i,k,jj)- vdmerd(i,k,jj))/dt
            ddtemp(i,k,jj)=  (ddtemp(i,k,jj)-ttp(i,k,jj))/dt
          enddo
        enddo
      enddo
          ! #region agent log
          call geps_dbg_ssq_grid(ddtemp, nxp, lev, my_max, 'ddtemp_post_dt')
          ! #endregion
!
!  combine non-linear grid point terms via gaussian quadrature
!
          ! #region agent log
          ! Baseline for the GPU probes 425/426/427. GPU measured
          !   cc_cg 0 -> 1.259e7 (joinrs), temten 0 -> 3.333e4 (tranrs)
          ! Expected GPU/CPU ratio is 32 for both (spectral, rank-split).
          ! Whichever side is already off names the broken stage.
          dbgcc0 = sum(real(cc,     kind=8)**2)
          dbgtm0 = sum(real(temten, kind=8)**2)
      call joinrs(cc,ddtemp,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
          dbgcc1 = sum(real(cc,     kind=8)**2)
          if (myrank .eq. 0) print *,'DBGCC joinrs scc_before=',dbgcc0,' scc_after=',dbgcc1
          ! #endregion
          ! #region agent log
          if (myrank .eq. 0) print *,'DBGMARK tendget_tranrs_2'
          ! #endregion
      call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc  &
                   ,temten,1,nsizey)
          ! #region agent log
          dbgtm1 = sum(real(temten, kind=8)**2)
          if (myrank .eq. 0) print *,'DBGCC tranrs stemten_before=',dbgtm0,' stemten_after=',dbgtm1
          ! #endregion
!!      call joinrs(cc,ddtemp,qvadv,dummy,dummy,nx,my_max,lev           &
!!                 ,jlistnum,2,ncld)
!!      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc         &
!!                 ,wss,1+ncld,nsize)
!!      call ujoinrs(wss,temten,qten,dummy,dummy,jtrun,jtmax,lev        &
!!                 ,mlistnum,2,ncld)
      call mpe2d_unify_nx(ww1,deldm)
      call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,ww1           &
                 ,plten,nsizey)
!
!CWB2021 for single precision test
      call rstrandz (jtrun,jtmax,nx,my,my_max,levp,vdmerd,vdzonl      &
                    ,weight,cim,onocos,poly,dpoly,divten,vorten,nsizey)
!
!  no tendency of zero mode
!
      do 96 m =1,mlistnum
        mf=mlist(m)
        if ( mf.eq.1) then
        do k = 1, levp
         vorten(k,1,1,m)= 0.
         vorten(k,2,1,m)= 0.
         divten(k,1,1,m)= 0.
         divten(k,2,1,m)= 0.
         temten(k,1,1,m)= 0.
         temten(k,2,1,m)= 0.
        enddo
        endif
 96   continue
!
!  compute phiten from temden and plten
!
         mlst=ilist(1)
         if(mlst .ne. 0) then
           plten(1,mlst,1)= 0.
           plten(1,mlst,2)= 0.
         endif
!
!       do 100 m =1,mlistnum
!        mf=mlist(m)
!       do 100 n  = mf,jtrun
!       do 99 k = 1, lev
!         phiten(k,1,n,m) = spalm(k)*plten(n,m,1)
!         phiten(k,2,n,m) = spalm(k)*plten(n,m,2)
!       do 98 l = 1, lev
!         phiten(k,1,n,m) = phiten(k,1,n,m)+arrhyd(k,l)*temten(l,1,n,m)
!         phiten(k,2,n,m) = phiten(k,2,n,m)+arrhyd(k,l)*temten(l,2,n,m)
!98     continue
!99     continue
!100    continue
!
!ch>
        call mpe2d_unify_lev(temten,temten1,lev,levp,jtrun,jtmax,mlistnum,nsizex,row_comm)
          ! #region agent log
          ! Baseline for the GPU probes 420/421/423 inside tendget_gpu.
          ! Expected ratios GPU/CPU:
          !   arrhyd, spalm  -> 1   (per-level constants, NOT split by rank)
          !   plten, temten1 -> 32  (spectral, split across the 32 ranks)
          ! Whichever ratio is off names the broken quantity.
          dbgtpl = sum(real(plten,   kind=8)**2)
          dbgtsp = sum(real(spalm,   kind=8)**2)
          dbgtar = sum(real(arrhyd,  kind=8)**2)
          dbgtt1 = sum(real(temten1, kind=8)**2)
          if (myrank .eq. 0) print *,'DBGTG in  splten=',dbgtpl,' sspalm=',dbgtsp
          if (myrank .eq. 0) print *,'DBGTG in  sarrhyd=',dbgtar,' stemten1=',dbgtt1
          ! #endregion

        do 100 m =1,mlistnum
         mf=mlist(m)
        do 100 n  = mf,jtrun
        do 99 k = 1, lev
          phiten1(k,1,n,m) = spalm(k)*plten(n,m,1)
          phiten1(k,2,n,m) = spalm(k)*plten(n,m,2)
        do 98 L = 1, lev
          phiten1(k,1,n,m) = phiten1(k,1,n,m)+arrhyd(k,L)*temten1(L,1,n,m)
          phiten1(k,2,n,m) = phiten1(k,2,n,m)+arrhyd(k,L)*temten1(L,2,n,m)
 98     continue
 99     continue
 100    continue
          ! #region agent log
          ! CPU fuses init and accumulate in one nest, so there is no
          ! separate 'after init' point like the GPU's probe 422.
          dbgtp1 = sum(real(phiten1, kind=8)**2)
          if (myrank .eq. 0) print *,'DBGTG out sphiten1=',dbgtp1
          ! #endregion

        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun
            do k = 1, levp
               kk=Llist(k)
               phiten(k,1,n,m) = phiten1(kk,1,n,m)
               phiten(k,2,n,m) = phiten1(kk,2,n,m)
            enddo
          enddo
        enddo
!ch<

      return
      end

! #region agent log
! Used-region, NaN-safe sum of squares of an (nxp, lev, my_max) grid field,
! reduced to a GLOBAL value over MPI_COMM_gfs so the GPU's single rank can be
! compared 1:1 (grid fields are latitude-split, so rank 0 x 32 would also be
! valid, but the global sum needs no such assumption). Only i <= nxdef_2d(j),
! jj <= jlistnum is summed: the rest is uninitialised padding (7.1 rule 2).
! MPI_REAL8 on purpose: MPI_DOUBLE_PRECISION is 16 bytes in this OpenMPI
! build (configured with -fdefault-real-8 only), see ROCM_PORT_HANDOFF.md.
! A 1-rank caller skips MPI entirely: the 1-rank GPU binary corrupted its
! heap on MPI_GATHER/ALLREDUCE on 2026-09-15 (cause not yet identified).
      subroutine geps_dbg_ssq_grid(a, n1, n2, n3, tag)
      use rank, only : myrank, MPI_COMM_gfs
      use index, only : jlistnum, jlist1, nxdef_2d
      use mpi
      implicit none
      integer, intent(in) :: n1, n2, n3
      real(kind=8), intent(in) :: a(n1, n2, n3)
      character(len=*), intent(in) :: tag
      real(kind=8) :: sl(2), sg(2), v
      integer :: i, k, jj, j, nxj, np, ierr, nb
      sl = 0.0d0
      nb = 0
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do k = 1, n2
            do i = 1, min(nxj, n1)
               v = a(i, k, jj)
               if (v /= v) then
                  nb = nb + 1
               else
                  sl(1) = sl(1) + v*v
               end if
            end do
         end do
      end do
      sl(2) = real(nb, kind=8)
      call MPI_COMM_SIZE(MPI_COMM_gfs, np, ierr)
      if (np .eq. 1) then
         sg = sl
      else
         call MPI_ALLREDUCE(sl, sg, 2, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
      end if
      if (myrank .eq. 0) print *,'DBGGRID ',tag,' ssq=',sg(1),' nonfinite=',nint(sg(2))
      end subroutine geps_dbg_ssq_grid
! #endregion

! #region agent log
! Same as geps_dbg_ssq_grid but for the full-longitude (_sl) layout, whose
! used region is i <= nxdef(j) (reduced grid) rather than the partial-lon
! nxdef_2d(j). Identical at NPEX=1, kept separate so the intent is explicit.
      subroutine geps_dbg_ssq_full(a, n1, n2, n3, tag)
      use rank, only : myrank, MPI_COMM_gfs
      use index, only : jlistnum, jlist1, nxdef
      use mpi
      implicit none
      integer, intent(in) :: n1, n2, n3
      real(kind=8), intent(in) :: a(n1, n2, n3)
      character(len=*), intent(in) :: tag
      real(kind=8) :: sl(2), sg(2), v
      integer :: i, k, jj, j, nxj, np, ierr, nb
      sl = 0.0d0
      nb = 0
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef(j)
         do k = 1, n2
            do i = 1, min(nxj, n1)
               v = a(i, k, jj)
               if (v /= v) then
                  nb = nb + 1
               else
                  sl(1) = sl(1) + v*v
               end if
            end do
         end do
      end do
      sl(2) = real(nb, kind=8)
      call MPI_COMM_SIZE(MPI_COMM_gfs, np, ierr)
      if (np .eq. 1) then
         sg = sl
      else
         call MPI_ALLREDUCE(sl, sg, 2, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
      end if
      if (myrank .eq. 0) print *,'DBGGRID ',tag,' ssq=',sg(1),' nonfinite=',nint(sg(2))
      end subroutine geps_dbg_ssq_full
! #endregion

! #region agent log
! Profile of an (n1, n2, n3) grid field: used-region sum of squares broken
! down by level k, by latitude band (16 bands over jlist order) and by
! relative longitude band (16 bands of i/nxj). Global over MPI_COMM_gfs.
! Purpose (2026-09-15): ddtemp after the GPU f2p transpose holds ~0.48x the
! CPU's energy; this says WHICH index range is missing.
      subroutine geps_dbg_ssq_prof(a, n1, n2, n3, tag)
      use rank, only : myrank, MPI_COMM_gfs
      use index, only : jlistnum, jlist1, nxdef_2d
      use mpi
      implicit none
      integer, intent(in) :: n1, n2, n3
      real(kind=8), intent(in) :: a(n1, n2, n3)
      character(len=*), intent(in) :: tag
      real(kind=8), allocatable :: pk(:), gk(:)
      real(kind=8) :: pj(16), gj(16), pi(16), gi(16), v
      integer :: i, k, jj, j, nxj, np, ierr, bj, bi
      allocate(pk(n2), gk(n2))
      pk = 0.0d0; pj = 0.0d0; pi = 0.0d0
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         bj = min(16, 1 + (16*(j - 1))/max(1, n3))
         do k = 1, n2
            do i = 1, min(nxj, n1)
               v = a(i, k, jj)
               if (v /= v) cycle
               v = v*v
               pk(k) = pk(k) + v
               pj(bj) = pj(bj) + v
               bi = min(16, 1 + (16*(i - 1))/max(1, nxj))
               pi(bi) = pi(bi) + v
            end do
         end do
      end do
      call MPI_COMM_SIZE(MPI_COMM_gfs, np, ierr)
      if (np .eq. 1) then
         gk = pk; gj = pj; gi = pi
      else
         call MPI_ALLREDUCE(pk, gk, n2, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
         call MPI_ALLREDUCE(pj, gj, 16, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
         call MPI_ALLREDUCE(pi, gi, 16, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
      end if
      if (myrank .eq. 0) then
         do k = 1, n2
            print *,'DBGPROF ',tag,' k=',k,' ssq=',gk(k)
         end do
         do bj = 1, 16
            print *,'DBGPROF ',tag,' jband=',bj,' ssq=',gj(bj)
         end do
         do bi = 1, 16
            print *,'DBGPROF ',tag,' iband=',bi,' ssq=',gi(bi)
         end do
      end if
      deallocate(pk, gk)
      end subroutine geps_dbg_ssq_prof
! #endregion

! #region agent log
! Spectral sum of squares over the triangular region n >= mlist(m) only (the
! part every routine writes; the rest is never initialised), reduced to a
! global value. `a` is viewed as (n1, jtrun, jtmax): pass n1 = lev*2 for the
! (lev, 2, jtrun, jtmax) arrays. geps_dbg_ssq_spec2 is the (jtrun, jtmax, n3)
! layout used by plten/plmid/plnow.
      subroutine geps_dbg_ssq_spec(a, n1, jtrun, jtmax, tag)
      use rank, only : myrank, MPI_COMM_gfs
      use index, only : mlistnum, mlist
      use mpi
      implicit none
      integer, intent(in) :: n1, jtrun, jtmax
      real(kind=8), intent(in) :: a(n1, jtrun, jtmax)
      character(len=*), intent(in) :: tag
      real(kind=8) :: sl(2), sg(2), v
      integer :: i, n, m, mf, np, ierr, nb
      sl = 0.0d0
      nb = 0
      do m = 1, mlistnum
         mf = mlist(m)
         do n = mf, jtrun
            do i = 1, n1
               v = a(i, n, m)
               if (v /= v) then
                  nb = nb + 1
               else
                  sl(1) = sl(1) + v*v
               end if
            end do
         end do
      end do
      sl(2) = real(nb, kind=8)
      call MPI_COMM_SIZE(MPI_COMM_gfs, np, ierr)
      if (np .eq. 1) then
         sg = sl
      else
         call MPI_ALLREDUCE(sl, sg, 2, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
      end if
      if (myrank .eq. 0) print *,'DBGSPEC ',tag,' ssq=',sg(1),' nonfinite=',nint(sg(2))
      end subroutine geps_dbg_ssq_spec

      subroutine geps_dbg_ssq_spec2(a, jtrun, jtmax, n3, tag)
      use rank, only : myrank, MPI_COMM_gfs
      use index, only : mlistnum, mlist
      use mpi
      implicit none
      integer, intent(in) :: jtrun, jtmax, n3
      real(kind=8), intent(in) :: a(jtrun, jtmax, n3)
      character(len=*), intent(in) :: tag
      real(kind=8) :: sl(2), sg(2), v
      integer :: k, n, m, mf, np, ierr, nb
      sl = 0.0d0
      nb = 0
      do k = 1, n3
         do m = 1, mlistnum
            mf = mlist(m)
            do n = mf, jtrun
               v = a(n, m, k)
               if (v /= v) then
                  nb = nb + 1
               else
                  sl(1) = sl(1) + v*v
               end if
            end do
         end do
      end do
      sl(2) = real(nb, kind=8)
      call MPI_COMM_SIZE(MPI_COMM_gfs, np, ierr)
      if (np .eq. 1) then
         sg = sl
      else
         call MPI_ALLREDUCE(sl, sg, 2, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
      end if
      if (myrank .eq. 0) print *,'DBGSPEC ',tag,' ssq=',sg(1),' nonfinite=',nint(sg(2))
      end subroutine geps_dbg_ssq_spec2
! #endregion

! #region agent log
! Column-wise accumulator for the CPU physics driver (diabat.f90's big
! `do 290 jj` loop): geps_dbg_acc adds the used-region sum of squares of one
! column (a(1:nxj, 1:n2)) under `tag`; geps_dbg_acc_flush (after the loop)
! reduces over ranks and prints in the DBGGRID format so cmp_ckpt.py can
! read it next to the GPU's whole-domain g_* checkpoints.
      module geps_dbg_acc_mod
      implicit none
      integer, parameter :: mxtag = 16
      character(len=32) :: tags(mxtag) = ''
      real(kind=8) :: acc(2, mxtag) = 0.0d0
      integer :: ntag = 0
      integer :: kinds(mxtag) = 0     ! 0 = sum of squares, 1 = count(>0) / sum
      end module geps_dbg_acc_mod

! integer companion: acc(1) counts entries > 0, acc(2) sums them
      subroutine geps_dbg_acc_int(ia, n1, nxj, tag)
      use geps_dbg_acc_mod
      implicit none
      integer, intent(in) :: n1, nxj
      integer, intent(in) :: ia(n1)
      character(len=*), intent(in) :: tag
      integer :: i, t
      do t = 1, ntag
         if (tags(t) == tag) exit
      end do
      if (t > ntag) then
         if (ntag >= mxtag) return
         ntag = ntag + 1; t = ntag; tags(t) = tag; kinds(t) = 1
      end if
      do i = 1, min(nxj, n1)
         if (ia(i) > 0) acc(1, t) = acc(1, t) + 1.0d0
         acc(2, t) = acc(2, t) + real(ia(i), kind=8)
      end do
      end subroutine geps_dbg_acc_int

! GPU-side whole-domain companion of geps_dbg_acc_int (same output line)
      subroutine geps_dbg_cnt_grid(ia, n1, n3, tag)
      use rank, only : myrank, MPI_COMM_gfs
      use index, only : jlistnum, jlist1, nxdef_2d
      use mpi
      implicit none
      integer, intent(in) :: n1, n3
      integer, intent(in) :: ia(n1, n3)
      character(len=*), intent(in) :: tag
      real(kind=8) :: sl(2), sg(2)
      integer :: i, jj, j, nxj, np, ierr
      sl = 0.0d0
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do i = 1, min(nxj, n1)
            if (ia(i, jj) > 0) sl(1) = sl(1) + 1.0d0
            sl(2) = sl(2) + real(ia(i, jj), kind=8)
         end do
      end do
      call MPI_COMM_SIZE(MPI_COMM_gfs, np, ierr)
      if (np .eq. 1) then
         sg = sl
      else
         call MPI_ALLREDUCE(sl, sg, 2, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
      end if
      if (myrank .eq. 0) print *,'DBGCNT ',tag,' count=',nint(sg(1)),' sum=',nint(sg(2))
      end subroutine geps_dbg_cnt_grid

      subroutine geps_dbg_acc(a, n1, nxj, n2, tag)
      use geps_dbg_acc_mod
      implicit none
      integer, intent(in) :: n1, nxj, n2
      real(kind=8), intent(in) :: a(n1, n2)
      character(len=*), intent(in) :: tag
      real(kind=8) :: v
      integer :: i, k, t
      do t = 1, ntag
         if (tags(t) == tag) exit
      end do
      if (t > ntag) then
         if (ntag >= mxtag) return
         ntag = ntag + 1; t = ntag; tags(t) = tag
      end if
      do k = 1, n2
         do i = 1, min(nxj, n1)
            v = a(i, k)
            if (v /= v) then
               acc(2, t) = acc(2, t) + 1.0d0
            else
               acc(1, t) = acc(1, t) + v*v
            end if
         end do
      end do
      end subroutine geps_dbg_acc

      subroutine geps_dbg_acc_flush()
      use geps_dbg_acc_mod
      use rank, only : myrank, MPI_COMM_gfs
      use mpi
      implicit none
      real(kind=8) :: sg(2)
      integer :: t, np, ierr
      call MPI_COMM_SIZE(MPI_COMM_gfs, np, ierr)
      do t = 1, ntag
         if (np .eq. 1) then
            sg = acc(:, t)
         else
            call MPI_ALLREDUCE(acc(:, t), sg, 2, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
         end if
         if (myrank .eq. 0) then
            if (kinds(t) == 1) then
               print *,'DBGCNT ',trim(tags(t)),' count=',nint(sg(1)),' sum=',nint(sg(2))
            else
               print *,'DBGGRID ',trim(tags(t)),' ssq=',sg(1),' nonfinite=',nint(sg(2))
            end if
         end if
      end do
      acc = 0.0d0
      end subroutine geps_dbg_acc_flush
! #endregion

! #region agent log
! count of .true. in the used region of a logical (n1, n3) flag array
      subroutine geps_dbg_cntl_grid(l, n1, n3, tag)
      use rank, only : myrank, MPI_COMM_gfs
      use index, only : jlistnum, jlist1, nxdef_2d
      use mpi
      implicit none
      integer, intent(in) :: n1, n3
      logical, intent(in) :: l(n1, n3)
      character(len=*), intent(in) :: tag
      real(kind=8) :: sl, sg
      integer :: i, jj, j, nxj, np, ierr
      sl = 0.0d0
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do i = 1, min(nxj, n1)
            if (l(i, jj)) sl = sl + 1.0d0
         end do
      end do
      call MPI_COMM_SIZE(MPI_COMM_gfs, np, ierr)
      if (np .eq. 1) then
         sg = sl
      else
         call MPI_ALLREDUCE(sl, sg, 1, MPI_REAL8, MPI_SUM, MPI_COMM_gfs, ierr)
      end if
      if (myrank .eq. 0) print *,'DBGCNT ',tag,' count=',nint(sg)
      end subroutine geps_dbg_cntl_grid
! #endregion
