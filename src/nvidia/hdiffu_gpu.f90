#define NCCLCHECK(ierr) call nccl_check_helper(ierr, __FILE__, __LINE__)

subroutine hdiffu_gpu(dta, my, my_max, nx, jtrun, jtmax, lev, ncld, amp &
                      , rad, cosl, ut, vt, vornow, divnow, temnow, eps4 &
                      , trefs)
   ! Present on device: cosl, ut, vt, vornow, divnow, temnow, eps4, trefs
   ! Present on device: jlist1, nxdef_2d, Llist, hdk2
   use index
   use mpe
   use rank
   use const, only: hdk1, hdk2, radsq, doskeb, onocos, wcfac, wdfac &
                    , poly, dpoly, hord, vd, factop, RTYPE, coslr, polyf, dpolyf
   use param, only: octahedral
   use openacc
   use cudafor
   use mpi
   use spec_cuda_graph, only: cc_cg, gwk1_cg, wcc_fk_cg, &
                              wc_cg, ws_cg, fj_weight_cg

   implicit none

   integer my, my_max, nx, jtrun, jtmax, lev, ncld
   real rad

   real(kind=RTYPE) dta
   real(kind=RTYPE) vordiss(levp, 2, jtrun, jtmax), divdiss(levp, 2, jtrun, jtmax)

   real(kind=RTYPE) vornow(levp, 2, jtrun, jtmax), divnow(levp, 2, jtrun, jtmax), &
      temnow(levp, 2, jtrun, jtmax), trefs(levp, 2, jtrun, jtmax), &
      ut(nxp, lev, my_max), vt(nxp, lev, my_max), eps4(jtrun, jtmax), &
      cosl(my)
!
!     parameter ( ktop=4, ktop2=ktop/2 ) ! top "ktop" levels are inhenced
!
   real wmax(lev), wmax_buf(lev)
   real windmax1, windmax2, windmax3

   integer jj, j, nxj, k, i, m, n, mf, nc, kk, KL
   real xx, facd, facv, fact, amp, ddiffu, vdiffu, tdiffu
   real hfilt, nf, dec, coefu, powd, kfac, dect, hfilt2, hfiltd, powdd
   real c1, c2, c3
   logical windchk
   real wt

   data windmax1/80./, windmax2/100./, windmax3/130./
!!      data      windmax1/70./, windmax2/100./, windmax3/130./
!
   integer async_id, istat, ierr
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   nf = jtrun - 1

   !$acc enter data create(wmax_buf, vordiss, divdiss, wmax) async(async_id)
   !$acc kernels async(async_id)
   wmax = 0.0
   !$acc end kernels
   !$acc parallel loop gang async(async_id) private(wt)
   do k = 1, lev
      wt = 0.0
      !$acc loop vector collapse(2) reduction(max:wt) private(j,nxj,xx)
      do jj = 1, jlistnum
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            if (i .le. nxj) then
               xx = rad/cosl(j)
               wt = max(wt, xx*sqrt(ut(i, k, jj)**2 + vt(i, k, jj)**2))
            end if
         end do
      end do
      wmax(k) = wt
   end do

   !$acc host_data use_device(wmax, wmax_buf)
   NCCLCHECK(ncclAllReduce(wmax, wmax_buf, lev, ncclFloat64, ncclMax, nccl_comm_gfs, stream))
   !$acc end host_data

   !$acc parallel loop async(async_id)
   do i = 1, lev
      wmax(i) = wmax_buf(i)
   end do

   powd = float(hord)/2.
   hfilt = (radsq/(nf*(nf + 1)))**powd
   hfilt2= radsq/(nf*(nf+1))
   coefu = factop/float(hdk2(1) - hdk1)
   if (octahedral) then
      hfilt = hfilt/(6.*dta)
      hfilt2= hfilt2/(6.*dta)
   else
      hfilt = hfilt/dta
      hfilt2= hfilt2/dta
   end if

   !$acc parallel loop gang private(KL, kfac, dect, dec, facd, facv, fact) async(async_id)
   do k = 1, levp
      KL = Llist(k)
      kfac = min(coefu*max(float(hdk2(1) - KL), 0.), factop)
      dect = float(min(max(hdk1 - KL, 1 - hdk1), hdk1 - 1))/float((hdk1 - 1))
      if (dect .ge. 0.) then
         dec = 0.5*(1.+dect**(1./3.))
      else
         dect = -1.*dect
         dec = 0.5*(1.-dect**(1./3.))
      end if
      kfac = kfac*(1.+vd*dec)
      facd = max(1., kfac)*amp
      facv = max(1., kfac)*amp
!      fact = max(1., kfac)*amp

      if ( KL .lt. hdk1 ) then
        hfiltd = hfilt2
        powdd  = 1.
      else
        hfiltd = hfilt
        powdd  = powd
      endif

      !$acc loop worker
      do m = 1, mlistnum
         mf = mlist(m)
         !$acc loop vector
         do n = mf, jtrun
            c1 = 1.+dta*facv*hfilt*eps4(n, m)**powd
!            c2 = 1.+dta*facd*hfilt*eps4(n, m)**powd
            c2 = 1.+dta*facd*hfiltd*eps4(n, m)**powdd
!            c3 = 1.+dta*fact*hfilt*eps4(n, m)**powd
            vordiss(k, 1, n, m) = (1.-1./c1)*vornow(k, 1, n, m)
            vordiss(k, 2, n, m) = (1.-1./c1)*vornow(k, 2, n, m)
            divdiss(k, 1, n, m) = (1.-1./c2)*divnow(k, 1, n, m)
            divdiss(k, 2, n, m) = (1.-1./c2)*divnow(k, 2, n, m)
            vornow(k, 1, n, m) = vornow(k, 1, n, m)/c1
            vornow(k, 2, n, m) = vornow(k, 2, n, m)/c1
            divnow(k, 1, n, m) = divnow(k, 1, n, m)/c2
            divnow(k, 2, n, m) = divnow(k, 2, n, m)/c2
!            temnow(k, 1, n, m) = (temnow(k, 1, n, m) + (c3 - 1.)*trefs(k, 1, n, m))/c3
!            temnow(k, 2, n, m) = (temnow(k, 2, n, m) + (c3 - 1.)*trefs(k, 2, n, m))/c3
         end do
      end do
   end do
!
!       estimate the dissipation of kinetic energy
   if (doskeb) then
      !$acc wait(async_id)
      ! transfer the dissipation to grid point from spectrum
      call tranuv_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, &
                                 coslr, wcfac, wdfac, polyf, dpolyf, &
                                 vordiss, divdiss, ut, vt, nsizey, &
                                 cc_cg, gwk1_cg, ws_cg(1, 1), ws_cg(1, 2), wc_cg, &
                                 wcc_fk_cg, fj_weight_cg(1, 1), fj_weight_cg(1, 2))
      !$acc wait(async_id)
   end if
   !$acc exit data copyout(wmax) delete(wmax_buf, vordiss, divdiss) async(async_id)
   !$acc wait(async_id)
!
!-------------------------------------------------------------------
!
!  filter layers near the upper bound if too strong wind speed
!  happens at the top layer
!
!  2003/10/7 :
!  sometimes wind speed greater than 100m/s happens at k=2
!
   windchk = .false.
   do k = 1, hdk1
      if (wmax(k) .gt. windmax3) windchk = .true.
   end do
   if (windchk) then
      ! Present on device: temnow, vornow, divnow, Llist, mlist
      call filter_top_gpu(jtrun, jtmax, levp, hdk1, ncld, temnow, vornow, divnow)
   end if
!--------------------------------------------------------------------
   return
end

subroutine filter_top_gpu(jtrun, jtmax, lev, ktop, ncld, temnow &
                          , vornow, divnow)
   ! Present on device: temnow, vornow, divnow, Llist, mlist
!
!  apply Lanczos filter to top "ktop" layers
!
   use index
   use mpe
   use const, only: RTYPE
   use openacc
   use cudafor
!
   implicit none

   integer ktop, ktopm1
   integer jtrun, jtmax, lev, ncld

   real(kind=RTYPE) temnow(lev, 2, jtrun, jtmax), &
      vornow(lev, 2, jtrun, jtmax), &
      divnow(lev, 2, jtrun, jtmax)

   real wvn_top(ktop + 1), djt

   integer k, mode, m, mf, n, nflt, KL
   real pi, flt, fac
   integer async_id, istat, ierr
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   wvn_top(1) = max(min(jtrun/3., 155.), 55.)
   wvn_top(ktop + 1) = jtrun
   do k = 2, ktop
      djt = (wvn_top(ktop + 1) - wvn_top(1))*exp(-0.5*(k - 1))
!        wvn_top(k) = (jtrun + wvn_top(k-1))*0.5
      wvn_top(k) = min(wvn_top(ktop + 1) - djt, float(jtrun))
   end do
!
   pi = 3.141596
!
!  mode = 0 : just truncate into assigned wavenumbers without
!             extra filtering
!  mode = 1 or other : add fitering along with truncating
!
   mode = 1
!
   if (mode .eq. 0) then
      do k = 1, lev
!2dMPI >
         KL = Llist(k)
         if (KL .le. ktop) then
!2dMPI <
            do m = 1, mlistnum
               mf = max(2, mlist(m))
               do n = mf, jtrun
                  nflt = int(wvn_top(k)/float(n))
                  flt = min(1.0, float(nflt))
                  vornow(k, 1, n, m) = vornow(k, 1, n, m)*flt
                  divnow(k, 1, n, m) = divnow(k, 1, n, m)*flt
                  temnow(k, 1, n, m) = temnow(k, 1, n, m)*flt
                  vornow(k, 2, n, m) = vornow(k, 2, n, m)*flt
                  divnow(k, 2, n, m) = divnow(k, 2, n, m)*flt
                  temnow(k, 2, n, m) = temnow(k, 2, n, m)*flt
               end do
            end do
!          do m = 1, mlistnum
!            mf=max(2,mlist(m))
!            do n = mf, jtrun
!              nflt = int ( wvn_top(k) / float(n) )
!              flt = min( 1.0, float(nflt) )
!              do nc=1,ncld
!                kk=(nc-1)*lev+k
!                qnow(kk,1,n,m)  = qnow(kk,1,n,m)*flt
!                qnow(kk,2,n,m)  = qnow(kk,2,n,m)*flt
!              enddo
!            enddo
!          enddo
!2dMPI >
         end if
!2dMPI <
      end do
   else
      !$acc enter data copyin(wvn_top) async(async_id)
      !$acc parallel loop gang async(async_id) private(KL)
      do k = 1, lev
!2dMPI >
         KL = Llist(k)
         if (KL .le. ktop) then
!2dMPI <
            !$acc loop worker private(mf)
            do m = 1, mlistnum
               mf = max(2, mlist(m))
               !$acc loop vector private(fac, flt)
               do n = mf, jtrun
                  fac = min(wvn_top(k), float(n - 1))*pi/wvn_top(k)
                  flt = sin(fac)/fac
                  vornow(k, 1, n, m) = vornow(k, 1, n, m)*flt
                  divnow(k, 1, n, m) = divnow(k, 1, n, m)*flt
                  temnow(k, 1, n, m) = temnow(k, 1, n, m)*flt
                  vornow(k, 2, n, m) = vornow(k, 2, n, m)*flt
                  divnow(k, 2, n, m) = divnow(k, 2, n, m)*flt
                  temnow(k, 2, n, m) = temnow(k, 2, n, m)*flt
               end do
            end do
!          do m = 1, mlistnum
!            mf=max(2,mlist(m))
!            do n = mf, jtrun
!              fac = min(wvn_top(k),float(n-1)) * pi / wvn_top(k)
!              flt = sin(fac)/fac
!              do nc=1,ncld
!                kk=(nc-1)*lev+k
!                qnow(kk,1,n,m)  = qnow(kk,1,n,m)*flt
!                qnow(kk,2,n,m)  = qnow(kk,2,n,m)*flt
!              enddo
!            enddo
!          enddo
!2dMPI >
         end if
!2dMPI <
      end do
      !$acc exit data delete(wvn_top) async(async_id)
   end if
!
   return
end

      subroutine filter_skeb_gpu(jtrun,jtmax,lev,dissest,wvn_top, async_id)
!
!  apply Lanczos filter to estimation of kinetic energy dissipation.
!
      use index
      use mpe
      use const, only : RTYPE
!
      implicit  none

!
      integer   jtrun,jtmax,lev,ncld

      real(kind=RTYPE) dissest(lev,2,jtrun,jtmax)
!
      real      wvn_top

      integer   k,mode,m,mf,n,nflt,IERR,KL, async_id
      real      pi,flt,fac
!

!2dMPI >
!     if(ktop.gt.levp)then
!        print *,'filter_top fatal: ktop greater than lev partial !'
!        call MPI_FINALIZE(IERR)
!        stop
!     endif
!2dMPI <

!
      pi = 3.141596
!
!  mode = 0 : just truncate into assigned wavenumbers without 
!             extra filtering
!  mode = 1 or other : add fitering along with truncating
!
      mode = 1
!
      if( mode .eq. 0 ) then
         !$acc parallel loop gang collapse(2) private(KL, mf) async(async_id)
         do k = 1, lev
!2dMPI >
            do m = 1, mlistnum
               KL=Llist(k)
               mf=max(2,mlist(m))
               !$acc loop vector private(nflt, flt)
               do n = mf, jtrun
                  nflt = int ( wvn_top / float(n) )
                  flt = min( 1.0, float(nflt) )
                  dissest(k,1,n,m)= dissest(k,1,n,m)*flt
                  dissest(k,2,n,m)= dissest(k,2,n,m)*flt
               enddo
            enddo
!2dMPI <
         enddo
      else
         !$acc parallel loop gang collapse(2) private(KL, mf) async(async_id)
         do k = 1, lev
!2dMPI >
            do m = 1, mlistnum
               KL=Llist(k)
               mf=max(2,mlist(m))
               !$acc loop vector private(fac, flt)
               do n = mf, jtrun
                  fac = min(wvn_top,float(n-1)) * pi / wvn_top
                  flt = sin(fac)/fac
                  dissest(k,1,n,m)= dissest(k,1,n,m)*flt
                  dissest(k,2,n,m)= dissest(k,2,n,m)*flt
               enddo
            enddo
!2dMPI <
         enddo
      endif
!  
      return
      end
