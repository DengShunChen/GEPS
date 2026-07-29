#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

      subroutine tranuv1_gpu(jtrun, jtmax, nx, my, my_max, lev, onocos, wcfac &
                             , wdfac, poly, dpoly, vor, div, ut, vt, nsize)
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
         use const, only: RTYPE
         use index
!     use paramt
         use fftcom
         use spec_cuda_graph, only: fft_cg => tranuv1_fft_cg, lt_cg => tranuv1_lt_cg, &
            cc => cc_cg2, gwk1 => gwk1_cg2, fj_ws3 => fj_ws3_cg, fj_wc => fj_wc_cg, fj_tcc => fj_tcc_cg
         use openacc
         use cudafor
         use cublas

         implicit none

         integer jtrun, jtmax, nx, my, my_max, lev, nsize
         integer myhalf, lev2, mlx, j, m, mf, l, k, nb, jchk
         integer jje, jlistnum_fj, l_fj, j_fj, llistnum_fj
         integer kk, ll, jj, jx, j2, i, jtrunj, mchk, mm, mp, mlst
         integer mm1, mp1, mlst1, mm2, mp2, mlst2, mm3, mp3, mlst3, nxj, ierr

         real sa00, sa10, sa20, sa30, dummy
!
         real(kind=RTYPE) onocos(my), wcfac(jtrun, jtmax), wdfac(jtrun, jtmax) &
            , poly(jtrun, my/2, jtmax), dpoly(jtrun, my/2, jtmax) &
            , vor(jtrun, jtmax, 2), div(jtrun, jtmax, 2) &
            , ut(nxp, my_max), vt(nxp, my_max)
!
         real(kind=RTYPE) wcc_fk(lev, 2, 2, jtmax, my_max*nsize)
         real(kind=RTYPE) twcc_fk(lev, 2, 2, jtmax*nsize, my_max)
         real(kind=RTYPE) tcc(lev, 2, 2, my, mlistnum)
         real(kind=RTYPE) ws3(lev, 2, 2, jtrun, mlistnum)
         real(kind=RTYPE) ws4(lev, 2, 2, jtrun, mlistnum)
!
         real(kind=RTYPE) tc2(lev, 2, 2, my, mlistnum)
         real(kind=RTYPE) wc(jtrun, my/2, mlistnum), wd(jtrun, my/2, mlistnum)
         integer(kind=cuda_stream_kind) :: stream, lt_cg_stream(jtmax)
         type(cudaEvent) :: spread_event, pack_event
         type(cublashandle) :: handle


!CWB2015
!     real      coslr(jm)
!     save coslr
         real(kind=RTYPE), dimension(my) ::  coslr

!
         logical lfirst
         data lfirst/.true./
         save lfirst
!
         integer jlist_fj(my/2, mlistnum)
         real fj_ws4(lev*2*2, jtrun, mlistnum)
         real fj_wd(jtrun, my/2, mlistnum)
         real fj_tc2(lev*2*2, my, mlistnum)
         integer :: async_id = 1, istat
!

!CWBinit
         stream = acc_get_cuda_stream(async_id)


         myhalf = my/2
         lev2 = lev*2
         mlx = (jtrun/2)*((jtrun + 1)/2)
         nb = 32
         jchk = iand(myhalf, 1)
         jje = myhalf - jchk
!
         !$acc data create(wcc_fk, wc, wd, ws3, ws4, tcc, &
         !$acc&     tc2, twcc_fk, coslr, jlist_fj) async(async_id)
         !$acc parallel loop async(async_id)
         do j = 1, my
            coslr(j) = 1./onocos(j)
         end do
         !$acc host_data use_device(wcc_fk, tcc, fj_ws3, fj_tcc, fj_wc)
         istat = cudaMemsetAsync(wcc_fk, real(0.0, RTYPE), size(wcc_fk), stream)
         istat = cudaMemsetAsync(tcc, real(0.0, RTYPE), size(tcc), stream)
         istat = cudaMemsetAsync(fj_ws3, real(0.0), size(fj_ws3), stream)
         istat = cudaMemsetAsync(fj_tcc, real(0.0), size(fj_tcc), stream)
         istat = cudaMemsetAsync(fj_wc, real(0.0), size(fj_wc), stream)
         !$acc end host_data


         !$acc parallel loop gang collapse(2) private(mf) async(async_id)
         do m = 1, mlistnum
            do j = 1, myhalf
               mf = mlist(m)
               if (mf .le. mtrundef(j)) then
                  !$acc loop vector
                  do l = mf, jtrun
                     wc(l, j, m) = wcfac(l, m)*poly(l, j, m)
                     wd(l, j, m) = wdfac(l, m)*dpoly(l, j, m)*coslr(j)
                  end do
               end if
            end do
         end do
!
         !$acc parallel loop collapse(3) private(mf) async(async_id)
         do m = 1, mlistnum
            do l = 1, jtrun
               do k = 1, lev
                  mf = mlist(m)
                  if (l .ge. mf) then
                     ws3(k, 1, 1, l, m) = +div(l, m, 2)
                     ws3(k, 2, 1, l, m) = -div(l, m, 1)
                     ws3(k, 1, 2, l, m) = +vor(l, m, 2)
                     ws3(k, 2, 2, l, m) = -vor(l, m, 1)
                     ws4(k, 1, 1, l, m) = +vor(l, m, 1)
                     ws4(k, 2, 1, l, m) = +vor(l, m, 2)
                     ws4(k, 1, 2, l, m) = -div(l, m, 1)
                     ws4(k, 2, 2, l, m) = -div(l, m, 2)
                  end if
               end do
            end do
         end do

         !$acc parallel loop private(mf, jlistnum_fj) async(async_id)
         do m = 1, mlistnum
            mf = mlist(m)
            jlistnum_fj = 0
            !$acc loop seq
            do j = 1, jje
               if (mf .le. mtrundef(j)) then
                  jlistnum_fj = jlistnum_fj + 1
                  jlist_fj(jlistnum_fj, m) = j
               end if
            end do
         end do
         !$acc parallel loop gang private(mf) async(async_id)
         do m = 1, mlistnum
            mf = mlist(m)
            !$acc loop vector private(l_fj)
            do l = mf, jtrun
               l_fj = l - mf + 1
               !$acc loop seq
               do k = 1, lev2*2
                  fj_ws3(k, l_fj, m) = ws3(k, 1, 1, l, m)
               end do
            end do
         end do

         !$acc parallel loop gang collapse(2) private(jlistnum_fj) async(async_id)
         do m = 1, mlistnum
            do l = 1, jtrun
               jlistnum_fj = tcolt_jlist(1, m)
               !$acc loop vector private(mf, j, l_fj)
               do j_fj = 1, jlistnum_fj
                  mf = mlist(m)
                  if (l .ge. mf) then
                     j = jlist_fj(j_fj, m)
                     l_fj = l - mf + 1
                     fj_wc(l_fj, j_fj, m) = wc(l, j, m)
                  end if
               end do
            end do
         end do
         if (.not. (lt_cg%created)) then
            istat = cudaEventCreate(pack_event)
            istat = cudaEventCreate(spread_event)
            handle = cublasGetHandle()
            do m = 1, mlistnum
               lt_cg_stream(m) = acc_get_cuda_stream(m + 1)
            end do

            CUDACHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal))
            istat = cudaEventRecord(spread_event, stream)
            do m = 1, mlistnum
               istat = cudaStreamWaitEvent(lt_cg_stream(m), spread_event, 0)
               mf = mlist(m)
               jlistnum_fj = tcolt_jlist(1, m)
               llistnum_fj = jtrun - mf + 1
               istat = cublasSetStream(handle, lt_cg_stream(m))
               !$acc host_data use_device(fj_ws3, fj_wc, fj_tcc)
               call dgemm('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, &
                        1.0d+0, fj_ws3(:,:,m), lev*2*2, fj_wc(:,:,m), &
                        jtrun, 1.0d+0, fj_tcc(:,:,m), lev*2*2)
               !$acc end host_data
               istat = cudaEventRecord(pack_event, lt_cg_stream(m))
               istat = cudaStreamWaitEvent(stream, pack_event, 0)
            end do
            CUDACHECK(cudaStreamEndCapture(stream, lt_cg%graph))
            CUDACHECK(cudaGraphInstantiate(lt_cg%graph_exec, lt_cg%graph, 0))
            lt_cg%created = .true.
            istat = cudaEventDestroy(pack_event)
            istat = cudaEventDestroy(spread_event)
         end if
         CUDACHECK(cudaGraphLaunch(lt_cg%graph_exec, stream))

         !$acc parallel loop collapse(2) gang private(jlistnum_fj) async(async_id)
         do m = 1, mlistnum
            do k = 1, lev2*2
               jlistnum_fj = tcolt_jlist(1, m)
               !$acc loop vector private(j)
               do j_fj = 1, jlistnum_fj
                  j = jlist_fj(j_fj, m)
                  tcc(k, 1, 1, j, m) = tcc(k, 1, 1, j, m) + fj_tcc(k, j_fj, m)
               end do
            end do
         end do

         !$acc parallel loop collapse(3) async(async_id)
         do m = 1, mlistnum
            do l = 1, jtrun
               do k = 1, lev2*2
                  mf = mlist(m)
                  if (l .ge. mf) then
                     l_fj = l - mf + 1
                     fj_ws3(k, l_fj, m) = ws4(k, 1, 1, l, m)
                  end if
               end do
            end do
         end do
         
         !$acc parallel loop gang collapse(2) private(mf, jlistnum_fj) async(async_id)
         do m = 1, mlistnum
            do l = 1, jtrun
               mf = mlist(m)
               if (l .ge. mf) then
                  jlistnum_fj = tcolt_jlist(1, m)
                  !$acc loop vector private(j, l_fj)
                  do j_fj = 1, jlistnum_fj
                     j = jlist_fj(j_fj, m)
                     l_fj = l - mf + 1
                     fj_wc(l_fj, j_fj, m) = wd(l, j, m)
                  end do
               end if
            end do
         end do
         !$acc host_data use_device(fj_tcc)
         istat = cudaMemsetAsync(fj_tcc, real(0.0), size(fj_tcc), stream)
         !$acc end host_data

         !do m = 1, mlistnum
         !   mf = mlist(m)
         !   jlistnum_fj = tcolt_jlist(1, m)
         !   llistnum_fj = jtrun - mf + 1
         !   istat = cublasSetStream(handle, lt_cg_stream(m))
         !   !$acc host_data use_device(fj_ws3, fj_wc, fj_tcc)
         !   call dgemm('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, &
         !              fj_ws3(:,:,m), lev*2*2, fj_wc(:,:,m), &
         !              jtrun, 1.0d+0, fj_tcc(:,:,m), lev*2*2)
         !   !$acc end host_data
         !end do
         CUDACHECK(cudaGraphLaunch(lt_cg%graph_exec, stream))

         !$acc parallel loop gang collapse(2) private(jlistnum_fj) async(async_id)
         do m = 1, mlistnum
            do k = 1, lev2*2
               jlistnum_fj = tcolt_jlist(1, m)
               !$acc loop vector private(j)
               do j_fj = 1, jlistnum_fj
                  j = jlist_fj(j_fj, m)
                  tcc(k, 1, 1, j, m) = tcc(k, 1, 1, j, m) + fj_tcc(k, j_fj, m)
               end do
            end do
         end do

!
! odd number
!
         if (jchk .eq. 1) then
            j = myhalf
            !$acc parallel loop gang collapse(2) private(mf) async(async_id)
            do m = 1, mlistnum
               do kk = 1, lev2*2, nb
                  mf = mlist(m)
                  if (mf .le. mtrundef(j)) then
                     !$acc loop vector private(sa00, sa10, sa20, sa30)
                     do k = kk, min(kk + nb - 1, lev2*2), 4
                        !$acc loop seq
                        do ll = mf, jtrun, nb
                           sa00 = tcc(k, 1, 1, j, m)
                           sa10 = tcc(k + 1, 1, 1, j, m)
                           sa20 = tcc(k + 2, 1, 1, j, m)
                           sa30 = tcc(k + 3, 1, 1, j, m)
                           !$acc loop seq
                           do l = ll, min(ll + nb - 1, jtrun)
                              sa00 = sa00 + ws3(k, 1, 1, l, m)*wc(l, j, m) + ws4(k, 1, 1, l, m)*wd(l, j, m)
                              sa10 = sa10 + ws3(k + 1, 1, 1, l, m)*wc(l, j, m) + ws4(k + 1, 1, 1, l, m)*wd(l, j, m)
                              sa20 = sa20 + ws3(k + 2, 1, 1, l, m)*wc(l, j, m) + ws4(k + 2, 1, 1, l, m)*wd(l, j, m)
                              sa30 = sa30 + ws3(k + 3, 1, 1, l, m)*wc(l, j, m) + ws4(k + 3, 1, 1, l, m)*wd(l, j, m)
                           end do
                           tcc(k, 1, 1, j, m) = sa00
                           tcc(k + 1, 1, 1, j, m) = sa10
                           tcc(k + 2, 1, 1, j, m) = sa20
                           tcc(k + 3, 1, 1, j, m) = sa30
                        end do
                     end do
                  end if
               end do
            end do
         end if
         !$acc parallel loop gang collapse(2) private(mf) async(async_id)
         do m = 1, mlistnum
            do k = 1, lev
               mf = mlist(m)
               !$acc loop vector
               do l = mf, jtrun, 2
                  ws3(k, 1, 1, l, m) = +div(l, m, 2)
                  ws3(k, 2, 1, l, m) = -div(l, m, 1)
                  ws3(k, 1, 2, l, m) = +vor(l, m, 2)
                  ws3(k, 2, 2, l, m) = -vor(l, m, 1)
                  ws4(k, 1, 1, l, m) = -vor(l, m, 1)
                  ws4(k, 2, 1, l, m) = -vor(l, m, 2)
                  ws4(k, 1, 2, l, m) = +div(l, m, 1)
                  ws4(k, 2, 2, l, m) = +div(l, m, 2)
               end do
            end do
         end do
         !$acc parallel loop gang collapse(2) private(mf) async(async_id)
         do m = 1, mlistnum
            do k = 1, lev
               mf = mlist(m)
               !$acc loop vector
               do l = mf + 1, jtrun, 2
                  ws3(k, 1, 1, l, m) = -div(l, m, 2)
                  ws3(k, 2, 1, l, m) = +div(l, m, 1)
                  ws3(k, 1, 2, l, m) = -vor(l, m, 2)
                  ws3(k, 2, 2, l, m) = +vor(l, m, 1)
                  ws4(k, 1, 1, l, m) = +vor(l, m, 1)
                  ws4(k, 2, 1, l, m) = +vor(l, m, 2)
                  ws4(k, 1, 2, l, m) = -div(l, m, 1)
                  ws4(k, 2, 2, l, m) = -div(l, m, 2)
               end do
            end do
         end do
!
         nb = 32
         jchk = iand(myhalf, 1)
         !$acc host_data use_device(tc2, fj_ws3, fj_tcc, fj_wc)
         istat = cudaMemsetAsync(tc2, real(0.0, RTYPE), size(tc2), stream)
         istat = cudaMemsetAsync(fj_ws3, real(0.0), size(fj_ws3), stream)
         istat = cudaMemsetAsync(fj_tcc, real(0.0), size(fj_tcc), stream)
         istat = cudaMemsetAsync(fj_wc, real(0.0), size(fj_wc), stream)
         !$acc end host_data
         
         !$acc parallel loop collapse(3) async(async_id)
         do m = 1, mlistnum
            do l = 1, jtrun
               do k = 1, lev2*2
                  mf = mlist(m)
                  if (l .ge. mf) then
                     l_fj = l - mf + 1
                     fj_ws3(k, l_fj, m) = ws3(k, 1, 1, l, m)
                  end if
               end do
            end do
         end do
         !$acc parallel loop collapse(2) private(mf, jlistnum_fj) async(async_id)
         do m = 1, mlistnum
            do l = 1, jtrun
               mf = mlist(m)
               jlistnum_fj = tcolt_jlist(1, m)
               !$acc loop vector private(j, l_fj)
               do j_fj = 1, jlistnum_fj
                  if (l .ge. mf) then
                     j = jlist_fj(j_fj, m)
                     l_fj = l - mf + 1
                     fj_wc(l_fj, j_fj, m) = wc(l, j, m)
                  end if
               end do
            end do
         end do
         !do m = 1, mlistnum
         !   mf = mlist(m)
         !   jlistnum_fj = tcolt_jlist(1, m)
         !   llistnum_fj = jtrun - mf + 1
         !   istat = cublasSetStream(handle, lt_cg_stream(m))
         !   !$acc host_data use_device(fj_ws3, fj_wc, fj_tcc)
         !   call dgemm('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, &
         !              fj_ws3(:,:,m), lev*2*2, fj_wc(:,:,m), &
         !              jtrun, 1.0d+0, fj_tcc(:,:,m), lev*2*2)
         !   !$acc end host_data
         !end do
         CUDACHECK(cudaGraphLaunch(lt_cg%graph_exec, stream))
         !$acc parallel loop gang collapse(2) private(jlistnum_fj) async(async_id)
         do m = 1, mlistnum
            do k = 1, lev2*2
               jlistnum_fj = tcolt_jlist(1, m)
               !$acc loop vector private(j)
               do j_fj = 1, jlistnum_fj
                  j = jlist_fj(j_fj, m)
                  tc2(k, 1, 1, j, m) = tc2(k, 1, 1, j, m) + fj_tcc(k, j_fj, m)
               end do
            end do
         end do
         !$acc parallel loop collapse(3) async(async_id)
         do m = 1, mlistnum
            do l = 1, jtrun
               do k = 1, lev2*2
                  mf = mlist(m)
                  if (l .ge. mf) then
                     l_fj = l - mf + 1
                     fj_ws3(k, l_fj, m) = ws4(k, 1, 1, l, m)
                  end if
               end do
            end do
         end do
         !$acc parallel loop gang collapse(2) private(mf, jlistnum_fj) async(async_id)
         do m = 1, mlistnum
            do l = 1, jtrun
               mf = mlist(m)
               if (l .ge. mf) then
                  jlistnum_fj = tcolt_jlist(1, m)
                  !$acc loop vector private(j, l_fj)
                  do j_fj = 1, jlistnum_fj
                     j = jlist_fj(j_fj, m)
                     l_fj = l - mf + 1
                     fj_wc(l_fj, j_fj, m) = wd(l, j, m)
                  end do
               end if
            end do
         end do
         !$acc host_data use_device(fj_tcc)
         istat = cudaMemsetAsync(fj_tcc, real(0.0), size(fj_tcc), stream)
         !$acc end host_data

         !do m = 1, mlistnum
         !   mf = mlist(m)
         !   jlistnum_fj = tcolt_jlist(1, m)
         !   llistnum_fj = jtrun - mf + 1
         !   istat = cublasSetStream(handle, lt_cg_stream(m))
         !   !$acc host_data use_device(fj_ws3, fj_wc, fj_tcc)
         !   call dgemm('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, &
         !              fj_ws3(:,:,m), lev*2*2, fj_wc(:,:,m), &
         !              jtrun, 1.0d+0, fj_tcc(:,:,m), lev*2*2)
         !   !$acc end host_data
         !end do
         CUDACHECK(cudaGraphLaunch(lt_cg%graph_exec, stream))
         
         !$acc parallel loop gang collapse(2) private(jlistnum_fj) async(async_id)
         do m = 1, mlistnum
            do k = 1, lev2*2
               jlistnum_fj = tcolt_jlist(1, m)
               !$acc loop vector private(j)
               do j_fj = 1, jlistnum_fj
                  j = jlist_fj(j_fj, m)
                  tc2(k, 1, 1, j, m) = tc2(k, 1, 1, j, m) + fj_tcc(k, j_fj, m)
               end do
            end do
         end do
         
         
!
! odd number
!
         if (jchk .eq. 1) then
            j = myhalf
            !$acc parallel loop gang collapse(2) private(mf) async(async_id)
            do m = 1, mlistnum
               do kk = 1, lev2*2, nb
                  mf = mlist(m)
                  if (mf .le. mtrundef(j)) then
                     !$acc loop vector private(sa00, sa10, sa20, sa30)
                     do k = kk, min(kk + nb - 1, lev2*2), 4
                        !$acc loop seq
                        do ll = mf, jtrun, nb
                           sa00 = tc2(k, 1, 1, j, m)
                           sa10 = tc2(k + 1, 1, 1, j, m)
                           sa20 = tc2(k + 2, 1, 1, j, m)
                           sa30 = tc2(k + 3, 1, 1, j, m)
                           !$acc loop seq
                           do l = ll, min(ll + nb - 1, jtrun)
                              sa00 = sa00 + ws3(k, 1, 1, l, m)*wc(l, j, m) + ws4(k, 1, 1, l, m)*wd(l, j, m)
                              sa10 = sa10 + ws3(k + 1, 1, 1, l, m)*wc(l, j, m) + ws4(k + 1, 1, 1, l, m)*wd(l, j, m)
                              sa20 = sa20 + ws3(k + 2, 1, 1, l, m)*wc(l, j, m) + ws4(k + 2, 1, 1, l, m)*wd(l, j, m)
                              sa30 = sa30 + ws3(k + 3, 1, 1, l, m)*wc(l, j, m) + ws4(k + 3, 1, 1, l, m)*wd(l, j, m)
                           end do
                           tc2(k, 1, 1, j, m) = sa00
                           tc2(k + 1, 1, 1, j, m) = sa10
                           tc2(k + 2, 1, 1, j, m) = sa20
                           tc2(k + 3, 1, 1, j, m) = sa30
                        end do
                     end do
                  end if
               end do
            end do
         end if
         !$acc parallel loop gang collapse(3) private(jj, jx, j2) async(async_id)
         do m = 1, mlistnum
            do j = 1, myhalf
               do k = 1, lev*2*2
                  jj = jlist2(j)
                  jx = my - j + 1
                  j2 = jlist2(jx)
                  wcc_fk(k, 1, 1, m, jj) = tcc(k, 1, 1, j, m)
                  wcc_fk(k, 1, 1, m, j2) = tc2(k, 1, 1, j, m)
               end do
            end do
         end do   ! end of big m loop

!ch   call mpe_transpose_sr(wcc_fk,twcc_fk,lev*2*2,jtmax,my_max,nsize)
!     call mpe_transpose_sr(wcc_fk,twcc_fk,lev*2*2,jtmax,my_max,nsize,col_comm)
         call mpe_transpose_sr_sp_gpu(wcc_fk, twcc_fk, lev*2*2, jtmax, my_max, nsize, nccl_col_comm)
         !$acc parallel loop collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, (nx + 2)*lev*2
               gwk1(i, 1, 1, jj) = 0.
            end do
         end do
!
         !$acc parallel loop gang collapse(2) private(j, jtrunj, mchk) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               j = jlist1(jj)
               jtrunj = mtrundef(j)
               mchk = iand(jtrunj, 3)
               !$acc loop vector private(mm, mp, mlst)
               do m = 1, mchk
                  mm = 2*m - 1
                  mp = mm + 1
                  mlst = nlist(m)
                  gwk1(mm, k, 1, jj) = twcc_fk(k, 1, 1, mlst, jj)
                  gwk1(mp, k, 1, jj) = twcc_fk(k, 2, 1, mlst, jj)
                  gwk1(mm, k, 2, jj) = twcc_fk(k, 1, 2, mlst, jj)
                  gwk1(mp, k, 2, jj) = twcc_fk(k, 2, 2, mlst, jj)
               end do
            end do
         end do

         !$acc parallel loop gang collapse(2) private(j, jtrunj, mchk) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               j = jlist1(jj)
               jtrunj = mtrundef(j)
               mchk = iand(jtrunj, 3)
               !$acc loop vector private(mm, mp, mlst, mm1, mp1, mlst1, mm2, mp2, &
               !$acc&     mlst2, mm3, mp3, mlst3)
               do m = mchk + 1, jtrunj, 4
                  mm = 2*m - 1
                  mp = mm + 1
                  mlst = nlist(m)
                  mm1 = 2*(m + 1) - 1
                  mp1 = mm1 + 1
                  mlst1 = nlist(m + 1)
                  mm2 = 2*(m + 2) - 1
                  mp2 = mm2 + 1
                  mlst2 = nlist(m + 2)
                  mm3 = 2*(m + 3) - 1
                  mp3 = mm3 + 1
                  mlst3 = nlist(m + 3)
                  gwk1(mm, k, 1, jj) = twcc_fk(k, 1, 1, mlst, jj)
                  gwk1(mp, k, 1, jj) = twcc_fk(k, 2, 1, mlst, jj)
                  gwk1(mm, k, 2, jj) = twcc_fk(k, 1, 2, mlst, jj)
                  gwk1(mp, k, 2, jj) = twcc_fk(k, 2, 2, mlst, jj)
                  gwk1(mm1, k, 1, jj) = twcc_fk(k, 1, 1, mlst1, jj)
                  gwk1(mp1, k, 1, jj) = twcc_fk(k, 2, 1, mlst1, jj)
                  gwk1(mm1, k, 2, jj) = twcc_fk(k, 1, 2, mlst1, jj)
                  gwk1(mp1, k, 2, jj) = twcc_fk(k, 2, 2, mlst1, jj)
                  gwk1(mm2, k, 1, jj) = twcc_fk(k, 1, 1, mlst2, jj)
                  gwk1(mp2, k, 1, jj) = twcc_fk(k, 2, 1, mlst2, jj)
                  gwk1(mm2, k, 2, jj) = twcc_fk(k, 1, 2, mlst2, jj)
                  gwk1(mp2, k, 2, jj) = twcc_fk(k, 2, 2, mlst2, jj)
                  gwk1(mm3, k, 1, jj) = twcc_fk(k, 1, 1, mlst3, jj)
                  gwk1(mp3, k, 1, jj) = twcc_fk(k, 2, 1, mlst3, jj)
                  gwk1(mm3, k, 2, jj) = twcc_fk(k, 1, 2, mlst3, jj)
                  gwk1(mp3, k, 2, jj) = twcc_fk(k, 2, 2, mlst3, jj)
               end do
            end do
         end do
!
         if (length_fft .eq. 0 .and. lreduce .eq. 0) then
            call rfftmlt_gpu(cc, gwk1, trigs, ifax, 1, nx + 2, nx, jlistnum*lev*2, 1)
         else

            if (.not. (fft_cg%created)) then
               call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, &
                  nxdef, jlistnum, nx + 2, lev*2, 1, fft_cg%graph)
               CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
               fft_cg%created = .true.
            end if
            CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
            CUDACHECK(cudaStreamSynchronize(stream))
            
         end if
!
         !$acc parallel loop collapse(2) private(j) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nxp
               j = jlist1(jj)
               if (i .le. nxjlen(j)) then
                  ut(i, jj) = cc(nxjstart(j) - 1 + i, 1, 1, jj)
                  vt(i, jj) = cc(nxjstart(j) - 1 + i, 1, 2, jj)
               end if
            end do
         end do
         !$acc end data
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
