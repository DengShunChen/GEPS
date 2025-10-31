      subroutine adjptqintp_gpu(ut, vt, tt, qt, qp, pnew, pten, nxjp, nxp, my_max, &
                                lev, ncld, dta)

         use const, only: RTYPE, dsigma, sigma, qmin
         use grid, only: ndslvvar
         use param, only: my
         use index, only: jlistnum, jlist1

         implicit none

         integer nxjp(my), nxp, my_max, lev, ncld, myim(my_max)
         real dta
         real(kind=RTYPE) ut(nxp, lev, my_max), vt(nxp, lev, my_max), tt(nxp, lev, my_max), &
            qt(nxp, lev*ncld, my_max), qp(nxp, lev*ncld, my_max), &
            pnew(nxp, my_max), pten(nxp, my_max), pold(nxp, my_max)

         integer i, k, n, kk, ki, jj
         integer kuu, kvv, ktt
         real(kind=RTYPE) plnew(nxp, lev + 1, my_max), plold(nxp, lev + 1, my_max)
         real(kind=RTYPE) rqda(nxp, lev, ndslvvar, my_max), rqnn(nxp, lev, ndslvvar, my_max)
         
         ! fixend_cell_plm_intp_gpu
         real(kind=RTYPE) hfdp(nxp, lev, my_max), rdsi(nxp, lev, my_max)
         real(kind=RTYPE) locf(nxp, lev + lev, my_max), df(nxp, lev + lev - 1, my_max), &
                          dt(nxp, lev + lev - 1, my_max)
         real(kind=RTYPE) dq(nxp, lev, ndslvvar, my_max), &
                          qmi, qpi
         real(kind=RTYPE) valf(nxp, lev + lev - 1, ndslvvar, my_max)
         real(kind=RTYPE) rdd(nxp, lev, my_max)
         real(kind=RTYPE) dd, check, ss, dfr, rddr, ss1, ss2, ss3
         integer js(nxp, lev + 1, my_max), i4j(nxp, lev + lev, my_max), jsr, jsr1
         real dsigp, qtot, qtota, odpondp, ptmp, rdsir, rdsip, rqdar, rqdam, &
            rqdap, poldr, dsigma1, dsigma2, hfdpr, hfdpr1
         integer j, ip, in, m
         integer :: async_id = 1

         kuu = 1
         kvv = kuu + 1
         ktt = kvv + 1
         !$acc data create(myim, pold, plnew, plold, rqda, locf, &
         !$acc&     i4j, js, hfdp, rdsi, df, rdd, dt, dq, &
         !$acc&     valf, rqnn) async(async_id)
         !$acc parallel loop async(async_id)
         do jj = 1, jlistnum
            j = jlist1(jj)
            myim(jj) = nxjp(j)
         end do

         
         !$acc parallel loop collapse(2) async(async_id) private(ptmp, ki, &
         !$acc&         dsigp, qtot, qtota, kk, poldr)
         do jj = 1, jlistnum
            do i = 1, nxp
               if (i .le. myim(jj)) then
                  poldr = pnew(i, jj)
                  ! new surface pressure
                  plnew(i, lev + 1, jj) = 0.
                  plold(i, lev + 1, jj) = 0.
                  ptmp = 0.
                  !$acc loop seq
                  do k = 1, lev
                     ki = lev - k + 1
                     dsigp = dsigma(k, 1)*poldr + dsigma(k, 2)
                     plold(i, ki, jj) = plold(i, ki + 1, jj) - dsigp
                     qtot = 0.
                     qtota = 0.
                     !$acc loop seq
                     do n = 1, ncld
                        kk = k + (n - 1)*lev
                        qtot = qtot + qp(i, kk, jj)
                        qtota = qtota + qt(i, kk, jj)
                     end do
                     ptmp = ptmp + dsigp*(1.-qtot + qtota)
                  end do
                  pnew(i, jj) = ptmp
                  pten(i, jj) = (pnew(i, jj) - pten(i, jj))/dta
                  pold(i, jj) = poldr
               end if
            end do
         end do
         !$acc parallel loop collapse(2) async(async_id) private(ki, dsigma1, &
         !$acc&         dsigma2, odpondp, kk)
         do jj = 1, jlistnum
            do i = 1, nxp
               if (i .le. myim(jj)) then
                  ! mass adjustment of all tracers and virtual potential temperature
                  !$acc loop seq
                  do k = 1, lev
                     ki = lev - k + 1
                     dsigma1 = dsigma(k, 1)
                     dsigma2 = dsigma(k, 2)
                     dsigp = dsigma1*pnew(i, jj) + dsigma2
                     plnew(i, ki, jj) = plnew(i, ki + 1, jj) - dsigp
                     odpondp = (dsigma1*pold(i, jj) + dsigma2)/dsigp
                     rqda(i, ki, kuu, jj) = ut(i, k, jj)*odpondp
                     rqda(i, ki, kvv, jj) = vt(i, k, jj)
                     rqda(i, ki, ktt, jj) = tt(i, k, jj)
                     !$acc loop seq
                     do n = 1, ncld
                        kk = k + (n - 1)*lev
                        rqda(i, ki, ktt + n, jj) = max(qt(i, kk, jj)*odpondp, qmin)
                     end do
                  end do
               end if
            end do
         end do
         
         !$acc parallel loop collapse(2) async(async_id) private(ip, in)
         do jj = 1, jlistnum
            do i = 1, nxp
               if (i .le. myim(jj)) then
                  !call fixend_cell_plm_intp_gpu(plold, rqda, plnew, rqnn, lev, ndslvvar)
                  locf(i, 1, jj) = plold(i, 1, jj)
                  i4j(i, 1, jj) = 1
                  js(i, 1, jj) = 1
                  ip = 2
                  in = 2
                  !$acc loop seq
                  do j = 2, lev + lev - 1
                     if (plold(i, ip, jj) .le. plnew(i, in, jj)) then
                        locf(i, j, jj) = plold(i, ip, jj)
                        ip = ip + 1
                     else
                        locf(i, j, jj) = plnew(i, in, jj)
                        js(i, in, jj) = j
                        in = in + 1
                     end if
                     i4j(i, j, jj) = ip - 1
                  end do
                  locf(i, lev + lev, jj) = plold(i, lev + 1, jj)
                  js(i, lev + 1, jj) = lev + lev
               end if
            end do
         end do
      !
      ! interpolation coefficient
      !
         !$acc parallel loop gang collapse(2) private(dd) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               !$acc loop vector
               do i = 1, myim(jj)
                  dd = 0.0
                  !$acc loop seq
                  do j = js(i, k, jj), js(i, k + 1, jj) - 1
                     dfr = locf(i, j + 1, jj) - locf(i, j, jj)
                     dd = dd + dfr
                     df(i, j, jj) = dfr
                  end do
                  rdd(i, k, jj) = 1.0/dd
               end do
            end do
         end do
         !$acc parallel loop gang collapse(2) async(async_id) private(hfdpr, hfdpr1)
         do jj = 1, jlistnum
            do k = 1, lev
               !$acc loop vector
               do i = 1, myim(jj)
                  hfdpr = 0.5*(plold(i, k + 1, jj) - plold(i, k, jj))
                  hfdpr1 = 0.5*(plold(i, k, jj) - plold(i, k - 1, jj))
                  hfdp(i, k, jj) = hfdpr
                  if (k .ge. 2) rdsi(i, k, jj) = 1.0/(hfdpr + hfdpr1)
               end do
            end do
         end do
         
         !$acc parallel loop gang collapse(2) private(k) async(async_id)
         do jj = 1, jlistnum
            do j = 1, lev + lev - 1
               !$acc loop vector
               do i = 1, myim(jj)
                  k = i4j(i, j, jj)
                  dt(i, j, jj) = (locf(i, j, jj) + 0.5*df(i, j, jj)) - (plold(i, k, jj) + hfdp(i, k, jj))
               end do
            end do
         end do
      !
      ! start interpolation by integral of ppm
      !
         !$acc parallel loop collapse(3) async(async_id) private(check, k, &
         !$acc&         qpi, qmi, rdsir, rdsip, rqdar, rqdam, rqdap)
         do jj = 1, jlistnum 
            do m = 1, ndslvvar
               do i = 1, nxp
                  if (i .le. myim(jj)) then
                     dq(i, 1, m, jj) = 0.0
                     dq(i, lev, m, jj) = 0.0
                     rqdar = rqda(i, 2, m, jj)
                     rqdam = rqda(i, 1, m, jj)
                     rdsir = rdsi(i, 2, jj)
                     !$acc loop seq
                     do k = 2, lev - 1
                        rqdap = rqda(i, k + 1, m, jj)
                        rdsip = rdsi(i, k + 1, jj)
                        qmi = (rqdar - rqdam)*rdsir
                        qpi = (rqdap - rqdar)*rdsip
                        rqdam = rqdar
                        rqdar = rqdap
                        rdsir = rdsip
                        check = qmi*qpi
                        if (check .lt. 0.0) then
                           dq(i, k, m, jj) = 0.0
                        else
                           dq(i, k, m, jj) = 0.5*(qmi + qpi)
                        end if
                     end do
                     !$acc loop seq
                     do j = 1, lev + lev - 1
                        k = i4j(i, j, jj)
                        valf(i, j, m, jj) = rqda(i, k, m, jj) + dq(i, k, m, jj)*dt(i, j, jj)
                     end do
                  end if
               end do
            end do
         end do
         !$acc parallel loop gang vector collapse(3) async(async_id) private(ss1, ss2, &
         !$acc&         ss3, dfr, rddr, kk, ki)
         do jj = 1, jlistnum 
            do k = 1, lev
               do i = 1, nxp
                  if (i .le. myim(jj)) then
                     ki = lev - k + 1
                     ss1 = 0.0
                     ss2 = 0.0
                     ss3 = 0.0
                     jsr = js(i, k, jj)
                     jsr1 = js(i, k + 1, jj)
                     !$acc loop seq
                     do j = jsr, jsr1 - 1
                        dfr = df(i, j, jj)
                        ss1 = ss1 + valf(i, j, kuu, jj)*dfr
                        ss2 = ss2 + valf(i, j, kvv, jj)*dfr
                        ss3 = ss3 + valf(i, j, ktt, jj)*dfr
                     end do
                     rddr = rdd(i, k, jj)
                     ut(i, ki, jj) = ss1*rddr
                     vt(i, ki, jj) = ss2*rddr
                     tt(i, ki, jj) = ss3*rddr
                  end if
               end do
            end do
         end do
         !$acc parallel loop gang collapse(2) async(async_id) private(ss, kk, ki, jsr, jsr1)
         do jj = 1, jlistnum 
            do k = 1, lev
               ki = lev - k + 1
               !$acc loop vector
               do i = 1, myim(jj)
                  !$acc loop seq
                  do n = 1, ncld
                     jsr = js(i, k, jj)
                     jsr1 = js(i, k + 1, jj)
                     ss = 0.0
                     !$acc loop seq
                     do j = jsr, jsr1 - 1
                        ss = ss + valf(i, j, n + ktt, jj)*df(i, j, jj)
                     end do
                     kk = ki + (n - 1)*lev
                     qt(i, kk, jj) = ss*rdd(i, k, jj)
                  end do
               end do
            end do
         end do
               ! end inlined fixend_cell_plm_intp_gpu

         !$acc end data

         return
      end

! -------------------------------------------------------------------------
