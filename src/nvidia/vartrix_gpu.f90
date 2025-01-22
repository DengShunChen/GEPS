      subroutine vartrix_gpu_old(vorten, divten, phiten, lev, x, no, &
                                 rad, omega, h, iflag)
!
!  purpose: Scatter between vorten, divten, phiten and x vector, and
!           convers between dimension and dimensionless forms
!
!-----------------------------------------------------------------------
!  ****  input ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!  x       :  vector matrix( when iflag=-2)
!  lev     :  total vertical levels
!  rad     :  radius of earth
!  omega   :  angular velocity of earth
!  h       :  coefficient array of (non)dimensionlize
!  no      :  leading dimension of x array
!  iflag   : = +2 : to nondimension form and variable tendency to x vector (5.43)
!              -2 : x vector(correction) to (vors,divs,phis) and to dimensional form
!
!  **** output ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!  x       : vector matrix(when iglag=+2)
!
!  **** local ****
!  nbig    : how many  wave number for initialization
!
!  **** present on device ****
!  vorten, divten, phiten, x
!  rad, omega, h
!  mlistnum, mlist, Lstart, Lend, nnmivm
!  jtrun, jtmax, lev
!-----------------------------------------------------------------------
         use param, only: jtrun, jtmax
         use const, only: RTYPE, nnmivm
         use index, only: mlistnum, mlist, Lstart, Lend

         implicit none

         real(kind=RTYPE), intent(inout):: vorten(lev, 2, jtrun, jtmax), &
                                           divten(lev, 2, jtrun, jtmax), &
                                           phiten(lev, 2, jtrun, jtmax)
         real, intent(inout):: x(no, 2, 2, jtmax, nnmivm)
         real, intent(in):: h(jtrun, jtmax, nnmivm)
         real, intent(in):: rad, omega
         integer, intent(in) :: lev, no, iflag

         real omega2, omega3r2, temp1, temp2, h_t
         real(kind=RTYPE) :: vorten_im, vorten_re, &
                             divten_im, divten_re, &
                             phiten_im, phiten_re

         integer k, kk, l, j, n, i, m, nbig, mf
         integer async_id
         async_id = 1
         !
         omega2 = omega*omega
         omega3r2 = omega**3*rad*rad
         if (iflag .eq. 2) then
            !
            !  construct variable vector
            !
            !$acc parallel loop collapse(3) private(k, mf, nbig, n, i, j, l, h_t, vorten_re, vorten_im, divten_re, divten_im, phiten_re, phiten_im, temp1, temp2) async(async_id)
            do kk = 1, nnmivm
               do m = 1, mlistnum
                  do i = 1, jtrun/2
                     if ((kk .ge. Lstart) .and. (kk .le. Lend)) then
                        k = kk - Lstart + 1
                        mf = mlist(m)
                        nbig = jtrun - mf + 1
                        n = nbig/2
                        if (i .le. n) then
                           j = 1 + (i - 1)*3
                           l = mf + (i - 1)*2

                           h_t = h(l, m, kk)

                           vorten_re = -h_t*vorten(k, 1, l, m)/omega2
                           vorten_im = -h_t*vorten(k, 2, l, m)/omega2

                           temp1 = divten(k, 2, l, m)/omega2
                           temp2 = divten(k, 1, l, m)/omega2

                           divten_re = -h_t*temp1
                           divten_im = h_t*temp2

                           phiten_re = phiten(k, 1, l, m)/omega3r2
                           phiten_im = phiten(k, 2, l, m)/omega3r2

                           x(j, 1, 1, m, kk) = phiten_re
                           x(j, 2, 1, m, kk) = phiten_im

                           x(j + 1, 1, 1, m, kk) = divten_re
                           x(j + 1, 2, 1, m, kk) = divten_im

                           x(j, 1, 2, m, kk) = vorten_re
                           x(j, 2, 2, m, kk) = vorten_im
                           ! --------------------
                           l = l + 1
                           h_t = h(l, m, kk)

                           vorten_re = -h_t*vorten(k, 1, l, m)/omega2
                           vorten_im = -h_t*vorten(k, 2, l, m)/omega2

                           temp1 = divten(k, 2, l, m)/omega2
                           temp2 = divten(k, 1, l, m)/omega2

                           divten_re = -h_t*temp1
                           divten_im = h_t*temp2

                           phiten_re = phiten(k, 1, l, m)/omega3r2
                           phiten_im = phiten(k, 2, l, m)/omega3r2

                           x(j + 2, 1, 1, m, kk) = vorten_re
                           x(j + 2, 2, 1, m, kk) = vorten_im

                           x(j + 1, 1, 2, m, kk) = phiten_re
                           x(j + 1, 2, 2, m, kk) = phiten_im

                           x(j + 2, 1, 2, m, kk) = divten_re
                           x(j + 2, 2, 2, m, kk) = divten_im
                        end if

                        if ((mod(nbig, 2) .ne. 0) .and. (i .eq. 1)) then
                           j = 1 + n*3
                           l = mf + n*2
                           h_t = h(l, m, kk)

                           vorten_re = -h_t*vorten(k, 1, l, m)/omega2
                           vorten_im = -h_t*vorten(k, 2, l, m)/omega2

                           temp1 = divten(k, 2, l, m)/omega2
                           temp2 = divten(k, 1, l, m)/omega2

                           divten_re = -h_t*temp1
                           divten_im = h_t*temp2

                           phiten_re = phiten(k, 1, l, m)/omega3r2
                           phiten_im = phiten(k, 2, l, m)/omega3r2
                           ! symmetric case
                           x(j, 1, 1, m, kk) = phiten_re
                           x(j, 2, 1, m, kk) = phiten_im

                           x(j + 1, 1, 1, m, kk) = divten_re
                           x(j + 1, 2, 1, m, kk) = divten_im

                           !  antiammetrix case
                           x(j, 1, 2, m, kk) = vorten_re
                           x(j, 2, 2, m, kk) = vorten_im
                        end if

                     end if
                  end do
               end do
            end do

         else if (iflag .eq. -2) then
            ! inverse transform
            do kk = 1, nnmivm
               if ((kk .ge. Lstart) .and. (kk .le. Lend)) then
                  k = kk - Lstart + 1
                  do m = 1, mlistnum
                     mf = mlist(m)
                     nbig = jtrun - mf + 1
                     n = nbig/2
                     do i = 1, n
                        j = 1 + (i - 1)*3
                        l = mf + (i - 1)*2

                        ! symmetric case
                        phiten(k, 1, l, m) = x(j, 1, 1, m, kk)
                        phiten(k, 2, l, m) = x(j, 2, 1, m, kk)

                        divten(k, 1, l, m) = x(j + 1, 1, 1, m, kk)
                        divten(k, 2, l, m) = x(j + 1, 2, 1, m, kk)

                        vorten(k, 1, l + 1, m) = x(j + 2, 1, 1, m, kk)
                        vorten(k, 2, l + 1, m) = x(j + 2, 2, 1, m, kk)

                        ! antiammetrix case
                        vorten(k, 1, l, m) = x(j, 1, 2, m, kk)
                        vorten(k, 2, l, m) = x(j, 2, 2, m, kk)

                        phiten(k, 1, l + 1, m) = x(j + 1, 1, 2, m, kk)
                        phiten(k, 2, l + 1, m) = x(j + 1, 2, 2, m, kk)

                        divten(k, 1, l + 1, m) = x(j + 2, 1, 2, m, kk)
                        divten(k, 2, l + 1, m) = x(j + 2, 2, 2, m, kk)
                     end do
                     if (mod(nbig, 2) .ne. 0) then
                        j = 1 + n*3
                        l = mf + n*2

                        ! symmetric case
                        phiten(k, 1, l, m) = x(j, 1, 1, m, kk)
                        phiten(k, 2, l, m) = x(j, 2, 1, m, kk)

                        divten(k, 1, l, m) = x(j + 1, 1, 1, m, kk)
                        divten(k, 2, l, m) = x(j + 1, 2, 1, m, kk)

                        ! antiammetrix case
                        vorten(k, 1, l, m) = x(j, 1, 2, m, kk)
                        vorten(k, 2, l, m) = x(j, 2, 2, m, kk)
                     end if

                  end do
               end if
            end do
         end if
         !
         return
      end subroutine vartrix_gpu_old

      subroutine vartrix_gpu(vorten, divten, phiten, lev, x, no, &
                             rad, omega, h, iflag)
!
!  purpose: Scatter between vorten, divten, phiten and x vector, and
!           convers between dimension and dimensionless forms
!
!-----------------------------------------------------------------------
!  ****  input ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!  x       :  vector matrix( when iflag=-2)
!  lev     :  total vertical levels
!  rad     :  radius of earth
!  omega   :  angular velocity of earth
!  h       :  coefficient array of (non)dimensionlize
!  no      :  leading dimension of x array
!  iflag   : = +2 : to nondimension form and variable tendency to x vector (5.43)
!              -2 : x vector(correction) to (vors,divs,phis) and to dimensional form
!
!  **** output ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!  x       : vector matrix(when iglag=+2)
!
!  **** local ****
!  nbig    : how many  wave number for initialization
!
!  **** present on device ****
!  vorten, divten, phiten, x
!  rad, omega, h
!  mlistnum, mlist, Lstart, Lend, nnmivm
!  jtrun, jtmax, lev
!-----------------------------------------------------------------------
         use param, only: jtrun, jtmax
         use const, only: RTYPE, nnmivm
         use index, only: mlistnum, mlist, Lstart, Lend

         implicit none
         !$acc routine(vartran_loc) seq

         real(kind=RTYPE), intent(inout):: vorten(lev, 2, jtrun, jtmax), &
                                           divten(lev, 2, jtrun, jtmax), &
                                           phiten(lev, 2, jtrun, jtmax)
         real, intent(inout):: x(no, 2, 2, jtmax, nnmivm)
         real, intent(in):: h(jtrun, jtmax, nnmivm)
         real, intent(in):: rad, omega
         integer, intent(in) :: lev, no, iflag

         real omega2, omega2r2, omega3r2, temp1, temp2, h_t
         real(kind=RTYPE) :: vor_im, vor_re, &
                             div_im, div_re, &
                             phi_im, phi_re

         integer k, kk, i, l, j, m, mf
         integer async_id
         async_id = 1
         !
         omega2 = omega*omega
         omega2r2 = omega2*rad*rad
         omega3r2 = omega**3*rad*rad
         if (iflag .eq. 2) then
            !
            !  construct variable vector
            !
            !$acc parallel loop collapse(3) private(k, mf, i, j, l, &
            !$acc& h_t, vor_re, vor_im, div_re, div_im, phi_re, phi_im, &
            !$acc& temp1, temp2) async(async_id)
            do kk = 1, nnmivm
               do m = 1, mlistnum
                  do l = 1, jtrun
                     mf = mlist(m)
                     if ((l .ge. mf) &
                         .and. (kk .ge. Lstart) .and. (kk .le. Lend)) then
                        k = kk - Lstart + 1
                        j = 1 + ((l - mf)/2)*3
                        i = l - mf + 1

                        h_t = h(l, m, kk)

                        vor_re = -h_t*vorten(k, 1, l, m)/omega2
                        vor_im = -h_t*vorten(k, 2, l, m)/omega2

                        temp1 = divten(k, 2, l, m)/omega2
                        temp2 = divten(k, 1, l, m)/omega2

                        div_re = -h_t*temp1
                        div_im = h_t*temp2

                        phi_re = phiten(k, 1, l, m)/omega3r2
                        phi_im = phiten(k, 2, l, m)/omega3r2

                        if (mod(i, 2) .eq. 1) then
                           x(j, 1, 1, m, kk) = phi_re
                           x(j, 2, 1, m, kk) = phi_im

                           x(j + 1, 1, 1, m, kk) = div_re
                           x(j + 1, 2, 1, m, kk) = div_im

                           x(j, 1, 2, m, kk) = vor_re
                           x(j, 2, 2, m, kk) = vor_im
                        else
                           x(j + 2, 1, 1, m, kk) = vor_re
                           x(j + 2, 2, 1, m, kk) = vor_im

                           x(j + 1, 1, 2, m, kk) = phi_re
                           x(j + 1, 2, 2, m, kk) = phi_im

                           x(j + 2, 1, 2, m, kk) = div_re
                           x(j + 2, 2, 2, m, kk) = div_im
                        end if
                     end if
                  end do
               end do
            end do

         else if (iflag .eq. -2) then
            ! inverse transform
            !$acc parallel loop collapse(3) private(k, mf, i, j, l, &
            !$acc& h_t, vor_re, vor_im, div_re, div_im, phi_re, phi_im, &
            !$acc& temp1, temp2) &
            !$acc& copyin(omega, omega2r2) async(async_id)
            do kk = 1, nnmivm
               do m = 1, mlistnum
                  do l = 1, jtrun
                     mf = mlist(m)
                     if ((l .ge. mf) &
                         .and. (kk .ge. Lstart) .and. (kk .le. Lend)) then
                        k = kk - Lstart + 1
                        j = 1 + ((l - mf)/2)*3
                        i = l - mf + 1

                        if (mod(i, 2) .eq. 1) then
                           phi_re = x(j, 1, 1, m, kk)
                           phi_im = x(j, 2, 1, m, kk)

                           div_re = x(j + 1, 1, 1, m, kk)
                           div_im = x(j + 1, 2, 1, m, kk)

                           vor_re = x(j, 1, 2, m, kk)
                           vor_im = x(j, 2, 2, m, kk)
                        else
                           vor_re = x(j + 2, 1, 1, m, kk)
                           vor_im = x(j + 2, 2, 1, m, kk)

                           phi_re = x(j + 1, 1, 2, m, kk)
                           phi_im = x(j + 1, 2, 2, m, kk)

                           div_re = x(j + 2, 1, 2, m, kk)
                           div_im = x(j + 2, 2, 2, m, kk)
                        end if
                        !
                        !  do dimensionlize which are correct term now
                        !
                        if (l .ne. 1) then
                           h_t = omega/h(l, m, kk)
                           vor_re = -vor_re*h_t
                           vor_im = -vor_im*h_t
                           temp1 = div_im*h_t
                           temp2 = -div_re*h_t
                           div_re = temp1
                           div_im = temp2
                        end if
                        vorten(k, 1, l, m) = vor_re
                        vorten(k, 2, l, m) = vor_im

                        divten(k, 1, l, m) = div_re
                        divten(k, 2, l, m) = div_im

                        phiten(k, 1, l, m) = phi_re*omega2r2
                        phiten(k, 2, l, m) = phi_im*omega2r2
                     end if

                  end do
               end do
            end do
         end if
         !
         return
      end subroutine vartrix_gpu
