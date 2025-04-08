      subroutine vartran_gpu(vorten, divten, phiten, nw, jtrun, jtmax, lev, &
                             rad, omega, h, iflag)
!
!    purpose : conversion between dimension and dimensionless forms
!--------------------------------------------------------------------
! *** input ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!  nw      :  total wavenumber index array
!  lev     :  total vertical levels
!  k       :  vetical level of (non)dimensionlize
!  rad     :  radius of earth
!  omega   :  angular velocity of earth
!  h       :  coefficient array of (non)dimensionlize
!  iflag   : +2 = to nondimension form
!          : -2 = to dimensional form
! **** output ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!************************************
         use index, only: levf, Lstart, Lend, mlist, mlistnum
         use const, only: RTYPE, nnmivm

         implicit none

         integer jtrun, jtmax, lev, k, iflag, m, l, mf, kk

         real(kind=RTYPE) vorten(lev, 2, jtrun, jtmax), &
            divten(lev, 2, jtrun, jtmax), &
            phiten(lev, 2, jtrun, jtmax)
         real h(jtrun, jtmax, nnmivm) ! 2dMPI
         integer nw(jtrun, jtmax)

         real omega2, omega, omga2r2, rad, temp1, temp2

         omega2 = omega*omega
         omga2r2 = omega2*rad*rad
!
!  do nondimensionlize which are tendency
!
         if (iflag .eq. 2) then

            do m = 1, mlistnum
               mf = mlist(m)
               do l = mf, jtrun
                  do kk = 1, nnmivm
                     if ((kk .ge. Lstart) .and. (kk .le. Lend)) then
                        k = kk - Lstart + 1
                        vorten(k, 1, l, m) = -h(l, m, kk)*vorten(k, 1, l, m)/omega2
                        vorten(k, 2, l, m) = -h(l, m, kk)*vorten(k, 2, l, m)/omega2

                        temp1 = divten(k, 2, l, m)/omega2
                        temp2 = divten(k, 1, l, m)/omega2

                        divten(k, 1, l, m) = -h(l, m, kk)*temp1
                        divten(k, 2, l, m) = h(l, m, kk)*temp2

                        phiten(k, 1, l, m) = phiten(k, 1, l, m)/omga2r2/omega
                        phiten(k, 2, l, m) = phiten(k, 2, l, m)/omga2r2/omega
                     end if
                  end do
               end do
            end do
!
!  do dimensionlize which are correct term now
!
         else if (iflag .eq. -2) then

            do m = 1, mlistnum
               mf = mlist(m)
               do l = mf, jtrun
                  if (l .ne. 1) then
                     do kk = 1, nnmivm
                        if ((kk .ge. Lstart) .and. (kk .le. Lend)) then
                           k = kk - Lstart + 1
                           vorten(k, 1, l, m) = -vorten(k, 1, l, m)*omega/h(l, m, kk)
                           vorten(k, 2, l, m) = -vorten(k, 2, l, m)*omega/h(l, m, kk)
                           temp1 = divten(k, 2, l, m)*omega/h(l, m, kk)
                           temp2 = -divten(k, 1, l, m)*omega/h(l, m, kk)
                           divten(k, 1, l, m) = temp1
                           divten(k, 2, l, m) = temp2
                        end if
                     end do
                  end if

                  do kk = 1, nnmivm
                     if ((kk .ge. Lstart) .and. (kk .le. Lend)) then
                        k = kk - Lstart + 1
                        phiten(k, 1, l, m) = phiten(k, 1, l, m)*omga2r2
                        phiten(k, 2, l, m) = phiten(k, 2, l, m)*omga2r2
                     end if
                  end do
               end do
            end do
            !
         end if
!
         return
      end subroutine vartran_gpu
