      subroutine couvtq ( nxj,mn,kk,ktpbl,dt,wktri,xkm,hgt,uu,beta &
                         ,dhgt,ro2)
!
!####################################################################
!
! 1. function description
!
!     prepare tri-diagonal matrix for backward time scheme of
!     the diffusion equation: du/dt = d(xkm*(du/dz))/dz
!
! 2. common block specification
!
! 3. parameter specification
!
!     mn   : number of tri-diagonal matrix (horizontal dimension)      i
!     kk   : vertical dimension                                        i
!     ktpbl: starting lveles to consider vertical mixing effects       i
!     dt   : time step                                                 i
!     wktri: output tri-diagonal matrix coef. + forcing      (mn,kk,4) o
!     xkm  : input eddy mixing coefficient at even levels    (mn,kk)   i
!     hgt  : input height at odd levels                      (mn,kk)   i
!     uu   : variable to be updated by backward time scheme  (mn,kk)   i
!     beta : parameter to control time scheme ;                        i
!            1.0 -> backward, 1/2. -> crank-nicholson, 0. -> forward.
!
! 4. local variable
!
! 5. calling modules
!
!    mixpbl
!
! 6. usage
!
!    call couvtq (mn,kk,ktpbl,dt,wktri,xkm,hgt,uu,beta)
!
! 7. modules called
!
! 8. limitation
!
! 9. date
!
!    created    nov. 1991
!    modified   nov. 1991
!    remodified mar. 1992
!
! 10. author
!       c.-s. liou/ f.-j. wang,
!       modify to f90 by C-H Lee and sort by River Chen in 2015
!
!#####################################################################
!
      use const, only: RTYPE
!
      implicit  none

      integer   nxj,mn,kk,ktpbl,kt,kt1,k,i,ik,iksm,iks,iksp
      real      alpha,beta,dt,xx,yy,zz

      real      wktri(mn,kk,4),xkm(mn,kk),hgt(mn,kk)
      real      dhgt(mn,kk),ro2(mn,kk)
      real(kind=RTYPE) uu(mn,kk)
!
      alpha = 1.0 - beta
      kt = kk - ktpbl + 1
      kt1 = kt - 1
!
      do 120 k = 1, kt1
      do 120 i = 1, nxj
      ik = (k-1)*mn + i
      iksm =(ktpbl-2)*mn + ik
      iks  =(ktpbl-1)*mn + ik
      iksp = ktpbl   *mn + ik
      xx = -dt*xkm(iksm,1)*ro2(iksm,1)/( dhgt(iks,1) *      &
                                 (hgt(iksm,1)-hgt(iks,1)) )
      yy = -dt*xkm(iks,1)*ro2(iks,1)/( dhgt(iks,1)*         &
                               (hgt(iks,1)-hgt(iksp,1)) )
      wktri(ik,1,1) = beta * xx
      wktri(ik,1,3) = beta * yy
      wktri(ik,1,2) = 1.0 - wktri(ik,1,1) - wktri(ik,1,3)
      wktri(ik,1,4) = uu(iks,1)   * ( 1.0 + alpha*(xx+yy) ) &
                   - uu(iksm,1) * alpha * xx                &
                   - uu(iksp,1) * alpha * yy
  120 continue
!
!     fix upper b.c. (move -uu(ktpbl-1)*wktri(1,1) to forcing term)
!
      do 200 i = 1, nxj
      wktri(i,1,4) = wktri(i,1,4) - uu(i,ktpbl-1)*wktri(i,1,1)
  200 continue
!
!     fix lower b.c.
!
      do 220 i = 1, nxj
      zz = - dt*xkm(i,kk-1)*ro2(i,kk-1)/( dhgt(i,kk) *       &
                                  (hgt(i,kk-1)-hgt(i,kk)) )
      wktri(i,kt,1) = beta * zz
      wktri(i,kt,3) = 0.0
      wktri(i,kt,2) = 1.0 - wktri(i,kt,1) - wktri(i,kt,3)
      wktri(i,kt,4) = uu(i,kk)   * ( 1.0 + alpha*zz )        &
                    - uu(i,kk-1) * alpha * zz
  220 continue
!
      return
      end
