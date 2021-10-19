      subroutine reduceintp(a,lonfd,lonf,latg)
!
! do  mass conserving interpolation from different grid at given latitude
!
! author: hann-ming henry juang 2008
!
      implicit none

!
      integer    lonf,latg,j,i,imp,lonfd(latg)
      real       a(lonf,latg)
!
      real      old(lonf),new(lonf)
      real      xpast(lonf+1),xnext(lonf+1)
      real      two_pi,dxp,dxf,hfdxp,hfdxf,sc,pi
!

! ..................................
! ..................................
       pi  = 4.0 * atan(1.0)
       two_pi = 2.0 * pi
       dxf = two_pi / lonf
       hfdxf = 0.5 * dxf
!$omp parallel do                                       &
!$omp private(j,i,imp,dxp,hfdxp,xpast,xnext,sc,old,new) &
!$omp schedule(dynamic)
      do j=1,latg
       imp=lonfd(j)
       if( imp.ne.lonf ) then
        dxp = two_pi / imp
        hfdxp = 0.5 * dxp


        do i=1,imp+1
          xpast(i) = (i-1) * dxp - hfdxp
        enddo

        do i=1,lonf+1
          xnext(i) = (i-1) * dxf - hfdxf
        enddo

        sc=two_pi

        old(1:imp)=a(1:imp,j)
!CWB2021 for ndsl single precision test
!       call cyclic_cell_ppm_intp(xpast,old,xnext,new,lonf,1,imp,lonf,sc)
        call cyclic_cell_ppm_intp_dp(xpast,old,xnext,new,lonf,1,imp,lonf,sc)

!        call cyclic_cell_plm_intp(xpast,old,xnext,new,lonf,1,imp,lonf,sc)

        a(1:lonf,j)=new(1:lonf)
       endif
      enddo
!$omp end parallel do
! .................

      return
      end subroutine reduceintp
