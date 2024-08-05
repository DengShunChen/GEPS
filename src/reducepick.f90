      subroutine reducepick(a,lonfd,lonf,latg)
!
! do  mass conserving interpolation from different grid at given
! latitude
!
! author: hann-ming henry juang 2008
!
      use const, only: RTYPE
!
      implicit none

!
      integer    lonf,latg,j,i,imp,lonfd(latg)
      real       a(lonf,latg)
!
      real(kind=RTYPE) old(lonf),new(lonf)
      real(kind=RTYPE) xpast(lonf+1),xnext(lonf+1)
      real(kind=RTYPE) sc
      real      two_pi,dxp,dxf,hfdxp,hfdxf,pi
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


        do i=1,lonf+1
          xpast(i) = (i-1) * dxf - hfdxf
        enddo

        do i=1,imp+1
          xnext(i) = (i-1) * dxp - hfdxp
        enddo

        sc=two_pi

        old(1:lonf)=a(1:lonf,j)
!CWB2021 for ndsl single precision test
        call cyclic_cell_ppm_intp(xpast,old,xnext,new,lonf,1,lonf,imp,sc)
!        call cyclic_cell_ppm_intp_dp(xpast,old,xnext,new,lonf,1,lonf,imp,sc)

!        call cyclic_cell_plm_intp(xpast,old,xnext,new,lonf,1,lonf,imp,sc)

        a(1:imp,j)=new(1:imp)
       endif
      enddo
!$omp end parallel do
! .................

      return
      end subroutine reducepick
!
      subroutine reducepicki(a,lonfd,lonf,latg)
!
! pick the reduce-grid value from the nearest regular grid
! then fill the tailing points to be the same as the last
! reduce-grid point.
!
      integer a(lonf,latg)
      dimension lonfd(latg)
      integer,allocatable:: tmp(:)
      allocate(tmp(lonf))
      dg=360./float(lonf)
      do j=1,latg
        dr=360./float(lonfd(j))
        do i=1,lonfd(j)
          ii=nint((i-1.)*dr/dg + 1.0)
          tmp(i)=a(ii,j)
        enddo
        ii=lonfd(j)
        do i=1,ii
          a(i,j)=tmp(i)
        enddo
        do i=ii+1,lonf
          a(i,j)=tmp(ii)
        enddo
      enddo
      deallocate(tmp)
      return
      end

      subroutine reducepickl(a,lonfd,lonf,latg)
!
! pick the reduce-grid value from the nearest regular grid
! then fill the tailing points to be the same as the last
! reduce-grid point.
!
      logical a(lonf,latg)
      dimension lonfd(latg)
      logical,allocatable:: tmp(:)
      allocate(tmp(lonf))
      dg=360./float(lonf)
      do j=1,latg
        dr=360./float(lonfd(j))
        do i=1,lonfd(j)
          ii=nint((i-1.)*dr/dg + 1.0)
          tmp(i)=a(ii,j)
        enddo
        ii=lonfd(j)
        do i=1,ii
          a(i,j)=tmp(i)
        enddo
        do i=ii+1,lonf
          a(i,j)=tmp(ii)
        enddo
      enddo
      deallocate(tmp)
      return
      end
!
      subroutine reducepickr(a,lonfd,lonf,latg)
!
! pick the reduce-grid value from the nearest regular grid
! then fill the tailing points to be the same as the last
! reduce-grid point.
!
      real a(lonf,latg)
      dimension lonfd(latg)
      real,allocatable:: tmp(:)
      allocate(tmp(lonf))
      dg=360./float(lonf)
      do j=1,latg
        dr=360./float(lonfd(j))
        do i=1,lonfd(j)
          ii=nint((i-1.)*dr/dg + 1.0)
          tmp(i)=a(ii,j)
        enddo
        ii=lonfd(j)
        do i=1,ii
          a(i,j)=tmp(i)
        enddo
        do i=ii+1,lonf
          a(i,j)=tmp(ii)
        enddo
      enddo
      deallocate(tmp)
      return
      end
!
      subroutine reducepickr_sp(a,lonfd,lonf,latg)
!
! pick the reduce-grid value from the nearest regular grid
! then fill the tailing points to be the same as the last
! reduce-grid point.
!
      use const, only: RTYPE
      real(kind=RTYPE) a(lonf,latg)
      dimension lonfd(latg)
      real,allocatable:: tmp(:)
      allocate(tmp(lonf))
      dg=360./float(lonf)
      do j=1,latg
        dr=360./float(lonfd(j))
        do i=1,lonfd(j)
          ii=nint((i-1.)*dr/dg + 1.0)
          tmp(i)=a(ii,j)
        enddo
        ii=lonfd(j)
        do i=1,ii
          a(i,j)=tmp(i)
        enddo
        do i=ii+1,lonf
          a(i,j)=tmp(ii)
        enddo
      enddo
      deallocate(tmp)
      return
      end
