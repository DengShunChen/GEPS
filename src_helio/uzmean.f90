      subroutine uzmean (nx,my,my_max,lev,ut,uzm)
!
!  subroutine to compute mean zonal wind speed used in implicit
!  zonal advection of vorticity and moisture
!
!  *** input ***
!
!  nx: e-w (zonal) dimension
!  my: n-s dimension
!  lev: vertical dimension
!  ut: 3-dimensional u component array
!
!  *** output ***
!
!  uzm: zonal average wind for each latitude
!
!  ************************************
!
      use index
      use mpe

      implicit  none

      integer   nx,my,my_max,lev,jj,j,k,i,nxj,kk
      real      umax,umin

!     real      ut(nx,lev,my_max),uzm(my,lev)
      real      ut(nxp,lev,my_max),uzm(my,lev),ut1(nx,levp,my_max)

      call mpe2d_transpose1_nxp_lev(ut,ut1,nx,nxp,lev,levp,my,my_max,jlistnum,jlen,nsizex,row_comm)

      do 10 jj =1, jlistnum
      j=jlist1(jj)
      nxj=nxdef(j)

!     do 20 k=1,lev
      do 20 k=1,levp
      umax= 1.0e-20
      umin= 1.0e+20
      do 30 i=1,nxj
!     umax= max(umax,ut(i,k,jj))
!     umin= min(umin,ut(i,k,jj))
      umax= max(umax,ut1(i,k,jj))
      umin= min(umin,ut1(i,k,jj))
   30 continue

!     uzm(j,k)= 0.5*(umax+umin)
      uzm(jj,k)= 0.5*(umax+umin)

   20 continue
   10 continue

!     call mpe_unify(uzm,my,lev,1,mpe_double)

      call mpe2d_unify_uzmean(uzm,my,lev)
!
      return
      end
