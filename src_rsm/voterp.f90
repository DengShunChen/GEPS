      subroutine voterp(nx,my,my_max,lev,lpout,pk,pklp,ff,flp,pkout, &
                        t,tensy)
!
!          a bicubic spline interpolator to interpolate from a grid
!          with constant i (first dimension) grid spacing and variable
!          k (second dimension) grid spacing to a grid with
!          variable grid spacing. all grids are assumed to have point
!          (i,1) in the upper left corner with i increasing to the right
!          and k increasing downward.
!
! **** input ****
!
!  nx: no. of points in e-w direction of input arrays
!  my: no. of points in n-s direction of input arrays
!  lev: no. of levels in input arrays
!  lpout: no. of levels in output arrays
!  pk: independent interpolation variable for input grid
!  pklp: indepentdent inter var appended to bottom of pk array
!  ff: dependent variable for input grid
!  flp: bottom boundary condition of f, defined at pklp
!  pkout: independent interpolation variable for output grid
!  tensy: cubic spline tension factor.
!
! **** output ****
!
!  t: dependent variable on output grid
!
      use mpe
      use index
      use const, only: RTYPE

      implicit  none

      integer   nx,my,my_max,lev,lpout

      integer   levp1,jym2,nxjym2,mn,ll,k,i,jj,j,nxj,ii

!byl      real      ff(nxp,lev,my_max),t(nx,my,lpout),pkout(lpout)    &
!byl      , pk(nxp,lev,my_max),tensy(lev+1),pklp(nx,my),flp(nx,my)
      real      ff(nxp,lev,my_max),pkout(lpout)                        &
      , pk(nxp,lev,my_max),tensy(lev+1),pklp(nxp,my_max)               &
      , flp(nxp,my_max),tp1(nxp,lpout,4),pout(nxp,lpout)               &
      , pkk(nxp,lev+1)
!
      real(kind=RTYPE) fxx(nxp,lev+1),fyy(nxp,lev+1),pjy(nxp,lpout,4)  &
      , f(nxp,lev+1),t(nxp,my_max,lpout)

      integer   ipt(nxp,lpout)
!
!      real      wk1(nx,my)
!
!          compute ipt and pjy
!
      levp1= lev+1
      jym2=levp1-2
!
!$omp parallel do                                                      &
!$omp private(j,i,ii,pkk,f,pjy,ipt,tp1,fxx,fyy,k,nxj,mn,pout) &
!$omp schedule(dynamic)
!
      do 200 jj =1, jlistnum
      j=jlist1(jj)
!1d   nxj=nxdef(j)
      nxj=nxdef_2d(j)
!!      nxj=nxjp(j)

!     nxjym2=nxj*jym2
      mn= nxj*lpout
!     ll= nxjym2+nxj
!
      do k=1,lpout
       do i=1,nxj
        pout(i,k)= pkout(k)
       enddo
      enddo

      do k=1,lev
       do i=1,nxj
        pkk(i,k)= pk(i,k,jj)
        f(i,k)= ff(i,k,jj)
       enddo
      enddo
!
      do i=1,nxj
        pkk(i,levp1)= pklp(i,jj)
        f(i,levp1)= flp(i,jj)
      enddo
!
      call setupv(pkk(1:nxj,1:levp1),pout(1:nxj,1:lpout),mn,nxj,levp1,   &
                  pjy(1:nxj,1:lpout,1:4),ipt(1:nxj,1:lpout),           &
                  tp1(1:nxj,1:lpout,1:4))
!
!          compute fyy
!
      do k=1,lev
       do i=1,nxj
        fxx(i,k+1)=pkk(i,k+1)-pkk(i,k)
       enddo
      enddo
!
      do k=1,jym2
       do i=1,nxj
        fyy(i,k+1)=(fxx(i,k+2)*(f(i,k)-f(i,k+1))+fxx(i,k+1)*(f(i,k+2)  &
                    -f(i,k+1)))/(fxx(i,k+2)*fxx(i,k+2)*fxx(i,k+1))
       enddo
      enddo

!!ocl serial
      do k=1,jym2
       do i=1,nxj
        fxx(i,k+1)= fxx(i,k+1)/fxx(i,k+2)
       enddo
      enddo
!
      call trdivv(nxj,jym2,fxx(1:nxj,2:lev),fyy(1:nxj,2:lev))
!
      do i=1,nxj
       fyy(i,1)= 0.0
       fyy(i,levp1)=0.0
      enddo
!
!  apply tension
!
      do k=1,levp1
       do i=1,nxj
        fyy(i,k)= fyy(i,k)*tensy(k)
       enddo
      enddo
!
      call gathv(mn,nxj,levp1,ipt(1:nxj,1:lpout),fyy(1:nxj,1:levp1),     &
                 f(1:nxj,1:levp1),tp1(1:nxj,1:lpout,1:4))
!
      do k=1,lpout
       do i=1,nxj
        t(i,jj,k)=tp1(i,k,1)*pjy(i,k,1)+tp1(i,k,2)*pjy(i,k,2)+          &
                    tp1(i,k,3)*pjy(i,k,3)+tp1(i,k,4)*pjy(i,k,4)
       enddo
      enddo
!
  200 continue
!$omp end parallel do
!
! do 'reduceintp' after votertical interplot
!1d   if( lreduce.eq.1 ) then
!1d     do k=1,lpout
!1d       call reduceintp(t(1,1,k),nxdef,nx,my)
!1d     enddo
!1d   endif

!!        do k=1,lpout
!!          call mpe_unify(t(1,1,k),nx,my,2,mpe_double)
!2d>
!!          if( lreduce.eq.1 ) then
!!           do jj =1, jlistnum
!!             j=jlist1(jj)
!!             call reduceintp(t(1,j,k),nxdef(j),nx,1)
!!           enddo
!!           call mpe_unify(t(1,1,k),nx,my,5,mpe_double) 
!!          endif
!2d<
!!        enddo
!

      return
      end
