 subroutine voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,ff,flp,pkout, &
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
      use openacc
      implicit  none
      integer   nx,my,my_max,lev,lpout 



!byl      real      ff(nxp,lev,my_max),t(nx,my,lpout),pkout(lpout)    &
!byl      , pk(nxp,lev,my_max),tensy(lev+1),pklp(nx,my),flp(nx,my)

      real     ff(nxp,lev,my_max)            &
              ,pk(nxp,lev,my_max)            &
              ,tensy(lev+1),   pkout(lpout)  &
              ,pklp(nxp,my_max)              &
              ,flp(nxp,my_max)
!
      real(kind=RTYPE) t(nxp,my_max,lpout)

!local variable
      integer levp1,jym2,nxjym2,mn,ll,k,i,jj,j,nxj,ii
      integer l,m,n,mn1,ix,nm,inx
      real:: eps = 1.0e-5
      real,allocatable:: f(:,:,:), pkk(:,:,:), pout(:,:,:),wrk(:,:,:)
      real,allocatable:: tp1(:,:,:,:)
      real(kind=RTYPE),allocatable:: pjy(:,:,:,:),fxx(:,:,:),fyy(:,:,:)
      integer,allocatable::  ipt(:,:,:)

      integer ::stream, async_id = 1
      stream = acc_get_cuda_stream(async_id)


!
!      real      wk1(nx,my)
!
!          compute ipt and pjy
!


      levp1= lev+1
      jym2 = levp1-1
      allocate(  wrk(nxp,lev+1  ,my_max)  )
      allocate(  pkk(nxp,lev+1  ,my_max) ,f(nxp,lev+1,my_max) )
      allocate( pout(nxp,lpout  ,my_max) )
      allocate(  ipt(nxp,lpout  ,my_max) )
      allocate(  tp1(nxp,lpout,4,my_max) )
      allocate(  pjy(nxp,lpout,4,my_max) )
      allocate(  fxx(nxp,lev+1  ,my_max), fyy(nxp,lev+1,my_max) )

!! -----input
!$acc wait(async_id)
!!!!!$acc enter data copyin(tensy)  async(async_id)
!! -----output
!!!!!$acc enter data create( t ) async(async_id)

!! ----- local
!$acc enter data create(pout,pkk,f,fxx,fyy,tp1,pjy,ipt,wrk) async(async_id)

!! -----
mn1 = nxp*lpout

mn  = nxj*lpout
!fyy = 0.0
!tp1 = 0.0
!pjy = 0.0


!$acc parallel loop collapse(3) private(nxj,j) async(async_id)
do jj =1, jlistnum
  do k=1,lpout
   do i=1,nxp
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
     pout(i,k,jj)= pkout(k)
    endif
    ipt(i,k,jj) = 0
   enddo
  enddo
enddo

!$acc parallel loop collapse(3) private(j,nxj)  async(async_id)
do jj =1, jlistnum
  do k=1,lev
    do i=1,nxp
     j=jlist1(jj)
     nxj=nxdef_2d(j)
     if(i<=nxj)then
      pkk(i,k,jj) = pk(i,k,jj)
      f  (i,k,jj) = ff(i,k,jj)
     endif
    enddo
  enddo
enddo

!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
do jj =1, jlistnum
  do i=1,nxp
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
      pkk(i,levp1,jj) = pklp(i,jj)
      f  (i,levp1,jj) = flp (i,jj)
    endif
  enddo
enddo
!
!!$acc wait(async_id)
!!$acc update self( pout,pkk,f ) async(async_id)
!!$acc wait(async_id)
!  call setupv(pkk(1:nxj,1:levp1),pout(1:nxj,1:lpout),mn,nxj,levp1, &
!              pjy(1:nxj,1:lpout,1:4),ipt(1:nxj,1:lpout),              &
!              tp1(1:nxj,1:lpout,1:4))
!!$acc wait(async_id)
!!$acc update device(pjy,ipt,tp1) async(async_id)
!!$acc wait(async_id)


  !n   = mn1/nxp

!$acc parallel loop collapse(3) private(ix,ii,j,nxj) async(async_id)
do jj =1, jlistnum
  do k=1,lpout
   do  i = 1,nxp !mn1
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
     ii=i +(k-1)*nxp
     ix = ii - ((ii-1)/nxp)*nxp
     if ( pout(i,k,jj).le.pkk(ix,1    ,jj) )  pout(i,k,jj) = pkk(ix,1    ,jj) + eps
     if ( pout(i,k,jj).ge.pkk(ix,levp1,jj) )  pout(i,k,jj) = pkk(ix,levp1,jj) - eps
    endif
   enddo
  enddo
enddo
!
!$acc parallel loop collapse(2) private(ii,k,j,nxj) async(async_id)
do jj =1, jlistnum
  do  i=1,nxp
   j=jlist1(jj)
   nxj=nxdef_2d(j)
   if(i<=nxj)then
    k=1
    do  m=1,lpout !n
      ii= nxp*(m-1)
    6 k= k+1
      if(pout(i+ii,1,jj).gt.pkk(i,k,jj)) go to 6
      ipt(i+ii,1,jj)= i+nxp*(k-1)
      tp1(i+ii,1,1,jj)= pout(i+ii,1,jj) -pkk(i,k-1,jj)
      tp1(i+ii,1,2,jj)= pkk (i,k   ,jj) -pkk(i,k-1,jj)
      k= k-1
    enddo
   endif
  enddo
enddo
!
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
do jj =1, jlistnum
  do k=1,lpout
   do i=1,nxp!mn1
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
     pjy(i,k,3,jj)= tp1(i,k,1,jj)/tp1(i,k,2,jj)
     pjy(i,k,4,jj)= 1.0-pjy(i,k,3,jj)
     pjy(i,k,1,jj)= pjy(i,k,3,jj)*pjy(i,k,3,jj)-1.0
     pjy(i,k,2,jj)= pjy(i,k,4,jj)*pjy(i,k,4,jj)-1.0
     pjy(i,k,1,jj)= pjy(i,k,1,jj)*tp1(i,k,1,jj)*tp1(i,k,2,jj)
     tp1(i,k,1,jj)= tp1(i,k,2,jj)-tp1(i,k,1,jj)
     pjy(i,k,2,jj)= pjy(i,k,2,jj)*tp1(i,k,1,jj)*tp1(i,k,2,jj)
    endif
   enddo
  enddo
enddo

!
!          compute fyy
!
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
do jj =1, jlistnum
  do k=2,lev+1
   do i=1,nxp
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
     fxx(i,k,jj)=pkk(i,k,jj)-pkk(i,k-1,jj)
    endif
   enddo
  enddo
enddo
!
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
do jj =1, jlistnum
  do k=2,jym2
   do i=1,nxp
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
    fyy(i,k,jj)= ( fxx(i,k+1,jj) * ( f(i,k-1,jj)-f(i,k,jj) )   + &
                   fxx(i,k  ,jj) * ( f(i,k+1,jj)-f(i,k,jj) ) ) / &
                 ( fxx(i,k+1,jj) * fxx(i,k+1,jj) * fxx(i,k,jj) )
    endif
   enddo
  enddo
enddo

!!ocl serial

!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
do jj =1, jlistnum
  do i=1,nxp
   j=jlist1(jj)
   nxj=nxdef_2d(j)
   if(i<=nxj)then
    do k=2,jym2
     fxx(i,k,jj)= fxx(i,k,jj)/fxx(i,k+1,jj)
    enddo
   endif
  enddo
enddo
!
!!$acc wait(async_id)
!!$acc update self( fxx,fyy ) async(async_id)
!!$acc wait(async_id)
!      call trdivv(nxj,jym2,fxx(1:nxj,2:lev),fyy(1:nxj,2:lev))
!!$acc wait(async_id)
!!$acc update device(fyy) async(async_id)
!!$acc wait(async_id)


!$acc parallel loop collapse(2) private(j,nxj,nm) async(async_id)
do jj =1, jlistnum
  do i=1,nxp
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
      nm = jym2 - 1

      wrk(i,2,jj) = 0.5 / (  1.0 + fxx(i,2,jj) )
      fyy(i,2,jj) = fyy(i,2,jj) * wrk(i,2,jj)
! gaussian elimination
      do m=3, nm 
        wrk(i,m  ,jj)= 1.0/(2.0+fxx(i,m  ,jj)*(2.0-wrk(i,m-1,jj)))
        fyy(i,m  ,jj) = (fyy(i,m  ,jj)-fxx(i,m,jj)*fyy(i,m-1,jj))*wrk(i,m,jj)
      enddo
      101 continue
! 
      fyy(i,jym2,jj)= (fyy(i,jym2,jj)-fxx(i,jym2,jj)*fyy(i,nm,jj)) /  & 
                   ( 2.0 + fxx(i,jym2,jj) * ( 2.0 - wrk(i,nm,jj) ) )
      202 continue
!
! backwards substitution
!
      do  k=nm,2,-1
        fyy(i,k,jj)= fyy(i,k,jj)-wrk(i,k,jj)*fyy(i,k+1,jj)
      enddo
      104 continue
    endif

  enddo
enddo
!=====end subroutine trdivv
!
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
do jj =1, jlistnum
  do i=1,nxp
   j=jlist1(jj)
   nxj=nxdef_2d(j)
   if(i<=nxj)then
    fyy(i,1    ,jj)= 0.0
    fyy(i,levp1,jj)= 0.0
   endif
  enddo
enddo
!
!  apply tension
!
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
do jj =1, jlistnum
  do k=1,levp1
   do i=1,nxp
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
     fyy(i,k,jj)= fyy(i,k,jj)*tensy(k)
    endif
   enddo
  enddo
enddo

!!$acc wait(async_id)
!!$acc update self( fyy ) async(async_id)
!!$acc wait(async_id)
!      call gathv(mn,nxj,levp1,ipt(1:nxj,1:lpout),fyy(1:nxj,1:levp1),     &
!                 f(1:nxj,1:levp1),tp1(1:nxj,1:lpout,1:4))
!!$acc wait(async_id)
!!$acc update device(tp1) async(async_id)
!!$acc wait(async_id)


!$acc parallel loop collapse(3) private(inx,j,nxj)  async(async_id) 
do jj =1, jlistnum
  do k=1,lpout
   do i=1,nxp!mn1
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
     tp1(i,k,1,jj) = fyy(ipt(i,k,jj),1,jj)
     tp1(i,k,3,jj) = f  (ipt(i,k,jj),1,jj)
     inx= ipt(i,k,jj) - nxp
     tp1(i,k,2,jj) = fyy(inx        ,1,jj)
     tp1(i,k,4,jj) = f  (inx        ,1,jj)
    endif
   enddo
  enddo
enddo

!
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
do jj =1, jlistnum
  do k=1,lpout
   do i=1,nxp
    j=jlist1(jj)
    nxj=nxdef_2d(j)
     if(i<=nxj)then
     t(i,jj,k)=tp1(i,k,1,jj)*pjy(i,k,1,jj)+tp1(i,k,2,jj)*pjy(i,k,2,jj)+  &
               tp1(i,k,3,jj)*pjy(i,k,3,jj)+tp1(i,k,4,jj)*pjy(i,k,4,jj)
    endif
   enddo
  enddo
enddo !  jj =1, jlistnum
  200 continue

!$acc wait(async_id)
!! -----local
!$acc exit data delete(pout,pkk,f,fxx,fyy,tp1,pjy,ipt,wrk) async(async_id)
!! -----input
!!!!!!$acc exit data delete( tensy)  async(async_id)
!! -----output
!!!!!!$acc update self( t ) async(async_id)
!$acc wait(async_id)
deallocate( wrk )
deallocate( pkk ,f  ,ipt ,pout)
deallocate( tp1 ,pjy,fxx ,fyy)


!!$omp end parallel do
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
