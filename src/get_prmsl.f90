subroutine get_prmsl(nx,my,my_max,lev,ncld &
          ,sgeo,pt,tt,qt,pk,pk2,plt,phi,pdiff,mslp)
!  compute mean sea level pressure and update pdiff
!  The method is based on one used by ecmwf, reseach manual 2 (1988)
!
!  llts layer's temperature is used to derive an alternative
!  surface skin temperature
!
!202604  make a subroutine to get pdiff and mean sea level pressure
!
!=====
use index, only:nxp,jlistnum,jlist1,nxdef_2d
use const , only:RTYPE, ptop, rgas, grav, cp
implicit none
!input
integer nx,my,my_max,lev,ncld
real(kind=RTYPE) :: tt(nxp,lev,my_max),qt(nxp,lev*ncld,my_max),pk(nxp,lev,my_max) &
                  ,pk2(nxp,lev,my_max),plt(nxp,lev,my_max)    ,phi(nxp,lev,my_max)
real(kind=RTYPE) :: sgeo(nxp,my_max),pt(nxp,my_max)
!output
real :: mslp(nxp,my_max), pdiff(nxp,my_max)
!local
integer i,j,jj,nxj,k
integer:: llts
real(kind=RTYPE):: alaps=0.0065 ,rdg,ttt,ttb,ttp,ttt1,ttt2,apha,anlslp
real(kind=RTYPE)::hld1(nxp,my_max), hld2(nxp,my_max)
integer, parameter:: async_id = 1
!=====
rdg = rgas/grav
llts= lev-5
!!=====get phi
!!
!!  hydrostatic equation
!!
#ifdef USE_CUDA
!$acc wait(async_id)
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
do jj =1,jlistnum
  do i=1,nxp
    j=jlist1(jj)
    nxj=nxdef_2d(j)
    if(i<=nxj)then
    phi(i,lev,jj)= cp*tt(i,lev,jj)*(pk2(i,lev,jj)-pk(i,lev,jj)) &
                 + sgeo(i,jj)
    do k=lev-1,1,-1
        phi(i,k,jj)= phi(i,k+1,jj)+cp*(tt(i,k,jj)*(pk2(i,k,jj)-pk(i,k,jj)) &
                   + tt(i,k+1,jj)*(pk(i,k+1,jj)-pk2(i,k,jj)))
    enddo
   endif
  enddo
enddo
!$acc wait(async_id)
#else
do jj =1,jlistnum
  j=jlist1(jj)
  nxj=nxdef_2d(j)
  do i=1,nxj
    phi(i,lev,jj)= cp*tt(i,lev,jj)*(pk2(i,lev,jj)-pk(i,lev,jj)) &
                 + sgeo(i,jj)
  enddo
  do k=lev-1,1,-1
    do i=1,nxj
      phi(i,k,jj)= phi(i,k+1,jj)+cp*(tt(i,k,jj)*(pk2(i,k,jj)-pk(i,k,jj)) &
                 + tt(i,k+1,jj)*(pk(i,k+1,jj)-pk2(i,k,jj)))
    enddo
  enddo
enddo
#endif
!=====
#ifdef USE_CUDA
!$acc wait(async_id)
!$acc enter data create(hld1,hld2) async(async_id)

!$acc parallel loop collapse(2) &
!$acc& private(j,nxj,ttt,ttb,ttp,ttt1,ttt2,apha,anlslp) async(async_id)
do jj = 1, jlistnum
  do i=1,nxp
   j=jlist1(jj)
   nxj=nxdef_2d(j)
   if(i<=nxj)then 
    ttb  = tt(i,lev ,jj)*pk(i,lev ,jj)/(1.0+0.608*qt(i,lev ,jj))
    ttp  = tt(i,llts,jj)*pk(i,llts,jj)/(1.0+0.608*qt(i,llts,jj))
    ttt1 = ttb + alaps*rdg*ttb*   &
           ((pt(i,jj)+ptop)/plt(i,lev,jj)-1.0)
    ttt2 = ttp + alaps*(phi(i,llts,jj)-sgeo(i,jj))/grav
    hld1(i,jj) = 0.25*ttt1 + 0.75*ttt2
    hld2(i,jj) = hld1(i,jj) + alaps*sgeo(i,jj)/grav
    if( sgeo(i,jj) .lt. 0.1 ) then
      anlslp = pt(i,jj) + ptop
    else if( hld1(i,jj) .le. 290.5 .and. hld2(i,jj) .gt. 290.5 ) then
      apha = rgas*(290.5-hld1(i,jj))/sgeo(i,jj)
      ttt  = sgeo(i,jj)/(rgas*hld1(i,jj))
      anlslp = (pt(i,jj)+ptop)*exp( ttt*(1.0-0.5*apha*ttt+0.333333*  &
               apha*ttt*apha*ttt) )
    else if( hld1(i,jj) .gt. 290.5 .and. hld2(i,jj) .gt. 290.5 ) then
      hld1(i,jj) = (hld1(i,jj)+290.5)*0.5
      anlslp = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,jj)) )
    else if( hld1(i,jj) .lt. 255.0 .and. hld2(i,jj) .lt. 255.0 ) then
      hld1(i,jj) = (hld1(i,jj)+255.0)*0.5
      anlslp = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,jj)) )
    else
      apha = alaps * rdg
      ttt  = sgeo(i,jj)/(rgas*hld1(i,jj))
      anlslp = (pt(i,jj)+ptop)*exp( ttt*(1.0-0.5*apha*ttt+0.333333*  &
               apha*ttt*apha*ttt) )
    endif
    mslp (i,jj) =  anlslp
    pdiff(i,jj) =  anlslp - pt(i,jj)
   endif
  enddo
enddo
#else
do jj = 1, jlistnum
  j=jlist1(jj)
  nxj=nxdef_2d(j)
  do i=1,nxj
    ttb  = tt(i,lev ,jj)*pk(i,lev ,jj)/(1.0+0.608*qt(i,lev ,jj))
    ttp  = tt(i,llts,jj)*pk(i,llts,jj)/(1.0+0.608*qt(i,llts,jj))
    ttt1 = ttb + alaps*rdg*ttb*   &
           ((pt(i,jj)+ptop)/plt(i,lev,jj)-1.0)
    ttt2 = ttp + alaps*(phi(i,llts,jj)-sgeo(i,jj))/grav
    hld1(i,jj) = 0.25*ttt1 + 0.75*ttt2
    hld2(i,jj) = hld1(i,jj) + alaps*sgeo(i,jj)/grav
    if( sgeo(i,jj) .lt. 0.1 ) then
      anlslp = pt(i,jj) + ptop
    else if( hld1(i,jj) .le. 290.5 .and. hld2(i,jj) .gt. 290.5 ) then
      apha = rgas*(290.5-hld1(i,jj))/sgeo(i,jj)
      ttt  = sgeo(i,jj)/(rgas*hld1(i,jj))
      anlslp = (pt(i,jj)+ptop)*exp( ttt*(1.0-0.5*apha*ttt+0.333333*  &
               apha*ttt*apha*ttt) )
    else if( hld1(i,jj) .gt. 290.5 .and. hld2(i,jj) .gt. 290.5 ) then
      hld1(i,jj) = (hld1(i,jj)+290.5)*0.5
      anlslp = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,jj)) )
    else if( hld1(i,jj) .lt. 255.0 .and. hld2(i,jj) .lt. 255.0 ) then
      hld1(i,jj) = (hld1(i,jj)+255.0)*0.5
      anlslp = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,jj)) )
    else
      apha = alaps * rdg
      ttt  = sgeo(i,jj)/(rgas*hld1(i,jj))
      anlslp = (pt(i,jj)+ptop)*exp( ttt*(1.0-0.5*apha*ttt+0.333333*  &
               apha*ttt*apha*ttt) )
    endif
    mslp (i,jj) =  anlslp
    pdiff(i,jj) =  anlslp - pt(i,jj)
  enddo
enddo
#endif

!$acc wait(async_id)
!$acc exit data delete(hld1,hld2) async(async_id)
return
end subroutine
