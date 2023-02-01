      subroutine rqkir ( nr,nx,lev,ign,ozon,coefkr,pa,ta,fkg )
!
!#####################################################################
!
!  compute the absorption coefficient at the temperature t for the
!  19 pressures by
!        ln k = a + b * ( t - tcon ) + c * ( t - tcon ) ** 2
!  and the absorption coefficient at conditions other than those 19
!  pressures is interpolated linearly with pressure (Fu, 1991).
!  ( note : tcon = 250 for ozone, otherwise, tcon = 245 )
!
!  ++ input variables :
!     ign    : total given g-value intervals
!     ozon   : logical constant for checking if it is ozone absorption
!     coefkr : coefficients to calculate the absorption coefficients
!              i.e. the above's a, b, and c 
!
!  ++ output variables :
!     fkg    : absorption coefficients
!
!######################################################################
!
      implicit  none

      integer   nr,nx,lev,ign

      logical   ozon

      real      coefkr(3,19,ign),pa(nx*lev),ta(nx*lev),fkg(nx*lev,ign)
!
!  local arrays
!
      real      taa(nx*lev),paa(nx*lev)
      real      stanp(19)
!
      data stanp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51,     &
                   3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,   &
                   100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /

      integer   ltop,i,j,ig,jd,jj
      real      tcon,t1,t2,x1,x2
!
      ltop = 9
!
      if ( ozon ) then 
        tcon = 250.
      else
	tcon = 245.
      end if
!
!  the below polynomial function valid only for 200 to 300 K
!
      do 50 i = 1, nr*lev
      taa(i) = min( 300., max( 200., ta(i) ) )
      paa(i) = min( 999.999, max( stanp(ltop)+0.001, pa(i) ) )
   50 continue
!
      jd = 1
      do 500 ig = 1, ign
       do 400 j = ltop, 19, jd
       jj = j + jd
        do 300 i = 1, nr*lev
        t1 = (taa(i)-tcon)
        t2 = t1 * t1
        if( paa(i).gt.stanp(j) .and. paa(i).le.stanp(jj) ) then
 	x1 = coefkr(1,j,ig)+coefkr(2,j,ig)*t1+coefkr(3,j,ig)*t2
 	x2 = coefkr(1,jj,ig)+coefkr(2,jj,ig)*t1+coefkr(3,jj,ig)*t2
	fkg(i,ig) = exp ( x1 + ( x2-x1 ) / ( stanp(jj)-stanp(j) )  &
                    * ( paa(i)-stanp(j) ) )
        end if
  300   continue
  400  continue
  500 continue
!
!ccc----
!     do 200 ig = 1, ign
!     do 100  i = 1, nr*lev
!     ll = int( ( log10(pa(i)/1000.)/0.2 )+19.0001 ) 
!     ll = min(ll,18)
!org  t1 = (ta(i)-tcon)
!     t1 = (taa(i)-tcon)
!     t2 = t1 * t1
!     if ( ll .lt. 1 ) then
!       x1 = exp ( coefkr(1,1,ig) + coefkr(2,1,ig) * t1
!    1       + coefkr(3,1,ig) * t2 )
!       fkg(i,ig) = x1 * pa(i) / stanp(1)
!     else
!       x1 = exp ( coefkr(1,ll,ig) + coefkr(2,ll,ig)
!    1       * t1 + coefkr(3,ll,ig) * t2 )
!       x2 = exp ( coefkr(1,ll+1,ig) + coefkr(2,ll+1,ig)
!    1       * t1 + coefkr(3,ll+1,ig) * t2 )
!       fkg(i,ig) = x1 + ( x2 - x1 ) / ( stanp(ll+1) - stanp(ll) )
!    1              * ( pa(i) - stanp(ll) )
!     end if 
! 100 continue
! 200 continue
!cc----
      return
      end
