      subroutine stemp(incrup,nx,lev,cp,grav,sgeo,pk,pk2,phi,tt)
!
!  stemp computes approximatly hydrostatically consistent thickness
!  tempratures from analysis heights.  algorithm prevents buckling
!  of the temperature soundings
!

      implicit  none

      integer   nx,lev,lm2,l,l1,l2,i,ll,l3,jj,jim

      real      phi(nx,lev),tt(nx,lev),pk(nx,lev),pk2(nx,lev)  &
      , sgeo(nx)
!
      logical incrup
!
      real      phidif(nx),that(nx,lev)
      real      pexp(nx),xlaps(14),tempps(nx,2*lev-1)
!
      data xlaps/.002,.002,.00175,.0015,.001625,.00175,.0018   &
      ,   .0019,.002,.0021,.0022,.0023,.0024,.0025/

      real      dkapa,capa,pratio,rgas,ocp,grav,cp,wrk,pp,xl1,xl2,phimax
!
!  thickness temps on half levels
!
!
!  pressure at full and half levels
!
      lm2= 2*lev-1
      dkapa= 3.5
      capa= 1.0/dkapa
      pratio=1.4142
      rgas=287.04
      ocp= 1.0/cp
!
      do 5 l=1,lev-1
      l1=2*l-1
      l2=l1+1
      do 55 i=1,nx
      that(i,l)= ocp*(phi(i,l)-phi(i,l+1))/    &
       (pk(i,l+1)-pk(i,l))
   55 continue
    5 continue
!
      do 201 i=1,nx
      wrk=that(i,1)*pk2(i,1)-that(i,2)*pk2(i,2)+0.5
      jj=int(wrk)+1
      jj= max(1,jj)
      jj= min(14,jj)
      pp= xlaps(jj)
      pexp(i)=rgas*pp/grav
  201 continue
!
!  begin interation
!
      do 100 jim=1,10
!
!  analysis temps match co-located thickness temps at half levels
!
      do 6 l=1,lev-1
      l1= 2*l
      do 66 i=1,nx
      tempps(i,l1)=that(i,l)
   66 continue
    6 continue
!
!  interpolate to full levels for analysis temps
!
      do 7 l=1,lev-2
      do 7 i=1,nx
      xl1=(pk(i,l+1)-pk2(i,l))/(pk2(i,l+1)-pk2(i,l))
      tempps(i,2*l+1)=that(i,l)+(that(i,l+1)-that(i,l))*xl1
    7 continue
!
!
!  temperature boundary conditions at top and bottom
!
!      if(incrup) then
      do 56 i=1,nx
      xl1=(pk(i,lev)-pk(i,lev-1))/(pk2(i,lev-1)-pk(i,lev-1))
      xl2=(pk(i,2)-pk(i,1))/(pk(i,2)-pk2(i,1))
      tempps(i,lm2)= tempps(i,lm2-2)+xl1*(tempps(i,lm2-1) &
       -tempps(i,lm2-2))
      tempps(i,1)= tempps(i,3)+xl2*(tempps(i,2)-tempps(i,3))
   56 continue
!
!      else
!
!      do 78 i=1,nx
!      tempps(i,1)=that(i,1)*(pk2(i,2)/pk(i,1))*pratio**pexp(i)
!      tempps(i,lm2)= ocp*(phi(i,lev)-sgeo(i))/(pk2(i,lev)-pk(i,lev))
!   78 continue
!      endif
!
!  adjust thickness temperatures to be consistent with both analysis
!  temperatures and geopotentials
!
      phimax= 0.0
      do 8 l=lev-1,1,-1
      l1=2*l-1
      l2=l1+1
      l3=l2+1
      do 88 i=1,nx
      xl1=cp*(pk2(i,l)-pk(i,l))
      xl2=cp*(pk(i,l+1)-pk2(i,l))
      phidif(i)= (tempps(i,l1)+tempps(i,l2))*xl1     &
        +(tempps(i,l2)+tempps(i,l3))*xl2 
      phidif(i)= 0.5*phidif(i)-phi(i,l)+phi(i,l+1)
      that(i,l)= that(i,l)-phidif(i)/(xl1+xl2)
      phimax= max(phimax,abs(phidif(i)))
   88 continue
    8 continue
      if(phimax.lt.1.0) go to 110
  100 continue
  110 continue
!
      do 13 l=1,lev
      ll= 2*l-1
      do 13 i=1,nx
      tt(i,l)= tempps(i,ll)*pk(i,l)
   13 continue
!
      return
      end
