      subroutine getphi(nx,lev,cp,tt,pk,pk2,sgeo,phi)

      implicit  none
      integer   nx,lev
      real      cp
      real      tt(nx,lev),pk(nx,lev),pk2(nx,lev),sgeo(nx), &
                phi(nx,lev)
      real      that(nx,lev)
      integer   k,i
!
      do 160 k=1,lev-1
      do 160 i=1,nx
      that(i,k+1)= tt(i,k)-(tt(i,k)-tt(i,k+1))  &
       *(pk(i,k+1)-pk2(i,k))/(pk(i,k+1)-pk(i,k))
      phi(i,k)= that(i,k+1)*(pk(i,k+1)-pk(i,k))
  160 continue
      do 100 i=1,nx
      phi(i,lev)= sgeo(i)+cp*tt(i,lev)*(pk2(i,lev)-pk(i,lev))
  100 continue
      do 105 k=lev-1,1,-1
      do 105 i=1,nx
      phi(i,k)= phi(i,k+1)+cp*phi(i,k)
  105 continue
!
      return
      end
