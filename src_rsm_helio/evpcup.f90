      subroutine evpcup ( nn,nx,lev,ktmin,prevap,ql,qls,pmassl &
                        , pprod,dqo )
!
!  subroutine to partially evaporate cumlus precipitation as
!  it falls.  evaporation is proportional to layer mass and
!  inversely proportional to layer relative humidity.
!
!  parameters:
!
!  input:
!
!     nn    : number of convective points in latitude rung
!     nx    : horizontal dimension
!     lev   : number of levels
!     ktmin : deepest cloud top index
!     prevap: fraction of precip evaporated
!     ql    : specific humidity`
!     qls   : saturation specific humidity
!     pmassl:  layer mass
!     pprod : precip mass generated in each layer penetrated
!            by clouds
!     dqo   : water mass budget
!
!  output:
!
!     dqo   : adjusted water mass budget
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use paramt

      implicit  none

      integer   nn,nx,lev,ktmin
      real      prevap

      real      ql(nx,lev),qls(nx,lev),pmassl(nx,lev),  &
                pprod(nx,lev,lev),dqo(nx,lev,lev)
!
      integer   k,i,l,ktype
      real      rhmass(im,lm),xx(im,lm),dqcup(im,lm)
!
!  (layer mass)/relative humidity
!
      do 20 k = ktmin, lev
      do 20 i = 1, nn
      rhmass(i,k)= qls(i,k)/ql(i,k)
   20 continue
!
      do 30 i = 1, nn
      xx(i,lev-1)= rhmass(i,lev-1)*pmassl(i,lev-1)
   30 continue
!
      do 40 k = lev-2, ktmin, -1
      do 40 i = 1, nn
      xx(i,k) = xx(i,k+1)+rhmass(i,k)*pmassl(i,k)
   40 continue
!
      do 50 k = ktmin, lev-1
      do 50 i = 1, nn
      xx(i,k) = prevap/xx(i,k)
   50 continue
!
      do 60 ktype = ktmin, lev-1
!
      do 70 l = ktype, lev-1
      do 70 i = 1, nn
      dqcup(i,l)= pprod(i,l,ktype)*xx(i,l)
   70 continue
!
      do 80 l = ktype, lev-1
      do 80 k = l, lev
      do 80 i = 1, nn
      dqo(i,k,ktype)= dqo(i,k,ktype)+dqcup(i,l)*rhmass(i,k)
   80 continue
   60 continue
!
      return
      end
