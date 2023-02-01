      module stochastic
 
      implicit none

      public

      integer, parameter :: num_domn=20,n_area=15

      integer n_domn,nx_st,my_st,ix_st(2,num_domn),jy_st(2,num_domn)
      common/stochasI/n_domn,nx_st,my_st,ix_st,jy_st

      real    fstoc,taurepst,size_st,position(2,num_domn), &
              area_lon(n_area),area_lat(n_area)
      common/stochasR/fstoc,taurepst,size_st,position,area_lon,area_lat

      logical dostoc
      common/stochasL/dostoc
 
      end module stochastic
