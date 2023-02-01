      subroutine sortml (jtrun,mlmax,msort,lsort,mlsort)
!
!  sortml builds pointer arrays for functional dependency between
!  1-d spectral index and zonal and total wavenumber indices.  this
!  subroutine reflects the coefficient storage strategy used in the
!  model
!
! ***input***
!
!  jtrun:  zonal and total wavenumber limit
!  mlmax:  total number of 1-d spectral index for triangular trunc
!
!  ***output***
!
!  msort:  zonal wavenumber as function of 1-d spectral index
!  lsort:  total wavenumber as function of 1-d spectral index
!  mlsort: total wavenumber index as function of zonal and total
!          wavenumber
!
! *******************************************************************
!
      implicit  none
      integer   jtrun,mlmax
      integer   msort(mlmax),lsort(mlmax),mlsort(jtrun,jtrun)
      integer   mlx,ml,k,m,mlp
!
      mlx= (jtrun/2)*((jtrun+1)/2)
      ml= 0
      do 1 k=1,jtrun-1,2
      do 1 m=1,jtrun-k
      ml= ml+1
      mlp= ml+mlx
      mlsort(m,m+k)= mlp
      mlsort(m,m+k-1)= ml
      msort(ml)= m
      lsort(ml)= m+k-1
      msort(mlp)= m
      lsort(mlp)= m+k
    1 continue
!
      ml= mlp
      do 2 m=2,jtrun,2
      ml= ml+1
      mlsort(m,jtrun)= ml
      msort(ml)= m
      lsort(ml)= jtrun
    2 continue
      return
      end
