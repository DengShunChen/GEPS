      subroutine srms (im,lm,jm,jm_max,rmsx,x)
!
!  routine to compute global mean rms of a variable on each model
!  level.
!
! **** input ****
!
!  im: no. of points in e-w direction of input grid x
!  lm: no. of vertical levels of x
!  jm: no. of points in n-s direction of input grid x
!  x: 3-d input grid
!
! **** output ****
!
!  rmsx: array of rms values
!
      use mpe
      use index

      implicit  none
      integer   im,lm,jm,jm_max
      real      rmsx(lm),x(im,lm,jm_max)
      integer   i,j,k,jj
      real      tem

      tem= 1.0/(im*jm)
!
      do 1 k=1,lm
      rmsx(k)= 0.0
    1 continue
      do 2 jj =1, jlistnum
      j=jlist1(jj)
      do 2 k=1,lm
      do 2 i=1,im
      rmsx(k)= rmsx(k)+x(i,k,jj)*x(i,k,jj)
    2 continue
!
      call mpe_global_sum(rmsx,lm,mpe_double)
!
      do 3 k=1,lm
      rmsx(k)= sqrt(rmsx(k)*tem)
    3 continue
!
      return
      end
