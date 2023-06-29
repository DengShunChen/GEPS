      subroutine roptir( nr,nx,lev,ib,dz,cld,clwc,ciwc,cre,cde,tauwc  &
                       , tauic,wwc,wic,wwc1,wwc2,wic1,wic2 )
!
!#####################################################################
!
!  (1)compute the optical properties of clouds for a given ir band.
!
!  (2)the optical properties for water clouds use
!  the mean single scattering properties of the eight drop
!  size distributions in each spectral band, the single scattering
!  properties of a water cloud with the given liquid water content
!  and effective radius are obtained by interpolating (Eqs. 4.25 -
!  4.27 of Fu, 1991).
!
!  (3)the optical properties for ice clouds follow the parameterization
!  of Fu et al. (1998).
!
!  ++ input variables :
!   ib  : ir band number ( 1 - 12 )
!   dz(nx*lev) : layer thickness in unit of km
!   clwc(nx*lev) : cloud liquid water content (g/m3)
!   cre(nx*lev)  : effective size of water cloud (um)
!   ciwc(nx*lev) : cloud ice water content (g/m3)
!   cde(nx*lev)  : effective size of ice cloud (um)
!
!  ++ output variables :
!   tauwc, tauic : optical depth for water and ice clouds
!   wwc, wic : single scattering albedo
!   wwc1~2, wic1~2 : first two expansion coeffs. of the phase function
!
!#####################################################################
!
      implicit  none

      integer   nr,nx,lev,ib

      real      dz(nx*lev),clwc(nx*lev),ciwc(nx*lev),cre(nx*lev)    &
               ,cde(nx*lev),tauwc(nx*lev),tauic(nx*lev),wwc(nx*lev) &
               ,wic(nx*lev),wwc1(nx*lev),wwc2(nx*lev),wic1(nx*lev)  &
               ,wic2(nx*lev),cld(nx*lev)
!
!  local arrays
!
      real      re(8),fl(8),bz(8,12),wz(8,12),gz(8,12)              &
               ,ap(3,12),bp(4,12),cp(4,12)

      integer   ibs(nx*lev),jre(nx*lev)
!
! bz, wz and gz are the extinction coefficient(1/km), single scattering
! albedo and asymmetry factor for the water clouds (St II, Sc I, St I,
! As, Ns, Sc II, Cu, and Cb) in different bands.   re is the effective
! radius and fl is the liquid water content (LWC).  See Tables 4.2-4.4
! of Fu (1991).
!
      data re / 4.18, 5.36, 5.89, 6.16, 9.27, 9.84, 12.10, 31.23 /
      data fl / 0.05, 0.14, 0.22, 0.28, 0.50, 0.47, 1.00, 2.50 /
      data bz /                                                         &
          22.44, 57.35, 84.41,103.50,103.49, 84.17,152.77,132.07        &
      ,   18.32, 52.69, 76.67,100.31,105.46, 92.86,157.82,133.03        &
      ,   17.27, 50.44, 74.18, 96.76,105.32, 95.25,158.07,134.48        &
      ,   13.73, 44.90, 67.70, 90.85,109.16,105.48,163.11,136.21        &
      ,   10.30, 36.28, 57.23, 76.43,106.45,104.90,161.73,136.62        &
      ,    7.16, 26.40, 43.51, 57.24, 92.55, 90.55,149.10,135.13        &
      ,    6.39, 21.00, 33.81, 43.36, 66.90, 63.58,113.83,125.65        &
      ,   10.33, 30.87, 47.63, 60.33, 79.54, 73.92,127.46,128.21        &
      ,   11.86, 35.64, 54.81, 69.85, 90.39, 84.16,142.49,135.25        &
      ,   10.27, 33.08, 51.81, 67.26, 93.24, 88.60,148.71,140.42        &
      ,    6.72, 24.09, 39.42, 51.68, 83.34, 80.72,140.14,143.57        &
      ,    3.92, 14.76, 25.32, 32.63, 60.85, 58.81,112.30,145.62 /
      data wz /                                                         &
          .9193, .8962, .8859, .8810, .8127, .7816, .7754, .6373        &
      ,   .8747, .8611, .8478, .8516, .7871, .7729, .7531, .6186        &
      ,   .7647, .7524, .7365, .7434, .6712, .6593, .6394, .5499        &
      ,   .8075, .8087, .7959, .8054, .7505, .7555, .7094, .5719        &
      ,   .7533, .7720, .7672, .7770, .7512, .7609, .7125, .5682        &
      ,   .6327, .6763, .6846, .6935, .7079, .7177, .6824, .5528        &
      ,   .2888, .3484, .3716, .3803, .4545, .4657, .4754, .4938        &
      ,   .2618, .3062, .3213, .3330, .3929, .4068, .4174, .4845        &
      ,   .2958, .3399, .3524, .3655, .4162, .4303, .4352, .4913        &
      ,   .3012, .3547, .3693, .3819, .4336, .4473, .4474, .4869        &
      ,   .2437, .3187, .3446, .3527, .4279, .4389, .4459, .4772        &
      ,   .1090, .1872, .2268, .2249, .3313, .3359, .3748, .4570 /
      data gz /                                                         &
          .818, .805, .824, .830, .815, .801, .820, .845                &
      ,   .810, .802, .826, .840, .829, .853, .840, .868                &
      ,   .774, .766, .799, .818, .815, .869, .834, .869                &
      ,   .734, .728, .767, .797, .796, .871, .818, .854                &
      ,   .693, .688, .736, .772, .780, .880, .808, .846                &
      ,   .643, .646, .698, .741, .759, .882, .793, .839                &
      ,   .564, .582, .637, .690, .719, .871, .764, .819                &
      ,   .466, .494, .546, .609, .651, .823, .701, .766                &
      ,   .375, .410, .455, .525, .583, .773, .637, .710                &
      ,   .262, .301, .334, .406, .485, .695, .545, .631                &
      ,   .144, .181, .200, .256, .352, .562, .413, .517                &
      ,   .060, .077, .088, .112, .181, .310, .222, .327 /
!
!  ap and bp are coefficients to calculate the extiction and
!  absorption coefficients (1/m) respectively for a randomly
!  oriented hexagonal ice crystal in the ir band.
!  cp is the coefficients to calculate the asymmetry factor.
!  They are based on the table 1 of Fu et al. (1998).
!  The units of mean effective size and ice water content
!  in these calculation are um and g/m3 respectively.
!
      data ap /                                                         &
       -2.3088e-03, 2.8140, 1.0722e+00,-2.4652e-03, 2.8331,-4.2275e-01  &
      ,-3.0345e-03, 2.9000,-1.8499e+00,-4.9366e-03, 3.0877,-3.8842e+00  &
      ,-8.1786e-03, 3.4012,-8.8128e+00,-8.3726e-03, 3.4550,-1.5166e+01  &
      ,-1.6916e-03, 2.7657,-8.3310e+00,-4.1594e-03, 3.0473,-5.0615e+00  &
      ,-9.5241e-03, 3.5877,-1.0688e+01,-1.3348e-02, 4.0438,-2.1710e+01  &
      , 3.3257e-03, 2.6013,-1.9096e+01, 4.9196e-03, 2.3277,-1.3908e+01 /
      data bp/                                                          &
          4.3464e-01,  1.7214e-02, -1.6232e-04,  5.5615e-07             &
      ,   7.4289e-01,  1.2796e-02, -1.3918e-04,  5.1801e-07             &
      ,   8.8624e-01,  1.2265e-02, -1.5230e-04,  6.0008e-07             &
      ,   7.1522e-01,  1.6217e-02, -1.8685e-04,  7.0787e-07             &
      ,   5.8743e-01,  1.8766e-02, -2.0458e-04,  7.5100e-07             &
      ,   5.4095e-01,  1.9496e-02, -2.0509e-04,  7.3646e-07             &
      ,   1.1955e+00,  3.3506e-03, -5.2669e-05,  2.2333e-07             &
      ,   1.4664e+00, -2.1292e-03, -1.3616e-05,  1.1936e-07             &
      ,   9.5514e-01,  1.3097e-02, -1.7936e-04,  7.3133e-07             &
      ,   3.0037e-01,  2.0515e-02, -1.9316e-04,  6.5830e-07             &
      ,   2.0055e-01,  2.1326e-02, -1.7510e-04,  5.3558e-07             &
      ,   8.8697e-01,  2.1184e-02, -2.7814e-04,  1.0945e-06 /
      data cp/                                                          &
          7.9627e-01,  3.0034e-03, -2.0823e-05,  5.3665e-08             &
      ,   8.4729e-01,  2.5599e-03, -2.1826e-05,  6.8799e-08             &
      ,   8.7416e-01,  2.4554e-03, -2.4569e-05,  8.6412e-08             &
      ,   8.5228e-01,  2.5236e-03, -2.1491e-05,  6.6850e-08             &
      ,   8.6096e-01,  2.2004e-03, -1.7481e-05,  5.1766e-08             &
      ,   8.9062e-01,  1.9032e-03, -1.7335e-05,  5.8550e-08             &
      ,   8.6633e-01,  2.7979e-03, -3.1870e-05,  1.2172e-07             &
      ,   7.9840e-01,  3.9771e-03, -4.4719e-05,  1.6949e-07             &
      ,   7.3634e-01,  4.7982e-03, -4.5132e-05,  1.5257e-07             &
      ,   7.2604e-01,  2.6643e-03, -1.2511e-05,  2.2433e-08             &
      ,   6.8914e-01,  6.1922e-03, -6.4595e-05,  2.4369e-07             &
      ,   4.9492e-01,  1.1861e-02, -1.2676e-04,  4.6035e-07 /
!

      integer i,j,ii,nt
      real    gg,x1,x2,fw1,fw2,fw3,bext,babs

      do 20 i = 1, nr*lev
      tauwc(i) = 0.0
      tauic(i) = 0.0
      wwc(i)   = 0.0
      wic(i)   = 0.0
      wwc1(i)  = 0.0
      wwc2(i)  = 0.0
      wic1(i)  = 0.0
      wic2(i)  = 0.0
   20 continue
!
      nt = 0
      do 50 i = 1, nr*lev
       if( cld(i) .gt. 0.01 .and. clwc(i) .gt. 1.0e-5 ) then
         cre(i) = max( re(1)+0.0001, min( re(8)-0.0001, cre(i) ) )
         nt = nt + 1
         ibs(nt) = i
       end if
   50 continue
!
      if( nt .eq. 0 ) go to 190
!
      do 60 j  = 1, 7
      do 60 ii = 1, nt
       i = ibs(ii)
       if(cre(i) .ge. re(j) .and. cre(i) .le. re(j+1)) then
         jre(ii) = j
       end if
   60 continue
!
!-----------------------------------
!  water cloud optical properties
!-----------------------------------

      do 100 ii = 1, nt
        i = ibs(ii)
        j = jre(ii)
!
! A cloud with the effective radius smaller than 4.18 um is assumed
! to have an effective radius of 4.18 um with respect to the single
! scattering properties.
!
! A cloud with the effective radius larger than 31.23 um is assumed
! to have an effective radius of 31.18 um with respect to the single
! scattering properties.
!
       tauwc(i) = dz(i) * clwc(i) * ( bz(j,ib) / fl(j) +                &
                ( bz(j+1,ib) / fl(j+1) - bz(j,ib) / fl(j) ) /           &
                ( 1.0 / re(j+1) - 1.0 / re(j) ) * ( 1.0 / cre(i)        &
                - 1.0 / re(j) ) )
         wwc(i) = wz(j,ib) + ( wz(j+1,ib) - wz(j,ib) ) /                &
                ( re(j+1) - re(j) ) * ( cre(i) - re(j) )
             gg = gz(j,ib) + ( gz(j+1,ib) - gz(j,ib) ) /                &
                ( re(j+1) - re(j) ) * ( cre(i) - re(j) )
             x1 = gg
             x2 = x1 * gg
        wwc1(i) = 3.0 * x1
        wwc2(i) = 5.0 * x2
  100 continue

  190 continue

!---------------------------------
!  ice cloud optical properties
!---------------------------------

      do 200  i = 1, nr*lev
        if ( cld(i) .gt. 0.01 .and. ciwc(i) .gt. 1.0e-5 ) then
          cde(i)= max( 10., cde(i) )
          fw1 = cde(i)
          fw2 = fw1 * cde(i)
          fw3 = fw2 * cde(i)
              bext = ap(1,ib) + ap(2,ib)/fw1 + ap(3,ib)/fw2
          tauic(i) = dz(i) * 1000.0 * ciwc(i) * bext
              babs = bp(1,ib)/fw1+bp(2,ib)+bp(3,ib)*fw1+bp(4,ib)*fw2
            wic(i) = 1.0 - babs/bext
                gg = cp(1,ib)+cp(2,ib)*fw1+cp(3,ib)*fw2+cp(4,ib)*fw3
                x1 = gg
                x2 = x1 * gg
           wic1(i) = 3.0 * x1
           wic2(i) = 5.0 * x2
        end if
  200 continue
!
      return
      end
