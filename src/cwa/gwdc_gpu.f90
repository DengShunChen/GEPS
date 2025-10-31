      subroutine gwdc_gpu(ix, iy, km, u1, v1, t1, q1, &
                         pmid1, pint1, dpmid1, &
                         ktop, kbot, kcnv, cldf, qmax, &
                         grav, cp, rd, fv, deltim, dlength, &
                         utgwc, vtgwc, tauctx, taucty)
!***********************************************************************
!        original code for parameterization of convectively forced
!        gravity wave drag from yonsei university, korea
!        based on the theory given by chun and baik (jas, 1998)
!        modified for implementation into the gfs/cfs by
!        ake johansson  --- aug 2005
!***********************************************************************
!lzl +add === #
         use rank
         use index
         use param, only: my, my_max, ncld
         use const, only: RTYPE
!lzl -end === #

!lzl      use machine , only : kind_phys
         implicit none

!---------------------------- arguments --------------------------------
!
!  input variables
!
!  u        : midpoint zonal wind
!  v        : midpoint meridional wind
!  t        : midpoint temperatures
!  pmid     : midpoint pressures
!  pint     : interface pressures
!  dpmid    : midpoint delta p ( pi(k)-pi(k-1) )
!  lat      : latitude index
!  qmax     : deep convective heating
!  kcldtop  : vertical level index for cloud top    ( mid level )
!  kcldbot  : vertical level index for cloud bottom ( mid level )
!  kcnv     : (0,1) dependent on whether convection occur or not
!
!  output variables
!
!  utgwc    : zonal wind tendency
!  vtgwc    : meridional wind tendency
!
!-----------------------------------------------------------------------

         integer ix, iy, km, lat, ipr, jj, iii, kkk, iinpt, myim(my_max), j, im
         integer, dimension(ix, my_max) :: ktop, kbot, kcnv

!lzl      real(kind=kind_phys) grav,cp,rd,fv,fhour,fhourpr
!lzl      real(kind=kind_phys), dimension(ix)      :: qmax
!lzl     &,                                           tauctx, taucty
!lzl      real(kind=kind_phys), dimension(im)      :: cldf,dlength
!lzl      real(kind=kind_phys), dimension(ix,km)   :: u1,v1,t1,q1,
!lzl     &                                            pmid1,dpmid1
!lzl!    &,                                           cumchr1
!lzl      real(kind=kind_phys), dimension(iy,km)   :: utgwc,vtgwc
!lzl      real(kind=kind_phys), dimension(ix,km+1) :: pint1
!lzl
!lzl      logical lprnt
!
!lzl add ---- ##
         real grav, cp, rd, fv, fhour, fhourpr, deltim
         real, dimension(ix, my_max)      :: qmax, tauctx, taucty
         real, dimension(ix, my_max)      :: cldf, dlength
         real, dimension(ix, km, my_max)   :: u1, v1, t1, pmid1, dpmid1
!    &,                                           cumchr1
         real, dimension(iy, km, my_max)   :: utgwc, vtgwc
         real, dimension(ix, km + 1, my_max) :: pint1
         real(kind=RTYPE), dimension(ix, km*ncld, my_max) :: q1

         logical lprnt
!
!------------------------- local workspace -----------------------------
!
!  i, k     : loop index
!  kk       : loop index
!  cldf     : deep convective cloud fraction at the cloud top.
!  ugwdc    : zonal wind after gwdc paramterization
!  vgwdc    : meridional wind after gwdc parameterization
!  plnmid   : log(pmid) ( mid level )
!  plnint   : log(pint) ( interface level )
!  dpint    : delta pmid ( interface level )
!  tauct    : wave stress at the cloud top calculated using basic-wind
!             parallel to the wind vector at the cloud top ( mid level )
!  tauctx   : wave stress at the cloud top projected in the east
!  taucty   : wave stress at the cloud top projected in the north
!  qmax     : maximum deep convective heating rate ( k s-1 ) in a
!             horizontal grid point calculated from cumulus para-
!             meterization. ( mid level )
!  wtgwc    : wind tendency in direction to the wind vector at the cloud top level
!             due to convectively generated gravity waves ( mid level )
!  utgwcl   : zonal wind tendency due to convectively generated
!             gravity waves ( mid level )
!  vtgwcl   : meridional wind tendency due to convectively generated
!             gravity waves ( mid level )
!  taugwci  : profile of wave stress calculated using basic-wind
!             parallel to the wind vector at the cloud top
!  taugwcxi : profile of zonal component of gravity wave stress
!  taugwcyi : profile of meridional component of gravity wave stress
!
!  taugwci, taugwcxi, and taugwcyi are defined at the interface level
!
!  bruni    : brunt-vaisala frequency ( interface level )
!  brunm    : brunt-vaisala frequency ( mid level )
!  rhoi     : air density ( interface level )
!  rhom     : air density ( mid level )
!  ti       : temperature ( interface level )
!  basicum  : basic-wind profile. basic-wind is parallel to the wind
!             vector at the cloud top level. (mid level)
!  basicui  : basic-wind profile. basic-wind is parallel to the wind
!             vector at the cloud top level. ( interface level )
!  riloc    : local richardson number ( interface level )
!  rimin    : minimum richardson number including both the basic-state
!             and gravity wave effects ( interface level )
!  gwdcloc  : horizontal location where the gwdc scheme is activated.
!  break    : horizontal location where wave breaking is occurred.
!  critic   : horizontal location where critical level filtering is
!             occurred.
!  dogwdc   : logical flag whether the gwdc parameterization is
!             calculated at a grid point or not.
!
!  dogwdc is used in order to lessen cpu time for gwdc calculation.
!
!-----------------------------------------------------------------------

         integer i, ii, k, kk, kb, ilev, npt(my_max), kcb, kcldm(my_max), npr
         integer, dimension(ix, my_max) :: ipt

!lzl      real(kind=kind_phys) tem, tem1,  tem2, qtem, wtgwc, tauct,
!lzl     &                     windcltop,  shear, nonlinct, nonlin, nonlins,
!lzl     &                     n2,   dtdp,  crit1, crit2, pi, p1, p2,
!lzl     &                     gsqr,  onebg
!lzl!    &                     taus, n2,   dtdp,  crit1, crit2, pi, p1, p2

!lzl +add -------- ##
         real tem, tem1, tem2, qtem, wtgwc, tauct, &
            !lzl c           windcltop,  shear, nonlinct, nonlin, nonlins,            &
            shear, nonlinct, nonlin, nonlins, &
            n2, dtdp, crit1, crit2, pi, p1, p2, &
            gsqr, onebg, wrk, rimin, tauctxl, tauctyl, windcltop
!lzl -end -------- ##
         integer :: kcldtop(ix, my_max), kcldbot(ix, my_max)
         logical :: do_gwc(ix, my_max)
!lzl      real(kind=kind_phys), allocatable :: tauctxl(:), tauctyl(:),
!lzl     &                                     gwdcloc(:), break(:),
!lzl     &                                     critic(:),
!lzl!    &                                     critic(:),  angle(:),
!lzl     &                                     cosphi(:),  sinphi(:),
!lzl     &                                     xstress(:), ystress(:),
!lzl     &                                     ucltop(:),  vcltop(:),
!lzl     &                                     wrk(:),
!lzl     &                                     dlen(:),       gqmcldlen(:)
!lzl!     real(kind=kind_phys), allocatable :: plnint(:,:),   dpint(:,:),
!lzl!    &                                     taugwci(:,:),  taugwcxi(:,:),
!lzl!    &                                     taugwcyi(:,:), bruni(:,:),
!lzl!    &                                     taugwcyi(:,:), bruni(:,:),
!lzl      real(kind=kind_phys), allocatable :: plnint(:,:),
!lzl     &                                     taugwci(:,:),  bruni(:,:),
!lzl     &                                     rhoi(:,:),     basicui(:,:),
!lzl     &                                     ti(:,:),       riloc(:,:),
!lzl     &                                     rimin(:,:),    pint(:,:)
!lzl!     real(kind=kind_phys), allocatable :: ugwdc(:,:),    vgwdc(:,:),
!lzl      real(kind=kind_phys), allocatable ::
!lzl!    &                                     plnmid(:,:),   wtgwc(:,:),
!lzl     &                                     plnmid(:,:),
!lzl     &                                     utgwcl(:,:),   vtgwcl(:,:),
!lzl     &                                     basicum(:,:),  u(:,:),v(:,:),
!lzl     &                                     t(:,:),        spfh(:,:),
!lzl     &                                     pmid(:,:),     dpmid(:,:),
!lzl!    &                                     pmid(:,:),     cumchr(:,:),
!lzl     &                                     brunm(:,:),    rhom(:,:)
!lzl!
!lzl +add -------- ##
         real, dimension(ix, my_max) :: cosphi, sinphi, ucltop, vcltop, &
                                xstress, ystress, dtfac
         !real, dimension(ix) :: gwdcloc, break, &
         !                       critic, &
         !                       !    &                                     critic(:),  angle(:),
         !                       tmpwindcltop
         real, dimension(ix, my_max) :: dlen, gqmcldlen
         real, dimension(ix, km, my_max) :: plnmid, utgwcl, vtgwcl, basicum, &
                                    u, v, t, spfh, pmid, dpmid, brunm, rhom, velco, taugw
         real, dimension(ix, km + 1, my_max) :: taugwci, bruni, rhoi, basicui, ti, riloc, &
                                        pint
         real, dimension(ix, 2:km + 1, my_max) :: plnint
!     real(kind=kind_phys), allocatable :: plnint(:,:),   dpint(:,:),
!    &                                     taugwci(:,:),  taugwcxi(:,:),
!    &                                     taugwcyi(:,:), bruni(:,:),
!    &                                     taugwcyi(:,:), bruni(:,:),
!     real(kind=kind_phys), allocatable :: ugwdc(:,:),    vgwdc(:,:),
!    &                                     plnmid(:,:),   wtgwc(:,:),
!    &                                     pmid(:,:),     cumchr(:,:),

!lzl -end -------- ##
!-----------------------------------------------------------------------
!
!  ucltop    : zonal wind at the cloud top ( mid level )
!  vcltop    : meridional wind at the cloud top ( mid level )
!  windcltop : wind speed at the cloud top ( mid level )
!  shear     : vertical shear of basic wind
!  cosphi    : cosine of angle of wind vector at the cloud top
!  sinphi    : sine   of angle of wind vector at the cloud top
!  c1        : tunable parameter
!  c2        : tunable parameter
!  dlength   : grid spacing in the direction of basic wind at the cloud top
!  nonlinct  : nonlinear parameter at the cloud top
!  nonlin    : nonlinear parameter above the cloud top
!  nonlins   : saturation nonlinear parameter
!  taus      : saturation gravity wave drag == taugwci(i,k)
!  n2        : square of brunt-vaisala frequency
!  dtdp      : dt/dp
!  xstress   : vertically integrated zonal momentum change due to gwdc
!  ystress   : vertically integrated meridional momentum change due to gwdc
!  crit1     : variable 1 for checking critical level
!  crit2     : variable 2 for checking critical level
!
!-----------------------------------------------------------------------

!lzl      real(kind=kind_phys), parameter ::
!lzl     &                      c1=1.41,          c2=-0.38,     ricrit=0.25
!lzl     &,                     n2min=1.e-32,     zero=0.0,     one=1.0
!lzl     &,                     taumin=1.0e-20,   tauctmax=-5.
!lzl     &,                     qmin=1.0e-10,     shmin=1.0e-20
!lzl     &,                     rimax=1.0e+20,    rimaxm=0.99e+20
!lzl     &,                     rimaxp=1.01e+20,  rilarge=0.9e+20
!lzl     &,                     riminx=-1.0e+20,  riminm=-1.01e+20
!lzl     &,                     riminp=-0.99e+20, rismall=-0.9e+20
!
!lzl +add -------- ##
         real, parameter :: &
            c1 = 1.41, c2 = -0.38, ricrit = 0.25, &
            n2min = 1.e-32, zero = 0.0, one = 1.0, &
            !           taumin=1.0e-20,   tauctmax=-5.,                   &
            taumin = 1.0e-20, tauctmax = -20., &
            qmin = 1.0e-10, shmin = 1.0e-20, &
            rimax = 1.0e+20, rimaxm = 0.99e+20, &
            rimaxp = 1.01e+20, rilarge = 0.9e+20, &
            riminx = -1.0e+20, riminm = -1.01e+20, &
            riminp = -0.99e+20, rismall = -0.9e+20
         
         
         ! for GPU porting
         integer :: async_id = 1
!
         iii = 1132
         kkk = 24
!
!       if(myrank.eq.0) print*,'inside gwdc '

         
         !$acc data create(myim, npt, ipt, u, v, t, spfh, rhom, &
         !$acc&     plnmid, utgwcl, vtgwcl, basicum, taugwci, &
         !$acc&     bruni, rhoi, ti, basicui, riloc, plnint, kcldtop, &
         !$acc&     kcldbot, gqmcldlen, kcldm, &
         !$acc&     cosphi, sinphi, do_gwc, velco, dtfac, &
         !$acc&     taugw, xstress, ystress) async(async_id)
         
         

         !$acc parallel loop async(async_id) gang private(j, npr)
         do jj = 1, jlistnum
            j = jlist1(jj)
            im = nxjp(j)
            myim(jj) = im
            npr = 0
            !$acc loop seq
            do i = 1, im
               ipt(i, jj) = 0
               if (kcnv(i, jj) /= 0 .and. qmax(i, jj) > zero) then
                  !!$acc atomic capture
                  npr = npr + 1
                  ipt(npr, jj) = i
                  !!$acc end atomic
!              if(i.eq.iii)then
!                print*,' npt = ',npt
!                iinpt = npt
!              endif
               end if
            end do
            npt(jj) = npr
         end do
!

         !$acc parallel loop async(async_id) gang collapse(2)
         do jj = 1, jlistnum
            do k = 1, km + 1
               !$acc loop vector
               do i = 1, myim(jj)
                  if (k .eq. km + 1) then
                     tauctx(i, jj) = 0.0
                     taucty(i, jj) = 0.0
                  else
                     utgwc(i, k, jj) = 0.0
                     vtgwc(i, k, jj) = 0.0
                  end if
!         brunm(i,k) = 0.0
!         rhom(i,k)  = 0.0
               end do
            end do
         end do

!
!      if (npt == 0) print *,' npt =0 no gwdc calculation done! '
            !if (npt == 0) return      ! no gwdc calculation done!

!***********************************************************************
!
!  begin gwdc
!
!***********************************************************************

!-----------------------------------------------------------------------
!        write out incoming variables
!-----------------------------------------------------------------------
            fhourpr = zero
!     if (lprnt) then
!       if (fhour >= fhourpr) then
!         print *,' '
!         write(*,*) 'inside gwdc raw input start print at fhour = ',
!    &               fhour
!-------- u1 v1 t1 ----------
!         write(*,9130)
!         do ilev=km,1,-1
!           write(*,9140) ilev,u1(ipr,ilev),v1(ipr,ilev),t1(ipr,ilev)
!         enddo
!
!         print *,' '
!         print *,' inside gwdc raw input end print'
!       endif
!     endif

!9100 format(//,14x,'pressure levels',//,
!    +' ilev',6x,'pint1',7x,'pmid1',6x,'dpmid1',/)
!9110 format(i4,2x,f10.3)
!9120 format(i4,12x,2(2x,f10.3))
!9130 format(//,' ilev',7x,'u1',10x,'v1',10x,'t1',/)
!9140 format(i4,3(2x,f10.3))
!
!     allocate local arrays

!-----------------------------------------------------------------------
!        create local arrays with reversed vertical indices
!        and initialize local variables
!-----------------------------------------------------------------------
            gsqr = grav*grav
            onebg = one/grav

!lzl c      if (lprnt) then
!lzl c        npr = 1
!lzl c        do i=1,npt
!lzl c          if (ipr == ipt(i))then
!lzl c            npr = i
!lzl c            exit
!lzl c          endif
!lzl c        enddo
!lzl c      endif
         !$acc parallel loop async(async_id) gang collapse(2) private(ii)
         do jj = 1, jlistnum
            do k = 1, km + 1
!byl        k1 = km - k + 1
               !$acc loop vector
               do i = 1, npt(jj)
                  ii = ipt(i, jj)
                  if (k .le. km) then
                     u(i, k, jj) = u1(ii, k, jj)
                     v(i, k, jj) = v1(ii, k, jj)
                     t(i, k, jj) = t1(ii, k, jj)
                     spfh(i, k, jj) = max(q1(ii, k, jj), qmin)
   !         cumchr(i,k)   = cumchr1(ii,k)

                     rhom(i, k, jj) = pmid1(ii, k, jj)/ &
                                      (rd*t(i, k, jj)*(1.0 + fv*spfh(i, k, jj)))
                     plnmid(i, k, jj) = log(pmid1(ii, k, jj))
   !         ugwdc(i,k)    = zero
   !         vgwdc(i,k)    = zero
                  end if
                  taugwci(i, k, jj) = zero
                  if (k .ge. 2) plnint(i, k, jj) = log(pint1(ii, k, jj))
               end do
            end do
         end do
         
         
         pi = 2.*asin(1.)
   !

!     if (lprnt) then
!       if (fhour.ge.fhourpr) then
!
!           write(*,9201) kcnv(i),kcldbot(i),kcldtop(i)
!         enddo
!      endif
!
!       endif
!     endif
!
!9200  format(//,'  inside gwdc local variables start print',//,     &
!       2x,'kcnv',2x,'kcldbot',2x,'kcldtop',//)
!9201  format(i4,2x,i5,4x,i5)

!***********************************************************************


!-----------------------------------------------------------------------
!
!                              pressure variables
!
!  interface 1 ======== pint(1)           *********
!  mid-level 1 --------          pmid(1)            dpmid(1)
!            2 ======== pint(2)           dpint(2)
!            2 --------          pmid(2)            dpmid(2)
!            3 ======== pint(3)           dpint(3)
!            3 --------          pmid(3)            dpmid(3)
!            4 ======== pint(4)           dpint(4)
!            4 --------          pmid(4)            dpmid(4)
!              ........
!           17 ======== pint(17)          dpint(17)
!           17 --------          pmid(17)           dpmid(17)
!           18 ======== pint(18)          dpint(18)
!           18 --------          pmid(18)           dpmid(18)
!           19 ======== pint(19)          *********
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!                              thermal variables
!
!  interface 1 ========       ti(1)           rhoi(1)            bruni(1)
!            1 -------- t(1)         rhom(1)           brunm(1)
!            2 ========       ti(2)           rhoi(2)            bruni(2)
!            2 -------- t(2)         rhom(2)           brunm(2)
!            3 ========       ti(3)           rhoi(3)            bruni(3)
!            3 -------- t(3)         rhom(3)           brunm(3)
!            4 ========       ti(4)           rhoi(4)            bruni(4)
!            4 -------- t(4)         rhom(4)           brunm(4)
!              ........
!           17 ========
!           17 -------- t(17)        rhom(17)          brunm(17)
!           18 ========       ti(18)          rhoi(18)           bruni(18)
!           18 -------- t(18)        rhom(18)          brunm(18)
!           19 ========       ti(19)          rhoi(19)           bruni(19)
!

!-----------------------------------------------------------------------
!
!  calculate interface level temperature, density, and brunt-vaisala
!  frequencies based on linear interpolation of temp in ln(pressure)
!
!-----------------------------------------------------------------------
         !$acc parallel loop async(async_id) gang collapse(2) private(tem1, &
         !$acc&         tem2, qtem, dtdp, n2)
         do jj = 1, jlistnum
            do k = 1, km + 1
               !$acc loop vector
               do i = 1, npt(jj)
                  ii = ipt(i, jj)
                  if (k .eq. 1) then
                     kcldtop(i, jj) = km - ktop(ii, jj) + 1
                     kcldbot(i, jj) = km - kbot(ii, jj) + 1
      !                                    (g*qmax(ii)*cldf(ii)*dlength(ii))
                     gqmcldlen(i, jj) = grav*qmax(ii, jj)*cldf(ii, jj)*dlength(ii, jj)
!
!  top interface temperature is calculated assuming an isothermal
!  atmosphere above the top mid level.

                     ti(i, 1, jj) = t(i, 1, jj)
                     rhoi(i, 1, jj) = pint1(ii, 1, jj)/(rd*ti(i, 1, jj))
                     bruni(i, 1, jj) = sqrt(gsqr/(cp*ti(i, 1, jj)))
                  elseif (k .eq. km + 1) then
!
!  bottom interface temperature is calculated assuming an isothermal
!  atmosphere below the bottom mid level

                     ti(i, km + 1, jj) = t(i, km, jj)
                     rhoi(i, km + 1, jj) = pint1(ii, km + 1, jj)/ &
                                           (rd*ti(i, km + 1, jj)*(1.0 + fv*spfh(i, km, jj)))
                     bruni(i, km + 1, jj) = sqrt(gsqr/(cp*ti(i, km + 1, jj)))
                  else

                     tem1 = (plnmid(i, k, jj) - plnint(i, k, jj))/ &
                            (plnmid(i, k, jj) - plnmid(i, k - 1, jj))
                     tem2 = one - tem1
                     ti(i, k, jj) = t(i, k - 1, jj)*tem1 + t(i, k, jj)*tem2
                     qtem = spfh(i, k - 1, jj)*tem1 + spfh(i, k, jj)*tem2
                     rhoi(i, k, jj) = pint1(ii, k, jj)/(rd*ti(i, k, jj)*(1.0 + fv*qtem))
                     dtdp = (t(i, k, jj) - t(i, k - 1, jj))/ &
                            (pmid1(ii, k, jj) - pmid1(ii, k - 1, jj))
                     n2 = gsqr/ti(i, k, jj)*(1./cp - rhoi(i, k, jj)*dtdp)
                     bruni(i, k, jj) = sqrt(max(n2min, n2))
                  end if
               end do
            end do
         end do

            !deallocate (spfh)
!-----------------------------------------------------------------------
!
!  determine the mid-level brunt-vaisala frequencies.
!             based on interpolated interface temperatures [ ti ]
!
!-----------------------------------------------------------------------

         
         !!$acc parallel loop async(async_id) gang collapse(2) private(dtdp, n2, brunm)
         !do jj = 1, jlistnum
         !   do k = 1, km
         !      !$acc loop vector
         !      do i = 1, npt(jj)
         !         dtdp = (ti(i, k + 1, jj) - ti(i, k, jj))/ &
         !                (pint(i, k + 1, jj) - pint(i, k, jj))
         !         n2 = gsqr/t(i, k, jj)*(1./cp - rhom(i, k, jj)*dtdp)
         !         brunm = sqrt(max(n2min, n2))
         !      end do
         !   end do
         !end do

!-----------------------------------------------------------------------
!        printout
!-----------------------------------------------------------------------

!     if (lprnt) then
!       if (fhour.ge.fhourpr) then

!-------- pressure levels ----------
!
!
!         write(*,9101)
!         do ilev=1,km
!           write(*,9111) ilev,(0.01*pint(ipr,ilev)),        &
!                              (0.01*dpint(ipr,ilev)),plnint(ipr,ilev)
!           write(*,9121) ilev,(0.01*pmid(ipr,ilev)),        &
!                              (0.01*dpmid(ipr,ilev)),plnmid(ipr,ilev)
!         enddo
!         ilev=km+1
!         write(*,9111) ilev,(0.01*pint(ipr,ilev)),          &
!                            (0.01*dpint(ipr,ilev)),plnint(ipr,ilev)

!                2
!-------- u v t n  ----------
            ipr = iii
!         if(j.eq.499) then
!         write(*,9102)
!         do ilev=1,km
!           print*,'ti bruni = ',ti(iii,kkk),bruni(iii,kkk)
!           print*,'u v t brunm =',u(iii,kkk),v(iii,kkk),t(iii,kkk),brunm(iii,kkk)
!           write(*,9112) ilev,ti(ipr,ilev),(100.*bruni(ipr,ilev))
!           write(*,9122) ilev,u(ipr,ilev),v(ipr,ilev),     &
!                         t(ipr,ilev),(100.*brunm(ipr,ilev))
!         enddo
!         ilev=km+1
!          print*,'k ti bruni =',ilev,ti(ipr,ilev),bruni(ipr,ilev)
!         write(*,9112) ilev,ti(ipr,ilev),(100.*bruni(ipr,ilev))
!
!        endif

!       endif
!     endif

!9101 format(//,14x,'pressure levels',//,       &
!      ' ilev',4x,'pint',4x,'pmid',4x,'dpint',3x,'dpmid',5x,'lnp',/)
!9111 format(i4,1x,f8.2,9x,f8.2,9x,f8.2)
!9121 format(i4,9x,f8.2,9x,f8.2,1x,f8.2)
!9102 format(//' ilev',5x,'u',7x,'v',5x,'ti',7x,'t',    &
!      5x,'bruni',3x,'brunm',//)
!9112 format(i4,16x,f8.2,8x,f8.3)
!9122 format(i4,2f8.2,8x,f8.2,8x,f8.3)

!***********************************************************************
!
!        big loop over grid points                    only done if kcnv=1
!
!***********************************************************************
         !$acc parallel loop async(async_id) gang  private(kk, kb, windcltop)
         do jj = 1, jlistnum
            kcldm(jj) = 1
            !$acc loop vector
            do i = 1, npt(jj)
               kk = kcldtop(i, jj)
               kb = kcldbot(i, jj)
               !$acc atomic
               kcldm(jj) = max(kcldm(jj), kk)

!-----------------------------------------------------------------------
!
!  determine cloud top wind component, direction, and speed.
!  here, ucltop, vcltop, and windcltop are wind components and
!  wind speed at mid-level cloud top index
!
!-----------------------------------------------------------------------

!       windcltop = sqrt( ucltop(i)*ucltop(i) + vcltop(i)*vcltop(i) )
!lzl tmp        windcltop = 1.0 / sqrt( ucltop(i)*ucltop(i)                    &
!lzl tmp                             + vcltop(i)*vcltop(i) )
!
!lzl tmp        cosphi(i) = ucltop(i)*windcltop
!lzl tmp        sinphi(i) = vcltop(i)*windcltop
!
!lzl tmp 2        windcltop(i) = 1.0 / sqrt( ucltop(i)*ucltop(i)                 &
!lzl tmp 2                             + vcltop(i)*vcltop(i) )
!
               windcltop = 1.0/sqrt(u(i, kk, jj)*u(i, kk, jj) + &
                                           v(i, kk, jj)*v(i, kk, jj))

               cosphi(i, jj) = u(i, kk, jj)*windcltop
               sinphi(i, jj) = v(i, kk, jj)*windcltop

!       angle(i)  = acos(cosphi)*180./pi
!
            end do
         end do
!
!-----------------------------------------------------------------------
!
!  calculate basic state wind projected in the direction of the cloud
!  top wind.
!  input u(i,k) and v(i,k) is defined at mid level
!
!-----------------------------------------------------------------------


         !$acc parallel loop async(async_id) gang collapse(2) 
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, npt(jj)
                  basicum(i, k, jj) = u(i, k, jj)*cosphi(i, jj) + &
                                      v(i, k, jj)*sinphi(i, jj)
               end do
            end do
         end do

!-----------------------------------------------------------------------
!
!  basic state wind at interface level is also calculated
!  based on linear interpolation in ln(pressure)
!
!  in the top and bottom boundaries, basic-state wind at interface level
!  is assumed to be vertically uniform.
!
!-----------------------------------------------------------------------

!
         !$acc parallel loop async(async_id) gang collapse(2) private(tem1, &
         !$acc&         tem2, shear, tem, ii) firstprivate(grav)
         do jj = 1, jlistnum
            do k = 1, km + 1
               !$acc loop vector
               do i = 1, npt(jj)
                  ii = ipt(i, jj)
                  if (k .eq. 1) then
                     basicui(i, 1, jj) = basicum(i, 1, jj)
                  elseif (k .eq. km + 1) then
                     basicui(i, km + 1, jj) = basicum(i, km, jj)
                  else
                     tem1 = (plnmid(i, k, jj) - plnint(i, k, jj))/ &
                            (plnmid(i, k, jj) - plnmid(i, k - 1, jj))
                     tem2 = one - tem1
                     basicui(i, k, jj) = basicum(i, k, jj)*tem2 + &
                                         basicum(i, k - 1, jj)*tem1
                     
                     shear = grav*rhoi(i, k, jj)*(basicum(i, k, jj) - &
                             basicum(i, k - 1, jj)) &
                             /(pmid1(ii, k, jj) - pmid1(ii, k - 1, jj))
                     if (abs(shear) < shmin) then
                        riloc(i, k, jj) = rimax
                     else
                        tem = bruni(i, k, jj)/shear
                        riloc(i, k, jj) = tem*tem
                        if (riloc(i, k, jj) >= rimax) riloc(i, k, jj) = rilarge
                     end if

!-----------------------------------------------------------------------
!
!  calculate local richardson number
!
!  basicum   : u at mid level
!  basicui   : ui at interface level
!
!  interface 1 ========       ui(1)            rhoi(1)  bruni(1)  riloc(1)
!  mid-level 1 -------- u(1)
!            2 ========       ui(2)  dpint(2)  rhoi(2)  bruni(2)  riloc(2)
!            2 -------- u(2)
!            3 ========       ui(3)  dpint(3)  rhoi(3)  bruni(3)  riloc(3)
!            3 -------- u(3)
!            4 ========       ui(4)  dpint(4)  rhoi(4)  bruni(4)  riloc(4)
!            4 -------- u(4)
!              ........
!           17 ========       ui(17) dpint(17) rhoi(17) bruni(17) riloc(17)
!           17 -------- u(17)
!           18 ========       ui(18) dpint(18) rhoi(18) bruni(18) riloc(18)
!           18 -------- u(18)
!           19 ========       ui(19)           rhoi(19) bruni(19) riloc(19)
!
!-----------------------------------------------------------------------
                  end if
               end do
            end do
         end do
         !$acc parallel loop async(async_id) gang
         do jj = 1, jlistnum
            !$acc loop vector
            do i = 1, npt(jj)
               riloc(i, 1, jj) = riloc(i, 2, jj)
               riloc(i, km + 1, jj) = riloc(i, km, jj)
            end do
         end do

!

!     if (lprnt.and.(i.eq.ipr)) then
!       if (fhour.ge.fhourpr) then

!         write(*,9104) ucltop,vcltop,windcltop,angle,kk
!         do ilev=1,km
!           write(*,9114) ilev,basicui(ipr,ilev),dpint(ipr,ilev),   &
!            rhoi(ipr,ilev),(100.*bruni(ipr,ilev)),riloc(ilev)
!           write(*,9124) ilev,(basicum(ipr,ilev))
!         enddo
            ilev = km + 1
!         write(*,9114) ilev,basicui(ipr,ilev),dpint(ipr,ilev),     &
!            rhoi(ipr,ilev),(100.*bruni(ipr,ilev)),riloc(ilev)

!       endif
!     endif

!9104 format(//,'wind vector at cloudtop = (',f6.2,' , ',f6.2,' ) = ', &
!      f6.2,' in direction ',f6.2,4x,'kk = ',i2,//,                    &
!      ' ilev',2x,'basicum',2x,'basicui',4x,'dpint',6x,'rhoi',5x,      &
!      'bruni',6x,'ri',/)
9114        format(i4, 10x, f8.2, 4(2x, f8.2))
9124        format(i4, 1x, f8.2)

!-----------------------------------------------------------------------
!
!  calculate gravity wave stress at the interface level cloud top
!
!  kcldtopi  : the interface level cloud top index
!  kcldtop   : the midlevel cloud top index
!  kcldbot   : the midlevel cloud bottom index
!
!  a : find deep convective heating rate maximum
!
!      if kcldtop(i) is less than kcldbot(i) in a horizontal grid point,
!      it can be thought that there is deep convective cloud. however,
!      deep convective heating between kcldbot and kcldtop is sometimes
!      zero in spite of kcldtop less than kcldbot. in this case,
!      maximum deep convective heating is assumed to be 1.e-30.
!
!  b : kk is the vertical index for interface level cloud top
!
!  c : total convective fractional cover (cldf) is used as the
!      convective cloud cover for gwdc calculation instead of
!      convective cloud cover in each layer (concld).
!                       a1 = cldf*dlength
!      you can see the difference between cldf(i) and concld(i)
!      in (4.a.2) in description of the ncar community climate
!      model (ccm3).
!      in ncar ccm3, cloud fractional cover in each layer in a deep
!      cumulus convection is determined assuming total convective
!      cloud cover is randomly overlapped in each layer in the
!      cumulus convection.
!
!  d : wave stress at cloud top is calculated when the atmosphere
!      is dynamically stable at the cloud top
!
!  e : cloud top wave stress and nonlinear parameter are calculated
!      using density, temperature, and wind that are defined at mid
!      level just below the interface level in which cloud top wave
!      stress is defined.
!      nonlinct is defined at the interface level.
!
!  f : if the atmosphere is dynamically unstable at the cloud top,
!      gwdc calculation in current horizontal grid is skipped.
!
!  g : if mean wind at the cloud top is less than zero, gwdc
!      calculation in current horizontal grid is skipped.
!
!  h : maximum cloud top stress, tauctmax =  -20 n m^(-2),
!  h : max stress -5 (*j*)5/2015 tauctmax =  - 5 n m^(-2),
!      in order to prevent numerical instability.
!
!-----------------------------------------------------------------------
!
         
         !$acc parallel loop async(async_id) collapse(2) private(kk, tem, tem1, &
         !$acc&         nonlinct, tem2, tauct, ii, rimin, tauctxl, tauctyl, &
         !$acc&         crit1, crit2, nonlin, nonlins)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. npt(jj)) then
                  kk = kcldtop(i, jj)
                  ii = ipt(i, jj)
                  rimin= 0.
                  tauctxl = zero
                  tauctyl = zero
                  if (abs(basicui(i, kk, jj)) > zero .and. riloc(i, kk, jj) > ricrit) then
   !
                     tem = basicum(i, kk, jj)
                     tem1 = tem*tem
                     nonlinct = gqmcldlen(i, jj)/(bruni(i, kk, jj)*t(i, kk, jj)*tem1)    ! mu
                     tem2 = c2*nonlinct
   !                                  rhou^3c1(c2mu)^2/ndx
                     tauct = -rhom(i, kk, jj)*tem*tem1*c1*tem2*tem2 &
                             /(bruni(i, kk, jj)*dlength(ii, jj))

                     tauct = max(tauctmax, tauct)
                     tauctxl = tauct*cosphi(i, jj)           ! x stress at cloud top
                     tauctyl = tauct*sinphi(i, jj)           ! y stress at cloud top
                     taugwci(i, kk, jj) = tauct                                    !  *1
                     do_gwc(i, jj) = .true.
   !
                  else
   !
                     do_gwc(i, jj) = .false.
                  end if
                  tauctx(ii, jj) = tauctxl
                  taucty(ii, jj) = tauctyl

!

!       if (lprnt.and.(i.eq.ipr)) then
!         if (fhour.ge.fhourpr) then
!lzl write out-------------#
!            write(*,9210) tauctx(ipr),taucty(ipr),tauct(ipr),angle,kk
!         endif
!       endif

9210        format(/, 5x, 'stress vector = ( ', f8.3, ' , ', f8.3, ' ) = ', f8.3, &
                    ' in direction ', f6.2, 4x, 'kk = ', i2,/)
!-----------------------------------------------------------------------
!
!  at this point, mean wind at the cloud top is larger than zero and
!  local ri at the cloud top is larger than ricrit (=0.25)
!
!  calculate minimum of richardson number including both basic-state
!  condition and wave effects.
!
!          g*q_0*alpha*dx                  ri_loc*(1 - mu*|c2|)
!  mu  =  ----------------  ri_min =  -----------------------------
!           c_p*n*t*u^2                (1 + mu*ri_loc^(0.5)*|c2|)^2
!
!  minimum ri is calculated for the following two cases
!
!  (1)   riloc < 1.e+20
!  (2)   riloc = 1.e+20  ----> vertically uniform basic-state wind
!
!  riloc cannot be smaller than zero because n^2 becomes 1.e-32 in the
!  case of n^2 < 0.. thus the sign of rinum is determined by
!  1 - nonlin*|c2|.
!
!-----------------------------------------------------------------------

                  !$acc loop seq
                  do k = kcldm(jj), 1, -1
                     if (do_gwc(i, jj)) then
                        if (k > kk) cycle
                        if (k /= 1) then
                           tem1 = (u(i, k, jj) + u(i, k - 1, jj))*0.5
                           tem2 = (v(i, k, jj) + v(i, k - 1, jj))*0.5
                           crit1 = u(i, kk, jj)*tem1
                           crit2 = v(i, kk, jj)*tem2
                           velco(i, k, jj) = tem1*cosphi(i, jj) + tem2*sinphi(i, jj)
                        else
                           crit1 = u(i, kk, jj)*u(i, 1, jj)
                           crit2 = v(i, kk, jj)*v(i, 1, jj)
                           velco(i, 1, jj) = u(i, 1, jj)*cosphi(i, jj) + &
                                         v(i, 1, jj)*sinphi(i, jj)
                        end if

                        if (abs(basicui(i, k, jj)) > zero .and. crit1 > zero &
                            .and. crit2 > zero) then
                           tem = basicui(i, k, jj)*basicui(i, k, jj)
                           nonlin = gqmcldlen(i, jj)/(bruni(i, k, jj)*ti(i, k, jj)*tem)
                           tem = nonlin*abs(c2)
                           if (riloc(i, k, jj) < rimaxm) then
                              tem1 = 1 + tem*sqrt(riloc(i, k, jj))
                              rimin = riloc(i, k, jj)*(1 - tem)/(tem1*tem1)
                           else if ((riloc(i, k, jj) > rimaxm) .and. &
                                    (riloc(i, k, jj) < rimaxp)) then
                              rimin = (1 - tem)/(tem*tem)
                           end if
                           if (rimin <= riminx) then
                              rimin = rismall
                           end if
                        else
                           rimin = riminx
                        end if

!-----------------------------------------------------------------------
!
!  if minimum ri at interface cloud top is less than or equal to 1/4,
!  gwdc calculation for current horizontal grid is skipped
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  calculate gravity wave stress profile using the wave saturation
!  hypothesis of lindzen (1981).
!
!  assuming kcldtop(i)=10 and kcldbot=16,
!
!                             taugwci  riloc  rimin   utgwc
!
!  interface 1 ========       - 0.001         -1.e20
!            1 --------                               0.000
!            2 ========       - 0.001         -1.e20
!            2 --------                               0.000
!            3 ========       - 0.001         -1.e20
!            3 --------                               -.xxx
!            4 ========       - 0.001  2.600  2.000
!            4 --------                               0.000
!            5 ========       - 0.001  2.500  2.000
!            5 --------                               0.000
!            6 ========       - 0.001  1.500  0.110
!            6 --------                               +.xxx
!            7 ========       - 0.005  2.000  3.000
!            7 --------                               0.000
!            8 ========       - 0.005  1.000  0.222
!            8 --------                               +.xxx
!            9 ========       - 0.010  1.000  2.000
!            9 --------                               0.000
! kcldtopi  10 ========  $$$  - 0.010
! kcldtop   10 --------  $$$                          yyyyy
!           11 ========  $$$  0
!           11 --------  $$$
!           12 ========  $$$  0
!           12 --------  $$$
!           13 ========  $$$  0
!           13 --------  $$$
!           14 ========  $$$  0
!           14 --------  $$$
!           15 ========  $$$  0
!           15 --------  $$$
!           16 ========  $$$  0
! kcldbot   16 --------  $$$
!           17 ========       0
!           17 --------
!           18 ========       0
!           18 --------
!           19 ========       0
!
!-----------------------------------------------------------------------
!
!   even though the cloud top level obtained in deep convective para-
!   meterization is defined in mid-level, the cloud top level for
!   the gwdc calculation is assumed to be the interface level just
!   above the mid-level cloud top vertical level index.
!
!-----------------------------------------------------------------------

                        if (k < kk .and. k > 1) then
                           if (abs(taugwci(i, k + 1, jj)) > taumin) then                  ! taugwci
                              if (riloc(i, k, jj) > ricrit) then                         ! riloc
                                 if (rimin > ricrit) then                       ! rimin
                                    taugwci(i, k, jj) = taugwci(i, k + 1, jj)
                                 elseif (rimin > riminp) then
                                    tem = 2.0 + 1.0/sqrt(riloc(i, k, jj))
                                    nonlins = (1.0/abs(c2))*(2.*sqrt(tem) - tem)
                                    tem1 = basicui(i, k, jj)
                                    tem2 = c2*nonlins*tem1
                                    taugwci(i, k, jj) = -rhoi(i, k, jj)*c1*tem1*tem2*tem2 &
                                                    /(bruni(i, k, jj)*dlength(ii, jj))
                                 elseif (rimin > riminm) then
                                    taugwci(i, k, jj) = zero
                                 end if                                              ! rimin
                              else

!!!!!!!!!! in the dynamically unstable environment, there is no gravity wave stress

                                 taugwci(i, k, jj) = zero
                              end if                                                ! riloc
                           else
                              taugwci(i, k, jj) = zero
                           end if                                                  ! taugwci

                           if ((basicum(i, k + 1, jj)*basicum(i, k, jj)) < 0.) then
                              taugwci(i, k + 1, jj) = zero
                              taugwci(i, k, jj) = zero
                           end if

                           if (abs(taugwci(i, k, jj)) > abs(taugwci(i, k + 1, jj))) then
                              taugwci(i, k, jj) = taugwci(i, k + 1, jj)
                           end if

                        elseif (k == 1) then

!!!!!! upper boundary condition - permit upward propagation of gravity wave energy

                           taugwci(i, 1, jj) = taugwci(i, 2, jj)
                        end if
                     end if
                  end do                     ! end of i=1,npt loop

                  dtfac(i, jj) = 1.0
                  xstress(i, jj) = zero
                  ystress(i, jj) = zero
               end if
            end do
         end do
         
         !$acc parallel loop async(async_id) gang private(kk, tem)
         do jj = 1, jlistnum
            !$acc loop vector
            do i = 1, npt(jj)
               !$acc loop seq
               do k = 1, km
                  ii = ipt(i, jj)
                  if (do_gwc(i, jj)) then
                     kk = kcldtop(i, jj)
                     if (k < kk) then
                        taugw(i, k, jj) = (taugwci(i, k + 1, jj) - taugwci(i, k, jj))/ &
                                      dpmid1(ii, k, jj)/onebg
                        if (taugw(i, k, jj) /= 0.0) then
                           tem = deltim*taugw(i, k, jj)
                           !!$acc atomic
                           dtfac(i, jj) = min(dtfac(i, jj), abs(velco(i, k, jj)/tem))
                        end if
                     else
                        taugw(i, k, jj) = 0.0
                     end if
                  else
                     taugw(i, k, jj) = 0.0
                  end if
               end do
            end do
         end do

!!!!!! Vertical differentiation
!!!!!!
!

         !$acc parallel loop async(async_id) gang private(kk, wtgwc)
         do jj = 1, jlistnum
            !$acc loop vector
            do i = 1, npt(jj)
               !$acc loop seq
               do k = 1, km
                  ii = ipt(i, jj)
                  utgwcl(i, k, jj) = zero
                  vtgwcl(i, k, jj) = zero
                  if (do_gwc(i, jj)) then
                     kk = kcldtop(i, jj)
                     if (k < kk) then
!              wtgwc       = (taugwci(i,k+1) - taugwci(i,k)) / dpmid(i,k)
                        wtgwc = taugw(i, k, jj)*dtfac(i, jj)
                        utgwcl(i, k, jj) = wtgwc*cosphi(i, jj)
                        vtgwcl(i, k, jj) = wtgwc*sinphi(i, jj)
                     end if
!-----------------------------------------------------------------------
!
!  calculate momentum flux = stress deposited above cloup top
!  apply equal amount with opposite sign within cloud
!
!-----------------------------------------------------------------------
                     if (k .le. kcldm(jj)) then
                        !!$acc atomic
                        xstress(i, jj) = xstress(i, jj) + utgwcl(i, k, jj)*dpmid1(ii, k, jj)*onebg
                        !!$acc atomic
                        ystress(i, jj) = ystress(i, jj) + vtgwcl(i, k, jj)*dpmid1(ii, k, jj)*onebg
                     end if
                  end if
               end do
            end do
         end do
!-----------------------------------------------------------------------
!        alt 1      only uppermost layer
!-----------------------------------------------------------------------

!     kk = kcldtop(i)
!     tem1 = g / dpmid(i,kk)
!     utgwc(i,kk) = - tem1 * xstress
!     vtgwc(i,kk) = - tem1 * ystress

!-----------------------------------------------------------------------
!        alt 2      sin(kt-kb)
!-----------------------------------------------------------------------

!
         !$acc parallel loop async(async_id) gang collapse(2) private(ii, wrk, kk, p1, p2, tem)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, npt(jj)
                  ii = ipt(i, jj)
                  if (do_gwc(i, jj)) then
                     kk = kcldtop(i, jj)
                     if (k >= kk .and. k <= kcldbot(i, jj)) then
                        wrk = 0.5*pi/(pint1(ii, kcldbot(i, jj) + 1, jj) - &
                                         pint1(ii, kk, jj))
                        p1 = sin(wrk*(pint1(ii, k, jj) - pint1(ii, kk, jj)))
                        p2 = sin(wrk*(pint1(ii, k + 1, jj) - pint1(ii, kk, jj)))
                        tem = -(p2 - p1)/dpmid1(ii, k, jj)/onebg
                        utgwcl(i, k, jj) = tem*xstress(i, jj)
                        vtgwcl(i, k, jj) = tem*ystress(i, jj)
                     end if
                  end if
                  utgwc(ii, k, jj) = utgwcl(i, k, jj)
                  vtgwc(ii, k, jj) = vtgwcl(i, k, jj)
               end do
            end do
         end do
         !$acc end data
!-----------------------------------------------------------------------
!        alt 3      from kt to kb  proportional to conv heating
!-----------------------------------------------------------------------

!     do k=kcldtop(i),kcldbot(i)
!     p1=cumchr(i,k)
!     p2=cumchr(i,k+1)
!     utgwcl(i,k) = - g*xstress*(p1-p2)/dpmid(i,k)
!     enddo
!-----------------------------------------------------------------------
!
!  the gwdc should accelerate the zonal and meridional wind in the
!  opposite direction of the previous zonal and meridional wind,
!  respectively
!
!-----------------------------------------------------------------------
!     do k=1,kcldtop(i)-1
!      if (utgwcl(i,k)*u(i,k) .gt. 0.0) then
!-------------------- x-component-------------------
!       write(6,'(a)')
!    +  '(gwdc) warning: the gwdc should accelerate the zonal wind '
!       write(6,'(a,a,i3,a,i3)')
!    +  'in the opposite direction of the previous zonal wind',
!    +  ' at i = ',i,' and j = ',lat
!       write(6,'(4(1x,e17.10))') u(i,kk),v(i,kk),u(i,k),v(i,k)
!       write(6,'(a,1x,e17.10))') 'vcld . v =',
!    +  u(i,kk)*u(i,k)+v(i,kk)*v(i,k)
!       if(u(i,kcldtop(i))*u(i,k)+v(i,kcldtop(i))*v(i,k).gt.0.0)then
!       do k1=1,km
!         write(6,'(i2,36x,2(1x,e17.10))')
!    +             k1,taugwcxi(i,k1),taugwci(i,k1)
!         write(6,'(i2,2(1x,e17.10))') k1,utgwcl(i,k1),u(i,k1)
!       end do
!       write(6,'(i2,36x,1x,e17.10)') (km+1),taugwcxi(i,km+1)
!       end if
!-------------------- along wind at cloud top -----
!       do k1=1,km
!         write(6,'(i2,36x,2(1x,e17.10))')
!    +             k1,taugwci(i,k1)
!         write(6,'(i2,2(1x,e17.10))') k1,wtgwc(i,k1),basicum(i,k1)
!       end do
!       write(6,'(i2,36x,1x,e17.10)') (km+1),taugwci(i,km+1)
!      end if
!      if (vtgwc(i,k)*v(i,k) .gt. 0.0) then
!       write(6,'(a)')
!    +  '(gwdc) warning: the gwdc should accelerate the meridional wind'
!       write(6,'(a,a,i3,a,i3)')
!    +  'in the opposite direction of the previous meridional wind',
!    +  ' at i = ',i,' and j = ',lat
!       write(6,'(4(1x,e17.10))') u(i,kcldtop(i)),v(i,kcldtop(i)),
!    +                            u(i,k),v(i,k)
!       write(6,'(a,1x,e17.10))') 'vcld . v =',
!    +                    u(i,kcldtop(i))*u(i,k)+v(i,kcldtop(i))*v(i,k)
!       if(u(i,kcldtop(i))*u(i,k)+v(i,kcldtop(i))*v(i,k).gt.0.0)then
!       do k1=1,km
!         write(6,'(i2,36x,2(1x,e17.10))')
!    +                        k1,taugwcyi(i,k1),taugwci(i,k1)
!         write(6,'(i2,2(1x,e17.10))') k1,vtgwc(i,k1),v(i,k1)
!       end do
!       write(6,'(i2,36x,1x,e17.10)') (km+1),taugwcyi(i,km+1)
!       end if
!      end if
!     enddo
!1000 continue
!***********************************************************************
!     if (lprnt) then
!       if (fhour.ge.fhourpr) then
!-------- utgwc vtgwc ----------
!         write(*,9220)
!         do ilev=1,km
!           write(*,9221) ilev,(86400.*utgwcl(ipr,ilev)),
!    +                         (86400.*vtgwcl(ipr,ilev))
!         enddo
!       endif
!     endif

!9220 format(//,14x,'tendency due to gwdc',//,
!    +' ilev',6x,'utgwc',7x,'vtgwc',/)
!9221 format(i4,2(2x,f10.3))

!-----------------------------------------------------------------------
!
!  for gwdc performance analysis
!
!-----------------------------------------------------------------------
!     do k = 1, kk-1
!       do i = 1, nct
!         kk = kcldtop(i)
!         if ( (abs(taugwci(i,kk)) > taumin) ) then
!           gwdcloc(i) = one
!        if ( abs(taugwci(i,k)-taugwci(i,kk)) > taumin ) then
!         break(i) = 1.0
!         go to 2000
!        endif
!       enddo
!2000   continue

!       do k = 1, kk-1
!        if ( ( abs(taugwci(i,k)).lt.taumin ) .and.
!    &        ( abs(taugwci(i,k+1)).gt.taumin ) .and.
!    &        ( basicum(i,k+1)*basicum(i,k) .lt. 0. ) ) then
!         critic(i) = 1.0
!         print *,i,k,' inside gwdc  taugwci(k) = ',taugwci(i,k)
!         print *,i,k+1,' inside gwdc  taugwci(k+1) = ',taugwci(i,k+1)
!         print *,i,k,' inside gwdc  basicum(k) = ',basicum(i,k)
!         print *,i,k+1,' inside gwdc  basicum(k+1) = ',basicum(i,k+1)
!         print *,i,' inside gwdc  critic = ',critic(i)
!         goto 2010
!        endif
!       enddo
!2010   continue
!      endif
!     enddo
!-----------------------------------------------------------------------
!        convert back local gwdc tendency arrays to gfs model vertical indices
!        outgoing (fu1,fv1)=(utgwc,vtgwc)
!-----------------------------------------------------------------------
!         brunm(ii,kk) = brunm(i,k)
!         brunm(i,k)  = tem
!         rhom(ii,kk) = rhom(i,k)
!
!byl        k1 = km - k + 1
! still keep: 1(top)-to-60(bot)

! ---  update the wind components with  gwdc tendencies
!
!!      do k =1,km
!!        do i = 1,im
!!          u1(i,k) = u1(i,k)+utgwc(i,k)*deltim
!!          v1(i,k) = v1(i,k)+vtgwc(i,k)*deltim
!!        enddo
!!      enddo
!
         
!
!     if (lprnt) then
!       if (fhour.ge.fhourpr) then
!-------- utgwc vtgwc ----------
!lzl +add turn on write out-----#
!         write(*,9225)
!         do ilev=km,1,-1
!           write(*,9226) ilev,(86400.*fu1(ipr,ilev)),
!     +                         (86400.*fv1(ipr,ilev))
!           write(*,9226) ilev,(utgwc(ipr,ilev)),
!     +                         (vtgwc(ipr,ilev))
!         enddo
!       endif
!     endif

!9225 format(//,14x,'tendency due to gwdc - to gbphys',//,
!     +' ilev',6x,'utgwc',7x,'vtgwc',/)
!9226 format(i4,2(2x,f10.3))
!

         return
      end
