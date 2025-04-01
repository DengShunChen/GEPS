      subroutine readaeroclx(nx,my,my_max,lev,naero,julian,&
                             itimestep,monsave,&
                             ggdef,aeroclx)
!
!  read climt data from data base
!
!----------------------------------------------------------------
!  input :
!         nx      : dimension of e-w direction
!         my      : dimension of s-n direction
!         julian  : julian day
!  output :
!         aeroclx : aerosol climatology from MERRA2
!
!...............................................
!  naero-species aerosol :
!    SO4 : Sulphate Mixing Ratio
!    DU1 : Dust Mixing Ratio (bin 001)
!    DU2 : Dust Mixing Ratio (bin 002)
!    DU3 : Dust Mixing Ratio (bin 003)
!    DU4 : Dust Mixing Ratio (bin 004)
!    DU5 : Dust Mixing Ratio (bin 005)
!    SS1 : Sea Salt Mixing Ratio (bin 001)
!    SS2 : Sea Salt Mixing Ratio (bin 002)
!    SS3 : Sea Salt Mixing Ratio (bin 003)
!    SS4 : Sea Salt Mixing Ratio (bin 004)
!    SS5 : Sea Salt Mixing Ratio (bin 005)
!    BLC : Hydrophilic Black Carbon
!    BBC : Hydrophobic Black Carbon
!    OLC : Hydrophilic Organic Carbon (Particulate Matter)
!    OBC : Hydrophobic Organic Carbon (Particulate Matter)
!    MSA : Methanesulphonic acid
!    DMS : Dimethylsulphide
!    SO2 : Sulphur dioxide
!...............................................
!  nlut-species aerosol :
!    SO4 : Sulphate
!    SOT : Soot
!    IOC : Insoluble OC
!    WOC : Water soluble OC
!    SAM : Sea salt accumulation mode
!    SCM : Sea salt coarse mode
!    DNM : Dust nucleation mode
!    DAM : Dust accumulation mode
!    DCM : Dust coarse mode
!-----------------------------------------------------------------
!
      use index
      use mpe
      use rank
      use const, only: ihdgi,ifilin_aero,RTYPE, &
                       naso4,nadu1,nadu2,nadu3,nadu4,nadu5, &
                       nass1,nass2,nass3,nass4,nass5,nablc, &
                       nabbc,naolc,naobc,namsa,nadms,naso2
#ifdef LUT_aero
      use module_gocart_coupling, only : nlut,nsuso,nsoot,ninso,   &
                                         nwaso,nssam,nsscm,nminm,  &
                                         nmiam,nmicm
#endif
      use phygrid, only : aerosave1, aerosave2
!
      implicit  none

      integer   nx,my,my_max,julian,lev,naero
      integer   itimestep,monsave,monnxt
      real(kind=RTYPE) aeroclx(nxp,naero*lev,my_max)

!
!  working array as climt data base
!
      real      glob(nx,my)
      real      wrk(nxp,my_max,2)

      character blnk*1,ggdef*4,typ*6
      character(len=3) aerokey(18)
      integer   mon(12),mondy(13)
      data mon/15,46,74,105,135,166,196,227,258,288,319,349/
      data mondy/0,31,59,90,120,151,181,212,243,273,304,334,365/
      data blnk/' '/

      integer   i,j,k,m,n,jul,nxj,mm,istat,monidex,lncrec,jj,ii,nm
      real      coef1,coef2
      real*8    tm_1, tm_2, tm_use, mpi_wtime

      tm_1 = mpi_wtime()

#ifdef LUT_aero
      if ( nsuso.le.naero ) aerokey(nsuso)  = 'SO4'
      if ( nsoot.le.naero ) aerokey(nsoot)  = 'SOT'
      if ( ninso.le.naero ) aerokey(ninso)  = 'IOC'
      if ( nwaso.le.naero ) aerokey(nwaso)  = 'WOC'
      if ( nssam.le.naero ) aerokey(nssam)  = 'SAM'
      if ( nsscm.le.naero ) aerokey(nsscm)  = 'SCM'
      if ( nminm.le.naero ) aerokey(nminm) = 'DNM'
      if ( nmiam.le.naero ) aerokey(nmiam) = 'DAM'
      if ( nmicm.le.naero ) aerokey(nmicm) = 'DCM'
#else
      if ( naso4.le.naero ) aerokey(naso4) = 'SO4'
      if ( nadu1.le.naero ) aerokey(nadu1) = 'DU1'
      if ( nadu2.le.naero ) aerokey(nadu2) = 'DU2'
      if ( nadu3.le.naero ) aerokey(nadu3) = 'DU3'
      if ( nadu4.le.naero ) aerokey(nadu4) = 'DU4'
      if ( nadu5.le.naero ) aerokey(nadu5) = 'DU5'
      if ( nass1.le.naero ) aerokey(nass1) = 'SS1'
      if ( nass2.le.naero ) aerokey(nass2) = 'SS2'
      if ( nass3.le.naero ) aerokey(nass3) = 'SS3'
      if ( nass4.le.naero ) aerokey(nass4) = 'SS4'
      if ( nass5.le.naero ) aerokey(nass5) = 'SS5'
      if ( nablc.le.naero ) aerokey(nablc) = 'BLC'
      if ( nabbc.le.naero ) aerokey(nabbc) = 'BBC'
      if ( naolc.le.naero ) aerokey(naolc) = 'OLC'
      if ( naobc.le.naero ) aerokey(naobc) = 'OBC'
      if ( namsa.le.naero ) aerokey(namsa) = 'MSA'
      if ( nadms.le.naero ) aerokey(nadms) = 'DMS'
      if ( naso2.le.naero ) aerokey(naso2) = 'SO2'
#endif

      lncrec = nx*my

      jul = julian
      if ( jul .ge. 366 ) jul = 365
!
!  what month is it ?
!
      do k = 1,12
        if ( jul .gt. mondy(k) .and. jul .le. mondy(k+1) ) monidex = k
      enddo

      ! initial time , let monsave = monidex
      if ( itimestep .lt. 1 ) monsave = monidex

      if ( myrank .eq. 0 ) then
        print *, 'before readaeroclx : monsave=',monsave,&
                 ' , monidex=',monidex
      endif

#ifdef I38K
  15  format(a6,'  gbck',a4,4x,i2.2,6x)
#else
  15  format(a6,'gbck',a4,4x,i2.2,6x)
#endif

!
!  to interpolat linearly based on julian day
!
!-- climat dataset has 12 months

      if ( jul .le. mon(1) ) jul = jul + 365

      do n = 1, naero

        if ( myrank .eq. 0 ) print *, 'read aerotype : ',aerokey(n)

        do m = 1, lev
          if ( lev .lt. 100 ) then
            write(typ,'("M",i2.2,a3)') m,aerokey(n)
          else
            write(typ,'("N",i2.2,a3)') mod(m,100),aerokey(n)
          endif

          wrk = 0.

          if ( jul .gt. mon(12) ) then
            monnxt = 1
            coef1 = float(jul - mon(12))/float(380 - mon(12))
            coef2 = 1. - coef1

            ! if initial monsave = monidex = 12 :
            ! read both mon=12 and mon=1 dmsfile -> aerosave1 and aerosave2
            if ( itimestep .lt. 1 ) then
              ! mon=12 :
              write(ihdgi,15) typ,ggdef,12
              call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
              call unify_reducepick(nx,my,my_max,glob,wrk(1,1,1))
              ! mon=1 :
              write(ihdgi,15) typ,ggdef,1
              call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
              call unify_reducepick(nx,my,my_max,glob,wrk(1,1,2))

              do j = 1, my_max
                do i = 1, nxp
                  aerosave1(i,m,j,n) = wrk(i,j,1)  !mon=12
                  aerosave2(i,m,j,n) = wrk(i,j,2)  !mon=1
                enddo
              enddo

            endif  !end of if itimestep<1

            ! if monsave=11 and monidex=12 :
            ! read only mon=12 dmsfile -> aerosave2
            if ( monsave .lt. monidex ) then
              ! mon=12 :
              write(ihdgi,15) typ,ggdef,12
              call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
              call unify_reducepick(nx,my,my_max,glob,wrk(1,1,2))

              do j = 1, my_max
                do i = 1, nxp
                  aerosave1(i,m,j,n) = aerosave2(i,m,j,n)  !mon=11
                  aerosave2(i,m,j,n) = wrk(i,j,2)          !mon=12
                enddo
              enddo

            endif  !end of if monidex>monsave

            ! interpolate :
            do jj = 1,jlistnum
              j = jlist1(jj)
              nxj = nxdef_2d(j)
              do i = 1,nxj
                nm = (n - 1)*lev + m
                aeroclx(i,nm,jj) = max(0.0,coef1*aerosave2(i,m,jj,n)+ &
                                   coef2*aerosave2(i,m,jj,n))
              enddo
            enddo

          endif   !end of if jul>mon(12)

          do k = 2 , 12
            if ( jul .gt. mon(k-1) .and. jul .le. mon(k) ) then
              monnxt = k
              coef1 = float(jul-mon(k-1))/float(mon(k)-mon(k-1))
              coef2 = 1. - coef1

              ! if initial monsave = monidex = k-1 :
              ! read both mon=k-1 and mon=k dmsfile -> aerosave1 and aerosave2
              if ( itimestep .lt. 1 ) then
                ! mon=k-1 :
                write(ihdgi,15) typ,ggdef,k-1
                call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
                call unify_reducepick(nx,my,my_max,glob,wrk(1,1,1))
                ! mon=k :
                write(ihdgi,15) typ,ggdef,k
                call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
                call unify_reducepick(nx,my,my_max,glob,wrk(1,1,2))

                do j = 1, my_max
                  do i = 1, nxp
                    aerosave1(i,m,j,n) = wrk(i,j,1)  !mon=k-1
                    aerosave2(i,m,j,n) = wrk(i,j,2)  !mon=k
                  enddo
                enddo

              endif  !end of if itimestep<1

              ! if monsave=k-2 and monidex=k-1 :
              ! read only mon=k dmsfile -> aerosave2
              if ( monsave .lt. monidex ) then
                ! mm=k :
                write(ihdgi,15) typ,ggdef,k
                call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
                call unify_reducepick(nx,my,my_max,glob,wrk(1,1,2))

                do j = 1, my_max
                  do i = 1, nxp
                    aerosave1(i,m,j,n) = aerosave2(i,m,j,n)  !mon=k-1
                    aerosave2(i,m,j,n) = wrk(i,j,2)          !mon=k
                  enddo
                enddo

              endif  !end of if monsave<monidex

              ! interpolate :
              do jj = 1,jlistnum
                j = jlist1(jj)
                nxj = nxdef_2d(j)
                do i = 1,nxj
                  nm = (n - 1)*lev + m
                  aeroclx(i,nm,jj) = max(0.0,coef1*aerosave2(i,m,jj,n)+ &
                                     coef2*aerosave2(i,m,jj,n))
                enddo
              enddo

            endif  !end of if jul>mon(k-1) & jul<=mon(k)

          enddo  !end of do k=2,12

        enddo  !end of do m=1,lev

      enddo  !end of do n=1,naero

      ! update monsave
      if ( monsave .lt. monidex ) monsave = monidex

      if ( myrank .eq. 0 ) then
        print *, 'after readaeroclx : monsave=',monsave,&
                 ' , monidex=',monidex
      endif
      tm_2 = mpi_wtime()
      tm_use = tm_2-tm_1
      if (myrank .eq. 0) print 265, tm_use
  265 format('readaeroclx Timing=',f8.3,' elapse seconds')

      return
      end
