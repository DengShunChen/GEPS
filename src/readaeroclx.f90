      subroutine readaeroclx(nx,my,my_max,lev,naero,julian,&
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
!    n=1 , SO4 : Sulphate Mixing Ratio
!    n=2 , DU1 : Dust Mixing Ratio (bin 001)
!    n=3 , DU2 : Dust Mixing Ratio (bin 002)
!    n=4 , DU3 : Dust Mixing Ratio (bin 003)
!    n=5 , DU4 : Dust Mixing Ratio (bin 004)
!    n=6 , DU5 : Dust Mixing Ratio (bin 005)
!    n=7 , SS1 : Sea Salt Mixing Ratio (bin 001)
!    n=8 , SS2 : Sea Salt Mixing Ratio (bin 002)
!    n=9 , SS3 : Sea Salt Mixing Ratio (bin 003)
!    n=10, SS4 : Sea Salt Mixing Ratio (bin 004)
!    n=11, SS5 : Sea Salt Mixing Ratio (bin 005)
!    n=12, BLC : Hydrophilic Black Carbon
!    n=13, BBC : Hydrophobic Black Carbon
!    n=14, OLC : Hydrophilic Organic Carbon (Particulate Matter)
!    n=15, OBC : Hydrophobic Organic Carbon (Particulate Matter)
!    n=16, MSA : Methanesulphonic acid
!    n=17, DMS : Dimethylsulphide
!    n=18, SO2 : Sulphur dioxide
!-----------------------------------------------------------------
!
      use index
      use mpe
      use rank
      use const, only: ihdgi,ifilin_aero
!
      implicit  none

      integer   nx,my,my_max,julian,lev,naero
      real      aeroclx(nxp,naero*lev,my_max)

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
      data aerokey/'SO4','DU1','DU2','DU3','DU4','DU5', &
                   'SS1','SS2','SS3','SS4','SS5','BLC', &
                   'BBC','OLC','OBC','MSA','DMS','SO2'/

      integer   i,j,k,m,n,jul,nxj,mm,istat,monidex,lncrec,jj,ii,nm
      real      coef1,coef2

      lncrec = nx*my

      jul = julian
      if ( jul .ge. 366 ) jul = 365
!
!  what month is it ?
!
      do k = 1,12
        if ( jul .gt. mondy(k) .and. jul .le. mondy(k+1) ) monidex = k
      enddo

#ifdef I38K
  15  format(a6,'  gbck',a4,4x,i2.2,6x)
#else
  15  format(a6,'gbck',a4,4x,i2.2,6x)
#endif

!
!  to interpolat linearly based on julian day
!
!-- climat dataset has 12 months

      do n = 1, naero

        if ( jul .le. mon(1) ) jul = jul + 365

        if ( jul .gt. mon(12) ) then
          coef1 = float(jul - mon(12))/float(380 - mon(12))
          coef2 = 1. - coef1
          do m = 1, lev
            ! read dmsfile :
            if ( lev .lt. 100 ) then
              write(typ,'("M",i2.2,a3)') m,aerokey(n)
            else
              write(typ,'("N",i2.2,a3)') mod(m,100),aerokey(n)
            endif
            ! mm=12 :
            write(ihdgi,15) typ,ggdef,12
            call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
            call unify_reducepick(nx,my,my_max,glob,wrk(1,1,1))
            ! mm=1 :
            write(ihdgi,15) typ,ggdef,1
            call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
            call unify_reducepick(nx,my,my_max,glob,wrk(1,1,2))

            ! interpolate :
            do jj = 1,jlistnum
              j = jlist1(jj)
              nxj = nxdef_2d(j)
              do i = 1,nxj
                nm = (n - 1)*lev + m
                aeroclx(i,nm,jj) = coef1*wrk(i,jj,2) + &
                                   coef2*wrk(i,jj,1)
              enddo
            enddo
          enddo   !end of do m=1,lev
        endif   !end of if jul>mon(12)

        do k = 2 , 12
          if ( jul .gt. mon(k-1) .and. jul .le. mon(k) ) then
            coef1 = float(jul-mon(k-1))/float(mon(k)-mon(k-1))
            coef2 = 1. - coef1
            do m = 1, lev
              wrk = 0.
              ! read dmsfile :
              if ( lev .lt. 100 ) then
                write(typ,'("M",i2.2,a3)') m,aerokey(n)
              else
                write(typ,'("N",i2.2,a3)') mod(m,100),aerokey(n)
              endif
              ! mm=k-1 :
              write(ihdgi,15) typ,ggdef,k-1
              call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
              call unify_reducepick(nx,my,my_max,glob,wrk(1,1,1))
!              call qmaxn3_r(glob,1,1,1,nx,my,1)
              ! mm=k :
              write(ihdgi,15) typ,ggdef,k
              call dmsread(nx,my,lncrec,'H',ifilin_aero,glob,istat)
              call unify_reducepick(nx,my,my_max,glob,wrk(1,1,2))
!              call qmaxn3_r(glob,1,1,1,nx,my,1)

              ! interpolate :
              nm = (n - 1)*lev + m
              do jj = 1,jlistnum
                j = jlist1(jj)
                nxj = nxdef_2d(j)
                do i = 1,nxj
                  aeroclx(i,nm,jj) = coef1*wrk(i,jj,2) + &
                                     coef2*wrk(i,jj,1)
                enddo
              enddo
            enddo  !end of do m
          endif  !end of if jul>mon(k-1) & jul<=mon(k)
        enddo  !end of do k=2,12

      enddo  !end of do n=1,naero

      return
      end
