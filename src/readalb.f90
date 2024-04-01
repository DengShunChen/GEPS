      subroutine readalb(nx,my,my_max,julian,            &
                         alvsf,alvwf,alnsf,alnwf,facsf,facwf)
!----------------------------------------------------------------
!
!  read climt data from data base
!
!----------------------------------------------------------------
!  input :
!         nx      : dimension of e-w direction
!         my      : dimension of s-n direction
!         julian  : julian day
!
!  -- 4 months data
!         alvsfcl : mean vis albedo with strong cosz dependency
!         alvwfcl : mean vis albedo with weak cosz dependency
!         alnsfcl : mean nir albedo with strong cosz dependency
!         alnwfcl : mean nir albedo with weak cosz dependency
!                 
!  output :
! -- time intepolation to julian day
!         alvsf  : mean vis albedo with strong cosz dependency
!         alvwf  : mean vis albedo with weak cosz dependency
!         alnsf  : mean nir albedo with strong cosz dependency
!         alnwf  : mean nir albedo with weak cosz dependency
!  -- annual mean
!         facsf   : fractional coverage with strong cosz dependency
!         facwf   : fractional coverage with weak cosz dependency
!
!  By Mei-Yu Chang, 2014/09/07
!----------------------------------------------------------------
! 
      use rank
      use mpe
      use index
      use physpara  ,only : ialbflg
      use const, only: ihdgi,bckfile,ggdef

      implicit none

      integer nx,my,my_max,julian,ii

      real  alvsf(nxp,my_max),alvwf(nxp,my_max), &
            alnsf(nxp,my_max),alnwf(nxp,my_max), &
            facsf(nxp,my_max),facwf(nxp,my_max)
! local working array
      integer i,j,jj,nxj,lncrec,mm,nn,k,jul,istat
      real    coef1,coef2
      real  alvsfcl(nxp,my_max,2),alvwfcl(nxp,my_max,2),   &
            alnsfcl(nxp,my_max,2),alnwfcl(nxp,my_max,2)
      real work(nx,my)
!
      character blnk*1
      integer   mon(12),mmse(2),mmax,mmt,mon1(12),mon2(12)
      data mon1/ 74,166,258,349,  0,  0,  0,  0,  0,  0,  0,  0/
      data mon2/ 15, 46, 74,105,135,166,196,227,258,288,319,349/

      if ( ialbflg .eq. 0 ) then
        mon  = mon1
        mmax = 4
        mmt  = 3
      endif
      
      if ( ialbflg .eq. 1 ) then
        mon  = mon2
        mmax = 12
        mmt  =  1
      endif
!
      data blnk/' '/
!
      lncrec=nx*my
!
      jul=julian
      if(jul .ge. 366)jul=365
!----------------------------------------------------------------
!  to interpolat linearly based on julian day
!
!-- climate dataset has 12 months

      if(jul .le. mon(1))jul=jul+365
!
      if(jul .gt. mon(mmax))then
        mmse(1)= mmax * mmt
        mmse(2)=    1 * mmt
        coef1=float(jul-mon(mmax))/float(365+mon(1)-mon(mmax))
        coef2=1.-coef1
      else
        do k=2,mmax
          if(jul .gt. mon(k-1) .and. jul .le. mon(k))then
            mmse(1)= (k-1) * mmt
            mmse(2)=    k  * mmt
            coef1=float(jul-mon(k-1))/float(mon(k)-mon(k-1))
            coef2=1.-coef1
          endif
        enddo
      endif
!----------------------------------------------------------------
! read albedo data
!----------------------------------------------------------------
      do nn=1,2
      mm=mmse(nn)


      if (ialbflg.eq.0) write(ihdgi,31)ggdef,mm
      if (ialbflg.eq.1) write(ihdgi,37)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,work,istat)
      call unify_reducepick(nx,my,my_max,work,alvsfcl(1,1,nn))
!     call qmax2d(work,1,1,nx,my)


      if (ialbflg.eq.0) write(ihdgi,32)ggdef,mm
      if (ialbflg.eq.1) write(ihdgi,38)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,work,istat)
      call unify_reducepick(nx,my,my_max,work,alvwfcl(1,1,nn))
!     call qmax2d(work,1,1,nx,my)

      if (ialbflg.eq.0) write(ihdgi,33)ggdef,mm
      if (ialbflg.eq.1) write(ihdgi,39)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,work,istat)
      call unify_reducepick(nx,my,my_max,work,alnsfcl(1,1,nn))
!     call qmax2d(work,1,1,nx,my)

      if (ialbflg.eq.0) write(ihdgi,34)ggdef,mm
      if (ialbflg.eq.1) write(ihdgi,40)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,work,istat)
      call unify_reducepick(nx,my,my_max,work,alnwfcl(1,1,nn))
!     call qmax2d(work,1,1,nx,my)

      enddo ! end of nn
!
!-- facsf (0-100)
!
      write(ihdgi,35)ggdef
!     write(*,*)'35, lrec=',lrec
      call dmsread(nx,my,lncrec,'H',bckfile,work,istat)
      call unify_reducepick(nx,my,my_max,work,facsf)
!     call qmax2d(work,1,1,nx,my)
!
!-- facwf (0-100)
!
      write(ihdgi,36)ggdef
!     write(*,*)'35, lrec=',lrec
      call dmsread(nx,my,lncrec,'H',bckfile,work,istat)
      call unify_reducepick(nx,my,my_max,work,facwf)
!     call qmax2d(work,1,1,nx,my)
!
!----------------------------------------------------------------
#ifdef I38K
  31  format('S0003A','  GBCK',a4,4x,i2.2,6x)  ! alvsfcl
  32  format('S0003B','  GBCK',a4,4x,i2.2,6x)  ! alvwfcl
  33  format('S0003C','  GBCK',a4,4x,i2.2,6x)  ! alnsfcl
  34  format('S0003D','  GBCK',a4,4x,i2.2,6x)  ! alnwfcl
  35  format('S0003E','  GBCK',a4,12x)         ! facsf
  36  format('S0003F','  GBCK',a4,12x)         ! facwf
  37  format('S00X3A','  GBCK',a4,4x,i2.2,6x)  ! alvsfcl
  38  format('S00X3B','  GBCK',a4,4x,i2.2,6x)  ! alvwfcl
  39  format('S00X3C','  GBCK',a4,4x,i2.2,6x)  ! alnsfcl
  40  format('S00X3D','  GBCK',a4,4x,i2.2,6x)  ! alnwfcl
#else
  31  format('S0003A','GBCK',a4,4x,i2.2,6x)  ! alvsfcl
  32  format('S0003B','GBCK',a4,4x,i2.2,6x)  ! alvwfcl
  33  format('S0003C','GBCK',a4,4x,i2.2,6x)  ! alnsfcl
  34  format('S0003D','GBCK',a4,4x,i2.2,6x)  ! alnwfcl
  35  format('S0003E','GBCK',a4,12x)         ! facsf
  36  format('S0003F','GBCK',a4,12x)         ! facwf
  37  format('S00X3A','GBCK',a4,4x,i2.2,6x)  ! alvsfcl
  38  format('S00X3B','GBCK',a4,4x,i2.2,6x)  ! alvwfcl
  39  format('S00X3C','GBCK',a4,4x,i2.2,6x)  ! alnsfcl
  40  format('S00X3D','GBCK',a4,4x,i2.2,6x)  ! alnwfcl
#endif
!----------------------------------------------------------------

! 
!---  time interpolation
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
      do i=1,nxj
         alvsf(i,jj)=coef1*alvsfcl(i,jj,1)+coef2*alvsfcl(i,jj,2)
         alvwf(i,jj)=coef1*alvwfcl(i,jj,1)+coef2*alvwfcl(i,jj,2)
         alnsf(i,jj)=coef1*alnsfcl(i,jj,1)+coef2*alnsfcl(i,jj,2)
         alnwf(i,jj)=coef1*alnwfcl(i,jj,1)+coef2*alnwfcl(i,jj,2)
      enddo
      enddo

! 
!---  unit change
!
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          alvsf(i,jj)=alvsf(i,jj)*0.01
          alvwf(i,jj)=alvwf(i,jj)*0.01
          alnsf(i,jj)=alnsf(i,jj)*0.01
          alnwf(i,jj)=alnwf(i,jj)*0.01
          facsf(i,jj)=facsf(i,jj)*0.01
          facwf(i,jj)=facwf(i,jj)*0.01
        enddo
      enddo

!     if (myrank .eq. 0) then
!         print *,'*** for readalb.f ***'
!         print *,'julian :',julian
!         print *,'*** alvsf :'
!         call qmax2d(alvsf,1,1,nx,my)
!         print *,'*** alvwf :'
!         call qmax2d(alvwf,1,1,nx,my)
!         print *,'*** alnsf :'
!         call qmax2d(alnsf,1,1,nx,my)
!         print *,'*** alnwf :'
!         call qmax2d(alnwf,1,1,nx,my)
!         print *,'*** facsf :'
!         call qmax2d(facsf,1,1,nx,my)
!         print *,'*** facwf :'
!         call qmax2d(facwf,1,1,nx,my)
!     endif
!
      return
      end
