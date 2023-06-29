       subroutine prerrtmg(                                       &     
!    for input :                                                  &
                nx,my,my_max,idtg,tau,dt,hours,frad,uprad,        &
                isubc_sw,isubc_lw,d2r,xlon,myrank,me,             &   
!    for output  
                idat,jdat,solhr,dtsw,dtlw,lsswr,lslwr,            &
                slag,sdec,cdec,solcon,                            &
                xlonr,ixseed)

!--------------------------------------------------------------------
!   Prepare for RRTMG scheme 
!--------------------------------------------------------------------
      use module_radiation_driver, only :radupdate

      use mersenne_twister, only : random_setseed, random_index, &        
                                   random_stat
      use index

      implicit none


      integer i,j,jj,nxj,nx,my,my_max,me,k
      integer numrdm(nx*my*2),ixseed(nx,my,2),ipsdlim
      type (random_stat)::stat
!
      integer idat(8),jdat(8)
      integer*8 idtg,idtg2
      character cdtg*12
      integer itau,isubc_sw,isubc_lw,myrank,ii,ipsd0,ipseed
      logical uprad,lsswr,lslwr
      real tau,slag,sdec,cdec,solcon,frad,dt,d2r
      real solhr,dtsw,dtlw,hours
      real xlon(nx,my_max),xlonr(nxp,my_max)
      parameter(ipsdlim=100000000)!upper limit for random seed

! CWB2016 
      ipsd0=0
      numrdm=0

      do i=1,8
         idat(i)=0
         jdat(i)=0
      enddO
!
! --- get initial time : idat
!
      write(cdtg,900)idtg
 900  format(i12.12)
      read(cdtg,'(i4,i2,i2,i2,i2)') idat(1),idat(2),idat(3),idat(5) &
                                   ,ii
!
! --- get forecast time : jdat
!
      itau=int(tau+0.001) ! tau from intgrt.f
!      if (myrank .eq. 0) print *,'+itau=',itau,' tau=',tau
!!      call dtgfix12(idtg,idtg2,itau)
!!      write(cdtg,900)idtg2
!!      read(cdtg,'(i4,i2,i2,i2,i2)') jdat(1),jdat(2),jdat(3),jdat(5) &
!!                                    ,ii
      call dtgfix12_new(idtg,jdat(1),jdat(2),jdat(3),jdat(5),ii,itau)
!
! --- set solhr = forecast hours if not at initial time
!
!        if (hours .gt. 0.0) then
           solhr = hours
!        else
!           solhr = idat(5)  !initial time
!        endif

!        if (myrank .eq. 0) print *,'solhr=',solhr
!
! --- set dtsw, dtlw, lsswr, lslwr
!
      dtsw  = 3600.0 * frad ! sw freq =3600s (frad=1 hr)
      dtlw  = 3600.0 * frad
      lsswr = uprad
      lslwr = uprad
!
!---------------------------------------------------------------------
!     uprad = T : run radiation scheme    
!---------------------------------------------------------------------
!      if (myrank .eq. 0) print *,'*** for uprad = T'
!
! --- run radupdate
!
      if ( uprad ) then
      call radupdate                                                 &
! --- inputs:                                                        &
        ( idat, jdat, dtsw, dt, lsswr, me, myrank,                   &
! --- outputs:                                                       &
         slag, sdec, cdec, solcon)
!
!      if (myrank .eq. 0) then 
!          print *,'cdtg=',cdtg
!          print *,'idtg=',idtg
!          print *,'idtg2=',idtg2
!          print *,'idat(1)=',idat(1)
!          print *,'idat(2)=',idat(2)
!          print *,'idat(3)=',idat(3)
!          print *,'idat(5)=',idat(5)
!          print *,'jdat(1)=',jdat(1)
!          print *,'jdat(2)=',jdat(2)
!          print *,'jdat(3)=',jdat(3)
!          print *,'jdat(5)=',jdat(5)
!      endif

!       if (myrank .eq.0 )print *,'*** radupdate ok ! ***'
!       if (myrank .eq.0 )print *,'slag=',slag,' solcon=',solcon
!       if (myrank .eq.0 )print *,'sdec=',sdec,' cdec=',cdec
!       if (myrank .eq.0 )print *,'uprad=',uprad
!       if (myrank .eq.0 )print *,'solhr=',solhr

!---------------------------------------------------------------------
!      generate initial permutation seed for random number generator
!---------------------------------------------------------------------
! --- set isubc_sw/isubc_lw for radn.h and block.f
!
      if ( (isubc_lw .eq. 2) .or. (isubc_sw .eq. 2) ) then
          ipsd0 = 17*idat(5) + 43*idat(2) + 37*idat(3) + 23*idat(1)
      endif

! --- set initial seed ipsd0 by initial time
!      if (myrank .eq. 0) then
!           print *,'  Radiation sub-cloud initial seed =',ipsd0       
!      endif


! --- generate 2-d random seeds array for sub-grid cloud-radiation
!
      if ( (isubc_lw .eq. 2) .or. (isubc_sw .eq. 2) ) then
      ipseed = mod(nint(100.0*sqrt(tau*3600)), ipsdlim) + 1 + ipsd0

      call random_setseed                                            &
!  ---  inputs:                                                      &
           ( ipseed,                                                 &
!  ---  outputs:
             stat)                                                   

        call random_index                                            &
!  ---  inputs:                                                      &
           ( ipsdlim,                                                &
!  ---  outputs:
             numrdm, stat)                                       

        do k = 1, 2
          do j = 1, my
            nxj=nxdef(j)
            do i = 1, nxj
!            do i = 1, nx
              ixseed(i,j,k) = numrdm(i+(j-1)*nx+(k-1)*my*nx)
!              ixseed(i,j,k) = numrdm(i+(j-1)*nxj+(k-1)*my*nxj)
            enddo
          enddo
        enddo

      endif ! for isub_lw=2 .or isub_sw=2
      endif ! if ( uprad )

!      if (myrank .eq. 0) print *,'random_index ok!'
     
!--------------------------------------------------------------------
!    setting xlonr for grrad : input 
!--------------------------------------------------------------------
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       ii=nxjstart(j)
      do i = 1, nxj
         xlonr(i,jj)=xlon(ii,jj)*d2r
         if (xlon(ii,jj).lt. 0) xlonr(i,jj)=(xlon(ii,jj)+360.)*d2r
         ii=ii+1
      enddo
      enddo
!      if (myrank .eq. 0) print *,'for xlonr setting ok!'
!--------------------------------------------------------------------
      return
      end

