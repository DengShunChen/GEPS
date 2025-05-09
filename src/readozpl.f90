      subroutine read_ozplin(myrank)
!######################################################################
! subroutine to read the coefficients of ozone chemistry parametrization 
! The forecast model includes a prognostic equation for the ozone mass 
! mixing ration (kg/kg) 
!   dO3/dt = C1 + C2(O3-O3mean) + C3(T-Tmean) + C4(O3F - O3Fmean)
! The parameterization assumes that chemiscal changes in ozone can be 
! described by a linear relaxation towards a phtochemical equilibrium.
! Here Ci are relaxation rates which comes from NCEP-CFSV2.
!  
! Modifie by : Mei-Yu Chang. 2014/6/20
! #####################################################################
 
      use machine, only : kind_io4, kind_phys
      use ozne_def
      integer n,k,i,j,istat,myrank,k1,me

      real (kind=kind_io4), allocatable :: pl_lat4(:), pl_pres4(:),       &
                            pl_time4(:)
      real (kind=kind_io4), allocatable :: tempin(:)
      real (kind=kind_io4), allocatable :: tempin4(:,:,:,:)
      real (kind=kind_phys), allocatable :: xmin(:,:), xmax(:,:)
      real (kind=kind_phys), allocatable :: rmin(:), rmax(:)
      
      character*80 filename

!      data pathname/'FIXDIR'/
!      data pathname/'GFSWRK'/
      data filename/'global_o3prdlos'/
      me=1
!       
    open(kozpl,file=filename,form='unformatted',convert='BIG_ENDIAN')
!
        rewind (kozpl)
        read (kozpl) pl_coeff, latsozp, levozp, timeoz

        allocate (pl_lat(latsozp), pl_pres(levozp),pl_time(timeoz+1))
        allocate (pl_lat4(latsozp), pl_pres4(levozp),pl_time4(timeoz+1))

        rewind (kozpl)
        read (kozpl,end=110) pl_coeff, latsozp, levozp, timeoz, pl_lat4,   &
                      pl_pres4, pl_time4

        pl_lat(:)  = pl_lat4(:)
        pl_time(:) = pl_time4(:)

        do k = 1, levozp
           k1=levozp-k+1
           pl_pres(k1) = pl_pres4(k)
        enddo
!
      allocate(ozplin(latsozp,levozp,pl_coeff,timeoz))!OZONE P-L coeffcients
      allocate(tempin4(latsozp,levozp,pl_coeff,timeoz))!OZONE P-L coeffcients
      allocate(tempin(latsozp))

      do i=1,timeoz
        do n=1,pl_coeff
          do k=1,levozp
             read(kozpl,end=110) tempin
             tempin4(:,k,n,i) = tempin(:)
          enddo
        enddo
      enddo

      do i=1,timeoz
        do n=1,pl_coeff
          do k=1,levozp    ! k  = 1 (p_bottom)     to k = levozp (p_top) 
             k1=levozp-k+1 ! k1 = levozp(p_bottom) to k1= 1 (p_top)
          do j=1,latsozp
              ozplin(j,k1,n,i) = tempin4(j,k,n,i)
          enddo
        enddo
      enddo
      enddo

      if (myrank.eq.0 .and. me .eq. 0) then
          print *,'*** in read_ozpl ****'
          print *,' pl_coeff=',pl_coeff
          print *,' latsozp=',latsozp,' levozp=',levozp,' timeoz=',timeoz
          write(*,111)pl_lat 
          write(*,112)pl_pres
          write(*,113)pl_time
!
      
      allocate(xmax(pl_coeff,timeoz))
      allocate(xmin(pl_coeff,timeoz))
      allocate(rmax(pl_coeff))
      allocate(rmin(pl_coeff))

      do i=1,pl_coeff
         rmax(i)=-1.0e25
         rmin(i)=1.0e25
      enddo
      do i=1,timeoz
      do n=1,pl_coeff
         xmax(n,i)=-1.0e25
         xmin(n,i)=1.0e25
      enddo
      enddo


!     print out max/min
      
      amax=-1.0e25
      amin=1.0e25
      do i=1,timeoz
        do n=1,pl_coeff
          do k=1,levozp    ! k  = 1 (p_bottom)     to k = levozp (p_top) 
          do j=1,latsozp
              if (ozplin(j,k,n,i) .ge. xmax(n,i)) xmax(n,i)=ozplin(j,k,n,i)
              if (ozplin(j,k,n,i) .le. xmin(n,i)) xmin(n,i)=ozplin(j,k,n,i)
              if (ozplin(j,k,n,i) .ge. rmax(n)) rmax(n)=ozplin(j,k,n,i)
              if (ozplin(j,k,n,i) .le. rmin(n)) rmin(n)=ozplin(j,k,n,i)
              if (ozplin(j,k,n,i) .ge. amax) amax=ozplin(j,k,n,i)
              if (ozplin(j,k,n,i) .le. amin) amin=ozplin(j,k,n,i)
          enddo
        enddo
      enddo
      enddo

      print *,'*** in read_ozplin ****' 
      print *,'*** ozplin: max=',amax,' amin=',amin
      do n=1,pl_coeff
         write(*,121)n,rmax(n),rmin(n)
      enddo

      do i=1,timeoz
      do n=1,pl_coeff
         write(*,122)i,n,xmax(n,i),xmin(n,i)
      enddo
      enddo

      endif

 110  continue
 111  format (' pl_lat=',10f6.1)
 112  format (' pl_pres=',6f11.4)
 113  format (' pl_time=',10f6.1)
 121  format(' pl_coeff=',i2,' ozplin-max=',e15.4,                 &
             ' ozplin-min=',e15.4)
 122  format(' month=',i2,' pl_coeff=',i2,' ozplin-max=',e15.4,    &
             ' ozplin-min=',e15.4)
      return
      end
