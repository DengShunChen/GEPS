module mod_split
      integer,parameter,private::NMAX=50
      real, allocatable:: fldsave(:,:,:)
      character:: ihdgo2save(NMAX)*28
      integer :: ptp1save(9,NMAX)
      integer :: oid
contains
      subroutine split_v3(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!
      use rank
      use const, only: RTYPE,kflag,ihdgo,ihdgo2,outgrb2 ,outdms
      use mod_grb2_param , only :wrt_grb2_v2 , wrt_grb2_accu_v2
!
      implicit none
!
      integer   nx,my,nc,istat,lenc,lev
      real(kind=RTYPE) glob(nx,my),mout(nx,my)
!      character*26 ihdg,ihdg2
!      character*80 ifilout
      integer::ptp0(9),ptp1(9)
      integer::itau_spl

      if (myrank .eq. nc) then
        mout   = glob
        ihdgo2 = ihdgo
        ptp1   = ptp0
      endif
        nc    = nc + 1
!
      if ( nc .eq. nsize .and. myrank .lt. nc ) then
        if(outdms.gt.0)then
          call dmswrit_split(nx,my,lenc,kflag,mout,istat)
        endif

        if(outgrb2 == 1 )then
#ifdef O38K
         read(ihdgo2(7:12),'(I6)')itau_spl
#else
         read(ihdgo2(7:10),'(I4)')itau_spl
#endif


!#ifdef USE_CUDA
!         oid = oid+1
!         if( oid > NMAX) stop
!         fldsave(:,:,oid) = mout(:,:)
!         ihdgo2save(oid)  = ihdgo2
!         ptp1save(:,oid)  = ptp1
!#else  
         if(ptp1(8)==-999)then
          call wrt_grb2_v2(itau_spl,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
              ,ptp1(6),ptp1(7),mout)
         else
          call wrt_grb2_accu_v2(itau_spl,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
              ,ptp1(6),ptp1(7),ptp1(8),ptp1(9),mout)
         endif
!#endif
        endif !outgrb2 == 1
        nc    = 0
      endif

      return
      end subroutine split_v3

      subroutine allocate_split()
        use param ,only:nx,my
        allocate( fldsave(nx,my,NMAX) )
        oid = 0
      end subroutine

      subroutine split_wrt
      use param ,only:nx,my
      use const, only: RTYPE,ihdgo2,outgrb2 ,ifilout_grb,idtg
      use mod_grb2_param , only :wrt_grb2_v2 , wrt_grb2_accu_v2,ofdir
      integer::ptp1(9)
      real(kind=RTYPE) mout(nx,my)
      integer::itau_spl


      if(outgrb2 == 1 )then
 134                    format( A  ,A ,I10.10 , i4.4       )
         
       do n = 1 , oid
         mout(:,:) = fldsave(:,:,n) 
         ihdgo2 = ihdgo2save(n)  
         ptp1 = ptp1save(:,n) 


#ifdef O38K
         read(ihdgo2(7:12),'(I6)')itau_spl
#else
         read(ihdgo2(7:10),'(I4)')itau_spl
#endif
         write(ofdir,134 )trim(ifilout_grb),'/',idtg/100 ,itau_spl
         if(ptp1(8)==-999)then
          call wrt_grb2_v2(itau_spl,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
              ,ptp1(6),ptp1(7),mout)
         else
          call wrt_grb2_accu_v2(itau_spl,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
              ,ptp1(6),ptp1(7),ptp1(8),ptp1(9),mout)
         endif
 
       enddo

      endif
      oid = 0

      end subroutine 

end module mod_split
