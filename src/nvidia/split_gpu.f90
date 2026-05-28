      subroutine split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!
      use rank
      use const, only: RTYPE,kflag,ihdgo,ihdgo2,outgrb2 ,outdms
      use mod_grb2_param , only :wrtgrb2_v2_gpu , wrtgrb2_accu_v2_gpu
!
      implicit none
!
      integer   i,nx,my,nc,istat,lenc,lev
      real(kind=RTYPE) glob(nx,my),mout(nx,my)
!      character*26 ihdg,ihdg2
!      character*80 ifilout
      integer::ptp0(9),ptp1(9)
      integer::itau_spl

      if (myrank .eq. nc) then
        !$acc parallel loop
        do i=1,nx*my
          mout(i,1)   = glob(i,1)
        enddo
        ihdgo2 = ihdgo
        ptp1   = ptp0
      endif
        nc    = nc + 1
!
      if ( nc .eq. nsize .and. myrank .lt. nc ) then
        if(outdms.gt.0)then
          !$acc update self(mout)
          call dmswrit_split(nx,my,lenc,kflag,mout,istat)
        endif
        if(outgrb2 == 1 )then
#ifdef O38K
         read(ihdgo2(7:12),'(I6)')itau_spl
#else
         read(ihdgo2(7:10),'(I4)')itau_spl
#endif
         if(ptp1(8)==-999)then
          call wrtgrb2_v2_gpu(itau_spl,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
              ,ptp1(6),ptp1(7),mout)
         else
          call wrtgrb2_accu_v2_gpu(itau_spl,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
              ,ptp1(6),ptp1(7),ptp1(8),ptp1(9),mout)
         endif
        endif
        nc    = 0
      endif
!
      return
      end subroutine split_v2_gpu

