#if defined(RSM) || defined(CWB_MPMD)
      subroutine send_idate(idate)
!
      use rank, only : root_rsm,myrank,itag
!
      integer idate
!
      if (myrank.eq.0) then

        itag=1
        call mpmd_send(idate,1,root_rsm,itag,'I')
        print*,'idate in send_idate =',idate
        print*,'send gfs idate itag=',itag
      endif
!
      return
      end
#endif
