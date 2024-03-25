#ifdef RSM
      subroutine wrte_idate(idate)
!
      use rank, only : myrank
!
      integer idate,nsig
      character cidtg*10
!
      if (myrank.eq.0) then
        nsig=51
        write(cidtg,'(I10.10)') idate
      open(nsig,file='rsm_idate_'//cidtg,status='unknown', &
          form='unformatted',iostat=ios)
      write(nsig) idate
        print*,'idate in wrte_idate =',idate
      close(nsig)
      endif
!
      return
      end
#endif
