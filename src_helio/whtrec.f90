      subroutine whtrec (label,ntau,taudir,whtlev,num)
!
      implicit none

      integer ntau,num,n,lev

      character*6 label,labx
      character*16 taudir(ntau)
!
      real  whtlev(100)
!
      num= 0
      do 10 n=1,ntau
      read(taudir(n),'(a6,1x,i4)') labx,lev
      if(label.eq.labx) then
      num= num+1
      whtlev(num)= lev
      endif
   10 continue
!
      return
      end
