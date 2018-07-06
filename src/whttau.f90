      subroutine whttau (itau,numout,outdir,ntau,taudir)
!
      implicit none

      integer  itau,numout,ntau,ii,ktau

      character*16 outdir(numout),taudir(numout)
!
      ntau= 0
      do 10 ii=1,numout
      if(outdir(ii)(1:8).eq.'nomodata')return
      read(outdir(ii),'(12x,i4)') ktau
      if(ktau.eq.itau) then
      ntau= ntau+1
      taudir(ntau)= outdir(ii)
      endif
   10 continue
!
      return
      end
