      module radn

      implicit none

      public

!--------------------------------------------------------------------------   
! for rad_initialize 
!--------------------------------------------------------------------------   
      integer, save ::  levr,ictm,isol,ico2,iaer,ialb,iems,         &
               ntoz,iovr_sw,iovr_lw,isubc_sw,isubc_lw,              &
               icliq_sw,icice_sw,icliq_lw,icice_lw,                 &
               iflip,me,irad,ntcw,ioutsigr,ntiw,ntrw,ntsw,          &
               ntgl,ntinc,ntrnc,nthl
                 
      integer, save :: idate(8) 
      logical, save :: sashal,crick_proof,ccnorm,norad_precip
      real,    save :: sdec,cdec,slag,solcon
      real,    save :: solhr,rsolhr

      end module radn
