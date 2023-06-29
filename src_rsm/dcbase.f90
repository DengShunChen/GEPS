      subroutine dcbase ( nx,im,lm,hltm,kcbase,zl2,ql,qls2,hl,hls2 &
                        , game,e,qcl,hcl,cu,hbase,qbase,cubase )
!
!***********************************************************************
!
!    1.   function description
!           perform calculation of cloud base index and cloud base
!                   condition for a-s cumulus parameterization
!
!    2.   common block specification
!            none
!
!    3.   parameter specification
!     input :
!       im        :   dimension of horizontal direction
!       nx        :   no of points in horizontal direction to process
!       lm        :   position index in the latitudinal direction
!       hltm      :   latent heat constant of water vapor (2.52e+6)
!       zl2(im,lm):   height above surface (meters)
!       ql (im,lm):   water vapor mixing ratio on model level(g/g)
!       qls2(im,lm):   saturation water vapor mixing ratio on evenl leve
!       hl (im,lm):   moist static energy on model level
!       hls2(im,lm):  saturation moist static energy on model level
!       game(im,lm)
!     output:
!       kcbase(im):   cloud base index
!       e(im,lm)  :   normalized mass flux
!       qcl(im,lm):   water vapor mixing ratio of cloud (kg/kg)
!       hcl(im,lm):   moist static energy of cloud
!       cu (im,lm):   saturation water vapor mixing ratio of cloud
!
!    4.  calling modules
!            cup92
!
!    5.  usage
!           call dcbase ( nx,im,lm,hltm,kcbase,zl2,ql,qls2,hl,hls2
!          1            , game,e,qcl,hcl,cu,hbase,qbase,cubase )
!
!    6.  modules called
!            none
!
!    7.  date
!          created    1992           by      c-s chen  ( cwb  )
!          modify to f90 by C-H Lee and sort by River Chen in 2015
!
!***********************************************************************
!
      implicit  none

      integer   nx,im,lm
      real      hltm

!     input arrays
!
      real      zl2(im,lm),ql(im,lm),hl(im,lm),hls2(im,lm),game(im,lm)
      real      qls2(im,lm)
!
!     output arrays
!
      real      e(im,lm),qcl(im,lm),hcl(im,lm),cu(im,lm)
      real      hbase(im),qbase(im),cubase(im)
      integer   kcbase(im)
!

      integer   i,l
      real      tx

      do 10 i = 1, nx
      kcbase(i)= 1
      hbase(i) = 0.0
      qcl(i,lm)= ql(i,lm)
      hcl(i,lm)= hl(i,lm)
   10 continue
!
! --  assume cloud base can not be higher than lm/2
!
      do 30 l = lm-1, lm/2, -1
      do 40 i = 1,nx
      if (kcbase(i).eq.1) then
      tx= 1.0-zl2(i,l+1)/zl2(i,l)
      qcl(i,l)= qcl(i,l+1)+tx*(ql(i,l+1)-qcl(i,l+1))
      hcl(i,l)= hcl(i,l+1)+tx*(hl(i,l+1)-hcl(i,l+1))
      cu(i,l) = qcl(i,l)-(qls2(i,l)+game(i,l)/((1.+game(i,l))*hltm) &
                                   *(hcl(i,l)-hls2(i,l)))
!
      if ( cu(i,l).gt.0.0 ) then
      kcbase(i) = l
      hbase(i)  = hcl(i,l)
      qbase(i)  = qcl(i,l)
      cubase(i) = cu(i,l)
      endif
      endif
  40  continue
  30  continue
!
      do 60 l = lm-1, lm/2, -1
      do 60 i = 1, nx
      if (l.ge.kcbase(i))  then
      e(i,l) = zl2(i,l)/zl2(i,kcbase(i))
      cu(i,l+1) = 0.0
      endif
   60 continue
!
      return
      end
