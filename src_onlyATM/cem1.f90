      subroutine cem1 ( nn,nx,lev,lt,cp,r,g,pcon,ccritl,etop,kcbase   &
                      , cadj,htop,hbase, hswo,ho,qo,qswo,hswe,qswe    &
                      , gamqo,gamqe,dzlc,dzl2c,e,exist,hc,qc,cu,rain  &
                      , zl,zl2 )
!
!***********************************************************************
!
!    1.   function description
!           cloud ensemble model for a-s scheme
!           normalized cloud mass flux is linear function of height
!
!    2.   common block specification
!            none
!
!    3.   parameter specification
!     input :
!       nn        :   no. of potential convective points in x-direction
!       nx        :   dimension of horizontal direction
!       lev       :   position index in the latitudinal direction
!       cp        :   heat specific constant of the air (=1004.)
!       r         :   gas constant (287.0)
!       g         :   gravity (9.806)
!       cadj      :   logical flag for water budget calculations
!       kcbase(nx):   cloud base index
!       htop(nx)  :   moist static energy on cloud top
!       hbase(nx) :   moist static energy on cloud base
!       ho (nx,lev):   moist static energy on model level
!       hswo(nx,lev):  saturation moist static energy on model level
!       hswe(nx,lev):  saturation moist static energy on even level
!       qo (nx,lev):   temperature on model level
!       qswo(nx,lev):  saturation mixing ratio on model level
!       qswe(nx,lev):  saturation mixing ratio on even level
!       gamqo(nx,lev)
!       gamqe(nx,lev)
!       dzlc(nx,lev))
!       dzl2c(nx,lev)
!     output:
!       e(nx,lev) : mass flux normalized with cloud base mass flux
!       hc(nx,lev): moist static energy for cloud
!       qc(nx,lev): mixing ratio for cloud
!       cu(nx,lev): cloud liquid water
!       rain(nx,lev): rain water
!     working array
!
!    4.  calling modules
!            cup92
!
!    5.  usage
!         call cem1 (nn,nx,lev,lt,cp,r,g,pcon,ccritl,etop,kcbase,cadj
!         1         ,htop,hbase,hswo,ho,qo,qswo,hswe,qswe,gamqo,gamqe
!         2         , dzlc,dzl2c,e,exist,hc,qc,cu,rain )
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
      use paramt

      implicit  none

      integer   nn,nx,lev,lt,i,l,lk,kklb

!     input array
      integer   kcbase(nx)
      real      hbase(nx),htop(nx,lev),ho(nx,lev),                  &
                dzlc(nx,lev),dzl2c(nx,lev),qo(nx,lev),qswo(nx,lev), &
                qswe(nx,lev),hswe(nx,lev),hswo(nx,lev),             &
                gamqo(nx,lev),gamqe(nx,lev)
      logical   exist(nx),cadj

!     output array
      real      e(nx,lev),hc(nx,lev),qc(nx,lev),cu(nx,lev),rain(nx,lev)
!
!  local work arrays
!
      real      alamda(im),addtop(im)
      real      alammax(im,lm),zl(im,lm),zl2(im,lm)
      real      etop,ccritl,pcon,g,r,cp,qctem,totliq
!
! subroutine cem1  calculates the thermaodynamic profile of
!  an updraft sub-ensemble according to a given lamda
!
! all the properties of updrafts are definded at even levels except
!   at the cloud tops where they are defined at odd levels
! mixing and its related processes take place at odd levels
! level indices for all properties of cumulus clouds are added by 1
!   the cloud base level is identified with level index (lb+1)
!   the cloud top levels are identified with level indices lt's
! ice phase is not considered in this version
!
!
!     constant
!
      do i=1,nn
        kklb=kcbase(i)
        alammax(i,lt) = 5.0/(zl(i,lt) - zl2(i,kklb))
      enddo
!
      do 50 i = 1,nn
      if (exist(i)) then
      addtop(i)= etop*(htop(i,lt)-ho(i,lt))
      alamda(i)= dzlc(i,lt)*(htop(i,lt)-ho(i,lt)+addtop(i))
      endif
   50 continue
!
      if(lt.lt.lev-1) then
      do 60 l=lt+1,lev-1
      do 60 i=1,nn
      if (exist(i).and.l.le.kcbase(i)) then
      alamda(i) = alamda(i) + dzl2c(i,l)*(htop(i,lt)-ho(i,l) &
                                        +addtop(i))
      endif
  60  continue
      endif
!
      do 70 i=1,nn
      if (exist(i)) then
      alamda(i) = (hbase(i) - htop(i,lt))/alamda(i)
      exist(i)= exist(i).and.alamda(i).ge.0.0
      exist(i)= exist(i).and.alamda(i).lt.alammax(i,lt)
      addtop(i)= dzlc(i,lt)
      endif
  70  continue
!
      do 110 lk= lev,lt+1,-1
      do 110 i = 1,nn
      if (exist(i).and.lk.le.kcbase(i))  then
      e(i,lk-1) = e(i,lk)+alamda(i)*dzl2c(i,lk)
      hc(i,lk-1)=(e(i,lk)*hc(i,lk)+alamda(i)*ho(i,lk)*dzl2c(i,lk)) &
                / e(i,lk-1)
      addtop(i)= addtop(i)+dzl2c(i,lk)
      else if (lk.le.kcbase(i)) then
      e(i,lk-1) = 0.0
      hc(i,lk-1)= 0.0
      endif
 110  continue
!
      do 120 i=1,nn
      if (exist(i))  then
      e(i,lt-1) = e(i,lt) +alamda(i)*(dzlc(i,lt)+etop*addtop(i))
      hc(i,lt-1)=(e(i,lt)*hc(i,lt)+(e(i,lt-1)-e(i,lt))*ho(i,lt)) &
                / e(i,lt-1)
      else
      e(i,lt-1) = 0.0
      hc(i,lt-1)= 0.0
      endif
  120 continue
!
!  return if water budgets not need (cem1 calls with adjusted soundings)
!
      if (cadj)  return
!
      do 100 lk= lev,lt+1,-1
      do 100 i=1,nn
      if (exist(i).and.lk.le.kcbase(i))  then
      qctem= qswe(i,lk-1)+(hc(i,lk-1)-hswe(i,lk-1))*gamqe(i,lk-1)
!
      qc(i,lk-1)=(e(i,lk)*qc(i,lk)+alamda(i)*qo(i,lk) &
               *dzl2c(i,lk))/e(i,lk-1)
      totliq= qc(i,lk-1)-qctem
!
      if (totliq .ge. ccritl)  then
!
      cu(i,lk-1) = (totliq+ccritl*dzl2c(i,lk)) &
                  /(1.+pcon*dzl2c(i,lk))
      rain(i,lk)=totliq-cu(i,lk-1)
      qc(i,lk-1)= qctem+cu(i,lk-1)
!
      else
!
        rain(i,lk) = 0.0
        cu(i,lk-1)= 0.0
      endif
!
      endif
 100  continue
!
!-----cloud top
!
      do 500 i=1,nn
      if (exist(i))  then
      qctem= qswo(i,lt)+(hc(i,lt-1)-hswo(i,lt))*gamqo(i,lt)
!
      qc(i,lt-1) =(e(i,lt)*qc(i,lt)+(e(i,lt-1)-e(i,lt))*qo(i,lt)) &
                /e(i,lt-1)
      totliq= qc(i,lt-1)-qctem
!
      if (totliq .ge. ccritl)  then
!
      cu(i,lt-1) = (totliq+ccritl*dzl2c(i,lt)) &
                  /(1.+pcon*dzlc(i,lt))
      rain(i,lt)=totliq-cu(i,lt-1)
      qc(i,lt-1)= qctem+cu(i,lt-1)
!
      else
!
        rain(i,lt) = 0.0
        cu(i,lt-1)= 0.0
      endif
!
      endif
 500  continue
      return
      end
