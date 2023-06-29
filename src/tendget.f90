      subroutine tendget (phiten)

!CWB2021 ndsl single precision test

!
!  purpose : compute tendency of vorticity,divergence and geopotential
!            using subroutine of forecast model
!-----------------------------------------------------------------------
!  **** input *****
!  the input are pass bye those '*.h' files (see below)
!
!  **** output *****
!  phiten   : geopotential tendency array
!  other output are pass by 'spec.h'
!-----------------------------------------------------------------------
      use param
      use mpe
      use rank
      use index
      use const
      use spec
      use grid

      implicit none

      integer  m,mf,j,jj,k,n,l,mlst,nxj,nk,kk,kl,i,ii,nl
      real     sqhaf

      real(kind=RTYPE) phiten(levp,2,jtrun,jtmax)
!byl      real tbar(lev),qbar(lev*ncld)
      real(kind=RTYPE) sdpbl(nxp,my_max)

! for Semi-Lagrangian
      real(kind=RTYPE) deldm(nxp,my_max)

!CWB2021 ndsl single precision test
      real(kind=RTYPE) uum_sl(nx,levp,my_max)              &
      ,vvm_sl(nx,levp,my_max)                              &
      ,ttm_sl(nx,levp,my_max)                              &
      ,pten_sl(nx,levp,my_max)                             &
      ,qm_sl(nx,levp*ncld,my_max)                          &
      ,pdot(nxp,lev+1,latpart)                             &
      ,vdmerdr(nxp,lev,my_max),vdzonlr(nxp,lev,my_max)     &
      ,vdmerd(nxp,lev,my_max),vdzonl(nxp,lev,my_max)       &
      ,ddtemp(nxp,lev,my_max),qvadv(nxp,lev*ncld,my_max)   &
      ,diveng(nxp,lev,my_max),pten(nxp,lev,my_max)

!byl      real cc(nx+2,levp,1+ncld,my_max)
      real(kind=RTYPE) cc(nx+2,levp,1,my_max),dummy
!byl,wss(levp,2,1+ncld,jtrun,jtmax)
!
      integer   ierr
      real(kind=RTYPE) dta
      real(kind=RTYPE) ww1(nx,my_max)

!for 2dMPI
      real(kind=RTYPE) temten1(lev,2,jtrun,jtmax)
      real(kind=RTYPE) phiten1(lev,2,jtrun,jtmax)
      logical forward
      forward = .false.
      phiten1=0.

!
!  global mean tempertures (tbar)
!
      sqhaf= sqrt(0.5)
!byl      tbar=0. ; qbar=0.
!     do 160 m = 1, mlistnum
!      mf=mlist(m)
!      if ( mf.eq.1) then
!       do 161 k=1,lev
!         tbar(k)= sqhaf*temnow(k,1,1,m)
! 161   continue
!       do 162 k=1,lev*ncld
!         qbar(k)= sqhaf*qnow(k,1,1,m)
! 162   continue
!      endif
! 160 continue
!

!2dMPI
!byl      do 160 m = 1, mlistnum
!byl       mf=mlist(m)
!byl       if ( mf.eq.1) then
!byl        do 161 k=1,levp
!byl          KK=Llist(k)
!byl          tbar(KK)= sqhaf*temnow(k,1,1,m)
!byl  161   continue
!byl        do 162 n=1,ncld
!byl           nk=(n-1)*levp
!byl           nL=(n-1)*lev 
!byl        do 162 k=1,levp
!byl           KK=nk+k
!byl           KL=nL+Llist(k)
!byl          qbar(KL)= sqhaf*qnow(KK,1,1,m)
!byl  162   continue
!byl       endif
!byl  160 continue

!byl      call mpe_global_sum(tbar,lev,mpe_double)
!byl      call mpe_global_sum(qbar,lev*ncld,mpe_double)
!
      do 170 m = 1, mlistnum
         mf=mlist(m)
      do 170 n = mf, jtrun
      plten(n,m,1) = 0.0
      plten(n,m,2) = 0.0
  170 continue
!
      do 180 m = 1, mlistnum
         mf=mlist(m)
      do 180 n = mf, jtrun
      do 180 k = 1, levp
      divten(k,1,n,m) = 0.0
      divten(k,2,n,m) = 0.0
      vorten(k,1,n,m) = 0.0
      vorten(k,2,n,m) = 0.0
      temten(k,1,n,m) = 0.0
      temten(k,2,n,m) = 0.0
!ncld qten(k,1,n,m) = 0.0
      hldten(k,1,n,m) = 0.0
      hldten(k,2,n,m) = 0.0
      phiten(k,1,n,m) = 0.0
      phiten(k,2,n,m) = 0.0
  180 continue
!
!!      do 183 m = 1, mlistnum
!!         mf=mlist(m)
!!      do 183 n = mf, jtrun
!!      do 183 k = 1, levp*ncld*2
!!      qten(k,1,n,m) = 0.0
!!  183 continue
!
       pdot=0.
       vdmerd=0.
       vdzonl=0.
       qvadv=0.
       ddtemp=0.
       pten=0.
!
       dta=0.5*dt
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do k = 1, lev
          do i = 1, nxj
            up(i,k,jj) = ut(i,k,jj)
            vp(i,k,jj) = vt(i,k,jj)
            ttp(i,k,jj)= tt(i,k,jj)
          enddo
        enddo
        do k = 1, lev*ncld
          do i = 1, nxj
            qm(i,k,jj) = qt(i,k,jj)
          enddo
        enddo
!!        do i = 1, nxj
!!          ptp(i,jj)= pt(i,jj)
!!        enddo
      enddo
!ch>
! transpose partial to full: ut -> ut_sl, vt -> vt_sl, ttm -> ttm_sl, qp -> qm_sl

!#ifdef MULTIPLE
!      call mpe2d_transpose_ndsl_p2f_multi(ut   ,vt   ,ttp   ,dummy,dummy,qp,    &
!                                    ut_sl,vt_sl,ttm_sl,dummy,dummy,qm_sl, &
!                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm,4)
!#else
      call mpe2d_transpose_ndsl_p2f(ut,ut_sl,    &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(vt,vt_sl,    &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(ttp,ttm_sl,  &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(qm,qm_sl,    &
                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!#endif

!!    call mpe2d_unify_nx(pt_sl,pt)
      uum_sl=ut_sl
      vvm_sl=vt_sl
!!    ptp_sl=pt_sl
!ch<

!
!     Semi-Lagrangian
!       Horizontal Advection

        call ndslfv_monoadvh2(ttm_sl,qm_sl,pten_sl,uum_sl,vvm_sl,nxdef &
                             ,dta,levp)

!ch>
! transpose full to partial: ttm_sl -> ddtemp,  pten_sl -> pten, uum_sl -> vdzonl
!                            vvm_sl -> vdmerd,  qm_sl -> qvadv

!#ifdef MULTIPLE
!      call mpe2d_transpose_ndsl_f2p_multi(ttm_sl,pten_sl,uum_sl,vvm_sl,qm_sl, &
!                                    ddtemp   ,pten   ,vdzonl   ,vdmerd   ,qvadv   , &
!                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!#else
      call mpe2d_transpose_ndsl_f2p(ttm_sl,ddtemp, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(pten_sl,pten, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(uum_sl,vdzonl, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(vvm_sl,vdmerd, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(qm_sl,qvadv,   &
                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
!#endif

!ch<

      do jj =1, jlistnum
      j=jlist1(jj)
!     nxj=nxdef(j)
!
!  p**capa quantities
!
       call prexp_hybrid_cwb (nxjp(j),nxp,lev,ptop,sigma,pt(1,jj),   &
                          pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj))
!
! Calculate Vertical velocity & Stream Functions
!
        call gridnl_hybrid_ndsl (nxjp(j),nxp,lev,ncld                  &
        , cp,radsq,ut(1,1,jj),vt(1,1,jj),rdiv(1,1,jj),tt(1,1,jj)       &
        , qt(1,1,jj),phi(1,1,jj),pt(1,jj),dtpl(1,jj),dlpl(1,jj),sinl(j)&
        , pk(1,1,jj),pk2(1,1,jj),dsigma,sigma,onocos(j),cor(j)         &
        , diveng(1,1,jj),vdmerdr(1,1,jj),vdzonlr(1,1,jj),pten(1,1,jj)  &
        , deldm(1,jj),sdpbl(1,jj),sd(1,1,jj),pdot(1,1,jj),vvel(1,1,jj) &
        , sgeo(1,jj) )

      enddo !jj = 1,jlistnum
!
!CWB2021 for single precision test
        call joinrs(cc,diveng,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
        call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc       &
                   ,hldten,1,nsizey)

        call trngra3(jtrun,jtmax,nx,levp,my,my_max,cim,poly,dpoly      &
                   ,hldten,dlphi,dtphi,nsizey)
!
!       update all horizontal informations
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            do k=1,lev
              vdmerdr(i,k,jj)=vdmerdr(i,k,jj)-dtphi(i,k,jj)/radsq/onocos(j)! &
!                             -ut(i,k,jj)*cor(j)-(ut(i,k,jj)*ut(i,k,jj)      &
!                             +vt(i,k,jj)*vt(i,k,jj))*onocos(j)*sinl(j) 

              vdzonlr(i,k,jj)=vdzonlr(i,k,jj)-dlphi(i,k,jj)/radsq!           &
!                             +vt(i,k,jj)*cor(j)
            enddo
          enddo
        enddo !jj = 1,jlistnum
      call ndslfv_update(nxjp,vdzonl,vdmerd,vdzonlr,vdmerdr,dta,forward)

!CWB2021 ndsl single precision test
!
!
!       Vertical Advection
!
      call ndslfv_monoadvv(ddtemp,qvadv,vdzonl,vdmerd,pdot,pt,nxjp,dta,forward)

!CWB2021 ndsl single precision test

!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do k = 1, lev
          do i = 1, nxj
            vdzonl(i,k,jj) = (vdzonl(i,k,jj)- up(i,k,jj))/dt
            vdmerd(i,k,jj) = (vp(i,k,jj)- vdmerd(i,k,jj))/dt
            ddtemp(i,k,jj)=  (ddtemp(i,k,jj)-ttp(i,k,jj))/dt
          enddo
        enddo
      enddo
!
!  combine non-linear grid point terms via gaussian quadrature
!
      call joinrs(cc,ddtemp,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc  &
                   ,temten,1,nsizey)
!!      call joinrs(cc,ddtemp,qvadv,dummy,dummy,nx,my_max,lev           &
!!                 ,jlistnum,2,ncld)
!!      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc         &
!!                 ,wss,1+ncld,nsize)
!!      call ujoinrs(wss,temten,qten,dummy,dummy,jtrun,jtmax,lev        &
!!                 ,mlistnum,2,ncld)
      call mpe2d_unify_nx(ww1,deldm)
      call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,ww1           &
                 ,plten,nsizey)
!
!CWB2021 for single precision test
      call rstrandz (jtrun,jtmax,nx,my,my_max,levp,vdmerd,vdzonl      &
                    ,weight,cim,onocos,poly,dpoly,divten,vorten,nsizey)
!
!  no tendency of zero mode
!
      do 96 m =1,mlistnum
        mf=mlist(m)
        if ( mf.eq.1) then
        do k = 1, levp
         vorten(k,1,1,m)= 0.
         vorten(k,2,1,m)= 0.
         divten(k,1,1,m)= 0.
         divten(k,2,1,m)= 0.
         temten(k,1,1,m)= 0.
         temten(k,2,1,m)= 0.
        enddo
        endif
 96   continue
!
!  compute phiten from temden and plten
!
         mlst=ilist(1)
         if(mlst .ne. 0) then
           plten(1,mlst,1)= 0.
           plten(1,mlst,2)= 0.
         endif
!
!       do 100 m =1,mlistnum
!        mf=mlist(m)
!       do 100 n  = mf,jtrun
!       do 99 k = 1, lev
!         phiten(k,1,n,m) = spalm(k)*plten(n,m,1)
!         phiten(k,2,n,m) = spalm(k)*plten(n,m,2)
!       do 98 l = 1, lev
!         phiten(k,1,n,m) = phiten(k,1,n,m)+arrhyd(k,l)*temten(l,1,n,m)
!         phiten(k,2,n,m) = phiten(k,2,n,m)+arrhyd(k,l)*temten(l,2,n,m)
!98     continue
!99     continue
!100    continue
!
!ch>
        call mpe2d_unify_lev(temten,temten1,lev,levp,jtrun,jtmax,mlistnum,nsizex,row_comm)

        do 100 m =1,mlistnum
         mf=mlist(m)
        do 100 n  = mf,jtrun
        do 99 k = 1, lev
          phiten1(k,1,n,m) = spalm(k)*plten(n,m,1)
          phiten1(k,2,n,m) = spalm(k)*plten(n,m,2)
        do 98 L = 1, lev
          phiten1(k,1,n,m) = phiten1(k,1,n,m)+arrhyd(k,L)*temten1(L,1,n,m)
          phiten1(k,2,n,m) = phiten1(k,2,n,m)+arrhyd(k,L)*temten1(L,2,n,m)
 98     continue
 99     continue
 100    continue

        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun
            do k = 1, levp
               kk=Llist(k)
               phiten(k,1,n,m) = phiten1(kk,1,n,m)
               phiten(k,2,n,m) = phiten1(kk,2,n,m)
            enddo
          enddo
        enddo
!ch<

      return
      end
