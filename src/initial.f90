      subroutine initial(no,jtrun,jtmax,lev,nx,my,my_max,mlmax)

!
!  purpose : do nonlinear normal mode initialization
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!-----------------------------------------------------------
      use param, only : ncld
      use mpe
      use rank
      use index
      use const, only : eigval,omega,rad,nnmiit,doincr,evecin,nnmivm,   &
                        cutfreq,evectr,pmcor,tmcor,poly,dpoly,cim,      &
                        wdfac,wcfac,onocos,RTYPE
      use spec,  only : temold,vorten,vorold,divten,divold,plnow,temnow,&
                        divnow,vornow,plold
!byl                        divnow,vornow,qold,qnow,plold
      use grid,  only : pt,tt,rdiv,rvor,dtpl,dlpl,ut,vt
      use fftcom

      implicit none
      integer  no,jtrun,jtmax,lev,nx,my,my_max,mlmax
!
!  working array
!
!     integer, parameter ::  no=2*((jtrun+1)/2)+(jtrun/2)+10
      real      eval(no),evec(no*no),epos(no*no),x(no*2)
!byl      real      a(jtrun,jtrun,lev),b(jtrun,jtrun,lev)
      real      a(jtrun,jtrun,nnmivm),b(jtrun,jtrun,nnmivm)
      real      mx(no*no),h(jtrun,jtmax,nnmivm),c(jtrun,jtrun,nnmivm)
      integer   nw(jtrun,jtmax)
      real(kind=RTYPE) phiten(levp,2,jtrun,jtmax),dummy
      real      wk(no*no),wc(no*2),wd(no*2),ew(no*no)
      real(kind=RTYPE) cc(nx+2,levp,1,my_max)
!byl      real      cc(nx+2,levp,3,my_max),wss(levp,2,3,jtrun,jtmax)
      real, allocatable :: bal_tmp(:)
      character lab*10,lrec*16

      integer   mlmax2,j,k,l,ic,m,mf,n,ns,na,nbig,kk,kL
      integer   jtrun_l,jtmax_l,lev_l,nx_l,my_l,my_max_l,no_l
      real      bal
      logical   nnmical
      integer   brank   !root rank of row_broadcast
      ! #region agent log
      integer(kind=8) :: dbg0, dbg1, dbg2
      ! #region agent log: CPU baseline for the GPU probes (handoff layer 8)
      real(kind=8) :: dbgx_in, dbgx_out, dbgev
      real(kind=8) :: dbgx2_tot
      real(kind=8) :: dbgcv, dbgcd, dbgcp
      real(kind=8) :: dbgipt, dbgiut, dbgivt, dbgitt, dbgird, dbgirv, dbgidt, dbgidl
      integer(kind=8) :: dbgnd_tot
      real(kind=8) :: dbgbal_prev, dbgbal_d, dbgev_s, dbgx2
      ! #endregion
      interface
        subroutine geps_dbg_log(hyp, locid, irank, p0, p1, p2)
          integer hyp, locid, irank
          integer(kind=8) p0, p1, p2
        end subroutine
      end interface
      ! #endregion
!
      phiten=0.
      mlmax2 = mlmax*2
      jtrun_l = jtrun
      jtmax_l = jtmax
      lev_l = lev
      nx_l = nx
      my_l = my
      my_max_l = my_max
      no_l = no
      allocate(bal_tmp(jtrun_l))
      bal_tmp = 0.
      if(myrank .eq. 0) print *,'jrtun=',jtrun_l,' mlmax=',mlmax
      lab='pt'
      call  check(pt,nx,my,my_max,lab)
!
!  define constants and comput coefficients for initializatin
!
      call  inicons (rad,omega,eigval,lev_l,jtrun_l,jtmax_l, &
                     nw,a,b,c,h,nnmivm)
!
!  begin to iterration, now doing 3 iterrations
!
      do 120 ic=1,nnmiit
!
        if(myrank .eq. 0) print *,'iteration=',ic
!
!  get tendency of vorticity,divergence,geopotential
!
          ! #region agent log
          ! tendget's INPUTS, same set and same point as the GPU probes
          ! 415/416/417. Decides whether tendget produces the bad phiten
          ! or merely inherits it. Same caveat as DBGCHAIN 412/413:
          ! rank 0's slice only, so compare GROWTH and run-to-run spread,
          ! not absolute values against the single-rank GPU run.
          dbgipt = sum(real(pt,   kind=8)**2)
          dbgiut = sum(real(ut,   kind=8)**2)
          dbgivt = sum(real(vt,   kind=8)**2)
          dbgitt = sum(real(tt,   kind=8)**2)
          dbgird = sum(real(rdiv, kind=8)**2)
          dbgirv = sum(real(rvor, kind=8)**2)
          dbgidt = sum(real(dtpl, kind=8)**2)
          dbgidl = sum(real(dlpl, kind=8)**2)
          if (myrank .eq. 0) print *,'DBGIN 415 spt=',dbgipt,' sut=',dbgiut,' svt=',dbgivt
          if (myrank .eq. 0) print *,'DBGIN 416 stt=',dbgitt,' srd=',dbgird,' srv=',dbgirv
          if (myrank .eq. 0) print *,'DBGIN 417 sdt=',dbgidt,' sdl=',dbgidl
          ! #endregion
        call tendget (phiten)
          ! #region agent log
          ! Sigma(.)^2 at the same two points as the GPU probes 412/413,
          ! to bisect where x's magnitude blows up. ABSOLUTE values are
          ! NOT comparable with the GPU (this is rank 0's slice, jtmax_l,
          ! while the GPU run is 1 rank holding everything) -- what IS
          ! comparable is the GROWTH from this point to the next.
          dbgcv = sum(real(vorten, kind=8)**2)
          dbgcd = sum(real(divten, kind=8)**2)
          dbgcp = sum(real(phiten, kind=8)**2)
          if (myrank .eq. 0) print *,'DBGCHAIN 412 after tendget', &
             ' svor=',dbgcv,' sdiv=',dbgcd,' sphi=',dbgcp
          ! #endregion
!
!--------------------------------------------------------
!  incremental initialization
!
        if(doincr)then
          do m=1,mlistnum
            mf=mlist(m)
            do n=mf,jtrun_l
              do j =1,2
                do k =1,levp
                  phiten(k,j,n,m)=phiten(k,j,n,m)-temold(k,j,n,m)
                  vorten(k,j,n,m)=vorten(k,j,n,m)-vorold(k,j,n,m)
                  divten(k,j,n,m)=divten(k,j,n,m)-divold(k,j,n,m)
                enddo
              enddo
            enddo
          enddo
        endif
!
!--------------------------------------------------------
!
!  do vertical transform
!
        call zx (evecin,vorten,divten,phiten,jtrun_l,jtmax_l,lev_l)
          ! #region agent log
          dbgcv = sum(real(vorten, kind=8)**2)
          dbgcd = sum(real(divten, kind=8)**2)
          dbgcp = sum(real(phiten, kind=8)**2)
          if (myrank .eq. 0) print *,'DBGCHAIN 413 after zx     ', &
             ' svor=',dbgcv,' sdiv=',dbgcd,' sphi=',dbgcp
          ! #endregion
!
!     conversion  structure of variables from forecast model to
!     initialization
!     model
!
!   begin initialize the first liz (now,liz=3) modes.
!

        do 110 L=1,nnmivm
          nnmical=.false.
          brank=0
          if( (L .ge. Lstart) .and. (L .le. Lend)) then
            k=L-Lstart+1
            nnmical=.true.
            brank=row_rank
          endif
!
          if(myrank .eq. 0) print *,'vertical mode l=',L
          ! #region agent log
          if (ic .eq. 1 .and. L .eq. 1) then
            dbg0 = loc(jtrun)
            dbg1 = loc(eval)
            dbg2 = loc(a)
            call geps_dbg_log(2, 104, myrank, dbg0, dbg1, dbg2)
          endif
          ! #endregion
!
! set convergence variable
!
          bal=0.
          ! #region agent log
          dbgx_in = 0.0d0; dbgx_out = 0.0d0; dbgev = 0.0d0
          dbgx2_tot = 0.0d0; dbgnd_tot = 0_8
          ! #endregion
!
!   do nondimensionalize of the variables
!         +2:nondimensionalize
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
          call vartran (vorten,divten,phiten,nw,jtrun_l,jtmax_l,levp,k, &
                       rad,omega,h(1,1,L),+2)
          ! #region agent log
          if (ic .eq. 1 .and. L .eq. 1) then
            dbg0 = loc(vorten)
            dbg1 = 0
            dbg2 = 0
            if (nnmical) dbg1 = 1
            call geps_dbg_log(5, 116, myrank, dbg0, dbg1, dbg2)
          endif
          ! #endregion
!
!   begin to initialize horizontal modes
!
          do 100 m=1,mlistnum
            mf=mlist(m)
!
!  calculate the size of symmetric and antisymmetric matrix which
!  include the gravity and rossby wave
!
            nbig = jtrun_l-mf+1
            ns   = 2*int((nbig+1)/2)+int(nbig/2)
            na   = 2*int(nbig/2)+int((nbig+1)/2)
!
!  do symmetric case
!
!  create variable vector
!       +1:symmetric, +2:creation
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call vartrix (vorten,divten,phiten,levp, &
                         k,x,ns,m,nbig,jtrun_l,jtmax_l,+1,+2)

          ! #region agent log
          if (ic .eq. 1 .and. L .eq. 1 .and. m .eq. 1) then
            dbg0 = loc(x)
            dbg1 = no
            dbg2 = brank
            call geps_dbg_log(3, 140, myrank, dbg0, dbg1, dbg2)
          endif
          ! #endregion
      call mpe2d_row_broadcast(x,no_l*2,brank)

!
!  construct coefficient matrix
!       +1:symmetric
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call coftrix (mf,L,mx,ns,a(1,1,L),b(1,1,L),c(1,1,L),jtrun_l,lev_l,+1)

      call mpe2d_row_broadcast(mx,no_l*no_l,brank)
!
!   fine the eigenvector and eigenvalues of the symmetric matrix
!
          ! #region agent log
          if (ic .eq. 1 .and. L .eq. 1 .and. m .eq. 1) then
            dbg0 = loc(mx)
            dbg1 = loc(eval)
            dbg2 = ns
            call geps_dbg_log(1, 154, myrank, dbg0, dbg1, dbg2)
          endif
          ! #endregion
            call eigen (mx,ns,eval,evec,epos,wk)
          ! #region agent log
          if (ic .eq. 1 .and. L .eq. 1 .and. m .eq. 1) then
            dbg0 = loc(evec)
            dbg1 = loc(epos)
            dbg2 = 1
            call geps_dbg_log(1, 155, myrank, dbg0, dbg1, dbg2)
          endif
          ! #endregion
!
!   perform nonlinear normal mode initialization
!
          ! #region agent log
          dbgx_in = max(dbgx_in, maxval(abs(x(1:ns*2))))
          dbgev   = max(dbgev,   maxval(abs(evec(1:ns*ns))))
          dbgbal_prev = bal
          dbgev_s = maxval(abs(eval(1:ns)))
          dbgx2   = sum(x(1:ns*2)**2)
          ! global total, to compare against the GPU's locid 407,
          ! which sums over ALL ind at once (nnlist holds ns and na)
          dbgx2_tot = dbgx2_tot + dbgx2
          dbgnd_tot = dbgnd_tot + int(ns*2, 8)
          ! #endregion
            call nnmi (no_l,x,epos,eval,evec,ns,wc,wd,bal,cutfreq)
          ! #region agent log
          dbgx_out = max(dbgx_out, maxval(abs(x(1:ns*2))))
          dbgbal_d = bal - dbgbal_prev
          ! per-(mf,L) distribution, not just the global maxima - `bal` is a SUM,
          ! so comparable maxima say nothing about comparable sums (handoff layer 9)
          if (myrank .eq. 0 .and. ic .eq. 1 .and. L .eq. 1) &
             print *,'DBGDIST mf=',mf,' ns=',ns,' maxeval=',dbgev_s, &
                     ' sumx2=',dbgx2,' dbal=',dbgbal_d
          ! #endregion
!
!   decompose the variable vector back to 3 individual variable
!       +1:symmetric , -2:decompose
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call vartrix (vorten,divten,phiten,levp, &
                          k,x,ns,m,nbig,jtrun_l,jtmax_l,+1,-2)
!
!   do antisymmtric case
!
!   construct variable vector from 3 individual variable
!       -1:antisymmetric , +2 : creation
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call vartrix (vorten,divten,phiten,levp, &
                          k,x,na,m,nbig,jtrun_l,jtmax_l,-1,+2)

      call mpe2d_row_broadcast(x,no_l*2,brank)
!
!  construct coefficient matrix
!       -1:antisymmetric
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call coftrix (mf,L,mx,na,a(1,1,L),b(1,1,L),c(1,1,L),jtrun_l,lev_l,-1)

      call mpe2d_row_broadcast(mx,no_l*no_l,brank)
!
!   fine the eigenvector and eigenvalues of the antisymmetric matrix
!
            call eigen (mx,na,eval,evec,epos,wk)
!
!   perform nonlinear normal mode initialization
!
          ! #region agent log
          dbgx_in = max(dbgx_in, maxval(abs(x(1:na*2))))
          dbgev   = max(dbgev,   maxval(abs(evec(1:na*na))))
          dbgx2_tot = dbgx2_tot + sum(x(1:na*2)**2)
          dbgnd_tot = dbgnd_tot + int(na*2, 8)
          ! #endregion
            call nnmi (no_l,x,epos,eval,evec,na,wc,wd,bal,cutfreq)
          ! #region agent log
          dbgx_out = max(dbgx_out, maxval(abs(x(1:na*2))))
          ! #endregion
!
!   decompose the variable vector back to 3 individual variable
!       +1:symmetric , -2:decompose
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call vartrix (vorten,divten,phiten,levp, &
                          k,x,na,m,nbig,jtrun_l,jtmax_l,-1,-2)
            bal_tmp(mf) = bal
            bal = 0.
!
 100      continue
!
          ! amdflang even at -O0 reloads &jtrun from a spill that
          ! mpe_unify/MPI clobbers (SEGV movslq (%rax) rax=NULL).
          ! #region agent log
          if (ic .eq. 1 .and. L .eq. 1) then
            dbg0 = loc(jtrun)
            dbg1 = jtrun_l
            dbg2 = loc(bal_tmp)
            call geps_dbg_log(17, 258, myrank, dbg0, dbg1, dbg2)
          endif
          ! #endregion
          call mpe_unify(bal_tmp(1),1,jtrun_l,3,mpe_double)
          ! #region agent log
          if (ic .eq. 1 .and. L .eq. 1) then
            dbg0 = loc(jtrun)
            dbg1 = jtrun_l
            dbg2 = loc(bal_tmp)
            call geps_dbg_log(17, 259, myrank, dbg0, dbg1, dbg2)
          endif
          ! #endregion
          do mf = 1, jtrun_l
            bal = bal + bal_tmp(mf)
          enddo
          ! #region agent log
          if (ic .eq. 1 .and. L .eq. 1) then
            dbg0 = loc(bal_tmp)
            dbg1 = jtrun_l
            dbg2 = 0
            call geps_dbg_log(17, 280, myrank, dbg0, dbg1, dbg2)
          endif
          ! #endregion
!
          if(myrank .eq. 0) print *,'bal=',bal
          ! #region agent log: baseline for the GPU probes 400/401 and the
          ! orthonormality check - magnitudes of the same quantities on CPU.
          if(myrank .eq. 0) print *,'DBGCPU maxx_in=',dbgx_in, &
                                    ' maxx_out=',dbgx_out,' maxevec=',dbgev
          if(myrank .eq. 0) print *,'DBGTOT sumx2_tot=',dbgx2_tot, &
                                    ' nused=',dbgnd_tot,' bal=',bal
          ! #endregion
!
!   dimensionlize the variables
!     -2:dimensionalize
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
          call vartran (vorten,divten,phiten,nw,jtrun_l,jtmax_l,levp,k, &
                        rad,omega,h(1,1,L),-2)
!
 110    continue
!
!   initial none initialized mode
!

!2dMPI>
!!         kL=1
!!         kk=Llist(nnmivm+1)
!!         if(kk .le.levp)then
!!            kL=nnmivm+1
!!         endif
!2dMPI<

         do m=1,mlistnum
           mf=mlist(m)
           do n=mf,jtrun_l
             do j = 1, 2
!              do k = nnmivm+1,lev
               do k = 1,levp
                 kk=Llist(k)
                 if ( kk .gt. nnmivm ) then
                   vorten(k,j,n,m)= 0.
                   divten(k,j,n,m)= 0.
                   phiten(k,j,n,m)= 0.
                 endif
               enddo
             enddo
           enddo
         enddo
!
!  conversion  structure of variables ( phiten,vorten,divten)
!
!   vertical transform back
!
        call zx (evectr,vorten,divten,phiten,jtrun_l,jtmax_l,lev_l)
!
!   add the correction to variables
!
        call correct (vorten,divten,phiten,tmcor,pmcor,vornow, &
                     divnow,temnow,plnow,jtrun_l,jtmax_l,lev_l)
!
        do m=1,mlistnum
          mf=mlist(m)
          if (mf.eq.1) then
            do k = 1, levp
              vornow(k,1,1,m)= 0.
              vornow(k,2,1,m)= 0.
              divnow(k,1,1,m)= 0.
              divnow(k,2,1,m)= 0.
            enddo
          endif
        enddo
!
!   transform back to phyical space
!
!!        call joinsr(wss,vornow,divnow,temnow,dummy,jtrun,jtmax,levp &
!!                   ,mlistnum,3,1)
!!        call transr(jtrun,jtmax,nx,my,my_max,levp,poly,wss,cc,3,nsizey)
!!        call ujoinsr(cc,rvor,rdiv,tt,dummy,nx,my_max,lev,jlistnum,3,1)
        call transr(jtrun_l,jtmax_l,nx_l,my_l,my_max_l,levp,poly,vornow,cc,1,nsizey)
        call ujoinsr(cc,rvor,dummy,dummy,dummy,nx_l,my_max_l,lev_l,jlistnum,1,1)
        call transr(jtrun_l,jtmax_l,nx_l,my_l,my_max_l,levp,poly,divnow,cc,1,nsizey)
        call ujoinsr(cc,rdiv,dummy,dummy,dummy,nx_l,my_max_l,lev_l,jlistnum,1,1)
        call transr(jtrun_l,jtmax_l,nx_l,my_l,my_max_l,levp,poly,temnow,cc,1,nsizey)
        call ujoinsr(cc,tt,dummy,dummy,dummy,nx_l,my_max_l,lev_l,jlistnum,1,1)
        call transr1(jtrun_l,jtmax_l,nx_l,my_l,my_max_l,poly,plnow,pt,nsizey)
!
!  compute zonal and meridional gradients of terrain pressure
!
        call trngra (jtrun_l,jtmax_l,nx_l,my_l,my_max_l,cim,poly,dpoly,plnow, &
                    dlpl,dtpl,nsizey)
!
!   transform spectrum vorticity , divergence to physical u , v
!
        call tranuv (jtrun_l,jtmax_l,nx_l,my_l,my_max_l,levp,onocos,wcfac,wdfac &
                   ,poly,dpoly,vornow,divnow,ut,vt,nsizey)
!
        lab='pt '
        call  check (pt,nx_l,my_l,my_max_l,lab)
 120  continue
!
      do m=1,mlistnum
        mf=mlist(m)
        do n=mf,jtrun_l
          do k= 1, levp
            vorold(k,1,n,m)= vornow(k,1,n,m)
            vorold(k,2,n,m)= vornow(k,2,n,m)
            divold(k,1,n,m)= divnow(k,1,n,m)
            divold(k,2,n,m)= divnow(k,2,n,m)
            temold(k,1,n,m)= temnow(k,1,n,m)
            temold(k,2,n,m)= temnow(k,2,n,m)
          enddo
        enddo
      enddo
!
!byl      do m=1,mlistnum
!        mf=mlist(m)
!        do n=mf,jtrun
!          do k= 1, levp*ncld
!            qold(k,1,n,m)= qnow(k,1,n,m)
!            qold(k,2,n,m)= qnow(k,2,n,m)
!          enddo
!        enddo
!byl      enddo
!
      do m=1,mlistnum
        mf=mlist(m)
        do n=mf,jtrun_l
          plold(n,m,1)= plnow(n,m,1)
          plold(n,m,2)= plnow(n,m,2)
        enddo
      enddo
!
      deallocate(bal_tmp)
      return
      end
