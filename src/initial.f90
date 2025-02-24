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
      real      bal_tmp(jtrun)
      character lab*10,lrec*16

      integer   mlmax2,j,k,l,ic,m,mf,n,ns,na,nbig,kk,kL
      real      bal
      logical   nnmical
      integer   brank   !root rank of row_broadcast
!
      phiten=0.
      mlmax2 = mlmax*2
      if(myrank .eq. 0) print *,'jrtun=',jtrun,' mlmax=',mlmax
      lab='pt'
      call  check(pt,nx,my,my_max,lab)
!
!  define constants and comput coefficients for initializatin
!
      call  inicons (rad,omega,eigval,lev,jtrun,jtmax, &
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
        call tendget (phiten)
!
!--------------------------------------------------------
!  incremental initialization
!
        if(doincr)then
          do m=1,mlistnum
            mf=mlist(m)
            do n=mf,jtrun
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
        call zx (evecin,vorten,divten,phiten,jtrun,jtmax,lev)
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
!
! set convergence variable
!
          bal=0.
!
!   do nondimensionalize of the variables
!         +2:nondimensionalize
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
          call vartran (vorten,divten,phiten,nw,jtrun,jtmax,levp,k, &
                       rad,omega,h(1,1,L),+2)
!
!   begin to initialize horizontal modes
!
          do 100 m=1,mlistnum
            mf=mlist(m)
!
!  calculate the size of symmetric and antisymmetric matrix which
!  include the gravity and rossby wave
!
            nbig = jtrun-mf+1
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
                         k,x,ns,m,nbig,jtrun,jtmax,+1,+2)

      call mpe2d_row_broadcast(x,no*2,brank)

!
!  construct coefficient matrix
!       +1:symmetric
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call coftrix (mf,L,mx,ns,a(1,1,L),b(1,1,L),c(1,1,L),jtrun,lev,+1)

      call mpe2d_row_broadcast(mx,no*no,brank)
!
!   fine the eigenvector and eigenvalues of the symmetric matrix
!
            call eigen (mx,ns,eval,evec,epos,wk)
!
!   perform nonlinear normal mode initialization
!
            call nnmi (no,x,epos,eval,evec,ns,wc,wd,bal,cutfreq)
!
!   decompose the variable vector back to 3 individual variable
!       +1:symmetric , -2:decompose
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call vartrix (vorten,divten,phiten,levp, &
                          k,x,ns,m,nbig,jtrun,jtmax,+1,-2)
!
!   do antisymmtric case
!
!   construct variable vector from 3 individual variable
!       -1:antisymmetric , +2 : creation
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call vartrix (vorten,divten,phiten,levp, &
                          k,x,na,m,nbig,jtrun,jtmax,-1,+2)

      call mpe2d_row_broadcast(x,no*2,brank)
!
!  construct coefficient matrix
!       -1:antisymmetric
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call coftrix (mf,L,mx,na,a(1,1,L),b(1,1,L),c(1,1,L),jtrun,lev,-1)

      call mpe2d_row_broadcast(mx,no*no,brank)
!
!   fine the eigenvector and eigenvalues of the antisymmetric matrix
!
            call eigen (mx,na,eval,evec,epos,wk)
!
!   perform nonlinear normal mode initialization
!
            call nnmi (no,x,epos,eval,evec,na,wc,wd,bal,cutfreq)
!
!   decompose the variable vector back to 3 individual variable
!       +1:symmetric , -2:decompose
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
            call vartrix (vorten,divten,phiten,levp, &
                          k,x,na,m,nbig,jtrun,jtmax,-1,-2)
            bal_tmp(mf) = bal
            bal = 0.
!
 100      continue
!
          call mpe_unify(bal_tmp,1,jtrun,3,mpe_double)
          do mf = 1, jtrun
            bal = bal + bal_tmp(mf)
          enddo
!
          if(myrank .eq. 0) print *,'bal=',bal
!
!   dimensionlize the variables
!     -2:dimensionalize
!
!      if(L.eq.Llist(L)) &
      if(nnmical) &
          call vartran (vorten,divten,phiten,nw,jtrun,jtmax,levp,k, &
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
           do n=mf,jtrun
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
        call zx (evectr,vorten,divten,phiten,jtrun,jtmax,lev)
!
!   add the correction to variables
!
        call correct (vorten,divten,phiten,tmcor,pmcor,vornow, &
                     divnow,temnow,plnow,jtrun,jtmax,lev)
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
        call transr(jtrun,jtmax,nx,my,my_max,levp,poly,vornow,cc,1,nsizey)
        call ujoinsr(cc,rvor,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
        call transr(jtrun,jtmax,nx,my,my_max,levp,poly,divnow,cc,1,nsizey)
        call ujoinsr(cc,rdiv,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
        call transr(jtrun,jtmax,nx,my,my_max,levp,poly,temnow,cc,1,nsizey)
        call ujoinsr(cc,tt,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
        call transr1(jtrun,jtmax,nx,my,my_max,poly,plnow,pt,nsizey)
!
!  compute zonal and meridional gradients of terrain pressure
!
        call trngra (jtrun,jtmax,nx,my,my_max,cim,poly,dpoly,plnow, &
                    dlpl,dtpl,nsizey)
!
!   transform spectrum vorticity , divergence to physical u , v
!
        call tranuv (jtrun,jtmax,nx,my,my_max,levp,onocos,wcfac,wdfac &
                   ,poly,dpoly,vornow,divnow,ut,vt,nsizey)
!
        lab='pt '
        call  check (pt,nx,my,my_max,lab)
 120  continue
!
      do m=1,mlistnum
        mf=mlist(m)
        do n=mf,jtrun
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
        do n=mf,jtrun
          plold(n,m,1)= plnow(n,m,1)
          plold(n,m,2)= plnow(n,m,2)
        enddo
      enddo
!
      return
      end
