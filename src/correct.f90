      subroutine correct(vorten,divten,phiten,tmcor,pmcor,vornow,  &
                         divnow,temnow,plnow,jtrun,jtmax,lev)
!
!  purpose :  add correction to variables for NNMI
!-----------------------------------------------------------------------
!  **** input ****
!  vorten  :  vorticity correction array in spectrum domaon
!  divten  :  divergence correction array in spectrum domain
!  phiten  :  geopotantial correction array in spectrum domain
!  pmcor   :  array for compute surface pressure correction
!  tmcor   :  array for compute temperature correction
!  lev     :  total vertical levels
!  **** output ****
!  vornow  : vorticity spectrum coefficient
!  divnow  : divergence spectrum coefficient
!  temnow  : temperature spectrum coefficient
!  plnow   : surface pressure spectrum coefficient
!
! modify to f90 in 2015 by C-H Lee
! modify 2015 by River Chen
!-------------------------------------------------------------------------
!
      use const, only : RTYPE
      use index
!
      real(kind=RTYPE) vorten(levp,2,jtrun,jtmax),  &
                       divten(levp,2,jtrun,jtmax),  &
                       phiten(levp,2,jtrun,jtmax),  &
                       vornow(levp,2,jtrun,jtmax),  &
                       divnow(levp,2,jtrun,jtmax),  &
                       temnow(levp,2,jtrun,jtmax),  &
                       temnow1(lev,2,jtrun,jtmax),  &
                       phiten1(lev,2,jtrun,jtmax),  &
                       plnow(jtrun,jtmax,2)
      real(kind=RTYPE) tmcor(lev,lev),pmcor(lev)

!2dMPI >


      call mpe2d_unify_spec_lev(phiten,phiten1,lev,levp,jtrun,jtmax,mlistnum,nsizex,row_comm)
      call mpe2d_unify_spec_lev(temnow,temnow1,lev,levp,jtrun,jtmax,mlistnum,nsizex,row_comm)
!2dMPI <

!
!  add vorticity and divergence correction to vorticity and divergence
!  spectrum array
!
      do 100 m=1,mlistnum
         mf=mlist(m)
      do 100 n=mf,jtrun
      do 100 k =1,levp
         vornow(k,1,n,m)=vornow(k,1,n,m)+vorten(k,1,n,m)
         vornow(k,2,n,m)=vornow(k,2,n,m)+vorten(k,2,n,m)
         divnow(k,1,n,m)=divnow(k,1,n,m)+divten(k,1,n,m)
         divnow(k,2,n,m)=divnow(k,2,n,m)+divten(k,2,n,m)
 100  continue
!
!
!  compute surface pressure correction from geopotential correction and
!  add to surface spectrum array
!
      do 110 m=1,mlistnum
         mf=mlist(m)
      do 110 n=mf,jtrun
      do 110 k =1,lev
!        plnow(n,m,1)=plnow(n,m,1)+pmcor(k)*phiten(k,1,n,m)
!        plnow(n,m,2)=plnow(n,m,2)+pmcor(k)*phiten(k,2,n,m)
         plnow(n,m,1)=plnow(n,m,1)+pmcor(k)*phiten1(k,1,n,m)
         plnow(n,m,2)=plnow(n,m,2)+pmcor(k)*phiten1(k,2,n,m)
 110  continue
!
!
!  compute temperature correction from geopotential correction and
!  add to temperature spectrum array
!
      do 140 m=1,mlistnum
         mf=mlist(m)
      do 140 n=mf,jtrun
      do 140 i=1,lev
      do 140 k=1,lev
!       temnow(k,1,n,m)=temnow(k,1,n,m)+tmcor(k,i)*phiten(i,1,n,m)
!       temnow(k,2,n,m)=temnow(k,2,n,m)+tmcor(k,i)*phiten(i,2,n,m)
        temnow1(k,1,n,m)=temnow1(k,1,n,m)+tmcor(k,i)*phiten1(i,1,n,m)
        temnow1(k,2,n,m)=temnow1(k,2,n,m)+tmcor(k,i)*phiten1(i,2,n,m)
 140  continue
!

!2dMPI >
      do m=1,mlistnum
         mf=mlist(m)
      do n=mf,jtrun
      do k=1,levp
         KK=Llist(k)
         temnow(k,1,n,m)=temnow1(KK,1,n,m)
         temnow(k,2,n,m)=temnow1(KK,2,n,m)
      enddo
      enddo
      enddo
!2dMPI <

      return
      end
