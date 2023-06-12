      subroutine siimpl( jtrun,jtmax,lev,dta,ptmean,dsigma,spalm,eps4    &
                , eigval,evecin,evectr,arrhyd,arsddt,temold,divold,plold &
                , temnow,divnow,plnow,temten,divten,plten,alpha )
!
!
!  computes corrections to explicit tendencies to convert model to a
!  semi-implicit model
!
!  *** input ***
!
!  dta: time step in seconds
!  ptmean: mean terrain pressure chosen for best stability properties
!  dsigma: thickness of sigma layers
!  spalm: mean energy conversion terms for each level
!  eps4: spherical harmonic laplacian operator
!  eigval: gravity mode phase speed eigenvalues
!  evecin: inverse of gravity mode eigenvector matrix
!  evectr: gravity mode eigenvector matrix
!  arrhyd: linearized hydrostatic matrix
!  arsddt: linerized vertical temperature advection matrix
!  temold: (t-dt) spectral temperature
!  divold: (t-dt) spectral divergence
!  plold: (t-dt) spectral terrain pressure
!  temnow: current time temperature
!  divnow: current time divergence
!  plnow: current time terrain pressure
!  temten: explicit temperature tendency
!  divten: explicit divergence tendency
!  plten: explicit terrain pressure tendency
!
! *** output ***
!
!  temten: semi-implicit temperature tendency
!  divten: semi-implicit divergence tendency
!  plten: semi-implicit terrain pressure tendency
!
! **************************************************
!

!CWB2017 2dMPI version

      use index
      use paramt
      use const,only : eps4L,RTYPE
      use spec ,only : plnowL,ploldL,pltenL,jtwvp

      implicit  none
      integer   jtrun,jtmax,lev
 
      real(kind=RTYPE) dsigma(lev,2),eps4(jtrun,jtmax),eigval(lev),evecin(lev,lev) &
      ,     evectr(lev,lev),arrhyd(lev,lev),arsddt(lev,lev),spalm(lev)
 
      real(kind=RTYPE) temold(levp,2,jtrun,jtmax),divold(levp,2,jtrun,jtmax) &
      ,                temnow(levp,2,jtrun,jtmax),divnow(levp,2,jtrun,jtmax) &
      ,                temten(levp,2,jtrun,jtmax),divten(levp,2,jtrun,jtmax)
      real(kind=RTYPE) plold(jtrun,jtmax,2),plnow(jtrun,jtmax,2),plten(jtrun,jtmax,2)
 
      real(kind=RTYPE) divavg(lev,2),phiave(lev,2,jtrun),eps4e(lev,jtrun)
 
      integer   m,mf,k,n,l,j
      real      alpha
      real(kind=RTYPE) dta,dd,odd,dd2,ptmean,tem,s1,s2,d1,d2,dp

      real(kind=RTYPE) wrk1(lev,2,jtp),wrk2(lev,2,jtp),wrk3(lev,2,jtp),&
                       wrk4(lev,2,jtp),wrk5(lev,2,jtp),wrk6(lev,2,jtp)

      dd = alpha*dta
      odd= 1.0/dd
      dd2= dd*dd

#ifdef MULTIPLE
      call mpe2d_reshape_pl_multi(plten, plnow, plold, pltenL, plnowL, ploldL)

      call mpe2d_transpose_siimpl_multi(temold,temnow,temten,divold,divnow,divten, &
                                   wrk1  ,wrk2  ,wrk3  ,wrk4  ,wrk5  ,wrk6  , &
                                   levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
#else
      call mpe2d_reshape_pl(plten, pltenL)
      call mpe2d_reshape_pl(plnow, plnowL)
      call mpe2d_reshape_pl(plold, ploldL)

      call mpe2d_transpose_siimpl(temold, &
                                   wrk1 , &
                                   levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
      call mpe2d_transpose_siimpl(temnow, &
                                   wrk2 , &
                                   levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
      call mpe2d_transpose_siimpl(temten, &
                                   wrk3 , &
                                   levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
      call mpe2d_transpose_siimpl(divold, &
                                   wrk4 , &
                                   levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
      call mpe2d_transpose_siimpl(divnow, &
                                   wrk5 , &
                                   levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
      call mpe2d_transpose_siimpl(divten, &
                                   wrk6 , &
                                   levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
#endif

      do 700 m = 1,jtlen
      n=jtwvp(m)

      do 100 k = 1, lev
      tem= dd2*eigval(k)
      eps4e(k,n)= 1.0/(1.0+tem*eps4L(m))
  100 continue

        if (n.ne.1) then
          do k=1,lev
             divavg(k,1)= wrk1(k,1,m)+dd*wrk3(k,1,m)-wrk2(k,1,m)
             divavg(k,2)= wrk1(k,2,m)+dd*wrk3(k,2,m)-wrk2(k,2,m)
          enddo
          do 130 k = 1, lev
            s1= spalm(k)*(ploldL(m,1)+dd*pltenL(m,1)-plnowL(m,1))
            s2= spalm(k)*(ploldL(m,2)+dd*pltenL(m,2)-plnowL(m,2))
            do 125 L = 1, lev
              s1= s1+arrhyd(k,L)*divavg(L,1)
              s2= s2+arrhyd(k,L)*divavg(L,2)
  125       continue
            phiave(k,1,n)= eps4L(m)*s1
            phiave(k,2,n)= eps4L(m)*s2
            wrk6(k,1,m)= dd*(phiave(k,1,n)+wrk6(k,1,m))  &
                           - wrk5(k,1,m)+wrk4(k,1,m)
            wrk6(k,2,m)= dd*(phiave(k,2,n)+wrk6(k,2,m))  &
                           - wrk5(k,2,m)+wrk4(k,2,m)
  130     continue
        endif
!
! transform time averaged divergence to eigenspace and compute
! semi-implicit values
!
        if (n.eq.1) then

          do 160 l = 1, lev
            wrk6(l,1,m)= 0.0
            wrk6(l,2,m)= 0.0
  160     continue

        else

          do 170 L = 1, lev
            d1= evecin(L,1)*wrk6(1,1,m)
            d2= evecin(L,1)*wrk6(1,2,m)
            do 162 k = 2, lev
              d1= d1+evecin(L,k)*wrk6(k,1,m)
              d2= d2+evecin(L,k)*wrk6(k,2,m)
  162       continue
            divavg(L,1)= d1*eps4e(L,n)
            divavg(L,2)= d2*eps4e(L,n)
  170     continue

          do 195 k = 1, lev
            s1= evectr(k,1)*divavg(1,1)
            s2= evectr(k,1)*divavg(1,2)
            do 193 L = 2, lev
              s1= s1+evectr(k,L)*divavg(L,1)
              s2= s2+evectr(k,L)*divavg(L,2)
  193       continue
            wrk6(k,1,m)= s1
            wrk6(k,2,m)= s2
  195     continue
        endif
!
! add contributions of time averaged divergence to temperature
! tendency
!
        if (n.ne.1) then
          do 205 k = 1, lev
            s1=wrk3(k,1,m)
            s2=wrk3(k,2,m)
            do 200 L = 1, lev
              s1= s1-arsddt(k,L)*wrk6(L,1,m)
              s2= s2-arsddt(k,L)*wrk6(L,2,m)
  200       continue
            wrk3(k,1,m)= s1
            wrk3(k,2,m)= s2
  205     continue
        endif
!
! add contribution of vertically integrated time-averaged divergence
! to surface pressure tendency.  convert time-averaged divergence to
! semi-implicit divergence tendency.
!
        if (n.ne.1) then
         do 230 j = 1, 2
         do 230 k = 1, lev
          dp= dsigma(k,1)*ptmean+dsigma(k,2)
          pltenL(m,j)= pltenL(m,j)-dp*wrk6(k,j,m)
          wrk6(k,j,m)=(wrk5(k,j,m)+wrk6(k,j,m)-wrk4(k,j,m))*odd
  230    continue
        endif

  700 continue   ! end of large m loop
 
      call mpe2d_reshape_pl_back(pltenL, plten)

#ifdef MULTIPLE
      call mpe2d_transpose_siimpl_back_multi(wrk3,wrk6,temten,divten,levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
#else
      call mpe2d_transpose_siimpl_back(wrk3,temten,levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
      call mpe2d_transpose_siimpl_back(wrk6,divten,levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
#endif

      return
      end
