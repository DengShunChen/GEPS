      subroutine zimadv (my,my_max,lev,ncld                     &
                       , jtrun,jtmax,dt,poly,onocos,w,uz        &
                       , s1,x1now,x1old,s2,x2now,x2old,nsize)
!
! this extension of the semi implicit scheme is based on the report
! "the design and performance of the new ecmwf operational model"
! by a.j. simons and m. jarraud in the european centre for medium
! range weather forecasts' seminar 1983, "numerical methods for
! weather predition", 5 - 9 september, volume 2, pp113-164
! (see section 4: the time-stepping scheme, pp120-130)
!
! ***input***
!
!  my: no. of gaussian latitudes
!  jtrun:  zonal wavenumber truncation limit
!  mlmax: total number of triangular truncation spherical harm. coeff.
!  dt: time step
!  poly: associated legendre coefficients
!  onocos:  1.0/cos(lat)**2
!  w: gaussian quadrature weights
!  uz: zonal mean u-comp for each latitude and level
!  s: spectral tendency array
!  xnow: current time level spectral dependent variable
!  xold: t-dt time level spectral dependent variable
!
!  ***output***
!
!  s: adjusted spectral tendency array
!
! ******************************************************************
!
      use index
!     use paramt
!
      implicit none

      integer  my,my_max,lev,ncld,jtrun,jtmax,nsize
      integer  myhalf,lev2,na2,na4,k,j,m,mf,l,n,n2,nk,kk,nb,jchk
      integer  jje,jlistnum_fj,j_fj,l_fj,llistnum_fj,jx0,jx,ll,j2,n3,n4,lchk,lle,jj

      real     sa00,sa10,sa20,sa30,sb00,sb10,sb20,sb30,dt,tem,workx,workdt,xx

!ch   real     uz(my,lev),w(my),onocos(my)
      real     uz(my,levf),w(my),onocos(my)
      real     x1now(lev,2,jtrun,jtmax), x1old(lev,2,jtrun,jtmax)
      real     x2now(lev*ncld,2,jtrun,jtmax),x2old(lev*ncld,2,jtrun,jtmax)
      real     poly(jtrun,my/2,jtmax)
      real     s1(lev,2,jtrun,jtmax)
      real     s2(lev*ncld,2,jtrun,jtmax)
      real     wcc2(lev,2,1+ncld,my)
!
      real     work(lev,my)
      real     wss(lev,2,2*(1+ncld),jtrun)
      real     wcc_fk(lev,2,2*(1+ncld),my_max*nsize)

      real     ws2(lev,2,2*(1+ncld),jtrun)
!
      integer  jlist_fj(my/2)
      real     fj_wcc_fk(lev*2*2*(1+ncld),my_max*nsize)
      real     fj_poly(jtrun,my/2),fj_poly_sum(my/2,jtrun),fj_poly_dif(my/2,jtrun)
      real     fj_wss(lev*2*2*(1+ncld),jtrun)
      real     fj_ws2(lev*2*2*(1+ncld),jtrun)
      real     fj_wccSUM(lev*2*(1+ncld),my), fj_wccDIF(lev*2*(1+ncld),my)
      real     fj_wss_sum(lev*2*2*(1+ncld),jtrun),fj_wss_dif(lev*2*2*(1+ncld),jtrun)
!
! calculate the complex fourier series for s = cc
! and for (xnow - xold) = bb
!
      myhalf=my/2
      lev2=lev*2
      na2=1+ncld
      na4=na2*2
!
      do k = 1, lev
         KK=Llist(k)
      do j = 1, my
!ch     work(k,j)= uz(j,k)*onocos(j)
        work(k,j)= uz(j,KK)*onocos(j)
      enddo
      enddo
!
!--- do start
      do m = 1, mlistnum
        mf=mlist(m)
!
        do l = mf, jtrun
        do k = 1, lev2
          wss(k,1,1,l)=s1(k,1,l,m)
          wss(k,1,2,l)=x1now(k,1,l,m)-x1old(k,1,l,m)
          wss(k,1,3,l)=s2(k,1,l,m)
          wss(k,1,4,l)=x2now(k,1,l,m)-x2old(k,1,l,m)
        enddo
        enddo
!
        do n = 1, ncld
         nk=(n-1)*lev
        do l = mf, jtrun
        do k = 1, lev
          kk=nk+k
          wss(k,1,2+n,l)=s2(kk,1,l,m)
          wss(k,2,2+n,l)=s2(kk,2,l,m)
          wss(k,1,2+ncld+n,l)=x2now(kk,1,l,m)-x2old(kk,1,l,m)
          wss(k,2,2+ncld+n,l)=x2now(kk,2,l,m)-x2old(kk,2,l,m)
        enddo
        enddo
        enddo

!
        do l = mf, jtrun, 2
        do k = 1, lev2*na4
         ws2(k,1,1,l) = wss(k,1,1,l)
        enddo
        enddo

        do l = mf+1, jtrun, 2
        do k = 1, lev2*na4
         ws2(k,1,1,l) = -wss(k,1,1,l)
        enddo
        enddo

!ibm---beg  .. using alternate unrolling method with cache block
!        .. the innermost loop is DDOT kernel

        do j = 1, my
        do k = 1, lev2*na4
          wcc_fk(k,1,1,j)= 0.0
        enddo
        enddo

!
        nb=32
        jchk= iand(myhalf, 1)
        jje = myhalf-jchk
!
        fj_wcc_fk = 0.0
        fj_poly = 0.0
        fj_wss = 0.0

        jlistnum_fj = 0
        do j = 1,jje
           if( mf.le.mtrundef(j) ) then
              jlistnum_fj = jlistnum_fj + 1
              jlist_fj(jlistnum_fj) = j
           endif
        enddo

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do l = mf,jtrun
              l_fj = l - mf + 1
              fj_poly(l_fj,j_fj) = poly(l,j,m)
           enddo
        enddo

        do l = mf,jtrun
           l_fj = l - mf + 1
           do k = 1,lev2*na4
              fj_wss(k,l_fj) = wss(k,1,1,l)
           enddo
        enddo

        llistnum_fj = jtrun - mf + 1
        call dgemm('n','n',lev2*na4,jlistnum_fj,llistnum_fj,1.0d+0,fj_wss,lev*2*2*(1+ncld), &
             fj_poly,jtrun,1.0d+0,fj_wcc_fk,lev*2*2*(1+ncld))

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           do k = 1,lev2*na4
              wcc_fk(k,1,1,j) = wcc_fk(k,1,1,j) + fj_wcc_fk(k,j_fj)
           enddo
        enddo

        fj_ws2 = 0.0
        fj_wcc_fk = 0.0

        do l = mf,jtrun
           l_fj = l - mf + 1
           do k = 1,lev2*na4
              fj_ws2(k,l_fj) = ws2(k,1,1,l)
           enddo
        enddo

        llistnum_fj = jtrun - mf + 1
        call dgemm('n','n',lev2*na4,jlistnum_fj,llistnum_fj,1.0d+0,fj_ws2,lev*2*2*(1+ncld), &
             fj_poly,jtrun,1.0d+0,fj_wcc_fk,lev*2*2*(1+ncld))

        do j_fj = 1,jlistnum_fj
           j = jlist_fj(j_fj)
           jx0=my-j+1
           do k = 1,lev2*na4
              wcc_fk(k,1,1,jx0) = wcc_fk(k,1,1,jx0) + fj_wcc_fk(k,j_fj)
           enddo
        enddo

!
! odd number
!
        if( jchk.eq.1 )then
          j = myhalf
          jx= my-j+1
        if( mf.le.mtrundef(j) ) then
        do kk=1,lev2*na4,nb
        do ll=mf,jtrun,nb
          do k = kk,min(kk+nb-1,lev2*na4),4
            sa00=wcc_fk(k,1,1,j)
            sa10=wcc_fk(k+1,1,1,j)
            sa20=wcc_fk(k+2,1,1,j)
            sa30=wcc_fk(k+3,1,1,j)
            sb00=wcc_fk(k,1,1,jx)
            sb10=wcc_fk(k+1,1,1,jx)
            sb20=wcc_fk(k+2,1,1,jx)
            sb30=wcc_fk(k+3,1,1,jx)
            do l = ll,min(ll+nb-1,jtrun)
              sa00=sa00+poly(l,j,  m)*wss(k,1,1,l)
              sa10=sa10+poly(l,j,  m)*wss(k+1,1,1,l)
              sa20=sa20+poly(l,j,  m)*wss(k+2,1,1,l)
              sa30=sa30+poly(l,j,  m)*wss(k+3,1,1,l)
              sb00=sb00+poly(l,j,  m)*ws2(k,1,1,l)
              sb10=sb10+poly(l,j,  m)*ws2(k+1,1,1,l)
              sb20=sb20+poly(l,j,  m)*ws2(k+2,1,1,l)
              sb30=sb30+poly(l,j,  m)*ws2(k+3,1,1,l)
            enddo
            wcc_fk(k,1,1,j  )=sa00
            wcc_fk(k+1,1,1,j  )=sa10
            wcc_fk(k+2,1,1,j  )=sa20
            wcc_fk(k+3,1,1,j  )=sa30
            wcc_fk(k,1,1,jx)=sb00
            wcc_fk(k+1,1,1,jx)=sb10
            wcc_fk(k+2,1,1,jx)=sb20
            wcc_fk(k+3,1,1,jx)=sb30
          enddo
        enddo
        enddo
        endif
        endif

!
! calculate c3 = the fourier coefficients of the new tendency
!
! for fully implicit calculation leave work alone and define
! bb below as (.5*work*bb)
!
      tem=2.0*(mf-1.0)
!
        do j = 1, my
        if( mf.le.mtrundef(j) ) then
        do k = 1, lev
         workx=work(k,j)*tem
         workdt=workx*dt
         xx= 1.0/(1.0+workdt*workdt)
         wcc2(k,1,1,j)=                                            &
                (wcc_fk(k,1,1,j)-(workx*wcc_fk(k,2,2,j))+workdt    &
               *(wcc_fk(k,2,1,j)+(workx*wcc_fk(k,1,2,j))))*xx
         wcc2(k,2,1,j)=                                            &
                (wcc_fk(k,2,1,j)+(workx*wcc_fk(k,1,2,j))-workdt    &
               *(wcc_fk(k,1,1,j)-(workx*wcc_fk(k,2,2,j))))*xx
        enddo

        do n = 1, ncld
         n2=1+n
         n3=2+n
         n4=2+ncld+n
        do k = 1, lev
         workx=work(k,j)*tem
         workdt=workx*dt
         xx= 1.0/(1.0+workdt*workdt)
         wcc2(k,1,n2,j)=                                           &
                (wcc_fk(k,1,n3,j)-(workx*wcc_fk(k,2,n4,j))+workdt  &
               *(wcc_fk(k,2,n3,j)+(workx*wcc_fk(k,1,n4,j))))*xx
         wcc2(k,2,n2,j)=                                           &
                (wcc_fk(k,2,n3,j)+(workx*wcc_fk(k,1,n4,j))-workdt  &
               *(wcc_fk(k,1,n3,j)-(workx*wcc_fk(k,2,n4,j))))*xx
        enddo
        enddo
        endif
        enddo
!
!
        do l = 1, jtrun
        do k = 1, lev2*na4
          wss(k,1,1,l)=0.
        enddo
        enddo

      lchk = iand(jtrun - mf + 1,1)
      lle = jtrun - lchk

      fj_wccSUM = 0.0
      fj_wccDIF = 0.0
      fj_poly_sum = 0.0
      fj_poly_dif = 0.0
      fj_wss_sum = 0.0
      fj_wss_dif = 0.0

      do j_fj = 1,jlistnum_fj
         j = jlist_fj(j_fj)
         jj=my-j+1
         do k = 1, lev2*na2
            fj_wccSUM(k,j_fj)=wcc2(k,1,1,j)+wcc2(k,1,1,jj)
            fj_wccDIF(k,j_fj)=wcc2(k,1,1,j)-wcc2(k,1,1,jj)
         enddo
      enddo

      if (jtrun .gt. mf) then
         do j_fj = 1,jlistnum_fj
            j = jlist_fj(j_fj)
            do l = mf,lle,2
               l_fj = (l - mf + 2) / 2
               fj_poly_sum(j_fj,l_fj) = poly(l,j,m)*w(j)
               fj_poly_dif(j_fj,l_fj) = poly(l+1,j,m)*w(j)
            enddo
         enddo

         llistnum_fj = l_fj

         call dgemm('n','n',lev2*na2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wccSUM,lev*2*(1+ncld), &
           fj_poly_sum,my/2,1.0d+0,fj_wss_sum,lev*2*2*(1+ncld))

         do l = mf,lle,2
            l_fj = (l - mf + 2) / 2
            do k = 1,lev2*na2
               wss(k,1,1,l) = wss(k,1,1,l) + fj_wss_sum(k,l_fj)
            enddo
         enddo

         call dgemm('n','n',lev2*na2,llistnum_fj,jlistnum_fj,1.0d+0,fj_wccDIF,lev*2*(1+ncld), &
              fj_poly_dif,my/2,1.0d+0,fj_wss_dif,lev*2*2*(1+ncld))

         do l = mf,lle,2
            l_fj = (l - mf + 2) / 2
            do k = 1,lev2*na2
               wss(k,1,1,l+1) = wss(k,1,1,l+1) + fj_wss_dif(k,l_fj)
            enddo
         enddo
      endif

      do j_fj = 1,jlistnum_fj
         j = jlist_fj(j_fj)
         fj_poly_sum(jtrun,j_fj) = poly(jtrun,j,m)*w(j)
      enddo

      do j_fj = 1,jlistnum_fj
         j = jlist_fj(j_fj)
         do k = 1, lev2*na2
            wss(k,1,1,jtrun  )=wss(k,1,1,jtrun  )+fj_poly_sum(j_fj,jtrun)*fj_wccSUM(k,j_fj)
         enddo
      enddo

        do l = mf, jtrun
        do k = 1, lev2
         s1(k,1,l,m)=wss(k,1,1,l)
        enddo

        do n=1,ncld
         nk=(n-1)*lev
        do k = 1, lev
         kk=nk+k
         s2(kk,1,l,m)=wss(k,1,1+n,l)
         s2(kk,2,l,m)=wss(k,2,1+n,l)
        enddo
        enddo
        enddo

      enddo
!

      return
      end
