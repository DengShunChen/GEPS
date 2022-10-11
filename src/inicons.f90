      subroutine inicons(rad,omega,eigval,lev,jtrun,jtmax, &
                         nw,a,b,c,h,nnmivm)
!
!  purpose:  define constants and calculate the coefficients for
!            initialization modual
!----------------------------------------------------------------------
!  **** input ****
!
!  rad          : radius of earth
!  omega        : angular velocity of earth
!  eigval       : vertical equivalent depth array
!  lev          : number of vertical levels
!  jtrun        : zonal wave number truncation limit
!
!  **** output ****
!
!  nw           : total wavenumber index array
!  a            : coefficient array of coefficient matrix
!  b            : coefficient array of coefficient matrix
!  c            : coefficient array of coefficient matrix
!  h            : coefficient array for (non)dimensionlize variables
!
!  modify to f90 bt C-H Lee and sort by River Chen in 2015
!------------------------------------------------------------------------
!
      use index
      use const, only : RTYPE

      implicit none
      integer  lev,jtrun,jtmax,nw(jtrun,jtmax),nnmivm
!byl      real     a(jtrun,jtrun,lev),b(jtrun,jtrun,lev),c(jtrun,jtrun,lev)
      real     a(jtrun,jtrun,nnmivm),b(jtrun,jtrun,nnmivm),c(jtrun,jtrun,nnmivm)
!byl      real     eigval(lev),h(jtrun,jtmax,lev)
      real     h(jtrun,jtmax,nnmivm)
      real(kind=RTYPE) eigval(lev)

      integer  k,l,m,n,m1,mf,mm,mn,nnp1
      real     omega,omega2,omga2r2,rad,rad2,tem,epxn

      rad2=rad*rad
      omega2=omega*omega
      omga2r2=omega2*rad2
!
!  compute total wavenumber index array
!
      mn=0
      do m=1,mlistnum
        mf=mlist(m)
        do n=mf,jtrun
          nw(n,m)=n*(n-1)
        enddo
      enddo
!
!  compute (non)dimensionlize array
!
!CWB2014 fixed the undefined value h in vartran loop
      h=0.
!byl      do k=1,lev
      do k=1,nnmivm
        tem=eigval(k)/omga2r2
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun
            if (n.ne.1) h(n,m,k)=(tem/nw(n,m))**0.5
          enddo
        enddo
      enddo
!
!  compute coefficients array of coefficient matrix
!
!byl      do k=1,lev
      do k=1,nnmivm
        a(1,1,k)=0.
        b(1,1,k)=0.
        c(1,1,k)=0.
        do m=1,jtrun
          m1=m-1
          mm=m
          if(m.eq.1)mm=2
          do l=mm,jtrun
            n=l-1
            nnp1=n*(n+1)
            epxn=((n*n-m1*m1)/(4.*n*n-1))**0.5
            c(m,l,k)=2.*m1/nnp1
            b(m,l,k)=(eigval(k)*nnp1/omga2r2)**0.5
            a(m,l,k)=-2.*epxn/n*((n+1)*(n-1))**0.5
          enddo
        enddo
      enddo
!
      return
      end
