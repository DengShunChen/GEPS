      subroutine lgndr (my2,jtrun,jtmax,sinl,poly,dpoly)
!
!  generate legendre polynomials and their derivatives on the
!  gaussian latitudes
!
! ***input***
!
!  my2: number of gaussian latitudes from south pole and equator
!  jtrun: zonal wavenumber truncation limit
!  jtmax: maximum amount of zonal wave located in each pe
!  sinl: sin of gaussian latitudes
!
!  ***output***
!
!  poly: associated legendre coefficients
!  dpoly: d(poly)/d(sinl)
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
! ******************************************************************
!
! ref= belousov, s. l., 1962= tables of normalized associated
!        legendre polynomials. pergamon press, new york
!
      use index
      use const, only: RTYPE
!
      implicit none

      integer  my2,jtrun,jtmax

      real(kind=RTYPE) poly(jtrun,my2,jtmax),dpoly(jtrun,my2,jtmax), &
                       sinl(my2)
!
      real     pnm(jtrun+1,jtrun+1),dpnm(jtrun+1,jtrun+1)
!
! sinl is sin(latitude) = cos(colatitude)
! pnm(np,mp) is legendre polynomial p(n,m) with np=n+1, mp=m+1
! pnm(mp,np+1) is x derivative of p(n,m) with np=n+1, mp=m+1
!
      integer  n,j,jtrunp,k,kp,m,mp,np,nps,mf,l,mlst
      real     xx,sn,sn2i,rt2,c1,c3,c4,c5,c6,c7,c8,cf,theta,fn,fn2,  &
               fn2s,ang,s1,s2,a,b,c,d,e,fk,fm,fm1,fm2,fm3,fms,fnp,fnp2
!
      jtrunp= jtrun+1
      pnm=0.0
      do 1001 j=1,my2
      xx= sinl(j)
      sn= sqrt(1.0-xx*xx)
	sn2i = 1.0/(1.0 - xx*xx)
      rt2= sqrt(2.0)
	c1 = rt2
!
      pnm(1,1) = 1.0/rt2
      theta=-atan(xx/sqrt(1.0-xx*xx))+2.0*atan(1.0)
!
      do 20 n=1,jtrun
	np = n + 1
      fn=n
	fn2 = fn + fn
	fn2s = fn2*fn2
! eq 22
      c1= c1*sqrt(1.0-1.0/fn2s)
      c3= c1/sqrt(fn*(fn+1.0))
	ang = fn*theta
	s1 = 0.0
	s2 = 0.0
	c4 = 1.0
	c5 = fn
	a = -1.0
	b = 0.0
!
      do 27 kp=1,np,2
	k = kp - 1
      s2= s2+c5*sin(ang)*c4
      if (k.eq.n) c4 = 0.5*c4
      s1= s1+c4*cos(ang)
	a = a + 2.0
	b = b + 1.0
      fk=k
	ang = theta*(fn - fk - 2.0)
	c4 = (a*(fn - b + 1.0)/(b*(fn2 - a)))*c4
	c5 = c5 - 2.0
   27 continue
! eq 19
	pnm(np,1) = s1*c1
! eq 21
	pnm(np,2) = s2*c3
   20 continue
!
      do 4 mp=3,jtrunp
	m = mp - 1
      fm= m
	fm1 = fm - 1.0
	fm2 = fm - 2.0
	fm3 = fm - 3.0
      c6= sqrt(1.0+1.0/(fm+fm))
! eq 23
	pnm(mp,mp) = c6*sn*pnm(m,m)
      if (mp - jtrunp) 3,4,4
    3 continue
	nps = mp + 1
!
      do 41 np=nps,jtrunp
	n = np - 1
      fn= n
	fn2 = fn + fn
	c7 = (fn2 + 1.0)/(fn2 - 1.0)
	c8 = (fm1 + fn)/((fm + fn)*(fm2 + fn))
      c= sqrt((fn2+1.0)*c8*(fm3+fn)/(fn2-3.0))
      d= -sqrt(c7*c8*(fn-fm1))
      e= sqrt(c7*(fn-fm)/(fn+fm))
! eq 17
	pnm(np,mp) = c*pnm(np-2,mp-2) &
                  + xx*(d*pnm(np-1,mp-2) + e*pnm(np - 1,mp))
   41 continue
    4 continue
!
      do 50 mp=1,jtrun
      fm= mp-1.0
	fms = fm*fm
      do 50 np=mp,jtrun
      fnp= np
	fnp2 = fnp + fnp
	cf = (fnp*fnp - fms)*(fnp2 - 1.0)/(fnp2 + 1.0)
      cf= sqrt(cf)
! der
      dpnm(np,mp)   = -sn2i*(cf*pnm(np+1,mp) - fnp*xx*pnm(np,mp))
   50 continue
!
      do 71 m=1,mlistnum
      mf=mlist(m)
      do 71 l=mf,jtrun
!ibm--- poly & dpoly: 2nd & 3rd dimension is transposed
!ibm     poly(l,m,j)= pnm(l,mf)
!ibm     dpoly(l,m,j)=dpnm(l,mf)
      poly(l,j,m)= pnm(l,mf)
      dpoly(l,j,m)=dpnm(l,mf)
   71 continue
      mlst=ilist(1)
!ibm     if(mlst .ne. 0) dpoly(1,mlst,j)= 0.0
      if(mlst .ne. 0) dpoly(1,j,mlst)= 0.0
!
 1001 continue
!
      return
      end
