      subroutine pnmy (jtrun,sinlj,pnm)
!
!  generate legendre polynomials and their derivatives on the
!  gaussian latitudes
!
! ***input***
!
!  jtrun: zonal wavenumber truncation limit
!  sinlj: sin of gaussian latitudes
!
!  ***output***
!
! ******************************************************************
!
!
      use const, only: RTYPE
!
      implicit  none
!
      integer   jtrun
      real(kind=RTYPE) sinlj
      real(kind=RTYPE) pnm(jtrun+1,jtrun+1)
!
! sinlj is sin(latitude) = cos(colatitude)
! pnm(np,mp) is legendre polynomial p(n,m) with np=n+1, mp=m+1
! pnm(mp,np+1) is x derivative of p(n,m) with np=n+1, mp=m+1
!
      real      xx,sn,sn2i,rt2,c1,theta,fn,fn2,fn2s,c3,ang
      real      s1,s2,c4,c5,c6,c7,c8,a,b,c,d,e,fk
      real      fm,fm1,fm2,fm3
      integer   jtrunp,n,np,kp,k,mp,m,nps

      jtrunp= jtrun+1
      pnm=0.0
!  
      xx= sinlj
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
	pnm(np,mp) = c*pnm(np-2,mp-2) + xx*(d*pnm(np-1,mp-2) + e*pnm(np - 1,mp))
   41 continue
    4 continue
!
      return
      end
