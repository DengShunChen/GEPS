      subroutine qprnt3 (a,t1,t2,ic,jc,kc,m,n,my_max,l,amul,add)
!
!  quick print a 3-d array a(m,l,n) at kc level starting at address
!  a(ic,kc,jc).  values are recomputed with add and amul, and then
!  are automatically scaled to allow 2 digit integer format printing.
!
!  ***input***
!
! a    : fwa of m x n array
! t1,t2: title (2*8 character)
! ic,jc: lower left corner coords to be printed
! kc   : level to be printed
! m,n,l: first,third,second dimension respectively
! amul : scaling multiplier for output
! add  : additive constant added before output
!
!  ***************************************************************
!
      use mpe
      use rank
      use index

      implicit  none
!
      real      a(m,l,my_max),amul,add
      character*8 t1, t2
      integer   ic,jc,kc,m,n,l,my_max

      integer   ix(41),jj,j,i,ie,jl,kp,kscale,j2,ii,jli
      real      haf,xm,af,x99,fk
!
      haf= 0.5
      ie = min0(ic+40,m)
      jl = n
!
      xm=0.
      do 110 jj =1,jlistnum
      j=jlist1(jj)
      if(j.ge.jc .and. j.le.jl) then
      do 100 i=ic,ie
      af= a(i,kc,jj)*amul + add
  100 xm= max(xm,abs(af))
      endif
  110 continue
!
      call mpe_global_max(xm,1,mpe_double)
!
      if(xm.lt.1.e-35) xm=99.
      x99= 99.
      xm = log10(x99/xm)
      kp = xm
      if (xm.lt.0.) kp = kp - 1
!
      kscale = - kp
      if(myrank .eq. 0) print 190, t1,t2,kscale,amul,add,(i,i=ic,ie,2)
  190 format(1h0,2a8,'  scale=',i3,' mult=',g12.5,' add=',g12.5, /1x,22i6)
!
      fk = 10.**kp
!
!     quickprint field
!
      do 200 j=jl,jc,-1
      do 205 jj =1,jlistnum
      j2=jlist1(jj)
      if(j2 .eq. j) then
        jli= jj
        ii = 0
        if (kp.eq.0) then
          do 300 i=ic,ie
          ii= ii+1
          af= a(i,kc,jli)*amul+add
          ix(ii) = af + sign(haf,af)
  300     continue
        else
          do 400 i=ic,ie
          ii= ii+1
          af= a(i,kc,jli)*amul+add
    7     ix(ii) = af*fk + sign(haf,af)
  400     continue
        endif
        call mpe_send_print(ix, 41, j, mpe_integer)
      endif
  205 continue
!
  200 continue

      if(myrank .eq. 0) then
      do j=jl,jc,-1
        call mpe_recv_print(ix, 41, j, mpe_integer)
        print 500, j, (ix(i),i=1,ii), j
  500   format (i4,42i3)
      enddo
      endif

      return
      end
