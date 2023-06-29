      subroutine invmtx (a,na,v,nv,n,d,ip,ier)
!
!
! dimension of           a(na,n),v(nv,n),ip(2*n)
! arguments
!
!
! purpose                invmtx calculates the inverse of the input
!                        matrix a using gaussian elimination with full
!                        pivoting.
!
!
! usage                  call invmtx (a,na,v,nv,n,d,ip,ier)
!
! arguments
!
! on input               a
!                          a two-dimensional variable with row dimension
!                          na and column dimension n.  on input, a
!                          contains the elements of the matrix to be
!                          inverted.
!
!                        na
!                          an integer input variable set equal to the
!                          row dimension of a as declared in the
!                          dimension statement of the calling program.
!
!                        nv
!                          an integer input variable set equal to the
!                          row dimension of v as declared in the
!                          dimension statement of the calling program.
!
!                        n
!                          an integer input variable set equal to the
!                          column dimension of a.
!
!                        ip
!                          an integer array used internally for working
!                          storage.  it must have dimension at least
!                          2*n.
!
! on output              v
!                          a two-dimensional variable with row dimension
!                          nv and column dimension n.  on output, v
!                          contains the inverse of a.
!                               if a is entered twice in the parameter
!                               list, replacing v, then on output the
!                               array a will contain the inverse matrix,
!                               and the original a will be destroyed.
!
!                        d
!                          a real variable which on output contains the
!                          determinant of a.
!
!                        ier
!                          an integer error flag.
!                            = 33  if the matrix a is singular.
!                            =  0  otherwise.
!
! entry points           invmtx
!
! common blocks          none
!
! i/o                    if the
!                        matrix a has a zero pivot element (i.e., a is
!                        singular), the message
!                            *matrix singular in invmtx*
!                        is printed.
!
! precision              single
!
! language               fortran
!
!
!
      use rank
      use const, only : RTYPE

!     implicit real (a-h,o-z)
      implicit none

!
!     include '../include/rank.h'
!
!
      integer  na,nv,n,ier
      real(kind=RTYPE) a(na,n),v(nv,n)
      real     d

!CWB2015 orig
!     dimension ip(1)
      integer   ip(2*n)

      integer   i,j,k,l,m
      real      vmax,vh,pvt,hold

!
      ier = 0
!
! store a in v
!
      do 102 j=1,n
	 ip (j) = 0
	 do 101 i=1,n
	    v(i,j) = a(i,j)
  101    continue
  102 continue
      d = 1.
      do 111 m=1,n
	 vmax = 0.
!fj
         k = m
         l = m
!fj
	 do 107 j=1,n
	    if (ip(j)) 107,103,107
!
! find maximum pivot element
!
  103       do 106 i=1,n
		if (ip(i)) 106,104,106
  104          vh = abs(v(i,j))
		if (vmax-vh) 105,106,106
  105          vmax = vh
		k = i
		l = j
  106       continue
  107    continue
!
	 ip(l) = k
	 ip(n+m) = l
	 d = d*v(k,l)
	 pvt = v(k,l)
	 if (pvt .eq. 0.) go to 114
	 v(k,l) = 1.
	 do 108 j=1,n
	    hold = v(k,j)
	    v(k,j) = v(l,j)
	    v(l,j) = hold/pvt
  108    continue
	 do 110 i=1,n
	    if (i .eq. l) go to 110
	    hold = v(i,l)
	    v(i,l) = 0.
	    do 109 j=1,n
		v(i,j) = v(i,j)-v(l,j)*hold
  109       continue
  110    continue
  111 continue
!
! permute final inverse matrix
!
      do 113 j=1,n
	 m = n-j+1
	 l = ip(n+m)
	 k = ip(l)
	 if (k .eq. l) go to 113
	 d = -d
	 do 112 i=1,n
	    hold = v(i,l)
	    v(i,l) = v(i,k)
	    v(i,k) = hold
  112    continue
  113 continue
      return
  114 ier = 33
!fj
            if(myrank .eq. 0) print 888
!fj
  888 format(1x,'matrix singualr in invmtx')
      return
      end
