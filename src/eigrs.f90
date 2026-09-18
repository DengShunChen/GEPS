      subroutine eigrs (rx,ir,idum,eval,evec,ix,wrk,ier)
      use rank
!
      dimension rx(ir,ir),eval(ir),evec(ir,ir),wrk(ir,ir)
!
!tom  double precision xx(500),eig(150),vec(150*150),work(500)
!cjh      dimension        xx(500),eig(150),vec(150*150),work(500)
!jh      dimension        xx(1000),eig(300),vec(300*300),work(1000)
!jh      dimension        xx(2000),eig(800),vec(800*800),work(2000)
!ch   dimension        xx(8000),eig(3200),vec(3200*3200),work(8000)
      dimension        xx(8000),eig(ir),vec(ir*ir),work(8000)
      ! #region agent log
      integer(kind=8) :: dbg0, dbg1, dbg2
      interface
        subroutine geps_dbg_log(hyp, locid, irank, p0, p1, p2)
          integer hyp, locid, irank
          integer(kind=8) p0, p1, p2
        end subroutine
      end interface
      dbg0 = loc(rx)
      dbg1 = ir
      dbg2 = loc(vec)
      call geps_dbg_log(14, 11, myrank, dbg0, dbg1, dbg2)
      ! #endregion
!
      if(ir.gt.2) then
!
      do 10 i=1,ir
      xx(i+2*ir)= dble(rx(i,i))
   10 continue
      do 20 i=2,ir
      xx(i+ir)= dble(rx(i-1,i))
   20 continue
      do 30 i=3,ir
      xx(i)= dble(rx(i-2,i))
   30 continue
      xx(1)= 0.d0
      xx(2)= 0.d0
      xx(1+ir)= 0.d0
!
!  find eigenvalues and eigenvectors for coefficient matrix
!
      ! #region agent log
      dbg0 = loc(vec)
      dbg1 = ir
      dbg2 = loc(xx)
      call geps_dbg_log(6, 29, myrank, dbg0, dbg1, dbg2)
      ! #endregion
      call rsb (ir,ir,3,xx,eig,1,vec,work,work(ir+1),ier)
      ! #region agent log
      dbg0 = loc(vec)
      dbg1 = ier
      dbg2 = ir
      call geps_dbg_log(6, 30, myrank, dbg0, dbg1, dbg2)
      ! #endregion
!
      do 40 i=1,ir*ir
      evec(i,1)= real(vec(i))
   40 continue
!
      do 50 i=1,ir
      eval(i)= real(eig(i))
   50 continue
!
      else
!
      do 35 i=1,ir
      eval(i)= 0.0
   35 continue
      do 37 i=1,ir*ir
      evec(i,1)= 0.0
   37 continue
      endif
!
      return
      end
