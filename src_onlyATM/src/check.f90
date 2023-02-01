!ch   subroutine check(pt,nx,my,my_max,lab)
      subroutine check(pt_2d,nx,my,my_max,lab)
!
!  purpose : write harmonic dia data to a file
!--------------------------------------------------------------------------
!  **** input ****
!  pt    : data array to be output
!  nx    : e-w direction dimension
!  my    : s-n direction dimension
!  lab   : label character
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!--------------------------------------------------------------------------
      use mpe
      use rank
      use index
      use fftcom
      use const, only: RTYPE

      implicit  none

      integer   nx,my,my_max,jj,latchk,i

!     real      holdit(2000),workx(2000), chold(2,1000)
!CWB20170801 fix for T853
      real      holdit(4000),workx(4000), chold(2,2000)
      real(kind=RTYPE) pt(nx,my_max)
      real(kind=RTYPE) pt_2d(nxp,my_max)

      equivalence (holdit,chold)
      character lab*10
      logical flag
!
      workx=0.

      flag = .false.
      latchk=60
!
      call mpe2d_unify_nx(pt,pt_2d)

      do 889 jj =1, jlistnum
      if(jlist1(jj) .eq. latchk) then
!
      do 888 i=1,nxdef(latchk)
      holdit(i)= pt(i,jj)
  888 continue
!
!     flag = .true.
      if(row_rank.eq.0) flag = .true.

      endif
 889  continue
!
      call mpe_broadcast(holdit,nx,flag,mpe_double)
!
      call rfftmlt(holdit,workx,trigsj(1,latchk),ifaxj(1,latchk),1, &
                   nx+2,nxdef(latchk),1,-1)
!
      if(myrank .eq. 0) write(*,887) lab,chold(1,7),chold(2,7)
  887 format(2x,a10,2x,f12.5,2x,f12.5)
!
      return
      end
