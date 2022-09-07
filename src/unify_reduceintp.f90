      subroutine unify_reduceintp(nx,my,my_max,fp8,ff8)
      use mpe
      use index
      use const, only: RTYPE

      implicit none

      integer   nx,my,my_max
      integer   i,j,jj,nxj
      real*8    fp8(nxp,my_max),ff8(nx,my)
      real(kind=RTYPE) fp(nxp,my_max),ff(nx,my),ffx(nx,my_max)



      fp=fp8
      call mpe2d_unify_nx(ffx,fp)
      if( lreduce.eq.1 ) then
        do jj =1, jlistnum
          j=jlist1(jj)
          call reduceintp(ffx(1,jj),nxdef(j),nx,1)
        enddo
      endif
      call mpe2d_unify_my(ff,ffx)
      ff8=ff
!
      return
      end subroutine unify_reduceintp
