      subroutine unify_reduceintp(nx,my,my_max,fp,ff)
      use mpe
      use index

      implicit none

      integer   nx,my,my_max
      integer   i,j,jj,nxj
      real      fp(nxp,my_max),ff(nx,my)

      do jj =1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
         ff(i,j)=fp(i,jj)
        enddo
      enddo

      call mpe_unify(ff,nx,my,2,mpe_double)

      if( lreduce.eq.1 ) then
        do jj =1, jlistnum
          j=jlist1(jj)
          call reduceintp(ff(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(ff,nx,my,5,mpe_double)
      endif
!
      return
      end subroutine unify_reduceintp
