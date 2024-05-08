      subroutine unify_reducepick(nx,my,my_max,ff,fp)
      use mpe
      use index
      ! pick the reduce-grid value from the nearest regular grid.
      ! this is only for changing low boundary data to reducegrid from
      ! regular grid.

      implicit none

      integer   nx,my,my_max
      integer   i,j,jj,nxj,ii
      real      fp(nxp,my_max),ff(nx,my)

      fp=0.

      do jj = 1, jlistnum
       j=jlist1(jj)
       ii=nxjstart(j)
       nxj=nxdef_2d(j)
       if( lreduce.eq.1 )call reducepickr (ff(1,j),nxdef(j),nx,1)
        do i = 1, nxj
          fp(i,jj) = ff(ii,j)
          ii=ii+1
        enddo
      enddo
  
!
      return
      end subroutine unify_reducepick
!
      subroutine unify_reducepicki(nx,my,my_max,ff,fp)
      use mpe
      use index
      ! pick the reduce-grid value from the nearest regular grid.
      ! this is only for changing low boundary data to reducegrid from
      ! regular grid.

      implicit none

      integer   nx,my,my_max
      integer   i,j,jj,nxj,ii
      integer   fp(nxp,my_max),ff(nx,my)

      fp=0.

      do jj = 1, jlistnum
       j=jlist1(jj)
       ii=nxjstart(j)
       nxj=nxdef_2d(j)
       if( lreduce.eq.1 )call reducepicki (ff(1,j),nxdef(j),nx,1)
        do i = 1, nxj
          fp(i,jj) = ff(ii,j)
          ii=ii+1
        enddo
      enddo
  
!
      return
      end subroutine unify_reducepicki
