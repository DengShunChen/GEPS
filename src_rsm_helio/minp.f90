      subroutine minp ( mn,aa,amn,imn )
!
!     find out the min value and location of the array aa
!
      implicit  none
      integer   mn,imn,i
      real      aa(mn),amn
!
      imn = 1
      amn = aa(1)
      do 100 i = 2, mn
      if ( aa(i) .lt. amn )  then
         imn = i
         amn = aa(i)
      endif
  100 continue
!
!ccc  imn = ismin ( mn,aa,1 )
!ccc  amn = aa(imn)
!
      return
      end

!fj add

      subroutine minpp2( nx,my,my_max,aa,amn,imn,jmn )
!     include '../include/index.h'
!     include '../include/mpe.h'
      use index
      use mpe

      implicit  none
      integer   nx,my,imn,jmn,jj,j,i,idummy,my_max,nxj
      real      aa(nx,my_max),amn
!
      amn = 1.797693134862316d+308
      imn=0 ; jmn=0
      do jj =1,jlistnum
      j=jlist1(jj)
      nxj=nxdef(j)
      do i=1,nxj
        if (aa(i,jj).lt.amn) then
          amn = aa(i,jj)
          imn = i
          jmn = j
        end if
      end do
      end do
!
      call mpe_global_minloc(amn,2,imn,jmn,idummy,idummy,mpe_double)
!
      return
      end

      subroutine minpp2l( lev,my,aa,amn,kmn,jmn )
!     include '../include/index.h'
!     include '../include/mpe.h'
      use index
      use mpe

      implicit  none
      integer   lev,my,kmn,jmn,k,j
      real      aa(my,lev),amn
!
      amn = 1.797693134862316d+308
      kmn=0 ; jmn=0
      do k =2,lev
      do j =1,my
        if (aa(j,k).lt.amn) then
          amn = aa(j,k)
          kmn = k
          jmn = j
        end if
      end do
      end do
!
!      call mpe_global_minloc(amn,2,kmn,jmn,idummy,idummy,mpe_double)
!
      return
      end
!------------------------------------------------------------------
!ch for 2dMPI with array(nx_partial,my_partial)

      subroutine minpp2_2d( aa,amn,imn,jmn )
      use index
      use param, only : my, my_max
      use mpe

      implicit  none
      integer   imn,jmn,jj,j,i,nxj,idummy
      real      aa(nxp,my_max),amn
 
      amn = 1.797693134862316d+308
      imn=0 ; jmn=0
      do jj =1,jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do i=1,nxj
        if (aa(i,jj).lt.amn) then
          amn = aa(i,jj)
!ch       imn = i
          imn = map2to1(i,j)
          jmn = j
        end if
      end do
      end do
 
      call mpe_global_minloc(amn,2,imn,jmn,idummy,idummy,mpe_double)
 
      return
      end
