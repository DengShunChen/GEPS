      subroutine maxp ( mn,aa,amx,imx )
!
!     find out the max value and location of the array aa
!
      implicit  none
      integer   mn,imx,i
      real      aa(mn),amx
!
      imx = 1
      amx = aa(1)
      do 100 i = 2, mn
      if ( aa(i) .gt. amx )  then
         imx = i
         amx = aa(i)
      endif
  100 continue
!
!ccc  imx = ismax ( mn,aa,1 )
!ccc  amx = aa(imx)
!
      return
      end
!

!fj add

      subroutine maxpp2( nx,my,my_max,aa,amx,imx,jmx )
!     include '../include/index.h'
!     include '../include/mpe.h'
      use index
      use mpe

      implicit  none
      integer   nx,my,my_max,imx,jmx,jj,j,i,nxj,idummy
      real      aa(nx,my_max),amx
!
      amx = -1.797693134862316d+308
      imx=0 ; jmx=0
      do jj =1,jlistnum
      j=jlist1(jj)
!ch   nxj=nxdef(j)
      nxj=nxdef_2d(j)
      do i=1,nxj
        if (aa(i,jj).gt.amx) then
          amx = aa(i,jj)
!ch       imx = i
          imx = map2to1(i,j)
          jmx = j
        end if
      end do
      end do
!
      call mpe_global_maxloc(amx,2,imx,jmx,idummy,idummy,mpe_double)
!
      return
      end

      subroutine maxpp2l( lev,my,aa,amx,kmx,jmx )
!     include '../include/index.h'
!     include '../include/mpe.h'
      use index
      use mpe

      implicit  none
      integer   lev,my,kmx,jmx,k,j,idummy
      real      aa(my,lev),amx
!
      amx = -1.797693134862316d+308
      kmx=0 ; jmx=0
      do k =2,lev
      do j =1,my
        if (aa(j,k).gt.amx) then
          amx = aa(j,k)
          jmx = j
          kmx = k
        end if
      end do
      end do
!
!      call mpe_global_maxloc(amx,2,kmx,jmx,idummy,idummy,mpe_double)
!
      return
      end

!------------------------------------------------------------------
!ch for 2dMPI with array(nx_partial,my_partial)

      subroutine maxpp2_2d( aa,amx,imx,jmx )
      use index
      use param, only : my, my_max
      use mpe

      implicit  none
      integer   imx,jmx,jj,j,i,nxj,idummy
      real      aa(nxp,my_max),amx
 
      amx = -1.797693134862316d+308
      imx=0 ; jmx=0
      do jj =1,jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do i=1,nxj
        if (aa(i,jj).gt.amx) then
          amx = aa(i,jj)
!ch       imx = i
          imx = map2to1(i,j)
          jmx = j
        end if
      end do
      end do
 
      call mpe_global_maxloc(amx,2,imx,jmx,idummy,idummy,mpe_double)
 
      return
      end
