      subroutine unify_reduceintp_gpu(nx,my,my_max,fp,ff)
      use mpe
      use index
      use const, only: RTYPE
      use mod_ndslfv_monoadv_gpu, only: cyclic_cell_ppm_intp_gpu

      implicit none

      integer   nx,my,my_max
      integer   i,j,jj,nxj
      real(kind=RTYPE) fp(nxp,my_max),ff(nx,my),ffx(nx,my_max)

      !variable on gpu for inline code
      real pi,dxp,dxf,hfdxp,hfdxf
      real(kind=RTYPE) xpast(nx+1,my_max),  old(nx,my_max), &
                       xnext(nx+1,my_max),  new(nx,my_max)
      real(kind=RTYPE) two_pi
      integer :: imp, imlst_inp(my_max), imlst_out(my_max)
      integer,parameter:: async_id = 1
      integer ierr
      pi  = 4.0 * atan(1.0)
      two_pi = 2.0 * pi
!$acc wait(async_id)
!$acc enter data create(ffx) async(async_id) 
!$acc enter data create(xpast,xnext,old, new ) async(async_id) 
!$acc enter data create(imlst_inp,imlst_out) async(async_id) 


!$acc wait(async_id)
      call mpe2d_unify_nx_gpu(ffx,fp)
!$acc wait(async_id)


      if( lreduce.eq.1 ) then
!$acc parallel loop private(j,imp ) async(async_id)
        do jj =1, jlistnum
          j=jlist1(jj)
          imp=nxdef(j)
          imlst_inp(jj) = imp
          imlst_out(jj) =  nx
        enddo

!$acc  parallel loop collapse(2) &
!$acc& private(j,imp,dxp ,hfdxp ,dxf ,hfdxf ) async(async_id)
        do jj =1, jlistnum
          do i=1,nx+1
            j=jlist1(jj)
            imp=nxdef(j)
            if (i <= imp+1 )then
              dxp = two_pi / imp
              hfdxp = 0.5 * dxp
              xpast(i,jj) = (i-1) * dxp - hfdxp
            endif
            dxf = two_pi / nx
            hfdxf = 0.5 * dxf
            xnext(i,jj) = (i-1) * dxf - hfdxf
            if(i<=imp)then
             old(i,jj)=ffx(i,jj)
            endif
          enddo
        enddo

        !$acc wait(async_id)
        call cyclic_cell_ppm_intp_gpu(xpast   ,old  &
                                     ,xnext   ,new  &
                                     ,imlst_inp ,imlst_out      &
                                     ,jlistnum ,nx, my_max, 1, 1, two_pi &
                                     ,async_id)

        !enddo

        !$acc parallel loop collapse(2)  async(async_id)
        do jj =1, jlistnum
          do i=1,nx
          ffx(i,jj)=new(i,jj)
          enddo
        enddo

      endif


!$acc wait(async_id)
      call mpe2d_unify_my_gpu (ff,ffx)
!
!$acc wait(async_id)
!!!!!!$acc update self(ff) async(async_id)
!$acc exit data delete( ffx ) async(async_id)
!$acc exit data delete(imlst_inp,imlst_out) async(async_id)
!$acc exit data delete( xpast,xnext,old,new ) async(async_id)
!$acc wait(async_id)


      return
      end subroutine unify_reduceintp_gpu

!helio>
      subroutine unify_reduceintp_idw_gpu(nx,my,my_max,fp,ff)
      use mpe
      use index
      use const, only: RTYPE

      implicit none

      integer   nx,my,my_max
      integer   i,j,jj,nxj
      real(kind=RTYPE) fp(nxp,my_max),ff(nx,my),ffx(nx,my_max)
      integer,parameter:: async_id = 1
!$acc wait(async_id)
!$acc update self(fp) async(async_id) 
!$acc wait(async_id)

      call mpe2d_unify_nx(ffx,fp)
      if( lreduce.eq.1 ) then
        do jj =1, jlistnum
          j=jlist1(jj)
          call reduceintp_idw(ffx(1,jj),nxdef(j),nx,my,jj)
        enddo
      endif
      call mpe2d_unify_my(ff,ffx)
!
!$acc wait(async_id)
!$acc update device(ff) async(async_id) 
!$acc wait(async_id)
      return
      end subroutine unify_reduceintp_idw_gpu
!helio<
