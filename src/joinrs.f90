      subroutine joinrs(cc,r1,r2,r3,r4,nx,my_max,lev,jlistnum,num,ncld)
!
      use const, only: RTYPE
!
      implicit  none
      real(kind=RTYPE) cc(*), r1(*), r2(*), r3(*), r4(*)
      integer   nx,my_max,lev,jlistnum,num,ncld
!
      if(num .eq. 1) call join1rs(cc,r1,nx,my_max,lev,jlistnum,ncld)
      if(num .eq. 2) call join2rs(cc,r1,r2,nx,my_max,lev,jlistnum,ncld)
      if(num .eq. 3) call join3rs(cc,r1,r2,r3,nx,my_max,lev,jlistnum,ncld)
!
      return
      end
!
      subroutine join1rs(cc,r1,nx,my_max,lev,jnum,ncld)
!
!     include '../include/index.h'
      use index
      use const, only: RTYPE

      implicit  none
      real(kind=RTYPE) cc(nx+2,levp,ncld,my_max)
      real(kind=RTYPE) r1(nxp,lev*ncld,my_max)
      real(kind=RTYPE) bufA(nxp,lev, ncld,my_max)
      real(kind=RTYPE) bufB(nx ,levp,ncld,my_max)
      integer   nx,my_max,lev,jnum,ncld
      integer   jj,j,nxj,k,i,n,nk,kk
!
!CWB 2019
      bufA=0.
      bufB=0.
!
      do jj =1, jlistnum
      do n=1,ncld
      nk=(n-1)*lev
      do k=1,lev
      kk=nk+k
      do i=1,nxp
         bufA(i,k,n,jj)=r1(i,kk,jj)
      enddo
      enddo
      enddo
      enddo
!
      call mpe2d_transpose_nxp_lev(bufA,bufB,nxp,nx,lev,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)

      do jj =1, jlistnum
      do n=1,ncld
      do k=1,levp
      do i=1,nx
         cc(i,k,n,jj)= bufB(i,k,n,jj)
      enddo
      enddo
      enddo
      enddo

      return
      end
!
      subroutine join2rs(cc,r1,r2,nx,my_max,lev,jnum,ncld)
!
!     include '../include/index.h'
      use index
      use const, only: RTYPE

      implicit  none
      real(kind=RTYPE) cc(nx+2,levp,1+ncld,my_max)
      real(kind=RTYPE) r1(nxp,lev,my_max)
      real(kind=RTYPE) r2(nxp,lev*ncld,my_max)
      real(kind=RTYPE) bufA(nxp,lev, 1+ncld,my_max)
      real(kind=RTYPE) bufB(nx ,levp,1+ncld,my_max)
      integer   nx,my_max,lev,jnum,ncld
      integer   jj,j,nxj,k,i,n,nk,kk
!
!CWB 2019
      bufA=0.
      bufB=0.
!
      do jj =1, jlistnum
      do k=1,lev
      do i=1,nxp
         bufA(i,k,1,jj)=r1(i,k,jj)
      enddo
      enddo
      enddo

      do jj =1, jlistnum
      do n=1,ncld
      nk=(n-1)*lev
      do k=1,lev
      kk=nk+k
      do i=1,nxp
         bufA(i,k,1+n,jj)=r2(i,kk,jj)
      enddo
      enddo
      enddo
      enddo

      call mpe2d_transpose_nxp_lev(bufA,bufB,nxp,nx,lev,levp,1+ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)

      do jj =1, jlistnum
      do n=1,1+ncld
      do k=1,levp
      do i=1,nx
         cc(i,k,n,jj)= bufB(i,k,n,jj)
      enddo
      enddo
      enddo
      enddo
!
      return
      end
!
      subroutine join3rs(cc,r1,r2,r3,nx,my_max,lev,jnum,ncld)
!
!     include '../include/index.h'
      use index
      use const, only: RTYPE

      implicit  none
      real(kind=RTYPE) cc(nx+2,levp,2+ncld,my_max)
      real(kind=RTYPE) r1(nxp,lev,my_max)
      real(kind=RTYPE) r2(nxp,lev,my_max)
      real(kind=RTYPE) r3(nxp,lev*ncld,my_max)
      real(kind=RTYPE) bufA(nxp,lev, 2+ncld,my_max)
      real(kind=RTYPE) bufB(nx ,levp,2+ncld,my_max)
      integer   nx,my_max,lev,jnum,ncld
      integer   jj,j,nxj,k,i,n,nk,kk
!
      bufA=0.
      bufB=0.
!
      do jj =1, jlistnum
      do k=1,lev
      do i=1,nxp
         bufA(i,k,1,jj)=r1(i,k,jj)
         bufA(i,k,2,jj)=r2(i,k,jj)
      enddo
      enddo
      enddo

      do jj =1, jlistnum
      do n=1,ncld
      nk=(n-1)*lev
      do k=1,lev
      kk=nk+k
      do i=1,nxp
         bufA(i,k,2+n,jj)=r3(i,kk,jj)
      enddo
      enddo
      enddo
      enddo

      call mpe2d_transpose_nxp_lev(bufA,bufB,nxp,nx,lev,levp,2+ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)

      do jj =1, jlistnum
      do n=1,2+ncld
      do k=1,levp
      do i=1,nx
         cc(i,k,n,jj)= bufB(i,k,n,jj)
      enddo
      enddo
      enddo
      enddo
!
      return
      end
!-----------------------------------------------------------------
!CWB2021 for single precision test

      subroutine joinrs_sp(cc,r1,r2,r3,r4,nx,my_max,lev,jlistnum,num,ncld)
!
      implicit  none
      real*4      cc(*), r1(*), r2(*), r3(*), r4(*)
      integer   nx,my_max,lev,jlistnum,num,ncld
!
      if(num .eq. 1) call join1rs_sp(cc,r1,nx,my_max,lev,jlistnum,ncld)
      if(num .eq. 2) call join2rs_sp(cc,r1,r2,nx,my_max,lev,jlistnum,ncld)
      if(num .eq. 3) call join3rs_sp(cc,r1,r2,r3,nx,my_max,lev,jlistnum,ncld)
!
      return
      end
!
      subroutine join1rs_sp(cc,r1,nx,my_max,lev,jnum,ncld)
!
!     include '../include/index.h'
      use index

      implicit  none
      real*4      cc(nx+2,levp,ncld,my_max)
      real*4      r1(nxp,lev*ncld,my_max)
      real*4      bufA(nxp,lev, ncld,my_max)
      real*4      bufB(nx ,levp,ncld,my_max)
      integer   nx,my_max,lev,jnum,ncld
      integer   jj,j,nxj,k,i,n,nk,kk
!
!CWB 2019
      bufB=0.
!
      do jj =1, jlistnum
      do n=1,ncld
      nk=(n-1)*lev
      do k=1,lev
      kk=nk+k
      do i=1,nxp
         bufA(i,k,n,jj)=r1(i,kk,jj)
      enddo
      enddo
      enddo
      enddo
!
      call mpe2d_transpose_nxp_lev_sp(bufA,bufB,nxp,nx,lev,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)

      do jj =1, jlistnum
      do n=1,ncld
      do k=1,levp
      do i=1,nx
         cc(i,k,n,jj)= bufB(i,k,n,jj)
      enddo
      enddo
      enddo
      enddo

      return
      end
!
      subroutine join2rs_sp(cc,r1,r2,nx,my_max,lev,jnum,ncld)
!
!     include '../include/index.h'
      use index

      implicit  none
      real*4      cc(nx+2,levp,1+ncld,my_max)
      real*4      r1(nxp,lev,my_max)
      real*4      r2(nxp,lev*ncld,my_max)
      real*4      bufA(nxp,lev, 1+ncld,my_max)
      real*4      bufB(nx ,levp,1+ncld,my_max)
      integer   nx,my_max,lev,jnum,ncld
      integer   jj,j,nxj,k,i,n,nk,kk
!
!CWB 2019
      bufB=0.
!
      do jj =1, jlistnum
      do k=1,lev
      do i=1,nxp
         bufA(i,k,1,jj)=r1(i,k,jj)
      enddo
      enddo
      enddo

      do jj =1, jlistnum
      do n=1,ncld
      nk=(n-1)*lev
      do k=1,lev
      kk=nk+k
      do i=1,nxp
         bufA(i,k,1+n,jj)=r2(i,kk,jj)
      enddo
      enddo
      enddo
      enddo

      call mpe2d_transpose_nxp_lev_sp(bufA,bufB,nxp,nx,lev,levp,1+ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)

      do jj =1, jlistnum
      do n=1,1+ncld
      do k=1,levp
      do i=1,nx
         cc(i,k,n,jj)= bufB(i,k,n,jj)
      enddo
      enddo
      enddo
      enddo
!
      return
      end
!
      subroutine join3rs_sp(cc,r1,r2,r3,nx,my_max,lev,jnum,ncld)
!
!     include '../include/index.h'
      use index

      implicit  none
      real*4      cc(nx+2,levp,2+ncld,my_max)
      real*4      r1(nxp,lev,my_max)
      real*4      r2(nxp,lev,my_max)
      real*4      r3(nxp,lev*ncld,my_max)
      real*4      bufA(nxp,lev, 2+ncld,my_max)
      real*4      bufB(nx ,levp,2+ncld,my_max)
      integer   nx,my_max,lev,jnum,ncld
      integer   jj,j,nxj,k,i,n,nk,kk
!
      do jj =1, jlistnum
      do k=1,lev
      do i=1,nxp
         bufA(i,k,1,jj)=r1(i,k,jj)
         bufA(i,k,2,jj)=r2(i,k,jj)
      enddo
      enddo
      enddo

      do jj =1, jlistnum
      do n=1,ncld
      nk=(n-1)*lev
      do k=1,lev
      kk=nk+k
      do i=1,nxp
         bufA(i,k,2+n,jj)=r3(i,kk,jj)
      enddo
      enddo
      enddo
      enddo

      call mpe2d_transpose_nxp_lev_sp(bufA,bufB,nxp,nx,lev,levp,2+ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)

      do jj =1, jlistnum
      do n=1,2+ncld
      do k=1,levp
      do i=1,nx
         cc(i,k,n,jj)= bufB(i,k,n,jj)
      enddo
      enddo
      enddo
      enddo
!
      return
      end
