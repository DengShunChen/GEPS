      module index

      implicit none

      public

!      integer mlistnum,mlist(1000),nlist(1000),ilist(1000)
!      common /mwork/mlistnum,mlist,nlist,ilist

!      integer jlistnum,jlist1(2560),jlist2(2560),jlistnum_tmp(2560)
!      common /jwork/jlistnum,jlist1,jlist2,jlistnum_tmp
!! for Semi-Lagrangian
!      integer jlistnum_sl(200),jlist1_sl(2560,200)

!! reduceg
!      integer lreduce,nxdef(2560),mtrundef(2560),mfj(2560)
!      common /reduceindx/lreduce,nxdef,mtrundef,mfj

!      end module index

!CWB2017 for high resolution and high number of MPI processes
      integer mlistnum
      common /mwork/mlistnum

      integer, dimension(:),allocatable,save :: mlist,nlist,ilist

      integer jlistnum
      common /jwork/jlistnum

      integer, dimension(:),  allocatable,save :: jlist1,jlist2,jlistnum_sl
      integer, dimension(:,:),allocatable,save :: jlist2_2d,jlist1_sl

! reduceg
      integer lreduce
      common /reduceindx/lreduce

      integer, dimension(:),  allocatable,save :: nxdef,mtrundef

!2dMPI
      integer nsizex,nsizey,mrow,ncol,row_comm,col_comm,row_rank,col_rank
      common/mpi2d_comm/nsizex,nsizey,mrow,ncol,row_comm,col_comm,row_rank,col_rank

      integer Lstart,Lend,Llen,nxp,levp,jlen,nxf,levf,myf,Llistnum,      &
              Lstart_ncld,Lend_ncld,jtf,jtp,jtstart,jtend,jtlen
      common /iwork/Lstart,Lend,Llen,nxp,levp,jlen,nxf,levf,myf,Llistnum,&
              Lstart_ncld,Lend_ncld,jtf,jtp,jtstart,jtend,jtlen

      integer, dimension(:),allocatable,save :: Llist,Llist_ncld,jtlen_all,    &
                            nxjp,nxjstart,nxjend,nxjlen,nxdef_2d

      integer, dimension(:,:),allocatable,save :: nxjstart_all, &
                              nxjend_all,nxjlen_all,map2to1

      contains 

         subroutine allocate_index_array

           use param
           use rank

           integer  ierr

           allocate (mlist(jtmax),nlist(jtmax*nsize),ilist(jtmax*nsize), &
                     jlist1(my_max),jlist2(my),jlist2_2d(nsizex,my),     &
                     jlistnum_sl(nsize), nxdef(my),mtrundef(my),         &
                     Llist(lev),Llist_ncld(lev*ncld),         &
                     jtlen_all(nsizex),nxjp(my),   &
                     nxjstart(my),nxjend(my),nxjlen(my),nxdef_2d(my),    &
                     nxjstart_all(nsizex,my),nxjend_all(nsizex,my),      &
                     nxjlen_all(nsizex,my),map2to1(nxp,my),stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_index : allocate fail 1'
               stop
           end if

           allocate (jlist1_sl(my_max,nsize),stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_index : allocate fail 2'
               stop
           end if

           return

         end subroutine

         subroutine deallocate_index_array

           deallocate (mlist,nlist,ilist,jlist1,jlist2,jlist2_2d,jlistnum_sl, &
                       jlist1_sl,nxdef,mtrundef,Llist,Llist_ncld,      &
                       jtlen_all,nxjp,nxjstart,nxjend,nxjlen,   &
                       nxdef_2d,nxjstart_all,nxjend_all,&
                       nxjlen_all,map2to1)

           return

         end subroutine

      end module index
