<<<<<<< HEAD
     subroutine read_mtnvar(nx,my,mtnv,hprime_b)
=======
      subroutine read_mtnvar(nx,my,mtnv,hprime_b)
>>>>>>> 64a9f40a91bbbba7e4796b7967f37f7109a84f91

!
      use index
      use param ,only : my_max
!
!      parameter (nx=3072, my=1536, mtnv=14)
!      parameter ( mtnv=14)
      implicit none
      integer i,j,v,jt,nrec,ii,jj,nxj,nx,my,mtnv
      real*4 hprime_a(nx,my)
      real hprime_b(nxp,mtnv,my_max),hprime_a8(nx,my)
!      real hprime_a(nx,my),hprime_aa(nx,my)
! 
!--------------------------------------------
      nrec=nx*my*4
<<<<<<< HEAD
!!
!open(22,file='/nwpr/gfs/p037/data/NCEPGFSGWDlzl/lzlwrkf2/terrain/'  &
!!      //'terr_640.2576.1280/global_mtnvar.t640.2576.1280.f77',
!&
=======
!!      open(22,file='/nwpr/gfs/p037/data/NCEPGFSGWDlzl/lzlwrkf2/terrain/'  &
!!      //'terr_640.2576.1280/global_mtnvar.t640.2576.1280.f77',           &
>>>>>>> 64a9f40a91bbbba7e4796b7967f37f7109a84f91
!!        form='unformatted',status='old',access='direct',recl=nrec )
      open(22,file='mtnvar.dat',                                      &
        form='unformatted',status='old',access='direct',recl=nrec )
!     &   form='binary',status='old' )
       do v =1,mtnv
         read(22,rec=v) hprime_a
         hprime_a8=hprime_a
!      if( myrank .eq. 0 ) &
<<<<<<< HEAD
!      print*,' in read_mtnvar hprime_a =
!      ',(hprime_a(1500,155,i),i=1,mtnv)
=======
!      print*,' in read_mtnvar hprime_a = ',(hprime_a(1500,155,i),i=1,mtnv)
>>>>>>> 64a9f40a91bbbba7e4796b7967f37f7109a84f91

!--------------------------------------------

        do jj = 1, jlistnum
          j=jlist1(jj)
          jt = my - j + 1
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
<<<<<<< HEAD
          if( lreduce.eq.1 ) &
             call reducepick(hprime_a8(1,jt),nxdef(j),nx,1)        
=======
          if( lreduce.eq.1 )call reducepick (hprime_a8(1,jt),nxdef(j),nx,1)        
>>>>>>> 64a9f40a91bbbba7e4796b7967f37f7109a84f91
          do i = 1, nxj
            hprime_b(i,v,jj) = hprime_a8(ii,jt)
            ii=ii+1
          enddo
        enddo
      enddo

!-- transpose even though glob 30" is from S to N and NCEP std is N to S


!      if( myrank .eq. 0 ) &
<<<<<<< HEAD
!      print*,'in read_mtnvar hprime_aa =
!      ',(hprime_aa(1500,my-155+1,i),i=1,mtnv)
=======
!      print*,'in read_mtnvar hprime_aa = ',(hprime_aa(1500,my-155+1,i),i=1,mtnv)
>>>>>>> 64a9f40a91bbbba7e4796b7967f37f7109a84f91
!
      return
      end
