      subroutine make_list
!
!  m        : start point of Lugendre number
!  mlist    : index of Fourier number on PE
!  mlistnum : index of Fourier number on PE
!
!
      use param
      use index
      use mpe
      use rank
      use spec,only : jtwv,jtwvp

      implicit  none
      integer   mm,m,ipe,n,jj,j,j1,j2,ii,i1,i2

      do mm=1,jtmax
         mlist(mm)=0
      enddo

      mlistnum=0
      do m=1,jtmax*nsize
         nlist(m)=0 ; ilist(m)=0
      enddo

!
! for spectral space, use folded cyclic allocation type
!

! 1dMPI
!     m=jtrun
!     do mm=1,jtmax+1,2
!        do ipe=nsize,1,-1
!           if(m.gt.0)then
!            if(ipe-1 .eq. myrank) then
!              mlistnum=mlistnum+1
!              mlist(mm)=m
!            endif
!              n=(ipe-1)*jtmax+mm
!              nlist(m)=n
!              if(ipe-1 .eq. myrank) ilist(m)=mm
!              m=m-1
!           endif
!        enddo
!        do ipe=1,nsize
!           if(m.gt.0)then
!            if(ipe-1 .eq. myrank) then
!              mlistnum=mlistnum+1
!              mlist(mm+1)=m
!            endif
!              n=(ipe-1)*jtmax+mm+1
!              nlist(m)=n
!              if(ipe-1 .eq. myrank) ilist(m)=mm+1
!              m=m-1
!           endif
!        enddo
!     enddo

      jlistnum=0
! for Semi-Lagrangian
!byl      do ipe=1,nsize
      do ipe=1,nsizey
         jlistnum_sl(ipe)=0
      enddo

      do jj=1,my_max
         jlist1(jj)=0
! for Semi-Lagrangian
!byl       do ipe=1,nsize
       do ipe=1,nsizey
          jlist1_sl(jj,ipe)=0
       enddo
      enddo

      do j=1,my
         jlist2(j)=0
      enddo

!
! for grid space, use cyclic allocation type
!

! 1dMPI
!     j=1
!     do jj=1,my_max
!        do ipe=1,nsize
!           if(j.le.my)then
!            if(ipe-1 .eq. myrank) then
!              jlistnum=jlistnum+1
!              jlist1(jj)=j
!            endif
!              j2=(ipe-1)*my_max+jj
!              jlist2(j)=j2
!              j=j+1
!           endif
!        enddo
!     enddo
!


! for Semi-Lagrangian
      j=1
      do jj=1,my_max
!ch      do ipe=1,nsize
         do ipe=1,nsizey
            if(j.le.my)then
               jlistnum_sl(ipe)=jlistnum_sl(ipe)+1
               jlist1_sl(jj,ipe)=j
               j=j+1
            endif
         enddo
      enddo
!

!for 2dMPI >>

      jlist2_2d=0

      j=1
      do jj=1,my_max
         do ipe=1,nsizey
            if(j.le.my)then
             if(ipe-1 .eq. col_rank) then
               jlistnum=jlistnum+1
               jlist1(jj)=j
             endif
               do ii=1,nsizex
                  j2=(ipe-1)*my_max*nsizex+(ii-1)*my_max+jj
                  jlist2_2d(ii,j)=j2
               enddo
               j2=(ipe-1)*my_max+jj
               jlist2(j)=j2
               j=j+1
            endif
         enddo
      enddo

      Llistnum=0

      do n=1,lev
         if(n .ge. Lstart .and. n .le. Lend) then
           Llistnum=Llistnum+1
           Llist(Llistnum)=n
         endif
      enddo

      do n=1,ncld
         i1=(n-1)*levp
         i2=(n-1)*lev
      do ii=1,levp
         mm=i1+ii
         jj=i2+Llist(ii)
         Llist_ncld(mm)=jj
      enddo
      enddo

      m=jtrun
      do mm=1,jtmax+1,2
         do ipe=nsizey,1,-1
            if(m.gt.0)then
             if(ipe-1 .eq. col_rank) then
               mlistnum=mlistnum+1
               mlist(mm)=m
             endif
               n=(ipe-1)*jtmax+mm
               nlist(m)=n
               if(ipe-1 .eq. col_rank) ilist(m)=mm
               m=m-1
            endif
         enddo
         do ipe=1,nsizey
            if(m.gt.0)then
             if(ipe-1 .eq. col_rank) then
               mlistnum=mlistnum+1
               mlist(mm+1)=m
             endif
               n=(ipe-1)*jtmax+mm+1
               nlist(m)=n
               if(ipe-1 .eq. col_rank) ilist(m)=mm+1
               m=m-1
            endif
         enddo
      enddo


! for siimpl

        i1=1
        jtf=0
        n=0

        do m=1,mlistnum
           mm=mlist(m)

           do j2=mm,jtrun
              n=n+1
              jtwv(n)=j2
           enddo

          i2=jtrun-mm+1
          jtf=jtf+i2

          i1=i1+i2

        enddo

        j1=jtf/nsizex
        j2=mod(jtf,nsizex)
        if(j2.eq.0)then
          jtp=j1
        else
!         jtp=(jtf/nsizex)+1
          jtp=(jtf-j2)/nsizex+1
        endif

!       jtstart=row_rank*j1+1+min(row_rank,j2)
!       jtend=jtstart+j1-1
!       if(j2> row_rank)jtend=jtend+1
!       jtlen=jtend-jtstart+1

        if(j2.eq.0)then
          jtstart=row_rank*jtp+1
          jtend=jtstart+jtp-1
        else
          jtstart=row_rank*jtp+1
          jtend=jtstart+jtp-1
          if(row_rank .eq. (nsizex-1)) jtend=jtend-(nsizex-j2)
        endif
        jtlen=jtend-jtstart+1

        j=0
        do m=jtstart,jtend
           j=j+1
           jtwvp(j)=jtwv(m)
        enddo

!       call get_all_jtlen   ! gather all jtlen of a row

        call mpe2d_get_allrow(1,jtlen,jtlen_all)   ! gather all jtlen of a row

!       print *,'total_wave jtf jtp jtstart jtend jtlen=',n,jtf,jtp,jtstart,jtend,jtlen
!       print *,'           jtf jtp jtstart jtend jtlen=',  jtf,jtp,jtstart,jtend,jtlen
        do j=jtstart,jtend
!          print *,j,jtwv(j)
        enddo
!ch <<

!        do ipe=1,nsizey
!              print *,ipe,jlistnum_sl(ipe)
!        enddo
!        print *,'   '
!        do ipe=1,nsizey
!        do jj=1,my_max
!              print *,col_rank,ipe,jj,jlist1_sl(jj,ipe)
!        enddo
!        enddo

!  Usage:
!        do ipe=1,npe
!        do mm=1,mlistnum(ipe)
!           mf=mlist(mm,ipe); m=nlist(mf)
!           ...
!        enddo
!        enddo

      return
      end

!-------------------------------------------------------
      subroutine make_list_nx

!
!  decompose nx for 2dMPI
!
      use param
      use mpe
      use rank
      use index

      implicit  none
      integer   i,ii,j,n1,n2,nxj,tmp(my,nsizex),lat

! nxp  : nx  partial
! nxjp : nxj partial

! nxjstart : nxj start
! nxjend   : nxj end

! map2to1  : 2d local index to 1d global index

         do j=1,my
            nxj=nxdef(j)

            n1=nxj/nsizex
            n2=mod(nxj,nsizex)

            nxjstart(j)=row_rank*n1+1+min(row_rank,n2)
            nxjend(j)=nxjstart(j)+n1-1
            if(row_rank < n2)nxjend(j)=nxjend(j)+1

            nxjlen(j)=nxjend(j)-nxjstart(j)+1
            nxjp(j)=nxjlen(j)
            nxdef_2d(j)=nxjlen(j)

            ii=nxjstart(j)
!ch         do i=1,nxjlen(j)
            do i=1,nxp
               map2to1(i,j)=ii
               ii=ii+1
            enddo
         enddo

      call mpe2d_get_allrow(my,nxjstart,tmp)  ! gather all nxjstart of a row
      nxjstart_all=transpose(tmp)
      call mpe2d_get_allrow(my,nxjend,tmp)    ! gather all nxjend of a row
      nxjend_all=transpose(tmp)
      call mpe2d_get_allrow(my,nxjlen,tmp)    ! gather all nxjlen of a row
      nxjlen_all=transpose(tmp)

      do j=1,my
         nxj=nxdef(j)
!        print 101, j,nxj,nxjstart(j),nxjend(j),nxjlen(j)
         do i=1,nsizex
!        print 102, j,i,nxjstart_all(i,j),nxjend_all(i,j),nxjlen_all(i,j)
         enddo
      enddo

      nxjp_acc = 0
      nxjp_acc(1) = 1
      do j = 1,jlistnum
         lat = jlist1(j)
         nxjp_acc(j+1) = nxjp_acc(j) + nxjp(lat)
      end do
      nxptot = nxjp_acc(jlistnum+1)-1

!     stop

101   format('j,nxj,nxjstart,nxjend,nxjlen=',5I6)
102   format('j,i,nxjstart,nxjend,nxjlen=',5I6)

      return
      end
