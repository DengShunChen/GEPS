      subroutine read_mtnvar(nx,my,mtnv,hprime_b)

!
      use index
      use param ,only : nco,my_max
!
!      parameter (nx=3072, my=1536, mtnv=14)
!      parameter ( mtnv=14)
      implicit none
      integer i,j,v,jt,nrec,ii,jj,nxj,nx,my,mtnv
      real*4 hprime_a(nx,my)
      real hprime_b(nxp,mtnv,my_max),hprime_a8(nx,my)
!      real hprime_a(nx,my),hprime_aa(nx,my)
      character rfile*40
! 
!--------------------------------------------
      nrec=nx*my*4

      if ( nco .gt. 999 ) then
        write(rfile,100) nco,nx,my
      else
        if ( nx .gt. 999 .and. my .gt. 999 ) write(rfile,101) nco,nx,my
        if ( nx .gt. 999 .and. my .le. 999 ) write(rfile,102) nco,nx,my
        if ( nx .le. 999 .and. my .le. 999 ) write(rfile,103) nco,nx,my
      endif
 100  format('global_mtnvar.t',i4.4,'.',i4.4,'.',i4.4,'.f77')
 101  format('global_mtnvar.t',i3.3,'.',i4.4,'.',i4.4,'.f77')
 102  format('global_mtnvar.t',i3.3,'.',i4.4,'.',i3.3,'.f77')
 103  format('global_mtnvar.t',i3.3,'.',i3.3,'.',i3.3,'.f77')

      open(22,file=rfile,form='unformatted',status='old'         &
          ,access='direct',recl=nrec )
!
       do v =1,mtnv
         read(22,rec=v) hprime_a
         hprime_a8=hprime_a
!      if( myrank .eq. 0 ) &
!      print*,' in read_mtnvar hprime_a = ',(hprime_a(1500,155,i),i=1,mtnv)

!--------------------------------------------

        do jj = 1, jlistnum
          j=jlist1(jj)
          jt = my - j + 1
          ii=nxjstart(jt)
          nxj=nxdef_2d(jt)
          if( lreduce.eq.1 ) &
             call reducepick(hprime_a8(1,jt),nxdef(jt),nx,1)        
          do i = 1, nxj
            hprime_b(i,v,jj) = hprime_a8(ii,jt)
            ii=ii+1
          enddo
        enddo
      enddo

!-- transpose even though glob 30" is from S to N and NCEP std is N to S


!      if( myrank .eq. 0 ) &
!      print*,'in read_mtnvar hprime_aa = ',(hprime_aa(1500,my-155+1,i),i=1,mtnv)
!
      return
      end
