   subroutine read_mtnvar(nx,my,mtnv,hprime_b)
      use index
      use param ,only : jtrun,my_max
!
!      parameter (nx=3072, my=1536, mtnv=14)
!      parameter ( mtnv=14)
      implicit none
      integer i,j,v,jt,nrec,ii,jj,nxj,nx,my,mtnv
      real*4 hprime_a(nx,my,mtnv)
      real hprime_b(nxp,mtnv,my_max),hprime_a8(nx,my)
      real work(nxp,my_max)
!      real hprime_a(nx,my),hprime_aa(nx,my)
      character rfile*40
! 
!--------------------------------------------
!!      nrec=nx*my*4
      work = 0.

!xb118
!      if ( jtrun .gt. 999 ) then
!       write(rfile,100) jtrun,nx,my
!      else
!       if ( nx .gt. 999 .and. my .gt. 999 ) write(rfile,101) jtrun,nx,my
!       if ( nx .gt. 999 .and. my .le. 999 ) write(rfile,102) jtrun,nx,my
!       if ( nx .le. 999 .and. my .le. 999 ) write(rfile,103) jtrun,nx,my
!      endif
       if ( jtrun .gt. 999 ) then
        write(rfile,105) jtrun,nx,my
       else
        if ( nx .gt. 999 .and. my .gt. 999 ) write(rfile,106) jtrun,nx,my
        if ( nx .gt. 999 .and. my .le. 999 ) write(rfile,107) jtrun,nx,my
        if ( nx .le. 999 .and. my .le. 999 ) write(rfile,108) jtrun,nx,my
       endif
!100  format('global_mtnvar.t',i4.4,'.',i4.4,'.',i4.4,'.f77')
!101  format('global_mtnvar.t',i3.3,'.',i4.4,'.',i4.4,'.f77')
!102  format('global_mtnvar.t',i3.3,'.',i4.4,'.',i3.3,'.f77')
!103  format('global_mtnvar.t',i3.3,'.',i3.3,'.',i3.3,'.f77')
 105  format('global_mtnvarw.t',i4.4,'.',i4.4,'.',i4.4,'.f77')
 106  format('global_mtnvarw.t',i3.3,'.',i4.4,'.',i4.4,'.f77')
 107  format('global_mtnvarw.t',i3.3,'.',i4.4,'.',i3.3,'.f77')
 108  format('global_mtnvarw.t',i3.3,'.',i3.3,'.',i3.3,'.f77')
!xb118

      open(22,file=rfile,form='unformatted',status='old' ,convert='BIG_ENDIAN')
!!      open(22,file=rfile,form='unformatted',status='old'         &
!!          ,access='direct',recl=nrec )
!
!!       do v =1,mtnv
!         read(22,rec=v) hprime_a
      read(22) hprime_a
      do v =1,mtnv
        do j = 1, my
          jt = my - j + 1
          hprime_a8(:,j)=hprime_a(:,jt,v)
        enddo   
!      if( myrank .eq. 0 ) &
!      print*,' in read_mtnvar hprime_a = ',(hprime_a(1500,155,i),i=1,mtnv)

!--------------------------------------------
        call unify_reducepick(nx,my,my_max,hprime_a8,work)
        hprime_b(:,v,:) = work(:,:)
      enddo
!-- transpose even though glob 30" is from S to N and NCEP std is N to S


!      if( myrank .eq. 0 ) &
!      print*,'in read_mtnvar hprime_aa = ',(hprime_aa(1500,my-155+1,i),i=1,mtnv)
!
      close(22)
    
    end subroutine read_mtnvar
