      subroutine reduceintp(a,lonfd,lonf,latg)
!
! do  mass conserving interpolation from different grid at given latitude
!
! author: hann-ming henry juang 2008
!
      use const, only: RTYPE
!
      implicit none

!
      integer    lonf,latg,j,i,imp,lonfd(latg)
      real(kind=RTYPE) a(lonf,latg)
!
      real(kind=RTYPE) old(lonf),new(lonf)
      real(kind=RTYPE) xpast(lonf+1),xnext(lonf+1)
      real(kind=RTYPE) sc
      real      two_pi,dxp,dxf,hfdxp,hfdxf,pi
!

! ..................................
! ..................................
       pi  = 4.0 * atan(1.0)
       two_pi = 2.0 * pi
       dxf = two_pi / lonf
       hfdxf = 0.5 * dxf
!$omp parallel do                                       &
!$omp private(j,i,imp,dxp,hfdxp,xpast,xnext,sc,old,new) &
!$omp schedule(dynamic)
      do j=1,latg
       imp=lonfd(j)
       if( imp.ne.lonf ) then
        dxp = two_pi / imp
        hfdxp = 0.5 * dxp


        do i=1,imp+1
          xpast(i) = (i-1) * dxp - hfdxp
        enddo

        do i=1,lonf+1
          xnext(i) = (i-1) * dxf - hfdxf
        enddo

        sc=two_pi

        old(1:imp)=a(1:imp,j)
!CWB2021 for ndsl single precision test
       call cyclic_cell_ppm_intp(xpast,old,xnext,new,lonf,1,imp,lonf,sc)
!        call cyclic_cell_ppm_intp_dp(xpast,old,xnext,new,lonf,1,imp,lonf,sc)

!        call cyclic_cell_plm_intp(xpast,old,xnext,new,lonf,1,imp,lonf,sc)

        a(1:lonf,j)=new(1:lonf)
       endif
      enddo
!$omp end parallel do
! .................

      return
      end subroutine reduceintp

!helio>
! --------------------------------------
! inverse distance interpolation
! by helio
! --------------------------------------
      subroutine reduceintp_idw(a,lonfd,lonf,my,jj)
      use rank
      use index
      use phygrid ,only : outp,ls_full,ls_redu
      use const   ,only : RTYPE

      implicit none

      integer   i,j,k,h,jj,my
      integer   lonf,imp,lonfd !lonf=nx
      real(kind=RTYPE) a(lonf)
      real      old(lonf),new(lonf)

!for idw
      integer  al1,al2,al3,al4
      real     tmplist(4)
      real     dlist(4),axb(4)
      real     eab,ea,d

!for lat avg
      integer  lonr(my),dlonr(my/2),head,tail
      real     tmp,suma,m
      integer  z,fdum,rdum

!for checking each point
      integer  chk

!===================================
! lon points for each lat
      do i = 1,my/2
         dlonr(i) = lonf-4*(i-1)
      end do
      lonr(my/2+1:my) = dlonr

! reverse 1d array
      head = 1
      tail = my/2
      do while (head .lt. tail)
         tmp = dlonr(head)
         dlonr(head) = dlonr(tail)
         dlonr(tail) = tmp
         head = head+1
         tail = tail-1
      end do
! reverse 1d array end

      lonr(1:my/2) = dlonr

!===================================
      new(:)      = 0.
      tmplist(:)  = 0.
      dlist(:)    = 0.
      axb(:)      = 0.
      j           = jlist1(jj)


      do i = 1,lonf
        suma = 0.
        m    = 0.

        dlist(1)   = outp(i,jj,1)
        dlist(2)   = outp(i,jj,2)
        dlist(3)   = outp(i,jj,3)
        dlist(4)   = outp(i,jj,4)
        al1 = int(outp(i,jj,5)+0.00001)
        al2 = int(outp(i,jj,6)+0.00001)
        al3 = int(outp(i,jj,7)+0.00001)
        al4 = int(outp(i,jj,8)+0.00001)

        h=0
        if (al1.gt.0) then
           h=1
           tmplist(1)=a(al1)
        end if
        if (al2.gt.0) then
           h=2
           tmplist(2)=a(al2)
        end if
        if (al3.gt.0) then
           h=3
           tmplist(3)=a(al3)
        end if
        if (al4.gt.0) then
           h=4
           tmplist(4)=a(al4)
        end if

        if (h.eq.0) then
!if no value, average latitude with land sea mask
           fdum = ls_full(i,jj)

           do z = 1,lonr(j)
              rdum = ls_redu(z,jj)
!             if (j.eq.118) print*,'fdum',fdum,'rdum',rdum
              if (fdum.eq.0) then
                 if (my.eq.1280.and.i.eq.1802.and.j.eq.106) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.1803.and.j.eq.106) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.1179.and.j.eq.86) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.1180.and.j.eq.86) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.1181.and.j.eq.86) then
                    suma = suma+a(z)
                    m = m+1.
                 end if
                 if (rdum.eq.0) then
                    suma = suma+a(z)
                    m = m+1.
                 end if
              else if (fdum.eq.1)then
                 if (my.eq.768.and.j.eq.118) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.1190.and.j.eq.280) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.495.and.j.eq.294) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.271.and.j.eq.307) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.2129.and.j.eq.192) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.1236.and.j.eq.352) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.1048.and.j.eq.363) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.2184.and.j.eq.392) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.2130.and.j.eq.193) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.2162.and.j.eq.200) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.2181.and.j.eq.205) then
                    suma = suma+a(z)
                    m = m+1.
                 else if (my.eq.1280.and.i.eq.2182.and.j.eq.205) then
                    suma = suma+a(z)
                    m = m+1.
                 end if
                 if (rdum.eq.1) then
                    suma = suma+a(z)
                    m = m+1.
                 end if
              end if
           end do

           new(i) = suma/m

!          print*,'helio AVG i',i,'j',j,'ls=',fdum,'new = ',new(i)
        else if (h.gt.0) then
           fdum = ls_full(i,jj)
           axb(1:h) = dlist(1:h)*tmplist(h:1:-1)

           eab = 0.
           ea  = 0.
           do k = 1,h
             eab = eab+axb(k)
             ea  = ea +dlist(k)
           end do
           new(i) = eab/ea

           if (dlist(1).le.0) then
              new(i) = tmplist(1)
           end if


!          print*,'helio IDW i',i,'j',j,'ls=',fdum,'new = ',new(i)

        end if

      end do

!check missing data
!     do i = 1,lonf
!      if ((new(i).gt.0)) then
!      else
!      print*,'helio i = ',i,'j = ',j,'new = ',new(i),'ls =',ls_full(i,jj)
!      end if
!     end do


      a(1:lonf)=new(1:lonf)

      return
      end subroutine reduceintp_idw

!helio<


