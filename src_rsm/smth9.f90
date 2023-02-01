!-----------------------------------------------------------
! 9-point smoother
!-----------------------------------------------------------
      subroutine smth9(im,jm,data,sdata,ismth)
      use const, only: RTYPE
!
      implicit   none
      integer                   :: im,jm,ismth
      real(kind=RTYPE),dimension(im,jm) :: data,sdata

      integer                   ::i,j,ip1,im1,jp1,jm1

      do j=1,jm
        jp1=j+1
        jm1=j-1
      do i=1,im
        ip1=i+1
        im1=i-1
        if(ip1 > im)ip1=ip1-im
        if(im1 < 1 )im1=im1+im

        if(j==1  .or. j==jm)then
! boundary by 3-point smoothing
          if(ismth==1)then
           sdata(i,j)=data(i,j)+(data(ip1,j)-2.*data(i,j)+data(im1,j))/4.
          else
           sdata(i,j)=data(i,j)
          endif
        else
! by 9-point smoothing
         sdata(i,j)=data(i,j)+(data(i,jp1)+data(i,jm1)+data(ip1,j)      &
                   +data(im1,j)-4.*data(i,j))/8.+(data(ip1,jp1)         &
                   +data(ip1,jm1)+data(im1,jp1)+data(im1,jm1)           &
                   -4.*data(i,j))/16.
        endif
      enddo
      enddo
      return
      end
