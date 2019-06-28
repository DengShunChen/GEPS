       subroutine get_phi(nxj,nx,lev,ptop,cp,r,g,pk,pk2,tt,qt,phii)

        integer nxj,nx,lev,i,k,kc
        real    ptop,cp,r,g,pk(nx,lev),pk2(nx,lev),         &
                tt(nx,lev),qt(nx,lev)
        real    pk2x(nx,lev),pkx(nx,lev),dhgtz(nx,lev),     &
                ppd,ppp,ppu,ttv,dhgt

        real    phii(nx,lev+1)
!
! geopotential height over model interface for WSM6
!
      do k=1,lev
      call vlog(pk2x(1,k),pk2(1,k),nxj)
      call vlog(pkx(1,k), pk(1,k), nxj)

      do i=1,nxj
        pk2x(i,k)=pk2x(i,k)*(cp/r)
        pkx(i,k) = pkx(i,k)*(cp/r)
      enddo
      call vexp(pk2x(1,k),pk2x(1,k),nxj)
      call vexp(pkx(1,k), pkx(1,k), nxj)
      enddo
!
      do 105 i = 1, nxj
      ppd = pk2x(i,1) * 1000.
      dhgt = ( ppd - ptop ) * 100. / g
      ppp = pkx(i,1) * 1000.
      ttv = tt(i,1)*(1.+0.608*qt(i,1))
!      ttv = tt(i,1)
      dhgtz(i,1)= dhgt * r * ttv / (100.*ppp)
  105 continue
!
      do 110 k=2, lev
      do 110 i=1, nxj
      ppu=pk2x(i,k-1) * 1000.
      ppd=pk2x(i,k) * 1000.
      dhgt = ( ppd - ppu ) * 100. / g
      ppp = pkx(i,k) * 1000.
      ttv = tt(i,k)*(1.+0.608*qt(i,k))
!      ttv = tt(i,k)
      dhgtz(i,k)= dhgt * r * ttv / (100.*ppp)
  110 continue
!
      do i=1,nxj
         phii(i,1)=0.
      enddo
       do k=lev,1,-1
         kc=lev-k+2
       do i = 1, nxj
          phii(i,kc)=phii(i,kc-1)+dhgtz(i,k)*g
      enddo
      enddo
!
      return
      end
