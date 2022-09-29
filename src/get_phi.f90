       subroutine get_phi(nxj,nx,lev,ptop,cp,r,g,sgeo,pk,pk2,tt,qt,  &
                          phii,phi)
!
        use const, only : RTYPE
!

        integer nxj,nx,lev,i,k,kc
        real    ptop,cp,r,g,pkd(nx),pk2d(nx)
        real    pk2x(nx,lev),pkx(nx,lev),dhgtz(nx,lev),     &
                ppd,ppp,ppu,ttv,dhgt,theda(nx,lev)

        real    phii(nx,lev+1)
        real(kind=RTYPE) tt(nx,lev),qt(nx,lev),phi(nx,lev), &
                         sgeo(nx),pk(nx,lev),pk2(nx,lev)
!
! geopotential height at model interface
!
      do k=1,lev
      pkd(:) =pk(:,k)
      pk2d(:)=pk2(:,k)
      call vlog(pk2x(1,k),pk2d,nxj)
      call vlog(pkx(1,k), pkd, nxj)

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
!
! geopotential height at model layer
!
      do k = 1, lev
        do i = 1, nxj
          theda(i,k) = tt(i,k)*(1.0+0.608*qt(i,k))/pk(i,k)
        enddo
      enddo
!
      do i = 1 ,nxj
        phi(i,lev) = sgeo(i)+cp*theda(i,lev)*(pk2(i,lev)-pk(i,lev))
      enddo
!
      do k = lev-1, 1, -1
        do i = 1,nxj
          phi(i,k) = phi(i,k+1) +cp*(theda(i,k)*(pk2(i,k)-pk(i,k)) &
                    +theda(i,k+1)*(pk(i,k+1)-pk2(i,k)))
        enddo
      enddo
!
      return
      end
