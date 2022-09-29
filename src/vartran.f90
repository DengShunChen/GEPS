      subroutine vartran(vorten,divten,phiten,nw,jtrun,jtmax,lev,k, &
                         rad,omega,h,iflag)
!
!    purpose : conversion between dimension and dimensionless forms
!--------------------------------------------------------------------
! *** input ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!  nw      :  total wavenumber index array
!  lev     :  total vertical levels
!  k       :  vetical level of (non)dimensionlize
!  rad     :  radius of earth
!  omega   :  angular velocity of earth
!  h       :  coefficient array of (non)dimensionlize
!  iflag   : +2 = to nondimension form
!          : -2 = to dimensional form
! **** output ****
!  vorten  :  vorticity tendency(correction) array in spectrum domaon
!  divten  :  divergence tendency(correction) array in spectrum domain
!  phiten  :  geopotantial tendency(correction) array in spectrum domain
!************************************
      use index
      use const, only: RTYPE

      implicit none

      integer  jtrun,jtmax,lev,k,iflag,m,l,mf

      real(kind=RTYPE) vorten(lev,2,jtrun,jtmax),   &
                       divten(lev,2,jtrun,jtmax),   &
                       phiten(lev,2,jtrun,jtmax)
!     real h(jtrun,jtmax,lev)
!byl      real h(jtrun,jtmax,levF) ! 2dMPI
      real h(jtrun,jtmax) ! 2dMPI
      integer nw(jtrun,jtmax)

      real omega2,omega,omga2r2,rad,temp1,temp2

      omega2=omega*omega
      omga2r2=omega2*rad*rad
!
!  do nondimensionlize which are tendency
!
      if(iflag.eq.2)then

         do 110 m = 1, mlistnum
          mf=mlist(m)
         do 110 l = mf, jtrun
!byl          vorten(k,1,l,m)=-h(l,m,k)*vorten(k,1,l,m)/omega2
!byl          vorten(k,2,l,m)=-h(l,m,k)*vorten(k,2,l,m)/omega2
          vorten(k,1,l,m)=-h(l,m)*vorten(k,1,l,m)/omega2
          vorten(k,2,l,m)=-h(l,m)*vorten(k,2,l,m)/omega2
          temp1=divten(k,2,l,m)/omega2
          temp2=divten(k,1,l,m)/omega2
!byl          divten(k,1,l,m)=-h(l,m,k)*temp1
!byl          divten(k,2,l,m)= h(l,m,k)*temp2
          divten(k,1,l,m)=-h(l,m)*temp1
          divten(k,2,l,m)= h(l,m)*temp2
          phiten(k,1,l,m)=phiten(k,1,l,m)/omga2r2/omega
          phiten(k,2,l,m)=phiten(k,2,l,m)/omga2r2/omega
 110    continue
!
!  do dimensionlize which are correct term now
!
      else if(iflag.eq.-2)then

         do 130 m = 1, mlistnum
          mf=mlist(m)
         do 130 l = mf, jtrun
	  if (l.ne.1) then
!byl            vorten(k,1,l,m)=-vorten(k,1,l,m)*omega/h(l,m,k)
!byl            vorten(k,2,l,m)=-vorten(k,2,l,m)*omega/h(l,m,k)
!byl            temp1= divten(k,2,l,m)*omega/h(l,m,k)
!byl            temp2=-divten(k,1,l,m)*omega/h(l,m,k)
            vorten(k,1,l,m)=-vorten(k,1,l,m)*omega/h(l,m)
            vorten(k,2,l,m)=-vorten(k,2,l,m)*omega/h(l,m)
            temp1= divten(k,2,l,m)*omega/h(l,m)
            temp2=-divten(k,1,l,m)*omega/h(l,m)
            divten(k,1,l,m)= temp1
            divten(k,2,l,m)= temp2
	  endif
 130     continue
!
         do 135 m = 1, mlistnum
          mf=mlist(m)
         do 135 l = mf, jtrun
          phiten(k,1,l,m)= phiten(k,1,l,m)*omga2r2
          phiten(k,2,l,m)= phiten(k,2,l,m)*omga2r2
 135     continue
      endif
!
      return
      end
