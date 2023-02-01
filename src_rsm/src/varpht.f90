      subroutine varpht (nx,my,my_max,lev,cp,dsig,pk2,pk,pt,sgeo,tt,phi &
      ,                  phicor,tcor)
!
!  varpht does a calculus of variations adjustment of input
!  geopotential and temperature in arrays phi and tt, yielding
!  hydrostatically consistent output geopotential and temperature
!  in arrays phi and tt.  see barker (1980) monthly weather
!  review no. 8 "solving for temperature using unnaturally
!  latticed hydrostatic equations" for more technical details.
!
!  **** input ****
!
!  cp: specific heat of air
!  dsig: sigma coordinate layer thicknesses
!  pk2: half level exner function
!  pk: full level exner function
!  pt: terrain pressure
!  sgeo: terrain geopotential
!  tt: temperatures
!  phi: geopotentials
!  phicor:  work array
!  tcor:  work array
!
! **** output ****
!
!  tt: adjusted temperatures
!  phi: adjusted geopotentials
!
!     include '../include/mpe.h'
!     include '../include/rank.h'
!     include '../include/index.h'
      use mpe
      use rank
      use index

      implicit  none

      integer   nx,my,my_max,lev
      integer   i,j,k,jj,kk
      real      cp,temx

      real      dsig(lev),tt(nx,lev,my_max)                          &
      , phi(nx,lev,my_max),pt(nx,my),sgeo(nx,my),pk(nx,lev,my_max)   &
      , pk2(nx,lev,my_max),phicor(nx,lev,my_max),tcor(nx,lev,my_max)
!
!sun  include '../include/paramt.h'.. change im,lm to nx,lev

      real      dp(nx,lev),elam(nx,lev),a(nx,lev),b(nx,lev),c(nx,lev)
      real      tweght(lev),rmst(lev),rmsf(lev),biast(lev),biasf(lev)
      character*16 crec
!
!  tweght is an array of coefficients used tt determine the
!  respective weights given tt the input geopotential and
!  temperature.  the larger the value of tweght, the more the
!  adjustment is biased to fitting the geopotential at the expense
!  of temperature fit.  values in excess of those used here will
!  cause unaceptable buckling of the temperature sounding.
!
!      data tweght/lev*0.0001/
      do 1 k=1,lev
        tweght(k)=0.0001
 1    continue
      tweght(lev)=   0.00001
      tweght(lev-1)= 0.0001
      tweght(1)=.0005
      tweght(2)=.00025
!
!  compute pressure variables in the sigma levels of the
!  forecast model.
!
      do 100 jj =1, jlistnum
      j=jlist1(jj)
      do 121 i=1,nx
      dp(i,1)= dsig(1)*pt(i,j)
      tcor(i,1,jj)= cp*pk2(i,1,jj)/pk(i,1,jj)-cp
  121 continue
!ocl nounroll
      do 111 k=1,lev-1
      do 113 i=1,nx
      dp(i,k+1)= dsig(k+1)*pt(i,j)
      tcor(i,k+1,jj)= cp*(pk2(i,k+1,jj)/pk(i,k+1,jj))-cp
      phicor(i,k,jj)= cp-cp*(pk2(i,k,jj)/pk(i,k+1,jj))
      elam(i,k)= phi(i,k,jj)-phi(i,k+1,jj)-(tcor(i,k,jj)*tt(i,k,jj) &
       +phicor(i,k,jj)*tt(i,k+1,jj))
      a(i,k+1)= (tcor(i,k+1,jj)*phicor(i,k,jj)*tweght(k+1)-1.0)     &
       *dp(i,k)/dp(i,k+1)
      b(i,k)=tcor(i,k,jj)*tcor(i,k,jj)*tweght(k)+1.0                &
       +dp(i,k)/dp(i,k+1)*(phicor(i,k,jj)*phicor(i,k,jj)            &
       *tweght(k+1)+1.0)
      c(i,k)=phicor(i,k,jj)*tcor(i,k+1,jj)*tweght(k+1)-1.0
  113 continue
  111 continue
!
      do 101 i=1,nx
      elam(i,lev)= phi(i,lev,jj)-sgeo(i,j)-tcor(i,lev,jj)*tt(i,lev,jj)
      b(i,lev)=tcor(i,lev,jj)*tcor(i,lev,jj)*tweght(lev)+1.0
  101 continue
!
      call trdiv2(nx,lev,a,b,c,elam)
!
!ocl nounroll
      do 400 k=2,lev
      do 400 i=1,nx
      tcor(i,k,jj)=(elam(i,k)*tcor(i,k,jj)*dp(i,k)                &
       +elam(i,k-1)*phicor(i,k-1,jj)*dp(i,k-1))*tweght(k)/dp(i,k)
      phicor(i,k,jj)=-(elam(i,k)*dp(i,k)-elam(i,k-1)              &
       *dp(i,k-1))/dp(i,k)
      phi(i,k,jj)= phi(i,k,jj)+phicor(i,k,jj)
      tt(i,k,jj)= tt(i,k,jj)+tcor(i,k,jj)
  400 continue
!
      do 510 i=1,nx
      tcor(i,1,jj)= elam(i,1)*tcor(i,1,jj)*tweght(1)
      phicor(i,1,jj)= -elam(i,1)
      phi(i,1,jj)= phi(i,1,jj)+phicor(i,1,jj)
      tt(i,1,jj)= tt(i,1,jj)+tcor(i,1,jj)
  510 continue
!
  100 continue
!
      do 600 kk=1,lev
      write(crec,800) kk
  800 format('var phi cng k=',i2)
!fj+
!     call qmaxn3 (phicor,crec(1:8),crec(9:16),1,1,kk,nx,my,lev)
      call qmaxn3p(phicor,crec(1:8),crec(9:16),1,1,kk,nx,my_max,lev)
!fj-
  600 continue
!
      do 610 kk=1,lev
      write(crec,900) kk
  900 format('var t chng k=',i2)
!fj+
!     call qmaxn3 (tcor,crec(1:8),crec(9:16),1,1,kk,nx,my,lev)
      call qmaxn3p(tcor,crec(1:8),crec(9:16),1,1,kk,nx,my_max,lev)
!fj-
  610 continue
!
!fj >>
!fj   call srms(nx,lev,my,rmst,tcor)
!fj   call srms(nx,lev,my,rmsf,phicor)
      call srms(nx,lev,my,my_max,rmst,tcor)
      call srms(nx,lev,my,my_max,rmsf,phicor)
!fj <<
!
      temx= 1.0/(nx*my)
      do 615 k=1,lev
      biast(k)= 0.0
      biasf(k)= 0.0
  615 continue
      do 625 jj =1, jlistnum
      j=jlist1(jj)
      do 625 k=1,lev
      do 625 i=1,nx
      biast(k)= biast(k)+tcor(i,k,jj)
      biasf(k)= biasf(k)+phicor(i,k,jj)
  625 continue
!
      call mpe_global_sum(biast,lev,mpe_double)
      call mpe_global_sum(biasf,lev,mpe_double)
!
      do 626 k=1,lev
      biast(k)= biast(k)*temx
      biasf(k)= biasf(k)*temx
  626 continue
!
!fj Change to print only RANK 0 PE >>
      if(myrank .eq. 0) then
      print 620
      print 630,(k,rmst(k),rmsf(k),biast(k),biasf(k),k=1,lev)
      endif
!fj <<
  620 format(' lev:    rms tcor:   rms phicor:  tcor bias: ',2x   &
      ,'phicor bias',/)
  630 format(i4,2x,f10.3,2x,f10.3,3x,f10.3,3x,f10.3)
!
      return
      end
