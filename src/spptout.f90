      subroutine spptout(nx,my,my_max,lev,sppt2d500,sppt2d1000  &
                        ,sppt2d2000,sppt3d,recn,tau)

      use mpe
      use rank
      use index

      implicit none

      integer      i,j,k,jj,nxj,ihead
      integer      nx,my,my_max,lev,recn,nxmy4
      real,intent(in):: sppt2d500(nxp,my_max),sppt2d1000(nxp,my_max) &
                       ,sppt2d2000(nxp,my_max)
      real         sppt3d(nxp,lev,my_max),tau
      real         glob(nx,my),temp(nxp,my_max)
      real*4       glob4(nx,my)

      ihead=15

      nxmy4=nx*my*4
      if ( myrank .eq. 0 )  &
        open(ihead,file='strucsppt3dv5.dat',access='direct'         &
            ,form='unformatted',recl=nxmy4,status='unknown')
!    
      call mpe2d_unify(glob,sppt2d500)
      do jj=1,jlistnum
        j=jlist1(jj)
        if( lreduce.eq.1 ) call reduceintp (glob(1,j),nxdef(j),nx,1)
      enddo
      call mpe_unify(glob,nx,my,5,mpe_double)
      if ( myrank .eq. 0 ) then
        glob4=glob
        write(ihead,rec=recn) glob4
        recn=recn+1
      endif
!
      call mpe2d_unify(glob,sppt2d1000)
      do jj=1,jlistnum
        j=jlist1(jj)
        if( lreduce.eq.1 ) call reduceintp (glob(1,j),nxdef(j),nx,1)
      enddo
      call mpe_unify(glob,nx,my,5,mpe_double)
      if ( myrank .eq. 0 ) then
        glob4=glob
        write(ihead,rec=recn) glob4
        recn=recn+1
      endif
!
      call mpe2d_unify(glob,sppt2d2000)
      do jj=1,jlistnum
        j=jlist1(jj)
        if( lreduce.eq.1 ) call reduceintp (glob(1,j),nxdef(j),nx,1)
      enddo
      call mpe_unify(glob,nx,my,5,mpe_double)
      if ( myrank .eq. 0 ) then
        glob4=glob
        write(ihead,rec=recn) glob4
        recn=recn+1
      endif
!
      do k=lev,1,-1
        do jj=1,jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
          temp(i,jj)=sppt3d(i,k,jj)
          enddo
        enddo
        call mpe2d_unify(glob,temp)
        do jj=1,jlistnum
          j=jlist1(jj)
         if( lreduce.eq.1 ) call reduceintp (glob(1,j),nxdef(j),nx,1)
        enddo
        call mpe_unify(glob,nx,my,5,mpe_double)
        if ( myrank.eq.0 ) then
          glob4=glob
          write(ihead,rec=recn) glob4
          recn=recn+1
        endif
      enddo

      if ( myrank.eq.0 ) close(ihead)
!creating ctl file
      call spptctl(nx,my,lev,tau)
!     

 
      return
      end
!----
      subroutine spptctl(nx,my,lev,tau)
!
      use const, only : idtg,sinl,sigma
      use rank, only : myrank
!
!
      integer nx,my,nxj,j,k,lev,itau,ch,iter,remd,js,je,kk,lev1
      real    pi,r2d,dlon,mlat(my),prsl(lev),tau
      character*30 forydef,forzdef
      
      character yy*4,dd*2,hh*2,mm*2,dtg*12
      character*3 mon(12)
      logical jrem

      data mon/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep' &
               ,'Oct','Nov','Dec'/

!
      itau=tau
      ch=10
      lev1=1
!
      pi=4.0*atan(1.0)
      r2d=180./pi
      dlon=360./float(nx)
!
      do j=1,my
        mlat(j)=asin(sinl(j))*r2d
      enddo

      do k=1,lev
        kk=lev-k+1
        prsl(kk)=sigma(k,2)+sigma(k+1,2)
        prsl(kk)=prsl(kk)+(sigma(k,1)+sigma(k+1,1))*1000.
        prsl(kk)=0.5*prsl(kk)/1000.
      enddo
!
        write(forydef,8) my
        write(forzdef,9) lev

        write(dtg,'(I12)') idtg
        read(dtg,'(A4,I2,A2,A2,A2)') yy,mn,dd,hh,mm

        if ( myrank .eq. 0 ) then
        OPEN(UNIT=ch, FILE='sppt.ctl', STATUS='UNKNOWN'             &
           , ACCESS='SEQUENTIAL')

        write(ch,'(A23)') 'dset ^strucsppt3dv5.dat'
        write(ch,'(A18)') 'options big_endian'
        write(ch,'(A12)') 'undef -999.0'
        write(ch,12) 'ydef' ,my, 'levels'
        iter=my/8
        remd=mod(my,8)
        jrem=(remd .eq. 0)
        write(forydef,8) remd
        do j=1,iter
         js=1+8*(j-1)
         je=js+7
         write(10,10) mlat(js:je)
        enddo
        if ( .not. jrem ) write(10,forydef) mlat(je+1:my)
         
        write(ch,13) 'xdef'  ,nx, 'linear 0.0',dlon
        write(ch,14) 'tdef',itau, 'linear',hh,'Z',dd,mon(mn),yy,'1hr'
        write(ch,12) 'zdef' ,lev, 'levels '
        iter=lev/8
        remd=mod(lev,8)
        jrem=(remd .eq. 0)
        write(forzdef,9) remd
        do j=1,iter
         js=1+8*(j-1)
         je=js+7
         write(10,11) prsl(js:je)
        enddo
        if ( .not. jrem ) write(10,forzdef) prsl(je+1:my)
        write(ch,'(A6)') 'vars 4'
        write(ch,15) 'sppt500' ,lev1,'99','500km perturbation'
        write(ch,15) 'sppt1000',lev1,'99','1000km perturbation'
        write(ch,15) 'sppt2000',lev1,'99','1000km perturbation'
        write(ch,15) 'tsppt'   , lev,'99','temp perturbation'
        write(ch,'(A7)') 'endvars'

        close(ch)
        endif

8       format("(",I4,"(2x,F11.7))")
9       format("(",I4,"(2x,F7.5))")
10      format(8(2x,F11.7))
11      format(8(2x,F7.5))
12      format(A4,1X,I4,1X,A6)
13      format(A4,1X,I4,1X,A10,1X,F10.7)
14      format(A4,1X,I4,1X,A6,1X,A2,A1,A2,A3,A4,1X,A3)
15      format(A9,1X,I3,1X,A2,1X,A20)

        return
        end
