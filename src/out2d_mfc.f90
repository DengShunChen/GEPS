      subroutine out2d_mfc (nx,lev,my,my_max,ifilout,itau,idtg,ntau    &
                            ,rain1,raintot,glob,t2,q2,rh2,rh10,u10,v10 &
                            ,tmax,tmin,td,rld,sld,ctot,slpty,ggdef )
!
      use rank
      use mpe
      use index
!
      implicit  none

      integer   nx,lev,my,my_max,itau,ntau,num
      parameter (num=14)

      real      raintot(nxp,my_max),t2(nxp,my_max),u10(nxp,my_max),   &
                v10(nxp,my_max),ctot(nxp,my_max),slpty(nxp,my_max)

      real rain1(nxp,my_max),q2(nxp,my_max),rh2(nxp,my_max),          &
           rh10(nxp,my_max),tmax(nxp,my_max),tmin(nxp,my_max),        &
           td(nxp,my_max),rld(nxp,my_max),sld(nxp,my_max)
!
      real mfcout(nxp,my_max,num)
!
      character*4 ggdef
      integer*8 idtg
      character*6 dmskey(num)
!
      real      glob(nx,my),mout(nx,my)
!
      character*80 ifilout
      character*26 ihdg,ihdg2
!
      integer   n,levz,lenc,lenc2,i,ia,kk,j,nxj,istat,jj
      real      tnshun
!
      data dmskey/'b00621','b0062t','b02100','b02500','b02510', &
                  'b10200','b10210','b02171','b02181','b02150', &
                  's003x0','s003u0','x00770','ssl010'/


!
      tnshun= 1.0
      lenc= nx*my
      lenc2= lev*my
!      
      mfcout(:,:,1)=rain1(:,:)
      mfcout(:,:,2)=raintot(:,:)
      mfcout(:,:,3)=t2(:,:)
      mfcout(:,:,4)=q2(:,:)
      mfcout(:,:,5)=rh2(:,:)
      mfcout(:,:,6)=u10(:,:)
      mfcout(:,:,7)=v10(:,:)
      mfcout(:,:,8)=tmax(:,:)
      mfcout(:,:,9)=tmin(:,:)
      mfcout(:,:,10)=td(:,:)
      mfcout(:,:,11)=rld(:,:)
      mfcout(:,:,12)=sld(:,:)
      mfcout(:,:,13)=ctot(:,:)
      mfcout(:,:,14)=slpty(:,:)
!
! Total Precp.
!byl      call mpe2d_unify(glob,raintot)
      do n=1,num
        call syslbl (dmskey(n),idtg,ntau,ggdef,ihdg)
        call unify_reduceintp(nx,my,my_max,mfcout(1,1,n),mout)
        if ( myrank .eq. n-1 ) then
          glob=mout
          ihdg2=ihdg
        endif
      enddo
!
      if (myrank .lt. num ) call dmswrit_split(nx,my,ihdg2,lenc,'H',ifilout,glob,istat)
!
!! rh10
!!byl      call mpe2d_unify(glob,rh10)
!      call syslbl ('b10510',idtg,ntau,ggdef,ihdg)
!      call unify_reduceintp(nx,my,my_max,rh10,glob)
!!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!!     call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
!      call dmswrit_mfc(nx,my,ihdg,lenc,'H',ifilout,glob,istat)


      return
      end

