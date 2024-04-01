      subroutine outflds_hp( itau,nx,my,my_max,idtg          &
                           ,raincu1,rainlp1, raincu3,rainlp3,ggdef)
      use mpe
      use rank
      use index
      use const ,only: outdms ,outgrb2 ,ifilout_grb ,RTYPE,kflag
      use mod_grb2_param , only :ofdir,wrt_grb2_v2,wrt_grb2_accu_v2

      implicit  none

      integer   itau,nx,my,my_max,istat,lenc

      real,intent(in):: raincu1(nxp,my_max),rainlp1(nxp,my_max)
      real              raincu3(nxp,my_max),rainlp3(nxp,my_max)

      character ggdef*4
      integer*8 idtg
!
! local work arrays
!
      real(kind=RTYPE) glob(nx,my),wrk(nxp,my_max)

      integer   jj,j,nxj,i
!
!accum
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i = 1,nxj
          rainlp3(i,jj)=rainlp3(i,jj)+rainlp1(i,jj)
          raincu3(i,jj)=raincu3(i,jj)+raincu1(i,jj)
        enddo
      enddo
      
      if( mod( itau , 3 ) == 0 )then
        if( outgrb2 == 1)then
 134                      format( A  ,A ,I10.10 , i4.4       )
             write(ofdir,134 )trim(ifilout_grb),'/',idtg/100 ,itau
             if(myrank==0) call system("mkdir -p "//trim(ofdir) )
        endif

      lenc = nx*my
!=======================================================================
      if(myrank .eq. 0) print*,'   in outflds_hp for tau= ',itau
!=======================================================================
!      call mpe_unify_1(glob,raincu3,nx,my,2,mpe_double)
!      call mpe_unify_1(glob1,rainlp3,nx,my,2,mpe_double)
      call syslbl_w ('b00632',idtg,itau,ggdef)
      wrk=raincu3
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      if(outdms.gt.0)call dmswrit(nx,my,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) &
                 call wrt_grb2_accu_v2(itau,0,1,10,2,103,0,0,1,3,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
!
      call syslbl_w ('b00642',idtg,itau,ggdef)
      wrk=rainlp3
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      if(outdms.gt.0)call dmswrit(nx,my,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) &
                  call wrt_grb2_accu_v2(itau,0,1,9,2,103,0,0,1,3,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
!
      call syslbl_w ('b00622',idtg,itau,ggdef)
      do 98 jj = 1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 98 i=1,nxj
       wrk(i,jj)=raincu3(i,jj)+rainlp3(i,jj)
 98   continue
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      if(outdms.gt.0)call dmswrit(nx,my,lenc,kflag,glob,istat)
      if(outgrb2==1.and.myrank==0) &
                  call wrt_grb2_accu_v2(itau,0,1,7,2,103,0,0,1,3,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)

      endif !( mod( itau , 3 ) == 0 )
!=======================================================================
      return
      end
