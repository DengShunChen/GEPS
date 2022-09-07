      subroutine outflds_hp( itau,nx,my,my_max,idtg,ifilout          &
                           , raincu3,rainlp3,ggdef)
      use mpe
      use rank
      use index
      use const, only: kflag

      implicit  none

      integer   itau,nx,my,my_max,istat,lenc

      real      raincu3(nxp,my_max),rainlp3(nxp,my_max)

      character ifilout*60, ggdef*4, ihdg*26
      integer*8 idtg
!
! local work arrays
!
      real      glob(nx,my),wrk(nxp,my_max) 

      integer   jj,j,nxj,i
!
      lenc = nx*my
!=======================================================================
      if(myrank .eq. 0) print*,'   in outflds_hp for tau= ',itau
!=======================================================================
!byl      call mpe2d_unify(glob,raincu3)
!byl      call mpe2d_unify(glob1,rainlp3)
!      call mpe_unify_1(glob,raincu3,nx,my,2,mpe_double)
!      call mpe_unify_1(glob1,rainlp3,nx,my,2,mpe_double)
      call syslbl ('b00632',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      call unify_reduceintp(nx,my,my_max,raincu3,glob)
      call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
!
      call syslbl ('b00642',idtg,itau,ggdef,ihdg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob1,nxdef,nx,my)
      call unify_reduceintp(nx,my,my_max,rainlp3,glob)
      call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
!
      call syslbl ('b00622',idtg,itau,ggdef,ihdg)
      do 98 jj = 1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 98 i=1,nxj
       wrk(i,jj)=raincu3(i,jj)+rainlp3(i,jj)
 98   continue
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)

!=======================================================================
      return
      end
