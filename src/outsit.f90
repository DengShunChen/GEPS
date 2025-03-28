      SUBROUTINE sitout(nx,my,my_max,itau,idtg,num  &
                       ,whtlev,ggdef)
      
      use rank
      use mpe
      use index
      use const,              only:RTYPE,kflag,outdms
      use mod_sitgrid
      use mod_sit_control,    only:outsitlev

      implicit none
      integer nx,my,my_max,itau,num
      real whtlev(num)
   

      integer*8 idtg
      character*4 ggdef
      integer lenc,n,k,jj,j,nxj,ii,istat
      character*6 lrec
      real(kind=RTYPE) wk1(nx,my),pout(nx,my),globp(nxp,my_max)
      integer ncnt
 

      lenc= nx*my
  
      ncnt=0
      do 10 k = 0, outsitlev+1
        write( lrec, '(i3.3,a3)' ) k,'SWT'         !!sit wt
        call syslbl_w (lrec,idtg,itau,ggdef)
        globp=sitwt(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,wk1)
        call split(nx,my,lenc,ncnt,wk1,pout)
   10 continue
      if(outdms.gt.0) then
      if(myrank .lt. ncnt) call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif

      ncnt=0
      do 20 k = 0, outsitlev+1
        write( lrec, '(i3.3,a3)' ) k,'OWT'        !!sit obswt
        call syslbl_w (lrec,idtg,itau,ggdef)
        globp=obswt(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,wk1)
        call split(nx,my,lenc,ncnt,wk1,pout)
   20 continue
      if(outdms.gt.0) then
      if(myrank .lt. ncnt) call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif

  
      end subroutine sitout


!--------------------------------------------------------------------------

      SUBROUTINE storesittau(nxj,jj,nxp,my_max,lkvl  &
                            ,wttau,wstau,wutau,wvtau,dtx)

      use rank
      use mpe
      use mod_sitgrid,       only:sitwttau,sitwstau,sitwutau,sitwvtau &
                                ,sitmask
      use mod_sit_control,only: xmissing

      implicit none

      integer nxj,jj,nxp,my_max,lkvl

      real  wttau(nxp,my_max,0:lkvl+1), wstau(nxp,my_max,0:lkvl+1)    &
          , wutau(nxp,my_max,0:lkvl+1), wvtau(nxp,my_max,0:lkvl+1)

      real dtx
      integer k,ii

       do ii = 1, nxj

        do k=0,lkvl+1
         if ((sitmask(ii,jj).eq. 1.) .and. (wttau(ii,jj,k) .ne. xmissing)) then
          sitwttau(ii,jj,k)=sitwttau(ii,jj,k)+wttau(ii,jj,k)*dtx
         else
          sitwttau(ii,jj,k)=0.
         endif

         if ((sitmask(ii,jj).eq. 1.) .and. (wstau(ii,jj,k) .ne. xmissing)) then
          sitwstau(ii,jj,k)=sitwstau(ii,jj,k)+wstau(ii,jj,k)*dtx
         else
          sitwstau(ii,jj,k)=0.
         endif

         if ((sitmask(ii,jj).eq. 1.) .and. (wutau(ii,jj,k) .ne. xmissing)) then
          sitwutau(ii,jj,k)=sitwutau(ii,jj,k)+wutau(ii,jj,k)*dtx
         else
          sitwutau(ii,jj,k)=0.
         endif

         if ((sitmask(ii,jj).eq. 1.) .and. (wvtau(ii,jj,k) .ne. xmissing)) then
          sitwvtau(ii,jj,k)=sitwvtau(ii,jj,k)+wvtau(ii,jj,k)*dtx
         else
          sitwvtau(ii,jj,k)=0.
         endif

        enddo
       enddo


      END SUBROUTINE storesittau


!--------------------------------------------------------------------------

      SUBROUTINE storesit24(nxj,jj,nxp,my_max,lkvl  &
                           ,wt24,ws24,wu24,wv24,dtx)

      use rank
      use mpe
      use mod_sitgrid,     only: sitwt24,sitws24,sitwu24,sitwv24 &
                                ,sitmask
      use mod_sit_control,only: xmissing

      implicit none

      integer nxj,jj,nxp,my_max,lkvl

      real  wt24(nxp,my_max,0:lkvl+1), ws24(nxp,my_max,0:lkvl+1)    &
          , wu24(nxp,my_max,0:lkvl+1), wv24(nxp,my_max,0:lkvl+1)  

      real dtx
      integer k,ii

       do ii = 1, nxj

        do k=0,lkvl+1
          if ((sitmask(ii,jj).eq. 1.) .and. (wt24(ii,jj,k) .ne. xmissing)) then
            sitwt24(ii,jj,k)=sitwt24(ii,jj,k)+wt24(ii,jj,k)*dtx
          else
            sitwt24(ii,jj,k)=0.
          endif

          if ((sitmask(ii,jj).eq. 1.) .and. (ws24(ii,jj,k) .ne. xmissing)) then
            sitws24(ii,jj,k)=sitws24(ii,jj,k)+ws24(ii,jj,k)*dtx
          else
            sitws24(ii,jj,k)=0.
          endif

          if ((sitmask(ii,jj).eq. 1.) .and. (wu24(ii,jj,k) .ne. xmissing)) then
            sitwu24(ii,jj,k)=sitwu24(ii,jj,k)+wu24(ii,jj,k)*dtx
          else
            sitwu24(ii,jj,k)=0.
          endif

          if ((sitmask(ii,jj).eq. 1.) .and. (wv24(ii,jj,k) .ne. xmissing)) then
            sitwv24(ii,jj,k)=sitwv24(ii,jj,k)+wv24(ii,jj,k)*dtx
          else
            sitwv24(ii,jj,k)=0.
          endif

        enddo
       enddo


      END SUBROUTINE storesit24


!--------------------------------------------------------------------------
      SUBROUTINE writesitmean(nx,my,my_max,lkvl,itau,idtg,ggdef)

      use rank
      use mpe
      use index
      use const,           only:RTYPE,kflag,outdms
      use mod_sitgrid,     only:sitwttau,sitwstau,sitwutau,sitwvtau &
                               ,dtsittau
      use mod_sit_control, only: xmissing,outsitlev

      implicit none

      integer nx,my,my_max,lkvl,itau


      integer*8 idtg
      character*4 ggdef
      integer lenc,k,jj,j,nxj,ii,istat
      character*6 lrec
      real(kind=RTYPE) wk1(nx,my),pout(nx,my),glob2d(nxp,my_max)
      integer ncnt

      lenc= nx*my


      ncnt=0
      do 10 k = 0, outsitlev+1
        write( lrec, '(i3.3,a3)' ) k,'WTT'
        call syslbl_w (lrec,idtg,itau,ggdef)
        if(dtsittau .ne. 0.) then
          glob2d(:,:)= sitwttau(:,:,k)/dtsittau
        else
          glob2d=xmissing
        endif
        call unify_reduceintp(nx,my,my_max,glob2d,wk1)
        call split(nx,my,lenc,ncnt,wk1,pout)
   10 continue

      if(outdms.gt.0)then
      if(myrank .lt. ncnt) call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif

      sitwttau=0.
      dtsittau=0.

      END SUBROUTINE writesitmean



!--------------------------------------------------------------------------
      SUBROUTINE outsit24(nx,my,my_max,lkvl,itau,idtg,ggdef)
      
      use rank
      use mpe
      use index
      use const,          only: RTYPE,kflag,outgrb2,outdms
      use mod_sitgrid,    only: sitwt24,sitws24,sitwu24,sitwv24 &
                               ,dtsit24,wtfn0,wsfn0,obswt,sitwt
      use mod_sit_control,only: xmissing,outsitlev
      use mod_grb2_param,only: wrt_grb2_v2
      implicit none

      integer nx,my,my_max,lkvl,itau
   
      integer*8 idtg
      character*4 ggdef
      integer lenc,k,jj,j,nxj,ii,istat
      character*6 lrec
      real(kind=RTYPE) tm1(nxp,my_max),tm2(nxp,my_max),tm3(nxp,my_max)
      real(kind=RTYPE) wk1(nx,my),pout(nx,my)
      integer ncnt
      integer:: ptp0(9),ptp1(9) !save grib info

      tm1=xmissing
      tm2=xmissing
      tm3=xmissing

      lenc= nx*my
  
      ncnt=0
      do 10 k = 0, outsitlev+1
        do jj =1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do ii = 1, nxj
            if(dtsit24 .ne. 0.) then
              tm1(ii,jj)  = sitwt24(ii,jj,k)/dtsit24
!              tm2(ii,j)  = wtfn0(ii,jj,k)/dtsit24
              if(obswt(ii,jj,k) .ge. 271. ) then
                tm2(ii,jj) = (obswt(ii,jj,k)-sitwt(ii,jj,k))/dtsit24
              else
                tm2(ii,jj) = 0.
              endif
              tm3(ii,jj)  = wsfn0(ii,jj,k)/dtsit24
            endif
          enddo
        enddo

        write( lrec, '(i3.3,a3)' ) k,'WTF'
        call syslbl_w (lrec,idtg,itau,ggdef)
        call unify_reduceintp(nx,my,my_max,tm1,wk1)
        ptp0=(/10,3,0,2,168,0,(k+1),-999,-999/)
        call split2(nx,my,lenc,ncnt,wk1,pout,ptp0,ptp1)
!        call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
   10 continue
      if(myrank .lt. ncnt)then
        if(outdms.gt.0)then
          call dmswrit_split(nx,my,lenc,kflag,pout,istat) 
        endif
        if(outgrb2==1)then
            call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
              ,ptp1(6),ptp1(7),pout)
        endif
      endif
      
     
      sitwt24=0.
      wtfn0=0.
      wsfn0=0.
      dtsit24=0. 
  
      end subroutine outsit24


!--------------------------------------------------------------------------
      SUBROUTINE outsitmon(nx,my,my_max,lkvl,itau,idtg,ggdef)
      
      use rank
      use mpe
      use index
      use const,              only: RTYPE,kflag,outdms
      use mod_sitgrid,        only: dtsitmon,wtfn,wtfns,wsfn,wsfns
      use mod_sit_control,    only: xmissing,outsitlev
  
      implicit none

      integer nx,my,my_max,lkvl,itau
  


      real dtmon
      integer*8 idtg
      character*4 ggdef
      integer lenc,i,j,k,ii,jj,nxj,istat
      character*6 lrec
      real(kind=RTYPE) wk1(nx,my),pout(nx,my),glob2d(nxp,my_max)
      integer ncnt

      lenc= nx*my
  
      ncnt=0
      do 10 k = 0, outsitlev+1
        write( lrec, '(i3.3,a3)' ) k,'TFM'
        call syslbl_w (lrec,idtg,itau,ggdef)
        glob2d(:,:)=wtfn(:,:,k)/dtsitmon
        call unify_reduceintp(nx,my,my_max,glob2d,wk1) 
        call split(nx,my,lenc,ncnt,wk1,pout)
   10 continue
      if(outdms.gt.0) then
      if(myrank .lt. ncnt) call dmswrit_split(nx,my,lenc,kflag,pout,istat) 
      endif

        write( lrec, '(i3.3,a3)' ) k,'TFS'
        call syslbl_w (lrec,idtg,itau,ggdef)
        glob2d(:,:)=wtfns(:,:)/dtsitmon
        call unify_reduceintp(nx,my,my_max,glob2d,wk1)
        if(outdms.gt.0) &
        call dmswrit(nx,my,lenc,kflag,wk1,istat)



!reset wtfn,wsfn,wtfns,wsfns
      dtsitmon=0.
      wtfn=0.
      wsfn=0.
      wtfns=0.
      wsfns=0.

      end subroutine outsitmon

!-------------------------------------------------------------	  
      subroutine rerun_sitgrid1(itau)

      use index
      use mpe
      use const
      use mod_sitgrid
      use const, only: RTYPE


      real(kind=RTYPE), dimension(:,:), allocatable ::     &
                   tmp1,tmp2,tmp3,tmp4,tmp5    &
                  ,tmp6,tmp7,tmp8,tmp9,tmp10   &
                  ,tmp11,tmp12,tmp13,tmp14
      real(kind=RTYPE), dimension(:,:,:), allocatable::    &
                   tm12, tm13, tm14
      real(kind=RTYPE) globp(nxp,my_max)


      integer itau,lphy
      character nfs*10, rfile*80
      integer i,j,k,ii,jj,nxj

      allocate(tmp1(nx,my),tmp2(nx,my),tmp3(nx,my),tmp4(nx,my)   &
               ,tmp5(nx,my),tmp6(nx,my),tmp7(nx,my),tmp8(nx,my)   &
               ,tmp9(nx,my),tmp10(nx,my),tmp11(nx,my))
      allocate(tmp12(nx,my),tmp13(nx,my),tmp14(nx,my))
      allocate(tm12(nx,0:1,my),tm13(nx,0:1,my),tm14(nx,0:3,my))




      do k = 0, 3  
        if(k .eq. 0) then
          globp=sitcc
          call unify_reduceintp(nx,my,my_max,globp,tmp1)
          globp=sithc
          call unify_reduceintp(nx,my,my_max,globp,tmp2)
          globp=engwac
          call unify_reduceintp(nx,my,my_max,globp,tmp3)
          globp=sc
          call unify_reduceintp(nx,my,my_max,globp,tmp4)
          globp=saltwac
          call unify_reduceintp(nx,my,my_max,globp,tmp5)
          globp=wtfns
          call unify_reduceintp(nx,my,my_max,globp,tmp6)
          globp=wsfns
          call unify_reduceintp(nx,my,my_max,globp,tmp7)
          globp=grndcapc
          call unify_reduceintp(nx,my,my_max,globp,tmp8)
          globp=grndhflx
          call unify_reduceintp(nx,my,my_max,globp,tmp9)
          globp=grndflux
          call unify_reduceintp(nx,my,my_max,globp,tmp10)
          globp=obswtb
          call unify_reduceintp(nx,my,my_max,globp,tmp11)
          globp=zsi(:,:,k)
          call unify_reduceintp(nx,my,my_max,globp,tmp12)
          globp=silw(:,:,k)
          call unify_reduceintp(nx,my,my_max,globp,tmp13)
          globp=tsnic(:,:,k)
          call unify_reduceintp(nx,my,my_max,globp,tmp14)
        endif
        if(k .eq. 1) then
          globp=zsi(:,:,k)
          call unify_reduceintp(nx,my,my_max,globp,tmp12)
          globp=silw(:,:,k)
          call unify_reduceintp(nx,my,my_max,globp,tmp13)
          globp=tsnic(:,:,k)
          call unify_reduceintp(nx,my,my_max,globp,tmp14)
        endif
        if(k .ge. 2) then
          globp=tsnic(:,:,k)
          call unify_reduceintp(nx,my,my_max,globp,tmp14)
        endif

        do j=1,my
          do i=1,nx
            if(k .le. 1) then
              tm12(i,k,j)=tmp12(i,j)
              tm13(i,k,j)=tmp13(i,j)
              tm14(i,k,j)=tmp14(i,j)
            elseif(k .ge. 2) then
              tm14(i,k,j)=tmp14(i,j)
            endif
          enddo
        enddo

      enddo

      call chlen (phyout,80,lphy)  

      if(myrank .eq. 0) then
        print *,'ready to write rerun_sitgrid1'
        write (nfs,801) itau,'_',1
  801   format('sit',i4.4,a1,i2.2)
        rfile = phyout(1:lphy)//nfs
        open (unit=10,file=rfile,form='unformatted')
        write(10) tmp1,tmp2,tmp3,tmp4,tmp5  &
                 ,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11 &
                 ,tm12,tm13,tm14
        call flush(10)
        close (10)
      endif

      deallocate(tmp1,tmp2,tmp3,tmp4,tmp5)
      deallocate(tmp6,tmp7,tmp8,tmp9,tmp10)
      deallocate(tmp11,tmp12,tmp13,tmp14)
      deallocate(tm12,tm13,tm14)
  
      return
      end subroutine rerun_sitgrid1

!-------------------------------------------------------------	  
      subroutine rerun_sitgrid2(itau)

      use index
      use mpe
      use const
      use mod_sitgrid

  
      real(kind=RTYPE), dimension(:,:), allocatable::     &
                  tmp11,tmp12,tmp13,tmp14,tmp15,tmp16

      real(kind=RTYPE), dimension(:,:,:), allocatable::   &
                tm11,tm12,tm13,tm14,tm15,tm16
      real(kind=RTYPE) globp(nxp,my_max)
   
      integer itau,lphy
      character nfs*10, rfile*80
      integer i,j,k,ii,jj,nxj

      allocate( tmp11(nx,my),tmp12(nx,my),tmp13(nx,my) &
               ,tmp14(nx,my),tmp15(nx,my),tmp16(nx,my) )
      allocate( tm11(nx,0:lkvl+1,my),tm12(nx,0:lkvl+1,my) &
               ,tm13(nx,0:lkvl+1,my),tm14(nx,0:lkvl+1,my) &
               ,tm15(nx,0:lkvl+1,my),tm16(nx,0:lkvl+1,my) )

  
  
      do k = 0, lkvl+1
        globp=sitwt(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,tmp11)
        globp=sitwu(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,tmp12)
        globp=sitwv(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,tmp13)
        globp=sitww(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,tmp14)
        globp=sitws(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,tmp15)
        globp=sitwtke(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,tmp16)

        do j=1,my
          do i=1,nx
            tm11(i,k,j)=tmp11(i,j)
            tm12(i,k,j)=tmp12(i,j)
            tm13(i,k,j)=tmp13(i,j)
            tm14(i,k,j)=tmp14(i,j)
            tm15(i,k,j)=tmp15(i,j)
            tm16(i,k,j)=tmp16(i,j)
          enddo
        enddo
      enddo



      call chlen (phyout,80,lphy)  

      if(myrank .eq. 0) then
        print *,'ready to write rerun_sitgrid2'
        write (nfs,801) itau,'_',2
  801   format('sit',i4.4,a1,i2.2)
        rfile = phyout(1:lphy)//nfs
        open (unit=10,file=rfile,form='unformatted')
        write(10) tm11,tm12,tm13,tm14,tm15,tm16
        call flush(10)
        close (10)
      endif

      deallocate(tmp11,tmp12,tmp13,tmp14,tmp15,tmp16) 
      deallocate( tm11, tm12, tm13, tm14, tm15, tm16) 
 
      return
      end subroutine rerun_sitgrid2

!-------------------------------------------------------------	  
      subroutine rerun_sitgrid3(itau)

      use index
      use mpe
      use const
      use mod_sitgrid

  
      real(kind=RTYPE), dimension(:,:), allocatable:: tmp11,tmp12
      real(kind=RTYPE), dimension(:,:,:), allocatable:: tm11,tm12
      real(kind=RTYPE) globp(nxp,my_max)

      integer itau,lphy
      character nfs*10, rfile*80
      integer i,j,k,ii,jj,nxj

      allocate( tmp11(nx,my),tmp12(nx,my) )
      allocate( tm11(nx,0:lkvl+1,my),tm12(nx,0:lkvl+1,my) )
  
  
      do k = 0, lkvl+1
        globp=wtfn(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,tmp11)
        globp=wsfn(:,:,k)
        call unify_reduceintp(nx,my,my_max,globp,tmp12)

        do j=1,my
          do i=1,nx
            tm11(i,k,j)=tmp11(i,j)
            tm12(i,k,j)=tmp12(i,j)
          enddo
        enddo
      enddo



      call chlen (phyout,80,lphy)  

      if(myrank .eq. 0) then
        print *,'ready to write rerun_sitgrid3'
        write (nfs,801) itau,'_',3
  801   format('sit',i4.4,a1,i2.2)
        rfile = phyout(1:lphy)//nfs
        open (unit=10,file=rfile,form='unformatted')
        write(10) tm11,tm12
        call flush(10)
        close (10)
      endif

      deallocate(tmp11,tmp12,tm11,tm12)

  
      return
      end subroutine rerun_sitgrid3


!-------------------------------------------------------------	  
      subroutine readrerun_sitgrid1(itau)

      use index
      use mpe
      use const
      use mod_sitgrid

      real(kind=RTYPE), dimension(:,:), allocatable:: &
                          tmp1,tmp2,tmp3,tmp4,tmp5    &
                         ,tmp6,tmp7,tmp8,tmp9,tmp10   &
                         ,tmp11,tmp12,tmp13,tmp14
      real(kind=RTYPE), dimension(:,:,:), allocatable:: tm12,tm13,tm14


      integer itau,lphy
      character nfs*10, rfile*80  
      logical flag
      integer i,j,k,ii,jj,nxj,mpe_typ
       

      allocate( tmp1(nx,my),tmp2(nx,my),tmp3(nx,my),tmp4(nx,my)   &
               ,tmp5(nx,my),tmp6(nx,my),tmp7(nx,my),tmp8(nx,my)   &
               ,tmp9(nx,my),tmp10(nx,my),tmp11(nx,my) )
      allocate( tmp12(nx,my),tmp13(nx,my),tmp14(nx,my) )
      allocate( tm12(nx,0:1,my),tm13(nx,0:1,my),tm14(nx,0:3,my) )

      call chlen (phyout,80,lphy)
  
      if(myrank .eq. 0) then
        print *,'ready to read rerun_sitgrid1'
        write (nfs,801) itau,'_',1
  801   format('sit',i4.4,a1,i2.2)
        rfile = phyout(1:lphy)//nfs
        open (unit=10,file=rfile,form='unformatted')
        read(10) tmp1,tmp2,tmp3,tmp4,tmp5  &
                 ,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11 &
                 ,tm12,tm13,tm14
        close (10)
      endif
  
      flag=.false.
      if(myrank .eq. 0) flag=.true.
#ifdef SP
      mpe_typ=mpe_single
#else 
      mpe_typ=mpe_double
#endif

      call mpe_bcast(tmp1,nx*my,0,mpe_typ)
      call mpe_bcast(tmp2,nx*my,0,mpe_typ)
      call mpe_bcast(tmp3,nx*my,0,mpe_typ)
      call mpe_bcast(tmp4,nx*my,0,mpe_typ)
      call mpe_bcast(tmp5,nx*my,0,mpe_typ)
      call mpe_bcast(tmp6,nx*my,0,mpe_typ)
      call mpe_bcast(tmp7,nx*my,0,mpe_typ)
      call mpe_bcast(tmp8,nx*my,0,mpe_typ)
      call mpe_bcast(tmp9,nx*my,0,mpe_typ)
      call mpe_bcast(tmp10,nx*my,0,mpe_typ)
      call mpe_bcast(tmp11,nx*my,0,mpe_typ)
      do k=0, 3
        if(k .le. 1)then
          call mpe_bcast(tm12(:,k,:),nx*my,0,mpe_typ)
          call mpe_bcast(tm13(:,k,:),nx*my,0,mpe_typ)
          call mpe_bcast(tm14(:,k,:),nx*my,0,mpe_typ)
        elseif(k .ge. 2) then
          call mpe_bcast(tm14(:,k,:),nx*my,0,mpe_typ)
        endif
      enddo

      do jj=1,jlistnum
        j=jlist1(jj)
        i=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) then
          call reducepickr_sp (tmp1(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp2(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp3(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp4(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp5(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp6(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp7(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp8(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp9(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp10(1,j),nxdef(j),nx,1)
          call reducepickr_sp (tmp11(1,j),nxdef(j),nx,1)
          do k=0, 3
            if(k .le. 1)then
              call reducepickr_sp (tm12(1,k,j),nxdef(j),nx,1)
              call reducepickr_sp (tm13(1,k,j),nxdef(j),nx,1)
              call reducepickr_sp (tm14(1,k,j),nxdef(j),nx,1)
            elseif(k .ge. 2) then
              call reducepickr_sp (tm14(1,k,j),nxdef(j),nx,1)
            endif
          enddo
        endif


        do ii=1,nxj
          i=nxjstart(j)+ii-1
          sitcc(ii,jj)   =tmp1(i,j)
          sithc(ii,jj)   =tmp2(i,j)
          engwac(ii,jj)  =tmp3(i,j)
          sc(ii,jj)      =tmp4(i,j)
          saltwac(ii,jj) =tmp5(i,j)
          wtfns(ii,jj)   =tmp6(i,j)
          wsfns(ii,jj)   =tmp7(i,j)
          grndcapc(ii,jj)=tmp8(i,j)
          grndhflx(ii,jj)=tmp9(i,j)
          grndflux(ii,jj)=tmp10(i,j)
          obswtb(ii,jj)  =tmp11(i,j)

          do k = 0, 3
            if(k .le. 1) then
              zsi(ii,jj,k)  =tm12(i,k,j)
              silw(ii,jj,k) =tm13(i,k,j)
              tsnic(ii,jj,k)=tm14(i,k,j)
            endif
            if(k .ge. 2) then
              tsnic(ii,jj,k)=tm14(i,k,j)
            endif
          enddo

        enddo
      enddo


      deallocate(tmp1,tmp2,tmp3,tmp4,tmp5)   
      deallocate(tmp6,tmp7,tmp8,tmp9,tmp10)   
      deallocate(tmp11,tmp12,tmp13,tmp14)   
      deallocate( tm12, tm13, tm14)   

      return
      end subroutine readrerun_sitgrid1

!-------------------------------------------------------------	  
      subroutine readrerun_sitgrid2(itau)

      use index
      use mpe
      use const
      use mod_sitgrid

      real, dimension(:,:,:), allocatable :: tm11,tm12,tm13
      real, dimension(:,:,:), allocatable :: tm14,tm15,tm16

      
      integer itau,lphy
      character nfs*10, rfile*80
      logical flag
      integer i,j,k,ii,jj,nxj

      allocate (tm11(nx,0:lkvl+1,my))
      allocate (tm12(nx,0:lkvl+1,my))
      allocate (tm13(nx,0:lkvl+1,my))
      allocate (tm14(nx,0:lkvl+1,my))
      allocate (tm15(nx,0:lkvl+1,my))
      allocate (tm16(nx,0:lkvl+1,my))

      call chlen (phyout,80,lphy)
  
      if(myrank .eq. 0) then
        print *,'ready to read rerun_sitgrid2'
        write (nfs,801) itau,'_',2
  801   format('sit',i4.4,a1,i2.2)
        rfile = phyout(1:lphy)//nfs
        open (unit=10,file=rfile,form='unformatted')
        read(10) tm11,tm12,tm13,tm14,tm15,tm16
        close (10)
      endif

      if(myrank .eq. 0) then
        print *,'read, tm11(914,:,265)=',tm11(914,:,265)
      endif
  
      flag=.false.
      if(myrank .eq. 0) flag=.true. 

      do k = 0, lkvl+1
        call mpe_bcast(tm11(:,k,:),nx*my,0,mpe_double)
        call mpe_bcast(tm12(:,k,:),nx*my,0,mpe_double)
        call mpe_bcast(tm13(:,k,:),nx*my,0,mpe_double)
        call mpe_bcast(tm14(:,k,:),nx*my,0,mpe_double)
        call mpe_bcast(tm15(:,k,:),nx*my,0,mpe_double)
        call mpe_bcast(tm16(:,k,:),nx*my,0,mpe_double)

        do jj=1,jlistnum
          j=jlist1(jj)
          i=nxjstart(j)
          nxj=nxdef_2d(j)
          if( lreduce.eq.1 ) then
            call reducepickr (tm11(1,k,j),nxdef(j),nx,1)
            call reducepickr (tm12(1,k,j),nxdef(j),nx,1)
            call reducepickr (tm13(1,k,j),nxdef(j),nx,1)
            call reducepickr (tm14(1,k,j),nxdef(j),nx,1)
            call reducepickr (tm15(1,k,j),nxdef(j),nx,1)
            call reducepickr (tm16(1,k,j),nxdef(j),nx,1)
          endif
  
          do ii = 1, nxj
            i=nxjstart(j)+ii-1
            sitwt(ii,jj,k)  = tm11(i,k,j)
            sitwu(ii,jj,k)  = tm12(i,k,j)
            sitwv(ii,jj,k)  = tm13(i,k,j)
            sitww(ii,jj,k)  = tm14(i,k,j)
            sitws(ii,jj,k)  = tm15(i,k,j)
            sitwtke(ii,jj,k)  = tm16(i,k,j)
          enddo

        enddo !end jj
      enddo   !end k

      deallocate  (tm11)
      deallocate  (tm12)
      deallocate  (tm13)
      deallocate  (tm14)
      deallocate  (tm15)
      deallocate  (tm16)

      end subroutine readrerun_sitgrid2

!-------------------------------------------------------------	  
      subroutine readrerun_sitgrid3(itau)

      use index
      use mpe
      use const
      use mod_sitgrid

      real, dimension(:,:,:), allocatable :: tm11,tm12
  
   
      integer itau,lphy
      character nfs*10, rfile*80
      logical flag
      integer i,j,k,ii,jj,nxj


      allocate (tm11(nx,0:lkvl+1,my))
      allocate (tm12(nx,0:lkvl+1,my))
      call chlen (phyout,80,lphy)
  
      if(myrank .eq. 0) then
        print *,'ready to read rerun_sitgrid3'
        write (nfs,801) itau,'_',3
  801   format('sit',i4.4,a1,i2.2)
        rfile = phyout(1:lphy)//nfs
        open (unit=10,file=rfile,form='unformatted')
        read(10) tm11,tm12
        close (10)
      endif
  
      flag=.false.
      if(myrank .eq. 0) flag=.true. 
      do k = 0, lkvl+1
        call mpe_bcast(tm11(:,k,:),nx*my,0,mpe_double)    
        call mpe_bcast(tm12(:,k,:),nx*my,0,mpe_double)    


        do jj =1, jlistnum
          j=jlist1(jj)
          i=nxjstart(j)
          nxj=nxdef_2d(j)
          if(lreduce .eq. 1) then 
            call reducepickr (tm11(1,k,j),nxdef(j),nx,1)
            call reducepickr (tm12(1,k,j),nxdef(j),nx,1)
          endif

          do ii = 1, nxj
            i=nxjstart(j)+ii-1
            wtfn(ii,jj,k) = tm11(i,k,j)
            wsfn(ii,jj,k) = tm12(i,k,j)
          enddo
        enddo

      enddo

      deallocate (tm11)
      deallocate (tm12)
  
      end subroutine readrerun_sitgrid3

!-------------------------------------------------------------	  
      subroutine readpre6hr_sit(nx,my,itaup,ifilin,idtg,ggdef)

      use index
      use mpe
      use rank
      use mod_sitgrid

      implicit none
      integer nx,my
      integer itaup,nxmy,i,j,k,ii,jj,istat,nxj
      character*60 ifilin
      character*4 ggdef
      integer*8 idtg,idtg2
      character*6 typ
  
  
      real, dimension(:,:), allocatable:: tm1,tm2,tm3   &
                                         ,tm4,tm5,tm6   

  
      allocate (tm1(nx,my))
      allocate (tm2(nx,my))
      allocate (tm3(nx,my))
      allocate (tm4(nx,my))
      allocate (tm5(nx,my))
      allocate (tm6(nx,my))
  
      call dtgfix12(idtg,idtg2,-itaup)
  
  
      nxmy=nx*my
      do k= 0, lkvl+1 
  
        write( typ, '(i3.3,a3)' ) k,'SWT'           !!sit wt
        call syslbl_r (typ,idtg2,itaup,ggdef)
        call dmsread(nx,my,nxmy,'H',ifilin,tm1,istat)

        write( typ, '(i3.3,a3)' ) k,'SWU'           !!sit wu
        call syslbl_r (typ,idtg2,itaup,ggdef)
        call dmsread(nx,my,nxmy,'H',ifilin,tm2,istat)

        write( typ, '(i3.3,a3)' ) k,'SWV'           !!sit wv
        call syslbl_r (typ,idtg2,itaup,ggdef)
        call dmsread(nx,my,nxmy,'H',ifilin,tm3,istat)

        write( typ, '(i3.3,a3)' ) k,'SWW'           !!sit ww
        call syslbl_r (typ,idtg2,itaup,ggdef)
        call dmsread(nx,my,nxmy,'H',ifilin,tm4,istat)

        write( typ, '(i3.3,a3)' ) k,'SWS'           !!sit ws
        call syslbl_r (typ,idtg2,itaup,ggdef)
        call dmsread(nx,my,nxmy,'H',ifilin,tm5,istat)

        write( typ, '(i3.3,a3)' ) k,'TKE'           !!sit wtke
        call syslbl_r (typ,idtg2,itaup,ggdef)
        call dmsread(nx,my,nxmy,'H',ifilin,tm6,istat)


        do jj =1, jlistnum
          j=jlist1(jj)
          i=nxjstart(j)
          nxj=nxdef_2d(j)
          if(lreduce.eq.1 ) then
            call reducepickr (tm1(1,j),nxdef(j),nx,1)
            call reducepickr (tm2(1,j),nxdef(j),nx,1)
            call reducepickr (tm3(1,j),nxdef(j),nx,1)
            call reducepickr (tm4(1,j),nxdef(j),nx,1)
            call reducepickr (tm5(1,j),nxdef(j),nx,1)
            call reducepickr (tm6(1,j),nxdef(j),nx,1)
          endif
          do ii = 1, nxj
            i=nxjstart(j)+ii-1
            sitwt(ii,jj,k)  = MERGE(tm1(i,j),sitwt(ii,jj,k),tm1(i,j).NE.xmissing)
            sitwu(ii,jj,k)  = MERGE(tm2(i,j),sitwu(ii,jj,k),tm2(i,j).NE.xmissing)
            sitwv(ii,jj,k)  = MERGE(tm3(i,j),sitwv(ii,jj,k),tm3(i,j).NE.xmissing)
            sitww(ii,jj,k)  = MERGE(tm4(i,j),sitww(ii,jj,k),tm4(i,j).NE.xmissing)
            sitws(ii,jj,k)  = MERGE(tm5(i,j),sitws(ii,jj,k),tm5(i,j).NE.xmissing)
            sitwtke(ii,jj,k)= MERGE(tm6(i,j),sitwtke(ii,jj,k),tm6(i,j).NE.xmissing)
          enddo
        enddo

      enddo   !end k

      deallocate  (tm1)
      deallocate  (tm2)
      deallocate  (tm3)
      deallocate  (tm4)
      deallocate  (tm5)
      deallocate  (tm6)

      end subroutine readpre6hr_sit


      SUBROUTINE outtseadiffSIT24(nx,my,my_max,ratioSIT,dt24,itau,idtg,ggdef)

      use mpe
      use index
      use const,             only: kflag,RTYPE,outdms,outgrb2,ihdgo,ihdgo2
      use mod_sitgrid,       only: tseadiffSIT24
      use mod_grb2_param,only: wrt_grb2_v2

      implicit none

      integer   nx,my,my_max,itau
      real      dt24
      real(kind=RTYPE) wrk(nxp,my_max),glob(nx,my),wrk2(nxp,my_max)
      real      ratioSIT(nxp,my_max)
      integer*8 idtg
      character*4  ggdef
      integer   imax,jmax,lenc,j,nxj,i,istat,jj

      imax=nx
      jmax=my
      lenc= imax*jmax
!
      do jj = 1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
         wrk(i,jj)=tseadiffSIT24(i,jj)!/dt24
         wrk2(i,jj)=ratioSIT(i,jj)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('w0002f',idtg,itau,ggdef)
      if(outdms.gt.0) &
      call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if( outgrb2==1.and.myrank==0 ) then
      ihdgo2 = ihdgo
      call wrt_grb2_v2(itau,10,3,199,2,168,0,1,glob)
      endif
      tseadiffSIT24=0.

      call unify_reduceintp(nx,my,my_max,wrk2,glob)
      call syslbl_w ('w00002',idtg,itau,ggdef)
      if(outdms.gt.0) &
      call dmswrit(imax,jmax,lenc,kflag,glob,istat)

      END SUBROUTINE outtseadiffSIT24
