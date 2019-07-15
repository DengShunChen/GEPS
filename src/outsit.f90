#define myrank_check 50
#define i_check 428
#define jj_check 2

      SUBROUTINE sitout(nx,my,itau,ifilout,idtg,num  &
                       ,whtlev,ggdef)
      
      use rank
      use mpe
      use index
      use mod_sitgrid 

      implicit none
      integer nx,my,itau,num
      real whtlev(num)
   

      character*60 ifilout
      integer*8 idtg
      character*4 ggdef
      integer lenc,n,k,jj,j,nxj,ii,istat
      character*26 lrec
      character*6 typ
      real, dimension(:,:), allocatable:: tm1,tm2,tm3   &
                                        ,tm4,tm5,tm6    &
                                        ,tm10,tm11,tm12

  
      allocate(tm1(nx,my))
      allocate(tm2(nx,my))
      allocate(tm3(nx,my))
      allocate(tm4(nx,my))
      allocate(tm5(nx,my))
      allocate(tm6(nx,my))
      allocate(tm10(nx,my))
      allocate(tm11(nx,my))
      allocate(tm12(nx,my))

      lenc= nx*my
  
      if (myrank.eq.0) print *,'in sitout, lkvl=',lkvl

      do 40 n=1, num
        do 30 k = 0, lkvl+1
            do jj =1, jlistnum
              j=jlist1(jj)
              nxj=nxdef_2d(j)
              do ii = 1, nxj
                tm1(ii,j)  = sitwt(ii,jj,k)
                tm2(ii,j)  = sitwu(ii,jj,k)
                tm3(ii,j)  = sitwv(ii,jj,k)
                tm4(ii,j)  = sitww(ii,jj,k)
                tm5(ii,j)  = sitws(ii,jj,k)
                tm6(ii,j)  = sitwtke(ii,jj,k)
                tm10(ii,j)  = obswt(ii,jj,k)
                tm11(ii,j)  = wkm(ii,jj,k)
                tm12(ii,j)  = wkh(ii,jj,k)
              enddo
            enddo
      
            call mpe_unify(tm1,nx,my,2,mpe_double)
            call mpe_unify(tm2,nx,my,2,mpe_double)
            call mpe_unify(tm3,nx,my,2,mpe_double)
            call mpe_unify(tm4,nx,my,2,mpe_double)
            call mpe_unify(tm5,nx,my,2,mpe_double)
            call mpe_unify(tm6,nx,my,2,mpe_double)
            call mpe_unify(tm10,nx,my,2,mpe_double)
            call mpe_unify(tm11,nx,my,2,mpe_double)
            call mpe_unify(tm12,nx,my,2,mpe_double)

            if( lreduce.eq.1 ) then
              do jj = 1, jlistnum
                j=jlist1(jj)
                call reduceintp (tm1(1,j),nxdef(j),nx,1)
                call reduceintp (tm2(1,j),nxdef(j),nx,1)
                call reduceintp (tm3(1,j),nxdef(j),nx,1)
                call reduceintp (tm4(1,j),nxdef(j),nx,1)
                call reduceintp (tm5(1,j),nxdef(j),nx,1)
                call reduceintp (tm6(1,j),nxdef(j),nx,1)
                call reduceintp (tm10(1,j),nxdef(j),nx,1)
                call reduceintp (tm11(1,j),nxdef(j),nx,1)
                call reduceintp (tm12(1,j),nxdef(j),nx,1)
              enddo
              call mpe_unify(tm1,nx,my,5,mpe_double)
              call mpe_unify(tm2,nx,my,5,mpe_double)
              call mpe_unify(tm3,nx,my,5,mpe_double)
              call mpe_unify(tm4,nx,my,5,mpe_double)
              call mpe_unify(tm5,nx,my,5,mpe_double)
              call mpe_unify(tm6,nx,my,5,mpe_double)
              call mpe_unify(tm10,nx,my,5,mpe_double)
              call mpe_unify(tm11,nx,my,5,mpe_double)
              call mpe_unify(tm12,nx,my,5,mpe_double)
            endif



            write( typ, '(i3.3,a3)' ) k,'SWT'         !!sit wt
            call syslbl (typ,idtg,itau,ggdef,lrec)
!            if( lreduce.eq.1 ) call reduceintp (tm1,nxdef,nx,my)
            call dmswrit(nx,my,lrec,lenc,'H',ifilout,tm1,istat)

            write( typ, '(i3.3,a3)' ) k,'SWU'         !!sit wu
            call syslbl (typ,idtg,itau,ggdef,lrec)
!            if( lreduce.eq.1 ) call reduceintp (tm2,nxdef,nx,my)
            call dmswrit(nx,my,lrec,lenc,'H',ifilout,tm2,istat)

            write( typ, '(i3.3,a3)' ) k,'SWV'         !!sit wv
            call syslbl (typ,idtg,itau,ggdef,lrec)
!            if( lreduce.eq.1 ) call reduceintp (tm3,nxdef,nx,my)
            call dmswrit(nx,my,lrec,lenc,'H',ifilout,tm3,istat)

            write( typ, '(i3.3,a3)' ) k,'SWW'         !!sit ww
            call syslbl (typ,idtg,itau,ggdef,lrec)
!            if( lreduce.eq.1 ) call reduceintp (tm4,nxdef,nx,my)
            call dmswrit(nx,my,lrec,lenc,'H',ifilout,tm4,istat)

            write( typ, '(i3.3,a3)' ) k,'SWS'         !!sit ws
            call syslbl (typ,idtg,itau,ggdef,lrec)
!            if( lreduce.eq.1 ) call reduceintp (tm5,nxdef,nx,my)
            call dmswrit(nx,my,lrec,lenc,'H',ifilout,tm5,istat)

            write( typ, '(i3.3,a3)' ) k,'TKE'         !!sit wtke
            call syslbl (typ,idtg,itau,ggdef,lrec)
!            if( lreduce.eq.1 ) call reduceintp (tm6,nxdef,nx,my)
            call dmswrit(nx,my,lrec,lenc,'H',ifilout,tm6,istat)

            write( typ, '(i3.3,a3)' ) k,'OWT'        !!sit obswt
            call syslbl (typ,idtg,itau,ggdef,lrec)
!            if( lreduce.eq.1 ) call reduceintp (tm10,nxdef,nx,my)
            call dmswrit(nx,my,lrec,lenc,'H',ifilout,tm10,istat)

            write( typ, '(i3.3,a3)' ) k,'WKM'
            call syslbl (typ,idtg,itau,ggdef,lrec)
!            if( lreduce.eq.1 ) call reduceintp (tm11,nxdef,nx,my)
            call dmswrit(nx,my,lrec,lenc,'H',ifilout,tm11,istat)
            
            write( typ, '(i3.3,a3)' ) k,'WKH'
            call syslbl (typ,idtg,itau,ggdef,lrec)
!            if( lreduce.eq.1 ) call reduceintp (tm12,nxdef,nx,my)
            call dmswrit(nx,my,lrec,lenc,'H',ifilout,tm12,istat)


   30   continue
   40 continue
      
      deallocate(tm1)
      deallocate(tm2)
      deallocate(tm3)
      deallocate(tm4)
      deallocate(tm5)
      deallocate(tm6)
      deallocate(tm10)
      deallocate(tm11)
      deallocate(tm12)
  
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
      SUBROUTINE writesitmean(nx,my,lkvl,ifilout,itau,idtg,ggdef)

      use rank
      use mpe
      use index
      use mod_sitgrid,     only:sitwttau,sitwstau,sitwutau,sitwvtau &
                               ,dtsittau
      use mod_sit_control,only: xmissing

      implicit none

      integer nx,lkvl,my,itau


      character*60 ifilout
      integer*8 idtg
      character*4 ggdef
      integer lenc,k,jj,j,nxj,ii,istat
      character*26 ihdg
      character*6 lrec
      real tm1(nx,my),tm2(nx,my),tm3(nx,my),tm4(nx,my)

      tm1=0.
!      tm2=0.
!      tm3=0.
!      tm4=0.

      lenc= nx*my


      if (myrank.eq.0) print *,'in writesitmean lkvl=',lkvl
      do k = 0, lkvl+1
        do jj =1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do ii = 1, nxj
            if(dtsittau .ne. 0.) then
              tm1(ii,j)  = sitwttau(ii,jj,k)/dtsittau
!              tm2(i,j)  = sitwstau(i,jj,k)/dtsittau(i,jj,k)
!              tm3(i,j)  = sitwutau(i,jj,k)/dtsittau(i,jj,k)
!              tm4(i,j)  = sitwvtau(i,jj,k)/dtsittau(i,jj,k)
            endif
           enddo
        enddo

        call mpe_unify(tm1,nx,my,2,mpe_double)
!        call mpe_unify(tm2,nx,my,2,mpe_double)
!        call mpe_unify(tm3,nx,my,2,mpe_double)
!        call mpe_unify(tm4,nx,my,2,mpe_double)

        if( lreduce.eq.1 ) then
          do jj = 1, jlistnum
            j=jlist1(jj)
            call reduceintp (tm1(1,j),nxdef(j),nx,1)
          enddo
          call mpe_unify(tm1,nx,my,5,mpe_double)
        endif


        write( lrec, '(i3.3,a3)' ) k,'WTT'
        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm1,nxdef,nx,my)
        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm1,istat)

!        write( lrec, '(i3.3,a3)' ) k,'WST'
!        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm2,nxdef,nx,my)
!        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm2,istat)

!        write( lrec, '(i3.3,a3)' ) k,'WUT'
!        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm3,nxdef,nx,my)
!        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm3,istat)

!        write( lrec, '(i3.3,a3)' ) k,'WVT'
!        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm4,nxdef,nx,my)
!        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm4,istat)
      enddo

      sitwttau=0.
!      sitwstau=0.
!      sitwutau=0.
!      sitwvtau=0.
      dtsittau=0.

      END SUBROUTINE writesitmean



!--------------------------------------------------------------------------
      SUBROUTINE outsit24(nx,my,lkvl,ifilout,itau,idtg,ggdef)
      
      use rank
      use mpe
      use index
      use mod_sitgrid,    only: sitwt24,sitws24,sitwu24,sitwv24 &
                               ,dtsit24,wtfn0,wsfn0,obswt,sitwt
      use mod_sit_control,only: xmissing

      implicit none

      integer nx,lkvl,my,itau
   
      character*60 ifilout
      integer*8 idtg
      character*4 ggdef
      integer lenc,k,jj,j,nxj,ii,istat
      character*26 ihdg
      character*6 lrec
      real tm1(nx,my),tm2(nx,my),tm3(nx,my),tm4(nx,my)

      tm1=xmissing
      tm2=xmissing
      tm3=xmissing
!      tm4=0.

      lenc= nx*my

  
      if (myrank.eq.0) print *,'in outsit24 lkvl=',lkvl

      do k = 0, lkvl+1
        do jj =1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do ii = 1, nxj
            if(dtsit24 .ne. 0.) then
              tm1(ii,j)  = sitwt24(ii,jj,k)/dtsit24
!              tm2(ii,j)  = wtfn0(ii,jj,k)/dtsit24
              if(obswt(ii,jj,k) .ge. 271. ) then
                tm2(ii,j) = (obswt(ii,jj,k)-sitwt(ii,jj,k))/dtsit24
              else
                tm2(ii,j) = 0.
              endif
              tm3(ii,j)  = wsfn0(ii,jj,k)/dtsit24
!              tm2(i,j)  = sitws24(i,jj,k)/dtsit24(i,jj,k)
!              tm3(i,j)  = sitwu24(i,jj,k)/dtsit24(i,jj,k)
!              tm4(i,j)  = sitwv24(i,jj,k)/dtsit24(i,jj,k)
            endif
          enddo
        enddo
      
        call mpe_unify(tm1,nx,my,2,mpe_double)
        call mpe_unify(tm2,nx,my,2,mpe_double)
        call mpe_unify(tm3,nx,my,2,mpe_double)
!        call mpe_unify(tm4,nx,my,2,mpe_double)
          
        if( lreduce.eq.1 ) then
          do jj = 1, jlistnum
            j=jlist1(jj)
            call reduceintp (tm1(1,j),nxdef_2d(j),nx,1)
            call reduceintp (tm2(1,j),nxdef_2d(j),nx,1)
            call reduceintp (tm3(1,j),nxdef_2d(j),nx,1)
          enddo
          call mpe_unify(tm1,nx,my,5,mpe_double)
          call mpe_unify(tm2,nx,my,5,mpe_double)
          call mpe_unify(tm3,nx,my,5,mpe_double)
        endif


        write( lrec, '(i3.3,a3)' ) k,'WTF'
        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm1,nxdef,nx,my)
        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm1,istat)

        write( lrec, '(i3.3,a3)' ) k,'TFN'
        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm2,nxdef,nx,my)
        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm2,istat)

        write( lrec, '(i3.3,a3)' ) k,'SFN'
        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm3,nxdef,nx,my)
        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm3,istat)

!        write( lrec, '(i3.3,a3)' ) k,'WSF'
!        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm2,nxdef,nx,my)
!        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm2,istat)

!        write( lrec, '(i3.3,a3)' ) k,'WUF'
!        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm3,nxdef,nx,my)
!        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm3,istat)

!        write( lrec, '(i3.3,a3)' ) k,'WVF'
!        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm4,nxdef,nx,my)
!        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm4,istat)
      enddo
     
      sitwt24=0.
      wtfn0=0.
      wsfn0=0.
!      sitws24=0.
!      sitwu24=0.
!      sitwv24=0.
      dtsit24=0. 
  
      end subroutine outsit24


!--------------------------------------------------------------------------
      SUBROUTINE outsitmon(nx,my,lkvl,ifilout,itau,idtg,ggdef)
      
      use rank
      use mpe
      use index
      use mod_sitgrid,        only: dtsitmon,wtfn,wtfns,wsfn,wsfns
      use mod_sit_control,    only: xmissing
  
      implicit none

      integer nx,my,lkvl,itau
  


      real dtmon
      character*60 ifilout
      integer*8 idtg
      character*4 ggdef
      integer lenc,i,j,k,ii,jj,nxj,istat
      character*26 ihdg
      character*6 lrec
      real tm1(nx,my),tm2(nx,my),tm3(nx,my),tm4(nx,my)

      tm1=xmissing
      tm2=xmissing
      tm3=xmissing
      tm4=xmissing

      lenc= nx*my

  
      if (myrank.eq.0) print *,'in outsitmon lkvl=',lkvl


      do k = 0, lkvl+1
        do jj =1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do ii = 1, nxj
            tm1(ii,j)  = wtfn(ii,jj,k)/dtsitmon
            tm2(ii,j)  = wsfn(ii,jj,k)/dtsitmon
            if(k .eq. 0) then
              tm3(ii,j)  = wtfns(ii,jj)/dtsitmon
              tm4(ii,j)  = wsfns(ii,jj)/dtsitmon
            endif
            if((ii .eq. 914) .AND. (jj .eq. 3) .AND. (myrank .eq. 72)) then
              print *, 'wtfn(',ii,',3,',k,')=',wtfn(ii,jj,k) &
                     , 'wsfn(',ii,',3,',k,')=',wsfn(ii,jj,k) &
                     , 'wtfns(',ii,',3,)=',wtfns(ii,jj) &
                     , 'wsfns(',ii,',3,)=',wsfns(ii,jj) &
                     , 'dtsitmon=',dtsitmon
            endif
          enddo
        enddo
      
        call mpe_unify(tm1,nx,my,2,mpe_double)
        call mpe_unify(tm2,nx,my,2,mpe_double)
        if(k.eq.0)then
          call mpe_unify(tm3,nx,my,2,mpe_double)
          call mpe_unify(tm4,nx,my,2,mpe_double)
        endif
          
        if( lreduce.eq.1 ) then
          do jj = 1, jlistnum
           j=jlist1(jj)
           call reduceintp (tm1(1,j),nxdef_2d(j),nx,1)
           call reduceintp (tm2(1,j),nxdef_2d(j),nx,1)
           if(k.eq.0)then
             call reduceintp (tm3(1,j),nxdef_2d(j),nx,1)
             call reduceintp (tm4(1,j),nxdef_2d(j),nx,1)
           endif
          enddo
          call mpe_unify(tm1,nx,my,5,mpe_double)
          call mpe_unify(tm2,nx,my,5,mpe_double)
          if(k .eq. 0)then
            call mpe_unify(tm3,nx,my,5,mpe_double)
            call mpe_unify(tm4,nx,my,5,mpe_double)
          endif
        endif

        write( lrec, '(i3.3,a3)' ) k,'TFM'
        call syslbl (lrec,idtg,itau,ggdef,ihdg)
        if(myrank .eq. 0) print *, 'idtg=',idtg,',itau=',itau,',ggdef=',ggdef
!        if( lreduce.eq.1 ) call reduceintp (tm1,nxdef,nx,my)
        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm1,istat)

        write( lrec, '(i3.3,a3)' ) k,'SFM'
        call syslbl (lrec,idtg,itau,ggdef,ihdg)
!        if( lreduce.eq.1 ) call reduceintp (tm2,nxdef,nx,my)
        call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm2,istat)

        if(k .eq. 0) then
          write( lrec, '(i3.3,a3)' ) k,'TFS'
          call syslbl (lrec,idtg,itau,ggdef,ihdg)
!         if( lreduce.eq.1 ) call reduceintp (tm3,nxdef,nx,my)
          call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm3,istat)

          write( lrec, '(i3.3,a3)' ) k,'SFS'
          call syslbl (lrec,idtg,itau,ggdef,ihdg)
!         if( lreduce.eq.1 ) call reduceintp (tm4,nxdef,nx,my)
          call dmswrit(nx,my,ihdg,lenc,'H',ifilout,tm4,istat)
        endif

      enddo


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


      real, dimension(:,:), allocatable ::     &
                   tmp1,tmp2,tmp3,tmp4,tmp5    &
                  ,tmp6,tmp7,tmp8,tmp9,tmp10   &
                  ,tmp11,tmp12,tmp13,tmp14
      real, dimension(:,:,:), allocatable::    &
                   tm12, tm13, tm14


      integer itau,lphy
      character nfs*10, rfile*80
      integer i,j,k,ii,jj,nxj

      allocate(tmp1(nx,my),tmp2(nx,my),tmp3(nx,my),tmp4(nx,my)   &
               ,tmp5(nx,my),tmp6(nx,my),tmp7(nx,my),tmp8(nx,my)   &
               ,tmp9(nx,my),tmp10(nx,my),tmp11(nx,my))
      allocate(tmp12(nx,my),tmp13(nx,my),tmp14(nx,my))
      allocate(tm12(nx,0:1,my),tm13(nx,0:1,my),tm14(nx,0:3,my))




      do k = 0, 3  
        do jj =1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do ii = 1, nxj
            if(k .eq. 0) then
              tmp1(ii,j)  = sitcc(ii,jj)
              tmp2(ii,j)  = sithc(ii,jj)
              tmp3(ii,j)  = engwac(ii,jj)
              tmp4(ii,j)  = sc(ii,jj)
              tmp5(ii,j)  = saltwac(ii,jj)
              tmp6(ii,j)  = wtfns(ii,jj)
              tmp7(ii,j)  = wsfns(ii,jj)
              tmp8(ii,j)  = grndcapc(ii,jj)
              tmp9(ii,j)  = grndhflx(ii,jj)
              tmp10(ii,j) = grndflux(ii,jj)
              tmp11(ii,j) = obswtb(ii,jj)
              tmp12(ii,j) = zsi(ii,jj,k)
              tmp13(ii,j) = silw(ii,jj,k)
              tmp14(ii,j) = tsnic(ii,jj,k)
            endif
            if(k .eq. 1) then
              tmp12(ii,j) = zsi(ii,jj,k)
              tmp13(ii,j) = silw(ii,jj,k)
              tmp14(ii,j) = tsnic(ii,jj,k)
            endif
            if(k .ge. 2) then
              tmp14(ii,j) = tsnic(ii,jj,k)
            endif

          enddo
        enddo

        if(k .eq. 0) then
          call mpe_unify(tmp1,nx,my,2,mpe_double)  
          call mpe_unify(tmp2,nx,my,2,mpe_double)  
          call mpe_unify(tmp3,nx,my,2,mpe_double)  
          call mpe_unify(tmp4,nx,my,2,mpe_double)  
          call mpe_unify(tmp5,nx,my,2,mpe_double)  
          call mpe_unify(tmp6,nx,my,2,mpe_double)  
          call mpe_unify(tmp7,nx,my,2,mpe_double)  
          call mpe_unify(tmp8,nx,my,2,mpe_double)  
          call mpe_unify(tmp9,nx,my,2,mpe_double)  
          call mpe_unify(tmp10,nx,my,2,mpe_double)  
          call mpe_unify(tmp11,nx,my,2,mpe_double)  
          call mpe_unify(tmp12,nx,my,2,mpe_double)  
          call mpe_unify(tmp13,nx,my,2,mpe_double)
          call mpe_unify(tmp14,nx,my,2,mpe_double)
        endif
        if(k .eq. 1) then
          call mpe_unify(tmp12,nx,my,2,mpe_double)  
          call mpe_unify(tmp13,nx,my,2,mpe_double)
          call mpe_unify(tmp14,nx,my,2,mpe_double)
        endif
        if(k .ge. 2) then
          call mpe_unify(tmp14,nx,my,2,mpe_double)
        endif


        if( lreduce.eq.1 ) then
          do jj = 1, jlistnum
            j=jlist1(jj)
            if(k .eq. 0) then
              call reduceintp (tmp1(1,j),nxdef(j),nx,1)
              call reduceintp (tmp2(1,j),nxdef(j),nx,1)
              call reduceintp (tmp3(1,j),nxdef(j),nx,1)
              call reduceintp (tmp4(1,j),nxdef(j),nx,1)
              call reduceintp (tmp5(1,j),nxdef(j),nx,1)
              call reduceintp (tmp6(1,j),nxdef(j),nx,1)
              call reduceintp (tmp7(1,j),nxdef(j),nx,1)
              call reduceintp (tmp8(1,j),nxdef(j),nx,1)
              call reduceintp (tmp9(1,j),nxdef(j),nx,1)
              call reduceintp (tmp10(1,j),nxdef(j),nx,1)
              call reduceintp (tmp11(1,j),nxdef(j),nx,1)
              call reduceintp (tmp12(1,j),nxdef(j),nx,1)
              call reduceintp (tmp13(1,j),nxdef(j),nx,1)
              call reduceintp (tmp14(1,j),nxdef(j),nx,1)
            endif
            if(k .eq. 1) then
              call reduceintp (tmp12(1,j),nxdef(j),nx,1)
              call reduceintp (tmp13(1,j),nxdef(j),nx,1)
              call reduceintp (tmp14(1,j),nxdef(j),nx,1)
            endif
            if(k .eq. 2) then
              call reduceintp (tmp14(1,j),nxdef(j),nx,1)
            endif
          enddo

          call mpe_unify(tmp1,nx,my,5,mpe_double)
          call mpe_unify(tmp2,nx,my,5,mpe_double)
          call mpe_unify(tmp3,nx,my,5,mpe_double)
          call mpe_unify(tmp4,nx,my,5,mpe_double)
          call mpe_unify(tmp5,nx,my,5,mpe_double)
          call mpe_unify(tmp6,nx,my,5,mpe_double)
          call mpe_unify(tmp7,nx,my,5,mpe_double)
          call mpe_unify(tmp8,nx,my,5,mpe_double)
          call mpe_unify(tmp9,nx,my,5,mpe_double)
          call mpe_unify(tmp10,nx,my,5,mpe_double)
          call mpe_unify(tmp11,nx,my,5,mpe_double)
          call mpe_unify(tmp12,nx,my,5,mpe_double)
          call mpe_unify(tmp13,nx,my,5,mpe_double)
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

  
      real, dimension(:,:), allocatable::       &
                  tmp11,tmp12,tmp13,tmp14,tmp15,tmp16

      real, dimension(:,:,:), allocatable::       &
                tm11,tm12,tm13,tm14,tm15,tm16
   
      integer itau,lphy
      character nfs*10, rfile*80
      integer i,j,k,ii,jj,nxj

      allocate( tmp11(nx,my),tmp12(nx,my),tmp13(nx,my) &
               ,tmp14(nx,my),tmp15(nx,my),tmp16(nx,my) )
      allocate( tm11(nx,0:lkvl+1,my),tm12(nx,0:lkvl+1,my) &
               ,tm13(nx,0:lkvl+1,my),tm14(nx,0:lkvl+1,my) &
               ,tm15(nx,0:lkvl+1,my),tm16(nx,0:lkvl+1,my) )

  
  
      do k = 0, lkvl+1
        do jj =1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do ii = 1, nxj
            tmp11(ii,j) = sitwt(ii,jj,k)
            tmp12(ii,j) = sitwu(ii,jj,k)
            tmp13(ii,j) = sitwv(ii,jj,k)
            tmp14(ii,j) = sitww(ii,jj,k)
            tmp15(ii,j) = sitws(ii,jj,k)
            tmp16(ii,j) = sitwtke(ii,jj,k)
          enddo
        enddo

        call mpe_unify(tmp11,nx,my,2,mpe_double)  
        call mpe_unify(tmp12,nx,my,2,mpe_double)  
        call mpe_unify(tmp13,nx,my,2,mpe_double)  
        call mpe_unify(tmp14,nx,my,2,mpe_double)
        call mpe_unify(tmp15,nx,my,2,mpe_double)  
        call mpe_unify(tmp16,nx,my,2,mpe_double)  

!        if( lreduce.eq.1 ) then
!          call reduceintp (tmp11,nxdef,nx,my)
!          call reduceintp (tmp12,nxdef,nx,my)
!          call reduceintp (tmp13,nxdef,nx,my)
!          call reduceintp (tmp14,nxdef,nx,my)
!          call reduceintp (tmp15,nxdef,nx,my)
!          call reduceintp (tmp16,nxdef,nx,my)
!        endif

        if( lreduce.eq.1 ) then
          do jj = 1, jlistnum
            j=jlist1(jj)
            call reduceintp (tmp11(1,j),nxdef(j),nx,1)
            call reduceintp (tmp12(1,j),nxdef(j),nx,1)
            call reduceintp (tmp13(1,j),nxdef(j),nx,1)
            call reduceintp (tmp14(1,j),nxdef(j),nx,1)
            call reduceintp (tmp15(1,j),nxdef(j),nx,1)
            call reduceintp (tmp16(1,j),nxdef(j),nx,1)
          enddo
          call mpe_unify(tmp11,nx,my,5,mpe_double)
          call mpe_unify(tmp12,nx,my,5,mpe_double)
          call mpe_unify(tmp13,nx,my,5,mpe_double)
          call mpe_unify(tmp14,nx,my,5,mpe_double)
          call mpe_unify(tmp15,nx,my,5,mpe_double)
          call mpe_unify(tmp16,nx,my,5,mpe_double)
        endif



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

  
      real, dimension(:,:), allocatable:: tmp11,tmp12
      real, dimension(:,:,:), allocatable:: tm11,tm12

      integer itau,lphy
      character nfs*10, rfile*80
      integer i,j,k,ii,jj,nxj

      allocate( tmp11(nx,my),tmp12(nx,my) )
      allocate( tm11(nx,0:lkvl+1,my),tm12(nx,0:lkvl+1,my) )
  
  
      do k = 0, lkvl+1
        do jj =1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do ii = 1, nxj
            tmp11(ii,j) = wtfn(ii,jj,k)
            tmp12(ii,j) = wsfn(ii,jj,k)
          enddo
        enddo

        call mpe_unify(tmp11,nx,my,2,mpe_double)  
        call mpe_unify(tmp12,nx,my,2,mpe_double)  

!        if( lreduce.eq.1 ) then
!          call reduceintp (tmp11,nxdef,nx,my)
!          call reduceintp (tmp12,nxdef,nx,my)
!        endif

        if( lreduce.eq.1 ) then
          do jj = 1, jlistnum
            j=jlist1(jj)
            call reduceintp (tmp11(1,j),nxdef(j),nx,1)
            call reduceintp (tmp12(1,j),nxdef(j),nx,1)
          enddo
          call mpe_unify(tmp11,nx,my,5,mpe_double)
          call mpe_unify(tmp12,nx,my,5,mpe_double)
        endif

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

      real, dimension(:,:), allocatable::             &
                          tmp1,tmp2,tmp3,tmp4,tmp5    &
                         ,tmp6,tmp7,tmp8,tmp9,tmp10   &
                         ,tmp11,tmp12,tmp13,tmp14
      real, dimension(:,:,:), allocatable:: tm12,tm13,tm14


      integer itau,lphy
      character nfs*10, rfile*80  
      logical flag
      integer i,j,k,ii,jj,nxj
       

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
!      call mpe_broadcast(tmp1,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp2,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp3,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp4,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp5,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp6,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp7,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp8,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp9,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp10,nx*my,flag,mpe_double)
!      call mpe_broadcast(tmp11,nx*my,flag,mpe_double)
!      if(lreduce .eq. 1) then
!        call reducepick (tmp1,nxdef,nx,my)
!        call reducepick (tmp2,nxdef,nx,my)
!        call reducepick (tmp3,nxdef,nx,my)
!        call reducepick (tmp4,nxdef,nx,my)
!        call reducepick (tmp5,nxdef,nx,my)
!        call reducepick (tmp6,nxdef,nx,my)
!        call reducepick (tmp7,nxdef,nx,my)
!        call reducepick (tmp8,nxdef,nx,my)
!        call reducepick (tmp9,nxdef,nx,my)
!        call reducepick (tmp10,nxdef,nx,my)
!        call reducepick (tmp11,nxdef,nx,my)
!      endif

      call mpe_bcast(tmp1,nx*my,0,mpe_double)
      call mpe_bcast(tmp2,nx*my,0,mpe_double)
      call mpe_bcast(tmp3,nx*my,0,mpe_double)
      call mpe_bcast(tmp4,nx*my,0,mpe_double)
      call mpe_bcast(tmp5,nx*my,0,mpe_double)
      call mpe_bcast(tmp6,nx*my,0,mpe_double)
      call mpe_bcast(tmp7,nx*my,0,mpe_double)
      call mpe_bcast(tmp8,nx*my,0,mpe_double)
      call mpe_bcast(tmp9,nx*my,0,mpe_double)
      call mpe_bcast(tmp10,nx*my,0,mpe_double)
      call mpe_bcast(tmp11,nx*my,0,mpe_double)
      do k=0, 3
        if(k .le. 1)then
          call mpe_bcast(tm12(:,k,:),nx*my,0,mpe_double)
          call mpe_bcast(tm13(:,k,:),nx*my,0,mpe_double)
          call mpe_bcast(tm14(:,k,:),nx*my,0,mpe_double)
        elseif(k .ge. 2) then
          call mpe_bcast(tm14(:,k,:),nx*my,0,mpe_double)
        endif
      enddo

      do jj=1,jlistnum
        j=jlist1(jj)
        i=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) then
          call reducepick (tmp1(1,j),nxdef(j),nx,1)
          call reducepick (tmp2(1,j),nxdef(j),nx,1)
          call reducepick (tmp3(1,j),nxdef(j),nx,1)
          call reducepick (tmp4(1,j),nxdef(j),nx,1)
          call reducepick (tmp5(1,j),nxdef(j),nx,1)
          call reducepick (tmp6(1,j),nxdef(j),nx,1)
          call reducepick (tmp7(1,j),nxdef(j),nx,1)
          call reducepick (tmp8(1,j),nxdef(j),nx,1)
          call reducepick (tmp9(1,j),nxdef(j),nx,1)
          call reducepick (tmp10(1,j),nxdef(j),nx,1)
          call reducepick (tmp11(1,j),nxdef(j),nx,1)
          do k=0, 3
            if(k .le. 1)then
              call reducepick (tm12(1,k,j),nxdef(j),nx,1)
              call reducepick (tm13(1,k,j),nxdef(j),nx,1)
              call reducepick (tm14(1,k,j),nxdef(j),nx,1)
            elseif(k .ge. 2) then
              call reducepick (tm14(1,k,j),nxdef(j),nx,1)
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
!        call mpe_broadcast(tm11(:,k,:),nx*my,flag,mpe_double)
!        call mpe_broadcast(tm12(:,k,:),nx*my,flag,mpe_double)
!        call mpe_broadcast(tm13(:,k,:),nx*my,flag,mpe_double)
!        call mpe_broadcast(tm14(:,k,:),nx*my,flag,mpe_double)
!        call mpe_broadcast(tm15(:,k,:),nx*my,flag,mpe_double)
!        call mpe_broadcast(tm16(:,k,:),nx*my,flag,mpe_double)

!        if(lreduce .eq. 1) then 
!          call reducepick (tm11(:,k,:),nxdef,nx,my)
!          call reducepick (tm12(:,k,:),nxdef,nx,my)
!          call reducepick (tm13(:,k,:),nxdef,nx,my)
!          call reducepick (tm14(:,k,:),nxdef,nx,my)
!          call reducepick (tm15(:,k,:),nxdef,nx,my)
!          call reducepick (tm16(:,k,:),nxdef,nx,my)
!        endif  

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
            call reducepick (tm11(1,k,j),nxdef(j),nx,1)
            call reducepick (tm12(1,k,j),nxdef(j),nx,1)
            call reducepick (tm13(1,k,j),nxdef(j),nx,1)
            call reducepick (tm14(1,k,j),nxdef(j),nx,1)
            call reducepick (tm15(1,k,j),nxdef(j),nx,1)
            call reducepick (tm16(1,k,j),nxdef(j),nx,1)
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
!        call mpe_broadcast(tm11(:,k,:),nx*my,flag,mpe_double)    
!        call mpe_broadcast(tm12(:,k,:),nx*my,flag,mpe_double)    
!        if(lreduce .eq. 1) then 
!          call reducepick (tm11(:,k,:),nxdef,nx,my)
!          call reducepick (tm12(:,k,:),nxdef,nx,my)
!        endif

        call mpe_bcast(tm11(:,k,:),nx*my,0,mpe_double)    
        call mpe_bcast(tm12(:,k,:),nx*my,0,mpe_double)    


        do jj =1, jlistnum
          j=jlist1(jj)
          i=nxjstart(j)
          nxj=nxdef_2d(j)
          if(lreduce .eq. 1) then 
            call reducepick (tm11(1,k,j),nxdef(j),nx,1)
            call reducepick (tm12(1,k,j),nxdef(j),nx,1)
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
      character*26 lrec
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
        call syslbl (typ,idtg2,itaup,ggdef,lrec)
        call dmsread(nx,my,lrec,nxmy,'H',ifilin,tm1,istat)
!        if( lreduce.eq.1 ) call reducepick (tm1,nxdef,nx,my)

        write( typ, '(i3.3,a3)' ) k,'SWU'           !!sit wu
        call syslbl (typ,idtg2,itaup,ggdef,lrec)
        call dmsread(nx,my,lrec,nxmy,'H',ifilin,tm2,istat)
!        if( lreduce.eq.1 ) call reducepick (tm2,nxdef,nx,my)

        write( typ, '(i3.3,a3)' ) k,'SWV'           !!sit wv
        call syslbl (typ,idtg2,itaup,ggdef,lrec)
        call dmsread(nx,my,lrec,nxmy,'H',ifilin,tm3,istat)
!        if( lreduce.eq.1 ) call reducepick (tm3,nxdef,nx,my)

        write( typ, '(i3.3,a3)' ) k,'SWW'           !!sit ww
        call syslbl (typ,idtg2,itaup,ggdef,lrec)
        call dmsread(nx,my,lrec,nxmy,'H',ifilin,tm4,istat)
!        if( lreduce.eq.1 ) call reducepick (tm4,nxdef,nx,my)

        write( typ, '(i3.3,a3)' ) k,'SWS'           !!sit ws
        call syslbl (typ,idtg2,itaup,ggdef,lrec)
        call dmsread(nx,my,lrec,nxmy,'H',ifilin,tm5,istat)
!        if( lreduce.eq.1 ) call reducepick (tm5,nxdef,nx,my)

        write( typ, '(i3.3,a3)' ) k,'TKE'           !!sit wtke
        call syslbl (typ,idtg2,itaup,ggdef,lrec)
        call dmsread(nx,my,lrec,nxmy,'H',ifilin,tm6,istat)
!        if( lreduce.eq.1 ) call reducepick (tm6,nxdef,nx,my)


        do jj =1, jlistnum
          j=jlist1(jj)
          i=nxjstart(j)
          nxj=nxdef_2d(j)
          if(lreduce.eq.1 ) then
            call reducepick (tm1(1,j),nxdef(j),nx,1)
            call reducepick (tm2(1,j),nxdef(j),nx,1)
            call reducepick (tm3(1,j),nxdef(j),nx,1)
            call reducepick (tm4(1,j),nxdef(j),nx,1)
            call reducepick (tm5(1,j),nxdef(j),nx,1)
            call reducepick (tm6(1,j),nxdef(j),nx,1)
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

