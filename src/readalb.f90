      subroutine readalb(bckfile,nx,my,my_max,julian,ggdef,            &
                         alvsf,alvwf,alnsf,alnwf,facsf,facwf)
!----------------------------------------------------------------
!
!  read climt data from data base
!
!----------------------------------------------------------------
!  input :
!         nx      : dimension of e-w direction
!         my      : dimension of s-n direction
!         julian  : julian day
!
!  -- 4 months data
!         alvsfcl : mean vis albedo with strong cosz dependency
!         alvwfcl : mean vis albedo with weak cosz dependency
!         alnsfcl : mean nir albedo with strong cosz dependency
!         alnwfcl : mean nir albedo with weak cosz dependency
!                 
!  output :
! -- time intepolation to julian day
!         alvsf  : mean vis albedo with strong cosz dependency
!         alvwf  : mean vis albedo with weak cosz dependency
!         alnsf  : mean nir albedo with strong cosz dependency
!         alnwf  : mean nir albedo with weak cosz dependency
!  -- annual mean
!         facsf   : fractional coverage with strong cosz dependency
!         facwf   : fractional coverage with weak cosz dependency
!
!  By Mei-Yu Chang, 2014/09/07
!----------------------------------------------------------------
! 
      use rank
      use mpe
      use index

      implicit none

      integer nx,my,my_max,julian,ii

      real  alvsf(nxp,my_max),alvwf(nxp,my_max), &
            alnsf(nxp,my_max),alnwf(nxp,my_max), &
            facsf(nxp,my_max),facwf(nxp,my_max)
! local working array
      integer i,j,jj,nxj,lncrec,mm,nn,k,jul,istat
      real    coef1,coef2
      real  alvsfcl(nxp,my_max,2),alvwfcl(nxp,my_max,2),   &
            alnsfcl(nxp,my_max,2),alnwfcl(nxp,my_max,2)
      real work(nx,my)
!
      character bckfile*80,lrec*26,blnk*1,ggdef*4
      integer   mon(12),mmse(2)
!      data mon/74,166,258,349/
      data mon/15,46,74,105,135,166,196,227,258,288,319,349/
!
      data blnk/' '/
!
      lncrec=nx*my
!
      jul=julian
      if(jul .ge. 366)jul=365
!----------------------------------------------------------------
!  to interpolat linearly based on julian day
!
!-- climate dataset has 12 months

      if(jul .le. mon(1))jul=jul+365
!
      if(jul .gt. mon(12))then
        mmse(1)=12
        mmse(2)=1
        coef1=float(jul-mon(12))/float(365+mon(1)-mon(12))
        coef2=1.-coef1
      else
        do k=2,12
          if(jul .gt. mon(k-1) .and. jul .le. mon(k))then
            mmse(1)=k-1
            mmse(2)=k
            coef1=float(jul-mon(k-1))/float(mon(k)-mon(k-1))
            coef2=1.-coef1
          endif
        enddo
      endif
!----------------------------------------------------------------
! read albedo data
!----------------------------------------------------------------
      do nn=1,2
      mm=mmse(nn)
      write(lrec,31)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,work,istat)
!     call qmax2d(work,1,1,nx,my)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
        do i=1,nxj
          alvsfcl(i,jj,nn)=work(ii,j)
          ii=ii+1
        enddo
      enddo

      write(lrec,32)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,work,istat)
!     call qmax2d(work,1,1,nx,my)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
        do i=1,nxj
          alvwfcl(i,jj,nn)=work(ii,j)
          ii=ii+1
        enddo
      enddo

      write(lrec,33)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,work,istat)
!     call qmax2d(work,1,1,nx,my)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
        do i=1,nxj
          alnsfcl(i,jj,nn)=work(ii,j)
          ii=ii+1
        enddo
      enddo

      write(lrec,34)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,work,istat)
!     call qmax2d(work,1,1,nx,my)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
        do i=1,nxj
          alnwfcl(i,jj,nn)=work(ii,j)
          ii=ii+1
        enddo
      enddo

      enddo ! end of nn
!
!-- facsf (0-100)
!
      write(lrec,35)ggdef
!     write(*,*)'35, lrec=',lrec
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,work,istat)
!     call qmax2d(work,1,1,nx,my)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
      do i=1,nxj
         facsf(i,jj)=work(ii,j)
         ii=ii+1
      enddo
      enddo
!
!-- facwf (0-100)
!
      write(lrec,36)ggdef
!     write(*,*)'35, lrec=',lrec
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,work,istat)
!     call qmax2d(work,1,1,nx,my)
!byl      if( lreduce.eq.1 ) call reducepick (work,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         if( lreduce.eq.1 ) call reducepick (work(1,j),nxdef(j),nx,1)
      do i=1,nxj
         facwf(i,jj)=work(ii,j)
         ii=ii+1
      enddo
      enddo
!----------------------------------------------------------------
  31  format('S0003A','GBCK',a4,4x,i2.2,6x)  ! alvsfcl
  32  format('S0003B','GBCK',a4,4x,i2.2,6x)  ! alvwfcl
  33  format('S0003C','GBCK',a4,4x,i2.2,6x)  ! alnsfcl
  34  format('S0003D','GBCK',a4,4x,i2.2,6x)  ! alnwfcl
  35  format('S0003E','GBCK',a4,12x)         ! facsf
  36  format('S0003F','GBCK',a4,12x)         ! facwf
!----------------------------------------------------------------

! 
!---  time interpolation
!
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
      do i=1,nxj
         alvsf(i,jj)=coef1*alvsfcl(i,jj,1)+coef2*alvsfcl(i,jj,2)
         alvwf(i,jj)=coef1*alvwfcl(i,jj,1)+coef2*alvwfcl(i,jj,2)
         alnsf(i,jj)=coef1*alnsfcl(i,jj,1)+coef2*alnsfcl(i,jj,2)
         alnwf(i,jj)=coef1*alnwfcl(i,jj,1)+coef2*alnwfcl(i,jj,2)
      enddo
      enddo


!     if (myrank .eq. 0) then
!         print *,'*** for readalb.f ***'
!         print *,'julian :',julian
!         print *,'*** alvsf :'
!         call qmax2d(alvsf,1,1,nx,my)
!         print *,'*** alvwf :'
!         call qmax2d(alvwf,1,1,nx,my)
!         print *,'*** alnsf :'
!         call qmax2d(alnsf,1,1,nx,my)
!         print *,'*** alnwf :'
!         call qmax2d(alnwf,1,1,nx,my)
!         print *,'*** facsf :'
!         call qmax2d(facsf,1,1,nx,my)
!         print *,'*** facwf :'
!         call qmax2d(facwf,1,1,nx,my)
!     endif
!
      return
      end
