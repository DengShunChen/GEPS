#ifdef RSM
      subroutine wrte_data(idtg,fhour,nx,my,lsoil,lev,ncld,temp_gfs,spfh_gfs   &
     &              ,clwr_gfs,rain_gfs,qice_gfs,snow_gfs,grpl_gfs    &
     &              ,ozon_gfs,geop_gfs,u_gfs,v_gfs    &
     &              ,tg_gfs,smc_gfs,snr_gfs,stc_gfs,cice_gfs   &
!yj2019     
     &              ,terr_gfs,slmsk_gfs                        &
     &              )
!
      use rank, only : myrank
!      
      integer(kind=8)  :: idtg
      integer nsig, kh, ndig
      integer lsoil
      integer ncld
      real fhour
      character cfhour*16,cform*40,cidtg*12

      dimension temp_gfs(nx*my,lev),spfh_gfs(nx*my,lev)      &
     & ,clwr_gfs(nx*my,lev),rain_gfs(nx*my,lev),qice_gfs(nx*my,lev) &
     & ,snow_gfs(nx*my,lev),grpl_gfs(nx*my,lev),ozon_gfs(nx*my,lev) &
     & ,geop_gfs(nx*my),u_gfs(nx*my,lev),v_gfs(nx*my,lev)    &
     & ,tg_gfs(nx*my),smc_gfs(nx*my,lsoil),snr_gfs(nx*my)        &
     & ,stc_gfs(nx*my,lsoil),cice_gfs(nx*my)                     &
!yj2019     
     & ,terr_gfs(nx*my),slmsk_gfs(nx*my)                     

!from call      if(myrank.eq.0)then
        print *,' start global wrte_data'
        nsig=52

        kh=nint(fhour)
        ndig=max(log10(kh+0.5)+1.,3.)
        write(cform,'("(i",i1,".",i1,")")') ndig,ndig
        write(cfhour,cform) kh
        write(cidtg,'(I12.12)') idtg

        open(nsig,file='rsm_data_'//cidtg//'.f'//cfhour,status='unknown', &
            form='unformatted',iostat=ios)

        write(nsig) fhour

        write(nsig) temp_gfs

        write(nsig) spfh_gfs

        write(nsig) clwr_gfs
        if ( ncld .ge. 7 ) then
          write(nsig) rain_gfs
          write(nsig) qice_gfs
          write(nsig) snow_gfs
          write(nsig) grpl_gfs
        endif

        write(nsig) ozon_gfs

        write(nsig) geop_gfs
!yj2019
        write(nsig) terr_gfs

        write(nsig) u_gfs

        write(nsig) v_gfs

        write(nsig) tg_gfs

        write(nsig) smc_gfs

        write(nsig) snr_gfs

        write(nsig) stc_gfs

        write(nsig) cice_gfs
!yj2019
        write(nsig) slmsk_gfs

        print *,' global wrte_data'
        close(nsig)
!from call      endif
!
      return
      end
#endif
