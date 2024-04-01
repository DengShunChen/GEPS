#if defined(RSM) && defined(CWB_MPMD)
      subroutine send_data(fhour,nx,my,lsoil,lev,ncld,temp_gfs,spfh_gfs &
     &              ,clwr_gfs,rain_gfs,qice_gfs,snow_gfs,grpl_gfs    &
     &              ,ozon_gfs,geop_gfs,u_gfs,v_gfs    &
     &              ,tg_gfs,smc_gfs,snr_gfs,stc_gfs,cice_gfs   &
!yj2019     
     &              ,terr_gfs,slmsk_gfs)
!
      use rank, only : root_rsm,myrank,itag
!      

      integer lsoil
      integer ncld
      real fhour
      dimension temp_gfs(nx*my,lev),spfh_gfs(nx*my,lev)      &
     & ,clwr_gfs(nx*my,lev),rain_gfs(nx*my,lev),qice_gfs(nx*my,lev) &
     & ,snow_gfs(nx*my,lev),grpl_gfs(nx*my,lev),ozon_gfs(nx*my,lev) &
     & ,geop_gfs(nx*my),u_gfs(nx*my,lev),v_gfs(nx*my,lev)    &
     & ,tg_gfs(nx*my),smc_gfs(nx*my,lsoil),snr_gfs(nx*my)        &
     & ,stc_gfs(nx*my,lsoil),cice_gfs(nx*my)                     &
!yj2019     
     & ,terr_gfs(nx*my),slmsk_gfs(nx*my)                      

!from call      if(myrank.eq.0)then
        nxmy=nx*my
        nxmyl=nx*my*lev
        print *,' start global send_data'

        itag=itag+1
        call mpmd_send(fhour,1,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(temp_gfs,nxmyl,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(spfh_gfs,nxmyl,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(clwr_gfs,nxmyl,root_rsm,itag,'R')

        if ( ncld .ge. 7 ) then
        itag=itag+1
        call mpmd_send(rain_gfs,nxmyl,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(qice_gfs,nxmyl,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(snow_gfs,nxmyl,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(grpl_gfs,nxmyl,root_rsm,itag,'R')
        endif

        itag=itag+1
        call mpmd_send(ozon_gfs,nxmyl,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(geop_gfs,nxmy,root_rsm,itag,'R')
!yj2019
        itag=itag+1
        call mpmd_send(terr_gfs,nxmy,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(u_gfs,nxmyl,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(v_gfs,nxmyl,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(tg_gfs,nxmy,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(smc_gfs,nxmy*lsoil,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(snr_gfs,nxmy,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(stc_gfs,nxmy*lsoil,root_rsm,itag,'R')

        itag=itag+1
        call mpmd_send(cice_gfs,nxmy,root_rsm,itag,'R')
!yj2019
        itag=itag+1
        call mpmd_send(slmsk_gfs,nxmy,root_rsm,itag,'R')

        print *,' global send_data: itag=',itag
!from call      endif
!
      return
      end
#endif
