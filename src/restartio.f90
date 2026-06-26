subroutine restartio( iaction )
!iaction = 0 : read
!        = 1 : write
use rank   ,only: myrank
use const  ,only: tau, hours, cwbout, phyout, hfiltx, alphax, pdryi,      &
                  itimestep
use spec   ,only: plnow , vornow, divnow, temnow,  &
                  plold , vorold, divold, temold, trefs
use phygrid,only: land,ice,ocean,alb,z0,tgclim,gwclim,ustar,tstar,qstar,  &
                  raincu,rainlp,totalp,raintot,curate,plcl,cumtop,snr,    &
                  tg,tg_diff,tg_ocn,gwr,gwet,fpsp,fpsp1,hflux,qflux,ss,   &
                  rs,ctot,tsflw,cosz,o3l,ftp,fqp,ftp1,fqp1,asl,atl,dtrad, &
                  deltaq,cnvwr,cnvcr,dtcup,ducup,dvcup,dtshl,dushl,dvshl, &
#ifdef TIMCOMCPL
                  ssu,ssv,                                                &
#endif
                  dtlsp
use grid   ,only: sgeo,std,tsave,qt,pdiff,xy
use noah   ,only: smc,stc,slc,sld,rld,sncover,sndepth,cice,xtice,zice,    &
                  sfemis,sfalb
use raddiag,only: asol,asold,olr,clds
use radn   ,only: sdec,cdec,slag,solcon,solhr,rsolhr
!=====
integer iaction ,ierr
!
character::rfile*255, ctau*7, ccore*4, cmdxx*256
integer::itau, itaui
logical::flg
integer::i,ios
!=====
if(iaction == 0 )then

        write (ctau,800) itaui
  800   format(i7.7)
!
        write(ccore,'(i4.4)') myrank
        rfile = trim(cwbout)//'cwbout_'//ctau//'/'//ccore
!
        if(myrank .eq. 0) then
          print*, 'read restart at tau, rfile = ', taui, trim(rfile)
        endif
        flg=.false.
        i=200+myrank
        open (i,file=trim(rfile),form='unformatted',status='old' &
             ,iostat=ios)
!
        if(ios .ne.0) then
          print *,'restartio : mpe_broadcast error !! , rfile=',trim(rfile)
          call mpe_finalize
          stop
        endif

        read(i) vornow
        read(i) divnow
        read(i) temnow
        read(i) vorold
        read(i) divold
        read(i) temold
        read(i) trefs
        read(i) plnow
        read(i) plold
        close(i)
!
        rfile = trim(phyout)//'phyout_'//ctau//'/'//ccore

        i=200+myrank
        open (i,file=trim(rfile),form='unformatted',status='old',iostat=ios)
        if (ios .ne. 0) then
          print *,'restartio : open rfile error !! , rfile=',trim(rfile)
          call mpe_finalize
          stop
        endif
  !
        read(i) land
        read(i) ocean
        read(i) ice
        read(i) alb
        read(i) z0
        read(i) tgclim
        read(i) gwclim
        read(i) sgeo
        read(i) canopy
        read(i) ustar
        read(i) tstar
        read(i) qstar
        read(i) raincu
        read(i) rainlp
        read(i) totalp
        read(i) raintot
        read(i) curate
        read(i) plcl
        read(i) cumtop
        read(i) snr
        read(i) sncover
        read(i) sndepth
        read(i) tg
        read(i) tg_diff
        read(i) tg_ocn
        read(i) gwr
        read(i) gwet
        read(i) cice
        read(i) xtice
        read(i) zice
        read(i) fpsp
        read(i) fpsp1
  !
        read(i) hflux
        read(i) qflux
        read(i) ss
        read(i) rs
        read(i) asol
        read(i) olr
        read(i) sld
        read(i) rld
        read(i) asold
  !jwhwu 202003
        read(i) sfemis
        read(i) sfalb
        read(i) std
  !jwhwu
        read(i) ctot
        read(i) clds
        read(i) pdiff
        read(i) tsave
  !jwhwu 201702 | add
        read(i) tsflw
  !jwhwu 201702 | add
        read(i) cosz
        read(i) slag,sdec,cdec,solcon,solhr,rsolhr, &
  ! overcome round off problem in restart
               tau,hours,itimestep,xy

        read(i) o3l
        read(i) ftp
        read(i) fqp
        read(i) ftp1
        read(i) fqp1
  !jwhwu 201702 ! add
        read(i) asl
        read(i) atl
        read(i) dtrad
  !jwhwu 201910 ! add
        read(i) deltaq
        read(i) cnvwr
        read(i) cnvcr
  !jwhwu 202208 ! add
        read(i) dtcup
        read(i) ducup
        read(i) dvcup
        read(i) dtshl
        read(i) dushl
        read(i) dvshl
        read(i) dtlsp
        read(i) hfiltx
        read(i) alphax
        read(i) pdryi
  !
        read(i) qt
  !       read(i) qp   ! not need
        read(i) smc
        read(i) stc
        read(i) slc
#ifdef TIMCOMCPL
        read(i) ssu
        read(i) ssv
#endif
        close(i)
!
elseif(iaction == 1 )then

        if(myrank.eq.0) write(*,*) 'Write out restart files in itau=',itau
        write (ctau, 801) itau
801     format(i7.7)
        write (ccore, '(i4.4)') myrank
!
        cmdxx =' '
        cmdxx = 'mkdir -p '//trim(cwbout)//'cwbout_'//ctau
        call system(trim(cmdxx))
        rfile = trim(cwbout)//'cwbout_'//ctau//'/'//ccore
        i = 200 + myrank
        open (i, file=rfile, form='unformatted')
        write(i) vornow
        write(i) divnow
        write(i) temnow
        write(i) vorold
        write(i) divold
        write(i) temold
        write(i) trefs
        write(i) plnow
        write(i) plold
        close(i)

        cmdxx = 'mkdir -p '//trim(phyout)//'phyout_'//ctau
        call system(trim(cmdxx))
        rfile = trim(phyout)//'phyout_'//ctau//'/'//ccore
!         if (myrank .eq. 0) &
!           open (unit=10,file=rfile,form='unformatted')
        i = 200 + myrank
        open (i, file=rfile, form='unformatted')
!jwhwu 202004 avoid undefined values.
        write (i) land
        write (i) ocean
        write (i) ice
        write (i) alb
        write (i) z0
        write (i) tgclim
        write (i) gwclim
        write (i) sgeo
        write (i) canopy
        write (i) ustar
        write (i) tstar
        write (i) qstar
        write (i) raincu
        write (i) rainlp
        write (i) totalp
        write (i) raintot
        write (i) curate
        write (i) plcl
        write (i) cumtop
        write (i) snr
        write (i) sncover
        write (i) sndepth
        write (i) tg
        write (i) tg_diff
        write (i) tg_ocn
        write (i) gwr
        write (i) gwet
        write (i) cice
        write (i) xtice
        write (i) zice
        write (i) fpsp
        write (i) fpsp1
!
        write (i) hflux
        write (i) qflux
        write (i) ss
        write (i) rs
        write (i) asol
        write (i) olr
        write (i) sld
        write (i) rld
        write (i) asold
!jwhwu 202003
        write (i) sfemis
        write (i) sfalb
        write (i) std
!jwhwu
        write (i) ctot
        write (i) clds
        write (i) pdiff
        write (i) tsave
!jwhwu 201702 | add
        write (i) tsflw
!jwhwu 201702 | add
        write (i) cosz
        write (i) slag, sdec, cdec, solcon, solhr, rsolhr, &
           ! overcome round off problem in restart
           tau, hours, itimestep, xy

        write (i) o3l
        write (i) ftp
        write (i) fqp
        write (i) ftp1
        write (i) fqp1
!jwhwu 201702 ! add
        write (i) asl
        write (i) atl
        write (i) dtrad
!jwhwu 201910 ! add
        write (i) deltaq
        write (i) cnvwr
        write (i) cnvcr
!jwhwu 202208 ! add
        write (i) dtcup
        write (i) ducup
        write (i) dvcup
        write (i) dtshl
        write (i) dushl
        write (i) dvshl
        write (i) dtlsp
        write (i) hfiltx
        write (i) alphax
        write (i) pdryi
!
        write (i) qt
!         write(i) qp  !  not need
        write (i) smc
        write (i) stc
        write (i) slc
#ifdef TIMCOMCPL
        write (i) ssu
        write (i) ssv
#endif
        close (i)

else
  if(myrank==0)print*,'In restartIO no any action. '
  continue
endif

return
end subroutine
