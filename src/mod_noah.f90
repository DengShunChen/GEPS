      module noah

! modify to f90 by C-H Lee and sort by River Chen in 2015

      use param
      use index

      implicit none

      public

      integer, parameter :: km_soil=4
      integer, parameter :: ntype=9
      integer, parameter :: ngrid=22

      integer isot,ivegsrc

      integer, allocatable,save ::                                   &
              slopetyp(:,:),istyp(:,:),ivegtyp(:,:)

      real, dimension(:,:),allocatable,save ::                       &
              canopy,runoff,rld,sigmaf,sld,gfx,                      &
              zice,cice,xtice, sncover,sndepth,                      &
              shdmax,shdmin,snoalb,sfalb,sfemis                                    

      real    ref(ntype),wlt(ntype),tsat(ntype),                     &
              dfkt(ngrid,ntype),                                     &
              xktk(ngrid,ntype),dfk(ngrid,ntype)

      common /soil/ ref,wlt,tsat, dfkt, xktk,dfk

      real, dimension(:,:,:),allocatable,save :: smc,stc,slc

      contains 

         subroutine allocate_noah_array

           integer  ierr
           integer,parameter :: async_id = 1

           allocate (smc(nxp,km_soil,my_max), stc(nxp,km_soil,my_max), &
                     slc(nxp,km_soil,my_max), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_noah : allocate fail 1 '
               stop
           end if
!
           stc=0.
           smc=0.
           slc=0.
!$acc enter data create(stc, smc, slc) async(async_id)
!
           allocate (canopy(nxp,my_max), runoff(nxp,my_max),rld(nxp,my_max), &
              sigmaf(nxp,my_max), sld(nxp,my_max),gfx(nxp,my_max),           &
              zice(nxp,my_max),cice(nxp,my_max),      &
              xtice(nxp,my_max),                          &
              sncover(nxp,my_max),sndepth(nxp,my_max),&
              shdmax(nxp,my_max),shdmin(nxp,my_max),  &
              snoalb(nxp,my_max),sfalb(nxp,my_max),sfemis(nxp,my_max), stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_noah : allocate fail 2 '
               stop
           end if
!CWB2016
           xtice  =0.0
           sfemis =0.0
           sfalb  =0.0
           zice   =0.0
           sncover=0.0
           gfx    =0.0
           canopy =0.0
           sndepth=0.0
           cice   =0.0
           sld    =0.0
           rld    =0.0
!$acc enter data create(canopy,runoff,rld,sigmaf,sld,gfx,zice &
!$acc& ,cice,xtice,sncover,sndepth,shdmax,shdmin,snoalb,sfalb &
!$acc& ,sfemis) async(async_id)

           allocate (slopetyp(nxp,my_max),istyp(nxp,my_max),&
                     ivegtyp(nxp,my_max), stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_noah : allocate fail 3 '
               stop
           end if
!$acc enter data create(slopetyp,istyp,ivegtyp) async(async_id)

           return

         end subroutine

         subroutine deallocate_noah_array

           deallocate (smc,stc,slc)

           deallocate (canopy,runoff,rld,sigmaf, sld,gfx,zice,cice,xtice,&
                       sncover,sndepth,shdmax,shdmin,snoalb,sfalb,sfemis)

           deallocate (slopetyp,istyp,ivegtyp)

           return

         end subroutine

      end module noah
