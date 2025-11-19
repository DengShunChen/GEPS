module gfs_cpl
contains
subroutine gfs_cpl_init(compid)
  use param,     only: nlon_glb=>nx, my_max, my
  use index,     only: nlon=>nxp, nlat=>jlen, & 
                       mylon=>row_rank, mylat=>col_rank, &
                       jlistnum, jlist1, nxdef_2d, nxjstart, nxdef
  use rank,      only: MPI_COMM_atm, myid=>myrank
  use cpl_rank,  only: id_gocn
  use cpl_gsmap, only: cpl_gsmap_init
  use cpl_attr,  only: cpl_attr_init
  use cpl_smat,  only: cpl_smat_init
  implicit none

  integer, intent(in) :: compid
  integer :: j, jj, rank_root, seg_strt(nlat-1), seg_leng(nlat-1), ierr

  rank_root = 0
  do jj = 1, jlistnum
    j = jlist1(jj)
    seg_strt(jj) = mylon*nlon + (j-1)*nlon_glb + 1 
    seg_leng(jj) = nlon
  enddo

  call cpl_gsmap_init(myid, nlat-1, seg_strt, seg_leng, rank_root, mpi_comm_atm, compid, "datm GSMap:")
  call cpl_smat_init(myid, rank_root, compid, id_gocn, mpi_comm_atm, "./rmp_tco2timcom.nc")
  call cpl_attr_init(compid, id_gocn, MPI_COMM_atm)
end subroutine gfs_cpl_init

subroutine gfs_cpl_send2gocn(compid, u10m, v10m, t02m, q02m, &
                                     pslv, swup, swdn, lwdn, &
                                     rain, snow, tgfs)
  use param,        only: nlon_glb=>nx, nlat_glb=>my
  use index,        only: nlon=>nxp, nlat=>jlen, jlistnum, jlist1, &
                          nxdef_2d, nxjstart, mylon=>row_rank
  use rank,         only: myid=>myrank
  use cpl_rank,     only: id_gocn
  use cpl_attr,     only: AttrVect_importRAttr, gfs_export_attr, gfsExport_rList
  use cpl_sendrecv, only: cpl_send
  use const       , only: RTYPE
  implicit none

  integer, intent(in) :: compid
  real(kind=8), dimension(nlon,nlat), intent(in)  :: u10m, v10m, &
                                                     t02m, q02m, &
                                                     pslv, swup, &
                                                     swdn, lwdn, &
                                                     rain, snow, tgfs

  real(kind=8), dimension(nlon_glb, nlat_glb,11) :: glob_var

  real(kind=RTYPE), dimension(nlon,nlat)             :: ocnwrk1
  real(kind=RTYPE), dimension(nlon_glb, nlat_glb,11) :: ocnwrk2

  integer :: cnt, i, ii, j, jj, nxj 

  ocnwrk1=0.
  ocnwrk2=0.
  glob_var=0.

  ocnwrk1 = u10m
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,1))
  glob_var(:,:,1) = ocnwrk2(:,:,1)

  ocnwrk1 = v10m
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,2))
  glob_var(:,:,2) = ocnwrk2(:,:,2)

  ocnwrk1 = t02m
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,3))
  glob_var(:,:,3) = ocnwrk2(:,:,3)

  ocnwrk1 = q02m
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,4))
  glob_var(:,:,4) = ocnwrk2(:,:,4)

  ocnwrk1 = pslv
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,5))
  glob_var(:,:,5) = ocnwrk2(:,:,5)

  ocnwrk1 = swup
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,6))
  glob_var(:,:,6) = ocnwrk2(:,:,6)

  ocnwrk1 = swdn
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,7))
  glob_var(:,:,7) = ocnwrk2(:,:,7)

  ocnwrk1 = lwdn
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,8))
  glob_var(:,:,8) = ocnwrk2(:,:,8)

  ocnwrk1 = rain
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,9))
  glob_var(:,:,9) = ocnwrk2(:,:,9)

  ocnwrk1 = snow
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,10))
  glob_var(:,:,10) = ocnwrk2(:,:,10)

  ocnwrk1 = tgfs
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, ocnwrk1, ocnwrk2(:,:,11))
  glob_var(:,:,11) = ocnwrk2(:,:,11)
  
  cnt = 0
  do jj = 1, jlistnum
    j = jlist1(jj)
    ii = mylon*nlon + 1 !nxjstart(j)
    nxj = nlon      !nxdef_2d(j)
    do i = ii, ii+nxj-1
      cnt = cnt + 1
      gfs_export_attr%rAttr(1,cnt) = glob_var(i,j,1)
      gfs_export_attr%rAttr(2,cnt) = glob_var(i,j,2)
      gfs_export_attr%rAttr(3,cnt) = glob_var(i,j,3)
      gfs_export_attr%rAttr(4,cnt) = glob_var(i,j,4)
      gfs_export_attr%rAttr(5,cnt) = glob_var(i,j,5)*100.d0
      gfs_export_attr%rAttr(6,cnt) = glob_var(i,j,6)
      gfs_export_attr%rAttr(7,cnt) = glob_var(i,j,7)
      gfs_export_attr%rAttr(8,cnt) = glob_var(i,j,8)
      gfs_export_attr%rAttr(9,cnt) = glob_var(i,j,9)
      gfs_export_attr%rAttr(10,cnt) = glob_var(i,j,10)
      gfs_export_attr%rAttr(11,cnt) = glob_var(i,j,11) - 273.15d0
    end do
  enddo

  call cpl_send(compid, id_gocn, gfs_export_attr, myid)

end subroutine gfs_cpl_send2gocn

!jwhwu 20241007
!subroutine gfs_cpl_recv4gocn(compid, mask_lnd, tgfs, ssufs, ssvfs)
  subroutine gfs_cpl_recv4gocn(compid, mask_lnd, mask_ice, tgfs, ssufs, ssvfs, ifrac, icedp, snodp)
!jwhwu
  use param,        only: nx, my, my_max
  use index,        only: nxp, jlistnum, jlist1, &
                          nxdef_2d, nxjstart, lreduce, nxdef
  use mpe,          only: mpe_double
  use rank,         only: myrank
  use cpl_rank,     only: id_gocn
  use cpl_attr,     only: recv_AV=>gfs_recv_gocn_Attr
  use cpl_sendrecv, only: cpl_recv
  use const       , only: RTYPE
  implicit none

  integer, intent(in) :: compid
!jwhwu 20241007
! logical, intent(inout) :: mask_lnd(nxp, my_max)
  logical, intent(inout) :: mask_lnd(nxp, my_max), mask_ice(nxp, my_max)
!jwhwu
  real, intent(inout) :: tgfs(nxp, my_max), ssufs(nxp,my_max), ssvfs(nxp,my_max), ifrac(nxp,my_max), icedp(nxp,my_max), snodp(nxp,my_max)
  real, dimension(nx, my) :: tg_glb, sst_glb, ssu_glb, ssv_glb, ifrac_glb, icedp_glb, snodp_glb
  real, dimension(nx, my_max) :: sst_nxj
  real, dimension(nxp, my_max) :: SST, SSU, SSV, CICE, ZICE, ZSNO
  real(kind=RTYPE) :: ocnwrk3(nxp, my_max)
  real(kind=RTYPE), dimension(nx, my) :: ocnwrk4
  real(kind=RTYPE), dimension(nx, my_max) :: ocnwrk5
  real(kind=RTYPE), dimension(nxp, my_max) :: ocnwrk6

  
  integer :: cnt, i, ii, j, jj, nxj

  call cpl_recv(compid, id_gocn, recv_AV, myrank)

  sst_nxj = 0.0
  sst_glb = 0.0
  cnt = 0
  SST=0.
  SSU=0.
  SSV=0.
  CICE=0.
  ZICE=0.
  ZSNO=0.
  ifrac_glb=0.           
  icedp_glb=0.           
  snodp_glb=0.           
  ssu_glb=0.
  ssv_glb=0.
  ssufs=0.
  ssvfs=0.
  ocnwrk3=0.
  ocnwrk4=0.
  ocnwrk5=0.
  ocnwrk6=0.

  do jj = 1, jlistnum
    j = jlist1(jj)
    nxj = nxp  !nxdef_2d(j)
    do i = 1, nxj
      cnt = cnt + 1
      SST(i,jj) = recv_AV%rAttr(1,cnt)
      SSU(i,jj) = recv_AV%rAttr(2,cnt)
      SSV(i,jj) = recv_AV%rAttr(3,cnt)
      CICE(i,jj) = recv_AV%rAttr(4,cnt)
      ZICE(i,jj) = recv_AV%rAttr(5,cnt)
      ZSNO(i,jj) = recv_AV%rAttr(6,cnt)
    end do
  end do
!      write(*,*) 'SSUmax=', maxval(SSU), 'myrank=', myrank
  ocnwrk3 = tgfs
  call unify_reduceintp(nx, my, my_max, ocnwrk3, ocnwrk4)
  tg_glb = ocnwrk4

  ocnwrk6 = SST
  call mpe2d_unify(ocnwrk4, ocnwrk6, .true.)
  sst_glb = ocnwrk4

  ocnwrk6 = SSU
  call mpe2d_unify(ocnwrk4, ocnwrk6, .true.)
  ssu_glb = ocnwrk4

  ocnwrk6 = SSV
  call mpe2d_unify(ocnwrk4, ocnwrk6, .true.)
  ssv_glb = ocnwrk4

  ocnwrk6 = CICE
  call mpe2d_unify(ocnwrk4, ocnwrk6, .true.)
  ifrac_glb = ocnwrk4

  ocnwrk6 = ZICE
  call mpe2d_unify(ocnwrk4, ocnwrk6, .true.)
  icedp_glb = ocnwrk4

  ocnwrk6 = ZSNO
  call mpe2d_unify(ocnwrk4, ocnwrk6, .true.)
  snodp_glb = ocnwrk4

!      if(myrank .eq. 0) write(*,*) 'ssu_glbmax=', maxval(ssu_glb)
  !call mpe2d_unify_nx(sst_nxj, SST)
  !call mpe2d_unify_my(sst_glb, sst_nxj)
  !do j = 1, my
  !  do i = 1, nx
  !    if(sst_glb(i,j).gt.271.0) tg_glb(i,j) = sst_glb(i,j) 
  !  end do
  !end do

  !call unify_reducepick(nx, my, my_max, tg_glb, tgfs)
  do jj = 1, jlistnum
    j   = jlist1(jj)
    ii  = nxjstart(j)
    nxj = nxdef_2d(j)
    if( lreduce.eq.1 ) then
      call reducepick(sst_glb(1,j), nxdef(j), nx, 1)
      call reducepick(ssu_glb(1,j), nxdef(j), nx, 1)
      call reducepick(ssv_glb(1,j), nxdef(j), nx, 1)
      call reducepick(ifrac_glb(1,j), nxdef(j), nx, 1)
      call reducepick(icedp_glb(1,j), nxdef(j), nx, 1)
      call reducepick(snodp_glb(1,j), nxdef(j), nx, 1)
    end if
      
    do i = 1, nxj
      if(.not.mask_lnd(i,jj).and.sst_glb(ii,j).gt.271.0) then
        tgfs(i,jj)  = sst_glb(ii,j) 
        ssufs(i,jj) = ssu_glb(ii,j)
        ssvfs(i,jj) = ssv_glb(ii,j)
      end if
      if(.not.mask_lnd(i,jj)) then
        ifrac(i,jj) = ifrac_glb(ii,j)
        icedp(i,jj) = icedp_glb(ii,j)
        snodp(i,jj) = snodp_glb(ii,j)
      endif
      ii = ii + 1 
    end do
  end do

end subroutine gfs_cpl_recv4gocn

end module gfs_cpl
