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
  call cpl_smat_init(myid, rank_root, compid, id_gocn, mpi_comm_atm, "/nwpr/gfs/xb124/TCo_TIMCOM/remap/TCo383_to_tcm1536_consv_modify.nc")
  call cpl_attr_init(compid, id_gocn, MPI_COMM_atm)
end subroutine gfs_cpl_init

subroutine gfs_cpl_send2gocn(compid, taux, tauy, lath, senh, &
                                     swnt, lwnt, evap, prec, &
                                     u10m, v10m, t10m, q10m, pslv, tgfs)
  use param,        only: nlon_glb=>nx, nlat_glb=>my
  use index,        only: nlon=>nxp, nlat=>jlen, jlistnum, jlist1, &
                          nxdef_2d, nxjstart, mylon=>row_rank
  use rank,         only: myid=>myrank
  use cpl_rank,     only: id_gocn
  use cpl_attr,     only: AttrVect_importRAttr, gfs_export_attr, gfsExport_rList
  use cpl_sendrecv, only: cpl_send
  implicit none

  integer, intent(in) :: compid
  real(kind=8), dimension(nlon,nlon), intent(in)  :: taux, tauy, &
                                                     lath, senh, &
                                                     swnt, lwnt, &
                                                     evap, prec, &
                                                     u10m, v10m, &
                                                     t10m, q10m, &
                                                     pslv, tgfs

  real(kind=8), dimension(nlon_glb, nlat_glb,12) :: glob_var
  integer :: cnt, i, ii, j, jj, nxj 

  call unify_reduceintp(nlon_glb, nlat_glb, nlat, tgfs, glob_var(:,:,1))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, u10m, glob_var(:,:,2))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, v10m, glob_var(:,:,3))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, t10m, glob_var(:,:,4))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, q10m, glob_var(:,:,5))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, pslv, glob_var(:,:,6))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, prec, glob_var(:,:,7))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, lath, glob_var(:,:,8))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, senh, glob_var(:,:,9))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, lwnt, glob_var(:,:,10))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, swnt, glob_var(:,:,11))
  call unify_reduceintp(nlon_glb, nlat_glb, nlat, evap, glob_var(:,:,12))
  
  cnt = 0
  do jj = 1, jlistnum
    j = jlist1(jj)
    ii = mylon*nlon + 1 !nxjstart(j)
    nxj = nlon      !nxdef_2d(j)
    do i = ii, ii+nxj-1
      cnt = cnt + 1
      gfs_export_attr%rAttr(1,cnt)  = glob_var(i,j,1) - 273.15d0
      gfs_export_attr%rAttr(2,cnt)  = glob_var(i,j,2)
      gfs_export_attr%rAttr(3,cnt)  = glob_var(i,j,3)
      gfs_export_attr%rAttr(4,cnt)  = glob_var(i,j,4)
      gfs_export_attr%rAttr(5,cnt)  = glob_var(i,j,5)
      gfs_export_attr%rAttr(6,cnt)  = glob_var(i,j,6)*100.d0
      gfs_export_attr%rAttr(7,cnt)  = glob_var(i,j,7)
      gfs_export_attr%rAttr(8,cnt)  =-glob_var(i,j,8)
      gfs_export_attr%rAttr(9,cnt)  =-glob_var(i,j,9)
      gfs_export_attr%rAttr(10,cnt) = glob_var(i,j,10)
      gfs_export_attr%rAttr(11,cnt) = glob_var(i,j,11)
      gfs_export_attr%rAttr(12,cnt) = glob_var(i,j,12)
    end do
  enddo

  call cpl_send(compid, id_gocn, gfs_export_attr, myid)

end subroutine gfs_cpl_send2gocn

subroutine gfs_cpl_recv4gocn(compid, mask_ocn, tgfs)
  use param,        only: nx, my, my_max
  use index,        only: nxp, jlistnum, jlist1, &
                          nxdef_2d, nxjstart, lreduce, nxdef
  use mpe,          only: mpe_double
  use rank,         only: myrank
  use cpl_rank,     only: id_gocn
  use cpl_attr,     only: recv_AV=>gfs_recv_gocn_Attr
  use cpl_sendrecv, only: cpl_recv
  implicit none

  integer, intent(in) :: compid
  logical, intent(inout) :: mask_ocn(nxp, my_max)
  real, intent(inout) :: tgfs(nxp, my_max)
  real, dimension(nx, my) :: tg_glb, sst_glb
  real, dimension(nx, my_max) :: sst_nxj
  real, dimension(nxp, my_max) :: SST, SSU, SSV
  
  integer :: cnt, i, ii, j, jj, nxj

  call cpl_recv(compid, id_gocn, recv_AV, myrank)

  sst_nxj = 0.0
  sst_glb = 0.0
  cnt = 0
  do jj = 1, jlistnum
    j = jlist1(jj)
    nxj = nxp  !nxdef_2d(j)
    do i = 1, nxj
      cnt = cnt + 1
      SST(i,jj) = recv_AV%rAttr(1,cnt)
      SSU(i,jj) = recv_AV%rAttr(2,cnt)
      SSV(i,jj) = recv_AV%rAttr(3,cnt)
    end do
  end do

  call unify_reduceintp(nx, my, my_max, tgfs, tg_glb)
  call mpe2d_unify(sst_glb, SST, .true.)
  !call mpe2d_unify_nx(sst_nxj, SST)
  !call mpe2d_unify_my(sst_glb, sst_nxj)
  do j = 1, my
    do i = 1, nx
      if(sst_glb(i,j).gt.271.0) tg_glb(i,j) = sst_glb(i,j)
    end do
  end do

  !call unify_reducepick(nx, my, my_max, tg_glb, tgfs)
  do jj = 1, jlistnum
    j   = jlist1(jj)
    ii  = nxjstart(j)
    nxj = nxdef_2d(j)
    if( lreduce.eq.1 ) call reducepick(tg_glb(1,j), nxdef(j), nx, 1)
    do i = 1, nxj
      tgfs(i,jj) = tg_glb(ii,j) 
      ii = ii + 1 
    end do
  end do

end subroutine gfs_cpl_recv4gocn

end module gfs_cpl
