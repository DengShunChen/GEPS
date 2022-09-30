      subroutine scatter_spec(work_io,vornow,divnow,temnow,qnow,plnow   &
        ,vorold,divold,temold,qold,plold,dsqgeo,spgeo,trefs             &
        ,lev,ncld,jtrun,jtmax,my,nsize)
!
      use const, only: RTYPE
!
      implicit  none

      integer   lev,ncld,jtrun,jtmax,my,nsize

      real(kind=RTYPE) work_io( ( ( (7+2*ncld)*2*                       &
                lev+8)*jtrun*jtmax+my*lev)*nsize)

      real(kind=RTYPE) vornow(lev,2,jtrun,jtmax)                        &
      , divnow(lev,2,jtrun,jtmax),trefs(lev,2,jtrun,jtmax)              &
      , temnow(lev,2,jtrun,jtmax),qnow(lev*ncld,2,jtrun,jtmax)          &
      , vorold(lev,2,jtrun,jtmax),divold(lev,2,jtrun,jtmax)             &
      , temold(lev,2,jtrun,jtmax),qold(lev*ncld,2,jtrun,jtmax) 
      real(kind=RTYPE) plnow(jtrun,jtmax,2), plold(jtrun,jtmax,2)       &
      , dsqgeo(jtrun,jtmax,2),    spgeo(jtrun,jtmax,2)
!
      integer   ns,len1,len2,len3,len1w,indx

      ns = nsize
      len1=2*lev*jtrun*jtmax
      len2=2*jtrun*jtmax
      len3=my*lev
      len1w=2*lev*ncld*jtrun*jtmax
!
      indx=0
      call mpe_scatter_io(work_io(1           ),vornow,len1,ns)
      indx=indx+len1
      call mpe_scatter_io(work_io(1+  indx *ns),divnow,len1,ns)
      indx=indx+len1
      call mpe_scatter_io(work_io(1+  indx *ns),temnow,len1,ns)
      indx=indx+len1
      call mpe_scatter_io(work_io(1+  indx *ns),qnow  ,len1w,ns)
      indx=indx+len1w
      call mpe_scatter_io(work_io(1+  indx *ns),plnow ,len2,ns)
      indx=indx+len2
      call mpe_scatter_io(work_io(1+  indx *ns),vorold,len1,ns)
      indx=indx+len1
      call mpe_scatter_io(work_io(1+  indx *ns),divold,len1,ns)
      indx=indx+len1
      call mpe_scatter_io(work_io(1+  indx *ns),temold,len1,ns)
      indx=indx+len1
      call mpe_scatter_io(work_io(1+  indx *ns),qold  ,len1w,ns)
      indx=indx+len1w
      call mpe_scatter_io(work_io(1+  indx *ns),plold ,len2,ns)
      indx=indx+len2
      call mpe_scatter_io(work_io(1+  indx *ns),dsqgeo,len2,ns)
      indx=indx+len2
      call mpe_scatter_io(work_io(1+  indx *ns),spgeo ,len2,ns)
      indx=indx+len2
      call mpe_scatter_io(work_io(1+  indx *ns),trefs ,len1,ns)
      indx=indx+len1
!!      call mpe_scatter_io(work_io(1+  indx *ns),qrefs ,len1w,ns)
!!      indx=indx+len1w
!
      return
      end
