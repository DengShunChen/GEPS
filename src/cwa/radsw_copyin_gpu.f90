module radsw_copyin_gpu
      implicit none

      private

      public copyin_radsw_datatb_gpu

contains

      subroutine copyin_radsw_datatb_gpu(async_id)
      implicit none
      integer, intent(in) :: async_id

      call copyin_radsw_sflux(async_id)
      call copyin_radsw_ref(async_id)
      call copyin_radsw_cldprtb(async_id)
      call copyin_radsw_kgb16(async_id)
      call copyin_radsw_kgb17(async_id)
      call copyin_radsw_kgb18(async_id)
      call copyin_radsw_kgb19(async_id)
      call copyin_radsw_kgb20(async_id)
      call copyin_radsw_kgb21(async_id)
      call copyin_radsw_kgb22(async_id)
      call copyin_radsw_kgb23(async_id)
      call copyin_radsw_kgb24(async_id)
      call copyin_radsw_kgb25(async_id)
      call copyin_radsw_kgb26(async_id)
      call copyin_radsw_kgb27(async_id)
      call copyin_radsw_kgb28(async_id)
      call copyin_radsw_kgb29(async_id)
      !$acc wait(async_id)

      end subroutine

      subroutine copyin_radsw_sflux(async_id)
      use module_radsw_sflux
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(strrat, specwt, layreffr, ix1, ix2, ibx, sfluxref01, &
      !$acc&      sfluxref02, sfluxref03) async(async_id)

      end subroutine


      subroutine copyin_radsw_ref(async_id)
      use module_radsw_ref
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(preflog, tref) async(async_id)

      end subroutine


      subroutine copyin_radsw_cldprtb(async_id)
      use module_radsw_cldprtb
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(extliq1, ssaliq1, asyliq1, extice2, ssaice2, &
      !$acc&      asyice2, extice3, ssaice3, asyice3, fdlice3, abari, bbari, &
      !$acc&      cbari, dbari, ebari, fbari, b0r, b0s, b1s, c0r, c0s) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb16(async_id)
      use module_radsw_kgb16
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb17(async_id)
      use module_radsw_kgb17
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb18(async_id)
      use module_radsw_kgb18
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb19(async_id)
      use module_radsw_kgb19
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb20(async_id)
      use module_radsw_kgb20
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      absch4) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb21(async_id)
      use module_radsw_kgb21
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, selfref, forref, absb) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb22(async_id)
      use module_radsw_kgb22
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb23(async_id)
      use module_radsw_kgb23
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, rayl, selfref, forref) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb24(async_id)
      use module_radsw_kgb24
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      abso3a, abso3b, rayla, raylb) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb25(async_id)
      use module_radsw_kgb25
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, rayl, abso3a, abso3b) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb26(async_id)
      use module_radsw_kgb26
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(rayl) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb27(async_id)
      use module_radsw_kgb27
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, rayl) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb28(async_id)
      use module_radsw_kgb28
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb) async(async_id)

      end subroutine


      subroutine copyin_radsw_kgb29(async_id)
      use module_radsw_kgb29
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      absh2o, absco2) async(async_id)

      end subroutine



end module radsw_copyin_gpu
