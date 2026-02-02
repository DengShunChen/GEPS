module radlw_copyin_gpu
      implicit none

      private

      public copyin_radlw_datatb_gpu

contains

      subroutine copyin_radlw_datatb_gpu(async_id)
      implicit none
      integer, intent(in) :: async_id

      call copyin_radlw_avplank(async_id)
      call copyin_radlw_ref(async_id)
      call copyin_radlw_cldprlw(async_id)
      call copyin_radlw_kgb01(async_id)
      call copyin_radlw_kgb02(async_id)
      call copyin_radlw_kgb03(async_id)
      call copyin_radlw_kgb04(async_id)
      call copyin_radlw_kgb05(async_id)
      call copyin_radlw_kgb06(async_id)
      call copyin_radlw_kgb07(async_id)
      call copyin_radlw_kgb08(async_id)
      call copyin_radlw_kgb09(async_id)
      call copyin_radlw_kgb10(async_id)
      call copyin_radlw_kgb11(async_id)
      call copyin_radlw_kgb12(async_id)
      call copyin_radlw_kgb13(async_id)
      call copyin_radlw_kgb14(async_id)
      call copyin_radlw_kgb15(async_id)
      call copyin_radlw_kgb16(async_id)

      end subroutine

      subroutine copyin_radlw_avplank(async_id)
      use module_radlw_avplank
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(totplnk) async(async_id)

      end subroutine


      subroutine copyin_radlw_ref(async_id)
      use module_radlw_ref
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(pref, preflog, tref, chi_mls) async(async_id)

      end subroutine


      subroutine copyin_radlw_cldprlw(async_id)
      use module_radlw_cldprlw
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(ipat, absliq1, absice0, absice1, &
      !$acc&      absice2, absice3) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb01(async_id)
      use module_radlw_kgb01
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb, ka_mn2, kb_mn2) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb02(async_id)
      use module_radlw_kgb02
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb03(async_id)
      use module_radlw_kgb03
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb, ka_mn2o,  kb_mn2o) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb04(async_id)
      use module_radlw_kgb04
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb05(async_id)
      use module_radlw_kgb05
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb, ka_mo3, ccl4) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb06(async_id)
      use module_radlw_kgb06
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, selfref, forref, &
      !$acc&      fracrefa, ka_mco2, cfc11adj, cfc12) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb07(async_id)
      use module_radlw_kgb07
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb, ka_mco2, kb_mco2) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb08(async_id)
      use module_radlw_kgb08
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb, ka_mo3, ka_mco2, kb_mco2, &
      !$acc&      cfc12, ka_mn2o, kb_mn2o, cfc22adj) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb09(async_id)
      use module_radlw_kgb09
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb, ka_mn2o, kb_mn2o) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb10(async_id)
      use module_radlw_kgb10
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb11(async_id)
      use module_radlw_kgb11
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb, ka_mo2, kb_mo2) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb12(async_id)
      use module_radlw_kgb12
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, selfref, forref, &
      !$acc&      fracrefa) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb13(async_id)
      use module_radlw_kgb13
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, selfref, forref, &
      !$acc&      fracrefa, fracrefb, ka_mco2, ka_mco, kb_mo3) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb14(async_id)
      use module_radlw_kgb14
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb15(async_id)
      use module_radlw_kgb15
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, selfref, forref, &
      !$acc&      fracrefa, ka_mn2) async(async_id)

      end subroutine


      subroutine copyin_radlw_kgb16(async_id)
      use module_radlw_kgb16
      implicit none

      integer, intent(in) :: async_id

      !$acc enter data copyin(absa, absb, selfref, forref, &
      !$acc&      fracrefa, fracrefb) async(async_id)

      end subroutine

end module radlw_copyin_gpu
