# 回覆 AMD（Micky Cheng，2026-10-16 來信）草稿

> 寄出前要決定並填入的三處以【】標示：原始碼存取方式、輸入資料傳遞方式、commit。

---

主旨：RE: GEPS/TCoGFS 完整 app build/run 說明（afar 24.3.0）

Hi Micky，您好：

感謝回覆。完整 app 的 build/run 說明已整理好，放在 GEPS repo 的 `qa/amd_app_repro/` 目錄，
入口是 `README.md`（英文，方便內部 team 直接使用），並附上照步驟執行的腳本。依您列的項目說明如下：

1. **原始碼取得方式**
   - 模式：`https://github.com/DengShunChen/GEPS`，branch `feature/multi-compiler`，commit 【commit hash】
   - 相依函式庫：`https://github.com/DengShunChen/GEPS_LIB`，branch `main`，commit `fe66918d557d0997ee0f317e9c4d35e9a49cebc6`
   - 兩個 repo 目前是 private。【請提供需要存取權的 GitHub 帳號，我們開 read 權限／或改以 git bundle 檔案提供】
   - 輸入資料（初始場 2022-08-15 00Z 與靜態場，CWA DMS 格式）約 10 GB，壓縮檔附 sha256。【傳遞方式】

2. **相依函式庫版本與編譯方式**：OpenMPI 4.1.6、netCDF-C 4.9.0／netCDF-Fortran 4.6.0（HDF5 1.13.2）、
   OpenBLAS 0.3.27、FFTW 3.3.5，以及 CWA 的 DMS 資料庫函式庫等，原始碼都在 GEPS_LIB 內。
   全部用 afar 重新編譯（Fortran `.mod` 綁編譯器版本），腳本 `02_build_deps.sh`，約 30 分鐘。

3. **編譯指令**：`03_build_app.sh`（CMake 參數、afar 編譯器路徑、HSA 配置快取一併建置），約 20 分鐘。

4. **執行方式與環境變數**：`04_run.sh`（1 小時）、`04_run.sh 24`（24 小時）。README §4 逐項列出我們使用的設定與原因，
   包括 `HSA_ENABLE_SDMA=0`、`LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0`、`LIBOMPTARGET_STACK_SIZE`、
   HSA 配置快取（`hsa_front`，上限 48 GB）、`ROCFFT_RTC_CACHE_PATH`。

5. **測試案例與驗證**：TCo383L72、2 MPI ranks（每張 MI300X 一個 rank，需要兩張 192 GB 的卡）。
   1 小時（6 步）約 3 分鐘，24 小時（144 步）約 10 分鐘。每步輸出的 `surf pres tend rms(GPU)` 是驗證值，
   `check_run.py` 會自動與我們的結果（`reference.json`）比對，並列出每步時間。
   我們的參考值：afar 24.3.0 每步 2.31 s（1 小時第 3–6 步平均）、24 小時 431.8 s；ROCm 7.2.2 為 2.21 s、415.7 s。

此外，README §6 列出先前回報的每個問題在 app 上的觸發方式（只改一個設定），例如：
關閉我們的配置快取看 libomptarget pool 的 OUT_OF_RESOURCES、移除 rocFFT 快取看 JIT 時間、
還原 sflx 修改並把快取調到 96 GB 重現 device malloc 在 VRAM 吃緊時變慢、以及 KFD eviction 的計數器觀察方式。

我們的測試環境：2 × MI300X（SPX/NPS1）、amdgpu/KFD 6.16.13、Ubuntu 24.04.5（kernel 6.8.0-138）、
therock-afar-24.3.0-multiarch-10.1.0-592954c。

若重現過程中遇到任何問題，歡迎直接聯絡，我們可以配合遠端排查或提供 log。

謝謝！
登舜
