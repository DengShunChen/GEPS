# GEPS on AMD ROCm — 接手筆記 (Handoff)

> 最後更新：2026-09-15（逐波數 1:1 對照定位到水平半拉格朗日平流；見 §6 第 9 層末尾）
> 目標：**讓 GEPS / TCoGFS 在 AMD MI300X (gfx942) 上跑起來**

---

## 0. 一分鐘摘要

| 項目 | 狀態 |
|---|---|
| ROCm CPU build (`build_rocm_cpu`, amdflang host-only OpenMP) | ✅ 可編、**可完整跑完 1 小時預報**（2026-08-27 01:52 的 log 有 `PROGRAM CWBGFS HAS ENDED`） |
| ROCm GPU build (`build_rocm_hip`, OpenMP target offload + HIP) | ✅ 可編；✅ **2026-09-10 首次跑完全部 3 輪 NNMI 並進入預報積分**（換成 2 張實體 MI300X / SPX 之後）；❌ 但**數值是錯的**：`iteration=2` 起 VRAM 耗盡，垂直平流輸入 `dd` 全為零，積分第一步就 `surf pres tend rms(GPU) = NaN`，最後仍以 `HSA_STATUS_ERROR_OUT_OF_RESOURCES` 中止 |
| **2026-09-18 進版控** | ✅ **跑通的版本已 commit 並 push 到 GitHub**：`feature/multi-compiler` `cefc44d0..4424cf83`（https://github.com/DengShunChen/GEPS/commit/4424cf83），55 個檔案、+11500/−108，工作樹乾淨、與 `origin` 同步。內容＝2026-09-17 18:55 驗證通過的狀態（acc2omp 修法、四類 stream 競態修法、`src/rocm/` 新 kernel、job 腳本、ckpt 比對工具、`qa/rocm_repro/`、本文件與兩份 ticket）。**探針還在 commit 裡**（commit message 有註明），拆探針要另開 commit。`.gitignore` 新增 `job/*.log.*` 與 `/compile_*.log`，約 70 個日期命名的 run log（~70 MB）留在本地。這台機器 push 的兩個坑見 §8「版控」。 |
| **2026-09-17 18:55 最新狀態** | ✅ **ROCm GPU 六步 `surf pres tend rms` 與 CPU 全部對到 1e-11～1e-14**（`job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_val`，tau=1 跑完）；深對流觸發數 29,009 = CPU；物理結束點 `qt` 差 8e-11。最後一片是 `samfdeepcnv_kh_gpu` 的 `val1/val2`（kernel 外初始化卻列 `private` ⇒ ROCm 讀垃圾 ⇒ `heso` NaN ⇒ 0 觸發），acc2omp 的宣告解析補上「無 `::` 的續行宣告」後自動改成 firstprivate。**移植數值上完成。速度（2026-09-18）：探針已拆（無影響）；每步 740 → 500（FFT plan set 8）→ 247（垂直平流 otile 8192）→ 58（D2D 複製改走 omp_target_memcpy；HIP 不認得 libomptarget 的指標，原本走 CPU large-BAR 26 MB/s）→ 28（`kernels` 區域序列化的 def_cfl/hh 改平行、pack/unpack tile-local）→ 11.8 s（水平平流 PPM 的 tile 8→64 列）→ **4.9 s（2026-09-19 00:41，`…_2gpu_hsapool`：`rocprofv3 --hsa-amd-trace` 量到每步 4.1 s 在 HSA alloc/free（threshold 0 讓每個 `data create` 都是裸 HSA 配置，free 不分大小 0.4–5 ms），自製 `LD_PRELOAD` 配置快取 `src/rocm/geps_hsa_pool.cc` 解掉；六步 `sptend` 逐位元不變；VRAM 158/206 GB，快取上限 48 GB/rank）→ 4.63 s（01:28，`…_2gpu_advx`：massadvx 三個序列區域展開＋wait_all 只等 dirty stream）→ 4.33 s（03:45，`…_2gpu_wmax2`：hdiffu/rayleifr wind-max 兩段 kernel、FFT 4 host thread）→ 4.06 s（04:20，`…_2gpu_dgemmq`：Legendre dgemm 迴圈單一 stream 排隊、迴圈後一次 wait）→ 3.49 s（05:05，`…_2gpu_spmd`：物理 96 個 kernel 的經度迴圈折進 collapse(3)、mpe2d 兩個單 thread 抄陣列 kernel 展開）→ 2.80 s（07:50，`…_2gpu_norpc`：把 device code 殘留的 `write/stop` 靜音，映像不再含 `__llvm_rpc_client`，libomptarget 不開 RPC server，kernel 完成通知不再晚 ~一個 kernel 長度）→ 2.51 s（14:25，`…_2gpu_tile2`：NDSL tile 64→128／8192→32768）→ 2.43 s（15:30，`…_2gpu_spmd2`：SPMD pass 擴到 diabat/gwd/rrtmg 的單層 owner）**；CPU 32 核 10.5 s ⇒ 已快過 CPU 4.3 倍；之後是 launch 延遲與 kernel 效率（§9 8f）；待辦：1-GPU、上游回報。** |
| **2026-09-17 05:00** | ✅ **ROCm GPU 首次完整跑完 1 小時預報**（2 張 MI300X，`job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_l18g`，`PROGRAM CWBGFS HAS ENDED`）。step-1 `sptend` 與 CPU 15 位相同；step 2–5 比 CPU 低 4–6%：逐常式對照（08:40）PBL 與 CPU 逐位元等價，**差異在深對流 `samfdeepcnv_kh_gpu`：ROCm 上 0 個 column 觸發（CPU 29,009）**——根因（17:00 量到）：`heso` 整條 NaN，因為 `val1/val2` 在 kernel 外初始化卻被列為 `private`（與 `snoexp` 同病），acc2omp 的 pass 漏掉「無 `::` 的續行宣告」；已修、重編驗證中（`…_2gpu_val`）。第 18 層根因：**runtime 大小的 private 陣列在裝置堆疊上互相覆寫**（`LIBOMPTARGET_STACK_SIZE=65536` 修，`job/regression.ksh`）＋ `adjptqintp_gpu` 的 `plold(i,0,jj)` 越界讀（acc2omp guard）。帶探針每步 ~750 s，拆探針後再量速度。 |
| **2026-09-16 21:17（第 17 層；第 18 層見 §6 末）** | ✅ **積分第 1 步在 2 張 MI300X 上與 CPU 對上**：`surf pres tend rms(GPU) = 0.43807396825059497`（CPU 0.4380739682505907）、NNMI `bal` 三輪 1.7619218910e-05 / 6.3945974578e-06 / 1.6451316149e-06、**105 個 step-1 檢查點全部相對差 < 1e-9**（`s1_hdiffu_vormid` 2e-15、`s2_pcorr` 1e-8）。**但第 1 步末的微物理 `fall_flux_gpu` 仍卡死**（rocgdb 實測：`__omp_offloading_…fall_flux_gpu_l625`，即 sedimentation `DO while (notlast)` kernel，916 個 wave 全在裡面轉；自 `mm0` 起每一次 step-1 正確的 run 都停在這裡，與 pool/競態無關，是獨立的 bug）。根因：**第 4 類 stream 競態——裸 `NCCLCHECK(ncclAllReduce…)` 之後沒有 wait**（`mpe2d_gpu` 12 處、`hdiffu_gpu`、`rayleifr_gpu`、`intgrt_gpu`），修在 `cmake/acc2omp.py`（§6 第 17 層）。待辦：跑完 tau=1、拆探針做乾淨對照。log：`job/TCo383L72_IC_sample_rocm.log.20260916_2gpu_l17`、CPU 基準 `job/TCo383L72_IC_sample_rocm_cpu.log.20260916_hd`。 |
| GPU 端是否已完成 tau=1 積分 | ❌ **還沒**，但 **NNMI 已於 2026-09-15 修好**：三輪 `bal`/`pt`/`qgini` 與 CPU/NVIDIA 對上 8-9 位（第 10 層 memset 競態 ＋ 第 7 層延伸到 6 個 LT 常式，都在 `cmake/acc2omp.py`）。積分第 1 步 `surf pres tend rms` = 132.6（應為 0.42，不再是 NaN），第 2 步死於 VRAM 耗盡（舊問題）。**2026-09-15 定位過程**：`tranrs`/FFT/LT/轉置整條鏈已證明**逐波數與 CPU 相同**（384 個 `mf` 全部 1.000）、`hldten` 一直是對的、垂直平流也清白；**壞的是水平半拉格朗日平流 `ndslfv_monoadvh2_gpu_refactor` 的輸出**（`ddtemp` 進垂直平流前 Σ² 只有 CPU 的 0.48 倍，`vdzonl` 0.54、`vdmerd` 1.38），其 kernel 正是 `src/rocm/cyclic_cell_*_gpu.f90` 的手寫 tiled 改寫。詳見 §6 第 9 層末尾「逐波數 1:1 對照」以下各節。以下為 2026-09-14 以前的舊摘要：崩潰類問題**都修好了**（第 1-3、6、7 層，`memory access fault` = 0），GPU 已能跑完三輪 NNMI 並進入積分。剩下的是**數值錯誤**：`bal` 比 CPU 大 8 個數量級（2624 vs 1.76e-05）、跑批間非決定性、`pt` 不變、積分 `NaN`。**已排除**：hipSOLVER、`use_device_addr` 指標傳遞、特徵向量正確性、以及第 9 層的「內層 private 被丟掉造成競態」（修法已 landed 並通過編譯、911 個指令受惠，但 `bal` 完全沒變）。**目前最強線索**：GPU 自己的數字矛盾——`bal` 推得 `\|wrk\| ~ 1`，`x_out` 推得 `\|wrk\| ~ 1e-3`，相差 1000 倍，指向 dgemm 寫的緩衝區與 OpenMP kernel 讀的不是同一塊（第 1、2 層那類 bug） |
| 最後一次 CPU 重跑（2026-08-28 06:14, NPEY=32） | ❌→✅ 已修好，見下 |
| CPU 回歸（2026-09-08 00:30, NPEY=32） | ✅ **已修復並重跑成功**，`tau=1.000` + `PROGRAM CWBGFS HAS ENDED`，`surf pres tend rms` 序列與 08-27 完全一致（見 §6） |
| 工作樹（2026-09-18 起） | ✅ 乾淨，全部在 `4424cf83`。以下為 2026-09-10 的舊描述：GPU **修法**：2 個檔案、112 行淨改動（`src/nvidia/mpe2d_gpu.f90` 77 行 + `src/nvidia/zx_gpu.f90` 35 行），全部 `#ifdef USE_HIP` guard、NVIDIA 路徑完全不受影響。另有 2 個檔案是 2026-09-10 新增的**除錯探針**（`cmake/acc2omp.py`、`src/rocm/vertical_cell_advect_gpu.f90`），**會改變執行時行為**（含 `hipDeviceSynchronize()`），詳見 §6「第 5 層根因」末。**全部尚未 commit**，branch = `feature/multi-compiler` |

**下一步最短路徑**：CPU 已綠燈。GPU 端第 1-3 層是真原始碼 bug（已修好）；**第 6 層**（RCCL 集合通訊後的競態）已修好並驗證，但乾淨對照顯示它對記憶體與 `NaN` 無影響；**2026-09-11 查到第 7 層，就是 `NaN` 的根因**：`nnmi_gpu` / `tendget_gpu` / `zx_gpu` 把計算包在 HIP graph capture 裡，而 acc2omp 丟掉 `async` 之後，capture 內的 OpenMP kernel 既在錄製期搶先執行、又沒被錄進 graph——**GPU 的 NNMI 因此是個完全的 no-op**（`bal` 精確為 0、`pt` 逐位元不變），模式帶著不平衡場進積分才炸成 `NaN`。**第 6 層與第 7 層是同一個病根（丟掉 `async`）的兩種表現。** 下一步是把這 3 個檔案的 ROCm 路徑改成 eager 執行（**注意 `cg_created` 陷阱，見 §6 第 7 層**），驗證只要跑到 `iteration= 1` 之後（約 25 分鐘）看 `bal` 與 `pt` 即可。剩下的「記憶體不歸還」是另一個獨立問題。

---

## 1. 這個專案是什麼

**GEPS / TCoGFS** 是中央氣象署（CWA，舊稱 CWB）的**全球數值天氣預報模式**：

- 主程式 `tcogfs.x`，Fortran 為主（`src/` 下 376 個檔），MPI + OpenMP，2D 分解（NPEX × NPEY）。
- 譜方法 + TCo（cubed-octahedral / triangular-cubic-octahedral）格點，解析度以 `JCAP` 指定（測試用 **TCo383L72** = JCAP 383、72 層）。
- 物理過程齊全：RRTMG 輻射、GCE 3ICE / Thompson 微物理、YSU 邊界層、對流、重力波拖曳、陸面、海冰、隨機物理（ensemble 用）。
- 資料 I/O 走 **DMS**（CWA 自有的 key-value 氣象資料庫，`rdmsdbcrt` / `rdmscrt` / `rdmspurge` 工具 + `~/.dmsrc` 指路徑）。輸出可選 DMS 或 GRIB2。
- 原生支援三個平台世代：Fujitsu FX1000（ARM，`build.sh`）、x86_64 + **NVIDIA A100（OpenACC + CUDA Fortran）**、以及本次新增的 **AMD ROCm**。
- 相依函式庫全部由**兄弟 repo `../GEPS_LIB`** 自己 build（每個 compiler 一棵 install tree）。

**倉庫位置**：`/mlsteam/workspace/data/geps/GEPS`（git 上游 `git.rdc.cwb/tco/tco639l72`，branch `feature/multi-compiler`）

---

## 2. 目錄地圖

```
/mlsteam/workspace/data/geps/
├── GEPS/                       ← 模式本體（本文件所在）
│   ├── CMakeLists.txt          ← 頂層；USE_CUDA / USE_HIP / USE_ACC / USE_OMP 開關
│   ├── compile                 ← 多編譯器 driver：./compile {nvidia|gnu|intel|fujitsu|rocm|rocm-cpu}
│   ├── cmake/
│   │   ├── compilers/ROCM.cmake ← AMD 工具鏈定義（flags / hipBLAS / RCCL / offload-arch）
│   │   └── acc2omp.py          ← ★ 核心：OpenACC + CUDA Fortran → OpenMP target 的原始碼轉譯器（~1600 行）
│   ├── src/                    ← CPU 主體
│   │   ├── nvidia/  (71 檔)    ← NVIDIA GPU kernel（OpenACC/CUDA Fortran）＝ AMD 的翻譯輸入
│   │   ├── cwa/                ← CWA 物理 GPU 版
│   │   └── rocm/               ← ★ AMD 專用 shim / 手寫替代 kernel（見 §4）
│   ├── job/
│   │   ├── regression.ksh      ← 共用跑批腳本（已加 rocm 分支）
│   │   ├── TCo383L72_IC_sample_rocm      ← GPU smoke（1 MPI，指向 build_rocm_hip）
│   │   ├── TCo383L72_IC_sample_rocm_cpu  ← CPU smoke（32 MPI，指向 build_rocm_cpu）
│   │   └── rocm_gpu_bind.sh    ← 只綁 PCI 0000:c5:00.0
│   ├── build_rocm/             ← ⚠️ 舊的（Aug 24），`./compile rocm` 會砍掉重建
│   ├── build_rocm_cpu/         ← CPU-only 樹（USE_HIP=OFF）
│   ├── build_rocm_hip/         ← ★ 真正的 GPU 樹（USE_HIP=ON），**手動 cmake 建的，別動**
│   ├── work_rocm/  work_rocm_cpu/  ← 跑批 workdir（regression.ksh 產生）
│   ├── fix/                    ← 靜態場 / 氣候場（orography、O3、aerosol…）
│   └── modulefiles/modulefile.tcogfs.rocm
├── GEPS_LIB/                   ← 相依函式庫；install/{gnu,rocm}/ + install/modulefiles/
└── dmsdb/                      ← 本機 DMS 資料庫（IC / 背景場 / 輸出）
    ├── ncep_ana_n1.ufs/TCo383l72_22081500/   ← 初始場（2022-08-15 00Z + 前 6h FG）
    ├── bckdms.ufs/BCK_TCo383_GI30S_xnew/     ← 背景/靜態場
    └── TCo383L72.ufs/                        ← 執行期 MASOPS / FCSTDMS
```

---

## 3. 建置路徑（Build matrix）

```sh
export GEPS_LIB_ROOT=/mlsteam/workspace/data/geps/GEPS_LIB
./compile rocm       # → build_rocm      (USE_HIP=ON)   ⚠️ 會 rm -rf build_rocm
./compile rocm-cpu   # → build_rocm_cpu  (USE_HIP=OFF)
```

| Compiler | 編譯器 | GPU | 數學庫 | 輸出目錄 |
|---|---|---|---|---|
| `nvidia` | NVHPC + OpenACC/CUDA | A100 | NVHPC BLAS + cuBLAS/cuFFT/NCCL | `build_nvidia` |
| `gnu` / `intel` / `fujitsu` | CPU only | – | OpenBLAS / MKL / SSL2 | `build_<name>` |
| **`rocm`** | amdclang / amdclang++ / **amdflang** | **MI300X gfx942** | OpenBLAS + **hipBLAS / hipFFT / hipSOLVER / hipSPARSE / RCCL** | `build_rocm` |
| **`rocm-cpu`** | 同上，host OpenMP | – | OpenBLAS | `build_rocm_cpu` |

> ⚠️ **重要落差**：`job/TCo383L72_IC_sample_rocm` 和 `regression.ksh` 都指向 **`build_rocm_hip`**，但 `./compile rocm` 產生的是 **`build_rocm`**。目前的 HIP 樹是**手動 `cmake -Bbuild_rocm_hip`** 建的，且 `CMAKE_Fortran_COMPILER=/opt/rocm/bin/amdflang`（不是 `mpifort`，MPI 靠 ROCM.cmake 手動塞 include dir + FindMPI）。**要嘛把 `compile` 的輸出目錄改成 `build_rocm_hip`，要嘛把 job 腳本改回 `build_rocm`** — 這是待辦第一項。

---

## 4. AMD 移植策略（核心設計，務必先讀懂這段）

原則：**不改 NVIDIA 原始碼**。`amdflang` 不能把 OpenACC offload 到 AMDGPU，所以：

### (a) `cmake/acc2omp.py` — 編譯期原始碼轉譯
建置時把 `src/nvidia/*.f90`、`src/cwa/*.f90` 等 GPU 檔複製到 `build_*/src/hip_src/` 並轉譯：
- `!$acc parallel loop` → `!$omp target teams distribute parallel do`
- `!$acc enter/exit data`、`update`、`host_data use_device` → 對應 `!$omp target ...`
- CUDA Fortran 擴充（`device`/`pinned` attribute、`<<<>>>`…）消毒
- 額外處理：`#ifdef USE_CUDA` 解析、重複宣告刪除、`intent(in)` 放寬、device 端 I/O 消音、`parameter` 不進 map clause 等
- `apply_file_fixups()` 針對特定檔案做點名修補，甚至用 `src/rocm/*.f90` **整支替換** subroutine

### (b) `src/rocm/` — 手寫的 AMD 側程式碼
| 檔案 | 作用 |
|---|---|
| `cudafor.f90` `openacc.f90` `cublas.f90` `cufft.f90` `cusolverDn.f90` `cusparse.f90` `nccl.f90` `nvtx.f90` | 假裝成 NVIDIA 的 Fortran module，API 轉呼叫 HIP |
| `hip_compat.cc` (666 行) | ★ C++ 實作層：hipStream/hipBLAS/hipFFT/hipSolver/hipSparse/RCCL handle 管理、`omp_get_mapped_ptr` 做 host→device 指標轉換 |
| `flang_compat.f90` `flang_device_rt.c` | 補 amdflang 缺的 runtime |
| `eigen_mx_gpu.f90` `ndslfv_para_gpu.f90` `cyclic_cell_*_gpu.f90` `vertical_cell_advect_gpu.f90` `ndslfv_monoadvv_gpu.f90` `mpe2d_row_broadcast_gpu.f90` | 「方言太重、翻譯器搞不定」的 kernel，**直接手寫 OpenMP target 版**取代 NVIDIA 檔（在 `src/CMakeLists.txt` 用 `REMOVE_ITEM` 換掉） |
| `stack_probe.c` | debug 用（stack 位址探針） |

### (c) 編譯細節（踩過的雷都寫在註解裡）
- **只有 GPU 檔加 `--offload-arch=gfx942`**。全域加會讓 device LTO 吃進上百個空的 amdgcn object，link 要好幾小時。
- `-Wl,--allow-multiple-definition`：`sbyte_` 同時在 `libw3.a` 和 `libg2.a`，lld 會擋。
- amdflang 不吃 gfortran 的 `-ffree-line-length-none` / `-fallow-argument-mismatch` / `-fno-range-check`。
- `mpi.mod` 在 OpenMPI 的 `lib/`（不是 `include/`），ROCM.cmake 手動加 include dir。

---

## 5. 環境事實（這台機器）

- **ROCm 7.2.2**（`/opt/rocm` → `/opt/rocm-7.2.2`），`amdflang` / `amdclang` / `hipcc` / `rocm-smi` 都在 `/usr/bin`
- **GPU：2 張實體 MI300X（2026-09-10 起，由系統管理者調整）**，兩張都是 **NPS1 / SPX**（未分割）：`rocm-smi` Device 0 = KFD node 4、Device 1 = KFD node 5，**是兩張真的卡，不再是同一張卡的 DPX 分區**。每張 304 CU（`simd_count=1216`, `num_xcc=8`）、VRAM 206,141,652,992 B（~192 GiB）。**這 2 張卡是專屬的，不與其他使用者共用**，所以跑批不必排隊、可以放心跑滿（單次 smoke 約 1 小時、峰值 135 GB VRAM）。目前 job 仍只綁 GPU0（`HIP_VISIBLE_DEVICES=0` / `ROCR_VISIBLE_DEVICES=0`），**GPU1 全程閒置**；多卡要等 RCCL 擴展（§9）。
  - **⚠️ 舊記錄的訂正**：2026-09-09 以前這台是 1 張 DPX 分區卡，當時 §6 第 4 層推測「DPX 壓縮了固定資源配額」。**換成 SPX 後實測，KFD 的 per-node 配額數字一個都沒變**（`num_sdma_engines=2`、`num_sdma_queues_per_engine=8`、`num_gws=64`、VRAM total 相同），所以那個推測在這些**可見屬性**上並不成立；真正的差異出在 **runtime 實際配置得到的量**（DPX 全程卡在 ~46 GB，SPX 可以長到 135 GB），機制尚未查明。
- **CPU：384 cores**
- MPI：`GEPS_LIB/install/rocm/openmpi-4.1.6`（用 amdflang 自己編的，`--enable-mpi-fortran=all`，`use mpi` 走 auto-detect `!DIR$ IGNORE_TKR`）。**⚠️ 它的 `MPI_DOUBLE_PRECISION` 是 16 bytes**（configure 時只給了 `-fdefault-real-8`，沒有 `-fdefault-double-8`），模式的 `double precision` 是 8 bytes ⇒ 用到它的呼叫都搬 2 倍位元組。目前 NPEX=1 的路徑上沒踩到，但 NPEX>1 / IO server / `sigful` 會。探針一律用 `MPI_REAL8`。詳見 §6 第 9 層末尾「附帶發現 1」。`env.sh` 會先找原生 ROCm MPI（`GEPS_ROCM_MPI_ROOT` / `ROCM_MPI_ROOT` / `MPI_HOME`），沒有才用自編的
- DMS 工具：`GEPS_LIB/install/rocm/dms-4.0.0/bin`，資料庫 `GEPS_DMSDB=/mlsteam/workspace/data/geps/dmsdb`
- 跑批以 root 執行 → `OMPI_ALLOW_RUN_AS_ROOT[_CONFIRM]=1`
- **⚠️ 容器 rootfs 不持久**：只有 `/mlsteam/workspace`（NFS）跨 session 保留；`/`、`/usr`、`/lib` 等系統套件每次容器重建都會被重置。曾經手動裝過、但不在底層 image 裡的套件（目前已知：`libtirpc3t64` `libhwloc15` `libevent-core-2.1-7t64` `libevent-pthreads-2.1-7t64`，供 rocm-cpu/rocm-hip 的 `tcogfs.x` 動態連結用；另外 **`libdw1t64` 是 `rocprofv3` 需要的**，缺了會以 `libdw.so.1: cannot open shared object file` 失敗）**每次新 session 都可能要重裝**。裝法：`apt-get install -y libtirpc3t64 libhwloc15 libevent-core-2.1-7t64 libevent-pthreads-2.1-7t64 libdw1t64 libtirpc-dev libcurl4-openssl-dev && ldconfig`。
  - **⚠️ 2026-09-10 新發現：還有兩個「編譯期」相依 `libtirpc-dev` / `libcurl4-openssl-dev`**。執行期的 `libtirpc3t64` 只給 `libtirpc.so.3`，**link 需要的是 dev 套件的 `libtirpc.so` symlink**。缺了會在 link 階段以 `clang-linker-wrapper: error: unable to find library -ltirpc`（接著 `-lcurl`）失敗——也就是說**容器即使能跑也可能不能重編**。更麻煩的是 **cmake 在 link 前會先刪掉既有的 `bin/tcogfs.x`，所以一次失敗的 link 會讓你暫時沒有可用的執行檔**，必須裝好套件重 link 才會回來。要確認能不能編，最快是看 `build_rocm_hip/src/CMakeFiles/tcogfs.x.dir/link.txt` 裡的每個 `-l` 是否都找得到對應的 `.so`/`.a`。。判斷是否需要重裝：`ldd build_rocm_cpu/bin/tcogfs.x`（先 module load）看有沒有 `not found`。
- Smoke 設定（`SMOKE=1`）：`taue=1.0, tauo=1.0, taup=6.0`（taup 必須留 6，getrdy 才讀得到 22081418 的 6h FG）；`outgrb2=0, outdms=1`

---

## 6. 已知問題 / 踩過的雷

### 已繞過（workaround 已在工作樹裡）
1. **amdflang `-O3` spill bug** — 在超大 subroutine 裡，spill 到 stack 的位址會被整數運算蓋掉，造成 `r8=NULL` 的 SIGSEGV（`intgrt` 的 `jtrun`→`hdiffu`、`initial` 的 `bal_tmp`）。
   → `src/CMakeLists.txt` 對 **`intgrt.f90` / `initial.f90` / `eigenlib.f90` / `mpe_unify.f90` 個別加 `-O0`**（只在 `NOT USE_HIP` 時）。
   → `initial.f90` 另外把 `jtrun/jtmax/lev/nx/my/my_max/no` 複製成 local `*_l`、`bal_tmp` 改成 `allocatable`。
2. **GRIB2 輸出 SIGSEGV** — `wrt_grb2`/`addgrid` 在 tau=0 掛掉（`imax=jmax=-1`）。→ smoke 一律 `outgrb2=0, outdms=1`。**這個 bug 還沒查。**
3. **1-rank GPU 記憶體不足** — `src/nvidia/mpe2d_gpu.f90` 的 `transpose_ndsl_p2f/f2p`，單 rank 時 `c1`/`c2` 各數百 MB，TCo383 在 NDSL 之後只剩 ~26 GB。→ 加了 `nsizex==1` / `proc==1` 的直通捷徑（**這是有改 NVIDIA 原始碼，未來要想辦法搬走**）。
4. `regression.ksh` 原本寫死 CWA 內網路徑（`/data/common/gfs/...`、`/package/x86_64/dms/...`、`/users/xb80/bin/Caldtg.ksh`）→ 全部加了 `rocm_lab=1` 分支。

### 已解決（2026-09-08）
- **2026-08-28 06:14 的 CPU 回歸失敗**（`fct model fail`）**根因不是模式或 job 腳本，是容器本身**：這個工作環境的 **容器 rootfs 不持久**（`/mlsteam/workspace` 是 NFS mount 會留著，但 `/`、`/usr`、`/lib` 這些系統套件會在容器重建時被重置回底層 image）。08-27 01:52 成功那次的容器裝過 `libtirpc3t64` / `libhwloc15` / `libevent-core-2.1-7t64` / `libevent-pthreads-2.1-7t64`（Ubuntu 24.04 套件），08-28 用的是**重建過的新容器**，這些套件不在裡面，`ldd` 才會對 `libtirpc.so.3` / `libopen-rte.so.40` / `libopen-pal.so.40` / `libhwloc.so.15` / `libevent_*` 回報 `not found`。`libopen-rte`/`libopen-pal` 本身是 GEPS_LIB 自編的（存在於 `install/rocm/openmpi-4.1.6/lib`），只要模組正確 `module load` 就找得到；真正遺失的系統套件是 **libtirpc / libhwloc / libevent**。
  → **修法**：`apt-get install -y libtirpc3t64 libhwloc15 libevent-core-2.1-7t64 libevent-pthreads-2.1-7t64`（套件是 apt cache 裡現成的，不需要外部網路重新整理 index）。裝完 `ldconfig`，`ldd build_rocm_cpu/bin/tcogfs.x`（配合 `module load`）就乾淨了。
  → **重跑驗證**：`job/TCo383L72_IC_sample_rocm_cpu.log.20260908` 完整跑完 `tau=1.000` + `PROGRAM CWBGFS HAS ENDED`，數值與 08-27 完全一致。
  → **⚠️ 這代表每次容器重建（新 session / 新 pod）都可能要重裝這 4 個套件**，才能跑 rocm-cpu / rocm-hip 的 binary。之後可以考慮寫一支 `job/ensure_system_deps.sh` 或在 job 腳本 preflight 加一段自動偵測 + 提示。
### GPU 端目前不會算完 tau=1 — 深入除錯記錄（2026-09-08 ~ 09-10，逐層剝洋蔥）

之前（§0/舊版）以為 08-27 那次 GPU 跑批「跑得起來、寫出 FCSTDMS 輸出」代表有算完。**重新看 log 才發現：那 509 筆 DMSRPUT key 全部是 `outflds_green for tau=0` 階段寫出的初始場回顯（echo of IC），不是預報結果**；真正的預報積分根本還沒開始就會崩潰。

這個 bug 的除錯過程是「剝洋蔥」式的：每修好一層，程式就跑得更遠，然後撞到**下一個、完全不同**的崩潰。已經剝掉 3 層，**iteration=1（第一次 NNMI 迭代，含 tendget_gpu + zx_gpu + vartrix_gpu + nnmi_gpu 全套）現在可以完整跑完**，卡在 iteration=2 開頭的第 4 層。以下完整記錄每一層，供下一次接手者不要重複繞路。

#### 第 1 層（已修好）：`mpe2d_unify_spec_lev_zx_gpu` 對 reshape 過的 zx_buf 做 memset

**症狀**：`iteration=1` 後立刻崩潰，`hsa_amd_memory_lock: HSA_STATUS_ERROR`，size 固定 510935040 bytes（= `lev*6*jtrun*jtmax*8`，即整個 `zx_buf`）。

**定位過程**：`#ifdef USE_CUDA`（ROCm 建置也定義這個巨集，見 §4a）下，`src/gfcst.f90:108` 呼叫的是 `initial_gpu.f90`（`src/nvidia/initial_gpu.f90`），**不是** CPU 版的 `src/initial.f90`。`initial_gpu.f90` 的 `do ic=1,nnmiit` 迴圈裡依序呼叫 `tendget_gpu` → `zx_gpu` → `vartrix_gpu` → `nnmi_gpu`。在 `zx_gpu.f90` 開頭，`CALL mpe2d_unify_spec_lev_zx_gpu(wrk, ...)` 把 `zx_buf`（1D，`initial_gpu.f90` 宣告）以 `wrk`（5D，`zx_gpu.f90` 宣告，`lev,2,3,jtrun,jtmax`）的姿態傳進去，再以 `aout`（5D，`mpe2d_unify_spec_lev_zx_gpu` 宣告）的姿態接住——**同一塊記憶體，三層 subroutine call 各自宣告成不同 shape**（Fortran sequence association，這份程式碼到處都是這招）。`mpe2d_unify_spec_lev_zx_gpu` 一開始用 `!$acc host_data use_device(aout, work)` + `cudaMemsetAsync` 把 `aout` 清零——**這就是崩潰點**：amdflang 的 `!$omp target data use_device_addr` 不認得「一個已存在陣列的不同 shape 重新詮釋」是同一塊記憶體，嘗試做一次全新的 host→device 拷貝，失敗。

**修法**（`src/nvidia/mpe2d_gpu.f90`，`#ifdef USE_HIP` guard，NVIDIA 不受影響）：把 `cudaMemsetAsync` 換成兩個 `!$acc parallel loop`（implicit-map kernel）逐元素清零。Implicit-map kernel 不像 `use_device_addr` 一樣需要「presence 檢查」，親測可以正確處理 reshape 過的陣列。

**已排除、不要重試的假設**（每個都花了一次完整 15 分鐘 build+run 驗證）：
1. ❌ GPU/driver 壞掉——獨立 HIP 小程式（`hipHostRegister`+`hipMalloc`+`hipMemcpy`，同樣 510935040 bytes）直接成功。
2. ❌ host 記憶體不足/ulimit——`free -h` 1.8 TiB 可用，`ulimit -l` unlimited。
3. ❌ 純粹的「陣列切片/reshape 定址」機制本身有問題——寫了 6 支獨立 amdflang repro（單一大陣列切片、累積 map 40 個 ~687MB 陣列到 27.5GB、跨編譯單元、與真實碼完全相同的 buffer size、implicit-map kernel 對 sequence-associated dummy）**全部成功，重現不出來**。純粹的 reshape/切片機制是好的，問題只在 `use_device_addr`/`host_data use_device` 這個特定 construct 上。
4. ❌ CUDA/HIP graph capture mode 太嚴格（`hipStreamCaptureModeGlobal=0`）——改 `hipStreamCaptureModeRelaxed=2` 完全沒用。
5. ❌ CUDA/HIP graph capture 本身不相容——把整個 capture 機制改成 no-op（eager 執行）依然在同一位置崩潰，只是晚 20 秒。
   （#4 #5 的實驗已完全 revert，`git diff` 乾淨。）
- `-g` debug info 加了也沒用，amdflang 目前不會把 source location 寫進 OpenMP target 的 metadata。
- `LIBOMPTARGET_INFO=-1` 全開一次能衝到 84GB（每個 map 事件都印整張表），過濾到只留關鍵字才可控。**若要重複這類追蹤，先用窄 grep pattern 或直接加 `write(0,*)`+`flush(0)` 探針，不要無腦開 `-1`。**

#### 第 2 層（已修好）：`zx_gpu` 讀取 `cudaMallocAsync` 出來的裸 device 指標

**症狀**：修好第 1 層後重跑，換了個位置崩潰——一樣的 `hsa_amd_memory_lock`，一樣的 size，但發生在 `zx_gpu.f90` 自己的最後一個 `!$acc parallel loop`（讀 `vars(idx)` 寫入 `vorten/divten/phiten`），而不是 `mpe2d_unify_spec_lev_zx_gpu`。

**定位過程**：`zx_gpu.f90` 宣告 `REAL(kind=RTYPE), dimension(:), allocatable, device :: vars`——`device` 是 CUDA Fortran 專屬屬性，NVHPC 原生支援，amdflang 不支援。acc2omp.py 把它翻譯成普通 `pointer`（`build_rocm_hip/src/hip_src/nvidia_zx_gpu.f90:38`），`vars` 透過 `cudaMallocAsync` 拿到一塊**純 device 記憶體、完全沒有 host 對應**。NVIDIA 的 OpenACC 編譯器原生認得 `device` 屬性陣列已經在裝置上，不需要任何映射；amdflang 對翻譯後的普通 `pointer` 沒有這個資訊，讀取 `vars` 的 implicit-map kernel 誤以為它是需要拷貝的 host 陣列，拿裝置位址當 host 位址去拷貝，當然失敗。

**修法**（`src/nvidia/zx_gpu.f90`）：把讀 `vars` 的 `!$acc parallel loop` 在 ROCm 路徑換成手寫的 `!$omp target teams distribute parallel do ... has_device_addr(vars)`（NVIDIA 路徑維持原本的 `!$acc parallel loop`）。`has_device_addr` 明確告訴 compiler「這個變數的位址已經是合法的 device 位址，不要映射」。
- 中途試過 `is_device_ptr(vars)`（OpenMP 對這種情境更常見的 clause），但 amdflang 要求該 clause 的引數必須是**字面上的 `type(c_ptr)`**，`vars` 是 Fortran pointer array 不符合，編譯期就報錯——改用 `has_device_addr` 才過。

#### 第 3 層（已修好）：`has_device_addr` 對「graph capture 內用 `cudaMallocAsync` 配置的記憶體」在 replay 時失效

**症狀**：修好第 1、2 層後，`iteration=1` **完整跑完**（包含 `vertical mode l=1/2/3` 的 `bal=` 收斂輸出），但一進 `iteration=2` 就以全新的錯誤崩潰：
```
OFFLOAD ERROR: memory access fault by GPU ... at virtual address 0x...  Reasons: Unknown (0)
```
（不是 `hsa_amd_memory_lock`，是真正的 GPU 記憶體存取錯誤/等同 device-side SIGSEGV。）錯誤前印出的「Last 8 kernels launched」顯示最後一個 kernel 正是 `zx_gpu_ @ 99`——也就是第 2 層剛修好的那個 `has_device_addr(vars)` kernel。

**根因**：`vars` 是在 `if (.not. cg_created)` 區塊裡、`accx_async_begin_capture` **之後**才呼叫 `cudaMallocAsync` 配置的——也就是說配置動作本身也被錄進了 capture 的 graph 裡。HIP/CUDA 的 graph 對這種「capture 內配置、capture 內使用」的記憶體有專門的「graph memory node」機制，會在每次 replay 時重新虛擬化/更新位址；但這個機制是給**原生 CUDA/HIP API 呼叫**設計的，`!$omp target ... has_device_addr` 的 kernel launch 很可能繞過了這個機制，把 capture 當下的（可能只是暫時性的）位址直接烤進 kernel 參數裡，到了下一次 replay（iteration=2）位址已经失效，就變成存取一個無效位址。

**修法**（`src/nvidia/zx_gpu.f90`，ROCm 路徑）：把 `cudaMallocAsync(vars,...)` 搬到 `accx_async_begin_capture` **之前**執行（配上一次 `cudaStreamSynchronize` 確保配置真的完成），讓 `vars` 在 capture 開始之前就拿到一個**真實、穩定、不會變的**位址；capture 只錄「使用」這個位址的動作，不錄「配置」。對應地拿掉 `cudaFreeAsync(vars, stream)`（ROCm 路徑永久不釋放——這是刻意的：因為這個位址被烤進了 graph，只要這個 graph 還會被 replay，這塊記憶體就必須一直有效；程式結束時由 OS 統一回收，不算真的洩漏）。

**驗證**：三層修好後重跑，`iteration=1` 完整跑完（`vertical mode l=1/2/3`、`bal=0.` ×3），成功進入 `iteration=2`，之後才撞到第 4 層（見下）。這是這個 GPU port 有史以來跑得最遠的一次。

#### 第 4 層（✅ **已於 2026-09-10 結案**——先確認是硬體/分區觸發，當天稍晚進一步查明**真正的根因是 GEPS 這邊每輪 NNMI 洩漏 ~150 GB VRAM**，與第 5 層同源，詳見下方「第 5 層根因」一節）

**症狀**：`iteration=2` 剛開始就崩潰，錯誤訊息**跟前 3 層都不一樣**：
```
"PluginInterface" error: Failure to allocate device memory: "unknown or internal error"
error in hsa_amd_memory_pool_allocate: HSA_STATUS_ERROR_OUT_OF_RESOURCES:
The runtime failed to allocate the necessary resources. ...
```
這是**真正的資源耗盡**（配置失敗），不是位址/映射問題。Log 見 `job/TCo383L72_IC_sample_rocm.log.20260908_iter2crash`。

**精確定位（2026-09-09，已用 print 探針確認，非猜測）**：在 `tendget_gpu.f90`/`mpe2d_gpu.f90` 裡對每個 `enter data`/`exit data`/`graph_launch` 呼叫都加了 `write(0,*)`+`flush(0)` 探針。結果：`mpe2d_unify_spec_lev_zx_gpu` 的 `enter data create(work, ain)`/`exit data delete(ain, work)` 在 `iteration=1` 內執行了兩次（`zx_gpu` 前後各一次），**兩次都成功**。真正失敗的是 **`tendget_gpu.f90` 自己最上面那段 `!$acc enter data create(sdpbl, deldm, ww1, ttm_sl, pten_sl, uum_sl, vvm_sl, qm_sl, pdot, vdmerdr, vdzonlr, vdmerd, vdzonl, ddtemp, qvadv, diveng, pten)`（17 個純 scratch 陣列）——探針顯示 `T1 before` 印出、`T1 after` 沒印出，崩潰精確發生在這個 enter_data 的**第二次呼叫**（`iteration=2` 進入 `tendget_gpu` 的第一件事）。

**已排除、不要重試的兩個修法**：
1. ❌ **把這段 enter_data 改成只做一次（比照 `vars` 的持久化修法，用 `cg_created` 保護）**——這是錯的：`sdpbl`/`deldm`/... 是**普通 Fortran automatic（stack）區域變數**，不是像 `vars` 那樣的持久 device 指標，它們的 host 位址**不保證**每次呼叫都一樣。跳過 enter_data 導致該位址範圍沒被正確追蹤，結果在 `iteration=1`（不是 2）就整個提早炸掉，錯誤變成完全不同的 `explicit extension not allowed: host address ... (510935040 bytes) ... device allocation maps to host at ... (687449088 bytes)`——這是 `zx_buf` 的 stack 位址跟某個 `mod_grid.f90` 陣列的位址發生碰撞，證實「跳過 enter/exit」會讓 stack 區域被其他變數重用，破壞正確性。**已完整 revert。**
2. ❌ **在 exit_data 後面加 `!$acc wait(async_id)`（強制同步，確保 async 的 free 真的做完才返回）**——完全沒用，跟原本一模一樣的位置、一模一樣的錯誤。**已完整 revert。**

**已排除的環境變數實驗（cheap，不需要 rebuild，直接改 env 重跑）**：
3. ❌ `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0`（完全關閉 libomptarget 內建的 free-list 記憶體池，所有配置走原始 `hipMalloc`/`hipFree`）——**讓情況變得更糟**：崩潰提早發生在 `iteration=1` 內部（`tranrs_gpu_cuda_graph_`/`mpe_transpose_rs_sp_gpu_`/`gridnl_hybrid_ndsl_gpu_refactor_` 一帶，完全不同的程式碼區域），而且是真正的 `OFFLOAD ERROR: memory access fault`（不是資源耗盡）。這證明**記憶體池預設是有在幫忙的**——很可能靠著 free-list 位址重用，意外掩蓋了程式碼裡其他潛在的、跟前三層同類型的「reshape/reuse 陣列」bug。
4. ❌ `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=4294967296`（4GB，把門檻拉高讓更多配置走池子而不是直接 hipMalloc）——**結果跟預設值一模一樣**，同一個錯誤、同一個位置。代表預設門檻早就已經大到能覆蓋 tendget_gpu 這些陣列了，調高沒有意義；問題不在「池子門檻」，而在池子/底層配置機制本身處理不了這種重複配置模式。

**用 `rocm-smi` 直接量測 VRAM（決定性證據，非猜測）**：在整個崩潰過程中每 0.5 秒輪詢一次 `rocm-smi --showmeminfo vram`，總量 206,141,652,992 bytes（~192 GiB），使用量穩定停在 ~49,114,161,152 bytes（~46 GiB），**從頭到尾幾乎不變**，還有 **~150GB 完全沒用到**。**這徹底排除「真的沒記憶體了」**——`HSA_STATUS_ERROR_OUT_OF_RESOURCES` 說的「resources」根本不是 VRAM 位元組，是某種**數量型**的資源（queue、signal、或類似的固定配額）。

**這台機器的 KFD topology 線索**：`/sys/class/kfd/kfd/topology/nodes/4/properties`（我們這張 MI300X）顯示 `num_sdma_engines=2`、`num_sdma_queues_per_engine=8`（合計固定 16 條 SDMA queue）、`num_gws=64`。這些都是**很小、固定**的配額。如果 libomptarget 每次 `enter data`/異步資料搬移都跟 SDMA queue 池要一條、且沒有正確歸還（queue 洩漏），累積一整輪 `iteration=1`（tendget_gpu ×1 + mpe2d_unify_spec_lev_zx_gpu ×2 + 其他好幾個 subroutine 各自的 async 資料搬移）下來，很可能就把這 16 條配額用光，導致 iteration=2 一開始下一個配置就失敗——**這是目前最有根據的假設，但還沒有像前 3 層那樣用獨立 repro 或更細的探針證實到「就是這個」的程度**。

**已排除（2026-09-09）：把 tendget_gpu 這 17 個陣列合併成更少次的 enter/exit data**——這是 §9 建議的「風險中等、最有機會透過原始碼修正」的方向，**已經動手做並測試過，結果是負面但很有資訊量的**：
- 逐一verify了這 17 個陣列會被傳進去的**每一個** downstream subroutine（`mpe2d_transpose_ndsl_p2f_gpu`/`f2p_gpu`、`ndslfv_monoadvh2_gpu_refactor`、`joinrs_gpu`/`join1rs_gpu`、`ndslfv_update_gpu`、`ndslfv_monoadvv_gpu`、`mpe2d_unify_nx_gpu`、`tranrs1_gpu`、`trandv_gpu_cuda_graph`），確認除了 `mpe2d_unify_nx_gpu`（直接對 `work`/`a` dummy 做 `host_data use_device`，且沒有 1-rank early-return guard）以外，其餘全部安全（不是有 1-rank/single-rank 的 early-return 提早跳過 host_data，就是 host_data 只碰它們自己的區域暫存變數）。
- 因此把 15 個陣列（排除 `deldm`、`ww1`，因為它們會進 `mpe2d_unify_nx_gpu`）合併成**一個**用 F2008 pointer bounds-remapping（`p(l:u) => buf(a:b)`）指到的暫存 buffer，只做**一次** `enter data create`/`exit data delete`（外加 `deldm`、`ww1` 各自維持原樣），從 17 個追蹤實體降到 3 個。amdflang **可以編過**這種 pointer remapping 語法。
- **結果：完全沒用，一模一樣的崩潰位置、一模一樣的錯誤**（`iteration=2` 一開始，`HSA_STATUS_ERROR_OUT_OF_RESOURCES`）。已完整 revert。
- **這個負面結果本身很關鍵**：徹底推翻了「resource 消耗量跟這個 enter_data 子句裡列出的陣列『數量』成正比」的假設——如果真的跟數量成正比，17→3 應該要有感的差異（撐更久、或錯誤位置移動），但完全沒有。**這代表問題根本不是 tendget_gpu 自己 17 個陣列的事**：更可能是 `iteration=1` 整輪下來，橫跨**所有**呼叫到的 subroutine（`mpe2d_transpose_ndsl_p2f/f2p_gpu` ×多次、`mpe2d_row_broadcast_gpu` ×3、`mpe2d_unify_spec_lev_zx_gpu` ×2、`tranrs1_gpu`、`trandv_gpu_cuda_graph`、`join1rs_gpu`、`ndslfv_*` 等）的 enter_data/host_data **總數**把某個固定配額耗盡了，而 `tendget_gpu` 在 `iteration=2` 開頭的呼叫只是「配額被耗盡後，第一個新提出請求的人」——換誰先請求都會炸，不是 `tendget_gpu` 本身的問題。

**結論**：這一層看起來**不是 GEPS 原始碼的邏輯錯誤**，也**不是靠改單一 subroutine 能修的**（前 3 層都是真的程式碼問題，這層目前所有證據——包含這次consolidation 實驗的負面結果——都指向底層 ROCm/HSA runtime 或這台**共享/分割過的 MI300X 分區**本身的固定資源配額，且是**整輪 iteration 累積**的效應，不是單一 subroutine 的事）。

**已試過、rocprof 工具鏈本身壞掉（2026-09-09）**：想直接量測 HSA queue/signal 用量來證實/推翻上面的假設，結果 ROCm 7.2.2 這台機器上兩代 profiling 工具都用不了：
- `rocprofv3`（新版）：**✅ 2026-09-10 在 SPX 新環境上重新確認，仍然完全重現**。不管開什麼 trace flag，只要載入就在**遠早於任何 GPU 活動之前** SIGSEGV，是它自己的 OMPT 掛鉤機制跟 amdflang libomp runtime 衝突的工具 bug。實際 stack（ROCm 7.2.2）：`SIGSEGV (@0x2a8) @ pthread_mutex_lock` ← `omp_get_num_devices` ← `ompt_post_init` ← `__kmp_do_middle_initialize()` ← `__kmp_middle_initialize` ← `__kmp_api_omp_get_num_places` ← `blas_get_cpu_number` ← `gotoblas_init`（OpenBLAS 初始化階段）。**注意 `rocprofv3` 還需要 `libdw1t64` 才載入得起來**，否則會先以 `libdw.so.1: cannot open shared object file` 失敗（見 §5 的容器套件問題）。
- `rocprof`（舊版 v1）：~~HSA API 呼叫一筆都沒記錄到，roctracer 攔截機制失效~~ **❌ 這個結論是錯的，2026-09-10 已推翻——見下方訂正。**
- ~~**結論：這條「直接量測 queue/signal」的路徑被工具本身擋死**~~ **❌ 已推翻，見下方訂正。**

**★ 訂正（2026-09-10）：`rocprof` v1 的 HSA API 攔截其實完全正常，「直接量測 queue/signal」這條路是通的。**

當初的「最小案例」選錯了程式：filter（`job/rocprof_hsa_filter.txt`）列的是 `hsa_queue_create` / `hsa_amd_memory_pool_allocate` / `hsa_signal_create` 等，而 **`rocminfo` 根本不會呼叫這些函式**——它只列舉 agent（`hsa_iterate_agents` / `hsa_agent_get_info` / `hsa_isa_get_info_alt`…），從不建 queue、不從記憶體池配置、不建 signal。**「零筆」是正確行為，不是攔截失效。**

用一支會真的建 queue 的最小 HIP 程式（`hipMalloc` + `hipMemcpy` + 8 次 kernel launch，`hipcc --offload-arch=gfx942`）重跑同一組 `rocprof --hsa-trace -i rocprof_hsa_filter.txt`，計數器全部抓得到：

| HSA API | Calls |
|---|---|
| `hsa_queue_create` | 1 |
| `hsa_amd_memory_pool_allocate` | 4 |
| `hsa_amd_memory_pool_free` | 1 |
| `hsa_signal_create` | 9 |
| `hsa_amd_signal_create` | 84 |
| `hsa_amd_memory_async_copy_on_engine` | 2 |

**⚠️ 這段訂正本身只對了一半——2026-09-10 下午再訂正一次：**

「攔截器沒壞」是對的（用正常結束的程式證實過）。但把 GEPS 跑批的零筆歸因於「最小案例選錯」是**錯的**——那個歸因只對 `rocminfo` 那個測試成立。**GEPS 跑批之所以零筆，真正的原因是它每次都以 `SIGABRT` 結束。**

實測（同一支 HIP 小程式，只改結束方式）：

| 結束方式 | HSA API 筆數 | KERNEL 筆數 | 有無 `hsa_stats.csv` |
|---|---|---|---|
| 正常 `return 0` | 40 | 1 | 有 |
| **`abort()`** | **0** | 1 | **完全沒有** |
| `SIGTERM` | 39 | 1 | 有 |
| `SIGINT` | 39 | 1 | 有 |

**`rocprof` v1 的 HSA API trace 是記憶體緩衝、只在正常結束（或可攔截的訊號）時寫出；`SIGABRT` 會讓整份遺失。** kernel/copy trace 因為是邊跑邊寫進 `rpl_data_*` 所以留得下來——這正好對應 GEPS 那次的產出（`trace.db` 只有 `KERN` 與 `COPY`，沒有 `HSA` 表）。`mpiexec`、`--timestamp on`、OpenMP offload 都已逐一排除，不是原因。

**還有第二個障礙讓 rocprof 不適用於這個洩漏**：未加 profiling 時 tau=0 輸出階段 VRAM 只用 ~36 GB，掛上 `rocprof --hsa-trace` 後同一階段就衝到 205.7 GB 全滿並提早 abort（連 `iteration=1` 都到不了）。**profiler 自身的配置量遠大於待測目標**，trace 會被它自己的開銷淹沒。

**影響**：
1. **AMD 工單原本的「Issue 3：舊版 rocprof HSA 攔截完全失效」仍然不該送出**——攔截器確實沒壞，但正確的說法是「trace 在 abort 時會遺失」，這比較像是可以順帶回報的使用性問題，不是攔截失效。
2. **查洩漏不要用 rocprof**。已改用 `LIBOMPTARGET_INFO=8`（見下方「洩漏追查結果」一節），開銷小、且能直接歸戶到原始檔。
3. 若將來真的需要 rocprof 的 HSA trace，**必須讓程式正常結束或用 SIGTERM/SIGINT 收尾**，不能讓它 abort。
4. `rocprofv3` 的 SIGSEGV 是**另一回事、仍然成立**（見下）。

**★ 但意外撿到一個有用的新線索**：用 rocprof v1（雖然抓不到 HSA API，但 kernel/copy trace 有效、確實加了 instrumentation overhead）重跑，**同一個 `HSA_STATUS_ERROR_OUT_OF_RESOURCES` 崩潰提早發生了非常多**——不是 `iteration=2` 開頭，而是**連 `iteration=1` 都還沒開始、初始 tau=0 輸出階段**就炸了。這個結果：
- 進一步排除「問題只跟 `tendget_gpu` 的第二次呼叫綁死」——同一種資源耗盡可以在完全不同的呼叫點發生。
- 把假設從「固定總配額用完」修正為更精確的版本：**很可能是配置與回收速率不平衡**（某種非同步資源在被系統實際歸還之前持續堆積，程式跑得越慢/overhead 越大，堆積得越快，門檻就提早觸發）。跟 VRAM 位元組無關、是數量型/時序型資源的大方向沒變，只是「固定總量」這個最簡單的版本可能不夠精確。

下一步：
1. ~~查 SDMA queue / HSA signal 實際使用量~~ **已嘗試，rocprof/rocprofv3 在這台機器上都用不了，見上**。如果還要走這條路，需要先解決 profiling 工具鏈本身的問題（可能要找 AMD 支援、換 ROCm 版本、或找到這台機器上真正能用的 HSA-level tracing 方式），不是本專案這邊能單獨排除的。
2. ~~consolidate tendget_gpu 的 17 個陣列~~ **已嘗試，無效，見上**——不要再花時間在「合併單一 subroutine 的陣列」這個方向，因為證據顯示問題是整輪迭代累積的，不是單一 subroutine 的。如果還要走「減少 enter_data 次數」這條路，必須是**橫跨整個 NNMI pipeline**的大規模重構（多個檔案、風險高很多），不是這次這種局部嘗試。
3. ~~確認這是否為這個 GPU 分區特有的限制~~ **✅ 已完成（2026-09-10），而且答案是「對」**：換成未分割的 SPX 完整卡後崩潰完全消失，詳見本節末「結案」。
4. 回報給 AMD ROCm：現在**有乾淨對照組了**（同一 binary，DPX 會炸、SPX 不會炸），值得回報的還有 `rocprof`/`rocprofv3` 在這台機器上兩代都失效這件事（見上）。注意回報內容要用訂正後的版本——不要再說「DPX 的 SDMA queue 配額比較小」，實測那些 KFD 數字兩邊一模一樣。


**★ 結案（2026-09-10）：換成未分割的 SPX 完整卡後，這個崩潰完全消失。**

系統管理者把這台機器從「1 張 DPX 分區卡」改成「**2 張實體 MI300X，兩張都是 NPS1 / SPX**」之後，用**完全相同的 binary**（`build_rocm_hip/bin/tcogfs.x`，2026-09-09 03:58 建置，比兩個修好的 source 都新，**沒有重新編譯**）、**完全相同的輸入與 job 腳本**重跑 `job/TCo383L72_IC_sample_rocm`：

- 跑批 `06:11:21 → 07:06:48`（約 55 分鐘），log：`job/TCo383L72_IC_sample_rocm.log.20260910_spx`
- **全跑批 `OUT_OF_RESOURCES` / `hsa_amd_memory_lock` / `memory access fault` / `omptarget error` 出現次數 = 0**
- log 第 2001 行 `iteration= 2` 印出後**沒有崩潰**——舊 log（`...log.20260908_iter2crash`）的第 2001 行同樣是 `iteration= 2`，但第 2002 行就是 `HSA_STATUS_ERROR_OUT_OF_RESOURCES`。**同一行號、同一個位置，這次直接跨過去了。**

**這是這個假設能拿到的最乾淨對照**：唯一變數是 GPU 分區模式，binary / 輸入 / 環境變數 / job 腳本全部沒動。因此第 4 層**確定不是 GEPS 原始碼的問題**，§6 這一節裡所有「已排除的修法」（持久化 enter_data、`!$acc wait`、17 陣列 consolidation、兩個 `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD` 實驗）之所以全部無效，是因為它們本來就修不到根因。

**但推測的機制要訂正**：原本猜「DPX 分區壓縮了固定資源配額（SDMA queue 之類）」。實測 SPX 之後，**KFD 的 per-node 配額數字跟 DPX 時期記錄的一模一樣**：

| 屬性 | DPX（舊記錄） | SPX（2026-09-10 實測） |
|---|---|---|
| `num_sdma_engines` | 2 | 2 |
| `num_sdma_queues_per_engine` | 8 | 8 |
| `num_gws` | 64 | 64 |
| VRAM total | 206,141,652,992 | 206,141,652,992 |

所以**不是這些可見屬性造成的**。真正呈現差異的是 **runtime 實際能配置到的量**，用 `rocm-smi` 每 2 秒輪詢 VRAM 得到：

| | DPX（舊） | SPX（2026-09-10） |
|---|---|---|
| 全程 VRAM | **穩定 ~46 GB，幾乎不變**，然後 iteration=2 配置失敗 | 36.4 → 49.8 → 77.1 → 83.3 → **峰值 135.3 GB** |
| 趨勢 | 平的 | **單調爬升**：每分鐘 ~40 MB 持續爬，每輪 NNMI 迭代跳 ~50 GB |

**新的觀察**：每輪 NNMI 迭代淨吃掉約 50 GB 且不歸還——這是真的洩漏，之前在 46 GB 就撞牆所以看不到曲線。`nnmiit=3`（`src/init_block.f90:58` 預設，`namlsts` 沒覆寫），照這個斜率第三次迭代會到 ~185 GB / 206 GB，**很緊但可能塞得下**。這次不是記憶體用完而死（死在第 5 層的 CFL），所以洩漏在 2 輪內還不致命，**但它仍然是個要修的真問題**，只是優先度排在第 5 層後面。

**對外回報的價值**：如果還要送 AMD 工單，現在有了「同一 binary，DPX 會炸、SPX 不會炸」這個乾淨對照組，比原本純推測強很多；`rocprof`/`rocprofv3` 工具鏈失效那兩件事（見下）仍然獨立成立、值得回報。

#### 第 5 層（2026-09-10 首次觸及）：`iteration=2` 垂直平流 deformation CFL 爆掉 —— ⚠️ **本節是初步判斷，其中「嫌疑是手寫 kernel 算錯」的推論已被下一節推翻，保留是為了留下推理過程；結論請直接看下一節「第 5 層根因」**

**症狀**（`job/TCo383L72_IC_sample_rocm.log.20260910_spx` 尾端）：

```
 iteration= 2
  max def_cfl   1.5440282E+04 needs *** steps in   0 rank in advv (gpu) processing
 → exit code 1
```

`***` 是 `i3` 格式溢位。中止機制在 `src/mod_ndslfv_monoadv_gpu.f90` 的 `def_cfl_step_gpu_type2`：`check_max = 1.544e4` → `nstep = int(check_max/0.99) + 1 = 15597`，`nstep > 10` 就 `call exit(1)`。**注意 `call exit(1)` 在 `print *, "Error in ..."` 之前，所以那行永遠印不出來**，也沒有 stack trace——log 尾端只有 mpiexec 的「exit code 1」，不要誤判成靜默崩潰。

**性質跟前 4 層完全不同**：不是崩潰、不是資源耗盡，是**模式自己偵測到數值不合理而主動中止**。`check = abs((dist(i+1,j) - dist(i,j))/del(i,j))` 是垂直方向的 deformation CFL，衝到 1.5e4 代表**出發點距離場 `dd` 已經是垃圾值**。

**已確認的三個佐證（都不是猜測）**：
1. **`def_cfl` 在過去所有 log 裡一次都沒出現過**——`log`、`log.20260908_0118`、`log.20260908_0141`、`log.20260908_iter2crash`、`log.rocprof4`、CPU 版 `log.20260908` 全部 grep 計數為 0。是第一次跑到這裡。
2. **CPU 版跑批完整跑完 `tau=1.000` 且從不觸發這個檢查**，所以不是模式設定或 IC 本身不穩定，是 **GPU 路徑算出來的值不對**。
3. **失敗的 `vertical_cell_advect_gpu` 正是 `src/rocm/` 手寫替換掉的兩支 kernel 之一**（`src/CMakeLists.txt` 的 `geps_acc2omp()` 對 `mod_ndslfv_monoadv_gpu` 額外掛了 `rocm/vertical_cell_advect_gpu.f90` + `rocm/ndslfv_monoadvv_gpu.f90` 當相依，由 `acc2omp.py` 的 `apply_file_fixups()` 整支替換）。已在產生出來的 `build_rocm_hip/src/hip_src/mod_ndslfv_monoadv_gpu.f90:1808` 確認用的就是手寫版（帶 `geps_dbg_vram` 探針呼叫），CFL 呼叫在該檔第 1914 行。**這是高風險嫌疑犯，本來就是「方言太重、翻譯器搞不定所以手寫」的檔案。**

**另一條相關線索**：`bal= 0.` 在 iteration=1 與 iteration=2 的三個 vertical mode **全部都是零**（handoff 舊版已標記為可疑）。NNMI 沒有做出任何修正量，代表送進積分的場是不平衡的——跟垂直平流 CFL 爆掉很可能是**同一個根因的兩個表現**，建議一起查。

**建議的查法（尚未動手）**：
1. 先確認 `bal=0.` 是真是假：CPU 版同一階段的 `bal` 是什麼數量級？（CPU log 有完整序列可比對）
2. 把 `vertical_cell_advect_gpu` 的輸入（`dd`, `ds`）在 iteration=1 / iteration=2 各 dump 一次，跟 CPU 版同點比對——**這是第一次可以做 CPU vs GPU 數值比對的機會**，因為 GPU 終於跑到這裡了。
3. 逐一檢查 `src/rocm/vertical_cell_advect_gpu.f90`（203 行）與被它取代的原版邏輯差異，特別是 reduction / 邊界處理。


#### ★ 第 5 層根因（2026-09-10 已查明）：`tendget_gpu` / `trngra_gpu` / `tranuv_gpu` 每輪 NNMI 洩漏 ~150 GB VRAM

**這一節推翻了上面「第 5 層是 `vertical_cell_advect_gpu` 數值 bug」的初步判斷，也把第 4 層一起解釋掉了。兩層是同一個根因。**

##### 怎麼查到的

在兩個地方加了探針，重編一次就同時回答了兩個問題（相關改動見本節末「工作樹狀態」）：

1. `cmake/acc2omp.py` 的 `initial_gpu` 注入器：在 NNMI 迴圈每個階段呼叫前插入 VRAM 探針（locid 300/310/311/320/321/330/340/350），補上原本 25 分鐘完全沒有儀器的空窗。
2. `src/rocm/vertical_cell_advect_gpu.f90`：在 `!$omp target update from(dd, ds)` 之後，**在 host 端用與裝置 kernel 完全相同的公式重算 CFL**（locid 257），report `max|dd|` 與 `check_max`。這樣就能分辨「輸入真的是壞值」與「裝置端映射/記憶體問題」。

跑批 `07:52:52 → 08:50`（約 57 分鐘），log：`job/TCo383L72_IC_sample_rocm.log.20260910_probe`

##### 洩漏定位（決定性數據）

```
iteration 1
  07:56:14  -> tendget_gpu    free = 157.54 GB
  08:08:31  -> zx_gpu #1      free =  61.89 GB   ← tendget_gpu 一次吃掉 95.65 GB
  08:08:51  -> vartrix_gpu    free =  60.31 GB   (-1.58)
  08:08:51  -> nnmi_gpu       free =  60.31 GB   (0)
  08:15:53  -> trngra_gpu     free =  23.76 GB   (-36.48)
  08:20:08  -> tranuv_gpu     free =   5.95 GB   (-17.81)
iteration 2 / 3
  全程停在 free = 5.92 GB，之後所有階段 delta 都是 0（已經配置不到任何東西）
```

- **每輪 NNMI 淨消耗約 150 GB，且回傳後一個位元組都沒歸還**（`-> zx_gpu #1` 這個探針位在 `tendget_gpu` 回傳之後，free 仍是 61.89 GB）。
- 三個元兇：**`tendget_gpu` −95.65 GB**、**`trngra_gpu` −36.48 GB**、**`tranuv_gpu` −17.81 GB**。
- `zx_gpu` / `vartrix_gpu` / `nnmi_gpu` 合計只有 1.6 GB，**不是問題**。
- `tendget_gpu` 內部再細分：`vram:130`（`ndslfv_pack_gpu`，水平平流）−47.45 GB、`vram:245`（transpose 到垂直平流之間）−41.62 GB。

##### CFL 爆掉的真正原因（`vram:257` 的 host vs 裝置對照）

| | iteration 1 | iteration 2 / 3 |
|---|---|---|
| free VRAM | 91.04 GB | **0.00 GB** |
| `max\|dd\|`（host 重算） | 1.1420 | **0.0000** |
| `check_max`（host 重算） | 0.0450 | 0.0000 |
| 裝置端 `nstep` | 1 | 1 |
| `ds(1,1)` / `dsmin` | 3.614149 / 0.636140 | **完全相同** |

**`dd` 在 iteration 2 之後整個變成零**，而同一個 subroutine 裡、由已駐留的 `ssi` 算出來的 `ds` **完全正常且與 iteration 1 逐位元相同**。所以不是除以零、不是 `ds` 壞掉，是 **`dd` 的來源（`wwi` = `pdot`）在記憶體歸零的情況下沒有進到裝置上**。

**最關鍵的一點：這個症狀不是決定性的。**
- 09-10 **未加探針**的那次跑批（`...log.20260910_spx`）：同樣的 0 free 狀態，裝置端算出 `check_max = 1.5440282E+04` → `nstep=15597 > 10` → `call exit(1)`。
- 09-10 **加了探針**這次：同樣的 0 free 狀態，同一段程式算出 `dd` 全零 → `check_max = 0` → `nstep=1` → **不崩潰，繼續往下算**。
- 唯一差別是我加的探針帶了 8 個 `hipDeviceSynchronize()`（`geps_dbg_gpu_log_` 走 `dbg_ndjson(..., do_sync=1)`）。**只是多了同步就換了一種垃圾值——這是 Heisenbug，算術錯誤不會有這種行為。**

**結論：`vertical_cell_advect_gpu` 沒有數值 bug。CFL 那個 1.5e4 是「在 0 可用 VRAM 下執行」的下游症狀。不要再去改那支手寫 kernel。**

##### 這次的結局，以及一個新的風險

不崩潰**不代表比較好**。這次跑批的實際結果：

```
iteration= 1 / 2 / 3        三輪 NNMI 全部「跑完」（但 2、3 輪是拿 dd=0 在算）
forcast begin tau= 0 to tau= 1
tau=    0.167,     Timing=   0.000 elapse seconds
surf pres tend rms(GPU) = NaN  mb/hrs          ← 積分結果是 NaN
"PluginInterface" error: ... HSA_STATUS_ERROR_OUT_OF_RESOURCES
omptarget error: Call to getTargetPointer returned null pointer
Signal: Aborted (6)  →  error occured: fct model fail !!
```

兩個要記住的點：
1. **`HSA_STATUS_ERROR_OUT_OF_RESOURCES` 回來了**——就是第 4 層那個錯誤，只是這次拖到積分階段才發作。**第 4 層與第 5 層確定是同一個根因。**第 4 層當初精確定位到 `tendget_gpu.f90` 頂部的 `enter data` 是對的，但**那一行是受害者不是元兇**：它是「洩漏把記憶體吃光後，第一個提出大額請求的人」。§6 第 4 層裡「合併 17 個陣列無效」的負面結果現在也完全說得通——問題從來不是那 17 個陣列。
2. **靜默算錯比崩潰危險**。加了同步之後模式不再中止，而是帶著零值垂直平流一路算到積分、產出 `NaN`。如果哪天 NaN 沒有出現，就會得到「看起來合理但完全錯誤」的預報場。**在洩漏修好之前，任何 GPU 跑出來的數值都不可信。**

##### 下一步（取代舊的第 5 層建議）

1. **查 `tendget_gpu` 的 95.65 GB 去哪了**。它有成對的 `enter data` / `exit data`（`nvidia_tendget_gpu.f90` 頂部 17 個 scratch 陣列 + 各 downstream subroutine 自己的）。要確認的是**哪些 `exit data` 沒有真的釋放**。
2. **`rocprof --hsa-trace` 現在可以用了**（見上面的訂正），用 `job/rocprof_hsa_filter.txt` 比對 `hsa_amd_memory_pool_allocate` 與 `hsa_amd_memory_pool_free` 的次數與大小差額，直接列出沒被釋放的配置。**這是目前最直接的手段**，而且不必再改原始碼。
3. 細分 `tendget_gpu` 內部：`ndslfv_pack_gpu`（−47.45 GB）與 transpose 段（−41.62 GB）各自再加探針。
4. `trngra_gpu`（−36.48 GB）與 `tranuv_gpu`（−17.81 GB）同樣要查，它們的量也足以單獨壓垮一輪迭代。
5. 修好之後才回頭確認 `bal= 0.`（三輪 NNMI 的三個 vertical mode 全為零）是不是真實結果——在目前的記憶體狀態下這個值同樣不可信。

##### 工作樹狀態（這次新增的除錯改動，**都還沒 commit**）

| 檔案 | 改動 | 是否該保留 |
|---|---|---|
| `cmake/acc2omp.py` | `initial_gpu` 注入器加了 per-anchor 出現次數計數器與 6 個新錨點（locid 300/310/311/320/321/330/340/350） | 除錯用；洩漏修好後可移除。既有 7 個探針（148/168/180/199/221/234/235）都完整保留 |
| `src/rocm/vertical_cell_advect_gpu.f90` | 新增 host 端 CFL 重算探針（locid 257）＋ 3 個區域變數 `ddmax`/`hostchk`/`chk` | 除錯用；很便宜，建議留到洩漏修好為止 |

> ⚠️ 這些探針**會改變程式行為**（`geps_dbg_gpu_log_` 內含 `hipDeviceSynchronize()`）。要重現「CFL=1.5e4 崩潰」必須把 `initial_gpu` 的新探針拿掉；要看到「靜默算成 NaN」則保留。做數值驗證時**兩者都要移除**。


#### ★★ 洩漏追查結果（2026-09-10 下午）：NNMI 的消耗不在 OpenMP 映射層（另一個「積分階段洩漏」的結論已自行推翻）

用 `LIBOMPTARGET_INFO=8`（只記錄映射表的 create/remove，每筆帶 `Size` 與 `Name=<原始檔>`）跑完整一輪，log 1.34 GB（**務必寫到 `/tmp` overlay，不要寫 NFS——workspace 只剩 ~57 GB**；分析工具在 `GEPS/qa/rocm_repro/`），結局與前次相同（`NaN` → `OUT_OF_RESOURCES` → abort）。

##### 決定性數據：未歸還的映射位元組 vs 模式階段

```
phase marker              map events   outstanding(GB)   entries
in outflds for tau               531        31.526          531
iteration= 1                    7995        44.979          595
iteration= 2                 1045583        44.979          595   ← 完全沒變
iteration= 3                 2083167        44.979          595   ← 完全沒變
beginning integration        3121123        77.781          795
forcast begin                3121132        78.479          804
(log 結束)                   5027022       153.656         1016
```

全程 `Creating new map entry` 2,514,019 次 / `Removing map entry` 2,513,003 次，差 1,016 筆。

##### 問題 A（**目前的 blocker**）：NNMI 階段的消耗**不是映射洩漏**

**三輪 NNMI 之間，未歸還的映射位元組一動也不動（45.0 GB / 595 筆），期間卻發生了約 200 萬次映射事件，而 HIP 可見的 free VRAM 從 ~157 GB 掉到 ~0。**

也就是說：`enter data` / `exit data` **有正確配對**，映射表是乾淨的，但**底下的裝置記憶體沒有還給 driver**。消耗發生在 libomptarget 的記憶體管理器（memory manager pool）或 HSA 配置器那一層，不在 GEPS 的 OpenMP 指令層。

這一口氣解釋了先前好幾個看不懂的現象：
- **為什麼把 `tendget_gpu` 的 17 個陣列合併成 1 個完全沒用**（§6 第 4 層）——映射「數量」從來不是問題。
- **為什麼 `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0` 會讓行為劇烈改變**（提早在別處炸）——那個開關正是繞過這一層。
- **為什麼 `rocm-smi` 與 `hipMemGetInfo` 對不上**（前者 135 GB、後者只剩 2 MB）——池子持有的記憶體兩邊的帳算法不同。
- **為什麼探針量到「`tendget_gpu` 吃掉 95.65 GB 且回傳後沒歸還」**——它的 `exit data` 確實執行了，記憶體只是回到池子沒回到 driver。

> ⚠️ **這代表問題 A 很可能不是 GEPS 原始碼能修的**，方向要轉往 ROCm/libomptarget 的記憶體管理行為。**不要再花時間調整 GEPS 這邊 `enter/exit data` 的寫法或數量——已經有三個獨立證據顯示那個方向無效。**

**下一步建議**：
1. 掃 `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD` 的中間值（先前只試過 `0` 與 `4GB`，兩個極端）。
2. 試 `LIBOMPTARGET_INFO=32`（資料搬移）或 libomptarget 的 debug build，觀察池子的成長。
3. 這是很適合回報給 AMD 的題目：**映射表平衡、`exit data` 正常，但裝置記憶體單調不歸還直到耗盡**。比原本那份 Issue 1 更具體。

##### 問題 B：~~積分階段有真的映射洩漏~~ ❌ **不存在，2026-09-10 稍晚推翻**

> 本小節原本主張 `diabat_gpu.f90`（create 145 / remove 3）、`cwa_rrtmg_gpu.f90`（59/2）等
> 在積分階段大量洩漏映射。**這個結論是錯的，已完整推翻。保留這段是為了不讓下一個人重踩。**

**推翻的方法**：把這些檔案的 create/remove 事件標上 log 行號，看它們發生在什麼時候。

```
lines 5028787-5028936   CRE x145   diabat_gpu.f90
lines 5028936-5028973   REM x3     diabat_gpu.f90
lines 5028973-5029031   CRE x58    cwa_rrtmg_gpu.f90
lines 5029031-5029033   REM/CRE/REM x1
(log 到 5029082 行結束 = abort)
```

**`diabat_gpu` 一輩子只被呼叫過一次，而且程式死在它裡面。** 那 145 筆是它頂部
`enter data` 建立的正常工作集；exit data 在 `src/diabat_gpu.f90:4129`，**根本還沒執行到**。
3 筆 remove 是內層 `ozplin`（1336→1368）與 `slope_data`（1721→1755）自帶的小配對。
`cwa_rrtmg_gpu` 是從 diabat 裡呼叫的，同樣死在半途。

其餘被列為嫌疑的檔案也一樣，全部不是洩漏：

| 檔案 | create/remove | 事件所在行號 | 實際是什麼 |
|---|---|---|---|
| `phygrid.mod` | 128 / 0 | 3122933–3123085 | 積分開始時一次性映射，長生命期 |
| `raddiag.mod` | 30 / 0 | 3122941–3123107 | 同上 |
| `nvidia_mod_spec_cuda_graph.f90` | 22 / 0 | 925–1037 | init 階段一次性映射 |
| `mod_grid.f90` | 54 / 8 | 首見 787 | init 階段一次性映射 |
| `cwa_grrad_gpu.f90` | 13 / 0 | 5029034–5029047 | 死亡當下執行中的工作集 |
| `diabat_gpu.f90` | 145 / 3 | 5028787–5028973 | 死亡當下執行中的工作集 |

**沒有任何一支呈現「反覆 create/remove、未歸還量持續成長」的洩漏樣態。**

**當初錯在哪**：`create ≫ remove` 這個指標，只有在該段程式**被呼叫多次且每次都正常結束**的前提下才代表洩漏。
`diabat_gpu` 只跑了一次而且死在中途，`145:3` 正是預期結果。
原始碼也佐證：`src/diabat_gpu.f90` 的 27 對 `enter/exit data` 完整配對（產生檔 26:26），
`cwa_rrtmg_gpu.f90` 則連一個 `enter/exit data` 都沒有（全是 `!$omp target` 的隱式映射），
**兩支都沒有漏寫 exit data，沒有東西可以修**。

**教訓（下次用這個 log 分析時務必記住）**：程式以 abort 結束時，
「死亡當下仍存在的映射」= 合法的長生命期映射 ＋ 執行中呼叫的工作集，
**兩者都不是洩漏**。要判定洩漏，必須看**同一段程式跨多次完整呼叫**時未歸還量是否單調成長
（`omptinfo_timeline.py` 就是為此而寫的）。

##### 目前唯一成立的結論

**只有問題 A 是真的**：三輪 NNMI 之間映射表完全持平（45.0 GB / 595 筆不變、約 200 萬次事件全部正確配對），
而 HIP 可見的 free VRAM 從 ~157 GB 掉到 ~0。消耗發生在 libomptarget 記憶體池／HSA 配置器那一層，
**不在 GEPS 的 OpenMP 指令層，很可能不是 GEPS 原始碼能修的**。


#### ★★★ 第 6 層（2026-09-10 查明並 ✅ **已修好**）：RCCL 集合通訊之後的競態 —— acc2omp 把 `async` 子句丟掉造成的

**這是目前最重要的發現：一個真實、可重現、有明確修法的 GEPS/toolchain bug，而且它一直被 libomptarget 的記憶體池掩蓋著。**

##### 怎麼找到的：掃 `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD`

| 門檻值 | 結果 |
|---|---|
| `0`（池全關，2026-09-09 試過） | ❌ `iteration=1` 內 memory access fault |
| **`65536`（64 KB）** | ❌ 同上，345 秒 |
| **`1048576`（1 MB）** | ❌ 同上，345 秒 |
| **`16777216`（16 MB）** | ❌ 同上，345 秒 |
| **`268435456`（256 MB）** | ❌ 同上，345 秒 |
| `4294967296`（4 GB） | 與預設完全相同（洩漏但不 fault） |
| 預設 | 洩漏但不 fault |

四個新測的值**全部在 345 秒、同一個位置、同一種錯誤**崩潰，崩潰前的「Last 8 kernels」也完全一致（`tranrs_gpu_cuda_graph_` / `mpe_transpose_rs_sp_gpu_` / `rfftmlt_loop_identical_cuda_graph_` / `join1rs_gpu_` / `gridnl_hybrid_ndsl_gpu_refactor_`）。**這不是隨機現象，是穩定可重現的潛伏 bug。**

##### `OFFLOAD_TRACK_ALLOCATION_TRACES=true` 直接指出兇手

用 16 MB 門檻（345 秒就觸發，比預設組態快十倍）搭配 `OFFLOAD_TRACK_ALLOCATION_TRACES=true` 與 `OFFLOAD_TRACK_NUM_KERNEL_LAUNCH_TRACES=8` 重跑：

```
OFFLOAD ERROR: memory access fault by GPU 2 at virtual address 0x1065134b2000
Device pointer 0x1065134b2000 points into prior host-issued allocation:
Last deallocation:
  ... __tgt_target_data_end_mapper + 547
Last allocation of size 341066880 for host pointer 0x7ff97d73a890
                                  -> device pointer 0x106505a00000:
  ... __tgt_target_data_begin_mapper + 547
```

出錯位址落在那塊 **341,066,880 bytes（325 MiB）** 的配置內（起點 `0x106505a00000`，偏移 229 MB），而**該配置已經被釋放**。用 `llvm-symbolizer` 解析堆疊：

```
0x120bd3f (配置點)  -> mpe_transpose_rs_sp_gpu_
0x120c05a (釋放點)  -> mpe_transpose_rs_sp_gpu_
0x1219634           -> tranrs_gpu_cuda_graph_
0x11f1afb           -> tendget_gpu_
0x10f9886           -> initial_gpu_
```

呼叫鏈：`initial_gpu` → `tendget_gpu` → `tranrs_gpu_cuda_graph` → **`mpe_transpose_rs_sp_gpu`**。

##### 根因：acc2omp 把 `async(async_id)` 全部丟掉，破壞了順序保證

`src/nvidia/mpe_transpose_rs_gpu.f90`（NVIDIA 原始碼，**正確**）：

```fortran
!$acc enter data create(swork) async(async_id)
!$acc parallel loop collapse(4) async(async_id)
   ... 填 swork ...
call nccl_alltoall(swork, len_tr*lev, rbuf, len_tr*lev, comm, nsize, async_id)
!$acc exit data delete(swork) async(async_id)
```

四個動作**全在同一條 async queue 上**，所以 `delete` 一定排在集合通訊完成之後。

`build_rocm_hip/src/hip_src/nvidia_mpe_transpose_rs_gpu.f90`（acc2omp 產生，**有 bug**）：

```fortran
!$omp target enter data map(alloc:swork)                 ← async 被拿掉
!$omp target teams distribute parallel do collapse(4)    ← async 被拿掉
   ... 填 swork ...
call nccl_alltoall(swork, ..., async_id)                 ← 仍然非同步！
!$omp target exit data map(delete:swork)                 ← async 被拿掉 → 立即執行
```

關鍵在於 **`nccl_alltoall` 不是 OpenMP 構造，acc2omp 不會動它**。它在 `src/nvidia/helper.f90:66` 用
`stream = acc_get_cuda_stream(async_id)` 拿到 stream，發出 `ncclSend`/`ncclRecv` 之後**直接返回、不做任何同步**。
於是翻譯後的 `exit data delete(swork)` 在 host 端**立刻**執行，把 RCCL 還在讀的 buffer 釋放掉 → **use-after-free**。

**為什麼平常看不到**：預設門檻下，釋放的 325 MiB 區塊留在 libomptarget 的記憶體池裡，
位址仍然映射有效，RCCL 讀到的是「已釋放但還沒還給 driver」的記憶體，不會 fault，
只是**資料可能已經被後續配置覆寫**。把門檻調低，區塊真的還給 driver，存取就爆。
**這正好解釋了 §6「洩漏追查結果」問題 A 的所有現象**——池子不歸還記憶體，
既是耗盡的原因，也是這個 bug 的遮罩。

##### 影響範圍（很收斂）

掃過所有呼叫 `nccl_*` 的檔案，**只有 `src/nvidia/mpe_transpose_rs_gpu.f90` 有「nccl 呼叫後 6 行內接 exit data」這個模式，共 3 處**
（`mpe_transpose_rs_sp_gpu` ×1、`mpe_transpose_rs_gpu` ×2）。
`src/nvidia/mpe_transpose_rs1_gpu.f90` 與 `mpe_transpose_sr_gpu.f90` 也呼叫 `nccl_alltoall`，
但 exit data 距離較遠，**仍需逐一確認**。

##### 建議修法（尚未實作）

**在 toolchain 層修，不要改 `src/nvidia/`**（符合 §7 的工作規矩）。已經有現成的
`src/rocm/geps_acc_wait.f90`（`geps_acc_wait(async_id)` → `geps_hip_wait`）可用。

方向：在 `cmake/acc2omp.py` 裡，**當它把 `!$acc exit data ... async(id)` 的 async 子句丟掉時，
在該 `exit data` 之前插入 `call geps_acc_wait(id)`**。這樣就把原本「排在同一 queue 上」的順序保證
換成一次顯式等待，語意等價且安全。

- 保守做法：只針對 `mpe_transpose_rs_gpu` 這個檔案做 `apply_file_fixups()` 點名修補（風險最低）。
- 一般做法：對所有「原本帶 async 的 exit data」都插入等待（更完整，但會多出一些不必要的同步，
  且 `src/nvidia/intgrt_gpu.f90` 單檔就有 32 處，效能影響要評估）。

**驗證方式**：修完用 `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=16777216` 重跑——
這個組態現在 345 秒就會 fault，如果修對了就應該能跨過去。**這是個很便宜的迴歸測試，不必等 55 分鐘。**


##### ✅ 修法（2026-09-10 已實作並驗證通過）

**改在 `cmake/acc2omp.py`（build 層），沒有動 `src/nvidia/`——符合 §7 的工作規矩。**

規則很簡單：**任何帶 `async_id` 引數的 `call nccl_*(...)` 之後，立刻插入 `call geps_acc_wait(async_id)`。**
（`geps_acc_wait` 來自 `src/rocm/geps_acc_wait.f90`，acc2omp 本來就用它翻譯 `!$acc wait`，不需要 interface。）
`nccl_check_helper`（NCCLCHECK 巨集）沒有 `async_id`，自然不會被插入。

實作位在 `apply_file_fixups()` 開頭，對**所有**翻譯檔案生效，共插入 12 處：

| 檔案 | nccl 呼叫 | 插入 wait |
|---|---|---|
| `mpe_transpose_rs_gpu.f90` | 4 | 4 |
| `mpe_transpose_sr_gpu.f90` | 1 | 1 |
| `mpe_transpose_rs1_gpu.f90` | 1 | 1 |
| `mpe2d_gpu.f90` | 6 | 6 |
| `helper.f90`（定義處，非呼叫） | 0 | 0 ✓ |

###### ⚠️ 第一版修法失敗的教訓（不要重蹈）

**第一版只在 `exit data` 前插等待，而且只針對 `mpe_transpose_rs_gpu`——測試沒過。**
重跑後 fault 仍在 332 秒發生，但**錯誤原因從 `Unknown (0)` 變成 `Write access to a read-only page`**，
最近的 kernel 也從 `gridnl_hybrid_ndsl_gpu_refactor_` 變成 `trngra3_gpu_cuda_graph_` / **`mpe_transpose_sr_sp_gpu_`**。

原因是問題比「exit data 太早」更廣。看 `src/nvidia/mpe_transpose_sr_gpu.f90`：

```fortran
!$acc enter data create(rwork) async(async_id)
call nccl_alltoall(sbuf, ..., rwork, ..., async_id)   ! 非同步「寫入」rwork
!$acc parallel loop collapse(4) async(async_id)        ! 「讀」rwork ← 翻譯後同步，搶跑
!$acc exit data delete(rwork) async(async_id)          ! 第一版只擋到這裡
```

**跟在 nccl 之後、碰到它 buffer 的任何構造都會搶跑**，不只是 `exit data`：
- 讀 recv buffer 的 kernel → 讀到 RCCL 還沒填完的資料
- 釋放 send buffer 的 `exit data` → use-after-free

`mpe_transpose_rs1_gpu.f90` 也是同樣模式（nccl 之後接 `parallel loop` 讀 `rwork`）。
**當初用「nccl 呼叫後 6 行內有 exit data」這個窗口去掃描，漏掉了這兩個檔案**——
掃描條件本身就不對，正確的判準是「nccl 之後的所有東西」。

###### 驗證結果

用 §6 那個便宜的迴歸測試（`LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=16777216`，
修法前**必定**在 332–345 秒 fault）：

| | 修法前 | 第一版修法 | **通用修法** |
|---|---|---|---|
| fault | 345s，`Unknown (0)` | 332s，`Write access to a read-only page` | **無（跑過 16 分鐘仍為 0）** |
| 進度 | 死在 `iteration= 1` | 死在 `iteration= 1` | **`iteration= 1` + `vertical mode l= 1/2/3` 全部通過** |

**這是第一次有東西真正修好了這個 fault。**

###### ⚠️ 但修好 fault ≠ 修好模式：`NaN` 依舊

讓這個組態繼續跑完整一輪（`elapsed=3658s`，約 61 分鐘），結局**跟修法前完全一樣**：

```
iteration= 1 / 2 / 3        三輪 NNMI 全部通過（修法前這個門檻連 iteration=1 都撐不過）
forcast begin tau= 0 to tau= 1
tau=    0.167
surf pres tend rms(GPU) = NaN  mb/hrs        ← CPU 版同一行是 .4380739682505907
"PluginInterface" error: ... HSA_STATUS_ERROR_OUT_OF_RESOURCES
Signal: Aborted (6)  →  error occured: fct model fail !!
```

**第 6 層是必要條件、不是充分條件。** 數值錯誤（`NaN`）與記憶體耗盡都還在，另有原因。

###### 記憶體行為：有改善，但歸因要小心

| 探針 | 基準（預設門檻、無修法） | 本次（16 MB 門檻 ＋ 修法） |
|---|---|---|
| `-> tendget_gpu` → `-> zx_gpu #1` | **−95.65 GB** | **−76.29 GB** |
| `-> trngra_gpu` | −36.48 GB | −36.14 GB |
| `-> tranuv_gpu` | −17.81 GB | −17.81 GB |
| iteration 1 結束時 free | **5.92 GB** | **26.65 GB** |

`tendget_gpu` 少吃 19.4 GB、iteration 1 結束多留 20.7 GB 餘裕。
**但這兩次跑批同時改了兩個變數（修法 ＋ 門檻），不是乾淨的對照組**——
門檻本身就會改變池子行為，所以**不能把改善單獨歸功於修法**。
而且「不歸還」的本質沒變：iteration 2、3 的各階段 delta 全是 0，一路停在 26.65 GB。

**要乾淨對照，下一步應該跑一次「預設門檻 ＋ 修法」。**

###### ✅ 乾淨對照已完成（2026-09-11）：修法對記憶體與 `NaN` **完全沒有影響**

跑了「預設門檻 ＋ 第 6 層修法」（`elapsed=3367s`，log：`job/TCo383L72_IC_sample_rocm.log.20260911_control`），
與 09-10 的基準線只差修法這一個變數。探針數字**逐位元相同**：

| 探針 | 基準（無修法） | 對照（有修法） |
|---|---|---|
| `-> tendget_gpu` → `-> zx_gpu #1` | −95.65 GB | **−95.65 GB** |
| `-> vartrix_gpu #1` | −1.58 GB | **−1.58 GB** |
| `-> trngra_gpu` | −36.48 GB | **−36.48 GB** |
| `-> tranuv_gpu` | −17.81 GB | **−17.81 GB** |
| iteration 1 結束 free | 5.92 GB | **5.92 GB** |
| iteration 2 / 3 各階段 delta | 全為 0 | **全為 0** |
| 積分第一步 | `NaN` | **`NaN`** |
| 結局 | `OUT_OF_RESOURCES` → abort | **同上** |

**結論：第 6 層的修法只做了一件事——消除低門檻下才會出現的 memory access fault。
它對記憶體耗盡與 `NaN` 毫無幫助。**

這同時推翻了先前那個誘人的推測（「第 6 層的競態可能就是第 4/5 層怪現象的來源」）：
**不是**。先前 16 MB 那次看到的「`tendget_gpu` 只吃 76.29 GB、結束多留 20.7 GB」，
**完全來自門檻變更，與修法無關**。

> 📌 這也是為什麼當時刻意不把改善歸功於修法——同時改兩個變數的實驗不能用來歸因。

> 📌 另一個觀察：帶 `OFFLOAD_TRACK_ALLOCATION_TRACES=true` 重跑第一版修法時，18 分鐘都沒有 fault，
> 但同一顆 binary 不帶該旗標時 332 秒就 fault。**allocation tracking 的開銷會掩蓋這個競態**，
> 跟先前探針造成的 Heisenbug 同類。**用它定位問題可以，但不能用它驗證修法。**

##### 順帶記下的工具

- `OFFLOAD_TRACK_ALLOCATION_TRACES=true`：fault 時印出該位址屬於哪次配置／釋放的堆疊。**查這類問題的首選。**
- `OFFLOAD_TRACK_NUM_KERNEL_LAUNCH_TRACES=8`：印出最後 8 個 kernel launch 的堆疊。
- `llvm-symbolizer`（`/opt/rocm/lib/llvm/bin/llvm-symbolizer`）：把上面那些裸位址解成函式名。
  用法：`echo 0x120bd3f | llvm-symbolizer --obj=build_rocm_hip/bin/tcogfs.x --functions=linkage`
  （目前 binary 沒有行號資訊，只能解到函式名；要行號需加 `-g`）。
- **用低門檻當快速重現手段**：`LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=16777216` 讓這個 bug 從
  「55 分鐘後表現為 NaN／耗盡」變成「345 秒後明確 fault」。


#### ★★★ 第 7 層（2026-09-11 查明）：`NaN` 的根因 —— graph capture 區塊內的 OpenMP kernel 沒有被錄進 graph

**`bal= 0.` 從一開始就被標為「可疑」，結果它就是主因。這一層完全不需要跑新的實驗就查出來了——兩份現成的 log 對比即可。**

##### 零成本的決定性對比

`job/TCo383L72_IC_sample_rocm.log.20260911_control`（GPU）vs `job/TCo383L72_IC_sample_rocm_cpu.log.20260908`（CPU）：

| | CPU | GPU |
|---|---|---|
| 初始 `pt`（max/min） | 0.92182 / -7.66227 | 0.92182 / -7.66227 ✅ 一致 |
| iteration 1 的三個 `bal` | 1.76e-05 → 6.39e-06 → 1.65e-06 | **0. / 0. / 0.** |
| iteration 1 **之後**的 `pt` | **0.91441 / -7.67721**（變了） | **0.92182 / -7.66227**（逐位元不變） |
| iteration 2 的 `bal` | 4.09e-07 … | 0. |
| iteration 3 的 `bal` | 4.93e-08 …（明確收斂） | 0. |
| iteration 3 之後的 `pt` | 0.91483 / -7.67639 | 0.92182 / -7.66227（仍不變） |

**GPU 的 NNMI 是個完全的 no-op。** 模式帶著未經正規模態初始化的不平衡場進入積分，
第一步就 `surf pres tend rms(GPU) = NaN`（CPU 同一行是 `.4380739682505907`）。

##### 根因：跟第 6 層同源——acc2omp 丟掉 `async`，OpenMP kernel 因此不參與 stream capture

`src/nvidia/nnmi_gpu.f90` 的結構是「**整段計算包在 HIP graph capture 裡，只在第一次呼叫時錄製，之後靠 replay 執行**」：

```fortran
if (.not. cg_created) then
   call accx_async_begin_capture(async_id)
      ... 第一組 dgemm（hipBLAS，在 subwrk_stream 上，靠 event fork/join 併入 capture）
      !$acc parallel loop ... async(async_id)   ← 算 bal，並把 wrk 做 1/e_t 變換
      ... 第二組 dgemm，把 wrk 轉回 x
   call accx_async_end_capture(async_id, cg_graph)
   cg_created = .true.
end if
call accx_graph_launch(cg_graph, async_id)
```

翻譯後（`build_rocm_hip/src/hip_src/nvidia_nnmi_gpu.f90`）：

```
67:  call accx_async_begin_capture(async_id)
82:  call dgemm('T','N',...)                           ← 錄進 graph，不執行
94:  !$omp target teams distribute parallel do ...     ← ❌ 沒有 nowait / depend
125:    bal(mf, L) = bal_t
140: call dgemm('N','N',...)                           ← 錄進 graph，不執行
149: call accx_async_end_capture(async_id, cg_graph)
154: call accx_graph_launch(cg_graph, async_id)
```

第 94 行的 kernel **沒有 `nowait`/`depend`，走的是 libomptarget 自己的 stream，不是被 capture 的那條**。後果有兩個，剛好解釋全部三個觀察：

1. **capture 期間它立刻執行**，此時 dgemm 只被錄製、還沒跑，`wrk` 是空的
   → `bal_t = 0` → **`bal(mf,L) = 0`**（觀察 1 ✅）
2. **它沒有被錄進 graph**，replay 時那個 `1/e_t` 變換整個消失
   → graph 只剩 `dgemm('T')` 後接 `dgemm('N')`，而 `evec` 是正交矩陣
   → `x ← evec·(evecᵀ·x) ≈ x`，修正量實質為零
   → **`pt` 逐位元不變**（觀察 2、3 ✅）

**這跟第 6 層是同一個病根的兩種表現**：acc2omp 把 `async(async_id)` 丟掉之後，
被翻譯的 OpenMP 構造就脫離了原本的 stream 語意——
第 6 層是「脫離順序保證 → 競態」，第 7 層是「脫離 capture → kernel 掉出 graph」。

##### 影響範圍（3 個檔案，都是核心路徑）

| 產生檔 | capture 區塊 | 區塊內的 OpenMP target kernel |
|---|---|---|
| `nvidia_nnmi_gpu.f90` | 1 | 1（算 `bal`＋`wrk` 的 `1/e_t` 變換） |
| `nvidia_tendget_gpu.f90` | 1 | 2 |
| `nvidia_zx_gpu.f90` | 1 | 2 |

`tendget_gpu` 與 `zx_gpu` 的 capture 內也各有 2 個 kernel，**同樣會掉出 graph**——
也就是說這兩支的結果目前也應該是錯的，只是還沒有像 `bal` 那樣明顯的指標暴露出來。

##### 建議修法（尚未實作）

**方向：ROCm 路徑不要用 graph capture，改成 eager 執行。** 但有一個一定會踩到的陷阱：

> 程式碼的結構是 `if (.not. cg_created) then ... 錄製 ... end if` ＋ 每次呼叫都 `graph_launch`。
> **如果只是把 capture/launch 改成 no-op，那段工作就只會在第一次呼叫時執行一次，之後完全不做事。**
> 必須讓 `cg_created` 永遠保持 `.false.`（即那個區塊每次都進去），eager 執行才會正確。

具體做法（在 `cmake/acc2omp.py`，維持不改 `src/nvidia/` 的原則）：
1. 對這 3 個檔案，把 `if (.not. cg_created) then` 改寫成每次都成立的條件（或直接刪掉守衛與對應 `end if`）；
2. 把 `accx_graph_launch(...)` 改成 no-op（或在 `src/rocm/openacc.f90` 讓它直接 return）；
3. `accx_async_begin_capture` / `end_capture` 同樣 no-op 化。

**驗證方式很便宜**：跑批只要跑到 `iteration= 1` 之後（約 25 分鐘）就能看 `bal` 是不是不再是 `0.`、
`pt` 有沒有跟著 CPU 版一起變。**不必等整輪 55 分鐘，更不必等積分。**
CPU 版的參考值就在 `job/TCo383L72_IC_sample_rocm_cpu.log.20260908:1564` 起。

> ⚠️ §6 第 1 層有記過「把整個 capture 機制改成 no-op、eager 執行」試過而且無效——
> **但那是 2026-09-08、為了排除當時的 `hsa_amd_memory_lock` 崩潰而做的，而且那次沒有處理 `cg_created` 的問題**，
> 跟這裡要解的問題不同。前 3 層與第 6 層都修好之後，這個方向值得重做一次。


##### 修法嘗試（2026-09-11）：診斷確認，但 eager 化還沒讓數值變對

**結論先講：第 7 層的診斷是對的——NNMI 確實原本是個完全的 no-op，改成 eager 之後它開始運作了。
但目前的實作只做到「會動」，還沒做到「算對」。這是部分進展，不是完成。**

###### 實作內容（都在 `cmake/acc2omp.py`，未動 `src/nvidia/`）

對 `nnmi_gpu` / `tendget_gpu` / `zx_gpu` 三個檔案：

1. `if (.not. cg_created) then` → **無條件**。必須這樣：`spread_event`/`pack_event`/
   `subwrk_stream`/`matmul_stream`/`handle` 都是**沒有 `SAVE` 的區域變數**，
   setup 只做一次的話第二次呼叫就是未定義值。重跑 setup 很便宜
   （`acc_get_cuda_stream` 在 shim 有快取、`cublasGetHandle` 冪等），
   `cudaEventDestroy` 保留所以 event 不會累積。
2. `accx_async_begin_capture` / `accx_async_end_capture` → no-op。
3. `accx_graph_launch` → 換成 `call geps_acc_wait(async_id)`，保留原本的同步語意。
4. 原 capture 區塊內每個 `!$omp target` 之前插入 `call geps_acc_wait(async_id)`。
5. **`zx_gpu` 的 `vars` 加上 `save`**（見下，這是踩到才發現的）。

###### 踩到的坑：`vars` 沒有 `save`（已修）

第一次驗證在 `iteration= 2` 崩在 **`0x143000`**——這種極小的位址是未初始化指標。
原因是 `zx_gpu` 的 `REAL(kind=RTYPE), dimension(:), pointer :: vars`：
原設計「配置一次 → 位址烤進 graph → 之後只 replay」，
**這個 routine 從來不需要在後續呼叫讀 `vars`**，所以它是普通區域指標無妨。
改成 eager 之後每次呼叫都要讀它 → 第二次就是垃圾位址。
加上 `save` 後這個崩潰消失。（已確認 `nnmi_gpu`/`tendget_gpu` 沒有同類的一次性配置。）

###### 驗證結果：有動，但算錯

log：`job/TCo383L72_IC_sample_rocm.log.20260911_layer7b`

| | CPU 參考 | 修法前 | 修法後 |
|---|---|---|---|
| iteration 1 的三個 `bal` | 1.76e-05 / 6.39e-06 / 1.65e-06 | **0. / 0. / 0.** | **2172.9 / 72.0 / 187.4** |
| iteration 1 後 `pt` | 0.91441 / -7.67721 | 0.92182 / -7.66227（不變） | 0.92182 / -7.66227（仍不變） |
| iteration 2 的 `bal` | 4.09e-07 … | 0. | **2.39e+25 / 7.76e+23 / 2.33e+22** |
| iteration 2 後 `pt` | 0.91410 / -7.67534 | 不變 | **`************`（格式溢位）** |

**最終結局**：iteration 3 的三個 `bal` 全部變成 `NaN`、`pt` 也是 `NaN NaN`，
積分第一步之後以 `fct model fail` 結束。**整個跑批 `memory access fault` = 0**
（第 6 層的修法與 `vars` 的 `save` 都有效），但數值從 iteration 2 起就發散了。

三件事同時成立：
- ✅ **`bal` 不再是精確 0**——kernel 現在真的在跑、也真的讀到 dgemm 的輸出，**診斷得到證實**
- ❌ **量級差約 8 個數量級**，而且 iteration 2 起直接發散
- ❌ **`pt` 在 iteration 1 後仍然沒變**（到 iteration 2 才變，然後溢位）

###### ⚠️ 關鍵線索：iteration 1 本身就是非決定性的

加 `save` 之前那次跑批，iteration 1 的 `bal` 是 **2787.0 / 88.7 / 240.0**；
加 `save` 之後是 **2172.9 / 72.0 / 187.4**。
**`save` 只影響第二次以後的呼叫，iteration 1 在兩個 binary 裡應該完全相同。**
數值卻不同 → **iteration 1 就在讀未初始化/競態的資料**。

###### 下一步最有根據的假設：`pack_event` 被重複使用

看 `nnmi_gpu`（`tendget_gpu`/`zx_gpu` 同構）的迴圈：

```fortran
do L = ...; do m = ...; do k = 1, 2
   ind = ...
   CUDACHECK(cudaStreamWaitEvent(subwrk_stream(ind), spread_event, 0))
   call dgemm(...)                                            ! 在 subwrk_stream(ind) 上
   CUDACHECK(cudaEventRecord(pack_event, subwrk_stream(ind)))  ! ← 單一 event，反覆覆寫
   CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))
end do; end do; end do
```

**`pack_event` 是一個 event，卻在迴圈裡被每個 `ind` 反覆 record。**
在原本的 capture 模式下這沒問題——每組 record/wait 會被錄成 graph 的一條依賴邊，
語意是「stream 等這一個 dgemm」。**改成 eager 之後，event 只保有最後一次 record 的狀態**，
`hipStreamSynchronize(stream)` 因此只保證等到最後一個 dgemm，前面的都可能還在飛。
這正好解釋 iteration 1 的非決定性。

**建議的修法**：eager 路徑不要依賴這組 event，改成在每個 dgemm 之後直接同步它自己的流
（`hipStreamSynchronize(subwrk_stream(ind))`），或替每個 `ind` 配一個獨立 event。
前者最簡單、也最符合「eager 就是放棄非同步」的取捨。
可以在 acc2omp 把 `cudaEventRecord(pack_event, subwrk_stream(ind))` ＋
`cudaStreamWaitEvent(stream, pack_event, 0)` 這一對改寫成一次
`cudaStreamSynchronize(subwrk_stream(ind))`。

**若修好後 `bal` 仍不對**，下一個要查的是 `pt` 為什麼在 iteration 1 後沒變——
修正量的傳遞路徑是
`nnmi_gpu → x → vartrix_gpu(-2) → zx_gpu(evectr,…) → correct_gpu(… vornow/divnow/temnow/plnow) → transr_gpu_cuda_graph → pt`，
其中 `correct_gpu`（`src/nvidia/correct_gpu.f90`）是真正把修正量加回譜係數的地方，值得先看。

###### 第二次嘗試（2026-09-11）：`pack_event` 改成直接同步 —— 無效

依上面的假設，把三個檔案裡的

```fortran
CUDACHECK(cudaEventRecord(pack_event, subwrk_stream(ind)))   ! nnmi_gpu ×2
CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))        ! tendget/zx ×1
```

改寫成 `CUDACHECK(cudaStreamSynchronize(<該 dgemm 所在的流>))`（已在 `acc2omp.py` 實作）。

**結果：非決定性完全沒有改善。** 三次跑批的 iteration 1 `bal`：

| 版本 | `bal`（三個 vertical mode） | iteration 1 後 `pt` |
|---|---|---|
| eager 化 | 2787.0 / 88.7 / 240.0 | 不變 |
| ＋ `vars` 加 `save` | 2172.9 / 72.0 / 187.4 | 不變 |
| ＋ `pack_event` 直接同步 | 5803.4 / 177.4 / 505.7 | 不變 |
| **CPU 參考** | **1.76e-05 / 6.39e-06 / 1.65e-06** | **0.91441 / -7.67721** |

**所以 `pack_event` 不是（唯一的）非決定性來源。** 這個改動本身是對的（eager 下單一 event
確實無法為迴圈排序），建議保留，但它沒解決問題。

後續與前一次跑批完全同型態：iteration 2 的 `bal` 發散到 1.26e+23、
`pt` 變成 162810.47 / 552050.03，iteration 3 全為 `NaN`，積分照樣走不下去。
**三次跑批的 `memory access fault` 都是 0**——崩潰類問題確實都修掉了，剩下的純粹是數值。
log：`job/TCo383L72_IC_sample_rocm.log.20260911_layer7c`。

###### ★ 下一步：一個能一次定位的實驗

目前最有指向性的矛盾是：

- `bal` **非零且很大**（5803），而 `bal` 是由 `nnmi_gpu` 內部的 `wrk` 算出來的
  → 第一次 dgemm（`evec^T · x` → `wrk`）**確實有在跑、有輸出**
- 但 `pt` 在 iteration 1 後**逐位元不變**
  → 傳到 `vorten` 的修正量幾乎是零

修正量的路徑是
`nnmi_gpu`（第二次 dgemm 寫 `x`）→ `vartrix_gpu(...,-2)` → `zx_gpu(evectr,…)` → `correct_gpu`。
已確認 `src/nvidia/correct_gpu.f90:134` 是 `vornow(k,i,n,m) = vornow(k,i,n,m) + vorten(k,i,n,m)`，
**增量確實有被加回去**，所以斷點在它上游。

**建議的實驗（一次 25 分鐘的跑批就能定位）**：在 `src/nvidia/initial_gpu.f90` 的 NNMI 迴圈裡，
用 §6「第 5 層根因」那個 host 端重算探針的手法，依序 dump：

1. `max|x|`　**在 `nnmi_gpu` 回傳之後**　→ 第二次 dgemm 有沒有寫出東西？
2. `max|vorten|`　**在 `vartrix_gpu(...,-2)` 之後**
3. `max|vorten|`　**在 `zx_gpu(evectr,…)` 之後**
4. `max|vornow|`　**在 `correct_gpu` 前後**

第一個變成零（或變成垃圾）的地方就是斷點。CPU 版可以用同樣的探針跑一次當對照
（CPU 跑批只要 5 分鐘）。

**另一條值得並行查的線索**：`evec` / `eval` 來自 `eigen_mx_gpu`（`src/rocm/` 的手寫替換，
走 hipSOLVER）。如果特徵分解本身就錯，`bal` 的量級差 8 個數量級就有了直接解釋。
`hip_compat.cc` 裡有 `hip_compat:dsyevj` 的 before/after logging，但**實測 0 筆**——
那段 logging 沒有觸發，要查得先確認它的條件或改用別的方式。

###### 已排除

- **精度不符**：懷疑過 `nnmi_gpu` 宣告 `real` 卻呼叫雙精度 `dgemm`，但 `cmake/compilers/ROCM.cmake:16`
  有 `-fdefault-real-8 -fdefault-double-8`，與 NVIDIA 的 `-r8`、GNU 的 `-fdefault-real-8` 等價。**不是這個原因。**


#### 第 8 層（2026-09-11）：~~`eigen_mx_gpu` 的特徵向量矩陣是垃圾~~ ⚠️ **這個結論是錯的，`max|mx| > 1e12` 是量測假象——推翻過程見本節末，保留是為了不讓下一個人重踩**

**這是數值錯誤的真正源頭。一次 43 分鐘的探針實驗就定位到了。**

##### 實驗設計

在 `cmake/acc2omp.py` 對 `initial_gpu` 注入一條完整的探針鏈，沿著修正量的傳遞路徑
逐段量 `max|.|`（host 端 `!$omp target update from` 之後 `maxval(abs(...))`），
一次跑批（跑到 `correct_gpu` 就收批，約 18 分鐘重編 ＋ 18 分鐘跑批）涵蓋全部：

| locid | 位置 | 量測 |
|---|---|---|
| 400 | `nnmi_gpu` 之前 | `x`(輸入)、`mx`(evec)、`eval` |
| 401 | `nnmi_gpu` 之後 | `x`(輸出) |
| 402 | `vartrix_gpu(-2)` 之後 | `vorten` / `divten` / `phiten` |
| 403 | `zx_gpu(evectr)` 之後 | 同上 |
| 404 | `correct_gpu` 之後 | `vornow` / `divnow` / `temnow` |

##### 結果

```
nnmi 之前      x(輸入) = 0.00071    mx(evec) = >1e12  ⚠️   eval = 256.273
nnmi 之後      x(輸出) = 4e-06
vartrix(-2)後  vorten = 0           divten = 0            phiten = 7.1e-05
zx(evectr)後   vorten = 0           divten = 6e-06        phiten = 0.002442
correct 之後   vornow = 2.1e-05     divnow = 6e-06        temnow = 3333.62  ⚠️
```

**`mx` 是特徵向量矩陣，元素應該是 O(1)**（對稱矩陣的特徵向量正交歸一化，`|元素| ≤ 1`），
**實際超過 1e12**。這一口氣解釋了先前所有看不懂的現象：

- `bal = Σ wrk²`，而 `wrk = evec^T · x`。`evec ~1e12` → `bal` 差 8 個數量級 ✅
- 每次跑批的垃圾值不同 → **iteration 1 就非決定性** ✅
- `temnow` 被污染成 3333.62（溫度場單位是 K，正常量級 ~200-300）✅
- 修正量傳下去是垃圾，積分因此 `NaN` ✅

##### 已排除的上游

- **不是 sequence association**：`initial_gpu.f90:124` 宣告 `mx(no*no, 2*jtmax*nnmivm)`，
  `src/rocm/eigen_mx_gpu.f90:11` 的 dummy 完全同形，不是第 1 層那種 reshape 問題。
- **不是 RCCL 廣播**：`mx` 進 `eigen_mx_gpu` 前會經過
  `mpe2d_row_broadcast_gpu`（`src/rocm/` 手寫檔，**不經 acc2omp，所以第 6 層的修法對它無效**），
  但 `src/rocm/mpe2d_row_broadcast_gpu.f90:20` 是 `if (nsize <= 1 .or. n <= 0) return`，
  單 rank smoke 下它直接返回，是 no-op。
  （⚠️ **多 rank 時這支仍有第 6 層同類的風險**，`ncclBroadcast` 之後沒有任何等待，
  將來跑多 rank 要回頭檢查。）

##### 嫌疑：`geps_rocm_dsyevj_dev` 的 hipSOLVER 呼叫（`src/rocm/hip_compat.cc:488`）

```cpp
int geps_rocm_dsyevj_dev(double *A, int n, double *W, int *devinfo) {
  ...
  hipsolverDnXsyevjSetTolerance(g_syevj, 1e-15);
  hipsolverDnXsyevjSetSortEig(g_syevj, 0);          // 不排序
  ...
  rc = hipsolverDnDsyevj(h, HIPSOLVER_EIG_MODE_VECTOR, HIPBLAS_FILL_MODE_UPPER,
                         n, A, n, W, work, lwork, devinfo, params);
  hipDeviceSynchronize();
  hipFree(work);
  return rc;                                        // devinfo 從來沒被檢查
}
```

**最可疑的兩點**：
1. **`devinfo` 完全沒有被檢查。** `hipsolverDnDsyevj` 用 `devinfo` 回報未收斂的數量，
   而 `rc` 在這種情況下仍可能是 0。`eigen_mx_gpu.f90:41` 只檢查 `ierr`（即 `rc`），
   所以 Jacobi 不收斂會**完全靜默**。
   `hipsolverDnXsyevjSetTolerance(1e-15)` 搭配預設 max sweeps 是很容易不收斂的設定。
2. **`hipsolverDnXsyevjSetSortEig(g_syevj, 0)`（不排序）**——要確認 NVIDIA 版
   （`cusolverDnXsyevjSetSortEig`）是否也設 0。若兩邊排序行為不同，
   `eval`/`evec` 的對應順序就會錯，而 `nnmi_gpu` 是用 `eval(i)` 去除 `wrk(i)` 的
   （`swap = -w2_t/e_t`），順序錯會直接產生巨大的值。

> 註：Jacobi 旋轉本身保持正交性，**即使不收斂也不該得到 1e12 的元素**。
> 所以更可能是「輸出根本沒被寫入、讀到未初始化的裝置記憶體」，
> 或是 `A`/`lda` 的傳遞有問題。這也與「每次跑批數值都不同」吻合。

##### ✅ 已排除：hipSOLVER `Dsyevj` 本身沒有問題（2026-09-11，最小 repro）

寫了一支獨立的小程式（`scratchpad/eigrepro2.cc`，把 `geps_rocm_dsyevj_dev` 的呼叫序列原樣內嵌，
避免拖進整個 compat 層的相依），拿一個 4×4 對稱矩陣測試，**幾秒鐘就跑完**：

```
bufferSize rc=0 lwork=16
rc=0  devinfo=0
eigenvalues: 0.254719 1.822717 3.177283 4.745281      ← 正確
max|eigenvector entry| = 0.777951   (must be <= 1)    ← 正確
column 0..3 norm^2 = 1.000000                         ← 全部正交歸一
```

**hipSOLVER 的 `Dsyevj`、`SetTolerance(1e-15)`、`SetSortEig(0)`、buffer 配置、`devinfo`
全部正常運作。** 所以上面「Jacobi 不收斂」「排序不一致」兩個嫌疑都**不成立**。

##### ★★ 修正後的假設：解算結果沒有寫進我們讀的那塊記憶體

`mx` 是 **`intent(inout)`**：
- **輸入** = `coftrix_refactor`（`initial_gpu.f90:214`）在主機算出的係數矩陣，
  **可以合法地很大**（它不是正交矩陣）
- **輸出** = 特徵向量矩陣，**必須 ≤ 1**

既然求解器本身是好的，`max|mx| > 1e12` 最可能代表 **`mx` 根本沒被寫回，我們讀到的還是輸入值**。
`eigen_mx_gpu.f90:31` 是

```fortran
!$omp target data use_device_addr(mx, eval, devinfo)
   ierr = geps_rocm_dsyevj_dev(c_loc(mx(1, j)), nn, c_loc(eval(1, j)), c_loc(devinfo(j)))
```

**`use_device_addr` 正是第 1 層踩過的坑**（amdflang 對它的處理很脆弱）。
外加 `hip_compat.cc:44` 的 `device_ptr()` 會在 `omp_target_is_present` 為假時
**直接把指標原樣傳回並假設它已經是裝置位址**——如果 `use_device_addr` 沒生效，
傳進去的是**主機位址**，hipSOLVER 就會寫到主機記憶體（或直接無效），
而裝置上的 `mx` 原封不動。這與觀察到的現象完全一致。

##### ❌ 推翻（2026-09-11 稍晚）：`max|mx| > 1e12` 是我的探針量錯了

加了兩個診斷之後，上面的結論站不住：

**(a) 指標是對的，`use_device_addr` 有生效。** 在 `geps_rocm_dsyevj_dev` 印出收到的指標：

```
[dsyevj] n=2 A=0x892fd800000 is_present=0 mapped=(nil) resolved=0x892fd800000
[dsyevj] n=1 A=0x892fda9eb20 is_present=0 mapped=(nil) resolved=0x892fda9eb20
[dsyevj] n=3 A=0x892fdd3d640 ...
```

`0x892f...` 落在本專案其他裝置配置的位址範圍內（對照：`0x1065...`、`0xa890...`）。
**`is_present=0` 對「裝置位址」而言是正常的**——`omp_target_is_present()` 接受的是主機位址，
拿裝置位址去問當然回答否；`device_ptr()` 原樣傳回是正確行為。
所以「`device_ptr()` 把主機位址靜默傳過去」這個假設**不成立**。

**(b) `max|mx|` 在 eigen 前後都是 `>1e12`，而 `max|eval|` 從 0 變成 256.273。**
求解器確實有執行、也確實寫進了裝置記憶體（`eval` 被寫了）。
既然 hipSOLVER 已用最小 repro 證實正確（見上），`mx` 的最大值前後不變只說明一件事：
**我量的不是特徵向量。**

**為什麼**——CPU 與 GPU 的資料結構根本不同：

| | CPU `initial.f90:29` | GPU `initial_gpu.f90:124` |
|---|---|---|
| `mx` | `mx(no*no)`，**單一矩陣**，每個 `(mf,L)` 迴圈內重複使用、算完立刻消費 | `mx(no*no, 2*jtmax*nnmivm)`，**一次擺下全部**，由 `eigen_mx_gpu` 批次求解 |

GPU 版每個 column `j` 放一個 `nn×nn` 矩陣（`nn = nnlist(j)`，**每個 j 都不一樣**），
但 column 的長度固定是 `no*no`。**每個 column 超出 `nn*nn` 的尾端從來沒有被寫過**
（`mx` 是區域自動陣列，`enter data create` 在裝置端也不清零）。
對整個陣列取 `maxval(abs(mx))` 量到的就是那些未初始化的尾巴。

**功能上這些尾端不會被讀**：hipSOLVER 用 `lda=nn`、`nnmi_gpu` 也只取 `evec(1,ind)` 的前 `nn*nn`。
所以**這不是 bug，是探針設計錯誤**。（它同時也解釋了為什麼這個數字每次跑批都不同。）

##### 目前還站得住的事實

- `bal` 比 CPU 大 8 個數量級，且**跑批之間非決定性** ← 仍未解釋
- `pt` 在 iteration 1 後逐位元不變 ← 仍未解釋
- `nnmi_gpu` 輸出的 `max|x| = 4e-06`（輸入是 0.00071）
- `eval = 256.273`（是否合理未經驗證，CPU 版沒有對應的輸出可比）
- hipSOLVER `Dsyevj` 正確；`use_device_addr` 有生效；指標正確

##### ★ 下一步：把探針改對

**必須只量實際使用的區域**，也就是對幾個有效的 `j`，量 `mx(1:nn*nn, j)`
（`nn = nnlist(j)`）而不是整個陣列。同樣的量測要在 `eigen_mx_gpu` 前後各做一次：

- **前**：係數矩陣，數值可以大
- **後**：特徵向量，**必須 ≤ 1，且每個 column 的平方和必須 = 1**

直接在 host 端檢查正交歸一性（`sum(V(:,c)**2)`）比看最大值更可靠。
若這項通過，特徵分解就徹底清白，下一個要查的是 `nnlist` / `eval` 的排序
（`hipsolverDnXsyevjSetSortEig(0)` 不排序，要確認 NVIDIA 的 `cusolverDnXsyevjSetSortEig`
是否也設 0；`nnmi_gpu` 是用 `eval(i)` 去除 `wrk(i)`，順序錯會直接產生巨大的值）。

##### ✅ 特徵分解徹底清白（2026-09-11，正確的探針）

把探針改成**只掃實際使用的 `nn*nn` 區域**，並直接驗證正交歸一性
（`src/rocm/eigen_mx_gpu.f90` 末端，`!$omp target update from(mx, eval)` 之後在 host 端算）：

```
取樣 64 個矩陣                                   ← ⚠️ 這個取樣是錯的，見下
  max|entry|          = 1.000000      (必須 <= 1)
  max|colnorm^2 - 1|  = 0.000e+00     (必須 ~ 0)
```

> ⚠️ **這次量測不算數（2026-09-11 稍後發現）**：`dbg_n >= 64` 的上限配上迴圈順序
> `(L, m, k)`，只掃到**小矩陣**（ns = 48、50、144、146…），而 `bal` 的質量所在的
> ns 到 530 / nn 到 576 **從來沒有被檢查過**。這是同一天第五次踩「取樣沒先驗證代表性」。
> 結論本身後來被無上限的重跑證實（見下），但當時的證據不足以支持它。

##### ✅ 特徵分解徹底清白（2026-09-11 重跑，**無取樣上限**）

拿掉 `dbg_n >= 64`，改掃**全部矩陣**，並額外記錄「最差的那個矩陣的 `nn`」：

```
特徵向量正交歸一性（全部矩陣，無取樣上限）
  max|entry|          = 1.000000        (必須 <= 1)   ✅
  max|colnorm^2 - 1|  = 0.000000e+00    (必須 ~ 0)    ✅
  最差的矩陣 nn       = 576
```

**代表性檢查（新規矩 §7.1）**：`worstnn = 576` 表示「最差 column norm」的紀錄保持者
是 576 階矩陣 ⇒ **大矩陣確實被掃到了**，這次的樣本涵蓋了 `bal` 的質量所在。

**兩個「太漂亮」的數字都有正當解釋，不是單位矩陣：**

| 數字 | 看起來的問題 | 實際原因 |
|---|---|---|
| `max\|colnorm²−1\|` 恰好 0 | 雙精度不可能精確為 0 | **編碼地板**：探針傳的是 `p1 = nint(worstnorm*1e9)`，解析度只到 1e-9，真值 ~1e-15 必然量化成 0。真正的結論是 **< 5e-10** |
| `max\|entry\|` 恰好 1.0 | 隨機正交矩陣不該剛好碰到 1 | 矩陣 `nn` 範圍是 **1..576**，**1×1 矩陣的正規化特徵向量本來就恰好是 ±1** |

> 驗證發出點唯一性：`locid 500` 全專案只有兩個發出點——
> `src/rocm/eigen_mx_gpu.f90:111`（hyp=6 → `E`）與
> `src/rocm/cyclic_cell_ppm_gpu.f90:57`（hyp=1 → `A`）。
> 本次跑批的 `vram:500` 共 2 列，`E` 只有 1 列（`p0=1000000, p1=0, p2=576`），
> 確定來自 eigen 探針。

**`eigen_mx_gpu` 產生的特徵向量在所有尺寸（含 nn=576）都完全正交歸一。**
加上先前的最小 repro（hipSOLVER 正確）與指標診斷（`use_device_addr` 有生效），
**整條特徵分解路徑可以結案了**。

> ⚠️ **踩到的坑**：我用了 `locid 500`，但 `src/rocm/cyclic_cell_ppm_gpu.f90:57` 已經在用
> 同一個編號。兩者會混在一起，**必須靠 `hypothesisId` 分辨**（我的是 `E`/hyp=6，
> 既有的是 `A`/hyp=1）。新增探針前先 grep 一下 locid 有沒有被佔用。

##### ⚠️ 一個我自己提出又自己推翻的「算術矛盾」（保留以免重踩）

一度主張：`|evec| ≤ 1`、`|x| ≤ 7.1e-4`、`nn` 只有 1~3 ⇒ `bal` 應該 ~1e-5，
所以 `bal = 5803` 在算術上不可能，一定是 `wrk` 讀到別的緩衝區。

**這個推論是錯的，錯在 `nn` 的估計。** `[dsyevj]` 診斷我設了 `seen < 4`，
**只印了前 4 次呼叫**（n = 2, 1, 3, 3），我從這 4 筆樣本外推成「`nn` 只有 1~3」。

實際上 `initial_gpu.f90:187-191`：

```fortran
nbig = jtrun - mf + 1                              ! jtrun = 384
ns = 2*int((nbig + 1)/2) + int(nbig/2)             ! mf=1 → nbig=384 → ns=576
na = 2*int(nbig/2) + int((nbig + 1)/2)
```

**`nn` 的範圍是 1 ~ 576**（前 4 次呼叫剛好都是 `mf` 很大、矩陣很小的那幾個）。
重算上界：`|wrk| ≤ nn·|evec|·|x| ≈ 576 × 7.1e-4 ≈ 0.41`，
`bal_t ≤ 2·nn·0.41² ≈ 194`，對 384 個 `mf` 求和可到 ~7e4。
**`bal = 5803` 完全落在合理範圍內，沒有矛盾。**

> **教訓**：用診斷樣本外推母體特性之前，先確認樣本有沒有代表性。
> 這個 `seen < 4` 的上限是我自己設的，卻忘了它會讓樣本偏向迴圈前段。

##### 目前確實還站得住的事實

- **特徵分解路徑完全清白**（三個獨立證據，見上）
- `bal` 比 CPU 大 8 個數量級（5803 vs 1.76e-05），且**跑批之間非決定性**
- `pt` 在 iteration 1 後逐位元不變
- `max|x|` 進 `nnmi_gpu` 前 = 7.1e-4、出來 = 4e-06
- **還不知道 CPU 版在同一點的 `x` 是多少**——沒有對照，就無法判斷 7.1e-4 是對是錯

##### ★ 下一步（建議）

**先建立 CPU 對照，再談 GPU 哪裡錯。** 目前所有 GPU 數字都缺乏基準：
`x`、`wrk`、`evec` 的量級在 CPU 版分別是多少？沒有這個，任何「GPU 的值不合理」
都只是猜測（今天已經因此誤判三次）。

具體作法：在 `src/initial.f90`（CPU 路徑）的對應位置加同樣的 `max|.|` 輸出，
跑一次 CPU smoke（**只要 5 分鐘**），得到 `x` / `wrk` / `eval` 的基準量級。
然後與 GPU 的探針逐項比對，**第一個量級對不上的地方就是真正的斷點**。

這比繼續在 GPU 側猜測有效率得多，而且 CPU 跑批便宜到可以反覆做。

##### 教訓##### 教訓

**對「一次擺下多個不同大小子矩陣」的批次陣列，`maxval(abs(整個陣列))` 沒有意義**——
未使用的填充區會主導結果，而且因為未初始化，每次跑批都不一樣。
這跟 §6「問題 B」那次的誤判是同一類錯誤：**指標選得不對，就會得到看似有力、實則無效的證據**。


#### ★★★ 第 9 層（2026-09-11 定位）：acc2omp 把內層 `!$acc loop` 的 `private(...)` 丟掉 → 資料競態

**這是數值錯誤的根因，而且是 `cmake/acc2omp.py` 的系統性缺陷，不是單點 bug。**

##### 先建立 CPU 基準（這一步是關鍵，早該做）

在 `src/initial.f90` 的對應位置加上同樣的 `max|.|` 輸出，跑一次 CPU smoke（**5 分鐘**），
得到 GPU 探針的對照基準：

| 量 | CPU（iteration 1） | GPU | 判讀 |
|---|---|---|---|
| `max\|evec\|` | 0.9999999492 | 1.000000 | ✅ 一致 |
| `max\|x\|` 進 `nnmi` | 9.82e-05 | 7.1e-4 | 同量級（GPU 是批次最大值） |
| `max\|x\|` 出 `nnmi` | 7.97e-07 ~ 1.44e-06 | 4e-06 | 同量級 |
| **`bal`** | **1.76e-05** | **5803** | ❌ **差 3×10⁸** |

**`x` 與 `evec` 兩邊都對得上，唯獨 `bal` 離譜。**
而且 `x_out ≈ x_in / eval`（7.1e-4 / 256 ≈ 2.8e-6，實測 4e-06）
→ **`wrk` 的 `1/e_t` 變換確實有正確執行**，dgemm 也沒問題。
問題被鎖死在「算 `bal` 的那個 kernel」。

##### 根因

`src/nvidia/nnmi_gpu.f90:98-105`（NVIDIA 原始碼，**正確**）：

```fortran
!$acc parallel loop collapse(2) private(bal_t, mf) async(async_id)
do L = 1, nnmivm
   do m = 1, mlistnum
      mf = mlist(m)
      bal_t = 0.
      !$acc loop vector collapse(2) reduction(+:bal_t) &
      !$acc&           private(ind, nn, e_t, w1_t, w2_t, swap)     ← 注意這個 private
      do k = 1, 2
         do i = 1, no
            ...
```

`build_rocm_hip/src/hip_src/nvidia_nnmi_gpu.f90`（acc2omp 產生，**有 bug**）：

```fortran
!$omp target teams distribute parallel do collapse(2) private(bal_t, mf)
do L = 1, nnmivm
   do m = 1, mlistnum
      mf = mlist(m)
      bal_t = 0.
      ! acc2omp: nested loop (serial on parent thread)     ← private(...) 整個不見了
      do k = 1, 2
         do i = 1, no
            ...
```

內層迴圈被攤平成「serial on parent thread」——**reduction 變序列累加在語意上沒問題**，
但 **`private(ind, nn, e_t, w1_t, w2_t, swap)` 被整個丟棄**。
外層的 `private` 只有 `bal_t, mf`，於是
**`ind`、`nn`、`e_t`、`w1_t`、`w2_t`、`swap`（以及序列 DO 的 `i`、`k`）
在整個 `parallel do` 裡變成共享變數，每個執行緒都在寫同一份**。

這一口氣解釋了全部：
- `bal` 是垃圾值（每個執行緒讀到別人寫的 `w1_t`/`w2_t`）✅
- **跑批之間非決定性**（競態的典型特徵）✅
- `wrk` 的變換「看起來大致對」——因為 `bal_t` 本身是 private、
  而且除以 `e_t` 之後量級被平均掉了 ✅

##### 影響範圍：系統性，不是單點

```
原始碼中帶 private(...) 的內層 !$acc loop：79 處
產生檔中 "nested loop (serial on parent thread)"：
  cwa_radlw_main_gpu.f90          52
  cwa_radiation_aerosols_gpu.f90  33
  cwa_adjptqintp_gpu.f90          15
  cwa_gwdc_gpu.f90                14
  cwa_nor_gwdp_gpu.f90            13
  cwa_ozphys_2015_gpu.f90         13
  cwa_mod_stochastic_physics_gpu.f90 11
  ... 等等
```

**不是每一處都會造成競態**（要看該內層迴圈是否落在平行區域內、以及那些變數是否真的被多執行緒寫），
但只要落在 `parallel do` 裡面就會。**這解釋了為什麼 GPU 路徑到處都是「跑得起來但數值不對」。**

##### ★ 修法方向

在 `cmake/acc2omp.py` 把內層 `!$acc loop` 攤平成序列迴圈時，
**把它的 `private(...)` 變數併入外層 OpenMP 指令的 `private(...)` 子句**，
而不是丟掉。序列 DO 的迴圈變數（此例的 `i`、`k`）同樣要併進去。

- 這是**翻譯器層級的修正**，符合 §7「編譯器問題在 toolchain 解」的規矩。
- 修完之後 `nnmi_gpu` 的外層應該變成
  `private(bal_t, mf, ind, nn, e_t, w1_t, w2_t, swap, i, k)`。
- **驗證很便宜**：跑到 iteration 1 的第一個 `bal=`（約 20 分鐘），
  看它是否掉到 **1.76e-05** 的量級（CPU 基準就在
  `job/TCo383L72_IC_sample_rocm_cpu.log.20260911_base`）。
- 修好之後**整個 `src/cwa/` 的物理過程也會一起受惠**——那些檔案有大量同樣的模式。

##### ⚠️ 修法已完成並通過編譯，但**沒有修好 `bal`** —— 競態假設被推翻

**修法本身是對的、也已經 landed**：`cmake/acc2omp.py` 新增 `_privatise_nested_loops()`，
在 `_merge_parallel_loop(_join_continuations(lines))` 之後對 joined 清單再跑一次，
把內層 `!$acc loop` 的 `private(...)`（含被攤平的序列 DO 的迴圈變數）併入外層指令。
**完整重編通過（rc=0），911 個指令被補上 private 子句。**

`nnmi_gpu` 產出的正是預期的樣子（`build_rocm_hip/src/hip_src/nvidia_nnmi_gpu.f90:97`）：

```fortran
!$omp target teams distribute parallel do collapse(2) &
!$omp&  private(bal_t, mf) private(ind, nn, e_t, w1_t, w2_t, swap, k, i)
```

**但 `bal` 完全沒變**：

| 版本 | iteration 1 的前兩個 `bal` |
|---|---|
| 第 7 層 eager 化 | 2787.0 / 88.7 |
| ＋ `vars` save | 2172.9 / 72.0 |
| ＋ `pack_event` 直接同步 | 5803.4 / 177.4 |
| **＋ 第 9 層 private 修法** | **2624.6 / 84.7** |
| **CPU 基準** | **1.76e-05 / 6.39e-06** |

**所以「內層 private 被丟掉造成競態」不是 `bal` 錯誤的原因。**
這個修法仍然是**語意上正確的**（那些變數本來就該是 private），
建議保留，但它解決的不是這個問題。

###### ❌ 已推翻：`wrk` 的緩衝區身分沒有問題（2026-09-11，直接量測）

曾主張 dgemm 寫的緩衝區與 OpenMP kernel 讀的不是同一塊。**直接比對位址，證明是同一塊。**

作法（照 `[dsyevj]` 診斷的模式）：
- `src/rocm/hip_compat.cc` 新增 `geps_dbg_mapped_`，用 Fortran 傳參考慣例收主機位址，
  印出 `omp_target_is_present` 與 `omp_get_mapped_ptr`（**不需要 `iso_c_binding`**）；
- `geps_blas_dgemm` 印出實際交給 hipBLAS 的 `A`/`B`/`C` 指標；
- 經 acc2omp 注入，在 `initial_gpu` 與 `nnmi_gpu` 內各呼叫一次 `geps_dbg_mapped`。

```
[mapped] tag=1 host=0x7ffac48b9850 is_present=1 device=0xc7e2ee00000   (initial_gpu 的 nnmi_buf)
[mapped] tag=2 host=0x7ffac48b9850 is_present=1 device=0xc7e2ee00000   (nnmi_gpu 內的 wrk)
[dgemm]  m=2 n=2 k=2  C=0xc7e2ee00000      ← 與上面完全相同
[dgemm]  m=1 n=2 k=1  C=0xc7e2ee024a0      ← base + 偏移
[dgemm]  m=3 n=2 k=3  C=0xc7e2ee04940
```

`ind=1` 的 `C` 指標就是 OpenMP 映射 `wrk` 的位址，後續呼叫是 base＋偏移，
完全符合 `wrk(1,1,ind)` 的排列。**`use_device_addr` 與隱式映射解析到同一塊記憶體。**

> ⚠️ **第一次量測是無效的，差點又誤判**：`[dgemm]` 的過濾條件原本只寫 `seen < 6`，
> 抓到的是程式中最早的 6 次 dgemm（`m=144, n=1..6, k≈394`），發生在 `nnmi_gpu` 之前、來自別處。
> `nnmi_gpu` 的形狀是 `dgemm(op, op, nn, 2, nn, ...)`（**n 固定為 2、m==k**），
> 必須這樣過濾才比得對。
> **這與稍早「用前 4 筆 `[dsyevj]` 外推 `nn` 範圍」是同一類錯誤——
> 診斷輸出設了取樣上限，卻沒確認取到的是不是要的那些。**

###### 目前還站得住的事實（2026-09-11 結束時）

- 特徵向量正交歸一 ✅（**已用無取樣上限的重跑證實，涵蓋到 nn=576**）、hipSOLVER 正確 ✅、指標傳遞正確 ✅、`wrk` 緩衝區身分正確 ✅
- 內層 `private` 已補上（911 個指令）✅ **但 `bal` 沒變**
- ~~`max|x|` 進/出 `nnmi_gpu` 與 CPU 同量級~~ ❌ **錯**：GPU `max|x|_in = 2.483`，CPU 全域 `9.82e-5`，差 2.5 萬倍
- **`bal` 比 CPU 大 8 個數量級，且跑批之間非決定性**（2717 / 2624 / 5831 / 2787 / 2173…）

###### ❌ 已推翻：`nnlist` 的裝置副本是正確的（2026-09-11）

在 `initial_gpu` 呼叫 `nnmi_gpu` 之前，先記主機值、再 `!$omp target update from(nnlist)` 取裝置值：

```
主機(update 前)        nnlist(1)=2  nnlist(2)=1
裝置(update from 後)   nnlist(1)=2  nnlist(2)=1
```

**完全相同。** 靜態上索引範圍也早已確認一致
（寫入 `j = 1+(m-1)*2+(L-1)*2*jtmax` 的 `j`/`j+1`；讀取 `ind = k+(m-1)*2+(L-1)*jtmax*2`、k=1,2）。

> 附帶發現：`nn=2` 與 `nn=1` **都不滿足 bal kernel 的守衛 `nn > 2`**，
> 所以 `ind=1,2` 對 `bal` 的貢獻本來就是零。查這一層時要注意，
> 不能只看前幾個 `ind` 就下結論。

###### 2026-09-11 結束時：已排除清單（全部經直接量測）

| 假設 | 結果 |
|---|---|
| hipSOLVER `Dsyevj` 有問題 | ❌ 最小 repro 證明正確（`qa/rocm_repro/hipsolver_dsyevj_check.cc`） |
| `use_device_addr` 沒給到裝置位址 | ❌ `[dsyevj]` 診斷證明位址正確 |
| 特徵向量矩陣是垃圾（`max\|mx\|>1e12`） | ❌ 量測假象（未初始化填充區）；**全部矩陣（nn 到 576）實測正交歸一** |
| 內層 `private` 被丟掉造成競態 | ❌ 修法已 landed、911 指令受惠，`bal` 完全沒變 |
| dgemm 與 kernel 用不同緩衝區 | ❌ 位址逐一比對，完全相同 |
| `nnlist` 裝置副本錯誤 | ⚠️ **撤回**：只比對了 2 個元素；且「主機==裝置」驗證的是傳輸、不是內容 |
| 算 `bal` 的 kernel 壞掉 | ❌ `bal`=7827.69 ≈ 輸入的 `Σx²`=8758（誤差 11%），忠實反映輸入 |
| `nnmi_gpu` 本身壞掉 | ❌ `Σx²` 8758 → 0.3921，行為合理 |

**仍然成立的事實**：`bal` 比 CPU 大 8 個數量級、**跑批之間非決定性**、`pt` 不變、積分 `NaN`；
而 `evec` 正交歸一、`x` 進出的量級與 CPU 相當、緩衝區與 `nnlist` 都正確。

###### ✅ CPU 基準做細了（2026-09-11）：`dbal ≈ Σx²`

在 `src/initial.f90` 的 `call nnmi(...)` 前後加 DBGDIST 探針（記 `bal` 的增量、
該次的 `max|eval|` 與 `Σ x²`），CPU 跑批 5 分鐘、rc=0：

```
DBGDIST mf=353 ns= 48  maxeval=255.60509448842006  sumx2=3.918820019895011E-09  dbal=3.4691620522791227E-09
DBGDIST mf=352 ns= 50  maxeval=256.2730304312415   sumx2=4.451082094292787E-09  dbal=4.182790364293406E-09
DBGDIST mf=288 ns=146  maxeval=256.27313883800167  sumx2=1.723874769439609E-08  dbal=1.5837489268097616E-08
DBGDIST mf= 96 ns=434  maxeval=256.2728108321279   sumx2=4.335448859350668E-08  dbal=3.240303087678746E-08
DBGDIST mf= 32 ns=530  maxeval=256.2724832243963   sumx2=9.96895955495776E-08   dbal=3.7273196999687265E-08
```

**兩個結論：**

1. **`maxeval` = 256.273，與 GPU 完全吻合** ⇒ `eval` 是對的（這是先前缺的對照）。
2. **`dbal ≈ sumx2`**——因為 `evec` 是正交矩陣，正交變換保範數。
   這給了一條**可以拿來檢驗 GPU 的硬關係**，不再只是比極值。

###### ‼️ 收緊後的矛盾：四個前提不可能同時成立

把上面那條關係套到 GPU：

| # | 前提 | 來源 | 狀態 |
|---|---|---|---|
| 1 | `bal ≈ Σ‖x‖²`（evec 正交 ⇒ 保範數） | CPU DBGDIST，純算術 | 站得住 |
| 2 | GPU 的 `evec` 正交歸一（含 nn=576） | 無取樣上限重跑 | ✅ **本次證實** |
| 3 | GPU `max\|x\| = 7.1e-4`，約 1.15e6 個元素 | 探針 400/401 | ❌ **後來證實是假的，見下** |
| 4 | GPU `bal = 2603`（跑批間 2173…5831 非決定性） | 跑批輸出 | 站得住 |

由 3 推得上界 `Σ‖x‖² ≤ 1.15e6 × (7.1e-4)² = 0.58`，
而 4 給的是 2603 —— **超出上界 4500 倍**。

先前唯一的逃生出口是「大矩陣的 `evec` 不正交」，**本次重跑已經把它堵死**。
所以四個前提裡必有一個錯，而**唯一沒有依新規矩（§7.1）複驗過代表性的是第 3 個**。

> ## ❌ 撤回：這個「4500 倍矛盾」不存在
>
> 複驗之後，**錯的就是第 3 個前提**：`max|x| = 7.1e-4` 這個數字在實際跑批裡
> **量不出來**。同一次跑批的探針 400 給的是 **`max|x| = 2.483`**，差 3500 倍。
> 那個 7.1e-4 是我從別處帶進來、沒有複驗就寫成前提的。
> **矛盾從頭到尾是我自己造出來的**，`evec`、`bal`、`nnmi_gpu` 都沒有問題。
>
> 這是同一天第六次同型錯誤，但形態不同：前五次是**取樣不具代表性**，
> 這一次是**把一個沒有出處的數字當成已量測的事實**。§7.1 因此補一條：
> **寫進前提表的每個數字，都要能指回「哪一次跑批、哪一個 locid」**。

###### ✅ 複驗完成（2026-09-11）：探針沒問題，前提才是錯的

先依 §7.1 把 400/401 本身查一遍，四題全過：

| 檢查 | 結果 |
|---|---|
| `x` 在裝置上？ | ✅ 生成碼 L166 `map(alloc:… x …)`，L720 才 delete，探針在資料區間內，`update from` 有效 |
| 量的是餵給 `bal` 的那個 `x`？ | ✅ 400 插在 `call nnmi_gpu` **之前** ＝ 輸入 `x`，與 CPU `sum(x(1:ns*2)**2)`（`src/initial.f90:230`）同一取樣點 |
| rank 涵蓋？ | ✅ 單 rank 跑批，`myrank` 有記錄 |
| 編碼邊界？ | ✅ 距 `-1`（NaN）/ `-2`（>1e12）都很遠 |

`maxval(abs(x))` 掃整個宣告範圍雖然和 `max|mx|` 是同一種寫法，但**對上界論證是保守的**
（掃到 padding 只會讓 max 更大、上界更鬆），所以這一點不構成問題。

**探針是乾淨的——錯的是我拿來當前提的那個數字。** 見上面的撤回框。

###### ✅ 決定性量測：`bal` 是好的，壞的是餵進去的 `x`（2026-09-11）

改印 CPU 基準的同一個量 `Σx²`（locid 407 = 進 `nnmi_gpu` 前、408 = 出來後），
而且**一次印兩個和**：整個宣告範圍 `dbgsa` 與只含使用區域 `dbgsu`
（`x(1:nn, 1:2, ind)`，`nn = nnlist(ind)`）。沒有第二個和，數字一大就分不清
是真值還是未初始化 padding——`max|mx|` 就是這樣誤判的。
編碼改用 `log10(·)*1e6`，因為要跨 1e-9 到 1e3 共 12 個數量級
（**先前 `nint(·*1e9)` 的地板正是讓我誤判「恰好為 0」的原因**）。

```
locid 407（進 nnmi_gpu 前）   Sigma x^2 = 8758      max|x| = 2.483
bal（第一次）                              = 7827.69
locid 408（出 nnmi_gpu 後）   Sigma x^2 = 0.3921    max|x| = 0.0244
使用區域元素數 = 1334076；整個宣告範圍的和 == 使用區域的和（padding 是零）
```

**兩個探針交叉檢驗自洽**（同一次跑批、同一陣列、同一時間點）：

- 由 400 的 `max` 推上界：`Σx² ≤ N·max² = 8.2e6` ⊇ 實測 8758 ✅
- 由 407 的 `sum` 推下界：`max ≥ √(S/N) = 0.081` ≤ 實測 2.483 ✅

**結論：`bal`（7827.69）≈ 它輸入的 `Σ‖x‖²`（8758），誤差 11%。**
**算 `bal` 的 kernel 是對的，`nnmi_gpu` 也是對的**——它忠實反映了輸入。
`x` 經過 `nnmi_gpu` 後從 8758 降到 0.3921，行為也合理。

**問題整個移到 `nnmi_gpu` 上游。** `x` 的產生鏈（`src/nvidia/initial_gpu.f90`）：

```
L254  tendget_gpu(phiten, zx_buf, ...)
L293  zx_gpu(evecin, vorten, divten, phiten, ...)
L298  vartrix_gpu(vorten, divten, phiten, levp, x, no, ...)   ← x 在這裡產生
L307  mpe2d_row_broadcast_gpu(x(1,1,j), ...)                  （單 rank，等同 no-op）
L313  nnmi_gpu(x, ...)                                        ← 407 在這裡量
```

CPU 對照（`job/TCo383L72_IC_sample_rocm_cpu.log.20260911_dist`）：

```
DBGCPU maxx_in = 9.821561488206484E-05     ← 全域極值，涵蓋所有 mf/L/ic
GPU   max|x|   = 2.483                     ← 大 2.5 萬倍
```

而且 CPU 的 per-mf `sumx2 ~1e-8` 與這個全域極值**量級自洽**
（`ns*2 ≈ 100..1060` 個元素，RMS ≈ 3e-6，max ~1e-4），
所以 DBGDIST 那 12 筆 `L==1, ic==1` 的取樣**在量級上確實有代表性**——這點依 §7.1 驗過了。

###### ✅ CPU 全域基準（2026-09-11）：`bal` 差 8.6 個數量級，實測

在 `src/initial.f90` 加全域累加器（兩個 `call nnmi` 呼叫點都累加，因為 GPU 的
`nnlist` 同時涵蓋 `ns` 與 `na`），CPU 跑批 rc=0：

```
DBGTOT sumx2_tot= 8.682365690507563E-07  nused= 13860
CPU  bal= 1.7619218910310565E-05   （後續 6.39e-6, 1.65e-6, 4.09e-7, ...）
GPU  bal= 7827.691521336932        （後續 239.43, 690.66；跑批間非決定性）
```

**能無假設直接對比的只有 `bal`**：`DBGTOT` 的 print 有 `if(myrank .eq. 0)` 守衛，
所以 `sumx2_tot` / `nused` 都**只是 rank 0 的值**（那 9 筆是 9 次 NNMI 呼叫，
不是 9 個 rank）。**不能拿 rank 0 乘以 32 去推全域總和**——各 rank 的模態分配不均。
而 `bal` 本來就是跨 rank reduce 過的，兩邊都是全域量。

| | CPU | GPU | 比值 |
|---|---|---|---|
| `bal`（第一次） | **1.76e-5** | **7827.69** | **4.4e8（8.6 個數量級）** |
| `Σx²` | 8.68e-7（僅 rank 0） | 8758（單 rank ＝全部） | 不可直接比 |
| `nused` | 13860（僅 rank 0） | 1334076 | 不可直接比 |

**`bal ≈ Σx²` 的關係兩邊都成立**（GPU：7827 vs 8758，誤差 11%；
CPU：1.76e-5 對 rank 0 的 8.68e-7，約 20 倍，與 32 個 rank 的量級相符），
再次確認 **`bal` 的機制本身健康，壞的是輸入**。

> 待查的次要異常：GPU 的 `nused = 1334076`。若各 rank 平均，CPU 全域約
> `32 × 13860 = 443520`，GPU 是它的 3.0 倍。**但這個 3 倍是用「各 rank 平均」
> 推出來的，不是量測**（§7.1），要確認得讓所有 rank 都印 `nused`。
> 無論如何 3 倍解釋不了 1e9 倍，不是主因。

###### ⚠️ 二分結果（2026-09-11）：一半無效，但挖到真正的線索

```
412  tendget_gpu 之後   Sigma phiten^2 = 133.6
413  zx_gpu 之後        Sigma phiten^2 = 0.1635
407  vartrix_gpu 之後   Sigma x^2      = 0.0005607    N = 1335248
```

**`Sigma vorten^2` / `Sigma divten^2` 顯示 `<=0` —— 這是真實量測，不是假象。**

> ### ❌ 撤回我自己稍早的「探針假象」判斷
>
> 我一度寫下「`vorten`/`divten` 不在任何 map 子句裡，所以 `update from` 是空操作」。
> **那是錯的**：我只查了 `initial_gpu` 自己的 map 子句，範圍太窄。
> 它們**有**被映射，在 `mod_spec.f90:75` 的**模組層級**：
>
> ```fortran
> vorten=0. ; divten=0. ; temten=0.        ! 主機端清零
> !$omp target enter data map(alloc:temten, vorten, divten, hldten, …)
> ```
>
> （在 `allocate_spec_array` 裡。）OpenMP 的 `enter data` 建立的映射存在於
> **裝置資料環境**中、不受語彙範圍限制，所以在 `initial_gpu` 裡 `update from` 有效。
> `tendget_gpu` 的輸入 `pt, ut, vt, tt, qt, rdiv, phi, dlpl, dtpl` 同理——
> 它們在 `mod_grid.f90:68/78/79` 的模組層級被映射，**在 `initial_gpu` 裡量得到**。
>
> **教訓**：查「有沒有被 map」要查**整個裝置資料環境**（含模組層級的
> `enter data`），不是只查當前子程式的 map 子句。§7.1 第 3 條照此修正。

所以實際狀況是：**裝置端的 `vorten` / `divten` 在 412 與 413 兩點都是零**，
而同一個 `tendget_gpu` 的第三個輸出 `phiten` 有值（1.885e8）。

**CPU 對照做完了（`src/initial.f90:98` / `:129`），這條線索作廢：**

```
412 after tendget   svor=4.67e-16   sdiv=1.65e-15   sphi=6.32e-02
413 after zx        svor=1.42e-17   sdiv=1.57e-17   sphi=1.44e-04
```

CPU 上 `vorten`/`divten` 也是數值雜訊——相對於 `sphi` 小 13~14 個數量級。
**GPU 的「恰好 0」與 CPU 的「1e-16」沒有實質差異。**

**但同一次量測給出了 `sphi` 的 CPU 基準，那才是重點：**

| 點 | CPU | GPU 跑批 B | GPU 跑批 A |
|---|---|---|---|
| 412 `tendget` 後 | **6.32e-2** | 1.885e+08 | 133.6 |
| 413 `zx` 後 | **1.44e-4** | 7.109e+05 | 0.1635 |

GPU/CPU 在 412 = **3.0e9 倍**（跑批 B）／2114 倍（跑批 A）。

> 必要的保留（§7.1）：CPU 是 32 rank、**只有 rank 0 印**、`phiten` 是
> rank 本地切片（`jtmax_l`），GPU 是單 rank 持有全部 ⇒ 有約 32 倍的結構性差距。
> **但 32 倍不是 3e9 倍。**
>
> 而且最關鍵的證據**不需要 CPU 對照**：GPU 自己兩次跑批在同一點差 1.4e6 倍。

衰減比率三者相近（CPU 438×、GPU 817× / 265×），再次確認 `zx_gpu` 是好的。

**`tendget_gpu` 一回來，`phiten` 就已經又大又不決定——目標鎖定在它身上。**

###### ‼️ `nnlist` 有未初始化的洞（2026-09-11，強線索）

兩件獨立的事對上了：

**一、`Σx²` 與 `N` 都在跑批之間改變**

| 跑批 | `Σx²`（locid 407） | `N = Σ 2·nnlist(ind)` |
|---|---|---|
| 15:14 | **8758** | 1334076 |
| 15:58 | **0.0005607** | **1335248** |

`Σx²` 差 1.6e7 倍。**而 `N` 也變了**——`N` 是在**主機端**用**主機的** `nnlist` 算的
（探針刻意沒有 `update from(nnlist)`），而填值迴圈是決定性的。
**唯一能讓 `N` 改變的，就是有元素從未被寫入、拿到的是殘留記憶體。**

**二、程式碼對得上**

`src/nvidia/initial_gpu.f90:127` 宣告 `integer nnlist(2*jtmax*nnmivm)`，
**從宣告到使用之間沒有任何初始化**（沒有 `nnlist = 0`）。填值迴圈是：

```fortran
do L = 1, nnmivm
   do m = 1, mlistnum                    ! ← 只到 mlistnum
      j = 1 + (m - 1)*2 + (L - 1)*2*jtmax   ! ← 每段 2*jtmax 長
      nnlist(j) = ns
      nnlist(j + 1) = na
   end do
end do
```

每個 `L` 的區塊長 `2*jtmax`，迴圈卻只寫其中 `2*mlistnum` 個。
**`mlistnum < jtmax` 就有 `2*(jtmax - mlistnum)` 個洞。**

算 `bal` 的 kernel 讀 `wrk(i, …)` 只到 `i <= nn`，`nn = nnlist(ind)`。
`nn` 是殘留記憶體 ⇒ kernel 走進 dgemm 從未寫過的裝置記憶體 ⇒
**同時解釋「量級離譜」與「每次跑批都不同」**。

CPU 端沒有這個陣列——CPU 在迴圈裡直接用當下的 `ns`/`na`，從不索引一個扁平的
`nnlist`，所以 CPU 免疫。這也解釋了為什麼只有 GPU 會這樣。

> ## ❌ 撤回：「`nnlist` 的裝置副本是正確的」
>
> 稍早我用 `nnlist(1)`、`nnlist(2)` **兩個元素**主機/裝置相同，就把
> 「`nnlist` 錯誤」這個假說標成推翻。陣列有 `2*jtmax*nnmivm` 個元素。
> **這是同一天第八次「取樣沒先驗證代表性」**，而且諷刺的是，
> 當時我自己在旁邊寫了「不能只看前幾個 `ind` 就下結論」。
>
> 更糟的是：主機/裝置**一致**這件事本來就與「值是否正確」無關——
> 未初始化的垃圾 `update device` 過去之後，兩邊當然一致。
> **我驗證的是傳輸，卻當成驗證了內容。**

###### ✅ `nnlist` 的洞證實存在、已修，但**不是主因**（2026-09-11）

```
nnlist 總長          = 2310   = 2*jtmax*nnmivm   = 2*385*3
填值迴圈能寫到的      = 2304   = 2*mlistnum*nnmivm = 2*384*3
未被寫入的洞          = 6      = 2*(jtmax-mlistnum)*nnmivm
實際落在 [1,no] 之外  = 6      （6 個全是垃圾，已強制歸零）
```

**事前的靜態預測精確命中**：`mpiexec -n 1` ⇒ `nsizey=1` ⇒ `jtmax = 384/1+1 = 385`，
`mlistnum = 384`，`nnmivm = 3` ⇒ 洞 = `2*(385-384)*3 = 6`。

量化佐證：歸零後 `N = 1330560`；先前兩次是 1334076、1335248，
多出來的 3516 / 4688 正落在「6 個洞各貢獻 `2*nn`、`nn ≤ no ≈ 1060`（上限 12720）」的範圍內。

**但 `bal` 完全沒有掉**：修法生效後仍是 `2225.40 / 73.93 / 193.07`。

> 依事前講死的判讀標準：**`p0 > p1` 但 `bal` 沒掉 ⇒ 洞是真的但不是主因。**
> 這個修法該留（未初始化讀取本來就是 bug，而且它確實造成了 `N` 的跑批間漂移），
> 但它解釋不了 8 個數量級。

###### ‼️ 非決定性在 `tendget_gpu` 就已經存在（2026-09-11，**目前最強的線索**）

把兩次跑批的二分欄並排：

| 探針 | 這次 | 上次 | 比值 |
|---|---|---|---|
| 412 `tendget_gpu` 之後 | **1.885e+08** | **133.6** | 1.4e6 |
| 413 `zx_gpu` 之後 | 7.109e+05 | 0.1635 | 4.3e6 |
| 407 `vartrix_gpu` 之後（`Σx²`） | 2492 | 0.0005607 | 4.4e6 |

**`Σphiten²` 在 `tendget_gpu` 剛回來時就已經差 1.4e6 倍。**
（`phiten` 確實在 map 清單裡，所以這個量測有效——不像 `vorten`/`divten`。）

而**各段的衰減比率兩次跑批幾乎一樣**：

```
412 → 413    817×（上次）   vs   265×（這次）
413 → 407    292×（上次）   vs   285×（這次）   ← 幾乎相同
```

⇒ **`zx_gpu`、`vartrix_gpu`、`nnmi_gpu` 都是好的**，它們只是忠實地傳遞一個
已經壞掉的輸入。**問題在 `tendget_gpu` 或它的上游**——比今天一整天查的位置都更前面。

###### ‼️ `tendget_gpu` 是**繼承**壞值，不是產生（2026-09-11）

在 `call tendget_gpu` **之前**量它的八個輸入（全部確認有模組層級映射）：

```
輸入 415   pt=5.905e+11    ut=0.00027     vt=3.572e-05
輸入 416   tt=1.775e+13    rdiv=0.05614   rvor=0.09951
輸入 417   dtpl=7.458e+12  dlpl=6.603e+11
------ call tendget_gpu ------
輸出 412   phiten=2.276e+08
zx 後 413  phiten=8.584e+05
407        Sigma x^2 = 3011
bal        = 2692.02 / 86.54 / 232.41
```

**呼叫之前輸入就已經是 1e11 ~ 1e13。`tendget_gpu` 本身沒有放大什麼**
（輸入 1e13 → 輸出 2.3e8，其實是縮小的）。

> ⚠️ **先不要宣告「`pt`/`tt`/`dtpl` 壞了」。** 這些是**有物理單位的場**，
> 不是傾向量：`pt` 是氣壓（~1e5 Pa）、`tt` 是溫度（~250 K），
> 對數十萬個格點求平方和，1e11 ~ 1e13 **完全可能是正常的**。
> 沒有 CPU 基準就判斷，等於重犯 §7.1 第 5 條（把沒有出處的數字當事實）。
>
> **反而另一個方向可疑**：`ut = 2.7e-4`、`vt = 3.6e-5`。
> 風場應該是 O(1~100 m/s)，平方和不該這麼小——
> 這比那幾個大數字更像異常（接近未初始化／從未被寫入）。

`qt` 與 `phi` **不在 `initial_gpu` 的 `use grid, only:` 清單裡**，作用域外，量不到；
若之後需要，得先改 use 清單。

###### ✅ CPU 基準做完 ⇒ **`tendget_gpu` 就是元凶**（2026-09-11）

CPU 是 32 rank、只有 rank 0 印、rank 本地切片；GPU 是單 rank 持有全部
⇒ **預期的結構性因子是 32**。

| 輸入 | CPU (rank 0) | GPU (單 rank) | 比值 |
|---|---|---|---|
| `pt` | 1.846e10 | 5.905e11 | **31.99** |
| `ut` | 8.424e-6 | 2.700e-4 | **32.05** |
| `vt` | 1.118e-6 | 3.572e-5 | **31.95** |
| `tt` | 5.547e11 | 1.775e13 | **32.00** |
| `rdiv` | 1.753e-3 | 5.614e-2 | **32.02** |
| `rvor` | 3.042e-3 | 9.951e-2 | 32.71 |
| `dlpl` | 2.156e10 | 6.603e11 | 30.62 |
| `dtpl` | 1.401e12 | 7.458e12 | **5.32** ← 唯一離群 |

**八個輸入有七個幾乎正好 32 倍 ⇒ 那七個是正確的。**

而輸出：

```
CPU  sphi   = 6.316e-2
GPU  phiten = 2.276e+08      比值 3.6e9  ÷32 = 1.1e8
```

> ## ❌ 撤回：「`tendget_gpu` 是繼承壞值，不是產生」
>
> 我稍早看到「輸入 `tt=1.8e13` → 輸出 `phiten=2.3e8`」就說
> 「`tendget_gpu` 沒有放大什麼，其實是縮小的」。
> **那個比較毫無意義**——`Σpt²` 和 `Σphiten²` 是兩個不同的物理量，
> 彼此的絕對大小不能推論任何事。
>
> **正確的比較永遠是「同一個量的 GPU vs CPU」**，而那個比值是 **1.1e8**。
> `tendget_gpu` 拿到幾乎正確的輸入，吐出一個大 1.1e8 倍的 `phiten`。
>
> §7.1 補第 6 條：**只比較同一個量的兩個版本，絕不跨物理量比大小。**
> 「A 比 B 大／小」在不同單位之間沒有意義，看起來卻很像論證。

**次要異常**：`dtpl` 比值 5.32 而非 32 ⇒ GPU 的 `dtpl` 小了 6 倍。
其餘七個的離散度只有 ±5%，這個 6 倍性質完全不同，**是真的異常**——
但 6 倍解釋不了 1e8 倍，且 `dtpl` 只是 `tendget_gpu` 的輸入之一。

###### `tendget_gpu` 內部結構（2026-09-11 靜態分析，零成本）

`src/nvidia/tendget_gpu.f90` 的計算鏈（L455-505）：

```fortran
init:   phiten1(k,1:2,n,m) = spalm(k)*plten(n,m,1:2)        只在 n >= mf
        cudaEventRecord(spread_event, stream)
        每個 matmul_stream(m) 先 cudaStreamWaitEvent(spread_event)
dgemm:  phiten1(1,1,mf,m) += arrhyd * temten1(1,1,mf,m)     beta = ONE，累加
        llistnum_fj = 2*(jtrun-mf+1) 欄
copy:   phiten(k,1:2,n,m) = phiten1(kk,1:2,n,m)             只在 n >= mf
```

**✅ 線索 2（`zx_buf` 切片重疊）排除——靜態即可確定：**

- `temten1` 與 `phiten1` 都是 `(lev, 2, jtrun, jtmax)`
- 呼叫端偏移 `1 + 2*lev*jtrun*jtmax` **正好等於 `temten1` 的大小** ⇒ 兩段不重疊
- `zx_buf(lev*6*jtrun*jtmax)` 容納 `2 × 2*lev*jtrun*jtmax` 綽綽有餘
- dgemm 的 `llistnum_fj = 2*(jtrun-mf+1)` 與 init / copy 的 `n >= mf` 範圍一致

**唯一的放大機制是 dgemm 的 `beta = ONE`（累加）。**
但 `tendget_gpu` 只有一個呼叫點（`initial_gpu.f90:254`），每次呼叫只累加一次，
**線性成長給不出 1e8**。所以只剩兩種可能：

1. init 沒有生效／沒有在 dgemm 之前落地 ⇒ `phiten1` 疊在前一次的結果上；
2. 輸入 `plten` / `temten1` / `arrhyd` 其中之一本來就已經很大。

> 附帶確認：`phiten` 在 L191-192 確實有清零（與 `vorten`、`divten`、`temten`、
> `hldten` 同一個迴圈）。這也解釋了為什麼量到 `vorten`/`divten` 是 0——
> 它們在這裡被清零，之後在這條路徑上沒有被填。

###### ✅ `tendget_gpu` 內部量測完成（2026-09-12）：放大在 dgemm，但 dgemm 沒錯

```
420  plten=<=0 (全零)   spalm=5.613e+05   arrhyd=3.386e+05
421  temten1=3098
422  phiten1=<=0 (全零)      ← init 之後仍是零
423  phiten1=2.466e+08       ← dgemm 之後
424  phiten =2.466e+08       ← copy 之後
```

**事先立下的作廢條件沒有成立。** 我曾寫明「422 若回報 `<=0`，就代表
`temten1`/`phiten1` 這種指進 `zx_buf` 的別名 `update from` 失效，該組數字作廢」。
422 確實是 0，**但 423/424 用同一個陣列、同一個別名機制卻讀回 2.466e8** ——
別名更新是有效的，**422 = 0 是真實量測**。

而且它自洽：`plten = 0` ⇒ `init: phiten1 = spalm*plten = 0` ⇒ 422 = 0。

**放大發生在 dgemm**：`arrhyd (3.386e5) × temten1 (3098)` → `2.466e8`。
粗估上界 `‖arrhyd‖²·‖temten1‖² ≈ 1.05e9`，實測 2.466e8 在其內
⇒ **dgemm 的輸出與它的輸入相稱，dgemm 本身沒有異常放大。**

**dgemm 的索引也核對過是對的**：`temten1(lev, 2, jtrun, jtmax)` 以 `ld = lev`
從 `(1,1,mf,m)` 取 `llistnum_fj = 2*(jtrun-mf+1)` 欄，行優先展開的順序是
`(c,n) = (1,mf),(2,mf),(1,mf+1),…`，正好對應 CPU 三重迴圈的
`n = mf..jtrun` × 兩個分量（`src/tendget.f90:344-355`）。
`mpe2d_unify_lev` 兩邊也都有（CPU L342、GPU L445）。

⇒ **問題在 `arrhyd` 或 `temten1` 的「值」，不在 `tendget_gpu` 的結構。**

###### ‼️ 找到了：`temten1` 大 3.6e6 倍（2026-09-12）

| 量 | CPU (rank 0) | GPU | 比值 | 預期 | |
|---|---|---|---|---|---|
| `spalm` | 5.613e+05 | 5.613e+05 | **1** | 1 | ✅ |
| `arrhyd` | 3.386e+05 | 3.386e+05 | **0.9999** | 1 | ✅ |
| `plten` | 1.306e-08 | **0（全零）** | — | 32 | ⚠️ 次要 |
| `temten1` | 2.685e-05 | **3098** | **1.154e+08** | 32 | ❌ **離群** |
| `phiten1` | 0.06316 | 2.466e+08 | 3.904e+09 | 32 | （下游結果） |

**鑑別器奏效**：常數類（`arrhyd`、`spalm`）比值精確為 1，譜係數類預期 32，
**唯一偏離的是 `temten1`**。

- `arrhyd` / `spalm` 精確到四位數 ⇒ **`inicons` 清白**，常數矩陣沒問題
- **`temten1` 超出預期 `1.154e8 / 32 = 3.6e6` 倍 ⇒ 這就是壞掉的量**
- `phiten1` 的超出量是 `3.9e9 / 32 = 1.2e8`，比 `temten1` 的 3.6e6 多約 34 倍。
  `Σ(A·B)²` 不會與 `ΣB²` 嚴格成正比（取決於誤差落在哪些模態、`arrhyd` 如何加權），
  所以這個差距本身不構成矛盾，但**若之後要精算，這一項要交代清楚**。

**次要發現**：`plten` 在 GPU 是零、CPU 是 1.306e-08。
它只進 init 項（`spalm*plten ≈ 7e-3`，佔 CPU `phiten1`(0.063) 約 11%），
**是真的差異，但量級上解釋不了 1e8**，別把它當主因。

###### ✅ `mpe2d_unify_lev_gpu` 清白，問題在 `temten`（2026-09-12）

GPU 跑批 `NPEX=1, NPEY=1` ⇒ `proc = nsizex = 1`；CPU 跑批 `NPEX=1, NPEY=32`
⇒ **`proc` 也是 1**（`NPEX` 才是 x 方向分解）。
兩邊的 unify 都退化成 `aout(k,n,m) = work(k,n,m,1)` 的**恆等複製**
⇒ `temten1 ≈ temten`，問題在 `temten`。

> **先排除掉一個對自己量測的疑慮**（差點推翻整個結論）：
> `temten1`/`phiten1` 是別名進共用暫存 `zx_buf`（`map(alloc:)`），而
> `_ssq` 用的是**整陣列和**，所以「未填區域帶著其他用途的舊值」是個真實風險。
> 查證結果：
> 1. `src/rocm/cudafor.f90:288` 的 `cudaMemsetAsync_r8` shim **有做元素→位元組
>    換算**（`int(n) * 8`），所以 `aout`(=`temten1`) 是整陣列清零後才填
>    `m = 1..mlistnum` ⇒ 未填的 slab 是乾淨的零，整陣列和 == 已填區域和。
>    （`_r4` 變體也乘 8，但它的 `dst` 宣告本來就是 `real(c_double)`，
>    差別只在 `val` 的精度，不是 bug。）
> 2. `phiten1` 更直接：**探針 422 讀回恰好 0**。若未填區域有舊值，422 不可能是 0。
>
> ⇒ 兩個量測都乾淨，`temten1` 的離群成立。

###### ‼️‼️ 根因：第 7 層的 bug 有 11 個常式中招，修法只涵蓋 3 個（2026-09-12）

`temten` 由 **`tranrs_gpu_cuda_graph`** 產生
（`src/nvidia/tendget_gpu.f90:385-388`，輸入是 `joinrs_gpu` 填的 `cc_cg`）。

它的 capture 區間 **L389 `cudaStreamBeginCapture` → L436 `cudaStreamEndCapture`
之間含有一個 `!$acc parallel loop`**（L391-392）。
acc2omp 把它翻成 `!$omp target teams distribute parallel do` 並**丟掉 `async`**
⇒ **正是第 7 層診斷出的 bug**：capture 期間 kernel 立即執行
（讀到 dgemm 只「錄製」而尚未寫入的緩衝區），且從未被錄進 graph，
重播時少了這一步。

**但第 7 層的修法認的是 `accx_async_begin_capture` 包裝＋三個常式名**：

```python
if any(k in name for k in ("nnmi_gpu", "tendget_gpu", "zx_gpu")):
```

**全面掃描（capture 區間內含 OpenMP kernel 的常式）：**

| 常式 | capture 方式 | 第 7 層修法 |
|---|---|---|
| `zx_gpu` (L83-152, 2 kernel) | `accx_*` 包裝 | ✅ 已涵蓋 |
| `tendget_gpu` (L456-508, 2 kernel) | `accx_*` 包裝 | ✅ 已涵蓋 |
| `nnmi_gpu` (L70-156, 1 kernel) | `accx_*` 包裝 | ✅ 已涵蓋 |
| **`tranrs_gpu_cuda_graph`** (L389-436) | 原始 CUDA API | ❌ **沒涵蓋** ← 寫 `temten` |
| `transr_gpu_cuda_graph` (L402-442) | 原始 CUDA API | ❌ 沒涵蓋 |
| `transr1_gpu_cuda_graph` (L431-480) | 原始 CUDA API | ❌ 沒涵蓋 |
| `trandv_gpu_cuda_graph` (L452-505) | 原始 CUDA API | ❌ 沒涵蓋 |
| `tranuv_gpu_cuda_graph` (L584-634) | 原始 CUDA API | ❌ 沒涵蓋 |
| `trngra3_gpu_cuda_graph` (L263-315) | 原始 CUDA API | ❌ 沒涵蓋 |
| `rstrandz_gpu_cuda_graph` (L416-469) | 原始 CUDA API | ❌ 沒涵蓋 |
| `cache_fft_plan` (L514-576) | 原始 CUDA API | ❌ 沒涵蓋 |

**11 個中招，修法只涵蓋 3 個。** 第 7 層的**診斷是對的，修法的覆蓋面錯了**。

這一條解釋了所有觀察到的現象：`temten` 壞 → `temten1` → `phiten1` → `phiten`
→ `x` → `bal` 大 8 個數量級；而 capture 期間讀到「錄製中但尚未寫入」的緩衝區，
內容取決於執行時序 ⇒ **跑批之間非決定性**。

###### ❌ 只修 `tranrs_gpu_cuda_graph` 無效（2026-09-12）

修法已套用並確認進入執行檔（生成檔 11:51:15 含 6 處 `ROCm layer-7` 標記，
執行檔 12:05:48 連結於其後 ⇒ **不是生成檔陳舊的陷阱**）。

| 量 | 修法前 | 修法後 | 預期（CPU×32） |
|---|---|---|---|
| `temten1` | 3098 | **3.329e+04** | 8.6e-4 |
| `phiten1` | 2.466e+08 | **2.861e+09** | 2.0 |
| `phiten` | 2.276e+08 | **2.861e+09** | 2.0 |

**沒有任何一項接近預期 ⇒ 依事先講死的判準，這條路不成立。**
`tranrs_gpu_cuda_graph` 的 capture 不是（唯一）原因。

> ⚠️ **不要說「修法讓事情變糟」。** 數字看起來大了約 11 倍，但 `phiten`
> 在修法前的三次跑批是 **133.6 / 1.885e8 / 2.276e8**，跨度 1.7e6 倍。
> **11 倍完全落在既有雜訊之內**，這次量測分辨不出修法是有效、無效還是有害。
> 唯一能確定的是「沒有掉到 CPU 量級」。

###### ‼️ 方法論問題：非決定性讓單次跑批失去鑑別力

`phiten` 的跑批間跨度是 **1e6 倍**。在這個雜訊地板上，**單次跑批只能分辨一種結果**：
「掉到 CPU 量級」（明確有效）。其餘任何變化都無法歸因。

後果：**「一次只改一個」這個紀律在這裡失效了**。
它能找出「那唯一的原因」，但找不出「必要的子集合」——
若真正的原因是多個常式同時中招，逐一測試會全部顯示為「無效」，
而每次要花 35 分鐘，8 個就是 4.7 小時，且結論仍是錯的。

**建議改成**：
1. 一次把修法套用到全部 8 個原始 CUDA API 的常式；
2. 若 `bal` 掉到 1e-5 ⇒ 假說成立，**再二分找出哪幾個是必要的**；
3. 若仍無效 ⇒ 第 7 層的 capture 假說整條要重新檢視，
   而不是繼續一個一個試。

###### ✅ 量到了：`temten` 由 `tranrs_gpu_cuda_graph` 寫入且大 3.9e7 倍（2026-09-12）

```
425  joinrs 前   cc_cg=0          temten=0
426  tranrs 前   cc_cg=1.259e+07  temten=0
427  tranrs 後   temten=3.333e+04
421  unify 後    temten1=3.333e+04
```

CPU 對照（`src/tendget.f90` 的 `DBGCC`）：

```
temten  tranrs 前=0   後=2.685e-05
GPU     tranrs 後=3.333e+04   比值=1.241e+09   預期 32  ⇒ 大 3.9e7 倍
```

**三件事確立：**

1. `temten` 在 tranrs **之前是 0、之後是 3.333e4** ⇒ 壞值確實由
   `tranrs_gpu_cuda_graph` 寫入，兩邊的「前值」都是 0，一致；
2. **`temten1` 與 `temten` 完全相等**（3.333e4 = 3.333e4）
   ⇒ `mpe2d_unify_lev_gpu` 在 `proc=1` 確實是恆等複製。
   **這先前只是推論，現在是量測**（§7.1 第 5 條）；
3. `cc_cg` 經 `joinrs_gpu` 從 0 變成 1.259e7。

###### ⚠️ `cc` 的比較無效：CPU 側整陣列和是 **NaN**

```
DBGCC joinrs scc_before= NaN  scc_after= NaN
```

`cc` 在 CPU 上有未初始化的 padding（堆疊垃圾含 NaN），GPU 的 `cc_cg` padding
可能是零 ⇒ **兩邊的整陣列和根本不可比**。

> **padding 陷阱第三次出現**，這次在 **CPU 側**。前兩次是 `max|mx|`（GPU 未初始化
> 填充區）與 `temten1`/`phiten1` 的疑慮（後經 memset shim 與探針 422 排除）。
> §7.1 第 2 條要補一句：**「量測範圍」兩邊都要查，不是只查 GPU。**

###### ‼️ 鑑別器結果：`tranrs_gpu_cuda_graph` 的**兩個輸出都錯**（2026-09-12）

`tendget` 裡有兩組 joinrs/tranrs 配對，共用同一份 `tranrs` 程式碼：

| tranrs 呼叫 | 輸出 | GPU | CPU (rank 0) | 應為 CPU×32 | 偏差 |
|---|---|---|---|---|---|
| 第一組 | `hldten` | **0（全零）** | 4.148e+06 | 1.33e+08 | **完全沒寫出來** |
| 第二組 | `temten` | 3.38e+04 | 2.685e-05 | 8.6e-04 | **大 3.9e7 倍** |

**同一份程式碼、不同輸入，兩次都算錯，而且方向相反。**
依事先講死的判準：**`tranrs_gpu_cuda_graph` 本身壞掉，與輸入無關**
⇒ capture 假說回到檯面，我的修法不完整。

> ### ⚠️ 我造成的證據缺口：沒有 `hldten` 的修法前基準
>
> 套用 `tranrs` 修法之前，我只量了 `temten`，沒量 `hldten`。
> 所以**無法斷定 `hldten = 0` 是既有的、還是我的修法造成的**。
>
> §7.1 補第 7 條：**套修法之前，把所有會被它影響的輸出都先量一次當基準**，
> 不是只量當下感興趣的那一個。否則修法之後就分不清因果。
>
> ~~有一個結構性理由懷疑 `hldten = 0` 是既有的：`lt_cg` 被兩個呼叫共用，
> 第二組等於在重播第一組的計算圖（不同的 `s` 輸出指標、不同的 `cc` 輸入）。~~
>
> **❌ 這個論證不成立（靜態查證後撤回）：**
> - 輸出 `s`（`hldten`/`temten`）的複製迴圈在 **capture 區間之外**，沒被錄進 graph
> - `lt_cg` 的 capture 區間內只有 `fj_wp` kernel 與 dgemm，
>   用的是 `wcc` / `fj_wp` / `wss`——**兩次呼叫是同一組模組緩衝區**
> - `cc` 只在 FFT 路徑用，屬於另一個 graph（`fft_cg`）
>
> **CUDA graph 重播的是「對同一批指標的操作」，不是「快取的結果」。**
> 緩衝區被 joinrs 重新填過，重播就會算出新值 ⇒ 共用 graph 在這裡是合法的。
>
> ⇒ **`hldten = 0` 的成因仍然未知**，而且因為缺修法前基準，可能永遠不可判。
> 真正在兩次呼叫之間會變的是 `wcc`，它由 capture 區間**之前**的迴圈
> 從 `wcc_fk` 重排而來（那段每次都執行，沒被錄進 graph）。

###### 📈 附帶觀察：修法後 `temten` 變得可重現

```
修法前   3098          （只有一個樣本）
修法後   3.329e+04
修法後   3.333e+04
修法後   3.382e+04
```

**修法後三次彼此只差 1.5%。** 而先前 `phiten`（直接由 `temten` 決定）
的跨度是 **1e6 倍**（133.6 / 1.885e8 / 2.276e8）。

⇒ 修法可能**解決了非決定性（競態），但沒解決量級**——這是兩個不同的問題。

> 保留：修法前 `temten` 只有一個樣本，嚴格說不能斷言「以前不穩」。
> 但 `phiten` 修法前的三個樣本確實跨 1e6 倍，推斷成立的機率不低。

**若這成立，單次跑批重新有了鑑別力**——先前「雜訊地板 1e6 倍讓單次跑批
無法歸因」的困境可能已經解除，後續驗證會容易很多。

###### ❌ 已排除：`tranrs` 的區域陣列**有**正確映射（2026-09-14）

起疑的理由：`wcc` / `fj_wp` / `wss` / `wcc_fk` 在 `tranrs_gpu_cuda_graph` 裡是
**區域自動陣列**（生成碼 L271-274），而整個常式只有 `twcc_fk` 有 `enter data`，
但 dgemm 那段卻用 `!$omp target data use_device_addr(fj_wp, wcc, wss)`——
`use_device_addr` 要求列出的變數已在裝置資料環境中。

**用存在性檢查（`geps_dbg_mapped`）而不是取值**，因為對沒映射的陣列下
`update from` 會安靜回傳主機副本（§7.1 第 3 條，`vorten` 陷阱）：

```
wcc      is_present=1  host=0x113279e4a010  device=0x1130fc600000
fj_wp    is_present=1  host=0x11321b017010  device=0x11309d200000
wss      is_present=1  host=0x113251499010  device=0x1130d3a00000
wcc_fk   is_present=1  host=0x1132cb1ab010  device=0x11314dc00000
twcc_fk  is_present=1  host=0x7ff694ee9e40  device=0x1127b8a00000   ← 陽性對照
```

**五個全部有映射，陽性對照也正常 ⇒ 假說推翻。**
`use_device_addr` 正常運作，dgemm 拿到的是真的裝置指標。

> **設計要點（值得重複用）**：探針裡放一個**已知有映射**的陣列當陽性對照。
> 若全部回報 0，沒有對照就分不清「映射真的壞了」與「探針本身壞了」。
>
> **成本**：260 秒就拿到答案（先前幾次都要 1000+ 秒），因為探針放在**常式入口**
> 而不是計算鏈末端。**能在上游回答的問題，不要放到下游問。**

> 附帶觀察：四個區域陣列的**主機**位址都在 `0x113...`，而 `twcc_fk` 是正常堆疊位址
> `0x7ff6...`。前者看來是 runtime 特殊配置的區域，這解釋了它們為何不需顯式
> `enter data` 就能被映射。

###### ‼️ 定位到 `wcc`：call 1 是零（2026-09-14）

`tranrs_gpu_cuda_graph` 在 `tendget` 裡跑兩次，探針 430 各記一筆：

| | `wcc` | `fj_wp` | `wss` | 該次的輸出 |
|---|---|---|---|---|
| call 1 | **0（全零）** | 237.2 | **0（全零）** | `hldten` = 0 |
| call 2 | 1.254e+07 | 237.2 | 3.398e+04 | `temten` = 3.382e+04 |

**自洽性佐證**：call 2 的 `wss = 3.398e4` 與獨立量到的 `temten = 3.382e4`
吻合（差 0.5%）⇒ 「log 順序 = 呼叫順序」的假定正確，標籤沒貼反。

**三件事確立：**

1. **`fj_wp` 兩次完全相同（237.2）且非零** ⇒ capture 區間內那個 kernel
   **有執行也有產出**。**這對「capture 吞掉 kernel」的假說是反證**
   （至少對 `lt_cg` 這個 capture 區間而言）。
2. **`wss = wcc × fj_wp` 兩次都成立**（0×237→0；1.25e7×237→3.4e4）
   ⇒ **dgemm 與 `beta=0.0` 的語意都正確**。
3. **唯一的異常是 `wcc` 在 call 1 為零。**

依事先講死的判準 ⇒ 問題在 `wcc_fk → wcc` 的重排迴圈，或上游的 **FFT 路徑**。

###### ‼️‼️ 根因確立：**FFT 的 capture 路徑產出全零**（2026-09-14）

```
call 1   fft_cg%created = 0  （capture 路徑）
         cc 前 = 1.56e+18     cc 後 = 0（全零）     gwk1 = 0（全零）
         twcc_fk = 0  →  wcc_fk = 0  →  wcc = 0  →  wss = 0  →  hldten = 0

call 2   fft_cg%created = 1  （replay 路徑）
         cc 前 = 1.28e+07     cc 後 = 1.28e+07     gwk1 = 9.108e+09
         twcc_fk = 1.264e+07 → wcc_fk = 1.264e+07 → wcc = 1.264e+07 → wss = 3.432e+04
```

**`fft_cg%created` 與失效完全對應**：
- `0`（capture）⇒ 全鏈歸零
- `1`（replay）⇒ 全鏈正常

**這先前是推論，現在是量測**（§7.1 第 7 條的要求已滿足）。

**壞的那一步是 FFT**：call 1 的 `cc` 經過
`rfftmlt_loop_identical_cuda_graph` 之後變成全零、`gwk1` 也是零；
後面的 pack kernel、NCCL transpose、fill 迴圈、dgemm **全都忠實地傳遞零**
——下游每一個環節都沒問題。

> ~~**要標明的不確定**：`cc 前 = 1.56e+18` 是整陣列和，很可能由未初始化的
> padding 主導。~~
>
> **✅ 後續追認：那個值是真的。** 2026-09-14 量到 `cc_cg` 在 joinrs **之前是全零**、
> 之後才變成 1.56e18，而 `diveng` 本身就是 1.56e18 ⇒ 它是真實資料，不是 padding。
> （當時的保留是對的——沒有證據排除 padding——但**必須回頭追認**。）

**這解釋了先前所有觀察：**

| 觀察 | 解釋 |
|---|---|
| `hldten = 0`（第一次 tranrs 的輸出） | 第一次呼叫走 capture 路徑 |
| `temten` 大 3.9e7 倍 | 它自己這次是好的，但吃到被污染的上游狀態 |
| `fj_wp` 兩次都正常 | 它在 `lt_cg` 區間內，與 FFT 無關 |
| dgemm、transpose、unify 都正確 | 它們只是傳遞零 |
| 我修 `lt_cg` 無效 | 修錯了 capture 區間 |
| 修 `lt_cg` 後 `temten` 變穩定 | `lt_cg` 不再 capture ⇒ 第一次呼叫不再特殊 |

###### ✅ 完整機制（2026-09-14，靜態確認＋量測吻合）

`rfftmlt_loop_identical_cuda_graph`（`src/nvidia/cufft_wrapper.f90`）的
capture 區間內，最後一步是 FFT 的**正規化與寫回**：

```fortran
CUDACHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal))
   …
   cufftExecD2Z(plan_id, cc(1,1,jj), gwk1(1,1,jj))     ! cc -> gwk1
   …
   !$acc parallel loop gang collapse(2) private(j, nxj) async(async_id)
   …
            cc(i, k, jj) = gwk1(i, k, jj)/float(nxj)   ! ← 正規化並寫回 cc
CUDACHECK(cudaStreamEndCapture(stream, graph))
```

acc2omp 翻譯時**丟掉 `async`**，於是：

**capture 呼叫（call 1）**

1. `cufftExecD2Z` 被**錄製**，尚未執行
2. 正規化 kernel **立即執行**，讀到 `gwk1 = 0`（剛被 `cudaMemsetAsync` 清零），寫 `cc = 0`
3. 它**沒有被錄進 graph**
4. launch 執行 graph 裡的 cufftExec，但 `cc` 已是 0 ⇒ `gwk1 = FFT(0) = 0`
5. **結果 `cc = 0, gwk1 = 0`** ✓ 與量測吻合

**replay 呼叫（call 2+）**

1. graph 裡**只有** cufftExec
2. 正規化與寫回 `cc` 的 kernel **完全不會執行**
3. **結果 `cc` 維持變換前的值，FFT 結果擱在 `gwk1` 裡沒人用**
   ✓ 與量測吻合（`cc 前 = cc 後 = 1.28e7`，`gwk1 = 9.108e9`）

> ## ‼️ 兩次呼叫**都是壞的**，只是壞法不同
>
> - call 1 → 全零
> - call 2+ → **未經變換、也未正規化**的譜資料
>
> **這就是 `temten` 大 3.9e7 倍的來源**：少了 `1/nxj` 正規化（`nxj` 是該緯圈的
> 經向點數，量級 1e3），又是未變換的原始資料。
> 先前以為「call 2 是好的」是錯的——它只是**非零**，不是正確。

> **更正先前的掃描**：常式清單上標為 `cache_fft_plan` 的那一項，
> 其實是 `rfftmlt_loop_identical_cuda_graph`。掃描用
> `re.match(r"\s*subroutine\s+(\w+)")` 比對，**誤把 `interface` 區塊裡的
> 宣告當成常式起點**。寫掃描工具時要排除 interface 區塊。

###### ✅ FFT 修法已套用：`bal` 改善 1000 倍（2026-09-14）

修法兩部分（`cmake/acc2omp.py`）：

1. **`cufft_wrapper.f90`**：拿掉 `cudaStreamBeginCapture`/`EndCapture`，
   每個 `!$omp target` 前插 `geps_acc_wait`，
   `pack_event` 的 record/wait 對換成 `cudaStreamSynchronize(plan_stream(jj))`
   （單一 event 在迴圈裡重複錄製，與 `tendget_gpu` 同一個坑）；
2. **18 個呼叫端／10 個檔案**：capture 拿掉後常式不再產生 graph，
   所以 `if (fft_cg%created)` → `if (.false.)`、
   `if (.not. (fft_cg%created))` → `if (.true.)`，
   `cudaGraphInstantiate(fft_cg%…)` 與 `cudaGraphLaunch(fft_cg%…)` 一律註解掉。
   **限定 `fft_cg%`**：各檔都把自己的 FFT graph 別名成這個區域名，
   而 LT graph 一律 `lt_cg` ⇒ 不會誤傷（已驗證 10 檔 FFT 殘留 0、`lt_cg` 保留）。

**結果：**

| 量 | 修法前 | 修法後 | 應為 |
|---|---|---|---|
| `bal` | 2225 / 7827 / 2692 | **2.86 / 0.092 / 0.160** | 1.76e-5 |
| `hldten` | **0（全零）** | 5.156e+12 | 1.33e+08 |
| `temten` | 3.4e+04 | 10.74 | 8.6e-04 |

**`bal` 改善約 1000 倍（三個數量級）**，`hldten` 從全零變成有值 ⇒
**機制上修法成功**。但距離 CPU 量級還差 1.6e5 倍，**修法有效但不完整**。

分段資料（同一次跑批）：

```
call 1（hldten）  cc FFT後=1.942e15  gwk1=1.631e21  twcc_fk=1.942e15 → wcc=1.942e15 → wss=5.156e12
call 2（temten）  cc FFT後=1.778e17  gwk1=8.979e09  twcc_fk=1.957e04 → wcc=1.957e04 → wss=36.58
```

**`cc` 不再是零、也不再等於 FFT 前的值** ⇒ 正規化 kernel 確實有跑（驗收指標達成）。

###### ‼️ 修法揭露了被零遮住的上游問題：call 1 的 `cc` 本來就是壞的

call 1 整條鏈是 1e15 量級。先前被歸零遮住，現在露出來。

**這個值是真的，不是 padding 假象**：先前量到 `cc 前 = 1.56e18` 時我標註
「可能由未初始化 padding 主導」而未採信——但 **`twcc_fk` 也是 1.942e15**，
而 `twcc_fk` 在每次呼叫開頭都被 `cudaMemsetAsync` 清零、且 pack kernel
只讀 `cc` 的使用區域（`m <= jtrunj` 有界），**所以它量到的是真實資料**。

`call 1` 的 `cc` 由 `joinrs_gpu(cc_cg, diveng, …)`（`tendget_gpu.f90:313`）填入
⇒ **嫌疑移到 `diveng` 或 `joinrs_gpu`**。

> **附帶觀察（尚未解釋）**：`hldten` 大 3.9e4 倍、`temten` 大 1.25e4 倍，
> **兩者都在 ~1e4 量級**。共同因子暗示某個系統性的正規化問題，
> 但這只是觀察，**還沒有量測支持**，不要當成結論。

###### ⚠️ 我的修法移除了自己探針的錨點

探針 431（印 `fft_cg%created` 與 FFT 前的 `cc`）錨定在
`if (fft_cg%created) then`——**而修法把那一行改寫成 `if (.false.) then`**，
於是 431 沒有被注入（生成碼中 locid 431 出現 0 次，430/432/433/434 各 1 次）。

資料並未遺失（432-434、430 都有），但我的分析器把整段輸出綁在 431 上，
導致第一次解讀時顯示「沒抓到記錄」。

**§7.1 第 7 條的新形態**：不只要在套修法前量好基準，
還要檢查**修法會不會改掉探針所依賴的錨點**。

###### ✅ `joinrs_gpu` 清白，`diveng` 進去就已經壞了（2026-09-14）

```
joinrs 前   diveng = 1.56e+18     cc_cg = 0（全零）
joinrs 後   cc_cg  = 1.56e+18
```

**`joinrs_gpu` 忠實地把 `diveng` 搬進 `cc_cg`（0 → 1.56e18），它沒有問題。**
嫌疑移到產生 `diveng` 的 **`gridnl_hybrid_ndsl_gpu_refactor`**
（`tendget_gpu.f90:306`）。

> **追認一個先前過度保留的判斷**：探針 431 量到「`cc` FFT 前 = 1.56e18」時，
> 我標註「可能由未初始化 padding 主導」而不採信。
> 現在看到 `cc_cg` 在 joinrs **之前是全零**、之後才變成 1.56e18
> ⇒ **那個值一直都是真實的 `diveng` 資料**，不是 padding。
>
> 保留是對的（當時沒有證據排除 padding），但要記得**回頭追認**——
> 否則一個真線索會被自己的保留標註永久擱置。

###### 📈 FFT 修法也解決了非決定性

```
FFT 修法前   bal = 2225 / 7827 / 2692       跨度大，跑批間不可重現
FFT 修法後   bal = 2.863（上次）/ 2.844（這次）   差 0.7%
```

**`bal` 現在可重現。** 先前「雜訊地板 1e6 倍讓單次跑批無法歸因」的困境
（見前面的方法論段落）**已經解除**——後續驗證只需單次跑批即可判讀。

這也追認了先前對 `lt_cg` 修法的觀察（`temten` 從跨 1e6 倍變成三次差 1.5%）：
非決定性的來源是 capture 期間 kernel 讀到「錄製中但尚未寫入」的緩衝區，
兩處 capture 修掉之後就穩定了。

###### ✅ `diveng` 是對的 ⇒ 上游全部清白（2026-09-14）

```
CPU  diveng = 4.8758e+16     （第二筆；第一筆是 NaN，見下）
GPU  diveng = 1.56e+18
比值 = 31.99                  預期 32  ✅
```

比值精確落在結構性因子上 ⇒ **`gridnl_hybrid_ndsl_gpu_refactor` 清白**，
`1.56e18` 就是這個量的正常量級。

> **我先前擔心它「看起來很大」是多慮**——這正是 §7.1 第 6 條要防的：
> 沒有基準就不能說大就是錯。`tt`（溫度平方和）在 CPU 上也有 5.5e11。

> **要標明的取樣問題**：CPU 第一筆 `sdiveng` 是 **NaN**
> （`diveng` 在 CPU 也是區域自動陣列，首次呼叫時 padding 未初始化），
> 用的是第二筆。嚴格說這不是與 GPU call 1 完全對應的比較，
> 但比值 31.99 的精確度讓結論可信。

**順帶驗算正規化**：`Σgwk1² / Σcc²` = `1.631e21 / 1.942e15` = `8.4e5` = `nxj²`
（`nxj ≈ 916`），與程式碼 `cc(i,k,jj) = gwk1(i,k,jj)/float(nxj)` **恰好一致**
⇒ **正規化本身正確，而且現在確實有執行**。

###### 目前的狀態（2026-09-14）

**已排除（全部經直接量測）：**

`bal` kernel、`nnmi_gpu`、`vartrix_gpu`、`zx_gpu`、`mpe2d_unify_lev_gpu`、
`arrhyd`/`spalm`（精確為 1）、dgemm 與 `beta=0.0` 語意、`use_device_addr` 映射、
`fj_wp`、`joinrs_gpu`、`gridnl_hybrid_ndsl_gpu_refactor`、`diveng`、FFT 正規化係數。

**已修（都保留）：**

| 修法 | 效果 |
|---|---|
| `nnlist` 未初始化的 6 個洞 | 真 bug，但不是主因 |
| `lt_cg` 的 layer-7 eager 化 | `temten` 從跨 1e6 倍變成差 1.5% |
| **FFT 路徑的 layer-7 eager 化** | **`bal` 改善 1000 倍＋非決定性消失** |

**剩餘落差：**

```
bal      2.844 / 2.863      應為 1.76e-5     約 1.6e5 倍
hldten   5.156e+12          應為 1.33e+08    約 3.9e4 倍
temten   10.28 / 10.74      應為 8.6e-04     約 1.2e4 倍
```

**三者都在 1e4~1e5 量級**，且上游輸入（`diveng`）已證實正確
⇒ **落差產生在 `tranrs` 的 FFT→LT 鏈內部**。

###### ‼️ 分岔點在 `twcc_fk`：FFT 階段引入約 767 倍振幅誤差（2026-09-15）

CPU `src/tranrs.f90` 內部探針（`DBGTR`）對照 GPU：

| 量 | CPU (rank 0) | GPU (call 1) | 比值 | 預期 |
|---|---|---|---|---|
| `twcc_fk`（轉置**前**） | 1.03e+08 | 1.942e+15 | **1.885e+07** | 32 |
| `wcc_fk`（轉置後） | 1131 | 1.942e+15 | 1.717e+12 | — |
| `wss`（LT 後） | 3.463 | 5.156e+12 | 1.489e+12 | — |

> ### ⚠️ 後兩項的「離群」是假象，我的分析器算錯了
>
> `mpe_transpose_rs_sp` 會把資料**重新分配到各 rank**：
> CPU 有 32 個 rank，rank 0 轉置後只拿到自己那一份
> （`twcc_fk` 1.03e8 → `wcc_fk` 1131，掉了 9e4 倍）；
> GPU 單 rank 則是恆等（1.942e15 → 1.942e15）。
>
> **所以轉置之後的量不適用「×32」的預期**，對它們硬算比值沒有意義。
> 分析器把三項都標成「離群」，其中兩項是它自己造成的。
> §7.1 第 6 條的延伸：**比較前要先確認兩邊的資料分佈方式相同**，
> 不只是物理量相同。

**唯一有效的比較是 `twcc_fk`（轉置之前）**：
GPU 應為 `32 × 1.03e8 = 3.3e9`，實測 `1.942e15` ⇒ **大 5.9e5 倍**。

換算成振幅：`√5.9e5 ≈ 767`。
而 `diveng`（FFT 的輸入）已證實正確（比值 31.99）
⇒ **FFT 階段本身引入約 767 倍的振幅誤差**。

###### 767 這個數字的候選解釋（**尚未證實**）

`nx` 是執行期變數（`src/mod_param.f90:20`），註解裡的 `nx=1536` 是被註解掉的舊設定。
由 `jtrun = 2*((1+(nx-1)/3)/2) = 384` **反推**得 `nx ≈ 1152`：

| 候選 | 值 | 與 767.4 的比 |
|---|---|---|
| `nx/2` | 576 | 1.33 ❌ |
| **`2 × jtrun`（= `⅔ × nx`）** | **768** | **0.999** ✅ |

> **這是從公式反推的，不是從跑批量到的**（§7.1 第 8 條）。
> 在直接印出執行期的 `nx` / `nxj` 之前，這個對應只是候選，不能當結論。

**先前已驗證的**：`Σgwk1² / Σcc² = 8.4e5 = nxj²`（`nxj ≈ 916`），
與程式碼 `cc = gwk1/float(nxj)` 一致 ⇒ **確實有除一次 `nxj`**。
所以問題不是「完全沒正規化」，而是**正規化的因子或次數與 CPU 的 `rfftmlt` 慣例不一致**。

###### ✅ FFT 慣例一致 ⇒ FFT 階段清白（2026-09-15）

**CPU**（`src/fftx.f90` 的 `rfftmlt`，`isign=-1`）：

```fortran
scale = 1.0/real(n)
dfftw_plan_many_dft_r2c(…)     ! FFTW r2c，未正規化
a(i,j) = pa(i,j)*scale         ! 除以 n
```

reduced-grid 呼叫端（`src/tranrs.f90:101`）傳 `1, nx+2, nxj, lev*num, -1`
⇒ **`n = nxj`**。

**GPU**：`cc(i,k,jj) = gwk1(i,k,jj)/float(nxj)` ⇒ **也是 `nxj`**。

FFTW 的 r2c 與 cuFFT 的 D2Z **都是未正規化的 forward transform**
⇒ **係數與次數完全對應，慣例一致**。

**Parseval 檢驗（`after ≈ before/nxj`）兩邊都通過：**

```
CPU:  8.42e10 / 916 = 9.19e7    實測 1.03e8    比 1.12
GPU:  1.56e18 / 916 = 1.70e15   實測 1.942e15  比 1.14
```

###### ‼️ 落差在 FFT **之前**就存在（2026-09-15）

用 NaN-safe 探針量 CPU 的 `cc`（只掃使用區域 `jj<=jlistnum`、`i<=nxj`，
非有限值另計）：

```
CPU cc_before = 8.42e+10   nonfinite = 0
CPU cc_after  = 1.03e+08   nonfinite = 0
```

| 點 | CPU | GPU | 比值 | 預期 |
|---|---|---|---|---|
| `cc` FFT **前** | 8.42e+10 | 1.56e+18 | **1.853e+07** | 32 |
| `cc` FFT **後** | 1.03e+08 | 1.942e+15 | 1.885e+07 | 32 |

**比值在 FFT 前後幾乎不變** ⇒ FFT 只是原樣傳遞，**落差在它上游**。

> ### ❌ 撤回：「FFT 階段引入約 767 倍振幅誤差」
>
> 上一輪我只有 FFT **後**的比值（1.885e7），就把它歸因於 FFT。
> 現在量到 FFT **前**的比值是 1.853e7——**幾乎相同**。
> 767 這個數字（以及 `2×jtrun = 768` 的「候選解釋」）**整個作廢**。
>
> 教訓：**要歸因於某一段，必須同時有那一段的入口與出口量測。**
> 只有出口就只能說「出口是錯的」，不能說「這一段弄錯了」。

###### ❌ 撤回：「`diveng` 是對的（比值 31.99）」

那個結論比較的是**整陣列和**：CPU 4.876e16 vs GPU 1.56e18。
但 **CPU 的 `diveng` 整陣列和第一筆就是 `NaN`**——它同樣被 padding 污染，
我卻取了第二筆並當成乾淨數據。

現在的證據指向相反方向：

- CPU `cc_before`（**使用區域，乾淨，nonfinite=0**）= **8.42e10**
- CPU `diveng`（**整陣列**）= 4.876e16
- 若 `joinrs` 是複製，兩者應該相當 ⇒ 差 5.8e5 倍
  ⇒ **4.876e16 多半是 padding 撐出來的**

GPU 那側可信：探針 425 量到 `cc_cg` 在 joinrs **之前是全零**
⇒ padding 為 0，整陣列和 == 使用區域和。

**所以真實情況很可能是：GPU 的 `diveng` 使用區域就比 CPU 大 5.8e5 倍。**

###### ❌ 撤回一整串：`tranrs` 的 CPU 量測抓到了**別的呼叫端**（2026-09-15）

`tranrs` 是**共用常式**：`getrdy.f90:779`、`intgrt.f90`（七處）、
`incrini.f90:186`，以及 `tendget` 的兩處都呼叫它。
一次跑批產生 **39 筆** `DBGCB` 記錄，而我的分析器用 `re.search` **取第一筆**
—— 抓到的是 **`getrdy` 的呼叫**（它跑在 `tendget` 之前），不是 `tendget` 的。

**作廢清單：**

| 先前結論 | 狀態 |
|---|---|
| 「落差在 FFT 之前就存在（比值 1.853e7）」 | ❌ 拿 `getrdy` 的 `cc` 比 `tendget` 的 GPU 值 |
| 「FFT 階段引入 767 倍」 | 早已撤回 |
| 「CPU 的 `joinrs` 降量級 5.8e5」 | ❌ 由 `DBGCJ` 直接推翻 |

> **這是同一個取樣錯誤在今天的第四種形態**：
> ① `seen < N` 上限、② 累積式 log 的全域計數、③ 兩側「使用區域」定義不一致、
> ④ **共用常式的多個呼叫端**。
> §7.1 第 1 條補充：**探針放在共用常式裡時，必須有辦法認出自己的呼叫端**
> ——在呼叫端前印 marker，分析器只取 marker 之後的記錄。

###### ✅ 加 marker 之後：`tranrs` 到 `twcc_fk` 為止兩邊完全一致

```
DBGMARK tranrs_1 = 3 次        DBGCB 總筆數 = 39
```

| 量 | CPU (rank 0) | GPU | 比值 | |
|---|---|---|---|---|
| `cc` FFT 前 | 4.876e+16 | 1.56e+18 | **31.99** | ✅ |
| `cc` FFT 後 | 6.037e+13 | 1.942e+15 | **32.17** | ✅ |
| `twcc_fk`（pack 後） | 6.037e+13 | 1.942e+15 | **32.17** | ✅ |

Parseval：CPU 1.13、GPU 1.14 ⇒ 一致。

**⇒ `diveng`、`joinrs`、FFT、pack **全部清白**（各自都有量測支持）。**

###### ‼️ 轉置之後的比較**方法上無效**（2026-09-15，重要）

```
wcc_fk  (transpose 後)  CPU=1.225e+09  GPU=1.942e+15   比值 1.585e+06
wss     (LT 後)         CPU=4.148e+06  GPU=5.156e+12   比值 1.243e+06
```

看起來像離群，**但這個比較不成立**。`mpe_transpose_rs_sp` 改變的是分解維度：

```
twcc_fk(…, jtmax*nsize, my_max)   本 rank 的緯度 × 全部波數
        ↓ 轉置
wcc_fk (…, jtmax, my_max*nsize)   本 rank 的波數 × 全部緯度
```

- **轉置前**按**緯度**分解 ⇒ 大致均勻 ⇒ rank 0 ≈ 全域/32 ⇒ **比值 32 成立**
- **轉置後**按**波數**分解 ⇒ 能量高度集中在低波數 ⇒ **極不均勻**
  ⇒ rank 0 與 GPU（單 rank＝全域）的比值是一個**取決於能量分佈的未知數**，不是 32

> ## ❌ 連帶撤回：`hldten 應為 1.33e8`、`temten 應為 8.6e-4`
>
> 這兩個目標值都是用「rank 0 × 32」算的，而**兩者都是轉置後的量**
> ⇒ 這個算法對它們無效，**我一直在拿錯誤的標準檢驗修法**。

###### ✅ 唯一仍然有效的全域對照：`bal`

`bal` 在 CPU 端是**跨 rank reduce 過**的全域量，GPU 單 rank 也是全域量
⇒ **兩者可直接比，不需要任何分解係數**：

```
CPU  bal = 1.76e-05
GPU  bal = 2.844 / 2.863      比值 ≈ 1.6e5
```

**落差是真的**，但**無法用 rank 0 的轉置後探針定位**。

###### ★ 下一步：CPU 探針改成跨 rank 全域加總

轉置後的量要比較，CPU 側必須 `MPI_ALLREDUCE` 成全域值，
然後與 GPU（單 rank ＝ 全域）**直接比 1:1**（不是 32）。

要改的探針：`wcc_fk`、`wss`、`hldten`、`temten`。
`src/tranrs.f90` 與 `src/tendget.f90` 都已 `use mpe` / 有 comm 可用。

- 全域 `wss` ≈ GPU `wss` ⇒ `tranrs` 整條清白，落差在更下游
- 全域 `wss` 仍差很多 ⇒ 落差在 LT 或轉置，且是真的

> **同時要重新檢視所有既有的「×32」比較**：只有**轉置前**（緯度分解）
> 的量適用。`diveng`、`cc`、`twcc_fk` 都在轉置前 ⇒ 那些 ✅ 仍然有效。

---

###### ✅ 逐波數 1:1 對照做完了：`tranrs_gpu_cuda_graph` **完全正確**，壞的是餵進第 2 次呼叫的 `ddtemp`（2026-09-15）

上一節的計畫已執行。探針設計（兩邊完全對稱）：

- CPU `src/tranrs.f90`：每個 `m` 算 `Σwcc_fk(:,:,:,m,:)²`（全緯度）與
  `Σwss(:,:,:,mf:jtrun,m)²`（**只取 `l ≥ mf`**，因為 GPU 的 dgemm 只寫這一塊、
  CPU 才是整個 `wss=0.` 清零，其餘區域不可比），`MPI_GATHER` 到 rank 0 後
  **按 `mf` 印出**（`DBGTRM call mf wcc wss`），另 `MPI_ALLREDUCE` 出全域總和（`DBGTRG`）。
- GPU（`cmake/acc2omp.py` 注入 `tranrs_gpu_cuda_graph` 的 `return` 之前）：同一組量、同一格式。
- 兩邊都在 `tendget` 的兩個呼叫點前印 `DBGMARK tendget_tranrs_N`；
  分析器（`job/cmp_mf.py`，用法 `python3 job/cmp_mf.py <cpu.log> <gpu.log>`）以**出現順序**而非印出的計數器配對
  （GPU 執行檔裡 `getrdy` 走的是 CPU `tranrs`，兩個計數器都會印 `call= 1`）。
- log：CPU `job/TCo383L72_IC_sample_rocm_cpu.log.20260915_pm`、
  GPU `job/TCo383L72_IC_sample_rocm.log.20260915_mf`。

**結果（全域對全域，1:1，不需要任何 ×32）：**

| tendget 呼叫 | 量 | CPU 全域 | GPU | GPU/CPU |
|---|---|---|---|---|
| 第 1 次（`hldten`） | `wcc_fk`（LT 輸入） | 1.9421e15 | 1.9422e15 | **1.000** ✅ |
| 第 1 次（`hldten`） | `wss` = `hldten` | 5.1559e12 | 5.1559e12 | **1.000** ✅（384 個 `mf` 逐一都是 1.000） |
| 第 2 次（`temten`） | `wcc_fk`（LT 輸入） | **0.1473** | **1.93e4** | 1.3e5 ❌ |
| 第 2 次（`temten`） | `wss` = `temten` | 3.77e-4 | 10.6 | 2.8e4 ❌ |

> ### ❌ 撤回：「`hldten` 應為 1.33e8、差 3.9e4 倍」與「`temten` 應為 8.6e-4」
>
> 上一節已經預告這兩個目標值算法無效，現在有了全域量：**`hldten` 在 GPU 上一直都是對的**
> （5.1559e12，與 CPU 全域值 5 位有效數字相同）。`temten` 確實錯，但錯的倍數是 2.8e4，
> 且原因不在 `tranrs`——同一段程式在第 1 次呼叫連 384 個波數都逐一吻合。

⇒ **FFT → pack → 轉置 → LT 整條鏈在 ROCm 上是正確的。** 錯的是第 2 次呼叫的**輸入**
`cc_cg`，它來自 `joinrs(ddtemp)`。

###### ✅ 垂直平流清白，`ddtemp` 進垂直平流**之前**就已經錯了（2026-09-15）

`ddtemp` 的來源鏈（`tendget_gpu.f90:288-375`）：
`ttp →(p2f 轉置)→ ttm_sl →(ndslfv_monoadvh2_gpu_refactor 水平半拉格朗日平流)→(f2p 轉置)→ ddtemp
→(ndslfv_monoadvv_gpu 垂直平流)→ ddtemp = (ddtemp − ttp)/dt → joinrs`。

新增共用探針 `geps_dbg_ssq_grid`（`src/tendget.f90` 末尾，兩個 build 都連進去）：
**只掃使用區域**（`i ≤ nxdef_2d(j)`、`jj ≤ jlistnum`）、NaN 另計、
`MPI_ALLREDUCE` 成全域；GPU 側由 acc2omp 在同一位置注入 `target update from` 後呼叫同一個常式。
log：CPU `…rocm_cpu.log.20260915_vadv`、GPU `…rocm.log.20260915_vadv`。

| 量（使用區域 Σ²，全域） | CPU | GPU | GPU/CPU |
|---|---|---|---|
| `pdot`、`pt`、`ttp`（垂直平流輸入） | 149.516 / 5.9054e11 / 1.7750e13 | 13 位數相同 | **1.000** ✅ |
| `ddtemp` 垂直平流**前** | 1.7747e13 | **8.5431e12** | **0.481** ❌ |
| `vdzonl` 垂直平流前 | 2.7001e-4 | 1.4613e-4 | 0.541 ❌ |
| `vdmerd` 垂直平流前 | 3.5828e-5 | 4.9322e-5 | 1.377 ❌ |
| `ddtemp` 垂直平流**後** | 1.7748e13 | 8.5404e12 | 0.481（兩邊都只動 <0.03%） |
| `ddtemp` 經 `(ddtemp−ttp)/dt` 後 | **561.06** | **1.2730e7** | 2.3e4 ❌ |

三點結論：

1. **垂直平流（`src/rocm/ndslfv_monoadvv_gpu.f90` + `vertical_cell_advect_gpu.f90` 的 tiled 改寫）行為與 CPU 一致**
   （進出比例相同），**不是**元凶。§9 第 15 項那個「host 端 CFL 重算」探針雖然該清，但與數值無關。
2. **`ddtemp`、`vdzonl`、`vdmerd` 到達垂直平流時已經是錯的**——三者正是水平平流
   `ndslfv_monoadvh2_gpu_refactor` 的三個輸出（經 f2p 轉置；單 rank 下轉置是恆等）。
3. 為什麼下游會放大到 1e4~1e5 倍：`ddtemp ≈ ttp ≈ T`（Σ² 都是 1.77e13），
   而餵給 `tranrs` 的是**兩者的差** `(ddtemp−ttp)/dt`（CPU Σ² 只有 561）。
   `ddtemp` 的任何誤差都會在這個相減裡被原樣保留、相對放大。
   這也解釋了 `bal` 對 build 極度敏感（09-14 為 2.84、09-15 同一份程式碼重新產生後為 1.0e20），
   但 `pt` 卻很穩定（−0.405/−4.547 vs −0.394/−4.542）：`bal` 已是雜訊主導的量，**看 `pt` 比看 `bal` 可靠**。

**目前的定位：`ndslfv_monoadvh2_gpu_refactor`（`src/nvidia/ndslfv_monoadvh_gpu.f90:24`）內部。**
它的 kernel 正好是 `src/rocm/` 手寫的 tiled 改寫：`cyclic_cell_intpx_gpu.f90`、
`cyclic_cell_massadvy_gpu.f90`、`cyclic_cell_ppm_gpu.f90`（acc2omp 用 `_strip_subroutine`
把原版 `cyclic_cell_ppm_intp_two_loops_gpu` / `cyclic_cell_intpx_jlist_gpu` /
`cyclic_cell_massadvy_mylonlen_gpu` 拿掉、換成這些）。**這是整個移植裡唯一大段手寫、且從未被獨立驗證過的區域。**

已加上外層 bracket（`ttm_sl`/`uum_sl`/`vvm_sl`/`pten_sl` 進出 `ndslfv_monoadvh2` 前後，
`geps_dbg_ssq_full`，使用區域 `i ≤ nxdef(j)`）；CPU 基準（`…rocm_cpu.log.20260915_hadv`）：

```
ttm_sl  pre 1.77500e13  post 1.77467e13     (post == ddtemp_pre_vadv，一致)
uum_sl  pre 2.6997e-4   post 2.6986e-4
vvm_sl  pre 3.5717e-5   post 3.5683e-5
pten_sl pre 248.408     post 248.408
```

GPU 對照見下一節。

###### ✅ 水平平流清白（GPU 對照，2026-09-15）→ 落差在它**之後**的 f2p 轉置之後出現

GPU（`…rocm.log.20260915_hadv`）：

```
ttm_sl  post_hadv  1.7746738244772e13   CPU 1.7746738244774e13   (13 位相同)
uum_sl  post_hadv  2.69859267040144e-4  CPU 2.698592670400616e-4 (13 位相同)
vvm_sl  post_hadv  3.56834554459697e-5  CPU 3.568345544597914e-5 (13 位相同)
pten_sl            0.                   CPU 248.4  ← 兩邊都沒人寫它（CPU 是堆疊殘值），與結果無關
```

⇒ **`ndslfv_monoadvh2_gpu_refactor`（含 `src/rocm/cyclic_cell_*` 手寫 kernel）完全正確**。§9 8a 對它的懷疑撤回。
但同一份 `ttm_sl` 經 `mpe2d_transpose_ndsl_f2p_gpu`（單 rank 是純複製＋k 反轉）之後的 `ddtemp` 只有 8.54e12。

#### ★★★ 第 10 層（2026-09-15 定位並修法）：`cudaMemsetAsync` 與 OpenMP kernel 的競態——第 6/7 層同一病根的第三種形態

##### 決定性證據：同一份程式碼，兩次跑批結果不同；加一個同步 kernel 就「修好」

為了看 `ddtemp` 少掉的是哪一半，我在 f2p 轉置**之前**加了一個同步的裝置端 kernel
（`DBGNXJP`：在裝置上加總 `nxjp(jlist1(j))`，與主機比對）與一組 profile 探針。
結果 **`ddtemp` 變成對的了**：

| 量（使用區域 Σ²） | 跑批 `vadv`（14:05 build） | 跑批 `f2p`（14:59 build，f2p 前多一個同步 kernel） | CPU |
|---|---|---|---|
| `ttm_sl` post_hadv | 1.7746738e13（=CPU 13 位） | 1.7747300e13（差 3e-5） | 1.7746738e13 |
| `ddtemp` post_f2p / pre_vadv | **8.5431e12** | **1.7747300e13** | 1.7746738e13 |
| `ddtemp` post `(−ttp)/dt` | 1.27e7 | **557.8** | **561.1** |
| `vdzonl` pre_vadv | 1.46e-4 | 2.22e-4 | 2.70e-4 |
| `vdmerd` pre_vadv | 4.93e-5 | 6.84e-5 | 3.58e-5 |
| `nxjp` 裝置 vs 主機 | — | 603648 = 603648 ✅ | — |

`ddtemp` 的 k-profile（72 層）與 i-profile（16 個 `i/nxj` 帶）在 `f2p` 跑批裡逐項與 CPU 相同到 4 位。
同一段程式、同一份輸入、只差一個「多花時間的同步 kernel」就從 0.48 倍變成 1.00 倍——**這是競態，不是邏輯錯誤。**
`vdzonl`/`vdmerd` 仍然不對、`ttm_sl` 這次反而差了 3e-5，說明同類的競態點不只一處。

##### 根因

`ndslfv_monoadvh2_gpu_refactor`（翻譯後 `nvidia_ndslfv_monoadvh_gpu.f90:76-84`）：

```fortran
   !$omp target data use_device_addr(uulon, vvlon, qqlon, rrlon)
   istat = cudaMemsetAsync(qqlon, 0, size(qqlon), stream)   ! HIP stream 1，立刻返回
   istat = cudaMemsetAsync(rrlon, ...)
   istat = cudaMemsetAsync(uulon, ...)
   istat = cudaMemsetAsync(vvlon, ...)
   !$omp end target data
   !$omp target teams distribute parallel do ...             ! libomptarget 自己的 queue
      qqlon(i, kk+3, lan) = ddtemp(i, k, lan)   ← 與 memset 同時在跑
```

原版 OpenACC 的 pack kernel 是 `async(async_id)`，與 memset **同一條 stream**，所以 memset 一定先完成。
acc2omp 丟掉 `async` 之後，OpenMP kernel 跑在 libomptarget 的 queue 上，與 HIP stream **沒有任何順序關係**：
memset 慢（這幾個 buffer 各數 GB），常常在 kernel 寫完之後才把一部分區域清成 0。
清掉多少取決於時序 ⇒ 0.48 倍、非決定性、`bal` 隨 build 從 2.84 跳到 1e20、加探針就變好——全部都對上了。

這與第 6 層（NCCL 集合通訊後沒等）、第 7 層（graph capture 內的 kernel 不在 stream 上）是**同一個病根的第三種表現**：
凡是「HIP stream 上的非同步操作 → 接著 OpenMP kernel 讀/寫同一塊 buffer」都會中。
第 6 層修了 `nccl_*`、第 7 層修了 FFT/dgemm，**`cudaMemsetAsync`/`cudaMemcpyAsync` 這一類一直沒修**：
`src/nvidia/` 裡有 **111 處**，`mpe2d_gpu.f90` 36 處、`tendget_gpu` 路徑上的 `ndslfv_monoadvh_gpu.f90` 10 處，
先前一個等待都沒有。

反方向（OpenMP kernel 寫 → HIP 非同步讀）是安全的：沒有 `nowait` 的 `!$omp target` 對主機是同步的。

##### 修法（`cmake/acc2omp.py`，接在第 6 層的 NCCL pass 之後；`src/rocm/hip_compat.cc` + `geps_acc_wait.f90`）

- 新增 `geps_hip_wait_all()`（`hipDeviceSynchronize`）與 Fortran 包裝 `geps_acc_wait_all()`。
  用 device-wide sync 而不是 stream wait，因為翻譯器不一定知道那個 `stream` 變數是哪個 id。
- 翻譯器對每一段連續的 `cudaMemsetAsync`/`cudaMemcpyAsync`（含 `CUDACHECK(...)` 形式、含續行）
  **在段落結束後立刻插入一行** `call geps_acc_wait_all()`——不管下一行是什麼
  （第一版把 `!$omp` 指令當註解跳過，會把等待插進 `!$omp target teams` 與它的迴圈之間，已修）。
- 未動 `src/nvidia/`。效能代價：每段 memset 一次 device sync，收斂後再與 §9 第 17 項一起評估。

##### ✅ 驗證 1（探針開著，`…rocm.log.20260915_l10`）：`ddtemp` 整條鏈變成精確

```
ttm_sl  post_hadv   1.7746738244772664e13   CPU 1.7746738244774144e13   (13 位)
ddtemp  post_f2p    1.7746738244773043e13   CPU 1.7746738244774152e13   (13 位)  ← 上一輪是 8.54e12
ddtemp  post_vadv   1.7747770501671777e13   CPU 1.7747770501670860e13   (13 位)
ddtemp  post (−ttp)/dt   561.0596876986699   CPU 561.0596876990455      (12 位)  ← 上一輪是 1.27e7
bal (iteration 1)   0.1157 / 9.5e-3 / 4.1e-3   CPU 1.76e-5 / 6.4e-6 / 1.6e-6   ← 從 2.84 改善 25 倍，仍差 6.6e3
vdzonl  pre_vadv    5.2307e-4   CPU 2.7001e-4   (1.94x) ❌
vdmerd  pre_vadv    6.8411e-5   CPU 3.5828e-5   (1.91x) ❌
```

`ddtemp` 到 `tranrs` 第 2 次呼叫的輸入為止**全部精確**。剩下 `vdzonl`/`vdmerd`：
它們的 f2p 來源 `uum_sl`/`vvm_sl` 精確（13 位），所以壞在 f2p 之後、垂直平流之前——
也就是 `vdmerdr/vdzonlr` 的傾向更新迴圈（吃 `trngra3` 的 `dlphi`/`dtphi`）與 `ndslfv_update_gpu`。

##### 假設（待驗證）：`trngra3_gpu_cuda_graph` 還在第 7 層的壞路徑上

第 7 層的註解寫得很清楚：「Seven other routines have the same defect (transr, transr1, trandv,
tranuv, trngra3, rstrandz, cache_fft_plan) — deliberately NOT touched yet」。
`trngra3_gpu_cuda_graph` 的 `lt_cg` 區塊與 `tranrs` 形狀完全相同（`pack_event` dgemm ＋ capture 內一個 unpack kernel）：
在 HIP 上第一次呼叫時 unpack kernel 在 dgemm 之前就跑了（dgemm 只是被錄下）、之後的呼叫則完全不跑
⇒ `dlphi`/`dtphi` 是垃圾。CPU 基準（`…rocm_cpu.log.20260915_upd`）顯示這裡又是一個**相消**點：
`vdzonlr` 進迴圈前 Σ² 9.4e-11，減掉 `dlphi/radsq` 之後只剩 3.0e-13 —— `dlphi` 的任何誤差都會被原樣保留。

**能定案的量測**（探針已就位、兩邊對稱）：`dlphi_pre_upd`、`dtphi_pre_upd`（CPU 1.5565e17 / 1.6343e18）、
`vdzonlr_post_loop`、`vdmerdr_post_loop`（CPU 2.995e-13 / 2.573e-13）、`vdzonl_post_upd`（CPU 2.7001e-4）。
`dlphi`/`dtphi` 對上 ⇒ 假設錯；對不上而修法後對上 ⇒ 假設成立。

**修法已實作（同一個 acc2omp pass，2026-09-15）**：把 `lt_cg` 的 eager 化從「只有 `subroutine tranrs_gpu_cuda_graph`」
放寬為「任何含 `lt_cg%created` 的檔案裡的 `*_gpu_cuda_graph` 常式」。乾跑確認 7 個檔案
（trngra3/transr/transr1/trandv/tranuv/rstrandz/tranrs）都是 eager=1、capture/instantiate/launch 殘留 0。
`tranrs1_gpu` 與 `cache_fft_plan` 沒有 `lt_cg` 區塊，不受影響。

##### ✅✅ 驗證 2（第 10 層 ＋ 第 7 層延伸，探針開著，`…rocm.log.20260915_l7x`）：**`bal` 與 CPU / NVIDIA 對上 8-9 位**

傾向更新 bracket（與 CPU `…rocm_cpu.log.20260915_upd` 對照，皆為使用區域 Σ²、全域）：

```
dlphi   pre_upd     1.5565389534636e17   CPU 1.5565389534645e17   (12 位)
dtphi   pre_upd     1.634289281429e18    CPU 1.634289281432e18    (12 位)
vdzonlr post_loop   2.995101636105e-13   CPU 2.995101636107e-13   (12 位)   ← 相消後的殘差也對
vdmerdr post_loop   2.573338218581e-13   CPU 2.573338218582e-13   (12 位)
vdzonl  post_upd    2.700071560167e-4    CPU 2.700071560169e-4    (12 位)   ← 上一輪 5.23e-4
vdmerd  post_upd    3.582817991806e-5    CPU 3.582817991805e-5    (12 位)   ← 上一輪 6.84e-5
```

**NNMI iteration 1：**

| | ROCm MI300X（本次） | NVIDIA（使用者提供） | CPU（本機 / 另一台） |
|---|---|---|---|
| `bal` l=1 | **1.7619218910e-05** | 1.7619218833e-05 | 1.7619218876e-05 |
| `bal` l=2 | **6.3945974578e-06** | 6.3945974572e-06 | 6.3945974575e-06 |
| `bal` l=3 | **1.6451316149e-06** | 1.6451316158e-06 | 1.6451316153e-06 |

前一輪（只有第 10 層）是 0.1157 / 9.5e-3 / 4.1e-3；再前一輪 2.84 / 0.092 / 0.16。
**8-9 位有效數字一致**，差異在 dgemm 加總順序的層級。

> 誠實標註：`dlphi`/`dtphi` 在套第 7 層延伸**之前**沒有量過（探針與修法同一個 build），
> 所以「`trngra3` 是元凶」是推論（輸入端除了它全部精確、修它之後全部精確），不是前後對照。
> 對結果沒有影響，但方法上不如第 10 層那次乾淨。

**同一次跑批的後續（`…rocm.log.20260915_l7x`）：**

| 階段 | ROCm | NVIDIA 參考 | 判定 |
|---|---|---|---|
| NNMI iteration 2 `bal` | 4.0898e-07 / 3.2222e-07 / 9.0687e-08 | 4.0898e-07 / 3.2222e-07 / 9.0687e-08 | ✅ |
| NNMI iteration 3 `bal` | 4.9302e-08 / 2.9515e-08 / 1.1525e-08 | 4.9302e-08 / 2.9515e-08 / 1.1525e-08 | ✅ |
| `pt` 三次更新 | 0.91441/−7.67721 → 0.91410/−7.67534 → 0.91483/−7.67639 | 完全相同 | ✅ |
| `qgini, thdai, tengi` | 30.7242371636 / 4789.7340632 / 3544561.6535514 | 30.7242371636 / 4789.7340632 / 3544561.6535514 | ✅ |
| 積分第 1 步 `surf pres tend rms(GPU)` | **132.62 mb/hrs** | 0.4195（CPU 0.4232） | ❌ 有限值但差 ~300 倍（以前是 NaN） |
| 積分第 2 步 | `Failure to allocate device memory` → abort | — | 舊問題（§9 第 13 項的 VRAM 不歸還），與數值無關 |

⇒ **NNMI 在 ROCm 上已經正確。** 剩下的數值問題搬到 `intgrt_gpu`（積分第一步就差 300 倍），
外加獨立的 VRAM 耗盡讓跑批活不過第 2 步。

##### 尚未做、但一定要做：**關掉探針**的確認跑批

這次驗證的 build 帶著 `tendget_gpu` 裡十幾個 `!$omp target update from(...)` 探針，
它們本身就是同步點，**理論上可能遮住殘留的競態**。要把 acc2omp 裡 `# #region agent log` 的注入拿掉
（或加一個環境變數開關）重編、重跑，`bal` 仍是 1.7619e-05 才算修法自己站得住。
CPU 側 `src/tendget.f90` / `src/tranrs.f90` 的探針不影響 GPU 路徑，可以晚點清。

#### ★★ 第 11 層（2026-09-16）：積分第 1 步——`intgrt_gpu` 全步的 `trandv` 呼叫傳錯 `poly`（graph capture 遮住的原始碼 bug）

##### 方法：同一套 checkpoint 表，兩邊 1:1

`cmake/intgrt_ckpts.py` 是一張表（CPU 錨點、GPU 錨點、第幾次出現、之前/之後、要量的陣列），
同時餵給 `src/intgrt.f90`（python 直接改）與 acc2omp（`intgrt_gpu` 區段），兩邊 tag 完全相同；
`job/cmp_ckpt.py cpu.log gpu.log` 取每個 tag 的第一次出現（= 第 1 步）比 GPU/CPU。
新增兩個頻譜 helper：`geps_dbg_ssq_spec`（`(n1, jtrun, jtmax)` 視角，只掃 `n ≥ mlist(m)`）與
`geps_dbg_ssq_spec2`（`plten` 的 `(jtrun, jtmax, n3)` 版），都在 `src/tendget.f90` 末尾。
（第一版注入器對「同一錨點有第 1 次與第 2 次兩筆」會重複計數，第 2 次永遠配不到——已修。）

##### 靜態發現（在量測之前就找到了）

`src/nvidia/intgrt_gpu.f90:1085`（全步的 `trandv_gpu_cuda_graph`）傳的是 **`poly, dpoly`**——
`mod_const` 裡 CPU 版面 `(jtrun, my/2, jtmax)` 的三維陣列——而該常式的 dummy 是 packed 的 `polyf` 版面
`((jtrun+nsize)*my/2*jtmax)`；半步的同一個呼叫（866 行）與其它所有 `*_cuda_graph` 呼叫點都傳 `polyf, dpolyf`。
`trandv` 在 capture 區塊**內**用 `poly(ind)` 算 `fj_wp/fj_wd`：

- NVIDIA：graph 錄的是**第一次呼叫的指標**（polyf），之後 replay 永遠用它 ⇒ 錯的參數從來沒被讀過，bug 被遮住。
- ROCm（第 7 層 eager 化之後）：每次呼叫真的去讀傳進來的陣列 ⇒ 全步的 `vorten/divten` 是垃圾。

這是「eager 化改變了 NVIDIA 程式碼依賴 graph 指標凍結的語意」的第一個實例；`grep` 過所有 `*_cuda_graph` 呼叫點，
只有這一處。修法放在 acc2omp（`intgrt_gpu` 區段把那一行改成 `polyf, dpolyf`），未動 `src/nvidia/`；
**上游應該直接改原始碼**，因為它在 NVIDIA 上也只是碰巧沒事。

##### 驗證（`…rocm.log.20260916_ig` vs `…rocm_cpu.log.20260916_ig`）

`surf pres tend rms(GPU)`：**132.6 → 1.987** mb/hrs（本機 CPU 0.438）。
Checkpoint（GPU/CPU，都是全域使用區域 Σ²）：半步全部 1.0000（`s1_fgnl` → `s1_gridnl` → `s1_trngra3` → `s1_vadv` →
`s1_tranrs1` → `s1_trandv` → `s1_siimpl` → `s1_hdiffu` → `s1_mid` 共 28 個量）；全步 `s2_hadv_*`、`s2_vadv_*`、`s2_pt` 也是 1.0000；
**`s2_ptend` 19.9 倍、`s2_ptend_plten` 20.6 倍** ❌。所以剩下的落差在 `pt` 之後、`ptend` 之前：
`(tt−ttp)/dta` 迴圈 → `tranrs`（temten）→ `trandv`（vorten/divten）→ `pcorr` → `siimpl` → `transr1`。
這次 build 少了 `s2_trandv`/`s2_siimpl` 四個點（注入器 bug），下一輪補上並加 `s2_tendin_*`、`s2_tranrs_temten`、
`s2_siimplin_*`、`s2_pcorr`。NNMI 三輪 `bal` 與昨天**逐位元相同**（跑批間決定性）。跑批仍在第 2 步死於 VRAM 耗盡。

##### 第 11 層之後的細 bracket（`…_ig2` 跑批）：落差在 `(ut−up)/dta`，`ut` 精確、`up` 不對

| checkpoint | GPU/CPU |
|---|---|
| `s2_gridnl`、`s2_trngra3`、`s2_hadv_*`、`s2_vadv_*`、`s2_pt`（狀態場，相對差 1e-13） | 1.0000 ✅ |
| `s2_tendin_ddtemp` = `(tt−ttp)/dta` | 1.0001 ✅ |
| **`s2_tendin_vdzonl` = `(ut−up)/dta`** | **3.1e3** ❌ |
| **`s2_tendin_vdmerd` = `(vt−vp)/dta`** | **4.3e2** ❌ |
| `s2_tranrs_temten`、`s2_pcorr`、`s2_siimplin_{temnow,divnow,plnow,plten}` | 1.0000 ✅ |
| `s2_trandv_vorten` / `divten`（吃上面兩個） | 34.6 / 6.2 ❌ |
| `s2_ptend` / `s2_ptend_plten` | 19.9 / 20.6 ❌ |

`ut`、`vt`（`s2_vadv_ut/vt`）與 CPU 相對差 1e-13，所以錯的是 `up`/`vp`——上一時間層的風。
兩邊都只在時間步開頭寫它一次：CPU `intgrt.f90:632` 的迴圈 `up = ut`，GPU `intgrt_gpu.f90:722`
的 `cudaMemcpyAsync(up, ut, …, DeviceToDevice)`（在 `host_data use_device` 裡，翻譯後是 `use_device_addr`）。

#### ★★★ 第 12 層（2026-09-16，已實作、待驗證）：acc2omp 把 OpenACC 的 `exit data delete` 翻成 OpenMP 的 `map(delete:)`——語意不同

- `mod_grid.f90:82` 在模組層級 `enter data create(… up, vp, ttp …)`（引用計數 1）。
- `initial_gpu.f90:167` 再 `enter data create(rvor, up, vp, ttp, qm, …)`（→ 2），結尾 453 行 `exit data delete(…)`。
- **OpenACC 的 `delete` 是「引用計數減 1，歸零才釋放」**（沒有 `finalize` 就不會強制）⇒ NVIDIA 上 `up` 之後仍然在裝置上。
- **OpenMP 的 `map(delete:)` 是「直接歸零並釋放」**；等價於 OpenACC `delete` 的是 **`map(release:)`**。
- acc2omp 第 86 行把所有 `delete(` 翻成 `map(delete:`（翻譯後 221 處）⇒ `initial_gpu` 結束時 `up/vp/ttp/rvor` 被整個拆掉。
- 之後 `intgrt_gpu` 裡：`use_device_addr(up)` 對不存在的變數給的是**主機位址** ⇒ `hipMemcpyAsync` D2D 失敗（`istat` 沒人檢查）
  ⇒ 裝置上的 `up` 根本沒被更新；而 `ttp`「看起來對」是因為 kernel 引用未映射的陣列會被**隱式逐 kernel 映射**
  （主機副本來回搬），`ttp` 走 kernel、`up` 走 HIP memcpy，命運不同。
  這也解釋 `(ut−up)` 是「小而真實」的差（NNMI 末的 `ut` vs 現在的 `ut`），而不是垃圾。

**修法**：acc2omp 第 86 行 `delete(` → `map(release:`。`src` 裡沒有任何 `finalize`，所以沒有地方真的需要 `delete`。
翻譯後 `initial_gpu` 的那行變成 `map(release:rvor, up, vp, ttp, qm, ut_sl, vt_sl)`。

> 誠實標註：`up` 本身沒有在修法前量過（`s0_copy_up`/`s2_tendin_up` 探針與修法同一個 build），
> 「`up` 未映射」是從翻譯語意靜態推出來的，不是用 `geps_dbg_mapped` 量到的。下一輪若 `s2_tendin_vdzonl` 變 1.000、
> `sptend` 變 0.438，這條鏈就成立；否則要回頭量 `up` 的映射狀態。

**可能的副作用（正面）**：以前被錯誤拆掉的大陣列（`up/vp/ttp/rvor/qm`…）在積分裡每個 kernel 都隱式 alloc+copy+free
一次（每個 ~700 MB），這很可能就是 §9 第 13 項「VRAM 不歸還」的來源之一——`release` 之後這些搬運全部消失。下一輪看 `vram:` 軌跡。

##### ✅✅ 驗證 3（第 12 層，`…rocm.log.20260916_rel` vs `…rocm_cpu.log.20260916_ig3`）：**積分第 1 步與 CPU 相同到 14 位**

```
surf pres tend rms   GPU 0.43807396825059497   CPU 0.4380739682505907
```

歷程：132.6（09-15）→ 1.987（第 11 層 trandv 的 `polyf`）→ **0.43807396825059**（第 12 層 `release`）。
57 個 checkpoint（`s0_copy_*`、`s1_*`、`s2_*`，含 `s2_tendin_vdzonl/vdmerd`、`s2_trandv_vorten/divten`、`s2_ptend`）
全部 GPU/CPU = 1.0000；唯二例外是兩邊都沒人寫的 `pten_sl`。NNMI 三輪 `bal` 仍逐位元相同。
`up` 未映射的推論由結果反證成立（`(ut−up)/dta` 從 3.1e3 倍變 1.0000）。

**跑批仍在第 2 步死於 `Failure to allocate device memory`**——`release` 沒有改善 VRAM（NNMI 階段軌跡與之前逐點相同）。
數值問題到此為止都清了；**VRAM 是現在唯一的 blocker**，見下一節。

##### VRAM 耗盡（§9 第 13 項）的新線索（2026-09-16，從 `vram:` 軌跡讀出）

以前一輪完整跑批的 `freeB` 軌跡，列出每次 ≥ 3 GiB 的變化：

| 時間 | free（GiB） | 區間 |
|---|---|---|
| 2.2 min | 191 → 159 → 147 | initial_gpu 的 enter data |
| 3.6 | 147 → 102 | initial_gpu:417 → vram:130（tendget 水平平流第一次進入） |
| 4.3 | 98 → 124 | vram:718 → 244（水平平流離開，**歸還 26**） |
| 5.8 | 122 → 103 | vram:436 → 432（`tranrs_gpu_cuda_graph` 第一次） |
| 8.2 | 103 → 85 | vram:428 → 245（trngra3 + 垂直平流第一次） |
| 19.9 | 85 → 49 | vram:430 → 427（tranrs 第 2 次呼叫 → unify） |
| 25.9 | 55 → 22 | initial_gpu:404 → 340（correct → transr×3/ujoinsr×3/transr1 第一次） |
| 28.5 | 22 → 6 | initial_gpu:340 → 350（**`trngra_gpu` 第一次，16.6 GiB**） |
| 28.6 起 | 0 ↔ 5.8 震盪 | NNMI 第 2、3 輪：**不再增加** |

**模式：每個常式的「第一次呼叫」永久吃掉 15–45 GiB，之後的呼叫不再消耗。** 所以不是每步洩漏，是**第一次呼叫建立的狀態從不釋放**，
而且 OpenMP 映射表持平（先前已量）⇒ 在 HIP 層。候選（都是 first-call、per-routine 的）：

1. **hipfft plan**：`rfftmlt_loop*` 對每個緯度建一個 plan（`jlistnum=768` 個），每種 `(n, batch, isign)` 一組；
   `fft_create_plan(plan, 1)` 讓 auto-allocation **開著**，rocFFT 對每個 plan 配自己的 work buffer。
   若每個 plan 的 work buffer 是幾十 MB，768 個就是幾十 GiB——**量級對得上**。已在 `hip_compat.cc:geps_fft_make_plan_many`
   加了 `[fftplan]` 計數（work_size 與 `hipMemGetInfo` 差值，每 256 個印一次），下一個 build 就能定案。
2. `acc_get_cuda_stream(m+1)`：tranrs/trngra3 等對每個波數各開一條 stream（385~768 條），加上 FFT 每緯度一條。
3. hipblas handle / workspace。

#### ★★★ 第 13 層（2026-09-16，VRAM 耗盡的根因，已量測）：每個緯度一個 hipfft plan，rocFFT 每個 plan 固定吃 ~20 MB

##### 量測（獨立小程式，scratchpad `fftplan*.cc`，在閒置的 GPU1 上）

| 情境 | 384 個 plan 的 VRAM | 每個 plan |
|---|---|---|
| 同一個 n=1552、D2Z、batch 72 | 9.52 GB | 25.4 MB |
| 同一個 n=1024（2 的冪） | 6.25 GB | 16.7 MB |
| 384 個不同 n、batch 72 | 7.33 GB | 19.6 MB |
| 384 個不同 n、batch 1 | 7.18 GB | 19.2 MB |
| 384 個不同 n、Z2D | 6.80 GB | 18.1 MB |

`hipfftMakePlanMany` 回報的 work_size 只有 0–0.9 MB；**實際每個 plan 固定 17–25 MB**，與 n 是否重複、batch、方向都無關，
`hipfftDestroy` 才會還。建 384 個 plan：冷 42 s（rocFFT 執行期編譯 kernel，之後有磁碟快取）、熱 4.4 s；銷毀 1.9 s。

##### 模式裡的用法

`rfftmlt_loop_identical_cuda_graph`（17 個呼叫點）對**每個緯度**建一個 plan（`jlistnum=768`），每種 `(jump, batch, isign)` 一組；
NVIDIA 上這樣做是為了讓每個緯度的 FFT 在自己的 stream 上並行（cuFFT 的 plan 很便宜）。
一組 = 768 × ~22 MB ≈ **16.6 GB**，正好是 `trngra_gpu` 第一次呼叫吃掉的量；NNMI ＋ 積分共 ~9 組 ⇒ ~150 GB。
這就是「每個常式第一次呼叫永久吃 15–45 GB」的模式。

##### 修法（兩階段，都在 `src/rocm/rfftmlt_loop_gpu.f90`，由 acc2omp 在 layer-7 pass 之後 `_replace_subroutine` 進去）

1. **每個不同的 `nxj` 一個 plan**（plan 只跟 `(jump, nxj, m, isign)` 有關，與緯度無關）。實測 reduced grid 有 **384 個不同長度**
   （南北對稱各一半），所以一組從 768 個減到 384 個、~8 GB；7 組 = 46.8 GB（`[fftplan]` 計數器實測）。
   跑批 `…rocm.log.20260916_fft`：NNMI 結束時 free 從 0 變 **67 GB**，但積分第 2 步（物理第一次呼叫）仍 OOM。
2. **plan-set LRU**：最多保留 `GEPS_FFT_PLAN_SETS` 組（預設 4），要新組時銷毀最久沒用的（印 `[rfftmlt_loop rocm] evicted set`）。
   每步 ~9 組輪流用 ⇒ 每步約 6 組要重建（熱 4.4 s ＋ 銷毀 1.9 s 各）⇒ **每步多 ~40 s**，先求跑得完；
   結構性的解（減少不同長度數、或一個 plan 服務多個長度）之後再做。狀態放在 subroutine 的 `save` 變數裡，不再用 NVIDIA 的 C plan cache。

`src/CMakeLists.txt` 已把 `rocm/rfftmlt_loop_gpu.f90` 加進 `cufft_wrapper` 的相依（否則產生檔不會重生，§6 陷阱 1）。
**改 `src/CMakeLists.txt` 會觸發 reconfigure，reconfigure 前一定要 `export COMP_MP=mpifort C_COMP_MP=mpicc CXX_COMP_MP=mpicxx`**（§8 有指令；2026-09-16 踩到，link 找不到 `mpi_*` 且執行檔被刪）。

##### 第 13 層之後（2026-09-16 下午）：單 rank 裝不下物理階段 → 改跑 2 GPU；物理的映射工作集實測

- `GEPS_FFT_PLAN_SETS=3`（`…_lru`）：NNMI 結束 free 96 GB，積分第 1 步後仍 OOM；`=1`（`…_lru1`）：free 78.6 GB，仍 OOM。
  上面 09-10 的映射表量測（積分 78 → 154 GB）已經說明單 rank 的物理工作集 > 75 GB，剩多少都不夠。
- **2 GPU**（`job/TCo383L72_IC_sample_rocm_2gpu`，NPEY=2，`regression.ksh` 新增 `GEPS_ROCM_DEVICES=0,1` 開關；
  `device_init` 本來就是 `mod(rank, ndevices)`）：第一次撞到 `src/rocm/ndslfv_para_gpu.f90` 的 stub
  「nsizey>1 needs RCCL works/workr」——當初只寫了單 rank 路徑。已把 NVIDIA 原版的多 rank 路徑
  （pack → `nccl_alltoallv_stride` → unpack，`works/workr` 暫存）以 OpenMP 手寫成 `para_we2ns_gpu_multi` / `para_ns2we_gpu_multi`，
  `nsizey==1` 仍走省記憶體的快路徑。
- **2 GPU 結果（`…_2gpu`）**：RCCL 多 rank 轉置全部正確——NNMI 三輪 `bal` 與單卡**逐位元相同**（1.7619218910310247E-05…）、
  `pt` 正確、積分第 1 步 `surf pres tend rms` = 0.43807396825059497（與 CPU 相同）。**但物理第一次呼叫仍 OOM**，
  且失敗時每張卡還有 ~100 GB free。
- 碎片化排除：獨立小程式棋盤式打碎 9189 個 20 MB 區塊後，1–64 GB 的 `hipMalloc` 全部成功。
- 新探針（acc2omp 全域 pass）：每個 `target enter data` 前印 `DBGMAP <file>:<n> GiB=`（`sizeof` 總和 > 1 GiB 才印）。
  2 GPU、單一 rank 的積分物理階段：**>1 GiB 的 enter data 合計 253 GB**——`radiation_aerosols_gpu` 129 GB（46 筆，其中
  `alon…ssaaer` 那組 2.9 GB 是每次呼叫就 exit data 的暫存 ×32）、`radsw_main_gpu` 61 GB、`diabat_gpu` 18、`moninedmf` 16、
  `mp_scheme` 9.7 …。失敗點在 `mp_scheme_gpu` 一個 1.6 GB 的映射 ⇒ 前面的東西（至少輻射的一部分）到那時還活著。
  下一輪探針加上 `exit data`（`DBGUNMAP`）與每筆的 `hipMemGetInfo` free，直接算活著的量。

**推論（待下一輪確認）**：物理階段的裝置工作集是 O(100+ GB)/rank，NVIDIA 生產環境顯然是 ≥ 4–8 張卡分攤。
這台只有 2 張 → 若確認是真的量，選項是 (a) 物理改成分緯度批次呼叫（大改），(b) 找出宣告過大的陣列（例如用 `nx`/`my` 而非 `nxp`/`my_max`、
或多了用不到的維度），(c) 更多卡。

#### ★★★ 第 14 層（2026-09-16 晚，§9 第 13 項「VRAM 不歸還」的根因，小程式重現）：libomptarget 的記憶體管理器把大塊釋放的 buffer 留在池子裡

##### 從模式 log 看到的訊號（`…_2gpu_lw`，256 MB 門檻的 `DBGMAP/DBGUNMAP` ＋ 每筆的 `hipMemGetInfo`）

`radiation_aerosols_gpu` 每次呼叫 `enter data create(alon…ssaaer)` 2.93 GB、用完 `exit data`：
exit 之前 free 63.41，exit 之後、下一次 enter 之前 63.18，再 enter 之後 63.18——**釋放不還、再配不扣**，
典型的「同尺寸池子」行為。而 `radsw`/`radlw` 每個緯度 block 的欄位數都不同 ⇒ 每次呼叫都是新尺寸 ⇒ 池子只長不縮：
radsw 每 block 掉 3–5 GB、radlw 每隔一次掉 7–10 GB，到 `mp_scheme` 時 free 只剩 1 GB。
（碎片化與隱式映射兩個假設先前已用小程式排除。）

##### 獨立重現（scratchpad `pool.f90`，注意 libomptarget 只認 `ROCR_VISIBLE_DEVICES`，不認 `HIP_VISIBLE_DEVICES`——第一版量到別張卡）

8 次 `enter data map(to:a)` → kernel → `exit data map(release:a)`，每次 `a` 尺寸不同（1.54 → 1.91 GB）：

| | free（GB）逐次 |
|---|---|
| 預設 | 191.3 → 189.3 → 187.7 → 186.0 → 184.3 …（**每個新尺寸永久累加**） |
| `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0` | 每次 exit 後回到 190.8 |

##### 修法

`job/regression.ksh` 的 rocm 分支 `export LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0`（可被環境覆蓋）。
09-10 時試 `=0` 會「提早在別處炸」——那是第 6 層 NCCL use-after-free 被池子遮住（釋放的 buffer 位址還有效）；
第 6/7/10 層修掉之後這個開關才變得安全。**這也解釋了 09-10「映射表持平、HIP free 卻單調下降」的全部觀察**
（NNMI 階段的那部分是第 13 層的 hipfft plan，積分階段的是這一層）。

##### 第 14 層驗證（`…_2gpu_mm0`）＋ 它揭露的第 15 層

- 池子關掉後：進積分 free **141 GB**（之前 107），物理一路走到 `mp_scheme` 還有 **89 GB** ⇒ 記憶體問題解決，2 卡夠用。
- **但** `surf pres tend rms` 從 0.43807396825059497 變成 0.4380625890241813（第 5 位不同）。
  記憶體池的開關不該改變數值 ⇒ 還有一個「釋放後仍被讀」的競態被池子遮著（跟第 6 層當年一模一樣：池子讓已釋放的位址暫時仍有效）。
- 接著在微物理（`module_mp_gsfcgce_3ice_nuwrf_gpu`）**卡死**：兩張卡 100%、51 分鐘沒有輸出。
  它的沉降子步迴圈 `DO while (notlast)` 用 `del_tv = min(del_tv, z1/vtrr)` 累加到 `dtb`——輸入含 NaN/垃圾就永遠不收斂。已 kill。

#### ★★ 第 15 層（2026-09-16 晚，已實作、待驗證）：graph 區塊之外的 hipBLAS `dgemm` / hipFFT exec 後面沒有等待

掃描產生檔（scratchpad 腳本：非同步 HIP 呼叫 → 下一個 `!$omp` 構造或釋放之間沒有 wait），扣掉沒人呼叫的常式後剩下：
`correct_gpu`（`dgemm` → 下一個 kernel 讀結果）、`tranuv1_gpu`（同）、`cublas_wrapper` 的 `dgemm_async`（20 個呼叫端）、
`rfftmlt_gpu` 的 `fft_exec_async`。hipBLAS 排在 handle 目前的 stream（被 `cublasSetStream` 設成任何一條 `lt_cg_stream(m)`）就返回，
OpenMP kernel 沒有順序保證。**修法**：不逐點補 wait，直接在 `src/rocm/hip_compat.cc` 的 `geps_blas_dgemm` 與四個 `geps_fft_exec_*`
之後 `hipDeviceSynchronize()`——所有 BLAS/FFT 呼叫變成同步。並行本來就被 eager 路徑序列化了，代價很小。
驗證標準：`sptend` 回到 0.43807396825059（與池子開著時、與 CPU 相同），且微物理不再卡死。

##### 第 15 層驗證（`…_2gpu_l15`）：BLAS/FFT 同步化**沒有**改變結果 ⇒ 競態在別處

`sptend` = 0.4380603838（上一輪 0.4380625890，池子開著 0.4380739683）——仍非決定性。
但 step-1 checkpoint 還在，直接定位：`s1_siimpl_*` 之前全部與 CPU 相對差 < 1e-11，
**`s1_hdiffu_vormid/divmid/temmid` 開始偏 9.7e-4 / 5.2e-3 / 1.8e-6**，之後全部繼承；另外 `s2_pcorr`（`ptotc_gpu`）相對差 4.9e14 = 垃圾。
微物理跟著再次卡死（同一個沉降迴圈）。

#### ★★★ 第 16 層（2026-09-16 深夜，已實作、待驗證）：acc2omp 沒有把「迴圈內被賦值的純量」私有化——OpenACC 預設 firstprivate、OpenMP 預設共享

`hdiffu_gpu.f90:98`：`!$acc parallel loop gang private(KL, kfac, dect, dec, facd, facv, fact)`，
但迴圈體裡還賦值 `hfiltd`、`powdd`、`c1`、`c2`（依 k 而異），原作者沒列——**OpenACC 對 `parallel` 區域裡未列的純量預設每個 gang 一份（firstprivate），
所以在 NVIDIA 上沒事**；翻成 OpenMP `target teams distribute parallel do` 之後這些純量是 target 層級的單一副本、各 team 共用 ⇒ 跨 team 競態，
正好給出「每層係數偶爾錯」這種 1e-3 級、看時序的偏差。第 9 層只補了內層 `!$acc loop` 的 private 與迴圈變數，沒補這一類。

**修法（`cmake/acc2omp.py` 新 pass `_privatise_assigned_scalars`，接在 `_privatise_nested_loops` 之後）**：
對每個 `parallel loop`，掃它的迴圈巢狀內的賦值目標（含 `if (…) x = …` 單行形式），取其中在同一程式單元宣告為**純量**的名字
（單行型別宣告、無 `dimension/allocatable/pointer/parameter`），扣掉已在任何子句裡的，加成 `firstprivate(...)`
（用 firstprivate 而非 private：忠實對應 OpenACC 語意，也涵蓋「進迴圈前有初值、迴圈裡條件式改寫」的情況）。
cpp 條件深度比照第 9 層（宣告或賦值在 owner 之外的 `#ifdef` 裡就不提升——`mr2mc` 陷阱）。
翻譯全部 GPU 原始碼：1055 個 `target teams` 指令、0 個子句重名；`hdiffu` 那條變成
`… private(KL, …) private(m) private(n) firstprivate(hfiltd, powdd, mf, c1, c2)`。`radlw_main` 多了 19 條、微物理 6 條、`diabat` 4 條。

驗證標準同第 15 層：`sptend` 回到 0.43807396825059、`s1_hdiffu_*` 相對差 < 1e-9、微物理不卡死。
`ptotc_gpu` 的 `pcorr` 垃圾要另外看（它是 n=1 項，不進 `sptend`，但進 `plten`）。

**第 16 層結果（2026-09-16 19:54，`job/TCo383L72_IC_sample_rocm.log.20260916_2gpu_l16`）：沒有改變任何數字。**
`sptend` 仍是 0.4380603838412094（與第 15 層逐位元相同、決定性），`s1_hdiffu_vormid` 仍偏 9.7e-4，
微物理仍卡死（sedimentation `DO while (notlast)`）。所以第 16 層是正確但無關的修正（它修的是理論上的競態，
這個 case 沒有觸發），真正的原因在別處。

#### ★★★ 第 17 層（2026-09-16 20:09 建好）：兩個殘留的 HIP-stream 競態——`NCCLCHECK(ncclAllReduce…)` 裸呼叫、`hdiffu` 讀 host 端 `wmax`

第 6 層只在 **`call nccl_*`（包裝副程式）** 之後補 `geps_acc_wait`；但翻譯後的程式裡還有 **15 處直接寫
`NCCLCHECK(ncclAllReduce(...))` / `ncclAllGather` 的裸呼叫**：`mpe2d_gpu.f90` 12 處（`mpe2d_unify_nx_gpu`、
`mpe2d_unify_my1d_gpu` 等——`ptotc_gpu` 算 `pdry`→`pcorr` 就走這裡）、`hdiffu_gpu.f90:78`（`wmax` 全域最大值）、
`rayleifr_gpu.f90:60`（同型）、`intgrt_gpu.f90:452`（`wk4` → `tengi/qgini` 診斷）。
它們排在 HIP `stream` 上，後面緊接著讀結果的 OpenMP kernel（`wmax(i) = wmax_buf(i)`）或同步 `exit data map(from:)`——
和第 6 層同一類：**collective 還沒跑完就被讀走**。這可以解釋：
- `hdiffu` 各層係數偶爾錯（`wmax` 沒 reduce 完就決定 `windchk`／`KL`）⇒ `s1_hdiffu_*` 9.7e-4；
- `s2_pcorr` 偏 8e-6（l16：8.912134e-10 vs CPU 8.912207e-10）；
- 而且錯法**看時序**，但同一個 binary 兩次跑出的 `sptend` 相同（stream 排程在同一台機器上大致可重現）。
（這是假設；驗證靠 l17 的 `s1_hdiffu_*`、`s2_pcorr`、`sptend` 三個量。）

`hdiffu_gpu.f90:156-169` 還有第二個問題：`!$acc exit data copyout(wmax) … async(async_id)` 之後
**沒有 wait 就在 host 上讀 `wmax(k)` 判斷 `windchk`**——這在 NVIDIA 上也是競態（上游 bug，見 §9），
只是 NVIDIA 的 copyout 通常趕得上。OpenMP 的 `exit data map(from:)` 是同步的，所以 ROCm 這邊其實沒事，
但為保險再補一條 `!$omp target update from(wmax)`，並在 CPU/GPU 兩側都印 `DBGHD wmax(1:hdk1) max=` 對照。

**修法（都在 `cmake/acc2omp.py`）**：
1. 新 pass：每一行 `NCCLCHECK(nccl(AllReduce|AllGather|Broadcast|Reduce|ReduceScatter|Bcast)(` 之後插
   `call geps_acc_wait_all()   ! layer 17: raw RCCL collective`（產生檔裡 42 個 `NCCLCHECK` 中 15 個被補：
   `mpe2d_gpu` 12、`hdiffu_gpu` 1、`rayleifr_gpu` 1、`intgrt_gpu` 1；
   其餘是 `ncclGroupStart/End`、`ncclSend/Recv`（在 `nccl_*` 包裝內，第 6 層已涵蓋）與 init/destroy）。
2. `hdiffu_gpu`：`windchk = .false.` 前插 `!$omp target update from(wmax)` + DBGHD print。
3. `src/hdiffu.f90`（CPU，探針）：同樣的 DBGHD print。**第一次放錯副程式**——檔案裡第一個是 `hdiffu_3tl`
   （`k=1..8`），CPU 路徑呼叫的是第二個 `hdiffu`（`k=1..hdk1`），已改到正確位置。

**CPU 基準（`job/TCo383L72_IC_sample_rocm_cpu.log.20260916_hd`，20:14–20:17，`PROGRAM CWBGFS HAS ENDED`）**：
6 步 `surf pres tend rms` = 0.4380739682505907 / 0.4462741042925812 / 0.42752721473604815 /
0.40398104459499434 / 0.39681622293890123 / …；`s1_hdiffu_vormid` step1 = 1.6851306945711851E-07。

⚠️ 編 CPU 樹時踩到 §6 的 reconfigure 陷阱：`src/CMakeLists.txt` 改過 ⇒ `cmake --build build_rocm_cpu` 自動 reconfigure，
沒有 `COMP_MP` 就變成 `amdflang` 直編、`mpi.mod` 找不到。修法：`export COMP_MP=… ; cmake -Bbuild_rocm_cpu -S. …-DUSE_HIP=OFF…` 再 build（§8 已補指令）。

**第 17 層結果（2026-09-16 20:18 起跑，`job/TCo383L72_IC_sample_rocm.log.20260916_2gpu_l17`，2 GPU、`GEPS_FFT_PLAN_SETS=3`、pool threshold=0）——假設成立：**

| 量 | CPU（`…_cpu.log.20260916_hd`） | GPU l17 | GPU l16（修前） |
|---|---|---|---|
| `bal` 三輪 | 1.7619e-05 / 6.3946e-06 / 1.6451e-06 | 1.7619218910e-05 / 6.3945974578e-06 / 1.6451316149e-06 | 同 |
| `DBGHD wmax(1:hdk1) max`, step 1 | （CPU 探針放錯副程式，這次沒印；已改正，待下次 CPU 跑） | 156.9612129162252（hdk1=6，兩 rank 相同） | — |
| `s1_hdiffu_vormid` | 1.6851306945711851E-07 | 1.6851306945711814E-07（2e-15） | 偏 9.7e-4 |
| `s2_pcorr` | 8.912206591E-10 | 8.912206687E-10（1e-8） | 8.912134E-10（8e-6） |
| `surf pres tend rms` step 1 | 0.4380739682505907 | **0.43807396825059497** | 0.4380603838412094 |
| 微物理（sedimentation `DO while`） | — | **仍卡死**（20:50 起 log 不再更新，GPU 100%） | 卡死 |

⚠️ 一開始誤判「微物理已通過」——log 最後一行是 step-1 末的 `DBGMAP module_mp…:4`，之後 30 分鐘沒有新輸出、GPU 100%。
用 `rocgdb -p <pid> -batch -ex "info dispatches"` 實測（不必猜）：唯一活著的 dispatch 是
`__omp_offloading_10003b_3a0e4__QMmodule_mp_gsfcgce_3ice_nuwrf_gpuPfall_flux_gpu_l625`
＝ 產生檔第 625 行 `!$omp target teams distribute parallel do collapse(3) private(t_del_tv, del_tv, notlast, …)`
（非 SL_sedi 分支的 sedimentation `DO while (notlast)`），`info threads` 有 **916 個 wave（229 個 workgroup）全部停在這個 kernel**
⇒ 不是單一格點的壞資料，是系統性的（整片 `del_tv ≤ 0` 或 NaN 之類）。`-O3` 無 debug symbol，`print del_tv` 讀不到。
回看 2026-09-16 所有 step-1 正確的 run（`mm0`、`l15`、`l16`、`l17`）**全部停在同一個地方**，與 pool threshold、第 15-17 層無關——這是下一個獨立的 bug（第 18 層）。

`job/cmp_ckpt.py` 對全部 105 個 step-1 檢查點：**除兩個探針本身的問題外全部 |GPU/CPU − 1| < 1e-9**。
那兩個例外：`pten_sl_pre/post_hadv` GPU 端印 0（GPU 的 `tendget` 在那個位置 `pten_sl` 尚未填值——探針位置不對應，
下游 `s2_*` 全部對上所以不是計算問題；**未逐行核對，列為探針疑點**）；`s2_siimpl_*`（siimpl 輸出）GPU 側沒注入。

**結論**：第 17 層是積分第 1 步 132.6 → 0.4380603 → 0.43807397 這條線的最後一片；
與第 6、7、10 層合起來，**四類 stream 競態都源自 acc2omp 丟掉 `async`**（見 memory `geps-rocm-stream-race-classes`）。
第 15（BLAS/FFT 後 sync）與第 16（純量私有化）層是無關但保守正確的修正，先留著。

#### ★★ 第 18 層（2026-09-16 21:40 起，進行中）：微物理 `fall_flux_gpu` 的 sedimentation `DO while (notlast)` 不終止

**量測（不是猜）**：rocgdb `info dispatches` 指出唯一活著的 kernel 是 `…fall_flux_gpu_l625`（非 SL_sedi 分支的 sedimentation 迴圈），
916 個 wave 全在裡面。無 debug symbol 讀不到變數，所以在 `cmake/acc2omp.py` 加探針（只對 `module_mp_gsfcgce_3ice_nuwrf_gpu`）：
迴圈計數超過 100000 就把該 column 的狀態寫進 `geps_dbgw`、強制 `notlast=.false.` 跳出，kernel 結束後印出。

**第一版探針結果（`job/TCo383L72_IC_sample_rocm.log.20260916_2gpu_l18`，22:1x）**——每個 rank、每次呼叫都有 ~91 個 thread 觸頂：

```
DBGMP stuck nhydro,i,j= 4. 265. 231.  del_tv,t_del_tv,dtb= 3.477E-243 9.494E-105 300.
  min_q,max_q,count= 3. 45. 91.  z1 min/max= 120.18 61855.9  qr max= 2.74E-04  vtr max= 12.
  rho,tz,dz8w(kts)= 4.489E-04  720509.8  112349.8
```

- 卡的是 **nhydro=4（冰）**；`del_tv ≈ 1e-243`、`t_del_tv ≈ 1e-105`、`dtb=300` ⇒ `del_tv = min(del_tv, z1/vti)` 被壓到 0，
  永遠加不到 300 s。`z1` 本身正常（120–61856 m）。
- **`rho(i,kts,j)=4.5e-4、tz(i,kts,j)=720509 K、dz8w(i,kts,j)=112350 m` 是垃圾**（地面應為 ~1.2 kg/m³、~300 K、~50 m）
  ⇒ 至少有一個「device 認為在 `i ≤ myim(j)` 內」的 column，其輸入 `rho/dz8w`（引數）與 `tz`（本常式第一個 kernel 算的）是未初始化值。
  ⚠️ 這版探針是 91 個 thread 競寫同一組欄位，(i,j) 與 rho/tz 不一定來自同一個 thread——**已改成 atomic 只讓第一個觸頂的 thread 寫、並傾印整個 column**（第二版，22:20 重編中）。
- `DBGMP pre z1 min/max= 1e30 -1e30`：host 迴圈 `do i = its, myim(j)` 一次都沒跑——因為 **`myim` 只在 device 上被賦值**
  （`mp_scheme_gpu.f90:271` 的 parallel loop，host 副本從未初始化）。探針不代表目標（§7.1），第二版先 `update from(myim)`。
  這也提醒：**任何在 host 端讀 `myim` 的程式碼在 GPU build 都是讀垃圾**（要查 `mp_scheme_gpu` 後段 `do i=1,myim(jj)` 的 host 迴圈是否存在）。

回看 2026-09-16 所有 step-1 正確的 run（`mm0`、`l15`、`l16`、`l17`）都停在同一處，與 pool threshold、第 15–17 層無關。
候選原因（**尚未驗證**）：(a) `rho/dz8w/z/th/pii` 由 `mp_scheme_gpu` 送進來時某些 column 沒填（NPEY=2 的 `jlistnum`/`my_max` 邊界？）；
(b) 第 16 層的 firstprivate 把某個「應共享」的純量私有化（要對照 `gsfcgce_3ice_nuwrf_gpu` 的 `firstprivate(...)` 清單）；
(c) 上游 `mp_scheme_gpu` 自己的 kernel 有第 5 類競態。第二版探針的整條 column 可以分辨 (a) 與其他。

**第二版探針 + 物理鏈檢查點（`job/TCo383L72_IC_sample_rocm.log.20260916_2gpu_l18b`，23:2x）——找到了：**

新增 `cmake/diabat_ckpts.py`（與 `intgrt_ckpts.py` 同格式；`inject()` 現在接受 `ck=` 並支援同一行多個 hit——順帶把
GPU 側原本漏掉的 `s2_siimpl_*` 補回來）：CPU `src/diabat.f90` 在 `do 290 jj` 前/後印 `p0_*`/`p1_*`；
GPU `diabat_gpu` 在 rrtmg/dcyc2t3/pbl_noah/gwdps/samfdeepcnv/lightning/gwdc/samfshalcnv 各呼叫之後印 `g_*_{tt,qt,ut,vt}`（used-region Σ² + 非有限值計數）：

| 檢查點 | tt nonfinite | qt nonfinite | ut nonfinite |
|---|---|---|---|
| `p0`（rrtmg 前） | 0 | 0 | 0 |
| `g_rad`、`g_dcyc` | 0 | 0 | 0 |
| **`g_gwdps`（pbl_noah_gpu + gwdps_gpu 之後）** | **2,128,788（Σ²=Inf）** | **779,470** | 0 |
| `g_deep` 以後 | 2.36M | 779k | 1.8M |

（`g_pbl` 沒印出來：注入器把它塞進 `pbl_noah_gpu` 呼叫尾巴的 `#ifdef TIMCOMCPL` 分支——注入器的續行掃描還不認 cpp 行，待修。）
所以**壞在 `pbl_noah_gpu`（或緊接的 `gwdps_gpu`）**。`pbl_noah_gpu` 的核心是 `moninedmf_gpu`，它用
**cuSPARSE `gtsvInterleavedBatch`（algo 0 = Thomas）** 解三對角系統，而且**同一組係數 `alg/adg/aug` 連解兩次**
（`moninedmf_gpu.f90:1078/1093`：先 t 再 q；`1434/1449`：先 u 再 v），兩次之間**只把上對角 `aug` 從 `aug_backup` 還原**——
也就是原作者知道 cuSPARSE 的 Thomas 會改寫 `du`，但 `dl`/`d` 不動。

**量測（scratchpad `gtsv_test.cc`，hipSPARSE 直接呼叫，m=72、batch=1000）**：

```
algo=0 (Thomas) changed: dl=0  d=72000/72000  du=71000/72000
algo=1 (LU)     changed: dl=0  d=0            du=0
algo=2 (QR)     changed: dl=70299 d=71299 du=25590
```

**rocSPARSE 的 Thomas 連主對角 `d` 也整個改寫**（rocsparse 文件也標成 `[inout]`）。所以第二次解用的是被改壞的 `adg` ⇒
v1/q1 出來 Inf/NaN ⇒ `tt`/`qt`/`phi` 壞 ⇒ `fall_flux_gpu` 的 column 裡 `z1<0`、`dz8w<0`、`qr=NaN`（第二版探針傾印的 column 就是這樣：
`z1(4)=-70.3`、`dz8w(3)=-299.7`、`qr(3)=NaN`、`del_tv=-101.7`、`t_del_tv=-1.0e7`——負的 `del_tv` 永遠加不到 `dtb=300`）。

~~修法：在 `hip_compat.cc` 保存/還原 `dl/d/du`~~ → **❌ 撤回（2026-09-17 00:05，使用者要求重新檢視根因時發現）**：
`moninedmf_gpu.f90` 裡 **所有 `cusparseDgtsvInterleavedBatch_async` 呼叫都在 `if (.false.) then` 裡（1011、1380 行）**，
執行時走的是 `else` 分支的手寫 `tridin_gpu`（t/q/tracer）與 `tridi2_gpu`（u/v），**cuSPARSE 根本沒被呼叫**。
上面關於 rocSPARSE 改寫 `d` 的量測本身是對的、`hip_compat.cc` 的保存/還原也留著（無害、若日後啟用那條路徑會需要），
但它**不是這個 hang 的原因**。犯的錯還是 §7.1 那條：拿一個「看起來能解釋一切」的機制當結論，沒先確認那段程式碼有沒有在跑
（一個 `grep -n "if (.false.)"` 就能排除）。而且它和量測也對不上：u/v 走 `tridi2_gpu`、t/q 走 `tridin_gpu`，兩者結構相同，
若是 solver 問題不會只有 tt/qt 壞而 ut/vt 好。

**目前確定的**（量測）：非有限值在 `g_dcyc`（0）與 `g_gwdps`（tt 2.1M / qt 779k）之間出現，中間只有 `pbl_noah_gpu`、`get_phi_gpu`、`gwdps_gpu`；
`ut/vt` 在同一點仍是有限的。**下一步（第三版探針，00:10 重編）**：`cmake/pbl_ckpts.py` 在 `pbl_noah_gpu` 內
每個階段（`sfc_diff/sfc_ocean/sfc_drv/sfc_sice/sfc_diag/moninedmf`）之後印 `t1/q1/u1/v1` 與地表通量 `heat/evap/stress/ustar/tg/cd/cdq` 的
used-region Σ² 與 nonfinite 計數（`b_*` 標籤），加上修好注入器後正確落點的 `g_pbl_*`，一次定位到是哪個階段、是 column 場還是地表通量先壞。
候選（未驗證）：(a) 某個 `sfc_*_gpu`（尤其 Noah `sfc_drv_gpu`，只影響陸地 column ⇒ 10% 的比例吻合）在部分 column 產生 NaN 通量；
(b) `fpvs_gpu` 用的 common block 查表 `tbpvs`（`common/fpvscom/`）在 target region 內是隱式映射——OpenMP 對 common block 成員的隱式 map 若不完整就會查到垃圾；
(c) `moninedmf_gpu` 內部 kernel 的私有化/競態。

**第三版探針結果（`job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_l18c`，00:35）——定位到 `sfc_drv_gpu`（Noah 陸面）：**

| `pbl_noah_gpu` 內的階段 | `tg`（地表溫度）nonfinite | 其他 |
|---|---|---|
| `b_in`、`b_diff`（sfc_diff）、`b_ocean` | 0 | stress/ustar/cd/cdq 正常 |
| **`b_drv`（`sfc_drv_gpu` = Noah LSM）** | **19,426**（iter 1）→ 29,470（iter 2） | 其他通量仍有限 |
| `b_sice`、`b_diag` | 29,470（沿用） | |
| `b_edmf`（moninedmf 後） | — | `t1` 2,123,388、`q1` 756,288 nonfinite（= 29,470 column × 72 層 ✓）、`hpbl` 13,376 |
| `g_pbl`（tt 更新前） | tt 0 | `g_gwdps`（tt += dtdtc·dt 後）tt 2,123,388 |

也就是：**Noah 在 ~3 萬個 column 給出 Inf 的地表溫度 `tsurf`→`tskin`→`tg`**，moninedmf 拿 `tg` 算通量把整條 column 的 t1/q1 弄成 NaN，
再經 `tt += dtdtc*dt` 進入 `tt/qt`，最後在微物理 sedimentation 變成負的 `del_tv` 而卡死。以上每一步都是量到的，不是推的。

**根因（在 `src/nvidia/sflx_gpu.f90`，讀碼確認）**：整個 Noah 是一個 9000 行的 `!$acc parallel loop gang vector collapse(2) private(…280 個名字…)` kernel。
private 清單裡有 **`snoexp`**，而 `snoexp` 是 **`data snoexp/2.0/` 初始化、迴圈裡從未賦值、只在雪蓋分支被讀**：
`t1 = tfreez*sncovr**snoexp + t12*(1 - sncovr**snoexp)`（`sflx_gpu.f90:6028`，snopac）。
OpenACC `private` 的初值是未定義的，但 **nvfortran 實作上會用外層值去起始 private 純量**（所以 NVIDIA 對）；
OpenMP `private(snoexp)` 在 ROCm 上就是真的垃圾 ⇒ `sncovr**垃圾` = Inf ⇒ 有雪的陸地 column（9 月的南極、格陵蘭、高山，~3 萬個）地表溫度 Inf。
比例、只影響陸地、只在 Noah 之後出現，三點都吻合。

**修法（`cmake/acc2omp.py`，第 16 層 pass `_privatise_assigned_scalars` 擴充）**：
1. `private(...)` 裡的純量若在整個迴圈巢狀內**從未被賦值**（沒有 `x =`、不是 DO 變數、沒出現在任何 `call` 引數裡）⇒ 從 `private` 移到 `firstprivate`
   （讀未定義值不可能是本意；firstprivate 是 nvfortran 實際給的語意）。`sflx_gpu` 被移的只有 `snoexp`、`tsnow`、`ii`（後兩個沒用到）。
2. 順便修了三個讓第 16 層幾乎失效的 bug：(a) `_DECL_RE` 的 lazy attrs 讓 `integer :: a, b`（`::` 前有空白）的名字全部漏掉、
   `real, dimension(n) :: a` 的 `a` 反而被當純量；(b) 續行的宣告（`real :: a, &`）不掃；(c) 巢狀深度不算 `do while` 與 `do 100 … 100 continue`，
   `sflx_gpu` 的 kernel 在第一個 `end do ! do_while_loop` 就被誤判結束（所以之前 sflx 一個 firstprivate 都沒加）。
   修完後第 16 層對全部 GPU 原始碼多加了約 60 條 `firstprivate`（例如 `moninedmf_gpu` 的 `accui`、`gwdps_gpu` 的 `eng0/eng1`），都是「迴圈內賦值的純量」，語意與 OpenACC 一致。
3. 靜態迴歸：104 檔翻譯 0 失敗；改變的指令逐一看過（`sflx` 16、微物理 8、`mod_ndslfv_monoadv` 7、其餘 ≤3）。

**驗證結果（`job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_l18d`，01:35）：❌ 沒修好。** `b_drv_tg` nonfinite = 19,300 → 29,304
（前一次 19,426 → 29,470：**數字每次跑都不一樣**），`t1` 2,111,832，微物理照樣卡。
所以「`snoexp` private 未初始化」是一個真實的缺陷（已修、保留），但**不是這批 Inf 的來源**——再次把「讀碼看到的合理機制」當成結論、沒等量測（§7.1）。
acc2omp 的三個解析修正與 firstprivate 擴充仍然有效且正確（NNMI `bal`、step-1 `sptend` 0.43807396825059497 都不變）。

**目前掌握的事實**：Inf 只出現在 `sfc_drv_gpu`（Noah）之後的 `tg`、只在陸地、**column 數目跑批間變動**（19,300 vs 19,426）⇒
仍然是「讀到未初始化／殘留的裝置記憶體或暫存器」這一類，只是不是 `snoexp`。候選：(a) sflx 9000 行 kernel 裡某個 private 純量在某條路徑上先讀後寫
（AMD 的 VGPR 起始值是上一個 wave 的殘留，NVIDIA 常是 0）；(b) `sfc_drv_gpu` 送進來的 `intent(inout)` 陣列（`tsurf/sneqv/snowh/sh2o/stc/cmc/ch/cm/z0`）
在某些 column 沒被填（`enter data create` 的裝置端不清零）；(c) `old_2d/old_3d` 還原邏輯。
**第四版探針（01:45 重編）**：進 column body 時把 sflx 全部 221 個 private 實數純量設成 NaN、17 個整數設成 -999999，
在 body 結尾若 `t1(i,jj)` 非有限就（atomic 只取第一個）把全部 private 純量 + 36 個 2D 輸入 + `stc/smc/sh2o(1:4)` 傾印（`DBGSFLX` 標籤）。
結尾仍是 NaN 的純量 = 這條路徑上從未被寫；輸入若是垃圾 = 問題在 `sfc_drv_gpu` 之前。一次就能分辨 (a)/(b)。

**第四版探針結果（`job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_l18e`，02:25）——兩個候選都不是，真正的原因是 runtime 的裝置堆疊：**

第一個壞掉的 column（rank 0，i=497, jj=210，myim=1416）：
- **輸入全部正常**：無雪陸地（sneqv=0、snowh=0、sncovr=0）、sfctmp 297.3 K、sfcprs 988 hPa、q2 0.0164、vegtyp 2、soiltyp 9、`smc/sh2o(1:4)=0.465`、tbot 300.2、cmc 1.4e-5、ch/cm/z0 正常。
  ⇒ 不是 (b)「`sfc_drv_gpu` 送進來的陣列沒填」。
- private 純量結尾仍是 NaN 的都是**這條路徑本來就不用的**（snopac/sfcdif/snow 那些：`bfac, esdc, t11/t12, zetalt, …`），
  沒有一個「被讀但沒寫」的 ⇒ 不是 (a)。
- **矛盾的值**：`zz1 = -Inf`，但它的公式 `zz1 = df2/(-0.5*zsoil(1)*rch*rr) + 1`（`sflx_gpu.f90:4534`）的四個輸入結尾都是有限非零
  （df2=0.262、rch=140.3、rr=1.042、`dsoil=-0.5*zsoil(1)=0.05` ⇒ zsoil(1)=-0.1）；`ddz_2 = 1/(-0.5*zsoil(2)) = Inf` ⇒ 算的時候 `zsoil(2)` 是 0；
  `stc(1:4)` 全 NaN、`ssoil/sheat/t1` NaN。也就是**private 陣列 `zsoil`（`dimension(nsoil)`，nsoil 是 dummy 引數 ⇒ 執行期大小）在 kernel 執行中被別人改掉了**。

**獨立重現（scratchpad `privarr/t2.f90`，amdflang -O3 -fopenmp --offload-arch=gfx942，1552×384 個 thread）**：
`target teams distribute parallel do private(a1..a20)`，20 個 `real(8) :: aN(nsoil)`（nsoil 為引數）每個 thread 寫自己的值、算一陣、再檢查：

| 設定 | 錯的元素數（共 47.7M） |
|---|---|
| 20 個 runtime 大小的 private 陣列，預設 | **45,146,640**、每次不同 |
| 8 個 | 16,871,410 |
| 20 個但宣告成固定大小 `a(4)` | 0 |
| runtime 大小 + `LIBOMPTARGET_STACK_SIZE=16384` / 32768 | 30.3M / 10.2M |
| runtime 大小 + **`LIBOMPTARGET_STACK_SIZE=49152`（或 65536）** | **0**，且結果與固定大小版逐位元相同 |

⇒ **執行期大小的 private 陣列放在 libomptarget 的裝置堆疊上，預設大小不夠時相鄰 lane 互相覆寫**（門檻 ~48 KB/wave；131072 以上會被 runtime 夾到 131056）。
`sflx_gpu` 的 private 清單正好有 12 個 `dimension(nsoil)` + 7 個固定大小陣列，每 lane 約 632 bytes，剛好在門檻附近 ⇒ 只有一部分 column 壞、數目每次不同。
Noah 之外，任何有 runtime 大小 private 陣列的 kernel 都可能中（微物理、`moninedmf` 等），所以這是 toolchain 層級的修正。

**修法（`job/regression.ksh`）**：`export LIBOMPTARGET_STACK_SIZE="${LIBOMPTARGET_STACK_SIZE:-65536}"`（與 `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD` 放一起）。
不改 src、不改 acc2omp；同一個 binary 直接重跑驗證（`_2gpu_l18f`，02:35 起跑）。

**驗證（`job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_l18f`，03:05）：✅ Noah 修好、微物理不再卡。**

| 量 | l18e（修前） | l18f（`LIBOMPTARGET_STACK_SIZE=65536`） |
|---|---|---|
| `b_drv_tg` nonfinite | 19,381 → 29,383 | **0 → 0**（Σ² 5.113e10，正常） |
| `b_edmf_t1` / `g_gwdps_tt` / `g_premp_tt` nonfinite | 2.1M / 2.1M / 2.3M | **0 / 0 / 0** |
| `DBGSFLX`（Noah 壞 column） | 14,747 + 15,119 | **無** |
| `DBGMP stuck`（sedimentation 卡死） | 424 / 391 / 1277 / 1308 | **無**；微物理走完 |
| NNMI `bal`、step-1 `sptend` | 不變 | 不變（1.7619218910e-05…、0.43807396825059497） |

**但緊接著在 `adjptqintp_gpu`（diabat 末尾的 pt/q 調整）兩個 rank 同時 `memory access fault`**（GPU 2 at 0x53a7d7fc000、GPU 3 at 0x5fd0bdfc000）。
`src/cwa/adjptqintp_gpu.f90:160`：`hfdpr1 = 0.5*(plold(i, k, jj) - plold(i, k - 1, jj))` 在 `do k = 1, lev` 內，`plold` 是 `(nxp, lev+1, my_max)`
⇒ k=1 時讀 `plold(i, 0, jj)`，i=1、jj=1 時讀到**陣列起點之前 8·nxp = 12,416 bytes**。`hfdpr1` 只在 `k ≥ 2` 用到，CPU/NVIDIA 上這個越界讀是無害的
（讀到前一塊記憶體）；ROCm 這邊 pool 關掉後每個緩衝區都從新的 page 開始，起點前 12,416 bytes 落在上一個 page（0x…fc000）——**兩張 GPU 的故障位址尾碼都是 `fc000`，與 base − 0x3080 落在 page 0x…fc000 完全吻合**。
上游 bug（CPU `adjptqintp.f90` 應該也有同樣的讀法，只是沒事）。
**修法（`cmake/acc2omp.py`，只對 `adjptqintp_gpu`）**：`hfdpr1 = 0.0; if (k .ge. 2) hfdpr1 = …`。03:10 重編中（改了 acc2omp ⇒ 全部重翻譯）。

#### ★★★★ 2026-09-17 05:00：ROCm GPU 首次跑完 tau=1（`job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_l18g`，`PROGRAM CWBGFS HAS ENDED`）

2 張 MI300X（NPEY=2）、binary 03:22（第 17/18 層全部修法 + 全部探針）、`LIBOMPTARGET_STACK_SIZE=65536`、pool threshold 0、`GEPS_FFT_PLAN_SETS=3`：

| 步 | tau | GPU `surf pres tend rms` | 本機 CPU（`…_cpu.log.20260916_p0`） | GPU/CPU |
|---|---|---|---|---|
| 1 | 0.167 | **0.43807396825059497** | 0.4380739682505907 | 1.000000（15 位） |
| 2 | 0.333 | 0.4284005468058377 | 0.4462741042925812 | 0.960 |
| 3 | 0.500 | 0.40145238371886277 | 0.42752721473604815 | 0.939 |
| 4 | 0.667 | 0.3866248341297478 | 0.40398104459499434 | 0.957 |
| 5 | 0.833 | 0.3785245054189599 | 0.39681622293890123 | 0.954 |
| 6 | 1.000 | 0.37468923721414743 | （CPU log 第 6 步略） | — |

- 沒有 NaN、沒有 memory fault、沒有 hang；VRAM 穩定（`free=` 一路 ~118–124 GB）。
- **每步 ~740–790 s** —— 這是探針的代價（每步幾十次 `target update from` + host 端 Σ²、DBGMAP/DBGFREE、微物理/Noah 的 NaN 播種與傾印），不是模式速度；NVIDIA 參考是 0.34 s/步。拆探針後再量。
- **第 2 步起 GPU 比 CPU 低 4–6%**：第 1 步的動力學逐位元相同，差異來自第 1 步末的物理（第 2 步 `sptend` 是第一個含物理效應的量）。
  參考：使用者的 NVIDIA log 對同機 CPU 也是逐步偏低 1–3%（0.4195/0.4122/0.4022 vs 0.4235/0.4232/0.4134），
  所以 GPU 物理與 CPU 有系統性差異本來就存在（GPU 版物理不是逐位元移植：`gpupbl/gpucup/…` 走不同的實作），但 ROCm 的 4–6% 比 NVIDIA 的 1–3% 大，還不能說已與 NVIDIA 等價。
  已有的線索：step-1 物理前後 `qt` 的 Σ² 增量 CPU +1.320、GPU +3.135（但 GPU `p1` 在 `mp_scheme` 後、CPU `p1` 在整個 290 迴圈後含 `adjptqintp`，**不可比**；
  `tt/ut/vt` 在 `p0/p1` 兩側單位不同也不可比）。
- **下一步**：兩側都加了可比的 `p2_*`（`adjptqintp` 之後 = 物理真正結束點；`cmake/diabat_ckpts.py` + `src/diabat.f90`），CPU 已重編並在跑（`…_cpu.log.20260917_p2`）；
  GPU 重編後比 `p2_qt`（同單位）與 `p2_tt`（需確認兩側此時的 tt 單位相同）。若 GPU 物理差異要追到底，路線是 §7.1 的老辦法：逐個物理常式（rrtmg → pbl → gwd → deep → shal → mp）在 CPU 側也做同樣的 used-region Σ²（CPU 的 big-loop 需改成累加後再印）。

**`p2` 比對結果（GPU `job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_p2`（binary 05:16）vs CPU `…_cpu.log.20260917_p2`，06:17）**：

`p2_qt`（`qt` 的 used-region Σ²，兩側同單位、同位置＝`adjptqintp` 之後）：

| 步 | CPU `p0_qt`（物理前） | CPU `p2_qt`（物理後） | CPU Δ | GPU `p2_qt` | GPU Δ（對 GPU `p0`） | GPU/CPU `p2` |
|---|---|---|---|---|---|---|
| 1 | 1146.2196982 | 1147.5400307 | **+1.320** | 1149.3578778 | **+3.138**（p0 = 1146.2196982，同 CPU） | 1.00158 |
| 2 | 1146.9040409 | 1148.0941234 | +1.190 | 1150.7080204 | — | 1.00228 |
| 3 | 1147.5062006 | 1148.5255797 | +1.019 | 1151.9087146 | — | 1.00295 |

- 物理前（`p0`）兩側 `qt` 相同到 8e-9；**物理後 GPU 的 `qt` Σ² 增量是 CPU 的 2.4 倍**（+3.14 vs +1.32），之後每步差距累積（0.16% → 0.23% → 0.30%）。
  這就是 step 2 起 `sptend` 低 4–6% 的來源：GPU 版物理（`gpupbl/gpucup/gpushl/gpulsp…` 走 `*_gpu` 實作）對水汽的改變與 CPU 版不同。
- GPU 的 `p0` 只在第 1 步印（`rrtmg_gpu` 不是每步呼叫），所以 GPU 側第 2 步起沒有 Δ；要每步都有，把 GPU 的 `p0` 錨點改到 `diabat_gpu` 進物理前一定執行的敘述。
- `p2_tt/ut/vt` **仍不可比**：CPU 在 `290 continue` 前已把 `tt` 轉成虛位溫（`diabat.f90:2262`），GPU 的對應轉換在 `adjptqintp_gpu` 之後的 kernel（`diabat_gpu.f90:3539`），
  `p2` 探針卡在兩者之間；`ut/vt` 兩側單位也不同（CPU 2.7e-4 vs GPU 1.96e10）。下次把 GPU `p2` 移到 3539 那個 kernel 之後即可比 `tt`。
- **決定性**：這次的 `sptend` 序列（0.43807396825059497 / 0.4284005468058377 / 0.40145238371886277）與 l18g 逐位元相同 ⇒ GPU 現在是 run-to-run 決定性的，第 18 層前的非決定性已消失。
- **這 2.4 倍是不是 ROCm 的問題，還不知道**：NVIDIA 參考 log 的 `sptend` 對其 CPU 也是系統性偏低，但沒有 NVIDIA 的 `p2_qt`。最省力的判定法：在使用者的 NVIDIA 機器上跑同一組 `p0/p2` 探針（`src/diabat.f90` 的 CPU 端已有；GPU 端要把 `cmake/diabat_ckpts.py` 的 `STATE` 探針以 OpenACC `!$acc update self` 形式放進 `src/diabat_gpu.f90`，或直接用這裡的 acc2omp 產生檔對照）。
  若 NVIDIA 也是 +3.1 ⇒ 是 GPU 物理實作本身的差異、不是移植問題；若 NVIDIA 是 +1.3 ⇒ 用 `g_*` 逐常式檢查點在 ROCm 上找哪個物理常式偏了（`g_pbl/g_gwdps/g_deep/g_light/g_gwdc/g_shal/g_premp` 的 `qt` 已在 log 裡，缺的是 CPU 側逐常式的對照）。

**逐常式對照（2026-09-17 07:30；CPU 端新增 column 累加探針 `geps_dbg_acc/_flush`（`src/tendget.f90` 末、`src/diabat.f90` 的 `do 290` 迴圈內），
與 GPU 的 `g_*` 同位置、同單位（此時兩側 `tt` 都是溫度）；CPU log `…_cpu.log.20260917_acc2`，GPU `…_2gpu_p2`，step 1）：**

| 階段 | CPU `qt` Σ² | GPU `qt` Σ² | 相對差 | CPU `tt` Σ² | GPU `tt` Σ² | 相對差 |
|---|---|---|---|---|---|---|
| 物理前 `p0` | 1146.2196982 | 1146.2196982 | 8e-11 | — | — | — |
| **PBL 趨勢套用後**（CPU `c_pbl` / GPU `g_gwdps`） | 1149.2591646466 | 1149.2591645539 | **8e-11** | 2694257430725.097 | 2694257430724.987 | **4e-14** |
| **深對流後 `qtc`**（CPU `c_deep_qtc` / GPU `g_deep_qtc`） | 1146.9588079（**−2.300**） | 1149.1794646（**−0.080**） | 1.9e-3 | 2694254950078.07（−2.5e6） | 2694155797608.85（**−1.0e8**） | 3.8e-5 |
| 淺對流後 `qtc`（CPU `c_shal_qtc` / GPU `g_shal_qtc`） | 1146.3473874（−0.611） | 1148.2952252（−0.884） | 1.7e-3 | 2694278260993.50（+2.3e7） | 2694189520476.24（+3.4e7） | 3.3e-5 |
| **微物理前**（CPU `c_premp` / GPU `g_premp`） | 1146.4312966 | 1148.3766645 | **1.7e-3** | 2694278260993.50 | 2694189520476.24 | **3.3e-5** |
| 物理後 `p2` | 1147.5400307 | 1149.3578778 | 1.6e-3 | （單位不同） | | |

- **PBL（`pbl_noah_gpu` + Noah + `moninedmf_gpu`）在 ROCm 上與 CPU 逐位元等價**（qt 8e-11、tt 4e-14）——第 18 層修完後 PBL 完全正確。
- **差異在對流段，而且定位到深對流 `samfdeepcnv_kh_gpu`**（GPU `…_2gpu_conv`，08:40，括號內是該階段的 Σ² 增量）：
  - 深對流：CPU 把 `qtc` Σ² 拉低 **−2.300**，GPU 只有 **−0.080**（差 29 倍）；但 `ttc` 的變化 GPU 反而是 CPU 的 **40 倍**（−1.0e8 vs −2.5e6）。
    「濕度幾乎不動、溫度動得很大」不是「同一套方案的實作差異」會有的樣子，比較像 GPU 版把趨勢加錯了對象／單位／層序（例如 `q` 的 tendency 加到了別的 tracer 或漏套用、`t` 用了錯的 `dt`），或是觸發判斷（`cnvflg`/`kuo`）只在少數 column 成立。
  - 淺對流：兩側同量級（CPU −0.611、GPU −0.884；`ttc` +2.3e7 vs +3.4e7；輸入狀態已不同，這種差距不足以下結論）。
  - **下一步（便宜的先做）**：在 `samfdeepcnv_kh(_gpu)` 之後兩側都印 (i) 觸發 column 數（CPU `kuo(1:nxj,jj)`／`kbot/ktop`；GPU 同名陣列）、(ii) `rcup`（對流降水）Σ²、(iii) `cldwrk` Σ²。
    觸發數相同而 `qtc` 差 29 倍 ⇒ 是趨勢套用的問題（去 `src/nvidia/samfdeepcnv_kh_gpu.f90` 找最後把 `q1/t1` 寫回的 kernel）；觸發數差很多 ⇒ 是觸發判斷（`cnvflg` 的 reduction、`private` 純量、或第 18 層那類 runtime 大小 private 陣列——先確認它有沒有 `dimension(km)` 的 private 陣列，`LIBOMPTARGET_STACK_SIZE` 現在是 65536，若不夠會出現同樣的隨機性；`sptend` 現在 run-to-run 逐位元相同，所以不像）。
  - 同時值得查：使用者 NVIDIA 參考 log 的逐步偏低是否也是這個深對流差異——若 NVIDIA 上 `g_deep_qtc` 也只有 −0.08，那是 `samfdeepcnv_kh_gpu` 本身（上游）的問題。
- 原先的概述（保留）：CPU 深＋淺對流把 `qt` Σ² 從 1149.259 拉到 1146.431（−2.83），GPU 只拉到 1148.377（−0.88）；`tt` 也在這段分岔（3.3e-5）。
  微物理＋`adjptqintp` 兩側增量相近（CPU +1.109、GPU +0.981）。
- 對流段 GPU 走 `samfdeepcnv_kh_gpu`／`samfshalcnv_kh_gpu`（`lightning_ec_gpu`、`gwdc_gpu` 夾在中間），輸出寫在 `qtc/ttc`，之後才套回 `qt`（所以 `g_deep/g_shal` 的 `qt` 沒變是正常的）。
  兩側的深／淺對流 `qtc/ttc` 檢查點都已量到（上表）：**是深對流**。
  `…_2gpu_conv` 也完整跑到 tau=1（`HAS ENDED`，第三次完整跑），6 步 `sptend` 與 l18g、p2 兩次**逐位元相同** ⇒ 三次獨立 run 決定性。

**深對流觸發數 / 降水 / 雲功比對（2026-09-17 12:10；CPU `…_cpu.log.20260917_acc3`，GPU `…_2gpu_deep`，step 1，深對流呼叫之後）：**

| 量 | CPU | GPU（ROCm） |
|---|---|---|
| 觸發 column 數（`kuo>0`） | **29,009** | **0** |
| `ktop>0` 數 / Σktop | 29,009 / 932,793 | 0 / 0 |
| `kbot>0` 數 / Σkbot | 603,648 / 42,282,775（平均 70.0） | 603,648 / 44,066,304（平均 73.0 = km+1，即全部是「無對流」預設值） |
| `rcup` Σ²（對流降水） | 1.111e-2 | **0** |
| `cldwrk` Σ²（雲功） | 5.645e9 | **0** |

⇒ **ROCm 上深對流 `samfdeepcnv_kh_gpu` 一個 column 都沒有觸發**（`cnvflg` 全被關掉），不是趨勢套用的問題。
（它仍讓 `qtc/ttc` 動了 −0.08／−1.0e8，是因為即使不觸發，該常式仍把中間轉換過的 `to/qo` 寫回；這也解釋了「濕度幾乎不動、溫度動很大」。）
`totflg` 提前返回在 GPU 版是註解掉的，所以是 `cnvflg` 沿途某個判斷把所有 column 關掉。
**下一步（12:15 重編中）**：在 `samfdeepcnv_kh_gpu` 內每個 kernel 之前 `update from(cnvflg)` 並印 `.true.` 的 column 數（`DBGCNT g_dc_<行號>`，26 個 kernel），
一次定位是哪個判斷（雲底 `kbcon` 上限、`pdot` 條件、`aa1<=0`、CIN、`xmb`…）在 ROCm 上把 29,009 變成 0；然後看那個判斷用的量是不是 reduction／private／`dimension(km)` 那類已知病。

**逐 kernel `cnvflg` 計數（`job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_cnv`，13:10，step 1，used region 共 603,648 column）：**

```
g_dc_00314/00328/00346   count=      0   （cnvflg 尚未初始化，忽略）
g_dc_00417 … g_dc_00696  count= 603648   （初始化為 .true. 之後，經 kmax/kbmax/kbm、zo、pfld/qeso/heo/heso、kb 各 kernel 都沒變）
g_dc_00966 … g_dc_02651  count=      0   ← 全部在 696–966 那個 kernel 內被關掉
g_deep_kuo               count=      0
```

⇒ **關掉所有 column 的是產生檔 696 行開始的 kernel**（`private(flgr, sigma_con, cinpcr, …)`），裡面只有兩個判斷：
(a) **無 LFC**：`kbconr == kmaxr`（在 `k ≤ kbmax` 內找不到 `heo(i,kb) > heso(i,k)`）→ `cnvflg=.false.`；
(b) **CIN**：`pfld(kb) − pfld(kbcon) > cinpcri` → `.false.`。
（第二個 kernel 內的 `dthk` 判斷在 966 之後，已經是 0 了，不相干。）
603,648 個 column 全部被 (a) 或 (b) 關掉——包括 CPU 上明明觸發的 29,009 個——所以是 (a)/(b) 用到的**輸入**系統性地錯，而不是判斷本身：
`heo/heso`（= `phil + cp·to + hvap·qo/qeso`，GPU 版用 `phil` 而 CPU 版用 `0.5·g·(zo(k)+zo(k+1))`）、`kb/kbmax/kmax`（由 `prsl/ps` 決定）、`pfld`、或 `gdx/sigma`→`cinpcri`。
候選（未驗證）：`phil` 的單位／方向（若比 CPU 大很多，`heo(kb) − heso(k)` 的 `phil(kb) − phil(k)` 項會壓過一切 ⇒ 全域無 LFC）；`fpvs_gpu` 查表（同一個 common block 在 Noah 裡是對的，機率低）；`kbmax` 太小。
**第三版探針（13:25 重編中，`…_2gpu_lfc`）**：在 LFC 判斷之後把 rank 0 的 column (i=100, jj=100) 傾印：`kmax/kb/kbmax/kbcon/heo(kb)/kbm/ps/gdx` 與整條 `heo/heso/phil/to/qo/pfld` 剖面（`DBGLFC`），
用物理量級直接判斷是哪個輸入壞了（`phil` 應為 0→~5e5 m²/s² 隨高度增、`to` 200–300 K、`qo` ≤ 0.02、`pfld` 單調遞減的 mb）。
第一次（`…_2gpu_lfc`，跑完 tau=1）**沒印出來**：print 用 `if (myrank .eq. 0)` 守門，而 `samfdeepcnv_kh_gpu` 的 `myrank` 是**呼叫端沒傳的 dummy 引數**——
`diabat_gpu.f90` 傳 29 個實引數、常式宣告 34 個（多出 `myrank, icheck, kcheck, jjcheck, rcheck`），沒有 explicit interface 所以編譯器不擋，
那 5 個 dummy 是垃圾（常式內只有註解掉的 print 用到，對計算無影響；NVIDIA 同源碼也一樣，列入上游回報）。改用模組的 `geps_hide_myrank` 守門，15:15 重編重跑（`…_2gpu_lfc2`）。

**剖面結果（`…_2gpu_lfc2`，17:00，rank 0 column (100,100)）——找到了：**

```
DBGLFC kmax,kb,kbmax,kbcon = 54, 21, 30, 54    ← kbcon == kmax ⇒「無 LFC」
 k   heo          heso   phil        to        qo         pfld
 1   2.888E+05    NaN    2.105E+02   276.27    4.44E-03   1002.4
 2   2.887E+05    NaN    6.561E+02   275.81    4.34E-03    996.8
 …（k ≤ kmax 全部 heso = NaN；heo/phil/to/qo/pfld 都是正常量級）
```

`heo` 正常、`heso` 整條 NaN ⇒ `qesor`（飽和比濕）壞。追到產生檔 607–680 行的 kernel（介面層的 `to/qo/heo/heso` 重算）：
`qesor = max(qesor, val1)`、`qor1 = max(qor1, val2)`，而 **`val1`/`val2` 是在 kernel 之前於 host 上設的（`val1=1e-8; val2=1e-10`），
原 OpenACC 指令卻把它們列在 `private(...)`** ⇒ OpenMP 上 kernel 內讀到的是未初始化的 private 值（跟 `snoexp` 一模一樣的病，nvfortran 會用外層值起始所以 NVIDIA 沒事）。
垃圾 `val1` 若很大 ⇒ `qesor` 巨大 ⇒ `heso` Inf/NaN ⇒ `heo(kb) > heso(k)` 永遠不成立 ⇒ **全部 column 無 LFC ⇒ 深對流 0 觸發**。全部串起來了。

**為什麼第 18 層的「private 但從未賦值 ⇒ firstprivate」pass 沒抓到**：`val1/val2` 宣告在 `real(kind=kind_phys) adw, aup, … &`（**沒有 `::` 的舊式續行宣告**，中間還夾一行註解掉的續行），
pass 只對含 `::` 的續行宣告做接合，所以這些名字不在「純量」清單裡。修：續行宣告不論有無 `::` 都接合（`cmake/acc2omp.py`）。
靜態迴歸（104 檔）新增/改動的 `firstprivate`：`moninedmf_gpu` 11、`samfshalcnv_kh_gpu` 8、`samfdeepcnv_kh_gpu` 5（含 `val1, val2`）、`gwdps_gpu` 4、`nor_gwdp_gpu` 3、`ozphys_2015_gpu` 2、`mfpbl_gpu` 2、`diabat_gpu` 2——
都是「迴圈內賦值的純量」或「private 但從未賦值」兩類，語意與 OpenACC 一致。**17:10 全部重翻譯重編中（`…_2gpu_val`）**，驗證：`g_deep_kuo` 應 ≈ 29,009、`c_deep_qtc` vs `g_deep_qtc`、step-2 `sptend` 應向 CPU 的 0.44627 靠近。

附帶：這一類（OpenACC `private` 純量在 kernel 外初始化、kernel 內只讀）在 `src/nvidia`、`src/cwa` 裡目前找到 `snoexp`、`val1/val2`、`tor`、`tsnow`、`ii` 等；
它們在 NVIDIA 上能跑純屬 nvfortran 的實作行為，上游應改成 `firstprivate` 或直接拿掉 `private`。

#### ★★★★★ 2026-09-17 18:55：ROCm GPU 與 CPU 六步全部對上（`job/TCo383L72_IC_sample_rocm.log.20260917_2gpu_val`，binary 17:19，`PROGRAM CWBGFS HAS ENDED`）

| 步 | CPU `surf pres tend rms`（`…_cpu.log.20260917_acc3`） | GPU（ROCm，2×MI300X） | 相對差 |
|---|---|---|---|
| 1 | 0.4380739682505907 | 0.43807396825059497 | 1e-15 |
| 2 | 0.4462741042925812 | 0.4462741042925732 | 2e-14 |
| 3 | 0.42752721473604815 | 0.4275272147313331 | 1e-11 |
| 4 | 0.40398104459499434 | 0.403981044588656 | 2e-11 |
| 5 | 0.39681622293890123 | 0.3968162229541496 | 4e-11 |
| 6 | 0.39255042411024554 | 0.39255042414676805 | 9e-11 |

物理逐常式（step 1）：深對流觸發 **29,009 = 29,009**、`g_deep_qtc` 6e-12、`g_deep_rcup` 1e-13、`g_shal_qtc` 6e-12、`p2_qt`（物理結束）1147.54003062 vs 1147.54003071（8e-11）；
`DBGLFC` 的 `heso` 全部有限（2.899e5）、`kbcon=22`。六步的 `p2_qt` 序列 GPU/CPU 都在 1e-10 內。
**結論：ROCm 移植的動力與物理現在都與 CPU 等價到 1e-10 量級**（比使用者的 NVIDIA 參考 log 對 CPU 的 1–3% 還近——那 1–3% 很可能就是 `val1/val2`/`snoexp` 這類在 NVIDIA 上「碰巧能跑」的未初始化 private；值得在 NVIDIA 上用同樣的 firstprivate 修法驗證）。

修法總表（全部在 toolchain，`src/` 只有探針）：第 6/7/10/17 層 stream 競態等待 → 第 12 層 `map(release)` → 第 13/14 層 FFT plan LRU、pool threshold 0 → 第 15 層 BLAS/FFT 同步 →
第 16/18 層純量私有化（賦值 ⇒ firstprivate；private 未賦值 ⇒ firstprivate；宣告解析修三處＋`do while`/標號 DO）→ 第 18 層 `LIBOMPTARGET_STACK_SIZE=65536`、`adjptqintp` 越界讀 guard。

**剩下的**：拆探針做乾淨對照＋量速度（現在 ~750 s/步全是探針）、1-GPU 路徑、上游回報清單（§9 8d）。

**探針開關（2026-09-18 00:10）**：`cmake/acc2omp.py` 所有探針 pass（DBGMAP/DBGFREE/DBGTRM/DBGHD/DBGMP/DBGSFLX/DBGLFC/cnvflg 計數、intgrt/diabat/pbl 檢查點表、nnmi/tendget/tranrs 值探針）
改成 **opt-in：`GEPS_ACC2OMP_PROBES=1`**，預設關；`src/rocm/hip_compat.cc` 的 ndjson VRAM 探針（cyclic_cell 每次平流呼叫 34 次 `hipMemGetInfo`+檔案追加）改成 `GEPS_DBG_NDJSON=1` 才開。
修法（競態等待、firstprivate、release、guard、FFT plan LRU、BLAS/FFT 同步）不受開關影響。`src/` 端的 CPU 探針（`intgrt/tendget/tranrs/hdiffu/diabat.f90` 的 `! #region agent log`）只進 CPU 執行檔，暫留。
乾淨翻譯後 104 檔產生碼裡沒有任何 `DBG`/`geps_dbg` 呼叫。

**乾淨 2-GPU 跑結果（`job/TCo383L72_IC_sample_rocm.log.20260918_2gpu_clean`，binary 00:15，01:51 `HAS ENDED`）**：
6 步 `sptend` 與帶探針的 `_2gpu_val` **逐位元相同**（0.43807396825059497 / 0.4462741042925732 / 0.4275272147313331 / 0.403981044588656 / 0.3968162229541496 / 0.39255042414676805）⇒ 探針對數值零影響。
**但每步仍要 715–770 s**（`Timing=` 768/715/770/740/740）——之前「~750 s 全是探針」的判斷是錯的，探針幾乎不花時間。

真正的原因在 log 裡（`[rfftmlt_loop rocm]` 的 LRU 訊息）：模式每步用到 **7 個不同的 FFT plan set**（`(jump,m,isign)` = 1554 × {1,2,72,144} × ±1），
而 `GEPS_FFT_PLAN_SETS=3` 讓 LRU **每步重建 44 個 set × 384 個 rocFFT plan ≈ 17k 次 plan 建立**（rocFFT 每個 plan 幾十 ms ⇒ 幾百秒）。
每個 set ≈ 384 × ~18 MB ≈ 7 GB；7 個 set ≈ 49 GB/GPU，2 張卡各有 ~125 GB 空閒，放得下。
**修：`src/rocm/rfftmlt_loop_gpu.f90` 預設 budget 4 → 8**（`GEPS_FFT_PLAN_SETS` 仍可覆蓋）；02:20 用同一個乾淨 binary 跑 `GEPS_FFT_PLAN_SETS=8`（`…_2gpu_fft8`）量每步秒數。
之後若還慢，下一個嫌疑依序是：pool threshold 0（每次 enter data 都 `hipMalloc/hipFree`）、第 15 層在每次 dgemm/FFT 後的 `hipDeviceSynchronize`、第 10/17 層的 `geps_acc_wait_all()`（整卡同步）。要量，不要猜：`rocprofv3 --kernel-trace --stats` 或在 `intgrt_gpu` 主要階段前後印 `omp_get_wtime()`。

**`GEPS_FFT_PLAN_SETS=8` 結果（`…_2gpu_fft8`，03:12，`HAS ENDED`）**：plan set 重建 0 次；每步 **514 / 502 / 502 / 496 / 497 s**（原 715–770）；6 步 `sptend` 不變。
⇒ FFT plan 重建佔 ~250 s/步，**還有 ~500 s/步在別處**（NVIDIA 0.34 s）。不再猜：03:15 起用 `rocprofv3 --kernel-trace --memory-copy-trace --hip-runtime-trace --stats` 跑完整一次（`…_2gpu_prof`）。
踩到：rocprofv3 的 OMPT 工具在 OpenBLAS 初始化（`gotoblas_init → omp_get_num_places → ompt_post_init → omp_get_num_devices`）時 segfault，**`OMP_TOOL=disabled` 繞過**；
`job/regression.ksh` 新增 `GEPS_PROF_WRAP` 鉤子（放在執行檔前面的 profiler 前綴）。

**rocprofv3 結果（`…_2gpu_prof`，04:43，整個 run：NNMI + 6 步，rank 0；統計檔在 scratchpad `prof/r_3982478_*_stats.csv`）**：

| 項目 | 次數 | 總時間 | 說明 |
|---|---|---|---|
| **`hipMemcpyAsync`** | 9,710 | **1,248 s**（avg 128 ms、max 88.6 s） | `intgrt_gpu` 每步 14 個 D2D 狀態複製（`up=ut, um=ut, …, qm=qt`）。從 trace 看每步：3×3 s + 7×13 s + **1×87 s**（qm）——**~26 MB/s**，不是 D2D 的速度。shim 用 `hipMemcpyDefault`，libomptarget 配的緩衝區不被 HIP 的指標分類認成 device 記憶體 ⇒ CLR 走 host-staged 路徑。 |
| `hipModuleLoadData` | 25,753 | 492 s | rocFFT plan 建立時載入 code object；budget=8 之後只在暖機期（NNMI）發生，每步 0。 |
| kernel `ndslfv_monoadvv_gpu_l440` / `_l471`（垂直平流的 pack/unpack） | 84,888 / 84,888 | 372 s / 262 s（avg 4.4 / 3.1 ms） | `src/rocm/ndslfv_monoadvv_gpu.f90` 的 tile **`otile=32`**（OOM 時代的值）⇒ 每次呼叫 ~9,400 個 tile × ~8 次 launch，而 pack kernel 每個 tile 都掃過**全部** column（O(nxptot²/otile)）。 |
| `hipLaunchKernel` | 462,313 | 149 s（avg 323 µs！） | launch 本身平均 0.3 ms——大半是排在上面那些 tile kernel 後面等 queue。 |
| `hipDeviceSynchronize` / `hipStreamSynchronize` | 351k / 311k | 52 s / 45 s | 第 10/15/17 層的等待；每次 ~150 µs，總量可接受，之後再細修。 |
| memcpy D2D（真正的複製） | 3.04M | 41.6 s（avg 13.7 µs） | 正常。 |
| `hipMemsetAsync`/`D32Async` | 681k | 8 s | 正常。 |

⇒ 每步 ~500 s 的組成（估）：**狀態複製 ~190 s + 垂直平流 tiling ~70–100 s + launch/等待 ~40 s + 其他**。

**修法（05:00）**：
1. `src/rocm/hip_compat.cc` `geps_hip_memcpy_async`：kind=DeviceToDevice 時改用 **`hipMemcpyDtoDAsync`**（不做指標分類）。重編 05:05，跑 `…_2gpu_dtod`。
2. `src/rocm/ndslfv_monoadvv_gpu.f90` 與 `vertical_cell_advect_gpu.f90` 的 **`otile` 32 → 8192**（temp 每個 < 100 MB；pool 關掉且 2 GPU 後記憶體充裕）。⚠️ tile 大小會改變 `def_cfl_step_gpu_type2` 取 nstep 的範圍（per-tile max），
   要用 6 步 `sptend` 確認數值不變。下一次 build 生效（`…_2gpu_tile`）。
3. 之後：pack kernel 改成 tile-local（用 column→(i,j) 對照表，去掉 O(nxptot) 掃描）；第 15 層每次 BLAS/FFT 後的 `hipDeviceSynchronize` 改成 stream 等待；kernel launch 數（46 萬/run）本身也要壓。

**`…_2gpu_dtod` 結果（06:10）**：每步 504/502/500/544/544 s——**DtoD 沒有改善**，`sptend` 不變。
⇒ `hipMemcpyAsync` 的 128 ms/88 s 不是複製本身慢，而是**API 在等前面排隊的 GPU 工作**（host 在這裡阻塞，時間其實屬於前一階段的 kernel）。
對照 kernel 統計：真正的時間在垂直平流 tiling 的 kernel（634 s/run）與 46 萬次 launch 的排隊（`hipLaunchKernel` avg 323 µs 也是同一種阻塞）。
DtoD 改法保留（正確且無害）。**06:15 重編 `otile=8192`（`…_2gpu_tile`）。**

**`…_2gpu_tile` 結果（07:01，`otile=8192`）**：每步 **260 / 247 / 247 / 247 / 247 s**（原 500）；6 步 `sptend` **逐位元不變**（per-tile nstep 沒改變結果）。
每步累計進度：740（探針/FFT 重建）→ 500（`GEPS_FFT_PLAN_SETS=8`）→ **247**（`otile=8192`）。距 NVIDIA 0.34 s 還差 ~700 倍。07:05 再 profile 一次（`…_2gpu_prof2`，統計在 scratchpad `prof2/`）找下一個。

**prof2（07:55，`otile=8192` 的 build）**：kernel 總時間降到 432 s/run（~30 s/步；最大的是 `vertical_cell_ppm_intp` 60 s、`def_cfl_step` 58 s、pack/unpack 53+39 s、`cyclic_cell_ppm_intp_two_loops` 4 個 kernel ~105 s），
但 HIP API 仍 1,770 s：**`hipMemcpyDtoDAsync` 84 次 = 1,127 s（每次 13.4 s）**。這次 kernel 很少了，所以不再是「等 queue」——**D2D 複製本身就是 ~26 MB/s**。
機制：libomptarget 用 HSA 直接配裝置記憶體，HIP 的 runtime 不認識這些指標，`hipMemcpy*` 退到 **CPU 走 large-BAR 映射逐字讀寫裝置記憶體**（uncached，正是幾十 MB/s 的量級）；kernel 看到的資料是對的，所以數值沒事。
**修：`geps_hip_memcpy_async` 的 D2D 改用 `omp_target_memcpy`（libomptarget 自己的 `hsa_amd_memory_async_copy`）**，08:05 重編跑 `…_2gpu_ompcpy`。預期每步 247 → ~60 s。
之後的順序：`hipModuleLoadData` 461 s/run 只在暖機（rocFFT 25k 次載入，可用 rocFFT plan cache 或減少 plan 數）；剩下的 `hipMemcpyAsync` 9.6k 次 102 s（來源待查，可能是 rocFFT/NCCL 內部）；垂直平流 kernel 本身（30 s/步）要靠 tile-local pack 與更大的 tile。

**`…_2gpu_ompcpy` 結果（08:13）**：每步 **69 / 58 / 58 / 58 / 58 s**（原 247）；6 步 `sptend` 逐位元不變。
累計：740 → 500（FFT plan set 8）→ 247（`otile=8192`）→ **58 s**（D2D 走 `omp_target_memcpy`）。08:15 再 profile（`…_2gpu_prof3`）。

**prof3（08:55）**：kernel 413 s/run（每步 ~40 s：`vertical_cell_ppm_intp` 60、`def_cfl_step` 60、pack/unpack 53+39、`cyclic_cell_ppm_intp_two_loops` 四個 kernel ~105、fgnl 15 s，都是 9 次呼叫的總和）；
HIP API：`hipModuleLoadData` 465 s（只在暖機）、**`hipMemcpyAsync` 9,626 次 102 s，其中 48 次 ≥0.1 s 共 94 s**（每步一組 0.2/3.2/9.8/9.8/3.2/3.2 s），又是 26 MB/s 那條 CPU 路徑，但不是我們的 shim 發的（前面緊接 `hipThreadExchangeStreamCaptureMode`）。

**LD_PRELOAD backtrace（`…_2gpu_bt`，scratchpad `interpose/libmemcpy_bt.so` 攔截 ≥16 MB 的 `hipMemcpyAsync` 印 `backtrace()`，09:00）**：

```
[memcpy_bt] hipMemcpyAsync 81 MB kind=4 (hipMemcpyDefault)
   librccl.so.1(+0x22053e78) … (+0x22011eaa)          ← RCCL 內部
   tcogfs.x [0x17320b9]  = _QMncclPncclallgather        ← src/rocm/nccl.f90 的 ncclAllGather 包裝
   tcogfs.x [0x116b3be]  = mpe2d_unify_lev_gpu_
   tcogfs.x [0x120b38e]  = tendget_gpu_
   tcogfs.x [0x110e447]  = initial_gpu_
```

⇒ **RCCL 的 `ncclAllGather` 在 out-of-place 時用 `hipMemcpyAsync(hipMemcpyDefault)` 把本地那一塊複製到 recvbuff**，HIP 不認識 libomptarget 的緩衝區 ⇒ CPU/large-BAR 路徑 ⇒ 每步 ~30 s。
只有這一個 call chain（其他集合通訊沒被抓到 ≥16 MB 的複製）。
**修（`src/rocm/hip_compat.cc`）**：`geps_nccl_allgather` 改成 **in-place**——先用 `omp_target_memcpy` 把 sendbuff 搬到 `recvbuff + rank*count`，再以該位址當 sendbuff 呼叫（NCCL 規定的 in-place 形式，RCCL 就不做本地複製）；
`geps_nccl_broadcast` 同理（root 先 `omp_target_memcpy` s→r，所有 rank 用 s=r）。09:05 重編、跑 `…_2gpu_inplace`。預期每步 58 → ~40 s；之後剩下的就是 kernel 本身（垂直平流 + cyclic_cell，各 ~15–25 s/步）與 launch 數。

**`…_2gpu_inplace` 結果（09:23）**：每步 68/57/57/57/57 s——**幾乎沒變**（那些 AllGather 大複製主要在 NNMI 的 `initial_gpu/tendget_gpu`，預報步裡少）；`sptend` 不變。修法保留。

**prof3 的每步 kernel 分解（trace 最後 60 s = 一步）**：46,638 次 launch、GPU 忙 44.1 s：
`def_cfl_step_gpu_type2` 7.95 s（74 次 × 107 ms）、`vertical_cell_ppm_intp l2100` 7.5 s（111 × 68 ms）、pack `l440` 5.9 s（37 × 159 ms）、unpack `l471` 4.4 s、
`cyclic_cell_ppm_intp_two_loops` 四個 kernel ~10 s（各 573 次）、fgnl pack/unpack 2.6 s、其餘 < 1 s；另外 57 − 44 ≈ 13 s 是 host/launch 開銷。
三個最大的都是**「`!$acc kernels` 被翻成單執行緒 `!$omp target`」或「每 tile 掃全部 column」**的問題，不是演算法：
- `def_cfl_step_gpu_type2`：`kernels` + `loop` → 序列 target（8192×72 單執行緒）→ acc2omp 改成 `teams distribute parallel do reduction(max:check_max)`（去掉 atomic）。
- `vertical_cell_ppm_intp_gpu` 的 `hh(k,i) = pp(k+1,i) − pp(k,i)`：同病 → `collapse(2)` 平行。
- pack/unpack（`src/rocm/ndslfv_monoadvv_gpu.f90` 兩個常式）：改成 tile-local——先建 column→(i,j) 對照表 `col_i/col_j(nxptot)`（每次呼叫一個 kernel），pack/unpack 直接 `do ot=1,nloc; do k` collapse(2)。
`src/` 裡共 30 個 `!$acc kernels` 區域（大多是單一陣列指定或幾條純量），只有這兩個在熱路徑上；其餘先不動（一般性的 kernels→parallel 轉換有相依性風險）。
09:35 重編、跑 `…_2gpu_pack`。預期每步 57 → ~35 s；再來就是 cyclic_cell 的四個 kernel（~10 s）、launch 數與同步。

（第一次 `…_2gpu_pack` 是舊 binary：acc2omp 的檔名守門 `"ndslfv_monoadv" in name` 也命中 `ndslfv_monoadvh_gpu.f90`，assert 讓翻譯失敗、build 沒產生新執行檔——job 腳本照跑舊的。改成 `mod_ndslfv_monoadv_gpu`。**教訓：跑之前看 `bin/tcogfs.x` 的時間戳。**）

**`…_2gpu_pack` 結果（binary 09:57，10:40）**：每步 **39 / 28 / 28 / 28 / 28 s**（原 57）；6 步 `sptend` 逐位元不變。
累計：740 → 500 → 247 → 58 → **28 s/步**。13:10 再 profile（`…_2gpu_prof4`）。

**prof4（13:50，每步 ~30 s 的分解，trace 最後 30 s）**：46,726 次 launch、**GPU 只忙 15.6 s**，其中 `cyclic_cell_ppm_intp_two_loops_gpu` 6 個 kernel **11.2 s**（每個 574 次 = 逐緯度 launch，l253 每次 6.9 ms）、`vertical_cell_advect` pack 1.0 s、其他 < 1 s；
launch 數：rocFFT 小 kernel ~25k（384 個逐緯度 plan × bluestein/pre/post）、OpenMP kernel ~7.5k、dgemm 1.5k。
**GPU 閒置 11.9 s，散在 2,608 個 > 2 ms 的空檔**（平均 4.5 ms、最大 30 ms，發生在物理各常式之間、kernel 與 kernel 之間，空檔內幾乎沒有 HIP API 呼叫）⇒ 是 host 端每個 target region 的 libomptarget 開銷（映射查表／隱式 map 的 `hipMalloc`+H2D+`hipFree`——threshold=0 時每次都真的配/釋放）或 host 端 Fortran/MPI。
**兩條路**：(a) 便宜的實驗：`LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=16777216`（小配置走 pool、大塊照樣釋放；當初設 0 是為了 VRAM 與遮罩競態，競態已修）→ `…_2gpu_mm16` 跑中；
(b) 結構性：`src/rocm/cyclic_cell_ppm_gpu.f90` 的逐緯度 kernel 改成整個 domain 一次 launch（574×6 → 6），並把 rocFFT 逐緯度 plan 換成少數幾個 batched plan（同長度的緯度合併）。

**`…_2gpu_mm16` 結果（13:37）**：每步 38/27/27/27/27 s——**只快 1 s**，`sptend` 不變 ⇒ 空檔不是 `hipMalloc/hipFree`。threshold 先維持 0（VRAM 保守）。
下一個量測：沒有 `perf`，用 rocgdb 對 rank 0 在預報階段取 40 次 host backtrace（scratchpad `sample_bt.sh` → `host_samples.txt`），看 host 在空檔裡卡在哪（libomptarget 映射／MPI 等待／host Fortran）。跑 `…_2gpu_bt2`。

**host backtrace 直方圖（`…_2gpu_bt2`，rank 0，預報第 2–4 步，rocgdb 取樣 30 次，scratchpad `host_samples2.txt`；第一版腳本 `pgrep -f` 抓到 mpiexec 而不是 tcogfs.x，改用 PID）**：

| 樣本 | host 在哪 |
|---|---|
| 29/30 | `libhsa-runtime64 → AMDGPUStreamTy::synchronize → __tgt_target_kernel`（**等 kernel 跑完**） |
| 1/30 | `DeviceTy::allocData → mmap64`（一次配置） |

等待的 kernel 全部在 `cyclic_cell_ppm_intp_two_loops_gpu`（`src/rocm/cyclic_cell_ppm_gpu.f90`），呼叫者分布：`cyclic_cell_intpx_jlist_gpu` 14、`cyclic_cell_massadvx_jlist_gpu` 8、`cyclic_cell_massadvy_mylonlen_gpu` 6（+1 在 massadvy 本身）；再上去都是 `ndslfv_monoadvh2_xy_gpu_refactor → ndslfv_monoadvh_gpu_refactor → intgrt_gpu`。
⇒ **每步 ~28 s 幾乎全部是水平平流（NDSL 半拉格朗日的 PPM 內插）**，host 沒有在做別的事，也不是 MPI 或 libomptarget 映射；prof4 看到的「GPU 閒置空檔」是這些 kernel 之間的 launch 延遲/尾巴（kernel 本身占用率極低）。
根因：這三個 src/rocm 檔（`cyclic_cell_ppm_gpu / intpx_gpu / massadvy_gpu`）都用 **`otile = 8`**（OOM 時代）：每個 kernel 只有 8 列 × 72 層 = **576 個 team**，i 方向（1552 個經度）在 team 內序列做 ⇒ 一次 launch 用不到 MI300X 的 1/50；每次呼叫 48 個 tile × 6 個 kernel、每步 3.4k 次 launch。
**修：三個檔的 `otile` 8 → 64**（4608 team/launch、launch 數 ÷8；PPM 暫存區 `dqmono/qi(3·lonn, nvars, 72, 64)` 各 ~1.5 GB，總計 ~4 GB，2 GPU 放得下）。14:15 重編、跑 `…_2gpu_ot64`。
之後若要再快：把 `n`（tracer）維度也 collapse 進去、或 otile 直接 = outer_size（記憶體 ~25 GB）；rocFFT 逐緯度 plan（每步 25k 個 3 µs 的小 kernel）合併成 batched plan。

**`…_2gpu_ot64` 結果（14:12）**：每步 **23 / 11.7 / 11.8 / 11.8 / 11.8 s**（原 28）；6 步 `sptend` 逐位元不變。累計 740 → 500 → 247 → 58 → 28 → **11.8 s/步**。
14:15 試 `otile = 384`（一個 tile 涵蓋全部列；暫存 ~25 GB）→ `…_2gpu_ot384`。
**`…_2gpu_ot384` 結果（14:36）**：每步 **15.5 s——比 64 慢**（每次呼叫要配/釋放 ~25 GB 的暫存，threshold=0 時是真的 hipMalloc/hipFree；或快取效應），`sptend` 不變。**保留 otile = 64**（已改回）。
14:40 在 11.8 s 的 build 上再 profile（`…_2gpu_prof5`）。

**prof5（15:00，每步 ~12.5 s 的分解）**：40.7k 次 launch，**GPU 只忙 6.0 s**（kernel 之間無重疊）：
`vertical_cell_advect` pack 1.0 s（74 × 14 ms，還是逐緯度掃描的那個）、`trngra` 0.7 s（2 × 357 ms）、`cyclic_cell_ppm` 6 個 kernel 共 1.8 s（每個 61 次，avg 2–10 ms）、RCCL 0.43 s、`massadvx` 0.2、`hdiffu` 0.35、其餘 < 0.15 s。
launch 組成：rocFFT 小 kernel ~25k（bluestein/pre/post，逐緯度 plan）、dgemm 1.5k、OpenMP kernel ~5k。host 端：`hipDeviceSynchronize` 16.5k 次 1.33 s（第 15 層每次 FFT/BLAS 後一次 + 第 10/17 層）、`hipModuleLaunchKernel` 32k 次 1.08 s。
⇒ 剩下 ~6.5 s 是 host/launch 開銷：**逐緯度 FFT（每步 ~3k 次 `hipfftExec`，每次 3–6 個 3 µs 的 kernel + 一次 device sync）** 與 OpenMP 每次 launch 的 libomptarget 成本。再取一次 host backtrace 直方圖（`…_2gpu_bt3`）確認比例後，下一個結構性修法：**把同長度的緯度合併成 batched rocFFT plan**（384 個 plan → ~20 個：distinct nxj 數）並移除逐次 sync（改成整批後一次 wait）。

**host backtrace 直方圖 #2（`…_2gpu_bt3`，11.8 s/步的 build，rank 0，40 個樣本，scratchpad `host_samples3.txt`）**：

| 類別 | 樣本 | 說明 |
|---|---|---|
| 等 GPU（`AsyncInfoTy::synchronize` → HSA signal wait） | **27/40** | 其中 `cyclic_cell_ppm_intp` 14（massadvy 7 / intpx 4 / massadvx 3）、`vertical_cell_advect` pack 4、`geps_fft_exec_z2d`（第 15 層每次 FFT 後的 `hipDeviceSynchronize`）4、`hdiffu` 1、dgemm 1 |
| HIP runtime 內部 | 4 | launch 路徑 |
| rocFFT host 端 | 3 | plan 執行的 host 開銷 |
| libomptarget 映射（`deleteData/submitData`） | 4 | `enter/exit data` 的查表與釋放 |
| `ioctl`/其他 | 2 | HSA queue |

對照 prof5 的每步 HIP API 計數：`hipModuleLaunchKernel` 19.9k（rocFFT）、`hipExtModuleLaunchKernel` 2.7k（hipBLAS）、`hipDeviceSynchronize` 9.8k、`hipStreamSynchronize` 4.1k、`hipSetDevice/GetDevice` 33k（rocFFT/hipBLAS 每次呼叫都做）、`hipMemsetAsync` 只有 247（45 萬次那是 NNMI）。
**解讀**：host 大部分時間在等 kernel，但 kernel 實際只跑 6 s ⇒ 差額是 **launch→開始的延遲乘上 4 萬次 launch**（每次幾十到一百多 µs），而不是某個 kernel 慢。所以方向是「減少 launch 數與同步數」：
- 逐緯度 FFT：因為 reduced grid 每條緯度 `nxj` 都不同（同一個 rank 內沒有兩條相同），**不能**合併成 batched plan（先前的想法作廢，除非補零改變長度——那會改變結果）。能做的是拿掉第 15 層每次 exec 後的 `hipDeviceSynchronize`，整個 384 次迴圈後一次 `geps_acc_wait(1)`（rocFFT 的 kernel 會在 stream 上排隊流水）。
- `vertical_cell_advect_gpu` 自己的 pack kernel（第 73 行那個逐緯度掃描，1.0 s + 等待）改成 tile-local（同 wrapper 的 `col_i/col_j` 做法）。
- `trngra_gpu l67`（2 × 357 ms）、`hdiffu_gpu`（2 × 90 ms）之後看。
- OpenMP 每次 launch 的 libomptarget 固定成本（~100–200 µs）×5k 只能靠合併 kernel（例如 cyclic_cell 的 6 個 kernel、pack/unpack 的多個 kernel）。
15:30 做前兩項，跑 `…_2gpu_fftw`。

**`…_2gpu_fftw` 結果（15:39）**：每步 **11.8 s——沒變**（FFT 逐次 sync 與 vertical pack 都不是瓶頸），`sptend` 不變。兩個改法保留（無害）。
推論再次失準 ⇒ 改成直接量：acc2omp 加 **`GEPS_ACC2OMP_TIMERS=1`**（build 時的環境變數）——在 `intgrt_gpu` 每個 host 層級的 `call xxx(...)` 前後取 `system_clock`，按被呼叫者累加，
在 `surf pres tend rms` 那行一起印 `TIMER <callee> <秒>`（`src/rocm/geps_acc_wait.f90` 加了 `geps_wtime`）。順便修了一個優先序 bug：`if PROBES and "intgrt_gpu" in name or "tendget_gpu" in name` 讓 tendget 的探針在「乾淨」build 裡其實還在（沒印東西但有 update from）。15:45 重編、跑 `…_2gpu_timer`。
（第一次 timer build 失敗：續行呼叫中夾了註解行，計時碼插到引數中間；injector 已改成跳過續行內的註解。）

**TIMER 結果（`…_2gpu_timer`，binary 16:05，16:15；rank 0，第 3 步 = 11.8 s；兩個 rank 都印所以 log 裡每項出現兩次）**：

| `intgrt_gpu` 內的呼叫 | 秒/步 | 占比 |
|---|---|---|
| `ndslfv_monoadvh_gpu_refactor`（水平半拉格朗日平流） | 2.69 | 23% |
| `ndslfv_monoadvh_fgnl_gpu_refactor` | 1.32 | 11% |
| `ndslfv_monoadvv_gpu`（垂直） | 1.25 | 11% |
| `trngra_gpu` | 0.78 | 7% |
| `hdiffu_gpu` | 0.47 | 4% |
| `ndslfv_monoadvv_fgnl_gpu` | 0.43 | 4% |
| `tranrs1_gpu` / `tranrs` / `transr` / `trandv` / `tranuv` / `trngra3` / `transr1`（譜轉換，各 0.1–0.3） | 1.9 | 16% |
| `siimpl_gpu` | 0.24 | 2% |
| **`intgrt_gpu` 內合計** | **8.6** | 73% |
| 不在 `intgrt_gpu` 內的（`diabat_gpu` 物理由主程式呼叫、輸出、MPI） | ~3.2 | 27% |

**結論**：沒有單一主宰項了——每個部件都比 NVIDIA 慢 10–30 倍，符合「大量小 kernel、launch 延遲主導」的型態（prof5：每步 4 萬次 launch、GPU 只忙 6 s）。
接下來是一般性的最佳化工程而不是找 bug：(1) 水平平流（3.9 s，仍是第一）——`cyclic_cell_ppm` 6 個 kernel 合併、把 tracer 維度 collapse 進去；(2) 譜轉換（1.9 s）——逐 mf 的小 kernel 合併；
(3) 全域：減少 `hipDeviceSynchronize`（第 10/17 層的 `geps_acc_wait_all` 改成 stream 級 wait）；(4) runtime 旋鈕：`LIBOMPTARGET_AMDGPU_STREAM_BUSYWAIT`（signal 等待改忙等，省每次 sync 的喚醒延遲；16:20 試跑 `…_2gpu_busywait`）。
記憶體：目前 2 GPU、pool threshold 0、`GEPS_FFT_PLAN_SETS=8`、otile 64/8192 下，`free=` 仍 > 100 GB/卡，1-GPU 路徑應該也放得下（待跑）。
**`…_2gpu_busywait` 結果（16:21）**：`LIBOMPTARGET_AMDGPU_STREAM_BUSYWAIT=2000000` 每步 11.8 s——**沒差**，`sptend` 不變。不採用。

**速度小結（2026-09-18 16:25）**：2×MI300X 每步 **11.8 s**（起點 740 s，全部在 toolchain/`src/rocm`，數值六步逐位元不變）：
FFT plan set budget 8 ／ 垂直平流 tile 32→8192 ／ D2D 複製走 `omp_target_memcpy` ／ RCCL in-place ／ `kernels` 序列區域改平行（def_cfl、hh）／ pack/unpack tile-local ／ 水平平流 tile 8→64。
距 NVIDIA 0.34 s 還有 ~35 倍，性質已從「某處有 bug」變成「4 萬次小 kernel 的 launch 延遲」，需要 kernel 合併等一般性最佳化（§9 8f）。

#### ★★★ 效能層 P1（2026-09-18 16:30，實作中）：acc2omp 把所有內層 `!$acc loop vector/worker` 序列化——這才是 10–30 倍慢的系統性原因

看 `trngra_gpu` 的產生碼才發現：`!$acc parallel loop gang` + 內層 `!$acc loop vector` 被翻成 `target teams distribute parallel do`（外層）+ **「`! acc2omp: nested loop (serial on parent thread)`」（內層序列）**，
也就是每個 kernel 的執行緒數 = 外層迭代數（384 條緯度／384 個波數），MI300X 20k 條 lane 只用 2%；全部 GPU 原始碼有 **302 個**這樣的內層 `loop vector/worker`。這解釋了為什麼每個部件都慢 10–30 倍、GPU 卻只忙 6 s。
（第 9 層當時把它序列化是為了正確性：內層 `parallel do` 需要正確的私有化，而且不能塞在已經是 `parallel do` 的外層裡。）

**新 pass `_vectorise_inner_loops`（`cmake/acc2omp.py`，預設開，`GEPS_ACC2OMP_VECTORISE=0` 關）**：
- 有內層 `loop` 的 `parallel loop` 外層改成 **`target teams distribute`**（每個 team 一個外層迭代），第一層內層 `loop [vector|worker]` 改成 **`parallel do`**（team 內平行）；更深的 `loop` 仍序列。
- 內層 `parallel do` 的 `private(...)` = 原 `loop` 的 private + 巢狀內被賦值的純量（與第 16/18 層同一套宣告解析）+ 更深迴圈的 DO 變數。
- **有 `reduction(` 或 `seq` 的內層迴圈不動**（平行化 reduction 會改變加法順序，六步 `sptend` 就不會逐位元相同）；外層 `private/firstprivate` 照舊。
- `!$acc end parallel loop` 不再輸出 `end target teams distribute parallel do`（外層可能已變成 `teams distribute`，OpenMP 的 end 指令本來就可省）。
靜態迴歸：104 檔 0 失敗，產生 277 個 `parallel do`（1,061 個外層）。**16:40 全部重編（`gpu_build61`），跑 `…_2gpu_vec` 驗證六步 `sptend` 逐位元相同＋每步秒數。**
同時 `src/rocm/cyclic_cell_ppm_gpu.f90` 的 6 個 kernel 也改成 (ot, inner, i) 三維平行（K2+K3 合併；i 方向沒有跨元素相依，`dqq/dpp` 的累加仍在單一執行緒內同序）。

**`…_2gpu_vec`（16:56，vectorise + PPM i-平行）**：NNMI `bal` 三輪與 step-1 `sptend`（0.43807396825059497）**逐位元相同**；step 1 **20.3 s**（原 11.8——變慢了？含第一步的初始化，要看第 2 步）；
但 step 2 開始時 **memory access fault**（兩個 rank，最後 8 個 kernel 都在 `cyclic_cell_ppm_intp_two_loops_gpu`/`massadvy`）。
二分：先只把 PPM 的 i-平行改寫退回（保留 vectorise pass），重編跑 `…_2gpu_vec2`（17:00）。
`…_2gpu_vec2`（17:18）：**同樣 fault** ⇒ 元凶是 vectorise pass（PPM 改寫先放一邊，備份在 scratchpad `cyclic_cell_ppm_gpu.f90.ipar`）。
rocgdb 下跑太慢（30 分鐘還在 NNMI，每次 launch 都被攔截）→ 改用 `OFFLOAD_TRACK_ALLOCATION_TRACES=true`（`…_2gpu_trk`，18:05）：
```
Device pointer 0x9a9eefff000 points into prior host-issued allocation … Last deallocation:
  cyclic_cell_ppm_intp_two_loops_gpu_ ← cyclic_cell_massadvy_mylonlen_gpu_ ← ndslfv_monoadvh2_fgnl_yx ← ndslfv_monoadvh_fgnl ← intgrt_gpu
Kernel 1: cyclic_cell_ppm_intp_two_loops_gpu @ 183
```
⇒ **use-after-free**：出錯的 kernel 正在讀一塊已被 PPM 常式自己的 `exit data map(delete:…)` 釋放的工作陣列。PPM 是 src/rocm 的純 OpenMP（沒被 vectorise 改到），所以最可能的機制是
**target region 非同步執行**：libomptarget 的 AMDGPU 外掛可以延後同步（`OMPX_FORCE_SYNC_REGIONS` 這個環境變數的存在就是證據），vectorise 之後 kernel 的執行時間分布變了，
`exit data` 的 `hipFree`（threshold=0 立刻釋放）跑在前一個 kernel 結束之前。18:10 用同一 binary 加 `OMPX_FORCE_SYNC_REGIONS=1` 重跑（`…_2gpu_sync`）驗證這個假設。
另外 pass 已加 **`lastprivate`**：內層迴圈裡賦值、迴圈後仍被讀的純量保留序列語意（靜態掃描 104 檔只有 `samfshalcnv_kh_gpu` 的 `tem1`）。
`…_2gpu_sync`（`OMPX_FORCE_SYNC_REGIONS=1`，18:04）：**同樣 fault、同一個 kernel（PPM l200）** ⇒ 不是非同步。
step 1 過、step 2 的第一個水平平流就死，而 step 2 走的是 y-first 路徑（`monoadvh2_fgnl_yx`）。
新增 `GEPS_ACC2OMP_VECTORISE_SKIP=<子字串,…>`（跳過檔案）。`…_2gpu_vec3`（跳過全部 `*ndslfv*`，18:28）：**仍 fault**。
⇒ 元凶不在 NDSL 的檔案裡，可能是物理（step 1 末）弄壞記憶體、或根本不是 vectorise pass 而是同一批 build 的另外兩個改動（`end parallel loop` 改成註解、device-I/O 靜音器修正）。
18:30 `GEPS_ACC2OMP_VECTORISE=0` 全部重翻譯（其他改動保留）跑 `…_2gpu_novec` 定案。
`…_2gpu_novec`（18:51）：**六步全過、`sptend` 逐位元相同、11.8 s/步** ⇒ 元凶確定是 vectorise pass 裡某個（非 ndslfv 的）檔案。
二分 #2（18:55，`…_2gpu_vecdyn`）：只 vectorise 動力（tran*/trngra*/hdiffu/ujoinsr/mpe*/intgrt/cufft/rstrandz/ndslfv*），跳過物理
（`diabat, module_mp, samf, moninedmf, gwd, rrtmg, stochastic, adjptqintp, radsw, prerrtmg, mfpbl, ozphys, prexp, mp_scheme, lightning, dcyc2`）。
`…_2gpu_vecdyn`（19:16，114 個 parallel do）：**六步全過、逐位元相同** ⇒ fault 在物理檔；**但每步仍 11.8 s**——TIMER 顯示有好有壞：
`trngra` 0.78 → <0.16 s、`hdiffu` 0.47 → 0.21 s（大贏），但 `monoadvh` 2.69 → 3.56 s、`monoadvv` 1.25 → 1.57 s（**變慢**）。
原因：外層已經 `collapse(3)`（幾千個 team）而內層迴圈很短時，team 內再開 `parallel do` 的 fork/barrier 開銷大於序列執行；外層只有一個迴圈（384 個 team）時才值得。
**啟發式：外層 `collapse(≥2)` 的不提升**。重新翻譯後提升數 277 → 82（集中在譜轉換 `tranrs1/tranrs/transr/trandv/rstrandz/transr1/trngra` 與 `diabat`、幾個物理常式）。
19:20 全檔（含物理）重編跑 `…_2gpu_vec4`：若過關就直接是結果；若仍 fault，元凶在剩下的 diabat/gwdc/gwdps/nor_gwdp/prerrtmg/stochastic 之中，再二分。
`…_2gpu_vec4`（19:41，82 個提升）：**六步全過、逐位元相同，但每步仍 11.6–11.8 s**——`trngra/hdiffu` 的收益被 `monoadvh` 變慢（2.6 → 3.5 s，`mod_ndslfv_monoadv` 裡剩下的 2 個提升）抵消。
**結論**：在 AMD 上 `teams distribute` + team 內 `parallel do`（generic 模式）的每區域開銷很高，只有「外層迭代少、內層迴圈長」（譜轉換那類）才划算；
NDSL 那種 collapse 過的 kernel 要靠改寫成 (outer, inner, i) 的 `collapse(3)`（SPMD）才會快——也就是我對 `cyclic_cell_ppm` 做的那種改法。
**設定**：vectorise 預設跳過 `ndslfv, mod_ndslfv` 與全部物理檔（物理裡有一個檔提升後會 fault，尚未二分到；`GEPS_ACC2OMP_VECTORISE_SKIP` 可覆蓋預設清單），只留譜轉換／hdiffu／mpe2d／tranuv1；
同時把 PPM 的 (ot, inner, i) 三維平行改寫放回去（先前的 fault 已證明與它無關）。19:45 重編跑 `…_2gpu_vec5`。
`…_2gpu_vec5`（20:06）：**逐位元相同；每步仍 11.8 s**。`monoadvh_fgnl` 1.32 → 0.85 s（PPM 改寫有效），但 `monoadvh` 反而 2.6 → 3.5 s、`monoadvv` 1.25 → 1.5 s，總和不動——
每一版總時間都停在 11.8 s 很可疑：像是**被另一個 rank 或某個同步點卡住**（rank 0 快了就在集合通訊處等 rank 1）。TIMER 沒標 rank，兩個 rank 的行混在一起。
20:10：TIMER 改印 rank（`geps_hide_myrank_t()`），`GEPS_ACC2OMP_TIMERS` 改成檔名清單（`intgrt_gpu,ndslfv_monoadvh_gpu` → monoadvh 內部的 intpx/massadvx/massadvy/we2ns/ns2we 每次呼叫也印），跑 `…_2gpu_timer2`。

**`…_2gpu_timer2`（20:40，每 rank，第 3 步，秒）**：

| 呼叫 | rank 0 | rank 1 |
|---|---|---|
| `ndslfv_monoadvh_gpu_refactor`（含下列） | 3.59 | 3.54 |
| ├ `cyclic_cell_massadvy_mylonlen_gpu` | 1.93 | 1.92 |
| ├ `cyclic_cell_intpx_jlist_gpu` | 1.48 | 1.28 |
| ├ `cyclic_cell_massadvx_jlist_gpu` | 0.40 | 0.39 |
| ├ `para_ns2we_gpu` / `para_we2ns_gpu`（NCCL 轉置） | 0.37 / 0.25 | 0.48 / 0.24 |
| `ndslfv_monoadvv_gpu` | 1.58 | 1.45 |
| `ndslfv_monoadvh_fgnl` / `monoadvv_fgnl` | 0.88 / 0.76 | 0.80 / 0.67 |
| 譜轉換合計（tranrs/transr/trandv/tranuv/trngra3/transr1/tranrs1/trngra） | ~1.5 | ~1.9 |
| `siimpl` / `hdiffu` | 0.21 / 0.21 | 0.24 / 0.20 |

兩個 rank 平衡（intgrt 內各 ~9 s + 物理 ~2.8 s = 11.8）⇒ **不是 MPI 等待**；11.8 s 是實打實的工作，只是分散在很多地方。最大的三塊是 NDSL 水平平流的 `massadvy`（1.9）、`intpx`（1.5）與垂直平流（1.6）。
20:45 對這個 build 再做 kernel 級 profile（`…_2gpu_prof6`）看這三塊裡各 kernel 的時間/次數，決定下一個改寫對象（massadvy 內部有 `target update from(qq(…))` 回 host、def_cfl 等）。

**prof6（21:00，PPM 三維平行的 build）**：每步 42k 次 launch，**GPU 只忙 2.45 s**（最大的 kernel：RCCL 0.28、`massadvx l79` 0.21、`hdiffu l52` 0.17、`siimpl 轉置` 0.13、PPM `l102` 0.11 s；PPM 其他 kernel 已降到 0.3–1.5 ms）。
⇒ 11.8 s 裡 **~9.4 s 是 launch/同步的延遲**：平均每次 launch ~220 µs（HSA 派送 ~20–40 µs 之外，還有 libomptarget 每次 launch 對幾十個引數查映射表、以及我們 shim 裡 1.6 萬次 `hipDeviceSynchronize`）。
kernel 本身已經不是問題；接下來是**純粹減少 launch 數與每次 launch 的固定成本**：
(a) ROCr 旋鈕 `HSA_ENABLE_INTERRUPT=0`（signal 等待改輪詢，省每次等待的中斷喚醒；`…_2gpu_noint`）——**結果（21:00）：無效**，每步 11.80/11.82/11.77/11.82 s，六步 `sptend` 與 clean 逐位元相同；ROCr 的 signal 等待不是那 220 µs 的來源；
(b) `geps_acc_wait_all()` 從 device 級 `hipDeviceSynchronize` 改成 `hipStreamSynchronize`（null stream + shim 建過的每條 stream，`g_all_streams`），dgemm 後的 sync 也改成只等 handle 的 stream——**結果（`…_2gpu_ssync`，binary 23:44，00:01）：無效**，每步 11.79/11.83/11.76/11.81 s，六步 `sptend` 逐位元相同。改法保留（無害）。
⇒ 和 prof6 的 per-second 分解一致：每步 12 s 裡，**譜轉換段（~6 s）每秒 1 萬次 launch、GPU 忙 0.25 s/s**；**NDSL 水平＋垂直平流段（~4 s）每秒只有 60–700 次 launch、HIP API 幾乎沒在動、GPU 忙 <0.1 s/s**——後者不是 launch 延遲，是 host 在做別的事（假設：每次呼叫 `cyclic_cell_ppm_intp` 的 `enter/exit data` 配置 ~2.5 GB 工作陣列，threshold 0 ⇒ 每次都走 HSA alloc/free；massadvy/intpx/vertical 的 wrapper 也各自每次呼叫 map 一批）。00:02 用 `rocprofv3 --hsa-core-trace --hsa-amd-trace --memory-allocation-trace` 直接量（`…_2gpu_prof7`，scratchpad `prof7/`），同時把 PPM 的工作陣列改成 module 級持久配置（`src/rocm/cyclic_cell_ppm_gpu.f90`：`module cyclic_cell_ppm_work`，同形狀只配一次）。

**prof7 結果（`…_2gpu_prof7`，00:17；rank 0 主執行緒、第 4 步的 12 s 視窗；`rocprofv3 --hsa-core-trace --hsa-amd-trace --memory-allocation-trace`，結果在 scratchpad `prof7/r_3582039_results.db`（sqlite，分析腳本 `prof7/an2.py`、`an3.py`）；追蹤下每步仍 11.8 s，`sptend` 六步不變）**：

| 主執行緒在 HSA API 裡的時間／步 | 次數 | 秒 | 說明 |
|---|---|---|---|
| `hsa_signal_wait_scacquire`（等 GPU） | 10.2k | **4.88** | 含 GPU 真忙的 ~2.45 s ＋ launch→完成的延遲 |
| **`hsa_amd_memory_pool_free`** | 3,643 | **2.32**（avg 0.64 ms） | 1–16 MB 的 free 也要 0.43 ms、128 MB–1 GB 1.35 ms、≥1 GB 4.9 ms |
| **`hsa_amd_memory_pool_allocate`** | 3,634 | **1.78** | ≥1 GB 76 次 1.11 s（同一 1637 MB 的配置有時 5.7 ms、有時 12 ms；一個 6.5 GB 的要 229 ms）；128 MB–1 GB 611 次 0.55 s；<128 MB 共 0.11 s |
| `hsa_signal_store_screlease`（doorbell） | 49k | 0.48 | 每次 launch ~10 µs |
| 其他（async_copy、pointer_info…） | — | 0.14 | |

⇒ **每步 4.1 s（35%）花在 HSA 配置／釋放**——這就是 NDSL 段「GPU 閒、HIP API 沒動」的 host 時間（那些秒裡 alloc+free 佔 0.75–1.1 s/s）。來源是 `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0`（§6 第 6/7 層時為了 VRAM 不回收而關掉 libomptarget 的池子）讓每個 `enter/exit data`／`!$acc data create` 都直接走 HSA，而原始碼（OpenACC 設計，nvhpc 有池子）每次呼叫都 create/delete 工作陣列。先前 `mm16` 實驗無效的原因也清楚了：<16 MB 的 alloc 只佔 0.08 s，貴的是 free（不分大小）與 ≥128 MB 的 alloc。
另外剩下 12 − 9.6 ≈ 2.4 s 是 HSA 之外的 host 時間（libomptarget 查表、rocFFT/hipBLAS host 端、MPI、Fortran host 迴圈）。
00:20 實驗：`LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=68719476736`（64 GB，全部走池子；binary 不變）跑 `…_2gpu_mmpool`，同時每 10 s 記 `rocm-smi` VRAM（scratchpad `vram_mmpool.log`）看池子是否單調成長。
**`…_2gpu_mmpool` 結果（00:26）：❌ 如第 6/7 層所記，libomptarget 的池子單調成長**——VRAM（`vram_mmpool2.log`）73 → 76 → 82 → 95 → 197 GB（滿），第 2 步 `HSA_STATUS_ERROR_OUT_OF_RESOURCES`（exit 134）；第 1 步 `sptend` 仍 0.43807396825059497。池子本身不能用。

**launch/alloc 微基準（scratchpad `launchbench/lb.f90`，單卡、job 的環境變數；`run.sh`）**：trivial kernel 39 µs/launch（launch＋同步）；`target update to` 8 bytes 54 µs；`enter/exit data` 37 KB 28 µs、**28 MB 242 µs、850 MB 285 µs**（模型裡因碎片化更貴）。`LIBOMPTARGET_STACK_SIZE`、`HSA_XNACK` 對這些數字沒影響。另一個附帶量測：collapse(2) 外層＋每 thread 串行走 `i` 的 kernel（翻譯碼常見型態）讀 8 個 3D 陣列要 785 µs，是 coalesced 頻寬的 ~10 倍——GPU 忙的 2.45 s 裡有多少是這種 pattern，之後再查。

**修法（00:40）：自製 HSA 配置快取 `src/rocm/geps_hsa_pool.cc` → `build_rocm_hip/lib/libgeps_hsa_pool.so`（`src/rocm/build_hsa_pool.sh` 編；刻意不進 CMake，避免 reconfigure）**，由 `job/regression.ksh` `LD_PRELOAD`（`GEPS_HSA_POOL=0` 關）。攔 `hsa_amd_memory_pool_allocate/free`：≥256 KB 的配置依 (pool, flags, 大小 class) 快取，大小進位到 3 個尾數位（最多 12.5% 浪費，讓大小略有不同的輻射 block 也能重用），總量上限 48 GB/rank（`GEPS_HSA_POOL_CAP_GB`），超過就真釋放最久沒用的；配置失敗先清空快取重試。threshold 仍是 0（libomptarget 每次都來叫我們）。微基準：28 MB alloc+free 242 → 28 µs、850 MB 285 → 34 µs。跑 `…_2gpu_hsapool`（binary 不變，只加 preload；`GEPS_HSA_POOL_STATS=1` 在結束時印命中數）。
**`…_2gpu_hsapool` 結果（00:41）：每步 6.96 / 4.98 / 4.98 / 4.91 / 4.93 s（原 11.8），六步 `sptend` 逐位元相同**（`.43807396825059497 .4462741042925732 .4275272147313331 .403981044588656 .3968162229541496 .39255042414676805`），`HAS ENDED`。快取統計（整個 run、每 rank）：alloc 60.4k（hit 23.6k、miss 5.2k、<256 KB 直通 31.6k）、free 39.8k（進快取 25.7k、淘汰 1.8k、flush 0）、結束時快取 45.8 GB。VRAM（`vram_hsapool.log`）：NNMI 期 69–77 GB，預報期 95 → 135 → 155 → 158 GB 後持平（live ~110 GB ＋ 快取 ≤48 GB；206 GB 卡）。NDSL 計時：massadvy 1.45 → 0.21 s、intpx 0.90 → 0.16 s（每次 monoadvh 呼叫）。
`src/rocm/cyclic_cell_ppm_gpu.f90` 的「持久工作陣列」改法因此**不需要**（已還原，快取讓那 9 個 alloc 變成 ~30 µs 的命中）。
⚠️ 未驗證：長預報（tau ≫ 1）下的 VRAM 是否持平（快取有上限，live set 由模式決定；建議跑一次 tau=3 看 `rocm-smi`）；快取回傳的記憶體不清零（裸 HSA 實務上也不保證），六步逐位元相同是目前的證據。
下一步：在 4.9 s 的狀態重新做 kernel/HIP API profile（prof8），看剩下的 ~2.4 s 非 GPU 時間在哪（launch 延遲 vs rocFFT/hipBLAS host 端 vs MPI）。

**prof8（`…_2gpu_prof8`，00:58；kernel＋HSA＋HIP trace，追蹤下每步 7.6 s；scratchpad `prof8/an8.py`）**，rank 0 主執行緒、一步：
- kernel 39.7k 次、GPU 忙 2.59 s：RCCL 56 次 0.43 s（兩卡之間是 **PCIe**，`rocm-smi --showtopo`）、`massadvx_jlist l79` 2×104 ms、`hdiffu l52` 2×83 ms、`siimpl 轉置` 12×11 ms、PPM `l102` 75×1.5 ms、微物理 saticel 各 2×15–38 ms、gwdc 34、diabat 30、深/淺對流 25/24 ms……
- 主執行緒：等 signal 3.64 s；**`hipStreamSynchronize` 26 萬次 0.78 s**——(b) 的 `geps_hip_wait_all` 對 shim 建過的每條 stream 都等（模式用了 ~350 個 async id ⇒ 350 條 stream）；`hipModuleLaunchKernel`（rocFFT）31k 次 0.72 s；HSA alloc/free 剩 0.09 s（快取生效）。
- 0.5 s 分箱的時間軸：**譜轉換（FFT＋dgemm）段佔 ~5.5 s（GPU 忙 30%）**、NDSL 段 1.5 s、物理 0.5 s。FFT 段每 0.5 s 有 2–5k 次 rocFFT kernel launch，是純 host launch 綁定。

修法（01:00，同一個 build）：
1. `geps_hip_wait_all` 只等「上次 wait 之後有派工作的 stream」（`use_stream()` 在每個吃 stream 的 shim 入口標記，wait 後清空）；
2. acc2omp 對 `ndslfv_pack_gpu` 加 fixup：`cyclic_cell_massadvx_jlist_gpu` 三個在裝置上序列跑的區域（`xreg_dup(:, k, :) = xreg` 是 72 個 thread 各抄 148k 元素、`outer_index(i, :) = lons_size`、`kernels` 的 `nstep_less = (nst .le. nstep)`）展成 collapse(2)/(3) 的顯式迴圈。
**`…_2gpu_advx` 結果（01:28）：每步 6.63 / 4.63 / 4.64 / 4.63 / 4.69 s（原 4.9），六步 `sptend` 逐位元相同。**

微基準補充（`launchbench/lb2.f90`）：同一個 8 陣列 stencil，`collapse(2)`＋thread 內串行 `i` 741 µs、`collapse(3)` 244 µs、`teams distribute`＋內層 `parallel do` 238 µs——翻譯碼「外層 collapse＋內層 seq」型態的 kernel 記憶體存取不 coalesced，慢 3 倍；GPU 忙的 2.6 s 裡這類 kernel 之後逐一處理。

下一步（01:30 編譯中）：FFT 段是 host launch 綁定 ⇒ `src/rocm/rfftmlt_loop_gpu.f90` 的逐緯度 exec 迴圈改成 **`GEPS_FFT_THREADS`（預設 4）個 host thread 各用自己的 stream 發**（每個 plan 有自己的 work area，plan 之間無共享；worker thread 要先 `hipSetDevice` 到 rank 的卡，HIP 的 current device 是 per-thread）；同樣的 kernel、同樣的資料 ⇒ 逐位元不變。跑 `…_2gpu_fftthr`。
**結果（01:55–02:05）**：`…_2gpu_fftthr` 因為預設值的 bug 實際跑 1 thread（4.58/4.68/4.63/4.64 s，= 基準）；`…_2gpu_fftthr4`（`GEPS_FFT_THREADS=4`）**4.44 / 4.54 / 4.53 / 4.58 s——只快 ~0.1 s**，六步 `sptend` 逐位元相同。FFT 迴圈本身不是那 5.5 s 的主體（或 HIP runtime 的 launch 路徑有全域鎖，多 thread 發也序列化）。預設改成 4（無害），下一步用 TIMER build 把譜轉換段按常式拆開再看（02:05 編譯：含 hdiffu/rayleifr 的 `reduction(max:)` 內層迴圈平行化——acc2omp 的 vectorise pass 原本跳過所有 reduction，max/min 是精確運算，順序無關，逐位元安全；hdiffu `l52` 是 72 個 thread 各掃 20 萬元素 = 83 ms×2）。

**❌ `…_2gpu_wmax` 結果（02:25）：`sptend` 變了**——`.4380482957086226 .44625022265080605 …`（第 1 步就差 3e-5，不是捨入）。每步 4.28–4.35 s 但數值錯，不算數。
**根因（02:35，standalone 證實，scratchpad `redmax/t.f90`、`s.f90`）：amdflang roc-7.2.2 的「`target teams distribute` 之下巢狀 `parallel do reduction(...)`」結果永遠是 0**——max 與 `+` 都一樣、`-O0` 一樣、有無 collapse 一樣、`teams`＋`distribute` 拆開也一樣；同一個 reduction 放在扁平的 `target teams distribute parallel do reduction` 上是對的。所以 hdiffu 的 `wmax` 全成 0 ⇒ 水平擴散係數變了。這是 **compiler bug**（要進 §9 上游回報清單），acc2omp 的 vectorise pass 對 reduction 一律不升級（`GEPS_ACC2OMP_VEC_REDMAX=1` 才開，預設關）。**教訓：任何要在 team 內做 reduction 的 kernel，在這個編譯器上都得手寫成扁平 kernel 或兩段式（先每 thread 寫部分結果到陣列，再一個 kernel 合併）。** 02:40 重編（gate 關）驗證回到逐位元相同。

TIMER（`…_2gpu_wmax`，intgrt 內各常式，rank 0 最後一步；數值雖錯、時間分佈仍有參考價值）：`read_mtnvar` **0.66**（！每步讀地形檔？）、`ndslfv_monoadvh_gpu_refactor` 0.63、`ndslfv_monoadvv_gpu` 0.42、`monoadvh_fgnl` 0.35、`tranrs` 0.22、`transr` 0.19、`siimpl` 0.19、`monoadvv_fgnl` 0.19、`trandv` 0.18、`tranuv` 0.16、`trngra3` 0.11、`transr1` 0.08、其餘 <0.05；合計 ~3.4 s（intgrt 內）＋物理（diabat，不在這張表）。
（`read_mtnvar` 只在第一步出現一次——一次性初始化，不是每步成本。）

**`…_2gpu_timer2`（02:50，gate 關、timers 在 intgrt＋diabat）：六步 `sptend` 回到逐位元相同；每步 4.49 / 4.57 / 4.50 / 4.59 s。** rank 0 每步（intgrt 呼叫＋diabat 呼叫）：NDSL 平流 1.63（monoadvh 0.67、monoadvv 0.42、monoadvh_fgnl 0.36、monoadvv_fgnl 0.18）、譜轉換 1.06（tranrs 0.23、trandv 0.19、transr 0.19、tranuv 0.15、trngra3 0.11、transr1 0.09、tranrs1 0.07、trngra 0.04）、hdiffu 0.21、siimpl 0.19、物理常式 ~0.45（深對流 0.11、adjptqintp 0.08、淺對流 0.07、gwdc 0.07、其餘 <0.02）；**`rrtmg_gpu` 1.96 s 但 6 步只叫 1 次**（第 1 步 6.5 s 與其他步 4.5 s 的差；長預報要看它的呼叫頻率）。callee 合計 3.6 s，其餘 ~0.9 s 是 intgrt/diabat 內嵌的 kernel 與 host 碼（diabat 整體的 timer 因為 injector 把累加行放進 `#ifdef Readaeroclx` 裡而沒算到——已修：續行走訪跳過 cpp 行並吞掉 `#else…#endif` 尾巴）。

03:00 編譯：hdiffu/rayleifr 的 wind-max 掃描改成兩段精確 kernel（acc2omp fixup：每 (k,jj) 部分 max → 每 k 合併；max 精確，逐位元安全；不用巢狀 reduction）。
**`…_2gpu_wmax2` 結果（03:45）：每步 6.39 / 4.29 / 4.36 / 4.33 / 4.36 s（原 4.5），六步 `sptend` 逐位元相同。** 修好的 diabat 計時：**`diabat_gpu` 每步 1.20 s**（其中 `mp_scheme_gpu` 微物理 0.48、深對流 0.11、adjptqintp 0.08、淺對流 0.07、gwdc 0.07，其餘 ~0.35 是 diabat 內嵌 kernel/host 碼）。每步 4.33 s 的組成（rank 0）：NDSL 平流 1.54、物理 1.20、譜轉換 ~1.05、siimpl 0.20、其他 ~0.3。
（附帶：自製 SIGUSR2 取樣 profiler（scratchpad `hostsamp/`）會讓模式卡在 iteration 1——signal 打斷 HSA 等待；mpiexec 也會把 signal 轉發給 rank。作廢，host 取樣仍用 rocgdb（`sample_bt2.sh`）。）

**host backtrace 直方圖 #3（`…_2gpu_bt4`，4.33 s/步的 build，rank 0，30 個樣本，scratchpad `host_samples4.txt`、`sample_bt3.sh`）**：21/30 在 `__tgt_target_kernel`（等 OpenMP kernel 跑完：微物理 `saticel_s` 4、PPM 3、深/淺對流各 2、輻射 2、其餘各 1）、dgemm 相關 4（`geps_blas_dgemm` 2、**`hipblasSetStream` 2**——Legendre 迴圈每個 m 都換 stream，rocBLAS 換 stream 不便宜）、rocFFT exec 2、libomptarget `deleteData` 1。⇒ host 端自己的計算幾乎不佔時間，**剩下的是 kernel 執行本身（物理的長 column kernel）＋每次 launch 的等待**。
修法（04:00，acc2omp `_defer_graph_loop_syncs`）：六個 `*_gpu_cuda_graph` 的 Legendre dgemm 迴圈（tranrs/transr/trandv/tranuv/trngra3/transr1，每步 ~5k 個 dgemm）原本每個 m：`cublasSetStream(lt_cg_stream(m))`＋`cudaStreamWaitEvent`＋`geps_acc_wait`＋dgemm（shim 內再 sync）＋`cudaStreamSynchronize`；改成迴圈前 `SetStream(handle, stream)` 一次、`geps_blas_defer_sync(1)` 關掉 shim 的逐 dgemm sync，迴圈內只剩 dgemm（同一條 stream 依序排隊），迴圈後 `defer_sync(0)`＋`geps_acc_wait_all()`。同樣的 kernel、同樣的順序 ⇒ 逐位元不變。跑 `…_2gpu_dgemmq`。
**`…_2gpu_dgemmq` 結果（04:20）：每步 6.13 / 4.04 / 4.12 / 4.03 / 4.08 s（原 4.33），六步 `sptend` 逐位元相同。** 譜轉換：tranrs 0.24→0.14、transr 0.19→0.13、trandv 0.19→0.11、tranuv 0.15→0.10、trngra3 0.11→0.08（合計 1.05 → 0.75 s）。

**下一個線索（prof8 的 kernel 派送參數）**：物理的 column kernel（微物理 `saticel_s` 各 15–38 ms、深/淺對流 25 ms、diabat 內嵌 kernel）都是 **grid = 27648 work-item = 108 個 256-thread workgroup**——trip count 27648（每次呼叫的 column 數）÷ 256 = 108 個 team，**MI300X 有 304 個 CU，七成 CU 沒事做**，而這些 kernel 每 thread 用 128 VGPR＋216 AGPR、scratch 72–296 B/lane（register 壓力大，每 SIMD 只能放 1–2 個 wave），所以更需要把 team 攤到所有 CU。試 `OMP_TEAMS_THREAD_LIMIT=64`（不用重編；432 個 64-thread team）。注意：翻譯碼裡有幾個 `reduction(+:…)` 的扁平 kernel（ubartmp/vbartmp/rolltmp/wk4_*/dutmp），team 數變了加總順序會變，`sptend` 可能不再逐位元相同——跑 `…_2gpu_tl64` 看。
**`…_2gpu_tl64` 結果（04:35）：整體變慢 4.06 → 4.35 s**（NDSL 0.61→0.81、rrtmg 1.96→2.37 變慢；物理只快 0.05）；`sptend` 仍逐位元相同（那些 `+` reduction 顯然不影響 sptend，或 runtime 的樹狀加總與 team 數無關）。不採用。

**真正的原因（04:40，看翻譯碼）**：物理 kernel 是 `collapse(2)` 跑 (緯度 jj, 層 k)，**經度迴圈 `do i = 1, myim(jj)` 在每個 thread 裡串行**（原 OpenACC 是 `gang collapse(2)` ＋ `loop vector` over i；vectorise pass 對物理檔一律跳過，而且「外層已 collapse 就不提升」的啟發式也擋掉）。所以 27648 = 384 緯度 × 72 層 個 thread，各自走 ~1000 個經度、跨 thread 的記憶體存取不 coalesced（微基準：慢 3 倍）。全模式有 96 個這種 kernel（module_mp 24、samfdeepcnv 13、moninedmf 11、samfshalcnv 10、gwdc 6、adjptqintp 3、…）。
修法（acc2omp 新 pass `_spmd_collapse_inner`，`GEPS_ACC2OMP_SPMD=0` 關、`GEPS_ACC2OMP_SPMD_SKIP` 跳檔）：對 `parallel loop collapse(2)` 底下**完美巢狀**的內層 `loop vector` `do i = lo, myim(jj)`，把 i 併進 collapse(3)：`do i = lo, <檔案的經度上限 ix/nxp/ite>` ＋ `if (i .le. myim(jj))` 守衛；`seq`/reduction/gang 或非完美巢狀的不動。每個 column 的算式完全一樣 ⇒ 逐位元不變。首批檔案（`_SPMD_EXT`）：samfdeepcnv/samfshalcnv/moninedmf/mfpbl/ozphys_2015/gwdc（`ix`）、adjptqintp（`nxp`）、module_mp_gsfcgce（`ite`）；同一個 build 也把 `mpe2d_gpu` 兩個「一個 thread 抄整段陣列」的 kernel（siimpl 轉置 12×11 ms、reshape_pl）展成 collapse 迴圈。04:45 編譯，跑 `…_2gpu_spmd`。
**`…_2gpu_spmd` 結果（05:05）：每步 5.60 / 3.49 / 3.54 / 3.46 / 3.48 s（原 4.06），六步 `sptend` 逐位元相同。** 物理 `diabat_gpu` 1.20 → 0.81（`mp_scheme` 0.48 → 0.20、深對流 0.11 → 0.064、adjptqintp 0.08 → 0.064、淺對流 0.07 → 0.036）、siimpl 0.18 → 0.035（轉置 kernel）。折進去的 column 迴圈：module_mp 24、samfdeepcnv 11、moninedmf ~11、samfshalcnv ~9、gwdc、adjptqintp 2、ozphys、mfpbl。剩下沒折的（非完美巢狀或不在清單）：rrtmg 3、nor_gwdp 3、samfdeepcnv 2、rozphys 2、dcyc2/lightning/ozphys/samfshalcnv 各 1、gwdps 1——都是每步 <0.02 s 的常式（rrtmg 除外，它 1.98 s/次但不是每步）。
現在每步 3.5 s 的組成（rank 0）：NDSL 平流 1.58（monoadvh 0.63、monoadvv 0.41、fgnl 0.35＋0.20）、物理 0.81、譜轉換 ~0.7、其他 ~0.4。

**prof9（`…_2gpu_prof9`，kernel trace only，05:20；scratchpad `prof9/`）**：一步 39.7k 次 launch，**GPU 只忙 1.64 s**（RCCL 0.30、ndsl-h 0.26、rocFFT 0.19、ndsl-v 0.16、diabat 0.15、微物理 0.11、dgemm 0.08…）；把「GPU 空檔」歸到前一個 kernel：**rocFFT kernel 之間 30k 個空檔 × 30 µs = 0.92 s**、dgemm 之後 2.9k × 73 µs = 0.21 s（hipblasDgemm 的 host 端 ~70 µs/次）、**PPM `l102` 之後 75 次平均 2.4 ms（696 次呼叫裡 288 次有 ~6 ms 的停頓；prof8 的 HSA trace 顯示那段主執行緒完全沒在 HSA/HIP API 裡——host 在做別的事或被擋住，原因未明）**、vertical PPM 各 kernel 之後 0.2–1 ms、def_cfl 之後 0.8 ms。
微基準補充（`launchbench/lb3–lb7`）：**每個 allocatable（有 descriptor）陣列引數讓一次 launch 多 ~10 µs**（2 個 35–45 µs、17 個 190 µs；flang 每次 launch 都把 descriptor 重新 `to` 一次，`LIBOMPTARGET_INFO` 可見「Copying data … Size=88」），explicit-shape 引數幾乎免費（17 個 24 µs）；`-fno-defer-desc-map`、`map(present,alloc:)`、`has_device_addr` 都救不了。模式的 kernel 大多用 explicit-shape dummy，module allocatable（`gglati`、`fa1..4`、`jlist1`…）引用要注意。
環境變數實驗（不重編）：`GPU_MAX_HW_QUEUES=16`＋`GEPS_FFT_THREADS=8` → **3.65 s，變慢**（不採用）。
`GEPS_FFT_THREADS=1` → 3.61 s（4 thread 值 0.12 s，預設維持 4）。`HSA_ENABLE_INTERRUPT=0` → 3.45 s（≈ 3.49，雜訊內）。
**PPM `l102` 之後的 ~6 ms 停頓（05:40，prof8 逐 API 時間軸）**：主執行緒在 doorbell 之後進入一次 `hsa_signal_wait_scacquire`，**kernel 本身 1.2 ms 就跑完（GPU 時戳），但 wait 7.2 ms 才回來**；同一步裡跑 <0.8 ms 的 `l102` 呼叫之後只有 25 µs 空檔，跑 >1 ms 的都停 4–6 ms（288/696 次）。前面 4-byte 的 `has_error` H2D copy 早就完成。`HSA_ENABLE_INTERRUPT=0`（改輪詢）沒改善 ⇒ 不是中斷喚醒延遲，像是 completion signal 本身晚到（或 libomptarget 的 stream 等待在 >1 ms 的 kernel 上換成別的等法）。05:50 用 `HSA_ENABLE_INTERRUPT=0` 再跑一次 kernel＋HSA core trace（`…_2gpu_prof10`）看輪詢下空檔是否還在。
**prof10 結果（06:05）：輪詢下空檔還在**（kernel 2–5 ms 的之後 62% 有 >2 ms 空檔，中位數 3.2 ms；主執行緒的 `hsa_signal_wait` 在 kernel 結束後 4.5 ms 才回）。統計（prof9，3 步）：**>2 ms 的空檔合計 4.0 s ⇒ 每步 ~1.3 s（步長的 40%）**，全部跟在 `private_segment_size = 0` 的 kernel 後面（PPM `l102` 97/225 次、vertical PPM `l2148/l2241` 各 45/333、def_cfl、adjptqintp `l242` 19 ms 的 kernel 之後停 25 ms…），沒有其他 kernel 同時在跑（全在 queue 1），RCCL 也沒在跑。微基準（`launchbench/lb8`、`lb9`）在單獨程序裡看不到：kernel 0.26–7.6 ms 時 wall = kernel 時間，不均勻的 workgroup 也沒事。
已排除（各跑一次，不重編）：`HSA_ENABLE_INTERRUPT=0`（3.45）、`LIBOMPTARGET_AMDGPU_STREAM_BUSYWAIT`（先前）、`HSA_SCRATCH_SINGLE_LIMIT=2GB`（3.48）、`GPU_MAX_HW_QUEUES=16`（3.65 更慢）、`OMP_TEAMS_THREAD_LIMIT=64`（更慢）；`HSA_NO_SCRATCH_RECLAIM=1`（單獨或加 8 GB single limit、或把快取上限降到 12 GB）都在第 2 步 `HSA_STATUS_ERROR_OUT_OF_RESOURCES`（VRAM 只用到 85–96 GB ⇒ 是 scratch 上限，不是 VRAM）——證明**有 dispatch 要的 scratch 超過 queue 的上限**（prof9：6 個 dispatch 的 private segment 是 65536 = `LIBOMPTARGET_STACK_SIZE`，65536×64 lane×全佔用 wave 數 = 數十 GB，ROCr 走 use-once 路徑），但那 6 個不是主要停頓來源。`LIBOMPTARGET_STACK_SIZE=1024` 會卡死在第 1 步（sflx 需要它）。
KFD queue 數（`/sys/class/kfd/kfd/proc/<host pid>/queues`）：每 rank 在自己的卡上 8 個 compute queue＋2 個 SDMA（另一張卡上 2 個，RCCL peer）——沒有超額訂閱（24 個硬體 slot）。
裝置映像裡有 `__llvm_rpc_client`／`__llvm_omp_emissary_rpc`（某處 device 端 I/O 沒被靜音器拿掉），libomptarget 因此啟用 RPC server thread——**假設**：kernel 等待路徑跟 RPC 有關（`HSA_ENABLE_INTERRUPT`／busywait 都無效與此相容），待驗證：找出還在 device code 裡的 print/stop 並移除，讓映像不含 RPC client。06:55 先試 `HSA_ENABLE_SCRATCH_ASYNC_RECLAIM=0`（`…_2gpu_noasyncrec`）。
`…_2gpu_noasyncrec`：3.50 s，無效。
**找到了（07:10）：裝置映像裡的 RPC client 就是停頓來源。** 用 `llvm-objcopy --dump-section=.llvm.offloading` 抽出裝置 ELF，解析 `s_getpc/s_add` 呼叫目標，找到還在 device code 裡的 I/O：`sflx_gpu l309` 的 `write (*, *)`＋`stop 333`（靜音器的 regex `write\s*\(\s*\*\b` 在 `*,` 前永遠不匹配）、`radiation_aerosols` 的 `stop`、`gocart mass2icn_gpu`（`!$acc routine` 裝置常式）的 `stop '…'`。這些把 `_FortranAio*`/`_FortranAStop*` → `__llvm_omp_emissary_rpc`/`__llvm_rpc_client` 連進映像，libomptarget 看到 `__llvm_rpc_client` 就開 RPC server thread，kernel 等待改走 RPC 路徑——長 kernel 的完成要晚 ~一個 kernel 長度才被看到。
修法（acc2omp `_silence_device_io`）：regex 改成 `write\s*\(\s*(\*|6)\s*[,)]`；`stop`/`error stop` 也註解掉，`if (…) stop '…'` 改成 `if (…) continue`；bare `!$omp declare target`（`!$acc routine` 來的）到 `end subroutine` 之間視為 device code；續行的 `&` 後面可以跟註解（sflx 的 `!%f`）。重編後映像沒有 rpc/Fortran I/O 符號。
**`…_2gpu_norpc` 結果（07:50）：每步 4.26 / 2.82 / 2.81 / 2.82 / 2.78 s（原 3.49），六步 `sptend` 逐位元相同。** NDSL 1.58 → 1.15（monoadvh 0.49、monoadvv 0.28、fgnl 0.25＋0.13）、物理 0.81 → 0.67、rrtmg 1.97 → 1.43（每 6 步一次）。
教訓（值得記）：**device code 裡任何 `print/write/stop` 不只會 Recursive I/O，只要連進映像就讓所有 kernel 的完成通知變慢**；檢查方法：抽裝置 ELF 看 `.dynsym` 有沒有 `__llvm_rpc_client`。

**prof11（`…_2gpu_prof11`，08:05，2.8 s build 的 kernel trace；追蹤下每步 3.57 s）**：一步 39.7k 次 launch，GPU 忙 1.60 s、空檔 1.97 s。分組（busy／之後的空檔）：rocFFT 31.6k 次 0.19／**0.87**、ndsl-h 859 次 0.26／0.24、ndsl-v 1566 次 0.16／0.30、RCCL 56 次 0.26／0.03、dgemm 5k 次 0.08／0.19、diabat 0.15／0.06、微物理 0.11、對流 0.07、gwd 0.07。**注意 rocprof 本身每個 kernel 加 ~20 µs**（追蹤 3.57 vs 不追蹤 2.80，差 0.77 s ≈ 31.6k×24 µs），所以 FFT 那 0.87 s 空檔大部分是追蹤造成的，不追蹤時譜轉換整段（含 dgemm）約 0.6 s（TIMER）。
不追蹤時每步 2.8 s 的組成：NDSL 1.15（monoadvh 0.49＋fgnl 0.25、monoadvv 0.28＋0.13；其中 RCCL 轉置 ~0.25、PPM kernel ~0.4、其餘是 ~2.4k 次 launch 的空檔）、物理 0.67、譜轉換 ~0.6、其他 ~0.4。
08:10 編譯：NDSL tile 加大（水平 PPM/intpx/massadvy `otile` 64 → 128、垂直 `otile` 8192 → 32768）把 launch 數砍半／砍四分之三；工作陣列大小成比例（PPM 的 dqmono/qi 各 3.3 GB、垂直 qmi/qpi 各 360 MB，有快取無妨）。跑 `…_2gpu_tile2`。
**`…_2gpu_tile2` 結果（14:25）：每步 4.30 / 2.49 / 2.52 / 2.49 / 2.54 s（原 2.80），六步 `sptend` 逐位元相同。** NDSL 1.15 → 0.86（monoadvh 0.49→0.40、monoadvv 0.28→0.19、fgnl 0.25→0.19、0.13→0.075）。再加大一級（水平 256、垂直 131072，看 VRAM）跑 `…_2gpu_tile3`。
`…_2gpu_tile3`（14:55）：2.55 / 2.55 / 2.54 / 3.11 s——**沒有更快**（水平 0.40→0.38，垂直 0.19→0.24 變慢，最後一步 3.1 s），`sptend` 不變。退回 128／32768。
**SPMD pass 第二批（15:15）**：owner 沒有 collapse（只有 `do jj`）也處理（→ collapse(2)），且允許外層 DO 與內層 `loop` 之間夾純量賦值前綴（`j = jlist1(jj); nxj = nxdef_2d(j)`，複製進折疊後的迴圈體；前綴變數在迴圈體裡不能再被賦值）；新增檔案 diabat（`nxp`，42 個）、gwdc 7、gwdps 2、rrtmg 6（`nx`）、nor_gwdp 2、dcyc2、lightning。**`…_2gpu_spmd2` 結果（15:30）：每步 4.11 / 2.41 / 2.43 / 2.44 / 2.46 s（原 2.51），六步 `sptend` 逐位元相同。** diabat 0.67 → 0.58。
(c) rocFFT 逐緯度 exec（每步 ~6k 次、25k 個 kernel）：Bluestein 長度每次 3–4 個 kernel——無法 batch（長度全不同），只能減少 sync；
(d) 長期：把常一起出現的小 kernel 合併（例如 `vertical_cell_ppm_intp` 的 6 個、PPM 的 7 個）以及讓 target region 非同步（`nowait` + depend）把延遲疊起來。

一般性教訓（值得記）：**在 ROCm 上，任何拿 libomptarget 映射的裝置指標去呼叫 `hipMemcpy*`（直接或經由 RCCL/rocFFT/hipBLAS 內部）都可能掉到 CPU 路徑**；D2D 一律用 `omp_target_memcpy`，集合通訊用 in-place 形式，並用 rocprofv3 的 `hipMemcpyAsync` avg 時間（正常 <1 ms）當健康指標。
- 這與 NVIDIA 參考 log「GPU 偏低」的方向一致，**可能就是 GPU 版 SAS 對流與 CPU 版的實作差異而非 ROCm 問題**——但在 ROCm 上先做到「GPU 深/淺對流各自和 CPU 差多少」再下結論；
  之後若要判定移植正確性，最乾淨的是拿 NVIDIA 跑同一組 `c_*`/`g_*`。

###### ⚠️ 附帶發現 1：自編 OpenMPI 的 `MPI_DOUBLE_PRECISION` 是 **16 bytes**（2026-09-15，已量測）

`GEPS_LIB/src/openmpi/compile.sh` 用 `FCFLAGS="-O1 -fPIC -fdefault-real-8"` 設定 OpenMPI，
**沒有 `-fdefault-double-8`**，所以 configure 時 amdflang 的 `DOUBLE PRECISION` 被提升成 16 bytes：
`opal_config.h` 裡 `OMPI_SIZEOF_FORTRAN_DOUBLE_PRECISION 16`。模式本身用
`-fdefault-real-8 -fdefault-double-8`，`double precision` 是 8 bytes。實測（scratchpad `mpisize.f90`）：

```
MPI_TYPE_SIZE: DOUBLE_PRECISION=16  REAL8=8  REAL=8  REAL4=4  INTEGER=4
```

後果：任何用 `MPI_DOUBLE_PRECISION` 的呼叫都搬 **2 倍** 的位元組——傳送端讀過界、接收端寫過界。
我的探針就是這樣死的（32 rank：gather 回來一半是垃圾、`deallocate` 時 `double free or corruption`；
4 rank 獨立重現：allreduce 得到 26 而非 52）。

**對模式的影響（目前設定下是潛伏的，不是今天 `bal` 的原因）**：`src/` 裡用到它的地方
（`mpe2d.f90` 的 `mpe2d_transpose_ndsl_*_multi`/`mpe2d_row_broadcast`/`mpe2d_unify_uzmean`/`mpe2d_unify_nx_lev*`、
`mpe_unify.f90` 的 `mpe_unify4_r`、`mod_mpe.f90` 的 `mpe_global_maxloc/minloc_r8`（`MPI_2DOUBLE_PRECISION`）、
`mpe_send_data.f90`/`mpe_recv_data.f90`）在 NPEX=1 時不是 size-1 通訊子（不搬資料）、就是沒被呼叫
（`uzmean`、`rstrantq`、`intgrt_3tl`）、或只是診斷用的 maxloc；GPU 端對應的路徑走 NCCL。
**但 NPEX>1、IO server、或 `sigful`（輸出時的 `mpe2d_unify_nx_lev_red`）都會踩到。**

修法在 toolchain（符合 §7 規矩）：`compile.sh` 的 FCFLAGS 加上 `-fdefault-double-8` 重編 OpenMPI。
在那之前，探針一律用 `MPI_REAL8`（已驗證 8 bytes）。

###### ⚠️ 附帶發現 2：GPU 執行檔裡單 rank 的 `MPI_GATHER`/`MPI_ALLREDUCE`（`MPI_REAL8`）會弄壞 heap（2026-09-15，未解）

同一段探針改成 `MPI_REAL8` 後，32-rank CPU 執行檔正常，**1-rank GPU 執行檔**在 `getrdy` 的
CPU `tranrs` 探針之後立刻 `Fatal glibc error: malloc.c: sysmalloc assertion`（top chunk 被寫壞）。
獨立小程式（同編譯旗標、1 與 4 rank、`jtmax=385`）**重現不出來**。
`np == 1` 時改成不呼叫 MPI 直接複製（`src/tranrs.f90`、`geps_dbg_ssq_grid`）就不再崩潰，
所以肇事者確定是那兩個 MPI 呼叫在 GPU 執行檔裡的行為，但機制未查
（候選：GPU 執行檔用 `-fopenmp --offload-arch` 連結時 `mpi` 模組介面的差異）。**只影響探針，不影響模式。**

###### 📎 使用者提供的兩份參考 log（2026-09-15）

- 另一台機器的 **CPU**：`bal = 1.7619e-05 / 6.3946e-06 / 1.6451e-06`，`surf pres tend rms ≈ 0.42 mb/hrs`，
  與本機 CPU build **完全一致** ⇒ 本機 CPU 基準可信。
- **NVIDIA GPU**：`bal = 1.7619218832540225E-005 / 6.3945974572339345E-006 / 1.6451316157579604E-006`，
  `pt 0.92182 −7.66227 → 0.91441 −7.67721`，**與 CPU 逐位元相同**，每步 0.34 s。
  ⇒ GPU 版程式碼本身是對的、`bal` 沒有任何「GPU 容許誤差」的空間；ROCm 版必須做到同一個數字。

###### ★ 下一步：把 CPU 基準做細一點

目前的 CPU 對照只有 `max|x|`、`max|evec|`、`bal` 三個量，**沒有 `eval` 的基準**，
也沒有**逐 (m, L) 的分布**——只有全域最大值。
而 `bal = Σ_mf Σ_{k,i≤nn} wrk²`，是**整個分布的和**，
最大值相當不代表和也相當（GPU 的 `max|x|` 只比 CPU 大 7 倍，`bal` 卻大 1e8 倍，
這本身就說明差異在分布而不在極值）。

建議的下一個量測（CPU 與 GPU 各跑一次，各約 5 / 17 分鐘）：
1. **`eval` 的量級**——CPU 版目前沒印，GPU 是 256.273，缺對照；
2. **逐 (m, L) 的 `bal_tmp(mf, L)`**——把幾個代表性的 `mf` 的值印出來，
   看是「所有項都偏大」還是「少數項爆掉」。後者會直接指出是哪個 `ind`，
   再回頭看那個 `ind` 的 `nn`/`eval`/`wrk`。
3. 同樣印 `Σ wrk²` 與 `Σ x²`（不是 max），兩邊比對。

**這比繼續提假設有效率**——今天在這一層提了六個假設、六個都被量測推翻，
共同原因是**只有全域極值、沒有分布**，極值相當就誤以為輸入相當。

###### （已作廢）下一個該查的：`nnlist` 在裝置上的值

**非決定性是目前最強的訊號**——它代表某處讀到未初始化的記憶體。
`wrk`（= `nnmi_buf`）是 `enter data create` 建立的，**裝置端不清零**。
算 `bal` 的 kernel 只讀 `i <= nn`，其中 `nn = nnlist(ind)`。
**如果裝置上的 `nnlist` 是錯的或過期的，`nn` 就會偏大，kernel 就會讀到 dgemm 從未寫過的區域**——
正好同時解釋「量級離譜」與「每次跑批都不同」。

`nnlist` 在 `initial_gpu` 是 `enter data create(nnlist)` ＋ `!$acc update device(mx, nnlist)`
（`initial_gpu.f90:221`）。值得查的點：
1. 那個 `update device` 有沒有真的生效（可用 `geps_dbg_mapped` 同樣的手法，或
   在 kernel 內把 `nn` 寫進 `bal_tmp` 之類的旁路輸出回來比對）；
2. CPU 版對應的 `ns`/`na` 值是多少（`initial.f90:190-191` 算出來的，可直接印）。

**驗證成本**：跑到第一個 `bal=` 約 17 分鐘（`job/TCo383L72_IC_sample_rocm.log.20260911_wrk2` 是最近一次）。

###### 方法論：今天在這一層誤判了五次 → 已訂成規矩，見 **§7.1 探針規矩**

問題 B（abort 當下的工作集當成洩漏）、`max|mx|`（未初始化填充區當成壞特徵向量）、
`nn` 範圍（前 4 筆樣本外推）、`wrk` 緩衝區（過濾條件抓錯呼叫）、
**正交歸一性檢查（`dbg_n >= 64` 上限只掃到小矩陣）**。
**四次的共同點都是「指標／樣本選得不對，就得到看似有力、實則無效的證據」。**
其中兩次還是靠「算術上界推不出觀測值」這種論證去推斷的——
**那類論證今天錯了兩次，不要再單獨依賴它**；
先把量測本身的有效性確認好（取樣是否有代表性、量到的是不是目標區域），再談推論。

###### （已作廢）GPU 自己的數字互相矛盾

- `bal = Σ wrk²` ≈ 2624 ⇒ 需要 `|wrk| ~ 1`
- `x_out = evec·(wrk/e_t)` = 4e-6，而 `|e_t| ≤ 256` ⇒ 只需要 `|wrk| ~ 1e-3`

**同一個 `wrk`，兩個讀取者推得的量級差約 1000 倍。**
`bal` 是 OpenMP kernel 讀的，`x_out` 是第二次 hipBLAS `dgemm` 讀的。
這與「**dgemm 透過 `use_device_addr` 寫的緩衝區，和 OpenMP kernel 透過一般映射讀的緩衝區，不是同一塊**」
高度吻合——也就是第 1、2 層那一類的 bug。

**下一步（直接驗證這個假設）**：在 `nnmi_gpu` 裡比對兩個位址——
在 `!$omp target data use_device_addr(evec, x, wrk)` 區塊內把 `c_loc(wrk(1,1,1))` 印出來，
在區塊外用 `omp_get_mapped_ptr(c_loc(wrk(1,1,1)), omp_get_default_device())` 取得
OpenMP 認定的裝置位址，**兩者不一致就確認了**。
`src/rocm/hip_compat.cc` 的 `[dsyevj]` 診斷（`geps_rocm_dsyevj_dev` 開頭那段）
就是現成可以照抄的作法，而且已經證實有效。
若確認，修法可比照**第 2 層**：給讀 `wrk` 的 kernel 加 `has_device_addr(wrk)`，
或改用 `omp_get_mapped_ptr` 明確取址後再交給 `dgemm`。

###### 實作這個修法時踩過的四個坑（都已解決，記錄供日後修改 acc2omp 參考）

1. **外層指令是多行接續**：第一版在「第一行」找 `private(...)`，而該行的
   `map(tofrom:loc_max, &` 括號未閉合，子句被插進中間 → `Unmatched ')'`。
   **在 joined 層做就沒有這個問題**（`_join_continuations` 已把續行接好）。
2. **去重沒處理巢狀括號**：`private(i, ..., dqmono(:3*lons))` 用
   `\(([^()]*)\)` 掃不到外層 → 重複加入 → `'i' appears in more than one data-sharing clause`。
   要用**括號配對掃描**。
3. **owner 的作用範圍**：`!$acc parallel loop` 只涵蓋緊接其後的迴圈巢狀。
   不追蹤 DO 深度就會把後面（甚至別的 subroutine）的 private 加上去 →
   `No explicit type declared for '<name>'`。
4. **不可跨條件編譯提升**（最難找的一個）：
   `module_mp_gsfcgce_3ice_nuwrf_gpu.f90` 裡 `mr2mc` 的宣告與內層 loop **都在
   `#ifdef Readaeroclx` 之內**（`CMakeLists.txt:46` 預設 **OFF**），
   而 owner 指令在外面。把 private 提升等於**把名稱從條件編譯區塊裡拉出來**，
   巨集關閉時就引用到不存在的東西。**必須比對 cpp 條件深度，不同就不提升。**

###### 靜態迴歸（這次建立，建議保留）

對 `src/nvidia/*.f90`、`src/*_gpu.f90`、`src/cwa/*.f90` 全部翻譯一次，檢查：

1. 每個 `!$omp` **邏輯指令**的括號平衡（續行接起來、跨 `#ifdef`）→ 抓坑 1
2. 同一指令內資料共享子句的變數不得重複 → 抓坑 2
3. private 名稱必須在同一 subroutine 內出現 → 抓坑 3

三項現在都是 104 檔 0 問題。**坑 4 這三項抓不到**（名稱在檔案裡確實存在），
只有編譯器會抓——所以改完 acc2omp 一定要真的編一次。

###### ⚠️ 改壞 acc2omp 之後一定要做的事

**cmake 的相依掃描會先讀既有的產生檔**，遇到語法/語意錯誤就中止，**不會重新翻譯**。
所以改壞一次之後即使修好了，build 仍會對著舊檔重複報一模一樣的錯
（症狀：build 只跑 10~20 秒就失敗）。

```bash
rm -f build_rocm_hip/src/hip_src/*.f90     # 強制全部重新翻譯
```

這次因此浪費了兩輪，以為修法沒生效。

##### 方法論上的教訓

這一層之所以卡了這麼久（2026-09-11 一整天、誤判三次），
是因為**一直在 GPU 側猜「哪個值不合理」，而沒有先建立 CPU 基準**。
CPU smoke 只要 5 分鐘、可以反覆跑；
一旦有了 `x`/`evec`/`bal` 的對照，「只有 `bal` 錯」這件事立刻就浮出來，
搜尋範圍從整條 NNMI 鏈縮到單一個 kernel。
**以後遇到數值不符，第一件事就是讓 CPU 版吐出同樣的量。**

#### 其他備忘
- 工作樹裡 `eigen.f90`/`eigenlib.f90`/`eigrs.f90`/`initial.f90`/`mpe_init.f90`/`mpe_unify.f90` 的 `! #region agent log` 除錯探針，**經確認是給 CPU 路徑的 amdflang -O3 spill bug 用的（見上面已繞過的問題 1），跟這整段 GPU omptarget 崩潰無關**（GPU build 走的是 `initial_gpu.f90` 等 `src/nvidia/` 檔案，不是那些有探針的 CPU 檔案）。
- `src/` 裡的這些 `! #region agent log` 探針（寫 NDJSON 到 `.cursor/debug-*.log`，其中一個已 49 MB）在 GPU port 收斂前可以獨立清掉，不影響 GPU 除錯。
- `compile` 的輸出目錄與 job 腳本不一致（見 §3 警告）。
- 目前 GPU 修法的**淨變動只有兩個檔案**：`src/nvidia/mpe2d_gpu.f90`（77 行）+ `src/nvidia/zx_gpu.f90`（35 行），全部 `#ifdef USE_HIP`，NVIDIA 路徑完全不受影響、字面上一行都没变。所有診斷用的 print/sync 探針、graph-capture-disable 實驗都已經清乾淨、revert 掉。
- 所有改動**都還沒 commit**。

---

## 7. 工作規矩（`.cursor/rules/prefer-no-src-changes.mdc`，仍然有效）

> **盡可能不要改 `GEPS/src`；編譯器 / MPI ABI 問題在 toolchain 解。**

**不要做**
- 為了 MPI 傳錯 pointer，在 GEPS 裡自幹 C trampoline / `mpe_allgather_*` wrapper
- 為了 amdflang VLA/descriptor，大規模把 automatic array 改 `allocate`、改 1D pack
- ❌ 跑 `./compile rocm`（會 `rm -rf build_rocm`）
- ❌ 動 `build_rocm_hip/`

**要做**
- MPI 炸掉：先查有沒有原生 ROCm MPI；沒有才在 **OpenMPI 編法**（`GEPS_LIB/src/openmpi/compile.sh`）解
- CPU 增量編譯：`cmake --build build_rocm_cpu --target tcogfs.x`（先 `module load geps/1.0`，`COMP_MP=mpifort C_COMP_MP=mpicc`）
- 真的必須改 src：只改與 ABI 無關的 bug，並寫清楚為什麼 toolchain 解不了

### 7.1 探針規矩（2026-09-11 加上，使用者明確要求「要記得錯誤 不要再犯了」）

> **在從探針下任何結論之前，先證明「這個量測代表我要推論的對象」。**

這一天所有被推翻的判斷（五次）都是同一個病因：**量測本身沒先驗證過代表性**，
但它看起來像證據、讀起來也像證據。

| # | 錯誤 | 真正的原因 |
|---|------|-----------|
| 1 | 把 `abort()` 之後的 working set 當成 leak | `abort()` 會丟掉 buffer 裡的 trace；`SIGTERM`/`SIGINT` 才保留 |
| 2 | `max\|mx\| > 1e12`，宣告特徵向量是垃圾 | 整段 `maxval` 量到**未初始化的 padding**：column 長 `no*no`，只用 `nn*nn` |
| 3 | 從 `nn=1~3` 推出「算術矛盾」 | 只看了前 4 筆（`seen < 4`）；實際 `nn = 1..576` |
| 4 | `[dgemm]` 印出的形狀對不上 | `seen < 6` 抓到最早的 dgemm（m=144, n=1..6），不是 nnmi 的。改用形狀過濾 `n == 2 && m == k` |
| 5 | 正交歸一性檢查「乾淨」 | `dbg_n >= 64` 上限 + 迴圈順序 `(L, m, k)` ⇒ 只掃到**小矩陣**（ns=48…146），`bal` 的質量所在的 ns 到 530 從沒被檢查過 |

**規矩（下結論前必須先答得出這兩題）**

1. **樣本代表性**：這個探針實際看到的是哪個子集合？為什麼那個子集合能代表結論的對象？
   任何 `seen < N` / `dbg_n >= N` 之類的上限，**先把迴圈順序追出來**再用。
2. **量測範圍**：量到的量是否覆蓋真正重要的區域？**兩邊都要查，不是只查 GPU。**
   2026-09-12 我對 `cc` 兩側取整陣列平方和，CPU 那側回報 **NaN**——
   CPU 的 `cc` 有未初始化的堆疊 padding。GPU 側看起來很正常（1.259e7），
   差點就拿一個有意義的數字去對一個無意義的數字。
   （只含已初始化的範圍、正確的呼叫點、活著的行程而不是屍體。）

兩題有一題答不出來 ⇒ **這個探針還不算證據**，先擴大再重跑，不要寫進 handoff。

3. **`.cursor/debug-a51a4c.log` 是跨跑批累積的——數記錄一定要從本次的 `RUN_START` marker 起算。**
   正確寫法是 `START=$(wc -l < $DBG)` 再 `tail -n +$START $DBG | grep -c …`。
   2026-09-14 我在等待條件裡用了全域 `grep -c '"vram:427"'`，而該 locid 在
   先前多次跑批就存在 ⇒ 條件在第一次檢查（t=10s）就成立，
   **跑批在什麼都還沒算出來時就被殺掉**，而輸出是「探針都到了 (10s)」——
   一個看起來很像成功的假象。先前幾次僥倖沒中，只因為每次都用了新的 locid。

4. **`target update from(v)` 之前，先確認 `v` 在裝置資料環境裡——查的是整個環境，不是當前子程式。**
   對沒被 map 的陣列下 `update from` **不報錯、沒有任何跡象**，只會安靜地給你主機舊值。
   但反過來也要小心：2026-09-11 我查了 `initial_gpu` 的 map 子句沒看到 `vorten`，
   就宣告量測無效——其實它在 `mod_spec.f90:75` 的**模組層級** `enter data` 裡，
   量測完全有效，我差點把一條真線索當成假象丟掉。
   **`enter data` 的映射是裝置資料環境的全域狀態，與宣告它的子程式無關。**

5. **「主機 == 裝置」驗證的是傳輸，不是內容。**
   未初始化的垃圾 `update device` 過去之後，兩邊當然一致。
   要驗證的是**值本身合不合理**（範圍、有沒有被寫過），不是兩邊一不一致。

6. **只比較同一個量的兩個版本，絕不跨物理量比大小。**
   2026-09-11 我看到「輸入 `Σtt²=1.8e13` → 輸出 `Σphiten²=2.3e8`」
   就斷定「這個常式沒有放大什麼」。溫度的平方和與位勢傾向的平方和
   **單位不同、本來就不可比**，這個推論是空的。
   同一個量的 GPU/CPU 比值才是證據——實際上是 **1.1e8 倍**，結論完全相反。

7. **套修法之前，把所有會被它影響的輸出都先量一次當基準。**
   2026-09-12 我對 `tranrs_gpu_cuda_graph` 套修法前只量了 `temten`，
   沒量同一常式另一個呼叫的輸出 `hldten`。修法後發現 `hldten = 0`（應為 1.33e8），
   卻**無法斷定是既有的還是修法造成的**。一個量測的缺口讓因果永久不可判。

8. **數字要有出處**：寫進「前提 / 事實」表的每個數字，都要能指回
   **哪一次跑批、哪一個 locid**。2026-09-11 我把一個沒有出處的 `max|x| = 7.1e-4`
   當成已量測的事實，據此推出一個「超出上界 4500 倍」的矛盾，
   花了一整輪去找不存在的洞——實測是 **2.483**。
   **憑空的前提比錯誤的量測更危險，因為它連被複驗的機會都沒有。**

另外：假說就標成假說，並寫下「什麼量測能定案」；量測跑完才升格成結論。
這一天有六個聽起來很有道理的理論被直接量測推翻
（hipSOLVER 正確性、`use_device_addr` 指標傳遞、「特徵向量是垃圾」、private/race、
「dgemm 和 kernel 用不同 buffer」——位址其實一模一樣 `0xc7e2ee00000`、`nnlist` device copy——host==device），
其中兩個已經寫進 handoff 才被迫撤回。撤回比慢一點得到答案傷得更重。

---

---

## 8. 快速指令

```bash
# 環境
export GEPS_LIB_ROOT=/mlsteam/workspace/data/geps/GEPS_LIB
export GEPS_DMSDB=/mlsteam/workspace/data/geps/dmsdb
export GEPS_COMPILER=rocm
. /usr/share/modules/init/bash && module purge
source $GEPS_LIB_ROOT/env.sh && load_geps_compiler
module use $GEPS_LIB_ROOT/install/modulefiles && module load geps/1.0
```

```bash
# CPU 增量編譯（安全）
cmake --build /mlsteam/workspace/data/geps/GEPS/build_rocm_cpu --target tcogfs.x -j 32
```

```bash
# GPU 增量編譯（不要 reconfigure）
cmake --build /mlsteam/workspace/data/geps/GEPS/build_rocm_hip --target tcogfs.x -j 32
```

```bash
# 萬一 reconfigure 了（改過 src/CMakeLists.txt 就會）：一定要先 export 這三個，否則 link 找不到 mpi_*（2026-09-16 踩到）
export COMP_MP=mpifort C_COMP_MP=mpicc CXX_COMP_MP=mpicxx
cmake -Bbuild_rocm_hip -S. -DCMAKE_BUILD_TYPE=Release -DUSE_RSM=OFF -DCWBSUM=OFF -DTIMCOMCPL=OFF -DUSE_aeroclx=OFF \
      -DGEPS_COMPILER=rocm -DUSE_CUDA=OFF -DUSE_HIP=ON -DUSE_ACC=OFF -DUSE_PAR=OFF -DUSE_OMP=ON -DUSE_NDMS=OFF -DAMD_GPU_ARCHS=gfx942
```

```bash
# CPU 樹萬一 reconfigure 了（同上，改過 src/CMakeLists.txt 就會；2026-09-16 20:10 踩到：mpi.mod 找不到）
export COMP_MP=mpifort C_COMP_MP=mpicc CXX_COMP_MP=mpicxx
cmake -Bbuild_rocm_cpu -S. -DCMAKE_BUILD_TYPE=Release -DUSE_RSM=OFF -DCWBSUM=OFF -DTIMCOMCPL=OFF -DUSE_aeroclx=OFF \
      -DGEPS_COMPILER=rocm -DUSE_CUDA=OFF -DUSE_HIP=OFF -DUSE_ACC=OFF -DUSE_PAR=OFF -DUSE_OMP=ON -DUSE_NDMS=OFF
cmake --build build_rocm_cpu --target tcogfs.x -j 48
```

```bash
# CPU smoke（1h 預報）
cd /mlsteam/workspace/data/geps/GEPS/job && ./TCo383L72_IC_sample_rocm_cpu 2>&1 | tee cpu_$(date +%Y%m%d_%H%M).log
```

```bash
# GPU smoke（1 MPI on MI300X）— 記得存 log
cd /mlsteam/workspace/data/geps/GEPS/job && ./TCo383L72_IC_sample_rocm 2>&1 | tee gpu_$(date +%Y%m%d_%H%M).log
```

```bash
# 判斷跑成功：log 尾巴要看到這兩行
#   tau=    1.000,     Timing= ...
#   PROGRAM CWBGFS HAS ENDED.
```

```bash
# 版控（2026-09-18 起）。remote = https://github.com/DengShunChen/GEPS.git，branch = feature/multi-compiler
# 坑 1：這台機器沒有 GitHub 憑證（無 gh、無 credential helper、無 SSH key），連 fetch 都要帳密。
#        不要把 token 存進 ~/.git-credentials（明文）；用一次性的 credential helper，token 只留在環境變數：
export GH_PAT=...   # classic PAT (ghp_…) 帶 repo scope；fine-grained 的話 Contents 要 Read and write，否則 403 "Write access not granted"
git -c credential.helper='!f() { echo username=x-access-token; echo password=$GH_PAT; }; f' push origin feature/multi-compiler
unset GH_PAT        # 用完到 GitHub revoke
# 坑 2：repo 的 pre-push / post-commit hook 需要 git-lfs（fix/*.f77、*.dat 走 LFS）。2026-09-18 已裝 git-lfs 3.4.1；
#        沒裝的話 push 會被 hook 擋，只有在 commit 沒碰 LFS 檔時才可以 git push --no-verify 繞過。
# 進版控前的檢查：git diff --cached --name-only | grep -i '\.log' 必須是空的（run log 不進 repo）。
```

---

## 9. 建議的下一步（依序）

1. ~~修 08-28 的 CPU 回歸~~ **已完成（2026-09-08）**：根因是容器重建後系統套件遺失（見 §5、§6）。
2. ~~GPU 端找崩潰點、修前 3 層~~ **已完成（2026-09-08）**：三個各自獨立的 root cause，`iteration=1` 從此能完整跑完。
3. ~~修第 4 層：`iteration=2` 開頭的 `HSA_STATUS_ERROR_OUT_OF_RESOURCES`~~ **✅ 已結案**：不是單一 subroutine 的問題，是記憶體耗盡的下游症狀。
4. ~~修 (B)：積分階段的映射洩漏~~ ❌ **不存在，2026-09-10 已推翻**：`diabat_gpu` / `cwa_rrtmg_gpu` 的 create≫remove 是 abort 當下執行中呼叫的工作集，不是洩漏。**不要再去找「漏掉的 exit data」。** 詳見 §6「問題 B」。

5. ~~修第 6 層（RCCL 集合通訊後的競態）~~ ✅ **已修好並驗證（2026-09-10）**：`cmake/acc2omp.py` 對每個帶 `async_id` 的 `nccl_*` 呼叫後插入 `call geps_acc_wait(async_id)`（4 檔 12 處，未動 `src/nvidia/`）。迴歸測試從「332s 必定 fault」變成「61 分鐘 0 fault、三輪 NNMI 全過」。**⚠️ 但只修掉 fault，`NaN` 與記憶體耗盡都還在**——第 6 層是必要非充分條件。第一版只在 `exit data` 前插等待、只改一個檔案，測試沒過，教訓見 §6「第 6 層」。
   `cmake/acc2omp.py` 把 `!$acc exit data delete(swork) async(async_id)` 的 `async` 子句丟掉，
   而 `nccl_alltoall`（`src/nvidia/helper.f90:66`）不是 OpenMP 構造、仍然把 `ncclSend`/`ncclRecv`
   丟進 stream 就返回，於是 325 MiB 的 `swork` 在 RCCL 還在讀時就被釋放。
   **修法**：在 acc2omp 對「原本帶 async 的 exit data」之前插入 `call geps_acc_wait(id)`
   （`src/rocm/geps_acc_wait.f90` 已存在）；保守起見可先只對 `mpe_transpose_rs_gpu` 點名修補
   （該檔 3 處，是唯一有「nccl 呼叫後緊接 exit data」模式的檔案）。
   **驗證很便宜**：`LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=16777216` 重跑，
   目前會在 **345 秒**明確 fault；修對了就該跨過去，不必等 55 分鐘。詳見 §6「第 6 層」。

6. ~~跑「預設門檻 ＋ 第 6 層修法」的乾淨對照~~ ✅ **已完成（2026-09-11），結果是負面的**：探針與基準線**逐位元相同**（`tendget_gpu` −95.65 GB、iteration 1 結束 5.92 GB、iteration 2/3 delta 全為 0、積分 `NaN`、結局 `OUT_OF_RESOURCES`）。**第 6 層的修法對記憶體問題零影響**；先前 16 MB 那次的「改善」完全來自門檻變更。log：`job/TCo383L72_IC_sample_rocm.log.20260911_control`。
   **在修好 #5 之前先不要把它當成獨立問題回報 AMD**——池子不歸還既是耗盡的原因，
   也是第 6 層 use-after-free 的遮罩，兩者高度相關；修好 #5 後要重新量一次，
   才知道還剩多少是真的 runtime 問題。已掃過的門檻值：`0` / 64KB / 1MB / 16MB / 256MB 都會觸發第 6 層的 fault，
   `4GB` 與預設相同。

7. ~~查 `NaN` 的來源~~ ✅ **已查明（2026-09-11）= 第 7 層**：GPU 的 NNMI 是完全的 no-op（`bal` 精確為 0、`pt` 逐位元不變），因為 graph capture 區塊內的 OpenMP kernel 不在被 capture 的 stream 上——錄製期搶先執行讀到空 `wrk`，又沒被錄進 graph。**與第 6 層同源。** 詳見 §6「第 7 層」。
8f. **速度（2026-09-18）**：探針已拆（無影響）；每步 740 → 11.8 s（見 §6 末的逐項）。下一步依序：(i) 1-GPU（NPEY=1）路徑跑通 tau=1；(ii) 水平平流 `cyclic_cell_ppm` 6 個 kernel 合併／tracer 維度 collapse；(iii) 譜轉換逐 mf 小 kernel 合併；(iv) `geps_acc_wait_all`（device 級）改 stream 級；(v) 上游回報清單。量測工具都在：`GEPS_ACC2OMP_TIMERS=1`（build 時）印每步各常式秒數、`GEPS_PROF_WRAP` 接 rocprofv3、scratchpad `sample_bt2.sh` 取 host backtrace。
8e. **✅ 2026-09-17 18:55：GPU 六步與 CPU 對到 1e-11**（深對流 `val1/val2` 修好）。優先序不變：拆探針＋量速度 → 1-GPU → 上游回報（新增：`samfdeepcnv_kh_gpu` 34 vs 29 個引數；`val1/val2`、`snoexp` 等 private-未初始化清單）。
8d. **✅ 2026-09-17 05:00：GPU 已能完整跑完 tau=1**（第 17 層 raw RCCL 等待、第 18 層 `LIBOMPTARGET_STACK_SIZE` + `adjptqintp` guard）。**現在的優先序**：
   (i) **拆探針做乾淨對照**：acc2omp 裡的 DBGMAP/DBGFREE/DBGTRM/DBGHD/DBGMP/DBGSFLX/intgrt_ckpts/diabat_ckpts/pbl_ckpts 全部關掉（建議加一個 `GEPS_ACC2OMP_PROBES=0` 開關，而不是刪碼），
   `src/` 裡的 `! #region agent log` 探針（`intgrt.f90`、`tendget.f90`、`tranrs.f90`、`hdiffu.f90`、`diabat.f90`）也移除，重跑 2-GPU 到 tau=1，確認 6 個 `sptend` 與 l18g 相同、並量真正的每步秒數（NVIDIA 0.34 s）。
   (ii) **量化 GPU 物理 vs CPU 的 4–6%**：用 `p2_*`（兩側已加）先看 `qt`；要追就把 CPU big-loop 的各物理常式輸出改成累加式 Σ² 與 GPU 的 `g_*` 對齊。先確認 NVIDIA 上同樣的檢查點差多少（使用者的 NVIDIA 機器）——若 NVIDIA 也是 4–6%，就不是 ROCm 的問題。
   (iii) 1 GPU（NPEY=1）也跑一次 tau=1，確認單卡路徑（VRAM 現在 pool 關掉後應該夠）。
   (iv-a) **✅ 2026-09-21：toolchain 缺陷的獨立 PoC 已整理在 `qa/rocm_toolchain_poc/`**（7 個目錄各附 `run.sh`，`run_all.sh` 全跑約 6 分鐘，`results/` 留有本機量測；README 是給 AMD 的英文報告）：巢狀 reduction 回 0、descriptor 每次 launch 重抄（assumed-shape 17 個 710 µs）、pool 只重用同大小（16 次 31 GB）＋關 pool 後 1 GB alloc/free 21 ms、`hipMemcpy` 對 libomptarget 指標 0.05 GB/s、裝置 I/O 帶進 `__llvm_rpc_client`（延遲效應獨立程式**重現不出**，README 有註明）、runtime 大小 private 陣列 89% 被蓋、`nowait`＋隱式 firstprivate 純量 NYI 編譯錯誤（新發現）。
   (iv) 上游回報：`hdiffu_gpu` host `wmax` 競態、`intgrt_gpu.f90:1085` `poly/dpoly` 引數、`adjptqintp_gpu.f90:160` 越界讀、`sflx_gpu` `private(snoexp)` 未初始化、OpenMPI `compile.sh` 缺 `-fdefault-double-8`。
8c. ~~**★★★ 積分第 1 步（2026-09-15 晚間起的最高優先）**~~ ✅ 已於 2026-09-16 結案（第 15–17 層；根因是第 17 層的裸 `NCCLCHECK` 集合通訊沒有等待）。原文：`surf pres tend rms(GPU)` = 132.6 vs 0.42。方法照抄今天的：CPU `src/intgrt.f90` 與 GPU `src/nvidia/intgrt_gpu.f90` 的呼叫順序對齊，用 `geps_dbg_ssq_grid`/`geps_dbg_ssq_full`（`src/tendget.f90` 末尾）與 `DBGTRM` 逐波數探針從兩端往中間夾。先做兩件便宜的：(i) 兩次跑同一個 binary 看第 1 步的數字是否相同（不同 ⇒ 還有競態，去 `hip_src/nvidia_intgrt_gpu.f90` 找沒被 wait 的 HIP 非同步操作）；(ii) 關掉探針重跑確認 NNMI 修法自己站得住。**VRAM 耗盡（第 13 項）現在直接擋在第 2 步前面，可能得先處理才有第 2 步可比。**
8a. ~~**★★★ 定位 `ndslfv_monoadvh2_gpu_refactor` 內部**~~ ✅ **已於 2026-09-15 結案**：水平平流本身精確（13 位），元凶是它前後的 memset 競態（第 10 層）與 `trngra3` 的 graph capture（第 7 層延伸）。原文：外層 bracket 已證明 `ddtemp`/`vdzonl`/`vdmerd` 在它的輸出端就錯了（0.48 / 0.54 / 1.38 倍），而輸入 `ttp`/`pt`/`pdot` 完全正確、下游 `tranrs` 整條鏈逐波數正確。下一步是把 bracket 推進去：CPU `src/ndslfv_monoadvh.f90:45` 與 GPU `src/nvidia/ndslfv_monoadvh_gpu.f90:24` 的呼叫順序完全相同（pack → `massadvx` → `intpx` → `we2ns` → `massadvy` → `ns2we` → `intpx` → `massadvx` → unpack），在每一段之後量 `qqlon`/`rrlon`/`qqlat` 的使用區域 Σ²，兩邊 1:1 比。**優先懷疑 `src/rocm/cyclic_cell_intpx_gpu.f90`、`cyclic_cell_massadvy_gpu.f90`、`cyclic_cell_ppm_gpu.f90`**（手寫 tiled 改寫、從未獨立驗證；0.48 ≈ 一半，`reducefactor=0.506`，「reduced grid ↔ full grid 的 intpx 只填了一半」是個可以一次量測定案的假設）。
8b. **修 OpenMPI 的 `MPI_DOUBLE_PRECISION`（toolchain）**：`GEPS_LIB/src/openmpi/compile.sh` 的 FCFLAGS 加 `-fdefault-double-8` 重編；重編前 `src/` 裡不要新增任何 `MPI_DOUBLE_PRECISION` 用法。
8. ~~**★★ 把 CPU 基準做細（最高優先）**~~ ✅ **已用另一種方式完成（2026-09-15）**：逐波數 + 全域 Σ²，見 §6。原文：目前 CPU 對照只有 `max|x|`、`max|evec|`、`bal` 三個**全域極值**，沒有 `eval` 基準、也沒有逐 (m,L) 分布。而 `bal` 是整個分布的和——GPU 的 `max|x|` 只比 CPU 大 7 倍、`bal` 卻大 1e8 倍，**差異必然在分布而非極值**。建議量：(a) `eval` 的量級（CPU 版目前沒印）；(b) 逐 (m,L) 的 `bal_tmp(mf,L)`，看是「所有項偏大」還是「少數項爆掉」——後者會直接指出是哪個 `ind`；(c) `Σwrk²` 與 `Σx²`（不是 max）。CPU 跑批只要 5 分鐘、可反覆做。**今天在這一層提了六個假設、六個都被量測推翻，共同原因就是只有極值沒有分布。**
9. ~~查 `nnlist` 在裝置上的值~~ ❌ **已推翻（2026-09-11）**：主機與裝置的 `nnlist(1)=2`、`nnlist(2)=1` 完全相同，索引範圍靜態上也一致。詳見 §6「第 9 層」。
10. ~~驗證 `wrk` 的緩衝區身分~~ ❌ **已推翻（2026-09-11）**：直接比對位址，`[dgemm]` 寫入的 `C` 指標與 `omp_get_mapped_ptr(wrk)` **完全相同**（`0xc7e2ee00000`），後續呼叫是 base＋偏移。dgemm 與 OpenMP kernel 用的是同一塊緩衝區。詳見 §6「第 9 層」。
11. ~~修第 9 層（內層 private 被丟掉）~~ ✅ **修法已 landed 並通過編譯**（`_privatise_nested_loops`，911 個指令受惠，靜態迴歸 104 檔 0 問題），**但 `bal` 完全沒變 → 競態不是原因**。修法語意正確、建議保留。實作時踩的四個坑與靜態迴歸腳本都記在 §6「第 9 層」。
12. **第 7 層（graph capture eager 化）**：診斷已證實、三項修法已實作（eager 化、`vars` 加 `save`、`pack_event` 直接同步），**崩潰類問題全部消除（`memory access fault` = 0）**。**2026-09-15 延伸到 trngra3/transr/transr1/trandv/tranuv/rstrandz（`lt_cg` 區塊），是 NNMI 修好的兩個必要條件之一。** 詳見 §6「第 7 層」與「第 10 層」。
   都是在 0 可用 VRAM 下產生的垃圾（`dd` 全零、積分 `NaN`）。修好後再跟 CPU 版
   `job/TCo383L72_IC_sample_rocm_cpu.log.20260908` 逐點比對，並回頭確認 `bal= 0.` 是不是真實結果。

13. **查記憶體不歸還（獨立問題，可能要找 AMD）**：映射表全程持平（45.0 GB / 595 筆、約 200 萬次事件全部正確配對），HIP 可見 free VRAM 卻從 157 GB 掉到 0。已排除 GEPS 的 `enter/exit data` 寫法與數量、第 6 層的競態。
14. **統一 build 目錄**：`compile` 與 job 腳本對齊到同一個 HIP 樹，讓 `./compile rocm` 成為唯一入口。

15. **清掉除錯探針**：`! #region agent log`（CPU 路徑 amdflang -O3 spill bug 用，與 GPU 無關），以及 2026-09-10 新增的兩處（`cmake/acc2omp.py` 的 initial_gpu 注入器、`src/rocm/vertical_cell_advect_gpu.f90` 的 host 端 CFL 重算）。⚠️ 後者**會改變執行時行為**（含 `hipDeviceSynchronize()`），做數值驗證前務必移除。**注意：`acc2omp.py` 裡的第 6 層修法（nccl 後插 wait）不是探針，要保留。**
   以及 2026-09-10 新增的兩處（`cmake/acc2omp.py` 的 initial_gpu 注入器、
   `src/rocm/vertical_cell_advect_gpu.f90` 的 host 端 CFL 重算）。
   ⚠️ 後者**會改變執行時行為**（含 `hipDeviceSynchronize()`），做數值驗證前務必移除。

16. ~~**收斂 commit**~~ ✅ **已完成（2026-09-18）**：全部 55 個檔案以單一 commit `4424cf83` 進 `feature/multi-compiler` 並 push（探針仍在，拆探針＝#15/#8d(i) 要另開 commit）。原文：GPU **修法**現在有 3 項：`src/nvidia/mpe2d_gpu.f90` + `src/nvidia/zx_gpu.f90`（112 行，第 1-3 層）＋ `cmake/acc2omp.py` 的第 6 層修法。另有除錯探針若干（見 #11）。加上其餘既有的 25+ 個檔案改動，全部未進版控。
    （`src/nvidia/mpe2d_gpu.f90` + `src/nvidia/zx_gpu.f90`，全部 `#ifdef USE_HIP`），
    加上其餘既有的 25+ 個檔案改動與上述探針，全部未進版控。

17. 之後才談效能：`-O0` 的四個檔要想辦法回到 `-O3`、`outgrb2` 的 SIGSEGV、多 GPU（RCCL）擴展（GPU1 目前全程閒置）、`vars` 永久不釋放的做法要不要統一清理。**注意：第 7 層若改成 eager 執行，等於放棄 graph replay 的效能優勢，收斂之後要回頭評估。**
    多 GPU（RCCL）擴展（GPU1 目前全程閒置）、`vars` 永久不釋放的做法要不要統一清理。

---

## 10. 給下一次的自己：讀這份文件的順序

1. §0 摘要 → **§6「★★★ 第 7 層」**（`NaN` 的根因，最高優先）→ **§6「★★★ 第 6 層」**（同源、已修好）→ §9 下一步
   - §6 的第 1-3 層是已修好的歷史；第 4 層、第 5 層兩節保留的是推理過程，**結論都在「第 5 層根因」那一節**，時間不夠就只讀它
2. §4（懂 acc2omp + src/rocm 的設計，這是整個移植的骨架）
3. §7 規矩（避免重蹈覆轍）
4. 動手前先跑 `git status` 和 `git log --oneline -3` 確認工作樹還是不是本文件描述的樣子（2026-09-18 基準：乾淨，HEAD = `4424cf83`）
