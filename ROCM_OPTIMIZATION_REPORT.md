# GEPS / TCoGFS ROCm 移植：效能優化報告

日期：2026-09-20　範圍：2026-09-18 ～ 2026-09-19 的加速工作（數值移植完成之後）
平台：2 × AMD Instinct MI300X（SPX、PCIe 互連）、ROCm 7.2.2、amdflang 22.0.0git（roc-7.2.2）、自編 OpenMPI 4.1.6
案例：TCo383L72 IC sample，tau=1（6 個積分步），2 MPI rank（NPEY=2，一卡一 rank）

---

## 1. 結果摘要

| 指標 | 起點（2026-09-18 00:00） | 現在（2026-09-19 15:30） |
|---|---|---|
| 每積分步 wall time（2 GPU） | **740 s** | **2.43 s**（第 1 步含初始化 4.1 s） |
| 數值 | 六步 `surf pres tend rms` 與 CPU 對到 1e-11～1e-14 | **每一次改動後六步 `sptend` 逐位元不變**（`.43807396825059497 .4462741042925732 .4275272147313331 .403981044588656 .3968162229541496 .39255042414676805`） |
| 對照 CPU（32 核，同一機器） | 10.5 s/步 | GPU 快 **4.3 倍** |
| 對照 NVIDIA 參考 log | 0.34 s/步 | 差距 305× → **7×** |
| VRAM（每卡 206 GB） | — | 峰值 ~158 GB（含 48 GB 配置快取） |

十三次改動、每次一個變數、每次跑完整 tau=1 驗證；所有 log 在 `job/TCo383L72_IC_sample_rocm.log.2026091[89]_2gpu_*`，逐項紀錄在 `ROCM_PORT_HANDOFF.md` §0/§6。

原始碼原則維持不變：改動全部在 `cmake/acc2omp.py`（OpenACC→OpenMP 翻譯器）、`src/rocm/`（ROCm 專用的 shim 與替換 kernel）與 `job/regression.ksh`；`src/`、`src/nvidia/` 未動。

---

## 2. 方法

**驗證準則**：每次改動後跑完整 tau=1，六個 `surf pres tend rms(GPU)` 必須與上一版逐位元相同（bit-identical），否則退回。此準則在本輪抓到一個編譯器 bug（§5.1）。

**一次一個變數**：每個 build/env 實驗只改一件事；同一 build 若含兩個改動會在報告中註明。

**先量再改**（§7.1 探針規則）：每一步都先用工具指出瓶頸，再改：

| 工具 | 用途 |
|---|---|
| `rocprofv3 --kernel-trace`（`GEPS_PROF_WRAP` 鉤子，需 `OMP_TOOL=disabled`） | kernel 數、GPU 忙碌時間、kernel 之間的空檔（歸到前一個 kernel） |
| `rocprofv3 --hsa-core-trace --hsa-amd-trace --memory-allocation-trace` | 主執行緒在 HSA API 裡的時間（找到 alloc/free、signal wait） |
| `rocprofv3 --hip-runtime-trace` | HIP API 次數（找到 26 萬次 `hipStreamSynchronize`） |
| `GEPS_ACC2OMP_TIMERS=<檔案清單>`（acc2omp 在 build 時插 `system_clock` 計時） | intgrt / diabat 內每個被呼叫常式的每步秒數 |
| rocgdb 取樣（`sample_bt3.sh`） | host 主執行緒在哪裡等 |
| 微基準（scratchpad `launchbench/`、`fftbench/`、`redmax/`） | 把假設從模式裡抽出來單獨驗證 |
| `llvm-objcopy --dump-section=.llvm.offloading` ＋ `llvm-readelf` | 檢查裝置映像的符號（找到 RPC client） |

---

## 3. 優化時間軸（每步秒數，2 × MI300X）

| # | 每步 | 根因（量到的） | 修法 | 位置 | log |
|---|---|---|---|---|---|
| 0 | 740 | 數值移植剛完成、探針剛拆 | — | — | `…0918_2gpu_clean` |
| 1 | 500 | rocFFT plan-set LRU 預算 4 < 每步 7 種 (jump,m,isign)，每步重建 ~8k 個 plan（`hipModuleLoadData` 492 s/run） | plan-set 預算 8 | `src/rocm/rfftmlt_loop_gpu.f90`（`GEPS_FFT_PLAN_SETS`） | `…_fft8` |
| 2 | 247 | 垂直平流 tile `otile=32`（OOM 時代）：每次呼叫 9.4k 個 tile × 8 launch，pack kernel 每 tile 掃全部 column | `otile` 32→8192；column map `col_i/col_j` 一次建好，pack tile-local | `src/rocm/ndslfv_monoadvv_gpu.f90`、`vertical_cell_advect_gpu.f90` | `…_ot8192` |
| 3 | 58 | **HIP 不認得 libomptarget 配的指標**：`hipMemcpyAsync` D2D 走 CPU large-BAR 路徑，26 MB/s（`qm` 一個 343 MB 陣列 87 s） | D2D 改 `omp_target_memcpy`；RCCL AllGather/Broadcast 改 in-place（RCCL 內部的本地複製同病） | `src/rocm/hip_compat.cc` | `…_d2d` |
| 4 | 28 | `!$acc kernels` 翻成序列 `!$omp target`（def_cfl_step_gpu_type2 一個 thread 掃 8192×72，107 ms × 74 次/步；vertical_cell_ppm_intp 的 hh 填值同病）；pack/unpack 逐 tile 掃全域 | acc2omp fixup 改成 `teams distribute parallel do reduction(max:)`／collapse(2)；tile-local pack | `cmake/acc2omp.py` | `…_defcfl` |
| 5 | 11.8 | 水平平流 PPM tile `otile=8`：576 個 team/launch，48 tile × 6 kernel/呼叫 | `otile` 8→64；PPM 六個 kernel 改成 (ot, inner, i) collapse(3) | `src/rocm/cyclic_cell_ppm_gpu.f90`、`cyclic_cell_intpx_gpu.f90`、`cyclic_cell_massadvy_gpu.f90` | `…_ot64` |
| 6 | 4.9 | **每步 4.1 s 在 HSA alloc/free**：`LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0`（必要，libomptarget 的池子只重用相同大小、輻射 block 會漲到 OOM——本輪重試 64 GB threshold 仍在第 2 步 OOM）讓每個 OpenACC `data create` 都是裸 HSA 配置；free 不分大小 0.4–5 ms、≥1 GB alloc 5–230 ms；每步 3.6k+3.6k 次 | **自製 `LD_PRELOAD` 配置快取** `libgeps_hsa_pool.so`：攔 `hsa_amd_memory_pool_allocate/free`，≥256 KB 依 (pool, flags, size class) 快取，大小進位到 3 個尾數位（≤12.5% 浪費），上限 48 GB/rank 淘汰最久未用，配置失敗先清快取重試 | `src/rocm/geps_hsa_pool.cc`、`build_hsa_pool.sh`、`job/regression.ksh`（`GEPS_HSA_POOL=0` 關） | `…_hsapool` |
| 7 | 4.63 | `cyclic_cell_massadvx_jlist` 三個序列區域（`xreg_dup(:,k,:) = xreg` 是 72 個 thread 各抄 148k 元素，104 ms×2）；`geps_acc_wait_all` 對 shim 建過的 350 條 stream 都 `hipStreamSynchronize`（26 萬次/步） | acc2omp fixup 展成 collapse 迴圈；wait_all 只等「上次 wait 後有派工作」的 stream | `cmake/acc2omp.py`、`src/rocm/hip_compat.cc` | `…_advx` |
| 8 | 4.33 | hdiffu/rayleifr 的 wind-max 掃描：`parallel loop gang` over k ＋ `loop vector reduction(max)` 被翻成 72 個 thread 各掃 20 萬點（83 ms×2）；FFT 逐緯度 exec 只有一條 host thread 在發 | 兩段精確 kernel（每 (k,jj) 部分 max → 每 k 合併；max 精確，逐位元安全）；FFT 用 `GEPS_FFT_THREADS`（4）條 host thread、各自 stream（每 thread 先 `hipSetDevice`） | `cmake/acc2omp.py`、`src/rocm/rfftmlt_loop_gpu.f90` | `…_wmax2` |
| 9 | 4.06 | 六個 `*_gpu_cuda_graph` 的 Legendre dgemm 迴圈：每個 m `hipblasSetStream`（rocgdb 取樣有 2/30 卡在裡面）＋`cudaStreamWaitEvent`＋`geps_acc_wait`＋dgemm（shim 內再 sync）＋`cudaStreamSynchronize`；每步 5k 個 dgemm | 迴圈前 `SetStream` 一次、`geps_blas_defer_sync(1)` 關掉 shim 逐 dgemm sync，迴圈內只剩 dgemm 在同一條 stream 依序排隊，迴圈後 `defer_sync(0)`＋`geps_acc_wait_all()` | `cmake/acc2omp.py`（`_defer_graph_loop_syncs`）、`hip_compat.cc`、`geps_acc_wait.f90` | `…_dgemmq` |
| 10 | 3.49 | **物理 kernel 的經度迴圈串行**：原 OpenACC `gang collapse(2)`(jj,k)＋`loop vector`(i)，翻譯後 27648 個 thread 各串行走 ~1000 個經度，跨 thread 存取不 coalesced（微基準慢 3 倍）、只用 108 個 team（304 CU） | acc2omp 新 pass `_spmd_collapse_inner`：完美巢狀的 `do i = lo, myim(jj)` 折進 collapse(3)，`do i = lo, <ix/nxp/ite>` ＋ `if (i .le. myim(jj))` 守衛；每個 column 算式不變。首批 96 個 kernel（module_mp 24、samfdeepcnv 11、moninedmf 11、samfshalcnv 10…）；同 build 展開 mpe2d 兩個單 thread 抄陣列的 kernel（siimpl 轉置 12×11 ms） | `cmake/acc2omp.py`（`_SPMD_EXT`、`GEPS_ACC2OMP_SPMD`） | `…_spmd` |
| 11 | 2.80 | **裝置映像含 `__llvm_rpc_client`**：`sflx_gpu` 的 `write (*, *)`＋`stop 333`、`radiation_aerosols` 的 `stop`、`gocart mass2icn_gpu`（`!$acc routine`）的 `stop '…'` 沒被靜音器拿掉（regex `\*\b` 在 `*,` 前不匹配；`stop` 從未處理），把 Fortran I/O runtime → emissary RPC 連進映像，libomptarget 因此開 RPC server thread，**每個 >1 ms 的 kernel 完成通知晚 ~一個 kernel 長度**（每步 ~1.3 s 空檔） | `_silence_device_io` 修 regex、處理 `stop`/`error stop`/`if (…) stop`、bare `declare target` 常式視為 device code、續行 `&` 後允許註解 | `cmake/acc2omp.py` | `…_norpc` |
| 12 | 2.51 | NDSL 是 launch 空檔綁定（prof11：ndsl-h 859 次、ndsl-v 1566 次 launch） | tile 加大：水平 64→128、垂直 8192→32768（再大一級：無效，退回） | `src/rocm/cyclic_cell_*`、`vertical_cell_advect_gpu.f90`、`ndslfv_monoadvv_gpu.f90` | `…_tile2` |
| 13 | **2.43** | diabat/gwdc/rrtmg 還有單層 `do jj` owner＋前綴純量（`j = jlist1(jj); nxj = nxdef_2d(j)`）的 column 迴圈 | SPMD pass 擴充：單層 owner → collapse(2)、純量前綴複製進迴圈體；再折 60 個 kernel | `cmake/acc2omp.py` | `…_spmd2` |

---

## 4. 根因分類

### 4.1 ROCm runtime 走錯路徑（三次，合計 740 → 4.9 s 的大部分）
1. **rocFFT plan 重建**（#1）：plan 每個 ~20 MB VRAM，384 個長度 × 7 組，LRU 預算太小。
2. **`hipMemcpy*` 拿到 libomptarget 指標走 CPU 路徑**（#3）：HIP runtime 的指標分類不認得 HSA 直接配的記憶體；結果正確、只是 26 MB/s。RCCL out-of-place 集合通訊內部同病。**規則：OpenMP 映射的指標一律用 `omp_target_memcpy` 或 in-place 集合通訊。**
3. **HSA alloc/free**（#6）：OpenACC 程式假設 runtime 有池子（nvhpc 有），libomptarget 的池子又不能用（exact-size 重用、單調成長）。自製快取是這一輪最大的單項（11.8 → 4.9）。
4. **RPC client**（#11）：device code 裡任何 `print/write/stop` 不只 Recursive I/O，連進映像就讓 kernel 完成通知變慢。

### 4.2 翻譯後的 kernel 結構（#2、#4、#5、#7、#8、#10、#13）
OpenACC 的 `kernels` 區域、`gang`＋`loop vector` 巢狀、陣列語法賦值、`loop vector reduction` 在直譯成 OpenMP 後都變成**單 thread 或少 thread 的序列 kernel**。acc2omp 現在有四類自動處理：
- `_spmd_collapse_inner`：column 迴圈折進 collapse（守衛式），適用「外層 collapse(1–2)、內層 `loop vector`、完美巢狀」；
- `_vectorise_inner_loops`：譜轉換那類「外層迭代少、內層長」的 → `teams distribute`＋團內 `parallel do`（NDSL、物理不用：team 內 fork 開銷大於收益）；
- per-file fixup：`kernels` 區域、陣列語法、wind-max 兩段化；
- tile 大小：src/rocm 的 wrapper 用 `otile` 控制每次 launch 的工作量（水平 128、垂直 32768）。

### 4.3 同步太多（#7、#9）
每個 dgemm 3 次 sync、`wait_all` 掃 350 條 stream。原則：**同一條 stream 上排隊的操作只在消費前等一次**；shim 記錄 dirty stream。

---

## 5. 過程中發現的 toolchain 缺陷（要回報上游／要記住）

1. **amdflang roc-7.2.2：`parallel do reduction(...)` 巢狀在 `target teams distribute` 之下結果永遠是 0**（max 與 `+` 都是、`-O0` 一樣、有無 collapse 一樣；扁平的 `target teams distribute parallel do reduction` 正確）。standalone 重現：scratchpad `redmax/t.f90`、`s.f90`。在模式裡表現為 hdiffu 的 `wmax` 全 0 → `sptend` 第一步差 3e-5。acc2omp 對 reduction 一律不升級（`GEPS_ACC2OMP_VEC_REDMAX=1` 才開）。
2. **flang 每次 launch 都重新 `to` allocatable 的 descriptor**：每個 allocatable 引數讓一次 launch 多 ~10 µs（2 個 35–45 µs、17 個 190 µs；explicit-shape 17 個只要 24 µs）。`-fno-defer-desc-map`、`map(present,alloc:)`、`has_device_addr` 都無效。模式的 kernel 大多是 explicit-shape dummy，影響有限，但 module allocatable（`gglati`、`fa1..4`、`jlist1`…）要注意。
3. **libomptarget 記憶體池只重用相同大小**，varying-size 的配置單調成長到 OOM（本輪再驗：threshold 64 GB 在第 2 步 `OUT_OF_RESOURCES`）。
4. **HIP runtime 對 libomptarget 指標的 `hipMemcpy*` 走 CPU 路徑**（§4.1）。
5. **裝置端 Fortran I/O/STOP 帶進 RPC client** 後 kernel 等待變慢（§4.1）。
6. 已在先前層次記錄：runtime 大小的 private 陣列需要 `LIBOMPTARGET_STACK_SIZE ≥ 49152`；OpenACC `private` 純量在 kernel 外初始化的語意差異；`plold(i,0,jj)` 越界讀。

---

## 6. 現在每步 2.43 s 的組成（rank 0，TIMER；GPU 忙碌約 1.6 s）

| 區塊 | 秒 | 內容 |
|---|---|---|
| NDSL 半拉格朗日平流 | 0.80 | 水平 `monoadvh` 0.40 ＋ `fgnl` 0.19、垂直 `monoadvv` 0.21 ＋ 0.07；其中 RCCL 轉置（PCIe）~0.25、PPM kernel ~0.3、其餘是 ~1.5k 次 launch 的空檔 |
| 物理 `diabat_gpu` | 0.58 | 微物理 `mp_scheme` 0.17、gwdc 0.06、pbl_noah 0.05、深對流 0.05、adjptqintp 0.04、淺對流 0.03、prerrtmg 0.04…；`rrtmg_gpu` 1.47 s 但每 6 步一次（第 1 步 4.1 s 的來源） |
| 譜轉換 | ~0.55 | tranrs 0.13、transr 0.11、trandv 0.10、transr1 0.09、tranuv 0.09、trngra3 0.06、tranrs1 0.05、trngra 0.04：~7k 次 rocFFT exec（31k 個 kernel）＋ ~5k 個 dgemm |
| 其他 | ~0.5 | hdiffu 0.05、siimpl 0.03、mpe 轉置、intgrt 內嵌 kernel（狀態複製、時間濾波） |

kernel 級（prof11，追蹤下）：每步 39.7k 次 launch；GPU 忙 1.60 s，空檔 1.97 s（其中 rocprof 自身每 kernel ~20 µs 佔 0.77 s）。最大的 kernel：RCCL 集合通訊 0.26 s（兩卡之間是 PCIe）、PPM `l102` 75 次 0.11 s、微物理 `saticel_s` 各 15–30 ms、深/淺對流 25 ms。

---

## 7. 試過但無效或退回的（都有 log）

| 實驗 | 結果 |
|---|---|
| libomptarget 池子 threshold 64 GB | 第 2 步 OOM（VRAM 73 → 197 GB） |
| `geps_acc_wait_all` device sync → 逐 stream sync（第一版） | 無差（後來發現變成 26 萬次 sync，改成 dirty list） |
| `HSA_ENABLE_INTERRUPT=0`（兩次） | 無差 |
| `LIBOMPTARGET_AMDGPU_STREAM_BUSYWAIT` | 無差 |
| `OMP_TEAMS_THREAD_LIMIT=64` | 變慢（NDSL、rrtmg） |
| `GPU_MAX_HW_QUEUES=16`＋8 條 FFT thread | 變慢 |
| `HSA_SCRATCH_SINGLE_LIMIT=2 GB` | 無差 |
| `HSA_NO_SCRATCH_RECLAIM=1`（單獨／＋8 GB／＋快取上限 12 GB） | 第 2 步 `OUT_OF_RESOURCES`（scratch 上限，非 VRAM） |
| `HSA_ENABLE_SCRATCH_ASYNC_RECLAIM=0` | 無差 |
| `LIBOMPTARGET_STACK_SIZE=1024` | 卡死在第 1 步（sflx 需要大堆疊） |
| `reduction(max:)` 內層迴圈升級 | **數值錯**（編譯器 bug，§5.1），退回 |
| vectorise pass 用在 NDSL／物理 | 變慢（team 內 fork 開銷），預設跳過 |
| NDSL tile 256／131072 | 無差／變慢，退回 128／32768 |
| PPM 工作陣列改 module 級持久配置 | 快取解掉後不需要，還原 |
| SIGUSR2 取樣 profiler | 讓模式卡死（打斷 HSA 等待），作廢 |

---

## 8. 剩餘機會

**保持逐位元不變的**：
1. NDSL launch 空檔：把 PPM 的 7 個 kernel、垂直 PPM 的 6 個合併（每步 ~1.5k launch × ~40 µs ≈ 0.06 s，收益有限）；PPM `l102` 的 `kstr/kend` 線性搜尋改二分（locs 單調，答案相同；0.1 s）。
2. RCCL 轉置 0.25 s：與計算重疊（目前發完就等）。
3. `has_error` 之類每次呼叫的小 H2D/D2H 合併。
4. rrtmg 1.47 s/次：輻射還有未折的 column 迴圈與 radsw/radlw 的 `pack_size(j2)` 型態。
5. 1-GPU（NPEY=1）路徑：沒有 RCCL 轉置，但單 rank 的 live set 約 2 倍（估 ~220 GB）放不進 192 GB。

**會改捨入（不再逐位元相同，但與 CPU 仍在 1e-11 量級；需要決定）**：
6. rocFFT 逐緯度 exec（每步 7k 次、31k 個 kernel，多為 Bluestein 4–5 個 kernel/次）：reduced grid 每條緯度長度不同、同一 rank 內沒有重複，無法 batch；自寫「一個 launch 處理所有緯度」的 FFT kernel 是唯一結構性解法。
7. Legendre 的 5k 個小 dgemm（每個 m 的形狀不同，無法 batched GEMM）：改自寫 kernel 或 grouped GEMM。

---

## 9. 交付物與重現方式

**程式**
- `cmake/acc2omp.py`：`_spmd_collapse_inner`（`_SPMD_EXT` 檔案→經度上限表、`GEPS_ACC2OMP_SPMD`、`GEPS_ACC2OMP_SPMD_SKIP`）、`_defer_graph_loop_syncs`、`_silence_device_io` 擴充、per-file fixup（massadvx、mpe2d、hdiffu/rayleifr、ndslfv_pack、mod_ndslfv_monoadv）、timer 注入（`GEPS_ACC2OMP_TIMERS`）、`_vectorise_inner_loops` 的 reduction gate。
- `src/rocm/geps_hsa_pool.cc`＋`build_hsa_pool.sh` → `build_rocm_hip/lib/libgeps_hsa_pool.so`（`GEPS_HSA_POOL`、`GEPS_HSA_POOL_CAP_GB`、`GEPS_HSA_POOL_MIN`、`GEPS_HSA_POOL_STATS`）。刻意不進 CMake（改 `src/CMakeLists.txt` 會觸發 reconfigure 陷阱）。
- `src/rocm/hip_compat.cc`：`omp_target_memcpy` D2D、in-place RCCL、dirty-stream `wait_all`、`geps_blas_defer_sync`、FFT exec 不逐次 sync。
- `src/rocm/rfftmlt_loop_gpu.f90`：plan-set LRU（`GEPS_FFT_PLAN_SETS`）、`GEPS_FFT_THREADS`。
- `src/rocm/cyclic_cell_*_gpu.f90`、`vertical_cell_advect_gpu.f90`、`ndslfv_monoadvv_gpu.f90`：tile 大小、tile-local pack、PPM collapse(3)。
- `job/regression.ksh`：`LIBOMPTARGET_STACK_SIZE=65536`、`LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0`、`LD_PRELOAD` 快取、`GEPS_PROF_WRAP`。

**重現**
```bash
# 環境（ROCM_PORT_HANDOFF.md §8）
export GEPS_LIB_ROOT=/mlsteam/workspace/data/geps/GEPS_LIB GEPS_DMSDB=/mlsteam/workspace/data/geps/dmsdb GEPS_COMPILER=rocm
. /usr/share/modules/init/bash && module purge && source $GEPS_LIB_ROOT/env.sh && load_geps_compiler
module use $GEPS_LIB_ROOT/install/modulefiles && module load geps/1.0
# 建置（不要 reconfigure）；acc2omp.py 改動會重翻譯全部檔案，連結約 13 分鐘
cmake --build /mlsteam/workspace/data/geps/GEPS/build_rocm_hip --target tcogfs.x -j 32
src/rocm/build_hsa_pool.sh
# 跑（2 GPU，約 12 分鐘：NNMI ~8 分鐘 ＋ 6 步）
cd job && GEPS_FFT_PLAN_SETS=8 ./TCo383L72_IC_sample_rocm_2gpu > TCo383L72_IC_sample_rocm.log.<tag> 2>&1
grep "surf pres tend rms(GPU)" TCo383L72_IC_sample_rocm.log.<tag>   # 六個值必須逐位元相同
grep "Timing=" TCo383L72_IC_sample_rocm.log.<tag>                  # 每步秒數
```

**注意事項**
- 快取上限 48 GB/rank 是以 6 步 VRAM 持平為依據，長預報（tau ≫ 1）請用 `rocm-smi --showmeminfo vram` 再確認一次。
- 任何新加的 device code 都要重新檢查映像沒有 `__llvm_rpc_client`（§5.5）。
- 任何 team 內的 reduction 都要寫成扁平 kernel 或兩段式（§5.1）。
- 每次 build 後先看 `bin/tcogfs.x` 的時間戳（有兩次建置失敗但舊 binary 照跑）。
