# ROCm 移植的除錯／量測小工具

`ROCM_PORT_HANDOFF.md` 與 `ROCM_OPTIMIZATION_REPORT.md` 引用的獨立小程式與腳本。原本放在
Claude 工作階段的 scratchpad（`/tmp`，容器重建會消失），2026-10-15 搬到這裡；內容未改。
只收原始碼與腳本；編譯產物、rocprofv3 trace 資料庫、取樣結果與備份檔沒有保存，文件裡對它們的
引用只是紀錄當時的量測。

多數程式是一次性的實驗：路徑、GPU 編號、PID 寫死在裡面，要重跑請先看過再改。
給 AMD 的正式重現程式在 `qa/rocm_toolchain_poc/`（整理過、各附 `run.sh`）。

| 檔案 | 用途 | 文件出處 |
|---|---|---|
| `redmax/t.f90`、`s.f90`、`t2.f90`、`t3.f90` | amdflang 7.2.2 巢狀 `parallel do reduction` 回 0 的獨立重現 | handoff「效能層 P1」（02:35 根因）、report §5 第 1 項 |
| `privarr/t*.f90` | runtime 大小 private 陣列在裝置堆疊上互相覆寫（`LIBOMPTARGET_STACK_SIZE`） | handoff 第 18 層（獨立重現） |
| `pool.f90` | libomptarget memory manager 只重用同大小區塊 | handoff 第 14 層 |
| `launchbench/lb*.f90`、`run.sh` | kernel launch／`target update`／`enter/exit data` 成本微基準（`lb5` 為映射引數數量） | handoff「效能層 P1」、report §2 |
| `fftbench/fb.cc` | rocFFT/hipFFT exec 的 host 端成本、1 vs N host thread | report §2 |
| `fftplan.cc`、`fftplan2.cc`、`fftplan3.cc` | rocFFT plan 建立成本（plan-set 預算） | handoff 第 13 層 |
| `gtsv_test.cc` | hipSPARSE gtsv 速度（m=72、batch=1000） | handoff 第 18 層 |
| `eigrepro2.cc` | `geps_rocm_dsyevj_dev` 呼叫序列的獨立重現（hipSOLVER） | handoff 第 8 層 |
| `mpisize.f90` | 自編 OpenMPI 的 `MPI_DOUBLE_PRECISION` 大小（16 bytes 問題） | handoff「效能層 P1」 |
| `interpose/memcpy_bt.cc` | `LD_PRELOAD` 攔截 ≥16 MB `hipMemcpyAsync` 並印 backtrace（編成 `libmemcpy_bt.so`） | handoff「2026-09-17 18:55」節 |
| `sample_bt.sh`、`sample_bt2.sh`、`sample_bt3.sh` | 預報階段用 rocgdb 對 rank 0 取 host backtrace（各版改 PID 抓法） | handoff「2026-09-17 18:55」節 host backtrace 直方圖（`sample_bt3.sh`：「效能層 P1」直方圖 #3） |
| `sample_bt4.sh` | 同上的通用版：`pgrep -x tcogfs.x` 跳過 zombie、樣本數參數 | handoff §9 (iv-b) rrtmg 尖峰 |
| `kfdmon.sh` | 每秒記錄 KFD `evicted_ms`／`faults`／`page_in`／`page_out`（目錄名是主機 PID；GPU id `29122`、`43483` 是這台機器的，換機要改） | handoff §9 (iv-b)、report §5 |
| `prof7/an1.py`–`an3.py`、`prof8/an8.py` | rocprofv3 sqlite 結果的分析腳本（每步 kernel／HSA／HIP 時間分解） | handoff「效能層 P1」prof7、prof8 |
