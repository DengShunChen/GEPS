# afar 24.3.0 安裝紀錄（2026-09-29）

AMD 要求用 afar 24.3.0 重測 toolchain PoC 並在 app 上驗證（結果見
`qa/rocm_toolchain_poc/AFAR_24.3.0_RETEST.md`、`ROCM_PORT_HANDOFF.md` §9 (iv-b)）。以下是當天實際執行的步驟，
從工作階段紀錄整理；指令照抄，只拿掉除錯用的輸出。全程約 1.5 小時（下載 10 分鐘、解壓 10 分鐘、OpenMPI 8 分鐘、
其餘函式庫 25 分鐘、app 樹）。

**⚠️ `/opt/afar` 在容器的 overlay rootfs 上，容器重建就會消失**（`/mlsteam/workspace` 才是持久的 NFS）。
`GEPS_LIB_afar/` 與 `GEPS/build_afar_hip/` 在 NFS 上會留著，但它們連結的是 `/opt/afar/current`，重建容器後要先做第 1–2 步。

## 1. 下載與解壓 SDK（`/opt/afar`）

套件是整包 TheRock SDK（flang 24.0.0git、HIP 7.16、HSA 1.21，自帶 libomptarget／rocFFT／rocBLAS／RCCL），
從 AMD 公開的 repo.radeon.com 下載，9.1 GB（9,680,852,441 bytes）：

```bash
mkdir -p /opt/afar && cd /opt/afar
curl -sSfL -o afar.tar.bz2 https://repo.radeon.com/rocm/misc/flang/therock-afar-24.3.0-multiarch-10.1.0-592954c.tar.bz2
```

單執行緒 `tar xjf` 太慢（光是 `tar tjf` 列目錄就 10 分鐘），改用平行的 lbzip2：

```bash
apt-get install -y lbzip2
cd /opt/afar && tar -I lbzip2 -xf afar.tar.bz2
ln -sfn /opt/afar/therock-afar-24.3.0-multiarch-10.1.0-592954c /opt/afar/current
/opt/afar/current/bin/amdflang --version
```

解出 `therock-afar-24.3.0-multiarch-10.1.0-592954c/`（`.info/version` = 10.2.0）。GEPS 的腳本一律用
`/opt/afar/current`。

## 2. 系統套件（不在 `init.sh` 裡）

GEPS_LIB 重編時依序撞到三個缺的工具：

| 失敗 | 原因 | 處理 |
|---|---|---|
| `libpng failed with exit code 127` | 建置腳本呼叫 `gtar`，容器沒有 | `GEPS_LIB_afar/.shim/gtar -> /usr/bin/tar`，建置時放進 PATH |
| `operlib failed with exit code 2` | `Makefile.gnu` 的 `SHELL` 是 ksh | `apt-get install -y ksh` |
| `netcdf failed with exit code 1` | configure 需要 m4 | `apt-get install -y m4` |

```bash
apt-get install -y lbzip2 ksh m4
```

## 3. 用 afar 重編 GEPS_LIB → `GEPS_LIB_afar/`

Fortran `.mod` 綁編譯器版本，所以整包重編，原本的 `GEPS_LIB/`（7.2.2）不動：

```bash
cd /mlsteam/workspace/data/geps
cp -a GEPS_LIB GEPS_LIB_afar && rm -rf GEPS_LIB_afar/install/rocm GEPS_LIB_afar/logs/*
cd GEPS_LIB_afar
```

OpenMPI 先編（`--enable-mpi-fortran=all`，約 8 分鐘；`env -u GEPS_LIB_ROOT` 避免指回原本的 GEPS_LIB）：

```bash
env -u GEPS_LIB_ROOT ROCM_PATH=/opt/afar/current GEPS_OPENMPI_DISTCLEAN=1 GEPS_OPENMPI_JOBS=16 \
  bash src/openmpi/compile.sh > logs/openmpi_afar.out 2>&1
```

`build_all.sh` 第一步會從原始碼 bootstrap cmake，很慢；cmake 與編譯器無關，直接從 7.2.2 那份複製：

```bash
rsync -a ../GEPS_LIB/install/rocm/cmake-3.21.4 install/rocm/
```

其餘函式庫（第 2 節的三個問題處理完後一次通過，約 15 分鐘）：

```bash
mkdir -p .shim && ln -sf /usr/bin/tar .shim/gtar
env -u GEPS_LIB_ROOT PATH=$PWD/.shim:$PATH ROCM_PATH=/opt/afar/current \
  GEPS_OPENMPI_MODULEFILES=$PWD/install/modulefiles \
  bash -c '. /usr/share/modules/init/bash; module purge; ./build_all.sh --compiler rocm \
    libpng w3lib operlib zlib netcdf jasper fftw openblas g2lib dms dms38key mct geps_coupler' \
  > logs/build_all_afar5.out 2>&1
grep -E "^\[OK\]|FATAL" logs/build_all_afar5.out
```

結果全部 `[OK]`：cmake、libpng、operlib、zlib、netcdf、jasper、fftw、openblas、g2lib、dms、dms38key、mct、
geps_coupler（`install/rocm/` 下另有 openmpi-4.1.6、w3lib-2.0.2）。各次嘗試的 log 在 `GEPS_LIB_afar/logs/`
（`build_all_afar{,2..5}.out`、`openmpi_afar.out`）。

GEPS_LIB_afar 的 git 狀態裡 jasper 的 `Makefile`／`config.*` 會顯示修改：那是 configure 把路徑寫成
`GEPS_LIB_afar`，不要 commit。`env.sh`、`src/openmpi/compile.sh` 的修改與 GEPS_LIB 相同（已在 GEPS_LIB commit）。

## 4. 建 GEPS app 樹 `build_afar_hip/`

```bash
cd /mlsteam/workspace/data/geps/GEPS
qa/afar_app/build_afar_hip.sh                 # 預設 build_afar_hip；腳本內含 configure 參數與 rank.mod 順序的繞法
ROCM_PATH=/opt/afar/current CXX=/opt/afar/current/llvm/bin/clang++ \
  src/rocm/build_hsa_pool.sh $PWD/build_afar_hip/lib   # 會另建 lib/hsa_front/（afar 的 libomptarget dlopen HSA）
```

原始碼要先有 `external split`（`src/outsigs.f90`、`src/outsit.f90`；flang 24 的 F2023 `SPLIT` intrinsic），
已在版控裡。

## 5. 跑

```bash
cd /mlsteam/workspace/data/geps/GEPS/job
ROCM_PATH=/opt/afar/current GEPS_LIB_ROOT=/mlsteam/workspace/data/geps/GEPS_LIB_afar \
  GEPS_ROCM_BUILD_DIR=build_afar_hip GEPS_FFT_PLAN_SETS=8 ./TCo383L72_IC_sample_rocm_2gpu
```

`regression.ksh` 偵測到 `lib/hsa_front/` 就自動設 `GEPS_HSA_POOL_CAP_GB=48`、`ROCFFT_RTC_CACHE_PATH`
（`qa/afar_app/rocfft_cache/afar_gfx942.db`，第一次跑約 4 分鐘產生 ~38 MB，之後正常啟動）與 `HSA_ENABLE_SDMA=0`。
長預報加 `GEPS_TAUE=<小時>`。PoC：`cd qa/rocm_toolchain_poc && ROCM_PATH=/opt/afar/current RESULTS=results_afar ./run_all.sh`。

## 容器重建後的最短恢復

`GEPS_LIB_afar/`、`build_afar_hip/`、rocFFT 快取都在 NFS 上，只要補回 `/opt/afar` 與系統套件：

```bash
apt-get install -y lbzip2 ksh m4
mkdir -p /opt/afar && cd /opt/afar
curl -sSfL -o afar.tar.bz2 https://repo.radeon.com/rocm/misc/flang/therock-afar-24.3.0-multiarch-10.1.0-592954c.tar.bz2
tar -I lbzip2 -xf afar.tar.bz2
ln -sfn /opt/afar/therock-afar-24.3.0-multiarch-10.1.0-592954c /opt/afar/current
```

把 `afar.tar.bz2` 留一份在 `/mlsteam/workspace`（NFS 剩約 64 GB）就不用重新下載（也不怕 repo 上的檔案被撤下）。
