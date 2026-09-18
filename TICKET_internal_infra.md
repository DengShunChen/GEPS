# 工單：GEPS on AMD MI300X — 機器設定回報與後續請求

日期：2026-09-10
提出：GEPS / TCoGFS ROCm 移植
對象：機房 / 機器管理單位

---

## 0. 先講結論：上次的請求已經解決問題，非常感謝

上一版工單請求「借一張未分割（SPX）的 MI300X 測試時段」，用來驗證我們卡了兩天的
崩潰是不是 DPX 分區造成的。

**貴單位已把機器調整成 2 張實體 MI300X（NPS1 / SPX，未分割），我們立刻重測，結論是：
問題確實出在分區設定，換卡之後崩潰完全消失。** 這個請求可以結案。

### 對照證據（唯一變數就是分區模式）

| | 調整前 | 調整後 |
|---|---|---|
| 設定 | 1 張 MI300X，DPX 分區 | 2 張實體 MI300X，NPS1 / SPX |
| 執行檔 | `build_rocm_hip/bin/tcogfs.x`（2026-09-09 03:58 建置） | **完全相同，沒有重新編譯** |
| 輸入 / job 腳本 / 環境變數 | — | **完全相同** |
| 結果 | 進入第 2 次 NNMI 迭代就 `HSA_STATUS_ERROR_OUT_OF_RESOURCES` | **完全沒有出現（整個跑批 0 次）**，順利跨過該點 |
| 跑批時間 | — | 06:11:21 → 07:06:48（約 55 分鐘） |

之後撞到的是另一個**我們自己程式的數值問題**，跟機器無關，由我們這邊處理。

---

## 1. 一個技術上的意外發現（供貴單位參考，非請求）

我們原本推測「DPX 分區把 SDMA queue 之類的固定配額壓小了」。**實測發現這個推測是錯的**
——SPX 之後，KFD 回報的 per-node 配額數字跟 DPX 時期**完全一樣**：

| `/sys/class/kfd/kfd/topology/nodes/N/properties` | DPX（舊記錄） | SPX（現在） |
|---|---|---|
| `num_sdma_engines` | 2 | 2 |
| `num_sdma_queues_per_engine` | 8 | 8 |
| `num_gws` | 64 | 64 |
| VRAM total | 206,141,652,992 B | 206,141,652,992 B |

真正有差的是**執行期實際配置得到的量**（用 `rocm-smi` 每 2 秒輪詢 VRAM）：

| | DPX（舊） | SPX（現在） |
|---|---|---|
| 全程 VRAM 用量 | **卡在 ~46 GB 不動**，然後配置失敗 | 36 GB → 77 GB → 83 GB → **峰值 135 GB** |

也就是說，**分區模式下有某個「可配置量」的天花板，但它不反映在 KFD 的可見屬性上**。
如果貴單位之後還要用分區模式跑 HPC / OpenMP offload 類的 workload，這點可能值得留意。
我們沒有進一步追查的權限與必要，僅供參考。

---

## 2. 請求

### 請求 1（高）：請維持目前的 SPX 設定，不要改回分區模式

我們的移植工作接下來會持續需要這個環境。這 2 張卡目前是**專屬給我們**使用，
單次跑批約 1 小時、峰值會用到 135 GB VRAM，維持專屬與 SPX 設定對我們的進度很重要。

如果之後因為其他使用者需求必須改回分區模式、或改成與他人共用，**請提前告知**，
我們需要重新評估與調整排程。

### 請求 2（中）：容器基礎 image 補上幾個系統套件

**問題**：容器 rootfs 不持久。只有 `/mlsteam/workspace`（NFS）跨 session 保留，
`/`、`/usr`、`/lib` 等會在容器重建時重置回底層 image。以下套件不在底層 image 裡，
每次新 session 都要重裝，否則我們的執行檔會因為動態連結失敗而跑不起來：

**(a) 執行期**（缺了模式跑不起來，動態連結失敗）：
```
libtirpc3t64
libhwloc15
libevent-core-2.1-7t64
libevent-pthreads-2.1-7t64
libdw1t64                  # ROCm 的 rocprofv3 需要，缺了會報 libdw.so.1 not found
```

**(b) 編譯期**（缺了**完全無法重新編譯**，link 階段就失敗）：
```
libtirpc-dev               # 提供 libtirpc.so，否則 link 報 unable to find library -ltirpc
libcurl4-openssl-dev       # 提供 libcurl.so，否則 link 報 unable to find library -lcurl
```

> (b) 是 2026-09-10 這次才發現的。執行期套件 `libtirpc3t64` 只提供 `libtirpc.so.3`，
> **連結需要的是 dev 套件提供的 `libtirpc.so` symlink**，兩者不同。這代表目前的容器
> 即使裝齊 (a) 也只能「跑」不能「重編」——我們這次要加除錯探針時才踩到，
> 而且失敗的 link 會先把既有的執行檔刪掉，等於暫時失去可用的 binary。

**請求**：把這 7 個套件加進基礎 image。它們都是 Ubuntu 24.04 官方套件、體積很小。

**目前的暫時作法**（每次新 session 手動執行）：
```bash
apt-get install -y libtirpc3t64 libhwloc15 libevent-core-2.1-7t64 \
                   libevent-pthreads-2.1-7t64 libdw1t64 \
                   libtirpc-dev libcurl4-openssl-dev && ldconfig
```

**附帶請求（低）**：如果容器 rootfs 有機會改成持久化，那是更根本的解法，
可以省掉這類反覆處理。若技術上有困難，補進 image 即可。

---

## 3. 不再需要的舊請求

- ~~借一張未分割（SPX）的 MI300X 測試時段~~ → **已達成，見第 0 節**
- ~~索取分區模式下的固定資源配額文件~~ → **已不需要**。實測顯示 KFD 的可見配額
  在兩種模式下相同（見第 1 節），所以那份文件回答不了我們原本的問題，就不麻煩了。
