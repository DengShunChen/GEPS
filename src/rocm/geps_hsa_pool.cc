// geps_hsa_pool: LD_PRELOAD cache in front of hsa_amd_memory_pool_allocate/free.
//
// 2026-09-19. With LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0 (needed because
// libomptarget's own pool only reuses exact sizes and grew to OOM on the
// radiation blocks, see ROCM_PORT_HANDOFF.md §6) every OpenACC-style
// `data create` in the model turns into a raw HSA allocation. rocprofv3
// (prof7) measured 3.6k allocs + 3.6k frees per step on rank 0, 4.1 s of an
// 11.8 s step: a free costs 0.4-5 ms whatever its size, a >=1 GB alloc
// 5-230 ms.
//
// This shim keeps freed blocks in size-class free lists (sizes rounded up to
// 3-mantissa-bit classes, so at most 12.5% slack) keyed by (pool, flags), and
// hands them back on the next allocation of the same class. Small blocks
// (< GEPS_HSA_POOL_MIN, default 256 KB) pass straight through. The cache is
// capped (GEPS_HSA_POOL_CAP_GB, default 48 GB per process); above the cap the
// least recently freed blocks are really released, and an allocation failure
// flushes the cache and retries. Memory handed back is not zeroed - the same
// as raw HSA in practice, and the model does not rely on it (six-step sptend
// stays bit-identical, see the handoff).
//
// GEPS_HSA_POOL=0 disables the cache (pass-through). GEPS_HSA_POOL_STATS=1
// prints counters at exit.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <list>
#include <map>
#include <mutex>
#include <unordered_map>

#include <hsa/hsa.h>
#include <hsa/hsa_ext_amd.h>

namespace {

using alloc_fn_t = hsa_status_t (*)(hsa_amd_memory_pool_t, size_t, uint32_t,
                                    void **);
using free_fn_t = hsa_status_t (*)(void *);

alloc_fn_t real_alloc = nullptr;
free_fn_t real_free = nullptr;

struct Key {
  uint64_t pool;
  uint32_t flags;
  size_t cls; // rounded size
  bool operator<(const Key &o) const {
    if (pool != o.pool)
      return pool < o.pool;
    if (flags != o.flags)
      return flags < o.flags;
    return cls < o.cls;
  }
};

struct Block {
  void *ptr;
  Key key;
};

struct Live {
  Key key;
};

std::mutex g_mu;
bool g_init = false;
bool g_enabled = true;
bool g_stats = false;
size_t g_min = 256 * 1024;
size_t g_cap = (size_t)48 << 30;
size_t g_cached = 0; // bytes in free lists
std::unordered_map<void *, Live> g_live; // every block we handed out (cached class)
// free lists: key -> list of ptrs; plus a global LRU of (ptr,key) for eviction
std::map<Key, std::list<void *>> g_free;
std::list<Block> g_lru; // front = most recently freed
std::unordered_map<void *, std::list<Block>::iterator> g_lru_pos;

// counters
unsigned long n_alloc = 0, n_hit = 0, n_miss = 0, n_free = 0, n_cached_free = 0,
              n_evict = 0, n_flush = 0, n_pass = 0;

void init_once() {
  if (g_init)
    return;
  g_init = true;
  real_alloc = (alloc_fn_t)dlsym(RTLD_NEXT, "hsa_amd_memory_pool_allocate");
  real_free = (free_fn_t)dlsym(RTLD_NEXT, "hsa_amd_memory_pool_free");
  const char *e = std::getenv("GEPS_HSA_POOL");
  if (e && e[0] == '0')
    g_enabled = false;
  if ((e = std::getenv("GEPS_HSA_POOL_MIN")))
    g_min = (size_t)std::atoll(e);
  if ((e = std::getenv("GEPS_HSA_POOL_CAP_GB")))
    g_cap = (size_t)std::atoll(e) << 30;
  if ((e = std::getenv("GEPS_HSA_POOL_STATS")) && e[0] == '1')
    g_stats = true;
}

// round up to 3 mantissa bits: 1.000, 1.125, ..., 1.875 x 2^k
size_t round_class(size_t n) {
  if (n <= g_min)
    return n;
  int k = 63 - __builtin_clzll((unsigned long long)n); // floor(log2 n)
  size_t unit = (size_t)1 << (k - 3);                  // 2^k / 8
  return (n + unit - 1) / unit * unit;
}

void lru_erase(void *p) {
  auto it = g_lru_pos.find(p);
  if (it == g_lru_pos.end())
    return;
  g_lru.erase(it->second);
  g_lru_pos.erase(it);
}

// really release one cached block (g_mu held)
void release_block(const Block &b) {
  auto fl = g_free.find(b.key);
  if (fl != g_free.end()) {
    fl->second.remove(b.ptr);
    if (fl->second.empty())
      g_free.erase(fl);
  }
  g_cached -= b.key.cls;
  g_live.erase(b.ptr);
  (void)real_free(b.ptr);
}

void evict_to(size_t target) {
  while (g_cached > target && !g_lru.empty()) {
    Block b = g_lru.back();
    g_lru.pop_back();
    g_lru_pos.erase(b.ptr);
    release_block(b);
    n_evict++;
  }
}

void flush_all() {
  n_flush++;
  evict_to(0);
}

void print_stats() {
  std::fprintf(stderr,
               "[geps_hsa_pool] alloc=%lu hit=%lu miss=%lu pass=%lu free=%lu "
               "cached_free=%lu evict=%lu flush=%lu cached_now=%.1f GB\n",
               n_alloc, n_hit, n_miss, n_pass, n_free, n_cached_free, n_evict,
               n_flush, g_cached / 1073741824.0);
}

struct AtExit {
  ~AtExit() {
    if (g_stats)
      print_stats();
  }
} g_atexit;

} // namespace

extern "C" {

hsa_status_t hsa_amd_memory_pool_allocate(hsa_amd_memory_pool_t pool,
                                          size_t size, uint32_t flags,
                                          void **ptr) {
  std::lock_guard<std::mutex> lock(g_mu);
  init_once();
  n_alloc++;
  if (!g_enabled || size < g_min) {
    n_pass++;
    return real_alloc(pool, size, flags, ptr);
  }
  Key key{pool.handle, flags, round_class(size)};
  auto fl = g_free.find(key);
  if (fl != g_free.end() && !fl->second.empty()) {
    void *p = fl->second.front();
    fl->second.pop_front();
    if (fl->second.empty())
      g_free.erase(fl);
    lru_erase(p);
    g_cached -= key.cls;
    n_hit++;
    *ptr = p;
    return HSA_STATUS_SUCCESS;
  }
  n_miss++;
  hsa_status_t st = real_alloc(pool, key.cls, flags, ptr);
  if (st != HSA_STATUS_SUCCESS && g_cached > 0) {
    flush_all();
    st = real_alloc(pool, key.cls, flags, ptr);
  }
  if (st == HSA_STATUS_SUCCESS)
    g_live[*ptr] = Live{key};
  return st;
}

hsa_status_t hsa_amd_memory_pool_free(void *p) {
  std::lock_guard<std::mutex> lock(g_mu);
  init_once();
  n_free++;
  auto it = g_live.find(p);
  if (it == g_live.end())
    return real_free(p); // pass-through block (or foreign)
  if (g_lru_pos.count(p))
    return HSA_STATUS_SUCCESS; // double free: ignore
  Key key = it->second.key;
  g_free[key].push_front(p);
  g_lru.push_front(Block{p, key});
  g_lru_pos[p] = g_lru.begin();
  g_cached += key.cls;
  n_cached_free++;
  if (g_cached > g_cap)
    evict_to(g_cap * 3 / 4);
  return HSA_STATUS_SUCCESS;
}

} // extern "C"
