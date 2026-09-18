/* amdflang does not emit stack-clash probes. Huge Fortran automatic
 * arrays then subtract RSP past the guard page → SIGSEGV even with
 * ulimit -s unlimited. Touch pages downward so the kernel maps them. */
#define _GNU_SOURCE
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <ucontext.h>
#include <mpi.h>

void geps_touch_stack(void) {
  const size_t bytes = 512ull << 20; /* 512 MiB */
  volatile char anchor = 0;
  volatile char *p = &anchor;
  size_t off;
  for (off = 4096; off < bytes; off += 4096) {
    p[-((ptrdiff_t)off)] = 0;
  }
}

/* #region agent log */
void geps_dbg_log_(int *hyp, int *locid, int *rank, long long *p0, long long *p1,
                   long long *p2) {
  FILE *f;
  volatile char anchor = 0;
  long long stk = (long long)(uintptr_t)&anchor;
  long long ms;
  const char *h;
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  ms = (long long)ts.tv_sec * 1000LL + ts.tv_nsec / 1000000LL;
  switch (*hyp) {
  case 1:
    h = "H1";
    break;
  case 2:
    h = "H2";
    break;
  case 3:
    h = "H3";
    break;
  case 5:
    h = "H5";
    break;
  case 6:
    h = "H6";
    break;
  case 14:
    h = "H14";
    break;
  case 15:
    h = "H15";
    break;
  case 16:
    h = "H16";
    break;
  case 17:
    h = "H17";
    break;
  case 18:
    h = "H18";
    break;
  case 21:
    h = "H21";
    break;
  default:
    h = "H?";
    break;
  }
  f = fopen("/mlsteam/workspace/data/geps/.cursor/debug-78f374.log", "a");
  if (!f)
    return;
  fprintf(f,
          "{\"sessionId\":\"78f374\",\"runId\":\"post-fix\",\"hypothesisId\":\"%s\","
          "\"location\":\"src:%d\",\"message\":\"probe\",\"data\":{\"rank\":%d,"
          "\"p0\":%lld,\"p1\":%lld,\"p2\":%lld,\"stack_addr\":%lld},"
          "\"timestamp\":%lld}\n",
          h, *locid, *rank, (long long)*p0, (long long)*p1, (long long)*p2, stk,
          ms);
  fflush(f);
  fclose(f);
}

static void geps_on_segv(int sig, siginfo_t *si, void *uctx) {
  ucontext_t *uc = (ucontext_t *)uctx;
  FILE *f;
  long long ms;
  struct timespec ts;
  (void)sig;
  clock_gettime(CLOCK_REALTIME, &ts);
  ms = (long long)ts.tv_sec * 1000LL + ts.tv_nsec / 1000000LL;
  f = fopen("/mlsteam/workspace/data/geps/.cursor/debug-78f374.log", "a");
  if (f) {
    fprintf(f,
            "{\"sessionId\":\"78f374\",\"runId\":\"post-fix\",\"hypothesisId\":\"H17\","
            "\"location\":\"src:0\",\"message\":\"segv\",\"data\":{\"addr\":%lld,"
            "\"rip\":%lld,\"rsi\":%lld,\"r8\":%lld,\"rsp\":%lld,\"rbp\":%lld},"
            "\"timestamp\":%lld}\n",
            (long long)(uintptr_t)si->si_addr,
            (long long)uc->uc_mcontext.gregs[REG_RIP],
            (long long)uc->uc_mcontext.gregs[REG_RSI],
            (long long)uc->uc_mcontext.gregs[REG_R8],
            (long long)uc->uc_mcontext.gregs[REG_RSP],
            (long long)uc->uc_mcontext.gregs[REG_RBP], ms);
    fflush(f);
    fclose(f);
  }
}

void geps_dbg_install_(void) {
  struct sigaction sa;
  memset(&sa, 0, sizeof(sa));
  sa.sa_sigaction = geps_on_segv;
  sa.sa_flags = SA_SIGINFO | SA_RESETHAND;
  sigaction(SIGSEGV, &sa, NULL);
}

/* Fortran ABI (all pointers). amdflang bind(C) VALUE on MPI_Allgather
 * made sendcount a pointer bit-pattern → malloc abort. */
void geps_mpi_allgather_dbl_(double *send, double *recv, int *count, int *fcomm,
                             int *ierr) {
  FILE *f;
  long long ms;
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  ms = (long long)ts.tv_sec * 1000LL + ts.tv_nsec / 1000000LL;
  f = fopen("/mlsteam/workspace/data/geps/.cursor/debug-78f374.log", "a");
  if (f) {
    fprintf(f,
            "{\"sessionId\":\"78f374\",\"runId\":\"post-fix\",\"hypothesisId\":\"H21\","
            "\"location\":\"src:263\",\"message\":\"c-allgather\","
            "\"data\":{\"rank\":-1,\"p0\":%d,\"p1\":%lld,\"p2\":%d,\"stack_addr\":%lld},"
            "\"timestamp\":%lld}\n",
            *count, (long long)(uintptr_t)send, *fcomm,
            (long long)(uintptr_t)recv, ms);
    fflush(f);
    fclose(f);
  }
  *ierr = MPI_Allgather(send, *count, MPI_DOUBLE, recv, *count, MPI_DOUBLE,
                        MPI_Comm_f2c(*fcomm));
}
/* #endregion */
