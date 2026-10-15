import sqlite3,collections,sys
c=sqlite3.connect('r_3582039_results.db'); cur=c.cursor()
t0=cur.execute("select min(start) from rocpd_kernel_dispatch").fetchone()[0]
# per-second kernel launch histogram to find step windows
bins=collections.Counter()
for (s,) in cur.execute("select start from rocpd_kernel_dispatch"): bins[(s-t0)//10**9]+=1
secs=sorted(bins)
# print seconds with >5000 launches and gaps to locate the 6 steps at the end
hi=[b for b in secs if bins[b]>3000]
print("busy seconds:", hi[:5], '...', hi[-40:])
