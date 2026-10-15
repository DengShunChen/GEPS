import sqlite3,collections,sys
c=sqlite3.connect('r_3582039_results.db'); cur=c.cursor()
t0=cur.execute("select min(start) from rocpd_kernel_dispatch").fetchone()[0]
kb=collections.Counter(); kbusy=collections.Counter()
for (s,e) in cur.execute("select start,end from rocpd_kernel_dispatch"):
    b=(s-t0)//10**9; kb[b]+=1; kbusy[b]+=(e-s)/1e9
# regions (HSA API) per second: total duration by name
rb=collections.defaultdict(collections.Counter); rc=collections.defaultdict(collections.Counter)
for (name,s,e) in cur.execute("select S.string, R.start, R.end from rocpd_region R join rocpd_string S on S.id=R.name_id"):
    b=(s-t0)//10**9
    if 700<=b<=745:
        rb[b][name]+=(e-s)/1e9; rc[b][name]+=1
for b in range(700,745):
    tot=sum(rb[b].values())
    print(f"=== sec {b}: launches {kb[b]:6d} busy {kbusy[b]:.2f} | api calls {sum(rc[b].values()):6d} api time {tot:.2f}")
    for k,v in rb[b].most_common(5): print(f"     {rc[b][k]:6d} {v:.3f}s {k}")
