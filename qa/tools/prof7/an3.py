import sqlite3,collections
c=sqlite3.connect('r_3582039_results.db'); cur=c.cursor()
t0=cur.execute("select min(start) from rocpd_kernel_dispatch").fetchone()[0]
lo,hi=t0+712*10**9, t0+724*10**9
# threads
tid_names={}
for (tid,name) in cur.execute("select id,name from rocpd_info_thread"): tid_names[tid]=name
print("threads:", len(tid_names))
agg=collections.defaultdict(lambda:[0,0.0]); bythr=collections.Counter()
for (name,tid,s,e) in cur.execute("select S.string, R.tid, R.start, R.end from rocpd_region R join rocpd_string S on S.id=R.name_id where R.start>=? and R.start<?", (lo,hi)):
    a=agg[name]; a[0]+=1; a[1]+=(e-s)/1e9; bythr[(tid,name)]+=(e-s)/1e9
print("== HSA API in one step (12 s window), by total time")
for k,(n,t) in sorted(agg.items(), key=lambda x:-x[1][1])[:16]: print(f"{n:8d} {t:8.3f}s avg {t/n*1e6:8.1f}us  {k}")
print("== signal wait by thread")
for (tid,name),t in sorted(bythr.items(), key=lambda x:-x[1])[:8]: print(f"  tid {tid} {tid_names.get(tid)} {t:.3f}s {name}")
# allocation sizes: memory_allocate table
print("== allocations in window by size class")
sz=collections.defaultdict(lambda:[0,0.0])
for (s,e,size,typ) in cur.execute("select start,end,size,type from rocpd_memory_allocate where start>=? and start<?", (lo,hi)):
    cls = '<1MB' if size<2**20 else ('1-16MB' if size<16*2**20 else ('16-128MB' if size<128*2**20 else ('128MB-1GB' if size<2**30 else '>=1GB')))
    a=sz[(typ,cls)]; a[0]+=1; a[1]+=(e-s)/1e9
for k,(n,t) in sorted(sz.items(), key=lambda x:-x[1][1]): print(f"  {k}: {n} calls {t:.3f}s")
print(cur.execute("select sql from sqlite_master where name='rocpd_memory_allocate'").fetchone()[0][:600])
