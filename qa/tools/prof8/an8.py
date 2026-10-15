import sqlite3,collections,re,sys
db=sys.argv[1] if len(sys.argv)>1 else 'r_3908792_results.db'
c=sqlite3.connect(db); cur=c.cursor()
# step markers: a kernel that runs once per step
marks=[s for (s,) in cur.execute("select K.start from rocpd_kernel_dispatch K join rocpd_info_kernel_symbol S on S.id=K.kernel_id where S.display_name like '%diabat_gpu__l2123%' order by K.start")]
print("markers:", len(marks))
lo,hi=marks[-3],marks[-2]
print("step window %.2f s"%((hi-lo)/1e9))
main_tid=cur.execute("select tid from rocpd_region R group by tid order by count(*) desc limit 1").fetchone()[0]
print("main tid", main_tid)
# kernels in window
kt=collections.defaultdict(lambda:[0,0.0]); busy=0.0; nk=0
for (name,s,e) in cur.execute("select S.display_name, K.start, K.end from rocpd_kernel_dispatch K join rocpd_info_kernel_symbol S on S.id=K.kernel_id where K.start>=? and K.start<?", (lo,hi)):
    m=re.search(r'__omp_offloading_[0-9a-f_]+_(.*?)(_l\d+)',name)
    key=(m.group(1)[-40:]+m.group(2)) if m else name[:50]
    a=kt[key]; a[0]+=1; a[1]+=(e-s)/1e9; busy+=(e-s)/1e9; nk+=1
print(f"kernels: {nk} launches, busy {busy:.2f} s")
for k,(n,t) in sorted(kt.items(), key=lambda x:-x[1][1])[:25]: print(f"   {n:6d} {t:.3f}s avg {t/n*1e6:7.1f}us {k}")
# regions by category on main thread
cat=collections.defaultdict(lambda:[0,0.0])
for (catname,name,s,e) in cur.execute("""select (select string from rocpd_string RS where RS.id=E.category_id), S.string, R.start, R.end
   from rocpd_region R join rocpd_event E on E.id=R.event_id join rocpd_string S on S.id=R.name_id
   where R.tid=? and R.start>=? and R.start<?""", (main_tid,lo,hi)):
    a=cat[(catname,name)]; a[0]+=1; a[1]+=(e-s)/1e9
tot=collections.Counter()
for (cn,nm),(n,t) in cat.items(): tot[cn]+=t
print("main-thread API time by category:", {k:round(v,2) for k,v in tot.items()})
for k,(n,t) in sorted(cat.items(), key=lambda x:-x[1][1])[:30]: print(f"   {n:7d} {t:.3f}s avg {t/n*1e6:7.1f}us {k[0]}:{k[1]}")
