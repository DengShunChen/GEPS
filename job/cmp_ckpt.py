#!/usr/bin/env python3
"""Compare DBGGRID/DBGSPEC checkpoints between a CPU and a GPU log.
Takes the FIRST occurrence of each tag on each side (= integration step 1;
tags repeat every step). Both sides print GLOBAL used-region sums, so the
expected ratio is 1. Usage: cmp_ckpt.py cpu.log gpu.log [tag-prefix]"""
import re, sys
def parse(path):
    d = {}
    for ln in open(path, errors='replace'):
        m = re.search(r'DBG(GRID|SPEC)\s+(\S+)\s+ssq=\s*(\S+)\s+nonfinite=\s*(\d+)', ln)
        if m and m.group(2) not in d:
            d[m.group(2)] = (float(m.group(3)), int(m.group(4)))
    return d
c, g = parse(sys.argv[1]), parse(sys.argv[2])
pre = sys.argv[3] if len(sys.argv) > 3 else ''
print(f"{'tag':26s} {'CPU':>14s} {'GPU':>14s} {'GPU/CPU':>10s}")
for t in [k for k in c if k.startswith(pre)]:
    if t not in g: print(f"{t:26s} {c[t][0]:14.6e} {'(missing)':>14s}"); continue
    r = g[t][0]/c[t][0] if c[t][0] else float('nan')
    flag = '' if abs(r-1) < 1e-3 else ('  <<' if abs(r-1) > 1e-2 else '  ~')
    nf = f"  nonfinite GPU={g[t][1]}" if g[t][1] else ''
    print(f"{t:26s} {c[t][0]:14.6e} {g[t][0]:14.6e} {r:10.4e}{flag}{nf}")
