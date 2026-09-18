#!/usr/bin/env python3
"""Per-wavenumber CPU-vs-GPU comparison of tranrs' post-transpose arrays.

For each `DBGMARK tendget_tranrs_N` the analyser takes the NEXT `DBGTRC call=K`
(so the record belongs to that call site, not to getrdy/intgrt which share
tranrs), then every `DBGTRM K mf wcc ws` line and the `DBGTRG call=K` line.
CPU records are gathered over all ranks by the probe itself, so both sides are
global and compared 1:1 (no x32).
"""
import re, sys, math
from collections import OrderedDict

def parse(path):
    """Records are keyed by the ORDER of DBGTRC lines, not the printed counter:
    in the GPU log the CPU tranrs (getrdy) and tranrs_gpu_cuda_graph each keep
    their own counter, so both print `call= 1`."""
    marks = {}          # marker -> block index
    per = {}            # block -> {mf: (wcc, wss)}
    glob = {}           # block -> (gwcc, gwss_lgem, gwss_all)
    pending = None
    blk = 0
    buf = ""
    for raw in open(path, errors="replace"):
        ln = buf + raw.rstrip("\n")
        buf = ""
        m = re.search(r"DBGMARK (tendget_tranrs_\d)", ln)
        if m:
            pending = m.group(1); continue
        if re.search(r"DBGTRC call=", ln):
            blk += 1
            if pending:
                marks[pending] = blk; pending = None
            continue
        m = re.search(r"DBGTRM\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)", ln)
        if m:
            per.setdefault(blk, {})[int(m.group(2))] = (float(m.group(3)), float(m.group(4)))
            continue
        if "DBGTRG call=" in ln and "gwss_all=" not in ln:
            buf = ln            # Fortran wrapped the record onto the next line
            continue
        m = re.search(r"DBGTRG call=\s*(\d+)\s+gwcc_fk=\s*(\S+)\s+gwss_lgem=\s*(\S+)\s+gwss\s*_?all=\s*(\S+)", ln.replace("gwss\n", "gwss"))
        if m:
            glob[blk] = tuple(float(x) for x in m.groups()[1:])
    return marks, per, glob

def ratio(a, b):
    if b == 0: return float('inf') if a else float('nan')
    return a / b

cpu, gpu = sys.argv[1], sys.argv[2]
cm, cp, cg = parse(cpu)
gm, gp, gg = parse(gpu)
print("CPU markers:", cm, " GPU markers:", gm)
for site in sorted(set(cm) & set(gm)):
    kc, kg = cm[site], gm[site]
    print(f"\n=== {site}: CPU call {kc}, GPU call {kg} ===")
    if kc in cg and kg in gg:
        for name, i in (("gwcc_fk", 0), ("gwss_lgem", 1), ("gwss_all", 2)):
            print(f"  {name:10s} CPU={cg[kc][i]:.4e}  GPU={gg[kg][i]:.4e}  GPU/CPU={ratio(gg[kg][i], cg[kc][i]):.4e}")
    c, g = cp.get(kc, {}), gp.get(kg, {})
    mfs = sorted(set(c) | set(g))
    print(f"  mf on CPU: {len(c)}, on GPU: {len(g)}, union {len(mfs)}")
    rows = []
    for mf in mfs:
        cw, cs = c.get(mf, (float('nan'),)*2)
        gw, gs = g.get(mf, (float('nan'),)*2)
        rows.append((mf, cw, gw, ratio(gw, cw), cs, gs, ratio(gs, cs)))
    # histogram of the wss ratio in decades
    hist = {}
    for r in rows:
        if r[6] == r[6] and r[6] not in (float('inf'),) and r[6] > 0:
            d = int(math.floor(math.log10(r[6])))
            hist[d] = hist.get(d, 0) + 1
    print("  wss ratio decades (log10 floor -> count):", dict(sorted(hist.items())))
    histw = {}
    for r in rows:
        if r[3] == r[3] and r[3] not in (float('inf'),) and r[3] > 0:
            d = int(math.floor(math.log10(r[3])))
            histw[d] = histw.get(d, 0) + 1
    print("  wcc ratio decades:", dict(sorted(histw.items())))
    print(f"  {'mf':>4} {'wcc CPU':>11} {'wcc GPU':>11} {'ratio':>9} | {'wss CPU':>11} {'wss GPU':>11} {'ratio':>9}")
    step = max(1, len(rows) // 40)
    show = rows[::step]
    # always include the worst wss ratios too
    worst = sorted((r for r in rows if r[6] == r[6]), key=lambda r: -abs(math.log10(r[6])) if r[6] > 0 else 0)[:8]
    for r in sorted(set(show) | set(worst), key=lambda r: r[0]):
        print(f"  {r[0]:4d} {r[1]:11.3e} {r[2]:11.3e} {r[3]:9.3e} | {r[4]:11.3e} {r[5]:11.3e} {r[6]:9.3e}")
