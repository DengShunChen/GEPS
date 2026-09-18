#!/usr/bin/env python3
"""Find leaked OpenMP device mappings from LIBOMPTARGET_INFO=8 output.

Pairs 'Creating new map entry' with 'Removing map entry' by TgtPtrBegin and
reports what is never removed, grouped by the source file libomptarget names.
Streams the file so it works on multi-GB logs.
"""
import re, sys, collections, os

LOG = sys.argv[1] if len(sys.argv) > 1 else "/tmp/omptinfo.log"
CRE = re.compile(r"Creating new map entry with .*?TgtPtrBegin=(0x[0-9a-f]+), Size=(\d+),.*?Name=(\S+)")
REM = re.compile(r"Removing map entry with .*?TgtPtrBegin=(0x[0-9a-f]+), Size=(\d+), Name=(\S+)")

live = {}                      # tgtptr -> (size, name, seq)
created = collections.Counter()      # name -> count
created_b = collections.Counter()    # name -> bytes
removed = collections.Counter()
removed_b = collections.Counter()
n_cre = n_rem = 0
seq = 0
# phase markers straight from the model's own stdout, interleaved in the same file
phase = "init"
phase_of = {}
PHASE = re.compile(r"^\s*(iteration=\s*\d+|forcast begin|in outflds_green for tau|beginning integration)")

with open(LOG, encoding="utf-8", errors="replace") as f:
    for line in f:
        p = PHASE.match(line)
        if p:
            phase = p.group(1).strip()
            continue
        m = CRE.search(line)
        if m:
            seq += 1
            n_cre += 1
            tgt, size, name = m.group(1), int(m.group(2)), os.path.basename(m.group(3))
            live[tgt] = (size, name, seq, phase)
            created[name] += 1
            created_b[name] += size
            continue
        m = REM.search(line)
        if m:
            n_rem += 1
            tgt, size, name = m.group(1), int(m.group(2)), os.path.basename(m.group(3))
            removed[name] += 1
            removed_b[name] += size
            live.pop(tgt, None)

print("log: %s (%.1f MB)" % (LOG, os.path.getsize(LOG)/1e6))
print("Creating new map entry : %d" % n_cre)
print("Removing  map entry    : %d" % n_rem)
print("never removed (live)   : %d entries, %.2f GB"
      % (len(live), sum(v[0] for v in live.values())/1e9))

print("\n=== 未被移除的映射，依原始檔歸戶 (top 15) ===")
by_name = collections.Counter()
cnt_name = collections.Counter()
for size, name, _, _ in live.values():
    by_name[name] += size
    cnt_name[name] += 1
print("  %-40s %10s %14s" % ("source file", "entries", "leaked(GB)"))
for name, b in by_name.most_common(15):
    print("  %-40s %10d %14.3f" % (name, cnt_name[name], b/1e9))

print("\n=== 依模式階段歸戶 ===")
by_phase = collections.Counter(); cnt_phase = collections.Counter()
for size, name, _, ph in live.values():
    by_phase[ph] += size; cnt_phase[ph] += 1
for ph, b in by_phase.most_common():
    print("  %-32s %8d entries  %10.3f GB" % (ph, cnt_phase[ph], b/1e9))

print("\n=== create/remove 收支不平衡 (依檔案，top 15) ===")
print("  %-40s %8s %8s %8s %12s" % ("source file", "create", "remove", "diff", "net(GB)"))
allnames = set(created) | set(removed)
rows = sorted(allnames, key=lambda n: -(created_b[n]-removed_b[n]))
for n in rows[:15]:
    print("  %-40s %8d %8d %8d %12.3f"
          % (n, created[n], removed[n], created[n]-removed[n],
             (created_b[n]-removed_b[n])/1e9))

print("\n=== 最大的 10 筆未移除映射 ===")
big = sorted(live.items(), key=lambda kv: -kv[1][0])[:10]
for tgt, (size, name, s, ph) in big:
    print("  %-18s %12.3f MB  %-28s phase=%s" % (tgt, size/1e6, name, ph))
