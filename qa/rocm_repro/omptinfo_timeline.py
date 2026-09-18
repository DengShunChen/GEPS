#!/usr/bin/env python3
"""Outstanding (created-but-not-removed) mapped bytes over the life of the run."""
import re, os, sys
LOG = sys.argv[1] if len(sys.argv) > 1 else "/tmp/omptinfo.log"
CRE = re.compile(r"Creating new map entry with .*?TgtPtrBegin=(0x[0-9a-f]+), Size=(\d+),")
REM = re.compile(r"Removing map entry with .*?TgtPtrBegin=(0x[0-9a-f]+), Size=(\d+),")
PHASE = re.compile(r"^\s*(iteration=\s*\d+|forcast begin|beginning integration|in outflds_green for tau|in outflds for tau)")
live = {}
outstanding = 0
phase = "init"
events = 0
marks = []
with open(LOG, encoding="utf-8", errors="replace") as f:
    for line in f:
        p = PHASE.match(line)
        if p:
            marks.append((p.group(1).strip(), events, outstanding, len(live)))
            phase = p.group(1).strip()
            continue
        m = CRE.search(line)
        if m:
            events += 1
            live[m.group(1)] = int(m.group(2))
            outstanding += int(m.group(2))
            continue
        m = REM.search(line)
        if m:
            events += 1
            s = live.pop(m.group(1), None)
            if s is not None:
                outstanding -= s
print("=== 每個模式階段開始時，未歸還的映射位元組 ===")
print("  %-30s %12s %14s %10s" % ("phase marker", "map events", "outstanding(GB)", "entries"))
for name, ev, out, n in marks:
    print("  %-30s %12d %14.3f %10d" % (name, ev, out/1e9, n))
print("  %-30s %12d %14.3f %10d" % ("(end of log)", events, outstanding/1e9, len(live)))
