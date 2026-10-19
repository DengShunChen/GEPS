#!/usr/bin/env python3
"""Summarise a GEPS run log and compare it with reference.json.

usage: check_run.py <log> [--ref afar-24.3.0|rocm-7.2.2]

Prints whether the run finished, per-step times, and the relative difference of
`surf pres tend rms(GPU)` from the reference. Exit status 1 if the run did not finish
or step 1 differs from the reference by more than 1e-10 relative.
"""
import json, re, statistics, sys, pathlib

def main():
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    ref_name = "afar-24.3.0"
    if "--ref" in args:
        i = args.index("--ref"); ref_name = args[i + 1]; del args[i:i + 2]
    text = pathlib.Path(args[0]).read_text(errors="replace")
    ref = json.loads((pathlib.Path(__file__).parent / "reference.json").read_text())[ref_name]

    sp = [float(m) for m in re.findall(r"surf pres tend rms\(GPU\)\s*=\s*([-+.\dEe]+)", text)]
    tm = [float(m) for m in re.findall(r"Timing=\s*([\d.]+)", text)]
    ended = "PROGRAM CWBGFS HAS ENDED" in text
    errors = sorted(set(re.findall(r"(HSA_STATUS_ERROR_\w+|memory access fault[^\n]*|NaN)", text)))[:5]

    print(f"log        : {args[0]}")
    print(f"finished   : {'yes' if ended else 'NO'}   steps with sptend: {len(sp)}")
    if errors:
        print(f"errors     : {', '.join(errors)}")
    # one Timing= line per step; step 1 prints 0.000, step 2 includes warm-up
    if len(tm) >= 6:
        print(f"s/step     : step 2 {tm[1]:.2f}, mean of steps 3-6 {statistics.mean(tm[2:6]):.2f}"
              f"  (reference {ref['sec_per_step_1h_mean_steps_3_6']})")
    if len(tm) > 6:
        print(f"           : {len(tm)} steps, sum of steps 3-{len(tm)} {sum(tm[2:]):.1f} s, median {statistics.median(tm[1:]):.2f} s,"
              f" max {max(tm):.2f} s  (reference 24 h: {ref['sec_24h_sum_steps_3_144']} s, median {ref['sec_24h_median']} s)")
        # the step after each 6-hourly history output (steps 37, 73, 109) takes ~22 s by design
        slow = [(i + 1, t) for i, t in enumerate(tm) if i > 1 and t > 6.0 and i % 36]
        if slow:
            print(f"           : steps over 6 s besides the history-output steps: {slow[:8]}")

    print(f"sptend vs {ref_name}:")
    for i, v in enumerate(sp[:6]):
        r = ref["sptend_1h"][i]
        print(f"  step {i + 1:3d}  {v:.17g}  ref {r:.17g}  rel {abs(v - r) / abs(r):.1e}")
    for k, r in ref["sptend_24h"].items():
        k = int(k)
        if len(sp) >= k:
            v = sp[k - 1]
            print(f"  step {k:3d}  {v:.17g}  ref {r:.17g}  rel {abs(v - r) / abs(r):.1e}")

    pool = [l for l in text.splitlines() if l.startswith("[geps_hsa_pool] alloc=")]
    if pool:
        print("HSA cache  : " + pool[-1].split("] ", 1)[1])
    ok = ended and sp and abs(sp[0] - ref["sptend_1h"][0]) / ref["sptend_1h"][0] < 1e-10
    print("RESULT     : " + ("OK" if ok else "CHECK"))
    sys.exit(0 if ok else 1)

main()
