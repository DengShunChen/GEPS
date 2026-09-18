#!/usr/bin/env python3
"""Translate NVIDIA OpenACC + CUDA Fortran extensions to OpenMP offload.

Used only for GEPS_COMPILER=rocm + USE_HIP. Original sources are not modified;
CMake writes translated copies into the build tree.
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

# Diagnostic probes (DBGMAP/DBGFREE/DBGTRM/DBGHD/DBGMP/DBGSFLX/DBGLFC, the
# intgrt/diabat/pbl checkpoint tables, nnmi/tendget/tranrs value probes) are
# opt-in: GEPS_ACC2OMP_PROBES=1. Default off = the clean production translation
# (2026-09-17; with them on a step costs ~750 s instead of seconds).
PROBES = os.environ.get("GEPS_ACC2OMP_PROBES", "0") == "1"

SENT = re.compile(r"^(\s*)!\$acc(\s|&)(.*)$", re.IGNORECASE)
CONT = re.compile(r"&\s*$")


def _join_continuations(lines: list[str]) -> list[tuple[str, bool]]:
    """Return (text, is_acc_directive) after joining !$acc continuations."""
    out: list[tuple[str, bool]] = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].rstrip("\n")
        m = SENT.match(line)
        if not m or line.lstrip().startswith("!!$"):
            out.append((line, False))
            i += 1
            continue
        buf = line
        while True:
            cont_self = bool(CONT.search(buf))
            nxt_is_acc = (
                i + 1 < n
                and SENT.match(lines[i + 1])
                and "&" in lines[i + 1][:12]
            )
            if not cont_self and not nxt_is_acc:
                break
            if i + 1 < n and lines[i + 1].lstrip().startswith("#"):
                # Keep `!$acc ... &` / `!$acc&` around `#ifdef`; never splice cpp into the directive.
                break
            i += 1
            nxt = lines[i].rstrip("\n")
            nm = SENT.match(nxt)
            if nm:
                rest = nm.group(3).lstrip()
                if rest.startswith("&"):
                    rest = rest[1:].lstrip()
                buf = CONT.sub("", buf).rstrip() + " " + rest
            else:
                buf = CONT.sub("", buf).rstrip() + " " + nxt.strip()
        out.append((buf, True))
        i += 1
    return out


def _strip_clause(text: str, name: str) -> str:
    # present(foo, bar(1:n))  /  async(id)  /  async
    pat = rf"\b{name}\s*(\([^()]*(?:\([^()]*\)[^()]*)*\))?"
    return re.sub(pat, " ", text, flags=re.IGNORECASE)


def _strip_paren_clause(text: str, name: str) -> str:
    """Strip `name(...)` only (do not eat prefixes of longer clause names)."""
    pat = rf"\b{name}\s*\([^()]*(?:\([^()]*\)[^()]*)*\)"
    return re.sub(pat, " ", text, flags=re.IGNORECASE)


def _rewrite_clauses(body: str, *, enter_data: bool = False) -> str:
    # Longer names first: present_or_copyin contains copyin, present_or_copy contains copy.
    # Structured target/data: copyin→tofrom so device writes are allowed
    # (LLVM AMDGPU maps `to` as read-only). enter data cannot use tofrom
    # (amdflang FIR) and must not remap already-present module arrays with
    # alloc+update (that GPU-faulted in initial_gpu). enter data copyin→map(to:).
    copyin_map = "map(to:" if enter_data else "map(tofrom:"
    copy_map = "map(to:" if enter_data else "map(tofrom:"
    for pat, repl in (
        (r"\bpresent_or_copyin\s*\(", copyin_map),
        (r"\bpresent_or_copyout\s*\(", "map(from:"),
        (r"\bpresent_or_copy\s*\(", copy_map),
        (r"\bpresent_or_create\s*\(", "map(alloc:"),
        (r"\bcopyin\s*\(", copyin_map),
        (r"\bcopyout\s*\(", "map(from:"),
        (r"\bcopy\s*\(", copy_map),
        (r"\bcreate\s*\(", "map(alloc:"),
        # OpenACC `exit data delete` DECREMENTS the dynamic reference count
        # (the buffer stays present while an outer enter data still holds it);
        # OpenMP `map(delete:)` forces the count to zero. The OpenMP equivalent
        # is `release`. With `delete`, initial_gpu's exit data unmapped the
        # module-level up/vp/ttp/rvor (mod_grid.f90:82), so intgrt_gpu's
        # use_device_addr(up) handed hipMemcpy a host pointer, the copy
        # failed silently and (ut-up)/dta came out 3100x (2026-09-16).
        # `finalize` is not used anywhere in src, so delete is never wanted.
        (r"\bdelete\s*\(", "map(release:"),
        (r"\bself\s*\(", "from("),
        (r"\bhost\s*\(", "from("),
        (r"\bdevice\s*\(", "to("),
        (r"\buse_device\s*\(", "use_device_addr("),
    ):
        body = re.sub(pat, repl, body, flags=re.IGNORECASE)
    body = re.sub(r"\bif_present\b", " ", body, flags=re.IGNORECASE)
    body = _strip_paren_clause(body, "present")
    body = _strip_clause(body, "async")
    body = _strip_clause(body, "deviceptr")
    body = _strip_paren_clause(body, "vector_length")
    body = _strip_paren_clause(body, "num_gangs")
    body = _strip_paren_clause(body, "num_workers")
    body = _strip_paren_clause(body, "tile")
    body = _strip_paren_clause(body, "wait")
    body = _strip_paren_clause(body, "cache")
    body = _strip_clause(body, "gang")
    body = _strip_clause(body, "vector")
    body = _strip_clause(body, "worker")
    body = re.sub(r"\b(independent|seq)\b", " ", body, flags=re.IGNORECASE)
    body = _rewrite_private_clause(body)
    body = re.sub(r"\s+", " ", body).strip()
    return body


def _rewrite_private_clause(body: str) -> str:
    """OpenMP private() cannot take array sections; keep the name only."""
    out: list[str] = []
    i = 0
    while True:
        m = re.search(r"\bprivate\s*\(", body[i:], flags=re.IGNORECASE)
        if not m:
            out.append(body[i:])
            break
        start = i + m.start()
        out.append(body[i:start])
        p0 = i + m.end() - 1
        depth = 0
        k = p0
        while k < len(body):
            if body[k] == "(":
                depth += 1
            elif body[k] == ")":
                depth -= 1
                if depth == 0:
                    break
            k += 1
        inner = body[p0 + 1 : k] if (k < len(body) and body[k] == ")") else body[p0 + 1 :]
        inner = re.sub(r"\b([A-Za-z_]\w*)\s*\([^()]*\)", r"\1", inner)
        inner = _dedup_csv(inner)
        if k < len(body) and body[k] == ")":
            out.append("private(" + inner + ")")
            i = k + 1
        else:
            # Unclosed private( continued across `#ifdef` / `!$acc&`.
            out.append("private(" + inner)
            i = len(body)
    return "".join(out)


def _split_top_commas(s: str) -> list[str]:
    items: list[str] = []
    buf: list[str] = []
    depth = 0
    for ch in s:
        if ch == "(":
            depth += 1
            buf.append(ch)
        elif ch == ")":
            depth -= 1
            buf.append(ch)
        elif ch == "," and depth == 0:
            items.append("".join(buf).strip())
            buf = []
        else:
            buf.append(ch)
    tail = "".join(buf).strip()
    if tail:
        items.append(tail)
    return items


def _dedup_csv(inner: str) -> str:
    seen: set[str] = set()
    kept: list[str] = []
    for part in inner.split(","):
        item = part.strip()
        if not item:
            continue
        name_m = re.match(r"([A-Za-z_]\w*)", item)
        key = name_m.group(1).lower() if name_m else item.lower()
        if key in seen:
            continue
        seen.add(key)
        kept.append(item)
    return ", ".join(kept)


def _wait_arg(body: str) -> str:
    m = re.search(r"wait\s*\((.+)\)", body, re.IGNORECASE)
    if m:
        return m.group(1).strip()
    return "1"


def _extract_paren_inners(text: str, names: tuple[str, ...]) -> tuple[list[str], str]:
    """Pull `name(a, b(1:n))` lists out of a clause string."""
    items: list[str] = []
    body = text
    for name in names:
        while True:
            m = re.search(rf"\b{name}\s*\(", body, re.IGNORECASE)
            if not m:
                break
            p0 = m.end() - 1
            depth = 0
            k = p0
            while k < len(body):
                if body[k] == "(":
                    depth += 1
                elif body[k] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                k += 1
            inner = body[p0 + 1 : k] if k < len(body) else body[p0 + 1 :]
            for part in _split_top_commas(inner):
                part = part.strip()
                if part:
                    items.append(part)
            body = (body[: m.start()] + body[k + 1 :]).strip() if k < len(body) else body[: m.start()]
    body = re.sub(r"\s+", " ", body).strip()
    return items, body


_WRITABLE_ENTER_COPYINS = frozenset({"has_error"})


def translate_directive(line: str) -> str:
    m = SENT.match(line)
    if not m:
        return line
    indent, _amp, rest = m.group(1), m.group(2), m.group(3).strip()
    is_cont = _amp == "&" or rest.startswith("&")
    if rest.startswith("&"):
        rest = rest[1:].strip()
    low = rest.lower()

    def emit(omp: str) -> str:
        # keep free-form line length reasonable
        if is_cont:
            text = f"{indent}!$omp& {omp}".rstrip()
        else:
            text = f"{indent}!$omp {omp}".rstrip()
        if len(text) <= 132:
            return text
        # break after a clause comma if possible
        chunks = []
        s = text
        while len(s) > 132:
            br = s.rfind(",", 0, 120)
            if br < 40:
                br = 120
            chunks.append(s[: br + 1] + " &")
            s = indent + "!$omp& " + s[br + 1 :].lstrip()
        chunks.append(s)
        return "\n".join(chunks)

    if is_cont:
        body = _rewrite_clauses(rest)
        if not body:
            return f"{indent}! acc2omp: stripped continuation"
        return emit(body)

    if low.startswith("cache"):
        return f"{indent}! acc2omp: dropped acc cache"

    if low.startswith("end serial"):
        return emit("end target")
    if low.startswith("end host_data"):
        return emit("end target data")
    if low.startswith("end loop"):
        # Combined `parallel loop` already closed the workshare; leftover
        # `end loop` is a nested-loop closer (serial on the parent thread).
        return f"{indent}! acc2omp: end nested loop"
    if low.startswith("end parallel loop") or low.startswith("end parallel"):
        return emit("end target teams distribute parallel do")
    if low.startswith("end kernels"):
        return emit("end target")
    if low.startswith("end data"):
        return emit("end target data")
    if low.startswith("wait"):
        arg = _wait_arg(rest)
        return f"{indent}call geps_acc_wait({arg})"

    if low.startswith("host_data"):
        body = _rewrite_clauses(rest[len("host_data") :])
        return emit(f"target data {body}".strip())
    if low.startswith("enter data"):
        raw = rest[len("enter data") :]
        # enter data may only map(to/alloc). copyin→tofrom is illegal here.
        # Default: map(to:) (OpenACC enter copyin). Writable device scalars
        # that kernels store into get alloc+update so the mapping is not RO.
        copyins, raw = _extract_paren_inners(
            raw, ("present_or_copyin", "present_or_copy", "copyin", "copy")
        )
        writable = [p for p in copyins if p.split("(", 1)[0].lower() in _WRITABLE_ENTER_COPYINS]
        readonly = [p for p in copyins if p.split("(", 1)[0].lower() not in _WRITABLE_ENTER_COPYINS]
        if readonly:
            raw = f"{raw} copyin({', '.join(readonly)})".strip()
        if writable:
            raw = f"{raw} create({', '.join(writable)})".strip()
        body = _rewrite_clauses(raw, enter_data=True)
        enter = emit(f"target enter data {body}".strip())
        if writable:
            return enter + "\n" + emit(f"target update to({', '.join(writable)})")
        return enter
    if low.startswith("exit data"):
        body = _rewrite_clauses(rest[len("exit data") :])
        return emit(f"target exit data {body}".strip())
    if low.startswith("update"):
        body = _rewrite_clauses(rest[len("update") :])
        return emit(f"target update {body}".strip())
    if low.startswith("atomic"):
        body = rest[len("atomic") :].strip()
        return emit(f"atomic {body}".strip())
    if low.startswith("routine"):
        body = re.sub(r"\bseq\b", " ", rest[len("routine") :], flags=re.IGNORECASE).strip()
        # `routine(name)` is a callee prototype; OpenMP wants declare target on the definition.
        if body.startswith("("):
            return f"{indent}! acc2omp: routine{body} (declare target lives on definition)"
        return emit("declare target")
    if low.startswith("kernels"):
        body = _rewrite_clauses(rest[len("kernels") :])
        extra = f" {body}" if body else ""
        return emit(f"target{extra}".strip())
    if low.startswith("serial"):
        body = _rewrite_clauses(rest[len("serial") :])
        extra = f" {body}" if body else ""
        return emit(f"target{extra}".strip())
    if low.startswith("data"):
        body = _rewrite_clauses(rest[len("data") :])
        return emit(f"target data {body}".strip())
    if low.startswith("parallel loop"):
        body = _rewrite_clauses(rest[len("parallel loop") :])
        extra = f" {body}" if body else ""
        return emit(f"target teams distribute parallel do{extra}")
    if low.startswith("parallel"):
        body = _rewrite_clauses(rest[len("parallel") :])
        extra = f" {body}" if body else ""
        return emit(f"target teams{extra}")
    if low.startswith("loop"):
        # Remaining `loop` after `_merge_parallel_loop` is nested inside
        # `parallel loop` / `kernels`. Emitting `parallel do` there starts a
        # nested OpenMP region on the GPU (wrong kkh/kkl, Recursive I/O).
        # Serial execution on the parent thread matches OpenACC gang+seq and
        # is valid for inner `loop vector` that contain calls / while.
        return f"{indent}! acc2omp: nested loop (serial on parent thread)"

    body = _rewrite_clauses(rest)
    return emit(body)


def sanitize_cuda_fortran(line: str) -> str:
    if line.lstrip().startswith("!"):
        return line
    raw = line
    # CUDA Fortran pointer-kind inquiry is not a valid kind expr on flang.
    line = re.sub(r"int_ptr_kind\s*\(\s*\)", "8", line, flags=re.IGNORECASE)
    # Drop CUDA Fortran data attributes (HIP shims use mapped host arrays / pointers).
    if re.search(r",\s*device\b", line, re.IGNORECASE) and re.search(
        r"\ballocatable\b", line, re.IGNORECASE
    ):
        line = re.sub(r"\ballocatable\b", "pointer", line, flags=re.IGNORECASE)
    line = re.sub(r",\s*device\b", "", line, flags=re.IGNORECASE)
    line = re.sub(r",\s*pinned\b", "", line, flags=re.IGNORECASE)
    line = re.sub(r",\s*managed\b", "", line, flags=re.IGNORECASE)
    line = re.sub(r"\s+,", ",", line)
    line = re.sub(r",\s+::", " ::", line)
    # NV Fortran allows LOGICAL .eq. .true.; amdflang requires .eqv. / .neqv.
    line = re.sub(r"\.eq\.\s*\.true\.", ".eqv. .true.", line, flags=re.IGNORECASE)
    line = re.sub(r"\.eq\.\s*\.false\.", ".eqv. .false.", line, flags=re.IGNORECASE)
    line = re.sub(r"\.ne\.\s*\.true\.", ".neqv. .true.", line, flags=re.IGNORECASE)
    line = re.sub(r"\.ne\.\s*\.false\.", ".neqv. .false.", line, flags=re.IGNORECASE)
    line = re.sub(r"==\s*\.true\.", ".eqv. .true.", line, flags=re.IGNORECASE)
    line = re.sub(r"==\s*\.false\.", ".eqv. .false.", line, flags=re.IGNORECASE)
    line = re.sub(r"/=\s*\.true\.", ".neqv. .true.", line, flags=re.IGNORECASE)
    line = re.sub(r"/=\s*\.false\.", ".neqv. .false.", line, flags=re.IGNORECASE)
    # Named INTERFACE wrapping two bind(C) procs with the same TKR is a generic;
    # amdflang rejects it. Anonymous INTERFACE still gives explicit interfaces.
    line = re.sub(r"\binterface\s+c_interface\b", "interface", line, flags=re.IGNORECASE)
    line = re.sub(r"\bend\s+interface\s+c_interface\b", "end interface", line, flags=re.IGNORECASE)
    for spec, gen in (
        ("dlog", "log"),
        ("dsqrt", "sqrt"),
        ("dexp", "exp"),
        ("dsin", "sin"),
        ("dcos", "cos"),
        ("dtan", "tan"),
        ("dabs", "abs"),
        ("dmax1", "max"),
        ("dmin1", "min"),
        ("dsign", "sign"),
        ("dmod", "mod"),
        ("alog10", "log10"),
        ("alog", "log"),
        ("dfloat", "real"),
        ("dble", "real"),
    ):
        line = re.sub(rf"\b{spec}\s*\(", f"{gen}(", line, flags=re.IGNORECASE)
    # F77-style: saturation function typed as a local REAL. amdflang then
    # sees a public func.func symbol. Keep the result name inside the function.
    if not re.search(r"\bfunction\s+fpvs", line, re.IGNORECASE):
        line = re.sub(r",\s*fpvs_gpu\b", "", line, flags=re.IGNORECASE)
        line = re.sub(r"\bfpvs_gpu\s*,\s*", "", line, flags=re.IGNORECASE)
        line = re.sub(r",\s*fpvs\b", "", line, flags=re.IGNORECASE)
        line = re.sub(r"\bfpvs\s*,\s*", "", line, flags=re.IGNORECASE)
    # -fdefault-real-8: D-exponents are REAL(16); use E so they stay default real.
    line = re.sub(r"\b(\d*\.\d+|\d+\.\d*|\d+)[dD]([-+]?\d+)", r"\1e\2", line)
    return line if line != raw else raw


def _acc_rest(line: str) -> str:
    m = SENT.match(line)
    if not m:
        return ""
    rest = m.group(3).strip()
    if rest.startswith("&"):
        rest = rest[1:].strip()
    return rest


def _merge_parallel_loop(joined: list[tuple[str, bool]]) -> list[tuple[str, bool]]:
    """Fold `!$acc parallel` + nested `!$acc loop` into `parallel loop`."""
    out: list[tuple[str, bool]] = []
    i = 0
    n = len(joined)
    while i < n:
        line, is_acc = joined[i]
        if is_acc:
            rest = _acc_rest(line)
            low = rest.lower()
            if (
                low.startswith("parallel")
                and not low.startswith("parallel loop")
                and not low.startswith("end")
            ):
                j = i + 1
                skipped_code = False
                while j < n:
                    lj, ja = joined[j]
                    if ja:
                        break
                    s = lj.strip()
                    if s and not s.startswith("!"):
                        skipped_code = True
                        break
                    j += 1
                if (
                    not skipped_code
                    and j < n
                    and joined[j][1]
                    and _acc_rest(joined[j][0]).lower().startswith("loop")
                ):
                    indent = SENT.match(line).group(1)
                    pbody = rest[len("parallel") :].strip()
                    lbody = _acc_rest(joined[j][0])[len("loop") :].strip()
                    merged = f"{indent}!$acc parallel loop {pbody} {lbody}".rstrip()
                    out.append((merged, True))
                    i = j + 1
                    continue
        out.append((line, is_acc))
        i += 1
    return out



_LOOP_DO = re.compile(r"^\s*(?:\w+\s*:\s*)?do\s+([A-Za-z_]\w*)\s*=", re.IGNORECASE)
# any DO statement: `do k = ...`, `do 100 k = ...`, `do while (...)`, bare `do`
# (layer 18: sflx_gpu's `do while` was not counted, so the nest tracker
# thought the 9000-line kernel ended at the first `end do ! do_while_loop`)
_DO_ANY = re.compile(r"^\s*(?:\w+\s*:\s*)?do\s*(?:(\d+)\s*,?\s*)?(?:while\b|[A-Za-z_]\w*\s*=|$)", re.IGNORECASE)
_LABEL_CONT = re.compile(r"^\s*(\d+)\s+continue\b", re.IGNORECASE)


def _do_start(code: str):
    """(True, label-or-None) if `code` opens a DO loop, else (False, None)."""
    m = _DO_ANY.match(code)
    if not m:
        return False, None
    return True, (int(m.group(1)) if m.group(1) else None)
_SCOPE_END = re.compile(r"^\s*(end\s+)?(subroutine|function|program|module|contains)\b", re.IGNORECASE)
_CPP_OPEN = re.compile(r"^\s*#\s*(if|ifdef|ifndef)\b")
_CPP_CLOSE = re.compile(r"^\s*#\s*endif\b")


def _clause_names(text: str) -> set:
    """Identifiers inside any clause, nesting included.

    A flat `\\(([^()]*)\\)` scan misses `private(i, dqmono(:3*lons))` - the outer
    group never matches, the names in it look unused, and adding them again
    trips "appears in more than one data-sharing clause".
    """
    names = set()
    i, n = 0, len(text)
    while i < n:
        if text[i] != "(":
            i += 1
            continue
        depth, j = 1, i + 1
        while j < n and depth:
            if text[j] == "(":
                depth += 1
            elif text[j] == ")":
                depth -= 1
            j += 1
        for tok in re.findall(r"[A-Za-z_]\w*", text[i + 1 : j - 1]):
            names.add(tok.lower())
        i = j
    return names


def _privatise_nested_loops(joined):
    """Give the enclosing parallel directive the private() of inner `!$acc loop`s.

    An inner `!$acc loop` that `_merge_parallel_loop` could not fold in (real
    code sits between it and the enclosing directive) is later flattened to a
    serial DO on the parent thread by the `loop` branch of translate_directive().
    Its clauses are dropped there, so its private() names - and the flattened
    DO's own iteration variables - stay SHARED across the threads of the
    enclosing `parallel loop` and every thread races on them. That is what made
    nnmi_gpu's `bal` garbage and different on every run (handoff layer 9).

    Runs on the joined list, so each directive is one logical string.

    Three things must be respected, each learned the hard way:
      * dedup against every name already in the directive (`_clause_names`),
        or flang rejects a variable in two data-sharing clauses;
      * a `parallel loop` owns exactly the loop nest that follows it, so close
        the region when that nest ends - otherwise later loops, even in another
        subroutine, get their privates appended to it;
      * never hoist across a cpp conditional. In
        module_mp_gsfcgce_3ice_nuwrf_gpu.f90 both `mr2mc`'s declaration and the
        inner loop live under `#ifdef Readaeroclx` (OFF by default), while the
        owner directive sits outside it - hoisting lifts the name out of the
        conditional and it no longer exists.
    """
    out = [list(t) for t in joined]
    owner = -1
    labels9: list = []
    depth = 0            # DO nesting under the owner
    armed = False        # its loop nest has started
    cpp = 0              # current cpp conditional depth
    owner_cpp = 0        # cpp depth where the owner was opened
    for i, (line, is_acc) in enumerate(joined):
        if not is_acc:
            code = line.strip()
            if _CPP_OPEN.match(code):
                cpp += 1
                continue
            if _CPP_CLOSE.match(code):
                cpp -= 1
                continue
            if owner < 0:
                continue
            low_c = code.lower()
            is_do, lab = _do_start(code)
            lc = _LABEL_CONT.match(code)
            if is_do:
                depth += 1
                labels9.append(lab)
                armed = True
            elif re.match(r"^end\s*do\b", low_c) or (lc and int(lc.group(1)) in labels9):
                if lc:
                    n_ = int(lc.group(1))
                    while labels9 and labels9[-1] == n_:
                        labels9.pop(); depth -= 1
                else:
                    if labels9: labels9.pop()
                    depth -= 1
                if armed and depth <= 0:
                    owner = -1
            elif _SCOPE_END.match(code):
                owner = -1
            continue

        rest = _acc_rest(line)
        low = rest.lower()
        if low.startswith("end"):
            if "parallel" in low or "kernels" in low:
                owner = -1
            continue
        if low.startswith("parallel"):
            # `parallel` / `parallel loop` spawn threads; `kernels` becomes a
            # single-threaded `!$omp target` and needs no privatisation.
            owner = -1 if line.rstrip().endswith("&") else i
            owner_cpp = cpp
            depth, armed = 0, False; labels9 = []
            continue
        if low.startswith(("kernels", "data", "serial", "routine", "declare")):
            owner = -1
            continue
        if not low.startswith("loop") or owner < 0:
            continue
        if cpp != owner_cpp:
            # inner loop is under a cpp conditional the owner is not in
            continue

        names = []
        mp = re.search(r"\bprivate\s*\(", rest, re.IGNORECASE)
        if mp:
            st = mp.end()
            d, k = 1, st
            while k < len(rest) and d:
                if rest[k] == "(":
                    d += 1
                elif rest[k] == ")":
                    d -= 1
                k += 1
            for item in re.split(r",(?![^()]*\))", rest[st : k - 1]):
                nm = re.match(r"\s*([A-Za-z_]\w*)", item)
                if nm:
                    names.append(nm.group(1))
        mc = re.search(r"\bcollapse\s*\(\s*(\d+)\s*\)", rest, re.IGNORECASE)
        want = int(mc.group(1)) if mc else 1
        seen = 0
        for j in range(i + 1, min(i + 40, len(joined))):
            if joined[j][1]:
                break
            dm = _LOOP_DO.match(joined[j][0])
            if dm:
                names.append(dm.group(1))
                seen += 1
                if seen >= want:
                    break
        if not names:
            continue
        target = out[owner][0]
        add = [n for n in names if n.lower() not in _clause_names(target)]
        if add:
            out[owner][0] = target.rstrip() + " private(" + ", ".join(add) + ")"
    return [tuple(t) for t in out]

_DECL_RE = re.compile(
    r"^\s*(real|integer|logical|double\s+precision|complex|character|type\s*\(\s*\w+\s*\))"
    r"(\s*\(\s*kind\s*=[^)]*\)|\s*\*\s*\d+|\s*\([^)]*\))?"
    r"(?P<attrs>[^:!]*?)(::)?(?P<ents>.*)$",
    re.IGNORECASE,
)
_ASSIGN_RE = re.compile(r"^\s*(?:if\s*\(.*\)\s*)?([A-Za-z_]\w*)\s*=(?!=)", re.IGNORECASE)


def _split_top(text: str) -> list[str]:
    out, depth, cur = [], 0, ""
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur); cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return out


_DATA_RE = re.compile(r"^\s*data\s+(.*)$", re.IGNORECASE)


def _data_names(code: str) -> list:
    """Names initialised by a DATA statement: `data a/2.0/, b, c/1., 2./`."""
    m = _DATA_RE.match(code)
    if not m:
        return []
    names = []
    parts = m.group(1).split("/")
    for k in range(0, len(parts) - 1, 2):      # even segments are name lists
        for tok in re.findall(r"[A-Za-z_]\w*", parts[k]):
            names.append(tok.lower())
    return names


def _private_names(text: str) -> set:
    names = set()
    for m in re.finditer(r"(?<![A-Za-z_])private\s*\(", text, re.IGNORECASE):
        j, depth, n = m.end(), 1, len(text)
        while j < n and depth:
            if text[j] == "(": depth += 1
            elif text[j] == ")": depth -= 1
            j += 1
        for e in _split_top(text[m.end(): j - 1]):
            mm = re.match(r"\s*([A-Za-z_]\w*)\s*(\(|$)", e.strip() + " ")
            if mm and mm.group(2) != "(":
                names.add(mm.group(1).lower())
    return names


def _strip_from_private(text: str, names: set) -> str:
    """Remove `names` from every private(...) clause (nesting-aware); drop a
    clause that becomes empty."""
    out, i, n = [], 0, len(text)
    while i < n:
        m = re.match(r"(?<![A-Za-z_])private\s*\(", text[i:], re.IGNORECASE)
        if not (m and text[i:i+1].lower() == "p"):
            out.append(text[i]); i += 1
            continue
        j, depth = i + m.end(), 1
        while j < n and depth:
            if text[j] == "(": depth += 1
            elif text[j] == ")": depth -= 1
            j += 1
        inner = text[i + m.end(): j - 1]
        keep = [e for e in _split_top(inner)
                if re.match(r"\s*([A-Za-z_]\w*)", e) and re.match(r"\s*([A-Za-z_]\w*)", e).group(1).lower() not in names]
        if keep:
            out.append("private(" + ", ".join(k.strip() for k in keep) + ")")
        else:
            # eat a trailing space so the directive stays tidy
            while j < n and text[j] == " ":
                j += 1
        i = j
    return "".join(out)


def _privatise_assigned_scalars(joined):
    """OpenACC makes every scalar in a `parallel` region firstprivate per gang
    unless a data clause says otherwise; OpenMP shares it across the target's
    teams. A scalar ASSIGNED inside the loop nest therefore races between
    teams after translation (handoff layer 16, 2026-09-16): hdiffu_gpu sets
    `hfiltd`/`powdd`/`c1`/`c2` per level k inside `parallel loop gang
    private(KL, kfac, ...)` and only some of them were listed - vormid/divmid
    came out 1e-3 off and run-dependent. Fix: add `firstprivate(...)` for every
    declared SCALAR of the program unit that is assigned inside the nest
    (arrays stay shared, as in OpenACC). Declarations are read per program
    unit from single-line type statements; a scalar declared on a continued
    declaration is simply not found, which leaves today's behaviour.
    """
    out = [list(t) for t in joined]
    scalars: dict = {}        # name -> cpp depth of its declaration
    datainit: set = set()     # DATA-initialised names of the program unit
    owner = -1
    owner_cpp = 0
    cpp = 0
    depth = 0
    armed = False
    assigned: list = []       # (name, cpp depth of the assignment)
    touched: set = set()      # names that appear in a `call` inside the nest
    pending_decl = ""         # continued declaration being joined for scanning
    labels16: list = []       # open DO loops' labels (None for end-do loops)
    loopvars: set = set()     # DO variables inside the nest

    def flush():
        nonlocal owner, assigned, touched
        if owner >= 0:
            target = out[owner][0]
            have = _clause_names(target)
            add = []
            for n, adepth in assigned:
                ln_ = n.lower()
                # never hoist a name out of a cpp conditional the owner is not
                # in (layer-9 trap: mr2mc under #ifdef Readaeroclx)
                if ln_ not in scalars or adepth != owner_cpp:
                    continue
                if scalars[ln_] not in (0, owner_cpp):
                    continue
                if ln_ not in have and ln_ not in [a.lower() for a in add]:
                    add.append(n)
            # Layer 18 (2026-09-17): a scalar listed in `private(...)` that is
            # NEVER assigned inside the nest is read uninitialised under
            # OpenMP; nvfortran happens to seed it with the outer value.
            # sflx_gpu lists the DATA-initialised `snoexp` (=2.0) as private
            # and uses it in sncovr**snoexp -> skin temperature Inf on ~30k
            # snow columns -> NaN through PBL -> microphysics spins forever.
            # Move such scalars from private to firstprivate.
            was_assigned = {n.lower() for n, _ in assigned}
            demote = set()
            for nm in _private_names(target):
                if nm in scalars and nm not in was_assigned and nm not in touched and nm not in loopvars:
                    demote.add(nm)
            if demote:
                target = _strip_from_private(target, demote)
                add += sorted(demote)
            if add:
                target = target.rstrip() + " firstprivate(" + ", ".join(add) + ")"
            out[owner][0] = target
        owner = -1
        assigned = []
        touched = set()
        loopvars.clear()

    for i, (line, is_acc) in enumerate(joined):
        code = line.strip()
        if not is_acc:
            if _CPP_OPEN.match(code):
                cpp += 1
                continue
            if _CPP_CLOSE.match(code):
                cpp -= 1
                continue
            if code.startswith("#") or code.startswith("!"):
                continue
            low_c = code.lower()
            if re.match(r"^(subroutine|function|program|module)\b", low_c) or re.match(r"^[a-z0-9_() ,=*]*\bfunction\s+\w+", low_c):
                flush(); scalars = {}; datainit = set()
                continue
            if _SCOPE_END.match(code):
                flush()
                continue
            if owner < 0:
                datainit.update(_data_names(code))
                # layer 18: continued declarations (`real :: a, b, &`) are
                # joined here for scanning only; the output keeps the lines.
                if pending_decl:
                    pending_decl += " " + code.lstrip("&").strip().rstrip("&").strip()
                    if code.rstrip().endswith("&"):
                        continue
                    code = pending_decl; pending_decl = ""
                    low_c = code.lower()
                elif _DECL_RE.match(code) and code.rstrip().endswith("&") and \
                        not re.search(r"\b(function|subroutine)\b", low_c):
                    # layer 18b: also the `::`-less form `real(kind=x) a, b, &`
                    # (samfdeepcnv_kh_gpu declares val1/val2 that way; they
                    # stayed private and were read uninitialised)
                    pending_decl = code.rstrip().rstrip("&").strip()
                    continue
                dm = _DECL_RE.match(code)
                if dm and "::" in code or (dm and not re.search(r"\b(function|subroutine)\b", low_c)):
                    if "::" in code:
                        # layer 18: `_DECL_RE`'s lazy attrs group leaves
                        # `, dimension(n) :: a, b` inside `ents`, so arrays
                        # looked like scalars. Split at `::` explicitly.
                        attrs, _, ents = code.partition("::")
                        attrs = attrs.lower()
                    else:
                        attrs = (dm.group("attrs") or "").lower()
                        ents = dm.group("ents")
                    if "dimension" in attrs or "allocatable" in attrs or "pointer" in attrs or "parameter" in attrs:
                        continue
                    for ent in _split_top(ents):
                        m = re.match(r"\s*([A-Za-z_]\w*)\s*(\(|=|$|\*)", ent.strip() + " ")
                        if m and m.group(2) != "(":
                            scalars.setdefault(m.group(1).lower(), cpp)
                continue
            is_do, lab = _do_start(code)
            if is_do:
                depth += 1; armed = True
                labels16.append(lab)
                lm = re.match(r"^\s*(?:\w+\s*:\s*)?do\s+(?:\d+\s*,?\s*)?([A-Za-z_]\w*)\s*=", code, re.IGNORECASE)
                if lm:
                    loopvars.add(lm.group(1).lower())   # predetermined private: never demote, never add
                continue
            lc = _LABEL_CONT.match(code)
            if re.match(r"^end\s*do\b", low_c) or (lc and int(lc.group(1)) in labels16):
                if lc:
                    n_ = int(lc.group(1))
                    while labels16 and labels16[-1] == n_:
                        labels16.pop(); depth -= 1
                else:
                    if labels16: labels16.pop()
                    depth -= 1
                if armed and depth <= 0:
                    flush()
                continue
            am = _ASSIGN_RE.match(code)
            if am:
                assigned.append((am.group(1), cpp))
            if owner >= 0 and re.match(r"^\s*(?:if\s*\(.*\)\s*)?call\b", low_c):
                touched.update(t.lower() for t in re.findall(r"[A-Za-z_]\w*", code))
            continue
        rest = _acc_rest(line)
        low = rest.lower()
        if low.startswith("end"):
            if "parallel" in low or "kernels" in low:
                flush()
            continue
        if low.startswith("parallel"):
            flush()
            owner = -1 if line.rstrip().endswith("&") else i
            owner_cpp = cpp
            depth, armed = 0, False; labels16 = []
            continue
        if low.startswith(("kernels", "data", "serial", "routine", "declare")):
            flush()
    flush()
    return [tuple(t) for t in out]


def _join_decl_continuations(lines: list[str]) -> list[str]:
    """Join `real :: a, &` declaration lists, but never across `#if` / `#else`."""
    out: list[str] = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        nxt_ok = (
            i + 1 < n
            and not lines[i + 1].lstrip().startswith("#")
            and not lines[i + 1].lstrip().startswith("!")
        )
        if (
            line.rstrip().endswith("&")
            and "::" in line
            and not line.lstrip().startswith("!")
            and nxt_ok
        ):
            buf = line.rstrip()[:-1].rstrip()
            i += 1
            while i < n:
                nxtline = lines[i]
                if nxtline.lstrip().startswith("#"):
                    break
                if nxtline.lstrip().startswith("!"):
                    i += 1
                    continue
                nxt = nxtline.strip()
                if nxt.startswith("&"):
                    nxt = nxt[1:].strip()
                ends = nxtline.rstrip().endswith("&")
                if ends:
                    if nxt.endswith("&"):
                        nxt = nxt[:-1].rstrip()
                    buf = buf + " " + nxt
                    i += 1
                else:
                    buf = buf + " " + nxt
                    i += 1
                    break
            out.append(buf)
            continue
        out.append(line)
        i += 1
    return out


def drop_duplicate_decls(lines: list[str]) -> list[str]:
    """amdflang rejects a name declared twice in the same subprogram (nvfortran does not).

    Preprocessor branches (`#ifdef` / `#else`) are treated as disjoint scopes so
    `RTYPE=4` / `RTYPE=8` alternatives are not collapsed.
    """
    decl = re.compile(
        r"^(\s*)((?:integer|real|logical|character|complex|double\s+precision)\b[^:\n]*::\s*)(.+)$",
        re.IGNORECASE,
    )
    declared: set[str] = set()
    bounded: set[str] = set()
    if_stack: list[tuple[set[str], set[str]]] = []
    proc_name = ""
    out: list[str] = []
    lines = _join_decl_continuations(lines)
    for line in lines:
        stripped = line.lstrip()
        low = stripped.lower()
        if stripped.startswith("#"):
            cmd = stripped[1:].split()[0].lower() if stripped[1:].split() else ""
            if cmd in ("if", "ifdef", "ifndef"):
                if_stack.append((set(declared), set(bounded)))
            elif cmd in ("else", "elif", "elseif"):
                if if_stack:
                    declared, bounded = set(if_stack[-1][0]), set(if_stack[-1][1])
            elif cmd == "endif":
                if if_stack:
                    if_stack.pop()
            out.append(line)
            continue
        if re.match(r"(end\s+)?(subroutine|function|module|program)\b", low):
            declared = set()
            bounded = set()
            if_stack = []
            if low.startswith("end"):
                proc_name = ""
            else:
                pm = re.search(r"\b(?:subroutine|function)\s+(\w+)", low)
                proc_name = pm.group(1) if pm else ""
            out.append(line)
            continue
        if stripped.startswith("!"):
            out.append(line)
            continue
        m = decl.match(line)
        if m and "&" not in m.group(3) and m.group(3).count("(") == m.group(3).count(")"):
            indent, typed, rest = m.group(1), m.group(2), m.group(3)
        else:
            old = re.match(
                r"^(\s*)((?:integer|real|logical|character|complex)\b(?:\s*\([^)]*\))?)\s+"
                r"(?!(?:function|subroutine|kind\b|\*))(.+)$",
                line,
                re.IGNORECASE,
            )
            if (
                old
                and "::" not in line
                and "&" not in old.group(3)
                and old.group(3).count("(") == old.group(3).count(")")
            ):
                indent, typed, rest = old.group(1), old.group(2) + " ", old.group(3)
            else:
                cm = re.match(
                    r"^(\s*)common\s*(/[^/]+/\s*)?(.*)$",
                    line,
                    re.IGNORECASE,
                )
                if cm and "::" not in line:
                    indent, blk, rest = cm.group(1), cm.group(2) or "", cm.group(3)
                    rest = rest.split("!")[0]
                    kept: list[str] = []
                    for item in _split_top_commas(rest):
                        item = item.strip()
                        if not item:
                            continue
                        nm = re.match(r"([A-Za-z_]\w*)(\s*\(.*\))?$", item)
                        if nm:
                            n = nm.group(1).lower()
                            if n in bounded and nm.group(2):
                                kept.append(nm.group(1))
                            else:
                                kept.append(item)
                            declared.add(n)
                            if nm.group(2):
                                bounded.add(n)
                        else:
                            kept.append(item)
                    out.append(indent + "common " + blk + ", ".join(kept))
                    continue
                out.append(line)
                continue
        kept: list[str] = []
        for item in _split_top_commas(rest):
            if not item:
                continue
            nm = re.match(r"([A-Za-z_]\w*)", item)
            if not nm:
                kept.append(item)
                continue
            n = nm.group(1).lower()
            if n in ("fpvs_gpu", "fpvs") and proc_name not in ("fpvs_gpu", "fpvs"):
                continue
            if n in declared:
                continue
            declared.add(n)
            if "(" in item:
                bounded.add(n)
            kept.append(item)
        if not kept:
            out.append(indent + "! acc2omp: duplicate decl removed")
        else:
            out.append(indent + typed + ", ".join(kept))
    return out


# Symbols that `use <mod>` without ONLY would import. Used to hide collisions
# with dummy arguments (nvfortran allows the shadowing; amdflang does not).
MODULE_PUBLIC: dict[str, set[str]] = {
    "index": {
        "mlistnum", "mlist", "nlist", "ilist", "jlistnum", "jlist1", "jlist2",
        "jlistnum_sl", "jlist2_2d", "jlist1_sl", "lreduce", "nxdef", "mtrundef",
        "nsizex", "nsizey", "mrow", "ncol", "row_comm", "col_comm", "row_rank",
        "col_rank", "nccl_row_comm", "nccl_col_comm", "tcolt_jlist", "poly_mlist",
        "lstart", "lend", "llen", "nxp", "levp", "jlen", "nxf", "levf", "myf",
        "llistnum", "lstart_ncld", "lend_ncld", "jtf", "jtp", "jtstart", "jtend",
        "jtlen", "nxptot", "llist", "llist_ncld", "jtlen_all", "nxjp", "nxjstart",
        "nxjend", "nxjlen", "nxdef_2d", "nxjp_acc", "nxjstart_all", "nxjend_all",
        "nxjlen_all", "map2to1",
    },
    "param": {
        "nco", "nx", "my", "lev", "jtrun", "mlmax", "octahedral", "ncld",
        "alpha", "mwhd", "nout", "io_quilting", "npe", "jtmax", "my_max",
        "nx_max", "npex", "npey",
    },
    "rank": {
        "nsize", "myrank", "nsize_all", "nsize_gfs", "nsize_io", "myrank_all",
        "myrank_gfs", "myrank_io", "root_gfs", "root_io", "mpi_comm_gfs_all",
        "root_rsm", "itag", "mpi_comm_gfs", "mpi_comm_io", "ntag", "ngfs", "nio",
        "mpi_comm_atm", "nccl_id", "nccl_comm_gfs",
    },
}


def _dummy_names_from(lines: list[str], start: int) -> tuple[set[str], int]:
    """Collect dummy argument names from a subroutine/function statement."""
    buf: list[str] = []
    i = start
    n = len(lines)
    depth = 0
    started = False
    while i < n:
        raw = lines[i]
        stripped = raw.lstrip()
        if stripped.startswith("#") or stripped.startswith("!"):
            i += 1
            continue
        for ch in raw.split("!")[0]:
            if ch == "(":
                depth += 1
                started = True
            elif ch == ")":
                depth -= 1
        buf.append(raw.split("!")[0])
        i += 1
        if started and depth <= 0:
            break
        if not started and not raw.rstrip().endswith("&"):
            break
    text = " ".join(buf)
    names: set[str] = set()
    m = re.search(r"\((.*)\)", text, re.DOTALL)
    if not m:
        return names, start + 1
    for item in _split_top_commas(m.group(1)):
        nm = re.match(r"([A-Za-z_]\w*)", item.strip())
        if nm:
            names.add(nm.group(1).lower())
    return names, i


def hide_use_collisions(lines: list[str]) -> list[str]:
    dummies: set[str] = set()
    out: list[str] = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        stripped = line.lstrip()
        low = stripped.lower()
        if re.match(r"(end\s+)?(subroutine|function|module|program)\b", low):
            if low.startswith("end"):
                dummies = set()
                out.append(line)
                i += 1
                continue
            dummies, nxt = _dummy_names_from(lines, i)
            out.extend(lines[i:nxt])
            i = nxt
            continue
        um = re.match(r"^(\s*)use\s+(\w+)\s*(!.*)?$", line, re.IGNORECASE)
        if um and "only" not in low and "=>" not in line:
            mod = um.group(2).lower()
            pub = MODULE_PUBLIC.get(mod, set())
            hits = sorted(dummies & pub)
            if hits:
                indent = um.group(1)
                renames = ", ".join(f"geps_hide_{n} => {n}" for n in hits)
                out.append(f"{indent}use {um.group(2)}, {renames}")
                i += 1
                continue
        out.append(line)
        i += 1
    return out


def inject_fpvs_external(lines: list[str]) -> list[str]:
    """implicit none + fpvs_gpu() needs an EXTERNAL, not a local REAL."""
    starts: list[int] = []
    for i, ln in enumerate(lines):
        low = ln.lstrip().lower()
        if re.match(r"(end\s+)?(subroutine|function)\b", low) and not low.startswith("end"):
            starts.append(i)
    if not starts:
        return lines
    out: list[str] = list(lines[: starts[0]])
    starts.append(len(lines))
    for s, e in zip(starts, starts[1:]):
        chunk = lines[s:e]
        text = "\n".join(chunk)
        n_call = len(re.findall(r"(?<![A-Za-z_])fpvs_gpu\s*\(", text, re.IGNORECASE))
        n_def = len(re.findall(r"\bfunction\s+fpvs_gpu\s*\(", text, re.IGNORECASE))
        if n_call > n_def and not re.search(
            r"\bexternal\b[^\n]*fpvs_gpu|\bfpvs_gpu\b[^\n]*\bexternal\b",
            text,
            re.IGNORECASE,
        ):
            injected = False
            new: list[str] = []
            for ln in chunk:
                new.append(ln)
                if not injected and re.match(r"\s*implicit\s+none\b", ln, re.IGNORECASE):
                    indent = re.match(r"^(\s*)", ln).group(1)
                    new.append(f"{indent}real, external :: fpvs_gpu")
                    injected = True
            if not injected:
                indent = re.match(r"^(\s*)", chunk[0]).group(1)
                new = [chunk[0], f"{indent}real, external :: fpvs_gpu", *chunk[1:]]
            out.extend(new)
        else:
            out.extend(chunk)
    return out


GPU_TRUE = {"USE_CUDA", "USE_HIP", "USE_GPU"}


def _ifdef_known(cmd: str, rest: str) -> bool | None:
    """True/False if this cpp condition is fully known for HIP builds; else None."""
    rest = rest.strip()
    if cmd == "ifdef":
        name = rest.split()[0] if rest.split() else ""
        if name.upper() in GPU_TRUE:
            return True
        return None
    if cmd == "ifndef":
        name = rest.split()[0] if rest.split() else ""
        if name.upper() in GPU_TRUE:
            return False
        return None
    if cmd == "if":
        m = re.fullmatch(r"defined\s*\(\s*(\w+)\s*\)", rest, re.IGNORECASE)
        if m and m.group(1).upper() in GPU_TRUE:
            return True
        m = re.fullmatch(r"!\s*defined\s*\(\s*(\w+)\s*\)", rest, re.IGNORECASE)
        if m and m.group(1).upper() in GPU_TRUE:
            return False
        if re.fullmatch(r"\w+", rest) and rest.upper() in GPU_TRUE:
            return True
        return None
    return None


def resolve_gpu_ifdefs(lines: list[str]) -> list[str]:
    """HIP always compiles with -DUSE_CUDA -DUSE_HIP; fold those branches for the scanner."""
    out: list[str] = []
    # stack: None = pass-through unknown; True/False = known keep
    stack: list[bool | None] = []
    keeping = True

    def active() -> bool:
        return all(s is not False for s in stack)

    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("#"):
            parts = stripped[1:].split(None, 1)
            cmd = parts[0].lower() if parts else ""
            rest = parts[1] if len(parts) > 1 else ""
            if cmd in ("if", "ifdef", "ifndef"):
                known = _ifdef_known(cmd, rest)
                if known is None:
                    stack.append(None)
                    if active():
                        out.append(line)
                else:
                    stack.append(known)
                continue
            if cmd in ("else", "elif", "elseif"):
                if not stack:
                    out.append(line)
                    continue
                top = stack[-1]
                if top is None:
                    if active():
                        out.append(line)
                else:
                    if cmd != "else":
                        # unknown elif on a known if: keep the directive
                        out.append(line)
                        stack[-1] = None
                    else:
                        stack[-1] = not top
                continue
            if cmd == "endif":
                if not stack:
                    out.append(line)
                    continue
                top = stack.pop()
                if top is None and active():
                    out.append(line)
                continue
        if active():
            out.append(line)
    return out


def translate_file(src: Path, dst: Path) -> None:
    text = src.read_text(errors="replace")
    lines = resolve_gpu_ifdefs(text.splitlines())
    joined = _privatise_nested_loops(_merge_parallel_loop(_join_continuations(lines)))
    joined = _privatise_assigned_scalars(joined)
    out: list[str] = [
        "! acc2omp: generated from " + src.as_posix(),
        "! Do not edit; regenerate via cmake/acc2omp.py",
    ]
    for line, is_acc in joined:
        if is_acc:
            # emit() / enter-data may return multiple physical lines joined by \n
            out.extend(translate_directive(line).splitlines())
        else:
            out.append(sanitize_cuda_fortran(line))
    out = drop_duplicate_decls(out)
    out = hide_use_collisions(out)
    out = inject_fpvs_external(out)
    out = apply_file_fixups(src.name.lower(), out)
    out = relax_intent_in_on_assignment(out)
    out = strip_parameter_maps(out)
    out = _silence_device_io(out)
    text_out = "\n".join(out) + "\n"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and dst.read_text(errors="replace") == text_out:
        return
    dst.write_text(text_out)


def _collect_parameters(lines: list[str]) -> set[str]:
    """Names declared `parameter ::` in this file (module or local)."""
    names: set[str] = set()
    pat = re.compile(r"\bparameter\b[^!\n]*::\s*(.+)$", re.IGNORECASE)
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("!") or stripped.startswith("#"):
            continue
        m = pat.search(line)
        if not m:
            continue
        rest = m.group(1).split("!")[0]
        for item in _split_top_commas(rest):
            nm = re.match(r"([A-Za-z_]\w*)", item.strip())
            if nm:
                names.add(nm.group(1).lower())
    return names


OMP_SENT = re.compile(r"^(\s*)!\$omp(\s|&)(.*)$", re.IGNORECASE)

_DEVICE_WRITE = re.compile(
    r"^(print\b|write\s*\(\s*(?:\*|6)\b)",
    re.IGNORECASE,
)


_DO_START = re.compile(r"^\s*do\b", re.IGNORECASE)
_DO_END = re.compile(r"^\s*end\s*do\b", re.IGNORECASE)


def _silence_device_io(lines: list[str]) -> list[str]:
    """Comment out host I/O inside OpenMP target compute regions.

    Device `print` / `write(6` aborts with Recursive I/O on unit 6 (seen in
    cyclic_cell_ppm_intp). Combined `target teams distribute parallel do` has
    no `end` directive; close it when the associated do-construct ends.
    """
    explicit = 0
    pending_do = 0
    assoc_depth: list[int] = []
    cont_io = False
    out: list[str] = []
    for line in lines:
        m = OMP_SENT.match(line)
        if m:
            rest = m.group(3).strip().lower()
            if rest.startswith("end target data"):
                pass
            elif rest.startswith("end target"):
                explicit = max(0, explicit - 1)
            elif rest.startswith(
                ("target data", "target enter", "target exit", "target update")
            ):
                pass
            elif "distribute parallel do" in rest:
                pending_do += 1
            elif rest.startswith("target"):
                explicit += 1
        else:
            stmt = line.split("!")[0]
            if pending_do and _DO_START.match(stmt):
                assoc_depth.append(1)
                pending_do -= 1
            elif assoc_depth and _DO_START.match(stmt):
                assoc_depth[-1] += 1
            elif assoc_depth and _DO_END.match(stmt):
                assoc_depth[-1] -= 1
                if assoc_depth[-1] <= 0:
                    assoc_depth.pop()
        in_device = explicit > 0 or bool(assoc_depth)
        code = line.lstrip()
        if cont_io or (in_device and _DEVICE_WRITE.match(code)):
            indent = line[: len(line) - len(code)]
            out.append(f"{indent}! acc2omp: device io {code}")
            cont_io = bool(re.search(r"&\s*$", line.rstrip()))
            continue
        out.append(line)
    return out


def _join_omp_continuations(lines: list[str]) -> list[str]:
    out: list[str] = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        m = OMP_SENT.match(line)
        if not m:
            out.append(line)
            i += 1
            continue
        buf = line
        while buf.rstrip().endswith("&") or (
            i + 1 < n and re.match(r"^\s*!\$omp&", lines[i + 1], re.IGNORECASE)
        ):
            if i + 1 < n and lines[i + 1].lstrip().startswith("#"):
                break
            i += 1
            nxt = lines[i]
            nm = OMP_SENT.match(nxt)
            if nm:
                rest = nm.group(3).lstrip()
                if rest.startswith("&"):
                    rest = rest[1:].lstrip()
                buf = re.sub(r"&\s*$", "", buf).rstrip() + " " + rest
            else:
                buf = re.sub(r"&\s*$", "", buf).rstrip() + " " + nxt.strip()
            if not buf.rstrip().endswith("&") and not (
                i + 1 < n and re.match(r"^\s*!\$omp&", lines[i + 1], re.IGNORECASE)
            ):
                break
        out.append(buf)
        i += 1
    return out


def _wrap_omp(line: str) -> list[str]:
    if len(line) <= 132 or not OMP_SENT.match(line):
        return [line]
    indent = OMP_SENT.match(line).group(1)
    chunks: list[str] = []
    s = line
    while len(s) > 132:
        br = s.rfind(",", 0, 120)
        if br < 40:
            br = 120
        chunks.append(s[: br + 1] + " &")
        s = indent + "!$omp& " + s[br + 1 :].lstrip()
    chunks.append(s)
    return chunks


def _rewrite_map_drop_params(line: str, params: set[str]) -> str:
    if "map(" not in line.lower() or not params:
        return line
    out: list[str] = []
    i = 0
    while True:
        m = re.search(r"\bmap\s*\(", line[i:], flags=re.IGNORECASE)
        if not m:
            out.append(line[i:])
            break
        start = i + m.start()
        out.append(line[i:start])
        p0 = i + m.end() - 1
        depth = 0
        k = p0
        while k < len(line):
            if line[k] == "(":
                depth += 1
            elif line[k] == ")":
                depth -= 1
                if depth == 0:
                    break
            k += 1
        inner = line[p0 + 1 : k]
        colon = inner.find(":")
        kind = inner[: colon + 1] if colon >= 0 else ""
        items = inner[colon + 1 :] if colon >= 0 else inner
        kept: list[str] = []
        for item in _split_top_commas(items):
            item = item.strip()
            if not item:
                continue
            nm = re.match(r"([A-Za-z_]\w*)", item)
            if nm and nm.group(1).lower() in params:
                continue
            kept.append(item)
        if kept:
            out.append("map(" + kind + ", ".join(kept) + ")")
        i = k + 1
    return "".join(out).rstrip()


def _rewrite_update_drop_params(line: str, params: set[str]) -> str:
    """PARAMETER names cannot appear in `target update to()/from()`."""
    if "target update" not in line.lower() or not params:
        return line
    out: list[str] = []
    i = 0
    while True:
        m = re.search(r"\b(to|from)\s*\(", line[i:], flags=re.IGNORECASE)
        if not m:
            out.append(line[i:])
            break
        start = i + m.start()
        out.append(line[i:start])
        p0 = i + m.end() - 1
        depth = 0
        k = p0
        while k < len(line):
            if line[k] == "(":
                depth += 1
            elif line[k] == ")":
                depth -= 1
                if depth == 0:
                    break
            k += 1
        inner = line[p0 + 1 : k]
        kept: list[str] = []
        for item in _split_top_commas(inner):
            item = item.strip()
            if not item:
                continue
            nm = re.match(r"([A-Za-z_]\w*)", item)
            if nm and nm.group(1).lower() in params:
                continue
            kept.append(item)
        if kept:
            out.append(m.group(1) + "(" + ", ".join(kept) + ")")
        i = k + 1 if k < len(line) else len(line)
    return "".join(out).rstrip()


def _drop_empty_omp_data(line: str) -> str:
    m = re.match(
        r"^(\s*)!\$omp\s+(target (?:enter|exit) data|target update|target data)\s*$",
        line,
        re.IGNORECASE,
    )
    if m:
        return f"{m.group(1)}! acc2omp: dropped empty {m.group(2)}"
    return line


def strip_parameter_maps(lines: list[str]) -> list[str]:
    """OpenMP map() requires variables; PARAMETER named constants are literals."""
    flat: list[str] = []
    for line in lines:
        flat.extend(line.splitlines())
    joined = _join_omp_continuations(flat)
    params = _collect_parameters(joined)
    out: list[str] = []
    for line in joined:
        if OMP_SENT.match(line):
            line = _rewrite_map_drop_params(line, params)
            line = _rewrite_update_drop_params(line, params)
            line = _drop_empty_omp_data(line)
            out.extend(_wrap_omp(line))
        else:
            out.append(line)
    return out


_DECL_START = re.compile(
    r"^(integer|real|logical|character|complex|double\s+precision|type|class|byte)\b",
    re.IGNORECASE,
)
_INTENT_IN_ONLY = re.compile(r"\bintent\s*\(\s*in\s*\)", re.IGNORECASE)


_SUB_START = re.compile(
    r"^\s*(?:(?:recursive|pure|elemental|module)\s+)*(?:subroutine|function)\b",
    re.IGNORECASE,
)
_SUB_END = re.compile(r"^\s*end\s*(?:subroutine|function)\b", re.IGNORECASE)
_IFACE_START = re.compile(r"^\s*interface\b", re.IGNORECASE)
_IFACE_END = re.compile(r"^\s*end\s*interface\b", re.IGNORECASE)


def _stmt_array_lhs_name(stmt: str) -> str | None:
    """Name of `arr(...) = ...` only. Scalars like `im = nx` are ignored."""
    s = stmt.split("!")[0].strip()
    if not s or s.startswith("#") or s.lower().startswith("!$"):
        return None
    low = s.lower()
    if low.startswith(("call ", "use ", "end", "do ", "if ", "else", "where ", "forall ")):
        return None
    if _DECL_START.match(s):
        return None
    m = re.match(r"([A-Za-z_]\w*)", s)
    if not m:
        return None
    name = m.group(1)
    rest = s[m.end() :].lstrip()
    if not rest.startswith("("):
        return None
    depth = 0
    for i, ch in enumerate(rest):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                after = rest[i + 1 :].lstrip()
                if after.startswith("=") and not after.startswith("=="):
                    return name.lower()
                return None
    return None


def _rewrite_intent_in_chunk(chunk: list[str]) -> list[str]:
    assigned: set[str] = set()
    for line in chunk:
        n = _stmt_array_lhs_name(line)
        if n:
            assigned.add(n)
    if not assigned:
        return chunk
    out: list[str] = []
    for line in chunk:
        stripped = line.lstrip()
        if stripped.startswith("!") or stripped.startswith("#") or "::" not in line:
            out.append(line)
            continue
        before, after = line.split("::", 1)
        if not _INTENT_IN_ONLY.search(before):
            out.append(line)
            continue
        names = []
        for item in _split_top_commas(after.split("!")[0]):
            nm = re.match(r"([A-Za-z_]\w*)", item.strip())
            if nm:
                names.append(nm.group(1).lower())
        if any(n in assigned for n in names):
            line = _INTENT_IN_ONLY.sub("intent(inout)", line, count=1)
        out.append(line)
    return out


def relax_intent_in_on_assignment(lines: list[str]) -> list[str]:
    """nvfortran lets you write INTENT(IN) arrays; amdflang does not.

    Scoped per subprogram, and only for array-section assignments, so
    dimension dummies like `im`/`levs` are not flipped to inout.
    """
    out: list[str] = []
    i = 0
    n = len(lines)
    iface = 0
    while i < n:
        s = lines[i].split("!")[0]
        if _IFACE_START.match(s):
            iface += 1
        elif _IFACE_END.match(s) and iface:
            iface -= 1
        if iface == 0 and _SUB_START.match(s) and not re.match(r"^\s*end\b", s, re.I):
            start = i
            depth = 1
            i += 1
            while i < n and depth:
                inner = lines[i].split("!")[0]
                if _SUB_START.match(inner) and not re.match(r"^\s*end\b", inner, re.I):
                    depth += 1
                elif _SUB_END.match(inner):
                    depth -= 1
                i += 1
            out.extend(_rewrite_intent_in_chunk(lines[start:i]))
            continue
        out.append(lines[i])
        i += 1
    return out


def _strip_subroutine(lines: list[str], name: str) -> list[str]:
    """Drop a whole Fortran subroutine so a handwritten src/rocm copy can own it."""
    start = re.compile(rf"^\s*subroutine\s+{re.escape(name)}\b", re.I)
    end = re.compile(rf"^\s*end\s+subroutine(\s+{re.escape(name)})?\s*$", re.I)
    out: list[str] = []
    skipping = False
    for ln in lines:
        if not skipping and start.search(ln):
            skipping = True
            continue
        if skipping:
            if end.search(ln):
                skipping = False
            continue
        out.append(ln)
    return out


def _extract_subroutine(lines: list[str], name: str) -> list[str]:
    start = re.compile(rf"^\s*subroutine\s+{re.escape(name)}\b", re.I)
    end_named = re.compile(rf"^\s*end\s+subroutine\s+{re.escape(name)}\s*$", re.I)
    iface = re.compile(r"^\s*interface\b", re.I)
    end_iface = re.compile(r"^\s*end\s+interface\b", re.I)
    got: list[str] = []
    grabbing = False
    iface_depth = 0
    for ln in lines:
        if not grabbing:
            if start.search(ln):
                grabbing = True
                got.append(ln)
            continue
        if iface.search(ln):
            iface_depth += 1
        elif end_iface.search(ln):
            iface_depth = max(0, iface_depth - 1)
        got.append(ln)
        if iface_depth == 0 and end_named.search(ln):
            return got
    return got


def _replace_subroutine(lines: list[str], name: str, new_src: Path) -> list[str]:
    """Swap a subroutine body for a handwritten OpenMP copy (keep module wrapper)."""
    new = _extract_subroutine(new_src.read_text().splitlines(), name)
    if not new:
        new = new_src.read_text().splitlines()
    start = re.compile(rf"^\s*subroutine\s+{re.escape(name)}\b", re.I)
    end_named = re.compile(rf"^\s*end\s+subroutine\s+{re.escape(name)}\s*$", re.I)
    iface = re.compile(r"^\s*interface\b", re.I)
    end_iface = re.compile(r"^\s*end\s+interface\b", re.I)
    out: list[str] = []
    skipping = False
    iface_depth = 0
    for ln in lines:
        if not skipping and start.search(ln):
            out.extend(new)
            skipping = True
            iface_depth = 0
            continue
        if skipping:
            if iface.search(ln):
                iface_depth += 1
            elif end_iface.search(ln):
                iface_depth = max(0, iface_depth - 1)
            if iface_depth == 0 and end_named.search(ln):
                skipping = False
            continue
        out.append(ln)
    return out


def apply_file_fixups(name: str, lines: list[str]) -> list[str]:
    """HIP/flang-only source tweaks for NVIDIA dialect that nvfortran accepts."""
    out = lines

    # ---- restore async ordering around RCCL collectives (2026-09-10) -------
    # The NVIDIA sources keep a whole sequence on one async queue, e.g.
    #     !$acc enter data create(rwork) async(async_id)
    #     call nccl_alltoall(sbuf, ..., rwork, ..., async_id)
    #     !$acc parallel loop collapse(4) async(async_id)   ! reads rwork
    #     !$acc exit data delete(rwork)   async(async_id)
    # so everything after the collective is ordered behind it.
    #
    # We drop `async` when translating the OpenACC constructs - harmless on its
    # own - but nccl_* is NOT an OpenACC/OpenMP construct and we do not touch
    # it. It queues ncclSend/ncclRecv on acc_get_cuda_stream(async_id) and
    # returns immediately (src/nvidia/helper.f90). Every translated construct
    # that follows therefore races against the in-flight collective:
    #   - a kernel reading the recv buffer reads it before RCCL fills it
    #   - `exit data delete` frees the send buffer while RCCL is reading it
    #
    # Both are normally invisible: libomptarget's memory manager keeps freed
    # blocks in its pool so addresses stay valid. Reproduce in ~345 s with
    # LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=16777216, which returns blocks to
    # the driver and turns the race into a hard memory access fault.
    #
    # Restore the ordering by waiting on the queue right after each collective.
    # Only calls that actually take async_id are collectives; nccl_check_helper
    # (the NCCLCHECK macro) has no queue and is left alone.
    _NCCL_CALL = re.compile(r"^call\s+nccl_\w+\s*\(", re.IGNORECASE)

    def _flush_nccl(buf: list[str], sink: list[str]) -> None:
        sink.extend(buf)
        if "async_id" in " ".join(buf).lower():
            indent = re.match(r"^(\s*)", buf[0]).group(1)
            sink.append(f"{indent}call geps_acc_wait(async_id)")

    _fixed: list[str] = []
    _buf: list[str] = []
    for _ln in out:
        _st = _ln.strip()
        if not _buf and _NCCL_CALL.match(_st):
            _buf = [_ln]
        elif _buf:
            _buf.append(_ln)
        else:
            _fixed.append(_ln)
            continue
        if not _st.endswith("&"):
            _flush_nccl(_buf, _fixed)
            _buf = []
    if _buf:
        _fixed.extend(_buf)
    out = _fixed
    # Same hazard, third form (2026-09-15): cudaMemsetAsync / cudaMemcpyAsync
    # are queued on the HIP stream and return at once, while the OpenMP
    # kernel that follows runs on libomptarget's queue. In the OpenACC
    # original that kernel was async(async_id) on the SAME stream, so the
    # memset was ordered before it; after translation the memset can finish
    # AFTER the kernel and wipe what it wrote. Measured in tendget_gpu:
    # ndslfv_monoadvh2's pack kernel follows four memsets of its own inputs,
    # and ddtemp came out with 0.48x the CPU's energy in one run and 1.00x in
    # the next (a synchronous probe kernel in between "fixed" it).
    # Insert a device-wide wait after each run of such calls (one wait per
    # run, not per line). The stream id is not always knowable here, hence
    # hipDeviceSynchronize rather than a stream wait.
    _ASYNC_MEM = re.compile(r"^(istat\s*=\s*|CUDACHECK\()?\s*cuda(Memset|Memcpy)Async\s*\(", re.IGNORECASE)
    _fixed = []
    _run = False
    _run_indent = ""
    _pending_cont = False
    for _ln in out:
        _st = _ln.strip()
        if _pending_cont:
            _fixed.append(_ln)
            _pending_cont = _st.endswith("&")
            continue
        if _ASYNC_MEM.match(_st):
            _fixed.append(_ln)
            _run = True
            _run_indent = re.match(r"^(\s*)", _ln).group(1)
            _pending_cont = _st.endswith("&")
            continue
        if _run:
            # first line that is not part of the run: emit the wait before
            # it, whatever it is (a directive counts - never let the wait
            # slip between `!$omp target ...` and its loop)
            _fixed.append(f"{_run_indent}call geps_acc_wait_all()")
            _run = False
        _fixed.append(_ln)
    if _run:
        _fixed.append(f"{_run_indent}call geps_acc_wait_all()")
    out = _fixed
    if "mpe2d_gpu" in name:
        out = [
            re.sub(
                r"\bain\(([^(),]+),([^(),]+),([^(),]+),([^(),]+),([^(),]+)\)",
                r"ain(\1,\2,\3,\4,\5, 1)",
                ln,
            )
            for ln in out
        ]
        # RCCL broadcast is handwritten in src/rocm/mpe2d_row_broadcast_gpu.f90
        # (explicit device ptr, no CUDA Fortran host_data).
        out = _strip_subroutine(out, "mpe2d_row_broadcast_gpu")
    if "moninedmf_gpu" in name:
        fixed = []
        for ln in out:
            # Only rewrite the host REAL declaration list. Wrapped
            # `!$omp ... private(..., kradr)` continuations have no `(`.
            if ln.lstrip().lower().startswith(("!$omp", "!$acc")):
                fixed.append(ln)
                continue
            if (
                re.search(r",\s*kradr\b", ln, re.IGNORECASE)
                and "(" not in ln
            ):
                ln = re.sub(r",\s*kradr\b", "", ln, flags=re.IGNORECASE)
                fixed.append(ln)
                indent = re.match(r"^(\s*)", ln).group(1)
                fixed.append(f"{indent}integer :: kradr")
                continue
            ln = _maybe_pad_rank(ln, "a2g", 1, 1)
            fixed.append(ln)
        out = fixed
    if PROBES and "nnmi_gpu" in name:
        # #region agent log: print the device address OpenMP mapped `wrk` to,
        # from inside nnmi_gpu itself, right before the dgemm loop. Compare it
        # with the C pointer [dgemm] reports - if they differ, hipBLAS and the
        # OpenMP kernel are using two different buffers, which would explain
        # bal implying |wrk|~1 while x_out implies |wrk|~1e-3 (handoff layer 9).
        fixed: list[str] = []
        done = False
        for ln in out:
            if (not done) and "call accx_async_begin_capture(" in ln.lower():
                indent = re.match(r"^(\s*)", ln).group(1)
                fixed.append(f"{indent}! #region agent log")
                fixed.append(f"{indent}call geps_dbg_mapped(wrk, 2)")
                fixed.append(f"{indent}! #endregion")
                done = True
            fixed.append(ln)
        if not done:
            # layer-7 fix removed begin_capture; fall back to the first dgemm
            fixed = []
            for ln in out:
                if (not done) and re.match(r"^\s*call\s+dgemm\s*\(", ln, re.I):
                    indent = re.match(r"^(\s*)", ln).group(1)
                    fixed.append(f"{indent}! #region agent log")
                    fixed.append(f"{indent}call geps_dbg_mapped(wrk, 2)")
                    fixed.append(f"{indent}! #endregion")
                    done = True
                fixed.append(ln)
        out = fixed
    if any(k in name for k in ("nnmi_gpu", "tendget_gpu", "zx_gpu")):
        # ---- layer 7: HIP graph capture drops our OpenMP kernels (2026-09-11)
        # These three routines record a computation into a HIP graph on the
        # first call and replay it afterwards:
        #     if (.not. cg_created) then
        #        <setup: cudaEventCreate, acc_get_cuda_stream, cublasGetHandle,
        #                and in zx_gpu a cudaMallocAsync for `vars`>
        #        call accx_async_begin_capture(async_id)
        #        <compute: hipBLAS dgemm on side streams + !$acc parallel loop>
        #        call accx_async_end_capture(async_id, cg_graph)
        #        cg_created = .true.
        #        <teardown: cudaEventDestroy>
        #     end if
        #     call accx_graph_launch(cg_graph, async_id)
        #
        # We drop `async` when translating, so the OpenMP kernels no longer run
        # on the captured stream. Two things go wrong at once:
        #   1. during capture they execute immediately, reading buffers the
        #      dgemms have only *recorded* writes to - i.e. garbage/zeros
        #   2. they are never recorded, so the replayed graph is missing them
        # In nnmi_gpu that removes the 1/e_t transform between the two dgemms,
        # leaving x <- evec*(evec^T*x) == x: NNMI becomes a no-op, `bal` stays
        # exactly 0 and `pt` never changes, and the integration produces NaN.
        #
        # Fix: drop the graph and run the whole block eagerly on every call.
        #   - the `cg_created` guard becomes unconditional. It has to: events,
        #     streams and the cuBLAS handle live in plain locals with no SAVE,
        #     so they would be undefined on the second call if setup ran once.
        #     Re-running setup is cheap - acc_get_cuda_stream is cached in the
        #     shim and cublasGetHandle is idempotent - and cudaEventDestroy at
        #     the end keeps events balanced.
        #   - zx_gpu is the one exception: it allocates `vars` with
        #     cudaMallocAsync and deliberately never frees it (layer-3 fix), so
        #     that single allocation keeps its own `cg_created` guard,
        #     otherwise it would leak ~510 MB per call.
        #   - graph_launch becomes a wait, keeping its synchronisation role.
        #   - a wait goes before each OpenMP construct in the former capture
        #     region so the async dgemms on the side streams have landed first.
        fixed: list[str] = []
        in_region = False
        pending_malloc_guard = 0
        drop_next_pack_wait = False
        for ln in out:
            low = ln.strip().lower()
            indent = re.match(r"^(\s*)", ln).group(1)
            if re.match(r"^if\s*\(\s*\.not\.\s*cg_created\s*\)\s*then$", low):
                fixed.append(f"{indent}! ROCm layer-7: eager path - this block re-runs on every call")
                fixed.append(f"{indent}if (.true.) then")
                continue
            if "cudamallocasync(vars" in low.replace(" ", ""):
                # keep this one allocation one-shot (layer-3: never freed)
                fixed.append(f"{indent}! ROCm layer-7: `vars` is never freed, so allocate it only once")
                fixed.append(f"{indent}if (.not. cg_created) then")
                fixed.append("   " + ln)
                pending_malloc_guard = 1
                continue
            if pending_malloc_guard and "cudastreamsynchronize" in low.replace(" ", ""):
                fixed.append("   " + ln)
                fixed.append(f"{indent}end if")
                pending_malloc_guard = 0
                continue
            if "call accx_async_begin_capture(" in low:
                in_region = True
                fixed.append(f"{indent}! ROCm layer-7: capture disabled, work runs eagerly")
                continue
            if "call accx_async_end_capture(" in low:
                in_region = False
                fixed.append(f"{indent}! ROCm layer-7: capture disabled, work runs eagerly")
                continue
            if "call accx_graph_launch(" in low:
                fixed.append(f"{indent}! ROCm layer-7: graph replay replaced by a queue wait")
                fixed.append(f"{indent}call geps_acc_wait(async_id)")
                continue
            # `pack_event` is a SINGLE event re-recorded for every ind/m in the
            # loop. Under capture each record/wait pair became its own graph
            # edge ("stream waits for THIS dgemm"), which is correct. Running
            # eagerly the event only ever holds the most recent record, so the
            # wait degenerates into "wait for the last dgemm" and the earlier
            # ones are still in flight - that is the nondeterminism we saw
            # (identical inputs gave bal=2787.0 on one build and 2172.9 on the
            # next). Replace the pair with a direct sync on the stream the
            # dgemm actually ran on.
            m_rec = re.match(
                r"^CUDACHECK\(cudaEventRecord\(pack_event,\s*(.+)\)\)\s*$", ln.strip(), re.I
            )
            if in_region and m_rec:
                fixed.append(
                    f"{indent}! ROCm layer-7: one shared pack_event cannot order a loop; sync directly"
                )
                fixed.append(f"{indent}CUDACHECK(cudaStreamSynchronize({m_rec.group(1)}))")
                drop_next_pack_wait = True
                continue
            if drop_next_pack_wait and "cudastreamwaitevent(stream, pack_event" in low.replace(" ", " "):
                drop_next_pack_wait = False
                continue
            drop_next_pack_wait = False
            if in_region and re.match(r"^!\$omp\s+target\b", low):
                fixed.append(f"{indent}call geps_acc_wait(async_id)")
            fixed.append(ln)
        out = fixed
    if "cufft_wrapper" in name:
        # ---- layer 7, part 1/2: the FFT capture itself (2026-09-14)
        #
        # MEASURED root cause. Inside rfftmlt_loop_identical_cuda_graph the
        # capture region ends with the FFT's normalise-and-write-back kernel:
        #     cc(i, k, jj) = gwk1(i, k, jj)/float(nxj)
        # acc2omp drops its `async`, so:
        #   capture call  - the kernel runs IMMEDIATELY, reading gwk1 while the
        #                   cufftExecs are only RECORDED (gwk1 still the memset
        #                   zeros) -> writes cc = 0; and it is never recorded,
        #                   so the replayed graph lacks it.
        #                   Measured: cc 1.56e18 -> 0, gwk1 = 0, whole chain 0.
        #   replay calls  - the graph holds only cufftExec, so cc is never
        #                   updated at all: the FFT result is stranded in gwk1.
        #                   Measured: cc before == cc after == 1.28e7,
        #                   gwk1 = 9.108e9.
        # BOTH calls are wrong; the second is merely non-zero, not correct.
        # That missing 1/nxj normalisation on un-transformed data is where
        # temten's 3.9e7 factor comes from.
        #
        # Fix: same eager treatment as the accx_* sites - no capture, wait
        # before each kernel so the cufftExecs have landed first.
        fixed: list[str] = []
        in_region = False
        drop_next_pack_wait = False
        for ln in out:
            low = ln.strip().lower()
            indent = re.match(r"^(\s*)", ln).group(1)
            if "cudastreambegincapture(" in low.replace(" ", ""):
                in_region = True
                fixed.append(f"{indent}! ROCm layer-7: FFT capture disabled, runs eagerly")
                continue
            if "cudastreamendcapture(" in low.replace(" ", ""):
                in_region = False
                fixed.append(f"{indent}! ROCm layer-7: FFT capture disabled, runs eagerly")
                continue
            # one shared pack_event re-recorded per jj: under capture each
            # record/wait was its own graph edge, but running eagerly the event
            # only holds the most recent record, so the wait degenerates into
            # "wait for the last cufftExec" and the earlier ones are in flight.
            m_rec = re.match(
                r"^CUDACHECK\(cudaEventRecord\(pack_event,\s*(.+)\)\)\s*$",
                ln.strip(), re.I,
            )
            if in_region and m_rec:
                fixed.append(
                    f"{indent}! ROCm layer-7: one shared pack_event cannot order a loop; sync directly"
                )
                fixed.append(f"{indent}CUDACHECK(cudaStreamSynchronize({m_rec.group(1)}))")
                drop_next_pack_wait = True
                continue
            if drop_next_pack_wait and "cudastreamwaitevent(stream,pack_event" in low.replace(" ", ""):
                drop_next_pack_wait = False
                continue
            drop_next_pack_wait = False
            if in_region and re.match(r"^!\$omp\s+target\b", low):
                fixed.append(f"{indent}call geps_acc_wait(async_id)")
            fixed.append(ln)
        out = fixed
    if "cufft_wrapper" in name:
        # ---- 2026-09-16: one hipfft plan per distinct reduced length, not
        # one per latitude. The per-latitude plan sets (768 x ~22 MB on
        # rocFFT, ~17 GB per (jump, batch, isign)) were the VRAM exhaustion;
        # see src/rocm/rfftmlt_loop_gpu.f90 for the measurement. Spliced
        # AFTER the layer-7 passes so the replacement text is final.
        rocm_fft = Path(__file__).resolve().parent.parent / "src" / "rocm" / "rfftmlt_loop_gpu.f90"
        if rocm_fft.is_file():
            out = _replace_subroutine(out, "rfftmlt_loop_identical_cuda_graph", rocm_fft)
    # ---- 2026-09-16 (layer 17): RAW RCCL collectives. The layer-6 rule waits
    # after `call nccl_*` helpers, but hdiffu_gpu (wmax allreduce),
    # intgrt_gpu:452 (wk4 allreduce -> pcorr) and mpe2d_gpu (allgather /
    # broadcast) call ncclAllReduce/ncclAllGather/ncclBroadcast directly via
    # NCCLCHECK(...). The OpenMP kernel that reads the result follows with no
    # wait. ncclSend/ncclRecv are excluded: they sit inside GroupStart/End
    # (helper.f90) and the caller-level rule already waits after the helper.
    _RAW_NCCL = re.compile(r"^\s*NCCLCHECK\(\s*nccl(AllReduce|AllGather|Broadcast|Reduce|ReduceScatter|Bcast)\s*\(", re.I)
    fixed = []
    pend = False
    for ln in out:
        fixed.append(ln)
        if pend:
            if not ln.rstrip().endswith("&"):
                indent = re.match(r"^(\s*)", ln).group(1)
                fixed.append(f"{indent}call geps_acc_wait_all()   ! layer 17: raw RCCL collective")
                pend = False
            continue
        if _RAW_NCCL.match(ln):
            if ln.rstrip().endswith("&"):
                pend = True
            else:
                indent = re.match(r"^(\s*)", ln).group(1)
                fixed.append(f"{indent}call geps_acc_wait_all()   ! layer 17: raw RCCL collective")
    out = fixed
    if PROBES and "hdiffu_gpu" in name:
        # hdiffu's windchk loop tests the HOST wmax, which is only ever written
        # on the device (wmax = 0 in a kernels region, the max in a parallel
        # loop, exit data delete) - so filter_top fires on stack garbage.
        # Original bug (NVIDIA has it too). CPU semantics: use the real max.
        fixed = []
        for ln in out:
            if re.match(r"^\s*windchk\s*=\s*\.false\.", ln, re.I):
                indent = re.match(r"^(\s*)", ln).group(1)
                fixed.append(f"{indent}!$omp target update from(wmax)")
                fixed.append(f"{indent}! #region agent log")
                fixed.append(f"{indent}print *, 'DBGHD wmax(1:hdk1) max=', maxval(wmax(1:hdk1)), ' hdk1=', hdk1")
                fixed.append(f"{indent}! #endregion")
            fixed.append(ln)
        out = fixed
    if PROBES and "module_mp_gsfcgce_3ice_nuwrf_gpu" in name:
        # ---- probe (2026-09-16, layer 18): fall_flux_gpu's non-SL_sedi
        # sedimentation `DO while (notlast)` kernel never terminates on ROCm
        # (rocgdb: every wave in ..fall_flux_gpu_l625).  Cap the loop, snapshot
        # the state of a stuck column into geps_dbgw (racy, any column will
        # do) and print the used-region range of z1/tz before the kernel.
        fixed = []
        in_ff = False
        for ln in out:
            low = ln.strip().lower()
            if re.match(r"^\s*subroutine\s+fall_flux_gpu\b", ln, re.I):
                in_ff = True
            if in_ff and re.match(r"^\s*end\s+subroutine\s+fall_flux_gpu\b", ln, re.I):
                in_ff = False
            indent = re.match(r"^(\s*)", ln).group(1)
            if in_ff and low.startswith("logical, intent(in) :: benchmark"):
                fixed.append(ln)
                fixed.append(f"{indent}real(8) :: geps_dbgw(24)   ! probe: stuck-column snapshot")
                fixed.append(f"{indent}real(8) :: geps_dbgcol(12, kts:kte)")
                fixed.append(f"{indent}integer :: geps_dbgit, geps_dbgk, geps_dbgflag, geps_dbgold")
                continue
            if in_ff and low == "dtb = dt":
                fixed.append(ln)
                fixed.append(f"{indent}geps_dbgw = 0.d0; geps_dbgcol = 0.d0; geps_dbgflag = 0")
                continue
            if in_ff and low.startswith("!$omp&") and "vtir, xlandr) private(k)" in ln:
                fixed.append(ln.rstrip() + " map(tofrom:geps_dbgw, geps_dbgcol, geps_dbgflag) private(geps_dbgit, geps_dbgk, geps_dbgold)")
                continue
            if in_ff and low == "notlast = .true.":
                fixed.append(ln)
                fixed.append(f"{indent}geps_dbgit = 0")
                continue
            if in_ff and low == "do while (notlast)":
                fixed.append(ln)
                fixed.append(f"{indent}   geps_dbgit = geps_dbgit + 1")
                fixed.append(f"{indent}   if (geps_dbgit .gt. 100000) then")
                fixed.append(f"{indent}      !$omp atomic capture")
                fixed.append(f"{indent}      geps_dbgold = geps_dbgflag")
                fixed.append(f"{indent}      geps_dbgflag = geps_dbgflag + 1")
                fixed.append(f"{indent}      !$omp end atomic")
                fixed.append(f"{indent}      if (geps_dbgold .eq. 0) then")
                fixed.append(f"{indent}      geps_dbgw(1) = nhydro; geps_dbgw(2) = i; geps_dbgw(3) = j")
                fixed.append(f"{indent}      geps_dbgw(4) = del_tv; geps_dbgw(5) = t_del_tv; geps_dbgw(6) = dtb")
                fixed.append(f"{indent}      geps_dbgw(7) = min_q; geps_dbgw(8) = max_q; geps_dbgw(9) = myim(j)")
                fixed.append(f"{indent}      geps_dbgw(10) = its; geps_dbgw(11) = ite; geps_dbgw(12) = jts; geps_dbgw(13) = jte")
                fixed.append(f"{indent}      geps_dbgw(14) = kts; geps_dbgw(15) = kte; geps_dbgw(16) = crmin")
                fixed.append(f"{indent}      do geps_dbgk = kts, kte")
                fixed.append(f"{indent}         geps_dbgcol(1, geps_dbgk) = z1(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(2, geps_dbgk) = rho(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(3, geps_dbgk) = tz(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(4, geps_dbgk) = dz8w(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(5, geps_dbgk) = qr(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(6, geps_dbgk) = qi(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(7, geps_dbgk) = qs(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(8, geps_dbgk) = qg(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(9, geps_dbgk) = vtr(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(10, geps_dbgk) = vts(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(11, geps_dbgk) = vtg(i, geps_dbgk, j)")
                fixed.append(f"{indent}         geps_dbgcol(12, geps_dbgk) = vti(i, geps_dbgk, j)")
                fixed.append(f"{indent}      end do")
                fixed.append(f"{indent}      end if")
                fixed.append(f"{indent}      notlast = .false.")
                fixed.append(f"{indent}      cycle")
                fixed.append(f"{indent}   end if")
                continue
            if in_ff and low.startswith("!$omp target exit data map(release:vtr, vts, vtg, vti, z1, pptall)"):
                fixed.append(f"{indent}if (geps_dbgflag .gt. 0) then")
                fixed.append(f"{indent}print *, 'DBGMP stuck count=', geps_dbgflag, ' nhydro,i,j=', geps_dbgw(1:3), ' del_tv,t_del_tv,dtb=', geps_dbgw(4:6), &")
                fixed.append(f"{indent}   ' min_q,max_q,myim(j)=', geps_dbgw(7:9), ' its,ite,jts,jte,kts,kte=', geps_dbgw(10:15), ' crmin=', geps_dbgw(16)")
                fixed.append(f"{indent}do geps_dbgk = kts, kte")
                fixed.append(f"{indent}   print '(a,i4,12(1x,es12.4))', 'DBGMPCOL k z1 rho tz dz8w qr qi qs qg vtr vts vtg vti', geps_dbgk, geps_dbgcol(1:12, geps_dbgk)")
                fixed.append(f"{indent}end do")
                fixed.append(f"{indent}end if")
                fixed.append(ln)
                continue
            fixed.append(ln)
        out = fixed
        # z1/tz used-region range right after the two init kernels
        fixed = []
        seen_z1_init = 0
        for ln in out:
            fixed.append(ln)
            if "z1(i, 1, j) = 0.9*(z(i, 1, j) - topo(i, j))" in ln:
                seen_z1_init = 1
                continue
            if seen_z1_init == 1 and ln.strip().lower() == "end do":
                seen_z1_init = 2   # closes `do i`
                continue
            if seen_z1_init == 2 and ln.strip().lower() == "end do":
                seen_z1_init = 3   # closes `do j`
                indent = re.match(r"^(\s*)", ln).group(1)
                fixed.append(f"{indent}!$omp target update from(z1, tz, myim)")
                fixed.append(f"{indent}print *, 'DBGMP dims its,ite,jts,jte,kts,kte=', its, ite, jts, jte, kts, kte, ' myim min/max=', minval(myim(jts:jte)), maxval(myim(jts:jte))")
                fixed.append(f"{indent}geps_dbgw(10) = 1.d30; geps_dbgw(11) = -1.d30; geps_dbgw(12) = 1.d30; geps_dbgw(13) = -1.d30; geps_dbgw(9) = 0.d0")
                fixed.append(f"{indent}do j = jts, jte; do k = kts, kte; do i = its, myim(j)")
                fixed.append(f"{indent}   geps_dbgw(10) = min(geps_dbgw(10), dble(z1(i, k, j))); geps_dbgw(11) = max(geps_dbgw(11), dble(z1(i, k, j)))")
                fixed.append(f"{indent}   geps_dbgw(12) = min(geps_dbgw(12), dble(tz(i, k, j))); geps_dbgw(13) = max(geps_dbgw(13), dble(tz(i, k, j)))")
                fixed.append(f"{indent}   if (z1(i, k, j) .ne. z1(i, k, j)) geps_dbgw(9) = geps_dbgw(9) + 1.d0")
                fixed.append(f"{indent}end do; end do; end do")
                fixed.append(f"{indent}print *, 'DBGMP pre z1 min/max=', geps_dbgw(10:11), ' tz min/max=', geps_dbgw(12:13), ' z1 NaN=', geps_dbgw(9), ' dtb=', dtb")
                fixed.append(f"{indent}geps_dbgw = 0.d0")
                continue
        out = fixed
    if any("fft_cg%" in ln for ln in out):
        # ---- layer 7, part 2/2: the FFT call sites (2026-09-14)
        #
        # With the capture removed above, rfftmlt_loop_identical_cuda_graph
        # computes eagerly and produces no graph, so every caller must stop
        # instantiating and replaying one. 18 call sites across 10 files, in
        # two shapes:
        #     if (fft_cg%created) then        <replay>  else <capture> end if
        #     if (.not. (fft_cg%created)) then <capture> end if ; <replay>
        # Forcing the guard makes the real call run every time in both shapes.
        #
        # Scoped to `fft_cg%` on purpose: every file aliases its own FFT graph
        # to that local name (fft_cg => tranrs_fft_cg, ...) while the LT graphs
        # are uniformly `lt_cg`, so this cannot touch them.
        fixed = []
        for ln in out:
            low = ln.strip().lower().replace(" ", "")
            indent = re.match(r"^(\s*)", ln).group(1)
            if low == "if(fft_cg%created)then":
                fixed.append(f"{indent}! ROCm layer-7: FFT runs eagerly; never replay a graph")
                fixed.append(f"{indent}if (.false.) then")
                continue
            if low == "if(.not.(fft_cg%created))then" or low == "if(.not.fft_cg%created)then":
                fixed.append(f"{indent}! ROCm layer-7: FFT runs eagerly on every call")
                fixed.append(f"{indent}if (.true.) then")
                continue
            if "cudagraphinstantiate(fft_cg%" in low:
                fixed.append(f"{indent}! ROCm layer-7: nothing captured, no FFT graph to instantiate")
                continue
            if "cudagraphlaunch(fft_cg%" in low:
                fixed.append(f"{indent}! ROCm layer-7: FFT already executed eagerly")
                continue
            fixed.append(ln)
        out = fixed
    if any("lt_cg%created" in ln for ln in out):
        # ---- layer 7, extended to the raw-CUDA-API capture sites (2026-09-12)
        #
        # 2026-09-15: gate widened from `tranrs_gpu` to every file with an
        # `lt_cg` block. trngra3 / transr / transr1 / trandv / tranuv /
        # rstrandz each hold exactly one block of the same shape (pack_event
        # dgemms + one OpenMP kernel inside the capture) and were left on the
        # broken path; measured consequence in tendget: dlphi/dtphi from
        # trngra3 feed vdzonlr/vdmerdr, and vdzonl/vdmerd came out 1.9x the
        # CPU's after every upstream stage had been shown exact.
        # Same bug as nnmi/tendget/zx, different spelling: this routine opens
        # its capture with cudaStreamBeginCapture instead of the accx_* wrapper,
        # so the original layer-7 pass (which keys on accx_async_begin_capture
        # and a three-name list) never saw it. Measured consequence: `temten`,
        # which tranrs_gpu_cuda_graph writes, comes out 3.6e6x too large and
        # nondeterministic, and that propagates
        #     temten -> temten1 -> phiten1 -> phiten -> x -> bal (8 orders).
        #
        # Scope is deliberately narrow:
        #   - only inside `subroutine tranrs_gpu_cuda_graph`
        #   - only the `lt_cg` graph. The other graph in this file (`fft_cg`,
        #     around L81/L334) has no OpenMP kernel inside its capture region,
        #     so it is left alone.
        # Seven other routines have the same defect (transr, transr1, trandv,
        # tranuv, trngra3, rstrandz, cache_fft_plan) - deliberately NOT touched
        # yet, so that this one change can be judged on its own.
        fixed: list[str] = []
        in_sub = False
        in_region = False
        drop_next_pack_wait = False
        for ln in out:
            low = ln.strip().lower()
            indent = re.match(r"^(\s*)", ln).group(1)
            if re.match(r"^\s*subroutine\s+\w+_gpu_cuda_graph\b", ln, re.I):
                in_sub = True
            elif re.match(r"^\s*end\s+subroutine\b", ln, re.I):
                in_sub = False
                in_region = False
            if in_sub:
                if re.match(r"^if\s*\(\s*\.not\.\s*\(?\s*lt_cg%created\s*\)?\s*\)\s*then$", low):
                    fixed.append(f"{indent}! ROCm layer-7: eager path - this block re-runs on every call")
                    fixed.append(f"{indent}if (.true.) then")
                    continue
                if "cudastreambegincapture(" in low.replace(" ", ""):
                    in_region = True
                    fixed.append(f"{indent}! ROCm layer-7: capture disabled, work runs eagerly")
                    continue
                if "cudastreamendcapture(" in low.replace(" ", ""):
                    in_region = False
                    fixed.append(f"{indent}! ROCm layer-7: capture disabled, work runs eagerly")
                    continue
                if "cudagraphinstantiate(lt_cg%" in low.replace(" ", ""):
                    # nothing was captured, so there is no graph to instantiate.
                    # Must match lt_cg explicitly: this same subroutine also
                    # instantiates fft_cg, whose graph IS still captured (by
                    # rfftmlt_loop_identical_cuda_graph) and must keep working.
                    fixed.append(f"{indent}! ROCm layer-7: nothing captured, no graph to instantiate")
                    continue
                if "cudagraphlaunch(lt_cg%graph_exec" in low.replace(" ", ""):
                    fixed.append(f"{indent}! ROCm layer-7: graph replay replaced by a queue wait")
                    fixed.append(f"{indent}call geps_acc_wait(async_id)")
                    continue
                # one shared pack_event cannot order a loop once we run eagerly
                # (see the comment in the accx_* layer-7 block above)
                m_rec = re.match(
                    r"^istat\s*=\s*cudaEventRecord\(pack_event,\s*(.+)\)\s*$",
                    ln.strip(), re.I,
                )
                if in_region and m_rec:
                    fixed.append(
                        f"{indent}! ROCm layer-7: one shared pack_event cannot order a loop; sync directly"
                    )
                    fixed.append(f"{indent}istat = cudaStreamSynchronize({m_rec.group(1)})")
                    drop_next_pack_wait = True
                    continue
                if drop_next_pack_wait and "cudastreamwaitevent(stream,pack_event" in low.replace(" ", ""):
                    drop_next_pack_wait = False
                    continue
                drop_next_pack_wait = False
                if in_region and re.match(r"^!\$omp\s+target\b", low):
                    fixed.append(f"{indent}call geps_acc_wait(async_id)")
            fixed.append(ln)
        out = fixed
    if PROBES and "tranrs_gpu" in name:
        # #region agent log: presence check, NOT a value probe.
        #
        # In tranrs_gpu_cuda_graph, wcc / fj_wp / wss / wcc_fk are LOCAL
        # automatic arrays (generated file L271-274) and the routine has an
        # `enter data` for twcc_fk ONLY. Yet the dgemm block does
        #     !$omp target data use_device_addr(fj_wp, wcc, wss)
        # and use_device_addr requires its list items to be present in the
        # device data environment. If they are not, the dgemms get host
        # addresses -- which would explain both broken outputs (hldten = 0,
        # temten 3.9e7 too large).
        #
        # Value probes would be useless here: `target update from(wcc)` on an
        # unmapped array silently returns the host copy (handoff 7.1 rule 3,
        # the vorten trap). So ask about presence first.
        #
        # twcc_fk is the POSITIVE CONTROL: it is explicitly mapped on the line
        # this probe is anchored to, so it must report is_present=1. If every
        # array reports 0, the probe itself is broken, not the mapping.
        # #region agent log: value probes for the dgemm's operands and result,
        # recorded once per CALL. tranrs_gpu_cuda_graph runs twice in tendget:
        #   call 1 -> writes hldten (measured as 0,      should be 1.33e8)
        #   call 2 -> writes temten (measured as 3.4e4,  should be 8.6e-4)
        # `wss` is written by the dgemm with beta = 0.0, i.e. wss = wcc x fj_wp,
        # so hldten = 0 means call 1's wss is zero, hence wcc or fj_wp is zero.
        #
        # Deliberately a GPU-vs-GPU comparison (call 1 vs call 2), not GPU-vs-CPU:
        #   - no CPU baseline run needed
        #   - immune to uninitialised padding, because it is the SAME array at
        #     the SAME point with the SAME untouched region both times, so the
        #     padding contribution cancels. That is exactly what was missing
        #     when the whole-array sum of `cc` came back NaN.
        # The two records appear in log order, first = call 1.
        def _ssq_tr(locid: int, varnames: list[str]) -> list[str]:
            body = ["      ! #region agent log"]
            for i, v in enumerate(varnames[:3]):
                body += [
                    f"      !$omp target update from({v})",
                    f"      geps_s = sum(real({v}, kind=8)**2)",
                    "      if (geps_s /= geps_s) then",
                    f"         geps_d{i} = -1_8",
                    "      else if (geps_s .le. 0.0d0) then",
                    f"         geps_d{i} = -2_8",
                    "      else",
                    f"         geps_d{i} = nint(log10(geps_s)*1.0d6, kind=8)",
                    "      end if",
                ]
            body += [
                f"      call geps_dbg_vram(6, {locid}, geps_d0, geps_d1, geps_d2)",
                "      ! #endregion",
            ]
            return body

        _TRIF = [
            "      ! #region agent log",
            "      integer(kind=8) :: geps_d0, geps_d1, geps_d2",
            "      real(kind=8) :: geps_s",
            "      integer, save :: geps_call = 0",
            "      integer :: geps_m, geps_mf",
            "      real(kind=8) :: geps_wc, geps_ws, geps_gwc, geps_gws",
            "      interface",
            "        subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)",
            "          integer hyp, locid",
            "          integer(kind=8) p0, p1, p2",
            "        end subroutine",
            "      end interface",
            "      ! #endregion",
        ]
        # #endregion
        fixed: list[str] = []
        in_cg_sub = False
        tr_iface = False
        for ln in out:
            if re.match(r"^\s*subroutine\s+tranrs_gpu_cuda_graph\b", ln, re.I):
                in_cg_sub = True
            elif re.match(r"^\s*end\s+subroutine\b", ln, re.I):
                in_cg_sub = False
            # declarations go in tranrs_gpu_cuda_graph's `implicit none`, NOT the
            # first subroutine's (this file has two)
            if in_cg_sub and (not tr_iface) and re.match(r"^\s*implicit\s+none", ln, re.I):
                fixed.append(ln)
                fixed.extend(_TRIF)
                tr_iface = True
                continue
            if in_cg_sub and re.match(r"^\s*lt_cg%created\s*=\s*\.true\.", ln, re.I):
                fixed.extend(_ssq_tr(430, ["wcc", "fj_wp", "wss"]))
            # #region agent log: stage-by-stage through the FFT -> pack ->
            # transpose chain that feeds wcc. Measured per CALL (two records
            # each, log order = call order).
            #   call 1 (-> hldten): wcc = 0, so something upstream yields zero
            #   call 2 (-> temten): wcc = 1.254e7, fine
            # The suspicion is that the failure tracks the capture path:
            # fft_cg%created is false on call 1 (capture) and true on call 2
            # (replay). 431 prints that flag so the correlation is measured,
            # not assumed.
            # Note the pack kernel and the NCCL transpose sit OUTSIDE the
            # `if (fft_cg%created)` block, so they run on every call -- wcc_fk
            # actually comes from mpe_transpose_rs_sp_gpu, not straight from
            # the FFT.
            if in_cg_sub and re.match(r"^\s*if\s*\(\s*fft_cg%created\s*\)\s*then", ln, re.I):
                fixed.append("      ! #region agent log")
                fixed.append("      geps_d0 = 0_8")
                fixed.append("      if (fft_cg%created) geps_d0 = 1_8")
                fixed.append("      !$omp target update from(cc)")
                fixed.append("      geps_s = sum(real(cc, kind=8)**2)")
                fixed.append("      if (geps_s /= geps_s) then")
                fixed.append("         geps_d1 = -1_8")
                fixed.append("      else if (geps_s .le. 0.0d0) then")
                fixed.append("         geps_d1 = -2_8")
                fixed.append("      else")
                fixed.append("         geps_d1 = nint(log10(geps_s)*1.0d6, kind=8)")
                fixed.append("      end if")
                fixed.append("      geps_d2 = 0_8")
                fixed.append("      call geps_dbg_vram(6, 431, geps_d0, geps_d1, geps_d2)")
                fixed.append("      ! #endregion")
            if in_cg_sub and "private(jj, jtrunj, mm, mp, mlst)" in ln:
                fixed.extend(_ssq_tr(432, ["cc", "gwk1"]))
            if in_cg_sub and "call mpe_transpose_rs_sp_gpu(" in ln.replace(" ", " "):
                fixed.extend(_ssq_tr(433, ["twcc_fk"]))
            fixed.append(ln)
            if in_cg_sub and "call mpe_transpose_rs_sp_gpu(" in ln:
                fixed.extend(_ssq_tr(434, ["wcc_fk"]))
            if "enter data map(alloc:twcc_fk)" in ln:
                ind = re.match(r"^(\s*)", ln).group(1)
                fixed.append(f"{ind}! #region agent log")
                for tag, v in ((11, "wcc"), (12, "fj_wp"), (13, "wss"),
                               (14, "wcc_fk"), (15, "twcc_fk")):
                    fixed.append(f"{ind}call geps_dbg_mapped({v}, {tag})")
                # collective-call counter, same as the CPU probe's dbg_call
                fixed.append(f"{ind}geps_call = geps_call + 1")
                fixed.append(f"{ind}print *,'DBGTRC call=',geps_call")
                fixed.append(f"{ind}! #endregion")
            # #region agent log: per-wavenumber sums at exit, mirroring the CPU
            # probe in src/tranrs.f90. Post-transpose arrays are split by
            # WAVENUMBER on the CPU, so the only valid comparison is per mf,
            # 1:1 against this single rank. wss is summed over l>=mf only:
            # the dgemm writes exactly that block (beta=0.0) and the CPU
            # zeroes the rest, so anything outside it is not comparable.
            # wcc_fk is padded with zeros on both sides (memset here,
            # twcc_fk=0. there), so the whole latitude range is summed.
            # Everything feeding these has already been synchronised: the
            # transpose is followed by geps_acc_wait, each dgemm's stream is
            # synchronised, and the s=wss kernel is a synchronous target.
            if in_cg_sub and re.match(r"^\s*return\s*$", ln, re.I):
                # this branch runs after fixed.append(ln): pull the `return`
                # back out so the block executes (first version was dead code)
                ret_ln = fixed.pop()
                fixed += [
                    "      ! #region agent log",
                    "      !$omp target update from(wcc_fk)",
                    "      !$omp target update from(wss)",
                    "      geps_gwc = 0.0d0",
                    "      geps_gws = 0.0d0",
                    "      do geps_m = 1, mlistnum",
                    "         geps_mf = mlist(geps_m)",
                    "         geps_wc = sum(real(wcc_fk(:, :, :, geps_m, :), kind=8)**2)",
                    "         geps_ws = sum(real(wss(:, :, :, geps_mf:jtrun, geps_m), kind=8)**2)",
                    "         geps_gwc = geps_gwc + geps_wc",
                    "         geps_gws = geps_gws + geps_ws",
                    "         print *,'DBGTRM',geps_call,geps_mf,geps_wc,geps_ws",
                    "      end do",
                    "      geps_s = sum(real(wss, kind=8)**2)",
                    "      print *,'DBGTRG call=',geps_call,' gwcc_fk=',geps_gwc, &",
                    "              ' gwss_lgem=',geps_gws,' gwss_all=',geps_s",
                    "      ! #endregion",
                    ret_ln,
                ]
            # #endregion
        out = fixed
        # #endregion
    if "zx_gpu" in name:
        # Sequence association: pointer element cannot be dummy array; use section.
        out = [re.sub(r"\bvars\(idx\)\s*,", "vars(idx:),", ln) for ln in out]
        # layer-7 follow-up: `vars` is allocated once (it is never freed, see the
        # layer-3 comment) but is a plain local pointer. That was fine while the
        # address lived inside a captured graph - the routine never had to read
        # `vars` again on later calls. Now that we run eagerly, every call reads
        # it, so it must survive between calls or the second call dereferences
        # garbage (observed: memory access fault at 0x143000).
        out = [
            re.sub(
                r"^(\s*REAL\(kind=RTYPE\), dimension\(:\), pointer)(\s*::\s*vars\s*)$",
                r"\1, save\2",
                ln,
            )
            for ln in out
        ]
    if "get_phi_gpu" in name:
        out = [
            re.sub(
                r"use param,\s*only:\s*my_max,\s*ncld",
                "use param, only: my, my_max, ncld",
                ln,
                flags=re.IGNORECASE,
            )
            for ln in out
        ]
    if PROBES and "initial_gpu" in name:
        # HIP-only probes after check(pt); last ckpt before abort = failing call.
        _DBG_IFACE = [
            "      ! #region agent log",
            "      integer(kind=8) :: dbg0, dbg1, dbg2",
            "      real(kind=8) :: dbgv",
            "      real(kind=8) :: dbgsa, dbgsu",
            "      integer :: dbgi, dbgk, dbgq, dbgnn, dbgnd",
            "      interface",
            "        subroutine geps_dbg_gpu_log(hyp, locid, irank, p0, p1, p2)",
            "          integer hyp, locid, irank",
            "          integer(kind=8) p0, p1, p2",
            "        end subroutine",
            "      end interface",
            "      ! #endregion",
        ]

        def _dbg(hyp: int, locid: int) -> list[str]:
            return [
                "      ! #region agent log",
                f"      dbg0=0_8; dbg1=0_8; dbg2={hyp}_8",
                f"      call geps_dbg_gpu_log({hyp}, {locid}, myrank, dbg0, dbg1, dbg2)",
                "      ! #endregion",
            ]

        # #region agent log: per-anchor occurrence counter so repeated call
        # sites (zx_gpu / vartrix_gpu appear twice) get distinct locids.
        _nseq: dict[str, int] = {}

        def _seq(key: str) -> int:
            n = _nseq.get(key, 0)
            _nseq[key] = n + 1
            return n

        def _amax(locid: int, varnames: list[str]) -> list[str]:
            # #region agent log: dump max|.| of one to three device arrays.
            # Copy back then reduce on the host - simple and immune to the
            # capture/stream problems we are chasing. Values are reported as
            # value*1e6; -1 means NaN, -2 means > 1e12 (would overflow int64).
            body = ["      ! #region agent log"]
            for i, v in enumerate(varnames[:3]):
                body += [
                    f"      !$omp target update from({v})",
                    f"      dbgv = maxval(abs({v}))",
                    "      if (dbgv /= dbgv) then",
                    f"         dbg{i} = -1_8",
                    "      else if (dbgv > 1.0d12) then",
                    f"         dbg{i} = -2_8",
                    "      else",
                    f"         dbg{i} = nint(dbgv*1.0d6, kind=8)",
                    "      end if",
                ]
            for i in range(len(varnames), 3):
                body.append(f"      dbg{i} = 0_8")
            body += [
                f"      call geps_dbg_gpu_log(7, {locid}, myrank, dbg0, dbg1, dbg2)",
                "      ! #endregion",
            ]
            return body

        def _sumx2(locid: int) -> list[str]:
            # #region agent log: Sigma x^2, the quantity the CPU baseline
            # (src/initial.f90 DBGDIST, `sum(x(1:ns*2)**2)`) reports, so the
            # two can be compared directly instead of comparing extremes.
            #
            # Two sums, deliberately:
            #   dbgsa  over the WHOLE declared x(no, 2, 2*jtmax*nnmivm)
            #   dbgsu  over the USED region only, x(1:nn, 1:2, ind),
            #          nn = nnlist(ind) -- the analogue of the CPU's ns*2
            # Without both, a large value cannot be told apart from
            # uninitialised padding; that is exactly the trap max|mx| fell
            # into (handoff layer 9 / rule 7.1).
            #
            # nnlist is read from the HOST copy: host and device were already
            # verified identical (handoff, "nnlist 的裝置副本是正確的").
            #
            # Reported as log10(value)*1e6 because the range to cover spans
            # 1e-9 (CPU) to 1e3 (GPU bal); -1 means NaN, -2 means <= 0.
            body = [
                "      ! #region agent log",
                "      !$omp target update from(x)",
                "      dbgsa = 0.0d0",
                "      dbgsu = 0.0d0",
                "      dbgnd = 0",
                "      do dbgq = 1, 2*jtmax*nnmivm",
                "         dbgnn = nnlist(dbgq)",
                "         if (dbgnn .lt. 0) dbgnn = 0",
                "         if (dbgnn .gt. no) dbgnn = no",
                "         do dbgk = 1, 2",
                "            do dbgi = 1, no",
                "               dbgsa = dbgsa + x(dbgi, dbgk, dbgq)**2",
                "            end do",
                "            do dbgi = 1, dbgnn",
                "               dbgsu = dbgsu + x(dbgi, dbgk, dbgq)**2",
                "            end do",
                "         end do",
                "         dbgnd = dbgnd + 2*dbgnn",
                "      end do",
            ]
            for i, v in enumerate(("dbgsa", "dbgsu")):
                body += [
                    f"      if ({v} /= {v}) then",
                    f"         dbg{i} = -1_8",
                    f"      else if ({v} .le. 0.0d0) then",
                    f"         dbg{i} = -2_8",
                    "      else",
                    f"         dbg{i} = nint(log10({v})*1.0d6, kind=8)",
                    "      end if",
                ]
            body += [
                "      dbg2 = int(dbgnd, 8)",
                f"      call geps_dbg_gpu_log(7, {locid}, myrank, dbg0, dbg1, dbg2)",
                "      ! #endregion",
            ]
            return body

        def _sumsq(locid: int, varnames: list[str]) -> list[str]:
            # #region agent log: Sigma(.)^2 for up to three device arrays, to
            # bisect where x's magnitude blows up between tendget_gpu and
            # vartrix_gpu. Same log10(.)*1e6 encoding as _sumx2 (the range to
            # cover spans ~10 orders); -1 means NaN, -2 means <= 0.
            #
            # Sums are over the WHOLE declared array. That is only meaningful
            # because the CPU probe sums the same declared extents, so the two
            # are compared like for like -- do not read an absolute value here
            # as "the physics", only the CPU/GPU ratio.
            #
            # Deliberately a SUM, not maxval: on this layer the extremes agree
            # to within 2.5e4 while the sums differ by ~9 orders, so extremes
            # say nothing about the distribution (handoff 7.1).
            body = ["      ! #region agent log"]
            for i, v in enumerate(varnames[:3]):
                body += [
                    f"      !$omp target update from({v})",
                    f"      dbgsa = sum(real({v}, kind=8)**2)",
                    "      if (dbgsa /= dbgsa) then",
                    f"         dbg{i} = -1_8",
                    "      else if (dbgsa .le. 0.0d0) then",
                    f"         dbg{i} = -2_8",
                    "      else",
                    f"         dbg{i} = nint(log10(dbgsa)*1.0d6, kind=8)",
                    "      end if",
                ]
            for i in range(len(varnames), 3):
                body.append(f"      dbg{i} = 0_8")
            body += [
                f"      call geps_dbg_gpu_log(7, {locid}, myrank, dbg0, dbg1, dbg2)",
                "      ! #endregion",
            ]
            return body

        def _nnlist_guard(locid: int) -> list[str]:
            # #region agent log: nnlist is declared nnlist(2*jtmax*nnmivm) and
            # NEVER initialised; the fill loop writes only
            #   j = 1 + (m-1)*2 + (L-1)*2*jtmax,  m = 1..mlistnum
            # i.e. 2*mlistnum entries out of every 2*jtmax-long L block. If
            # mlistnum < jtmax the rest keep whatever was in memory.
            #
            # That matches the observation that N = sum(2*nnlist(ind)) changed
            # between two otherwise identical runs (1334076 -> 1335248), and
            # gives the bal kernel a garbage `nn` to walk wrk(i,...) with --
            # explaining both the magnitude and the run-to-run nondeterminism.
            #
            # This block is BOTH the diagnostic and the fix:
            #   p0 = total entries, p1 = entries the loop can write,
            #   p2 = entries found outside [1, no]  (then forced to 0, which
            #        the bal kernel skips because it guards nn > 2).
            # p0 > p1 alone proves the hole, whatever the garbage happens to be.
            return [
                "      ! #region agent log",
                "      dbg0 = int(2*jtmax*nnmivm, 8)",
                "      dbg1 = int(2*mlistnum*nnmivm, 8)",
                "      dbgnd = 0",
                "      do dbgq = 1, 2*jtmax*nnmivm",
                "         if (nnlist(dbgq) .lt. 1 .or. nnlist(dbgq) .gt. no) then",
                "            dbgnd = dbgnd + 1",
                "            nnlist(dbgq) = 0",
                "         end if",
                "      end do",
                "      dbg2 = int(dbgnd, 8)",
                f"      call geps_dbg_gpu_log(7, {locid}, myrank, dbg0, dbg1, dbg2)",
                "      ! #endregion",
            ]

        def _before(ln_: str, hyp: int, locid: int) -> None:
            # insert probe BEFORE the call (calls are continued with '&',
            # appending after the first line would break the continuation)
            fixed.pop()
            fixed.extend(_dbg(hyp, locid))
            fixed.append(ln_)
        # #endregion

        fixed: list[str] = []
        iface = False
        for ln in out:
            if (not iface) and re.match(r"^\s*integer\s+async_id\s*$", ln, re.I):
                fixed.append(ln)
                fixed.extend(_DBG_IFACE)
                iface = True
                continue
            fixed.append(ln)
            low = ln.lower().replace(" ", "")
            if re.match(r"^\s*async_id\s*=\s*1\s*$", ln, re.I):
                fixed.extend(_dbg(5, 148))
            elif "callcheck(pt,nx,my,my_max,lab)" in low:
                fixed.extend(_dbg(5, 180))
            elif "targetupdateto(h,nnlist)" in low:
                fixed.extend(_dbg(1, 199))
            elif "targetupdateto(mx)" in low:
                fixed.extend(_dbg(1, 221))
            elif "callmpe2d_row_broadcast_gpu(" in low:
                pass  # nccl wrapper logs hyp B
            elif "calleigen_mx_gpu(" in low:
                fixed.pop()
                # mx is intent(inout): the coefficient matrix goes IN (may be
                # large), eigenvectors come OUT (must be <= 1). Measuring both
                # sides says whether the solve wrote anything at all.
                fixed.extend(_amax(405, ["mx", "eval"]))
                fixed.append(ln)
                fixed.extend(_amax(406, ["mx", "eval"]))
            elif low.startswith("calltendget_gpu("):
                # #region agent log: tendget_gpu's inputs, measured BEFORE the
                # call. phiten comes back from tendget_gpu already too large
                # and nondeterministic (1.885e8 vs 133.6 across two runs, CPU
                # baseline 6.32e-2), so this decides in ONE run whether
                # tendget_gpu PRODUCES that or merely INHERITS it.
                #
                # All eight are mapped at module scope in mod_grid.f90
                # (L68 ut/vt/rdiv/rvor/tt, L78 pt, L79 dtpl/dlpl) -- verified,
                # because `update from` on an unmapped array silently returns
                # the stale host copy (handoff 7.1 rule 3).
                #
                # qt and phi are NOT in initial_gpu's `use grid, only:` list,
                # so they are out of scope here and cannot be measured.
                # pt is the control: earlier work found it bit-identical
                # across runs, so pt stable + others unstable narrows it fast.
                fixed.pop()
                fixed.extend(_dbg(5, 300 + _seq("tendget")))
                fixed.extend(_sumsq(415, ["pt", "ut", "vt"]))
                fixed.extend(_sumsq(416, ["tt", "rdiv", "rvor"]))
                fixed.extend(_sumsq(417, ["dtpl", "dlpl"]))
                fixed.append(ln)
            elif low.startswith("callzx_gpu("):
                n_zx = _seq("zx")
                fixed.pop()
                fixed.extend(_dbg(5, 310 + n_zx))
                if n_zx == 0:
                    # forward transform; nothing but comments sits between
                    # tendget_gpu and here, so this measures tendget_gpu's
                    # output -- the first bisect point for x's blow-up
                    # phiten ONLY: vorten/divten come from `use spec` and are
                    # never in any map clause, so `target update from` on them
                    # is a no-op and the host copy reads as all-zero. The
                    # "Sigma vorten^2 <= 0" from the first bisect run was that
                    # artifact, not a finding.
                    fixed.extend(_sumsq(412, ["phiten"]))
                if n_zx == 1:
                    # backward transform; measure what vartrix(-2) produced
                    fixed.extend(_amax(402, ["vorten", "divten", "phiten"]))
                fixed.append(ln)
            elif low.startswith("callvartrix_gpu("):
                n_vt = _seq("vartrix")
                fixed.pop()
                fixed.extend(_dbg(5, 320 + n_vt))
                if n_vt == 0:
                    # this is the +2 call that BUILDS x; measuring here gives
                    # zx_gpu's output -- the second bisect point
                    fixed.extend(_sumsq(413, ["phiten"]))
                if n_vt == 1:
                    # this is the -2 call, i.e. right after nnmi_gpu returned
                    fixed.extend(_amax(401, ["x"]))
                    # same Sigma x^2, now on the way OUT of nnmi_gpu
                    fixed.extend(_sumx2(408))
                fixed.append(ln)
            elif low.startswith("callinicons("):
                fixed.pop()
                fixed.extend(_nnlist_guard(414))
                fixed.append(ln)
            elif low.startswith("callnnmi_gpu("):
                fixed.pop()
                fixed.extend(_dbg(5, 330 + _seq("nnmi")))
                # eigen inputs + the vector going in
                fixed.extend(_amax(400, ["x", "mx", "eval"]))
                fixed.extend(_sumx2(407))
                # #region agent log: print the device address OpenMP mapped
                # nnmi_buf (= nnmi_gpu's `wrk`) to, so it can be compared with
                # the C pointer hipBLAS dgemm receives via use_device_addr.
                fixed.extend([
                    "      ! #region agent log",
                    "      call geps_dbg_mapped(nnmi_buf, 1)",
                    "      ! #endregion",
                    # Does the device copy of nnlist match the host one? The
                    # bal kernel reads wrk(i,...) only for i <= nnlist(ind);
                    # if the device nnlist is larger or stale the kernel walks
                    # into memory the dgemm never wrote - `enter data create`
                    # does not zero it - which would explain both the absurd
                    # magnitude and the run-to-run variation (handoff layer 9).
                    "      ! #region agent log",
                    "      dbg0 = int(nnlist(1), 8)",
                    "      dbg1 = int(nnlist(2), 8)",
                    "      dbg2 = 0_8",
                    "      call geps_dbg_gpu_log(8, 410, myrank, dbg0, dbg1, dbg2)",
                    "      !$omp target update from(nnlist)",
                    "      dbg0 = int(nnlist(1), 8)",
                    "      dbg1 = int(nnlist(2), 8)",
                    "      dbg2 = 1_8",
                    "      call geps_dbg_gpu_log(8, 411, myrank, dbg0, dbg1, dbg2)",
                    "      ! #endregion",
                ])
                fixed.append(ln)
            elif low.startswith("callcorrect_gpu("):
                fixed.pop()
                # after zx_gpu(evectr,...): the correction as spectral tendency
                fixed.extend(_amax(403, ["vorten", "divten", "phiten"]))
                fixed.append(ln)
            elif low.startswith("calltransr_gpu_cuda_graph(") and _seq("transr") == 0:
                fixed.pop()
                # after correct_gpu: did the increment reach the spectral field?
                fixed.extend(_amax(404, ["vornow", "divnow", "temnow"]))
                fixed.append(ln)
            elif low.startswith("calltrngra_gpu("):
                _before(ln, 5, 340 + _seq("trngra"))
            elif low.startswith("calltranuv_gpu_cuda_graph("):
                _before(ln, 5, 350 + _seq("tranuv"))
            elif "map(alloc:rvor" in low:
                fixed.extend(_dbg(5, 168))
        out = fixed
    if "ndslfv_pack_gpu" in name:
        # Fat PPM / intpx temps live in src/rocm (OTILE maps).
        out = _strip_subroutine(out, "cyclic_cell_ppm_intp_two_loops_gpu")
        out = _strip_subroutine(out, "cyclic_cell_intpx_jlist_gpu")
        out = _strip_subroutine(out, "cyclic_cell_massadvy_mylonlen_gpu")
        _VRAM_IFACE = [
            "      ! #region agent log",
            "      integer(kind=8) :: dbg0, dbg1, dbg2",
            "      interface",
            "        subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)",
            "          integer hyp, locid",
            "          integer(kind=8) p0, p1, p2",
            "        end subroutine",
            "      end interface",
            "      ! #endregion",
        ]

        def _vram(hyp: int, locid: int, a: str, b: str, c: str) -> list[str]:
            return [
                "      ! #region agent log",
                f"      dbg0=int({a},8); dbg1=int({b},8); dbg2=int({c},8)",
                f"      call geps_dbg_vram({hyp}, {locid}, dbg0, dbg1, dbg2)",
                "      ! #endregion",
            ]

        want = {
            "cyclic_cell_massadvx_jlist_gpu",
            "cyclic_cell_massadvy_mylonlen_gpu",
        }
        fixed: list[str] = []
        cur = ""
        iface_done: set[str] = set()
        for ln in out:
            msub = re.match(r"^\s*subroutine\s+(\w+)", ln, re.I)
            if msub:
                cur = msub.group(1).lower()
            low = ln.lower().replace(" ", "")
            if (
                cur in want
                and cur not in iface_done
                and re.match(r"^\s*integer\s*::\s*async_id\s*$", ln, re.I)
            ):
                fixed.append(ln)
                fixed.extend(_VRAM_IFACE)
                iface_done.add(cur)
                continue
            if "map(alloc:ds,xreg,dist,step,nstep,xpast" in low:
                fixed.extend(_vram(2, 130, "lonfull", "nvars", "levs"))
                fixed.append(ln)
                fixed.extend(_vram(2, 131, "lonfull", "nvars", "jlistnum"))
                continue
            if "map(alloc:gglati_dup" in low:
                fixed.extend(_vram(2, 240, "latfull", "nvars", "mylonlen"))
                fixed.append(ln)
                continue
            if "map(alloc:has_error)" in low:
                fixed.extend(_vram(1, 500, "lonn", "nvars", "outer_size"))
                fixed.append(ln)
                continue
            if "map(alloc:hh,dqmono,qi)" in low:
                fixed.extend(_vram(1, 663, "lonn", "nvars", "inner_size"))
                fixed.append(ln)
                fixed.extend(_vram(1, 664, "lonn", "nvars", "outer_size"))
                continue
            if "map(alloc:dql_array,dqq_array)" in low:
                fixed.extend(_vram(1, 718, "lonn", "nvars", "outer_size"))
                fixed.append(ln)
                continue
            fixed.append(ln)
        out = fixed
    if "ndslfv_monoadv" in name:
        # Full-londim vertical temps OOM on 1-rank HIP; splice tiled OpenMP.
        root = Path(__file__).resolve().parent.parent / "src" / "rocm"
        vert = root / "vertical_cell_advect_gpu.f90"
        pack = root / "ndslfv_monoadvv_gpu.f90"
        if vert.is_file():
            out = _replace_subroutine(out, "vertical_cell_advect_gpu", vert)
        if pack.is_file():
            out = _replace_subroutine(out, "ndslfv_monoadvv_gpu", pack)
            out = _replace_subroutine(out, "ndslfv_monoadvv_fgnl_gpu", pack)
    if "intgrt_gpu" in name:
        # ---- 2026-09-16: the full-step trandv call passes the UNPACKED
        # poly/dpoly (CPU layout (jtrun,my/2,jtmax)) where the routine wants
        # the packed polyf layout; the half-step call (and every other
        # *_cuda_graph call site) passes polyf. On NVIDIA this was masked:
        # fj_wp is computed inside the captured graph, which replays the
        # pointer captured on the FIRST call (polyf). Running eagerly (layer
        # 7) the second call really reads `poly`. Upstream should fix
        # src/nvidia/intgrt_gpu.f90:1085; patched here to keep src/nvidia
        # untouched.
        fixed = []
        n_fix = 0
        for ln in out:
            if "onocos, poly, dpoly, vorten, divten" in ln:
                # note goes before the `call` line (the previous element),
                # not between continuation lines
                fixed.insert(len(fixed) - 1,
                             "      ! ROCm: packed polyf/dpolyf, see acc2omp (trandv reads them eagerly)")
                ln = ln.replace("onocos, poly, dpoly, vorten, divten", "onocos, polyf, dpolyf, vorten, divten")
                n_fix += 1
            fixed.append(ln)
        assert n_fix == 1, f"trandv poly/dpoly call site: expected 1, found {n_fix}"
        out = fixed
        # #region agent log: integration step-1 checkpoints, shared table with
        # src/intgrt.f90 (cmake/intgrt_ckpts.py) so tags match 1:1.
        if PROBES:
            import importlib.util as _ilu
            _spec = _ilu.spec_from_file_location("intgrt_ckpts", str(Path(__file__).resolve().parent / "intgrt_ckpts.py"))
            _ck = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_ck)
            out, _counts = _ck.inject(out, 1, gpu=True)
            for _e in _ck.CK:
                assert _counts.get(_e[1], 0) >= _e[2], f"intgrt_gpu checkpoint anchor not found: {_e[1]} #{_e[2]}"
        # #endregion
    if PROBES and "diabat_gpu" in name:
        # #region agent log: physics-chain checkpoints (layer 18), table in
        # cmake/diabat_ckpts.py; CPU side src/diabat.f90 carries p0/p1 only.
        import importlib.util as _ilu
        _spec = _ilu.spec_from_file_location("intgrt_ckpts", str(Path(__file__).resolve().parent / "intgrt_ckpts.py"))
        _ck = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_ck)
        _spec2 = _ilu.spec_from_file_location("diabat_ckpts", str(Path(__file__).resolve().parent / "diabat_ckpts.py"))
        _dk = _ilu.module_from_spec(_spec2); _spec2.loader.exec_module(_dk)
        out, _counts = _ck.inject(out, 1, gpu=True, ck=_dk.CK)
        for _e in _dk.CK:
            assert _counts.get(_e[1], 0) >= _e[2], f"diabat_gpu checkpoint anchor not found: {_e[1]} #{_e[2]}"
        # #endregion
    if PROBES and "pbl_noah_gpu" in name:
        # #region agent log: layer-18 stage checkpoints inside pbl_noah_gpu
        import importlib.util as _ilu
        _spec = _ilu.spec_from_file_location("intgrt_ckpts", str(Path(__file__).resolve().parent / "intgrt_ckpts.py"))
        _ck = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_ck)
        _spec3 = _ilu.spec_from_file_location("pbl_ckpts", str(Path(__file__).resolve().parent / "pbl_ckpts.py"))
        _pk = _ilu.module_from_spec(_spec3); _spec3.loader.exec_module(_pk)
        out, _counts = _ck.inject(out, 1, gpu=True, ck=_pk.CK)
        for _e in _pk.CK:
            assert _counts.get(_e[1], 0) >= _e[2], f"pbl_noah_gpu checkpoint anchor not found: {_e[1]} #{_e[2]}"
        # #endregion
    if PROBES and "sflx_gpu" in name:
        # #region agent log (layer 18, 2026-09-17): Noah gives Inf skin temp on
        # ~30k columns, count varies run to run -> a read-before-write of some
        # private scalar (garbage registers on AMD).  Seed every private real
        # scalar with NaN / every private integer with -999999 at the top of
        # the column body, and dump all of them plus the column's inputs for
        # the first column whose t1 comes out non-finite.  A scalar still NaN
        # at the end was never written on that path.
        src_lines = out
        # 1. private scalar names, typed from the declaration section
        di = next(i for i, l in enumerate(src_lines) if "!$omp target teams" in l and "private(bexp" in l)
        dend = di
        while src_lines[dend].rstrip().endswith("&"):
            dend += 1
        directive = " ".join(l.split("!$omp&")[-1] if k else l for k, l in enumerate(src_lines[di:dend + 1]))
        pnames = _private_names(directive)
        decls = {}
        k = 0
        while k < di:
            l = src_lines[k].split("!")[0].rstrip()
            if re.match(r"^\s*(real|integer|logical)\b", l, re.I) and "::" in l:
                full = l
                while full.rstrip().endswith("&"):
                    k += 1
                    full = full.rstrip().rstrip("&") + " " + src_lines[k].split("!")[0].strip()
                head, _, ents = full.partition("::")
                if "dimension" in head.lower():
                    k += 1; continue
                ty = re.match(r"^\s*(real|integer|logical)", head, re.I).group(1).lower()
                for e in _split_top(ents):
                    m = re.match(r"\s*([A-Za-z_]\w*)\s*(\(|$)", e.strip() + " ")
                    if m and m.group(2) != "(":
                        decls[m.group(1).lower()] = ty
            k += 1
        reals = [n for n in sorted(pnames) if decls.get(n) == "real"]
        ints = [n for n in sorted(pnames) if decls.get(n) == "integer" and n not in ("i", "jj")]
        in2d = ["ffrozp", "zlvl", "lwdn", "swdn", "swnet", "sfcprs", "sfctmp", "sfcems", "sfcspd", "prcp", "q2",
                "q2sat", "dqsdt2", "th2", "shdmin", "alb", "snoalb", "tbot", "cmc", "t1", "sneqv", "ch", "cm", "z0",
                "shdfac", "snowh", "sncovr", "etp", "ssoil", "beta", "eta", "sheat", "flx1", "flx2", "flx3", "snomlt"]
        in2i = ["couple", "icein", "vegtyp", "soiltyp", "slopetyp", "nroot"]
        nr = len(reals) + len(in2d) + 12
        ni = len(ints) + len(in2i) + 3
        fixed = []
        for i, l in enumerate(src_lines):
            st = l.strip()
            if st.startswith("data snoexp/2.0/"):
                fixed.append(f"         real(kind=kind_phys) :: geps_snan, geps_dbgs({nr})")
                fixed.append(f"         integer :: geps_dbgi({ni}), geps_dbgflag, geps_dbgold")
                fixed.append(l); continue
            if i == di:
                fixed.append("         geps_snan = transfer(9221120237041090560_8, 1.0_8)")
                fixed.append("         geps_dbgs = 0.0_8; geps_dbgi = 0; geps_dbgflag = 0")
                fixed.append(l); continue
            if i == dend:
                fixed.append(l.rstrip() + " map(tofrom: geps_dbgs, geps_dbgi, geps_dbgflag) private(geps_dbgold)")
                continue
            if st == "if (flag(i, jj) .and. flag_iter(i, jj)) then" and i > di:
                fixed.append(l)
                for n in reals:
                    fixed.append(f"         {n} = geps_snan")
                for n in ints:
                    fixed.append(f"         {n} = -999999")
                continue
            if st == "soilw(i, jj) = soilww/soilwm":
                fixed.append(l)
                fixed.append("         if (t1(i, jj) /= t1(i, jj) .or. abs(t1(i, jj)) > 1.0d10) then")
                fixed.append("            !$omp atomic capture")
                fixed.append("            geps_dbgold = geps_dbgflag")
                fixed.append("            geps_dbgflag = geps_dbgflag + 1")
                fixed.append("            !$omp end atomic")
                fixed.append("            if (geps_dbgold .eq. 0) then")
                kk = 0
                for n in reals:
                    kk += 1; fixed.append(f"               geps_dbgs({kk}) = {n}")
                for n in in2d:
                    kk += 1; fixed.append(f"               geps_dbgs({kk}) = {n}(i, jj)")
                for n in ("stc", "smc", "sh2o"):
                    for kz in range(1, 5):
                        kk += 1; fixed.append(f"               geps_dbgs({kk}) = {n}(i, {kz}, jj)")
                kk = 0
                for n in ints:
                    kk += 1; fixed.append(f"               geps_dbgi({kk}) = {n}")
                for n in in2i:
                    kk += 1; fixed.append(f"               geps_dbgi({kk}) = {n}(i, jj)")
                fixed.append(f"               geps_dbgi({kk+1}) = i; geps_dbgi({kk+2}) = jj; geps_dbgi({kk+3}) = myim(jj)")
                fixed.append("            end if")
                fixed.append("         end if")
                continue
            if st.startswith("return") and i > dend and not any("DBGSFLX" in x for x in fixed):
                fixed.append("         if (geps_dbgflag .gt. 0) then")
                fixed.append("            print *, 'DBGSFLX count=', geps_dbgflag")
                labels_r = reals + [n + "(i,jj)" for n in in2d] + [f"{n}({kz})" for n in ("stc", "smc", "sh2o") for kz in range(1, 5)]
                labels_i = ints + [n + "(i,jj)" for n in in2i] + ["i", "jj", "myim(jj)"]
                for kk, lab in enumerate(labels_r, 1):
                    fixed.append(f"            print '(a,a,es16.6)', 'DBGSFLX r ', '{lab}=', geps_dbgs({kk})")
                for kk, lab in enumerate(labels_i, 1):
                    fixed.append(f"            print '(a,a,i12)', 'DBGSFLX i ', '{lab}=', geps_dbgi({kk})")
                fixed.append("         end if")
                fixed.append(l); continue
            fixed.append(l)
        out = fixed
        # #endregion
    if "adjptqintp_gpu" in name:
        # ---- layer 18b (2026-09-17): out-of-bounds read plold(i, 0, jj) at
        # k = 1 (plold is (nxp, lev+1, my_max)); hfdpr1 is only used for
        # k >= 2. Harmless on NVIDIA/CPU (reads the previous column), but with
        # the pool off every buffer starts on a fresh page and both ranks
        # faulted at <base - 8*nxp> (address ...fc000). Guard the read.
        fixed = []
        n_fix = 0
        for ln in out:
            if ln.strip() == "hfdpr1 = 0.5*(plold(i, k, jj) - plold(i, k - 1, jj))":
                indent = re.match(r"^(\s*)", ln).group(1)
                fixed.append(f"{indent}hfdpr1 = 0.0")
                fixed.append(f"{indent}if (k .ge. 2) hfdpr1 = 0.5*(plold(i, k, jj) - plold(i, k - 1, jj))   ! acc2omp: no plold(i,0,jj) read")
                n_fix += 1
                continue
            fixed.append(ln)
        assert n_fix == 1, f"adjptqintp_gpu hfdpr1 site: expected 1, found {n_fix}"
        out = fixed
    if PROBES and "samfdeepcnv_kh_gpu" in name:
        # #region agent log (2026-09-17): deep convection triggers in 0 columns
        # on ROCm (CPU: 29009). Count cnvflg before every kernel to see which
        # test switches every column off. Tag = generated line of the kernel.
        fixed = []
        seen_first = False
        for ln in out:
            if ln.lstrip().startswith("!$omp target teams"):
                if seen_first:
                    indent = re.match(r"^(\s*)", ln).group(1)
                    fixed.append(f"{indent}!$omp target update from(cnvflg)")
                    fixed.append(f"{indent}call geps_dbg_cntl_grid(cnvflg, ix, my_max, 'g_dc_{len(fixed)+1:05d}')")
                seen_first = True
            fixed.append(ln)
        out = fixed
        # Second probe: the LFC-search kernel switches every column off.
        # Dump column (100,100) of rank 0 right after the LFC test: the
        # indices and the heo/heso/phil/to/qo/pfld profiles.
        fixed = []
        for ln in out:
            st = ln.strip()
            if st.startswith("real(kind=kind_phys) cinpcr,  cinpcrmx,  cinpcrmn, &"):
                fixed.append("      real(kind=kind_phys) :: geps_dbgp(6, km), geps_dbgi(8)   ! probe")
                fixed.append(ln); continue
            if ln.lstrip().startswith("!$omp target teams") and "private(flgr, sigma_con, cinpcr" in ln:
                fixed.append("      geps_dbgp = 0.0_8; geps_dbgi = 0.0_8")
                fixed.append(ln.rstrip() + (" &" if not ln.rstrip().endswith("&") else ""))
                continue
            if st == "if(kbconr == kmaxr) cnvflg(i, jj) = .false.":
                fixed.append(ln)
                fixed.append("               if (i .eq. 100 .and. jj .eq. 100) then")
                fixed.append("                  geps_dbgi(1) = kmaxr; geps_dbgi(2) = kbr; geps_dbgi(3) = kbmaxr; geps_dbgi(4) = kbconr")
                fixed.append("                  geps_dbgi(5) = heo(i, kbr, jj); geps_dbgi(6) = kbm(i, jj); geps_dbgi(7) = ps(i, jj); geps_dbgi(8) = gdx(i, jj)")
                fixed.append("                  do k = 1, km")
                fixed.append("                     geps_dbgp(1, k) = heo(i, k, jj); geps_dbgp(2, k) = heso(i, k, jj); geps_dbgp(3, k) = phil(i, k, jj)")
                fixed.append("                     geps_dbgp(4, k) = to(i, k, jj); geps_dbgp(5, k) = qo(i, k, jj); geps_dbgp(6, k) = pfld(i, k, jj)")
                fixed.append("                  end do")
                fixed.append("               end if")
                continue
            if st.startswith("call geps_dbg_cntl_grid(cnvflg, ix, my_max, 'g_dc_00966')"):
                fixed.append(ln)
                # NB: `myrank` here is a dummy argument the caller does NOT pass
                # (34 dummies vs 29 actuals, no interface) -> garbage; use the
                # module rank, which this routine renames to geps_hide_myrank.
                fixed.append("      if (geps_hide_myrank .eq. 0) then")
                fixed.append("         print *, 'DBGLFC kmax,kb,kbmax,kbcon,heo(kb),kbm,ps,gdx=', geps_dbgi")
                fixed.append("         do k = 1, km")
                fixed.append("            print '(a,i3,6(1x,es14.6))', 'DBGLFC k heo heso phil to qo pfld ', k, geps_dbgp(1:6, k)")
                fixed.append("         end do")
                fixed.append("      end if")
                continue
            fixed.append(ln)
        out = fixed
        # add the map clause on the directive's last continuation line: the
        # directive was extended with a trailing '&' above, so append a line
        fixed = []
        pend = False
        for ln in out:
            if pend and ln.lstrip().startswith("!$omp&"):
                fixed.append(ln)
                if not ln.rstrip().endswith("&"):
                    fixed[-1] = ln.rstrip() + " &"
                    fixed.append("      !$omp& map(tofrom: geps_dbgp, geps_dbgi)")
                    pend = False
                continue
            if ln.lstrip().startswith("!$omp target teams") and "private(flgr, sigma_con, cinpcr" in ln:
                pend = True
            fixed.append(ln)
        out = fixed
        # #endregion
    if PROBES and "intgrt_gpu" in name or "tendget_gpu" in name:
        _VIF = [
            "      ! #region agent log",
            "      integer(kind=8) :: geps_d0, geps_d1, geps_d2",
            "      real(kind=8) :: geps_s",
            "      interface",
            "        subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)",
            "          integer hyp, locid",
            "          integer(kind=8) p0, p1, p2",
            "        end subroutine",
            "      end interface",
            "      ! #endregion",
        ]

        def _ssq(locid: int, varnames: list[str]) -> list[str]:
            # #region agent log: Sigma(.)^2 along tendget_gpu's compute chain.
            # CPU/GPU comparison put the 1.1e8 amplification inside this
            # routine: its eight inputs are all ~32x CPU rank 0 (the expected
            # structural factor, i.e. correct) while phiten comes out 1.1e8
            # too large. These probes say WHICH step does it.
            #
            # log10(.)*1e6 encoding, same as _sumsq in the initial_gpu section
            # (range to cover is ~12 orders); -1 = NaN, -2 = <= 0.
            # Sums are whole-array, so only GPU-vs-CPU ratios of the SAME
            # quantity mean anything -- never compare two different arrays'
            # magnitudes against each other (handoff 7.1 rule 5).
            body = ["      ! #region agent log"]
            for i, v in enumerate(varnames[:3]):
                body += [
                    f"      !$omp target update from({v})",
                    f"      geps_s = sum(real({v}, kind=8)**2)",
                    "      if (geps_s /= geps_s) then",
                    f"         geps_d{i} = -1_8",
                    "      else if (geps_s .le. 0.0d0) then",
                    f"         geps_d{i} = -2_8",
                    "      else",
                    f"         geps_d{i} = nint(log10(geps_s)*1.0d6, kind=8)",
                    "      end if",
                ]
            for i in range(len(varnames), 3):
                body.append(f"      geps_d{i} = 0_8")
            body += [
                f"      call geps_dbg_vram(6, {locid}, geps_d0, geps_d1, geps_d2)",
                "      ! #endregion",
            ]
            return body

        def _vv(locid: int) -> list[str]:
            return [
                "      ! #region agent log",
                f"      geps_d0=0_8; geps_d1=0_8; geps_d2={locid}_8",
                f"      call geps_dbg_vram(4, {locid}, geps_d0, geps_d1, geps_d2)",
                "      ! #endregion",
            ]

        fixed: list[str] = []
        iface = False
        pending_244 = False
        # #region agent log: occurrence counters. `call joinrs_gpu(` and
        # `call tranrs_gpu_cuda_graph(` each appear TWICE in tendget_gpu:
        #   1st pair (L313-316): joinrs(cc_cg, diveng) -> tranrs(cc_cg, hldten)
        #   2nd pair (L371-375): joinrs(cc_cg, ddtemp) -> tranrs(cc_cg, temten)
        # Only the 2nd writes `temten`, the quantity under investigation, and
        # cc_cg is REUSED between the two pairs.
        n_joinrs = 0
        n_tranrs = 0
        _pending436 = [False]
        _pending_dt = [False]
        _pending_hadv = [False]
        _f2p_cont = [False]
        _after_trngra3 = [False]
        _hadv_cont = [False]
        # #endregion
        for ln in out:
            if (not iface) and re.match(r"^\s*implicit\s+none", ln, re.I):
                fixed.append(ln)
                fixed.extend(_VIF)
                iface = True
                continue
            low = ln.lower().replace(" ", "")
            # #region agent log: four probes along tendget_gpu's compute chain.
            # Anchors are POST-layer-7 text (this section runs after the
            # eager-execution rewrite) and each was verified unique in the
            # generated file.
            #   420 before cudaEventCreate   inputs: plten, spalm, arrhyd
            #   421 same point               input:  temten1
            #   422 before cudaEventRecord   phiten1 AFTER the init loop
            #                                (eager => the loop has landed)
            #   423 before the copy loop     phiten1 AFTER the dgemms
            #                                <- if it jumps 1e8 here, the
            #                                   beta=ONE accumulation is it
            #   424 before cg_created=.true. phiten AFTER the copy
            if "tendget_gpu" in name and n_joinrs == 1 and _pending436[0]:
                _pending436[0] = False
                fixed.extend(_ssq(436, ["cc_cg"]))
            if "tendget_gpu" in name:
                # #region agent log: is `temten` already bad when tranrs gets
                # it, or does tranrs make it bad? Measuring cc_cg on both sides
                # of joinrs and temten on both sides of tranrs answers that in
                # one run. ddtemp is deliberately NOT measured: its enter data
                # lives in intgrt_gpu, which runs AFTER initial_gpu, so it may
                # not be mapped here and `update from` would silently hand back
                # the host copy (handoff 7.1 rule 3, the vorten trap).
                if "calljoinrs_gpu(" in low:
                    n_joinrs += 1
                    if _pending_dt[0]:
                        # ddtemp after (ddtemp - ttp)/dt, i.e. what joinrs reads
                        _pending_dt[0] = False
                        fixed += [
                            "      ! #region agent log",
                            "      !$omp target update from(ddtemp)",
                            "      call geps_dbg_ssq_grid(ddtemp, nxp, lev, my_max, 'ddtemp_post_dt')",
                            "      ! #endregion",
                        ]
                    # #region agent log: the FIRST joinrs feeds tranrs call 1,
                    # whose whole chain came out at 1e15 once the FFT fix
                    # stopped zeroing it. diveng is its input, produced by
                    # gridnl_hybrid_ndsl_gpu_refactor. diveng is mapped locally
                    # (tendget_gpu's own enter data, L122) so it is measurable
                    # here -- unlike ddtemp, whose enter data lives in
                    # intgrt_gpu and has not run yet.
                    if n_joinrs == 1:
                        fixed.extend(_ssq(435, ["diveng", "cc_cg"]))
                        _pending436[0] = True
                    # #endregion
                    if n_joinrs == 2:
                        fixed.extend(_ssq(425, ["cc_cg", "temten"]))
                if "calltranrs_gpu_cuda_graph(" in low:
                    n_tranrs += 1
                    if n_tranrs == 2:
                        fixed.extend(_ssq(426, ["cc_cg", "temten"]))
                    # same marker text as src/tendget.f90 so one analyser
                    # serves both logs (tranrs is a shared routine)
                    fixed.append(f"      print *,'DBGMARK tendget_tranrs_{n_tranrs}'")
                if "calltrngra3_gpu_cuda_graph(" in low:
                    # hldten is the FIRST tranrs pair's output. Same tranrs
                    # code, different input. If hldten is off by the same
                    # ~3.9e7 as temten, tranrs itself is broken; if only
                    # temten is off, the problem is on its input side.
                    # Measured here, before trngra3 consumes it.
                    fixed.extend(_ssq(428, ["hldten"]))
                if "callmpe2d_unify_lev_gpu(" in low:
                    fixed.extend(_ssq(427, ["temten"]))
                # #endregion
                if "cudaeventcreate(spread_event)" in low:
                    fixed.extend(_ssq(420, ["plten", "spalm", "arrhyd"]))
                    fixed.extend(_ssq(421, ["temten1"]))
                elif "cudaeventrecord(spread_event,stream)" in low:
                    fixed.extend(_ssq(422, ["phiten1"]))
                elif "private(mf,kk)" in low:
                    fixed.extend(_ssq(423, ["phiten1"]))
                elif low.startswith("cg_created=.true."):
                    fixed.extend(_ssq(424, ["phiten"]))
            # #endregion
            if "tendget_gpu" in name and "calltrngra3_gpu_cuda_graph(" in low:
                _after_trngra3[0] = True
            if ("tendget_gpu" in name and _after_trngra3[0]
                    and low.startswith("!$omptargetteamsdistributeparalleldocollapse(3)private(j,nxj)")):
                # #region agent log: the vdmerdr/vdzonlr tendency loop is the
                # first such loop after trngra3. Bracket the update (see
                # src/tendget.f90 for the reasoning).
                _after_trngra3[0] = False
                for v, tag in (("vdzonl", "vdzonl_post_f2p"), ("vdmerd", "vdmerd_post_f2p"),
                               ("vdzonlr", "vdzonlr_pre_upd"), ("vdmerdr", "vdmerdr_pre_upd"),
                               ("dlphi", "dlphi_pre_upd"), ("dtphi", "dtphi_pre_upd")):
                    fixed += [
                        "      ! #region agent log",
                        f"      !$omp target update from({v})",
                        f"      call geps_dbg_ssq_grid({v}, nxp, lev, my_max, '{tag}')",
                        "      ! #endregion",
                    ]
                # #endregion
            if "tendget_gpu" in name and "callndslfv_update_gpu(" in low:
                for v, tag in (("vdzonlr", "vdzonlr_post_loop"), ("vdmerdr", "vdmerdr_post_loop")):
                    fixed += [
                        "      ! #region agent log",
                        f"      !$omp target update from({v})",
                        f"      call geps_dbg_ssq_grid({v}, nxp, lev, my_max, '{tag}')",
                        "      ! #endregion",
                    ]
                fixed.append(ln)
                for v, tag in (("vdzonl", "vdzonl_post_upd"), ("vdmerd", "vdmerd_post_upd")):
                    fixed += [
                        "      ! #region agent log",
                        f"      !$omp target update from({v})",
                        f"      call geps_dbg_ssq_grid({v}, nxp, lev, my_max, '{tag}')",
                        "      ! #endregion",
                    ]
                continue
            if "tendget_gpu" in name and "callmpe2d_transpose_ndsl_f2p_gpu(ttm_sl,ddtemp," in low:
                # #region agent log: the first f2p transpose is the stage
                # that loses ~half of ddtemp's energy (its input ttm_sl is
                # exact). The call spans two lines; probe after the
                # continuation closes.
                # Device-vs-host nxjp/jlist1 (the kernel's only data
                # dependent bound): if the device sum differs, the copy-in
                # in gfcst is stale/absent and that alone explains the loss.
                fixed += [
                    "      ! #region agent log",
                    "      geps_d0 = 0_8",
                    "      !$omp target teams distribute parallel do reduction(+:geps_d0)",
                    "      do j = 1, jlistnum",
                    "         geps_d0 = geps_d0 + int(nxjp(jlist1(j)), 8)",
                    "      end do",
                    "      geps_d1 = 0_8",
                    "      do j = 1, jlistnum",
                    "         geps_d1 = geps_d1 + int(nxjp(jlist1(j)), 8)",
                    "      end do",
                    "      print *,'DBGNXJP sum_nxjp device=',geps_d0,' host=',geps_d1, &",
                    "              ' jlistnum=',jlistnum,' nxp=',nxp,' lev=',lev,' levp=',levp,' levf=',levf",
                    "      ! #endregion",
                ]
                fixed.append(ln)
                _f2p_cont[0] = ln.rstrip().endswith("&")
                if not _f2p_cont[0]:
                    fixed += [
                        "      ! #region agent log",
                        "      !$omp target update from(ddtemp)",
                        "      call geps_dbg_ssq_grid(ddtemp, nxp, lev, my_max, 'ddtemp_post_f2p')",
                        "      call geps_dbg_ssq_prof(ddtemp, nxp, lev, my_max, 'ddtemp_post_f2p')",
                        "      ! #endregion",
                    ]
                continue
            if "tendget_gpu" in name and _f2p_cont[0]:
                fixed.append(ln)
                if not ln.rstrip().endswith("&"):
                    _f2p_cont[0] = False
                    fixed += [
                        "      ! #region agent log",
                        "      !$omp target update from(ddtemp)",
                        "      call geps_dbg_ssq_grid(ddtemp, nxp, lev, my_max, 'ddtemp_post_f2p')",
                        "      call geps_dbg_ssq_prof(ddtemp, nxp, lev, my_max, 'ddtemp_post_f2p')",
                        "      ! #endregion",
                    ]
                continue
            if "tendget_gpu" in name and "callndslfv_monoadvh2_gpu_refactor(" in low:
                # #region agent log: bracket the horizontal SL advection (see
                # src/tendget.f90 for the reasoning). pten_sl is written by
                # this call, so its "pre" value is only a padding check.
                def _full(v, tag):
                    return [
                        "      ! #region agent log",
                        f"      !$omp target update from({v})",
                        f"      call geps_dbg_ssq_full({v}, nx, levp, my_max, '{tag}')",
                        "      ! #endregion",
                    ]
                for v in ("ttm_sl", "uum_sl", "vvm_sl", "pten_sl"):
                    fixed += _full(v, f"{v}_pre_hadv")
                fixed.append(ln)
                _pending_hadv[0] = not ln.rstrip().endswith("&")
                _hadv_cont[0] = ln.rstrip().endswith("&")
                if _pending_hadv[0]:
                    for v in ("ttm_sl", "uum_sl", "vvm_sl", "pten_sl"):
                        fixed += _full(v, f"{v}_post_hadv")
                    _pending_hadv[0] = False
                    fixed.extend(_vv(244))   # keep the pre-existing marker
                continue
            if "tendget_gpu" in name and _hadv_cont[0]:
                # continuation line(s) of the monoadvh2 call
                fixed.append(ln)
                if not ln.rstrip().endswith("&"):
                    _hadv_cont[0] = False
                    def _full2(v, tag):
                        return [
                            "      ! #region agent log",
                            f"      !$omp target update from({v})",
                            f"      call geps_dbg_ssq_full({v}, nx, levp, my_max, '{tag}')",
                            "      ! #endregion",
                        ]
                    for v in ("ttm_sl", "uum_sl", "vvm_sl", "pten_sl"):
                        fixed += _full2(v, f"{v}_post_hadv")
                    fixed.extend(_vv(244))   # keep the pre-existing marker
                continue
            if "callndslfv_monoadvh_gpu_refactor(" in low or "callndslfv_monoadvh2_gpu_refactor(" in low:
                fixed.append(ln)
                if not ln.rstrip().endswith("&"):
                    fixed.extend(_vv(244))
                else:
                    pending_244 = True
                continue
            if pending_244:
                fixed.append(ln)
                if not ln.rstrip().endswith("&"):
                    fixed.extend(_vv(244))
                    pending_244 = False
                continue
            if "callndslfv_monoadvv_gpu(" in low:
                fixed.extend(_vv(245))
                # #region agent log: same three-point bracket as src/tendget.f90
                # (geps_dbg_ssq_grid lives there and is linked into this
                # binary). Every array is mapped here: ddtemp/vdzonl/vdmerd/
                # pdot by tendget_gpu's own enter data, ttp by initial_gpu,
                # pt by the module-level enter data. Sums are used-region and
                # NaN-safe, so they compare 1:1 with the CPU's global value.
                def _grid(v, dims, tag):
                    return [
                        "      ! #region agent log",
                        f"      !$omp target update from({v})",
                        f"      call geps_dbg_ssq_grid({v}, {dims}, '{tag}')",
                        "      ! #endregion",
                    ]
                if "tendget_gpu" in name:
                    fixed += _grid("ddtemp", "nxp, lev, my_max", "ddtemp_pre_vadv")
                    fixed += _grid("vdzonl", "nxp, lev, my_max", "vdzonl_pre_vadv")
                    fixed += _grid("vdmerd", "nxp, lev, my_max", "vdmerd_pre_vadv")
                    fixed += _grid("pdot", "nxp, lev+1, latpart", "pdot_pre_vadv")
                    fixed += _grid("pt", "nxp, 1, my_max", "pt_pre_vadv")
                    fixed += _grid("ttp", "nxp, lev, my_max", "ttp_pre_vadv")
                fixed.append(ln)
                if "tendget_gpu" in name:
                    fixed += _grid("ddtemp", "nxp, lev, my_max", "ddtemp_post_vadv")
                    fixed += _grid("vdzonl", "nxp, lev, my_max", "vdzonl_post_vadv")
                    fixed += _grid("vdmerd", "nxp, lev, my_max", "vdmerd_post_vadv")
                    _pending_dt[0] = True
                # #endregion
                continue
            fixed.append(ln)
        out = fixed
    if PROBES:
        # #region agent log: size of every `target enter data` (2026-09-16).
        # The physics' first call dies with "getTargetPointer returned null" while
        # 109 GB are still free on each of 2 GPUs, i.e. one single mapping is
        # absurdly large. Before each enter data, print the byte total of its
        # items when it exceeds 1 GiB (sizeof works for arrays, sections and
        # scalars in flang). Continuation lines (`!$omp&`) are joined for the
        # clause parse; the directive itself is left untouched.
        fixed = []
        i = 0
        nline = 0
        used_dbgmem = [False]
        while i < len(out):
            ln = out[i]
            mdir = re.match(r"^\s*!\$omp\s+target\s+(enter|exit)\s+data\b", ln, re.I)
            if mdir:
                tagname = "DBGMAP" if mdir.group(1).lower() == "enter" else "DBGUNMAP"
                j = i
                text = ln
                while text.rstrip().endswith("&") and j + 1 < len(out):
                    j += 1
                    text = text.rstrip()[:-1] + re.sub(r"^\s*!\$omp&?\s*", " ", out[j])
                items: list[str] = []
                for mm in re.finditer(r"map\(\s*\w+\s*:([^()]*(?:\([^()]*\)[^()]*)*)\)", text, re.I):
                    body = mm.group(1)
                    depth = 0
                    cur = ""
                    for ch in body:
                        if ch == "(":
                            depth += 1
                        elif ch == ")":
                            depth -= 1
                        if ch == "," and depth == 0:
                            items.append(cur.strip()); cur = ""
                        else:
                            cur += ch
                    if cur.strip():
                        items.append(cur.strip())
                if items:
                    indent = re.match(r"^(\s*)", ln).group(1)
                    used_dbgmem[0] = True
                    expr = " + ".join(f"int(sizeof({v}), 8)" for v in items)
                    nline += 1
                    fixed.append(f"{indent}! #region agent log")
                    fixed.append(f"{indent}if ({expr} > 268435456_8) call geps_dbgmap('{tagname} {name}:{nline}', &")
                    fixed.append(f"{indent}   real({expr}, 8)/1073741824.0d0)")
                    fixed.append(f"{indent}! #endregion")
                fixed.extend(out[i:j + 1])
                i = j + 1
                continue
            fixed.append(ln)
            i += 1
        out = fixed
        # #endregion
        # #region agent log: free-VRAM trace for the longwave radiation path
        # (2026-09-16): 7-10 GB vanish between radlw_main's enter data #9 and
        # exit data #17 on every other call. Print free VRAM before every
        # `!$omp` construct in these files so the drop points at the construct
        # just before it. Directive continuation lines (`!$omp&`) are skipped.
    if PROBES and any(k in name for k in ("radlw_main_gpu", "grrad_gpu")):
        fixed = []
        nfree = 0
        for ln in out:
            indent = re.match(r"^(\s*)", ln).group(1)
            if re.match(r"^\s*!\$omp\s+target\b", ln, re.I):
                nfree += 1
                fixed.append(f"{indent}call geps_dbgmap('DBGFREE {name}:{nfree}', 0.0d0)")
                fixed.append(ln)
            elif re.match(r"^\s*!\$omp\s+end\s+target\b", ln, re.I):
                # AFTER the end, never inside the region (device link failed
                # on an in-region call the first time)
                fixed.append(ln)
                nfree += 1
                fixed.append(f"{indent}call geps_dbgmap('DBGFREE {name}:{nfree}', 0.0d0)")
            else:
                fixed.append(ln)
        out = fixed
    # #endregion
    return out


def _maybe_pad_rank(line: str, name: str, expect_commas: int, pad: int) -> str:
    """If `name(a, b)` has expect_commas commas at depth 1, append `, pad`."""
    out: list[str] = []
    i = 0
    pat = re.compile(rf"\b{name}\s*\(", re.IGNORECASE)
    while True:
        m = pat.search(line[i:])
        if not m:
            out.append(line[i:])
            break
        start = i + m.start()
        p0 = i + m.end() - 1
        depth = 0
        k = p0
        commas = 0
        while k < len(line):
            if line[k] == "(":
                depth += 1
            elif line[k] == ")":
                depth -= 1
                if depth == 0:
                    break
            elif line[k] == "," and depth == 1:
                commas += 1
            k += 1
        inner = line[p0 + 1 : k]
        out.append(line[i:start])
        if commas == expect_commas:
            out.append(f"{name}({inner}, {pad})")
        else:
            out.append(line[start : k + 1])
        i = k + 1
    return "".join(out)


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: acc2omp.py IN.f90 OUT.f90", file=sys.stderr)
        return 2
    translate_file(Path(argv[1]), Path(argv[2]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
