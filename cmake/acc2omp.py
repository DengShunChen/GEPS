#!/usr/bin/env python3
"""Translate NVIDIA OpenACC + CUDA Fortran extensions to OpenMP offload.

Used only for GEPS_COMPILER=rocm + USE_HIP. Original sources are not modified;
CMake writes translated copies into the build tree.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

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


def _rewrite_clauses(body: str) -> str:
    # Longer names first: present_or_copyin contains copyin, present_or_copy contains copy.
    for pat, repl in (
        (r"\bpresent_or_copyin\s*\(", "map(to:"),
        (r"\bpresent_or_copyout\s*\(", "map(from:"),
        (r"\bpresent_or_copy\s*\(", "map(tofrom:"),
        (r"\bpresent_or_create\s*\(", "map(alloc:"),
        (r"\bcopyin\s*\(", "map(to:"),
        (r"\bcopyout\s*\(", "map(from:"),
        (r"\bcopy\s*\(", "map(tofrom:"),
        (r"\bcreate\s*\(", "map(alloc:"),
        (r"\bdelete\s*\(", "map(delete:"),
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
        body = _rewrite_clauses(rest[len("enter data") :])
        return emit(f"target enter data {body}".strip())
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
        is_seq = bool(re.search(r"\bseq\b", rest[len("loop") :], re.IGNORECASE))
        body = _rewrite_clauses(rest[len("loop") :])
        extra = f" {body}" if body else ""
        if is_seq:
            return f"{indent}! acc2omp: loop seq (host-order on device)"
        return emit(f"parallel do{extra}")

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
    joined = _merge_parallel_loop(_join_continuations(lines))
    out: list[str] = [
        "! acc2omp: generated from " + src.as_posix(),
        "! Do not edit; regenerate via cmake/acc2omp.py",
    ]
    for line, is_acc in joined:
        if is_acc:
            out.append(translate_directive(line))
        else:
            out.append(sanitize_cuda_fortran(line))
    out = drop_duplicate_decls(out)
    out = hide_use_collisions(out)
    out = inject_fpvs_external(out)
    out = apply_file_fixups(src.name.lower(), out)
    out = relax_intent_in_on_assignment(out)
    out = strip_parameter_maps(out)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(out) + "\n")


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
    joined = _join_omp_continuations(lines)
    params = _collect_parameters(joined)
    out: list[str] = []
    for line in joined:
        if OMP_SENT.match(line):
            line = _rewrite_map_drop_params(line, params)
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


def apply_file_fixups(name: str, lines: list[str]) -> list[str]:
    """HIP/flang-only source tweaks for NVIDIA dialect that nvfortran accepts."""
    out = lines
    if "mpe2d_gpu" in name:
        out = [
            re.sub(
                r"\bain\(([^(),]+),([^(),]+),([^(),]+),([^(),]+),([^(),]+)\)",
                r"ain(\1,\2,\3,\4,\5, 1)",
                ln,
            )
            for ln in out
        ]
    if "moninedmf_gpu" in name:
        fixed = []
        for ln in out:
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
    if "zx_gpu" in name:
        # Sequence association: pointer element cannot be dummy array; use section.
        out = [re.sub(r"\bvars\(idx\)\s*,", "vars(idx:),", ln) for ln in out]
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
