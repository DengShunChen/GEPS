# Checkpoints along integration step 1. Each entry:
#   (cpu_anchor, gpu_anchor, nth, where, probes)
#   where = 'after' (after the call statement ends) | 'before'
#   probes = list of (kind, var, dims, tag); kind in grid|full|spec|spec2
G = lambda v, d, t: ('grid', v, d, t)
S = lambda v, d, t: ('spec', v, d, t)
S2 = lambda v, d, t: ('spec2', v, d, t)
SC = lambda v, t: ('scal', v, '', t)
CK = [
 # right after the n-level copies (CPU: the do-loop that ends with `rdivm = rdiv`
 # -- anchored on the statement AFTER it; GPU: the memcpy block, anchored on
 # its last memcpy). 'after' on a non-call line works the same way.
 ("forward = .true.", "itt = min(itimestep, 2)", 1, 'before',
   [G('up','nxp, lev, my_max','s0_copy_up'), G('vp','nxp, lev, my_max','s0_copy_vp'), G('ttp','nxp, lev, my_max','s0_copy_ttp'),
    G('ut','nxp, lev, my_max','s0_copy_ut'), G('vt','nxp, lev, my_max','s0_copy_vt'), G('um','nxp, lev, my_max','s0_copy_um'),
    G('rdivm','nxp, lev, my_max','s0_copy_rdivm'), G('qm','nxp, lev*ncld, my_max','s0_copy_qm'), G('ptm','nxp, 1, my_max','s0_copy_ptm')]),
 ("call joinrs(cc, ddtemp,", "call joinrs_gpu(cc_cg, ddtemp,", 2, 'before',
   [G('vdzonl','nxp, lev, my_max','s2_tendin_vdzonl'), G('vdmerd','nxp, lev, my_max','s2_tendin_vdmerd'), G('ddtemp','nxp, lev, my_max','s2_tendin_ddtemp'), SC('pcorr','s2_pcorr'),
    G('up','nxp, lev, my_max','s2_tendin_up'), G('vp','nxp, lev, my_max','s2_tendin_vp'), G('ttp','nxp, lev, my_max','s2_tendin_ttp'), G('ut','nxp, lev, my_max','s2_tendin_ut'), G('vt','nxp, lev, my_max','s2_tendin_vt')]),
 ("call tranrs(", "call tranrs_gpu_cuda_graph(", 4, 'after',
   [S('temten','levp*2, jtrun, jtmax','s2_tranrs_temten')]),
 ("call siimpl(", "call siimpl_gpu(", 2, 'before',
   [S('temnow','levp*2, jtrun, jtmax','s2_siimplin_temnow'), S('divnow','levp*2, jtrun, jtmax','s2_siimplin_divnow'), S2('plnow','jtrun, jtmax, 2','s2_siimplin_plnow'),
    S2('plten','jtrun, jtmax, 2','s2_siimplin_plten'), S('temten','levp*2, jtrun, jtmax','s2_siimplin_temten'), S('divten','levp*2, jtrun, jtmax','s2_siimplin_divten')]),
 ("call mpe2d_transpose_ndsl_f2p(ttm_sl, ddtemp,", "call mpe2d_transpose_ndsl_f2p_gpu(ttm_sl, ddtemp,", 1, 'after',
   [G('vdzonl','nxp, lev, my_max','s1_fgnl_vdzonl'), G('vdmerd','nxp, lev, my_max','s1_fgnl_vdmerd'), G('ddtemp','nxp, lev, my_max','s1_fgnl_ddtemp')]),
 ("call joinrs(cc, diveng,", "call joinrs_gpu(cc_cg, diveng,", 1, 'before',
   [G('diveng','nxp, lev, my_max','s1_gridnl_diveng'), G('pdot','nxp, lev+1, latpart','s1_gridnl_pdot')]),
 ("call trngra3(", "call trngra3_gpu_cuda_graph(", 1, 'after',
   [S('hldten','levp*2, jtrun, jtmax','s1_trngra3_hldten'), G('dlphi','nxp, lev, my_max','s1_trngra3_dlphi'), G('dtphi','nxp, lev, my_max','s1_trngra3_dtphi')]),
 ("call ndslfv_monoadvv_fgnl(", "call ndslfv_monoadvv_fgnl_gpu(", 1, 'after',
   [G('vdzonl','nxp, lev, my_max','s1_vadv_vdzonl'), G('vdmerd','nxp, lev, my_max','s1_vadv_vdmerd'), G('ddtemp','nxp, lev, my_max','s1_vadv_ddtemp')]),
 ("call tranrs1(", "call tranrs1_gpu(", 1, 'after',
   [S2('plten','jtrun, jtmax, 2','s1_tranrs1_plten')]),
 ("call rstrandz(", "call trandv_gpu_cuda_graph(", 1, 'after',
   [S('temten','levp*2, jtrun, jtmax','s1_trandv_temten'), S('vorten','levp*2, jtrun, jtmax','s1_trandv_vorten'), S('divten','levp*2, jtrun, jtmax','s1_trandv_divten')]),
 ("call siimpl(", "call siimpl_gpu(", 1, 'after',
   [S('temmid','levp*2, jtrun, jtmax','s1_siimpl_temmid'), S('divmid','levp*2, jtrun, jtmax','s1_siimpl_divmid'), S2('plmid','jtrun, jtmax, 2','s1_siimpl_plmid')]),
 ("call hdiffu(", "call hdiffu_gpu(", 1, 'after',
   [S('vormid','levp*2, jtrun, jtmax','s1_hdiffu_vormid'), S('divmid','levp*2, jtrun, jtmax','s1_hdiffu_divmid'), S('temmid','levp*2, jtrun, jtmax','s1_hdiffu_temmid')]),
 ("call trngra(", "call trngra_gpu(", 1, 'after',
   [G('um','nxp, lev, my_max','s1_mid_um'), G('vm','nxp, lev, my_max','s1_mid_vm'), G('rdivm','nxp, lev, my_max','s1_mid_rdivm'), G('tm','nxp, lev, my_max','s1_mid_tm'),
    G('ptm','nxp, 1, my_max','s1_mid_ptm'), G('dlpl','nxp, 1, my_max','s1_mid_dlpl'), G('dtpl','nxp, 1, my_max','s1_mid_dtpl')]),
 ("call joinrs(cc, diveng,", "call joinrs_gpu(cc_cg, diveng,", 2, 'before',
   [G('diveng','nxp, lev, my_max','s2_gridnl_diveng'), G('pdot','nxp, lev+1, latpart','s2_gridnl_pdot')]),
 ("call trngra3(", "call trngra3_gpu_cuda_graph(", 2, 'after',
   [S('hldten','levp*2, jtrun, jtmax','s2_trngra3_hldten'), G('dlphi','nxp, lev, my_max','s2_trngra3_dlphi'), G('dtphi','nxp, lev, my_max','s2_trngra3_dtphi')]),
 ("call mpe2d_transpose_ndsl_f2p(qm_sl, qt,", "call mpe2d_transpose_ndsl_f2p_gpu(qm_sl, qt,", 1, 'after',
   [G('tt','nxp, lev, my_max','s2_hadv_tt'), G('ut','nxp, lev, my_max','s2_hadv_ut'), G('vt','nxp, lev, my_max','s2_hadv_vt'), G('qt','nxp, lev*ncld, my_max','s2_hadv_qt')]),
 ("call ndslfv_monoadvv(tt,", "call ndslfv_monoadvv_gpu(tt,", 1, 'after',
   [G('tt','nxp, lev, my_max','s2_vadv_tt'), G('ut','nxp, lev, my_max','s2_vadv_ut'), G('vt','nxp, lev, my_max','s2_vadv_vt'), G('qt','nxp, lev*ncld, my_max','s2_vadv_qt')]),
 ("call transr1(jtrun, jtmax, nx, my, my_max, poly, pltemp, pt,", "call transr1_gpu(jtrun, jtmax, nx, my, my_max, polyf, pltemp, pt,", 1, 'after',
   [G('pt','nxp, 1, my_max','s2_pt')]),
 ("call rstrandz(", "call trandv_gpu_cuda_graph(", 2, 'after',
   [S('temten','levp*2, jtrun, jtmax','s2_trandv_temten'), S('vorten','levp*2, jtrun, jtmax','s2_trandv_vorten'), S('divten','levp*2, jtrun, jtmax','s2_trandv_divten')]),
 ("call siimpl(", "call siimpl_gpu(", 2, 'after',
   [S('temten','levp*2, jtrun, jtmax','s2_siimpl_temten'), S('divten','levp*2, jtrun, jtmax','s2_siimpl_divten'), S2('plten','jtrun, jtmax, 2','s2_siimpl_plten')]),
 ("call transr1(jtrun, jtmax, nx, my, my_max, poly, plten, ptend,", "call transr1_gpu(jtrun, jtmax, nx, my, my_max, polyf, plten, ptend,", 1, 'after',
   [G('ptend','nxp, 1, my_max','s2_ptend'), S2('plten','jtrun, jtmax, 2','s2_ptend_plten')]),
]
def probe_lines(probes, gpu):
    out = ["      ! #region agent log"]
    for kind, v, d, t in probes:
        if kind == 'scal':
            # host scalar; same record format so cmp_ckpt.py can read it
            out.append(f"      if (myrank .eq. 0) print *,'DBGGRID {t} ssq=',real({v},kind=8)**2,' nonfinite= 0'")
            continue
        if gpu:
            out.append(f"      !$omp target update from({v})")
        fn = {'grid':'geps_dbg_ssq_grid','full':'geps_dbg_ssq_full','spec':'geps_dbg_ssq_spec','spec2':'geps_dbg_ssq_spec2','cnt':'geps_dbg_cnt_grid'}[kind]
        out.append(f"      call {fn}({v}, {d}, '{t}')")
    out.append("      ! #endregion")
    return out
def _cont(ln):
    """True if this source line continues (trailing '&', ignoring a trailing comment)."""
    code, q = [], None
    for ch in ln:
        if q:
            if ch == q: q = None
        elif ch in ("'", '"'): q = ch
        elif ch == '!': break
        code.append(ch)
    return ''.join(code).rstrip().endswith('&')
def inject(lines, anchors_key, gpu, ck=None):
    """anchors_key: index 0 (cpu) or 1 (gpu) into CK entries. Continuation-aware."""
    ck = CK if ck is None else ck
    counts = {}
    out = []
    i = 0
    while i < len(lines):
        ln = lines[i]
        st = ln.strip()
        hits = []
        seen_anchor = set()
        for ent in ck:
            a = ent[anchors_key]
            if a is not None and st.startswith(a):
                if a not in seen_anchor:          # count once per line
                    seen_anchor.add(a)
                    counts[a] = counts.get(a, 0) + 1
                if counts[a] == ent[2]:
                    hits.append(ent)
        if not hits:
            out.append(ln); i += 1; continue
        for hit in hits:
            if hit[3] == 'before':
                out += probe_lines(hit[4], gpu)
        out.append(ln)
        while _cont(lines[i]):
            i += 1
            out.append(lines[i])
            # comment / cpp lines inside a continued statement don't end it
            while (lines[i].lstrip().startswith('!') or lines[i].startswith('#')) \
                    and not _cont(lines[i]):
                i += 1
                out.append(lines[i])
        # statement tail split by `#ifdef X ... #else ... #endif`: we are on the
        # last line of the first branch; swallow the alternative branch too.
        if i + 1 < len(lines) and lines[i + 1].startswith('#else'):
            while not lines[i].startswith('#endif'):
                i += 1
                out.append(lines[i])
        for hit in hits:
            if hit[3] == 'after':
                out += probe_lines(hit[4], gpu)
        i += 1
    return out, counts
