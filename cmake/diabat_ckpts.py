# Checkpoints through the physics driver (diabat / diabat_gpu), same format as
# intgrt_ckpts.CK.  CPU side is one big column loop (do 290 jj), so the CPU
# only gets p0 (before the loop) and p1 (after it); the GPU side is a chain of
# whole-domain kernels, so it also gets g_* after each physics call.  Purpose
# (layer 18, 2026-09-16): find which physics kernel leaves garbage in tt/qt
# before the microphysics (fall_flux_gpu spins forever on rho/tz garbage).
G = lambda v, d, t: ('grid', v, d, t)
STATE = lambda pre: [G('tt','nxp, lev, my_max',pre+'_tt'), G('qt','nxp, lev*ncld, my_max',pre+'_qt'),
                     G('ut','nxp, lev, my_max',pre+'_ut'), G('vt','nxp, lev, my_max',pre+'_vt')]
CK = [
 ("do 290 jj =1, jlistnum", "call rrtmg_gpu", 1, 'before', STATE('p0')),
 ("290 continue", "call mp_scheme_gpu", 1, 'after', STATE('p1')),
 (None, "call rrtmg_gpu", 1, 'after', STATE('g_rad')),
 (None, "call dcyc2t3_gpu", 1, 'after', STATE('g_dcyc')),
 (None, "call pbl_noah_gpu(", 1, 'after', STATE('g_pbl')),
 (None, "call gwdps_gpu(", 1, 'after', STATE('g_gwdps')),
 (None, "call samfdeepcnv_kh_gpu(", 1, 'after', STATE('g_deep') + [G('qtc','nxp, lev, my_max','g_deep_qtc'), G('ttc','nxp, lev, my_max','g_deep_ttc'),
   G('rcup','nxp, 1, my_max','g_deep_rcup'), G('cldwrk','nxp, 1, my_max','g_deep_cldwrk'),
   ('cnt','kuo','nxp, my_max','g_deep_kuo'), ('cnt','kbot','nxp, my_max','g_deep_kbot'), ('cnt','ktop','nxp, my_max','g_deep_ktop')]),
 (None, "call lightning_ec_gpu(", 1, 'after', STATE('g_light')),
 (None, "call gwdc_gpu(", 1, 'after', STATE('g_gwdc')),
 (None, "call samfshalcnv_kh_gpu(", 1, 'after', STATE('g_shal') + [G('qtc','nxp, lev, my_max','g_shal_qtc'), G('ttc','nxp, lev, my_max','g_shal_ttc')]),
 (None, "call mp_scheme_gpu", 1, 'before', STATE('g_premp')),
 # end of physics on both sides: CPU after the big loop (adjptqintp is inside
 # it), GPU after adjptqintp_gpu.  p1 is NOT comparable (GPU p1 = after
 # mp_scheme only); p2 is.
 ("290 continue", "call adjptqintp_gpu(", 1, 'after', STATE('p2')),
]
CK_CPU = [e for e in CK if e[0] is not None]
