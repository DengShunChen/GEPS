# Layer 18 (2026-09-16): tt/qt come out of pbl_noah_gpu with ~10% non-finite
# points on ROCm.  GPU-only checkpoints after each stage inside pbl_noah_gpu;
# the signal is the `nonfinite=` count, so no CPU counterpart is needed.
G = lambda v, d, t: ('grid', v, d, t)
COL = lambda pre: [G('t1','nx, lev, my_max',pre+'_t1'), G('q1','nx, lev*ntrac, my_max',pre+'_q1'),
                   G('u1','nx, lev, my_max',pre+'_u1'), G('v1','nx, lev, my_max',pre+'_v1')]
SFC = lambda pre: [G('heat','nx, 1, my_max',pre+'_heat'), G('evap','nx, 1, my_max',pre+'_evap'),
                   G('stress','nx, 1, my_max',pre+'_stress'), G('ustar','nx, 1, my_max',pre+'_ustar'),
                   G('tg','nx, 1, my_max',pre+'_tg'), G('cd','nx, 1, my_max',pre+'_cd'), G('cdq','nx, 1, my_max',pre+'_cdq')]
CK = [
 (None, "call sfc_diff_gpu(", 1, 'before', COL('b_in') + SFC('b_in')),
 (None, "call sfc_diff_gpu(", 1, 'after', SFC('b_diff')),
 (None, "call sfc_ocean_gpu(", 1, 'after', SFC('b_ocean')),
 (None, "call sfc_drv_gpu(", 1, 'after', SFC('b_drv')),
 (None, "call sfc_sice_gpu(", 1, 'after', SFC('b_sice')),
 (None, "call sfc_diag_gpu(", 1, 'after', SFC('b_diag') + [G('hpbl','nx, 1, my_max','b_diag_hpbl')]),
 (None, "call moninedmf_gpu(", 1, 'after', COL('b_edmf') + [G('hpbl','nx, 1, my_max','b_edmf_hpbl')]),
]
