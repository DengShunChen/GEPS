#!/bin/bash
# Bind to GPU0 only.
#
# 2026-09-10: this lab now has TWO PHYSICAL MI300X cards, both NPS1/SPX (unpartitioned):
#   GPU[0] = KFD node 4, PCI 0000:46:00.0
#   GPU[1] = KFD node 5, PCI 0000:66:00.0
# GPU[1] is a REAL SECOND CARD, not a DPX partition (that was the pre-09-10 layout,
# where GPU[1] was the c5:00.1 partition of a single card and had to be avoided).
# We still pin to GPU0 because the model runs 1 MPI rank; multi-GPU needs RCCL work.
#
# GEPS_ROCM_PCI is informational only - nothing reads it. Kept for log breadcrumbs.
export GEPS_ROCM_PCI="${GEPS_ROCM_PCI:-0000:46:00.0}"
export HIP_VISIBLE_DEVICES=0
export ROCR_VISIBLE_DEVICES=0
export GPU_DEVICE_ORDINAL=0
export OMP_DEFAULT_DEVICE=0
exec "$@"
