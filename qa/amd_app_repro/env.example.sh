# Copy to env.sh and edit. Every script in this directory sources ./env.sh.
#
# Layout under $WORK (the three directories must be siblings; the job script finds
# GEPS_LIB as ../GEPS_LIB relative to GEPS/ unless GEPS_LIB_ROOT is set):
#   $WORK/GEPS       model source (git clone, see README §1)
#   $WORK/GEPS_LIB   dependency libraries (git clone, see README §2)
#   $WORK/dmsdb      input data (tarball, see README §1); must be writable, the run adds TCo383L72.ufs/
export WORK=${WORK:-$HOME/geps_repro}

export ROCM_PATH=${ROCM_PATH:-/opt/afar/current}     # afar 24.3.0 SDK (01_install_afar.sh)
export GEPS_SRC=$WORK/GEPS
export GEPS_LIB_ROOT=$WORK/GEPS_LIB
export GEPS_DMSDB=$WORK/dmsdb
export GEPS_ROCM_BUILD_DIR=build_afar_hip            # build tree inside $GEPS_SRC
export GEPS_ROCM_DEVICES=${GEPS_ROCM_DEVICES:-0,1}   # two MI300X, one MPI rank each
export NJ=${NJ:-32}                                   # parallel make jobs
