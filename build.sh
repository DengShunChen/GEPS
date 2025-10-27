#!/bin/bash
#------------------------------------------------------------------------------
#
#  build GFS forecast model
#
# Usage:
#   ./build.sh [MACHINE] [OPTION]
#     MACHINE: fx10, fx100, pcc, fx1000
#     OPTION : 2cpl_ice, 2cpl, 4cpl(optional, for coupling build)
#
#------------------------------------------------------------------------------
#
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
  echo "Usage: $0 [MACHINE] [OPTION: 2cpl 4cpl 2cpl_CICE]"
  exit 1
fi

export MACHINE=$1
OPTION=$2

machines='fx1000 fx100 fx10 pcc'
[[ $machines =~ (^|[[:space:]])$MACHINE($|[[:space:]]) ]] && known='True' || known='False'
if [ "${known}" == 'False' ]; then
  echo "Fatal Error: Unknown machine --> ${MACHINE}"
  exit 1
fi

if [ "${MACHINE}" == 'fx1000' ]; then
  if ! [[ $HOSTNAME =~ h6ln[0-9][0-9] ]]; then
    echo "Fatal Error: Please build on login node of HPC Gen6 (h6ln??)"
    exit 1
  fi
fi

set -x

# load libs
export MDIR=$(pwd)
. /usr/share/Modules/init/bash
module purge
module use  ${MDIR}/modulefiles
module av

# Decide module to load
case "$OPTION" in
  2cpl)
    module load modulefile.tcogfs.${MACHINE}_2cpl
    export TIMCOMCPL=TRUE
    export clpath=/data/common/gfs/GEPSv3_lib/coupler/fx1000
    ;;
  4cpl)
    module load modulefile.tcogfs.${MACHINE}_4cpl
    export TIMCOMCPL=TRUE
    export clpath=/data/common/gfs/GEPSv3_lib/coupler/fx1000
    ;;
  "")
    module load modulefile.tcogfs.${MACHINE}
    export TIMCOMCPL=FALSE
    ;;
  *)
    echo "Fatal Error: Unknown option --> $OPTION"
    echo "Usage: $0 [MACHINE] [OPTION: 2cpl_CICE | 2cpl | 4cpl]"
    exit 1
    ;;
esac

module list
module unuse ${MDIR}/modulefiles

# compile
cd src/
make clean
make -j24 TIMCOMCPL=$TIMCOMCPL

if [[ $? -ne 0 ]]; then
  echo "Build failed for MACHINE=$MACHINE OPTION=$OPTION"
  exit 1
fi

# Echo current build info
echo "✅ Build completed for:"
echo "   MACHINE : $MACHINE"
if [[ "$OPTION" == "2cpl_CICE" || "$OPTION" == "4cpl" || "$OPTION" == "2cpl" ]]
then
  echo "    MODE    : Coupled (${OPTION})"
else
  echo "    MODE    : Atmosphere-only"
fi
