#!/bin/bash

disp(){
 echo "$0 [fx1000|a100]" ; exit 1
}

if [[ "$1" =~ \-h|help ]] ; then
  disp 
fi

set -x

MDIR=$(pwd)
echo "MDIR=${MDIR}"
cd ${MDIR}

machine=${1:-fx1000}

# submit with sum group
group=sum

if [ $(id -gn) != $group ]; then
  exec sg $group "$0 $*"
fi

JID=$(pjsub -z jid TCo383L72_IC_sample_${machine} -x CMAKE_BUILD=1,GITLAB_CICD=1,machine=$machine -g ${group} )

# wait job finish
echo "Job ID :  $JID"
OUT=$(/usr/bin/pjwait ${JID})
echo "pjwait : ${OUT}"

PJM_CODE=$(echo $OUT | cut -d' ' -f2 )
EXIT_CODE=$(echo $OUT | cut -d' ' -f3 )
SIGNAL=$(echo $OUT | cut -d' ' -f4 )

if [ "$PJM_CODE" -eq 0 ] && [ "$EXIT_CODE" -eq 0 ] && [ "$SIGNAL" -eq 0 ]; then
  # Continue with the rest of your script if all variables are zero
  echo "All variables are zero. Continue with the rest of the script." && exit 0
else
  echo "One of the variables is not zero. Exiting with code 9." && exit 9
fi

