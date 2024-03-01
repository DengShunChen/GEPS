#!/bin/bash -x

JID=$(pjsub -z jid TCo383L72_IC_sample_fx1000 -x CMAKE_BUILD=1,GITLAB_CICD=1 -g sum)
OUT=$(/usr/bin/pjwait ${JID})

PJM_CODE=$(echo $OUT | cut -d' ' -f2 )
EXIT_CODE=$(echo $OUT | cut -d' ' -f3 )
SIGNAL=$(echo $OUT | cut -d' ' -f4 )

if [ "$PJM_CODE" -ne 0 ] || [ "$EXIT_CODE" -ne 0 ] || [ "$SIGNAL" -ne 0 ]; then
  echo "One of the variables is not zero. Exiting with code 9." && exit 9
else
  # Continue with the rest of your script if all variables are zero
  echo "All variables are zero. Continue with the rest of the script." && exit 0
fi

