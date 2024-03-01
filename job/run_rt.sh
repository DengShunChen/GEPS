#!/bin/bash -x

JID=$(pjsub -z jid TCo383L72_IC_sample_fx1000 -x CMAKE_BUILD=1,GITLAB_CICD=1 -g sum)
EXIT_CODE=$(/usr/bin/pjwait ${JID} | cut -d' ' -f3 )
exit ${EXIT_CODE}
