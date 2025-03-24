#!/bin/sh

 tautemp=0
 endtau=168
 updatehr=6
 while [ ${tautemp} -le ${endtau} ]
 do
   tau=$(printf "%04d" ${tautemp})

   echo "${tau}" >> gfsctl.${endtau}.${updatehr}
   
   tautemp=$((tautemp+updatehr))
   
  done
