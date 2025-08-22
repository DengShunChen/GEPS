#!/bin/sh

 tautemp=24
 endtau=168
 updatehr=24
 while [ ${tautemp} -le ${endtau} ]
 do
   tau=$(printf "%08d" ${tautemp})

   echo "${tau}" >> timectl.${endtau}.${updatehr}
   
   tautemp=$((tautemp+updatehr))
   
  done
