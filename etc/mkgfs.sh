#!/bin/sh

hh=0 
htot=$1
outpath=$2
intvh=$3
file=gfsctl_${htot}
cd $outpath

if [ -f $file ]; then
  cp $file ${file}.bak
  rm -rf ${file}
fi

while [ $hh -le $htot ]
do

if [ ${hh} -lt 10 ]; then
   hht=00${hh}
elif [ ${hh} -lt 100 ]; then
   hht=0${hh}
else
   hht=${hh}
fi
  echo "${hht}" >> $file



 let hh=${hh}+${intvh}
done

#hh=6
#while [ $hh -le 6312 ]
#do

#if [ ${hh} -lt 10 ]; then
#   hht=00${hh}
#elif [ ${hh} -lt 100 ]; then
#   hht=0${hh}
#else
#   hht=${hh}
#fi
#  echo "W00100      ${hht}" >> $file
#  echo "s00320      ${hht}" >> $file
#  echo "s004d0      ${hht}" >> $file
#  echo "s004f0      ${hht}" >> $file
#  echo "S00100      ${hht}" >> $file
# let hh=${hh}+6
#done
#  echo "nomodata         "  >> $file
