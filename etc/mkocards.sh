#!/bin/sh

hh=0
htot=$1
outpath=$2
cd ${outpath}
file=ocards_${htot}

if [ -f $file ]; then
  cp $file ${file}.bak
  rm -rf ${file}
fi

while [ $hh -le ${htot} ]
do

if [ ${hh} -lt 10 ]; then
   hht=00${hh}
elif [ ${hh} -lt 100 ]; then
   hht=0${hh}
else
   hht=${hh}
fi
  echo "SSL010      ${hht}" >> $file
  echo "B00010      ${hht}" >> $file
  echo "B00100      ${hht}" >> $file
  echo "B00200      ${hht}" >> $file
  echo "B00210      ${hht}" >> $file
  echo "S00300      ${hht}" >> $file
  echo "B00510      ${hht}" >> $file
  echo "B10510      ${hht}" >> $file
  echo "B00620      ${hht}" >> $file
  echo "B00626      ${hht}" >> $file
  echo "B00630      ${hht}" >> $file
  echo "B00640      ${hht}" >> $file
  echo "B00650      ${hht}" >> $file
  echo "B00651      ${hht}" >> $file
  echo "B006A0      ${hht}" >> $file
  echo "B02100      ${hht}" >> $file
  echo "B02510      ${hht}" >> $file
  echo "B10200      ${hht}" >> $file
  echo "B10210      ${hht}" >> $file
  echo "S00420      ${hht}" >> $file
  echo "S00430      ${hht}" >> $file
  echo "S00440      ${hht}" >> $file
  echo "S00450      ${hht}" >> $file
  echo "S00460      ${hht}" >> $file
  echo "S00040      ${hht}" >> $file
  echo "X00330      ${hht}" >> $file
  echo "X00340      ${hht}" >> $file
  echo "X00590      ${hht}" >> $file
  echo "S00030      ${hht}" >> $file
  echo "h2o       1 ${hht}" >> $file
  echo "h2o       2 ${hht}" >> $file
  echo "h2o       3 ${hht}" >> $file
  echo "h2o       5 ${hht}" >> $file
  echo "h2o       7 ${hht}" >> $file
  echo "h2o      10 ${hht}" >> $file
  echo "h2o      20 ${hht}" >> $file
  echo "h2o      30 ${hht}" >> $file
  echo "h2o      50 ${hht}" >> $file
  echo "h2o      70 ${hht}" >> $file
  echo "h2o     100 ${hht}" >> $file
  echo "h2o     150 ${hht}" >> $file
  echo "h2o     200 ${hht}" >> $file
  echo "h2o     250 ${hht}" >> $file
  echo "h2o     300 ${hht}" >> $file
  echo "h2o     350 ${hht}" >> $file
  echo "h2o     400 ${hht}" >> $file
  echo "h2o     450 ${hht}" >> $file
  echo "h2o     500 ${hht}" >> $file
  echo "h2o     550 ${hht}" >> $file
  echo "h2o     600 ${hht}" >> $file
  echo "h2o     650 ${hht}" >> $file
  echo "h2o     700 ${hht}" >> $file
  echo "h2o     750 ${hht}" >> $file
  echo "h2o     800 ${hht}" >> $file
  echo "h2o     850 ${hht}" >> $file
  echo "h2o     900 ${hht}" >> $file
  echo "h2o     925 ${hht}" >> $file
  echo "h2o     950 ${hht}" >> $file
  echo "h2o     975 ${hht}" >> $file
  echo "h2o    1000 ${hht}" >> $file
  echo "wnd       1 ${hht}" >> $file
  echo "wnd       2 ${hht}" >> $file
  echo "wnd       3 ${hht}" >> $file
  echo "wnd       5 ${hht}" >> $file
  echo "wnd       7 ${hht}" >> $file
  echo "wnd      10 ${hht}" >> $file
  echo "wnd      20 ${hht}" >> $file
  echo "wnd      30 ${hht}" >> $file
  echo "wnd      50 ${hht}" >> $file
  echo "wnd      70 ${hht}" >> $file
  echo "wnd     100 ${hht}" >> $file
  echo "wnd     150 ${hht}" >> $file
  echo "wnd     200 ${hht}" >> $file
  echo "wnd     250 ${hht}" >> $file
  echo "wnd     300 ${hht}" >> $file
  echo "wnd     350 ${hht}" >> $file
  echo "wnd     400 ${hht}" >> $file
  echo "wnd     450 ${hht}" >> $file
  echo "wnd     500 ${hht}" >> $file
  echo "wnd     550 ${hht}" >> $file
  echo "wnd     600 ${hht}" >> $file
  echo "wnd     650 ${hht}" >> $file
  echo "wnd     700 ${hht}" >> $file
  echo "wnd     750 ${hht}" >> $file
  echo "wnd     800 ${hht}" >> $file
  echo "wnd     850 ${hht}" >> $file
  echo "wnd     900 ${hht}" >> $file
  echo "wnd     925 ${hht}" >> $file
  echo "wnd     950 ${hht}" >> $file
  echo "wnd     975 ${hht}" >> $file
  echo "wnd    1000 ${hht}" >> $file
  echo "tmp       1 ${hht}" >> $file
  echo "tmp       2 ${hht}" >> $file
  echo "tmp       3 ${hht}" >> $file
  echo "tmp       5 ${hht}" >> $file
  echo "tmp       7 ${hht}" >> $file
  echo "tmp      10 ${hht}" >> $file
  echo "tmp      20 ${hht}" >> $file
  echo "tmp      30 ${hht}" >> $file
  echo "tmp      50 ${hht}" >> $file
  echo "tmp      70 ${hht}" >> $file
  echo "tmp     100 ${hht}" >> $file
  echo "tmp     150 ${hht}" >> $file
  echo "tmp     200 ${hht}" >> $file
  echo "tmp     250 ${hht}" >> $file
  echo "tmp     300 ${hht}" >> $file
  echo "tmp     350 ${hht}" >> $file
  echo "tmp     400 ${hht}" >> $file
  echo "tmp     450 ${hht}" >> $file
  echo "tmp     500 ${hht}" >> $file
  echo "tmp     550 ${hht}" >> $file
  echo "tmp     600 ${hht}" >> $file
  echo "tmp     650 ${hht}" >> $file
  echo "tmp     700 ${hht}" >> $file
  echo "tmp     750 ${hht}" >> $file
  echo "tmp     800 ${hht}" >> $file
  echo "tmp     850 ${hht}" >> $file
  echo "tmp     900 ${hht}" >> $file
  echo "tmp     925 ${hht}" >> $file
  echo "tmp     950 ${hht}" >> $file
  echo "tmp     975 ${hht}" >> $file
  echo "tmp    1000 ${hht}" >> $file
  echo "phi       1 ${hht}" >> $file
  echo "phi       2 ${hht}" >> $file
  echo "phi       3 ${hht}" >> $file
  echo "phi       5 ${hht}" >> $file
  echo "phi       7 ${hht}" >> $file
  echo "phi      10 ${hht}" >> $file
  echo "phi      20 ${hht}" >> $file
  echo "phi      30 ${hht}" >> $file
  echo "phi      50 ${hht}" >> $file
  echo "phi      70 ${hht}" >> $file
  echo "phi     100 ${hht}" >> $file
  echo "phi     150 ${hht}" >> $file
  echo "phi     200 ${hht}" >> $file
  echo "phi     250 ${hht}" >> $file
  echo "phi     300 ${hht}" >> $file
  echo "phi     350 ${hht}" >> $file
  echo "phi     400 ${hht}" >> $file
  echo "phi     450 ${hht}" >> $file
  echo "phi     500 ${hht}" >> $file
  echo "phi     550 ${hht}" >> $file
  echo "phi     600 ${hht}" >> $file
  echo "phi     650 ${hht}" >> $file
  echo "phi     700 ${hht}" >> $file
  echo "phi     750 ${hht}" >> $file
  echo "phi     800 ${hht}" >> $file
  echo "phi     850 ${hht}" >> $file
  echo "phi     900 ${hht}" >> $file
  echo "phi     925 ${hht}" >> $file
  echo "phi     950 ${hht}" >> $file
  echo "phi     975 ${hht}" >> $file
  echo "phi    1000 ${hht}" >> $file
  echo "x00730      ${hht}" >> $file
  echo "x00740      ${hht}" >> $file
  echo "x00750      ${hht}" >> $file
  echo "x00760      ${hht}" >> $file
  echo "x00770      ${hht}" >> $file
  echo "s00310      ${hht}" >> $file
  echo "s00320      ${hht}" >> $file
  echo "S004D0      ${hht}" >> $file
  echo "S004E0      ${hht}" >> $file
  echo "S00100      ${hht}" >> $file
  echo "M01200      ${hht}" >> $file
  echo "M02200      ${hht}" >> $file
  echo "M01210      ${hht}" >> $file
  echo "M02210      ${hht}" >> $file
  echo "W00100      ${hht}" >> $file
  echo "S005A0      ${hht}" >> $file
  echo "S005A1      ${hht}" >> $file
  echo "S003X0      ${hht}" >> $file
  echo "S003U0      ${hht}" >> $file
  echo "S005C0      ${hht}" >> $file
  echo "sa15b0      ${hht}" >> $file
  echo "sa25b0      ${hht}" >> $file
  echo "sa35b0      ${hht}" >> $file
  echo "sa45b0      ${hht}" >> $file
  echo "s015b0      ${hht}" >> $file
  echo "s02100      ${hht}" >> $file
  echo "s025b0      ${hht}" >> $file
  echo "sa15b1      ${hht}" >> $file
  echo "sa25b1      ${hht}" >> $file
  echo "sa35b1      ${hht}" >> $file
  echo "sa45b1      ${hht}" >> $file
  echo "sa1100      ${hht}" >> $file
  echo "sa2100      ${hht}" >> $file
  echo "sa3100      ${hht}" >> $file
  echo "sa4100      ${hht}" >> $file
  echo "B02170      ${hht}" >> $file
  echo "B02180      ${hht}" >> $file
  echo "W00092      ${hht}" >> $file
  echo "W00093      ${hht}" >> $file
  echo "W00094      ${hht}" >> $file
  echo "S00070      ${hht}" >> $file



 let hh=${hh}+6
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
  echo "nomodata         "  >> $file
