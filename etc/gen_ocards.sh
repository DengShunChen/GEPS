#!/bin/bash

#set -x
ofn='ocards_new'

tau_end=1080

rm -f ${ofn}

#=========================
tau_dt=6
dmslist="B02100 B02500 B02510 B10200 B10210 SSL010 X00340 X00590"
for varkey in  ${dmslist}
do
echo "varkey = ${varkey}"
for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
  ctau=`printf "%3.3i" ${tau}`
  echo "${varkey}      ${ctau}" >> ${ofn} 
done
done

#=========================
#  0-384/6  396-1080/12
dmslist="S00100 B00100 B00200 B00210 B00010"
for varkey in  ${dmslist}
do
echo "varkey = ${varkey}"
tau_dt=6
for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
  ctau=`printf "%3.3i" ${tau}`
  if [[ ${tau} -ge 384 ]];then tau_dt=12;fi
  echo "${varkey}      ${ctau}" >> ${ofn}
done
done

#=========================precipitation
tau_dt=12
dmslist="B00620 B00630 B00640"
for varkey in  ${dmslist}
do
echo "varkey = ${varkey}"
for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
  ctau=`printf "%3.3i" ${tau}`
  echo "${varkey}      ${ctau}" >> ${ofn} 
done
done

#=========================snow depth
tau_dt=24
dmslist="B00650"
for varkey in  ${dmslist}
do
echo "varkey = ${varkey}"
for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
  ctau=`printf "%3.3i" ${tau}`
  echo "${varkey}      ${ctau}" >> ${ofn} 
done
done
#=========================
varkey='phi'
echo "varkey = ${varkey}"
for plevel in 1000 925 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10 
do 
  tau_dt=6
  for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
    ctau=`printf "%3.3i" ${tau}`
    if [[ ${tau} -ge 384 ]];then tau_dt=12;fi
    if [[ ${plevel} = '850'  ]];then tau_dt=3 ;fi
    if [[ ${plevel} = '500'  ]];then tau_dt=3 ;fi
    if [[ ${plevel} -lt 200 ]];then tau_dt=24 ;fi
    clevel=`printf "%4i" ${plevel} `
    echo "${varkey}    ${clevel} ${ctau}" >> ${ofn}
  done
done
#=========================
varkey='tmp'
echo "varkey = ${varkey}"
for plevel in 1000 925 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10
do
  tau_dt=6
  for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
    ctau=`printf "%3.3i" ${tau}`
    if [[ ${tau} -ge 384 ]];then tau_dt=12;fi
    if [[ ${plevel} = '1000' ]];then tau_dt=1 ;fi
    if [[ ${plevel} = '925'  ]];then tau_dt=1 ;fi
    if [[ ${plevel} = '850'  ]];then tau_dt=1 ;fi
    if [[ ${plevel} = '300'  ]];then tau_dt=6 ;fi
    if [[ ${plevel} -lt 200 ]];then tau_dt=24 ;fi
    clevel=`printf "%4i" ${plevel} `
    echo "${varkey}    ${clevel} ${ctau}" >> ${ofn}
  done
done
#=========================
varkey='wnd'
echo "varkey = ${varkey}"
for plevel in 1000 925 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10
do
  tau_dt=6
  for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
    ctau=`printf "%3.3i" ${tau}`
    if [[ ${tau} -ge 384 ]];then tau_dt=12;fi
    if [[ ${plevel} = '850'  ]];then tau_dt=6 ;fi
    if [[ ${plevel} = '300'  ]];then tau_dt=6 ;fi
    if [[ ${plevel} -lt 200 ]];then tau_dt=24 ;fi
    clevel=`printf "%4i" ${plevel} `
    echo "${varkey}    ${clevel} ${ctau}" >> ${ofn}
  done
done
#=========================
varkey='h2o'
echo "varkey = ${varkey}"
for plevel in 1000 925 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10
do
  tau_dt=6
  for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
    ctau=`printf "%3.3i" ${tau}`
    if [[ ${tau} -ge 384 ]];then tau_dt=12;fi
    if [[ ${plevel} = '1000' ]];then tau_dt=1 ;fi
    if [[ ${plevel} = '925'  ]];then tau_dt=1 ;fi
    if [[ ${plevel} = '850'  ]];then tau_dt=1 ;fi
    if [[ ${plevel} = '300'  ]];then tau_dt=6 ;fi
    if [[ ${plevel} -lt 200 ]];then tau_dt=24 ;fi
    clevel=`printf "%4i" ${plevel} `
    echo "${varkey}    ${clevel} ${ctau}" >> ${ofn}
  done
done
#=========================
varkey='vor'
echo "varkey = ${varkey}"
for plevel in 850 700 500
do
  tau_end=1080
  tau_dt=6
  for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
    ctau=`printf "%3.3i" ${tau}`
    if [[ ${tau} -ge 384 ]];then tau_dt=12;fi
    if [[ ${plevel} = '850'  ]];then tau_dt=3 ;fi
    if [[ ${plevel} = '700'  ]];then tau_dt=3 ;fi
    clevel=`printf "%4i" ${plevel} `
    echo "${varkey}    ${clevel} ${ctau}" >> ${ofn}
  done
done
#=========================
varkey='div'
echo "varkey = ${varkey}"
for plevel in 700 300 200
do
  tau_end=1080
  tau_dt=6
  for (( tau=0 ;tau<=${tau_end} ;tau=${tau}+${tau_dt} ));do
    ctau=`printf "%3.3i" ${tau}`
    if [[ ${tau} -ge 384 ]];then tau_dt=12;fi
    if [[ ${plevel} = '300'  ]];then tau_dt=6 ;fi
    clevel=`printf "%4i" ${plevel} `
    echo "${varkey}    ${clevel} ${ctau}" >> ${ofn}
  done
done
#=========================
echo "nomodata" >> ${ofn}
exit 0
