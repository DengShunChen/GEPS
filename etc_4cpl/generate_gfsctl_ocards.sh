#!/bin/bash

. /nwpr/gfs/xb169/code_update/op_work_4cpl_t1/00_setDate.sh

dmskeys='B00100
         B00200
         B00210
         B00510
         B10200
         B10210
         B10100
         B10500
         B02100
         B02500
         B02510
         S00100
         SSL010'

vars='h2o 
      phi 
      tmp 
      vor 
      wnd'

levs=' 200    
       500       
       850      
      1000'

levs2='  10
         20
         30
         50
         70
        100
        150
        250
        300
        400
        700
        850
        925
       1000'


ft_ini=0
#ft_end=8784
ft_end=288 #${FCSTHR}
#ft_end=1095
ft_gap=6

ft_day=$(( ${ft_end}/24 ))

if [ -f ocards_${ft_day}d ]; then
  rm -f ocards_${ft_day}d
fi

for dmskey in $dmskeys
do
  for ft in $(seq ${ft_ini} ${ft_gap} ${ft_end})
  do
    if [ $ft -lt 10 ]; then
      ft='000'$ft
    else if [ $ft -lt 100 ]; then
      ft='00'$ft
    else if [ $ft -lt 1000 ]; then       
      ft='0'$ft
    fi
    fi
    fi
    echo $dmskey'     ' $ft >> ocards_${ft_day}d    
  done
done

for var in $vars
do
  for lv in $levs
  do
    if [ $lv -lt 100 ]; then
      varlv=$var'      '$lv
    else if [ $lv -lt 1000 ]; then
      varlv=$var'     '$lv
    else
      varlv=$var'    '$lv
    fi
    fi

    for ft in $(seq ${ft_ini} ${ft_gap} ${ft_end})
    do     
      if [ $ft -lt 10 ]; then
        ft='000'$ft
      else if [ $ft -lt 100 ]; then
        ft='00'$ft
      else if [ $ft -lt 1000 ]; then
        ft='0'$ft
      fi
      fi
      fi
      echo "${varlv}" $ft >> ocards_${ft_day}d
    done
  done

  for lv in $levs2
  do
    if [ $lv -lt 100 ]; then
      varlv=$var'      '$lv
    else if [ $lv -lt 1000 ]; then
      varlv=$var'     '$lv
    else
      varlv=$var'    '$lv
    fi
    fi

    for ft in $(seq ${ft_ini} ${ft_gap} ${ft_end})
    do
      if [ $ft -lt 10 ]; then
        ft='000'$ft
      else if [ $ft -lt 100 ]; then
        ft='00'$ft
      else if [ $ft -lt 1000 ]; then
        ft='0'$ft
      fi
      fi
      fi
      echo "${varlv}" $ft >> ocards_${ft_day}d
    done
  done


done

echo nomodata >> ocards_${ft_day}d

if [ -f gfsctl ]; then
  rm -f gfsctl
fi
for ft in $(seq ${ft_ini} 3 ${ft_end})
do
  if [ $ft -lt 10 ]; then
    ft='000'$ft
  else if [ $ft -lt 100 ]; then
    ft='00'$ft
  else if [ $ft -lt 1000 ]; then
    ft='0'$ft
  fi
  fi
  fi
  echo $ft >> gfsctl
done



