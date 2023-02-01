#!/bin/bash

dmskeys='S00000
         S00070
         B10200
         B10210
         S00420
         S0042F
         S00430
         S0043F
         S00310
         S0031F
         S00320
         S0032F
         S003X0
         X00330
         X0033F
         X00340
         X0034F
         X00360
         B10100
         B10500
         B02100
         B02500
         B02510
         S00100
         SSL010
         B0062F
         B00623
         B00630
         B00640'

vars='h2o 
      phi 
      tmp 
      vor 
      wnd'

levs2='10 20 30 50 70 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 925 950 975 1000'


ft_ini=0
#ft_end=8784
ft_end=1080
#ft_end=1095
ft_gap=6

if [ -f ocards ]; then
  rm -f ocards
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
    echo $dmskey'     ' $ft >> ocards    
  done
done

for var in $vars
do
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

    for ft in $(seq ${ft_ini} 24 840)
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
      echo "${varlv}" $ft >> ocards
    done
  done


done

echo nomodata >> ocards

if [ -f gfsctl ]; then
  rm -f gfsctl
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
  echo $ft >> gfsctl
done



