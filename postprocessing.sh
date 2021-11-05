#!/bin/bash

touch g_av_mins_and_maxes.dat
echo "period: 4.13566781E-03" > g_av_mins_and_maxes.dat

for filename in ./*
  do
    if  [[ $filename == ./g_av_vs_tau_* ]] && [[ $filename == *.dat ]] ;
    then
      echo "$filename" >> g_av_mins_and_maxes.dat
      tail -n 2 $filename >> g_av_mins_and_maxes.dat
    fi
  done
