#!/bin/bash

gfortran coherence_gain.f90 g_av_vs_tau.f90 *.f

./a.out 0.00 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.00 20.0 70.0 0
./a.out 0.10 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.10 20.0 70.0 0
./a.out 0.20 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.20 20.0 70.0 0
./a.out 0.30 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.30 20.0 70.0 0
./a.out 0.40 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.40 20.0 70.0 0
./a.out 0.50 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.50 20.0 70.0 0
./a.out 0.60 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.60 20.0 70.0 0
./a.out 0.70 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.70 20.0 70.0 0
./a.out 0.80 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.80 20.0 70.0 0
./a.out 0.90 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 0.90 20.0 70.0 0

./a.out 1.00 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.00 20.0 70.0 0
./a.out 1.10 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.10 20.0 70.0 0
./a.out 1.20 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.20 20.0 70.0 0
./a.out 1.30 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.30 20.0 70.0 0
./a.out 1.40 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.40 20.0 70.0 0
./a.out 1.50 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.50 20.0 70.0 0
./a.out 1.60 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.60 20.0 70.0 0
./a.out 1.70 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.70 20.0 70.0 0
./a.out 1.80 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.80 20.0 70.0 0
./a.out 1.90 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 1.90 20.0 70.0 0

./a.out 2.00 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.00 20.0 70.0 0
./a.out 2.10 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.10 20.0 70.0 0
./a.out 2.20 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.20 20.0 70.0 0
./a.out 2.30 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.30 20.0 70.0 0
./a.out 2.40 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.40 20.0 70.0 0
./a.out 2.50 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.50 20.0 70.0 0
./a.out 2.60 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.60 20.0 70.0 0
./a.out 2.70 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.70 20.0 70.0 0
./a.out 2.80 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.80 20.0 70.0 0
./a.out 2.90 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 2.90 20.0 70.0 0

./a.out 3.00 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.00 20.0 70.0 0
./a.out 3.10 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.10 20.0 70.0 0
./a.out 3.20 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.20 20.0 70.0 0
./a.out 3.30 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.30 20.0 70.0 0
./a.out 3.40 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.40 20.0 70.0 0
./a.out 3.50 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.50 20.0 70.0 0
./a.out 3.60 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.60 20.0 70.0 0
./a.out 3.70 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.70 20.0 70.0 0
./a.out 3.80 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.80 20.0 70.0 0
./a.out 3.90 20.0 70.0 10000
python3 -W ignore plot_vs_tau.py 3.90 20.0 70.0 0
