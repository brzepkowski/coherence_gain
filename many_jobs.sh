#!/bin/bash

# gfortran g_av_vs_t.f90 coherence_gain.f90 *.f -fallow-argument-mismatch -o g_av_vs_tau.out

# ./g_av_vs_tau.out 0.00 20.0 34.0 400000
# ./g_av_vs_tau.out 4.00 20.0 34.0 1000
#
# ./g_av_vs_tau.out 0.00 20.0 70.0 400000
# ./g_av_vs_tau.out 4.00 20.0 70.0 1000


# gfortran g_av_vs_t.f90 coherence_gain.f90 *.f -fallow-argument-mismatch -o g_av_vs_t.out

./g_av_vs_t.out MAX 4.00021982 0.00 34.0 1000000
./g_av_vs_t.out MIN 4.0012598 0.00 34.0 1000000
#
./g_av_vs_t.out MAX 4.00021982 0.00 70.0 1000000
./g_av_vs_t.out MIN 4.0012598 0.00 70.0 1000000

./g_av_vs_t.out MAX 0.424919993 0.00 34.0 1000000
./g_av_vs_t.out MIN 0.423879981 0.00 34.0 1000000

./g_av_vs_t.out MAX 0.424919993 0.00 70.0 1000000
./g_av_vs_t.out MIN 0.423879981 0.00 70.0 1000000
