#!/bin/bash

gfortran coherence_gain.f90 g_av_vs_t.f90 *.f

./a.out MAX 0.199560001 0.199560001 70.0 60000
./a.out MAX 0.298790008 0.298790008 70.0 60000
./a.out MAX 0.346660018 0.346660018 70.0 60000
./a.out MAX 0.499359995 0.499359995 70.0 60000
./a.out MAX 0.598609984 0.598609984 70.0 60000
./a.out MAX 0.699940026 0.699940026 70.0 60000
./a.out MAX 0.751639962 0.751639962 70.0 60000
./a.out MAX 0.898450017 0.898450017 70.0 60000
./a.out MAX 0.962559998 0.962559998 70.0 60000
./a.out MAX 1.09904003 1.09904003 70.0 60000
./a.out MAX 1.1983 1.1983 70.0 60000
./a.out MAX 1.29962003 1.29962003 70.0 60000
./a.out MAX 1.36785996 1.36785996 70.0 60000
./a.out MAX 1.49813998 1.49813998 70.0 60000
./a.out MAX 1.59738994 1.59738994 70.0 60000
./a.out MAX 1.69665003 1.69665003 70.0 60000
./a.out MAX 1.79798007 1.79798007 70.0 60000
./a.out MAX 1.89929998 1.89929998 70.0 60000
python3 plot_vs_t.py MAX 70.0
