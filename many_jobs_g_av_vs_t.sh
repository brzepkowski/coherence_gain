#!/bin/bash

gfortran coherence_gain.f90 g_av_vs_t.f90 *.f

# ./a.out MIN 0.0992499962 0.0992499962 70.0 500000
# ./a.out MAX 0.0983599946 0.0983599946 70.0 500000 # Max occured before min
# python3 plot_vs_t.py 0.0992499962 0.0983599946 70.0
#
# ./a.out MIN 0.198489994 0.198489994 70.0 500000
# ./a.out MAX 0.199550003 0.199550003 70.0 500000 # Max occured after min
# python3 plot_vs_t.py 0.198489994 0.199550003 70.0
#
# ./a.out MIN 0.299809992 0.299809992 70.0 500000
# ./a.out MAX 0.298790008 0.298790008 70.0 500000 # Max occured before min
# python3 plot_vs_t.py 0.299809992 0.298790008 70.0
#
# ./a.out MIN 0.399070024 0.399070024 70.0 500000
# ./a.out MAX 0.398039997 0.398039997 70.0 500000 # Max occured before min
# python3 plot_vs_t.py 0.399070024 0.398039997 70.0
#
# ./a.out MIN 0.599650025 0.599650025 70.0 500000
# ./a.out MAX 0.598609984 0.598609984 70.0 500000 # Max occured before min
# python3 plot_vs_t.py 0.599650025 0.598609984 70.0
#
# ./a.out MIN 0.694770038 0.694770038 70.0 500000
# ./a.out MAX 0.699940026 0.699940026 70.0 500000 # Max occured after min
# python3 plot_vs_t.py 0.694770038 0.699940026 70.0
#
# ./a.out MIN 0.794030011 0.794030011 70.0 500000
# ./a.out MAX 0.795059979 0.795059979 70.0 500000 # Max occured after min
# python3 plot_vs_t.py 0.794030011 0.795059979 70.0

./a.out MIN 0.893280029 0.893280029 70.0 500000
./a.out MAX 0.894320011 0.894320011 70.0 500000 # Max occured after min
# python3 plot_vs_t.py 0.893280029 0.894320011 70.0

./a.out MIN 0.998749971 0.998749971 70.0 500000
./a.out MAX 0.999779999 0.999779999 70.0 500000 # Max occured after min
# python3 plot_vs_t.py 0.998749971 0.999779999 70.0

./a.out MIN 1.49916995 1.49916995 70.0 500000
./a.out MAX 1.49813998 1.49813998 70.0 500000 # Max occured before min
# python3 plot_vs_t.py 1.49916995 1.49813998 70.0

./a.out MIN 1.99958992 1.99958992 70.0 500000
./a.out MAX 1.99855995 1.99855995 70.0 500000 # Max occured before min
# python3 plot_vs_t.py 1.99958992 1.99855995 70.0

./a.out MIN 2.49794006 2.49794006 70.0 500000
./a.out MAX 2.49898005 2.49898005 70.0 500000 # Max occured after min
# python3 plot_vs_t.py 2.49794006 2.49898005 70.0

./a.out MIN 2.99836016 2.99836016 70.0 500000
./a.out MAX 2.99939013 2.99939013 70.0 500000 # Max occured after min
# python3 plot_vs_t.py 2.99836016 2.99939013 70.0

./a.out MIN 3.49877 3.49877 70.0 500000
./a.out MAX 3.49980998 3.49980998 70.0 500000 # Max occured after min
# python3 plot_vs_t.py 3.49877 3.49980998 70.0
