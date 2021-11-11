# coherence_gain
Remember, that we are storing a constant PERIOD in find_mins_and_maxes.py! This PERIOD should be changed before generating final plots, in which we include the sum in a formula for E = 1000 + SUM(|g_k|^2)!

WORKFLOW:
- $ gfortran g_av_vs_tau.f90 coherence_gain.f90 *.f
- $ ./many_jobs_g_av_vs_tau.sh
- $ python3 find_mins_and_maxes.py
- Manually find minima and maxima of interest (lying close to each other in given range).
- $ gfortran g_av_vs_t.f90 coherence_gain.f90 *.f
- $ ./many_jobs_g_av_vs_t.sh
