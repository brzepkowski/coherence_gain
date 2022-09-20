# coherence_gain

## WORKFLOW:
- $ gfortran g_av_vs_tau.f90 coherence_gain.f90 *.f
- $ ./many_jobs_g_av_vs_tau.sh
- $ python3 find_mins_and_maxes.py
- Manually find minima and maxima of interest (lying close to each other in given range).
- $ gfortran g_av_vs_t.f90 coherence_gain.f90 *.f
- $ ./many_jobs_g_av_vs_t.sh

### Note: In newer version of fortran it might be necessary to compile with the flag "-fallow-argument-mismatch", e.g.:
gfortran g_av_vs_tau.f90 coherence_gain.f90 *.f -fallow-argument-mismatch
