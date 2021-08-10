#!/bin/bash

T=$1
t=$2
tau_min=$3
tau_max=$4

NAME=`echo "coherence_${T}_${t}_${tau_min}_${tau_max}"`

PBS="#!/bin/bash\n\
#PBS -N ${NAME}\n\
#PBS -l walltime=2:00:00\n\
#PBS -l select=1:ncpus=1:mem=256MB\n\
#PBS -l software=generate_data_1D_t_finite_detailed.py\n\
#PBS -m n\n\
cd \$PBS_O_WORKDIR\n\

python3 generate_data_1D_t_finite_detailed.py ${T} ${t} ${tau_min} ${tau_max}"

# Echo the string PBS to the function qsub, which submits it as a cluster job for you
# A small delay is included to avoid overloading the submission process

echo -e ${PBS} | qsub
#echo %{$PBS}
sleep 0.5
echo "done."
