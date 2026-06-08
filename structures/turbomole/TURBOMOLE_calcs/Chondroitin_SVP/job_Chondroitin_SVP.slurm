#!/bin/bash
#SBATCH --job-name=Chondroit_SVP
#SBATCH --partition=lusi2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28000M
#SBATCH --time=48:00:00
#SBATCH --output=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Chondroitin_SVP/slurm_%j.out

export TURBODIR=~/MASTER_SALUD_2026/TmoleX2024/TURBOMOLE
export TURBOMOLE_SYSNAME=x86_64-unknown-linux-gnu_smp
export PATH=$TURBODIR/scripts:$TURBODIR/bin/x86_64-unknown-linux-gnu_smp:$PATH
export PARA_ARCH=SMP
export PARNODES=16
export OMP_NUM_THREADS=16

cd ~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Chondroitin_SVP
jobex -ri -c 1000 > jobex_output.log 2>&1
[ -f GEO_OPT_CONVERGED ] && echo "CONVERGED" || echo "FAILED"
