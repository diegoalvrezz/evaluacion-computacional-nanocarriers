#!/bin/bash
#SBATCH --job-name=Chondr_SF
#SBATCH --partition=lusi2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28000M
#SBATCH --time=72:00:00
#SBATCH --output=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Chondroitin_SVP/sp_freq_%j.out

export TURBODIR=~/MASTER_SALUD_2026/TmoleX2024/TURBOMOLE
export TURBOMOLE_SYSNAME=x86_64-unknown-linux-gnu_smp
export PATH=$TURBODIR/scripts:$TURBODIR/bin/x86_64-unknown-linux-gnu_smp:$PATH
export PARA_ARCH=SMP
export PARNODES=16
export OMP_NUM_THREADS=16

DIR=~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Chondroitin_SVP
cd $DIR

mkdir -p SP
cp coord basis auxbasis mos SP/
cd SP
cp ../control ./control
sed -i '/^\$optimize/d;/^\$statpt/d;/^\$cosmo/d;/^\$cosmo_out/d;/^\$moldenfile/d;/^\$actual/d' control
sed -i 's/\$end/\$cosmo\n  epsilon=infinity\n\$cosmo_out file=Chondroitin_SP.cosmo\n\$scfconv 7\n\$end/' control
ridft > ridft_SP.log 2>&1
[ -f "Chondroitin_SP.cosmo" ] && echo "SP: OK" || echo "SP: FALLO"
cd $DIR

mkdir -p FREQ
cp coord basis auxbasis control mos FREQ/
cd FREQ
sed -i '/^\$optimize/d;/^\$statpt/d;/^\$actual/d' control
ridft > ridft_FREQ.log 2>&1
aoforce > aoforce.log 2>&1
[ -f "vibspectrum" ] && echo "FREQ: OK" || echo "FREQ: FALLO"
cd $DIR

echo -e "\nY\nY\n" | tm2molden > /dev/null 2>&1
[ -f "molden.input" ] && cp molden.input ~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/molden_files/Chondroitin.molden && echo "MOLDEN: OK" || echo "MOLDEN: FALLO"
