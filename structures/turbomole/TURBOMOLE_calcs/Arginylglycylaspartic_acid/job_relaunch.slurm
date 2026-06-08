#!/bin/bash
#SBATCH --job-name=Arginylglycylas
#SBATCH --partition=lusi2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28000M
#SBATCH --time=48:00:00
#SBATCH --output=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Arginylglycylaspartic_acid/slurm_relaunch_%j.out
#SBATCH --error=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Arginylglycylaspartic_acid/slurm_relaunch_%j.err

export TURBODIR=~/MASTER_SALUD_2026/TmoleX2024/TURBOMOLE
export TURBOMOLE_SYSNAME=x86_64-unknown-linux-gnu_smp
export PATH=$TURBODIR/scripts:$TURBODIR/bin/x86_64-unknown-linux-gnu_smp:$PATH
export PARA_ARCH=SMP
export PARNODES=16
export OMP_NUM_THREADS=16

cd /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Arginylglycylaspartic_acid

echo "=== Inicio: $(date) — Arginylglycylaspartic_acid ==="

# Regenerar orbitales iniciales
actual > /dev/null 2>&1 || true

# Etapa 1
jobex -ri -c 300 > jobex_output.log 2>&1

if [ -f "GEO_OPT_CONVERGED" ]; then
    echo "ETAPA 1: CONVERGIDO"
    
    # Etapa 2
    mkdir -p SP
    cp coord basis auxbasis mos SP/ 2>/dev/null || true
    cd SP
    cp ../control ./control
    python3 -c "
with open('control') as f:
    content = f.read()
lines = [l for l in content.split('\n')
         if not any(l.strip().startswith(k) for k in
         ['$optimize','$statpt','$cosmo','$cosmo_out','$moldenfile','$scfdamp'])]
content = '\n'.join(lines)
additions = '$cosmo\n  epsilon=infinity\n$cosmo_out file=Arginylglycylaspartic_acid_SP.cosmo\n$scfconv 7\n'
content = content.replace('$end', additions + '$end')
with open('control', 'w') as f:
    f.write(content)
"
    ridft > ridft_SP.log 2>&1
    cd /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Arginylglycylaspartic_acid

    # Etapa 3
    mkdir -p FREQ
    cp coord basis auxbasis control mos FREQ/ 2>/dev/null || true
    cd FREQ
    python3 -c "
with open('control') as f:
    content = f.read()
lines = [l for l in content.split('\n')
         if not any(l.strip().startswith(k) for k in ['$optimize','$statpt','$scfdamp'])]
content = '\n'.join(lines)
with open('control', 'w') as f:
    f.write(content)
"
    ridft > ridft_FREQ.log 2>&1
    aoforce > aoforce.log 2>&1
    cd /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Arginylglycylaspartic_acid
else
    echo "ETAPA 1: NO CONVERGIÓ"
fi

echo "=== Fin: $(date) ==="
