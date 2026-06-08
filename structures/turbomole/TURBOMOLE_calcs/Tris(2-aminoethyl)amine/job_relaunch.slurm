#!/bin/bash
#SBATCH --job-name=Tris(2-aminoeth
#SBATCH --partition=lusi2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28000M
#SBATCH --time=12:00:00
#SBATCH --output=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Tris(2-aminoethyl)amine/slurm_relaunch_%j.out
#SBATCH --error=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Tris(2-aminoethyl)amine/slurm_relaunch_%j.err

export TURBODIR=~/MASTER_SALUD_2026/TmoleX2024/TURBOMOLE
export TURBOMOLE_SYSNAME=x86_64-unknown-linux-gnu_smp
export PATH=$TURBODIR/scripts:$TURBODIR/bin/x86_64-unknown-linux-gnu_smp:$PATH
export PARA_ARCH=SMP
export PARNODES=16
export OMP_NUM_THREADS=16

cd /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Tris(2-aminoethyl)amine

# ETAPA 2 — Single Point COSMO-RS
mkdir -p SP
if [ ! -f "SP/Tris(2-aminoethyl)amine_SP.cosmo" ]; then
    echo "--- ETAPA 2: Single Point COSMO-RS ---"
    cp coord basis auxbasis mos SP/ 2>/dev/null || true
    cd SP
    cp ../control ./control
    python3 -c "
with open('control') as f:
    content = f.read()
lines = [l for l in content.split('\n')
         if not any(l.strip().startswith(k) for k in
         ['$optimize','$statpt','$cosmo','$cosmo_out','$moldenfile'])]
content = '\n'.join(lines)
additions = '$cosmo\n  epsilon=infinity\n$cosmo_out file=Tris(2-aminoethyl)amine_SP.cosmo\n$scfconv 7\n'
content = content.replace('$end', additions + '$end')
with open('control', 'w') as f:
    f.write(content)
"
    ridft > ridft_SP.log 2>&1
    if [ -f "Tris(2-aminoethyl)amine_SP.cosmo" ]; then
        echo "ETAPA 2: OK"
    else
        echo "ETAPA 2: FALLO"
    fi
    cd /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Tris(2-aminoethyl)amine
fi

# ETAPA 3 — Frecuencias IR
mkdir -p FREQ
if [ ! -f "FREQ/vibspectrum" ]; then
    echo "--- ETAPA 3: Frecuencias IR ---"
    cp coord basis auxbasis control mos FREQ/ 2>/dev/null || true
    cd FREQ
    python3 -c "
with open('control') as f:
    content = f.read()
lines = [l for l in content.split('\n')
         if not any(l.strip().startswith(k) for k in ['$optimize','$statpt'])]
content = '\n'.join(lines)
with open('control', 'w') as f:
    f.write(content)
"
    ridft > ridft_FREQ.log 2>&1
    aoforce > aoforce.log 2>&1
    if [ -f "vibspectrum" ]; then
        echo "ETAPA 3: OK"
    else
        echo "ETAPA 3: FALLO"
    fi
    cd /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Tris(2-aminoethyl)amine
fi
