#!/bin/bash
#SBATCH --job-name=Chitotriose
#SBATCH --partition=lusi2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28000M
#SBATCH --time=48:00:00
#SBATCH --output=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Chitotriose/slurm_%j.out
#SBATCH --error=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Chitotriose/slurm_%j.err

# ── Entorno TURBOMOLE ──────────────────────────────────────
export TURBODIR=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/TmoleX2024/TURBOMOLE
export TURBOMOLE_SYSNAME=x86_64-unknown-linux-gnu_smp
export PATH=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/TmoleX2024/TURBOMOLE/scripts:/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/TmoleX2024/TURBOMOLE/bin/x86_64-unknown-linux-gnu_smp:$PATH
export PARA_ARCH=SMP
export PARNODES=16
export OMP_NUM_THREADS=16

echo "========================================"
echo " Molécula : Chitotriose"
echo " Nodo     : $HOSTNAME"
echo " Inicio   : $(date)"
echo "========================================"

cd /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Chitotriose

# ════════════════════════════════════════════
# ETAPA 1 — Optimización B3LYP-D4/def2-TZVP
# ════════════════════════════════════════════
if [ ! -f "GEO_OPT_CONVERGED" ]; then
    echo "--- ETAPA 1: Optimización geometría ---"
    jobex -ri -c 300 > jobex_output.log 2>&1
    if [ -f "GEO_OPT_CONVERGED" ]; then
        echo "ETAPA 1: CONVERGIDO $(date)"
    else
        echo "ETAPA 1: NO CONVERGIÓ — abortando"
        exit 1
    fi
else
    echo "--- ETAPA 1: Ya convergida, saltando ---"
fi

# ════════════════════════════════════════════
# ETAPA 2 — Single Point COSMO-RS
# ════════════════════════════════════════════
mkdir -p SP
if [ ! -f "SP/Chitotriose_SP.cosmo" ]; then
    echo "--- ETAPA 2: Single Point COSMO-RS ---"

    # Copiar geometría optimizada y archivos de base
    cp coord     SP/
    cp basis     SP/
    cp auxbasis  SP/
    cp mos       SP/  2>/dev/null || true

    # Generar control para SP con COSMO
    cd SP
    # Usar control de la optimización como base y modificar
    cp ../control ./control_base
    python3 -c "
import sys
with open('control_base') as f:
    content = f.read()
# Limpiar directivas de optimización
lines = [l for l in content.split('\n')
         if not any(l.strip().startswith(k) for k in
         ['\$optimize','\$statpt','\$cosmo','\$cosmo_out','\$moldenfile'])]
content = '\n'.join(lines)
# Añadir COSMO-RS y export
additions = '\$cosmo\n  epsilon=infinity\n\$cosmo_out file=Chitotriose_SP.cosmo\n\$scfconv 7\n'
content = content.replace('\$end', additions + '\$end')
with open('control', 'w') as f:
    f.write(content)
print('control SP generado')
"
    # Lanzar single point
    ridft > ridft_SP.log 2>&1
    if [ -f "Chitotriose_SP.cosmo" ]; then
        echo "ETAPA 2: SP COSMO-RS completado $(date)"
    else
        echo "ETAPA 2: SP falló — revisar ridft_SP.log"
    fi
    cd /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Chitotriose
else
    echo "--- ETAPA 2: Ya completada, saltando ---"
fi

# ════════════════════════════════════════════
# ETAPA 3 — Frecuencias IR / Vibracionales
# ════════════════════════════════════════════
mkdir -p FREQ
if [ ! -f "FREQ/vibrational_spectrum" ]; then
    echo "--- ETAPA 3: Frecuencias IR (aoforce) ---"

    # Copiar geometría optimizada y archivos necesarios
    cp coord    FREQ/
    cp basis    FREQ/
    cp auxbasis FREQ/
    cp control  FREQ/
    cp mos      FREQ/ 2>/dev/null || true

    cd FREQ
    # Asegurarse que el control no tiene directivas de optimización
    python3 -c "
with open('control') as f:
    content = f.read()
lines = [l for l in content.split('\n')
         if not any(l.strip().startswith(k) for k in
         ['\$optimize','\$statpt'])]
content = '\n'.join(lines)
with open('control', 'w') as f:
    f.write(content)
"
    # Primero un single point para obtener los orbitales actualizados
    ridft > ridft_FREQ.log 2>&1
    # Luego calcular frecuencias
    aoforce > aoforce.log 2>&1

    if [ -f "vibrational_spectrum" ]; then
        echo "ETAPA 3: FRECUENCIAS completadas $(date)"
    else
        echo "ETAPA 3: aoforce falló — revisar aoforce.log"
    fi
    cd /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs/Chitotriose
else
    echo "--- ETAPA 3: Ya completada, saltando ---"
fi

echo "========================================"
echo " FIN: $(date)"
echo "========================================"
