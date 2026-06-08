#!/bin/bash
#SBATCH --job-name=1CX8_Tetraethyl
#SBATCH --partition=lusi2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4000M
#SBATCH --time=4:00:00
#SBATCH --output=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/1CX8/Tetraethylenepentamine_slurm_%j.out
#SBATCH --error=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/1CX8/Tetraethylenepentamine_slurm_%j.err

/lusitania_homes/ubu-iccram/ubu-iccram01/bin/vina \
    --receptor /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/receptores/1CX8.pdbqt \
    --ligand   /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/ligandos/Tetraethylenepentamine.pdbqt \
    --center_x -28.4 \
    --center_y -38.9 \
    --center_z 147.7 \
    --size_x   30.0 \
    --size_y   30.0 \
    --size_z   30.0 \
    --exhaustiveness 32 \
    --num_modes 9 \
    --energy_range 3     --cpu 8 \
    --out /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/1CX8/Tetraethylenepentamine_docked.pdbqt \


echo "DONE: 1CX8 x Tetraethylenepentamine"
