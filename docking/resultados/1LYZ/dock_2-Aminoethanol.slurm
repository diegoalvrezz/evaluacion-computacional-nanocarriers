#!/bin/bash
#SBATCH --job-name=1LYZ_2-Aminoeth
#SBATCH --partition=lusi2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4000M
#SBATCH --time=4:00:00
#SBATCH --output=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/1LYZ/2-Aminoethanol_slurm_%j.out
#SBATCH --error=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/1LYZ/2-Aminoethanol_slurm_%j.err

/lusitania_homes/ubu-iccram/ubu-iccram01/bin/vina \
    --receptor /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/receptores/1LYZ.pdbqt \
    --ligand   /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/ligandos/2-Aminoethanol.pdbqt \
    --center_x 17.0 \
    --center_y 10.0 \
    --center_z 25.0 \
    --size_x   20.0 \
    --size_y   20.0 \
    --size_z   20.0 \
    --exhaustiveness 32 \
    --num_modes 9 \
    --energy_range 3     --cpu 8 \
    --out /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/1LYZ/2-Aminoethanol_docked.pdbqt \


echo "DONE: 1LYZ x 2-Aminoethanol"
