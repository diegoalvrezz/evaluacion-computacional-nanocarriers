#!/bin/bash
#SBATCH --job-name=1TQN_N-Acetylgl
#SBATCH --partition=lusi2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4000M
#SBATCH --time=4:00:00
#SBATCH --output=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/1TQN/N-Acetylglucosamine_slurm_%j.out
#SBATCH --error=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/1TQN/N-Acetylglucosamine_slurm_%j.err

/lusitania_homes/ubu-iccram/ubu-iccram01/bin/vina \
    --receptor /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/receptores/1TQN.pdbqt \
    --ligand   /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/ligandos/N-Acetylglucosamine.pdbqt \
    --center_x -15.5 \
    --center_y -22.0 \
    --center_z -11.0 \
    --size_x   25.0 \
    --size_y   25.0 \
    --size_z   25.0 \
    --exhaustiveness 32 \
    --num_modes 9 \
    --energy_range 3     --cpu 8 \
    --out /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/1TQN/N-Acetylglucosamine_docked.pdbqt \


echo "DONE: 1TQN x N-Acetylglucosamine"
