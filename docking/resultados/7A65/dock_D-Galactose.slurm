#!/bin/bash
#SBATCH --job-name=7A65_D-Galactos
#SBATCH --partition=lusi2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4000M
#SBATCH --time=4:00:00
#SBATCH --output=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/7A65/D-Galactose_slurm_%j.out
#SBATCH --error=/lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/7A65/D-Galactose_slurm_%j.err

/lusitania_homes/ubu-iccram/ubu-iccram01/bin/vina \
    --receptor /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/receptores/7A65.pdbqt \
    --ligand   /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/ligandos/D-Galactose.pdbqt \
    --center_x 160.0 \
    --center_y 177.0 \
    --center_z 144.0 \
    --size_x   60.0 \
    --size_y   60.0 \
    --size_z   60.0 \
    --exhaustiveness 32 \
    --num_modes 9 \
    --energy_range 3     --cpu 8 \
    --out /lusitania_homes/ubu-iccram/ubu-iccram01/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados/7A65/D-Galactose_docked.pdbqt \


echo "DONE: 7A65 x D-Galactose"
