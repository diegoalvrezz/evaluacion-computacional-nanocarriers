#!/usr/bin/env python3
"""
lanzar_docking.py
=================
Lanza docking AutoDock Vina de 58 moléculas × 6 proteínas = 348 dockings
en el cluster HPC mediante jobs SLURM independientes.

Grid boxes basados en sitios activos de la literatura:
  7A65  (P-gp)    → sitio de unión a sustrato (Drug-binding pocket)
  1TQN  (CYP3A4)  → sitio activo hemo
  1CX8  (TfR1)    → sitio de unión apotransferrina
  4LRH  (FRα)     → sitio de unión folato
  1LYZ  (Lisozima)→ sitio activo (Asp52-Glu35)
  1AO6  (HSA)     → Sudlow site I (sitio IIA)

Uso:
  python3 lanzar_docking.py
"""

import os
import subprocess

# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURACIÓN
# ═══════════════════════════════════════════════════════════════════════════

BASE_DIR     = os.path.expanduser("~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo")
LIGANDOS_DIR = os.path.join(BASE_DIR, "docking/ligandos")
RECEPT_DIR   = os.path.join(BASE_DIR, "docking/receptores")
OUTPUT_DIR   = os.path.join(BASE_DIR, "docking/resultados")
VINA         = os.path.expanduser("~/bin/vina")
PARTITION    = "lusi2"
NCORES       = 8
MAX_HOURS    = 4

# ── Grid boxes por proteína ───────────────────────────────────────────────
# Formato: (center_x, center_y, center_z, size_x, size_y, size_z)
# Unidades: Angstrom
GRID_BOXES = {
    "7A65": {
        "name": "P-glicoproteína (P-gp/MDR1)",
        "site": "Drug-binding pocket (transmembrane)",
        # Centro del sitio de unión a sustrato en P-gp humana
        "cx": 160.0, "cy": 177.0, "cz": 144.0,
        "sx": 60.0,  "sy": 60.0,   "sz": 60.0,
    },
    "1TQN": {
        "name": "CYP3A4",
        "site": "Sitio activo hemo (F304-R212)",
        # Centro del sitio activo sobre el grupo hemo
        "cx": -15.5, "cy": -22.0, "cz": -11.0,
        "sx": 25.0,  "sy": 25.0,  "sz": 25.0,
    },
    "1CX8": {
        "name": "Receptor de transferrina (TfR1)",
        "site": "Sitio de unión apotransferrina",
        # Centro del sitio de unión en el dominio apical
        "cx": -28.4, "cy": -38.9, "cz": 147.7,
        "sx": 30.0,  "sy": 30.0,  "sz": 30.0,
    },
    "4LRH": {
        "name": "Receptor de folato α (FRα)",
        "site": "Sitio de unión folato",
        # Centro del bolsillo de unión al folato
        "cx": 5.2,   "cy": 25.8,  "cz": -1.4,
        "sx": 22.0,  "sy": 22.0,  "sz": 22.0,
    },
    "1LYZ": {
        "name": "Lisozima humana",
        "site": "Sitio activo (Asp52-Glu35)",
        # Centro del sitio activo de lisozima
        "cx": 17.0,  "cy": 10.0,  "cz": 25.0,
        "sx": 20.0,  "sy": 20.0,  "sz": 20.0,
    },
    "1AO6": {
        "name": "HSA (Albúmina sérica humana)",
        "site": "Sudlow site I (subdominio IIA)",
        # Centro del Sudlow site I
        "cx": 7.0,   "cy": 10.0,  "cz": 18.0,
        "sx": 25.0,  "sy": 25.0,  "sz": 25.0,
    },
}

# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    ligandos = sorted([f for f in os.listdir(LIGANDOS_DIR) if f.endswith(".pdbqt")])
    proteinas = sorted(GRID_BOXES.keys())

    print(f"Ligandos:  {len(ligandos)}")
    print(f"Proteínas: {len(proteinas)}")
    print(f"Total dockings: {len(ligandos) * len(proteinas)}\n")

    submitted = 0

    for prot_id, grid in GRID_BOXES.items():
        receptor = os.path.join(RECEPT_DIR, f"{prot_id}.pdbqt")
        if not os.path.exists(receptor):
            print(f"⚠️  Receptor no encontrado: {receptor}")
            continue

        prot_out = os.path.join(OUTPUT_DIR, prot_id)
        os.makedirs(prot_out, exist_ok=True)

        for lig_file in ligandos:
            mol_name = lig_file.replace(".pdbqt", "")
            ligando  = os.path.join(LIGANDOS_DIR, lig_file)
            out_file = os.path.join(prot_out, f"{mol_name}_docked.pdbqt")
            log_file = os.path.join(prot_out, f"{mol_name}_vina.log")

            # Saltar si ya está hecho
            if os.path.exists(out_file):
                print(f"  ⊘ {prot_id} × {mol_name} → ya existe")
                continue

            # Crear script SLURM
            slurm_path = os.path.join(prot_out, f"dock_{mol_name[:20]}.slurm")
            slurm_content = f"""#!/bin/bash
#SBATCH --job-name={prot_id}_{mol_name[:10]}
#SBATCH --partition={PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={NCORES}
#SBATCH --mem=4000M
#SBATCH --time={MAX_HOURS}:00:00
#SBATCH --output={prot_out}/{mol_name}_slurm_%j.out
#SBATCH --error={prot_out}/{mol_name}_slurm_%j.err

{VINA} \\
    --receptor {receptor} \\
    --ligand   {ligando} \\
    --center_x {grid['cx']} \\
    --center_y {grid['cy']} \\
    --center_z {grid['cz']} \\
    --size_x   {grid['sx']} \\
    --size_y   {grid['sy']} \\
    --size_z   {grid['sz']} \\
    --exhaustiveness 32 \\
    --num_modes 9 \\
    --energy_range 3 \
    --cpu 8 \\
    --out {out_file} \\


echo "DONE: {prot_id} x {mol_name}"
"""
            with open(slurm_path, "w") as f:
                f.write(slurm_content)

            # Lanzar
            result = subprocess.run(
                ["sbatch", slurm_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
            )
            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1]
                print(f"  ✓ {prot_id} × {mol_name} → job {job_id}")
                submitted += 1
            else:
                print(f"  ✗ {prot_id} × {mol_name} → ERROR: {result.stderr.strip()}")

    print(f"\n{'='*50}")
    print(f"Jobs lanzados: {submitted}")
    print(f"Monitoriza: squeue -u ubu-iccram01 | wc -l")
    print(f"Resultados: {OUTPUT_DIR}")
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
