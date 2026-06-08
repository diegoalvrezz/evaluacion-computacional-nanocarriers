"""
extraer_docking_energias.py
===========================
Extrae las energías de afinidad ΔG (kcal/mol) de los resultados de docking
AutoDock Vina y genera una matriz 58×6 (moléculas × proteínas).

Uso:
    python extraer_docking_energias.py

Salida:
    docking_energias.csv  → matriz con ΔG modo 1 (mejor pose) por proteína
    docking_todas_poses.csv → todas las poses (modos 1-9) por molécula y proteína
"""

import os
import glob
import pandas as pd
import numpy as np

# ── Configuración ─────────────────────────────────────────────
BASE_DOCKING = os.path.expanduser(
    "~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking/resultados"
)
OUTPUT_MATRIZ = os.path.expanduser(
    "~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking_energias.csv"
)
OUTPUT_POSES = os.path.expanduser(
    "~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/docking_todas_poses.csv"
)

PROTEINAS = ["7A65", "1TQN", "1CX8", "4LRH", "1LYZ", "1AO6"]

PROTEINA_NOMBRES = {
    "7A65": "P-gp_MDR1",
    "1TQN": "CYP3A4",
    "1CX8": "TfR1",
    "4LRH": "FRalpha",
    "1LYZ": "Lisozima",
    "1AO6": "HSA",
}

# ── Función para parsear PDBQT ────────────────────────────────
def extraer_energias_pdbqt(filepath):
    """
    Extrae todas las energías VINA RESULT de un archivo _docked.pdbqt.
    Devuelve lista de (modo, energia_kcal_mol).
    """
    energias = []
    modo = 0
    with open(filepath) as f:
        for line in f:
            if line.startswith("MODEL"):
                modo += 1
            elif "VINA RESULT:" in line:
                parts = line.split()
                try:
                    energia = float(parts[3])
                    energias.append((modo, energia))
                except (ValueError, IndexError):
                    continue
    return energias


# ── Main ──────────────────────────────────────────────────────
def main():
    resultados_matriz = {}   # {mol: {prot: dG_modo1}}
    resultados_poses  = []   # lista de filas con todas las poses

    for prot in PROTEINAS:
        prot_dir = os.path.join(os.path.expanduser(BASE_DOCKING), prot)
        if not os.path.exists(prot_dir):
            print(f"⚠️  Carpeta no encontrada: {prot_dir}")
            continue

        pdbqt_files = glob.glob(os.path.join(prot_dir, "*_docked.pdbqt"))
        print(f"\n{prot} ({PROTEINA_NOMBRES[prot]}): {len(pdbqt_files)} archivos")

        for filepath in sorted(pdbqt_files):
            filename = os.path.basename(filepath)
            mol_name = filename.replace("_docked.pdbqt", "")

            energias = extraer_energias_pdbqt(filepath)

            if not energias:
                print(f"  ⚠️  Sin energías: {mol_name}")
                continue

            # Modo 1 = mejor pose
            dG_modo1 = energias[0][1]

            # Guardar en matriz
            if mol_name not in resultados_matriz:
                resultados_matriz[mol_name] = {}
            resultados_matriz[mol_name][prot] = dG_modo1

            # Guardar todas las poses
            for modo, energia in energias:
                resultados_poses.append({
                    "molecula": mol_name,
                    "proteina": prot,
                    "proteina_nombre": PROTEINA_NOMBRES[prot],
                    "modo": modo,
                    "dG_kcal_mol": energia,
                })

            print(f"  ✓ {mol_name}: ΔG={dG_modo1:.3f} kcal/mol ({len(energias)} poses)")

    # ── Generar matriz 58×6 ───────────────────────────────────
    mols = sorted(resultados_matriz.keys())
    rows = []
    for mol in mols:
        row = {"molecula": mol}
        for prot in PROTEINAS:
            col = f"dG_{PROTEINA_NOMBRES[prot]}"
            row[col] = resultados_matriz[mol].get(prot, np.nan)
        rows.append(row)

    df_matriz = pd.DataFrame(rows)
    df_matriz.to_csv(OUTPUT_MATRIZ, index=False)

    # ── Generar CSV todas las poses ───────────────────────────
    df_poses = pd.DataFrame(resultados_poses)
    df_poses.to_csv(OUTPUT_POSES, index=False)

    # ── Resumen ───────────────────────────────────────────────
    print(f"\n{'='*55}")
    print(f"  ✓ Moléculas procesadas: {len(df_matriz)}")
    print(f"  ✓ Proteínas: {len(PROTEINAS)}")
    print(f"  ✓ Dockings con resultado: {df_matriz.notna().sum().sum() - len(df_matriz)}")
    print(f"  ✓ Matriz guardada: {OUTPUT_MATRIZ}")
    print(f"  ✓ Todas las poses: {OUTPUT_POSES}")
    print(f"\n  Estadísticas ΔG modo 1 (kcal/mol):")
    for prot in PROTEINAS:
        col = f"dG_{PROTEINA_NOMBRES[prot]}"
        vals = df_matriz[col].dropna()
        print(f"    {PROTEINA_NOMBRES[prot]:15s}: media={vals.mean():.2f}, "
              f"min={vals.min():.2f}, max={vals.max():.2f}")
    print(f"{'='*55}")


if __name__ == "__main__":
    main()
