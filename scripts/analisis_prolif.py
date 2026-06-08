"""
analisis_prolif.py
==================
Analiza interacciones proteína-ligando para los top 5 candidatos ITI
usando ProLIF. Genera tablas de interacciones y gráficas para la memoria.

Uso:
    python3 analisis_prolif.py
"""

import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
import prolif as plf
from pathlib import Path

BASE     = Path.home() / "Desktop/TFM_DiegoVallina/docking_poses"
OUT_DIR  = Path.home() / "Desktop/TFM_DiegoVallina/placeholder/resultados_ML"
OUT_DIR.mkdir(exist_ok=True)

PROTEINAS = {
    "7A65": "P-gp MDR1",
    "1TQN": "CYP3A4",
    "1CX8": "TfR1",
    "4LRH": "FRalpha",
    "1LYZ": "Lisozima",
    "1AO6": "HSA",
}

MOLECULAS = [
    "Triethylene_glycol",
    "Tetraethylenepentamine",
    "Pentaethylenehexamine",
    "Chitotriose",
    "Chitobiose",
]

MOL_LABELS = {
    "Triethylene_glycol":     "Triethylene glycol",
    "Tetraethylenepentamine": "Tetraethylenepentamine",
    "Pentaethylenehexamine":  "Pentaethylenehexamine",
    "Chitotriose":            "Chitotriose",
    "Chitobiose":             "Chitobiose",
}

INTERACTION_COLORS = {
    "HBDonor":      "#2196F3",
    "HBAcceptor":   "#4CAF50",
    "Hydrophobic":  "#FF9800",
    "Cationic":     "#F44336",
    "Anionic":      "#9C27B0",
    "PiStacking":   "#00BCD4",
    "PiCation":     "#795548",
    "VdWContact":   "#9E9E9E",
}

def parse_pdbqt_to_pdb(pdbqt_path):
    """Convierte PDBQT a PDB en memoria filtrando líneas ATOM/HETATM del modelo 1."""
    lines = []
    in_model1 = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if "1" in line.split()[1]:
                    in_model1 = True
                else:
                    in_model1 = False
            if in_model1 and (line.startswith("ATOM") or line.startswith("HETATM")):
                lines.append(line[:66] + "\n")
            if line.startswith("ENDMDL") and in_model1:
                break
    return lines

resultados = []

print("Analizando interacciones con ProLIF...")

for prot_code, prot_name in PROTEINAS.items():
    pdb_path = BASE / f"{prot_code}.pdb"
    if not pdb_path.exists():
        print(f"  ✗ {prot_code}.pdb no encontrado")
        continue

    for mol_name in MOLECULAS:
        pdbqt_path = BASE / prot_code / f"{mol_name}_docked.pdbqt"
        if not pdbqt_path.exists():
            print(f"  ✗ {mol_name} en {prot_code} no encontrado")
            continue

        try:
            # Cargar proteína
            u_prot = mda.Universe(str(pdb_path))

            # Convertir ligando PDBQT → PDB temporal
            lig_lines = parse_pdbqt_to_pdb(pdbqt_path)
            tmp_lig = str(OUT_DIR / "tmp_lig.pdb")
            with open(tmp_lig, "w") as f:
                f.write("".join(lig_lines))
                f.write("END\n")

            u_lig = mda.Universe(tmp_lig)

            # ProLIF fingerprint
            fp = plf.Fingerprint(count=False)
            fp.run_from_atomgroups(
                u_lig.atoms,
                u_prot.atoms,
                residues="all"
            )

            df_fp = fp.to_dataframe()

            if df_fp.empty:
                print(f"  ⚠ Sin interacciones: {mol_name} - {prot_code}")
                continue

            # Extraer interacciones
            for col in df_fp.columns:
                if df_fp[col].any():
                    residuo = col[0] if isinstance(col, tuple) else str(col)
                    tipo    = col[1] if isinstance(col, tuple) and len(col) > 1 else "Unknown"
                    resultados.append({
                        "proteina":    prot_name,
                        "prot_code":   prot_code,
                        "molecula":    MOL_LABELS[mol_name],
                        "residuo":     str(residuo),
                        "interaccion": str(tipo),
                    })

            n_int = df_fp.any().sum()
            print(f"  ✓ {mol_name} - {prot_code}: {n_int} interacciones")

        except Exception as e:
            print(f"  ✗ Error {mol_name} - {prot_code}: {e}")

if not resultados:
    print("No se obtuvieron interacciones. Verifica los archivos de entrada.")
    exit(1)

df = pd.DataFrame(resultados)
df.to_csv(OUT_DIR / "prolif_interacciones.csv", index=False)
print(f"\nTotal interacciones: {len(df)}")

# ── Gráfica: heatmap de tipos de interacción por molécula×proteína ──
print("\nGenerando visualizaciones...")

pivot = df.groupby(["molecula", "proteina", "interaccion"]).size().reset_index(name="count")
pivot_total = pivot.groupby(["molecula", "proteina"])["count"].sum().reset_index()
matrix = pivot_total.pivot(index="molecula", columns="proteina", values="count").fillna(0)

fig, ax = plt.subplots(figsize=(12, 6))
sns.heatmap(matrix, annot=True, fmt=".0f", cmap="YlOrRd",
            ax=ax, linewidths=0.5,
            cbar_kws={"label": "Número de interacciones"})
ax.set_title("Número de interacciones proteína-ligando\n(Top 5 candidatos ITI, modo 1 docking)",
             fontsize=13, fontweight="bold")
ax.set_xlabel("")
ax.set_ylabel("")
plt.tight_layout()
plt.savefig(OUT_DIR / "prolif_heatmap_interacciones.png", dpi=150, bbox_inches="tight")
plt.close()
print("  ✓ prolif_heatmap_interacciones.png")

# ── Gráfica: tipos de interacción por proteína (barplot apilado) ──
tipo_count = df.groupby(["proteina", "interaccion"]).size().reset_index(name="count")
tipos = tipo_count["interaccion"].unique()
proteinas = list(PROTEINAS.values())

fig, ax = plt.subplots(figsize=(12, 6))
bottom = np.zeros(len(proteinas))

for tipo in tipos:
    vals = []
    for prot in proteinas:
        v = tipo_count[(tipo_count["proteina"]==prot) &
                       (tipo_count["interaccion"]==tipo)]["count"].sum()
        vals.append(v)
    color = INTERACTION_COLORS.get(tipo, "#888888")
    ax.bar(proteinas, vals, bottom=bottom, label=tipo, color=color, alpha=0.85)
    bottom += np.array(vals)

ax.set_xlabel("")
ax.set_ylabel("Número de interacciones", fontsize=12)
ax.set_title("Tipos de interacciones proteína-ligando por proteína\n(Top 5 candidatos ITI)",
             fontsize=13, fontweight="bold")
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
ax.tick_params(axis="x", rotation=15)
plt.tight_layout()
plt.savefig(OUT_DIR / "prolif_tipos_interaccion.png", dpi=150, bbox_inches="tight")
plt.close()
print("  ✓ prolif_tipos_interaccion.png")

# ── Tabla resumen por molécula ──
resumen = df.groupby(["molecula", "proteina"])["interaccion"].apply(
    lambda x: ", ".join(sorted(x.unique()))
).reset_index()
resumen.columns = ["Molécula", "Proteína", "Tipos de interacción"]
resumen.to_csv(OUT_DIR / "prolif_resumen_interacciones.csv", index=False)

print(f"\n✓ Análisis completado")
print(f"  - prolif_interacciones.csv")
print(f"  - prolif_heatmap_interacciones.png")
print(f"  - prolif_tipos_interaccion.png")
print(f"  - prolif_resumen_interacciones.csv")
