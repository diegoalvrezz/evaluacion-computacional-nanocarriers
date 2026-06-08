"""
combinar_descriptores.py
========================
Combina todos los datasets de descriptores en un único DataFrame
listo para los modelos de Machine Learning.

Datasets:
    - descriptores_2D.csv          (60 mol × 38 desc)
    - sigma_profiles.csv           (58 mol × 100 desc)
    - fingerprints_maccs.csv       (60 mol × 167 desc)
    - fingerprints_morgan.csv      (60 mol × 2048 desc)
    - descriptores_3D_mordred.csv  (56 mol × 1826 desc)
    - docking_energias.csv         (57 mol × 6 desc)
    - base_molecular_pubchem.csv   (metadata: grupo)

Salida:
    dataset_completo.csv     → todos los descriptores
    dataset_ML.csv           → solo moléculas con todos los datos + limpieza
    dataset_por_grupos.csv   → con columna de grupo molecular
"""

import os
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.expanduser("~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo")

# ── Cargar datasets ───────────────────────────────────────────
print("Cargando datasets...")

df_2d     = pd.read_csv(os.path.join(BASE, "descriptores_2D.csv"))
df_sigma  = pd.read_csv(os.path.join(BASE, "sigma_profiles.csv"))
df_maccs  = pd.read_csv(os.path.join(BASE, "fingerprints_maccs.csv"))
df_morgan = pd.read_csv(os.path.join(BASE, "fingerprints_morgan.csv"))
df_3d     = pd.read_csv(os.path.join(BASE, "descriptores_3D_mordred.csv"))
df_dock   = pd.read_csv(os.path.join(BASE, "docking_energias.csv"))
df_meta   = pd.read_csv(os.path.join(BASE, "base_molecular_pubchem.csv"))

print(f"  2D RDKit:      {df_2d.shape}")
print(f"  Sigma profiles:{df_sigma.shape}")
print(f"  MACCS:         {df_maccs.shape}")
print(f"  Morgan:        {df_morgan.shape}")
print(f"  3D Mordred:    {df_3d.shape}")
print(f"  Docking ΔG:    {df_dock.shape}")

# ── Normalizar nombres ────────────────────────────────────────
# Usar 'nombre' como clave en todos
def norm(df, col):
    df = df.copy()
    df[col] = df[col].str.strip()
    return df

df_2d    = norm(df_2d, "nombre")
df_sigma = norm(df_sigma, "nombre")
df_maccs = norm(df_maccs, "nombre")
df_morgan= norm(df_morgan, "nombre")
df_3d    = norm(df_3d, "nombre")
df_dock  = norm(df_dock, "molecula").rename(columns={"molecula": "nombre"})
df_meta  = df_meta[["nombre_entrada", "grupo"]].rename(
    columns={"nombre_entrada": "nombre"})
df_meta  = norm(df_meta, "nombre")

# ── Merge progresivo ──────────────────────────────────────────
print("\nCombinando datasets...")

# Base: 2D descriptors
df = df_2d.copy()

# Añadir grupo desde metadata
df = df.merge(df_meta[["nombre", "grupo"]], on="nombre", how="left", suffixes=("", "_meta"))
if "grupo_meta" in df.columns:
    df["grupo"] = df["grupo"].fillna(df["grupo_meta"])
    df.drop(columns=["grupo_meta"], inplace=True)

# Sigma profiles
df_sigma_clean = df_sigma.drop(columns=["nombre"], errors="ignore").rename(
    columns={c: c for c in df_sigma.columns})
df = df.merge(df_sigma.rename(columns={"nombre": "nombre"}),
              on="nombre", how="left", suffixes=("", "_sigma"))

# MACCS (quitar columna grupo duplicada)
df_maccs_clean = df_maccs.drop(columns=["grupo"], errors="ignore")
df = df.merge(df_maccs_clean, on="nombre", how="left", suffixes=("", "_maccs"))

# Morgan (quitar columna grupo duplicada)
df_morgan_clean = df_morgan.drop(columns=["grupo"], errors="ignore")
df = df.merge(df_morgan_clean, on="nombre", how="left", suffixes=("", "_morgan"))

# 3D Mordred
df = df.merge(df_3d, on="nombre", how="left", suffixes=("", "_3d"))

# Docking energías
df = df.merge(df_dock, on="nombre", how="left", suffixes=("", "_dock"))

print(f"  Dataset completo: {df.shape}")

# ── Guardar dataset completo ──────────────────────────────────
df.to_csv(os.path.join(BASE, "dataset_completo.csv"), index=False)
print(f"  ✓ dataset_completo.csv guardado")

# ── Dataset ML (solo moléculas con docking completo) ─────────
docking_cols = [c for c in df.columns if c.startswith("dG_")]
df_ml = df.dropna(subset=docking_cols).copy()

# Eliminar columnas con >50% NaN
threshold = len(df_ml) * 0.5
df_ml = df_ml.dropna(axis=1, thresh=threshold)

# Eliminar columnas de varianza cero
numeric_cols = df_ml.select_dtypes(include=[np.number]).columns
df_ml = df_ml.loc[:, (df_ml[numeric_cols].std() != 0) | ~df_ml.columns.isin(numeric_cols)]

# Convertir columnas numéricas con errores a NaN
for col in df_ml.select_dtypes(include=[object]).columns:
    if col not in ["nombre", "grupo", "smiles"]:
        df_ml[col] = pd.to_numeric(df_ml[col], errors="coerce")

df_ml.to_csv(os.path.join(BASE, "dataset_ML.csv"), index=False)

print(f"\n{'='*60}")
print(f"  ✓ Moléculas en dataset completo: {len(df)}")
print(f"  ✓ Moléculas en dataset ML:       {len(df_ml)}")
print(f"  ✓ Features totales (ML):         {len(df_ml.columns)-3}")  # -nombre, grupo, smiles
print(f"\n  Columnas de docking:")
for col in docking_cols:
    vals = df_ml[col].dropna()
    print(f"    {col}: {len(vals)} valores, media={vals.mean():.2f}")
print(f"\n  ✓ Archivos guardados:")
print(f"    - dataset_completo.csv")
print(f"    - dataset_ML.csv")
print(f"{'='*60}")
