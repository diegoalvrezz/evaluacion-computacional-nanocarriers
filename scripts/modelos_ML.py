"""
modelos_ML.py
=============
Modelos de Machine Learning para predicción de afinidad de docking
y clasificación de moléculas para drug delivery.

Modelos:
    1. Regresión de afinidad ΔG por proteína (RF, XGBoost, SVR, MLP)
    2. Índice de Transporte (ITI) - métrica compuesta
    3. Clasificador favorable/desfavorable para drug delivery

Uso:
    python modelos_ML.py

Salida:
    resultados_ML/
        ├── modelo1_regresion_resultados.csv
        ├── modelo2_ITI.csv
        ├── modelo3_clasificador_resultados.csv
        └── resumen_ML.txt
"""

import os
import warnings
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier, GradientBoostingRegressor
from sklearn.svm import SVR, SVC
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.model_selection import cross_val_score, KFold, StratifiedKFold
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_regression, f_classif
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier

warnings.filterwarnings('ignore')

# ── Configuración ─────────────────────────────────────────────
BASE      = os.path.expanduser("~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo")
INPUT_CSV = os.path.join(BASE, "dataset_ML.csv")
OUT_DIR   = os.path.join(BASE, "resultados_ML")
os.makedirs(OUT_DIR, exist_ok=True)

DOCKING_COLS = ["dG_P-gp_MDR1", "dG_CYP3A4", "dG_TfR1", "dG_FRalpha", "dG_Lisozima", "dG_HSA"]
META_COLS    = ["nombre", "grupo", "smiles", "nombre_norm"]
SEED         = 42
N_FOLDS      = 5

# ── Cargar datos ──────────────────────────────────────────────
print("="*60)
print("MODELOS ML - TFM Diego Vallina")
print("="*60)

df = pd.read_csv(INPUT_CSV)
print(f"\nDataset: {df.shape[0]} moléculas × {df.shape[1]} columnas")

# Features numéricas (excluir meta y docking como features)
feature_cols = [c for c in df.columns 
                if c not in META_COLS + DOCKING_COLS
                and df[c].dtype in [np.float64, np.int64, float, int]]

X_raw = df[feature_cols].values
print(f"Features disponibles: {len(feature_cols)}")

# ── Pipeline de preprocesado ──────────────────────────────────
def make_pipeline(model, k_features=200):
    return Pipeline([
        ("scaler",   RobustScaler()),
        ("selector", SelectKBest(f_regression, k=min(k_features, len(feature_cols)))),
        ("model",    model),
    ])

def make_pipeline_clf(model, k_features=200):
    return Pipeline([
        ("scaler",   RobustScaler()),
        ("selector", SelectKBest(f_classif, k=min(k_features, len(feature_cols)))),
        ("model",    model),
    ])

# ── MODELO 1: Regresión de afinidad ΔG ───────────────────────
print("\n" + "="*60)
print("MODELO 1: Regresión de afinidad ΔG por proteína")
print("="*60)

modelos_reg = {
    "Random Forest":  RandomForestRegressor(n_estimators=200, random_state=SEED, n_jobs=-1),
    "GradientBoosting":        GradientBoostingRegressor(n_estimators=200, random_state=SEED),
    "SVR":            SVR(kernel="rbf", C=10, epsilon=0.1),
    "MLP":            MLPRegressor(hidden_layer_sizes=(256,128,64), max_iter=500,
                                   random_state=SEED, early_stopping=True),
}

kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
resultados_reg = []

for prot in DOCKING_COLS:
    y = df[prot].values
    print(f"\n  Proteína: {prot}")
    
    for nombre_modelo, modelo in modelos_reg.items():
        pipe = make_pipeline(modelo)
        
        r2_scores  = cross_val_score(pipe, X_raw, y, cv=kf, scoring="r2")
        mae_scores = -cross_val_score(pipe, X_raw, y, cv=kf, scoring="neg_mean_absolute_error")
        rmse_scores= np.sqrt(-cross_val_score(pipe, X_raw, y, cv=kf, scoring="neg_mean_squared_error"))
        
        r2_mean   = r2_scores.mean()
        mae_mean  = mae_scores.mean()
        rmse_mean = rmse_scores.mean()
        
        resultados_reg.append({
            "proteina":   prot,
            "modelo":     nombre_modelo,
            "R2_mean":    round(r2_mean, 4),
            "R2_std":     round(r2_scores.std(), 4),
            "MAE_mean":   round(mae_mean, 4),
            "RMSE_mean":  round(rmse_mean, 4),
        })
        
        print(f"    {nombre_modelo:20s} R²={r2_mean:.3f}±{r2_scores.std():.3f}  "
              f"MAE={mae_mean:.3f}  RMSE={rmse_mean:.3f}")

df_reg = pd.DataFrame(resultados_reg)
df_reg.to_csv(os.path.join(OUT_DIR, "modelo1_regresion_resultados.csv"), index=False)
print(f"\n  ✓ Resultados guardados: modelo1_regresion_resultados.csv")

# ── MODELO 2: Índice de Transporte (ITI) ──────────────────────
print("\n" + "="*60)
print("MODELO 2: Índice de Transporte (ITI)")
print("="*60)

# ITI = (afinidad_targeting - afinidad_eflujo) normalizado
# Alta afinidad TfR1 y FRα (valores más negativos = mejor)
# Baja afinidad P-gp y CYP3A4 (valores menos negativos = mejor para evitar)

df_iti = df[["nombre", "grupo"] + DOCKING_COLS].copy()

# Normalizar cada columna a [0,1]
for col in DOCKING_COLS:
    min_val = df_iti[col].min()
    max_val = df_iti[col].max()
    df_iti[col + "_norm"] = (df_iti[col] - min_val) / (max_val - min_val)

# ITI: favorece alta afinidad targeting (TfR1, FRα) y baja afinidad eflujo (P-gp, CYP3A4)
# Más negativo = más afinidad → valor_norm más alto = más afinidad
targeting  = ["dG_TfR1_norm", "dG_FRalpha_norm"]
eflujo     = ["dG_P-gp_MDR1_norm", "dG_CYP3A4_norm"]
transporte = ["dG_HSA_norm"]
biocompat  = ["dG_Lisozima_norm"]

df_iti["score_targeting"]  = df_iti[targeting].mean(axis=1)
df_iti["score_eflujo"]     = df_iti[eflujo].mean(axis=1)
df_iti["score_transporte"] = df_iti[transporte].mean(axis=1)
df_iti["score_biocompat"]  = df_iti[biocompat].mean(axis=1)

# ITI compuesto: targeting y transporte suman, eflujo resta
df_iti["ITI"] = (
    0.40 * df_iti["score_targeting"] +
    0.25 * df_iti["score_transporte"] +
    0.20 * df_iti["score_biocompat"] -
    0.15 * df_iti["score_eflujo"]
)

# Normalizar ITI a [0,100]
iti_min = df_iti["ITI"].min()
iti_max = df_iti["ITI"].max()
df_iti["ITI_score"] = ((df_iti["ITI"] - iti_min) / (iti_max - iti_min) * 100).round(2)

# Clasificar
df_iti["perfil"] = pd.cut(df_iti["ITI_score"],
                           bins=[0, 33, 66, 100],
                           labels=["Desfavorable", "Moderado", "Favorable"],
                           include_lowest=True)

df_iti_out = df_iti[["nombre", "grupo"] + DOCKING_COLS + 
                     ["score_targeting", "score_eflujo", "score_transporte",
                      "score_biocompat", "ITI_score", "perfil"]].sort_values(
                     "ITI_score", ascending=False)

df_iti_out.to_csv(os.path.join(OUT_DIR, "modelo2_ITI.csv"), index=False)

print(f"\n  Top 10 moléculas por ITI:")
print(df_iti_out[["nombre", "grupo", "ITI_score", "perfil"]].head(10).to_string(index=False))
print(f"\n  Distribución de perfiles:")
print(df_iti["perfil"].value_counts().to_string())
print(f"\n  ✓ Resultados guardados: modelo2_ITI.csv")

# ── MODELO 3: Clasificador favorable/desfavorable ─────────────
print("\n" + "="*60)
print("MODELO 3: Clasificador favorable/desfavorable")
print("="*60)

# Target: Favorable (ITI >= 50) vs Desfavorable (ITI < 50)
y_clf = (df_iti["ITI_score"] >= 50).astype(int).values
print(f"\n  Distribución clases: Favorable={y_clf.sum()} | Desfavorable={len(y_clf)-y_clf.sum()}")

modelos_clf = {
    "Random Forest": RandomForestClassifier(n_estimators=200, random_state=SEED, n_jobs=-1),
    "GradientBoosting":       GradientBoostingClassifier(n_estimators=200, random_state=SEED),
    "SVC":           SVC(kernel="rbf", C=10, probability=True, random_state=SEED),
    "MLP":           MLPClassifier(hidden_layer_sizes=(256,128), max_iter=500,
                                   random_state=SEED, early_stopping=True),
}

skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
resultados_clf = []

for nombre_modelo, modelo in modelos_clf.items():
    pipe = make_pipeline_clf(modelo)
    
    acc_scores = cross_val_score(pipe, X_raw, y_clf, cv=skf, scoring="accuracy")
    auc_scores = cross_val_score(pipe, X_raw, y_clf, cv=skf, scoring="roc_auc")
    f1_scores  = cross_val_score(pipe, X_raw, y_clf, cv=skf, scoring="f1")
    
    resultados_clf.append({
        "modelo":      nombre_modelo,
        "Accuracy":    round(acc_scores.mean(), 4),
        "Acc_std":     round(acc_scores.std(), 4),
        "ROC_AUC":     round(auc_scores.mean(), 4),
        "AUC_std":     round(auc_scores.std(), 4),
        "F1":          round(f1_scores.mean(), 4),
        "F1_std":      round(f1_scores.std(), 4),
    })
    
    print(f"  {nombre_modelo:20s} Acc={acc_scores.mean():.3f}±{acc_scores.std():.3f}  "
          f"AUC={auc_scores.mean():.3f}  F1={f1_scores.mean():.3f}")

df_clf = pd.DataFrame(resultados_clf)
df_clf.to_csv(os.path.join(OUT_DIR, "modelo3_clasificador_resultados.csv"), index=False)
print(f"\n  ✓ Resultados guardados: modelo3_clasificador_resultados.csv")

# ── Resumen final ─────────────────────────────────────────────
print("\n" + "="*60)
print("RESUMEN FINAL")
print("="*60)

# Mejor modelo de regresión por proteína
print("\n  Mejor modelo regresión por proteína:")
for prot in DOCKING_COLS:
    best = df_reg[df_reg["proteina"]==prot].sort_values("R2_mean", ascending=False).iloc[0]
    print(f"    {prot:20s}: {best['modelo']:20s} R²={best['R2_mean']:.3f}")

# Mejor clasificador
best_clf = df_clf.sort_values("ROC_AUC", ascending=False).iloc[0]
print(f"\n  Mejor clasificador: {best_clf['modelo']} (AUC={best_clf['ROC_AUC']:.3f})")

# Top 5 moléculas ITI
print(f"\n  Top 5 moléculas para drug delivery:")
print(df_iti_out[["nombre","grupo","ITI_score","perfil"]].head(5).to_string(index=False))

# Guardar resumen
with open(os.path.join(OUT_DIR, "resumen_ML.txt"), "w") as f:
    f.write("RESUMEN MODELOS ML - TFM Diego Vallina\n")
    f.write("="*60 + "\n\n")
    f.write(f"Dataset: {df.shape[0]} moléculas, {len(feature_cols)} features\n\n")
    f.write("MODELO 1 - Regresión ΔG:\n")
    f.write(df_reg.sort_values(["proteina","R2_mean"], ascending=[True,False]).to_string(index=False))
    f.write("\n\nMODELO 2 - ITI Top 10:\n")
    f.write(df_iti_out[["nombre","grupo","ITI_score","perfil"]].head(10).to_string(index=False))
    f.write("\n\nMODELO 3 - Clasificador:\n")
    f.write(df_clf.sort_values("ROC_AUC", ascending=False).to_string(index=False))

print(f"\n  ✓ Resumen guardado: resumen_ML.txt")
print(f"  ✓ Todos los resultados en: {OUT_DIR}")
print("="*60)
