"""
modelos_ML_v2.py
================
Versión mejorada con:
- ITI corregido (penaliza moléculas pequeñas, ratio targeting/eflujo)
- SHAP values para interpretabilidad
- Gráficas para la memoria del TFM

Uso:
    python modelos_ML_v2.py
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier
from sklearn.svm import SVR, SVC
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.preprocessing import RobustScaler
from sklearn.model_selection import cross_val_score, KFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_regression, f_classif
import shap

warnings.filterwarnings('ignore')

# ── Configuración ─────────────────────────────────────────────
BASE    = os.path.expanduser("~/Desktop/TFM_DiegoVallina/placeholder")
OUT_DIR = os.path.join(BASE, "resultados_ML")
os.makedirs(OUT_DIR, exist_ok=True)

DOCKING_COLS = ["dG_P-gp_MDR1", "dG_CYP3A4", "dG_TfR1", "dG_FRalpha", "dG_Lisozima", "dG_HSA"]
META_COLS    = ["nombre", "grupo", "smiles", "nombre_norm"]
SEED = 42
N_FOLDS = 5
PALETTE = {"monomeros_degradacion": "#4CAF50", "oligomeros_unidades": "#2196F3",
           "farmacos_modelo": "#F44336", "ligandos_targeting": "#FF9800"}

# ── Cargar datos ──────────────────────────────────────────────
df = pd.read_csv(os.path.join(BASE, "dataset_ML.csv"))
feature_cols = [c for c in df.columns
                if c not in META_COLS + DOCKING_COLS
                and df[c].dtype in [np.float64, np.int64, float, int]]
X_raw = df[feature_cols].values
print(f"Dataset: {df.shape[0]} moléculas × {len(feature_cols)} features")

# ── ITI CORREGIDO ─────────────────────────────────────────────
print("\n=== MODELO 2: ITI CORREGIDO ===")

df_iti = df[["nombre", "grupo", "MW"] + DOCKING_COLS].copy()

# Normalizar docking a [0,1] (más negativo = más afinidad = valor más alto)
for col in DOCKING_COLS:
    mn, mx = df_iti[col].min(), df_iti[col].max()
    df_iti[col + "_n"] = (df_iti[col] - mn) / (mx - mn)

# Score targeting: TfR1 + FRα (alta afinidad es bueno)
df_iti["score_targeting"] = df_iti[["dG_TfR1_n", "dG_FRalpha_n"]].mean(axis=1)

# Score eflujo: P-gp + CYP3A4 (alta afinidad es MALO → invertir)
df_iti["score_eflujo"] = df_iti[["dG_P-gp_MDR1_n", "dG_CYP3A4_n"]].mean(axis=1)

# Score transporte: HSA (moderado es bueno)
df_iti["score_transporte"] = df_iti["dG_HSA_n"]

# Score biocompat: Lisozima
df_iti["score_biocompat"] = df_iti["dG_Lisozima_n"]

# Ratio targeting/eflujo (penaliza moléculas que también se unen a P-gp)
df_iti["ratio_TE"] = df_iti["score_targeting"] / (df_iti["score_eflujo"] + 0.01)

# Factor MW: penaliza moléculas muy pequeñas (<100 Da) y muy grandes (>1000 Da)
def mw_factor(mw):
    if mw < 50:   return 0.1
    if mw < 100:  return 0.4
    if mw < 150:  return 0.7
    if mw > 900:  return 0.8
    return 1.0

df_iti["mw_factor"] = df_iti["MW"].apply(mw_factor)

# ITI compuesto corregido
df_iti["ITI_raw"] = (
    0.35 * df_iti["score_targeting"] +
    0.25 * df_iti["ratio_TE"] / df_iti["ratio_TE"].max() +
    0.20 * df_iti["score_transporte"] +
    0.10 * df_iti["score_biocompat"] -
    0.10 * df_iti["score_eflujo"]
) * df_iti["mw_factor"]

# Normalizar a [0,100]
mn, mx = df_iti["ITI_raw"].min(), df_iti["ITI_raw"].max()
df_iti["ITI_score"] = ((df_iti["ITI_raw"] - mn) / (mx - mn) * 100).round(2)

df_iti["perfil"] = pd.cut(df_iti["ITI_score"],
                           bins=[0, 33, 66, 100],
                           labels=["Desfavorable", "Moderado", "Favorable"],
                           include_lowest=True)

df_iti_sorted = df_iti.sort_values("ITI_score", ascending=False)
print("\nTop 15 moléculas por ITI corregido:")
print(df_iti_sorted[["nombre", "grupo", "MW", "ITI_score", "perfil"]].head(15).to_string(index=False))

df_iti_sorted.to_csv(os.path.join(OUT_DIR, "modelo2_ITI_corregido.csv"), index=False)

# ── MODELO 1: Regresión (RF y GB) ────────────────────────────
print("\n=== MODELO 1: Regresión ΔG ===")

def make_pipe(model, k=200):
    return Pipeline([
        ("scaler",   RobustScaler()),
        ("selector", SelectKBest(f_regression, k=min(k, len(feature_cols)))),
        ("model",    model),
    ])

modelos_reg = {
    "Random Forest":     RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1),
    "GradientBoosting":  GradientBoostingRegressor(n_estimators=200, random_state=SEED),
    "SVR":               SVR(kernel="rbf", C=10, epsilon=0.1),
}

kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
res_reg = []

for prot in DOCKING_COLS:
    y = df[prot].values
    print(f"\n  {prot}:")
    for nm, mod in modelos_reg.items():
        pipe = make_pipe(mod)
        r2  = cross_val_score(pipe, X_raw, y, cv=kf, scoring="r2")
        mae = -cross_val_score(pipe, X_raw, y, cv=kf, scoring="neg_mean_absolute_error")
        rmse= np.sqrt(-cross_val_score(pipe, X_raw, y, cv=kf, scoring="neg_mean_squared_error"))
        res_reg.append({"proteina": prot, "modelo": nm,
                        "R2": round(r2.mean(),4), "R2_std": round(r2.std(),4),
                        "MAE": round(mae.mean(),4), "RMSE": round(rmse.mean(),4)})
        print(f"    {nm:20s} R²={r2.mean():.3f}±{r2.std():.3f}  MAE={mae.mean():.3f}")

df_reg = pd.DataFrame(res_reg)
df_reg.to_csv(os.path.join(OUT_DIR, "modelo1_regresion.csv"), index=False)

# ── MODELO 3: Clasificador ────────────────────────────────────
print("\n=== MODELO 3: Clasificador ===")

y_clf = (df_iti["ITI_score"] >= 50).astype(int).values
print(f"Favorable: {y_clf.sum()} | Desfavorable: {len(y_clf)-y_clf.sum()}")

def make_pipe_clf(model, k=200):
    return Pipeline([
        ("scaler",   RobustScaler()),
        ("selector", SelectKBest(f_classif, k=min(k, len(feature_cols)))),
        ("model",    model),
    ])

modelos_clf = {
    "Random Forest":    RandomForestClassifier(n_estimators=300, random_state=SEED, n_jobs=-1),
    "GradientBoosting": GradientBoostingClassifier(n_estimators=200, random_state=SEED),
    "SVC":              SVC(kernel="rbf", C=10, probability=True, random_state=SEED),
}

skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
res_clf = []
for nm, mod in modelos_clf.items():
    pipe = make_pipe_clf(mod)
    acc = cross_val_score(pipe, X_raw, y_clf, cv=skf, scoring="accuracy")
    auc = cross_val_score(pipe, X_raw, y_clf, cv=skf, scoring="roc_auc")
    f1  = cross_val_score(pipe, X_raw, y_clf, cv=skf, scoring="f1")
    res_clf.append({"modelo": nm, "Accuracy": round(acc.mean(),4),
                    "ROC_AUC": round(auc.mean(),4), "F1": round(f1.mean(),4)})
    print(f"  {nm:20s} Acc={acc.mean():.3f}  AUC={auc.mean():.3f}  F1={f1.mean():.3f}")

df_clf = pd.DataFrame(res_clf)
df_clf.to_csv(os.path.join(OUT_DIR, "modelo3_clasificador.csv"), index=False)

# ── SHAP VALUES ───────────────────────────────────────────────
print("\n=== SHAP VALUES ===")

# Entrenar RF en P-gp (mejor R²) para SHAP
scaler   = RobustScaler()
selector = SelectKBest(f_regression, k=50)
y_pgp    = df["dG_P-gp_MDR1"].values

X_scaled   = scaler.fit_transform(X_raw)
X_selected = selector.fit_transform(X_scaled, y_pgp)
feat_names = np.array(feature_cols)[selector.get_support()]

rf_shap = RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1)
rf_shap.fit(X_selected, y_pgp)

explainer   = shap.TreeExplainer(rf_shap)
shap_values = explainer.shap_values(X_selected)

# Top 20 features por importancia SHAP
shap_importance = pd.DataFrame({
    "feature": feat_names,
    "shap_mean_abs": np.abs(shap_values).mean(axis=0)
}).sort_values("shap_mean_abs", ascending=False).head(20)

shap_importance.to_csv(os.path.join(OUT_DIR, "shap_importancia_Pgp.csv"), index=False)
print("Top 10 features SHAP (P-gp):")
print(shap_importance.head(10).to_string(index=False))

# ── VISUALIZACIONES ───────────────────────────────────────────
print("\n=== GENERANDO VISUALIZACIONES ===")
plt.style.use('seaborn-v0_8-whitegrid')
COLORS = ['#2196F3','#4CAF50','#F44336','#FF9800','#9C27B0','#00BCD4']

# 1. Heatmap matriz de afinidad
fig, ax = plt.subplots(figsize=(14, 10))
heatmap_data = df.set_index("nombre")[DOCKING_COLS].sort_values("dG_FRalpha")
col_labels = ["P-gp\nMDR1", "CYP3A4", "TfR1", "FRα", "Lisozima", "HSA"]
sns.heatmap(heatmap_data, cmap="RdYlGn_r", ax=ax,
            xticklabels=col_labels, yticklabels=True,
            cbar_kws={"label": "ΔG (kcal/mol)"}, linewidths=0.3)
ax.set_title("Matriz de Afinidad de Docking\n(ΔG modo 1, kcal/mol)", fontsize=14, fontweight='bold')
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=7)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "heatmap_docking.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  ✓ heatmap_docking.png")

# 2. R² por modelo y proteína
fig, ax = plt.subplots(figsize=(12, 6))
prot_labels = ["P-gp MDR1", "CYP3A4", "TfR1", "FRα", "Lisozima", "HSA"]
x = np.arange(len(DOCKING_COLS))
width = 0.28
modelos_names = df_reg["modelo"].unique()
for i, mod in enumerate(modelos_names):
    vals = [df_reg[(df_reg["proteina"]==p) & (df_reg["modelo"]==mod)]["R2"].values[0]
            for p in DOCKING_COLS]
    errs = [df_reg[(df_reg["proteina"]==p) & (df_reg["modelo"]==mod)]["R2_std"].values[0]
            for p in DOCKING_COLS]
    ax.bar(x + i*width, vals, width, label=mod, color=COLORS[i], alpha=0.85,
           yerr=errs, capsize=3, error_kw={"elinewidth":1})

ax.axhline(0.7, color='red', linestyle='--', linewidth=1.5, label='Umbral R²=0.7')
ax.set_xticks(x + width)
ax.set_xticklabels(prot_labels, fontsize=11)
ax.set_ylabel("R² (5-fold CV)", fontsize=12)
ax.set_title("Rendimiento Modelos de Regresión por Proteína", fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.set_ylim(-0.1, 1.1)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "r2_por_proteina.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  ✓ r2_por_proteina.png")

# 3. ITI score barplot top 20
fig, ax = plt.subplots(figsize=(12, 8))
top20 = df_iti_sorted.head(20)
colors_bar = [PALETTE.get(g, "#888") for g in top20["grupo"]]
bars = ax.barh(range(len(top20)), top20["ITI_score"], color=colors_bar, alpha=0.85, edgecolor='white')
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20["nombre"], fontsize=9)
ax.set_xlabel("ITI Score (0-100)", fontsize=12)
ax.set_title("Top 20 Moléculas por Índice de Transporte (ITI)", fontsize=13, fontweight='bold')
ax.axvline(66, color='green', linestyle='--', linewidth=1.5, alpha=0.7, label='Favorable (>66)')
ax.axvline(33, color='orange', linestyle='--', linewidth=1.5, alpha=0.7, label='Moderado (>33)')
ax.invert_yaxis()

# Leyenda grupos
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=v, label=k.replace("_"," ").title())
                   for k, v in PALETTE.items()]
legend_elements.append(plt.Line2D([0],[0], color='green', linestyle='--', label='Favorable (>66)'))
ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "iti_top20.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  ✓ iti_top20.png")

# 4. SHAP beeswarm / bar
fig, ax = plt.subplots(figsize=(10, 7))
top_feats = shap_importance.head(15)
ax.barh(range(len(top_feats)), top_feats["shap_mean_abs"],
        color=COLORS[0], alpha=0.85, edgecolor='white')
ax.set_yticks(range(len(top_feats)))
ax.set_yticklabels(top_feats["feature"], fontsize=9)
ax.set_xlabel("SHAP mean |value|", fontsize=12)
ax.set_title("Top 15 Features por Importancia SHAP\n(Modelo P-gp MDR1, Random Forest)", fontsize=12, fontweight='bold')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "shap_importancia.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  ✓ shap_importancia.png")

# 5. Scatter docking por grupos
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()
for idx, prot in enumerate(DOCKING_COLS):
    ax = axes[idx]
    for grupo, color in PALETTE.items():
        mask = df["grupo"] == grupo
        ax.scatter(df[mask]["MW"], df[mask][prot],
                   c=color, label=grupo.replace("_"," ").title(),
                   alpha=0.8, s=60, edgecolors='white', linewidth=0.5)
    ax.set_xlabel("MW (Da)", fontsize=10)
    ax.set_ylabel("ΔG (kcal/mol)", fontsize=10)
    ax.set_title(prot.replace("dG_","").replace("_"," "), fontsize=11, fontweight='bold')

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=4, fontsize=9,
           bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Energía de Docking vs Peso Molecular por Grupo Molecular",
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "scatter_docking_MW.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  ✓ scatter_docking_MW.png")

# 6. Radar chart top 5 moléculas
from matplotlib.patches import FancyArrowPatch
categories = ["P-gp\n(evitar)", "CYP3A4\n(evitar)", "TfR1\n(targeting)",
              "FRα\n(targeting)", "Lisozima\n(biocompat)", "HSA\n(transporte)"]
N = len(categories)
angles = [n / float(N) * 2 * np.pi for n in range(N)]
angles += angles[:1]

fig, axes = plt.subplots(1, 5, figsize=(18, 4), subplot_kw=dict(polar=True))
top5 = df_iti_sorted.head(5)

for i, (_, row) in enumerate(top5.iterrows()):
    ax = axes[i]
    mol_data = df[df["nombre"] == row["nombre"]][DOCKING_COLS].values[0]
    # Normalizar: más negativo = más afinidad = valor más alto en radar
    vals_norm = []
    for j, col in enumerate(DOCKING_COLS):
        mn = df[col].min()
        mx = df[col].max()
        vals_norm.append((mol_data[j] - mn) / (mx - mn))
    vals_norm += vals_norm[:1]

    ax.plot(angles, vals_norm, 'o-', linewidth=2,
            color=PALETTE.get(row["grupo"], "#888"))
    ax.fill(angles, vals_norm, alpha=0.25,
            color=PALETTE.get(row["grupo"], "#888"))
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, size=7)
    ax.set_ylim(0, 1)
    ax.set_title(row["nombre"][:20], size=8, fontweight='bold', pad=10)
    ax.set_yticklabels([])

fig.suptitle("Perfil de Afinidad Top 5 Moléculas (Radar Chart)",
             fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "radar_top5.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  ✓ radar_top5.png")

# ── Resumen final ─────────────────────────────────────────────
print("\n" + "="*60)
print("RESUMEN FINAL")
print("="*60)
print("\nMejor modelo regresión por proteína:")
for prot in DOCKING_COLS:
    best = df_reg[df_reg["proteina"]==prot].sort_values("R2", ascending=False).iloc[0]
    print(f"  {prot:20s}: {best['modelo']:20s} R²={best['R2']:.3f}")

print(f"\nMejor clasificador:")
best_clf = df_clf.sort_values("ROC_AUC", ascending=False).iloc[0]
print(f"  {best_clf['modelo']} AUC={best_clf['ROC_AUC']:.3f} Acc={best_clf['Accuracy']:.3f}")

print(f"\nTop 5 ITI corregido:")
print(df_iti_sorted[["nombre","grupo","MW","ITI_score","perfil"]].head(5).to_string(index=False))

print(f"\n✓ Todos los resultados en: {OUT_DIR}")
