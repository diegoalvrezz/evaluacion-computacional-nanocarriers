"""
graficas_pro.py
===============
Regenera TODAS las gráficas del TFM con estilo publicación científica.
"""
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from matplotlib import rcParams
from pathlib import Path

# ── Estilo publicación ────────────────────────────────────────
rcParams.update({
    'font.family':       'DejaVu Sans',
    'font.size':         11,
    'axes.linewidth':    0.8,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.grid':         True,
    'axes.grid.axis':    'y',
    'grid.alpha':        0.3,
    'grid.linewidth':    0.5,
    'grid.color':        '#CCCCCC',
    'xtick.direction':   'out',
    'ytick.direction':   'out',
    'xtick.major.size':  4,
    'ytick.major.size':  4,
    'figure.dpi':        150,
    'savefig.dpi':       200,
    'savefig.bbox':      'tight',
    'savefig.facecolor': 'white',
    'legend.framealpha': 0.9,
    'legend.edgecolor':  '#CCCCCC',
    'legend.fontsize':   9,
})

# Paleta sofisticada
PALETTE = {
    "monomeros_degradacion": "#2E86AB",
    "oligomeros_unidades":   "#A23B72",
    "farmacos_modelo":       "#E84855",
    "ligandos_targeting":    "#F18F01",
}
LABELS = {
    "monomeros_degradacion": "Monómeros degradación",
    "oligomeros_unidades":   "Oligómeros/Unidades",
    "farmacos_modelo":       "Fármacos modelo",
    "ligandos_targeting":    "Ligandos targeting",
}
PROT_COLORS = {
    "P-gp MDR1":     "#E84855",
    "CYP3A4":        "#F18F01",
    "TfR1":          "#2E86AB",
    "FRalpha":       "#A23B72",
    "Lisozima":      "#3BB273",
    "HSA":           "#7B2D8B",
}
INT_COLORS = {
    "H-Bond":        "#2E86AB",
    "Hidrofóbica":   "#F18F01",
    "Puente salino": "#E84855",
}

BASE    = Path.home() / "Desktop/TFM_DiegoVallina/placeholder"
PLIP    = Path.home() / "Desktop/TFM_DiegoVallina/plip_resultados"
OUT     = BASE / "resultados_ML"

grupos  = list(PALETTE.keys())
labels_g = [LABELS[g] for g in grupos]
colors_g = [PALETTE[g] for g in grupos]

# ── Cargar datos ──────────────────────────────────────────────
df_ml    = pd.read_csv(BASE / "dataset_ML.csv")
df_iti   = pd.read_csv(OUT / "modelo2_ITI_corregido.csv")
df_reg   = pd.read_csv(OUT / "modelo1_regresion.csv")
df_admet = pd.read_csv(BASE / "admet_swissadme.csv")
df_plip  = pd.read_csv(PLIP / "plip_parsed.csv")
df_shap  = pd.read_csv(OUT / "shap_importancia_Pgp.csv")

DOCKING_COLS = ["dG_P-gp_MDR1","dG_CYP3A4","dG_TfR1","dG_FRalpha","dG_Lisozima","dG_HSA"]
PROT_LABELS  = ["P-gp MDR1","CYP3A4","TfR1","FRα","Lisozima","HSA"]

def save(fig, name):
    fig.savefig(OUT / name, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  ✓ {name}")

# ═══════════════════════════════════════════════════════════════
# 1. HEATMAP DOCKING
# ═══════════════════════════════════════════════════════════════
import seaborn as sns

fig, ax = plt.subplots(figsize=(13, 11))
hdata = df_ml.set_index("nombre")[DOCKING_COLS].sort_values("dG_FRalpha")
hdata.columns = PROT_LABELS
sns.heatmap(hdata, cmap="RdYlGn_r", ax=ax,
            xticklabels=True, yticklabels=True,
            cbar_kws={"label":"ΔG (kcal/mol)", "shrink":0.6},
            linewidths=0.2, linecolor='#F0F0F0')
ax.set_title("Matriz de Afinidad de Docking", fontsize=14, fontweight='bold', pad=15)
ax.tick_params(axis='y', labelsize=7)
ax.tick_params(axis='x', labelsize=10)
save(fig, "heatmap_docking.png")

# ═══════════════════════════════════════════════════════════════
# 2. R² POR PROTEÍNA
# ═══════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 6))
modelos = df_reg["modelo"].unique()
MOD_COLORS = {"Random Forest":"#2E86AB", "GradientBoosting":"#A23B72", "SVR":"#F18F01"}
x = np.arange(len(DOCKING_COLS))
w = 0.26
for i, mod in enumerate(modelos):
    vals = [df_reg[(df_reg["proteina"]==p)&(df_reg["modelo"]==mod)]["R2"].values[0] for p in DOCKING_COLS]
    errs = [df_reg[(df_reg["proteina"]==p)&(df_reg["modelo"]==mod)]["R2_std"].values[0] for p in DOCKING_COLS]
    ax.bar(x+i*w, vals, w, label=mod, color=MOD_COLORS.get(mod,"#888"),
           alpha=0.85, yerr=errs, capsize=3, error_kw={"elinewidth":1,"ecolor":"#555"})
ax.axhline(0.7, color='#E84855', linestyle='--', linewidth=1.5, label='Umbral R²=0.7', zorder=5)
ax.set_xticks(x+w); ax.set_xticklabels(PROT_LABELS, fontsize=11)
ax.set_ylabel("R² (validación cruzada 5-fold)", fontsize=11)
ax.set_title("Rendimiento de los Modelos de Regresión por Proteína", fontsize=13, fontweight='bold')
ax.legend(fontsize=10, loc='upper right')
ax.set_ylim(-0.05, 1.05)
save(fig, "r2_por_proteina.png")

# ═══════════════════════════════════════════════════════════════
# 3. ITI TOP 20
# ═══════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(11, 8))
top20 = df_iti.sort_values("ITI_score", ascending=False).head(20)
colors_bar = [PALETTE.get(g,"#888") for g in top20["grupo"]]
ax.barh(range(len(top20)), top20["ITI_score"], color=colors_bar, alpha=0.85,
        edgecolor='white', linewidth=0.5)
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20["nombre"], fontsize=9)
ax.set_xlabel("ITI Score (0–100)", fontsize=11)
ax.set_title("Top 20 Moléculas — Índice de Eficiencia de Transporte (ITI)", fontsize=13, fontweight='bold')
ax.axvline(66, color='#3BB273', linestyle='--', linewidth=1.5, alpha=0.8, label='Umbral favorable (>66)')
ax.axvline(33, color='#F18F01', linestyle='--', linewidth=1.5, alpha=0.8, label='Umbral moderado (>33)')
ax.invert_yaxis()
from matplotlib.patches import Patch
legend_el = [Patch(facecolor=v, label=LABELS[k]) for k,v in PALETTE.items()
             if k in top20["grupo"].values]
legend_el += [plt.Line2D([0],[0],color='#3BB273',linestyle='--',label='Favorable (>66)'),
              plt.Line2D([0],[0],color='#F18F01',linestyle='--',label='Moderado (>33)')]
ax.legend(handles=legend_el, loc='lower right', fontsize=8)
save(fig, "iti_top20.png")

# ═══════════════════════════════════════════════════════════════
# 4. SHAP SUMMARY
# ═══════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(9, 7))
top15 = df_shap.head(15).sort_values("shap_mean_abs")
cmap_vals = plt.cm.Blues(np.linspace(0.35, 0.9, len(top15)))
ax.barh(range(len(top15)), top15["shap_mean_abs"], color=cmap_vals, edgecolor='white')
ax.set_yticks(range(len(top15)))
ax.set_yticklabels(top15["feature"], fontsize=9)
ax.set_xlabel("SHAP mean |value|", fontsize=11)
ax.set_title("Importancia de Features — Análisis SHAP\n(Modelo P-gp MDR1, Random Forest)", fontsize=12, fontweight='bold')
save(fig, "shap_importancia.png")

# ═══════════════════════════════════════════════════════════════
# 5. SCATTER DOCKING VS MW
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()
for idx, (prot, label) in enumerate(zip(DOCKING_COLS, PROT_LABELS)):
    ax = axes[idx]
    for grupo, color in PALETTE.items():
        mask = df_ml["grupo"] == grupo
        ax.scatter(df_ml[mask]["MW"], df_ml[mask][prot],
                   c=color, alpha=0.8, s=55, edgecolors='white', linewidth=0.4, zorder=3)
    ax.set_xlabel("MW (Da)", fontsize=9)
    ax.set_ylabel("ΔG (kcal/mol)", fontsize=9)
    ax.set_title(label, fontsize=11, fontweight='bold')
from matplotlib.patches import Patch
handles = [Patch(facecolor=v, label=LABELS[k]) for k,v in PALETTE.items()]
fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=9, bbox_to_anchor=(0.5,-0.02))
fig.suptitle("Energía de Docking vs Peso Molecular", fontsize=13, fontweight='bold')
plt.tight_layout()
save(fig, "scatter_docking_MW.png")

# ═══════════════════════════════════════════════════════════════
# 6. DISTRIBUCIÓN ITI
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
lip_data = [df_iti[df_iti["grupo"]==g]["ITI_score"].dropna().values for g in grupos]
bp = axes[0].boxplot(lip_data, patch_artist=True, medianprops=dict(color='#222',linewidth=2),
                     whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1.5))
for patch, color in zip(bp['boxes'], colors_g):
    patch.set_facecolor(color); patch.set_alpha(0.7)
axes[0].set_xticklabels([LABELS[g].replace(" ","\n") for g in grupos], fontsize=9)
axes[0].set_ylabel("ITI Score", fontsize=11)
axes[0].set_title("Distribución ITI por Grupo Molecular", fontsize=12, fontweight='bold')
axes[0].axhline(66, color='#3BB273', linestyle='--', alpha=0.7)
axes[0].axhline(33, color='#F18F01', linestyle='--', alpha=0.7)

for i, (grupo, color) in enumerate(zip(grupos, colors_g)):
    data = df_iti[df_iti["grupo"]==grupo]["ITI_score"].values
    parts = axes[1].violinplot(data, positions=[i], showmedians=True, showextrema=True)
    for pc in parts['bodies']:
        pc.set_facecolor(color); pc.set_alpha(0.6)
    parts['cmedians'].set_color('#222')
    jitter = np.random.normal(i, 0.05, len(data))
    axes[1].scatter(jitter, data, c=color, alpha=0.75, s=35, zorder=3, edgecolors='white', linewidth=0.3)
axes[1].set_xticks(range(len(grupos)))
axes[1].set_xticklabels([LABELS[g].replace(" ","\n") for g in grupos], fontsize=9)
axes[1].set_ylabel("ITI Score", fontsize=11)
axes[1].set_title("Violin Plot ITI por Grupo Molecular", fontsize=12, fontweight='bold')
axes[1].axhline(66, color='#3BB273', linestyle='--', alpha=0.7)
axes[1].axhline(33, color='#F18F01', linestyle='--', alpha=0.7)
plt.tight_layout()
save(fig, "distribucion_ITI_grupos.png")

# ═══════════════════════════════════════════════════════════════
# 7. CORRELACIÓN PROTEÍNAS
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
corr = df_ml[DOCKING_COLS].corr()
corr.index = corr.columns = PROT_LABELS
mask = np.triu(np.ones_like(corr, dtype=bool))
sns.heatmap(corr, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r',
            center=0, ax=axes[0], square=True,
            cbar_kws={"label":"r Pearson","shrink":0.8},
            annot_kws={"size":10}, linewidths=0.3)
axes[0].set_title("Correlación entre Energías de Docking", fontsize=12, fontweight='bold')

colors_s = [PALETTE.get(g,"#888") for g in df_ml["grupo"]]
axes[1].scatter(df_ml["dG_P-gp_MDR1"], df_ml["dG_TfR1"],
                c=colors_s, alpha=0.8, s=65, edgecolors='white', linewidth=0.4)
z = np.polyfit(df_ml["dG_P-gp_MDR1"].dropna(), df_ml["dG_TfR1"].dropna(), 1)
xl = np.linspace(df_ml["dG_P-gp_MDR1"].min(), df_ml["dG_P-gp_MDR1"].max(), 100)
axes[1].plot(xl, np.poly1d(z)(xl), 'k--', linewidth=1.5, alpha=0.6)
r = df_ml[["dG_P-gp_MDR1","dG_TfR1"]].corr().iloc[0,1]
axes[1].set_xlabel("ΔG P-gp MDR1 (kcal/mol)", fontsize=11)
axes[1].set_ylabel("ΔG TfR1 (kcal/mol)", fontsize=11)
axes[1].set_title(f"P-gp MDR1 vs TfR1 (r = {r:.3f})", fontsize=12, fontweight='bold')
handles = [Patch(facecolor=v, label=LABELS[k]) for k,v in PALETTE.items()]
axes[1].legend(handles=handles, fontsize=8, loc='lower right')
plt.tight_layout()
save(fig, "correlacion_proteinas.png")

# ═══════════════════════════════════════════════════════════════
# 8. PLIP INTERACCIONES
# ═══════════════════════════════════════════════════════════════
PROTS_PLIP = ["TfR1 (1CX8)","FRalpha (4LRH)","P-gp (7A65)"]
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

tipo_count = df_plip.groupby(["proteina","tipo"]).size().reset_index(name="n")
tipos = list(INT_COLORS.keys())
bottom = np.zeros(len(PROTS_PLIP))
for tipo in tipos:
    vals = [tipo_count[(tipo_count["proteina"]==p)&(tipo_count["tipo"]==tipo)]["n"].sum()
            for p in PROTS_PLIP]
    axes[0].bar(PROTS_PLIP, vals, bottom=bottom, label=tipo,
                color=INT_COLORS[tipo], alpha=0.85, edgecolor='white')
    bottom += np.array(vals, dtype=float)
axes[0].set_ylabel("Número de interacciones", fontsize=11)
axes[0].set_title("Interacciones proteína-ligando por tipo\n(Top 5 candidatos ITI, PLIP)", fontsize=12, fontweight='bold')
axes[0].legend(fontsize=10)
axes[0].tick_params(axis='x', rotation=10)

top_res = df_plip.groupby(["proteina","residuo","tipo"]).size().reset_index(name="n")
prot_plot = "FRalpha (4LRH)"
sub = top_res[top_res["proteina"]==prot_plot].sort_values("n",ascending=False).head(12)
colors_bar2 = [INT_COLORS.get(t,"#888") for t in sub["tipo"]]
axes[1].barh(range(len(sub)), sub["n"], color=colors_bar2, alpha=0.85, edgecolor='white')
axes[1].set_yticks(range(len(sub)))
axes[1].set_yticklabels(sub["residuo"], fontsize=9)
axes[1].invert_yaxis()
axes[1].set_xlabel("Nº interacciones", fontsize=11)
axes[1].set_title(f"Top residuos de interacción\n{prot_plot}", fontsize=12, fontweight='bold')
from matplotlib.patches import Patch
leg2 = [Patch(facecolor=v, label=k) for k,v in INT_COLORS.items()]
axes[1].legend(handles=leg2, fontsize=9)
plt.tight_layout()
save(fig, "plip_interacciones_pro.png")

# ═══════════════════════════════════════════════════════════════
# 9. ADMET
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Lipinski
lip_data = [df_admet[df_admet["grupo"]==g]["Lipinski #violations"].dropna().values for g in grupos]
bp = axes[0].boxplot(lip_data, patch_artist=True, medianprops=dict(color='#222',linewidth=2))
for patch, color in zip(bp['boxes'], colors_g):
    patch.set_facecolor(color); patch.set_alpha(0.7)
axes[0].set_xticklabels([LABELS[g].replace(" ","\n") for g in grupos], fontsize=8)
axes[0].set_ylabel("Nº violaciones", fontsize=11)
axes[0].set_title("Regla de Lipinski\npor grupo molecular", fontsize=11, fontweight='bold')
axes[0].axhline(0, color='#3BB273', linestyle='--', alpha=0.6, label='Drug-like')
axes[0].legend(fontsize=9)

# P-gp sustrato
pgp_pct = [(df_admet[df_admet["grupo"]==g]["Pgp substrate"]=="Yes").mean()*100 for g in grupos]
bars = axes[1].bar(labels_g, pgp_pct, color=colors_g, alpha=0.85, edgecolor='white')
axes[1].set_ylabel("% moléculas", fontsize=11)
axes[1].set_title("% Sustrato P-gp\npor grupo molecular", fontsize=11, fontweight='bold')
axes[1].set_ylim(0, 105)
axes[1].tick_params(axis='x', rotation=15)
for bar, val in zip(bars, pgp_pct):
    axes[1].text(bar.get_x()+bar.get_width()/2, bar.get_height()+1,
                f"{val:.0f}%", ha='center', fontsize=10, fontweight='bold')

# Bioavailability scatter
for grupo, color in PALETTE.items():
    sub = df_admet[df_admet["grupo"]==grupo]
    axes[2].scatter(sub["Consensus Log P"], sub["Bioavailability Score"],
                    c=color, label=LABELS[grupo], alpha=0.8, s=65,
                    edgecolors='white', linewidth=0.4)
axes[2].axvline(5, color='#E84855', linestyle='--', linewidth=1.5, alpha=0.7, label='LogP=5')
axes[2].axhline(0.55, color='#3BB273', linestyle='--', linewidth=1.5, alpha=0.7, label='BA=0.55')
axes[2].set_xlabel("Consensus Log P", fontsize=11)
axes[2].set_ylabel("Bioavailability Score", fontsize=11)
axes[2].set_title("Biodisponibilidad vs Lipofilia", fontsize=11, fontweight='bold')
axes[2].legend(fontsize=7, loc='lower right')

plt.tight_layout()
save(fig, "admet_completo.png")

print("\n✓ Todas las gráficas regeneradas con estilo publicación")
