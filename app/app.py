# ============================================================
# app.py  -  Drug Delivery Carrier Recommender
# TFM · Diego Vallina Álvarez · Universidad de Burgos · 2025-2026
# ============================================================

import warnings
warnings.filterwarnings('ignore')

import streamlit as st
import pandas as pd
import numpy as np
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

# Draw se importa de forma segura para entornos sin display
try:
    from rdkit.Chem import Draw
    DRAW_AVAILABLE = True
except ImportError:
    DRAW_AVAILABLE = False

import py3Dmol
import streamlit.components.v1 as components

# ── Rutas ─────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).parent
MODELS_DIR = BASE_DIR / "resultados_ML" / "models"
DATA_DIR   = BASE_DIR

# ── Constantes ────────────────────────────────────────────────
PROT_NAMES = ["P-gp MDR1", "CYP3A4", "TfR1", "FRα", "Lisozima", "HSA"]
PROT_KEYS  = ["dG_P_gp_MDR1", "dG_CYP3A4", "dG_TfR1",
              "dG_FRalpha", "dG_Lisozima", "dG_HSA"]
PROT_ROLES = {
    "P-gp MDR1": "eflujo",
    "CYP3A4":    "metabolismo",
    "TfR1":      "targeting",
    "FRα":       "targeting",
    "Lisozima":  "biocompat.",
    "HSA":       "transporte"
}

ITI_WEIGHTS = {
    "dG_P_gp_MDR1": -1.5,
    "dG_CYP3A4":    -1.0,
    "dG_TfR1":       1.2,
    "dG_FRalpha":    1.2,
    "dG_Lisozima":   0.8,
    "dG_HSA":        0.6,
}

GRUPO_LABELS = {
    "monomeros_degradacion": "Monomeros degradacion",
    "oligomeros_unidades":   "Oligomeros/Unidades",
    "farmacos_modelo":       "Farmacos modelo",
    "ligandos_targeting":    "Ligandos targeting"
}

GRUPO_COLORS = {
    "Monomeros degradacion": "#2E86AB",
    "Oligomeros/Unidades":   "#A23B72",
    "Farmacos modelo":       "#E84855",
    "Ligandos targeting":    "#F18F01"
}

# ── Carga de modelos ──────────────────────────────────────────
@st.cache_resource
def load_models():
    models = {}
    for key in PROT_KEYS:
        path = MODELS_DIR / f"rf_{key}.joblib"
        if path.exists():
            models[key] = joblib.load(path)
    models["classifier"] = joblib.load(MODELS_DIR / "gb_classifier.joblib")
    models["features"]   = joblib.load(MODELS_DIR / "feature_cols.joblib")
    return models

@st.cache_data
def load_dataset():
    df     = pd.read_csv(DATA_DIR / "dataset_ML.csv")
    df_iti = pd.read_csv(DATA_DIR / "resultados_ML" / "modelo2_ITI_corregido.csv")
    df     = df.merge(df_iti[["nombre", "ITI_score"]], on="nombre", how="left")
    df["grupo_label"] = df["grupo"].map(GRUPO_LABELS).fillna(df["grupo"])
    return df

@st.cache_data
def load_base_molecular():
    return pd.read_csv(DATA_DIR / "base_molecular_pubchem.csv")

# ── Descriptores ──────────────────────────────────────────────
def calc_descriptors_2d(mol):
    desc = {}
    desc["MW"]             = Descriptors.MolWt(mol)
    desc["ExactMW"]        = Descriptors.ExactMolWt(mol)
    desc["LogP"]           = Descriptors.MolLogP(mol)
    desc["TPSA"]           = Descriptors.TPSA(mol)
    desc["HBD"]            = rdMolDescriptors.CalcNumHBD(mol)
    desc["HBA"]            = rdMolDescriptors.CalcNumHBA(mol)
    desc["RotBonds"]       = rdMolDescriptors.CalcNumRotatableBonds(mol)
    desc["AromaticRings"]  = rdMolDescriptors.CalcNumAromaticRings(mol)
    desc["HeavyAtomCount"] = mol.GetNumHeavyAtoms()
    desc["NumAtoms"]       = mol.GetNumAtoms()
    desc["NumBonds"]       = mol.GetNumBonds()
    desc["BertzCT"]        = Descriptors.BertzCT(mol)
    desc["Chi0"]           = Descriptors.Chi0(mol)
    desc["Chi1"]           = Descriptors.Chi1(mol)
    desc["Kappa1"]         = Descriptors.Kappa1(mol)
    desc["Kappa2"]         = Descriptors.Kappa2(mol)
    desc["FractionCSP3"]   = rdMolDescriptors.CalcFractionCSP3(mol)
    desc["MolMR"]          = Descriptors.MolMR(mol)
    return desc

def calc_fingerprints(mol):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp  = gen.GetFingerprintAsNumPy(mol)
    return {f"Morgan_{i}": int(fp[i]) for i in range(2048)}

def build_feature_vector(mol, feature_cols):
    d2        = calc_descriptors_2d(mol)
    fps       = calc_fingerprints(mol)
    all_feats = {**d2, **fps}
    vec       = np.array([all_feats.get(f, 0.0) for f in feature_cols], dtype=float)
    return vec

def calc_esol_solubility(mol):
    # Modelo ESOL (Delaney 2004) para estimar solubilidad acuosa
    mw      = Descriptors.MolWt(mol)
    logp    = Descriptors.MolLogP(mol)
    rb      = rdMolDescriptors.CalcNumRotatableBonds(mol)
    ap      = sum(1 for a in mol.GetAromaticAtoms())
    # Formula ESOL: log(S) = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RB - 0.74*AP
    log_s   = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ap
    sol_mol = 10 ** log_s
    sol_mg  = sol_mol * mw
    if log_s > -1:
        clase = "Muy soluble"
    elif log_s > -2:
        clase = "Soluble"
    elif log_s > -4:
        clase = "Moderadamente soluble"
    elif log_s > -6:
        clase = "Poco soluble"
    else:
        clase = "Insoluble"
    return round(log_s, 2), round(sol_mg, 3), clase

# ── ITI y clasificación ───────────────────────────────────────
def calc_iti(dg_dict):
    vals  = [dg_dict[k] for k in ITI_WEIGHTS if k in dg_dict]
    if not vals:
        return 50.0
    mn, mx = min(vals), max(vals)
    rng    = mx - mn if mx != mn else 1.0
    score, total_w = 0.0, 0.0
    for key, w in ITI_WEIGHTS.items():
        if key in dg_dict:
            norm    = (dg_dict[key] - mn) / rng
            score  += w * norm
            total_w += abs(w)
    raw = score / total_w if total_w > 0 else 0.0
    iti = (raw + 1) / 2 * 100
    return round(max(0.0, min(100.0, iti)), 1)

def get_recommendation(iti, dg_dict):
    if iti >= 66:
        label = "Favorable"
        color = "#3BB273"
        if dg_dict.get("dG_TfR1", 0) < -5 or dg_dict.get("dG_FRalpha", 0) < -5:
            rec = "Candidato para ligando de targeting activo (TfR1/FRα)."
        elif dg_dict.get("dG_P_gp_MDR1", 0) > -3:
            rec = "Buen componente de recubrimiento PEGilado. Baja afinidad P-gp."
        else:
            rec = "Candidato favorable como componente de nanocarrier polimerico."
    elif iti >= 33:
        label = "Moderado"
        color = "#F18F01"
        rec   = "Perfil moderado. Evaluar con ensayos adicionales antes de uso en DDS."
    else:
        label = "Desfavorable"
        color = "#E84855"
        rec   = "Alta afinidad por P-gp o CYP3A4. No recomendado como componente de DDS."
    return label, color, rec

def check_lipinski(desc):
    rules = {
        "MW <= 500":  desc["MW"] <= 500,
        "LogP <= 5":  desc["LogP"] <= 5,
        "HBD <= 5":   desc["HBD"] <= 5,
        "HBA <= 10":  desc["HBA"] <= 10
    }
    violations = sum(1 for v in rules.values() if not v)
    return rules, violations

# ── Similitud Tanimoto ────────────────────────────────────────
def find_most_similar(mol, df_base):
    # Calcula similitud Tanimoto entre la molécula y todas las del dataset
    gen    = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp_mol = gen.GetFingerprint(mol)
    best_sim, best_name, best_smi = 0.0, "", ""
    for _, row in df_base.iterrows():
        ref = Chem.MolFromSmiles(row["smiles_pubchem"])
        if ref is None:
            continue
        fp_ref = gen.GetFingerprint(ref)
        sim    = DataStructs.TanimotoSimilarity(fp_mol, fp_ref)
        if sim > best_sim:
            best_sim  = sim
            best_name = row["nombre_entrada"]
            best_smi  = row["smiles_pubchem"]
    return best_name, best_smi, round(best_sim, 3)

# ── Visualización 3D ──────────────────────────────────────────
def show_3d_molecule(mol, style="stick"):
    # Genera coordenadas 3D y muestra la molécula con py3Dmol via HTML
    try:
        mol_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol_h)
        mol_block = Chem.MolToMolBlock(mol_h)

        viewer = py3Dmol.view(width=400, height=350)
        viewer.addModel(mol_block, "mol")
        if style == "stick":
            viewer.setStyle({"stick": {}})
        elif style == "sphere":
            viewer.setStyle({"sphere": {"scale": 0.3}, "stick": {"radius": 0.1}})
        elif style == "surface":
            viewer.setStyle({"stick": {}})
            viewer.addSurface(py3Dmol.VDW, {"opacity": 0.7})
        viewer.setBackgroundColor("white")
        viewer.zoomTo()
        html_str = viewer._make_html()
        components.html(html_str, height=350, width=400)
    except Exception as e:
        st.info(f"No se pudo generar la estructura 3D: {e}")

# ── Gráficas ──────────────────────────────────────────────────
def make_radar(dg_dict, mol_name):
    values    = [dg_dict.get(k, 0) for k in PROT_KEYS]
    mn, mx    = min(values), max(values)
    rng       = mx - mn if mx != mn else 1.0
    norm_vals = [(v - mn) / rng for v in values]
    norm_vals += norm_vals[:1]
    N         = len(PROT_NAMES)
    angles    = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles   += angles[:1]

    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(polar=True))
    ax.plot(angles, norm_vals, 'o-', linewidth=2.5, color="#2E86AB")
    ax.fill(angles, norm_vals, alpha=0.25, color="#2E86AB")
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(
        [f"{n}\n({PROT_ROLES[n]})" for n in PROT_NAMES], size=8, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75])
    ax.set_yticklabels(["0.25", "0.5", "0.75"], size=7, color='gray')
    ax.grid(color='gray', alpha=0.3)
    ax.set_title(f"Perfil de afinidad\n{mol_name[:30]}", size=10,
                 fontweight='bold', pad=20)
    plt.tight_layout()
    return fig

def make_comparison_chart(molecules_data):
    n_mols  = len(molecules_data)
    x       = np.arange(len(PROT_NAMES))
    width   = 0.8 / n_mols
    colors  = ["#2E86AB", "#A23B72", "#E84855", "#F18F01"]

    fig, ax = plt.subplots(figsize=(12, 5))
    for i, (name, dg_dict) in enumerate(molecules_data.items()):
        vals = [dg_dict.get(k, 0) for k in PROT_KEYS]
        ax.bar(x + i * width, vals, width, label=name[:20],
               color=colors[i % len(colors)], alpha=0.85, edgecolor='white')
    ax.set_xticks(x + width * (n_mols - 1) / 2)
    ax.set_xticklabels(PROT_NAMES, fontsize=10)
    ax.set_ylabel("ΔG predicho (kcal/mol)", fontsize=11)
    ax.set_title("Comparacion de perfiles de afinidad", fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    return fig

def make_session_heatmap(historial):
    # Heatmap de ΔG para todas las moléculas del historial
    nombres = [e["nombre"] for e in historial]
    matrix  = []
    for entry in historial:
        fila = [entry["dg_dict"].get(k, 0) for k in PROT_KEYS]
        matrix.append(fila)
    df_heat = pd.DataFrame(matrix, index=nombres, columns=PROT_NAMES)

    fig, ax = plt.subplots(figsize=(10, max(3, len(historial) * 0.5)))
    sns.heatmap(df_heat, cmap="RdYlGn_r", ax=ax, annot=True, fmt=".2f",
                linewidths=0.3, cbar_kws={"label": "ΔG (kcal/mol)", "shrink": 0.6})
    ax.set_title("Matriz de afinidad - sesion actual", fontsize=12, fontweight='bold')
    ax.set_ylabel("")
    plt.tight_layout()
    return fig

def make_iti_bar(iti):
    # Barra de progreso coloreada para el ITI
    if iti >= 66:
        color = "#3BB273"
    elif iti >= 33:
        color = "#F18F01"
    else:
        color = "#E84855"
    fig, ax = plt.subplots(figsize=(6, 0.6))
    ax.barh([0], [iti], color=color, height=0.5)
    ax.barh([0], [100], color="#EEEEEE", height=0.5, zorder=0)
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.axvline(66, color='#3BB273', linestyle='--', linewidth=1.2, alpha=0.7)
    ax.axvline(33, color='#F18F01', linestyle='--', linewidth=1.2, alpha=0.7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel("ITI Score", fontsize=9)
    plt.tight_layout()
    return fig

def results_to_csv(name, dg_dict, iti, label):
    rows = []
    for key, prot in zip(PROT_KEYS, PROT_NAMES):
        rows.append({
            "Molecula":      name,
            "Proteina":      prot,
            "DeltaG_pred":   dg_dict.get(key, ""),
            "ITI":           iti,
            "Clasificacion": label
        })
    buf = BytesIO()
    pd.DataFrame(rows).to_csv(buf, index=False)
    return buf.getvalue()

# ── Configuración de la página ────────────────────────────────
st.set_page_config(
    page_title="Drug Delivery Carrier Recommender",
    layout="wide"
)

# CSS personalizado para estilo cientifico
st.markdown("""
<style>
/* Fuente principal */
html, body, [class*="css"] {
    font-family: 'Inter', 'Helvetica Neue', Arial, sans-serif;
}

/* Cabecera principal */
h1 {
    font-size: 1.8rem;
    font-weight: 700;
    letter-spacing: -0.5px;
    color: #1a1a2e;
    border-bottom: 2px solid #2E86AB;
    padding-bottom: 0.4rem;
    margin-bottom: 0.5rem;
}

/* Subtítulos */
h2, h3 {
    font-weight: 600;
    color: #1a1a2e;
    letter-spacing: -0.3px;
}

/* Métricas */
[data-testid="metric-container"] {
    background: #f8f9fb;
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    padding: 12px 16px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04);
}
[data-testid="metric-container"] label {
    font-size: 0.75rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    color: #64748b;
}
[data-testid="metric-container"] [data-testid="stMetricValue"] {
    font-size: 1.5rem;
    font-weight: 700;
    color: #1a1a2e;
}

/* Botón primario */
[data-testid="baseButton-primary"] {
    background: #2E86AB;
    border: none;
    border-radius: 6px;
    font-weight: 600;
    letter-spacing: 0.3px;
}
[data-testid="baseButton-primary"]:hover {
    background: #1a6a8a;
}

/* Botones secundarios */
[data-testid="baseButton-secondary"] {
    border: 1px solid #cbd5e1;
    border-radius: 6px;
    color: #334155;
    font-weight: 500;
}

/* Pestañas */
[data-baseweb="tab-list"] {
    border-bottom: 2px solid #e2e8f0;
    gap: 4px;
}
[data-baseweb="tab"] {
    font-weight: 500;
    font-size: 0.9rem;
    color: #64748b;
    padding: 8px 16px;
}
[aria-selected="true"] {
    color: #2E86AB;
    border-bottom: 2px solid #2E86AB;
    font-weight: 600;
}

/* Cajas de resultados */
.resultado-box {
    background: #f8f9fb;
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    padding: 16px;
    margin: 8px 0;
}

/* Inputs */
[data-baseweb="input"] {
    border-radius: 6px;
    border: 1px solid #cbd5e1;
    font-size: 0.9rem;
}

/* Tablas */
[data-testid="stDataFrame"] {
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    overflow: hidden;
}

/* Sidebar y expanders */
[data-testid="stExpander"] {
    border: 1px solid #e2e8f0;
    border-radius: 8px;
}

/* Footer */
footer { visibility: hidden; }
</style>
""", unsafe_allow_html=True)

models  = load_models()
df_data = load_dataset()
df_base = load_base_molecular()

if "historial" not in st.session_state:
    st.session_state["historial"] = []
if "smiles_input" not in st.session_state:
    st.session_state["smiles_input"] = ""
if "mol_name_input" not in st.session_state:
    st.session_state["mol_name_input"] = "Molecula nueva"

# ── Cabecera ──────────────────────────────────────────────────
st.markdown("""
<div style="margin-bottom: 1rem;">
    <h1>Drug Delivery Carrier Recommender</h1>
    <p style="color: #64748b; font-size: 0.95rem; margin-top: 0.3rem;">
        Plataforma computacional para la evaluacion del perfil de afinidad
        de moleculas frente a proteinas de barrera biologica.
        Basada en modelos Random Forest entrenados con DFT, docking molecular
        y 2908 descriptores moleculares.
        <span style="color: #2E86AB; font-weight: 600;">
        TFM · Universidad de Burgos · 2025-2026
        </span>
    </p>
</div>
""", unsafe_allow_html=True)

tab1, tab2, tab3, tab4 = st.tabs([
    "Analizar molecula",
    "Ranking ITI del dataset",
    "Comparador",
    "Historial de sesion"
])

# ════════════════════════════════════════════════════════════════
# PESTANA 1 — Analizar molecula
# ════════════════════════════════════════════════════════════════
with tab1:
    st.subheader("Analizar una molecula nueva")

    input_method = st.radio(
        "Introduce la molecula:",
        ["Escribir SMILES", "Seleccionar del dataset del TFM"],
        horizontal=True
    )

    smiles_input   = st.session_state["smiles_input"]
    mol_name_input = st.session_state["mol_name_input"]

    if input_method == "Escribir SMILES":
        col_inp, col_ex = st.columns([2, 1])
        with col_inp:
            smiles_input   = st.text_input("SMILES",
                value=st.session_state["smiles_input"],
                placeholder="Ej: C(COCCOCCO)O  (Triethylene glycol)")
            mol_name_input = st.text_input("Nombre (opcional)",
                value=st.session_state["mol_name_input"])
            st.session_state["smiles_input"]   = smiles_input
            st.session_state["mol_name_input"] = mol_name_input
        with col_ex:
            st.markdown("**Ejemplos rapidos:**")
            ejemplos = {
                "Triethylene glycol": "C(COCCOCCO)O",
                "Acido folico":       "C1=CC(=CC=C1C(=O)NC(CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N",
                "Doxorrubicina":      "CC1C(C(CC(O1)OC2CC(CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)(C(=O)CO)O)N)O",
                "Chitobiose":         "C(C1C(C(C(C(O1)OC(C(CO)O)C(C(C=O)N)O)N)O)O)O",
            }
            for name, smi in ejemplos.items():
                if st.button(name, use_container_width=True):
                    st.session_state["smiles_input"]   = smi
                    st.session_state["mol_name_input"] = name
                    st.rerun()
    else:
        grupos_disp = ["Todos"] + list(df_data["grupo_label"].dropna().unique())
        grupo_sel   = st.selectbox("Filtrar por grupo:", grupos_disp)
        df_filtrado = df_data if grupo_sel == "Todos" else df_data[df_data["grupo_label"] == grupo_sel]
        mol_sel     = st.selectbox("Selecciona una molecula:", df_filtrado["nombre"].tolist())
        fila_base   = df_base[df_base["nombre_entrada"] == mol_sel]
        if not fila_base.empty:
            if smiles_input != fila_base["smiles_pubchem"].values[0]:
                st.session_state["smiles_input"]   = fila_base["smiles_pubchem"].values[0]
                st.session_state["mol_name_input"] = mol_sel
                smiles_input   = st.session_state["smiles_input"]
                mol_name_input = st.session_state["mol_name_input"]
            st.code(smiles_input, language=None)

    # Vista previa de la molécula
    if smiles_input:
        mol_preview = Chem.MolFromSmiles(smiles_input)
        if mol_preview:
            desc_prev   = calc_descriptors_2d(mol_preview)
            rules, viols = check_lipinski(desc_prev)
            log_s, sol_mg, sol_clase = calc_esol_solubility(mol_preview)

            col_2d, col_3d, col_props = st.columns([1, 1, 2])

            with col_2d:
                st.markdown("**Estructura 2D**")
                if DRAW_AVAILABLE:
                    img = Draw.MolToImage(mol_preview, size=(220, 180))
                    st.image(img)
                else:
                    st.info("Visualización 2D no disponible en este entorno.")

            with col_3d:
                st.markdown("**Estructura 3D interactiva**")
                style_3d = st.selectbox("Estilo:",
                    ["stick", "sphere", "surface"], key="style_preview")
                try:
                    show_3d_molecule(mol_preview, style=style_3d)
                except Exception:
                    st.info("No se pudo generar la estructura 3D.")

            with col_props:
                st.markdown("**Propiedades**")
                c1, c2, c3 = st.columns(3)
                c1.metric("MW (Da)",  f"{desc_prev['MW']:.1f}")
                c2.metric("LogP",     f"{desc_prev['LogP']:.2f}")
                c3.metric("TPSA (Å²)",f"{desc_prev['TPSA']:.1f}")
                c1.metric("HBD", int(desc_prev['HBD']))
                c2.metric("HBA", int(desc_prev['HBA']))
                c3.metric("Rot. bonds", int(desc_prev['RotBonds']))

                # Solubilidad ESOL
                sol_color = "normal" if sol_clase in ["Muy soluble","Soluble"] else "off"
                st.metric("Solubilidad estimada (ESOL)",
                          f"{sol_mg:.2f} mg/mL",
                          delta=sol_clase,
                          delta_color=sol_color)
                st.caption(f"log S = {log_s:.2f}")

                # Lipinski
                lip_str = "Cumple Lipinski" if viols == 0 else f"{viols} violacion(es) Lipinski"
                for rule, ok in rules.items():
                    st.markdown(f"{'OK' if ok else 'X'}  {rule}")

            # Molécula más similar del dataset
            sim_name, sim_smi, sim_score = find_most_similar(mol_preview, df_base)
            st.markdown(f"**Molecula mas similar en el dataset:** {sim_name} "
                        f"(Tanimoto = {sim_score:.3f})")
            if sim_smi:
                mol_sim = Chem.MolFromSmiles(sim_smi)
                if mol_sim:
                    if DRAW_AVAILABLE:
                        img_sim = Draw.MolToImage(mol_sim, size=(180, 140))
                        st.image(img_sim, caption=sim_name)
        else:
            st.error("SMILES no valido.")
            smiles_input = ""

    if smiles_input and st.button("Analizar", type="primary", use_container_width=True):
        mol = Chem.MolFromSmiles(smiles_input)
        if mol is None:
            st.error("SMILES no valido.")
            st.stop()

        with st.spinner("Calculando descriptores y predicciones..."):
            feat_vec = build_feature_vector(mol, models["features"])
            X_input  = feat_vec.reshape(1, -1)
            dg_preds = {}
            for key in PROT_KEYS:
                if key in models:
                    dg_preds[key] = round(float(models[key].predict(X_input)[0]), 3)
            iti              = calc_iti(dg_preds)
            label, color, rec = get_recommendation(iti, dg_preds)
            proba            = float(models["classifier"].predict_proba(X_input)[0][1])

        st.markdown("---")
        st.subheader("Resultados")

        c1, c2, c3 = st.columns(3)
        c1.metric("ITI Score", f"{iti:.1f} / 100")
        c2.metric("Clasificacion", label)
        c3.metric("Probabilidad favorable", f"{proba*100:.1f}%")

        # Barra ITI
        fig_bar = make_iti_bar(iti)
        st.pyplot(fig_bar, use_container_width=True)
        plt.close()

        col_radar, col_tabla = st.columns([1, 1])
        with col_radar:
            fig_radar = make_radar(dg_preds, mol_name_input)
            st.pyplot(fig_radar)
            plt.close()
        with col_tabla:
            rows = []
            for key, name in zip(PROT_KEYS, PROT_NAMES):
                dg    = dg_preds.get(key, 0)
                interp = "Alta afinidad" if dg < -6 else "Media" if dg < -4 else "Baja afinidad"
                rows.append({"Proteina": name, "Rol": PROT_ROLES[name],
                             "DG pred (kcal/mol)": dg, "Interpretacion": interp})
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        st.markdown(f"""
        <div style='background:{color}22; border-left:4px solid {color};
        padding:12px 16px; border-radius:4px; margin-top:12px'>
        <b>Recomendacion:</b> {rec}
        </div>
        """, unsafe_allow_html=True)

        with st.expander("Ver descriptores moleculares"):
            d2 = calc_descriptors_2d(mol)
            st.dataframe(pd.DataFrame([d2]).T.rename(columns={0: "Valor"}),
                         use_container_width=True)

        csv_data = results_to_csv(mol_name_input, dg_preds, iti, label)
        st.download_button("Descargar resultados (CSV)", data=csv_data,
            file_name=f"resultados_{mol_name_input.replace(' ','_')}.csv",
            mime="text/csv")

        st.session_state["historial"].append({
            "nombre":  mol_name_input,
            "smiles":  smiles_input,
            "iti":     iti,
            "label":   label,
            "dg_dict": dg_preds.copy()
        })
        st.success("Resultado guardado en el historial de sesion.")

# ════════════════════════════════════════════════════════════════
# PESTANA 2 — Ranking ITI
# ════════════════════════════════════════════════════════════════
with tab2:
    st.subheader("Ranking ITI del dataset del TFM")
    col_f1, col_f2 = st.columns(2)
    with col_f1:
        grupo_filtro = st.selectbox("Filtrar por grupo:",
            ["Todos"] + list(df_data["grupo_label"].dropna().unique()),
            key="filtro_ranking")
    with col_f2:
        n_mostrar = st.slider("Moleculas a mostrar:", 5, 60, 20)

    df_rank = df_data if grupo_filtro == "Todos" else df_data[df_data["grupo_label"] == grupo_filtro]
    df_rank = df_rank.dropna(subset=["ITI_score"]).sort_values("ITI_score", ascending=False).head(n_mostrar)

    fig_rank, ax_rank = plt.subplots(figsize=(10, max(4, len(df_rank) * 0.35)))
    colors_bar = [GRUPO_COLORS.get(g, "#888") for g in df_rank["grupo_label"]]
    ax_rank.barh(range(len(df_rank)), df_rank["ITI_score"],
                 color=colors_bar, alpha=0.85, edgecolor='white')
    ax_rank.set_yticks(range(len(df_rank)))
    ax_rank.set_yticklabels(df_rank["nombre"], fontsize=8)
    ax_rank.invert_yaxis()
    ax_rank.set_xlabel("ITI Score (0-100)", fontsize=11)
    ax_rank.axvline(66, color='#3BB273', linestyle='--', linewidth=1.5, alpha=0.8, label='Favorable')
    ax_rank.axvline(33, color='#F18F01', linestyle='--', linewidth=1.5, alpha=0.8, label='Moderado')
    ax_rank.legend(fontsize=9)
    ax_rank.spines['top'].set_visible(False)
    ax_rank.spines['right'].set_visible(False)
    plt.tight_layout()
    st.pyplot(fig_rank)
    plt.close()

    cols_mostrar = ["nombre", "grupo_label", "ITI_score"]
    docking_disp = [c for c in df_data.columns if c.startswith("dG_") and c in df_rank.columns]
    df_tabla = df_rank[cols_mostrar + docking_disp].rename(
        columns={"nombre": "Nombre", "grupo_label": "Grupo", "ITI_score": "ITI Score"})
    st.dataframe(df_tabla.reset_index(drop=True), use_container_width=True, hide_index=True)
    st.download_button("Descargar tabla (CSV)",
        data=df_tabla.to_csv(index=False).encode(),
        file_name="ranking_ITI.csv", mime="text/csv")

# ════════════════════════════════════════════════════════════════
# PESTANA 3 — Comparador
# ════════════════════════════════════════════════════════════════
with tab3:
    st.subheader("Comparar moleculas")
    default_mols = df_data.sort_values("ITI_score", ascending=False)["nombre"].head(3).tolist() \
                   if "ITI_score" in df_data.columns else []
    mols_comparar = st.multiselect("Selecciona hasta 4 moleculas:",
        df_data["nombre"].tolist(), default=default_mols)

    if len(mols_comparar) < 2:
        st.info("Selecciona al menos 2 moleculas.")
    elif len(mols_comparar) > 4:
        st.warning("Maximo 4 moleculas.")
    else:
        comp_data, comp_iti = {}, {}
        for mol_name in mols_comparar:
            fila = df_base[df_base["nombre_entrada"] == mol_name]
            if fila.empty:
                continue
            mol = Chem.MolFromSmiles(fila["smiles_pubchem"].values[0])
            if mol is None:
                continue
            X_inp = build_feature_vector(mol, models["features"]).reshape(1, -1)
            dg    = {k: round(float(models[k].predict(X_inp)[0]), 3) for k in PROT_KEYS if k in models}
            comp_data[mol_name] = dg
            comp_iti[mol_name]  = calc_iti(dg)

        fig_comp = make_comparison_chart(comp_data)
        st.pyplot(fig_comp)
        plt.close()

        filas_iti = [{"Molecula": n, "ITI Score": v,
                      "Clasificacion": get_recommendation(v, comp_data[n])[0]}
                     for n, v in comp_iti.items()]
        st.dataframe(pd.DataFrame(filas_iti), use_container_width=True, hide_index=True)

        cols_r = st.columns(len(mols_comparar))
        for col_r, (name, dg_dict) in zip(cols_r, comp_data.items()):
            fig_r = make_radar(dg_dict, name)
            col_r.pyplot(fig_r)
            plt.close()

# ════════════════════════════════════════════════════════════════
# PESTANA 4 — Historial
# ════════════════════════════════════════════════════════════════
with tab4:
    st.subheader("Historial de analisis de esta sesion")

    if not st.session_state["historial"]:
        st.info("Analiza alguna molecula en la primera pestana para ver el historial.")
    else:
        # Tabla resumen
        filas = [{"Nombre": e["nombre"], "ITI Score": e["iti"],
                  "Clasificacion": e["label"],
                  "SMILES": e["smiles"][:40] + "..." if len(e["smiles"]) > 40 else e["smiles"]}
                 for e in st.session_state["historial"]]
        st.dataframe(pd.DataFrame(filas), use_container_width=True, hide_index=True)

        # Heatmap de la sesión
        if len(st.session_state["historial"]) >= 2:
            st.markdown("**Mapa de calor de la sesion:**")
            fig_heat = make_session_heatmap(st.session_state["historial"])
            st.pyplot(fig_heat)
            plt.close()

            # Comparador rápido desde historial
            st.markdown("**Comparar moleculas del historial:**")
            names_hist = [e["nombre"] for e in st.session_state["historial"]]
            sel_hist   = st.multiselect("Selecciona:", names_hist,
                                        default=names_hist[:2], key="hist_comp")
            if len(sel_hist) >= 2:
                comp_hist = {e["nombre"]: e["dg_dict"]
                             for e in st.session_state["historial"]
                             if e["nombre"] in sel_hist}
                fig_hist = make_comparison_chart(comp_hist)
                st.pyplot(fig_hist)
                plt.close()

        if st.button("Limpiar historial"):
            st.session_state["historial"] = []
            st.rerun()

        filas_csv = []
        for entry in st.session_state["historial"]:
            for key, name in zip(PROT_KEYS, PROT_NAMES):
                filas_csv.append({"Nombre": entry["nombre"], "Proteina": name,
                                   "DG": entry["dg_dict"].get(key, ""),
                                   "ITI": entry["iti"], "Clasif.": entry["label"]})
        st.download_button("Descargar historial (CSV)",
            data=pd.DataFrame(filas_csv).to_csv(index=False).encode(),
            file_name="historial_analisis.csv", mime="text/csv")

st.markdown("""
<div style="margin-top: 2rem; padding-top: 1rem;
border-top: 1px solid #e2e8f0; color: #94a3b8; font-size: 0.8rem;">
    TFM · Evaluacion computacional de nanoparticulas polimericas y dendrimeros ·
    Diego Vallina Alvarez · Universidad de Burgos · 2025-2026
</div>
""", unsafe_allow_html=True)
