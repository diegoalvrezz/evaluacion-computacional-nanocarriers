"""
calcular_descriptores_2D.py
===========================
Calcula descriptores 2D de RDKit para todas las moléculas del CSV.
Genera un archivo descriptores_2D.csv con todos los descriptores.

Requisitos:
    pip install rdkit pandas

Uso:
    python calcular_descriptores_2D.py
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, QED
from rdkit.Chem import Lipinski, Crippen, rdFreeSASA
import warnings
warnings.filterwarnings('ignore')

# ── Configuración ────────────────────────────────────────────
INPUT_CSV  = "base_molecular_pubchem.csv"
OUTPUT_CSV = "descriptores_2D.csv"

# ── Descriptores a calcular ──────────────────────────────────
def calcular_descriptores(mol, nombre, smiles):
    """Calcula todos los descriptores 2D relevantes para drug delivery."""
    if mol is None:
        return None

    d = {}
    d["nombre"]  = nombre
    d["smiles"]  = smiles

    # ── Básicos ──────────────────────────────────────────────
    d["MW"]       = round(Descriptors.MolWt(mol), 3)
    d["ExactMW"]  = round(Descriptors.ExactMolWt(mol), 3)
    d["LogP"]     = round(Crippen.MolLogP(mol), 3)
    d["TPSA"]     = round(Descriptors.TPSA(mol), 3)
    d["HBD"]      = Lipinski.NumHDonors(mol)
    d["HBA"]      = Lipinski.NumHAcceptors(mol)
    d["RotBonds"] = Lipinski.NumRotatableBonds(mol)
    d["NumRings"] = rdMolDescriptors.CalcNumRings(mol)
    d["NumArRings"]= rdMolDescriptors.CalcNumAromaticRings(mol)
    d["NumHetAt"] = rdMolDescriptors.CalcNumHeteroatoms(mol)
    d["NumHeavyAt"]= mol.GetNumHeavyAtoms()
    d["NumAtoms"] = mol.GetNumAtoms()
    d["NumBonds"] = mol.GetNumBonds()
    d["Charge"]   = Chem.GetFormalCharge(mol)

    # ── Regla de Lipinski ────────────────────────────────────
    d["Lipinski_MW"]  = 1 if d["MW"] <= 500 else 0
    d["Lipinski_LogP"]= 1 if d["LogP"] <= 5 else 0
    d["Lipinski_HBD"] = 1 if d["HBD"] <= 5 else 0
    d["Lipinski_HBA"] = 1 if d["HBA"] <= 10 else 0
    d["Lipinski_OK"]  = 1 if sum([d["Lipinski_MW"], d["Lipinski_LogP"],
                                   d["Lipinski_HBD"], d["Lipinski_HBA"]]) >= 3 else 0

    # ── Drug-likeness ────────────────────────────────────────
    d["QED"] = round(QED.qed(mol), 4)

    # ── Complejidad y forma ──────────────────────────────────
    d["FractionCSP3"]   = round(rdMolDescriptors.CalcFractionCSP3(mol), 4)
    d["NumStereocenters"]= rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    d["NumAmideBonds"]  = rdMolDescriptors.CalcNumAmideBonds(mol)

    # ── Relevantes para drug delivery ────────────────────────
    d["MR"]     = round(Crippen.MolMR(mol), 3)       # Refracción molar
    d["LabuteASA"] = round(rdMolDescriptors.CalcLabuteASA(mol), 3)  # Área superficial

    # ── Índices topológicos ──────────────────────────────────
    d["BertzCT"]  = round(Descriptors.BertzCT(mol), 3)
    d["Chi0v"]    = round(Descriptors.Chi0v(mol), 4)
    d["Chi1v"]    = round(Descriptors.Chi1v(mol), 4)
    d["Kappa1"]   = round(Descriptors.Kappa1(mol), 4)
    d["Kappa2"]   = round(Descriptors.Kappa2(mol), 4)
    d["Kappa3"]   = round(Descriptors.Kappa3(mol), 4)

    # ── Conteo de grupos funcionales relevantes ──────────────
    d["NumOH"]    = sum(1 for a in mol.GetAtoms()
                        if a.GetAtomicNum() == 8 and a.GetTotalNumHs() > 0)
    d["NumNH"]    = sum(1 for a in mol.GetAtoms()
                        if a.GetAtomicNum() == 7 and a.GetTotalNumHs() > 0)
    d["NumCOOH"]  = len(mol.GetSubstructMatches(
                        Chem.MolFromSmarts('[CX3](=O)[OX2H1]')))
    d["NumEster"] = len(mol.GetSubstructMatches(
                        Chem.MolFromSmarts('[CX3](=O)[OX2][CX4]')))
    d["NumAmide"] = len(mol.GetSubstructMatches(
                        Chem.MolFromSmarts('[CX3](=O)[NX3]')))
    d["NumAmine"] = len(mol.GetSubstructMatches(
                        Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')))

    # ── Biodegradabilidad (enlaces hidrolizables) ────────────
    d["NumHydrolyzable"] = d["NumEster"] + d["NumAmide"] + \
                           len(mol.GetSubstructMatches(
                               Chem.MolFromSmarts('[CX3](=O)[OX2][CX3](=O)')))

    return d


# ── Main ─────────────────────────────────────────────────────
def main():
    df = pd.read_csv(INPUT_CSV)
    print(f"Moléculas en CSV: {len(df)}")

    resultados = []
    fallos = []

    for _, row in df.iterrows():
        nombre = row["nombre_entrada"]
        smiles = row["smiles_pubchem"]
        grupo  = row["grupo"]
        status = row["status"]

        if status != "ok":
            print(f"  [SKIP] {nombre} (status={status})")
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  [ERROR] SMILES inválido: {nombre}")
            fallos.append(nombre)
            continue

        desc = calcular_descriptores(mol, nombre, smiles)
        if desc:
            desc["grupo"] = grupo
            resultados.append(desc)
            print(f"  ✓ {nombre}")
        else:
            fallos.append(nombre)

    # Guardar CSV
    df_out = pd.DataFrame(resultados)
    # Reordenar columnas: nombre, grupo, smiles primero
    cols = ["nombre", "grupo", "smiles"] + \
           [c for c in df_out.columns if c not in ["nombre", "grupo", "smiles"]]
    df_out = df_out[cols]
    df_out.to_csv(OUTPUT_CSV, index=False)

    print(f"\n{'='*50}")
    print(f"  ✓ Descriptores calculados: {len(resultados)}")
    print(f"  ✗ Fallos: {len(fallos)}")
    print(f"  Columnas (descriptores): {len(df_out.columns) - 3}")
    print(f"  Archivo guardado: {OUTPUT_CSV}")
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
