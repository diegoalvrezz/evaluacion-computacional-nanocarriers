"""
calcular_fingerprints.py
========================
Calcula fingerprints moleculares Morgan/ECFP4 y MACCS keys
para todas las moléculas del CSV usando RDKit.

Uso:
    python calcular_fingerprints.py

Salida:
    fingerprints_morgan.csv  → 2048 bits ECFP4 (Morgan radio=2)
    fingerprints_maccs.csv   → 167 bits MACCS keys
    fingerprints_combined.csv → Morgan + MACCS combinados
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
import warnings
warnings.filterwarnings('ignore')

INPUT_CSV = "base_molecular_pubchem.csv"
OUT_MORGAN  = "fingerprints_morgan.csv"
OUT_MACCS   = "fingerprints_maccs.csv"
OUT_COMBINED = "fingerprints_combined.csv"

def main():
    df = pd.read_csv(INPUT_CSV)
    df = df[df["status"] == "ok"].reset_index(drop=True)
    print(f"Moléculas: {len(df)}")

    morgan_rows = []
    maccs_rows  = []

    for _, row in df.iterrows():
        nombre = row["nombre_entrada"]
        smiles = row["smiles_pubchem"]
        grupo  = row["grupo"]

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  ✗ {nombre} SMILES inválido")
            continue

        # Morgan / ECFP4 (radio=2, 2048 bits)
        fp_morgan = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        morgan_bits = list(fp_morgan.ToBitString())
        morgan_row = {"nombre": nombre, "grupo": grupo}
        for i, bit in enumerate(morgan_bits):
            morgan_row[f"Morgan_{i}"] = int(bit)
        morgan_rows.append(morgan_row)

        # MACCS keys (167 bits)
        fp_maccs = MACCSkeys.GenMACCSKeys(mol)
        maccs_bits = list(fp_maccs.ToBitString())
        maccs_row = {"nombre": nombre, "grupo": grupo}
        for i, bit in enumerate(maccs_bits):
            maccs_row[f"MACCS_{i}"] = int(bit)
        maccs_rows.append(maccs_row)

        print(f"  ✓ {nombre}")

    # Guardar Morgan
    df_morgan = pd.DataFrame(morgan_rows)
    df_morgan.to_csv(OUT_MORGAN, index=False)

    # Guardar MACCS
    df_maccs = pd.DataFrame(maccs_rows)
    df_maccs.to_csv(OUT_MACCS, index=False)

    # Combinar
    df_combined = pd.merge(
        df_maccs,
        df_morgan.drop(columns=["grupo"]),
        on="nombre"
    )
    df_combined.to_csv(OUT_COMBINED, index=False)

    print(f"\n{'='*55}")
    print(f"  ✓ Moléculas procesadas: {len(morgan_rows)}")
    print(f"  ✓ Morgan/ECFP4: {len(df_morgan.columns)-2} bits")
    print(f"  ✓ MACCS keys:   {len(df_maccs.columns)-2} bits")
    print(f"  ✓ Guardados: {OUT_MORGAN}, {OUT_MACCS}, {OUT_COMBINED}")
    print(f"{'='*55}")

if __name__ == "__main__":
    main()
