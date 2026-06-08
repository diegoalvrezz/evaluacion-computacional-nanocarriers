"""
generar_xyz.py
==============
Genera archivos .xyz con geometría 3D optimizada (MMFF94) para cada molécula
del CSV base_molecular_pubchem.csv, listos para importar en TmoleX.

Requisitos:
    pip install rdkit

Uso:
    python generar_xyz.py

    Por defecto busca 'base_molecular_pubchem.csv' en la misma carpeta.
    Los .xyz se guardan en la subcarpeta  xyz_molecules/

Notas:
    - Se eliminan duplicados por CID (ej. 6-Hydroxyhexanoic acid y 6-Hydroxycaproic
      acid tienen el mismo CID 14490 → solo se genera uno).
    - Cisplatino (Pt) se genera igualmente; TmoleX/TURBOMOLE asigna def2-TZVP+ECP
      automáticamente para metales pesados.
    - Si la optimización MMFF94 falla (raro), se usa UFF como fallback.
    - Si ambas fallan, se guarda la geometría 2D embebida (aviso en consola).
"""

import os
import sys
import csv

# ── Importar RDKit ──────────────────────────────────────────────────────────
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
except ImportError:
    print("ERROR: RDKit no está instalado.")
    print("Instálalo con:  pip install rdkit")
    sys.exit(1)


# ── Configuración ────────────────────────────────────────────────────────────
CSV_FILE   = "base_molecular_pubchem.csv"
OUTPUT_DIR = "xyz_molecules"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ── Función principal ────────────────────────────────────────────────────────
def smiles_to_xyz(smiles: str, mol_name: str, output_path: str) -> bool:
    """
    Convierte un SMILES a geometría 3D y guarda en formato .xyz.
    Devuelve True si tuvo éxito, False si falló.
    """
    # 1. Parsear SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  [ERROR] SMILES inválido: {smiles}")
        return False

    # 2. Añadir hidrógenos explícitos
    mol = Chem.AddHs(mol)

    # 3. Generar coordenadas 3D (embedding)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        print(f"  [WARN] ETKDGv3 falló para {mol_name}, intentando método clásico...")
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result == -1:
        print(f"  [ERROR] No se pudo generar geometría 3D para {mol_name}")
        return False

    # 4. Optimizar geometría con MMFF94 (fallback a UFF)
    ff_result = AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94", maxIters=2000)
    if ff_result != 0:
        print(f"  [WARN] MMFF94 no convergió para {mol_name}, probando UFF...")
        ff_result = AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        if ff_result != 0:
            print(f"  [WARN] UFF tampoco convergió para {mol_name}, usando geometría sin optimizar")

    # 5. Extraer coordenadas
    conf = mol.GetConformer()
    atoms = mol.GetAtoms()
    positions = conf.GetPositions()

    # 6. Escribir archivo .xyz
    num_atoms = mol.GetNumAtoms()
    with open(output_path, "w") as f:
        f.write(f"{num_atoms}\n")
        f.write(f"{mol_name}\n")
        for atom, pos in zip(atoms, positions):
            symbol = atom.GetSymbol()
            f.write(f"{symbol:<4s}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}\n")

    return True


# ── Leer CSV y procesar moléculas ────────────────────────────────────────────
def main():
    if not os.path.exists(CSV_FILE):
        print(f"ERROR: No se encuentra '{CSV_FILE}'")
        print("Asegúrate de que el CSV está en la misma carpeta que este script.")
        sys.exit(1)

    seen_cids   = set()       # para eliminar duplicados por CID
    ok_count    = 0
    fail_count  = 0
    skip_count  = 0
    fail_list   = []

    with open(CSV_FILE, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            grupo  = row["grupo"].strip()
            nombre = row["nombre_entrada"].strip()
            smiles = row["smiles_pubchem"].strip()
            cid    = row["cid"].strip()
            status = row["status"].strip()

            # Saltar filas con status != ok
            if status != "ok":
                print(f"  [SKIP] {nombre} (status={status})")
                skip_count += 1
                continue

            # Eliminar duplicados por CID
            if cid in seen_cids:
                print(f"  [DUP]  {nombre} (CID {cid} ya procesado, omitido)")
                skip_count += 1
                continue
            seen_cids.add(cid)

            # Nombre de archivo seguro (sin caracteres especiales)
            safe_name = nombre.replace(" ", "_").replace(",", "").replace("'", "").replace("/", "-")
            out_file  = os.path.join(OUTPUT_DIR, f"{safe_name}.xyz")

            print(f"  Procesando [{grupo}]  {nombre}  (CID {cid})...")

            success = smiles_to_xyz(smiles, nombre, out_file)
            if success:
                ok_count += 1
            else:
                fail_count += 1
                fail_list.append(nombre)

    # ── Resumen final ────────────────────────────────────────────────────────
    print("\n" + "="*55)
    print(f"  RESUMEN")
    print("="*55)
    print(f"  ✓ Generados con éxito : {ok_count}")
    print(f"  ✗ Fallidos            : {fail_count}")
    print(f"  ⊘ Omitidos/duplicados : {skip_count}")
    print(f"  Archivos guardados en : ./{OUTPUT_DIR}/")
    if fail_list:
        print(f"\n  Moléculas que fallaron:")
        for m in fail_list:
            print(f"    - {m}")
    print("="*55)


if __name__ == "__main__":
    main()