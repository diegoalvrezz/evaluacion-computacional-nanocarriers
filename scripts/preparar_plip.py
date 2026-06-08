"""
preparar_plip.py
================
Combina proteína PDB + ligando PDBQT en un único PDB listo para PLIP online.
Solo para las 3 proteínas más relevantes × 5 moléculas top ITI = 15 archivos.

Uso:
    python3 preparar_plip.py
"""
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path

BASE    = Path.home() / "Desktop/TFM_DiegoVallina/docking_poses"
OUT_DIR = Path.home() / "Desktop/TFM_DiegoVallina/plip_inputs"
OUT_DIR.mkdir(exist_ok=True)

PROTEINAS = ["7A65", "1CX8", "4LRH"]
MOLECULAS = [
    "Triethylene_glycol",
    "Tetraethylenepentamine",
    "Pentaethylenehexamine",
    "Chitotriose",
    "Chitobiose",
]

def pdbqt_to_pdb_lines(pdbqt_path, res_name="LIG", chain="Z"):
    """Extrae modo 1 del PDBQT y formatea como HETATM."""
    lines = []
    in_m1 = False
    atom_num = 1
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                in_m1 = line.split()[1].strip() == "1"
            if in_m1 and (line.startswith("ATOM") or line.startswith("HETATM")):
                # Reformatear como HETATM con resname LIG chain Z
                atom_name = line[12:16].strip()
                x = line[30:38]
                y = line[38:46]
                z = line[46:54]
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                new_line = (f"HETATM{atom_num:5d}  {atom_name:<3s} {res_name:3s} "
                           f"{chain}   1    {x}{y}{z}  1.00  0.00          {element:>2s}\n")
                lines.append(new_line)
                atom_num += 1
            if line.startswith("ENDMDL") and in_m1:
                break
    if not lines:
        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_name = line[12:16].strip()
                    x = line[30:38]
                    y = line[38:46]
                    z = line[46:54]
                    element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                    new_line = (f"HETATM{atom_num:5d}  {atom_name:<3s} {res_name:3s} "
                               f"{chain}   1    {x}{y}{z}  1.00  0.00          {element:>2s}\n")
                    lines.append(new_line)
                    atom_num += 1
    return lines

print("Preparando archivos para PLIP online...")
count = 0

for prot_code in PROTEINAS:
    pdb_path = BASE / f"{prot_code}.pdb"
    if not pdb_path.exists():
        print(f"  ✗ {prot_code}.pdb no encontrado")
        continue

    with open(pdb_path) as f:
        prot_lines = [l for l in f.readlines()
                      if l.startswith("ATOM") or l.startswith("HETATM") or
                      l.startswith("TER") or l.startswith("END")]

    for mol_name in MOLECULAS:
        pdbqt_path = BASE / prot_code / f"{mol_name}_docked.pdbqt"
        if not pdbqt_path.exists():
            print(f"  ✗ {mol_name} en {prot_code}")
            continue

        lig_lines = pdbqt_to_pdb_lines(pdbqt_path)
        if not lig_lines:
            print(f"  ✗ Sin coordenadas: {mol_name}")
            continue

        out_name = f"{prot_code}_{mol_name}.pdb"
        out_path = OUT_DIR / out_name

        with open(out_path, "w") as f:
            f.write(f"REMARK  Proteína: {prot_code}  Ligando: {mol_name}\n")
            f.write(f"REMARK  Preparado para PLIP online\n")
            for line in prot_lines:
                if not line.startswith("END"):
                    f.write(line)
            f.write("TER\n")
            for line in lig_lines:
                f.write(line)
            f.write("END\n")

        print(f"  ✓ {out_name}")
        count += 1

print(f"\n✓ {count} archivos generados en: {OUT_DIR}")
print("\nInstrucciones:")
print("1. Ve a https://plip-tool.biotec.tu-dresden.de/plip/session")
print("2. Sube cada archivo .pdb")
print("3. PLIP detectará automáticamente el ligando (cadena Z, resname LIG)")
print("4. Descarga el informe y anota los residuos de interacción")
