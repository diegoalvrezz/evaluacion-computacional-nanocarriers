"""
extraer_sigma_profiles.py
=========================
Extrae el perfil sigma (σ-profile) de archivos .cosmo de TURBOMOLE.
Genera un CSV con ~100 descriptores por molécula (muestreo cada 0.01 e/Å²).

El perfil σ se obtiene:
1. Leyendo los segmentos del archivo .cosmo (carga y área de cada segmento)
2. Calculando σ = carga / área para cada segmento
3. Histogramando en intervalos de 0.01 entre -0.025 y 0.025 e/Å²

Uso:
    python extraer_sigma_profiles.py

Salida:
    sigma_profiles.csv
"""

import os
import numpy as np
import pandas as pd
import glob

# ── Configuración ─────────────────────────────────────────────
COSMO_BASE = os.path.expanduser(
    "~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/TURBOMOLE_calcs"
)
OUTPUT_CSV = os.path.expanduser(
    "~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo/sigma_profiles.csv"
)

# Rango y paso del perfil σ
SIGMA_MIN  = -0.025   # e/Å²
SIGMA_MAX  =  0.025   # e/Å²
SIGMA_STEP =  0.01    # e/Å² → ~50 bins (de -0.025 a 0.025)
# O si queremos ~100 valores, usamos paso 0.001 entre -0.05 y 0.05
# Ajustamos según lo que pide el profesor: "cada 0.01 generando ~100"
SIGMA_MIN  = -0.05
SIGMA_MAX  =  0.05
SIGMA_STEP =  0.001   # 100 bins entre -0.05 y 0.05

# ── Función para parsear .cosmo ───────────────────────────────
def parse_cosmo(filepath):
    """
    Lee un archivo .cosmo de TURBOMOLE y devuelve
    listas de sigma (carga/área) y área de cada segmento.
    """
    sigmas = []
    areas  = []
    
    in_segments = False
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('$segment_information'):
                in_segments = True
                continue
            if in_segments and line.startswith('$'):
                in_segments = False
                continue
            if in_segments and line.startswith('#'):
                continue
            if in_segments and line:
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        charge = float(parts[3])  # columna 4: carga
                        area   = float(parts[4])  # columna 5: área
                        if area > 0:
                            sigma = charge / area
                            sigmas.append(sigma)
                            areas.append(area)
                    except (ValueError, IndexError):
                        continue
    
    return np.array(sigmas), np.array(areas)


def sigma_profile(sigmas, areas, sigma_min, sigma_max, step):
    """
    Construye el perfil σ como histograma ponderado por área.
    """
    bins = np.arange(sigma_min, sigma_max + step, step)
    profile = np.zeros(len(bins) - 1)
    
    for sigma, area in zip(sigmas, areas):
        idx = int((sigma - sigma_min) / step)
        if 0 <= idx < len(profile):
            profile[idx] += area
    
    return profile, bins[:-1]


# ── Main ──────────────────────────────────────────────────────
def main():
    # Buscar todos los archivos SP.cosmo
    cosmo_files = glob.glob(
        os.path.join(COSMO_BASE, "*", "SP", "*_SP.cosmo")
    )
    # También buscar en Cisplatin_SVP y Chondroitin_SVP
    cosmo_files += glob.glob(
        os.path.join(COSMO_BASE, "*_SVP", "SP", "*_SP.cosmo")
    )
    
    cosmo_files = sorted(set(cosmo_files))
    print(f"Archivos .cosmo encontrados: {len(cosmo_files)}")
    
    bins = np.arange(SIGMA_MIN, SIGMA_MAX + SIGMA_STEP, SIGMA_STEP)
    bin_centers = bins[:-1] + SIGMA_STEP / 2
    col_names = [f"sigma_{b:.4f}" for b in bin_centers]
    
    resultados = []
    
    for filepath in cosmo_files:
        # Nombre de la molécula desde la ruta
        parts = filepath.split(os.sep)
        mol_dir = parts[-3]  # carpeta de la molécula
        mol_name = mol_dir.replace("_SVP", "")
        
        sigmas, areas = parse_cosmo(filepath)
        
        if len(sigmas) == 0:
            print(f"  [ERROR] Sin segmentos: {mol_name}")
            continue
        
        profile, _ = sigma_profile(sigmas, areas, SIGMA_MIN, SIGMA_MAX, SIGMA_STEP)
        
        # Normalizar por área total
        total_area = areas.sum()
        profile_norm = profile / total_area if total_area > 0 else profile
        
        row = {"nombre": mol_name}
        for col, val in zip(col_names, profile_norm):
            row[col] = round(float(val), 6)
        
        resultados.append(row)
        print(f"  ✓ {mol_name} ({len(sigmas)} segmentos, área total={total_area:.2f} Å²)")
    
    df = pd.DataFrame(resultados)
    cols = ["nombre"] + col_names
    df = df[cols]
    df.to_csv(OUTPUT_CSV, index=False)
    
    print(f"\n{'='*55}")
    print(f"  ✓ Moléculas procesadas: {len(resultados)}")
    print(f"  ✓ Descriptores σ por molécula: {len(col_names)}")
    print(f"  ✓ Rango σ: [{SIGMA_MIN}, {SIGMA_MAX}] e/Å²")
    print(f"  ✓ Paso: {SIGMA_STEP} e/Å²")
    print(f"  ✓ Archivo guardado: {OUTPUT_CSV}")
    print(f"{'='*55}")


if __name__ == "__main__":
    main()
