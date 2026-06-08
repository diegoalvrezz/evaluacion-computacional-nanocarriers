# Computational Evaluation of Polymeric Nanoparticles and Dendrimers as Drug Delivery Systems

**TFM — Máster Universitario en Ingeniería Biomédica**  
Universidad de Burgos · Curso académico 2025–2026

**Autor:** Diego Vallina Álvarez  
**Directores:** Prof. Santiago Aparicio Martínez · Prof. Pedro Ángel Marcos Villa  
**Grupo:** Química-Física Computacional · Departamento de Química, UBU

---

## Descripción

Plataforma computacional para la evaluación sistemática de 58 moléculas representativas de nanopartículas poliméricas y dendrímeros frente a seis proteínas de barrera biológica que determinan la biodistribución y el metabolismo de los fármacos:

| PDB  | Proteína   | Función                          |
|------|------------|----------------------------------|
| 7A65 | P-gp/MDR1  | Bomba de expulsión de fármacos   |
| 1TQN | CYP3A4     | Metabolismo hepático             |
| 1CX8 | TfR1       | Targeting tumoral (endocitosis)  |
| 4LRH | FRα        | Targeting tumoral (folato)       |
| 1LYZ | Lisozima   | Biocompatibilidad tisular        |
| 1AO6 | HSA        | Transporte plasmático            |

---

## Metodología

1. **Cálculos DFT** (TURBOMOLE 7.8.1, B3LYP-D4/def2-TZVP) — optimización de geometrías y perfiles sigma COSMO-RS para las 58 moléculas en clúster HPC CENITS
2. **Docking molecular** (AutoDock Vina 1.2.5) — 342 cálculos sistemáticos
3. **Análisis de interacciones** (PLIP 2025) — 708 contactos no covalentes en 15 complejos
4. **Descriptores moleculares** (RDKit, Mordred) — 2908 características por molécula
5. **Propiedades ADMET** (SwissADME) — 49 propiedades farmacocinéticas
6. **Modelos ML** (scikit-learn) — regresión de afinidad (RF: R²=0.817 en P-gp), ITI y clasificador (AUC=0.983)
7. **Herramienta web** (Streamlit) — recomendación de carrier óptimo a partir de SMILES

---

## Estructura del repositorio

```
├── app/                    # Aplicación Streamlit + modelos entrenados
│   ├── app.py
│   └── models/             # Modelos RF serializados (.joblib)
├── data/                   # Datasets y descriptores
│   ├── base_molecular_pubchem.csv
│   ├── dataset_ML.csv
│   ├── descriptores_2D.csv
│   ├── descriptores_3D_mordred.csv
│   ├── docking_energias.csv
│   ├── sigma_profiles.csv
│   ├── admet_swissadme.csv
│   ├── plip_parsed.csv
│   └── admet_raw/          # CSVs originales de SwissADME por tanda
├── docking/                # Archivos de docking
│   ├── receptores/         # PDB de las 6 proteínas
│   ├── ligandos/           # PDBQT de las 58 moléculas
│   ├── resultados/         # Poses de docking (.pdbqt)
│   ├── plip_inputs/        # PDB combinados proteína+ligando para PLIP
│   └── plip_resultados/    # Imágenes y XMLs de PLIP
├── docs/                   # Memoria y anexos del TFM
│   ├── memoria.pdf
│   ├── anexos.pdf
│   ├── memoria.qmd
│   ├── anexos.qmd
│   └── qmd/                # Capítulos en formato Quarto
├── results/                # Resultados y figuras
│   ├── figures/            # Todas las figuras del TFM
│   ├── ml_results/         # CSVs de resultados ML
│   └── models/             # Modelos entrenados
├── scripts/                # Scripts Python del proyecto
│   ├── calcular_descriptores_2D.py
│   ├── calcular_fingerprints.py
│   ├── extraer_sigma_profiles.py
│   ├── combinar_descriptores.py
│   ├── lanzar_docking.py
│   ├── extraer_docking_energias.py
│   ├── modelos_ML_v2.py
│   ├── preparar_plip.py
│   └── graficas_pro.py
└── structures/             # Estructuras moleculares
    ├── xyz/                # Geometrías DFT optimizadas (.xyz)
    ├── sdf/                # Estructuras 3D (.sdf)
    ├── molden/             # Archivos de visualización (.molden)
    └── turbomole/          # Cálculos TURBOMOLE completos
```

---

## Resultados principales

| Métrica | Valor |
|---------|-------|
| Moléculas con DFT completo | 58/58 |
| Dockings válidos | 342/348 |
| R² Random Forest (P-gp) | 0.817 |
| R² Random Forest (TfR1) | 0.753 |
| AUC clasificador | 0.983 |
| Interacciones PLIP | 708 |

**Top candidatos ITI:** Triethylene glycol (100), Tetraethylenepentamine (87.3), Pentaethylenehexamine (85.0), Chitotriose (84.7), Chitobiose (82.7)

---

## Uso de la aplicación Streamlit

```bash
conda activate streamlit-env
cd app
streamlit run app.py
```

La app se abre en `http://localhost:8501`. Introduce el SMILES de cualquier molécula para obtener su perfil de afinidad, ITI y recomendación de carrier.

---

## Dependencias principales

```
rdkit, mordred, scikit-learn, shap, pandas, numpy
matplotlib, seaborn, streamlit, py3Dmol, stmol, joblib
```

---

## Cita

Si utilizas este trabajo, por favor cita:

> Vallina Álvarez, D. (2026). *Evaluación computacional de nanopartículas poliméricas y dendrímeros como sistemas de liberación de fármacos*. TFM, Universidad de Burgos.

---

## Licencia

MIT License — libre uso con atribución.
