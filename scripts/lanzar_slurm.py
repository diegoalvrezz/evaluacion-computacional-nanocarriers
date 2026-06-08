#!/usr/bin/env python3
"""
lanzar_slurm.py  — COMPLETO (3 etapas)
=======================================
Etapa 1: Optimización geometría  → B3LYP-D4/def2-TZVP
Etapa 2: Single point COSMO-RS   → B3LYP-D4/def2-TZVP + COSMO(epsilon=inf)
Etapa 3: Frecuencias IR          → B3LYP-D4/def2-TZVP + aoforce

Flujo por molécula:
  xyz → define → jobex (opt) → define SP → ridft (SP) → aoforce (freq)

Exports:
  - COSMOtherm: .cosmo  (Etapa 1 y 2)
  - Molden:     .molden (Etapa 1)
  - IR:         vibrational_spectrum (Etapa 3)

Uso:
  python3 lanzar_slurm.py [--etapa 1|2|3|all]

  --etapa 1   → solo optimización
  --etapa 2   → solo single point COSMO-RS (requiere Etapa 1 hecha)
  --etapa 3   → solo frecuencias IR (requiere Etapa 1 hecha)
  --etapa all → las 3 etapas secuenciales en el mismo job (por defecto)
"""

import os
import sys
import subprocess
import argparse

# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURACIÓN
# ═══════════════════════════════════════════════════════════════════════════

BASE_DIR    = os.path.expanduser("~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo")
XYZ_DIR     = os.path.join(BASE_DIR, "xyz_molecules")
OUTPUT_DIR  = os.path.join(BASE_DIR, "TURBOMOLE_calcs")
TURBODIR    = os.path.expanduser("~/MASTER_SALUD_2026/TmoleX2024/TURBOMOLE")
BIN_DIR     = os.path.join(TURBODIR, "bin/x86_64-unknown-linux-gnu_smp")
SCRIPTS_DIR = os.path.join(TURBODIR, "scripts")

# SLURM
PARTITION  = "lusi2"
NCORES     = 16
MEMORY_MB  = 28000
MAX_HOURS  = 72       # 72h para moléculas grandes
MAX_CYCLES = 300

# ═══════════════════════════════════════════════════════════════════════════
# UTILIDADES
# ═══════════════════════════════════════════════════════════════════════════

def get_env():
    env = os.environ.copy()
    env["TURBODIR"]          = TURBODIR
    env["TURBOMOLE_SYSNAME"] = "x86_64-unknown-linux-gnu_smp"
    env["PATH"]              = f"{SCRIPTS_DIR}:{BIN_DIR}:{env.get('PATH','')}"
    env["PARA_ARCH"]         = "SMP"
    env["PARNODES"]          = str(NCORES)
    env["OMP_NUM_THREADS"]   = str(NCORES)
    return env


def count_atoms(xyz_path):
    with open(xyz_path) as f:
        return int(f.readline().strip())


def xyz_to_coord(xyz_path, coord_path):
    ANGSTROM_TO_BOHR = 1.88972612462577
    with open(xyz_path) as f:
        lines = f.readlines()
    with open(coord_path, "w") as f:
        f.write("$coord\n")
        for line in lines[2:]:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            symbol = parts[0].lower()
            x = float(parts[1]) * ANGSTROM_TO_BOHR
            y = float(parts[2]) * ANGSTROM_TO_BOHR
            z = float(parts[3]) * ANGSTROM_TO_BOHR
            f.write(f"  {x:20.14f}  {y:20.14f}  {z:20.14f}  {symbol}\n")
        f.write("$end\n")


def get_charge(mol_name):
    exceptions = {"Chondroitin_sulfate_disaccharide": -2}
    return exceptions.get(mol_name, 0)


def is_converged(job_dir):
    return os.path.exists(os.path.join(job_dir, "GEO_OPT_CONVERGED"))


def sp_done(job_dir, mol_name):
    """Single point COSMO-RS completado si existe el .cosmo de SP."""
    return os.path.exists(os.path.join(job_dir, "SP", f"{mol_name}_SP.cosmo"))


def freq_done(job_dir):
    """Frecuencias completadas si existe vibrational_spectrum."""
    return os.path.exists(os.path.join(job_dir, "FREQ", "vibrational_spectrum"))


def run_define_opt(job_dir, charge):
    """Define para optimización: B3LYP-D4/def2-TZVP + RI + exports."""
    define_input = (
        "\n\na coord\n*\nno\n"
        "b all def2-TZVP\n*\n"
        f"eht\n{charge}\n\n\n"
        "dft\non\nfunc b3-lyp\ngrid m4\n\n"
        "ri\non\nm 8000\n\n"
        "disp\n4\n\n"
        "*\n"
    )
    return subprocess.run(
        ["define"], input=define_input, cwd=job_dir,
        capture_output=True, text=True, timeout=120, env=get_env()
    )


def run_define_sp(sp_dir, opt_dir, mol_name, charge):
    """
    Define para single point COSMO-RS.
    Usa las coordenadas optimizadas y activa COSMO con epsilon=infinity.
    """
    # Copiar coord optimizado al directorio SP
    opt_coord = os.path.join(opt_dir, "coord")
    sp_coord  = os.path.join(sp_dir, "coord")

    # Extraer geometría final del archivo GEO_OPT_CONVERGED o coord
    import shutil
    shutil.copy(opt_coord, sp_coord)

    define_input = (
        "\n\na coord\n*\nno\n"
        "b all def2-TZVP\n*\n"
        f"eht\n{charge}\n\n\n"
        "dft\non\nfunc b3-lyp\ngrid m4\n\n"
        "ri\non\nm 8000\n\n"
        "disp\n4\n\n"
        "*\n"
    )
    result = subprocess.run(
        ["define"], input=define_input, cwd=sp_dir,
        capture_output=True, text=True, timeout=120, env=get_env()
    )
    return result


def patch_control_opt(job_dir, mol_name):
    """Añade exports COSMOtherm y Molden al control de optimización."""
    control_path = os.path.join(job_dir, "control")
    if not os.path.exists(control_path):
        return False
    with open(control_path) as f:
        content = f.read()
    additions = ""
    if "$cosmo_out" not in content:
        additions += f"$cosmo_out file={mol_name}.cosmo\n"
    if "$moldenfile" not in content:
        additions += f"$moldenfile file={mol_name}.molden\n"
    if "$scfiterlimit" not in content:
        additions += "$scfiterlimit 300\n"
    if additions:
        content = content.replace("$end", additions + "$end")
        with open(control_path, "w") as f:
            f.write(content)
    return True


def patch_control_sp(sp_dir, mol_name):
    """
    Modifica el control para single point COSMO-RS:
    - Elimina $optimize (no queremos optimizar)
    - Añade $cosmo con epsilon=infinity
    - Añade export .cosmo
    """
    control_path = os.path.join(sp_dir, "control")
    if not os.path.exists(control_path):
        return False
    with open(control_path) as f:
        content = f.read()

    # Quitar bloque de optimización si existe
    lines = content.split("\n")
    lines_clean = [l for l in lines if not l.strip().startswith("$optimize")
                   and not l.strip().startswith("$statpt")]
    content = "\n".join(lines_clean)

    additions = ""
    # COSMO con epsilon infinito para COSMO-RS
    if "$cosmo" not in content:
        additions += "$cosmo\n  epsilon=infinity\n"
    if "$cosmo_out" not in content:
        additions += f"$cosmo_out file={mol_name}_SP.cosmo\n"
    if "$scfiterlimit" not in content:
        additions += "$scfiterlimit 300\n"

    content = content.replace("$end", additions + "$end")
    with open(control_path, "w") as f:
        f.write(content)
    return True


def write_slurm_all(job_dir, mol_name, num_atoms):
    """
    Script SLURM que ejecuta las 3 etapas secuencialmente:
    1. Optimización (jobex -ri)
    2. Single point COSMO-RS (ridft en subcarpeta SP/)
    3. Frecuencias IR (aoforce en subcarpeta FREQ/)
    """
    slurm_path = os.path.join(job_dir, f"job_{mol_name[:20]}.slurm")

    # Tiempo estimado según tamaño
    if num_atoms > 80:
        hours = 72
    elif num_atoms > 50:
        hours = 48
    elif num_atoms > 30:
        hours = 24
    else:
        hours = 12

    script = f"""#!/bin/bash
#SBATCH --job-name={mol_name[:15]}
#SBATCH --partition={PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={NCORES}
#SBATCH --mem={MEMORY_MB}M
#SBATCH --time={hours}:00:00
#SBATCH --output={job_dir}/slurm_%j.out
#SBATCH --error={job_dir}/slurm_%j.err

# ── Entorno TURBOMOLE ──────────────────────────────────────
export TURBODIR={TURBODIR}
export TURBOMOLE_SYSNAME=x86_64-unknown-linux-gnu_smp
export PATH={SCRIPTS_DIR}:{BIN_DIR}:$PATH
export PARA_ARCH=SMP
export PARNODES={NCORES}
export OMP_NUM_THREADS={NCORES}

echo "========================================"
echo " Molécula : {mol_name}"
echo " Nodo     : $HOSTNAME"
echo " Inicio   : $(date)"
echo "========================================"

cd {job_dir}

# ════════════════════════════════════════════
# ETAPA 1 — Optimización B3LYP-D4/def2-TZVP
# ════════════════════════════════════════════
if [ ! -f "GEO_OPT_CONVERGED" ]; then
    echo "--- ETAPA 1: Optimización geometría ---"
    jobex -ri -c {MAX_CYCLES} > jobex_output.log 2>&1
    if [ -f "GEO_OPT_CONVERGED" ]; then
        echo "ETAPA 1: CONVERGIDO $(date)"
    else
        echo "ETAPA 1: NO CONVERGIÓ — abortando"
        exit 1
    fi
else
    echo "--- ETAPA 1: Ya convergida, saltando ---"
fi

# ════════════════════════════════════════════
# ETAPA 2 — Single Point COSMO-RS
# ════════════════════════════════════════════
mkdir -p SP
if [ ! -f "SP/{mol_name}_SP.cosmo" ]; then
    echo "--- ETAPA 2: Single Point COSMO-RS ---"

    # Copiar geometría optimizada y archivos de base
    cp coord     SP/
    cp basis     SP/
    cp auxbasis  SP/
    cp mos       SP/  2>/dev/null || true

    # Generar control para SP con COSMO
    cd SP
    # Usar control de la optimización como base y modificar
    cp ../control ./control_base
    python3 -c "
import sys
with open('control_base') as f:
    content = f.read()
# Limpiar directivas de optimización
lines = [l for l in content.split('\\n')
         if not any(l.strip().startswith(k) for k in
         ['\$optimize','\$statpt','\$cosmo','\$cosmo_out','\$moldenfile'])]
content = '\\n'.join(lines)
# Añadir COSMO-RS y export
additions = '\$cosmo\\n  epsilon=infinity\\n\$cosmo_out file={mol_name}_SP.cosmo\\n\$scfconv 7\\n'
content = content.replace('\$end', additions + '\$end')
with open('control', 'w') as f:
    f.write(content)
print('control SP generado')
"
    # Lanzar single point
    ridft > ridft_SP.log 2>&1
    if [ -f "{mol_name}_SP.cosmo" ]; then
        echo "ETAPA 2: SP COSMO-RS completado $(date)"
    else
        echo "ETAPA 2: SP falló — revisar ridft_SP.log"
    fi
    cd {job_dir}
else
    echo "--- ETAPA 2: Ya completada, saltando ---"
fi

# ════════════════════════════════════════════
# ETAPA 3 — Frecuencias IR / Vibracionales
# ════════════════════════════════════════════
mkdir -p FREQ
if [ ! -f "FREQ/vibrational_spectrum" ]; then
    echo "--- ETAPA 3: Frecuencias IR (aoforce) ---"

    # Copiar geometría optimizada y archivos necesarios
    cp coord    FREQ/
    cp basis    FREQ/
    cp auxbasis FREQ/
    cp control  FREQ/
    cp mos      FREQ/ 2>/dev/null || true

    cd FREQ
    # Asegurarse que el control no tiene directivas de optimización
    python3 -c "
with open('control') as f:
    content = f.read()
lines = [l for l in content.split('\\n')
         if not any(l.strip().startswith(k) for k in
         ['\$optimize','\$statpt'])]
content = '\\n'.join(lines)
with open('control', 'w') as f:
    f.write(content)
"
    # Primero un single point para obtener los orbitales actualizados
    ridft > ridft_FREQ.log 2>&1
    # Luego calcular frecuencias
    aoforce > aoforce.log 2>&1

    if [ -f "vibrational_spectrum" ]; then
        echo "ETAPA 3: FRECUENCIAS completadas $(date)"
    else
        echo "ETAPA 3: aoforce falló — revisar aoforce.log"
    fi
    cd {job_dir}
else
    echo "--- ETAPA 3: Ya completada, saltando ---"
fi

echo "========================================"
echo " FIN: $(date)"
echo "========================================"
"""
    with open(slurm_path, "w") as f:
        f.write(script)
    # Hacer ejecutable
    os.chmod(slurm_path, 0o755)
    return slurm_path


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--etapa", default="all",
                        choices=["1","2","3","all"],
                        help="Etapa a ejecutar (default: all)")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    xyz_files = [f for f in os.listdir(XYZ_DIR) if f.endswith(".xyz")]
    if not xyz_files:
        print(f"ERROR: No se encontraron .xyz en {XYZ_DIR}")
        sys.exit(1)

    xyz_files = sorted(xyz_files,
                       key=lambda f: count_atoms(os.path.join(XYZ_DIR, f)))

    submitted = []
    skipped   = []
    failed    = []

    print(f"{'='*55}")
    print(f"  TURBOMOLE B3LYP-D4/def2-TZVP — Cluster SLURM")
    print(f"  Etapa: {args.etapa} | Partición: {PARTITION} | Núcleos: {NCORES}")
    print(f"{'='*55}\n")

    for i, xyz_file in enumerate(xyz_files, 1):
        mol_name = xyz_file.replace(".xyz", "")
        xyz_path = os.path.join(XYZ_DIR, xyz_file)
        job_dir  = os.path.join(OUTPUT_DIR, mol_name)
        num_atoms = count_atoms(xyz_path)

        # Verificar si todo ya está hecho
        if (is_converged(job_dir) and
            sp_done(job_dir, mol_name) and
            freq_done(job_dir)):
            print(f"[{i:02d}/{len(xyz_files)}] {mol_name} ({num_atoms} át.) → ⊘ TODO COMPLETO")
            skipped.append(mol_name)
            continue

        print(f"[{i:02d}/{len(xyz_files)}] {mol_name} ({num_atoms} át.)")
        os.makedirs(job_dir, exist_ok=True)

        try:
            # Si no tiene coord, generarlo desde xyz
            coord_path = os.path.join(job_dir, "coord")
            if not os.path.exists(coord_path):
                xyz_to_coord(xyz_path, coord_path)

            charge = get_charge(mol_name)

            # Si no tiene control (no se ha ejecutado define aún)
            if not os.path.exists(os.path.join(job_dir, "control")):
                print(f"         ejecutando define...")
                def_result = run_define_opt(job_dir, charge)
                with open(os.path.join(job_dir, "define_output.log"), "w") as lf:
                    lf.write(def_result.stdout)
                    lf.write(def_result.stderr)

                if not os.path.exists(os.path.join(job_dir, "control")):
                    print(f"         ✗ define falló")
                    failed.append(mol_name)
                    continue

                patch_control_opt(job_dir, mol_name)
                print(f"         define OK")
            else:
                print(f"         control ya existe, usando el actual")

            # Generar script SLURM con las 3 etapas
            slurm_path = write_slurm_all(job_dir, mol_name, num_atoms)

            # Lanzar job
            result = subprocess.run(
                ["sbatch", slurm_path],
                capture_output=True, text=True
            )

            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1]
                print(f"         ✓ Job lanzado → ID {job_id}")
                submitted.append((mol_name, job_id))
            else:
                print(f"         ✗ sbatch error: {result.stderr.strip()}")
                failed.append(mol_name)

        except Exception as e:
            print(f"         ✗ ERROR: {e}")
            failed.append(mol_name)

    # ── Resumen ─────────────────────────────────────────────────────────
    print(f"\n{'='*55}")
    print("RESUMEN")
    print(f"{'='*55}")
    print(f"  ✓ Jobs lanzados  : {len(submitted)}")
    print(f"  ⊘ Ya completos   : {len(skipped)}")
    print(f"  ✗ Fallidos       : {len(failed)}")
    print(f"\n  Monitoriza con:")
    print(f"    squeue -u ubu-iccram01")
    print(f"    squeue -u ubu-iccram01 -o '%.10i %.20j %.8T %.10M %.6D'")
    if failed:
        print(f"\n  Fallidos:")
        for m in failed:
            print(f"    - {m}")
    print(f"{'='*55}")


if __name__ == "__main__":
    main()
