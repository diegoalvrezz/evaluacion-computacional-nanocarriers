"""
lanzar_optimizaciones.py  (v3 - OPTIMIZADO)
============================================
Mejoras respecto a v2:
  - Salta moléculas ya convergidas (reanudable)
  - Ordena moléculas de menor a mayor número de átomos
  - Grid m3 en lugar de m4 (más rápido, calidad suficiente para TFM)
  - Criterios de convergencia ligeramente más laxos (gcart 3, energy 6)
  - Mueve resultados de test si existen
  - Estadísticas de tiempo por molécula

Uso:
  python lanzar_optimizaciones.py
"""

import os
import sys
import subprocess
import time
import shutil

# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURACIÓN
# ═══════════════════════════════════════════════════════════════════════════

XYZ_DIR      = r"C:\Users\Diego\Desktop\TFM\Archivos_Trabajo\xyz_molecules"
OUTPUT_DIR   = r"C:\Users\Diego\Desktop\TFM\Archivos_Trabajo\TURBOMOLE_calcs"
TEST_DIR     = r"C:\Users\Diego\Desktop\TFM\Archivos_Trabajo\TURBOMOLE_test"
TURBODIR     = r"D:\TmoleX2024\TURBOMOLE"
BIN_DIR      = r"D:\TmoleX2024\TURBOMOLE\bin\winnt_smp"
SCRIPTS_DIR  = r"D:\TmoleX2024\TURBOMOLE\scripts"
MAX_CYCLES   = 300      # aumentado por si acaso
LOG_FILE     = os.path.join(OUTPUT_DIR, "resumen_optimizaciones.log")

# ═══════════════════════════════════════════════════════════════════════════
# ENTORNO
# ═══════════════════════════════════════════════════════════════════════════

def get_env():
    env = os.environ.copy()
    env["TURBODIR"]          = TURBODIR
    env["TURBOMOLE_SYSNAME"] = "winnt_smp"
    env["PATH"]              = SCRIPTS_DIR + ";" + BIN_DIR + ";" + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = "4"
    env["PARNODES"]        = "4"
    env["PARA_ARCH"]       = "SMP"
    return env

# ═══════════════════════════════════════════════════════════════════════════
# UTILIDADES
# ═══════════════════════════════════════════════════════════════════════════

def log(msg):
    timestamp = time.strftime("%H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def count_atoms_and_elements(xyz_path):
    elements = set()
    with open(xyz_path, "r") as f:
        lines = f.readlines()
    num_atoms = int(lines[0].strip())
    for line in lines[2:]:
        parts = line.strip().split()
        if parts:
            elements.add(parts[0].capitalize())
    return num_atoms, elements


def xyz_to_coord(xyz_path, coord_path):
    ANGSTROM_TO_BOHR = 1.88972612462577
    with open(xyz_path, "r") as f:
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


def get_charge_multiplicity(mol_name):
    exceptions = {
        "Chondroitin_sulfate_disaccharide": (-2, 1),
    }
    return exceptions.get(mol_name, (0, 1))


def is_already_converged(job_dir):
    """Devuelve True si la molécula ya está calculada y convergida."""
    return os.path.exists(os.path.join(job_dir, "GEO_OPT_CONVERGED"))


def move_test_results(test_dir, output_dir):
    """
    Si existe TURBOMOLE_test/L-Lactic_acid con GEO_OPT_CONVERGED,
    lo mueve a TURBOMOLE_calcs para no recalcular.
    """
    if not os.path.exists(test_dir):
        return
    moved = []
    for mol_dir in os.listdir(test_dir):
        src = os.path.join(test_dir, mol_dir)
        dst = os.path.join(output_dir, mol_dir)
        if os.path.isdir(src) and is_already_converged(src):
            if not os.path.exists(dst):
                shutil.copytree(src, dst)
                moved.append(mol_dir)
    if moved:
        log(f"  Resultados de test movidos a TURBOMOLE_calcs: {', '.join(moved)}")


def run_define(job_dir, charge, env):
    """Ejecuta define de forma no interactiva. Grid m3 para mayor velocidad."""
    define_input = (
        "\n"
        "\n"
        "a coord\n"
        "*\n"
        "no\n"
        "b all def2-TZVP\n"
        "*\n"
        "eht\n"
        f"{charge}\n"
        "\n"
        "\n"
        "dft\n"
        "on\n"
        "func b3-lyp\n"
        "grid m3\n"      # ← m3 en lugar de m4, más rápido
        "\n"
        "ri\n"
        "on\n"
        "m 5000\n"
        "\n"
        "disp\n"
        "4\n"
        "\n"
        "*\n"
    )
    result = subprocess.run(
        ["define"],
        input=define_input,
        cwd=job_dir,
        env=env,
        capture_output=True,
        text=True,
        timeout=120
    )
    return result


def patch_control(job_dir, mol_name):
    """Añade exports COSMO y Molden, y aumenta límite SCF."""
    control_path = os.path.join(job_dir, "control")
    if not os.path.exists(control_path):
        return False
    with open(control_path, "r") as f:
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


def run_jobex(job_dir, env):
    result = subprocess.run(
        ["jobex", "-ri", "-c", str(MAX_CYCLES)],
        cwd=job_dir,
        env=env,
        capture_output=True,
        text=True,
        timeout=3600 * 12   # 12 horas máximo por molécula
    )
    return result


def check_convergence(job_dir, result):
    if os.path.exists(os.path.join(job_dir, "GEO_OPT_CONVERGED")):
        return True
    combined = result.stdout + result.stderr
    keywords = [
        "GEOMETRY OPTIMIZATION CONVERGED",
        "convergence criteria satisfied",
        "CONVERGENCE REACHED",
    ]
    for kw in keywords:
        if kw.lower() in combined.lower():
            return True
    return False


def format_elapsed(seconds):
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    if h > 0:
        return f"{h}h {m}m {s}s"
    return f"{m}m {s}s"


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with open(LOG_FILE, "w", encoding="utf-8") as f:
        f.write(f"{'='*60}\n")
        f.write(f"TURBOMOLE B3LYP-D4/def2-TZVP — Optimización geométrica v3\n")
        f.write(f"Inicio: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"{'='*60}\n\n")

    # Mover resultados del test si existen
    move_test_results(TEST_DIR, OUTPUT_DIR)

    # Leer todos los .xyz y ordenar por número de átomos (menor → mayor)
    xyz_files = [f for f in os.listdir(XYZ_DIR) if f.endswith(".xyz")]
    if not xyz_files:
        print(f"ERROR: No se encontraron .xyz en {XYZ_DIR}")
        sys.exit(1)

    def get_natoms(f):
        with open(os.path.join(XYZ_DIR, f)) as fh:
            return int(fh.readline().strip())

    xyz_files = sorted(xyz_files, key=get_natoms)
 

    env = get_env()
    ok_list      = []
    fail_list    = []
    skipped_list = []
    total_start  = time.time()

    log(f"Moléculas a procesar: {len(xyz_files)} (ordenadas menor→mayor tamaño)\n")

    for i, xyz_file in enumerate(xyz_files, 1):
        mol_name = xyz_file.replace(".xyz", "")
        xyz_path = os.path.join(XYZ_DIR, xyz_file)
        job_dir  = os.path.join(OUTPUT_DIR, mol_name)
        num_atoms, _ = count_atoms_and_elements(xyz_path)

        # ── Saltar si ya convergió ──────────────────────────────────────
        if is_already_converged(job_dir):
            log(f"[{i:02d}/{len(xyz_files)}] {mol_name} ({num_atoms} átomos) → ⊘ YA CONVERGIDO, saltando")
            skipped_list.append(mol_name)
            ok_list.append(mol_name)
            continue

        log(f"[{i:02d}/{len(xyz_files)}] {mol_name} ({num_atoms} átomos)")
        mol_start = time.time()

        os.makedirs(job_dir, exist_ok=True)

        try:
            # PASO 1 — coord
            xyz_to_coord(xyz_path, os.path.join(job_dir, "coord"))

            # PASO 2 — carga/multiplicidad
            charge, multiplicity = get_charge_multiplicity(mol_name)

            # PASO 3 — define
            log(f"         ejecutando define...")
            def_result = run_define(job_dir, charge, env)

            with open(os.path.join(job_dir, "define_output.log"), "w") as lf:
                lf.write("=== STDOUT ===\n" + def_result.stdout)
                lf.write("=== STDERR ===\n" + def_result.stderr)

            if not os.path.exists(os.path.join(job_dir, "control")):
                log(f"         ✗ define no generó 'control' — ver define_output.log")
                fail_list.append(mol_name)
                continue

            # PASO 4 — parchear control
            patch_control(job_dir, mol_name)

            # PASO 5 — jobex
            log(f"         lanzando jobex...")
            job_result = run_jobex(job_dir, env)

            with open(os.path.join(job_dir, "jobex_output.log"), "w") as lf:
                lf.write("=== STDOUT ===\n" + job_result.stdout)
                lf.write("=== STDERR ===\n" + job_result.stderr)

            elapsed = time.time() - mol_start

            # PASO 6 — convergencia
            if check_convergence(job_dir, job_result):
                log(f"         ✓ CONVERGIDO en {format_elapsed(elapsed)}")
                ok_list.append(mol_name)
            else:
                log(f"         ✗ NO convergió en {format_elapsed(elapsed)} (rc={job_result.returncode})")
                fail_list.append(mol_name)

        except subprocess.TimeoutExpired:
            log(f"         ✗ TIMEOUT (>12h)")
            fail_list.append(mol_name)
        except Exception as e:
            log(f"         ✗ ERROR: {e}")
            fail_list.append(mol_name)

    # ── Resumen ──────────────────────────────────────────────────────────
    total_elapsed = time.time() - total_start
    log("\n" + "="*60)
    log("RESUMEN FINAL")
    log("="*60)
    log(f"  ✓ Convergidos      : {len(ok_list)}")
    log(f"  ⊘ Ya existían      : {len(skipped_list)}")
    log(f"  ✗ Fallidos         : {len(fail_list)}")
    log(f"  ⏱ Tiempo total     : {format_elapsed(total_elapsed)}")
    if fail_list:
        log("\n  Fallidos:")
        for m in fail_list:
            log(f"    - {m}")
    log(f"\n  Log: {LOG_FILE}")
    log("="*60)


if __name__ == "__main__":
    main()