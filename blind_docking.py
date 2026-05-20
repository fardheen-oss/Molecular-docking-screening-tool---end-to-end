# =========================================================
# blind_docking.py
# PRODUCTION VERSION — Discovery Studio Compatible
# =========================================================
#
# Grid strategy:
#   The search box covers the ENTIRE protein (bounding box
#   + 10 Å padding, capped at 126 Å per Vina limit).
#   The centre (cx, cy, cz) is the geometric centre of all
#   receptor atoms — identical for every ligand in a run.
#
# Per-ligand outputs:
#   {name}_blind.pdbqt   — Vina output (all poses)
#   {name}_blind.pdb     — Best pose only (MODEL 1, clean PDB)
#   {name}_complex.pdb   — Receptor + ligand (DS-ready)
#
# CSV columns (blind_docking_results_exhaustiveness_N.csv):
#   Ligand, Binding_Affinity,
#   Center_X, Center_Y, Center_Z,
#   Size_X, Size_Y, Size_Z
#   (no file-path columns)
# =========================================================

import os
import subprocess
import pandas as pd


# =========================================================
# WHOLE-PROTEIN GRID
# =========================================================

def get_protein_grid(pdbqt_file):
    """
    Compute bounding-box centre and grid size from receptor PDBQT.

    Because this covers the entire protein surface, Center_X/Y/Z
    is the same number for every ligand docked against the same
    receptor — they all share the same search space.

    Size is capped at 126 Å on each axis (AutoDock Vina hard limit).
    Returns (center_list, size_list).
    """

    coords = []

    with open(pdbqt_file, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append((
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ))
                except ValueError:
                    pass

    if not coords:
        raise RuntimeError(
            f"No ATOM/HETATM coordinates found in: {pdbqt_file}"
        )

    mins   = [min(c[i] for c in coords) for i in range(3)]
    maxs   = [max(c[i] for c in coords) for i in range(3)]
    center = [round((mins[i] + maxs[i]) / 2, 3) for i in range(3)]
    size   = [min(126, round((maxs[i] - mins[i]) + 30, 3)) for i in range(3)]

    return center, size


# =========================================================
# SCORE EXTRACTOR
# =========================================================

def extract_vina_score(pdbqt_file):
    """Return best binding affinity (float) from Vina PDBQT, or None."""

    if not os.path.exists(pdbqt_file):
        return None

    with open(pdbqt_file, "r") as f:
        for line in f:
            if "REMARK VINA RESULT:" in line:
                try:
                    return float(line.split()[3])
                except (IndexError, ValueError):
                    return None
    return None


# =========================================================
# POSE CENTRE  (used by active refinement / P2Rank fallback)
# =========================================================

def get_pose_center(pdbqt_file):
    """
    Geometric centre of MODEL 1 atoms in a Vina PDBQT.
    Returns (cx, cy, cz) or (0, 0, 0) on failure.
    """

    coords    = []
    in_model1 = False

    if not os.path.exists(pdbqt_file):
        return (0.0, 0.0, 0.0)

    with open(pdbqt_file, "r") as f:
        for line in f:
            if line.startswith("MODEL") and "1" in line:
                in_model1 = True
                continue
            if in_model1 and line.startswith("ENDMDL"):
                break
            if in_model1 and line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append((
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ))
                except ValueError:
                    pass

    if not coords:
        return (0.0, 0.0, 0.0)

    n = len(coords)
    return (
        round(sum(c[0] for c in coords) / n, 3),
        round(sum(c[1] for c in coords) / n, 3),
        round(sum(c[2] for c in coords) / n, 3),
    )


# =========================================================
# BEST POSE → PDB  (MODEL 1 only)
# =========================================================

def extract_best_pose_to_pdb(pdbqt_file, output_pdb):
    """
    Extract MODEL 1 from a Vina PDBQT and convert to PDB via obabel.
    Returns True on success.
    """

    temp      = output_pdb.replace(".pdb", "_pose1_tmp.pdbqt")
    in_model1 = False

    with open(pdbqt_file, "r") as fin, open(temp, "w") as fout:
        for line in fin:
            if line.startswith("MODEL") and "1" in line:
                in_model1 = True
                continue
            if in_model1 and line.startswith("ENDMDL"):
                break
            if in_model1:
                fout.write(line)

    if not os.path.exists(temp) or os.path.getsize(temp) == 0:
        print(f"  [ERROR] MODEL 1 extraction failed from {pdbqt_file}")
        return False

    result = subprocess.run(
        ["obabel", temp, "-O", output_pdb],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    os.remove(temp)

    if not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
        print(f"  [ERROR] obabel conversion failed: "
              f"{result.stderr.decode().strip()[:120]}")
        return False

    return True


# =========================================================
# COMPLEX BUILDER  (Discovery Studio compatible)
# =========================================================

def build_complex(receptor_pdb, ligand_pdb, output_complex):
    """
    Write receptor + ligand as a strict 80-column PDB complex.
    Ligand atoms → HETATM, chain X, residue LIG 999.
    Returns True on success.
    """

    for fpath, label in ((receptor_pdb, "Receptor"), (ligand_pdb, "Ligand")):
        if not os.path.exists(fpath):
            print(f"  [ERROR] {label} PDB missing: {fpath}")
            return False

    serial = 9001

    with open(output_complex, "w") as out:

        # Receptor
        with open(receptor_pdb, "r") as rec:
            for line in rec:
                if line.startswith(("ATOM", "HETATM")):
                    out.write(line if line.endswith("\n") else line + "\n")
        out.write("TER\n")

        # Ligand — rewritten as HETATM with strict column layout
        with open(ligand_pdb, "r") as lig:
            for line in lig:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                try:
                    x    = float(line[30:38])
                    y    = float(line[38:46])
                    z    = float(line[46:54])
                    occ  = float(line[54:60]) if len(line) > 54 else 1.00
                    bfac = float(line[60:66]) if len(line) > 60 else 0.00
                except ValueError:
                    x, y, z, occ, bfac = 0.0, 0.0, 0.0, 1.00, 0.00

                raw  = line[12:16].strip()
                name = f" {raw:<3s}" if len(raw) < 4 else raw[:4]
                elem = line[76:78].strip() if len(line) >= 78 else raw[:1]

                out.write(
                    f"HETATM{serial:5d} {name} LIG X{999:4d}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"{occ:6.2f}{bfac:6.2f}          {elem:>2s}\n"
                )
                serial += 1

        out.write("END\n")

    if not os.path.exists(output_complex) or os.path.getsize(output_complex) == 0:
        print(f"  [ERROR] Complex file empty: {output_complex}")
        return False

    return True


# =========================================================
# MAIN: BLIND DOCKING  (standalone entry point)
# =========================================================

def run_blind_docking(
    receptor_pdb,
    ligand_folder,
    output_folder,
    exhaustiveness=2,
):
    """
    Blind docking of all PDBQT ligands in ligand_folder.
    Uses whole-protein bounding box — same grid for every ligand.

    Returns path to results CSV.
    """

    os.makedirs(output_folder, exist_ok=True)

    # Receptor PDBQT
    receptor_pdbqt = receptor_pdb.replace(".pdb", ".pdbqt")
    if not os.path.exists(receptor_pdbqt):
        print("Converting receptor PDB → PDBQT ...")
        result = subprocess.run(
            ["obabel", receptor_pdb, "-O", receptor_pdbqt,
             "-xr", "--partialcharge", "gasteiger"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        if not os.path.exists(receptor_pdbqt):
            raise RuntimeError(
                f"Receptor conversion failed: {result.stderr.decode().strip()}"
            )

    center, size = get_protein_grid(receptor_pdbqt)
    print(f"Whole-protein grid centre : {center}")
    print(f"Grid size (Å)             : {size}")
    print(f"(All ligands share this same search box)\n")

    ligands = sorted(f for f in os.listdir(ligand_folder) if f.endswith(".pdbqt"))
    if not ligands:
        raise RuntimeError(f"No PDBQT ligands found in: {ligand_folder}")

    results = []

    for fname in ligands:

        name   = fname.replace(".pdbqt", "")
        lpath  = os.path.join(ligand_folder, fname)
        print(f"Docking: {name}")

        out_pdbqt   = os.path.join(output_folder, f"{name}_blind.pdbqt")
        out_pdb     = os.path.join(output_folder, f"{name}_blind.pdb")
        out_complex = os.path.join(output_folder, f"{name}_complex.pdb")

        # Vina
        vina_result = subprocess.run(
            [
                "vina",
                "--receptor", receptor_pdbqt,
                "--ligand",   lpath,
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x",   str(size[0]),
                "--size_y",   str(size[1]),
                "--size_z",   str(size[2]),
                "--exhaustiveness", str(exhaustiveness),
                "--cpu", "1",
                "--out", out_pdbqt,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if not os.path.exists(out_pdbqt):
            print(f"  [FAIL] Vina error: {vina_result.stderr.decode().strip()[:120]}")
            continue

        score = extract_vina_score(out_pdbqt)
        print(f"  Score: {score} kcal/mol")

        if not extract_best_pose_to_pdb(out_pdbqt, out_pdb):
            print(f"  [FAIL] PDB extraction failed for {name}")
            continue
        print(f"  PDB     → {out_pdb}")

        if not build_complex(receptor_pdb, out_pdb, out_complex):
            print(f"  [FAIL] Complex build failed for {name}")
            continue
        print(f"  Complex → {out_complex}")

        # No file paths in CSV
        results.append({
            "Ligand":           name,
            "Binding_Affinity": score,
            "Center_X":         center[0],
            "Center_Y":         center[1],
            "Center_Z":         center[2],
            "Size_X":           size[0],
            "Size_Y":           size[1],
            "Size_Z":           size[2],
        })

    csv_path = os.path.join(
        output_folder,
        f"blind_docking_results_exhaustiveness_{exhaustiveness}.csv",
    )
    pd.DataFrame(results).to_csv(csv_path, index=False)

    print("\n" + "=" * 50)
    print("BLIND DOCKING COMPLETE")
    print("=" * 50)
    print(f"Results → {csv_path}")

    return csv_path