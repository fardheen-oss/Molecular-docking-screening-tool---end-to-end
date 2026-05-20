# =========================================================
# active_site_refine.py
# PRODUCTION VERSION — P2Rank + Discovery Studio Compatible
# =========================================================

import os
import subprocess
import pandas as pd


def run_vina(receptor_pdbqt, ligand_pdbqt, output_pdbqt, center, size=(20, 20, 20), exhaustiveness=8, num_modes=9):
    cmd = [
        "vina",
        "--receptor", receptor_pdbqt,
        "--ligand",   ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x",   str(size[0]),
        "--size_y",   str(size[1]),
        "--size_z",   str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--cpu", "1",
        "--out", output_pdbqt,
    ]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"  [VINA ERROR] {result.stderr.decode().strip()}")
        return False
    if not os.path.exists(output_pdbqt):
        print(f"  [ERROR] Vina produced no output: {output_pdbqt}")
        return False
    return True


def extract_vina_score(pdbqt_file):
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


def get_pose_center_from_pdbqt(pdbqt_file):
    coords = []
    in_model1 = False
    with open(pdbqt_file, "r") as f:
        for line in f:
            if line.startswith("MODEL") and "1" in line:
                in_model1 = True
                continue
            if in_model1 and line.startswith("ENDMDL"):
                break
            if in_model1 and line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    pass
    if not coords:
        return None
    n = len(coords)
    return (
        round(sum(c[0] for c in coords) / n, 3),
        round(sum(c[1] for c in coords) / n, 3),
        round(sum(c[2] for c in coords) / n, 3),
    )


def extract_best_pose_to_pdb(pdbqt_file, output_pdb):
    temp = output_pdb.replace(".pdb", "_pose1_tmp.pdbqt")
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
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    os.remove(temp)
    if not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
        print(f"  [ERROR] obabel PDB conversion failed: {result.stderr.decode().strip()[:120]}")
        return False
    return True


def build_complex(receptor_pdb, ligand_pdb, output_complex):
    for fpath, label in ((receptor_pdb, "Receptor"), (ligand_pdb, "Ligand")):
        if not os.path.exists(fpath):
            print(f"  [ERROR] {label} PDB missing: {fpath}")
            return False
    serial = 9001
    with open(output_complex, "w") as out:
        with open(receptor_pdb, "r") as rec:
            for line in rec:
                if line.startswith(("ATOM", "HETATM")):
                    out.write(line if line.endswith("\n") else line + "\n")
        out.write("TER\n")
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
        print(f"  [ERROR] Complex file is empty: {output_complex}")
        return False
    return True


def run_active_refinement(
    receptor_pdbqt,
    receptor_pdb,
    top_hits_csv,
    prepared_ligands_folder,
    blind_docking_folder,
    output_folder,
    exhaustiveness=8,
    p2rank_center=None,
    p2rank_size=None,
):
    """
    Docking centre priority:
      1. p2rank_center  — protein-based pocket (best, all ligands share same box)
      2. Geometric centre of MODEL 1 from blind docking PDBQT
      3. Center_X/Y/Z from top_hits.csv (whole-protein fallback)

    CSV columns: Rank, Ligand, Binding_Affinity, Center_X, Center_Y, Center_Z, Grid_Size
    No file-path columns in CSV.
    """

    os.makedirs(output_folder, exist_ok=True)

    for fpath in (receptor_pdbqt, receptor_pdb, top_hits_csv):
        if not os.path.exists(fpath):
            raise FileNotFoundError(f"Required file not found: {fpath}")

    top_df = pd.read_csv(top_hits_csv)
    if top_df.empty:
        print("[WARNING] top_hits.csv is empty — nothing to refine.")
        return None

    print(f"\nActive-site refinement: {len(top_df)} ligand(s)")
    print("=" * 50)

    if p2rank_center is not None:
        print(f"  Mode     : P2Rank pocket (all ligands share same centre)")
        print(f"  Centre   : {p2rank_center}")
        print(f"  Grid     : {p2rank_size}")
    else:
        print(f"  Mode     : Per-ligand blind-pose centre (P2Rank fallback)")

    refined_results = []

    for _, row in top_df.iterrows():

        name = str(row["Ligand"]).strip()
        rank = int(row.get("Rank", 0))

        print(f"\n[Rank {rank}] {name}")

        ligand_pdbqt = os.path.join(prepared_ligands_folder, f"{name}.pdbqt")
        if not os.path.exists(ligand_pdbqt):
            print(f"  [SKIP] Prepared ligand not found: {ligand_pdbqt}")
            continue

        # ── Centre and size selection ──────────────────────
        if p2rank_center is not None:
            center = p2rank_center
            size   = p2rank_size if p2rank_size else (20, 20, 20)
            print(f"  [INFO] Using P2Rank centre: {center}")
        else:
            blind_out = os.path.join(blind_docking_folder, f"{name}_blind.pdbqt")
            center = get_pose_center_from_pdbqt(blind_out) if os.path.exists(blind_out) else None
            size   = (20, 20, 20)
            if center is None:
                try:
                    center = (
                        float(row["Center_X"]),
                        float(row["Center_Y"]),
                        float(row["Center_Z"]),
                    )
                    print("  [INFO] Fallback: whole-protein centre from CSV")
                except (KeyError, ValueError):
                    print(f"  [SKIP] Cannot determine docking centre for {name}")
                    continue
            else:
                print(f"  [INFO] Centre from best blind pose: {center}")

        # ── Output paths ──────────────────────────────────
        out_pdbqt   = os.path.join(output_folder, f"{name}_active.pdbqt")
        out_pdb     = os.path.join(output_folder, f"{name}_active.pdb")
        out_complex = os.path.join(output_folder, f"{name}_complex.pdb")

        # ── Docking ───────────────────────────────────────
        print(f"  Docking → {out_pdbqt}")
        if not run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt,
                        center=center, size=size, exhaustiveness=exhaustiveness):
            print(f"  [FAIL] Docking failed for {name}")
            continue

        score = extract_vina_score(out_pdbqt)
        print(f"  Score   : {score} kcal/mol")

        # ── Best pose → PDB ───────────────────────────────
        if not extract_best_pose_to_pdb(out_pdbqt, out_pdb):
            print(f"  [FAIL] PDB extraction failed for {name}")
            continue
        print(f"  PDB     → {out_pdb}")

        # ── Complex ───────────────────────────────────────
        if not build_complex(receptor_pdb, out_pdb, out_complex):
            print(f"  [FAIL] Complex build failed for {name}")
            continue
        print(f"  Complex → {out_complex}")

        refined_results.append({
            "Rank":             rank,
            "Ligand":           name,
            "Binding_Affinity": score,
            "Center_X":         center[0],
            "Center_Y":         center[1],
            "Center_Z":         center[2],
            "Grid_Size":        size[0],
        })

    out_csv = os.path.join(output_folder, "active_refinement_results.csv")
    if refined_results:
        pd.DataFrame(refined_results).to_csv(out_csv, index=False)
        print(f"\nResults → {out_csv}")
    else:
        print("\n[WARNING] No ligands completed active-site refinement.")

    print("\n" + "=" * 50)
    print("ACTIVE REFINEMENT COMPLETE")
    print("=" * 50)

    return out_csv if refined_results else None