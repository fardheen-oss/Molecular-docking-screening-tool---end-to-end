# =========================================================
# main_pipeline.py
# PRODUCTION VERSION — P2Rank Active Site Pipeline
# =========================================================

import os
import argparse
import pandas as pd

from protein_prep import (
    clean_protein,
    protein_to_pdbqt,
)
from admet_ht import (
    run_admet_pipeline,
)
from blind_docking import (
    get_protein_grid,
    extract_vina_score,
    get_pose_center,
    extract_best_pose_to_pdb,
    build_complex,
)
from docking import (
    run_vina,
)
from site_predictor import (
    predict_binding_sites,
    pocket_grid_size,
    save_site_summary,
)
from active_site_refine import (
    run_active_refinement,
)


def setup_folders(base_dir):
    folders = {
        "admet":  os.path.join(base_dir, "admet"),
        "blind":  os.path.join(base_dir, "blind_docking"),
        "p2rank": os.path.join(base_dir, "p2rank"),
        "active": os.path.join(base_dir, "active_refinement"),
    }
    for path in folders.values():
        os.makedirs(path, exist_ok=True)
    return folders


def main():

    parser = argparse.ArgumentParser(
        description="Molecular docking pipeline with P2Rank site prediction"
    )
    parser.add_argument("--protein", required=True, help="Input protein PDB file")
    parser.add_argument("--ligands", required=True, help="Ligand CSV with columns: Ligand, SMILES")
    parser.add_argument(
        "--p2rank", default=None,
        help=(
            "Path to P2Rank install directory "
            "(e.g. C:\\tools\\p2rank or ~/tools/p2rank). "
            "Falls back to P2RANK_HOME env variable. "
            "If neither set, active refinement uses blind-pose centres."
        )
    )
    args = parser.parse_args()

    if not os.path.exists(args.protein):
        raise FileNotFoundError(f"Protein file not found: {args.protein}")
    if not os.path.exists(args.ligands):
        raise FileNotFoundError(f"Ligand CSV not found: {args.ligands}")

    # ── Exhaustiveness ─────────────────────────────────────
    print("\n===================================")
    print("DOCKING MODE")
    print("===================================")
    print("1 = Faster   (Exhaustiveness 2)")
    print("2 = Accurate (Exhaustiveness 8)")
    print("===================================")
    choice               = input("Select option (1 or 2): ").strip()
    blind_exhaustiveness = 8 if choice == "2" else 2
    active_exhaustiveness = 8

    result_dir = "results"
    folders    = setup_folders(result_dir)

    # =========================================================
    # STAGE 1: PROTEIN PREPARATION
    # =========================================================

    print("\n===================================")
    print("STAGE 1: PROTEIN PREPARATION")
    print("===================================\n")

    clean_pdb      = os.path.join(result_dir, "clean_protein.pdb")
    receptor_pdbqt = os.path.join(result_dir, "receptor.pdbqt")

    clean_protein(args.protein, clean_pdb)
    protein_to_pdbqt(clean_pdb, receptor_pdbqt)

    for f in (clean_pdb, receptor_pdbqt):
        if not os.path.exists(f):
            raise RuntimeError(f"Protein preparation failed — missing: {f}")

    print(f"  Clean PDB → {clean_pdb}")
    print(f"  PDBQT     → {receptor_pdbqt}")

    # =========================================================
    # STAGE 2: ADMET SCREENING
    # =========================================================

    print("\n===================================")
    print("STAGE 2: ADMET SCREENING")
    print("===================================\n")

    passed_ligands = run_admet_pipeline(
        ligand_csv=args.ligands,
        output_dir=folders["admet"],
    )
    prepared_ligands_folder = os.path.join(folders["admet"], "prepared_ligands")

    if not passed_ligands:
        print("\n[STOP] No ligands passed ADMET filters.\n")
        return

    print(f"\n  {len(passed_ligands)} ligand(s) passed ADMET.")

    # =========================================================
    # STAGE 3: BLIND DOCKING
    # Whole-protein grid — identical centre and size for ALL ligands
    # Per ligand: _blind.pdbqt  _blind.pdb  _complex.pdb
    # CSV columns: Ligand, Binding_Affinity, Center_X/Y/Z, Size_X/Y/Z
    # =========================================================

    print("\n===================================")
    print("STAGE 3: BLIND DOCKING")
    print("===================================\n")

    center, size = get_protein_grid(receptor_pdbqt)
    print(f"  Whole-protein grid centre : {center}")
    print(f"  Grid size (Å)             : {size}")
    print(f"  (All {len(passed_ligands)} ligands share this same search box)\n")

    blind_results = []

    for ligand_path in passed_ligands:

        name = os.path.basename(ligand_path).replace(".pdbqt", "")
        print(f"  Docking: {name}")

        out_pdbqt   = os.path.join(folders["blind"], f"{name}_blind.pdbqt")
        out_pdb     = os.path.join(folders["blind"], f"{name}_blind.pdb")
        out_complex = os.path.join(folders["blind"], f"{name}_complex.pdb")

        run_vina(
            receptor_pdbqt, ligand_path, out_pdbqt,
            center=center, size=size,
            exhaustiveness=blind_exhaustiveness,
        )

        if not os.path.exists(out_pdbqt):
            print(f"    [FAIL] Vina produced no output for {name}")
            continue

        score = extract_vina_score(out_pdbqt)
        if score is None:
            print(f"    [WARNING] No score parsed for {name} — skipped")
            continue
        print(f"    Score: {score} kcal/mol")

        if not extract_best_pose_to_pdb(out_pdbqt, out_pdb):
            print(f"    [FAIL] PDB extraction failed for {name}")
            continue
        print(f"    PDB     → {out_pdb}")

        if not build_complex(clean_pdb, out_pdb, out_complex):
            print(f"    [FAIL] Complex build failed for {name}")
            continue
        print(f"    Complex → {out_complex}")

        blind_results.append({
            "Ligand":           name,
            "Binding_Affinity": score,
            "Center_X":         center[0],
            "Center_Y":         center[1],
            "Center_Z":         center[2],
            "Size_X":           size[0],
            "Size_Y":           size[1],
            "Size_Z":           size[2],
        })

    if not blind_results:
        print("\n[STOP] No ligands docked successfully.\n")
        return

    blind_csv = os.path.join(
        folders["blind"],
        f"blind_docking_exhaustiveness_{blind_exhaustiveness}.csv",
    )
    pd.DataFrame(blind_results).to_csv(blind_csv, index=False)
    print(f"\n  Blind docking results → {blind_csv}")

    # =========================================================
    # STAGE 4: TOP-3 SELECTION
    # =========================================================

    print("\n===================================")
    print("STAGE 4: TOP-3 SELECTION")
    print("===================================\n")

    sorted_hits  = sorted(blind_results, key=lambda x: x["Binding_Affinity"])
    top_3        = sorted_hits[:3]
    top_hits_csv = os.path.join(result_dir, "top_hits.csv")
    top_df       = pd.DataFrame(top_3)
    top_df.insert(0, "Rank", range(1, len(top_df) + 1))
    top_df.to_csv(top_hits_csv, index=False)

    print(top_df[["Rank", "Ligand", "Binding_Affinity"]].to_string(index=False))
    print(f"\n  Top hits saved → {top_hits_csv}")

    # =========================================================
    # STAGE 5: P2RANK BINDING SITE PREDICTION
    #
    # P2Rank scores every point on the protein surface using a
    # trained Random Forest model — no ligand info used at all.
    # This gives a protein-based pocket that is scientifically
    # independent of which ligands you are testing.
    #
    # If P2Rank is not installed the pipeline falls back to the
    # geometric centre of each ligand's best blind pose.
    #
    # Install P2Rank:
    #   https://github.com/rdk/p2rank/releases
    #   Then run with:  --p2rank C:\tools\p2rank
    #   or set env var: P2RANK_HOME=C:\tools\p2rank
    # =========================================================

    print("\n===================================")
    print("STAGE 5: P2RANK SITE PREDICTION")
    print("===================================\n")

    pockets = predict_binding_sites(
        receptor_pdb=clean_pdb,
        output_folder=folders["p2rank"],
        p2rank_dir=args.p2rank,
        top_n=3,
    )

    p2rank_available = len(pockets) > 0

    if p2rank_available:
        save_site_summary(pockets, folders["p2rank"])
        best_pocket   = pockets[0]
        active_center = (best_pocket["center_x"], best_pocket["center_y"], best_pocket["center_z"])
        gs            = pocket_grid_size(best_pocket["score"])
        active_size   = gs

        print(f"\n  Best pocket  : Pocket {best_pocket['pocket']}")
        print(f"  P2Rank score : {best_pocket['score']}")
        print(f"  Centre       : {active_center}")
        print(f"  Grid size    : {active_size}")
        print(f"  Lining res.  : {best_pocket['residues'][:80]}")
        print("\n  All top-3 ligands will refine into this pocket.")

        # Store pocket info in top_hits.csv for reference
        top_df["Pocket_Center_X"] = active_center[0]
        top_df["Pocket_Center_Y"] = active_center[1]
        top_df["Pocket_Center_Z"] = active_center[2]
        top_df["Pocket_Grid"]     = active_size[0]
        top_df.to_csv(top_hits_csv, index=False)

    else:
        print(
            "\n  [FALLBACK] P2Rank not available.\n"
            "  Active refinement will use each ligand's own\n"
            "  blind-docking best-pose centre instead.\n"
            "\n  To enable P2Rank:\n"
            "    1. Download: https://github.com/rdk/p2rank/releases\n"
            "    2. Run with: python main_pipeline.py --protein X.pdb "
            "--ligands Y.csv --p2rank C:\\tools\\p2rank\n"
            "    OR set:      set P2RANK_HOME=C:\\tools\\p2rank\n"
        )
        active_center = None
        active_size   = None

    # =========================================================
    # STAGE 6: ACTIVE-SITE REFINEMENT
    # Per top-3 ligand: _active.pdbqt  _active.pdb  _complex.pdb
    # CSV: Rank, Ligand, Binding_Affinity, Center_X/Y/Z, Grid_Size
    # =========================================================

    print("\n===================================")
    print("STAGE 6: ACTIVE-SITE REFINEMENT")
    print("===================================")

    run_active_refinement(
        receptor_pdbqt=receptor_pdbqt,
        receptor_pdb=clean_pdb,
        top_hits_csv=top_hits_csv,
        prepared_ligands_folder=prepared_ligands_folder,
        blind_docking_folder=folders["blind"],
        output_folder=folders["active"],
        exhaustiveness=active_exhaustiveness,
        p2rank_center=active_center,
        p2rank_size=active_size,
    )

    # =========================================================
    # DONE
    # =========================================================

    print("\n===================================")
    print("PIPELINE COMPLETE")
    print("===================================\n")
    print(f"All results in: {result_dir}/\n")
    print("  admet/")
    print("    admet_results.csv")
    print("    prepared_ligands/  {name}.pdbqt")
    print()
    print("  blind_docking/   — all ADMET-passing ligands")
    print("    {name}_blind.pdbqt    Vina all poses")
    print("    {name}_blind.pdb      Best pose PDB")
    print("    {name}_complex.pdb    DS-ready complex")
    print(f"    blind_docking_exhaustiveness_{blind_exhaustiveness}.csv")
    print()
    if p2rank_available:
        print("  p2rank/")
        print("    predicted_sites.csv   Pocket list with scores and centres")
        print()
    print("  active_refinement/   — top-3 ligands only")
    print("    {name}_active.pdbqt   Vina all poses")
    print("    {name}_active.pdb     Best pose PDB")
    print("    {name}_complex.pdb    DS-ready complex")
    print("    active_refinement_results.csv")
    print()
    if p2rank_available:
        print(f"  Active site: P2Rank Pocket 1 @ {active_center} (score {best_pocket['score']})")
    else:
        print("  Active site: per-ligand blind-pose centre (P2Rank not used)")


if __name__ == "__main__":
    main()