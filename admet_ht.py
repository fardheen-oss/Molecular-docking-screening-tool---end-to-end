# =========================================================
# admet_ht.py
# PRODUCTION VERSION — SwissADME + ADMETlab 3.0 Style
# =========================================================
#
# Descriptor engines:
#   - RDKit  : MW, logP, TPSA, HBA, HBD, RotBonds, QED,
#              PAINS, Lipinski, Veber, Egan, Muegge rules
#   - ADMET-AI: ML predictions trained on TDC/ChEMBL
#               BBB, AMES, hERG, DILI, CYP panel,
#               Bioavailability, Caco2, Solubility
#
# Filtering logic mirrors:
#   SwissADME  — Lipinski Ro5, Veber, Egan, PAINS
#   ADMETlab 3.0 — toxicity thresholds, absorption flags
#
# Output CSV columns:
#   Identity   : Ligand, SMILES
#   Physicochemical: MW, logP, TPSA, HBA, HBD,
#                    Rotatable, QED, Lipinski_Violations
#   Drug-likeness rules: Lipinski, Veber, Egan, Muegge
#   PAINS      : PAINS flag
#   Absorption : Caco2, Bioavailability
#   Distribution: BBB
#   Metabolism : CYP1A2, CYP2C19, CYP2C9, CYP2D6, CYP3A4
#   Toxicity   : AMES, hERG, DILI, Solubility
#   Decision   : Status, Fail_Reasons
#
# =========================================================

import os
import subprocess
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Crippen,
    Lipinski,
    QED,
    rdMolDescriptors,
    AllChem,
)
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from admet_ai import ADMETModel


# =========================================================
# LOAD ADMET-AI MODEL  (once at import time)
# =========================================================

print("\nLoading ADMET-AI model...")
_MODEL = ADMETModel()
print("ADMET-AI loaded.\n")


# =========================================================
# PAINS CATALOG  (once at import time)
# =========================================================

def _build_pains():
    p = FilterCatalogParams()
    p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(p)

_PAINS = _build_pains()


# =========================================================
# SAFE VALUE GETTER
# =========================================================

def _safe(preds, key, default=np.nan):
    """
    Safely extract a float prediction from ADMET-AI dict.
    Returns rounded float or default on any failure.
    """
    try:
        v = preds.get(key, default)
        return default if v is None else round(float(v), 4)
    except Exception:
        return default


# =========================================================
# PHYSICOCHEMICAL DESCRIPTORS
# =========================================================

def calculate_descriptors(mol):
    """
    Calculate all physicochemical descriptors from an RDKit mol.
    Returns dict or None if mol is invalid.
    """

    if mol is None:
        return None

    mw   = round(Descriptors.MolWt(mol), 3)
    logp = round(Crippen.MolLogP(mol), 3)
    tpsa = round(rdMolDescriptors.CalcTPSA(mol), 3)
    hba  = Lipinski.NumHAcceptors(mol)
    hbd  = Lipinski.NumHDonors(mol)
    rot  = Lipinski.NumRotatableBonds(mol)
    qed  = round(QED.qed(mol), 3)

    # ── PAINS alert ──────────────────────────────────────
    pains = int(_PAINS.HasMatch(mol))

    # ── Lipinski Ro5 (oral drug-likeness) ────────────────
    # Rule: ≤1 violation allowed
    lip_v = sum([
        mw   > 500,
        logp > 5,
        hba  > 10,
        hbd  > 5,
    ])

    # ── Veber rules (oral bioavailability proxy) ──────────
    # Rule: rotatable ≤10 AND TPSA ≤140
    veber_ok = int(rot <= 10 and tpsa <= 140)

    # ── Egan rules (passive intestinal absorption) ────────
    # Rule: logP ≤5.88 AND TPSA ≤131.6
    egan_ok = int(logp <= 5.88 and tpsa <= 131.6)

    # ── Muegge rules (lead-likeness) ─────────────────────
    # MW 200-600, logP -2 to 5, HBD ≤5, HBA ≤10,
    # TPSA ≤150, rot ≤15, rings ≤7, no fused rings >4
    rings = rdMolDescriptors.CalcNumRings(mol)
    muegge_ok = int(
        200  <= mw   <= 600  and
        -2   <= logp <= 5    and
        hbd  <= 5            and
        hba  <= 10           and
        tpsa <= 150          and
        rot  <= 15           and
        rings <= 7
    )

    return {
        "MW":                  mw,
        "logP":                logp,
        "TPSA":                tpsa,
        "HBA":                 hba,
        "HBD":                 hbd,
        "Rotatable":           rot,
        "QED":                 qed,
        "PAINS":               pains,
        "Lipinski_Violations": lip_v,
        "Veber_Pass":          veber_ok,
        "Egan_Pass":           egan_ok,
        "Muegge_Pass":         muegge_ok,
    }


# =========================================================
# LIGAND PREPARATION  (SMILES → PDBQT)
# =========================================================

def prepare_ligand(smiles, name, output_dir):
    """
    3D-embed → MMFF optimise → write SDF → convert to PDBQT.
    Returns path to PDBQT or None on failure.
    """

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol = Chem.AddHs(mol)

        res = AllChem.EmbedMolecule(mol, randomSeed=42)
        if res == -1:
            # Try ETKDG with more attempts
            params = AllChem.EmbedParameters()
            params.randomSeed = 42
            params.maxIterations = 1000
            AllChem.EmbedMolecule(mol, params)

        AllChem.MMFFOptimizeMolecule(mol)

        sdf_path   = os.path.join(output_dir, f"{name}.sdf")
        pdbqt_path = os.path.join(output_dir, f"{name}.pdbqt")

        with Chem.SDWriter(sdf_path) as w:
            w.write(mol)

        result = subprocess.run(
            ["obabel", sdf_path, "-O", pdbqt_path,
             "--partialcharge", "gasteiger"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if not os.path.exists(pdbqt_path) or os.path.getsize(pdbqt_path) == 0:
            print(f"    [WARN] obabel failed for {name}: "
                  f"{result.stderr.decode().strip()[:120]}")
            return None

        return pdbqt_path

    except Exception as e:
        print(f"    [WARN] Ligand prep exception for {name}: {e}")
        return None


# =========================================================
# PASS / FAIL DECISION ENGINE
# =========================================================
#
# Thresholds are drawn from:
#   SwissADME  — Lipinski, Veber, Egan, PAINS
#   ADMETlab 3.0 — recommended probability cutoffs
#       hERG  inhibition > 0.5  → cardiotoxic risk
#       AMES  mutagenicity > 0.5 → mutagenic
#       DILI  > 0.7             → liver injury risk
#       BBB   note: low BBB is fine for non-CNS targets
#       Bioavailability_Ma > 0.5 → orally bioavailable
#       Solubility < -6         → practically insoluble
# =========================================================

def _decide(r):
    """
    Apply multi-filter decision.
    Returns (status, [reasons]).
    """

    reasons = []

    # ── Hard physicochemical filters ──────────────────────
    if r["MW"] > 500:
        reasons.append("MW > 500")
    if r["logP"] > 5:
        reasons.append("logP > 5")
    if r["TPSA"] > 140:
        reasons.append("TPSA > 140")
    if r["HBD"] > 5:
        reasons.append("HBD > 5")
    if r["HBA"] > 10:
        reasons.append("HBA > 10")
    if r["Rotatable"] > 10:
        reasons.append("Rotatable > 10")

    # ── PAINS structural alert ────────────────────────────
    if r["PAINS"] > 0:
        reasons.append("PAINS Alert")

    # ── ML toxicity filters ───────────────────────────────
    ames = r.get("AMES", np.nan)
    if not np.isnan(ames) and ames > 0.5:
        reasons.append(f"AMES mutagenic ({ames:.2f})")

    herg = r.get("hERG", np.nan)
    if not np.isnan(herg) and herg > 0.5:
        reasons.append(f"hERG cardiotoxic ({herg:.2f})")

    dili = r.get("DILI", np.nan)
    if not np.isnan(dili) and dili > 0.7:
        reasons.append(f"DILI hepatotoxic ({dili:.2f})")

    # ── Solubility filter ─────────────────────────────────
    sol = r.get("Solubility", np.nan)
    if not np.isnan(sol) and sol < -6:
        reasons.append(f"Practically insoluble (logS {sol:.2f})")

    status = "FAIL" if reasons else "PASS"
    return status, reasons


# =========================================================
# MAIN ADMET PIPELINE
# =========================================================

def run_admet(ligand_csv, output_folder):
    """
    Full ADMET screening pipeline.
    Returns list of PDBQT paths for passing ligands.
    """

    os.makedirs(output_folder, exist_ok=True)

    ligands_folder = os.path.join(output_folder, "prepared_ligands")
    os.makedirs(ligands_folder, exist_ok=True)

    output_csv = os.path.join(output_folder, "admet_results.csv")

    df = pd.read_csv(ligand_csv)

    final_results  = []
    passed_ligands = []

    for _, row in df.iterrows():

        name   = str(row["Ligand"]).strip()
        smiles = str(row["SMILES"]).strip()

        print(f"  Processing: {name}")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"    [SKIP] Invalid SMILES for {name}")
            continue

        # ── Physicochemical descriptors ───────────────────
        desc = calculate_descriptors(mol)
        if desc is None:
            print(f"    [SKIP] Descriptor calculation failed for {name}")
            continue

        # ── ADMET-AI ML predictions ───────────────────────
        try:
            preds = _MODEL.predict(smiles=smiles)
        except Exception as e:
            print(f"    [WARN] ADMET-AI prediction failed for {name}: {e}")
            preds = {}

        # ── Assemble result row ───────────────────────────
        result = {

            # Identity
            "Ligand":  name,
            "SMILES":  smiles,

            # Physicochemical
            "MW":                  desc["MW"],
            "logP":                desc["logP"],
            "TPSA":                desc["TPSA"],
            "HBA":                 desc["HBA"],
            "HBD":                 desc["HBD"],
            "Rotatable":           desc["Rotatable"],
            "QED":                 desc["QED"],

            # Drug-likeness rules
            "Lipinski_Violations": desc["Lipinski_Violations"],
            "Veber_Pass":          desc["Veber_Pass"],
            "Egan_Pass":           desc["Egan_Pass"],
            "Muegge_Pass":         desc["Muegge_Pass"],

            # Structural alerts
            "PAINS":               desc["PAINS"],

            # ── ADMET-AI ML predictions ───────────────────
            # Absorption
            "Caco2_Permeability":  _safe(preds, "Caco2_Wang"),
            "Bioavailability":     _safe(preds, "Bioavailability_Ma"),

            # Distribution
            "BBB_Penetration":     _safe(preds, "BBB_Martins"),

            # Metabolism — CYP inhibition probability
            # (>0.5 = likely inhibitor of that isoform)
            "CYP1A2_Inhibitor":    _safe(preds, "CYP1A2_Veith"),
            "CYP2C19_Inhibitor":   _safe(preds, "CYP2C19_Veith"),
            "CYP2C9_Inhibitor":    _safe(preds, "CYP2C9_Veith"),
            "CYP2D6_Inhibitor":    _safe(preds, "CYP2D6_Veith"),
            "CYP3A4_Inhibitor":    _safe(preds, "CYP3A4_Veith"),

            # Toxicity
            "AMES_Mutagenicity":   _safe(preds, "AMES"),
            "hERG_Cardiotoxicity": _safe(preds, "hERG"),
            "DILI_Hepatotoxicity": _safe(preds, "DILI"),
            "Solubility_logS":     _safe(preds, "Solubility_AqSolDB"),
        }

        # ── Pass/Fail decision ────────────────────────────
        status, reasons = _decide(result)
        result["Status"]       = status
        result["Fail_Reasons"] = "; ".join(reasons) if reasons else "—"

        final_results.append(result)

        # ── Prepare PDBQT for passing ligands ─────────────
        if status == "PASS":
            pdbqt = prepare_ligand(smiles, name, ligands_folder)
            if pdbqt:
                passed_ligands.append(pdbqt)
            else:
                print(f"    [WARN] PDBQT preparation failed for {name} "
                      f"— excluded from docking")

    # ── Save CSV ──────────────────────────────────────────
    out_df = pd.DataFrame(final_results)
    out_df.to_csv(output_csv, index=False)

    # ── Summary ───────────────────────────────────────────
    n_pass = sum(1 for r in final_results if r["Status"] == "PASS")
    n_fail = len(final_results) - n_pass

    print("\n" + "=" * 50)
    print("ADMET SCREENING COMPLETE")
    print("=" * 50)
    print(f"  Total   : {len(final_results)}")
    print(f"  PASS    : {n_pass}")
    print(f"  FAIL    : {n_fail}")
    print(f"  Results : {output_csv}")
    print("=" * 50 + "\n")

    return passed_ligands


# =========================================================
# WRAPPER (called by main_pipeline.py)
# =========================================================

def run_admet_pipeline(ligand_csv, output_dir):
    return run_admet(ligand_csv, output_dir)


# =========================================================
# CLI
# =========================================================

if __name__ == "__main__":

    import sys

    if len(sys.argv) != 3:
        print(
            "\nUsage:\n"
            "  python admet_ht.py ligands.csv output_folder\n\n"
            "Example:\n"
            "  python admet_ht.py ligands.csv results/admet\n"
        )
        sys.exit(1)

    run_admet(sys.argv[1], sys.argv[2])