# =========================================================
# protein_prep.py
# FINAL STABLE VERSION
# =========================================================

import os
import subprocess


# =========================================================
# CLEAN PROTEIN
# =========================================================

def clean_protein(
    input_pdb,
    output_pdb
):

    with open(input_pdb, "r") as inp, open(output_pdb, "w") as out:

        for line in inp:

            if line.startswith(("ATOM", "HETATM")):

                out.write(line)

    print(f"Protein cleaned: {output_pdb}")

    return output_pdb


# =========================================================
# CONVERT TO PDBQT
# =========================================================

def protein_to_pdbqt(
    input_pdb,
    output_pdbqt
):

    subprocess.run([

        "obabel",
        input_pdb,
        "-O",
        output_pdbqt,
        "-xr",
        "--partialcharge",
        "gasteiger"

    ],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL)

    print(f"PDBQT created: {output_pdbqt}")

    return output_pdbqt