import os
import subprocess


def extract_best_pose(pdbqt_file, output_pdb):
    """
    Extract ONLY MODEL 1 from vina output
    and convert it cleanly into PDB.
    """

    temp_pose = output_pdb.replace(".pdb", "_pose1.pdbqt")

    writing = False

    with open(pdbqt_file, "r") as fin, open(temp_pose, "w") as fout:

        for line in fin:

            if line.startswith("MODEL 1"):
                writing = True
                continue

            if writing:

                if line.startswith("ENDMDL"):
                    break

                fout.write(line)

    subprocess.run([
        "obabel",
        temp_pose,
        "-O",
        output_pdb
    ],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL)

    os.remove(temp_pose)

    return output_pdb


def build_best_complex(
    receptor_pdb,
    docked_pdbqt,
    output_complex
):
    """
    Create Discovery-Studio-compatible complex.
    """

    ligand_pdb = output_complex.replace(".pdb", "_ligand.pdb")

    extract_best_pose(
        docked_pdbqt,
        ligand_pdb
    )

    with open(output_complex, "w") as out:

        # =========================
        # RECEPTOR
        # =========================

        with open(receptor_pdb, "r") as rec:

            for line in rec:

                if line.startswith(("ATOM", "HETATM")):
                    out.write(line)

        out.write("TER\n")

        # =========================
        # LIGAND
        # =========================

        atom_id = 9000

        with open(ligand_pdb, "r") as lig:

            for line in lig:

                if line.startswith(("ATOM", "HETATM")):

                    newline = (
                        f"HETATM{atom_id:5d} "
                        f"{line[12:16]}"
                        f"LIG A9999    "
                        f"{line[30:]}"
                    )

                    out.write(newline)

                    atom_id += 1

        out.write("END\n")

    os.remove(ligand_pdb)

    print(f"Complex created: {output_complex}")