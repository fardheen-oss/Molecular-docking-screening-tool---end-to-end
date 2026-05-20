import os
import subprocess


# =========================================================
# RUN VINA
# =========================================================

def run_vina(
    receptor,
    ligand,
    out,
    center,
    size,
    exhaustiveness=2
):

    cmd = [

        "vina",

        "--receptor", receptor,
        "--ligand", ligand,

        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),

        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),

        "--exhaustiveness", str(exhaustiveness),

        # SAFE FOR LOW-END LAPTOPS
        "--cpu", "1",

        "--out", out
    ]

    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )


# =========================================================
# EXTRACT VINA SCORE
# =========================================================

def extract_vina_score(pdbqt_file):

    if not os.path.exists(pdbqt_file):
        return None

    with open(pdbqt_file, "r") as f:

        for line in f:

            if "REMARK VINA RESULT:" in line:

                try:
                    return float(line.split()[3])

                except:
                    return None

    return None


# =========================================================
# GET WHOLE PROTEIN GRID
# =========================================================

def get_protein_grid(pdbqt_file):

    coords = []

    with open(pdbqt_file, "r") as f:

        for line in f:

            if line.startswith(("ATOM", "HETATM")):

                try:

                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    coords.append((x, y, z))

                except:
                    pass

    if not coords:
        raise Exception("No coordinates found in receptor")

    mins = [min(c[i] for c in coords) for i in range(3)]
    maxs = [max(c[i] for c in coords) for i in range(3)]

    center = [

        round((mins[i] + maxs[i]) / 2, 3)

        for i in range(3)
    ]

    size = [

        min(126, round((maxs[i] - mins[i]) + 10, 3))

        for i in range(3)
    ]

    return center, size


# =========================================================
# GET POSE CENTER
# =========================================================

def get_pose_center(pdbqt_file):

    coords = []

    if not os.path.exists(pdbqt_file):
        return (0, 0, 0)

    with open(pdbqt_file, "r") as f:

        for line in f:

            if line.startswith(("ATOM", "HETATM")):

                try:

                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    coords.append((x, y, z))

                except:
                    pass

            if "ENDMDL" in line:
                break

    if not coords:
        return (0, 0, 0)

    return (

        round(sum(c[0] for c in coords) / len(coords), 3),

        round(sum(c[1] for c in coords) / len(coords), 3),

        round(sum(c[2] for c in coords) / len(coords), 3)
    )