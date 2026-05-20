import os
os.environ["JAVA_HOME"] = r"E:\script\java"
os.environ["P2RANK_HOME"] = r"E:\script\p2rank"
os.environ["PATH"] += os.pathsep + os.path.join(r"E:\script\java", "bin")
import subprocess
import platform
import pandas as pd


def predict_binding_sites(
    receptor_pdb,
    output_folder,
    p2rank_dir=None,
    top_n=3,
):
    if p2rank_dir is None:
        p2rank_dir = os.environ.get("P2RANK_HOME", "")

    if not p2rank_dir or not os.path.isdir(p2rank_dir):
        print(
            "\n[P2Rank] P2RANK_HOME not set or directory not found.\n"
            "  Download P2Rank from:\n"
            "    https://github.com/rdk/p2rank/releases\n"
            "  Then either:\n"
            "    set P2RANK_HOME=C:\\tools\\p2rank   (Windows)\n"
            "    export P2RANK_HOME=~/tools/p2rank   (Linux/Mac)\n"
            "  OR pass p2rank_dir= to predict_binding_sites().\n"
        )
        return []

    if platform.system() == "Windows":
        p2rank_bin = os.path.join(p2rank_dir, "prank.bat")
    else:
        p2rank_bin = os.path.join(p2rank_dir, "prank")

    if not os.path.exists(p2rank_bin):
        print(f"[P2Rank] Executable not found: {p2rank_bin}")
        return []

    os.makedirs(output_folder, exist_ok=True)

    print(f"\n[P2Rank] Predicting binding sites on: {receptor_pdb}")

    cmd = [
        p2rank_bin,
        "predict",
        "-f", os.path.abspath(receptor_pdb),
        "-o", os.path.abspath(output_folder),
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=p2rank_dir,
        shell=True,
    )

    if result.returncode != 0:
        print(f"[P2Rank] Error:\n{result.stderr.decode().strip()}")
        return []

   pdb_stem     = os.path.splitext(os.path.basename(receptor_pdb))[0]
pdb_basename = os.path.basename(receptor_pdb)

pred_csv = os.path.join(
    os.path.abspath(output_folder),
    f"{pdb_basename}_predictions.csv",
)

if not os.path.exists(pred_csv):
    pred_csv = os.path.join(
        os.path.abspath(output_folder),
        f"predict_{pdb_stem}",
        f"{pdb_stem}_predictions.csv",
    )

    if not os.path.exists(pred_csv):
        print(f"[P2Rank] Prediction CSV not found: {pred_csv}")
        print(f"  P2Rank stdout: {result.stdout.decode().strip()[:300]}")
        return []

    df = pd.read_csv(pred_csv)
    df.columns = [c.strip() for c in df.columns]

    required = {"name", "score", "center_x", "center_y", "center_z"}
    if not required.issubset(set(df.columns)):
        df = df.rename(columns=lambda c: c.strip())

    pockets = []

    for i, row in df.head(top_n).iterrows():
        pocket = {
            "pocket":   int(i + 1),
            "score":    round(float(row.get("score", 0)), 3),
            "center_x": round(float(row.get("center_x", 0)), 3),
            "center_y": round(float(row.get("center_y", 0)), 3),
            "center_z": round(float(row.get("center_z", 0)), 3),
            "residues": str(row.get("residue_ids", "")).strip(),
        }
        pockets.append(pocket)
        print(
            f"  Pocket {pocket['pocket']:2d} | "
            f"Score {pocket['score']:6.2f} | "
            f"Centre ({pocket['center_x']:.2f}, "
            f"{pocket['center_y']:.2f}, "
            f"{pocket['center_z']:.2f})"
        )

    print(f"[P2Rank] {len(pockets)} pocket(s) identified.")
    return pockets


def pocket_grid_size(pocket_score, base=20, max_size=30):
    fraction = min(pocket_score / 60.0, 1.0)
    size     = round(base + fraction * (max_size - base))
    return (size, size, size)


def save_site_summary(pockets, output_folder):
    rows = []
    for p in pockets:
        sx, sy, sz = pocket_grid_size(p["score"])
        rows.append({
            "Pocket":    p["pocket"],
            "Score":     p["score"],
            "Center_X":  p["center_x"],
            "Center_Y":  p["center_y"],
            "Center_Z":  p["center_z"],
            "Grid_Size": sx,
            "Residues":  p["residues"],
        })

    csv_path = os.path.join(output_folder, "predicted_sites.csv")
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    print(f"[P2Rank] Site summary → {csv_path}")
    return csv_path