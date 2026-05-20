# NatDock — Natural Compound Molecular Docking Pipeline

A fully automated, production-grade molecular docking pipeline for screening natural phytochemical compounds against protein targets. Built with open-source tools and machine learning-based ADMET prediction.

---

## What This Pipeline Does

NatDock takes a protein structure and a library of compounds, then automatically:

1. **Prepares the protein** — cleans PDB, removes heteroatoms, converts to PDBQT format
2. **Screens compounds for drug-likeness (ADMET)** — filters using Lipinski Ro5, Veber, Egan, Muegge rules and ML-based toxicity prediction (AMES mutagenicity, hERG cardiotoxicity, DILI hepatotoxicity, CYP inhibition, solubility, bioavailability)
3. **Performs blind docking** — docks all ADMET-passing compounds against the whole protein using AutoDock Vina
4. **Selects Top 3 hits** — ranked by binding affinity (kcal/mol)
5. **Predicts binding pockets** — uses P2Rank machine learning to identify the most druggable cavity on the protein surface, independent of any ligand
6. **Performs active-site refinement** — re-docks Top 3 compounds into the predicted binding pocket with higher exhaustiveness for accurate pose prediction
7. **Generates Discovery Studio-compatible complexes** — produces receptor-ligand complex PDB files with correct HETATM formatting for ligand interaction analysis, hydrogen bond detection, and 2D interaction diagrams

---

## Output Structure
results/
├── clean_protein.pdb
├── receptor.pdbqt
├── top_hits.csv
├── admet/
│   ├── admet_results.csv
│   └── prepared_ligands/
│       └── {compound}.pdbqt
├── blind_docking/
│   ├── {compound}_blind.pdbqt
│   ├── {compound}_blind.pdb
│   ├── {compound}_complex.pdb
│   └── blind_docking_exhaustiveness_N.csv
├── p2rank/
│   └── predicted_sites.csv
└── active_refinement/
├── {compound}_active.pdbqt
├── {compound}_active.pdb
├── {compound}_complex.pdb
└── active_refinement_results.csv
---

## Requirements

### Python packages
Install with:
```bash
pip install rdkit pandas numpy admet-ai biopython
```

### External tools (must be installed separately)

| Tool | Version | Download |
|---|---|---|
| AutoDock Vina | 1.2+ | https://vina.scripps.edu/downloads/ |
| Open Babel | 3.1+ | https://openbabel.org/wiki/Get_Open_Babel |
| P2Rank | 2.5.1 | https://github.com/rdk/p2rank/releases |
| Java (JDK) | 17+ | https://adoptium.net |

### Verify installations
```bash
vina --version
obabel --version
java -version
```

---

## Input Format

### Protein
Any standard `.pdb` file. Example: protein.pdb
### Ligand CSV
A CSV file with exactly these two columns:
```csv
Ligand,SMILES
quercetin,C1=CC(=C(C=C1)O)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O
apigenin,C1=CC(=CC=C1)C2=CC(=O)C3=C(O2)C=C(C=C3O)O
podophyllotoxin,CC1=CC(=CC(=C1OC)OC)OC2C3C(COC3=O)C(C4=CC5=C(C=C42)OCO5)O
```

---

## Usage

### Basic run
```bash
python main_pipeline.py --protein protein.pdb --ligands ligands.csv
```

### With P2Rank (recommended for accurate active site)
```bash
python main_pipeline.py --protein protein.pdb --ligands ligands.csv --p2rank E:\script\p2rank
```

### P2Rank environment setup (Windows)
If P2Rank and Java are not globally installed, set these at the top of `site_predictor.py`:
```python
os.environ["JAVA_HOME"] = r"E:\script\java"
os.environ["P2RANK_HOME"] = r"E:\script\p2rank"
os.environ["PATH"] += os.pathsep + os.path.join(r"E:\script\java", "bin")
```

---

## ADMET Filtering Criteria

| Property | Threshold | Source |
|---|---|---|
| Molecular Weight | ≤ 500 Da | Lipinski Ro5 |
| logP | ≤ 5 | Lipinski Ro5 |
| TPSA | ≤ 140 Å² | Veber / Egan |
| HBD | ≤ 5 | Lipinski Ro5 |
| HBA | ≤ 10 | Lipinski Ro5 |
| Rotatable Bonds | ≤ 10 | Veber |
| PAINS | None | RDKit FilterCatalog |
| AMES Mutagenicity | < 0.5 | ADMET-AI ML |
| hERG Cardiotoxicity | < 0.5 | ADMET-AI ML |
| DILI Hepatotoxicity | < 0.7 | ADMET-AI ML |
| Solubility (logS) | > −6 | ADMET-AI ML |

---

## Pipeline Modules

| File | Purpose |
|---|---|
| `main_pipeline.py` | Master orchestration script |
| `protein_prep.py` | PDB cleaning and PDBQT conversion |
| `admet_ht.py` | ADMET screening and ligand preparation |
| `blind_docking.py` | Whole-protein blind docking |
| `docking.py` | AutoDock Vina wrapper functions |
| `site_predictor.py` | P2Rank binding site prediction |
| `active_site_refine.py` | Active-site refinement docking |

---

## Example Results

| Rank | Ligand | Blind Affinity | Refined Affinity |
|---|---|---|---|
| 1 | podophyllotoxin | -8.190 kcal/mol | -8.565 kcal/mol |
| 2 | quercetin | -7.529 kcal/mol | -7.552 kcal/mol |
| 3 | apigenin | -7.496 kcal/mol | -7.477 kcal/mol |

---

## Visualization

Complex PDB files are compatible with:
- **BIOVIA Discovery Studio** — ligand interaction diagrams, hydrogen bond analysis, 2D interaction maps
- **PyMOL** — structure visualization
- **ChimeraX** — P2Rank generates `.cxc` scripts for direct pocket visualization

---

## Citation

If you use this pipeline in your research, please cite:

- **AutoDock Vina**: Trott O, Olson AJ. *J Comput Chem.* 2010
- **P2Rank**: Krivak R, Hoksza D. *J Cheminform.* 2018
- **ADMET-AI**: Swanson K et al. *Bioinformatics.* 2024
- **RDKit**: https://www.rdkit.org

---

## License

MIT License — free to use, modify, and distribute with attribution.

---

## Author

Developed as part of a computational drug discovery study targeting natural phytochemical compounds.