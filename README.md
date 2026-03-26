# conformer-shape

**Phase 24 — 3D Conformer Generation + Shape Similarity**
Track 1 — Cheminformatics Core

Generates RDKit ETKDGv3 conformers for a compound library, aligns each compound to a
reference via Open3D Alignment (O3A), and ranks compounds by 3D shape similarity —
revealing scaffold hops invisible to 2D fingerprints.

## Why 3D shape similarity finds scaffold hops

2D Tanimoto (ECFP4) measures atom-environment overlap. Two compounds from different
scaffold families can have low 2D similarity yet occupy the same 3D binding pocket
volume — a scaffold hop. Shape Tanimoto captures this: compounds in the top-left of
the `shape_vs_2d.png` scatter (high 3D / low 2D similarity) are scaffold hop candidates
that a 2D-only VS would miss entirely.

## Conformer generation: ETKDGv3

| Step | Detail |
|---|---|
| `Chem.AddHs(mol)` | Explicit Hs required before 3D embedding |
| `AllChem.ETKDGv3()` | Experimental Torsion Knowledge Distance Geometry v3; CSD-trained torsion library |
| `params.randomSeed` | Deterministic conformers (seed=42) |
| `AllChem.MMFFOptimizeMolecule` | MMFF94 strain minimization; preferred for drug-like organics |

## O3A alignment

`AllChem.GetO3A(mobile, ref)` maximizes pharmacophoric feature + shape overlap between
mobile and reference. **Alignment must be called before computing shape distances** —
`ShapeTanimotoDist` operates on current conformer positions.

```python
o3a = AllChem.GetO3A(mobile_mol_h, ref_mol_h)
o3a.Align()   # modifies mobile in-place
score = o3a.Score()
```

## Shape metrics

| Metric | Formula | Range | Interpretation |
|---|---|---|---|
| Shape Tanimoto | `1 - ShapeTanimotoDist(ref, mob)` | [0, 1] | 1 = identical shape |
| Protrude dist | `ShapeProtrudeDist(ref, mob)` | [0, 1] | Fraction of ref volume not covered by mob |

## Reading the scatter plot (`shape_vs_2d.png`)

- **Top-left quadrant** (high 3D / low 2D): shape-similar but fingerprint-dissimilar → scaffold hop candidates
- **Diagonal** (x ≈ y): compounds where 2D and 3D similarity track together (same scaffold family)
- **Bottom-right** (high 2D / low 3D): same 2D features but different 3D shape (e.g., different stereochemistry or conformational preference)

## Quickstart

```bash
python -m venv .venv
.venv\Scripts\pip install -r requirements.txt

PYTHONUTF8=1 .venv\Scripts\python main.py --input data/compounds.csv --reference benz_001_F
```

On Windows PowerShell:

```powershell
$env:PYTHONUTF8=1
.venv\Scripts\python main.py --input data/compounds.csv --reference benz_001_F
```

> **RDKit install note:** `pip install rdkit` works on most systems. If it fails,
> use conda: `conda install -c conda-forge rdkit`

## CLI flags

| Flag | Default | Description |
|---|---|---|
| `--input` | required | Compounds CSV (compound_name, smiles, pic50) |
| `--reference` | `benz_001_F` | compound_name of the reference for alignment |
| `--top-n` | `20` | Top hits for bar chart; top-10 for structure grid |
| `--output-dir` | `output` | Output directory |
| `--seed` | `42` | ETKDG random seed |

## Outputs

| File | Description |
|---|---|
| `output/shape_results.csv` | Per-compound: conformer_ok, o3a_score, shape_tanimoto, protrude_dist, tanimoto_2d |
| `output/shape_bar.png` | Top-N by shape Tanimoto; scaffold-family colored; annotated values |
| `output/shape_vs_2d.png` | 3D shape vs 2D Tanimoto scatter; quadrant annotations |
| `output/top_hits_grid.png` | 2D structure grid of top-10 hits with shape similarity labels |

## Implementation notes

- `matplotlib.use("Agg")` called before `import matplotlib.pyplot` (headless rendering)
- `RDLogger.DisableLog("rdApp.*")` suppresses RDKit warnings
- `MolsToGridImage` returns a PIL Image — save with `img.save(path)`, not `plt.savefig`
- If the reference conformer fails: exits with error (cannot rank without a reference)
- Scaffold family parsed from `compound_name` prefix before first `_`
