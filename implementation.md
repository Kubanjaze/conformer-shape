# Phase 24 — 3D Conformer Generation + Shape Similarity

**Version:** 1.1 (final as-built)
**Author:** Kerwyn Medrano
**Date:** 2026-03-26
**Track:** Track 1 — Cheminformatics Core
**Tier:** Small (3–5 hrs)
**API Cost:** $0.00 — pure RDKit + pandas + matplotlib

---

## 1. Project Overview

### Goal

Generate 3D conformers for a compound library using RDKit's ETKDGv3 force field,
align each compound to a reference conformer using Open3D Alignment (O3A), and
rank compounds by 3D shape similarity. Produce a bar chart of top hits, a scatter
comparing 3D shape vs 2D Tanimoto similarity, and a 2D structure grid of top hits
labeled with shape scores.

```bash
python main.py --input data/compounds.csv --reference benz_001_F
```

Outputs:
- `output/shape_results.csv` — per-compound: conformer status, O3A score, shape Tanimoto, protrude distance
- `output/shape_bar.png` — top-20 compounds ranked by shape Tanimoto, colored by scaffold family
- `output/shape_vs_2d.png` — scatter: 2D Tanimoto (x) vs 3D shape Tanimoto (y); highlight reference + near-misses
- `output/top_hits_grid.png` — 2D structure grid of top-10 hits labeled with shape similarity score

### What This Phase Teaches

| Concept | Detail |
|---|---|
| 3D conformer generation | `AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())` + `AllChem.MMFFOptimizeMolecule` |
| Hydrogen handling | `Chem.AddHs(mol)` before embedding; `Chem.RemoveHs` after for display |
| O3A alignment | `AllChem.GetO3A(mobile, ref)` — optimal shape+pharmacophore alignment; `.Align()` + `.Score()` |
| Shape Tanimoto | `rdShapeHelpers.ShapeTanimotoDist(ref, mob)` after alignment; similarity = 1 − dist |
| Protrude distance | `rdShapeHelpers.ShapeProtrudeDist(ref, mob)` — asymmetric: fraction of ref not covered by mob |
| 3D vs 2D comparison | ECFP4 Tanimoto (2D) vs shape Tanimoto (3D); divergence reveals shape-similar scaffold hops |

### Domain Context

3D shape similarity is used in ligand-based virtual screening to find scaffold hops —
compounds that are 2D-dissimilar but 3D-similar to a known active (same binding pocket
volume, different scaffold). ETKDG (Experimental Torsion Knowledge Distance Geometry)
generates low-energy conformers consistent with the Cambridge Structural Database.
O3A alignment maximizes overlap of pharmacophoric features + shape before computing
ShapeTanimoto. This phase builds on Phase 22 (VS pipeline, Tanimoto ECFP4 scoring)
to introduce the 3D shape dimension of VS.

---

## 2. Architecture

```
conformer-shape/
├── main.py
├── requirements.txt
├── README.md
├── .gitignore
├── data/
│   └── compounds.csv          — 45-compound dataset (reused from Phase 21/23)
└── output/
    ├── shape_results.csv
    ├── shape_bar.png
    ├── shape_vs_2d.png
    └── top_hits_grid.png
```

---

## 3. Input Format

### `data/compounds.csv`
```csv
compound_name,smiles,pic50
benz_001_F,C=CC(=O)Nc1ccc(F)cc1,7.25
...
```
- Reuse the 45-compound Phase 21/23 dataset verbatim
- `--reference`: compound_name of the reference compound (default: benz_001_F)

---

## 4. Pipeline

### Step 1 — Generate conformers (ETKDG)

```python
from rdkit.Chem import AllChem

mol_h = Chem.AddHs(mol)
params = AllChem.ETKDGv3()
params.randomSeed = 42
result = AllChem.EmbedMolecule(mol_h, params)
if result == -1:
    # Embedding failed — try with random coords
    params.useRandomCoords = True
    result = AllChem.EmbedMolecule(mol_h, params)
AllChem.MMFFOptimizeMolecule(mol_h)
```

- `AddHs` required before 3D embedding; explicit Hs improve conformer quality
- `MMFFOptimizeMolecule` minimizes strain energy; returns 0 on success, 1 if not converged
- Returns None (not-None mol) on success; `result == -1` means embedding failed
- Failed compounds: log count; exclude from shape calculations; include in CSV with NaN scores

### Step 2 — O3A alignment

```python
from rdkit.Chem import AllChem

o3a = AllChem.GetO3A(mobile_mol_h, ref_mol_h)
o3a_score = o3a.Score()          # higher = better overlap (pharmacophore + shape)
o3a.Align()                       # modifies mobile_mol_h in-place to align to ref
```

- O3A = Open3D Alignment; maximizes pharmacophoric + shape overlap
- Score is NOT normalized — use shape Tanimoto as the normalized similarity metric
- Must call `.Align()` BEFORE computing shape distances (distances require aligned conformers)

### Step 3 — Shape similarity

```python
from rdkit.Chem import rdShapeHelpers

shape_tanimoto = 1.0 - rdShapeHelpers.ShapeTanimotoDist(ref_mol_h, mobile_mol_h)
protrude_dist  = rdShapeHelpers.ShapeProtrudeDist(ref_mol_h, mobile_mol_h)
```

- `ShapeTanimotoDist` ∈ [0, 1]; 0 = identical shape; convert to similarity = 1 − dist
- `ShapeProtrudeDist` ∈ [0, 1]; fraction of reference volume not covered by mobile; asymmetric
- `shape_tanimoto = 1.0` means identical shape; reference compound scores 1.0 against itself

### Step 4 — 2D Tanimoto (for scatter plot)

```python
from rdkit.Chem import AllChem
from rdkit import DataStructs

ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol_2d, radius=2, nBits=2048)
fp = AllChem.GetMorganFingerprintAsBitVect(mol_2d, radius=2, nBits=2048)
tani_2d = DataStructs.TanimotoSimilarity(ref_fp, fp)
```

- Use original 2D mol (no Hs) for fingerprints — consistent with Phase 21/22

---

## 5. Module Specification

### `load_compounds(path)` → pd.DataFrame
- Standard loader: compound_name, smiles, pic50; skip invalid SMILES
- Return df with mol column (2D, no Hs)

### `generate_conformer(mol_2d, seed)` → Optional[Chem.Mol]
- Add Hs → embed with ETKDGv3 + seed → MMFF optimize → return mol_h (with Hs, 3D)
- Return None on embedding failure; log failure count

### `align_and_score(mobile_mol_h, ref_mol_h)` → dict
- O3A alignment (modifies mobile in-place) → score → ShapeTanimoto → ProtroudeDist
- Return: `{o3a_score, shape_tanimoto, protrude_dist}`

### `compute_2d_tanimoto(mol_2d, ref_mol_2d, radius, nbits)` → float
- ECFP4 Morgan fingerprint Tanimoto to reference (same params as Phase 22)

### `build_results(df, ref_name, seed)` → pd.DataFrame
- For each compound: generate conformer → align → score → compute 2D Tanimoto
- Return results df with columns: compound_name, smiles, pic50, conformer_ok (bool),
  o3a_score, shape_tanimoto, protrude_dist, tanimoto_2d
- Sort by shape_tanimoto descending

### `plot_shape_bar(results_df, top_n, output_path)`
- Horizontal bar chart: top-N by shape_tanimoto, x-axis = shape_tanimoto (0–1)
- Color by scaffold family (parsed from compound_name prefix: benz_, naph_, ind_, etc.)
- Annotate each bar with shape_tanimoto value
- Mark reference compound (score = 1.0 by definition) distinctly
- Save 150 dpi

### `plot_shape_vs_2d(results_df, output_path)`
- Scatter: tanimoto_2d (x) vs shape_tanimoto (y)
- Color by scaffold family
- Annotate quadrants: top-left = "3D similar / 2D dissimilar" (scaffold hop candidates)
- Mark reference compound with a star marker
- Draw x=y dashed line (compounds where 2D ≈ 3D similarity)
- Save 150 dpi

### `plot_top_hits_grid(results_df, df, top_n, output_path)`
- `rdkit.Chem.Draw.MolsToGridImage` of top-N compounds
- Legend per structure: "{compound_name}\nshape={shape_tanimoto:.3f}"
- mol_size=(200, 150), subImgSize=(250, 200)
- Save as PNG via PIL

### `main()`
- `--input` (required): compounds CSV
- `--reference` (default: benz_001_F): compound_name of reference
- `--top-n` (default: 20): top hits for bar chart; top-10 for grid
- `--output-dir` (default: output)
- `--seed` (default: 42): ETKDG random seed
- Print: N compounds, N conformers OK/failed, ranked top-10 table

---

## 6. Seed Data Design

Reuse `data/compounds.csv` verbatim from Phase 21/23 (45 compounds).

Reference compound: `benz_001_F` (`C=CC(=O)Nc1ccc(F)cc1`)

Expected shape similarity ranking intuition:
- Other benz_ variants (same acrylamide-aniline scaffold, different para-substituent) → highest shape Tanimoto
- naph_ variants (larger fused ring, same acrylamide group) → moderate shape Tanimoto
- ind_, quin_ (different heterocycle, no acrylamide) → lower shape Tanimoto
- pyr_, bzim_ (smaller scaffolds) → lowest shape Tanimoto

Expected scatter (shape_tanimoto vs tanimoto_2d):
- benz_ compounds: top-right cluster (high 2D + high 3D)
- naph_ compounds: moderate 2D, moderate-high 3D (similar shape, different ring size)
- ind_/quin_ compounds: low 2D, low-moderate 3D
- Any compound in top-left quadrant (low 2D / high 3D) = scaffold hop candidate

Expected conformer failures: 0 (all compounds are drug-like, ETKDG handles them well)

---

## 7. Verification Checklist

```bash
python main.py --input data/compounds.csv --reference benz_001_F

# Expected:
# - shape_results.csv: 45 rows; reference compound shape_tanimoto = 1.000
# - shape_bar.png: benz_ compounds dominate top-20; naph_ compounds mid-range
# - shape_vs_2d.png: positive correlation 2D vs 3D Tanimoto; scaffold clusters visible
# - top_hits_grid.png: top-10 structures with shape similarity labels
# - Console: reference compound ranks first; top-10 dominated by benz_ variants
```

---

## 8. Risks / Assumptions / Next Step

**Risks:**
- `AllChem.EmbedMolecule` may fail for strained or unusual SMILES → fallback with
  `useRandomCoords=True`; if still fails, exclude from ranking (log count + names)
- `AllChem.GetO3A` requires both mols to have 3D conformers — gate with conformer_ok check
- `rdShapeHelpers.ShapeTanimotoDist` values can exceed 1.0 for very different shapes in
  edge cases — clamp to [0, 1] before computing similarity
- `MolsToGridImage` returns a PIL Image, not a matplotlib figure — save with `.save(path)` not `plt.savefig`
- PIL/Pillow must be installed separately (add to requirements.txt)

**Assumptions:**
- ETKDGv3 with seed=42; MMFF optimization (not UFF) — MMFF is preferred for drug-like organics
- Single conformer per compound (not multi-conformer ensemble); lowest-energy MMFF conformer
- O3A alignment is sufficient for shape comparison; full ROCS-style alignment not needed
- 2D Tanimoto uses ECFP4 radius=2, nBits=2048 — consistent with prior phases
- `compound_name` prefix before `_` (e.g., `benz`, `naph`, `ind`) used as scaffold family label for coloring

**Next step:** Phase 25 — Pharmacophore Feature Mapping. Extract RDKit pharmacophore
features (donors, acceptors, aromatic, hydrophobic, positive/negative ionizable) from
the acrylamide-aniline series, visualize feature distributions across scaffold families,
and highlight shared vs family-specific features.

---

## 9. Actual Results (v1.1)

### Counts

| Metric | Value |
|---|---|
| Compounds loaded (valid) | 45 / 45 |
| Conformers generated (ETKDGv3) | 45 / 45 (0 failures) |
| Reference: benz_001_F | shape_tanimoto = 1.000 (as expected) |

### Top-10 hits (reference: benz_001_F)

| Rank | compound_name | shape_tanimoto | protrude_dist | tanimoto_2d | pic50 |
|---|---|---|---|---|---|
| 1 | benz_001_F (ref) | 1.000 | 0.000 | 1.000 | 7.25 |
| 2 | benz_009_OH | 0.989 | 0.002 | 0.679 | 6.30 |
| 3 | benz_007_Me | 0.961 | 0.003 | 0.655 | 6.60 |
| 4 | benz_011_difluoro | 0.943 | 0.001 | 0.576 | 7.40 |
| 5 | benz_006_NO2 | 0.872 | 0.000 | 0.576 | 7.45 |
| 6 | pyr_006_acetamide | 0.817 | 0.005 | 0.171 | 6.35 |
| 7 | benz_002_Cl | 0.784 | 0.092 | 0.655 | 7.65 |
| 8 | quin_007_CN | 0.784 | 0.107 | 0.064 | 8.10 |
| 9 | quin_005_OMe | 0.781 | 0.118 | 0.064 | 6.90 |
| 10 | benz_008_OMe | 0.772 | 0.081 | 0.594 | 6.35 |

### Key insights

- **All 45 conformers generated successfully** — drug-like compounds with ETKDGv3 (no failures).
- **benz_ family dominates** top shape hits (ranks 1–5, 7, 10) as expected — same scaffold as reference.
- **Scaffold hop candidate confirmed**: `quin_007_CN` ranks 8th by shape (0.784) but has only 0.064 2D Tanimoto — a quinoline CN compound that looks like the reference acrylamide-aniline in 3D despite being 2D-dissimilar. This is exactly the scaffold hop signature that 3D shape VS recovers vs 2D screening.
- **pyr_006_acetamide** at rank 6 (shape=0.817, 2D=0.171) — pyrimidine with acetamide group has high 3D overlap due to the amide carbonyl, another scaffold hop signal.
- `protrude_dist` near 0 for top benz_ compounds confirms they fit tightly inside the reference shape volume.

### Deviations from plan

| Item | Plan | Actual | Note |
|---|---|---|---|
| Conformer failures | 0 expected | 0 actual | All 45 drug-like compounds embedded cleanly |
| naph_ mid-range shape | ~0.6–0.8 expected | naph_ at 0.4–0.7 | Larger fused ring extends beyond reference's smaller envelope |
| quin_007_CN scaffold hop | predicted as candidate | confirmed at rank 8 (shape=0.784, 2D=0.064) | Most striking scaffold hop in dataset | Extract RDKit pharmacophore
features (donors, acceptors, aromatic, hydrophobic, positive/negative ionizable) from
the acrylamide-aniline series, visualize feature distributions across scaffold families,
and highlight shared vs family-specific features.
