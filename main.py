import sys

if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse
import math
import os
import warnings
from typing import Optional

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Draw, rdMolAlign, rdShapeHelpers

RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
# Scaffold family palette
# ---------------------------------------------------------------------------

FAMILY_COLORS = {
    "benz":  "#4C72B0",
    "naph":  "#DD8452",
    "ind":   "#55A868",
    "quin":  "#C44E52",
    "pyr":   "#8172B2",
    "bzim":  "#937860",
    "other": "#808080",
}


def scaffold_family(compound_name: str) -> str:
    """Extract prefix before first '_'; return 'other' if none or unrecognised."""
    prefix = compound_name.split("_")[0] if "_" in compound_name else compound_name
    return prefix if prefix in FAMILY_COLORS else "other"


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def load_compounds(path: str) -> pd.DataFrame:
    """Load CSV; validate columns; coerce pic50; skip invalid SMILES."""
    df = pd.read_csv(path)
    for col in ("compound_name", "smiles", "pic50"):
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col!r}")

    df["pic50"] = pd.to_numeric(df["pic50"], errors="coerce")
    n_nan_pic50 = int(df["pic50"].isna().sum())
    df = df.dropna(subset=["pic50"]).copy()

    records, n_invalid = [], 0
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is None:
            n_invalid += 1
            continue
        records.append({
            "compound_name": str(row["compound_name"]),
            "smiles":        str(row["smiles"]),
            "pic50":         float(row["pic50"]),
            "mol":           mol,
        })

    print(f"  loaded_rows={len(df) + n_nan_pic50}  valid={len(records)}"
          f"  invalid_smiles={n_invalid}  pic50_nan_dropped={n_nan_pic50}")
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Conformer generation
# ---------------------------------------------------------------------------

def generate_conformer(mol_2d: Chem.Mol, seed: int = 42) -> Optional[Chem.Mol]:
    """ETKDGv3 + MMFF; fallback useRandomCoords=True; return mol_h (3D) or None."""
    mol_h = Chem.AddHs(mol_2d)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    result = AllChem.EmbedMolecule(mol_h, params)
    if result == -1:
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol_h, params)
    if result == -1:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol_h)
    except Exception:
        pass  # keep embedded geometry if MMFF fails
    return mol_h


# ---------------------------------------------------------------------------
# 2D Tanimoto
# ---------------------------------------------------------------------------

def compute_2d_tanimoto(mol_2d: Chem.Mol, ref_mol_2d: Chem.Mol,
                        radius: int = 2, nbits: int = 2048) -> float:
    """ECFP4 Morgan Tanimoto similarity to reference (2D mols, no Hs)."""
    ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol_2d, radius, nBits=nbits,
                                                    useChirality=False)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol_2d, radius, nBits=nbits,
                                                useChirality=False)
    return float(DataStructs.TanimotoSimilarity(ref_fp, fp))


# ---------------------------------------------------------------------------
# O3A alignment + shape scoring
# ---------------------------------------------------------------------------

def align_and_score(mobile_mol_h: Chem.Mol, ref_mol_h: Chem.Mol) -> dict:
    """
    O3A align mobile to ref (modifies mobile in-place).
    Returns {o3a_score, shape_tanimoto, protrude_dist}; NaN values on failure.
    """
    NAN = float("nan")
    try:
        o3a = AllChem.GetO3A(mobile_mol_h, ref_mol_h)
        o3a_score = float(o3a.Score())
        o3a.Align()
    except Exception:
        # Fallback: RMSD-based alignment (works on molecules of different size via MCS)
        o3a_score = NAN
        try:
            rdMolAlign.AlignMol(mobile_mol_h, ref_mol_h)
        except Exception:
            pass  # proceed with unaligned conformer; shape scores will be low

    try:
        raw_tani = rdShapeHelpers.ShapeTanimotoDist(ref_mol_h, mobile_mol_h)
        raw_prot = rdShapeHelpers.ShapeProtrudeDist(ref_mol_h, mobile_mol_h)
        shape_tanimoto = max(0.0, min(1.0, 1.0 - float(raw_tani)))
        protrude_dist  = max(0.0, min(1.0, float(raw_prot)))
    except Exception:
        shape_tanimoto = NAN
        protrude_dist  = NAN

    return {"o3a_score": o3a_score, "shape_tanimoto": shape_tanimoto,
            "protrude_dist": protrude_dist}


# ---------------------------------------------------------------------------
# Build results table
# ---------------------------------------------------------------------------

def build_results(df: pd.DataFrame, ref_name: str, seed: int) -> pd.DataFrame:
    """
    For each compound: generate conformer -> O3A align -> shape scores -> 2D tanimoto.
    Returns results DataFrame sorted by shape_tanimoto desc.
    """
    # Locate reference
    ref_rows = df[df["compound_name"] == ref_name]
    if ref_rows.empty:
        raise ValueError(f"Reference compound {ref_name!r} not found in input CSV.")
    ref_mol_2d = ref_rows.iloc[0]["mol"]
    ref_conf = generate_conformer(ref_mol_2d, seed=seed)
    if ref_conf is None:
        raise RuntimeError(f"Could not generate conformer for reference {ref_name!r}. Exiting.")

    # Generate all conformers
    conf_map: dict[str, Optional[Chem.Mol]] = {}
    failed = []
    for _, row in df.iterrows():
        name = row["compound_name"]
        mol_3d = generate_conformer(row["mol"], seed=seed)
        conf_map[name] = mol_3d
        if mol_3d is None:
            failed.append(name)

    n_ok   = sum(1 for v in conf_map.values() if v is not None)
    n_fail = len(failed)
    print(f"  conformer_ok={n_ok}  conformer_fail={n_fail}")
    if failed:
        preview = ", ".join(failed[:10])
        suffix  = f" ... (+{len(failed)-10} more)" if len(failed) > 10 else ""
        print(f"  [warn] Failed: {preview}{suffix}", file=sys.stderr)

    records = []
    for _, row in df.iterrows():
        name   = row["compound_name"]
        mol_2d = row["mol"]
        mol_3d = conf_map.get(name)

        tani_2d = compute_2d_tanimoto(mol_2d, ref_mol_2d)

        if name == ref_name:
            # Reference scores 1.0 against itself by definition
            records.append({
                "compound_name":  name,
                "smiles":         row["smiles"],
                "pic50":          row["pic50"],
                "conformer_ok":   mol_3d is not None,
                "o3a_score":      0.0,   # align-to-self score is trivially 0
                "shape_tanimoto": 1.0,
                "protrude_dist":  0.0,
                "tanimoto_2d":    tani_2d,
            })
            continue

        if mol_3d is None:
            records.append({
                "compound_name":  name,
                "smiles":         row["smiles"],
                "pic50":          row["pic50"],
                "conformer_ok":   False,
                "o3a_score":      float("nan"),
                "shape_tanimoto": float("nan"),
                "protrude_dist":  float("nan"),
                "tanimoto_2d":    tani_2d,
            })
            continue

        scores = align_and_score(mol_3d, ref_conf)
        records.append({
            "compound_name":  name,
            "smiles":         row["smiles"],
            "pic50":          row["pic50"],
            "conformer_ok":   True,
            "o3a_score":      scores["o3a_score"],
            "shape_tanimoto": scores["shape_tanimoto"],
            "protrude_dist":  scores["protrude_dist"],
            "tanimoto_2d":    tani_2d,
        })

    result_df = pd.DataFrame(records)
    result_df = result_df.sort_values(
        by=["shape_tanimoto", "tanimoto_2d", "compound_name"],
        ascending=[False, False, True],
        na_position="last",
    ).reset_index(drop=True)
    return result_df


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_shape_bar(results_df: pd.DataFrame, ref_name: str, top_n: int,
                   output_path: str) -> None:
    """Horizontal bar chart: top-N by shape_tanimoto; family colors; annotated values."""
    ok_df = results_df[results_df["conformer_ok"]].dropna(subset=["shape_tanimoto"])
    if ok_df.empty:
        print("[warn] No valid shape scores — skipping shape_bar.", file=sys.stderr)
        return

    plot_df = ok_df.head(top_n).iloc[::-1]  # reverse for bottom-up display

    families = [scaffold_family(n) for n in plot_df["compound_name"]]
    colors   = [FAMILY_COLORS.get(f, FAMILY_COLORS["other"]) for f in families]
    labels   = [
        f"{n} *" if n == ref_name else n
        for n in plot_df["compound_name"]
    ]

    fig, ax = plt.subplots(figsize=(10, max(4, len(plot_df) * 0.35)))
    bars = ax.barh(range(len(plot_df)), plot_df["shape_tanimoto"].values,
                   color=colors, edgecolor="white", height=0.7)

    for i, (bar, val) in enumerate(zip(bars, plot_df["shape_tanimoto"].values)):
        ax.text(val + 0.01, bar.get_y() + bar.get_height() / 2,
                f"{val:.2f}", va="center", ha="left", fontsize=8)

    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("3D Shape Tanimoto Similarity", fontsize=11)
    ax.set_title(f"Shape Similarity to {ref_name} (top {len(plot_df)})",
                 fontsize=13, fontweight="bold")

    # Legend: only families present
    seen_families = dict.fromkeys(families)
    legend_handles = [
        mpatches.Patch(facecolor=FAMILY_COLORS.get(f, FAMILY_COLORS["other"]),
                       edgecolor="white", label=f)
        for f in seen_families
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=9, framealpha=0.85)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_shape_vs_2d(results_df: pd.DataFrame, ref_name: str,
                     output_path: str) -> None:
    """Scatter: 2D Tanimoto (x) vs 3D Shape Tanimoto (y); family colors; quadrant labels."""
    ok_df = results_df[results_df["conformer_ok"]].dropna(subset=["shape_tanimoto"])
    if ok_df.empty:
        print("[warn] No valid shape scores — skipping shape_vs_2d.", file=sys.stderr)
        return

    fig, ax = plt.subplots(figsize=(8, 7))

    # Plot by family
    seen_families = {}
    for _, row in ok_df.iterrows():
        fam   = scaffold_family(row["compound_name"])
        color = FAMILY_COLORS.get(fam, FAMILY_COLORS["other"])
        is_ref = row["compound_name"] == ref_name
        if is_ref:
            ax.scatter(row["tanimoto_2d"], row["shape_tanimoto"],
                       s=200, color="black", marker="*", zorder=5,
                       edgecolors="black", linewidths=0.5)
        else:
            ax.scatter(row["tanimoto_2d"], row["shape_tanimoto"],
                       s=60, color=color, alpha=0.8, edgecolors="white",
                       linewidths=0.3, zorder=3)
        if fam not in seen_families:
            seen_families[fam] = color

    # Diagonal x=y
    lims = [0, 1.05]
    ax.plot(lims, lims, "--", color="lightgray", linewidth=1, zorder=1)

    # Quadrant annotations
    ax.text(0.05, 0.92, "3D similar\n2D dissimilar\n(scaffold hop)",
            transform=ax.transAxes, fontsize=8, color="steelblue",
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="aliceblue",
                      edgecolor="steelblue", alpha=0.7))
    ax.text(0.95, 0.08, "2D similar\n3D dissimilar",
            transform=ax.transAxes, fontsize=8, color="firebrick",
            va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="mistyrose",
                      edgecolor="firebrick", alpha=0.7))

    # Legend
    legend_handles = [
        mpatches.Patch(facecolor=c, edgecolor="white", label=f)
        for f, c in seen_families.items()
    ]
    legend_handles.append(
        plt.Line2D([0], [0], marker="*", color="w", markerfacecolor="black",
                   markersize=12, label=f"{ref_name} (ref)")
    )
    ax.legend(handles=legend_handles, fontsize=9, framealpha=0.85,
              loc="upper left", bbox_to_anchor=(0.0, 0.75))

    ax.set_xlim(-0.05, 1.1)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("2D Tanimoto (ECFP4)", fontsize=11)
    ax.set_ylabel("3D Shape Tanimoto", fontsize=11)
    ax.set_title("3D Shape vs 2D Fingerprint Similarity", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_top_hits_grid(results_df: pd.DataFrame, top_n: int,
                       output_path: str) -> None:
    """2D structure grid of top-min(10,top_n) hits; labeled with shape similarity."""
    grid_n = min(10, top_n)

    # Choose ranking column
    has_shape = results_df["conformer_ok"].any() and \
                not results_df["shape_tanimoto"].isna().all()

    if has_shape:
        top_rows = (results_df
                    .dropna(subset=["shape_tanimoto"])
                    [results_df["conformer_ok"]]
                    .head(grid_n))
        label_key = "shape_tanimoto"
        label_fmt = lambda r: f"{r['compound_name']}\nshape={r[label_key]:.3f}"
    else:
        top_rows = (results_df
                    .sort_values("tanimoto_2d", ascending=False)
                    .head(grid_n))
        label_key = "tanimoto_2d"
        label_fmt = lambda r: f"{r['compound_name']}\n2D={r[label_key]:.3f}"

    mols    = [Chem.MolFromSmiles(smi) for smi in top_rows["smiles"]]
    legends = [label_fmt(r) for _, r in top_rows.iterrows()]

    # Filter out any None mols
    valid = [(m, l) for m, l in zip(mols, legends) if m is not None]
    if not valid:
        print("[warn] No valid mols for grid — skipping top_hits_grid.", file=sys.stderr)
        return
    mols, legends = zip(*valid)

    img = Draw.MolsToGridImage(
        list(mols), molsPerRow=5, subImgSize=(250, 200),
        legends=list(legends), useSVG=False,
    )
    img.save(output_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="3D conformer generation + shape similarity ranking.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input",     required=True,           help="Compounds CSV")
    parser.add_argument("--reference", default="benz_001_F",    help="Reference compound_name")
    parser.add_argument("--top-n",     type=int, default=20,    help="Top hits for bar chart")
    parser.add_argument("--output-dir",default="output",        help="Output directory")
    parser.add_argument("--seed",      type=int, default=42,    help="ETKDG random seed")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load ---
    print(f"\nLoading: {args.input}")
    df = load_compounds(args.input)

    # --- Build results ---
    print(f"\nGenerating conformers + aligning to reference: {args.reference!r}")
    try:
        results_df = build_results(df, args.reference, args.seed)
    except (ValueError, RuntimeError) as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)

    # --- Save CSV ---
    csv_path = os.path.join(args.output_dir, "shape_results.csv")
    out_cols  = ["compound_name", "smiles", "pic50", "conformer_ok",
                 "o3a_score", "shape_tanimoto", "protrude_dist", "tanimoto_2d"]
    results_df[out_cols].to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # Decide if we have usable 3D scores
    has_3d = results_df["conformer_ok"].any() and \
             not results_df["shape_tanimoto"].isna().all()

    # --- Plots ---
    if has_3d:
        bar_path = os.path.join(args.output_dir, "shape_bar.png")
        plot_shape_bar(results_df, args.reference, args.top_n, bar_path)
        print(f"Saved: {bar_path}")

        scatter_path = os.path.join(args.output_dir, "shape_vs_2d.png")
        plot_shape_vs_2d(results_df, args.reference, scatter_path)
        print(f"Saved: {scatter_path}")
    else:
        print("[warn] No valid 3D scores — skipping shape_bar and shape_vs_2d.",
              file=sys.stderr)

    grid_path = os.path.join(args.output_dir, "top_hits_grid.png")
    plot_top_hits_grid(results_df, args.top_n, grid_path)
    print(f"Saved: {grid_path}")

    # --- Console: top-10 table ---
    print(f"\n--- Top-10 by shape similarity ---")
    display_cols = ["compound_name", "shape_tanimoto", "protrude_dist",
                    "tanimoto_2d", "pic50"]
    top10 = results_df.dropna(subset=["shape_tanimoto"]).head(10)
    if not top10.empty:
        print(top10[display_cols].to_string(index=False))
    else:
        print("(no valid 3D scores)")

    print("\nDone.")


if __name__ == "__main__":
    main()
