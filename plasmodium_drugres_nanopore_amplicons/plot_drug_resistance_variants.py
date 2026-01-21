# %%
"""
Create stacked bar chart of genotype counts by gene -> protein_change,
with side-by-side bars for 2024 and 2025, from a CSV file.

- Excludes synonymous variants (type == 'synonymous_variant')
- Colors: reference(0)=blue, alternate(1)=orange, NA=grey
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

# %% set working directory
os.chdir("/mnt/storage11/sophie/plasmodium_cape_verde/CV_PF_DR_amplicons3/pf_drug_resistance/pf_drugres_nanopore_fastq_files/downsampled_bams_and_output")
print(os.getcwd())

CSV_PATH = "Pf_CV_drug_resistance_collated_results.variants_MPCorrected.csv"

# Load
df = pd.read_csv(CSV_PATH)

# %% Basic column checks (raises clear errors if missing)
required_cols = {"gene_name", "protein_change", "year",
    "genotype", "type"}

# %% 1) Keep only the missense variants 
df = df[df["type"] == "missense_variant"].copy()

# %% 2) Make genotype clean strings BEFORE grouping/pivoting
df["genotype"] = df["genotype"].map({1.0: "1", 0.0: "0"}).fillna("NA")

# %% 3) Group and pivot to counts
grp = (
    df.groupby(
        ["gene_name", "protein_change", "year", "genotype"]
    )
    .size()
    .reset_index(name="count")
)

# %%
pivot = (
    grp.pivot_table(
        index=["gene_name", "protein_change", "year"],
        columns="genotype",
        values="count",
        fill_value=0,
        aggfunc="sum",
    )
    .reset_index()
)

df["genotype"].unique()

# %% 
for col in ["0", "1", "NA"]:
    if col not in pivot.columns:
        pivot[col] = 0

# %% 4) Plot: one subplot per gene

color_map = {"0": "purple", "1": "darkorange", "NA": "rosybrown"}
years_order = [2024, 2025]

genes_order = list(pivot["gene_name"].drop_duplicates())
n_genes = len(genes_order)

# Create and save one figure per gene (separate plots)
for gene in genes_order:
    gsub = pivot[pivot["gene_name"] == gene]

    # Determine protein_changes for this gene, excluding those that are all NA
    all_proteins = list(gsub["protein_change"].drop_duplicates())
    proteins_keep = []
    for prot in all_proteins:
        psub = gsub[gsub["protein_change"] == prot]
        if psub[["0", "1"]].to_numpy().sum() > 0:
            proteins_keep.append(prot)

    # If nothing to plot for this gene, skip
    if len(proteins_keep) == 0:
        continue

    # --- Auto width per gene: scale with number of bars (2 per protein_change) ---
    n_ticks = len(proteins_keep) * 2
    # Base width + per-tick expansion (tweak factors if you want a little more/less space)
    fig_width = max(12, min(36, 6 + 0.38 * n_ticks))
    # -----------------------------------------------------------------------------

    fig, ax = plt.subplots(1, 1, figsize=(fig_width, 4), dpi=600)

    x_labels = []
    x_positions = []
    pos = 0.0
    bar_width = 0.35
    protein_gap = 0.8
    year_gap = 0.15   # small gap between 2024 & 2025 bars

    for prot in proteins_keep:
        psub = gsub[gsub["protein_change"] == prot]

        for i, yr in enumerate(years_order):
            row = psub[psub["year"] == yr]
            if not row.empty:
                counts = [int(row["0"].iloc[0]), int(row["1"].iloc[0]), int(row["NA"].iloc[0])]
            else:
                counts = [0, 0, 0]

            # x position of this year's bar
            bar_x = pos + i * (bar_width + year_gap)

            bottom = 0
            for label, ct in zip(["0", "1", "NA"], counts):
                bar = ax.bar(bar_x, ct, width=bar_width, bottom=bottom)
                bar[0].set_facecolor(color_map[label])

                if ct > 0:
                    ax.text(
                        bar_x,
                        bottom + ct/2,
                        str(ct),
                        ha="center", va="center",
                        fontsize=8, color="black", fontweight="bold"
                    )
                bottom += ct

            x_labels.append(f"{prot}\n{yr}")
            x_positions.append(bar_x)

        # after both year bars, move to next protein_change
        pos += (2 * bar_width + year_gap) + protein_gap

    # Formatting
    ax.set_title(f"Gene: {gene}", fontweight="bold")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_ylabel("Count")
    ax.tick_params(axis="x", labelsize=8)

    # Legend
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=color_map["0"], label="Reference (0)"),
        Patch(facecolor=color_map["1"], label="Alternate (1)"),
        Patch(facecolor=color_map["NA"], label="Missing (NA)"),
    ]
    ax.legend(handles=legend_handles, title="Genotype", loc="center left", bbox_to_anchor=(1.02, 0.5))

    plt.tight_layout()

    # Show interactively and save
    plt.show()
    safe_gene = gene.replace("/", "_").replace(" ", "_")
    plt.savefig(f"genotype_stacked_bars_{safe_gene}.png", dpi=600)
# %%
