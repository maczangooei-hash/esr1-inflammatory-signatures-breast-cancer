import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parents[1]
TABLE = BASE / "results/manuscript_tables/Table2_primary_adjusted_associations.tsv"
OUT = BASE / "figures"

with TABLE.open(encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

required = {
    "Platform", "Outcome", "Inflammatory signature",
    "Beta", "95% CI", "N"
}
if not rows or not required.issubset(rows[0]):
    raise ValueError("Table 2 is empty or has unexpected columns.")

signatures = [
    "Inflammatory response",
    "IL6/JAK/STAT3",
    "TNF-alpha/NF-kappa B",
]
outcomes = [
    ("Estrogen response, early", "Early estrogen response"),
    ("Estrogen response, late", "Late estrogen response"),
]
platforms = ["GPL11154", "GPL18573"]

def parse_ci(value):
    low, high = value.split(" to ")
    return float(low), float(high)

fig, axes = plt.subplots(
    1, 2, figsize=(12, 6.8), sharey=True,
    gridspec_kw={"wspace": 0.08}
)

# Two outcome groups, each with three signatures.
positions = {}
y = 0
for outcome, _ in outcomes:
    for sig in signatures:
        positions[(outcome, sig)] = y
        y += 1
    y += 0.8

for ax, platform in zip(axes, platforms):
    for outcome, outcome_label in outcomes:
        for sig in signatures:
            match = [
                r for r in rows
                if r["Platform"] == platform
                and r["Outcome"] == outcome
                and r["Inflammatory signature"] == sig
            ]
            if len(match) != 1:
                raise ValueError(
                    f"Expected exactly one row for {platform}, {outcome}, {sig}; "
                    f"found {len(match)}"
                )

            r = match[0]
            beta = float(r["Beta"])
            low, high = parse_ci(r["95% CI"])
            ypos = positions[(outcome, sig)]

            ax.errorbar(
                beta, ypos,
                xerr=[[beta - low], [high - beta]],
                fmt="o", color="black",
                ecolor="black", elinewidth=1.3,
                capsize=3, markersize=5
            )

    ax.axvline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_title(f"{platform} (n = " +
                 str(next(r["N"] for r in rows if r["Platform"] == platform)) +
                 ")")
    ax.set_xlabel(
        "Adjusted coefficient (95% CI)\n"
        "per 1-SD increase in inflammatory signature"
    )
    ax.grid(axis="x", linestyle=":", alpha=0.35)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

yticks = []
yticklabels = []
for outcome, outcome_label in outcomes:
    for sig in signatures:
        yticks.append(positions[(outcome, sig)])
        yticklabels.append(sig)

axes[0].set_yticks(yticks)
axes[0].set_yticklabels(yticklabels)
axes[1].tick_params(axis="y", left=False, labelleft=False)

# Outcome group labels and separation.
for ax in axes:
    ax.text(
        0.02, 1.02, "Early response (top) · Late response (bottom)",
        transform=ax.transAxes, fontsize=8, va="bottom"
    )
    ax.set_ylim(max(positions.values()) + 0.8, -0.8)

fig.suptitle(
    "Inflammatory signatures and estrogen-response scores",
    fontsize=13, fontweight="bold", y=1.04
)
fig.text(
    0.5, -0.015,
    "Models adjusted for ESR1 expression, age, tumor size, and PAM50 subtype. "
    "Outcomes are mean gene-wise z-scores; coefficients are not standardized "
    "outcome coefficients.",
    ha="center", fontsize=8, wrap=True
)

fig.tight_layout()

pdf = OUT / "Figure1_primary_adjusted_associations_revised.pdf"
tiff = OUT / "Figure1_primary_adjusted_associations_revised.tiff"
fig.savefig(pdf, bbox_inches="tight")
fig.savefig(tiff, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

print("Created:")
print(pdf)
print(tiff)
print(f"Plotted {len(rows)} estimates from Table 2.")
