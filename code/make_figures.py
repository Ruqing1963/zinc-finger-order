# -*- coding: utf-8 -*-
"""Figure 2: finger-order ladders for ZNF91 and ZNF721, from results/case_*.json.  Needs matplotlib.
Writes figures/fig_ladder.pdf."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import spectral_fibres as SF

FIG = os.path.join(SF.ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
INK, INK2, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 8.5, "pdf.fonttype": 42,
                     "axes.titlesize": 9, "axes.titleweight": "bold", "axes.titlelocation": "left",
                     "axes.titlecolor": INK})


def match(fS, fT):
    """for each finger of S (array position i), the array position of the same finger (same helix) in T;
    identical helices are matched in order"""
    out = [i if fT[i][2] == a[2] else None for i, a in enumerate(fS)]   # fingers that stay in place
    used = {i for i in out if i is not None}
    for i, a in enumerate(fS):
        if out[i] is None:
            cands = [j for j, b in enumerate(fT) if b[2] == a[2] and j not in used]
            j = min(cands, key=lambda j: abs(j - i))
            used.add(j)
            out[i] = j
    return out


genes = [g for g in ("ZNF91", "ZNF721") if os.path.exists(os.path.join(SF.RESULTS, f"case_{g}.json"))]
fig, axes = plt.subplots(len(genes), 1, figsize=(7.2, 2.25 * len(genes)))
for ax, g in zip(axes, genes):
    C = json.load(open(os.path.join(SF.RESULTS, f"case_{g}.json")))
    fS, fT = C["fS"], C["fT"]
    n = len(fS)
    m = match(fS, fT)
    yS, yT, w, h = 1.0, 0.0, 0.8, 0.34
    for i, a in enumerate(fS):
        moved = fT[i][3] != a[3]
        for y, lab in ((yS, a[3]), (yT, fT[i][3])):
            ax.add_patch(FancyBboxPatch((i - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.0,rounding_size=0.08",
                                        fc=("#fbe3d8" if moved else "#eeeeea"), ec=(ORANGE if moved else BASE), lw=0.7))
            ax.text(i, y, lab, ha="center", va="center", fontsize=4.6, family="monospace", color=INK, rotation=90)
    for i, j in enumerate(m):
        if i != j:
            ax.plot([i, j], [yS - h / 2, yT + h / 2], color=ORANGE, lw=0.7, alpha=0.8)
    ax.text(-1.0, yS, "wild type", ha="right", va="center", fontsize=7.5, color=INK2)
    ax.text(-1.0, yT, "fibre-mate", ha="right", va="center", fontsize=7.5, color=INK2)
    moved = sum(fT[i][3] != fS[i][3] for i in range(n))
    ax.set_title(f"{g}: {n} fingers, {moved} array positions with different recognition residues "
                 f"(same 25-spectrum)")
    ax.set_xlim(-4.5, n)
    ax.set_ylim(-0.45, 1.45)
    ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig_ladder.pdf"), bbox_inches="tight")
print("figure written")
