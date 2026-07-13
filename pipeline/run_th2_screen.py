#!/usr/bin/env python
"""
run_th2_screen.py
=================
Rerunnable pipeline: genome-scale CD4+ T cell Perturb-seq phenocopy screen for
intracellular Th2-collapse master regulators.

INPUT  : GWCD4i.DE_stats.h5ad  (genome-wide DESeq2 differential-expression object,
         Marson lab 2025 CD4+ T cell Perturb-seq atlas; 33,983 perturbation x
         condition rows, 10,282 measured genes; layers include per-gene z-scores)
OUTPUT : ranked shortlist CSV + per-axis annotation tables + figures.

Dependencies: h5py, numpy, pandas, scipy, matplotlib
External DB annotation (degradability, GWAS) is a separate step using the
Open Targets Platform GraphQL API and the GWAS Catalog REST API (see methods).

Usage:
    python run_th2_screen.py --h5ad GWCD4i.DE_stats.h5ad --out ./out
"""
import argparse, os
import numpy as np, pandas as pd, h5py

STIM = ("Stim8hr", "Stim48hr")


# ---------- load ----------
def load_de(h5ad):
    f = h5py.File(h5ad, "r")
    dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
    gene_name = dec(f["var"]["gene_name"][:])
    gene_ids  = dec(f["var"]["gene_ids"][:])

    def col(g):
        if isinstance(g, h5py.Group) and "categories" in g and "codes" in g:
            return pd.Categorical.from_codes(g["codes"][:], categories=dec(g["categories"][:]))
        a = g[:]
        return dec(a) if (a.dtype.kind == "S" or a.dtype == object) else a

    obs = pd.DataFrame({k: col(f["obs"][k]) for k in f["obs"].keys() if k != "index"})
    obs.index = dec(f["obs"]["index"][:])

    Z = f["layers"]["zscore"]; n, g = Z.shape
    zscore = np.empty((n, g), np.float32)
    for i in range(0, n, 2000):
        zscore[i:i+2000] = Z[i:i+2000].astype(np.float32)
    # adj_p loaded lazily per-anchor to build DE subspaces
    return f, obs, gene_name, gene_ids, zscore


# ---------- eligibility ----------
def eligible_mask(obs):
    return (obs["ontarget_significant"].values
            & ~obs["neighboring_gene_KD"].values
            & ~obs["distal_offtarget_flag"].values
            & ~obs["low_target_gex"].values)


def ridx(obs, gene, cond):
    m = np.where((obs.target_contrast_gene_name.values == gene)
                 & (obs.culture_condition.values == cond))[0]
    return int(m[0]) if len(m) else None


# ---------- axis ----------
def build_axis(f, obs, zscore, anchor, fdr=0.10):
    """Consensus stim z-vector + informative (anchor-DE) subspace mask."""
    rows = {c: ridx(obs, anchor, c) for c in STIM}
    rows = {c: i for c, i in rows.items() if i is not None}
    axis_z = np.vstack([zscore[i] for i in rows.values()]).mean(0)
    sig = np.zeros(zscore.shape[1], bool)
    for c, i in rows.items():
        sig |= (f["layers"]["adj_p_value"][i].astype(np.float32) < fdr)
    return axis_z, sig


def phenocopy(query_z, axis_z, sig):
    q = np.atleast_2d(query_z)[:, sig].astype(np.float64)
    a = axis_z[sig].astype(np.float64)
    cos = (q / (np.linalg.norm(q, axis=1, keepdims=True) + 1e-9)) @ (a / (np.linalg.norm(a) + 1e-9))
    return cos


# ---------- screen ----------
def run(h5ad, out):
    os.makedirs(out, exist_ok=True)
    f, obs, gname, gids, Z = load_de(h5ad)
    obs["eligible"] = eligible_mask(obs)

    axes = {a: build_axis(f, obs, Z, a) for a in ("GATA3", "STAT6")}

    def screen(anchor, cond):
        az, sig = axes[anchor]
        sel = np.where((obs.culture_condition.values == cond) & obs.eligible.values)[0]
        cos = phenocopy(Z[sel], az, sig)
        return pd.DataFrame({"gene": obs.target_contrast_gene_name.values[sel],
                             "cosine": cos,
                             "n_downstream": obs.n_downstream.values[sel]})

    results = {}
    for anchor in ("GATA3", "STAT6"):
        d8, d48 = screen(anchor, "Stim8hr"), screen(anchor, "Stim48hr")
        m = d8.merge(d48, on="gene", suffixes=("_8", "_48"))
        m["cosine_mean"] = m[["cosine_8", "cosine_48"]].mean(axis=1)
        m["cosine_min"] = m[["cosine_8", "cosine_48"]].min(axis=1)   # cross-condition robustness
        m["nd_min"] = m[["n_downstream_8", "n_downstream_48"]].min(axis=1)
        m = m[m.nd_min >= 20]                                    # signature-strength floor
        m = m.sort_values("cosine_mean", ascending=False).reset_index(drop=True)
        m["rank"] = np.arange(1, len(m) + 1)
        m.to_csv(f"{out}/screen_{anchor}_ranked.csv", index=False)
        results[anchor] = m
        # positive-control check
        r = m.index[m.gene == anchor]
        print(f"[{anchor}] self rank {int(m['rank'][r[0]]) if len(r) else 'NA'} "
              f"(cos {float(m['cosine_mean'][r[0]]):.3f})" if len(r) else f"[{anchor}] self filtered")
    f.close()
    print(f"Screens written to {out}/. Downstream axes (context, direction, "
          f"essentiality, degradability, GWAS) — see methods writeup.")
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out", default="./out")
    a = ap.parse_args()
    run(a.h5ad, a.out)
