"""
th2_collapse_scorer.py
=======================
Reusable Th2-collapse phenocopy scorer for the genome-scale CD4+ T cell
Perturb-seq atlas (Marson lab, 2025; GWCD4i.DE_stats.h5ad).

A "collapse axis" is the differential-expression z-score signature of an
anchor knockout (GATA3 or STAT6) in the stimulated CD4+ T cell conditions,
restricted to that anchor's own significantly DE genes (its informative
subspace). Any query perturbation is scored by how strongly its own
z-score vector aligns with that axis on the subspace.

Key design decisions (see methods writeup):
  * Similarity is computed ONLY on the anchor's DE genes. The full
    10,282-gene vector is ~93% near-zero noise; scoring on it drowns the
    signal (STAT6 vs GATA3 full-vector r ~ -0.05, uninformative).
  * A consensus axis averages the two stimulated timepoints (Stim8hr,
    Stim48hr). Screening is run per-timepoint and combined (mean + min)
    so only cross-condition-robust hits rank high.
  * A signature-strength floor (min trans-effect genes in both timepoints)
    removes sparse-vector artifacts where a tiny program aligns by chance.

The scorer is anchor-agnostic: pass any reference z-vector + subspace mask.
"""
import numpy as np


def build_axis(zscore, obs, adj_p_getter, anchor, stim_conditions=("Stim8hr", "Stim48hr"),
               fdr=0.10):
    """Build a consensus collapse axis for an anchor gene.

    Parameters
    ----------
    zscore : (n_obs, n_genes) float array  -- the DE z-score layer
    obs : DataFrame with columns target_contrast_gene_name, culture_condition
    adj_p_getter : callable(row_index) -> (n_genes,) FDR-adjusted p-values
    anchor : str  -- anchor gene symbol (e.g. "GATA3", "STAT6")
    Returns
    -------
    dict {axis_z: (n_genes,), sig_mask: bool (n_genes,), n_informative: int}
    """
    rows = {}
    for c in stim_conditions:
        m = np.where((obs.target_contrast_gene_name.values == anchor) &
                     (obs.culture_condition.values == c))[0]
        if len(m):
            rows[c] = int(m[0])
    if not rows:
        raise ValueError(f"anchor {anchor} not found in any of {stim_conditions}")
    axis_z = np.vstack([zscore[i] for i in rows.values()]).mean(axis=0)
    sig = np.zeros(zscore.shape[1], dtype=bool)
    for c, i in rows.items():
        sig |= (adj_p_getter(i) < fdr)
    return {"axis_z": axis_z, "sig_mask": sig, "n_informative": int(sig.sum())}


def phenocopy_score(query_z, axis_z, sig_mask):
    """Cosine and Pearson similarity of query z-vector(s) to an axis on its subspace.

    query_z : (n_genes,) or (n_pert, n_genes)
    Returns (cosine, pearson) arrays of length n_pert.
    """
    q = np.atleast_2d(query_z)[:, sig_mask].astype(np.float64)
    a = axis_z[sig_mask].astype(np.float64)
    qn = q / (np.linalg.norm(q, axis=1, keepdims=True) + 1e-9)
    an = a / (np.linalg.norm(a) + 1e-9)
    cosine = qn @ an
    qc = q - q.mean(axis=1, keepdims=True)
    ac = a - a.mean()
    pearson = (qc @ ac) / (np.linalg.norm(qc, axis=1) * np.linalg.norm(ac) + 1e-9)
    return cosine, pearson


def score_perturbation(query_z, axis):
    """Convenience: score a single perturbation z-vector against a built axis dict.
    Returns dict {cosine, pearson}."""
    cos, pear = phenocopy_score(query_z, axis["axis_z"], axis["sig_mask"])
    return {"cosine": float(cos[0]), "pearson": float(pear[0])}
