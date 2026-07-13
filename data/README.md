# Data — precomputed checkpoints

The analysis runs on ONE object from the released atlas:

`GWCD4i.DE_stats.h5ad` (16.8 GB) — genome-wide DESeq2 differential-expression
statistics (33,983 perturbation×condition rows × 10,282 genes).

**Data citation.** Zhu R, Dann E, Yan J, et al. "Genome-scale perturb-seq in primary
human CD4+ T cells maps context-specific regulators of T cell programs and human
immune traits." bioRxiv (2025), doi:10.64898/2025.12.23.696273 (senior authors
J.K. Pritchard & A. Marson). Analysis code: github.com/emdann/GWT_perturbseq_analysis_2025.
Dataset page: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq

## Two ways to reproduce

**A. Full pipeline from raw data (authoritative).**
Download `GWCD4i.DE_stats.h5ad` from the dataset page above, place it here
(`data/GWCD4i.DE_stats.h5ad`), then:

    python pipeline/run_th2_screen.py

Expected tail of stdout (the built-in positive control):

    [GATA3] self rank 1 (cos 0.870)
    [STAT6] self rank 1 (cos 0.860)

**B. Fast highlights from checkpoints (no 16 GB download).**
The notebook `notebooks/highlights.ipynb` runs off two precomputed checkpoints
that are NOT committed here because of size (~1.3 GB total):

    data/de_zscore_matrix.npz   (~1.30 GB)  z-score matrix + gene names
    data/de_obs.parquet         (~2.4 MB)   per-perturbation metadata

Both are derived deterministically from the h5ad by `run_th2_screen.py`, or are
available as Claude Science artifacts (version IDs in the notebook's provenance
table). Drop them in this folder and the notebook reproduces every headline
number in ~15 s.
