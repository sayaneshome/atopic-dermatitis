# An intracellular master-regulator map of atopic dermatitis

*A genome-scale CD4⁺ T-cell Perturb-seq phenocopy screen for druggable nodes beyond STAT6.*

*Submitted to [Built with Claude — Life Sciences](https://cerebralvalley.ai/e/built-with-claude-life-sciences) (Cerebral Valley), a global hackathon of ~500 selected participants. Built end-to-end with Claude Science.*

> **Researcher Track — Build From the Bench.** We started from a bench question —
> *which intracellular regulators beyond STAT6 would collapse the allergic program if
> degraded?* — found the dataset built to answer it (the Marson lab's genome-scale
> CD4⁺ T-cell Perturb-seq atlas), and used Claude Science to turn it into a discrete,
> reproducible finding plus a druggability-annotated target map. One command
> re-derives the positive controls; every claim links to a figure or a primary paper.

---

## Summary (100–200 words)

Every approved atopic dermatitis (AD) drug blocks the type-2 pathway from the
**outside** — injectable antibodies against extracellular cytokines. The field is
moving **inward**, to the intracellular transcription factors that run the allergic
program, because degrading a master node promises deeper, oral suppression. STAT6 is
the proof of concept. But no one had systematically asked which *other* intracellular
regulators would collapse the program if shut down.

We screened a genome-scale CD4⁺ T-cell Perturb-seq atlas (~34,000 knockout signatures),
scoring every gene knockout for how strongly it phenocopies the STAT6/GATA3
"Th2-collapse" signature. The pipeline recovers both anchors at **rank 1** — a built-in
positive control — and surfaced an unexpected finding: **in circulating cells STAT6 and
GATA3 knockouts anti-correlate**, so GATA3 maintains the type-2 program
STAT6-independently. We expanded to every AD immune axis a T cell can represent
(JAK-STAT, Th22/AHR, Th1, activation), cross-referenced the known AD drug universe, added
literature, degrader-tractability, human-genetics and lesional-skin validation layers, and
ranked repurposing candidates against ClinicalTrials.gov. Deliverables: a
druggability-annotated target map, a reusable Th2-collapse classifier, and a fully
rerunnable pipeline.

---

## Reproduce the positive control in one command

```bash
python pipeline/run_th2_screen.py    # loads the released 16 GB h5ad, runs the screen
# expected tail of stdout:
#   [GATA3] self rank 1 (cos 0.870)
#   [STAT6] self rank 1 (cos 0.860)
```



---

## The finding, in one figure

**STAT6 ≠ GATA3.** Knocking out GATA3 collapses the type-2 cytokines (IL13/IL5/IL4);
knocking out STAT6 does not — it de-represses Th1 genes. In circulating committed CD4⁺
cells, GATA3 maintains the allergic program through a STAT6-independent autoregulatory
loop. Three independent methods agree (cosine screen, a supervised classifier at ROC-AUC
0.94 that ranks GATA3 #4/6923 and STAT6 #4809, and an autoencoder latent space), and the
result reproduces classical 2000–2004 immunology the pipeline was never given
(see `literature_validation.csv`).

---

## Repository map

```
README.md                    ← you are here (summary + quickstart)
METHODS_and_RESULTS.md       ← full 16-section methods & results write-up
SUBMISSION_researcher.md     ← written summary, 3-min video script, criteria map
LICENSE
notebooks/
  highlights.ipynb           ← fast reproduction off checkpoints (no 16 GB download)
pipeline/
  run_th2_screen.py          ← end-to-end, rerunnable from the 16 GB h5ad
  th2_collapse_scorer.py     ← reusable scorer (build_axis / phenocopy_score)
results/
  figures/                   ← 12 publication-grade figures (Nature-style, 300 dpi)
    fig_control_validation.png     ← START HERE — pipeline recovers its controls
    fig_network_topology.png       ← the STAT6≠GATA3 two-arm network
    fig_multiaxis_target_map.png   ← the 8-axis expansion
    fig_classifier_validation.png  ← ML confirmation (AUC 0.94)
    fig_lesional_validation.png    ← independent validation in AD skin (GSE147424)
    fig_translational_priority.png ← which axis is most translational
    fig_repurposing_candidates.png ← approved oral drugs on nominated nodes
    ... (+ latent UMAP, redundancy, TYK2 drill-down, shortlist, landscape)
  tables/                    ← 15 ranked / annotated CSVs
    SHORTLIST_master_ranked.csv        ← type-2 axis, ranked + annotated
    MULTIAXIS_integrated_shortlist.csv ← all 8 axes, ranked + annotated
    literature_validation.csv          ← 7 core claims vs. primary literature
    ... (screens, degradability, GWAS, repurposing, trials, landscape)
data/
  README.md  ← how to get the 16 GB h5ad + the checkpoints (not committed) due to large-size
```

## Data

`GWCD4i.DE_stats.h5ad` (16.8 GB) — genome-wide DESeq2 differential-expression
statistics from the released genome-scale CD4⁺ T-cell Perturb-seq atlas. 33,983
(perturbation × condition) × 10,282 genes. Not committed to the repo; download from the
[dataset page](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq)
and point `run_th2_screen.py` at it.

> **Cite the data:** Zhu R, Dann E, Yan J, *et al.* "Genome-scale perturb-seq in primary
> human CD4+ T cells maps context-specific regulators of T cell programs and human immune
> traits." *bioRxiv* (2025), doi:10.64898/2025.12.23.696273 (senior authors J.K. Pritchard
> & A. Marson). Analysis code: `github.com/emdann/GWT_perturbseq_analysis_2025`.

## Scope note

The primary screen profiles **circulating** CD4⁺ T cells, not skin-resident pathogenic
Th2 cells. We executed the named validation step — top hits confirmed in independent
lesional AD skin scRNA-seq (GSE147424); 8/18, including the GATA3 anchor, are enriched.
Barrier, itch-neuron, and epithelial arms of AD are out of scope for a T-cell atlas.

## Databases used (all via Claude Science connectors)

Open Targets (tractability, AD target/drug universe) · GWAS Catalog (allergic-disease
loci) · OpenAlex (literature depth) · ClinicalTrials.gov (trial reality-check) · NIH
RePORTER + PubMed (competitive landscape) · GEO (lesional-skin validation).

*Literature citations throughout are via PubMed; DOIs are in `literature_validation.csv`.*
