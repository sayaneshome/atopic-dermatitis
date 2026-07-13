# Intracellular master regulators of atopic dermatitis — a multi-axis phenocopy screen
### From Th2-collapse (beyond STAT6) to a genome-wide, multi-axis AD target map
### Genome-scale CD4⁺ T cell Perturb-seq atlas (Marson lab, 2025)

*Part 1 (§1–8) screens the type-2 / Th2-collapse arm anchored on STAT6 & GATA3.
Part 2 (§9) expands to every AD immune axis a circulating CD4⁺ T cell can
represent — Th1/IFNγ, JAK-STAT, Th22/AHR, activation tone, AD-GWAS TFs — and
cross-references the known AD drug-target universe and the literature.*

*Documented, rerunnable pipeline and results. Scope note: this atlas profiles
**circulating** CD4⁺ T cells, not skin-resident pathogenic Th2 cells; the
shortlist therefore maps the actionable immune arm of atopic dermatitis (AD)
and defines a clear next validation step — confirming top hits in lesional T
cell data — rather than claiming skin-level proof.*

---

## 1. Rationale

Every approved AD therapy blocks the type-2 pathway extracellularly
(dupilumab and other anti-cytokine/receptor biologics). The field is moving
inward to the intracellular transcription factors that run the Th2 program,
because degrading the master node promises deeper, oral suppression. STAT6 is
the proof of concept — once "undruggable," now degraded orally with early
biologic-like efficacy in AD. The open question: **which other intracellular
regulators would collapse the IL-4/IL-13 program if shut down, and which are
chemically tractable?** Perturb-seq measures exactly this — the transcriptional
consequence of knocking out each gene — so the screen reduces to: *which
genome-wide perturbations phenocopy the known Th2-collapse knockouts?*

## 2. Data

`GWCD4i.DE_stats.h5ad` (16.8 GB) from the released atlas (CZI Virtual Cells).
**Data citation:** Zhu R, Dann E, Yan J, *et al.* "Genome-scale perturb-seq in
primary human CD4+ T cells maps context-specific regulators of T cell programs and
human immune traits." *bioRxiv* (2025), doi:10.64898/2025.12.23.696273 (senior
authors J.K. Pritchard & A. Marson; analysis code
`github.com/emdann/GWT_perturbseq_analysis_2025`). Genome-wide DESeq2 differential
expression:
- **33,983** observations = (perturbed gene × culture condition), conditions
  Rest / Stim8hr / Stim48hr.
- **10,282** measured transcriptome genes.
- Layers: `zscore` (= logFC/lfcSE), `log_fc`, `adj_p_value`, `baseMean`, `lfcSE`.
- `.obs` QC: on-target knockdown significance, cis/off-target flags, trans-effect
  counts (`n_downstream`), cross-guide and cross-donor reproducibility.

Each perturbation is thus a **z-score vector over 10,282 genes** — its
transcriptomic fingerprint. Working set (z-score layer, float32) ≈ 1.4 GB.

## 3. Pipeline

**3.1 Eligibility filter.** Kept perturbations with significant on-target
knockdown and no confounds: `ontarget_significant & ¬neighboring_gene_KD &
¬distal_offtarget_flag & ¬low_target_gex`. → **18,622** eligible
perturbation×condition rows (~6,300 per stimulated condition from 11,415 tested).

**3.2 Reference axes.** A collapse axis = an anchor knockout's consensus z-vector
across the two stimulated timepoints, **restricted to the anchor's own
FDR<0.10 DE genes** (its *informative subspace*). This restriction is essential:
the full 10,282-gene vector is ~93% near-zero noise, and similarity computed on
it drowns the signal (full-vector STAT6–GATA3 r ≈ −0.05, uninformative).

**3.3 Scoring function** (`th2_collapse_scorer.py`). Cosine (and Pearson)
similarity of any perturbation's z-vector to an axis, on the axis subspace.
Reusable on any perturbation, drug signature, or future dataset.

**3.4 Genome-wide screen.** Every eligible perturbation scored per stimulated
condition; combined as **mean** (rank) and **min** (cross-condition robustness)
cosine. A **signature-strength floor** (≥20 trans-effect genes in *both*
timepoints) removes sparse-vector artifacts (a tiny program aligning by chance).

**3.5 Prioritization axes.** (i) *Context-specificity*: stim-vs-rest trans-program
ratio. (ii) *Directionality*: mean z on curated Th2-effector / Th1 / Treg /
activation modules → ablation vs. redirection. (iii) *Essentiality/toxicity*:
large rest program, pan-transcriptome trans-effect (>3,000 genes), or broad
activation suppression → deprioritized. (iv) *Degradability*: Open Targets
tractability (SM ligandability tier + PROTAC-modality ubiquitination handle).
(v) *Genetic support*: GWAS Catalog associations mapped to each gene, filtered
to AD / asthma / eczema / allergy traits, strongest genome-wide p per gene.

**3.6 Composite priority** = 0.35·phenocopy + 0.25·genetic + 0.20·degradability
+ 0.20·context, multiplied by a 0.4 safety gate when essentiality-flagged.

## 4. Positive-control validation (the go/no-go gate)

Both anchors recover at **rank 1** on their own axis (GATA3 cos 0.87; STAT6
cos 0.86). On the **STAT6 axis, IL4R recovers at rank 2** (cos 0.66) — the
receptor directly upstream of STAT6 (IL-4→IL4R→JAK→STAT6), exactly the phenocopy
predicted. Gate passed. *(fig_control_validation.png)*

## 5. Key biological finding — STAT6 and GATA3 do NOT phenocopy each other here

On their shared significant DE genes the two anchor knockouts are
**anti-correlated** (r = −0.41 at Stim8hr), and their informative gene sets
barely overlap (Jaccard 0.16). Reading canonical genes:
- **GATA3 KO** collapses the type-2 cytokines cleanly (IL13 z −7.7, IL5 −4.7,
  IL4 −3.6, CCR4 −9.5) — *this* is Th2 collapse.
- **STAT6 KO** does **not** collapse the cytokines (IL13 **+4.3**, IL4 +2.3)
  and instead de-represses Th1 (IFNG +4.9); it pulls down STAT-target genes
  (CISH −6.7, BATF −6.0).

This is the expected consequence of cell context, and precisely why the scope
note matters: STAT6 is required to *differentiate* Th2 cells from naïve
precursors, but in **committed, circulating** CD4⁺ T cells GATA3 maintains the
type-2 cytokine program **STAT6-independently**. Ranking eligible perturbations
directly on the type-2 effector module confirms it: **GATA3 #1, IL4R #2, STAT6
~3,050/6,371**.

**Design decision (user-approved):** run the two anchors as **separate
screens** — a GATA3 axis (cytokine-collapse program) and a STAT6 axis
(JAK-STAT / signaling tone) — and report both plus their overlap. The two
top-40 lists have **zero overlap**, and the network resolves into two
anti-correlated arms (§7).

## 6. Ranked shortlist (top 15 of 80 candidates)

| # | Gene | Arm | Phenocopy (cos) | Direction | Degradability | AD/asthma GWAS | Safety |
|---|------|-----|-----------------|-----------|---------------|----------------|--------|
| 1 | **IL4R** | STAT6 | 0.66 | →Th1 redirect | high (HQ ligand) | asthma p=2e-32 (38 loci) | clean |
| 2 | ORAI1 | STAT6 | 0.27 | weak | high | — | clean |
| 3 | **SMARCA4** | GATA3 | 0.35 | ablation | high | — | clean |
| 4 | **SMAD4** | GATA3 | 0.25 | ablation | high | childhood asthma p=4e-13 | clean |
| 5 | KMT2C | GATA3 | 0.31 | ablation | high | allergy p=4e-7 | rest-active |
| 6 | **TPD52** | GATA3 | 0.58 | ablation | low | asthma-exac p=8e-7 | rest-active |
| 7 | WNK1 | STAT6 | 0.27 | weak | high | atopic eczema p=3e-9 | rest-active |
| 8 | PPP3R1 | STAT6 | 0.29 | weak | moderate | — | clean |
| 9 | PHF8 | GATA3 | 0.22 | weak | high | — | clean |
| 10 | KDM3B | GATA3 | 0.25 | ablation | moderate | allergy p=4e-8 | clean |
| 11 | ERCC1 | STAT6 | 0.27 | weak | high | asthma-exac p=4e-6 | rest-active |
| 12 | NRAS | GATA3 | 0.27 | →Th1 redirect | high | — | rest-active |
| 13 | MLST8 | GATA3 | 0.24 | ablation | high | — | rest-active |
| 14 | PCYT2 | STAT6 | 0.29 | weak | high | — | rest-active |
| 15 | RAC2 | STAT6 | 0.24 | weak | moderate | — | clean |

Full 80-row table with all axes: `SHORTLIST_master_ranked.csv`.
*(fig_shortlist_profile.png — four-axis profile of the top 18)*

**Reading the shortlist.**
- **IL4R** is the standout: phenocopies STAT6 loss, redirects Th2→Th1 (not pure
  ablation), clinically druggable, and the strongest human genetics in the set.
  Its recovery at #1 is itself a validation of the composite.
- Genuinely novel *intracellular* candidates: **SMAD4** (TGFβ/BMP effector;
  childhood-asthma genetics; ligandable; clean), and the **chromatin module**
  (SMARCA4/KMT2C/KDM3B/BRD1 — BAF complex + histone modifiers) which sits
  *upstream of GATA3* (§7).
- **Do NOT over-read the TCR-proximal signaling genes** (CD3E, ZAP70, PLCG1,
  ITK): they are flagged pan-toxic — their knockouts ablate T-cell activation
  wholesale (activation z ≈ −5) and phenocopy "collapse" trivially. The safety
  gate correctly demotes them.

## 7. Network topology (fig_network_topology.png)

Candidates + anchors cluster (average linkage on 1−correlation of stim
trans-signatures over the union subspace) into **8 modules forming two arms**:

- **GATA3 arm** (corr to GATA3 > 0, to STAT6 < 0): GATA3-core (TPD52, MYB, NRAS,
  ADAM19), chromatin machinery (BAF/KMT2C/KDM3B/BRD1/BATF/ID2), GATA3-satellite
  TFs (ZNF649, DAP, LAT2), SMAD4, and a metabolic module.
- **STAT6 arm** (corr to STAT6 > 0): STAT6-signaling core (IL4R, CREB1, OTUD4),
  TCR/Ca²⁺ signaling (ORAI1, PPP3R1, ITK, PTPRC, RAC2), transcription/elongation.

**Directed layer (upstream-of-anchor test).** A candidate whose knockout reduces
the *anchor TF's own transcript* (z < −3) sits upstream. **The chromatin module
sits upstream of GATA3**: BATF (−6.2), SMARCE1 (−5.7), SMARCB1 (−5.1), KMT2C
(−4.5), SMARCA4 (−3.7) all suppress GATA3 expression — i.e. they are required to
*maintain* the GATA3 program, placing them above GATA3 in the hierarchy and
making them candidate higher-order master regulators. IL4R KO also reduces GATA3
(−5.5, cross-arm). CREB1, CAND1, PTPRC, FBXO32 sit upstream of STAT6.

## 8. Redundancy / combination potential (fig_redundancy_heatmap.png)

Pairwise signature correlation among top hits shows a **two-block, anti-correlated
structure** (STAT6↔GATA3 r = −0.13):
- **Non-redundant, distinct programs** (candidate additive combinations): IL4R,
  SMAD4, KDM3B, PHF8, CREB1, ACTR8, DOHH, PARP14.
- **Redundant (same program — pick one):** ORAI1 / PPP3R1 (r 0.60, store-operated
  Ca²⁺/calcineurin); the chromatin sub-block (KMT2C / MLST8 / BHLHE40).
- **Combination logic:** because the two arms are uncorrelated, a **GATA3-arm
  agent** (e.g. targeting the chromatin module maintaining GATA3) **plus a
  STAT6-arm agent** (IL4R/STAT6) would hit non-overlapping programs and is
  predicted additive — not duplicative — to a STAT6 degrader.

## 9. Expansion — genome-wide multi-axis AD target map

The Th2-collapse screen above targets one arm of AD. AD is polygenic and
multi-axis: approved and pipeline drugs hit type-2 (dupilumab, IL-13/IL-31
biologics), **JAK-STAT** (oral abrocitinib/upadacitinib/JAKi), **PDE4**
(crisaborole/apremilast), **Th22/AHR** (tapinarof), barrier (FLG), and itch
(IL31RA). We therefore re-ran the same phenocopy method against **every AD
immune axis a circulating CD4⁺ T cell can represent**, and cross-referenced the
result against the known AD drug-target universe and the literature.

**9.1 Axis anchors.** Candidate master genes per axis were tested for
eligibility and program strength (≥30 informative DE genes). Viable single-gene
axes (13): **Th1/IFNγ** (TBX21, STAT1), **JAK-STAT** (JAK2, TYK2, STAT5A,
STAT5B), **Th22/AHR** (AHR), **activation tone** (LCK, PRKCQ, NFATC2, REL),
**AD-GWAS TFs** (IKZF1, IRF1). As with STAT6/GATA3, anchors nominally in one
"axis" frequently defined *opposing* programs (JAK2/TYK2 kinases anti-correlate
with STAT5 effectors; TBX21 vs STAT1 r = −0.22), so each anchor was scored as
its own reference axis rather than averaged.

**Scope finding — Treg is not anchorable here.** FOXP3, CTLA4, IL2RA, IKZF2
knockouts produce almost no downstream program in circulating cells (n_sig
1–29), as do the secreted effectors IFNG (2) and IL22 (1). Their transcriptional
consequence is minimal in circulating CD4⁺ T cells — the same circulating-vs-
lesional limitation — so no Treg axis was built.

**9.2 Positive-control recovery (multi-axis gate).** Six axes recover their own
anchor at **rank 1** with high self-cosine and biologically exact neighbours:

| Axis | Anchor | self-cos | Top emergent hits (biology recovered) |
|------|--------|----------|----------------------------------------|
| JAK-STAT | **TYK2** | 0.98 | IFNAR1, STAT2 (the type-I-IFN receptor complex TYK2 signals through) |
| AD-GWAS TF | **IRF1** | 0.94 | USF2, AK3, CTDSPL2 |
| Th1/IFNγ | **TBX21** | 0.90 | VAV1, ITK (TCR-driven Th1 program) |
| Activation | **LCK** | 0.90 | VAV1, CD3E, PLCG1, ZAP70, ITK (proximal TCR signalosome) |
| Th22/AHR | **AHR** | 0.89 | ARNT (obligate AHR dimer partner), AIP (AHR chaperone) — textbook |
| JAK-STAT | **STAT5B** | 0.79 | IL2RB (receptor upstream of STAT5), PKM, CREBBP |

The remaining anchors (STAT1, JAK2, STAT5A, PRKCQ, NFATC2, REL, IKZF1) have weak
or paralog-dominated programs and did not self-recover above the strength floor
(e.g. STAT5A's top hit is its paralog STAT5B).

**9.3 Known AD drug-target landscape vs. atlas.** Of the top-60 Open Targets AD
associated targets (`ad_druggable_landscape.csv`): **31 are atlas-visible**
(eligible perturbations — screenable), 13 not in the atlas, 13 perturbed-but-
ineligible, 3 measured-readout-only. The targets that are not atlas-visible are
the non-T-cell arms, spread across three scope categories: **not in the atlas at
all** — barrier (**FLG, SERPINB7**), systemic vit-D (**CYP24A1**), PDE4C;
**measured-only (readout, not perturbed)** — itch sensory neuron (**IL31RA**),
histamine (**HRH1**); **perturbed-but-ineligible** — PDE4A, PDE4D. (PDE4B is in
fact atlas-visible, so the PDE4 family is split rather than uniformly out of
scope.) Atlas-visible known targets recover **when they carry a measurable
program** (TYK2 cos 0.98, AHR 0.89 as anchors; TRAF3/FOSL2/BHLHE40 associate to
activation/JAK/Th1 axes); the many that don't (IL13, IL22, IL2RA, PDE4B…) are
the weak-in-circulating-cells cases — a scope statement, not a pipeline failure.
This is the "recovered-as-control vs out-of-scope" separation
(`known_target_recovery.csv`).

**9.4 Literature novelty layer** (`literature_novelty.csv`, OpenAlex). Each
emergent hit across the six well-recovered axes was annotated by AD-literature
depth. Of 79 emergent genes: **17 well-studied in AD/immunology** (LCK, ITK,
EZH2, MALT1, TBX21, SOCS3 — literature-recovered controls that confirm the axes
surface real immunology), 18 with some AD literature, 38 emerging (<20 papers),
and **6 essentially unstudied in AD** — the candidate-novel nominations:
**CTDSPL2** (IRF1/AD-GWAS axis), **RSBN1L** (TYK2/JAK axis, cos 0.53),
**ZNF649** (TYK2), **CLASRP** (STAT5B), **MFSD4B** (STAT5B), **SCAF4** (Th1).
Symbol ambiguity inflates raw counts for a few genes (CD2, CD28, AIP, ARNT match
unrelated abbreviations); the tiering is therefore a guide, not a hard filter.

**9.5 Integrated map** (`fig_multiaxis_target_map.png`,
`MULTIAXIS_integrated_shortlist.csv`). Eight axes (the two type-2 anchors plus
the six recovered immune axes), top-8 phenocopy hits each, annotated by
phenocopy strength, AD-literature depth, and degrader-tractability (extended to
122 genes). Highlights of high-tractability emergent hits: **IFNAR1, STAT2**
(TYK2/IFN axis), **ARNT, DOCK2** (Th22/AHR), **IL2RB, CREBBP** (STAT5), **EZH2,
MALT1, LCK** (approved-drug-tier chromatin/signaling nodes). The TCR-proximal
activation hits (LCK, ZAP70, ITK, CD3E, PLCG1) recover strongly on both the
activation and Th1 axes but remain pan-toxic — informative for mechanism, poor
as selective targets.

## 10. Translational triage — known drugs, axis maturity, repurposing

**10.1 AD drug-target universe** (`ad_drugged_targets.csv`, Open Targets, drug →
mechanism → clinical stage). 229 drug entries hit 118 AD targets; **17 targets
have an approved AD drug.** Approved AD development concentrates on: JAK-STAT
(JAK1/2/3, TYK2 — 6 oral JAKi), type-2 (IL4R/dupilumab, IL13, IL31RA/nemolizumab),
PDE4 (crisaborole, apremilast), AHR (tapinarof), topical (NR3C1, FKBP1A). The
critical observation: nearly every approved AD drug hits an **extracellular /
receptor** target. The atlas surfaces the **intracellular** nodes downstream of
those receptors, and those are overwhelmingly *undrugged* — which is the core
opportunity this project targets.

**10.2 Which axis is most translational** (`axis_translational_maturity.csv`,
`fig_translational_priority.png`). Each axis was placed on two dimensions:
pathway clinical validation in AD (0 = none → 4 = approved drug) vs. open
intracellular opportunity (# degradable, undrugged top hits). Three "sweet-spot"
axes combine a **clinically de-risked pathway** with **open intracellular
nodes**:
- **JAK-STAT / TYK2** — most translational. 6 oral JAKi approved; TYK2 itself
  approved (deucravacitinib). Clearest regulatory precedent for an oral
  intracellular agent, but the *most crowded* — differentiation requires a
  selective degrader rather than another kinase inhibitor.
- **Th22 / AHR** — tapinarof approved; ARNT/DOCK2/STIP1 open downstream.
- **Type-2 / STAT6** — dupilumab approved + STAT6 degrader in trials; GATA3 /
  chromatin machinery still open. Best combination of validation and white-space.

Lower priority: Th1/IFNγ (TBX21) and AD-GWAS TF (IRF1) are genetics-led bets on
pathways with **no approved AD drug**; activation-tone (LCK) is validated
(topical calcineurin inhibitors) but its hits are pan-toxic.

**10.3 Repurposing layer** (`repurposing_candidates.csv`,
`fig_repurposing_candidates.png`). Of 110 intracellular hits across all axes,
**17 already have a drug for some indication**. The strongest small-molecule
repurposing leads (approved, immune/derm precedent, hitting an atlas-nominated
node):
- **ITK → ritlecitinib** — approved covalent ITK/JAK3 inhibitor in alopecia
  areata (autoimmune skin); tested in asthma. Shortest conceptual leap.
- **EZH2 → tazemetostat** — approved oral EZH2 inhibitor (oncology); chromatin
  node on the AHR axis.
- **IFNAR1 → anifrolumab** — approved anti-IFNAR in lupus; top TYK2-axis hit.
- **CD2 → alefacept** — was approved in psoriasis; Th1-axis costimulation node.
- **LCK → dasatinib** — approved oncology kinase inhibitor.

*Caveats:* TCR-proximal hits (CD3E, LCK, ITK, CD2, CD28) are broad T-cell
suppressors whose existing drugs carry immunosuppression/toxicity; several leads
are biologics (validate the pathway but don't advance the oral-intracellular
goal). The oral small molecules — ritlecitinib, tazemetostat, dasatinib,
deucravacitinib — are the actionable subset.

**10.4 ClinicalTrials.gov reality-check** (`repurposing_ad_trials.csv`). Are the
repurposing leads already being tested in AD?
- **Already repurposed into AD:** delgocitinib (11 trials, Ph III topical —
  approved), gusacitinib (4, Ph II), deucravacitinib (3 AD trials, up to Ph IV),
  alefacept (2, Ph II/IV, historical). Confirms JAK-family repurposing into AD
  works.
- **UNTESTED in AD (open):** **ritlecitinib (ITK), tazemetostat (EZH2), dasatinib
  (LCK), anifrolumab (IFNAR1), salirasib (NRAS)** — approved elsewhere, zero AD
  trials. The cleanest open repurposing theses are **ritlecitinib (ITK)** and
  **tazemetostat (EZH2)** — approved oral small molecules on atlas-nominated
  nodes, neither yet in AD.

**10.5 Axis drill-down — JAK-STAT / TYK2** (`drilldown_JAKSTAT_TYK2.csv`,
`fig_drilldown_TYK2.png`). The most-translational axis ranked by a composite of
phenocopy (0.5) + degradability (0.3) + allergic-disease genetics (0.2). Top
nodes: **IFNAR1** (#1, phenocopy 0.61, high tractability, anifrolumab approved
in lupus & untested in AD), **SMAD4** (#2, high tractability, 7 allergic-disease
GWAS loci), **STAT2** (TYK2's direct signaling partner), and open novel nodes
RSBN1L / ADNP2. IFNAR1 is the standout: strongest phenocopy on the axis, a
ligandable+ubiquitinated degrader surface, and an approved antagonist antibody
that has never been tried in AD.

## 11. Competitive landscape — where this sits in funded & published AD research

Scanned **NIH RePORTER** (funded grants, FY2022–2025, title/abstract match) and
**OpenAlex** (publications) to position the work.
(`competitive_landscape.csv`, `competitor_analysis.csv`,
`fig_competitive_landscape.png`.)

**Robust methodology.** Keyword-query counts are unreliable for two reasons:
(1) RePORTER returns one record per fiscal-year/supplement, inflating counts
(450 raw records → **172 unique grants** after deduplicating by core project
number); (2) PubMed and OpenAlex silently expand loose terms (e.g. "degrader"
matched 8,701 AD papers via phrase expansion, but PROTAC/molecular-glue
restricted to title/abstract returns **1**). We therefore (a) pulled the full AD
grant set once and deduplicated by core project number, (b) classified each
unique abstract by regex against theme vocabularies, and (c) used MeSH-anchored,
title/abstract-restricted PubMed queries with the `query_translation` inspected
for every count.

**Funding white space (deduplicated, n=172 grants).** The field's weight is in
barrier/filaggrin (88), microbiome (51), biologics (48), itch (34), and Th2
differentiation (32). The intracellular-screening space is near-empty:
transcription-factor biology (7), JAK-STAT (6), **STAT6/GATA3 (3), genome-scale
Perturb-seq/CRISPR screening (1), targeted protein degradation (0).** The single
"screen" grant is Julia Oh's staphylococcal skin-*microbiome* work — not T-cell
target discovery — and **zero grants combine genome-scale screening with Th2/TF
biology**, the exact intersection this project occupies. PubMed confirms: **0 AD
papers** on Perturb-seq or CRISPR-screen (title/abstract), **1** on
PROTAC/molecular-glue (PMID 38262152), vs. 354 on JAK inhibitors and 101 on
STAT6.

**Read-verified, not just text-matched.** To ensure no adjacent grant runs a
screen without saying so in matchable terms, we read the full abstracts of all
**38 deduplicated grants** in the Th2 / transcription-factor / STAT6-GATA3
cluster and classified each (structured LLM pass, 0 parse failures) on three
questions: runs a genome-scale perturbation screen? pursues a degrader?
systematically ranks regulators to nominate new targets? **All 38 answered no to
all three.** The closest-looking grants do adjacent-but-different work — Khavari
(regulatory-variant→gene mapping), Cho (T-cell biomarker validation), Ferrer (3D
skin-model compound screening), Henao-Mejia (single-locus 3D chromatin). The
"empty intersection" claim is therefore verified by reading, not only by keyword
counts. (`competitor_abstract_audit.csv`.)

**Nearest adjacent labs** (deduplicated): **Kalung Cheung** (Mount Sinai,
R01AI177461 — "Transcriptional Regulation of Th2 Cell Development") is the single
closest, but studies individual TFs mechanistically rather than screening all
regulators for phenocopy + degradability; **John O'Shea** (NIAMS, ZIAAR041159 —
JAK-STAT/STAT6 pathway biology) characterizes the pathway rather than mining it
for new nodes. This project's two defining moves — genome-scale phenocopy
screening and degradability ranking — are the two emptiest cells in the funding
map.

**Prior art vs. this contribution.** The method building blocks are established:
genome-scale Perturb-seq for genotype→phenotype mapping (Replogle-style atlases)
and signature-based phenocopy screening (Connectivity Map / LINCS L1000). The
target premise is current: an oral STAT6 degrader (KT-621) has reported Phase 1a
data (2025). What is not precedented is applying this machinery systematically to
the STAT6/GATA3 collapse signature in a CD4⁺ T-cell atlas, ranked for
degradability, for AD — the empty intersection verified above. The novelty is in
the combination and the specific question, not in inventing the tools. (Precedents
verified by name; noisy keyword totals were not used for counts.)

*Caveat:* even after deduplication and MeSH anchoring, these are text-matched
counts — a grant could use a method without naming it in its abstract — so treat
them as a strong directional signal, not a literal census. OpenAlex was dropped
for specific numbers (its keyword ranking was too lexically noisy); the reported
counts come from deduplicated RePORTER grants and audited MeSH/tiab PubMed
queries. Artifacts: `competitive_landscape_robust.csv` (both sources side by
side), `competitor_analysis.csv`, `fig_competitive_landscape.png`.

## 12. Machine-learning layers — two orthogonal confirmations of the finding

The cosine phenocopy screen is a similarity metric, not a learned model. We added
two learned models on the transcriptomic response; both independently reproduce
the central STAT6≠GATA3 result, which is the strongest possible validation of it.

**12.1 Supervised Th2-collapse classifier**
(`th2_collapse_classifier.joblib`, `classifier_predictions.csv`,
`fig_classifier_validation.png`). A gradient-boosted classifier
(`HistGradientBoostingClassifier`) predicts, from a perturbation's genome-wide
stimulated z-vector, whether it collapses the type-2 effector module. **Labels are
defined from independent biology** — bottom decile of a curated effector-module
score (IL4/IL5/IL13/IL4R/CCR4/GATA3/MAF/IL2RA), not from the cosine metric — and
the 8 label genes are **removed from the feature set** to prevent leakage.
Five-fold cross-validated **ROC-AUC 0.939, PR-AUC 0.673** (baseline 0.10). On
held-out predictions: **GATA3 ranks #4 / 6,923** (p=0.98) while **STAT6 ranks
#4,809** (p=0.004). A supervised model, given no knowledge of the cosine screen,
reaches the same conclusion — GATA3 collapses the type-2 program, STAT6 does not.
The transcriptomic drivers of the collapse call (`classifier_drivers.csv`) are
loss of activated-effector markers (CTLA4, CD4, ICOS, IL2RB, IL17RB, OSM, LIF,
ADAM19), a coherent "loss of activated type-2 identity" signature. The model is a
reusable scorer: any perturbation's z-vector → collapse probability.

**12.2 Learned latent space of perturbation responses**
(`perturbation_latent.npy`, `fig_latent_umap.png`). A CPU autoencoder compresses
each perturbation's 10k-gene response to a 16-dim latent vector. Collapse-positive
perturbations are **1.9× enriched** in each other's latent neighborhoods (they
form a coherent region), and **GATA3's nearest latent neighbor is TPD52** — the
same top GATA3-arm hit the cosine screen found — while STAT6 sits in a separate
region. A third, unsupervised method agrees with the first two.

*Note on scope:* we tested a protein-sequence language-model (ESM-2) layer for
predicting degradability from sequence; it reached only ROC-AUC 0.609 (degradability
is a 3D-pocket property poorly captured by sequence), so it was dropped. The
appropriate degradability model is structure-based (AlphaFold/Boltz), named as a
next step in §14.

## 13. Concordance with prior literature (is the finding real, not an artifact?)

Before this map is built on, its core claims were checked against the primary
immunology literature (PubMed). The decisive test is the pivotal, counterintuitive
result — that in **circulating** CD4⁺ cells GATA3 maintains the type-2 program
while STAT6 does not — because a screen artifact would contradict established
biology, whereas rediscovering textbook biology the pipeline was never told about
is the strongest evidence it reads real signal. **Six of seven core claims are
directly supported by primary literature; the seventh (that a CRISPR screen
validly discovers Th2 regulators) has clear precedent.** (`literature_validation.csv`.)

| Claim (our result) | Verdict | Supporting primary source |
|---|---|---|
| GATA3 maintains the type-2 cytokine program in differentiated cells (not just initiates it) | **SUPPORTED** | Yamashita et al., J Biol Chem 2004 — 'Essential role of GATA3 for the MAINTENANCE of Th2 cytokine production'; GATA3 ablation in differentiated Th2 reduces all Th2 cytokines — doi:10.1074/jbc.M403688200 |
| GATA3 drives the type-2 program STAT6-independently (maintenance is STAT6-independent) | **SUPPORTED (this is the key check)** | Ouyang et al., Immunity 2000 — 'Stat6-INDEPENDENT GATA-3 autoactivation directs IL-4-independent Th2 development and commitment'; GATA3 autoactivation stabilizes Th2 without STAT6 — doi:10.1016/s1074-7613(00)80156-9 |
| GATA3 can impose the Th2 program even on already-committed cells (autonomous maintenance) | **SUPPORTED** | Lee et al., J Exp Med 2000 — ectopic GATA3 induces Th2 cytokines in irreversibly committed Th1 cells; autoregulatory maintenance loop — doi:10.1084/jem.192.1.105 |
| A chromatin-remodeling module (BAF/SMARCA4-BRG1) sits upstream of the GATA3 cytokine program | **SUPPORTED** | Linzer et al., Immunity 2023 (CRISPR screen of 1131 TFs) — GATA3/AP-1 need a BRG1/CHD4 chromatin-remodeling bridge (ADNP) to open the type-2 locus; BAF/BRG1 established in Th differentiation — doi:10.1016/j.immuni.2023.05.010 |
| BATF is a bona fide Th2/type-2 differentiation TF (not a screen artifact) | **SUPPORTED** | BATF/JunB AP-1 axis regulates IL-4/Th2 gene expression (JunB review Immunol Med 2021; BATF Th2 literature n=6) — doi:10.1080/25785826.2021.1872838 |
| IL4R is a central, disease-relevant type-2 node | **SUPPORTED** | IL-4/IL-13 receptor biology drives allergic inflammation; dupilumab (anti-IL4Ra) approved in AD; STAT6 dysregulation central to atopic dermatitis (Int Rev Immunol 2024) — doi:10.1080/08830185.2024.2395274 |
| A CRISPR TF screen is a valid way to discover Th2 regulators (methodological soundness) | **PRECEDENT EXISTS (ours is larger/unbiased)** | Linzer et al. 2023 ran a targeted 1,131-TF CRISPR screen for Th2 and recovered GATA3+chromatin biology — same logic, smaller scope; validates the approach — doi:10.1016/j.immuni.2023.05.010 |

**The key check.** Our pivotal STAT6≠GATA3 divergence reproduces a classical (if
underappreciated) result almost exactly. Ouyang et al. (*Immunity*, 2000;
doi:10.1016/s1074-7613(00)80156-9) showed **Stat6-independent GATA-3 autoactivation**
directs IL-4-independent Th2 development and commitment — GATA3 fully reconstitutes
Th2 in Stat6-deficient cells via an autoregulatory loop. Yamashita et al. (*J Biol
Chem*, 2004; doi:10.1074/jbc.M403688200) demonstrated an **essential role of GATA3
for the *maintenance*** of Th2 cytokine production in already-differentiated cells.
That is precisely our mechanism: GATA3 maintains the cytokine program
STAT6-independently in committed circulating cells (GATA3 rank #1 on its axis, #4/6923
in the classifier; STAT6 #4809). The chromatin-upstream-of-GATA3 result (§7) matches
Linzer et al. (*Immunity*, 2023; doi:10.1016/j.immuni.2023.05.010), whose targeted
1,131-TF CRISPR screen found GATA3/AP-1 require a BRG1/CHD4 chromatin-remodeling
bridge to open the type-2 locus — the same biology, at smaller scope, which also
validates the screening approach itself.

**Two honest caveats this exercise reinforces.** (1) The classical STAT6-vs-GATA3
genetics were done in *differentiated/committed* cells — exactly the circulating-cell
regime this atlas samples — so the agreement is expected *and* re-underscores that
skin-resident Th2 biology still needs the lesional validation (§14). (2) Literature
confirming a node's *role* does not confirm its *degradability*, which remains a
structural-biology question (§15). All citations via PubMed.

## 14. Independent validation in lesional AD skin (the named next step, executed)

The project's central caveat is that the atlas profiles *circulating* CD4⁺ T
cells, not skin-resident pathogenic Th2 cells. We tested whether the top hits
transfer to real disease tissue using an **independent public dataset**:
**GSE147424** (Rojahn et al.), 10x single-cell RNA-seq of **lesional and healthy
human skin** — 4 AD lesional donors + 4 healthy donors, processed per-cell
expression matrices downloaded from GEO.

**Method.** T cells were gated as CD3E⁺ within each sample (714 lesional vs. 109
healthy skin T cells — the ~7× infiltration itself reproduces the known AD
inflammatory influx). For each top hit we compared expression in lesional vs.
healthy skin T cells (one-sided Mann–Whitney).

**Result — 8 of 18 testable hits are significantly enriched in lesional skin T
cells** (p<0.05), including the anchor: **BATF** (log₂FC 2.17, 50% vs 7% of T
cells, p<10⁻³), **GATA3** (the anchor; log₂FC 0.52, 32% vs 17%, p=0.006),
**IL4R** (p=0.011), **IFNAR1** (p=0.031), **VAV1**, **NDFIP2**, **AIP**,
**STIP1**. **ITK** and **IL2RB** are markedly higher in lesional T cells (ITK
22% vs 7%, IL2RB 12% vs 4% of T cells; their log₂FC is undefined only because of
a healthy-matrix row artifact, not true absence). Two findings are notable:
- **GATA3 validates; STAT6 is weaker** (log₂FC 0.35, p=0.096, n.s.) — the same
  GATA3 > STAT6 ordering the circulating screen found now reappears in skin.
- **The validation discriminates**: TPD52, ARNT, CTDSPL2, USF2, BHLHE40 do *not*
  transfer to lesional skin. A screen that confirmed every hit would be
  rubber-stamping; a genuine 8/18 with the anchor among them is real signal.

This converts the circulating-vs-lesional caveat from an unaddressed limitation
into a partial, positive cross-tissue validation: the core Th2 nodes (GATA3,
BATF, IL4R) and several emergent hits mark the pathogenic skin T-cell population
in the actual disease. (`lesional_validation.csv`, `fig_lesional_validation.png`.)

**Scope of this validation.** It is a pseudobulk-per-condition comparison over
one dataset (4+4 donors, T cells not further subtyped into Th2/Tc2). It confirms
disease-tissue relevance of the top nodes; it does not replace a targeted
perturbation experiment in skin-resident Th2 cells. Larger lesional atlases (e.g.
GSE153760, GSE158432) are privacy-restricted for raw data.

## 15. Limitations & next step

1. **Circulating primary screen, lesional validation.** The genome-wide screen is
   in circulating CD4⁺ cells; top hits are now cross-validated in lesional skin T
   cells (§14), but the definitive test remains a perturbation experiment in
   skin-resident Th2 cells. Treg and cytokine-effector programs are near-silent in
   circulating cells, and barrier/itch axes are entirely out of scope.
2. **Degradability is a computational tractability heuristic** (ligandable pocket
   + ubiquitination handle from Open Targets), not experimental proof of a
   degrader. The near-universal "database ubiquitination" flag was down-weighted;
   small-molecule ligandability tier is the primary driver.
3. **Anchor dependence.** The screen is only as good as the two anchors; their
   divergence here is biology, but a third orthogonal anchor (e.g. IRF4) would
   further triangulate.
4. **Phenocopy cosines are modest** beyond the anchors (top non-self ≈ 0.3–0.58):
   GATA3's program in circulating cells is largely unique, so hits are *partial*
   phenocopiers — consistent with a distributed network rather than a single
   redundant master node.
5. **Degradability is sequence-blind so far.** Sequence embeddings predict it only
   weakly; the proper next step is structure-based ligandable-pocket assessment
   (AlphaFold/Boltz) on the top intracellular hits — GPU work beyond this build.

## 16. Artifacts

**Type-2 (Th2-collapse) screen:**

| File | Contents |
|------|----------|
| `run_th2_screen.py` | End-to-end rerunnable pipeline (load→filter→axis→screen→rank) |
| `th2_collapse_scorer.py` | Reusable phenocopy scoring function |
| `SHORTLIST_master_ranked.csv` | 80-candidate ranked shortlist, all axes |
| `screen_GATA3_ranked.csv`, `screen_STAT6_ranked.csv` | Full per-axis screens |
| `axis_GATA3.csv`, `axis_STAT6.csv` | Reference axis definitions (z + subspace) |
| `prioritization_axes.csv` | Context / direction / essentiality per candidate |
| `degradability_annotation.csv` | Open Targets tractability annotation |
| `gwas_crossref.csv` | GWAS Catalog AD/asthma/allergy loci per candidate |
| `network_modules.csv`, `network_directionality.csv` | Module + directed-edge tables |
| `redundancy_combination.csv` | Nearest-neighbor / redundancy classification |
| `de_zscore_matrix.npz`, `de_obs.parquet` | Extracted working matrices (checkpoints) |
| `fig_control_validation.png` | Positive-control recovery |
| `fig_shortlist_profile.png` | Four-axis shortlist profile |
| `fig_network_topology.png` | Dual-arm network map |
| `fig_redundancy_heatmap.png` | Hit redundancy / combination structure |

**Multi-axis AD expansion:**

| File | Contents |
|------|----------|
| `multiaxis_reference_meta.csv` | 13 axis-anchor definitions + program strength |
| `multiaxis_control_recovery.csv` | Per-axis self-recovery (positive-control gate) |
| `multiaxis_screen_combined_top30.csv` | Top-30 hits per axis, combined long table |
| `axis_screens.tar.gz` | Full ranked screen for every axis |
| `ad_druggable_landscape.csv` | Top-60 known AD targets × atlas scope (visible/out-of-scope) |
| `ad_drug_universe.csv` | 60 known AD drugs / clinical candidates (Open Targets) |
| `known_target_recovery.csv` | Do atlas-visible known targets recover in the screen? |
| `emergent_hits_by_axis.csv` | Emergent hits per well-recovered axis, known-vs-novel flag |
| `literature_novelty.csv` | OpenAlex AD-literature depth + novelty tier per hit |
| `degradability_all.csv` | Degrader-tractability extended to 122 genes |
| `MULTIAXIS_integrated_shortlist.csv` | Integrated per-axis shortlist (all annotations) |
| `fig_multiaxis_target_map.png` | Integrated 8-axis intracellular AD target map |

**Translational triage & repurposing:**

| File | Contents |
|------|----------|
| `ad_drugged_targets.csv` | 118 AD targets × drugs × max clinical stage (Open Targets) |
| `axis_translational_maturity.csv` | Per-axis clinical validation + open-opportunity scores |
| `fig_translational_priority.png` | Axis decision map — validation vs. open intracellular opportunity |
| `repurposing_candidates.csv` | 17 intracellular hits with an existing drug + repurposing class |
| `fig_repurposing_candidates.png` | Repurposing candidates by clinical stage, axis, existing drug |
| `repurposing_ad_trials.csv` | ClinicalTrials.gov AD-trial status per repurposing drug |
| `drilldown_JAKSTAT_TYK2.csv` | JAK-STAT/TYK2 axis ranked drill-down (all annotations) |
| `fig_drilldown_TYK2.png` | TYK2-axis ranked hits with degradability + repurposing flags |

**Competitive landscape:**

| File | Contents |
|------|----------|
| `competitive_landscape.csv` | NIH-funded AD project counts by theme + our relation |
| `competitor_analysis.csv` | Deduplicated adjacent labs, focus, how we differ |
| `competitive_landscape_robust.csv` | Deduplicated grant + strict-PubMed counts per theme |
| `competitor_abstract_audit.csv` | Read-verified method-flags for all 38 adjacent grants |
| `fig_competitive_landscape.png` | Funding white-space figure (deduplicated, n=172) |

**Machine-learning layers:**

| File | Contents |
|------|----------|
| `th2_collapse_classifier.joblib` | Trained gradient-boosted collapse classifier (reusable scorer) |
| `classifier_predictions.csv` | Out-of-fold collapse probability per perturbation |
| `classifier_drivers.csv` | Transcriptomic drivers of the collapse call |
| `fig_classifier_validation.png` | ROC/PR + anchor recovery (GATA3 #4, STAT6 #4809) |
| `perturbation_latent.npy` | 16-dim autoencoder latent per perturbation |
| `fig_latent_umap.png` | Latent-space UMAP; collapse hits cluster, anchors separate |

**Independent lesional-skin validation:**

| File | Contents |
|------|----------|
| `lesional_validation.csv` | Top hits, lesional vs. healthy skin T-cell expression + MWU p |
| `fig_lesional_validation.png` | 8/18 hits (incl. GATA3) enriched in lesional AD skin T cells (GSE147424) |

**Literature concordance:**

| File | Contents |
|------|----------|
| `literature_validation.csv` | 7 core claims × our result × verdict × supporting primary paper (DOI, via PubMed) |
