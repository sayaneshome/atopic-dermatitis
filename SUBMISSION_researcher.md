# Researcher-track submission — Intracellular master-regulator map of atopic dermatitis

## 1. Written summary (100–200 words)

**Every approved atopic dermatitis (AD) drug blocks the type-2 pathway from the
outside — injectable antibodies against extracellular cytokines. The field is
moving inward, to the intracellular transcription factors that run the allergic
program, because degrading a master node promises deeper, oral suppression.
STAT6 is the proof of concept. But no one had systematically asked which *other*
intracellular regulators would collapse the program if shut down.**

We screened a genome-scale CD4⁺ T-cell Perturb-seq atlas (~34,000 perturbation
signatures), scoring every gene knockout for how strongly it phenocopies the
STAT6 and GATA3 "Th2-collapse" signature. The pipeline recovers both anchors at
rank 1 — a built-in positive control — and surfaced a finding we didn't expect:
**in circulating cells, STAT6 and GATA3 knockouts anti-correlate**, meaning
GATA3 maintains the type-2 program STAT6-independently. We then expanded to every
AD immune axis a T cell can represent (JAK-STAT, Th22/AHR, Th1, activation),
cross-referenced the known AD drug universe, added literature and
degrader-tractability layers, and ranked repurposing candidates against
ClinicalTrials.gov. Deliverables: a druggability-annotated target map, a reusable
Th2-collapse classifier, and a fully rerunnable pipeline.

*(word count ≈ 190)*

**Track fit — Build From the Bench.** We started from a bench question (which
intracellular regulators beyond STAT6 collapse the allergic program if degraded?),
found the dataset built to answer it — the Marson lab's genome-scale CD4⁺ T-cell
Perturb-seq atlas (Zhu, Dann, Yan et al., bioRxiv 2025, doi:10.64898/2025.12.23.696273;
one of the track's suggested starting points) — and used Claude Science end-to-end to
deliver a discrete, reproducible finding plus a target map others can build on.

---

## 2. Three-minute demo video script (storyboard + voiceover)

**Total: 180 s. Screen-record the artifacts; talk over them. Do not show code
scrolling — show results.**

### [0:00–0:25] The gap (hook)
*On screen:* a slide — a cartoon cell with antibodies bouncing off the outside,
then an arrow pointing inward to "STAT6". *Optional cutaway to*
`fig_competitive_landscape.png` *on the word "no one".*
> "Every approved atopic dermatitis drug works from outside the cell. The frontier
> is inside — the transcription factors running the allergic program. STAT6 just
> became the first, degraded orally. But across 172 NIH-funded AD projects, no one
> is screening these intracellular nodes at genome scale, and no one is pursuing
> degraders. Our question: which *other* intracellular nodes would collapse this
> disease if you shut them down — and which are druggable?"

### [0:25–0:55] The method (one clear idea)
*On screen:* `fig_control_validation.png` — the two rank-vs-cosine panels.
> "We took a genome-scale Perturb-seq atlas — 34,000 gene-knockout signatures in
> CD4 T cells — and scored every knockout for how much it phenocopies the STAT6
> and GATA3 collapse signature. The validation: both anchors come back at rank 1.
> The pipeline re-derives its own positive controls."

### [0:55–1:30] The finding (the thing judges remember)
*On screen:* the canonical-gene table (IL13/IL4/IFNG effects) → then
`fig_network_topology.png`.
> "Then the surprise. STAT6 and GATA3 knockouts *anti-correlate*. Knocking out
> GATA3 collapses the type-2 cytokines. Knocking out STAT6 doesn't — it actually
> de-represses them. In circulating T cells, GATA3 maintains the program
> independently of STAT6. That's a real, testable biological claim the atlas
> handed us — and it reframes what 'beyond STAT6' even means."

### [1:20–1:30] Two more methods agree (rigor)
*On screen:* `fig_classifier_validation.png` (ROC 0.94) → `fig_latent_umap.png`.
> "And it's not one method. A supervised classifier — AUC 0.94, trained with the
> label genes removed — ranks GATA3 fourth and STAT6 four-thousandth. An
> autoencoder's latent space puts GATA3 next to its top hit and STAT6 across the
> map. Three independent models, same answer. And it's not just internally
> consistent — that STAT6-independent GATA3 maintenance loop is exactly what the
> classical immunology literature reported twenty years ago. The screen
> rediscovered textbook biology it was never told."

### [1:30–2:20] The scale-up (depth past first idea)
*On screen:* `fig_multiaxis_target_map.png` scrolling across the 8 axes → then
`fig_translational_priority.png`.
> "We didn't stop at one axis. We rebuilt the screen for every AD immune axis a T
> cell can represent — JAK-STAT, Th22, Th1, activation — recovering textbook
> biology each time: knock out TYK2, you get its receptor IFNAR1; knock out AHR,
> you get its dimer partner ARNT. Then we asked which axis is most translational —
> proven pathway *and* open intracellular nodes. JAK-STAT/TYK2, Th22/AHR, and the
> original type-2 axis win."

### [2:20–2:55] Translation (impact + Claude use)
*On screen:* `fig_repurposing_candidates.png` → highlight ITK/ritlecitinib row.
> "Finally, repurposing. Cross-referencing every hit against approved drugs and
> ClinicalTrials.gov, two oral drugs jump out — ritlecitinib and tazemetostat —
> approved elsewhere, hitting nodes our screen nominates, and *never tested in
> atopic dermatitis*. Claude drove all of this: querying Open Targets, GWAS
> Catalog, OpenAlex, and ClinicalTrials.gov, building every figure, and writing a
> pipeline anyone can rerun on the released data."

### [2:55–3:00] Close
*On screen:* the title slide + repo URL.
> "A data-driven, druggable map of the intracellular allergic program. All
> reproducible. Thank you."

**Recording tips:** lead with the anti-correlation finding if you have to cut for
time — it is the single most memorable, defensible result. Keep every on-screen
asset one of your saved figures (judges trust what they can see). Narrate
findings, never tool mechanics.

---

## 3. Repository / write-up structure

Organize the submitted repo so a judge can verify the finding in under five
minutes:

```
/README.md              ← the written summary + "run this" quickstart
/METHODS_and_RESULTS.md ← full 12-section methods (already written)
/pipeline/
   run_th2_screen.py       ← end-to-end, rerunnable from the 16 GB h5ad
   th2_collapse_scorer.py  ← the reusable classifier (build_axis / score)
/results/
   fig_control_validation.png      ← START HERE: proves the pipeline works
   fig_network_topology.png        ← the STAT6≠GATA3 finding
   fig_multiaxis_target_map.png    ← the 8-axis expansion
   fig_translational_priority.png  ← which axis to pursue
   fig_repurposing_candidates.png  ← the drug leads
   *.csv                            ← every ranked table
/data/
   README.md  ← link to the released AnnData object (do not commit 16 GB)
```

**README quickstart** must show the one command that reproduces the positive
control (`python pipeline/run_th2_screen.py`) and prints
`[GATA3] self rank 1` / `[STAT6] self rank 1`. A judge who sees that line trusts
everything downstream.

---

## 4. How this maps to the judging criteria

| Criterion | Weight | Your strongest evidence |
|-----------|--------|--------------------------|
| **Impact** | 25% | A finding others can build on: STAT6≠GATA3 divergence reframes "beyond STAT6"; a ranked, druggable target map + two untested oral repurposing leads (ritlecitinib, tazemetostat). Fits the problem statement (intracellular AD targets) exactly. |
| **Claude Use** | 25% | Claude orchestrated a full research pipeline: 4 scientific databases (Open Targets, GWAS Catalog, OpenAlex, ClinicalTrials.gov), a 16 GB Perturb-seq matrix, a custom phenocopy classifier, and every figure — end to end, not a single API call. |
| **Depth & Execution** | 20% | Pushed well past the first idea: single axis → 8 axes → translational triage → repurposing → trial check. Three independent methods (cosine screen, supervised classifier AUC 0.94, autoencoder latent) converge on the same STAT6≠GATA3 finding. Tested an ESM-2 sequence layer and *removed* it when it underperformed (AUC 0.61) — rigor over padding. Caught and corrected its own errors (anchor divergence, GWAS trait mislabeling, grant deduplication). Validated the finding against the primary literature — 6/7 core claims directly supported, the pivotal STAT6≠GATA3 result matching classical 2000–2004 immunology (Ouyang, Yamashita) the pipeline was never given. |
| **Demo** | 30% | Findings you trust: the pipeline re-derives its own controls; every claim has a figure; the anti-correlation result is visually unmistakable and biologically testable. |

**The one-sentence pitch for judges:** *"We turned a genome-scale knockout atlas
into a druggable map of the intracellular allergic program — and the screen
re-derives its own positive controls, so you can trust what it found."*

---

## 5. Competitive positioning (for the Impact score)

We checked the field before claiming novelty — and did it robustly. NIH
RePORTER returns duplicate records per fiscal year, so we deduplicated by core
project number (450 raw records → **172 unique AD grants**, FY2022–2025) and
classified each abstract by theme. The money goes to barrier/filaggrin (88),
microbiome (51), biologics (48), itch (34), and Th2 differentiation (32).
**Genome-scale perturbation screening for AD target discovery: 1 grant** — and
it's a skin-*microbiome* study, not target discovery. **Targeted protein
degradation in AD: 0 grants.** PubMed (MeSH-anchored, title/abstract) confirms:
**0 AD papers** on Perturb-seq/CRISPR-screen, **1** on PROTAC/molecular-glue,
vs. 354 on JAK inhibitors. The nearest adjacent lab (Kalung Cheung, Mount Sinai —
"Transcriptional Regulation of Th2 Cell Development") studies TFs one at a time;
we screen all of them at once and rank for degradability.

We didn't stop at keyword counts. We **read the full abstracts of all 38
deduplicated grants** in the Th2 / transcription-factor / STAT6-GATA3 cluster and
classified each for whether it runs a genome-scale screen, pursues a degrader, or
systematically ranks regulators to find new targets. **All 38: no on all three.**
The empty intersection this project occupies is verified by reading, not just by
text-matching.

**Judge-facing line:** *"This isn't a crowded space we're joining — it's an empty
one we're mapping. Zero funded NIH grants combine genome-scale T-cell screening
with degradability for AD; those are our two defining moves."*

*Honesty caveat to state:* even deduplicated and MeSH-anchored, these are
text-matched counts (a grant could use a method without naming it), so they're a
strong directional signal, not a literal census. We dropped OpenAlex for specific
numbers (too lexically noisy) and rely on deduplicated RePORTER grants + audited
PubMed queries. `fig_competitive_landscape.png` is the supporting slide.

### 5a. Prior art vs. our gap (preempt the "has this been done?" question)

Be precise, not sweeping — a judge knows the field, so claim the exact gap:

- **The method building blocks are established, and we rest on them.**
  Genome-scale Perturb-seq for genotype→phenotype mapping is proven (Replogle-style
  atlases); signature-based phenocopy screening has a decade of precedent
  (Connectivity Map / LINCS L1000). Our pipeline stands on validated foundations —
  that's a credibility strength, not a novelty debt.
- **The target premise is real and current.** Oral STAT6 degradation for AD is in
  first-in-human trials (KT-621, Phase 1a, 2025); a 2025 review covers intracellular
  small-molecule approaches in AD. The "move inward" thesis is happening now.
- **What no one has done — our contribution.** No prior work applies this machinery
  systematically to the STAT6/GATA3 collapse signature in a CD4⁺ T-cell atlas,
  ranked for degradability, for AD. Verified three ways: 0 / 172 deduplicated NIH
  grants, 0 AD Perturb-seq papers (PubMed), all 38 adjacent abstracts read.

**Judge-facing line:** *"The tools are proven — genome-scale Perturb-seq and
Connectivity-Map-style phenocopy screening — and oral STAT6 degraders are already
in trials. What's new is applying that machinery to ask which* other *intracellular
nodes collapse the allergic program and which are degradable. We're first to occupy
that intersection, and we verified it's empty."*

*(Sources verified by name, not by count — the KT-621, Connectivity Map, and
Perturb-seq precedents are real; we did not rely on noisy keyword totals.)*

### 5b. Concordance with prior literature (proves it's not an artifact)

A skeptical judge will ask: is the surprising finding real, or a pipeline
artifact? We checked the core claims against the primary immunology literature
(PubMed). **6 of 7 are directly supported; the 7th has clear precedent.** The
decisive one: our pivotal STAT6≠GATA3 divergence reproduces classical genetics
almost exactly — Ouyang et al. (*Immunity*, 2000) showed **Stat6-independent
GATA-3 autoactivation** directs Th2 commitment, and Yamashita et al. (*J Biol
Chem*, 2004) showed GATA3 is essential for the **maintenance** of Th2 cytokines
in differentiated cells. The chromatin-upstream-of-GATA3 result matches a 2023
targeted TF-CRISPR screen (Linzer et al., *Immunity*), which also validates the
screening approach. The pipeline rediscovered textbook biology it was never
told — the strongest evidence it reads real signal (`literature_validation.csv`;
all citations via PubMed).

**Judge-facing line:** *"The most surprising result isn't a fluke — it's a
twenty-year-old immunology finding the screen rederived from scratch, with the
answer genes hidden from it."*

## 6. Honest scope note (include it — it builds trust)

State plainly: the primary screen profiles **circulating** CD4⁺ T cells, not
skin-resident pathogenic Th2 cells. **We then executed the named validation
step** — testing the top hits in an independent lesional-skin scRNA-seq dataset
(GSE147424, 4 AD lesional + 4 healthy donors). **8 of 18 hits, including the
anchor GATA3 (and BATF, IL4R, IFNAR1), are significantly enriched in lesional AD
skin T cells vs. healthy** (`fig_lesional_validation.png`), and the same GATA3 >
STAT6 ordering the circulating screen found reappears in skin. The validation
discriminates — several hits do *not* transfer — so this is real cross-tissue
signal, not rubber-stamping. Remaining honest limits: it's a pseudobulk
comparison over one dataset, not a perturbation experiment in skin Th2 cells;
Treg/cytokine and barrier/itch axes stay out of scope. Judges reward a team that
tests its own caveat over one that only names it.

## 7. Conclusion (the pitch, in three sentences)

Atopic dermatitis is treated entirely from outside the cell; the frontier is the
intracellular network running the allergic program, and STAT6 is the first proof
that those nodes can be drugged orally. We turned a genome-scale CD4⁺ T-cell
Perturb-seq atlas into a druggable map of that network — a pipeline that
re-derives its own positive controls (STAT6 and GATA3 at rank 1), surfaced a
genuine finding (in circulating cells STAT6 and GATA3 knockouts anti-correlate,
so GATA3 maintains the type-2 program STAT6-independently), scaled to every AD
immune axis a T cell can represent, and layered on degradability, human genetics,
literature, and repurposing — pointing to two approved oral drugs (ritlecitinib,
tazemetostat) that hit atlas-nominated nodes and have never been tried in AD. And
we verified, by reading all 38 adjacent NIH grants, that no funded work occupies
this intersection: genome-scale screening plus degradability for AD is empty
space, and this is a reproducible, database-grounded map of it.
