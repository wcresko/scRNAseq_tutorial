# Comparison: Hope Healey's two scRNA-seq lectures vs. our workshop lectures

Review note (not part of the rendered site). Compares Hope's `Hope_scRNAseq_lecture_1`
(sequencing technology + Cell Ranger processing) and `Hope_scRNAseq_lecture_2`
(QC → analysis → advanced methods) against our `Lecture_Folder/` deck.

Three things were requested: (1) find content of Hope's that we're missing,
(2) incorporate her single-cell-platform descriptions (10x vs Parse) into the
appropriate initial lecture, (3) where Hope's images are clearer, add duplicate
slides so the old/new can be compared.

---

## Task 2 — platform descriptions: DONE (Lecture 00)

Added a platform deep-dive to `Lecture_00_Bulk_vs_scRNAseq.qmd` (the tech-overview
lecture), written in our own words / house style:

- **Two ways into the Chromium: polyA vs. probe** — 3′/5′ polyA GEX (any species,
  live) vs. Flex/probe capture (fixed incl. FFPE, human & mouse only, microarray-like).
- **A polyA-capture catch: 3′ UTR annotation** — see gap #1 below (this was also a
  genuine content gap, not just platform framing).
- **Cells vs. nuclei, and multi-omic readouts** — cytoplasmic vs nuclear RNA;
  Multiome (RNA+ATAC), CITE-seq protein, 5′ V(D)J, CRISPR/Perturb-seq; sample
  multiplexing (CellPlex/hashing/genotype demux).
- **Beyond 10x: combinatorial-barcoding platforms** — SPLiT-seq → Parse Evercode
  (split-pool, no instrument, fixed cells/nuclei, scales to millions of cells / many
  samples, more bench labor, no ATAC); Fluent PIP-seq mentioned.
- **10x Chromium vs. Parse Evercode — at a glance** — head-to-head comparison table.
- **Platform choice changes results, not just cost** — cost/replicate/depth design
  guidance + the Filippov et al. 2024 point (described, not figure-copied) that
  platform choice can change which populations you recover.
- Enriched the **timeline** bullets with named milestones (Tang 2009 → Streets 2014 →
  Drop-seq 2015 → Chromium 2015 → SPLiT-seq/Parse 2018 → PIP-seq 2023).

---

## Task 1 — content Hope has that we did not

### Incorporated now (committed to Lecture 00)
1. **3′ UTR annotation problem for polyA capture** *(Hope L1 s25–26; she cites Healey
   et al. 2022).* 3′ capture concentrates reads at the 3′ UTR, so quantification
   depends on good 3′ UTR annotation; missing/truncated UTRs make genes vanish — a
   **major issue in non-model systems** and relevant to this lab specifically. Added
   as its own slide. **High value, license-clean (written in our words).**
2. **Platform-choice-affects-biology** (Filippov 2024) and **cost/replicate design**
   numbers — added (see Task 2).
3. **SPLiT-seq / PIP-seq / Parse milestones** — added to timeline.

### Recommended additions (NOT yet written — your call on scope)
Ordered by my sense of value for this audience.

- **A. Trajectory / pseudotime / RNA velocity** *(Hope L2 s59–61: Velocyto/scVelo,
  CytoTRACE, Monocle).* We have **no lecture** on this — it's an entire analysis
  category we only mention in passing (UMAP caveats in Lec 02). Biggest single
  topical gap. Options: a short section in an existing "downstream" lecture, or a new
  bonus lecture paralleling WGCNA/Spatial.
- **B. Cell–cell communication / ligand–receptor** *(Hope L2 s62–63: CellChat,
  CellPhoneDB).* Also absent from our lectures. Good candidate for a short
  "downstream analyses" slide or bonus section.
- **C. More even-handed DE framing.** Our Lec 06 strongly defaults to pseudobulk and
  frames per-cell tests as "wrong" for condition contrasts. Hope presents **both
  camps** — pseudobulk (edgeR/limma/DESeq2) *and* cell-level mixed-model methods
  (**MAST, nebula, glmmTMB**), noting the literature disagrees (Gagnon 2022 etc.). I'd
  keep our pseudobulk-first stance but add one slide acknowledging the model-based
  cell-level methods and when they're defensible. (Low effort, improves balance.)
- **D. Quantitative integration metrics** *(Hope L2 s28: iLISI/LISI; kBET).* Our Lec 05
  teaches over/under-correction diagnostics qualitatively (ISG-score gap) but not a
  named mixing metric. A short "how do we *measure* integration success" slide (iLISI/
  kBET) would round it out.
- **E. Cell Ranger internals** *(Hope L1 s54–69).* We name Cell Ranger/STARsolo/alevin
  and teach the knee plot + EmptyDrops, but not the algorithmic steps (TSO/polyA
  trimming, STAR + MAPQ 255 logic, intron retention, CB correction by Hamming distance,
  UMI collapsing, read assignment, then **OrdMag + EmptyDrops** two-step cell calling).
  Optional "what Cell Ranger does under the hood" reference slide(s) — useful for the
  Talapas/raw-data track.
- **F. LLM/AI chatbots for annotation, with caution** *(Hope L2 s44–45).* Topical; a
  one-slide caution in Lec 04 (reference annotation) noting LLMs can suggest plausible-
  but-wrong markers. Low effort.
- **G. scWGCNA / pseudocells for networks** *(Hope L2 s55, s73).* Check whether Lec 10
  (WGCNA) already mentions scWGCNA + pseudocell aggregation; if not, add a bullet.
- **H. Long-read single-cell / isoforms (scISO-seq)** *(Hope L2 s68–69).* Advanced;
  ties to the 3′-annotation point. Optional bonus mention.
- **I. Experimental-design "spike-in" doublet estimation** *(Hope L2 s9: mix
  species/genotypes to measure multiplet rate).* One bullet in Lec 01's doublet slide.
- **J. Seurat object anatomy diagram** *(Hope L2 s29).* We show code but not an
  object-structure picture; minor.

### Things WE cover that Hope's lectures do not (coverage is broadly strong)
FAIR / data sharing (Lec 12); a full pseudobulk **DESeq2** walkthrough; **differential
abundance** with miloR as its own lecture; deeper **reference annotation** (SingleR +
Azimuth + CellTypist reconciliation, confidence scores); HVG mean–variance detail;
barcode-rank/knee plot; LFC shrinkage nuance; statistical-foundations appendix.

---

## Task 3 — image upgrades (shortlist + an IP note)

**Why I did not bulk-copy Hope's images into the repo:** almost all of her figures are
third-party copyrighted material — 10x Genomics diagrams and journal figures (Yu 2020,
Filippov 2024, Hafemeister & Satija 2019, Aran 2019, Manno 2018, Gulati 2020, Wolock
2019, Young & Behjati 2020, etc.). Committing/pushing them into the workshop repo (which
renders to a public site) would be redistributing copyrighted figures. So instead of
embedding ~40 such images, here is the shortlist of where Hope's figure is plausibly
clearer than ours, with her cited source. For any you want to keep, the clean paths are
(a) recreate it as an original house-style SVG (what our decks already do — I'm happy to
do the top picks), or (b) keep her figure with permission/attribution.

Candidates (Hope slide → our slide → source she cites):

| Hope figure | Our corresponding slide | Source (per Hope) | Why possibly clearer |
|---|---|---|---|
| GEM formation / chip / beads-oligos / barcoded library (L1 s17–23) | Lec 00 "Droplet generation", "Chromium library prep", "barcodes + UMIs" | 10x Genomics | Official 10x step-by-step diagrams; very legible |
| 10x molecule structure, 3′ kit (L1 s23) | Lec 00 "barcodes + UMIs" | 10x Genomics | Shows read1/read2/UMI/insert layout concretely |
| 3′ UTR alignment cartoon (L1 s25–26) | NEW 3′-UTR slide (Lec 00) | Healey et al. 2022 | We added the concept as text; a figure would help |
| Cell Ranger pipeline overview (L1 s50–57) | Lec 00 "Raw data processing" | 10x Genomics | Linear step diagram maps to the prose |
| Dropout / capture cartoon (L1 s28–36) | Lec 00 "Sparsity" simulated panel | Yu 2020 | Conceptual cartoon vs our simulated heatmap |
| NormalizeData vs SCTransform (L2 s17–19) | Lec 01 normalization slides | Hafemeister & Satija 2019 | The paper's before/after UMAPs are crisp |
| CCA integration schematic (L2 s25) | Lec 05 integration | Stuart & Butler 2019 | Clean anchor/MNN schematic |
| Harmony schematic + iLISI (L2 s26, s28) | Lec 05 integration | Korsunsky 2019 | If we add integration metrics |
| PCs 30 vs 40 → extra cluster (L2 s35) | Lec 02 "choosing resolution/PCs" | Filippov 2024 | Concrete "more PCs found a real cell type" example |
| RNA velocity / CytoTRACE (L2 s60–61) | (no slide yet) | Manno 2018; Gulati 2020 | Only if we add trajectory (gap A) |
| CellChat circle plot (L2 s63) | (no slide yet) | Efremova 2020 | Only if we add cell–cell comm (gap B) |

If you'd like, point me at the 3–5 you care about most and I'll recreate them as
original SVGs in our style (license-clean), or wire in your figures with attribution.
