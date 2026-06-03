# Learning-Outcomes Audit — scRNA-seq Workshop

**Prepared:** 2026-06-03 · **Scope:** the entire course as it will be shared with
the class — lectures (00–13), chapters/readings, laptop tutorials (ifnb),
Talapas advanced track + R scripts, appendices, glossary, and the top-level
site (index, setup, schedule, navbar).

**What this is:** a *prioritized* list of changes that would most improve student
learning outcomes, ranked **Must do** / **Should do** / **Might do**. This is a
review document only — nothing in the course has been changed. It is deliberately
oriented to *pedagogy and runnability* (can a student follow it, run it, and learn
the intended thing?), which is distinct from the recent spelling/cross-ref QA pass.

**How it was produced:** six parallel deep reads — one each for lectures,
chapters, laptop tutorials, the Talapas track, appendices + glossary, and the
top-level site — plus a few claims I verified by hand. Where a reviewer flagged
something I checked and found to be *fine*, I dropped it (noted at the end).

A "Must" is something that will **block a student, mislead them, or break the
1:1:1:1 promise** the course makes (lecture ↔ chapter ↔ laptop tutorial ↔ Talapas
script). A "Should" materially improves comprehension or trust. A "Might" is
polish.

---

## TL;DR — the highest-leverage fixes

If you only do five things before sharing:

1. **Fix the Talapas track so it runs end-to-end.** Three independent breakages
   (object directory split, a wrong column name in step 06, and missing `logs/`
   dirs) mean a student following the advanced track hits a hard stop they cannot
   debug. (M1–M3)
2. **Reconcile "teach Leiden / run Louvain."** Lecture 02 teaches Leiden
   (`algorithm = 4`); the tutorials run the default Louvain. Students will be
   confused about what they actually did. (M4)
3. **Fix the donor-metadata story in laptop Tutorials 06 & 08.** They still say
   ifnb "strips per-cell donor IDs" and fall back to synthetic pseudo-donors —
   but Tutorial 01 loads ifnb via `muscData::Kang18_8vs8()`, which *does* carry
   real donors. The pseudobulk lesson is taught on fake replicates when real ones
   exist. (M5)
4. **State the course's learning outcomes once, up front, and make the schedule
   realistic.** `index.qmd` is too thin and Day 3 is over-stuffed (Tut07 + Lec08 +
   Tut08 + closing in ~90 min). (M6, M7)
5. **Add expected outputs (cell/cluster counts, key numbers) to the laptop
   tutorials.** They render with `eval: false`, so students currently have *no*
   ground truth to check their run against. (M8)

---

## MUST DO

> Blocks a student, actively misleads, or breaks the 1:1:1:1 mapping.

### M1 — Talapas track: unify the object directory (`../data/` vs `../objects/`)
**Where:** `Tutorial_Talapas_03/04/05` (save/load `../data/`) vs
`Tutorial_Talapas_06–10` and **every** `scripts/*.R` (`OBJ_DIR`, default
`../objects/`).
**Problem:** Tut 03–05 write derived `.rds` to `../data/`; Tut 06+ and all scripts
read them from `../objects/`. A student following the chunks gets **file-not-found
at step 06** and the pipeline halts. This is the single highest-impact fix in the
whole course.
**Fix:** pick one location everywhere (recommend `../objects/` for derived objects,
reserving `../data/` for the downloaded matrix).

### M2 — Talapas track: broken manual-label column in reference annotation
**Where:** `scripts/06_reference_annotate.R` (and `Tutorial_Talapas_06`).
**Problem:** step 05 writes the manual label to `seu$celltype_manual`, but step 06
looks for `"celltype"` or `"seurat_annotations"` — neither exists — so
reconciliation indexes a non-existent column and errors / fills `celltype_final`
with NA. Step 06 fails for every student.
**Fix:** make the lookup prefer `celltype_manual`.

### M3 — Talapas track: jobs write to `logs/` but nothing creates it
**Where:** `scripts/run_rscript.sbatch`, `scripts/01_slurm_demo.sbatch`,
`Tutorial_Talapas_01` §6.
**Problem:** `--output=logs/...` with no prior `mkdir -p logs` → SLURM cannot open
the output file and the job **dies immediately, often with no `.out` to explain
why** — maddening for a first-time HPC user. (Tut 02 does it right; Tut 01 and the
wrapper don't.)
**Fix:** add `mkdir -p logs` before the first `sbatch`.

### M4 — "Teach Leiden, run Louvain" clustering mismatch
**Where:** Lecture 02 teaches `FindClusters(..., algorithm = 4)` (Leiden); the
laptop tutorials and `scripts/04_dimred_cluster.R` run `FindClusters()` with the
default (Louvain), and the script comment even *says* "Leiden."
**Problem:** students are taught one algorithm and run another, with no note on the
difference or on Leiden's extra `leidenalg`/`reticulate` dependency (a real HPC
gotcha). Undermines the core clustering lesson.
**Fix:** either run Leiden consistently (and document the dependency) or teach the
default Louvain and add one slide on why/when you'd switch to Leiden. Fix the
script comment either way.

### M5 — Donor metadata: ifnb *does* carry real donors (Tut 06 & 08 say otherwise)
**Where:** `Exercise_Folder/Tutorial_06_DESeq2_DE.qmd` (lines ~288/298) and
`Tutorial_08`, vs `Tutorial_01_QC_Preprocessing.qmd` (line 56).
**Problem:** Tut 01 loads ifnb via `muscData::Kang18_8vs8()`, which carries the real
donor `ind`. But Tut 06/08 still say the "SeuratData ifnb strips per-cell donor
IDs" and default to **synthetic pseudo-donors**. The entire pseudobulk lesson — the
"why biological replicates matter" punchline — is taught on fabricated replicates
when 8-vs-8 real donors are right there.
**Fix:** use the real `ind` donors for pseudobulk; remove the synthetic-donor
fallback framing (or demote it to a clearly-labeled "if you loaded ifnb the other
way" aside).

### M6 — No explicit, up-front statement of course learning outcomes
**Where:** `index.qmd` (too thin) and the top of the lecture sequence.
**Problem:** students arrive without a one-screen answer to "by the end of this
week, I will be able to ___." Orientation material is buried.
**Fix:** add a short, measurable "Learning outcomes" block to `index.qmd` (5–7
"you will be able to…" bullets) and move the orientation up.

### M7 — The schedule is unrealistic in two places
**Where:** Schedule / day plan.
**Problem:** Day 1 morning packs Welcome + Lec 00 + Lec 01; Day 3 working-lunch
crams **Tut 07 + Lec 08 + Tut 08 + closing into ~90 minutes**. Students fall behind
and the last day's material gets shortchanged.
**Fix:** rebalance — move at least one Day-3 item earlier, and give Lec 00/01
separate slots.

### M8 — Laptop tutorials give no expected outputs to check against
**Where:** all laptop tutorials (rendered with `eval: false`).
**Problem:** because code is not executed in the rendered site, students see **no
cell counts, no cluster counts, no numbers** — they have no ground truth to know
whether their own run worked.
**Fix:** add short "you should see roughly N cells / K clusters / this UMAP shape"
expectations at the key checkpoints (after QC, after clustering, after annotation).

### M9 — Chapter 12 (Trajectory) is a structural/depth outlier
**Where:** the Trajectory chapter (reading).
**Problem:** every other chapter has numbered sections, a worked example,
equations, and "Common errors" / "Going further" boxes; Ch 12 has almost none.
Trajectory is conceptually hard *and* the matching Tutorial 12 has the heaviest
dependencies — the thinnest reading sits under the hardest topic.
**Fix:** bring Ch 12 up to the structural template of its neighbors.

### M10 — Laptop Tutorial 12 (Trajectory) is thin and unscaffolded under the
heaviest dependencies
**Where:** `Exercise_Folder/Tutorial_12_*`.
**Problem:** heaviest install burden (Slingshot/Monocle3/scVelo/CytoTRACE family)
with the least scaffolding and no fallback. High risk students simply can't run it.
**Fix:** add a dependency/install preflight, a "if X won't install, do Y" fallback,
and the same Think-about-it scaffolding as earlier tutorials.

### M11 — Glossary is missing ~25 core analytical terms
**Where:** `Glossary`.
**Problem:** the glossary omits the terms students will actually look up — GEM,
ambient RNA, HVG, SNN, Leiden, Louvain, UMAP, Harmony, anchors, pseudobulk,
dispersion, shrinkage, FDR, SingleR, Azimuth, SCTransform, FAIR, and more.
**Fix:** add concise entries for the core vocabulary, cross-linked from first use.

### M12 — Appendix H (Talapas partitions) is stale and contradicts the tutorials
**Where:** Appendix H vs the SLURM tutorial and `Resources_Folder/VSCode_Talapas.qmd`.
**Problem:** Appendix H lists old partitions/hardware (`short`/`long`/`fat`,
K80/A100) that contradict the current docs (`compute`/`computelong`/`memory`/
`interactive`, A40–H100, `/scratch`). A student copying partition names from the
appendix submits jobs to queues that no longer exist.
**Fix:** align Appendix H with `VSCode_Talapas.qmd`'s partition table and storage
list (it's already correct there).

### M13 — `Software_Setup` verification step doesn't actually verify
**Where:** `Software_Setup` (and `index`/setup flow).
**Problem:** the "did it install?" check uses `requireNamespace(...)`, which returns
TRUE even when a package's *compiled* dependency is broken (it only checks the
package is findable, not loadable). Students get a green check, then fail mid-class.
**Fix:** verify with `library(pkg)` (which actually loads), and pre-warm the large
reference downloads (Azimuth/celldex) before the relevant session.

---

## SHOULD DO

> Materially improves comprehension or trust.

### Lectures
- **S1 — Add active-learning checkpoints.** Lectures 03/04/06/07/09/10/12/13 are
  near-continuous exposition. Add a "predict / discuss / try" prompt every ~10
  slides. (Lec 04 is especially thin and has no figures — give it a worked figure.)
- **S2 — Lecture 02 front-loads the math.** Move the intuition before the linear
  algebra; gate the derivations behind a "for the curious" aside.
- **S3 — Lecture 06 has nine goals.** Trim to 4–5 measurable outcomes; the rest
  become "we'll also touch on…".
- **S4 — Lecture 09 uses `algorithm = 3` with no explanation.** Either explain the
  choice or align with the clustering decision in M4.

### Chapters / readings
- **S5 — No chapter has a "Key takeaways" box.** Add a 3–5 bullet recap to each
  chapter; it's the single cheapest comprehension win across the readings.
- **S6 — Chapters 3 and 4 overlap.** Tighten the redundancy between
  markers/manual annotation and reference annotation so students see a clean
  hand-off, not a repeat.
- **S7 — Glossary linking is inconsistent.** Once M11 adds the terms, link first
  use consistently from chapters and tutorials.

### Laptop tutorials
- **S8 — Flag heavy downloads and give fallbacks.** Tut 04 (Azimuth/celldex
  references) and Tut 09 (~1.8 GB ATAC download) should warn about size up front
  and offer an Azimuth fallback / cached-data path.
- **S9 — Add time estimates to Tutorials 09–13.** Students need to know which
  steps are coffee-break long.
- **S10 — Tut 03 manual labeling has no self-check.** Add a "your cluster→label
  table should look roughly like this" checkpoint.
- **S11 — Add a data-flow / branch map.** Make explicit that Tut 04 and Tut 12 are
  side-branches off the main ifnb spine, so students don't think they broke the
  through-line.
- **S12 — Tut 13 has a known Seurat-v5 failure.** `SaveH5Seurat` fails on v5 and
  `glue` is an undeclared dependency — fix the save path and declare `glue`.
- **S13 — Make Think-about-it scaffolding consistent.** Some tutorials have rich
  prompts, others none; standardize.

### Talapas track
- **S14 — QC thresholds diverge between script and tutorial with no explanation.**
  `Tutorial_Talapas_03` filters `nFeature < 2500 & mt < 5`; `scripts/03` uses
  `nFeature < 6000 & nCount < 50000 & mt < 15`. The "mirror of the script" framing
  is broken when they disagree. Align them (the looser tumor-appropriate script
  thresholds are arguably the right canonical choice — and explain *why* tumor
  needs `mt < 15`).
- **S15 — `OUT_DIR` vs `OBJ_DIR` env-var naming.** Script 03 reads `OUT_DIR`;
  scripts 04–10 read `OBJ_DIR`. A student who exports `OBJ_DIR` to relocate to
  `/scratch/` finds script 03 ignores it. Standardize the variable name.
- **S16 — Step 06's annotated object is silently discarded.** `scripts/07` reads
  `nsclc_annotated.rds`, not the `..._final.rds` that step 06 produced — so the
  SingleR work never flows downstream. Defensible as a design, but it *looks*
  broken to a learner. Make the rationale loud or carry the `_final` object forward.
- **S17 — Random pseudo-donors → DE/DA results are pure noise.** Scripts 07/08/10
  assign donors/conditions by coin-flip, so pseudobulk DE and milo DA test random
  groups and (correctly) return ~no hits. Add an explicit callout: *"Expect
  near-zero significant results — the groups are random; this validates the
  plumbing, not biology,"* and point to CMO demultiplexing as the real prerequisite.
- **S18 — The track switches execution models midway with no bridge.** Tut 03–05
  read like "run interactively in RStudio"; Tut 06–10 assume SLURM submission via
  `run_rscript.sbatch` + `OBJ_DIR`. Add a bridging note that the canonical path is
  the R scripts via the wrapper and the `.qmd` chunks are the readable mirror.
- **S19 — `N_PCS` differs (script 30 vs tutorial 15)** in dimred/clustering — will
  change results between the two "equivalent" paths. Align.
- **S20 — Pin R version consistently and fail fast.** `Tutorial_Talapas_01` §6.3
  uses `R/4.3.2`; `run_rscript.sbatch` uses `R/4.4.0`; the wrapper's
  `module load R/4.4.0 || true` *swallows* a failed load. Pick one verified module
  string, use it everywhere, and drop `|| true` so failures surface at the real
  point.
- **S21 — Cell Ranger onboarding fallback.** Add `module spider cellranger` as the
  first step and warn that the signed download URL must be quoted and expires in
  minutes (the most common real Cell Ranger onboarding failure).
- **S22 — Fix the Tut 10 / Tut 09 code mirrors** so copy-pasting the rendered chunk
  works: Tut 10's milo mirror omits the `condition` synthesis (`condition not
  found`); Tut 09's GSEA mirror omits the gene-ID dedupe (`gseGO` errors on
  duplicates).

### Top-level site
- **S23 — Reconcile conflicting project-directory conventions** across
  `Software_Setup`, `Datasets`, and `FAQ` (they tell students to put data in
  different places).
- **S24 — Required-vs-bonus dataset checklist.** `stxBrain` is in the required
  install list but is bonus material; mark required vs optional so students don't
  burn disk/time on extras.
- **S25 — Add an end-to-end smoke test** students can run after setup to confirm
  the stack works before Day 1.
- **S26 — Surface useful pages in the navbar.** Datasets, FAQ,
  Compute_Environments, Accessibility, and Code_of_Conduct exist but aren't linked
  from the navbar.
- **S27 — Define the "chat channel"** (it's referenced but never named), and fix
  the circular recording-policy cross-references.
- **S28 — Standardize "1:1:1" vs "1:1:1:1" wording.** Now that the Talapas scripts
  are a fourth rail, the mapping is 1:1:1:1 — make the language consistent
  everywhere.

### Appendices / glossary
- **S29 — Appendix D: add the integration/abundance statistics** students meet in
  the tutorials (SpatialFDR, LISI/kBET) so the stats appendix matches the pipeline.
- **S30 — Link the orphaned appendices.** Appendices B/E/F/G/H/J aren't linked from
  the places students need them; add pointers from the relevant tutorials/lectures.
- **S31 — De-duplicate.** Appendix C §11's mini-glossary duplicates the main
  Glossary; Appendices B/F repeat FASTA examples and carry no scRNA content;
  Appendices I/J overlap on cell×gene vs gene×cell. Consolidate.

---

## MIGHT DO

> Polish; do if time permits.

- **G1 — Lecture 00:** open with a biological hook and make the goals measurable.
- **G2 — Remove leftover "(dup1/dup2/dup3)" image placeholder slides** (notably in
  Lec 00).
- **G3 — Add a single laptop↔Talapas mapping table** (Talapas 03 = laptop 01 … the
  "+2 offset"); it's stated per-tutorial but never collected in one place.
- **G4 — Talapas right-sizing:** Tut 09's submit line over-requests CPUs (wrapper
  default 8, step needs 4); show the full `--cpus-per-task` override consistently.
- **G5 — Modernize Seurat-v5 idioms:** `AggregateExpression(..., slot=)` →
  `layer=` in `scripts/08` (codebase is mid-migration).
- **G6 — Per-step "did it work?" checklist** on the Talapas track: expected `.out`
  sentinel line, expected output file, and `seff <jobid>` → `State: COMPLETED` (Tut
  01 teaches `seff`/`sacct` well — just connect it to the analysis steps).
- **G7 — Concrete footprint guidance:** recommend running the Talapas track under
  `/projects/<PIRG>/<duckid>/...` or `/scratch/` and note the ~50–80 GB total
  (refs + FASTQs + objects) so students don't blow a `/home` quota.
- **G8 — Accessibility page** makes an absolute claim that conflicts with the
  known-gaps list; soften to "working toward."
- **G9 — Disk budget mismatch** (setup says 10 GB; realistic total ~12 GB+).
- **G10 — `.qmd.txt` download rename stumble** — smooth the instruction so students
  know to strip the `.txt`.
- **G11 — Module version drift** in appendices (R/4.3.3 vs 4.4.1) — align once a
  canonical version is chosen (ties to S20).
- **G12 — Appendix E** OS-internals goes deeper than the audience needs; trim or
  mark as optional.
- **G13 — No `renv.lock` exists** though `Compute_Environments`/`Software_Setup`
  reference one. Either add a real lockfile (best reproducibility win) or stop
  referencing it. *(Listed here as low-urgency, but it's the cleanest path to true
  reproducibility if you ever want it.)*

---

## Reviewer claims I checked and dropped (no action needed)

- **Tutorial 09 ATAC link "mismatch."** A reviewer flagged the landing-page link
  vs the `wget` URL as hg19/hg38 inconsistent. I verified both (line 52 and line
  70): both are the v1/hg19 `atac_v1_pbmc_10k` dataset. **Consistent — false
  positive.**

---

*End of audit. Nothing in the course has been modified; this is a planning
document. Recommend tackling the MUST items (especially M1–M5) before sharing the
site with the class, since those either break a student's run or teach the wrong
thing.*
