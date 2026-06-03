# SeqDef project — context

## What this repo is
Analysis/manuscript repo for **SeqDef**, a phylogenetically weighted genomic-prioritization
statistic and R package. Companion R-package source lives in the sibling repo `SeqDef/`.

Manuscript: *"SeqDef: An R package for phylogenetically weighted genomic prioritization
across the tree of life"* — submitted to the **Journal of Evolutionary Biology**
(ID **JEB-2026-00140**), currently in **Major Revision**.

## Layout
- `function.R` — local copy of the `SeqDef()` and `calc_priority()` functions.
- `analysis.R` — Chondrichthyes pipeline: IUCN + NCBI + 100 posterior trees → SeqDef across
  trees → median scores → `calc_priority` (SeqDef × 2^GE) → circular ggtree (Fig 4).
  Also contains the commented toy-tree block used for Fig 1.
- `validation.R` — sensitivity (variance-vs-λ curves, Fig 2A/B) and topological-robustness
  (top target across 100 trees, Fig 3).
- `data/` — `chondrichthyes.nex` (100 posterior trees, VertLife/Stein et al. 2018),
  `iucn_assessment_data.rds`, `ncbi_assembly_data.rds`.
- `figures/` — fig1–4 PDFs.
- `seqdef.csv` — toy 10-taxon availability example.
- `revision/` — revision deliverables (see below).

## Key facts about the code (for the revision)
- **S** is any numeric availability score on [0,1]; `SeqDef()` does not require binary input.
  The case study uses binary (NCBI assembly presence/absence) by choice.
- Only the **exponential kernel** is implemented: `exp(-lambda * d / tree_depth)`.
- `lambda` options: `"auto_max"` (variance max, scans `seq(1,50,0.1)`, 10% drop tolerance),
  `"by_genus"` (half-life = median intra-genus distance), or numeric.
- Core cost is **O(n²)** (cophenetic matrix + weighted matrix–vector product); `auto_max`
  multiplies that by the λ-grid length.
- Two issues to fix before any public archive: **real API keys are present in `analysis.R`**
  (NCBI + IUCN, also in git history) and **`SeqDef/DESCRIPTION` has placeholder author info**.

## Revision work (this session — 2026-06-02)
Produced two deliverables in `revision/` responding to the JEB reviewer comments:
- `SeqDef_revision_strategy.docx` — per-point strategy + triage table + new-analyses list +
  pre-submission checklist + 8-week timeline.
- `SeqDef_response_to_reviewers.docx` — draft point-by-point rebuttal (R1 1.1–1.5,
  R2 2.1–2.6, editor) with `[bracketed]` placeholders for results/line numbers.

Reviewer-suggested references identified: **Soares et al. 2023, *Biological Conservation***
(line-40 cite; reviewer said "Conservation Biology" — wrong journal) and **EDGE2 = Gumbs
et al. 2023, *PLoS Biology***.

### Still to do (planned, not yet run)
> Full build spec for the coding agent: **`REVISION_ANALYSES.md`** (analyses 1–6, outputs, acceptance criteria).

1. λ rank-stability + auto_max-vs-by_genus head-to-head.
2. Kernel comparison (add `kernel=` arg: exponential/Gaussian/linear).
3. EDGE / EDGE2 benchmark on the Chondrichthyes data.
4. Runtime + memory benchmarks (n = 100–10,000) + complexity statement.
5. Regenerate Fig 1 from a real `SeqDef()` run with λ annotated; supplementary continuous-S example.
6. Notation pass (Eq 4 roman type; **restore missing λ in the kernel equation**).
7. Data Availability Statement + CRediT; rotate/scrub API keys; fill DESCRIPTION/citation; add tests.
