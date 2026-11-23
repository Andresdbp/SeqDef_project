# SeqDef project — context

## Current state (2026-11-23): JEB second revision (R2) prepared

- **Paper:** "SeqDef: An R package for phylogenetically weighted prioritization of genome-sequencing
  targets" (title narrowed in R2; was "… genomic prioritization across the tree of life").
  Journal of Evolutionary Biology, JEB-2026-00140.R1 → major revision (decision 2026-11-05;
  resubmission due within 8 weeks, i.e. by ~2026-12-31).
- **Repositories (public):** package github.com/Andresdbp/SeqDef (v0.2.0); analyses
  github.com/Andresdbp/SeqDef_project (this repo). Archiving: GitHub only for now; a DOI snapshot can
  be added if the editor asks.
- **Manuscript master:** Google Doc "manuscript". The R2 edits were applied to a .docx export as
  tracked changes, kept outside this public repo in `~/Documents/Lab/seqdef/R2_submission/`:
  `R2_manuscript_tracked.docx` (tracked) and `R2_manuscript_clean.docx` / `.pdf` (all changes
  accepted; response-letter line numbers refer to this PDF). Re-import into Google Docs (File > Open > Upload) to keep working there; the
  tracked changes become suggestions.
- **Response letter:** `revision/R2_response_to_reviewers.md` / `.docx`.
- **Still for the authors:** CRediT statement and Data Availability Statement on the title page;
  AI-use disclosure required by JEB (cover letter + Acknowledgements); verify line numbers after the
  final file is built; rotate the NCBI and IUCN API keys that were exposed in the old git history.

## Key R2 results (all 100 posterior trees; median [2.5%, 97.5%])

- SeqDef assigns one score to all species of a clade without sequenced members (822 unsequenced
  species in 178 [168, 189] clades; within-clade range ≤ 2e-8).
- vs nearest-sequenced-relative distance: ρ = 0.75 [0.63, 0.83] (clade level 0.64); vs family lacks
  assembly: ρ = 0.43; top-1 target first-in-family in 99/100 trees; top-10/25/50: 90/40/24 %;
  NSD is the large-λ limit (ρ 0.97 at λ = 10, 0.995 at 25).
- Priority composition: ρ(Priority, GE) = 0.86, ρ(Priority, SeqDef) = 0.28; 95 % of cross-category
  species pairs follow the IUCN order.
- Stability across λ ∈ [1, 25]: SeqDef ρ 0.82 (min 0.78), Priority 0.94 (min 0.89); top-10 Jaccard
  0.25 (SeqDef) / 0.72 (Priority); C. atromarginatus top-1 in 69–95 % of trees at every λ.
- Fig. 3: 85 trees C. atromarginatus alone; 14 trees 8-way tie of CR Rhinobatos + Acroteriobatus;
  1 tree 6-way tie of CR Squatina.
- Brownian kernel = linear kernel at λ = 0.5 on ultrametric trees (max SeqDef difference 1.6e-8).
- Provenance: 1,192 tips → 980 exact IUCN matches → minus 102 DD and 1 EX → 877; 498 molecular /
  379 imputed (Stein et al. 2018 Suppl. Table 1); 55/877 with an assembly (NCBI Assembly, any record
  released ≤ 2025-12-30; reproduced 2026-11-23 for 876/877); IUCN Red List 2025-2, accessed 2025-12-10.
- Full-tree sensitivity (all 1,192 tips, 60 with an assembly, 1,005 assessed): Priority ρ vs main =
  0.985 [0.945, 0.994]; C. atromarginatus unique top-1 in 97/100 trees.

## Analyses

Run from the repository root: `Rscript analyses/run_all.R` (see README for the script list).
`00_setup.R` loads cached data; it never calls external APIs. `analysis.R` is the original download
script and reads API tokens from environment variables only — never write a token into a file.

## Environment

Apple M1 Pro laptop (16 GB, R 4.5.2) for all analyses; TAMU Grace bigmem node (3 TB, R 4.4.2) for the
large-n dense-kernel benchmark in Fig. S3 (`grace_dense_kernels.R`, run with SeqDef v0.1.0).
