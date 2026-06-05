# SeqDef project — context (for a fresh chat)

## What this is
Analysis/manuscript repo for **SeqDef**, a phylogenetically weighted genomic-prioritization statistic
+ R package. Manuscript *"SeqDef: An R package for phylogenetically weighted genomic prioritization
across the tree of life"* — **Journal of Evolutionary Biology, JEB-2026-00140, Major Revision**.

Two local repos (both on `main`, both pushed to GitHub, **public**):
- Package: `/Users/andresbarboza/Documents/GitHub/SeqDef` → github.com/Andresdbp/SeqDef
- Project: `/Users/andresbarboza/Documents/GitHub/SeqDef_project` → github.com/Andresdbp/SeqDef_project

Manuscript + rebuttal live on **Google Drive** (read via the Drive connector; can't edit in place):
- "manuscript" Doc id `1DF6OUwHqvZFZgv40ngG5BTNBHvt7WZP773rATRwoYzY`
- "Answer to reviewers" Doc id `11-lz-Sv1BZSPetPTl99MzyfvryxcNLXj_StHhHsldfI` (currently holds only the
  editor/reviewer letter — our point-by-point response needs to be added there).

## Status: the revision is essentially built; the user is hand-applying edits to the Drive manuscript.

### Package (done, committed, pushed)
- `SeqDef(tree, df, data.col, invert, scale, lambda, kernel)` with **kernel = exponential (default) /
  gaussian / linear / brownian**.
  - **exponential** → exact **O(n)** two-pass tree traversal (`.seqdef_exp_avail`), no distance matrix.
  - **gaussian, linear** → dense O(n²) cophenetic path.
  - **brownian** → parameter-free phylogenetic correlation (`.seqdef_bm_avail`), also **O(n)**; ignores λ.
  - All verified identical to the dense computation (≤1e-15). Default behaviour unchanged.
- `lambda` = `"auto_max"` (variance-max scan 1–50) / `"by_genus"` (half-life = median intra-genus dist) / numeric.
- `calc_priority(seqdef_res, trait, model="exponential", base=2)` = SeqDef × base^trait.
- DESCRIPTION author = **Andres Barboza <andresdbp00@gmail.com>** (renders "Barboza, A.").
  README credits **Dr. Andres Barboza** + a citation line.
- NOTE: DESCRIPTION says `License: MIT + file LICENSE` but **no LICENSE file exists** (held per user;
  R CMD check will warn). `function.R` (project) mirrors the package SeqDef + calc_priority.

### Analyses (project `analyses/`, run from project root; `00_setup.R` shared; NEVER source analysis.R = keys)
01 λ rank-stability + method head-to-head · 02 kernel comparison · 03 EDGE/EDGE2 benchmark ·
04 runtime+memory (local) · 05 Fig 1 regen · 06 continuous-S (**DROPPED from manuscript**) ·
07 traversal O(n) scaling · 08 brownian. Grace HPC scripts: `grace_runtime_compute.R`, `grace_runtime.slurm`,
`recover_grace_runtime.R`. Outputs in `results/` (+ `SUMMARY.md` = all headline numbers).

### Headline results (reproduced; in `results/SUMMARY.md`)
- **n = 877** species (phylo∩IUCN; 55 have an NCBI assembly). *C. atromarginatus* top priority in **85/100** trees.
- λ: auto_max 3.40 [2.65,5.00], by_genus 3.19; two methods ρ=**0.998**, agree on top target **100/100**;
  Priority ranking ρ vs auto_max ≥0.88 (median 0.93) across λ∈[1,25].
- Kernels: Priority ρ vs exponential = 0.88 (gaussian), 0.94 (linear); top target preserved.
- EDGE: **ρ(EDGE, SeqDef)=0.07** (orthogonal); ρ(EDGE, Priority)=0.83; top-25 overlap 9/25. Divergence:
  *Sphyrna lewini* high-EDGE/low-SeqDef (congener *S. mokarran* sequenced); *C. atromarginatus* high on both.
- Efficiency: dense **O(n²)** (447 GB at n=100k on Grace 3TB node); **exponential O(n) traversal** does
  n=10⁶ in <1 s / ~0.4 GB on a laptop. 877-tip run 0.02 s (fixed λ) / 0.21 s (auto_max); 100-tree posterior ~40 s.
- Brownian = flat low-λ limit of exponential (ρ≈0.98 at λ≈0.5); less discriminating; different top target.

### Figures
- Fig 1 regenerated (λ-annotated; `figures/fig1.pdf`).
- Supplementary **renumbered S1–S5**: S1 `figS_lambda_stability` · S2 `figS_kernel` ·
  S3 `figS_edge_vs_seqdef` · S4 `figS_runtime` (dense O(n²) vs traversal O(n)) · S5 `figS_brownian`.
  `figS_continuousS.pdf` exists but is **dropped from the submission**.

### Key decisions
- **n = 877** (the 850 in the original draft was an error).
- **Exponential is the default**, justified four ways: (i) grounded in sequence evolution (shared identity
  decays exponentially with divergence; Jukes & Cantor 1969); (ii) tunable half-life horizon; (iii) standard
  in distance-decay/sequence-novelty work; (iv) the only kernel that is exact **O(n)** *and* tunable
  (Brownian is O(n) but parameter-free; gaussian/linear are O(n²)).
- **Continuous-S demo dropped** to avoid bloat — binary presence/absence is a legitimate, simple metric, and
  the Methods worked-S constructions already answer "examples beyond 0/1." Keep figures/analyses lean.
- **Brownian kept** (one supp figure) because it sharpens the "why exponential" argument.

## What remains
- **`revision/manuscript_changes.md`** now lists **only the pending edits** (equations/λ glyphs, Priority eq,
  doubling sentence, Sensitivity rank-stability, EDGE Methods+Results, Computational-benchmarks subsection,
  Brownian Discussion fix, rtrees paragraph, Conclusions, Fig 1 caption + supp captions, refs, compliance).
  The user is applying these to the Drive manuscript by hand.
- **`revision/SeqDef_response_to_reviewers_FINAL.md`/.docx** = the point-by-point response to add to the
  "Answer to reviewers" Drive Doc (user will paste/upload; do NOT auto-create a Drive doc).
- New refs via Paperpile: **Jukes & Cantor 1969, Gumbs et al. 2023, Kembel et al. 2010** (Soares 2023 already cited).
- User action items: insert figures into the Doc; CRediT; Data Availability + mint Zenodo DOI.

## ⚠️ Security
`analysis.R` had real NCBI + IUCN API keys. The working copy is scrubbed, **but the keys are still in the
public git history** (commit 126814d). **The user must rotate/revoke them at NCBI and IUCN.** Never echo or
re-commit keys; never source `analysis.R` in analyses (load the cached `.rds`/`.nex` directly).

## Git
Both repos on `main` (pushed). Branches: `main`, `revision-jeb`, `linear-time-exp`, `brownian-kernel`.
End commit messages with: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## Env
R 4.5.2 local (Apple Silicon, 16 GB). Grace HPC via `ssh grace` (TAMU; module `gfbf/2024a R/4.4.2`,
`bigmem`=3TB/80c; user keeps an SSH ControlMaster open). Pandoc available for md→docx.
`REVISION_ANALYSES.md` has the original analysis build spec.
