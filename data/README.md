# Data sources

| File | Content | Source | Terms |
| --- | --- | --- | --- |
| `chondrichthyes.nex` | 100 trees from the posterior distribution of the taxon-complete chondrichthyan phylogeny (1,192 species) | Stein et al. (2018) *Nat. Ecol. Evol.* 2: 288–298; VertLife (vertlife.org/data/sharks) | Cite Stein et al. (2018) |
| `stein2018_data_type.csv`, `stein2018_genetic_species.csv` | Whether each species was placed with molecular data (610) or by taxonomic constraints (582) | Stein et al. (2018), Supplementary Table 1 ("Data type" column) | Cite Stein et al. (2018) |
| `ncbi_assembly_data.rds` | Assembly counts and presence/absence for the 877 analysed species (table used in the paper; compiled 30 December 2025) | NCBI Assembly database | Public domain (NCBI) |
| `ncbi_assembly_species_summary.tsv` | All chondrichthyan assembly records summarised per NCBI species, with release-date snapshot counts (≤ 2025-12-30) and NCBI synonyms; query run 23 November 2026 | NCBI E-utilities, `analyses/ncbi_query_all_tips.R` (term `Chondrichthyes[Organism]`) | Public domain (NCBI) |
| `ncbi_assembly_all_tips.csv` | The snapshot mapped onto all 1,192 tree tips (exact name, NCBI synonym, or IUCN synonym) | `analyses/ncbi_map_all_tips.R` | Public domain (NCBI) |
| `iucn_assessment_data.rds` | **Not in this repository.** Latest global Red List assessments for class Chondrichthyes (1,266 records), Red List version 2025-2, accessed 10 December 2025 | IUCN Red List API v4 via the `iucnredlist` R package | IUCN Red List Terms of Use prohibit redistribution of Red List data |

## Re-creating the IUCN file

The IUCN Red List terms of use do not allow Red List data to be reposted or redistributed, so
the assessment file is not included. To reproduce it, obtain a free API token from
https://api.iucnredlist.org, set it as an environment variable (never write it into a script),
and run Step 1 of `analysis.R`:

```r
Sys.setenv(IUCN_REDLIST_KEY = "<your token>")   # or put it in ~/.Renviron
library(iucnredlist)
api <- init_api(Sys.getenv("IUCN_REDLIST_KEY"))
ids <- assessments_by_taxonomy(api, level = "class", name = "chondrichthyes",
                               latest = TRUE, scope_code = 1, wait_time = 2)
a_data <- assessment_data_many(api, ids$assessment_id, wait_time = 1)
saveRDS(a_data, "data/iucn_assessment_data.rds")
```

Assessments change between Red List versions, so a new download may differ from the
2025-2 version used in the paper.
