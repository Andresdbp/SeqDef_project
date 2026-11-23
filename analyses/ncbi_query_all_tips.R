# =============================================================================
# ncbi_query_all_tips.R   (needs internet + rentrez)
# Documents and reproduces the NCBI Assembly query behind the case study
# (Reviewer 2, Major 4).
#
# Query (NCBI E-utilities, database "assembly"):
#     esearch  term = "Chondrichthyes[Organism]"     -> every chondrichthyan
#              genome-assembly record (all assembly levels, GenBank and RefSeq,
#              any status: latest, replaced or suppressed)
#     esummary on those ids                         -> species, taxid, assembly
#              level, GenBank release date, status flags
#     efetch   db = "taxonomy" on the species taxids -> NCBI synonyms
# Snapshot rule: an assembly counts if its GenBank release date is on or before
# SNAPSHOT (2025-12-30, the date of the R1 data table). A species is coded
# S = 1 if it has >= 1 such record ("any genome assembly in NCBI").
# This reproduces the R1 table (data/ncbi_assembly_data.rds) for all 55
# sequenced species (presence/absence; see ncbi_map_all_tips.R).
# The query used for the R2 revision was run on 2026-11-23 (222 records).
#
# Output: data/ncbi_assembly_species_summary.tsv (one row per NCBI species)
# =============================================================================
suppressPackageStartupMessages({ library(rentrez); library(xml2) })
SNAPSHOT <- as.Date("2025-12-30")
if (nzchar(Sys.getenv("ENTREZ_KEY"))) set_entrez_key(Sys.getenv("ENTREZ_KEY"))

ids <- entrez_search(db = "assembly", term = "Chondrichthyes[Organism]", retmax = 10000)$ids
message(length(ids), " assembly records")
sm <- unlist(lapply(split(ids, ceiling(seq_along(ids) / 150)), function(b) {
  Sys.sleep(0.4); entrez_summary(db = "assembly", id = b, always_return_list = TRUE) }), recursive = FALSE)
rec <- do.call(rbind, lapply(sm, function(x) data.frame(
  accession = x$assemblyaccession, species = x$speciesname, taxid = x$speciestaxid,
  level = x$assemblystatus, release = as.Date(substr(x$asmreleasedate_genbank, 1, 10), "%Y/%m/%d"),
  latest = "latest" %in% unlist(x$propertylist), stringsAsFactors = FALSE)))
lvl <- c(Contig = 1, Scaffold = 2, Chromosome = 3, "Complete Genome" = 4)
out <- do.call(rbind, lapply(split(rec, rec$species), function(d) {
  k <- !is.na(d$release) & d$release <= SNAPSHOT
  data.frame(ncbi_species = d$species[1], taxid = d$taxid[1],
             n_assemblies_le_snapshot = sum(k), n_latest_le_snapshot = sum(k & d$latest),
             best_level_le_snapshot = if (any(k)) names(lvl)[max(lvl[d$level[k]])] else "",
             first_release = if (any(k)) format(min(d$release[k]), "%Y/%m/%d") else "",
             n_assemblies_at_query = nrow(d))
}))
tx <- read_xml(entrez_fetch(db = "taxonomy", id = unique(out$taxid), rettype = "xml"))
syn <- sapply(xml_find_all(tx, "/TaxaSet/Taxon"), function(t) {
  nm <- xml_text(xml_find_all(t, "./OtherNames/Synonym | ./OtherNames/Name/DispName"))
  paste(unique(grep("^[A-Z][a-z]+ [a-z-]+$", nm, value = TRUE)), collapse = ";") })
names(syn) <- xml_text(xml_find_all(tx, "/TaxaSet/Taxon/TaxId"))
out$ncbi_synonyms <- unname(syn[as.character(out$taxid)]); out$ncbi_synonyms[is.na(out$ncbi_synonyms)] <- ""
write.table(out[order(out$ncbi_species), ], "data/ncbi_assembly_species_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d species; %d with >= 1 assembly released on or before %s\n", nrow(out), sum(out$n_assemblies_le_snapshot > 0), SNAPSHOT))
