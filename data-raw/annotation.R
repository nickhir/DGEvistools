###################################################################
# Get information about the genes we have counts for form Ensembl #
###################################################################
library(biomaRt)

# for creation see counts.R
data(counts)

gene_ids <- rownames(counts)


ensembl = useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  verbose = T,
  mirror = "useast"
)

# In our case, we have ensembl gene ids and want to map them to the gene symbol and the
gene_attributes = c("ensembl_gene_id",
                    "mgi_symbol",
                    "description")

annotation <- getBM(
  attributes = gene_attributes,
  values = gene_ids,
  filters = "ensembl_gene_id",
  mart = ensembl
)

# reformat the annotations we got
annotation$description <- str_remove_all(annotation$description, " \\[.*\\]")

colnames(annotation) <- c("ensembl_id", "symbol", "description")

# for 21 genes, we did not find any additional information, so they were removed. We add them again
annotation <- dplyr::bind_rows(annotation, data.frame("ensembl_id"=gene_ids [!gene_ids %in% annotation$ensembl_id]))
data(counts)
# sort rows between counts and annotation match
annotation <- annotation[match(rownames(counts),annotation$ensembl_id),]

stopifnot(all(annotation$ensembl_id == rownames(counts)))

rownames(annotation) <- annotation$ensembl_id


usethis::use_data(annotation, overwrite = TRUE)

