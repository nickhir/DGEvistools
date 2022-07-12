##########################################
# PREPROCESS 20 MONTH OLD FRONTAL CORTEX #
##########################################
# load in the raw count data, for the 20 month old frontal cortex tissue. This was downloaded from the GEO side (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112575)

fc_20month <- data.table::fread("C:/Users/Nick/Rprojects/DGEvistools/data-raw/GSE112575_20month_cortex_raw_counts.txt.gz")

# only keep the columns we are interested in, i.e. the ENSEMBL ID and the actual count data
fc_20month_subset <- fc_20month %>%
  dplyr::filter(Description != "No description") %>%
  dplyr::select(ID, contains("20month")) %>%
  dplyr::select(-contains("Q331K/+")) %>% #remove heterozygous samples
  dplyr::distinct(., ID, .keep_all = T) %>%
  tibble::column_to_rownames("ID")



#########################################
# PREPROCESS 5 MONTH OLD FRONTAL CORTEX #
#########################################


# load in the raw count data, for the 5 month old frontal cortex tissue. This was downloaded from the GEO side (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99353)
fc_5month <- data.table::fread("C:/Users/Nick/Rprojects/DGEvistools/data-raw/GSE99353_frontal_cortex_raw_counts.txt.gz")

# only keep the columns we are interested in, i.e. the ENSEMBL ID and the actual count data
fc_5month_subset <- fc_5month %>%
  dplyr::filter(Description != "No description") %>%
  dplyr::filter(ID %in% rownames(fc_20month_subset)) %>%
  dplyr::select(ID, contains("5month")) %>%
  dplyr::select(-contains("Q331K_+")) %>% #remove heterozygous samples
  dplyr::distinct(., ID, .keep_all = T) %>%
  tibble::column_to_rownames("ID")


# combine both counts df into one big matrix
counts <- fc_5month_subset %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(., fc_20month_subset %>% rownames_to_column("ID")) %>%
  tibble::column_to_rownames("ID") %>%
  as.matrix()

# make col names more readable and parseable for R
colnames(counts) <- str_replace_all(colnames(counts),"5month-frontal-cortex_\\+_\\+","FC_5months_WT") %>%
  str_replace_all(.,"5month-frontal-cortex_Q331K_Q331K","FC_5months_MT") %>%
  str_replace_all(.,"20month-frontal-cortex\\[\\+/\\+\\]","FC_20months_WT_") %>%
  str_replace_all(.,"20month-frontal-cortex\\[Q331K/Q331K\\]","FC_20months_MT_")

# this is only supposed to be an example dataset. We reduce the size of it so it runs quicker.
# For the actual data set these steps might be unnecessary
counts <- counts[rowSums(counts)>50,]

# pick the 5k most variable genes. This is not the correct way to do this, because it operates on raw counts instead of normalized ones.
# but for demonstration purposes this is fine.
counts <- counts[order(matrixStats::rowVars(counts), decreasing = T)[1:5000], ]

usethis::use_data(counts, overwrite = T)




