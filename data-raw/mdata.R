##########################################
# PREPROCESS 20 MONTH OLD FRONTAL CORTEX #
##########################################

# load in sample metainfomartion from GEO
gset <- GEOquery::getGEO("GSE99354", getGPL = T, AnnotGPL = T)


# we only keep the samples which are also in the gene expression matrix.
# i.e. we remove heterozygous samples
FC_20month <- pData(gset$`GSE99354-GPL15103_series_matrix.txt.gz`) %>%
  dplyr::filter(!grepl("Q331K_\\+", title)) %>%
  dplyr::select(title, geo_accession, instrument_model, `age:ch1`,
                `genotype:ch1`, `strain background:ch1`, `tissue:ch1`)




# reformat the title column, so it matches the column names of the counts data frame
FC_20month$title <- str_replace_all(FC_20month$title,"20month-frontal-cortex\\[\\+_\\+\\]","FC_20months_WT_") %>%
  str_replace_all(.,"20month-frontal-cortex\\[Q331K_Q331K\\]","FC_20months_MT_")

# reformat column names
colnames(FC_20month) <- str_remove_all(colnames(FC_20month),":ch1") %>%
  str_replace_all(., " ","_")








##########################################
# PREPROCESS 5 MONTH OLD FRONTAL CORTEX #
##########################################


# we only keep the samples which are also in the gene expression matrix.
# i.e. we remove heterozygous samples
FC_5month <- pData(gset$`GSE99354-GPL17021_series_matrix.txt.gz`) %>%
  dplyr::filter(`tissue:ch1` == "Frontal cortex") %>%
  dplyr::filter(!grepl("Q331K/\\+", title)) %>%
  dplyr::select(title, geo_accession, instrument_model, `age:ch1`,
                `genotype/variation:ch1`, `strain background:ch1`, `tissue:ch1`)




# reformat the title column, so it matches the column names of the counts data frame
FC_5month$title <- str_replace_all(FC_5month$title,"5month-frontal-cortex\\[\\+/\\+\\]","FC_5months_WT_") %>%
  str_replace_all(.,"5month-frontal-cortex\\[Q331K/Q331K\\]","FC_5months_MT_")

# reformat column names
colnames(FC_5month) <- str_remove_all(colnames(FC_5month),":ch1") %>%
  str_replace_all(., " ","_") %>%
  str_replace_all(.,"genotype/variation", "genotype")


########################
# COMBINE THE METADATA #
########################
mdata <- rbind(FC_5month, FC_20month)

# lastly reorder the rows so they correspond to the column names of the counts matrix
data(counts)
mdata <- mdata[match(colnames(counts),mdata$title),]

stopifnot(all(mdata$title == colnames(counts)))

usethis::use_data(mdata, overwrite = TRUE)
