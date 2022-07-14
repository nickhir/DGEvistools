## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
devtools::load_all()
library(tidyverse)
# load in the example data. For generation of it, check scripts in data-raw directory.
# load raw counts
data(counts)

# load metadata which contains additional information about the sequenced samples, for example different batches, RNA integrity, ...
data(mdata)

# load annotation data, which contains additional information about the genes for which we have expression data
data(annotation)

# each row in the metadata corresponds to exactly one column in the counts data
stopifnot(all(colnames(counts) == rownames(mdata)))

# each row in the annotation table corresponds to one gene (row) in the counts table
stopifnot(all(rownames(counts) == rownames(annotation)))


## ---- warning=F, message=F----------------------------------------------------
library(DESeq2)

# create a `SummarizedExperiemnt`
se <- SummarizedExperiment(assays = list("counts"=counts),
                           colData=mdata,
                           rowData=annotation) 

# Construct DESeqDataSet. Because we are interested in the difference between genotypes we specify the design accordingly
ddsSE <- DESeqDataSet(se, design = ~ genotype)
ddsSE$genotype <- relevel(ddsSE$genotype, ref="WT")


# perform VST transformation for explorative analysis
vsd <- vst(ddsSE)

# plot PCA
pcaData <- plotPCA(vsd, intgroup=c("genotype", "age"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype, shape=age)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

## ---- warning=F, message=F----------------------------------------------------
# as determined above, we perform DGE separately for 5month and 20month old cortex
DGE_result <- list()

for (age in unique(colData(ddsSE)$age)){
  dds_tmp <- ddsSE[,ddsSE$age == age]
  dds_tmp <- DESeq(dds_tmp)
  res <- results(dds_tmp)
  res <- res[order(res$pvalue),]
  
  # add the gene symbol and description
  res <- res %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_id") %>%
    left_join(., data.frame(rowData(ddsSE)), by="ensembl_id")
  
  DGE_result[[age]]<- res
}




## ----volcano plot-------------------------------------------------------------
volcano_plot(de_res = DGE_result$`5 months`, title = "5-month old frontal cortex", 
             subtitle = "6 WT TDP43, 8 MT TDP43",  logFC_threshold = 0.25,
             annotate_by = c("Col4a2", "Fgf12", "Meaf6"), xlim = c(-2,2))

## ----correlation plot---------------------------------------------------------
plot_de_correlation(result_df1 = DGE_result$`5 months`,
                    result_df2 = DGE_result$`20 months`,
                    title_result_1 = "5-month old frontal cortex",
                    title_result_2 = "20-month old frontal cortex",
                    col="stat") 

plot_de_correlation(result_df1 = DGE_result$`5 months`,
                    result_df2 = DGE_result$`20 months`,
                    title_result_1 = "5-month old frontal cortex",
                    title_result_2 = "20-month old frontal cortex",
                    col="baseMean") 

## ----boxplots-----------------------------------------------------------------
# For this, we first transform our raw counts to log2(CPM).
assays(se)$cpm_counts <- edgeR::cpm(se, log=T)

# Simple boxplot using a wilcoxon test to asses significance. 
create_expression_boxplot(gene = "Gng3", SE = se, 
                          intgroup = "genotype")

# # If we have the limma results for the DGE, we can also use the p value that was determined by limma. 
# # This p value is more accurate and should always be used.
# # Make sure that limma was used to compare the same conditions as the boxplot.
# create_expression_boxplot(gene = "IGFBP5", eset = tissue_esets$cervical, 
#                           x_axis_groups = "sample.status", 
#                           test_comparison = list(c("diseased","healthy")), # optional
#                           colors = c("firebrick","steelblue"), #optional 
#                           x_lab = "Sample Phenotype", #optional
#                           limma_results = cervical) #optional
# 
# 
# # It is also very easy to visualize a varity of different genes and save them all in a pdf document:
# plots <- list()
# gene_list <- c("SEC62", "IGFBP5", "KRAS") # any genes you want to visualize.Automatically skips genes that are not found.
# for(gene in gene_list){
#   print(gene)
#   p <- create_expression_boxplot(gene = gene, eset = tissue_esets$cervical, 
#                           x_axis_groups = "sample.status", 
#                           test_comparison = list(c("diseased","healthy")), # optional
#                           colors = c("firebrick","steelblue"), #optional 
#                           x_lab = "Sample Phenotype", #optional
#                           limma_results = cervical) #optional
# 
#   
#   plots[[gene]] <- p
# }
# 
# # save them all as one big pdf to disk.
# pdf_plot <- cowplot::plot_grid(plotlist = plots, ncol=2)
# 
# # here we save it as a random temp file. 
# out_file = tempfile(fileext = ".pdf")
# cowplot::save_plot(out_file, pdf_plot,  base_height = 4.7, ncol = 2,
#           nrow=ceiling(length(plots)/2),
#           limitsize = FALSE)

