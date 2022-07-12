#####
# basic theme
#####
theme_jh <- function () {
  theme_bw() %+replace%
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(),
      strip.text = element_text(color = "black"),
      axis.text = element_text(colour = "black", size=14),
      axis.title = element_text(colour = "black", size=18),
      panel.background  = element_blank(),
      plot.background = element_rect(fill="white", colour=NA),
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      axis.line.y = element_line(), strip.text.x = element_text(face = "bold", margin = margin(t = 2,r = 0,b = 2,l=0))
    )
}


#' @title Plot a pretty volcano plot
#' @description This function plots a pretty volcano plot using the results of the output of \code{msdbulkseqtools::toptable_for_all}
#'
#' @param de_res A dataframe created using the \code{msdbulkseqtools::toptable_for_all} command. Alternatively a data frame
#' containing the following columns: 'symbol' (corresponding to gene name), 'padj', 'log2FoldChange'
#' @param title The title of the plot
#' @param subtitle The subtitle of the plot
#' @param annotate_by A vector with gene names which will be annotated in the volcano plot.
#' Gene names have to occur in the 'symbol' column of 'de_res'
#' @param annotation_size The font size of the annotation.
#' @param padj_threshold Threshold determining which genes should be counted as differentially expressed
#' @param logFC_threshold Threshold determining which genes should be counted as "strongly" differentially expressed
#' @param ymax Maximum of the y axis
#' @param xlim The x axis limits.
#' @return A ggplot of a volcano plot.
#'
#' @examples
#' \dontrun{
#' Also see vignette `Standard bulk RNA seq showcase`
#' res <- limmaVoomFit(v, design, cmat)
#' res_table <- msdbulkseqtools::toptable_for_all(res)
#' volcano_plot(res_table, annotate_by=c("SOD1"))
#' }
#'
#' @import ggplot2
#' @import grid
#' @importFrom dplyr group_by tally
#'
#' @export
volcano_plot <- function(de_res, title = NULL, subtitle = NULL, annotate_by = NULL, annotation_size=5.2,
                         padj_threshold=0.05, logFC_threshold=1, ymax=NULL, xlim=c(-3,3), res=300){

  # check if all required columns exist.
  if (!all(c("padj","log2FoldChange","symbol") %in% colnames(de_res))){
    stop("ERROR: `de_res` requires the following columns: `padj`,`log2FoldChange`,`symbol`")
  }

  # add some extra annotation about the magnitude of differential expression using log2FoldChange
  de_res <-
    dplyr::mutate(de_res,
                  sig = case_when(
                    padj >= padj_threshold ~ "non_sig",
                    padj < padj_threshold & abs(log2FoldChange) < logFC_threshold ~ "sig",
                    padj < padj_threshold & abs(log2FoldChange) >= logFC_threshold ~ "sig - strong"
                  )) %>%
    dplyr::mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
    dplyr::mutate(class = paste(sig, direction)) %>%
    dplyr::mutate(log2FoldChange = case_when(
      log2FoldChange > 3 ~ Inf,
      log2FoldChange < -3 ~ -Inf,
      TRUE ~ log2FoldChange
    ))

  # get number of differentially expressed genes and determine position where they should be plotted.
  de_tally <- de_res %>% group_by(sig, direction, class) %>% tally() %>%
    dplyr::filter(sig != "non_sig") %>%
    dplyr::mutate(position = ifelse(sig == "sig", 0.5, 2) ) %>%
    dplyr::mutate(position = ifelse( direction == "down", -1 * position, position)) %>%
    dplyr::mutate(n = formatC(n, format="f", big.mark=",", digits=0))

  # determine a suitable ymax if none was specified
  if (is.null(ymax)){
    ymax <- ifelse(min(de_res$padj, na.rm=T) < 1e-16, 18, -log10(min(de_res$padj, na.rm=T))+2)
  }

  # generate the volcano plot
  plot <- de_res %>%
    dplyr::mutate(padj = ifelse(padj < 1e-16, 1e-16, padj)) %>% #threshold at 1e16
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    #geom_point(aes(colour = class ), size = 0.5) +
    ggrastr::rasterise(geom_point(aes(colour = class), size = 0.8), dpi = res) + # otherwise if we plot thousands of individual points takes to long
    scale_colour_manual(values = c("non_sig up" = "gray",
                                   "non_sig down" = "gray",
                                   "sig up" = "#EB7F56",
                                   "sig - strong up" = "#B61927",
                                   "sig down" = "#75A4EF",
                                   "sig - strong down" = "#063A8C"
    )) +
    labs(y = expression(-log[10]~adj~P~value), x = expression(log[2]~"(fold change)"), title = title, subtitle = subtitle) +
    guides(colour = "none") +
    scale_y_continuous(expand = c(0,0), limits = c(0,ymax)) +
    theme_jh() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 19),
      plot.subtitle = element_text(hjust = 0.5, size=14),
      panel.border = element_blank(),
      axis.ticks = element_line(colour = "black")
    ) +
    scale_x_continuous(limits = xlim, breaks=scales::pretty_breaks(n = sum(abs(xlim))))

  # if we found any differentially regulated genes, add them to the plot
  if(nrow(de_tally)>0){
    plot <- plot +
      geom_text(fontface = "bold", data = de_tally, aes(x = position, y = ymax - 0.7, label = n, colour = class), size = 5.8 )
  }

  # add annotation if requested.
  if(!is.null(annotate_by)){
    # Check if all the genes we want to annotate occure in the column
    if (sum(de_res$symbol %in% annotate_by) != length(annotate_by)){
      stop("Specified gene names for annotation do not occure in the symbol column of `de_res`.")
    }
    plot <- plot +
      ggrepel::geom_text_repel(
        fontface = "italic",
        data = dplyr::filter(de_res, symbol %in% annotate_by),
        aes(x = log2FoldChange, y = -log10(padj), label = symbol),
        min.segment.length = unit(0, "lines"),
        size = annotation_size) +
      geom_point(
        data = dplyr::filter(de_res, symbol %in% annotate_by), size = 1.6, colour = "black"
      ) +
      geom_point(aes(colour = class ),
                 data = dplyr::filter(de_res, symbol %in% annotate_by), size = 1
      )

  }
  return(plot)

}





#' Plot correlation between two differential gene expression analyses
#' @description This function visualizes the correlation between two different differential gene expression analyses.
#' Might be used to compare two different tissues or two different timepoints and see how similar they behave.

#' @param result_df1  The first results dataframe. Usually created using the \code{msdbulkseqtools::toptable_for_all} command.
#' @param result_df2 The second results dataframe. Usually created using the \code{msdbulkseqtools::toptable_for_all} command.
#' @param title_result_1 Name of the first differential gene expression analysis. Will be the title on the x-axis.
#' @param title_result_2  Name of the second differential gene expression analysis. Will be the title on the y-axis.
#' @param col Name of the column which will be used to perform the correlation analysis.
#' For example \code{log2FoldChange} Has to occure in both \code{result_df1} and \code{result_df2}
#' @examples
#' \dontrun{
#' Also see vignette `Standard bulk RNA seq showcase`
#' res1 <- limmaVoomFit(v, design, cmat)
#' res2 <- limmaVoomFit(v, design, cmat)
#' res_table1 <- msdbulkseqtools::toptable_for_all(res1)
#' res_table2 <- msdbulkseqtools::toptable_for_all(res2)
#' plot_de_correlation(res_table1, res_table2, "first experiment", "second experiment")
#' }
#'
#' @export

plot_de_correlation <- function(result_df1, result_df2, title_result_1, title_result_2, col = "log2FoldChange"){
  # combine dataframes

  res <- dplyr::left_join(result_df1, result_df2,  by = "gene_id" , suffix = c(".1", ".2") )

  # when joining we add a suffix, so we also have to do it here
  x_string <- paste0(col, ".1")
  y_string <- paste0(col, ".2")

  # run linear regression to determine the slope and intercept.
  correlation_eq <- broom::tidy(lm(as.formula(paste0(y_string, " ~ ", x_string)), res))
  intercept <- correlation_eq[1,2] %>% pull()
  slope <- correlation_eq[2,2] %>% pull()


  plot <- res %>%
    ggplot(aes_string(x = x_string, y = y_string )) +
    labs(x = title_result_1, y = title_result_2) +
    theme_bw() +
    ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..) ) +
    geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
    geom_abline(slope = slope , intercept = intercept, linetype = 3, color="red") +
    geom_bin2d(bins = 300) +
    scale_fill_continuous(type = "viridis") +
    guides(fill = "none") +
    theme_jh()

  return(plot)
}



#' Generate a heatmap of GSEA results
#' @description This function can be used to compare GSEA results performed for different tissues, different timepoints or conditions.
#' It will generate a heatmap, where rows represent different gene sets and columns the different comparisons.
#'
#' @param gsea_df A dataframe which contains the results of \code{msdbulkseqtools::glattCameraPR}. For creation see \code{examples}
#' @param comparisons_column The name of the column which contains the different comparisons that you want to make.
#' @param comparisons A vector of names which correspond to the columns of the heatmap (i.e. the conditions you want to compare).
#' @param ordering By which condition should the p values of the displayed gene sets be ranked
#' @param top_n The number of gene sets you want to display.
#' @param data_src The data source for which you want to display the top enriched terms.
#' @param color A colormap for the heatmap
#' @param legend_title The title of the legend
#' @param title The title of the heatmap
#' @param remove_string A string that should be removed from all rownames of the heatmap.
#' @param ... Arguments are directly passed to \code{ComplexHeatmap::Heatmap}
#'
#' @examples
#' \dontrun{
#'Also see vignette `Standard bulk RNA seq showcase`
#' gsrc <- c("go_component","KEGG")
#' gsea <- msdbulkseqtools::glattCameraPR(g_src = gsrc,
#'                                       efit_list = res_list,
#'                                       specie = "Hs",
#'                                       g_cor_list = tcor_list,
#'                                       gs_low = 5,
#'                                       gs_high = 500)
#'
#' gsea_df <- lapply(gsea, function(Ontology){
#'   lapply(names(Ontology), function(Tissue){
#'     lapply(names(Ontology[[Tissue]]),function(Contrast){
#'       Ontology[[Tissue]][[Contrast]] %>% mutate(contrast = Contrast)
#'       }) %>% dplyr::bind_rows() %>% mutate(tissue = Tissue)
#'    }) %>% bind_rows()
#'  }) %>% bind_rows()
#'
#'enrichment_heatmap(gsea_df = gsea_df, comparisons = c("cervical","lumbar"),
#'                   top_n = 20, data_src = "go_component", comparisons_column = "tissue",
#'                   title = "Tissue comparison")
#' }
#' @export
enrichment_heatmap <- function(gsea_df, comparisons_column, comparisons=NULL,
                               ordering = NULL,top_n=10, data_src="KEGG",
                               color=viridis::plasma(50,end=0.75),
                               legend_title="-log10(FDR)", title="",
                               remove_string=NULL, ...){


  # this gets called by `cell_fun` of ComplexHeatmap and returns asteriscs according to p value.
  annotate_sig <- function(j, i, x, y, w, h, fill) {
    # case_when thorws an error and i dont know why.. this works however
    if(is.na(wide_matrix[i, j])) {
      grid::grid.text("", x, y, gp=grid::gpar(col="white", fontsize=25))
    } else if(wide_matrix[i, j] > -log10(1e-4)) {
      grid::grid.text("***", x, y, gp=grid::gpar(col="white", fontsize=25))
    } else if(wide_matrix[i, j] > -log10(1e-3)) {
      grid::grid.text("**", x, y, gp=grid::gpar(col="white", fontsize=25))
    } else if(wide_matrix[i, j] > -log10(1e-2)) {
      grid::grid.text("*", x, y, gp=grid::gpar(col="white", fontsize=25))
    }
  }


  # if no comparisons selected, take all values in the comparisons column
  if(all(is.null(comparisons))){
    comparisons <- unique(gsea_df[comparisons_column])
  }

  if (is.null(ordering)){
    ordering <- comparisons[1]
  }

  if (sum(gsea_df$src == data_src) ==0){
    stop("The 'data_src' you specified does not occure in 'gsea_df'")
  }

  subset <- gsea_df %>%
    dplyr::filter(src == data_src)

  # get the top n enriched terms for a selected tissue
  ontologies_of_interest <- subset %>%
    dplyr::filter(get(comparisons_column)==ordering)%>%
    dplyr::arrange(FDR) %>%
    slice_head(n=top_n) %>%
    pull(SUID)

  # get the p values of the selected ontologies for all other comparisons.
  wide_df <- subset %>%
    dplyr::filter(SUID %in% ontologies_of_interest) %>%
    dplyr::select(c("SUID","Description", "neglog10FDR",comparisons_column)) %>%
    tidyr::pivot_wider(names_from = comparisons_column,
                       values_from = neglog10FDR)

  # this matrix is necessary to annotate significance.
  wide_matrix <- wide_df %>%
    distinct(Description,.keep_all = T) %>%
    tibble::column_to_rownames("Description") %>%
    dplyr::select(all_of(comparisons)) %>%
    dplyr::arrange(desc(get(ordering))) %>%
    as.matrix()

  if (!is.null(remove_string)){
    rownames(wide_matrix) <- str_remove(rownames(wide_matrix), remove_string)
  }


  # Reformat the description, to make it more readable.
  wide_matrix <- wide_matrix %>%
    data.frame() %>%
    tibble::rownames_to_column("name") %>%
    dplyr::mutate(name = str_replace_all(name, "_"," ")) %>%
    dplyr::mutate(name = str_wrap(name, 30)) %>%
    tibble::column_to_rownames("name") %>%
    as.matrix()


  # calculate an appropriate width of the cells
  names_length <- max(nchar(colnames(wide_matrix)))
  cellwidth = names_length*3
  p <- ComplexHeatmap::Heatmap(wide_matrix, name=legend_title, cluster_rows = F,
                               cluster_columns = F,
                               border=T, width = ncol(wide_matrix)*unit(cellwidth, "mm"),
                               col = color, column_title=title,
                               column_names_rot = 45, cell_fun = annotate_sig,
                               ...)

  return(p)

}




#' Heatmap of counts after running GWENA
#' @description After running GWENA (WGCNA) this function can be used to visualize the results in a heatmap.
#' Here, the rows get split into the different modules and columns into either predefined columns clusters or they get clustered using \code{hclust}.
#'
#' @param scaled_counts A matrix of z-transformed counts. Columns should correspond to samples and rows to genes.
#' @param mdata A dataframe with information about the subjects.
#' @param info_col A column in the \code{mdata} object, which has entries that match the column names of \code{scaled_counts}.
#' They have to be in the same order, otherwise an error will occure.
#' @param column_cluster A named vector of subjects and corresponding cluster. Use only if you already determined how to cluster patients.
#' If none is specified, \code{hclust} is used to perform the clustering.
#' @param GWENA_modules Absolute path to the identified co-expression modules. Generated using \code{GWENA::detect_modules()} followed by \code{saveRDS()}.
#' The whole \code{GWENA} workflow has to be performed using the same subjects and genes that are in the \code{scaled_counts}.
#' @param heatmap_annotation A \href{https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html}{HeatmapAnnotation object}.
#' @param ... Other arguments that are passed directly to the \code{ComplexHeatmap::Heatmap} function. Might be something like \code{show_row_names=T}
#' @return A \code{ComplexHeatmap} which shows the GWENA results.
#' @examples
#' \dontrun{
#' See vignette `Standard bulk RNA seq showcase`
#' scaled_counts <- vst_counts %>% t() %>% scale() %>% t() # first have to vst normalize counts
#' module_save_path <- "/path/to/file/with/detected/modules" # saved output of GWENA::detect_modules. Same genes that are in scaled counts
#' mdata <- fread("path/to/file/with/metadata")
#' gwena_heatmap(scaled_counts = scaled_counts,
#' mdata = mdata, info_col = "sample.ID",
#' GWENA_modules = module_save_path,
#' heatmap_annotation = NULL)
#' }
#' @export
#'
gwena_heatmap <- function(scaled_counts, mdata, info_col,
                          column_cluster=NULL, GWENA_modules,
                          heatmap_annotation=NULL, ...){


  # only keep metadata for the samples that we are interested in
  mdata <- mdata %>%
    dplyr::filter(get(info_col) %in% colnames(scaled_counts))

  if (!all(dplyr::pull(mdata, info_col) == colnames(scaled_counts))){
    stop("Rows of `mdata` do not match to columns of `scaled_counts`")
  }


  # determine row split which is determined by the identified GWENA modules
  # read in precalculated identified modules
  modules <- readRDS(GWENA_modules)

  # reorder rows of the count matrix so it matches the genes in module
  module_genes <- modules$modules %>%
    unlist()

  scaled_counts <- scaled_counts[module_genes,] %>% as.matrix()

  # create a factor of modules, will be needed while splitting heatmap by module
  module_membership <- lapply(1:length(modules$modules),function(i){
    rep(names(modules$modules)[[i]],length(unlist(modules$modules[i])))
  }) %>% unlist() %>% as.factor()


  # check if we predetermined column clusters.
  if (is.null(column_cluster)){
    heatmap <-  Heatmap(scaled_counts,
                        show_row_names = F,
                        show_row_dend = F,
                        cluster_rows = F,
                        row_split = factor(module_membership, levels = unique(module_membership)),
                        row_gap = unit(2, "mm"),
                        show_column_dend = T,
                        show_column_names = F,
                        name="Zscore",
                        column_gap = unit(1.5, "mm"),
                        border = T,
                        top_annotation=heatmap_annotation)
    return(heatmap)
  } else {
    heatmap <-  Heatmap(scaled_counts,
                        show_row_names = F,
                        show_row_dend = F,
                        cluster_rows = F,
                        column_split = column_cluster,
                        row_split = factor(module_membership, levels = unique(module_membership)),
                        row_gap = unit(2, "mm"),
                        show_column_dend = T,
                        show_column_names = F,
                        name="Zscore",
                        column_gap = unit(1.5, "mm"),
                        border = T,
                        top_annotation=heatmap_annotation)

    return(heatmap)
  }


}


#' Create an expression boxplot to compare two conditions.
#' @description Given a properly formatted \code{ExpresisonSet}, this function can quickly plot boxplots to investigate the expression of specific genes between two conditions.
#' The \code{ExpresisonSet} can be easily generated from the \code{DGElist} which is used by \code{limma}. See example.
#' @param gene The gene symbol of the gene for which the expression should be plotted. Must occur in the \code{symbol} column of \code{fData} of the \code{ExpressionSet}.
#' @param eset The \code{ExpresisonSet}. For generation, see example.
#' @param x_axis_groups Specify which two conditions should be compared
#' @param test_comparison Optionally, you can also specify if and how conditions should be compared. Should be a list which contains the comparisons that you are interested in.
#' See examples and \href{https://cran.r-project.org/web/packages/ggsignif/vignettes/intro.html}{here}.
#' @param limma_results By default, if you specify \code{test_comparison} a \code{wilcox.test} is performed.
#' This might be fine to get an intuition if expression differences between two conditions
#'  are significantly different, but the proper way is to use the the p value that
#'  was determined using \code{limma}. For that, you can simply specify an object that was created \code{msdbulkseqtools::toptable_for_all} that contains the correct contrast/comparison.
#' @param colors Specify the colors of the boxplots. Number of colors must match number of factors in \code{x_axis_groups}.
#' @param x_lab X axis label. Defaults to \code{x_axis_groups}
#' @param y_lab Y axis label. Defaults to \code{log(CPM)}
#' @param ymax Max value of the Y axis. Sometimes useful to pick a different value than the default.
#' @param ymin Min value of the Y axis. Sometimes useful to pick a different value than the default.
#'
#' @return A ggplot of a boxplot
#'
#' @examples
#' \dontrun{
#' See vignette `Standard bulk RNA seq showcase`
#'
#' eset <- ExpressionSet(assayData = counts, # Normalized counts
#'        phenoData = AnnotatedDataFrame(mdata), # information about the samples
#'        featureData = AnnotatedDataFrame(ann_uniq)) # information about the genes. Must have a \code{symbol} column which corresponds to gene name
#'
#'
#' create_expression_boxplot(gene = "IGFBP5", eset = expressionSet,
#' x_axis_groups = "sample.status",
#' test_comparison = list(c("diseased","healthy")), # optional
#' colors = c("firebrick","steelblue"), #optional
#' x_lab = "Sample Phenotype", #optional
#' ymax=11, ymin=5) #optional
#' }
#'
#' @export
#'
create_expression_boxplot <- function(gene, eset, x_axis_groups,
                                      test_comparison=NULL, limma_results=NULL,
                                      colors=NULL, x_lab=ggplot2::waiver(),
                                      y_lab=expression(log[2](CPM)),
                                      ymax=NULL, ymin=NULL){

  keep <- fData(eset)$symbol == gene
  if (sum(keep) == 0){
    print("Gene symbol was not found in eset, returning empty plot")
    return()
  }

  # check if the specified `test_comparisons` exist in the eset
  for (comparison in test_comparison %>% unlist()){
    if (!(comparison %in% pData(eset)[x_axis_groups][,1])){
      stop(paste0("The specified `test_comparison` (", comparison,") does not occure in `x_axis_groups`"))
    }
  }

  # only keep information for the gene we are interested in.
  subset_eset <- eset[keep,]

  # turn eset into a long df. Eset was basically required to perform effective selection of gene of interest.
  counts <- t(exprs(subset_eset)) %>%
    set_colnames(gene)

  # Do some sanity checks.
  pdata <- pData(subset_eset)
  stopifnot(rownames(counts)==rownames(pdata))
  rownames(counts) <- NULL
  rownames(pdata) <- NULL
  fdata <- fData(subset_eset)
  rownames(fdata) <- NULL

  # create a long df which contains information about the expression as well as additional information such as description and meta data.
  plot_df <- cbind(fdata,counts, pdata)

  # in very few edgecases:
  if (n_distinct(colnames(plot_df))!=ncol(plot_df)){
    colnames(plot_df) <- make.names(colnames(plot_df),unique = T)
  }

  # get gene description for the plot
  gene_descr <- plot_df %>%
    pull("description") %>%
    head(1) %>%
    str_wrap(40)


  # make sure that the `x_axis_groups` are factors.
  plot_df[,x_axis_groups] <- factor(plot_df[,x_axis_groups])

  # If user did not specify anything try and determine appropriate  axis limits.
  # The thing is that we have to leave some space for the annotation (number of samples and significance),
  # so the default values ggplot picks might be unsuitable.
  if(is.null(ymin)){
    ymin <- plot_df %>%
      slice_min(get(gene),n=1) %>%
      pull(gene) - 1.5
  } else{ymin=ymin}

  # this can clearly be optimized but seem to work fine for now. Also user can just pick what looks good...
  if (is.null(ymax)){
    ymax <- plot_df %>%
      slice_max(get(gene),n=1) %>%
      pull(gene) + 0.5

    ymax <- case_when(
      ymax > 5 ~ ymax + 1,
      ymax > 10 ~ ymax + 2,
      ymax > 20 ~ ymax + 3,
      ymax > 30 ~ ymax + 5,
      T ~ ymax
    )
  } else {ymax=ymax}

  p<-ggplot(plot_df, aes_string(x=x_axis_groups, y=gene, fill=x_axis_groups))+
    geom_boxplot(width=0.3)+
    geom_jitter(width=0.01,size=0.7)+
    ggtitle(paste0(gene," (",gene_descr,")"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10,face="bold"),
          plot.title = element_text(hjust=0.5, face="bold"),
          strip.text.x = element_text(size = 10,face="bold"),
          legend.position = "none") +
    ylab(y_lab)+
    xlab(x_lab)+
    ggsignif::geom_signif(
      comparisons = test_comparison,
      map_signif_level = F,
      step_increase = 0.12
    )+
    ylim(ymin[1],ymax[1])+
    stat_summary(fun.data = give.n, geom = "text", fun.args = list(ypos=ymin+0.8),size=6)

  # if user specified colors, add them
  if (!is.null(colors)){p<-p+scale_fill_manual(values=colors)}

  # if we do not have limma results, simply return the p value.
  if (is.null(limma_results)){return(p)}

  # if we did specify limma results check if gene is found.
  if (nrow(limma_results %>% dplyr::filter(symbol==gene))==0){
    print("Gene symbol was not found in limma result, returning plot with p value from wilcox.test")
    return(p)
  }


  # Do a work around to display the limma pvalue instead of the ggsignif pvalue.
  # deconstruct the plot
  p_deconstructed <- ggplot_build(p)

  # get adjusted p value from limma
  adj_p_vals <- limma_results %>%
    dplyr::filter(symbol==gene) %>%
    pull(padj) %>%
    signif(., digits=3)


  # now replace the ggsignif pvalue by the limma pvalue.
  # Something very similar can be used also if you have multiple panels.
  p_deconstructed$data[[3]] <- p_deconstructed$data[[3]]  %>%
    dplyr::mutate(annotation = case_when(
      PANEL == 1 ~ adj_p_vals[1]
    ))

  ## Reconstruct plot
  myplot3 <- ggpubr::as_ggplot(ggplot_gtable(myplot2))

  return(myplot3)
}

# function to indicate number of samples per boxplot
give.n <- function(x,ypos){
  return(data.frame(y = ypos, label = paste0("n=",length(x))))
}





#' Association between eigengene expression and phenotype
#' @description The function visualizes the assocation between a modules eigengene expression and a phenotype.
#' For continuous phenotypes (such as age) the correlation will be plotted. For categorical phenotypes (such as gender) a boxplot will be plotted
#'
#' @param mdata_eigengene A dataframe which contains both the module eigengene expression and the phenotype
#' for which you want to evaluate the association
#' @param phenotype Column of \code{mdata_eigengene} which contains the phenotype of interest.
#' Make sure that continuous phenotypes are of class \code{numeric} or \code{integer} and that
#' categorical variables are of class \code{factor} or \code{characters}.
#' @param module_eigengene Column of \code{mdata_eigengene} which contains the eigengene expression of the module of interest
#' @param xlab The x-axis title
#' @param ylab The y-axis title
#' @param label_x_pos X-axis position to write the calculated correlation coefficient. Only relevant for continuous phenotypes.
#' @param label_y_pos Y-axis position to write the calculated correlation coefficient. Only relevant for continuous phenotypes.
#' @param label_size Font size of the correlation coefficient. Only relevant for continuous phenotypes.
#' @param ymin Specify minimum of the Y-axis range.
#' @param ymax Specify maximum of the Y-axis range.
#' @param comparison Specify the type of statistical comparison you want to perform for categorical phenotypes.
#' For examples how to specify comparisons, check the \href{https://cran.r-project.org/web/packages/ggsignif/vignettes/intro.html}{homepage} of \code{ggsignif}
#' or the vignette. By default, all factors will be compared. If no comparisons should be performed, specify \code{NULL}.
#' @param colors Specify the colors of the boxplots. Number of colors must match number of factors in \code{phenotype}.
#' @examples
#' \dontrun{
#' library(ggpubr)
#' modules <- readRDS(module_save_path) # generated using GWENA::detect_modules()
#' mdata <- fread("path/to/file/with/metadata")
#' #add the eigengene expression to the metadata
#' mdata_eigengene <- modules$modules_eigengenes %>%
#'   rownames_to_column("sample.ID") %>%
#'   left_join(., mdata, by="sample.ID")
#'
#' module_association(mdata_eigengene = mdata_eigengene,
#' phenotype = "sample.status",
#' module_eigengene = "ME2",
#' ymin=-0.5, ymax=0.4,
#' xlab = "Sample status",
#' ylab= "Module 2 eigengene epxression",
#' colors=c("firebrick","steelblue"))
#' }
#'
#'
#'
#' @export
#'
module_association <- function(mdata_eigengene, phenotype, module_eigengene,
                               xlab=ggplot2::waiver(), ylab=ggplot2::waiver(),
                               label_x_pos=ggplot2::waiver(),
                               label_y_pos=ggplot2::waiver(), label_size=7,
                               ymin =-0.6, ymax=0.6,
                               comparison="all", colors=NULL){

  # determine if the phenotype of interest is continuous. If yes use correlation to assess association
  if (is.numeric(mdata_eigengene[phenotype][,1])){
    print("Assuming phenotype is continuous")
    print(paste(sum(is.na(mdata_eigengene[phenotype][,1])), "entries contain NA. Removing them"))

    mdata_eigengene <- mdata_eigengene %>%
      dplyr::filter(!is.na(get(phenotype)))

    p <- ggscatter(mdata_eigengene, x = phenotype, y = module_eigengene, add = "reg.line") +
      stat_cor(label.x=label_x_pos, label.y = label_y_pos, size = label_size)+
      xlab(xlab)+
      ylab(ylab)+
      theme(legend.position = "none",
            axis.text = element_text(size=16),
            axis.title = element_text(size=18))
    return(p)
  } else {

    # turn into factor
    mdata_eigengene[phenotype][,1] <- factor(mdata_eigengene[phenotype][,1])


    # for statistical comparison:
    if (is.null(comparison)){
      print("No statistical comparisons performed")
    } else if (comparison=="all"){
      # determine all the comparisons
      combinations <- t(combn(levels(mdata_eigengene[phenotype][,1]), 2))

      # split them into a list that is accepted by ggsignif
      comparison = split(combinations, seq(nrow(combinations)))

    }
    else {
      # check if the specified `test_comparisons` exist in the eset
      for (c in comparison %>% unlist()){
        if (!(c %in% mdata_eigengene[phenotype][,1])){
          stop(paste0("The specified `comparison` (", c,") does not occure in `phenotype`"))
        }
      }
    }

    p<-ggplot(mdata_eigengene %>%
                dplyr::filter(!is.na(get(phenotype))), #remove NAs
              aes_string(x = phenotype, y = module_eigengene, fill=phenotype))+
      geom_boxplot(width=0.4)+
      xlab(xlab)+
      ylab(ylab)+
      theme_bw()+
      theme(legend.position = "none",
            axis.text = element_text(size=16),
            axis.title = element_text(size=18))+
      ggsignif::geom_signif(comparisons = comparison,
                            map_signif_level = F,
                            test="wilcox.test",
                            step_increase = 0.1)+
      ylim(ymin,ymax) +
      stat_summary(fun.data = give.n, geom = "text", fun.args = list(ypos=ymin+0.1),size=6)

    if (!is.null(colors)){p<-p+scale_fill_manual(values=colors)}


    return(p)

  }
}










