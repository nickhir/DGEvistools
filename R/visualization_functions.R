#####
# Plotting themes
#####
theme_volcano <- function () {
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

# this theme is taken from https://github.com/koundy/geom_roc_plot/blob/master/Theme_Publication.R
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}


#' @title Plot a pretty volcano plot
#' @description This function plots a pretty volcano plot using the results of the output of \code{DESeq2::results()}. Code was adapted from [here](https://github.com/jackhump/ALS_SpinalCord_QTLs).
#'
#' @param de_res A dataframe created using the \code{msdbulkseqtools::toptable_for_all} command. Alternatively a data frame
#' containing the following columns: 'padj', 'log2FoldChange'
#' @param title The title of the plot
#' @param subtitle The subtitle of the plot
#' @param annotate_by A vector with gene names which will be annotated in the volcano plot.
#' Gene names have to occur in the 'symbol' column of 'de_res'
#' @param annotation_size The font size of the annotation.
#' @param padj_threshold Threshold determining which genes should be counted as differentially expressed
#' @param logFC_threshold Threshold determining which genes should be counted as "strongly" differentially expressed
#' @param ymax Maximum of the y axis
#' @param xlim The x axis limits.
#' @param res Resolution of the plot. The higher the resolution, the longer takes the plot.
#' @param min_pval_cutoff All p values that are smaller than \code{min_pval_cutoff} will be set to \code{min_pval_cutoff} for the visualization
#' @return A ggplot of a volcano plot.
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import tibble
#' @import SummarizedExperiment
#' @import grid
#' @import ggthemes
#'
#' @export
volcano_plot <- function(de_res, title = NULL, subtitle = NULL, annotate_by = NULL, annotation_size=5.2,
                         padj_threshold=0.05, logFC_threshold=1, ymax=NULL, xlim=c(-3,3), res=300,
                         min_pval_cutoff=1e-16){

  # check if all required columns exist.
  if (!all(c("padj","log2FoldChange") %in% colnames(de_res))){
    stop("ERROR: `de_res` requires the following columns: `padj`,`log2FoldChange`")
  }

  # add some extra annotation about the magnitude of differential expression using log2FoldChange
  de_res <- de_res %>%
    dplyr::filter(dplyr::if_all(c(padj, log2FoldChange), ~ !is.na(.))) %>%
    dplyr::mutate(sig = case_when(
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
  de_tally <- de_res %>%
    group_by(sig, direction, class) %>% tally() %>%
    dplyr::filter(sig != "non_sig") %>%
    dplyr::mutate(position = ifelse(sig == "sig", 0.5, 2) ) %>%
    dplyr::mutate(position = ifelse( direction == "down", -1 * position, position)) %>%
    dplyr::mutate(n = formatC(n, format="f", big.mark=",", digits=0))

  # determine a suitable ymax if none was specified
  if (is.null(ymax)){
    ymax <- ifelse(min(de_res$padj, na.rm=T) < min_pval_cutoff,
                   min_pval_cutoff+0.05,
                   -log10(min(de_res$padj, na.rm=T))+2)
  }

  # generate the volcano plot
  plot <- de_res %>%
    dplyr::mutate(padj = ifelse(padj < min_pval_cutoff, min_pval_cutoff, padj)) %>% #threshold at 1e16
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
    theme_volcano() +
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
    if (! "symbol" %in% colnames(de_res)){
      stop("Dataframe doesnt contain a `symbol` column")
    }
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
        data = dplyr::filter(de_res, symbol %in% annotate_by), size = 2.5, colour = "black"
      ) +
      geom_point(aes(colour = class ),
                 data = dplyr::filter(de_res, symbol %in% annotate_by), size = 1.3
      )

  }
  return(plot)

}





#' Plot correlation between two differential gene expression analyses
#' @description This function visualizes the correlation between two different differential gene expression analyses.
#' Might be used to compare two different tissues or two different timepoints and see how similar they behave.

#' @param result_df1  The first results dataframe. Usually created using the \code{DESeq2::results()} command.
#' @param result_df2 The second results dataframe. Usually created using the \code{DESeq2::results()} command.
#' @param title_result_1 Name of the first differential gene expression analysis. Will be the title on the x-axis.
#' @param title_result_2  Name of the second differential gene expression analysis. Will be the title on the y-axis.
#' @param col Name of the column which will be used to perform the correlation analysis.
#' @param remove_outlier If true, only values between the 1st and 99th percentile are included in analysis.
#' For example \code{log2FoldChange} Has to occur in both \code{result_df1} and \code{result_df2}

#' @export

plot_de_correlation <- function(result_df1, result_df2,
                                title_result_1, title_result_2,
                                col = "log2FoldChange",
                                remove_outlier=T){
  # combine dataframes
  join_col_1 <- colnames(result_df1)[1]
  join_col_2 <- colnames(result_df2)[1]
  if (join_col_1 != join_col_2){
    stop("The first column name of result_df1 and result_df2 should contain the same information and column name.")
  }

  res <- dplyr::left_join(result_df1,
                          result_df2,
                          by = join_col_1 ,
                          suffix = c(".1", ".2"))

  # when joining we add a suffix, so we also have to do it here
  x_string <- paste0(col, ".1")
  y_string <- paste0(col, ".2")

  # sometimes we have extreme outliers which distort the results.
  # we only keep values between the 1st and 99th percentile
  if (remove_outlier){
    before = nrow(res)
    Q_x <- quantile(res[,x_string], probs=c(0.01, 0.99))
    Q_y <- quantile(res[,y_string], probs=c(0.01, 0.99))

    res <- res %>%
      dplyr::filter(get(x_string)>Q_x[1]) %>%
      dplyr::filter(get(x_string)<Q_x[2]) %>%
      dplyr::filter(get(y_string)>Q_y[1]) %>%
      dplyr::filter(get(y_string)<Q_y[2])
    removed = before - nrow(res)
    message(paste(removed, "entries were removed, because we are only plotting values between the 1st and 99th percentile. Specify `remove_outlier=F` to disable."))
  }




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
    theme_Publication() %+%
    theme(axis.title = element_text(size=15), axis.text = element_text(size=12))

  return(plot)
}



#' Generate a heatmap of GSEA results
#' @description This function can be used to compare GSEA results performed for different tissues, different timepoints or conditions.
#' It will generate a heatmap, where rows represent different gene sets and columns the different comparisons.
#'
#' @param gsea_df A dataframe which contains the results of \code{clusterProfiler::GSEA}.
#' @param comparisons_column The name of the column which contains the different comparisons that you want to make.
#' @param ordering By which condition should the p values of the displayed gene sets be ranked
#' @param top_n The number of gene sets you want to display.
#' @param color A colormap for the heatmap
#' @param legend_title The title of the legend
#' @param title The title of the heatmap
#' @param remove_string A string that should be removed from all rownames of the heatmap.
#' @param ... Arguments are directly passed to \code{ComplexHeatmap::Heatmap}
#'
#' @export
enrichment_heatmap <- function(gsea_df, comparisons_column,
                               ordering = NULL,top_n=10,
                               color=viridis::plasma(50,end=0.75),
                               legend_title="-log10(FDR)",
                               title="",
                               remove_string=NULL, ...){


  # this gets called by `cell_fun` of ComplexHeatmap and returns asteriscs according to p value.
  annotate_sig <- function(j, i, x, y, w, h, fill) {
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


  comparisons <- unique(gsea_df[[comparisons_column]])

  if (is.null(ordering)){
    ordering <- comparisons[1]
  }

  # get the top n enriched terms for a selected tissue
  ontologies_of_interest <- gsea_df %>%
    dplyr::filter(get(comparisons_column)==ordering)%>%
    dplyr::arrange(p.adjust) %>%
    slice_head(n=top_n) %>%
    pull(ID)

  # get the p values of the selected ontologies for all other comparisons.
  wide_df <- gsea_df %>%
    dplyr::filter(ID %in% ontologies_of_interest) %>%
    dplyr::mutate(neglog10FDR=-log10(p.adjust))%>%
    dplyr::select(c("ID","Description", "neglog10FDR",comparisons_column)) %>%
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
    as.matrix

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
#' @param mdata A dataframe with information about the samples
#' @param info_col A column in the \code{mdata} object, which has entries that match the column names of \code{scaled_counts}.
#' They have to be in the same order, otherwise an error will occure.
#' @param column_cluster A named vector of subjects and corresponding cluster. Use only if you already determined how to cluster patients.
#' If none is specified, \code{hclust} is used to perform the clustering.
#' @param GWENA_modules Absolute path to the identified co-expression modules. Generated using \code{GWENA::detect_modules()} followed by \code{saveRDS()}.
#' The whole \code{GWENA} workflow has to be performed using the same subjects and genes that are in the \code{scaled_counts}.
#' @param heatmap_annotation A \href{https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html}{HeatmapAnnotation object}.
#' @param ... Other arguments that are passed directly to the \code{ComplexHeatmap::Heatmap} function. Might be something like \code{show_row_names=T}
#' @return A \code{ComplexHeatmap} which shows the GWENA results.
#'
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
#' @description Given a properly formatted \code{SummarizedExperiment.}, this function can quickly plot boxplots to investigate the expression of specific genes between two conditions.
#' By default, CPM transformed counts are expected
#'
#' @param gene The gene symbol of the gene for which the expression should be plotted. Must occur in the \code{symbol} column of \code{rowData()} of the \code{SummarizedExperiment}.
#' @param SE A SummarizedExperiment. Usually the same which was used to generate the DESEq2 object.
#' @param intgroup Counts are grouped by this variable. Has to occur in \code{colData()} of the \code{SummarizedExperiment}.
#' @param test_comparison Optionally, you can also specify if and how conditions should be compared. Should be a list which contains the comparisons that you are interested in.
#' See examples and \href{https://cran.r-project.org/web/packages/ggsignif/vignettes/intro.html}{here}.
#' @param DESEq_res By default, if you specify \code{test_comparison} a \code{wilcox.test} is performed.
#'  This might be fine to get an intuition if expression differences between two conditions
#'  are significantly different, but the proper way is to use the the p value that
#'  was determined using \code{DESeq}. For that, you can simply specify an object that was created \code{DESeq2::results()} that contains the correct contrast/comparison.
#'  Also has to contain column \code{symbol}.
#' @param colors Specify the colors of the boxplots. Number of colors must match number of factors in \code{intgroup}.
#' @param x_lab X axis label. Defaults to \code{intgroup}
#' @param y_lab Y axis label. Defaults to \code{log(CPM)}
#' @param ymax Max value of the Y axis. Sometimes useful to pick a different value than the default.
#' @param ymin Min value of the Y axis. Sometimes useful to pick a different value than the default.
#'
#' @return A ggplot of a boxplot
#'
#' @examples
#' \dontrun{
#' See vignette `rna_seq_workflow`
#' }
#'
#' @export
#'
create_expression_boxplot <- function(gene, SE, intgroup,
                                      test_comparison=NULL, DESEq_res=NULL,
                                      colors=NULL, x_lab=ggplot2::waiver(),
                                      y_lab=expression(log[2](CPM)),
                                      ymax=NULL, ymin=NULL){

  # subset the SummarizedExperiment to only include the gene of interest.

  keep <- rowData(SE)$symbol == gene
  keep[is.na(keep)] <- F
  if (sum(keep) == 0){
    print("Gene symbol was not found in SummarizedExperiment, returning empty plot")
    return()
  }

  # check if the specified `test_comparisons` exist in the SummarizedExperiment
  for (comparison in test_comparison %>% unlist()){
    if (!(comparison %in% colData(SE)[intgroup][,1])){
      stop(paste0("The specified `test_comparison` (", comparison,") does not occure in `intgroup`"))
    }
  }

  # only keep information for the gene we are interested in.
  subset_SE <- SE[keep,]

  # turn SE into a long df. SE was basically required to perform effective selection of gene of interest.
  if (!"cpm_counts" %in% names(assays(SE))){
    stop("SummarizedExperiment should contain an assay called `cpm_counts`, which contains log2 CPM counts.")
  }

  counts <- assays(subset_SE)$cpm_counts %>%
    as.matrix() %>%
    t() %>%
    magrittr::set_colnames(gene)

  # create a long df which contains information about the expression as well as additional information such as description and meta data.
  plot_df <- counts %>%
    data.frame() %>%
    rownames_to_column("ID") %>%
    left_join(.,colData(subset_SE) %>% data.frame %>% rownames_to_column("ID"), by="ID")


  # in very few edgecases:
  if (n_distinct(colnames(plot_df))!=ncol(plot_df)){
    colnames(plot_df) <- make.names(colnames(plot_df),unique = T)
  }

  # get gene description for the plot
  if ("description" %in% colnames(rowData(subset_SE))){
    gene_descr <- rowData(subset_SE) %>%
      data.frame() %>%
      pull("description") %>%
      str_wrap(30)
  } else {
    gene_descr=""
  }



  # make sure that the `intgroup` are factors.
  plot_df[,intgroup] <- factor(plot_df[,intgroup])

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

  p<-ggplot(plot_df, aes_string(x=intgroup, y=gene, fill=intgroup))+
    geom_boxplot(width=0.3)+
    geom_jitter(width=0.01,size=0.7)+
    ggtitle(paste0(gene," (",gene_descr,")"))+
    theme_Publication()+
    ylab(y_lab)+
    xlab(x_lab)+
    ylim(ymin[1],ymax[1])+
    stat_summary(fun.data = give.n, geom = "text", fun.args = list(ypos=ymin+0.3),size=6)

  if (!is.null(test_comparison)){
    p <- p +
      ggsignif::geom_signif(
      comparisons = test_comparison,
      map_signif_level = F,
      step_increase = 0.12)
  }

  # if user specified colors, add them
  if (!is.null(colors)){p<-p+scale_fill_manual(values=colors)}

  # if we do not have DESeq results, simply return the p value.
  if (is.null(DESEq_res)){return(p)}

  # if we did specify DESeq results check if gene is found.
  if (nrow(DESEq_res %>% dplyr::filter(symbol==gene))==0){
    print("Gene symbol was not found in DESeq result, returning plot with p value from wilcox.test")
    return(p)
  }

  # Do a work around to display the DESeq pvalue instead of the ggsignif pvalue.
  # deconstruct the plot
  p_deconstructed <- ggplot_build(p)

  # get adjusted p value from DESeq
  adj_p_vals <- DESEq_res %>%
    dplyr::filter(symbol==gene) %>%
    pull(padj) %>%
    signif(., digits=3)


  # now replace the ggsignif pvalue by the deseq2 pvalue.
  # Something very similar can be used also if you have multiple panels.
  if ("annotation" %in% colnames( p_deconstructed$data[[4]])){
    p_deconstructed$data[[4]] <- p_deconstructed$data[[4]]  %>%
      dplyr::mutate(annotation = case_when(
        PANEL == 1 ~ adj_p_vals[1]
      ))
  } else if ("annotation" %in% colnames( p_deconstructed$data[[3]])){
    p_deconstructed$data[[3]] <- p_deconstructed$data[[3]]  %>%
      dplyr::mutate(annotation = case_when(
        PANEL == 1 ~ adj_p_vals[1]
      ))
  } else {
    message("Make sure you have ggplot version 3.3.6 installed. Returning wilcoxon p value.")
  }

  ## Reconstruct plot
  final_plot <- ggpubr::as_ggplot(ggplot_gtable(p_deconstructed))

  return(final_plot)
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
      theme_Publication()
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
      theme_Publication()+
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










