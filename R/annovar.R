# annovar.R
#     Copyright (C) 2019  Rahul Dhodapkar
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.



#' Helper function to concatenate the unique identifiers for a variant
#'
#' @param df dataframe read from annovar output file
#'
#' @return A vector containing all unique variant identifiers present
#'      in an annovar variant file.
#'
GenerateVaridFromAnnovarVFF <- function(df) {
    colname2ix <- hashmap(
        c("CHR", "POS_START", "POS_END", "REF", "ALT"),
        c(    3,           4,         5,     6,     7)
    )

    return(paste(
            as.character(df[,colname2ix[["CHR"]]]),
            as.character(df[,colname2ix[["POS_START"]]]),
            as.character(df[,colname2ix[["POS_END"]]]),
            as.character(df[,colname2ix[["REF"]]]),
            as.character(df[,colname2ix[["ALT"]]]),
            sep="."))
}

#' Helper function to group variants by gene.
#'
#' @param df dataframe read from annovar output file
#' @param read.depth.norm.col integer, column to use for normalization of variant
#'      call rates.
#'
#' @return A hashmap from gene name -> variant counts
#'
#' @importFrom stats aggregate
#'
GroupByGene <- function(df, read.depth.norm.col) {
    GENE_COL <- 2
    df$clean.gene.col <- gsub("[\\(\\,\\;].*","", df[,GENE_COL])
    norm.data <- aggregate(df[,read.depth.norm.col],
        by=list(Gene=df$clean.gene.col), FUN=sum)
    return(list(
        norm=hashmap(as.character(norm.data[,1]),norm.data[,2])
    ))
}

#' Load data from annovar-annotated variant calls.
#'
#' \code{LoadAnnovar} loads Annovar data for the provided cell to file mapping.
#' By default, does not perfom any quality control on variant call files.
#' \code{LoadAnnovar} should be provided \code{*.variant_function} files
#' for each cell. If additional exonic functional pruning is required,
#' \code{LoadAnnovar} will automatically extract the location of the corresponding
#' \code{*.exonic_variant_function} file locally.
#'
#' If a local \{*.exonic_variant_function} file is not found, \code{LoadAnnovar}
#' will fail to load a data matrix.
#'
#' @param cells character vector containing cell IDs
#' @param annovar.filenames character vector containing annovar file paths.
#' @param reduction reduction to use from annovar to Seurat assay.
#'     Currently 'by_gene'
#'
#' @param germline.filename file path for annovar file containing germline 
#'     variants. Used as a positive control for single-cell somatic 
#'     variant analysis.
#' @param spike.in.regex regex selection for 'CHR' vcf field inidcating
#'     synthetic spike-in. Used as a negative control for single-cell
#'     somatic variant analysis. Default '^ERCC-'
#' @param somatic.only boolean, default TRUE, produce counts based only
#'     on predicted somatic variants (subtact germline)
#' @param exonic.only boolean, default TRUE, produce counts based only
#'     on variants predicted to occur within exons.
#' @param read.depth.norm.col default 10, col with read depth to normalize
#'     variant calls by for single cell analysis. If using the default
#'     Annovar conversion tools from VCF, this will most likely be column 10.
#' @param unique.id.cols set of columns that uniquely define a variant in the
#'     *.variant_function file. Will be used to generate row names for the
#'     final aggregated data set.
#'
#' @return A matrix containing the reduced counts of each ENSG merged from
#'     all cells supplied in a Seurat-importable format.
#'
#' @examples
#' \dontrun{
#' 
#' LoadAnnovar(c('1', '2'),
#'     c('1.anno.variant_function', '2.anno.variant_function')
#' )
#' 
#' }
#'
#' @importFrom hashmap hashmap
#' @importFrom utils read.table
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stringr str_detect
#' @importFrom Seurat CreateAssayObject
#'
#' @rdname LoadAnnovar
#' @export LoadAnnovar
#'
LoadAnnovar <- function(cells, annovar.filenames, germline.filename,
                            reduction="by_gene", somatic.only=TRUE,
                            exonic.only=TRUE, read.depth.norm.col=10,
                            spike.in.regex="^ERCC-", unique.id.cols=c(3,4,5,6,7)) {

    annovar.gene.col <- 2
    annovar.exonic.variant.function.filenames <- 
        gsub("\\.variant_function$", "\\.exonic_variant_function", annovar.filenames)
    exonic.variant.function.join.col <- 1
    exonic.variant.mutation.type.col <- 2

    if (! is.null(germline.filename) ) {
        germline_df <- read.table(germline.filename, 
                                    comment.char = '', header = FALSE, sep="\t")
        germline_varstrings <- GenerateVaridFromAnnovarVFF(germline_df)
    } else {
        germline_varstrings = c()
    }
    germ.map <- hashmap(
        germline_varstrings, rep(1, length(germline_varstrings))
    )

    pb <- txtProgressBar(
        label = "Loading Annovar data files",
        min = 0,
        max = length(cells),
        style = 3)

    grouped_count_maps <- c()
    germline_counts <- c()
    total_counts <- c()
    spike_in_counts <- c()
    somatic_counts <- c()

    varstring2mutation.type <- hashmap(c(""), c(""))
    varstring2gene <- hashmap(c(""), c(""))

    for (i in 1:length(cells)) {
        if (file.size(annovar.filenames[[i]]) == 0) {
            grouped_count_maps <- c(grouped_count_maps, hashmap())
            next
        }

        temp_df <- read.table(annovar.filenames[[i]], 
                                comment.char = '', header = FALSE, sep = "\t")

        if (file.size(annovar.exonic.variant.function.filenames[[i]]) == 0) {
            # to allow code to continue gracefully
            temp_df$MutationType <- rep(NA, nrow(temp_df))
        }
        else {
            temp_exon_anno_df <- read.table(annovar.exonic.variant.function.filenames[[i]],
                              comment.char = '', header = FALSE, sep = "\t")
            temp_exon_anno_df[,exonic.variant.function.join.col] <- 
                gsub("line", "", temp_exon_anno_df[,exonic.variant.function.join.col])

            # add Mutation Type for exonic variants
            temp_df[,"MutationType"] <- rep(NA, nrow(temp_df))
            temp_df[as.numeric(temp_exon_anno_df[,exonic.variant.function.join.col]), "MutationType"] <-
                as.character(temp_exon_anno_df[,exonic.variant.mutation.type.col])
        }

        temp_df$temp_varstrings <- GenerateVaridFromAnnovarVFF(temp_df)

        germline_hits <- ! is.na(germ.map[[temp_df$temp_varstrings]])
        num_germline_hits <- sum(germline_hits, na.rm = TRUE)
        num_spike_in_hits <- sum(grepl(spike.in.regex, temp_df$temp_varstrings), na.rm=TRUE)
        num_total_hits <- nrow(temp_df)
        num_somatic_hits <- sum(is.na(germ.map[[temp_df$temp_varstrings]]), na.rm = TRUE)

        if (somatic.only) {
            temp_df <- temp_df[is.na(germ.map[[temp_df$temp_varstrings]]),]
        }

        germline_counts <- c(germline_counts, num_germline_hits)
        total_counts <- c(total_counts, num_total_hits)
        spike_in_counts <- c(spike_in_counts, num_spike_in_hits)
        somatic_counts <- c(somatic_counts, num_somatic_hits)

        if (exonic.only) {
            temp_df <- temp_df[str_detect(temp_df[,1], "exonic"),]
        }

        temp_counts_map <- hashmap(
            temp_df$temp_varstrings, 
            temp_df[,read.depth.norm.col])

        grouped_count_maps <- c(grouped_count_maps, temp_counts_map)

        exonic_df <- temp_df[!is.na(temp_df$MutationType),]
        varstring2mutation.type$insert(exonic_df$temp_varstrings, as.character(exonic_df$MutationType))
        varstring2gene$insert(temp_df$temp_varstrings, temp_df[,annovar.gene.col])

        setTxtProgressBar(pb, i)
    }
    close(pb)

    cell2germ_ct <- hashmap(cells, germline_counts)
    cell2total_ct <- hashmap(cells, total_counts)
    cell2spike_ct <- hashmap(cells, spike_in_counts)
    cell2som_ct <- hashmap(cells, somatic_counts)

    nonzero_keys <- unique(
        c(unlist(sapply(grouped_count_maps, function(a) a$keys())))
    )

    M <- matrix(0, nrow = length(nonzero_keys), ncol = length(cells))
    rownames(M) <- nonzero_keys
    colnames(M) <- cells


    print("Collating and generating variation matrix")
    pb <- txtProgressBar(min = 0,
        label = "Collating and generating variation matrix",
        max = length(cells),
        style = 3)

    for (i in 1:length(cells)) {
        colvals <- grouped_count_maps[[i]][[nonzero_keys]]
        colvals[is.na(colvals)] <- 0
        M[,i] = colvals

        setTxtProgressBar(pb, i)
    }
    close(pb)

    annovar.melange <- Melange(M)
    annovar.melange@meta.data$germline.calls <- cell2germ_ct[[cells]]
    annovar.melange@meta.data$total.calls <- cell2total_ct[[cells]]
    annovar.melange@meta.data$spike.in.calls <- cell2spike_ct[[cells]]
    annovar.melange@meta.data$somatic.calls <- cell2som_ct[[cells]]

    annovar.melange@meta.features$MutationType <- varstring2mutation.type[[rownames(M)]]
    annovar.melange@meta.features$Gene <- varstring2gene[[rownames(M)]]

    return(annovar.melange)
}




