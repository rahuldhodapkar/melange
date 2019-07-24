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
generate.varid.from.annovar.vff <- function(df) {
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
#'
#' @return A hashmap from gene name -> variant counts
#'
group_by_gene <- function(df) {
    GENE_COL <- 2
    df$clean.gene.col <- gsub("[\\(\\,\\;].*","", df[,GENE_COL])
    freq.table <- as.data.frame(table(df$clean.gene.col))
    return(hashmap(as.character(freq.table[,1]),freq.table[,2]))
}

#' Load data from annovar-annotated variant calls.
#'
#' \code{load.annovar} loads Annovar data for the provided cell to file mapping.
#' By default, does not perfom any quality control on variant call files.
#' \code{load.annovar} should be provided \code{*.variant_function} files
#' for each cell. If additional exonic functional pruning is required,
#' \code{load.annovar} will automatically extract the location of the corresponding
#' \code{*.exonic_variant_function} file locally.
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
#' 
#' @return A matrix containing the reduced counts of each ENSG merged from
#'     all cells supplied in a Seurat-importable format.
#'
#' @examples
#' \dontrun{
#' 
#' load.annovar(c('1', '2'), 
#'     c('1.anno.variant_function', '2.anno.variant_function')
#' )
#' 
#' }
#'
#' @importFrom hashmap hashmap
#' @importFrom utils read.table
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom methods new
#' @importClassesFrom Seurat Assay
#'
#' @rdname load.annovar
#' @export load.annovar
#'
load.annovar <- function(cells, annovar.filenames, germline.filename,
                            reduction="by_gene", somatic.only=TRUE,
                            spike.in.regex="^ERCC-") {

    if (! is.null(germline.filename) ) {
        germline_df <- read.table(germline.filename, 
                                    comment.char = '', header = FALSE, sep="\t")
        germline_varstrings <- generate.varid.from.annovar.vff(germline_df)
    } else {
        germline_varstrings = c()
    }
    germ.map <- hashmap(
        germline_varstrings, rep(1, length(germline_varstrings))
    )

    pb <- txtProgressBar(min = 0,
        max = length(cells),
        style = 3)

    grouped_count_maps <- c()
    germline_counts <- c()
    total_counts <- c()
    spike_in_counts <- c()
    somatic_counts <- c()

    for (i in 1:length(cells)) {
        temp_df <- read.table(annovar.filenames[[i]], 
                                comment.char = '', header = FALSE, sep = "\t")
        temp_varstrings <- generate.varid.from.annovar.vff(temp_df)

        germline_hits <- ! is.na(germ.map[[temp_varstrings]])
        num_germline_hits <- sum(germline_hits, na.rm = TRUE)
        num_spike_in_hits <- sum(grepl(spike.in.regex, temp_varstrings), na.rm=TRUE)
        num_total_hits <- nrow(temp_df)
        num_somatic_hits <- sum(is.na(germ.map[[temp_varstrings]]), na.rm = TRUE)

        if (somatic.only) {
            temp_df <- temp_df[is.na(germ.map[[temp_varstrings]]),]
        }

        germline_counts <- c(germline_counts, num_germline_hits)
        total_counts <- c(total_counts, num_total_hits)
        spike_in_counts <- c(spike_in_counts, num_spike_in_hits)
        somatic_counts <- c(somatic_counts, num_somatic_hits)

        temp_counts_map <- group_by_gene(temp_df)
        grouped_count_maps <- c(grouped_count_maps, temp_counts_map)

        setTxtProgressBar(pb, i)
    }

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

    for (i in 1:length(cells)) {
        colvals <- grouped_count_maps[[i]][[nonzero_keys]]
        colvals[is.na(colvals)] <- 0
        M[,i] = colvals
    }
    close(pb)

    annovar_assay <- new(
        Class = "Assay",
        counts = M,
        data = M,
        scale.data = M,
        key = 'annovar_',
        misc = list(
            meta.data = data.frame(row.names = colnames(x = M))
        )
    )

    annovar_assay@misc$meta.data$germline.calls <- cell2germ_ct[[cells]]
    annovar_assay@misc$meta.data$total.calls <- cell2total_ct[[cells]]
    annovar_assay@misc$meta.data$spike.in.calls <- cell2spike_ct[[cells]]
    annovar_assay@misc$meta.data$somatic.calls <- cell2som_ct[[cells]]

    return(annovar_assay)
}




