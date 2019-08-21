# rsem.R
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

#' Internal method to extract expression data
#'
#' @param df data frame containing individual RSEM file output
#' @param scaling.method switch containing scaling options
#'
#' @return a scaled vector of counts
ExtractQuantitationData <- function(df, scaling.method) {
    if (scaling.method == "count") {
        return(df$expected_count)
    } else if (scaling.method == "TPM") {
        return(df$TPM)
    } else if (scaling.method == "lengthScaledTPM") {
        return(df$TPM / df$effective_length)
    }
}

#' Load Data from RSEM transcript quantitation.
#'
#' \code{LoadRSEM} loads rsem data for the provided cell to file mapping
#'
#' @param cells character vector containing cell IDs
#' @param rsem.filenames character vector containing RSEM file paths.
#' @param quantitation.method String, default 'count', returns expected
#'     counts from RSEM files. 'TPM' also supported.
#' 
#' @return A matrix containing the TPM counts of each ENSG merged from
#'     all cells supplied in a Seurat-importable format.
#'
#' @examples
#'
#' \dontrun{
#' LoadRSEM(c('1', '2', '3'), c('1.rsem', '2.rsem', '3.rsem'))
#' }
#'
#' @importFrom hashmap hashmap
#' @importFrom utils read.table
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom Seurat CreateAssayObject
#'
#' @rdname LoadRSEM
#' @export LoadRSEM
#'
LoadRSEM <- function(cells, rsem.filenames, quantitation.method='TPM') {
    # create progress bar
    pb <- txtProgressBar(min = 0,
        label = "Reading Expression from RSEM Files",
        max = length(cells),
        style = 3)

    nonzero_cell_maps <- c()

    for (i in 1:length(cells)) {
        temp_df <- read.table(rsem.filenames[[i]], header = TRUE, sep = "\t")

        temp_df$scaled <- ExtractQuantitationData(temp_df, quantitation.method)
        temp_df <- temp_df[temp_df$scaled > 0,]

        temp_map <- hashmap(
            as.character(temp_df$gene_id),
            temp_df$scaled
        )

        nonzero_cell_maps <- c(nonzero_cell_maps, temp_map)

        setTxtProgressBar(pb, i)
    }

    nonzero_genes <- unique(
        c(unlist(sapply(nonzero_cell_maps, function(a) a$keys())))
    )

    M <- matrix(0, nrow = length(nonzero_genes), ncol = length(cells))
    rownames(M) <- nonzero_genes
    colnames(M) <- cells

    for (i in 1:length(cells)) {
        colvals <- nonzero_cell_maps[[i]][[nonzero_genes]]
        colvals[is.na(colvals)] <- 0
        M[,i] = colvals
    }

    close(pb)

    rsem.melange <- Melange(M)

    return(rsem.melange)
}
