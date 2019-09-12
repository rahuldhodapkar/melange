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
#'
ExtractQuantitationData <- function(df, scaling.method) {
    if (scaling.method == "count") {
        return(df$expected_count)
    } else if (scaling.method == "TPM") {
        return(df$TPM)
    }
}

#' Load Data from RSEM transcript quantitation.
#'
#' \code{LoadRSEM} loads rsem data for the provided cell to file mapping
#'
#' @param cells character vector containing cell IDs
#' @param rsem.filenames character vector containing RSEM file paths.
#' @param quantitation.method String, default 'TPM', returns expected
#'     counts from RSEM files. 'count' also supported.
#' @param min.quant.value numeric, containing the minimum quantitation
#'     a gene must have to be included in the output \code{@data} matrix
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
LoadRSEM <- function(cells, rsem.filenames,
                            quantitation.method='TPM',
                            min.quant.value=0) {
    # create progress bar
    print("Reading Expression from RSEM Files")
    pb <- txtProgressBar(min = 0,
        label = "Reading Expression from RSEM Files",
        max = length(cells),
        style = 3)

    above.threshold.cell.maps <- c()

    for (i in 1:length(cells)) {
        temp.df <- read.table(rsem.filenames[[i]], header = TRUE, sep = "\t")

        temp.df$scaled <- ExtractQuantitationData(temp.df, quantitation.method)
        temp.df <- temp.df[temp.df$scaled > min.quant.value,]

        temp.map <- hashmap(
            as.character(temp.df$gene_id),
            temp.df$scaled
        )

        above.threshold.cell.maps <- c(above.threshold.cell.maps, temp.map)

        setTxtProgressBar(pb, i)
    }

    close(pb)

    print("Collating and generating expression matrix")
    pb <- txtProgressBar(min = 0,
        label = "Collating and generating expression matrix",
        max = length(cells),
        style = 3)

    above.threshold.genes <- unique(
        c(unlist(sapply(above.threshold.cell.maps, function(a) a$keys())))
    )

    M <- matrix(0, nrow = length(above.threshold.genes), ncol = length(cells))
    rownames(M) <- above.threshold.genes
    colnames(M) <- cells

    for (i in 1:length(cells)) {
        setTxtProgressBar(pb, i)
        colvals <- above.threshold.cell.maps[[i]][[above.threshold.genes]]
        colvals[is.na(colvals)] <- 0
        M[,i] = colvals
    }
    close(pb)

    rsem.melange <- Melange(M)

    return(rsem.melange)
}
