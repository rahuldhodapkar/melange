# mixcr.R
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

#' Load data from MiXCR output files
#'
#' Given a particular TCR chain file set (e.g. TCRB / TCRA), load data into
#' a confidence score matrix where nucleotide sequence of CDR3 is the primary key.
#'
#' @param cells character vector containing cell IDs
#' @param mixcr.filenames character vector containing MiXCR file paths.
#' 
#' @return a melange object with all called TCR sequences for the provided cells
#'
#' @importFrom hashmap hashmap
#' @importFrom utils read.table
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @rdname LoadMiXCR
#' @export LoadMiXCR
#' 
LoadMiXCR <- function(cells, mixcr.filenames) {
    # create progress bar
    pb <- txtProgressBar(min = 0,
        label = "Reading TCR from MiXCR files",
        max = length(cells),
        style = 3)

    n2aa <- hashmap(character(), character())
    n2v <- hashmap(character(), character())
    n2d <- hashmap(character(), character())
    n2j <- hashmap(character(), character())

    nonzero_cell_maps <- c()

    for (i in 1:length(cells)) {
        temp_df <- read.table(mixcr.filenames[[i]], header = TRUE, sep = "\t")

        temp_map <- hashmap(
            as.character(temp_df$nSeqCDR3),
            temp_df$cloneFraction
        )

        n2aa[[temp_df$nSeqCDR3]] <- temp_df$aaSeqCDR3
        n2v[[temp_df$nSeqCDR3]] <- temp_df$allVHitsWithScore
        n2d[[temp_df$nSeqCDR3]] <- temp_df$allDHitsWithScore
        n2j[[temp_df$nSeqCDR3]] <- temp_df$allJHitsWithScore

        nonzero_cell_maps <- c(nonzero_cell_maps, temp_map)

        setTxtProgressBar(pb, i)
    }

    nonzero_clones <- unique(
        c(unlist(sapply(nonzero_cell_maps, function(a) a$keys())))
    )

    M <- matrix(0, nrow = length(nonzero_clones), ncol = length(cells))
    rownames(M) <- nonzero_clones
    colnames(M) <- cells

    for (i in 1:length(cells)) {
        colvals <- nonzero_cell_maps[[i]][[nonzero_clones]]
        colvals[is.na(colvals)] <- 0
        M[,i] = colvals
    }

    close(pb)

    mixcr.melange <- Melange(M)
    mixcr.melange@meta.features$aaSeqCDR3 <- n2aa[[rownames(M)]]
    mixcr.melange@meta.features$allVHitsWithScore <- n2v[[rownames(M)]]
    mixcr.melange@meta.features$allDHitsWithScore <- n2d[[rownames(M)]]
    mixcr.melange@meta.features$allJHitsWithScore <- n2j[[rownames(M)]]

    return(mixcr.melange)         
}