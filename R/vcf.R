# vcf.R
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

#' Load data from VCF output files
#'
#' Combine variants by position only, to use as input for describing positions of
#' mutational frequency.
#'
#' @param cells character vector containing cell IDs
#' @param vcf.filenames character vector containing VCF file paths
#' @param germline.filename single filename for germline vcf
#' @param consolidation.factor the factor to use for dividing 
#' 
#' @return a melange object with consolidated variants by consolidation.factor
#'         rows will be named with the form: _<chr_name>:<position>
#'
#' @importFrom hashmap hashmap
#' @importFrom utils read.table
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @rdname LoadVCF
#' @export LoadVCF
#' 
LoadVCF <- function(cells, vcf.filenames, germline.filename, consolidation.factor=10000) {

    VCF.CHROMOSOME.COL <- 1
    VCF.POS.COL <- 2

    empty.map <- hashmap(c(""), c(0))
    empty.map$clear()

    print("Reading germline VCF file")
    germline.df <- read.table(germline.filename, sep='\t', comment.char='#')
    code2germhits <- hashmap(
        paste0("_",
            germline.df[,VCF.CHROMOSOME.COL],
             ":",
            germline.df[,VCF.POS.COL]),
        rep(1, nrow(germline.df))
    )

    # create progress bar
    print("Reading variants from VCF files")
    pb <- txtProgressBar(min = 0,
        label = "Reading variants from VCF files",
        max = length(cells),
        style = 3)

    nonzero.cell.maps <- c()
    for (i in 1:length(cells)) {
        
        temp.df <- try(read.table(vcf.filenames[[i]], sep='\t', comment.char='#'))
        if (inherits(temp.df, 'try-error')) {
            nonzero.cell.maps <- c(nonzero.cell.maps, empty.map)
            next
        }

        # isolate somatic variants only
        temp.df <- temp.df[
            !is.na(
                code2germhits[[
                    paste0("_",
                    temp.df[,VCF.CHROMOSOME.COL],
                     ":",
                    temp.df[,VCF.POS.COL] / consolidation.factor)
                ]]
            )
        ,]
        
        grouped.names <- paste0("_",
            temp.df[,VCF.CHROMOSOME.COL],
             ":",
            as.integer(temp.df[,VCF.POS.COL] / consolidation.factor))
        names.table <- table(grouped.names)

        temp.map <- hashmap(
            names(names.table),
            as.vector(names.table)
        )

        nonzero.cell.maps <- c(nonzero.cell.maps, temp.map)
        setTxtProgressBar(pb, i)
    }
    close(pb)

    nonzero.locs <- unique(
        c(unlist(sapply(nonzero.cell.maps, function(a) a$keys())))
    )

    M <- matrix(0, nrow = length(nonzero.locs), ncol = length(cells))
    rownames(M) <- nonzero.locs
    colnames(M) <- cells

    print("Collating and generating variation matrix")
    pb <- txtProgressBar(min = 0,
        label = "Collating and generating variation matrix",
        max = length(cells),
        style = 3)

    for (i in 1:length(cells)) {
        colvals <- nonzero.cell.maps[[i]][[nonzero.locs]]
        colvals[is.na(colvals)] <- 0
        M[,i] = colvals
        setTxtProgressBar(pb, i)
    }

    close(pb)

    vcf.melange <- Melange(M)

    return(vcf.melange)         
}