# test_seurat_conversion.R
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

context("ConvertToSeurat")
library(melange)

dname = "../testdata/rsem"

cells <- c("cell1", "cell2", "cell3")
rsem.files <- paste(cells, "rsem_quant.genes.results", sep=".")
rsem.files.fullpath <- paste(dname, rsem.files, sep="/")

rsem.assay.lengthScaledTPM <- SeuratAssayFromMelange(LoadRSEM(
    cells = cells,
    rsem.filenames = rsem.files.fullpath,
    quantitation.method = 'lengthScaledTPM'
))

rsem.assay.TPM <- SeuratAssayFromMelange(LoadRSEM(
    cells = cells,
    rsem.filenames = rsem.files.fullpath,
    quantitation.method = 'TPM'
))

rsem.assay.count <- SeuratAssayFromMelange(LoadRSEM(
    cells = cells,
    rsem.filenames = rsem.files.fullpath,
    quantitation.method = 'count'
))

test_that("RSEM files to Seurat Assay from Melange", {
  expect_is(rsem.assay.lengthScaledTPM, "Assay")
  expect_is(rsem.assay.TPM, "Assay")
})
