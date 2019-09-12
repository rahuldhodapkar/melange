# test_load_rsem.R
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

context("LoadRSEM")
library(melange)

dname = "../testdata/rsem"

cells <- c("cell1", "cell2", "cell3")
rsem.files <- paste(cells, "rsem_quant.genes.results", sep=".")
rsem.files.fullpath <- paste(dname, rsem.files, sep="/")

rsem.assay.TPM <- LoadRSEM(
    cells = cells,
    rsem.filenames = rsem.files.fullpath,
    quantitation.method = 'TPM'
)

rsem.assay.count <- LoadRSEM(
    cells = cells,
    rsem.filenames = rsem.files.fullpath,
    quantitation.method = 'count'
)

test_that("RSEM files to Melange", {
  expect_is(rsem.assay.TPM, "Melange")
})
test_that("TPM generates correct counts", {
  expect_equal(rsem.assay.TPM@data["ERCC-00136","cell1"], 663.68)
  expect_equal(rsem.assay.TPM@data["ERCC-00136","cell2"], 0)
})
test_that("count generates correct counts", {
  expect_equal(rsem.assay.count@data["ERCC-00136","cell1"], 83.0)
  expect_equal(rsem.assay.count@data["ERCC-00136","cell2"], 0)
})


rsem.assay.min.quant.value <- LoadRSEM(
    cells = cells,
    rsem.filenames = rsem.files.fullpath,
    quantitation.method = 'TPM',
    min.quant.value=-1
)

test_that("include all genes", {
  expect_true(! "ERCC-00164" %in% rownames(rsem.assay.TPM))
  expect_true("ERCC-00164" %in% rownames(rsem.assay.min.quant.value@data))
})
