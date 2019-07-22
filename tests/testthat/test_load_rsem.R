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
rsem_files <- paste(cells, "rsem_quant.genes.results", sep=".")
rsem_files_fullpath <- paste(dname, rsem_files, sep="/")

rsem_assay <- load.rsem(cells = cells, rsem.filenames = rsem_files_fullpath)

test_that("RSEM files to Seurat Assay", {
  expect_is(rsem_assay, "Assay")
  expect_equal(rsem_assay@counts["ERCC-00136","cell1"], 663.68 / 777.31)
  expect_equal(rsem_assay@counts["ERCC-00136","cell2"], 0)
})
