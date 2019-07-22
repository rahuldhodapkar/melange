# test_load_annovar.R
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

context("LoadAnnovar")
library(melange)

dname = "../testdata/annovar"

cells <- c("cell1", "cell2", "cell3")
annovar_files <- paste(cells, "anno.variant_function", sep=".")
annovar_files_fullpath <- paste(dname, annovar_files, sep="/")

germline.filename <- paste(dname, 'germline.anno.variant_function', sep='/')

annovar_matrix <- load.annovar(
    cells=cells,
    annovar.filenames=annovar_files_fullpath,
    germline.filename=germline.filename
)

test_that("RSEM files to Assay", {
  expect_is(annovar_matrix, "Assay")
})