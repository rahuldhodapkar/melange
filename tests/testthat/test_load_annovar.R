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

cells <- c("cell1", "cell2", "cell3", "cell4", "cell_no_exons")
annovar.files <- paste(cells, "anno.variant_function", sep=".")
annovar.files.fullpath <- paste(dname, annovar.files, sep="/")

germline.filename <- paste(dname, 'germline.anno.variant_function', sep='/')

annovar.assay <- LoadAnnovar(
    cells=cells,
    annovar.filenames=annovar.files.fullpath,
    germline.filename=germline.filename
)

test_that("Annovar files to Melange", {
  expect_is(annovar.assay, "Melange")
})