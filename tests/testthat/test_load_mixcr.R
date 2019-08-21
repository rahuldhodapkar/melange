# test_load_mixcr.R
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

context("LoadMiXCR")
library(melange)

dname = "../testdata/mixcr"

cells <- c("cell1", "cell2", "cell3", "cell4", "cell5")
mixcr.files <- paste(cells, "_clones_TRB.txt", sep="")
mixcr.files.fullpath <- paste(dname, mixcr.files, sep="/")

mixcr.melange <- LoadMiXCR(
    cells=cells,
    mixcr.filenames=mixcr.files.fullpath
)

test_that("MiXCR files to Melange", {
  expect_is(mixcr.melange, "Melange")
})

test_that("Appropriate keys used", {
    expect_equal(mixcr.melange@data["TGTGCCAGCAGTGGAGACAGGGGGCTTGGAAACACCATATATTTT", "cell5"], 0.0526, tolerance = 1e-2)
})