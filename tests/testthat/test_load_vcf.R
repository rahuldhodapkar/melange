# test_load_vcf.R
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

context("LoadVCF")
library(melange)

dname = "../testdata/vcf"

cells <- c("cell1", "cell2", "cell3")
vcf.files <- paste(cells, "vcf", sep=".")
vcf.files.fullpath <- paste(dname, vcf.files, sep="/")
germline.filename <- paste(dname, 'germline.vcf', sep='/')


vcf.melange <- LoadVCF(
    cells=cells,
    vcf.filenames=vcf.files.fullpath,
    germline.filename=germline.filename
)

test_that("Annovar files to Melange", {
  expect_is(vcf.melange, "Melange")
})