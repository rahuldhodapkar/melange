# seurat_wrapper.R
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


#' Convert data from Melange to Seurat Assay
#'
#' Melange operations return a melange data object for interoperability
#' with different analysis pipelines, but many users may prefer to do
#' single-cell analysis within Seurat. As such, the \code{SeuratAssayFromMelange}
#' function provides an easy way to convert between the two formats.
#' Mostly, this is done by re-naming and re-organizing some data elements.
#'
#' @param object melange object to be converted into Seurat Assay form for
#'     easier downstream processing and data integration.
#'
#' @return A Seurat Assay object containing the data collected in the Melange
#'     object.
#'
#' @examples
#' \dontrun{
#' 
#' SeuratAssayFromMelange(melange.object)
#' 
#' }
#'
#' @importFrom Seurat CreateAssayObject
#'
#' @rdname SeuratAssayFromMelange
#' @export SeuratAssayFromMelange
#'
SeuratAssayFromMelange <- function(object) {
    seurat.assay.object <- CreateAssayObject(data=object@data)

    seurat.assay.object@misc$meta.data <- object@meta.data
    seurat.assay.object@meta.features <- object@meta.features

    return(seurat.assay.object)
}
