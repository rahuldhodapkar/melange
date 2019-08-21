# classes.R
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Melange Class
#'
#' The Melange object stores the cell-specific information aggregated from whichever
#' tool the user has chosen to import from. The number of rows (features) as well
#' as the datatype of the core data matrix depends on the particular input file
#' type, 
#'
#' @slot data           Raw counts data for the loaded data source. Can be any 
#'                      matrix.
#' @slot source         Unambiguous name for the source data in the Melange object. 
#' @slot meta.features  Feature-level metadata
#' @slot meta.data      Cell-level metadata
#' @slot misc           Utility slot for storing additional data
#'
Melange <- setClass(
  Class = 'Melange',
  slots = c(
    data = 'AnyMatrix',
    source = 'character',
    meta.features = 'data.frame',
    meta.data = 'data.frame',
    misc = 'ANY'
  )
)

#' Constructor function for the melange class 
#'
#' @param data.matrix  Matrix with data for the specific data type interrogated
#'
#' @name Melange
#' @rdname Melange
#' @export Melange
#'
#' @importFrom methods new
#'
Melange <- function(data.matrix) {

    init.meta.features <- data.frame(row.names = rownames(x = data.matrix))
    init.meta.data <- data.frame(row.names = colnames(x = data.matrix))

    melange <- new(
        Class = 'Melange',
        data = data.matrix,
        meta.features = init.meta.features,
        meta.data = init.meta.data
    )

    return(melange)
}

