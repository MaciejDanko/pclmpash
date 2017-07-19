# Correct the problem when the both vectors are indentical, but has only numerical fluctuactions coming from the double type coding

#' Rounding a number to the machine precission to avoid numerical mistakes during comparison
#'
#' @keywords internal
machine.round <- function(x) round(x,floor(-log10(.Machine$double.eps^0.8)))

#' Test if Age Vector Matches nx Vector
#'
#' @keywords internal
TestnxMatchx <- function (x, nx) {
  if(any(machine.round(diff(x)) != machine.round(nx[-length(nx)]))) { stop("Provided Age and nx vectors don't match.", call. = FALSE) }
}
