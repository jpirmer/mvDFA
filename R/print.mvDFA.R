#' print object of class mvDFA
#' @param x object of class DFA to print.
#' @param ... further parameters.
#' @exportS3Method print mvDFA
#' @export

print.mvDFA <- function(x, ...)
{
     obj <- unclass(x)
     obj <- lapply(obj[1:5], function(o) round(o, 3))
     print(obj)
     invisible(x)
}
