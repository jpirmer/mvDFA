#' print object of class DFA
#' @param x object of class DFA to print.
#' @param ... further parameters.
#' @exportS3Method print DFA
#' @return Truncates the output printed into the console of objects of class \code{DFA}, but does not change object itself.
#' @export

print.DFA <- function(x, ...)
{
     obj <- unclass(x)
     obj <- lapply(obj[1:2], function(o) round(o, 3))
     print(obj[1:2])
     invisible(x)
}
