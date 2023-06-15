#' print object of class mvDFA
#' @param x object of class DFA to print.
#' @param ... further parameters.
#' @exportS3Method print mvDFA
#' @return Truncates the output printed into the console of objects of class \code{mvDFA}, but does not change object itself.
#' @export

print.mvDFA <- function(x, ...)
{
     obj <- unclass(x)
     obj <- lapply(obj[1:5], function(o) round(o, 3))
     print(obj)
     invisible(x)
}
