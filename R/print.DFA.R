#' print object of class DFA
#' @param x object of class DFA to print.
#' @param ... further parameters.
#' @exportS3Method print DFA
#' @export

print.DFA <- function(x, ...)
{
     obj <- unclass(x)
     print(obj[1:2])
     invisible(x)
}
