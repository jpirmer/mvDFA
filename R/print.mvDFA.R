#' print object of class DFA
#' @param x object of class DFA to print.
#' @param ... further parameters.
#' @export

print.mvDFA <- function(x, ...)
{
     obj <- unclass(x)
     print(obj[1:10])
     invisible(x)
}
