mvDFAStartupMessage <- function()
{
     msg <- c(paste0("This is a beta version. Please report any bugs! mvDFA version ",
          packageVersion("mvDFA")),
          "\nType 'citation(\"mvDFA\")' for citing this R package in publications.")
     return(msg)
}

.onAttach <- function(lib, pkg)
{
     # unlock .mclust variable allowing its modification
     #unlockBinding(".mvDFA", asNamespace("mvDFA"))
     # startup message
     msg <- mvDFAStartupMessage()
     if(!interactive())
          msg[1] <- paste("Package 'mvDFA' version", packageVersion("mvDFA"))
     packageStartupMessage(msg)
     invisible()
}
