#Startup code
#register emmeans methods dynamically
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("emmeans", quietly = TRUE)) {
    if (utils::packageVersion("emmeans") < "1.4") {
      warning("please install a newer version of emmeans (>=1.4)")
      return(NULL)
    }
    emmeans::.emm_register(c("combTMB"), pkgname)
  }
}


.onUnload <- function (libpath) {
   library.dynam.unload("combTMB", libpath)
  library.dynam.unload("combTMB_TMBExports", libpath)
 }


#http://patorjk.com/software/taag/#p=testall&h=0&v=1&f=Blocks&t=combTMB
combTMBStartupMessage <- function()
{
  msg <- c(paste0("
                  _   _____  _____  _____
 ___  ___  _____ | |_|_   _||     || __  |
|  _|| . ||     || . | | |  | | | || __ -|
|___||___||_|_|_||___| |_|  |_|_|_||_____| version",
    utils::packageVersion("combTMB")),
    "\nType 'citation(\"combTMB\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- combTMBStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'combTMB' version",  utils::packageVersion("combTMB"))
  packageStartupMessage(msg)
  invisible()
}


