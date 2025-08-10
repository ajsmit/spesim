# Hard refresh package build

.rs.restartR() # if using RStudio; otherwise quit R and reopen
unlink("NAMESPACE") # let roxygen rebuild it
devtools::clean_dll() # clear compiled artefacts if any
devtools::document() # rebuild NAMESPACE + man
devtools::install(quick = FALSE) # force full reinstall
