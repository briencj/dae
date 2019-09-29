.onAttach <- function(...)
{ 

  if (!interactive() || sample.int(2, 1) == 1) 
    return()
  tips <- c("Need help? Enter help(package = 'dae') and click on 'User guides, package vignettes and other docs'.", 
            "Find out what has changed in dae: enter news(package = 'dae').",
            "Need help to produce randomized designs? Enter vignette('DesignNotes', package = 'dae').", 
            "Need help to do the canonical analysis of a design? Enter vignette('DesignNotes', package = 'dae').", 
            "Use suppressPackageStartupMessages() to eliminate all package startup messages.", 
            "To see all the intermittent, randomly-presented, startup tips enter ?daeTips.",
            "For versions between CRAN releases (and more) go to http://chris.brien.name/rpackages.")
  tip <- sample(tips, 1)
  packageStartupMessage(paste(strwrap(tip), collapse = "\n"))
}

