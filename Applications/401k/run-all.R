#-----------------------------------------------------------------------------
# This script run all scripts to replicate the empirical application  
#-----------------------------------------------------------------------------
# Library
library(here)
#-----------------------------------------------------------------------------
# Source files
source(here("Applications/401k/01-401k_itt_application.R"))
source(here("Applications/401k/02-401K_itt_plots.R"))
source(here("Applications/401k/03-401k_lte_application.R"))
source(here("Applications/401k/04-401K_lte_plots.R"))
#-----------------------------------------------------------------------------