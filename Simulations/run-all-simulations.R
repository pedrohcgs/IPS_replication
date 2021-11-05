#-----------------------------------------------------------------------------
# This script run all scripts to replicate all Monte Carlo Simulations
#-----------------------------------------------------------------------------
# Library
library(here)
#-----------------------------------------------------------------------------
# Stylized simulations
#-----------------------------------------------------------------------------
# Unconfoundedness setup
source(here("Simulations/Stylized/main_exog.R"))
# LTE setup
source(here("Simulations/Stylized/main_endog.R"))
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Empiricaly motivated simulations
#-----------------------------------------------------------------------------
# Unconfoundedness setup
source(here("Simulations/401k/main_401k_itt.R"))
# LTE setup
source(here("Simulations/401k/main_401k_lte.R"))
#-----------------------------------------------------------------------------