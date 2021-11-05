# Covariate Distribution Balance via Propensity Scores: <br> Replication files for Sant'Anna, Song and Xu (2021)

Contain codes and data for application and Monte Carlo simulations in  Sant'Anna, Song and Xu (2021), "Covariate Distribution Balance via Propensity Scores"

## Empirical Aplication
Our empirical application builds on [Benjamin (2003)](https://www.sciencedirect.com/science/article/abs/pii/S0047272701001670), [Abadie (2003)](https://www.sciencedirect.com/science/article/abs/pii/S0304407602002014), [Chernozhukov and Hansen (2004)](https://www.jstor.org/stable/3211794). We assess the causal effect of 401(k) eligibility and participation on different measures of wealth. We use 1991  Survey of Income and Program Participation (SIPP) data, exactly as in [Benjamin (2003)](https://www.sciencedirect.com/science/article/abs/pii/S0047272701001670), and [Chernozhukov and Hansen (2004)](https://www.jstor.org/stable/3211794). The folder `Application/401/data` contains the data in `Stata` and `R` format.


All codes to replicate the results of our empirical application are available in the `Applications/401k` folder. The `run-all.R` file will run all codes and produce all the results. The generated plots are saved in the `Applications/401k/plots` folder.

## Monte Carlo Simulations
All files to replicate our Monte Carlo results are available in the `Simulations` folder. The R file `run-all-simulations.R` source all scripts.

## Problems
If you have any trouble running the code, or find any errors, please file an issue on this repo and we will look into it.
