## landscape-sim

This repository contains files to run movement simulations and analyze resulting data for the PNAS 2018 manuscript: L.A. White, J.D. Forester & M.E. Craft. "Disease outbreak thresholds emerge from interactions between movement behavior, landscape structure, and epidemiology"

www.pnas.org/cgi/doi/10.1073/pnas.1801383115

### Scripts:
1) *sim_parallel_functions.R*- contains IBM_loop_function and sub functions needed to run simulations
2) *move_sim1.R*: Sample script to run IBM_loop_function in parallel via clusters
3) *Manuscript_figures_2017_12_12.Rmd*- Analyzes IBMsummary.csv simulation results to produce figures for manuscript
4) *RandomForestParallelLogit.R*- random forest analysis to test- is an outbreak successful?
5) *RandomForestParallelLogitDur.R*- given a successful outbreak- what determines outbreak duration?
*RandomForestParallelLogitPrev.R*- given a successful outbreak- what determines maximum prevalence?
6) *SimpleSIR.Rmd*- compare spatially-explicit simulations to comparable stochastic model that assumes homogeneous mixing (produces Figure S7)
7) *PartyRFLogit.R*- random forest analysis (using `cforest` function from `party` package) to test- is an outbreak successful?
8) *PartyRFLogitPrev.R*-random forest analysis (using `cforest` function from `party` package) to test- is an outbreak successful?
9) *PartyRFLogitDur.R*- random forest analysis (using `cforest` function from `party` package) to test-  what determines outbreak duration?

### Data/results:
1) *IBMsummary.zip*- contains *IBMsummary.csv* with all simulation results; headers: X= parameter set/combination; k= size of landscape (2^k+1); density= conspecific density (0.25 or 0.5); rec_rate=recovery rate; p=proportion available habitat; H=Hurst exponent; beta1 =strength of selection for resources; beta2= strength of selection for conspecifics; beta3= squared term in RSF for strength of selection of conspecifics; percep= perceptual range (1,2, or 3 cells out from current cell); duration= duration of outbreak; max_I= maximum number of individuals infected at any given time during outbreak; max_prev= maximum prevalence; betas= character string of betas used for parameter set; pbyH= string of values used for p and H to create landscape structure
2) *imp_logit.csv* and *imp_logitSD.csv*- variable importance and standard deviation (SD) results from *RandomForestParallelLogit.R*
3) *logit_imp_dur.csv* and *imp_logitdurSD.csv*- variable importance and standard deviation (SD) results  from *RandomForestParallelLogitDur.R*
4) *logit_imp_maxprev.csv* and *imp_logitprevSD.csv*- variable importance and standard deviation (SD) results from *RandomForestParallelLogitPrev.R*
5) *SIRdata_merged.csv*- simplistic, stochastic, homogeneous SIR model results produced from *SimpleSIR.Rmd*
6) *partyRF_logit1000.csv*- variable importance reuslts from *PartyRFLogit.R* using 1,000 trees
7) *partyRF_logitdur1000.csv* - variable importance reuslts from *PartyRFLogitDur.R* using 1,000 trees
8) *partyRF_logitprev1000.csv* - variable importance reuslts from *PartyRFLogitPrev.R* using 1,000 trees
