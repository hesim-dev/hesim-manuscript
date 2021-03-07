#' ---
#' title: "hesim: Health Economic Simulation Modeling and Decision Analysis"
#' author: "Devin Incerti, Jeroen P Jansen"
#' date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#' ---

#' This is a replication file for the manuscript. It can be used to produce
#' an HTML file that displays the results without the text with:
# rmarkdown::render("replication.R", output_format = "html_document")
#' Replication was performed by the authors on macOS Mojave 10.14.6. See the
#' [Session Information](#sessionInfo) below for further details.

#+ message = FALSE, warning = FALSE
#' # Settings and dependencies
#' Package versioning is managed using the `renv` package. Packages can be 
#' installed by cloning the GitHub repository available at 
#' https://github.com/hesim-dev/hesim-manuscript and running `renv::restore()`. 
library("data.table")
library("flexsurv")
library("ggplot2")
library("hesim")
library("knitr")
library("magrittr")
library("msm")
opts_chunk$set(
  tidy = FALSE, eval = TRUE, cache = TRUE,
  fig.height = 4, fig.width = 6
)
options(prompt = "> ", continue = "+  ", width = 120, useFancyQuotes = TRUE)
ggplot2::theme_set(ggplot2::theme_bw())
set.seed(102)

#' # Illustrative example
#' ## Individual continuous time state transition model
#' ### Setup model
#' #### Target population, treatment strategies, and model structure
# Define transitions
tmat <- rbind(
  c(NA, 1, 2),
  c(NA, NA, 3),
  c(NA, NA, NA)
)
colnames(tmat) <- rownames(tmat) <- c("Stable", "Progression", "Dead")
print(tmat)

# hesim data
n_patients <- 1000
patients <- data.table(
  patient_id = 1:n_patients,
  age = rnorm(n_patients, mean = 45, sd = 7),
  female = rbinom(n_patients, size = 1, prob = .51)
)

states <- data.table(
  state_id = c(1, 2),
  state_name = c("Stable", "Progression") # Non-death health states
)

strategies <- data.table(
  strategy_id = 1:3,
  strategy_name = c("SOC", "New 1", "New 2"),
  soc = c(1, 0, 0),
  new1 = c(0, 1, 0),
  new2 = c(0, 0, 1)
)

hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients,
  states = states
)
print(hesim_dat)

#' #### Labels for ID variables
labs_indiv <- get_labels(hesim_dat)
print(labs_indiv)

#' ### Parameterization
#' #### Transition model
# Estimation data
onc3[patient_id %in% c(1, 2)] %>%
  kable()

# Fit multi-state model
n_trans <- max(tmat, na.rm = TRUE) 
wei_fits <- vector(length = n_trans, mode = "list")
f <- as.formula(Surv(time, status) ~ factor(strategy_name) + female + age)
for (i in 1:length(wei_fits)){
  if (i == 3) f <- update(f, .~.-factor(strategy_name)) 
  wei_fits[[i]] <- flexsurvreg(f, data = onc3, 
                               subset = (transition_id == i),
                               dist = "weibull")
}
wei_fits <- flexsurvreg_list(wei_fits)

#' #### Utility model
utility_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(.8, .6),
             se = c(0.02, .05)
  ),
  dist = "beta")
print(utility_tbl)

#' #### Cost models
# Medical costs
medcost_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(2000, 9500),
             se = c(2000, 9500)
  ),
  dist = "gamma")
print(medcost_tbl)

# Drug costs
n_times <- 2
n_states <- nrow(states)
n_strategies <- nrow(strategies)
n_i <- n_states * n_times
drugcost_dt <-  data.table(
    strategy_id = rep(strategies$strategy_id, each = n_i),
    state_id = rep(rep(states$state_id, each = n_times), n_strategies),
    time_start = rep(c(0, 3/12), n_states * n_strategies),
    est = c(2000, 2000, 1500, 1200,
           12000, 12000, 1500, 1200,
           15000, 15000, 1500, 1200)
) 

drugcost_tbl <- stateval_tbl(
  drugcost_dt,
  dist = "fixed")
print(drugcost_tbl)


#' ### Simulation
#' #### Constructing the economic model
n_samples <- 1000 # Number of samples for PSA

#' ##### Transition model
transmod_data <- expand(hesim_dat,
                        by = c("strategies", "patients"))
head(transmod_data)
transmod <- create_IndivCtstmTrans(wei_fits, transmod_data,
                                   trans_mat = tmat, n = n_samples,
                                   uncertainty = "normal",
                                   clock = "reset",
                                   start_age = patients$age)

#' ##### Utility and cost models
# Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples, 
                               hesim_data = hesim_dat)

# Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples,
                                time_reset = TRUE, hesim_data = hesim_dat)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples,
                               hesim_data = hesim_dat)
costmods <- list(Drug = drugcostmod,
                 Medical = medcostmod)

#' ##### Economic model
ictstm <- IndivCtstm$new(trans_model = transmod,
                         utility_model = utilitymod,
                         cost_models = costmods)

#' #### Simulating outcomes
#' ##### Individual trajectories through model
ictstm$sim_disease(max_age = 100)
head(ictstm$disprog_)

#' ##### State probabilities
ictstm$sim_stateprobs(t = seq(0, 30 , 1/12))
head(ictstm$stateprobs_)
autoplot(ictstm$stateprobs_, labels = labs_indiv,
         ci = FALSE) 

#' ##### QALYs and costs
# QALYs
ictstm$sim_qalys(dr = c(0,.03))
head(ictstm$qalys_)

# Costs
ictstm$sim_costs(dr = .03)
head(ictstm$costs_)

#' ##### Summarize QALYs and costs
ce_sim_ictstm <- ictstm$summarize()
summary(ce_sim_ictstm, labels = labs_indiv) %>%
  format() %>%
  kable()

#' ## Partitioned survival model
#' ### Setup model
#' #### Target population for cohort model
hesim_dat$patients <- data.table(
  patient_id = 1:4,
  grp_id = c(1, 1, 2, 2),
  patient_wt = rep(1/4, 4),
  age = c(55, 65, 55, 65),
  female = c(1, 1, 0, 0),
  grp_name = rep(c("Female", "Male"), 
                 each = 2)
)

#' #### Labels for cohort model
labs_cohort <- get_labels(hesim_dat)

#' ### Parameterization
#' #### Survival models
#' Parameters are estimated using both individual patient data and summary data:
#' a model is fit for SOC and relative treatment effects are obtained from the 
#' literature.
# Estimation data
onc_pfs_os <- as_pfs_os(onc3, patient_vars = c("patient_id", "female", "age", 
                                               "strategy_name"))
onc_pfs_os[patient_id %in% c(1, 2)]

# Parameters
## Fit models for SOC
onc_pfs_os_soc <- onc_pfs_os[strategy_name == "SOC"]
fit_pfs_soc_wei <- flexsurv::flexsurvreg(
  Surv(pfs_time, pfs_status) ~ female,
  data = onc_pfs_os_soc,
  dist = "weibull")

fit_os_soc_wei <- flexsurvreg(
  Surv(os_time, os_status) ~  female,
  data = onc_pfs_os_soc,
  dist = "weibull")

psmfit_soc_wei <- partsurvfit(
  flexsurvreg_list(pfs = fit_pfs_soc_wei, os = fit_os_soc_wei),
  data = onc_pfs_os_soc
)

## Obtain relative treatment effects
## Fit using the "onc_pfs_os" dataset but assume the estimates are
## obtained from the literature for illustrative purposes
### Create dummies
onc_pfs_os[, new1 := ifelse(strategy_name == "New 1", 1, 0)]
onc_pfs_os[, new2 := ifelse(strategy_name == "New 2", 1, 0)]

### Fit models
fit_pfs_wei <- flexsurv::flexsurvreg(
  Surv(pfs_time, pfs_status) ~ female + new1 + new2,
  data = onc_pfs_os,
  dist = "weibull")

fit_os_wei <- flexsurvreg(
  Surv(os_time, os_status) ~ female + new1 + new2,
  data = onc_pfs_os,
  dist = "weibull")

### Extract coefficients (assumed taken from literature)
params_pfs_wei <- create_params(fit_pfs_wei, n = n_samples)
params_os_wei <- create_params(fit_os_wei, n = n_samples)
coef_pfs_wei <- params_pfs_wei$coefs$scale[, c("new1", "new2")]
coef_os_wei <- params_os_wei$coefs$scale[, c("new1", "new2")]
head(coef_pfs_wei)

#' #### Costs
# Update drug costs to work with cohort model
drugcost_tbl <- stateval_tbl(
  drugcost_dt[time_start ==  0][, time_start := NULL],
  dist = "fixed")
print(drugcost_tbl)

#' ### Simulation
#' #### Constructing the economic model
#' ##### Survival models
# Input data
survmods_data <- expand(hesim_dat, by = c("strategies", "patients"))
head(survmods_data)

# Parameters
## Parameters for SOC. We obtain a warning because the model did not
## converge during one of the bootstrap samples; this is allowed because
## max_errors = 5 >= 1.
survmods_params <- create_params(psmfit_soc_wei, n = n_samples, 
                                 uncertainty = "bootstrap", 
                                 max_errors = 5, silent = TRUE)

## "Manually" add relative treatment effects
survmods_params$pfs$coefs$scale <- cbind(survmods_params$pfs$coefs$scale,
                                         coef_pfs_wei)
survmods_params$os$coefs$scale <- cbind(survmods_params$os$coefs$scale,
                                        coef_os_wei)
survmods_params$pfs$coefs$scale[1:2, ]

# Combine input data and parameters to create the survival models
survmods_data[, ("(Intercept)") := 1]
survmods <- create_PsmCurves(survmods_params, 
                             input_data = survmods_data)

#' ##### Utility and cost models
utilitymod <- create_StateVals(utility_tbl, n = n_samples, 
                               hesim_data = hesim_dat)
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples, 
                                hesim_data = hesim_dat)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples, 
                               hesim_data = hesim_dat)
costmods <- list(Drug = drugcostmod, Medical = medcostmod)

#' ##### Economic model
psm <- Psm$new(survival_models = survmods,
               utility_model = utilitymod,
               cost_models = costmods)

#' #### Simulating outcomes
#' ##### Survival
times <- seq(0, 30, by = .1)
psm$sim_survival(t = times)
head(psm$survival_)

autoplot(psm$survival_, 
         labels = c(labs_cohort,
                    list(curve = c("PFS" = 1, "OS" = 2))
                    ),
         ci = TRUE, ci_style = "ribbon")


#' ##### State probabilities
psm$sim_stateprobs()
psm$stateprobs_[sample == 1 & patient_id == 1 & state_id == 2 & t == 12]


#' ##### QALYs and costs
psm$sim_qalys(dr = .03) 
psm$sim_costs(dr = .03)

#' ##### Summarize QALY and costs
ce_sim_psm <- psm$summarize(by_grp = TRUE)
summary(ce_sim_psm, labels = labs_cohort) %>%
  format() %>%
  kable()

#' ## Cohort discrete time state transition model
#' ### Parameterization
#' #### Transition model

#' We will first obtain the underlying parameter estimates. A model is fit
#' for SOC and relative treatment effects are obtained from the literature.

# Estimation data
head(onc3p)

# Fit multi-state model for SOC with msm
qinit <- matrix(0, nrow = 3, ncol = 3)
qinit[1, 2] <-  0.28; qinit[1, 3] <-  0.013; qinit[2, 3] <-  0.10
msm_fit <- msm(state_id ~ time, subject = patient_id, 
               data = onc3p[strategy_name == "SOC"],
               exacttimes = FALSE,
               covariates = ~ age + female,
               qmatrix = qinit, gen.inits = FALSE)

# Compare model fit with msm using exact times to an exponential model
## msm with exact transition times
tmat2 <- tmat; tmat2[is.na(tmat2)] <- 0
msm_efit <- msm(state_id ~ time, subject = patient_id, 
                data = onc3p[strategy_name == "SOC"],
                qmatrix = tmat2, exacttimes = TRUE)

## flexsurv with exponential model
exp_fits <- vector(length = n_trans, mode = "list")
for (i in 1:length(wei_fits)){
  exp_fits[[i]] <- flexsurvreg(Surv(time, status) ~ 1, data = onc3, 
                               subset = (transition_id == i & 
                                         strategy_name == "SOC"),
                               dist = "exponential")
}

## Compare the coefficients
msm_efit
exp_fits[[1]]$res
exp_fits[[2]]$res
exp_fits[[3]]$res

# Relative treatment effects (i.e., relative risk)
# Assume estimated from external source and simulate distribution
# for PSA
params_rr <- list(
  lrr_12_est = c(soc = log(1), new1 = log(.80), new2 = log(.71)),
  lrr_12_se = c(soc = 0, new1 = .03, new2 = .04),
  lrr_13_est = c(soc = log(1), new1 = log(.90), new2 = log(.85)),
  lrr_13_se = c(soc = 0, new1 = .02, new2 = .03)
)  

params_rr_def <- define_rng({
  list(
    rr_12 = lognormal_rng(lrr_12_est, lrr_12_se),
    rr_13 = lognormal_rng(lrr_13_est, lrr_13_se)
  )}, n = n_samples)
params_rr_rng <- eval_rng(params_rr_def, params_rr)

#' We now use the underlying parameter estimates and input data to generate
#' the "transformed" parameters.
# Start by predicting transition intensity matrices for each patient
# with SOC. Use a multivariate normal distribution to simulate samples
# for the PSA
transmod_data <- survmods_data
qmat_soc <- qmatrix(msm_fit, 
                    newdata = transmod_data[strategy_name == "SOC"],
                    uncertainty = "normal", n = n_samples)
dim(qmat_soc)
qmat_soc[,, 1]

# Take the matrix exponential to convert the transition intensity matrices to
# transition probability matrices
cycle_len <- 1/4
pmat_soc <- expmat(qmat_soc, t = cycle_len)
pmat_soc[,, 1]

# Align the transition matrices with the ID variables
tpmat_id <- tpmatrix_id(transmod_data, n_samples)
head(tpmat_id)

# Compute relative risks for each row in the input data and each parameter 
# sample for the PSA
xbeta <- function(x, beta) c(x %*% t(beta))
x_rr <- as.matrix(transmod_data[, .(soc, new1, new2)])
rr <- cbind(xbeta(x_rr, params_rr_rng$rr_12),
            xbeta(x_rr, params_rr_rng$rr_13))
head(rr) # Corresponds to 1st 6 rows in tpmat_id

# Apply the relative risks to SOC to get transition probabilities for each
# parameter sample, treatment strategy, and representative patient
pmat <- apply_rr(pmat_soc, rr = rr,
                 index = list(c(1, 2), c(1, 3)))
tprobs <- tparams_transprobs(pmat, tpmat_id)

#' ### Simulation
#' #### Constructing the economic model
#' ##### Transition model
transmod <- CohortDtstmTrans$new(params = tprobs, 
                                 cycle_length = cycle_len)

#' ##### Economic model: use the utility and cost models from the PSM
cdtstm <- CohortDtstm$new(trans_model = transmod,
                          utility_model = psm$utility_model,
                          cost_models = psm$cost_models)

#' #### Simulating outcomes
#' ##### State probabilities
cdtstm$sim_stateprobs(n_cycles = 30/cycle_len)

#' ##### QALYs and costs
cdtstm$sim_qalys(dr = .03)
cdtstm$sim_costs(dr = .03)

#' ##### Summarize QALYs and costs
ce_sim_cdtstm <- cdtstm$summarize()
summary(ce_sim_cdtstm, labels = labs_cohort) %>%
  format() %>%
  kable()

#' # Cost-effectiveness analysis
# "ce" object to use for analysis
ce_sim_ictstm

#' ## Perform cost-effectiveness analysis
wtp <- seq(0, 250000, 500) # Willingness to pay per QALY
cea_ictstm <- cea(ce_sim_ictstm, dr_costs = .03, dr_qalys = .03, k = wtp)
cea_pw_ictstm <- cea_pw(ce_sim_ictstm, comparator = 1,
                        dr_qalys = .03, dr_costs = .03,
                        k = wtp)

#' ## Incremental cost-effectiveness ratio
icer(cea_pw_ictstm, k = 100000, labels = labs_indiv) %>%
  format() %>%
  kable()

#' ## Representing decision uncertainty
#' ### Cost-effectiveness plane
plot_ceplane(cea_pw_ictstm, labels = labs_indiv)

#' ### Cost-effectiveness acceptability curve
#' #### Simultaneous comparison
plot_ceac(cea_ictstm, labels = labs_indiv)

#' #### Pairwise comparison
plot_ceac(cea_pw_ictstm, labels = labs_indiv)

#' ### Cost-effectiveness acceptability frontier
plot_ceaf(cea_ictstm, labels = labs_indiv)

#' ### Expected value of perfect information
plot_evpi(cea_ictstm, labels = labs_indiv)

#' # Session information {#sessionInfo}
sessionInfo()

