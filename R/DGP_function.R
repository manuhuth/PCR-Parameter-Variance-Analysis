files <- c('PCA_PropVar', 'PCA', 'PCR_cv', 'PCR_predict', 'PCR', 'random_discreteVariables', 'random_VCV') #names of files to read
for (i in 1:length(files)) { #loop to read all files
  source(paste('R/', files[i], '.R', sep = ''))
}


dgp_model <- function(n = 1000, var_err = 1,
                      mean_parent_educ = 13.342, variance_parent_educ = 21.215,
                      var_ability = 1, gamma_ability = 2,  gamma_parent_educ = 2,
                      breaks_test7_m = c(0, 0.141, 0.158, 0.185, 0.190, 0.212),  breaks_test11_m = c(0, 0.122,0.152, 0.157, 0.179, 0.199),
                      breaks_test7_r =  c(0, 0.166, 0.179, 0.188, 0.187, 0.165),  breaks_test11_r = c(0, 0.132, 0.163, 0.163, 0.176, 0.176), test_cat = TRUE,
                      mean_numberSiblings = 1.692, var_numberSiblings = 1.7^2,  max_yearsSchooling = 29, min_yearsSchooling = 0, prob_numbSiblings = 0.1,
                      mean_schooling = 13.342, variance_schooling = 21.215, q = 0.85, age_min = 33, age_max = 68,
                      probs_gap = c(0.59, 0.11,0.7,0.04, 0.03), gap_years = c(0,1,2,3,4), age_school_count = 4,
                      bound_constant_low = -4, bound_constant_up = 4, beta_min = c(0.03, 0.01, -0.06, -10, 0.01), beta_max = c(0.06, 0.06, -0.03, 10, 10),
                      variance_error_wage = 0.1,  mean_wage = 2.040, variance_wage = 1.5,  tau = 0.5,
                      starting_values = c(3.5, 0.046, 0.026, -0.037, 0, 0.05), mu = 1e-4, outer.it = 500, outer.eps = 1e-10

) {

  #parents' education
  mean_parent_educ <- mean_parent_educ
  eta_parent_educ <- 1 - sqrt(mean_parent_educ)/sqrt(variance_parent_educ)
  lambda_parent_educ <- mean_parent_educ*(1-eta_parent_educ)
  parent_educ_randomDraw <- rgenPois(n = n, lambda = lambda_parent_educ, eta = eta_parent_educ )
  parent_educ <- parent_educ_randomDraw$random_value
  mean_parent_educ <- mean_parent_educ

  parent_educ_unif <- parent_educ_randomDraw$distribution_value

  #Ability
  ability <- rnorm(n = n, mean = 0, sd = var_ability^0.5)

  #test scores -> modified version of Heckman, Muller (effect of schooling and ability on average test scores)
  err_test7_m <- rnorm(n = n, mean = 0, sd = var_err^0.5)
  test7_m_latent <- gamma_ability*ability + gamma_parent_educ*qnorm(parent_educ_unif, mean = 0, sd = var_ability^0.5) + err_test7_m
  test7_m <- pnorm(test7_m_latent, mean = 0, sd = ((gamma_ability^2 + gamma_parent_educ^2)*var_ability+var_err)^0.5) # uniform distributed

  err_test11_m <- rnorm(n = n, mean = 0, sd =var_err^0.5)
  test11_m_latent <- test7_m_latent + err_test11_m
  test11_m <- pnorm(test7_m_latent, mean = 0, sd = ((gamma_ability^2 + gamma_parent_educ^2)*var_ability + 2*var_err)^0.5) # uniform distributed

  err_test7_r <- rnorm(n = n, mean = 0, sd = var_err^0.5)
  test7_r_latent <- gamma_ability*ability +  gamma_parent_educ*qnorm(parent_educ_unif, mean = 0, sd = var_ability^0.5) + err_test7_r
  test7_r <- pnorm(test7_r_latent, mean = 0, sd = ((gamma_ability^2 + gamma_parent_educ^2)*var_ability+var_err)^0.5)# uniform distributed
  err_test11_r <- rnorm(n = n, mean = 0, sd = var_err^0.5)
  test11_r_latent <- test7_r_latent + err_test11_r
  test11_r <-pnorm(test7_r_latent, mean = 0, sd = ((gamma_ability^2 + gamma_parent_educ^2)*var_ability+ 2*var_err)^0.5)# uniform distributed


  breaks_test7_m = breaks_test7_m / sum(breaks_test7_m)
  breaks_test11_m = breaks_test11_m / sum(breaks_test11_m)
  breaks_test7_r = breaks_test7_r / sum(breaks_test7_r)
  breaks_test11_r = breaks_test11_r / sum(breaks_test11_r)
  if (isTRUE(test_cat)) {
    test7_m = as.numeric(cut(test7_m, breaks = cumsum(breaks_test7_m))) #categorize the test results and make them numeric (cut generates factors)
    test11_m = as.numeric(cut(test11_m, breaks = cumsum(breaks_test11_m)))
    test7_r = as.numeric(cut(test7_r, breaks = cumsum(breaks_test7_r)))
    test11_r = as.numeric(cut(test11_r, breaks = cumsum(breaks_test11_r)))
  }

  #number of siblings
  transformed_parentEduc <- max_yearsSchooling - parent_educ #to make it feasible with theory
  transformed_parentEduc <- replace(transformed_parentEduc, transformed_parentEduc < 0, 0)
  expec_var <- (max_yearsSchooling - mean_parent_educ)*prob_numbSiblings*(1-prob_numbSiblings)
  var_expec <- prob_numbSiblings^2 * variance_parent_educ
  eta_err_numbSiblings <- 1 - sqrt(mean_numberSiblings - (max_yearsSchooling - mean_parent_educ)*prob_numbSiblings) / sqrt(var_numberSiblings - expec_var - var_expec)
  lambda_err_numbSiblings <- (mean_numberSiblings - (max_yearsSchooling - mean_parent_educ)*prob_numbSiblings)*(1-eta_err_numbSiblings)
  err_numbSiblings <- rgenPois(n, lambda = lambda_err_numbSiblings, eta = eta_err_numbSiblings)$random_value
  numb_Siblings_binomSum <- sapply(transformed_parentEduc, rbinom_wrapper, p = prob_numbSiblings)
  numb_Siblings <- numb_Siblings_binomSum + err_numbSiblings


  #number of years of schooling -> generalized Poisson process to be in line with parents, last hard part, maybe incorporate channel from parents educ to test scores (through ability)
  k <- 2*(mean_schooling- mean_parent_educ *q)^(3/2) / ( (variance_schooling - mean_parent_educ * q* (1-q) - variance_parent_educ * q - 1/3 * (mean_schooling - mean_parent_educ*q)^2)^0.5)

  if (k < 0) {
    k = - k #adjust for two solution case
  }

  eta_err_schooling <- 1 - k/(2*(mean_schooling - mean_parent_educ*q))
  lambda_sch <- k*pnorm(ability, 0, var_ability^0.5)

  err_schooling <- sapply(lambda_sch, rgenPois, n = 1, eta = eta_err_schooling, length_cumDensity = 1000)
  err_schooling <- unlist(err_schooling)
  err_schooling <- err_schooling[c(TRUE, FALSE)]

  #err_schooling <- as.numeric(err_schooling[err_schooling != 'cumDensity'])
  schooling_binomSum <- sapply(parent_educ, rbinom_wrapper, p = q)
  schooling <- schooling_binomSum + err_schooling
  schooling <- replace(schooling, schooling > max_yearsSchooling, max_yearsSchooling)
  schooling <- replace(schooling, schooling < min_yearsSchooling, min_yearsSchooling)

  #age and working
  UK_pop <- as.matrix(read_excel("UK Population Estimates 1838-2015.xls",  sheet = "GB SYOA 1961-2015", range = "X2:X94"))
  UK_pop <- UK_pop[-c(1:(age_min+1), (age_max+1):90)]
  freq <- UK_pop/sum(UK_pop)
  age <- rdiscrete(n, freq, 33:70)
  gap <- rdiscrete(n, probs_gap, gap_years)
  working <- age - schooling - age_school_count - gap
  working <- replace(working, working < 0, 0)

  #logarithmic wage

  working_sqr <- working^2/100
  covariates <- data.frame(ability, schooling, working, working_sqr, numb_Siblings, parent_educ)
  VCV_cov <- var(covariates)
  error_wage <- rnorm(n, 0, variance_error_wage^0.5)

  #constraints
  ui <- rbind(diag( rep(1, (ncol(covariates))) ), diag(rep(-1, (ncol(covariates)))))
  ci <- c(bound_constant_low, beta_min, -bound_constant_up, -beta_max)

  wage_optim <- function(par) {
    beta <- par[-1]
    a <- par[1]
    beta_calc <- c(1, beta) #add ability for covariance

    objective <- tau*(mean_wage - a - colMeans(covariates) %*% beta_calc)^2 + (1-tau)*(variance_wage^0.5 - (t(beta_calc) %*% VCV_cov %*% beta_calc - variance_error_wage)^0.5 )^2
    return(objective)
  }

  #constrained optimization process
  fit <- constrOptim(theta = starting_values , f = wage_optim, ui = ui, ci = ci, mu = mu,
                     method = "Nelder-Mead",
                     outer.iterations = outer.it, outer.eps = outer.eps,
                     hessian = FALSE)


  #append to data frame
  data <- data.frame(ability, test7_m, test11_m, test7_r, test11_r, parent_educ, schooling, numb_Siblings, working)
  matrix_to_create_data <- data.frame(rep(1,n), schooling, working, working_sqr, numb_Siblings, parent_educ)

  logwage <- ability + as.matrix(matrix_to_create_data) %*% fit$par + error_wage

  data <- cbind(logwage, data, age)
return(data)
}

