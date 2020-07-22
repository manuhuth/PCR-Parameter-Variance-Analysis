files <- c('PCA_PropVar', 'PCA', 'PCR_cv', 'PCR_predict', 'PCR', 'random_discreteVariables', 'random_VCV') #names of files to read
for (i in 1:length(files)) { #loop to read all files
  source(paste('R/', files[i], '.R', sep = ''))
}


#DGP
set.seed(123)
n <- 2000
var_err <- 1

#parents' education
mean_parent_educ <- 7.2 #lambda/(1-eta) #function input
variance_parent_educ <- 4.827^2 # lambda/(1-eta)^3 #function input from A Generalization of the Poisson Distribution (Consul, Jain) (1973)
eta_parent_educ <- 1 - sqrt(mean_parent_educ)/variance_parent_educ^0.5
lambda_parent_educ <- mean_parent_educ*(1-eta_parent_educ)
parent_educ_randomDraw <- rgenPois(n = n, lambda = lambda_parent_educ, eta = eta_parent_educ )
parent_educ <- parent_educ_randomDraw$rnumb

#Ability
sd_ability <- 1 #function input
ability <- rnorm(n = n, mean = 0, sd = sd_ability)

#test scores -> modified version of Heckman, Muller (effect of schooling and ability on average test scores)
err_test7_m <- rnorm(n = n, mean = 0, sd = var_err^0.5) #if changed -> test7_m must be changed
gamma_ability <- 4 #function input
test7_m_latent <- gamma_ability*ability + err_test7_m
test7_m <- pnorm(test7_m_latent, mean = 0, sd = (gamma_ability^2*var_err+1)^0.5) # uniform distributed -> fits to data
err_test11_m <- rnorm(n = n, mean = 0, sd =var_err^0.5)
test11_m_latent <- test7_m_latent + err_test11_m
test11_m <- pnorm(test7_m_latent, mean = 0, sd = (gamma_ability^2*var_err + 1 + 1)^0.5) # uniform distributed -> fits to data

err_test7_r <- rnorm(n = n, mean = 0, sd = var_err^0.5)
test7_r_latent <- gamma_ability*ability + err_test7_r
test7_r <- pnorm(test7_r_latent, mean = 0, sd = (gamma_ability^2*var_err+1)^0.5) # uniform distributed -> fits to data
err_test11_r <- rnorm(n = n, mean = 0, sd = var_err^0.5)
test11_r_latent <- test7_r_latent + err_test11_r
test11_r <- pnorm(test7_r_latent, mean = 0, sd = (gamma_ability^2*var_err+ 1 + 1)^0.5) # uniform distributed -> fits to data

breaks_test7_m = c(0, 0.141, 0.158, 0.185, 0.190, 0.212) #function input
breaks_test7_m = breaks_test7_m / sum(breaks_test7_m)
breaks_test11_m = c(0, 0.122,0.152, 0.157, 0.179, 0.199) #function input
breaks_test11_m = breaks_test11_m / sum(breaks_test11_m)
breaks_test7_r = c(0, 0.166, 0.179, 0.188, 0.187, 0.165) #function input
breaks_test7_r = breaks_test7_r / sum(breaks_test7_r)
breaks_test11_r = c(0, 0.132, 0.163, 0.163, 0.176, 0.176) #function input
breaks_test11_r = breaks_test11_r / sum(breaks_test11_r)
test_cat = TRUE
if (isTRUE(test_cat)) {
  test7_m = as.numeric(cut(test7_m, breaks = cumsum(breaks_test7_m))) #categorize the test results and make them numeric (cut generates factors)
  test11_m = as.numeric(cut(test11_m, breaks = cumsum(breaks_test11_m)))
  test7_r = as.numeric(cut(test7_r, breaks = cumsum(breaks_test7_r)))
  test11_r = as.numeric(cut(test11_r, breaks = cumsum(breaks_test11_r)))
}

#number of siblings ->
gamma_parent <- 2 #function input
max_yearsSchooling <- 20 #maybe find convenient alternative
prob_numbSiblings <- 0.1 #func input: additional year of schooling reduces the expected number of children by 0.1
transformed_parentEduc <- max_yearsSchooling - parent_educ #to make it feasible with theory
transformed_parentEduc <- replace(transformed_parentEduc, transformed_parentEduc < 0, 0)
mean_numberSiblings <- 1.692 #func input
var_numberSiblings <- 2.5 #func input
expec_var <- (max_yearsSchooling - mean_parent_educ)*prob_numbSiblings*(1-prob_numbSiblings)
var_expec <- prob_numbSiblings^2 * variance_parent_educ
eta_err_numbSiblings <- 1 - sqrt(mean_numberSiblings - (max_yearsSchooling - mean_parent_educ)*prob_numbSiblings) / sqrt(var_numberSiblings - expec_var - var_expec)
lambda_err_numbSiblings <- (mean_numberSiblings - (max_yearsSchooling - mean_parent_educ)*prob_numbSiblings)*(1-eta_err_numbSiblings)
err_numbSiblings <- rgenPois(n, lambda = lambda_err_numbSiblings, eta = eta_err_numbSiblings)$rnumb
numb_Siblings_first <- sapply(transformed_parentEduc, rbinom_wrapper, n = 1, p = prob_numbSiblings)
numb_Siblings <- numb_Siblings_first + err_numbSiblings
var(numb_Siblings)
