library(dagitty)
causal_graph <- dagitty("dag {
Ability -> TestScores7thGrade Ability Ability -> TestScores11thGrade Ability -> YearsSchooling Ability -> Wage
TestScores7thGrade -> TestScores11thGrade
ParentEducation -> TestScores7thGrade ParentEducation -> TestScores11thGrade ParentEducation -> NumberSiblings ParentEducation -> YearsSchooling
YearsSchooling -> WorkExperience YearsSchooling -> Wage
NumberSiblings -> Wage ParentEducation -> Wage
WorkExperience -> Wage
Age -> Wage
}")
plot(graphLayout(causal_graph))


library(gridExtra)
library(ggplot2)
library(readxl)
library(tidyr)
library(matrixStats)
setwd('C:/Users/Mhuth/Desktop/PCRPVA')

files <- c('PCA_PropVar', 'PCA', 'PCR_cv', 'PCR_predict', 'PCR', 'random_discreteVariables', 'random_VCV', 'DGP_function', 'plot_functions','coef_var_analysis', 'OLS_estimate') #names of files to read
for (i in 1:length(files)) { #loop to read all files
  source(paste('R/', files[i], '.R', sep = ''))
}



trans <- 'normalized'
#trans <- 'standardized'

load('SimData/population_230000.Rda')

set.seed(123)
#population <- dgp_model(n = 230000)
N <- nrow(population)
X <- population[c('test7_m', 'test11_m', 'test7_r', 'test11_r', 'parent_educ', 'schooling', 'numb_Siblings', 'working')]
X <- cbind(X, X$working^2/100)
colnames(X) <- c('test7_m', 'test11_m', 'test7_r', 'test11_r', 'parent_educ', 'schooling', 'numb_Siblings', 'working', 'working_squ')

#histograms
population_names <- population
colnames(population_names) <- c('Logwage', 'Ability', '7th Grade Math', '11th Grade Math', '7th Grade Reading', '11th Grade Reading',
                                   'Parent Education', 'Years of Schooling', 'Number of Siblings', 'Years of Working', 'Age')
population_names <- population
colnames(population_names) <- c('Logwage', 'Ability', '7th Grade Math', '11th Grade Math', '7th Grade Reading', '11th Grade Reading',
                                'Parent Education', 'Years of Schooling', 'Number of Siblings', 'Years of Working', 'Age')
gathered_population <- gather(population_names)
ggplot((gathered_population), aes(x = value, y = ..density..) ) +
  geom_histogram(binwidth = 1, color = 'steelblue3', fill = "steelblue3" ) +  facet_wrap(~key, scales = 'free')

#deviations
imposed_means <- c('logwage' = 2.040, 'Parent Ed.' = 13.342, 'Scholing' = 13.342, 'Siblings' = 1.692)  #vector of imposed means
imposed_variances <- c('logwage' = 1.5,'Parent Ed.' = 21.215, 'Scholing' = 21.215, 'Siblings' = 2.89) #vector of imposed variances

check_moments <- population[c('logwage', 'parent_educ','schooling', 'numb_Siblings')] #get variables for which a moment was imposed
observed_means <- colMeans(check_moments) #calculate simulatd means
observed_variances <- colVars(as.matrix(check_moments)) # calculate simulated variances
'Relative Deviations from the Imposed Means'
(deviations_mean <- round(t(imposed_means - observed_means)/imposed_means, 5))
'Relative Deviations from the Imposed Variances'
(deviations_vars <- round(t(imposed_variances - observed_variances)/imposed_variances, 5))



#transform data
X_stand <- scale(X, center = TRUE, scale = FALSE)
VCV <- t(X_stand)%*%X_stand/N
H <- diag(diag(VCV)^(-0.5))
colnames(H) <-   c('test7_m', 'test11_m', 'test7_r', 'test11_r', 'parent_educ', 'schooling', 'numb_Siblings', 'working', 'working_squ')
rownames(H) <-   c('test7_m', 'test11_m', 'test7_r', 'test11_r', 'parent_educ', 'schooling', 'numb_Siblings', 'working', 'working_squ')
corr <-  H %*% VCV %*% H # H %*% \Sigma %*% H, as in equation a.23

phi <- eigen(VCV)$vectors
lambda <- eigen(VCV)$values



  #4.4 Variance of the estimated coefficients

############################################################################################################
#with confidence intervals
sample_size <- c(c(50,80,100),150, seq(200, 5000, by = 200))
set.seed(1234)
population_Y <- cbind(population$logwage, X)
M <- 9 # do it for all since it does not change the other
numb_it <- 400
#variance_analysis_results <- lapply(sample_size, coef_var_analysis, population = population_Y, true_phi = phi, M = M,
                                   # transform = trans, iterations = numb_it, numb_iterations = numb_it)
load("SimData/var_9M.Rda")
type_CI <- 'quantile' #either 'mean' or 'empirical'
line <- 'mean' #either 'mean' or 'empirical'
variance_plots_data_CI <- lapply(c(1:M), prepare_variances_plots_CI, list_variances = variance_analysis_results,
                                 sample_size = sample_size, confidence = 0.99, type_CI = type_CI, line = line, numb_it = numb_it)
xlab <- 'Sample Size'
ylab <- 'Mean of Variance Estimate'
numb_columns <- 3
col_prac_formula <- 'orange'
col_prac <- 'blue'
col_theo_formula <- 'green'
col_theo <- 'red'
xlab_string <- xlab
ylab_string <- ylab

#1
plots_variance_CI(variance_plots_data = variance_plots_data_CI, first_method = 'beta_prac_formula', second_method = 'beta_prac',
                  colour_1 = col_prac_formula, colour_2 = col_prac, xlab_string = xlab, ylab_string = ylab,
                  legend_df1 = 'Mean of Variance Stochastic Phi',
                  legend_df2 = 'Mean of Variance with Stochastic Phi (Empirical)',
                  numb_columns = numb_columns, sample_size = sample_size, subtitle = 'beta')

#2
plots_variance_CI(variance_plots_data = variance_plots_data_CI, first_method = 'beta_theo_formula', second_method = 'beta_theo',
                  colour_1 = col_theo_formula, colour_2 = col_theo, xlab_string = xlab, ylab_string = ylab,
                  legend_df1 = 'Mean of Variance with True Phi',
                  legend_df2 = 'Mean of Variance with True Phi (Empirical Distribution)',
                  numb_columns = numb_columns, sample_size = sample_size, subtitle = 'beta')

#3
plots_variance_CI(variance_plots_data = variance_plots_data_CI, first_method = 'beta_theo_formula', second_method = 'beta_prac_formula',
                  colour_1 = col_theo_formula, colour_2 = col_prac_formula, xlab_string = xlab, ylab_string = ylab,
                  legend_df1 = 'Mean of Variance with True Phi',
                  legend_df2 = 'Mean of Variance with Stochastic Phi',
                  numb_columns = numb_columns, sample_size = sample_size, subtitle = 'beta')
##################
#Variance of Y_hat
confidence <- 0.99

variances_Y_prac <- prepare_Y_variances(meth_interest = 'variances_Y_prac', list_variances = variance_analysis_results, confidence = confidence,
                                        type_CI = type_CI, line = line, numb_it = numb_it)
variances_Y_theo <- prepare_Y_variances(meth_interest = 'variances_Y_theo', list_variances = variance_analysis_results, confidence = confidence,
                                        type_CI = type_CI, line = line, numb_it = numb_it)
list_var_Y <- list(list('prac' = variances_Y_prac, 'theo' = variances_Y_theo))

plots_variance_CI(variance_plots_data = list_var_Y, first_method = 'theo', second_method = 'prac',
                  colour_1 = 'orange', colour_2 = 'blue', xlab_string = 'Sample Size', ylab_string = 'Variance of Y Estimator',
                  legend_df1 = 'Mean of Variance with Non-Stochastic (True) Phi',
                  legend_df2 = 'Mean of Variance with Stochastic Phi',
                  numb_columns = 1, sample_size = sample_size, subtitle = 'Y', make_subtitle = FALSE)


########################check robustness. 6 and 3 principal components

