population <- load("C:/Users/Mhuth/Desktop/PCRPVA/SimData/population.RData")

X <- population[c('test7_m', 'test11_m', 'test7_r', 'test11_r', 'parent_educ', 'schooling', 'working')]
X <- cbind(X, X$working^2/100)
colnames(X) <- c(c('test7_m', 'test11_m', 'test7_r', 'test11_r', 'parent_educ', 'schooling','working', 'working_squ'))
#note that scale uses (n-1) degrees of freedom, even though there are actually n in this case. However, dividing by (n-1) in the subsequent line cancels this out.
population_stand <- scale(population, center = TRUE, scale = TRUE)
VCV <- t(population_stand)%*%population_stand/(n-1) #since all variables are standardized, this equals the Pearson-Bravais correlations.

#Table note
VCV

#compute true eigenvectors
phi <- eigen(VCV)$vectors
lambda <- eigen(VCV)$values

#4.3 Convergence of eigenvectors in two dimensions -> since nice to illustrate

