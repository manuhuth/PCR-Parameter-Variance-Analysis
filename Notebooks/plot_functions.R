plot_2D_data <- function(sample_size, true_phi, true_X, legend, transform) {
  #input: - sample_size: vector of sample sizes for which the samples should be drawn
  #       - true_phi: matrix of the true eigenvectors of the whole population
  #       - true_X: matrix of the true values
  #       - legend: TRUE or FALSE indicates whether each plot should have a legend
  #output: - plots showing the deviations of the eigenvectors from the true eigenvectors

  X_sample_2D <- X_2D[sample(nrow(true_X), sample_size), ] #draw sample
  pc_sample_2D <- PCA(X_sample_2D, transform = transform) #compute pc of sample
  phi_sample_2D <- pc_sample_2D$phi #get eigenvectors of sample
  if (isTRUE(legend)) {
    ggplot(as.data.frame(pc_sample_2D$X_transformed), aes(parent_educ, schooling)) + geom_point() +
      geom_abline(aes(intercept = 0, slope = true_phi[1,1]/true_phi[2,1], colour = 'true eigenvectors')) +
      geom_abline(aes(intercept = 0, slope = true_phi[1,2]/true_phi[2,2])) +
      geom_abline(aes(intercept = 0, slope = phi_sample_2D[1,1]/phi_sample_2D[2,1], colour = 'first eigenvector')) +
      geom_abline(aes(intercept = 0, slope = phi_sample_2D[1,2]/phi_sample_2D[2,2], colour = 'second eigenvector')) +
      labs(colour = '') +
      scale_colour_manual(values=c('red','green','black')) + ggtitle(paste('N = ', sample_size, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
  } else{
    ggplot(as.data.frame(pc_sample_2D$X_transformed), aes(parent_educ, schooling)) + geom_point() +
      geom_abline(intercept = 0, slope = true_phi[1,1]/true_phi[2,1], colour = 'black') +
      geom_abline(intercept = 0, slope = true_phi[1,2]/true_phi[2,2], colour = 'black') +
      geom_abline(intercept = 0, slope = phi_sample_2D[1,1]/phi_sample_2D[2,1], colour = 'green') +
      geom_abline(intercept = 0, slope = phi_sample_2D[1,2]/phi_sample_2D[2,2], colour = 'red') +
      ggtitle(paste('N = ', sample_size, sep = ""))} +
    theme(plot.title = element_text(hjust = 0.5))
}
