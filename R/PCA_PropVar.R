PCA_PropVar <- function(pca, numb_components = NULL, barplot = FALSE, main_plot = NULL) {
  #input: - pca: object of class PCA, created by the PCA function
  #       - numb_components: number for principal components for which the proportional variance should be computed
  #       - barplot: boolean, if TRUE a barplot of the proportional variance is plotted
  #       - main_plot: main title of the barplot.
  #output: A list containing
  #       - var_proportion: vector containing the proportional variances
  if ( isFALSE(is(pca, 'PCA')) ){
    stop('The input `pca` must be an object which is an output from the PCA function.')
  }
  
  if (is.null(numb_components)) {
    numb_components = length(pca$eigenval)
  }
  
  var_proportion <- pca$eigenval[1:numb_components]/sum(pca$eigenval[1:numb_components])
  if (isTRUE(barplot)) {
    if (is.null(main_plot)) {
      main_plot = 'Proportion of the Variance that is captured by Z_m'
    }
    
    barplot(var_proportion, main = main_plot, xlab = 'Z_m', ylab = 'Proportion',
            ylim = c(0, max(var_proportion*1.2, 1) ), col = 'steelblue' )
  }
  
  return(var_proportion)
}

