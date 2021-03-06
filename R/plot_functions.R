plot_2D_data <- function(sample_size, true_phi, true_X, legend, transform) {
  #input: - sample_size: vector of sample sizes for which the samples should be drawn
  #       - true_phi: matrix of the true eigenvectors of the whole population
  #       - true_X: matrix of the true values
  #       - legend: TRUE or FALSE indicates whether each plot should have a legend
  #output: - plots showing the deviations of the eigenvectors from the true eigenvectors

  X_sample_2D <- true_X[sample(nrow(true_X), sample_size), ] #draw sample
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

plots_variance <- function(df_1, legend_df1, colour_1 = 'orange', df_2, legend_df2, colour_2 = 'blue',
                           sample_size, xlab_string, ylab_string, numb_columns, numb_plots = NULL, height = c(1, 0.2)) {
  #input: - df_1/2: data frame that consists of the variances to be plotted in one column. One sample size per column
  #       - legend_df1/2: name displayed in the legend
  #       - colour_1/2: colour of line plot
  #       - sample_size: vector of sample sizes in the columns of df_1/2
  #       - x/ylab_string: name of the xlab/ylab
  #       - numb_columns: Number of plots in each row in the final plot
  #       - numb_plots: Number of columns for which the columns should be plotted
  #output: plots of each columns of df_1 and df_2 in one graph

  plot_list <- list()
  if (is.null(numb_plots)) {
    numb_plots <- ncol(df_1)
  }

  for (index in 1:numb_plots) {
    df <- as.data.frame(cbind(df_1[, index], df_2[,index], sample_size))
    colnames(df) <- c('first', 'second', 'sample_size')
    plot_list[[index]] <- ggplot(df, aes(x=sample_size, y=first, colour = legend_df1)) + geom_line() +
      geom_line(aes(x=sample_size, y=second, colour = legend_df2)) + ggtitle(paste('beta', index, sep = "")) +
      xlab(xlab_string) + ylab(ylab_string) + labs(colour = '') +
      scale_colour_manual(values=c(colour_2,colour_1)) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
  }

  legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "bottom")) #cowplot::otherwise cowplot's ggplot theme is loaded
  p_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = numb_columns)
  cowplot::plot_grid(p_grid, legend, ncol = 1, rel_heights = height)
}

#Plots with CI
prepare_variances_plots_CI <- function(list_variances, sample_size, beta_column, confidence, type_CI = 'mean', line = 'mean', numb_it) {
  #input: - list_variances: object created by prepare_variance_plots function
  #       - sample_size: vector containing the sample sizes
  #       - beta_column: index for which beta the plot should be prepared
  #       - confidence: confidence level in which the variance estimate should be
  #       - type_CI/line: if 'mean' the confidence intervals/plotted line of the mean are plotted. If != 'mean' coverage probabilities are plotted
  #       - numb_it: umber of iterations used to compute the means of the standard errors
  #output: - 4 data frames containing the lower/upper bounds and the means of the variances
    low <- (1-confidence)/2
    up <- (1-confidence)/2 + confidence
    store_beta_prac <- c()
    store_beta_prac_formula <- c()
    store_beta_theo <- c()
    store_beta_theo_formula <- c()

    for (index in 1:length(list_variances)) {
      prac_mean <- mean(list_variances[[index]]$variances_beta_prac[,beta_column])
      prac_line <- list_variances[[index]]$variances_beta_prac[1,beta_column]
      if (type_CI == 'mean') {
        prac_sd <- var(list_variances[[index]]$variances_beta_prac[,beta_column])^0.5 / numb_it^0.5
        prac_up <- prac_mean + qt(up, numb_it-1)*prac_sd
        prac_low <- prac_mean - qt(up, numb_it -1)*prac_sd
      } else{
        prac_up <- quantile(list_variances[[index]]$variances_beta_prac[,beta_column], up)
        prac_low <- quantile(list_variances[[index]]$variances_beta_prac[,beta_column], low)
      }


      prac_formula_mean <- mean(list_variances[[index]]$variances_beta_prac_formula[,beta_column])
      prac_formula_line <- list_variances[[index]]$variances_beta_prac_formula[1,beta_column]
      if (type_CI == 'mean') {
        prac_formula_sd <- var(list_variances[[index]]$variances_beta_prac_formula[,beta_column])^0.5 /numb_it^0.5
        prac_formula_up <- prac_formula_mean + qt(up, numb_it-1)*prac_formula_sd
        prac_formula_low <- prac_formula_mean - qt(up, numb_it-1)*prac_formula_sd
      } else{
      prac_formula_up <-quantile(list_variances[[index]]$variances_beta_prac_formula[,beta_column], up)
      prac_formula_low <- quantile(list_variances[[index]]$variances_beta_prac_formula[,beta_column], low)
      }


      theo_mean <- mean(list_variances[[index]]$variances_beta_theo[,beta_column])
      theo_line <- list_variances[[index]]$variances_beta_theo[1,beta_column]
      if (type_CI == 'mean') {
        theo_sd <- var(list_variances[[index]]$variances_beta_theo[,beta_column])^0.5 /numb_it^0.5
        theo_up <- theo_mean + qt(up, numb_it-1)*theo_sd
        theo_low <- theo_mean - qt(up, numb_it-1)*theo_sd
      } else{
        theo_up <- quantile(list_variances[[index]]$variances_beta_theo[,beta_column], up)
        theo_low <- quantile(list_variances[[index]]$variances_beta_theo[,beta_column], low)
      }

      theo_formula_mean <- mean(list_variances[[index]]$variances_beta_theo_formula[,beta_column])
      theo_formula_line <- list_variances[[index]]$variances_beta_theo_formula[1,beta_column]
      if (type_CI == 'mean') {
        theo_formula_sd <- var(list_variances[[index]]$variances_beta_theo_formula[,beta_column])^0.5 /numb_it^0.5
        theo_formula_up <- theo_formula_mean + qt(up, numb_it-1)*theo_formula_sd
        theo_formula_low <- theo_formula_mean - qt(up, numb_it-1)*theo_formula_sd
      } else{
      theo_formula_up <-quantile(list_variances[[index]]$variances_beta_theo_formula[,beta_column], up)
      theo_formula_low <- quantile(list_variances[[index]]$variances_beta_theo_formula[,beta_column], low)
      }

      plot_line_prac <- prac_mean
      plot_line_prac_formula <- prac_formula_mean
      plot_line_theo <- theo_mean
      plot_line_theo_formula <- theo_formula_mean

      if (line != 'mean') {
        plot_line_prac <- prac_line
        plot_line_prac_formula <- prac_formula_line
        plot_line_theo <- theo_line
        plot_line_theo_formula <- theo_formula_line
      }

      if (line == 'none') {
        plot_line_prac <-  rep(NaN, length(prac_mean))
        plot_line_prac_formula <- rep(NaN, length(prac_mean))
        plot_line_theo <- rep(NaN, length(prac_mean))
        plot_line_theo_formula <- rep(NaN, length(prac_mean))
      }

      store_beta_prac <- rbind(store_beta_prac, cbind(prac_up, plot_line_prac, prac_low))
      store_beta_prac_formula <- rbind(store_beta_prac_formula, cbind(prac_formula_up, plot_line_prac_formula, prac_formula_low))
      store_beta_theo <- rbind(store_beta_theo, cbind(theo_up,plot_line_theo, theo_low))
      store_beta_theo_formula <- rbind(store_beta_theo_formula, cbind(theo_formula_up,plot_line_theo_formula, theo_formula_low))

    }
    colnames(store_beta_prac_formula) <- c('up','mean', 'low')
    colnames(store_beta_prac) <- c('up','mean', 'low')
    colnames(store_beta_theo_formula) <- c('up','mean', 'low')
    colnames(store_beta_theo) <- c('up','mean', 'low')
    list_return <- list('beta_prac_formula' = store_beta_prac_formula, 'beta_prac' = store_beta_prac,
                        'beta_theo_formula' = store_beta_theo_formula, 'beta_theo' = store_beta_theo)
    return(list_return)
}

plots_variance_CI <- function(variance_plots_data, first_method, second_method, colour_1, colour_2,
                              xlab_string, ylab_string, legend_df1, legend_df2, height = c(1,0.2),
                              numb_columns, numb_plots = NULL, sample_size, position_ylab = 0.6, subtitle, make_subtitle = TRUE) {
  #input: - variance_plots_data: list created by prepare_variances_plots_CI
  #       - first/second_method: either c(beta_prac, beta_prac_formula, beta_theo, beta_theo_formula)
  #       - colour_1/2: colour of line plot
  #       - sample_size: vector of sample sizes in the columns of df_1/2
  #       - x/ylab_string: name of the xlab/ylab
  #       - legend_df1/2: name displayed in the legend
  #       - numb_columns: Number of plots in each row in the final plot
  #       - numb_plots: Number of columns for which the columns should be plotted
  #       - height: scaling of the plot window
  #       - sample_size: vector containing the sample size
  #output: plots of each columns of the data in variance_plots_data
  plot_list <- list()

  for (index in 1:length(variance_plots_data)) {
    df1 <- as.data.frame(variance_plots_data[[index]][first_method])
    colnames(df1) <- c('up', 'mean', 'low')
    df2 <- as.data.frame(variance_plots_data[[index]][second_method])
    colnames(df2) <- c('up', 'mean', 'low')
    df <- as.data.frame(cbind(df1[,'up'],df1[,'mean'], df1[,'low'], df2[,'up'],df2[,'mean'], df2[,'low'], sample_size))

    colnames(df) <- c('up1', 'mean1', 'low1', 'up2', 'mean2', 'low2', 'sample_size')
    subt <- ggtitle(" ")
    if (isTRUE(make_subtitle)) {
      subt <- ggtitle(paste(subtitle, index, sep = ""))
    }
    plot_list[[index]] <- ggplot(data = df, aes(sample_size, up1)) +
      geom_ribbon(data=df,aes(ymin=low1,ymax=up1),alpha=0.3, fill = colour_1) +
      geom_line(data = df, aes(sample_size, mean1,colour = legend_df1)) +
      geom_ribbon(data=df,aes(ymin=low2,ymax=up2),alpha=0.3, fill = colour_2) +
      geom_line(data = df, aes(sample_size, mean2,colour = legend_df2)) +
      subt +
      xlab(xlab_string) + ylab(' ') + labs(colour = '') +
      scale_colour_manual(values=c(colour_1, colour_2)) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
  }
  legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "bottom")) #cowplot::otherwise cowplot's ggplot theme is loaded
  p_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = numb_columns)
  cowplot::plot_grid(p_grid, legend, ncol = 1, rel_heights = height) +
    cowplot::draw_label(ylab_string, x=  0, y=position_ylab, vjust= 1, angle=90)

}

prepare_Y_variances <- function(meth_interest, list_variances, confidence, numb_it, type_CI = 'mean',line = 'mean') {
  #input: - meth_interest: method to call from 'list variances'
  #       - list_variances: object created by prepare_variance_plots function
  #       - confidence: confidence level in which the variance estimate should be
  #       - type_CI/line: if 'mean' the confidence intervals/plotted line of the mean are plotted. If != 'mean' quantiles are plotted
  #output: - a matrix containing the lower/upper bounds and the means of the variances
  store_line <- c()
  store_mean <- c()
  store_up <- c()
  store_low <- c()
  low <- (1-confidence)/2
  up <- (1-confidence)/2 + confidence
  for (index in 1:length(list_variances)) {
    help_object <- list_variances[[index]][[meth_interest]]
    sd <- var(help_object)^0.5/numb_it^0.5
    mean <-  mean(help_object)
    store_mean <- rbind(store_mean, mean)
    store_line <- rbind(store_line, help_object[1])
    if (type_CI == 'mean') {
      store_up <- rbind(store_up, mean + qt(up, numb_it-1)*sd )
      store_low <- rbind(store_low, mean - qt(up, numb_it-1)*sd )
    } else {
      store_up <- rbind(store_up, quantile(help_object, up))
      store_low <- rbind(store_low, quantile(help_object, low))
    } #end if/else
  } #end index

  plot_line <- store_mean
  if (line != 'mean') {
    plot_line <- store_line
  }

  if (line == 'none') {
    plot_line <- rep(NaN, length(store_mean))
  }

  matrix_return <- cbind(store_up, plot_line, store_low)
  colnames(matrix_return) <- c('up', 'mean', 'low')
  return(matrix_return)
}
