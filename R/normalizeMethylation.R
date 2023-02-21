#' Title Normalize Methylation Values
#'
#' @param df a data frame containing methylation, position, and coverage values.
#' @param alpha
#' @param normalize_position
#' @param rescale
#' @param plots
#'
#' @return
#' @export
#'
#' @examples
normalizeMethylation <- function(df,
                                 alpha = 0.001,
                                 normalize_position = TRUE,
                                 rescale = FALSE,
                                 plots = TRUE) {
  # TODO: Change this to act on one sample at once, then write a second function
  #       as a wrapper to do multiple samples at once.
  # df must have columns 'Position' and coverage starting with 'cov' and betas
  # starting with 'methyl'
  df_noNA <- df %>% tidyr::drop_na()
  coverage_df <- df_noNA %>%
    dplyr::select(starts_with('cov'))
  if (normalize_position == TRUE) {
    for (i in 1:ncol(coverage_df)) {
      coverage_df[[i]] <- detrendCoverage(coverage = coverage_df[[i]],
                                          position = df_noNA[[2]])
    }
  }
  methyl_colnames <- df_noNA %>%
    dplyr::select(starts_with('methyl')) %>%
    colnames()
  methyl_df <- df_noNA %>%
    dplyr::select(starts_with('methyl')) %>%
    as.matrix()
  # Convert methylation beta values to number of methylated reads
  methyl_df_reads <- round(as.matrix(coverage_df)*methyl_df*100, digits = -2)/100
  # Normalize coverage by quantile for each sample
  coverage_df_quantnormalized <- coverage_df %>%
    as.matrix() %>%
    preprocessCore::normalize.quantiles(x = .)
  colnames(coverage_df_quantnormalized) <- colnames(coverage_df)
  colnames(methyl_df_reads) <- methyl_colnames
  # For each sample, fit lm to normalized coverage ~ original coverage,
  # then extract the intercept and slope from the lm object,
  # then use formula y = a + bx to normalize each sample's # methylated reads
  transformabx <- function(a, b, x) {a + b*x}
  coverage_df_normalized <- matrix(data = numeric(length = 1),
                                   nrow = nrow(coverage_df),
                                   ncol = ncol(coverage_df))
  colnames(coverage_df_normalized) <- colnames(coverage_df)
  methyl_df_normalized <- matrix(data = numeric(length = 1),
                                 nrow = nrow(methyl_df),
                                 ncol = ncol(methyl_df))
  colnames(methyl_df_normalized) <- methyl_colnames
  models <- vector('list', ncol(methyl_df))
  if (plots == TRUE) {
    model_plots <- vector('list', ncol(methyl_df))
    cov_plots <- vector('list', ncol(methyl_df))
    methyl_plots <- vector('list', ncol(methyl_df))
    methyl_read_plots <- vector('list', ncol(methyl_df))
  }
  r_squared_cov <- vector('list', ncol(methyl_df))
  r_squared_methyl <- vector('list', ncol(methyl_df))
  # Shows RMSE between beta==norm_beta (theoretical "perfect" transformation)
  # and actual best fit line
  rmse_list <- vector('list', ncol(methyl_df))

  for (i in 1:ncol(methyl_df)) {
    fit <- lm(coverage_df_quantnormalized[,i] ~ as.matrix(coverage_df)[,i])
    models[[i]] <- fit
    r_squared_cov[[i]] <- summary(fit)$r.squared
    a <- fit$coefficients[1]
    b <- fit$coefficients[2]
    coverage_df_normalized[,i] <- transformabx(a, b, as.matrix(coverage_df)[,i])
    methyl_df_normalized[,i] <-
      transformabx(a, b, methyl_df_reads[,i]) /
      (coverage_df_normalized[,i] + alpha)
    rmse_fit <- lm(methyl_df_normalized[,i] ~ methyl_df[,i])
    rmse_pred <- seq(0, 1, 1 / (length(methyl_df[,i]) - 1))
    rmse_list[[i]] <- rmse(rmse_fit$coefficients[2]*rmse_pred + rmse_fit$coefficients[1],
                           rmse_pred)
    r_squared_methyl[[i]] <- summary(rmse_fit)$r.squared
    if (rescale == TRUE) {
      methyl_df_normalized[,i] <- (methyl_df_normalized[,i] - rmse_fit$coefficients[1]) / rmse_fit$coefficients[2]
      methyl_df_normalized[methyl_df_normalized[,i] < 0, i] <- 0
      methyl_df_normalized[methyl_df_normalized[,i] > 1, i] <- 1
    }
    if (plots == TRUE) {
      model_plots[[i]] <- autoplot(fit)
      cov_plots[[i]] <- ggplot(data = data.frame(coverage = as.matrix(coverage_df)[,i],
                                                 norm_coverage = coverage_df_normalized[,i]),
                               aes(x = coverage,
                                   y = norm_coverage)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_smooth(method = 'lm', color = 'red')
      methyl_plots[[i]] <- ggplot(data = data.frame(beta = methyl_df[,i],
                                                    norm_beta = methyl_df_normalized[,i]),
                                  aes(x = beta,
                                      y = norm_beta)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_smooth(method = 'lm', color = 'red')
      methyl_read_plots[[i]] <- ggplot(data = data.frame(methyl_reads = methyl_df_reads[,i],
                                                         norm_methyl_reads = transformabx(a, b, methyl_df_reads[,i])),
                                       aes(x = methyl_reads,
                                           y = norm_methyl_reads)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_smooth(method = 'lm', color = 'red')
    }
  }
  suppressMessages(
    out_df <- bind_cols(df_noNA$Position, as.data.frame(coverage_df_normalized), as.data.frame(methyl_df_normalized))
  )
  colnames(out_df)[1] <- 'Position'
  if (plots == TRUE) {
    return(list(data = out_df,
                models = models,
                model_plots = model_plots,
                cov_plots = cov_plots,
                r_squared_cov = r_squared_cov,
                r_squared_methyl = r_squared_methyl,
                methyl_plots = methyl_plots,
                methyl_read_plots = methyl_read_plots,
                RMSE = rmse_list))
  } else {
    return(list(data = out_df,
                models = models,
                r_squared_cov = r_squared_cov,
                r_squared_methyl = r_squared_methyl,
                RMSE = rmse_list))
  }

}
