# utils.R
r_summary <- function(data, lev = NULL, model = NULL) {
    pred_sd <- sd(data$pred, na.rm = TRUE)
    obs_sd <- sd(data$obs, na.rm = TRUE)

    # Return NA if predictions or observations have issues
    if (any(is.na(data$pred)) || is.na(pred_sd) || pred_sd == 0 || is.na(obs_sd) || obs_sd == 0) {
      return(c(r = NA, RMSE = NA, Rsquared = NA, MAE = NA, rho = NA))
    }
    
    final_metrics <- c(
      RMSE = sqrt(mean((data$pred - data$obs)^2)),
      MAE = mean(abs(data$pred - data$obs)),
      r = cor(data$pred, data$obs, method = "pearson"),
      rho = cor(data$pred, data$obs, method = "spearman"),
      Rsquared_ss = 1 - sum((data$obs - data$pred)^2) /
        sum((data$obs - mean(data$obs))^2))

    final_metrics["Rsquared"] <- final_metrics["r"]^2
   
    return(final_metrics)
}