assumptions <- function(x) {
  descr <- as.data.frame(matrix(NA, nrow = 7, ncol = ncol(x)))
  rownames(descr) <- c("Number_of_Observations",
                       "Number_of_missing_values",
                       "min_value", "max_value",
                       "median_value","_skewness_", "_kurtosis_")
  descriptives <- pastecs::stat.desc(x) #calculate descriptive statistics.
  
  ##This function taken from https://www.r-bloggers.com/computing-the-mode-in-r/
  #Mode = function(x) {
  #  ta = table(x)
  #  tam = max(ta)
  #  if (all(ta == tam))
  #    mod = NA
  #  else
  #    if(is.numeric(x))
  #      mod = as.numeric(names(ta)[ta == tam])
  #  else
  #    mod = names(ta)[ta == tam]
  #  return(mod)
  #}
  
  #calculate modes
  descr[1:4, ] <- descriptives[c(1, 3, 4, 5), ]
  descr[5, ] <- descriptives[8, ]
  descr[6, ] <- moments::skewness(x, na.rm = T)
  descr[7, ] <- moments::kurtosis(x, na.rm = T)-3
  
  
  #Calculate VIF and Tolerance values
  #To obtain IF and TV values describe the model
  x_new <- x
  x_new$rn <- 1:nrow(x)
  model_for_collinearity <- lm(
    as.formula(paste(colnames(x_new)[ncol(x_new)], "~",
                     paste(colnames(x_new)[1:(ncol(x_new)-1)], collapse = "+"),
                     sep = ""
    )), data = x_new)
  mc_VIF_TOL <- as.data.frame(mctest::mctest(model_for_collinearity,
                                             type = "i")$idiags[,1:2]) #calculate VIF and Tollerance values
  #Calculate Condition Index
  mc_CI <- mctest::eigprop(mod = model_for_collinearity)$ci
  #A data frame for summary of multicollinearity
  mc_control <- data.frame(min_VIF = min(mc_VIF_TOL$VIF),
                           max_VIF = max(mc_VIF_TOL$VIF),
                           min_TOL = min(mc_VIF_TOL$TOL),
                           max_TOL = max(mc_VIF_TOL$TOL),
                           min_CI = min(mc_CI),
                           max_CI = max(mc_CI)) #giving a summary of multicollinearity
  
  #Mahalanobis Distance Calculation
  #To calculate mahalanobis distance, missing values are not accepted.
  distance <- as.matrix(mahalanobis(x, colMeans(x), cov = cov(x)))
  #Those with Mahalanobis Distance p values bigger than 0.001 were considered as outliers.
  Mah_significant <- x %>%
    transmute(row_number = 1:nrow(x),
              Mahalanobis_distance = distance,
              Mah_p_value = pchisq(distance, df = ncol(x), lower.tail = F)) %>%
    filter(Mah_p_value <= 0.001)
  #Calculate Mardia's kurtosis value for multivariate normality
  mardia_kurt <- psych::mardia(x, plot = F)
  
  #Return a list consist of descriptive statistics, multicollinearity, multivariate normality and Mahalanobis distance for multivariate outliers
  return(list(descriptives = round(descr, 2),
              multicollinearity = round(mc_control, 2),
              Mah_significant = Mah_significant,
              n_outlier = nrow(Mah_significant),
              Mardia_Kurtosis = mardia_kurt$kurtosis,
              Mardia_Kurtosis_p_value = mardia_kurt$p.kurt )) }
