#' Estimating the group average potential outcome
#' 
#' IPW estimator of the group average potential outcome.
#'
#' @param dta Data frame including treatment, outcome and covariates.
#' @param cov_cols The indices including the covariates of the ps model.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the ps, and the second one is the random effect variance.
#' If the second element is 0, the propensity score excludes random effects.
#' @param gamma_numer The coefficients of the ps model in the numerator.
#' If left NULL and estimand is 1, the coefficients in phi_hat will be used
#' instead.
#' @param alpha The values of alpha for which we want to estimate the group
#' average potential outcome.
#' @param neigh_ind List. i^{th} element is a vector with the row indices of
#' dta that are in cluster i. Can be left NULL.
#' @param trt_col If the treatment is not named 'V' in dta, specify the
#' treatment column index.
#' @param out_col If the outcome is not named 'Y', specify the outcome column
#' index.
#' @param alpha_re_bound The lower and upper end of the values for bi we will
#' look at. Defaults to 10, meaning we will look between - 10 and 10.
#' @param integral_bound The number of standard deviations of the random effect
#' that will be used as the lower and upper limit.
#' @param keep_re_alpha Logical. If set to TRUE the "random" effect that makes
#' the average probability of treatment equal to alpha will be returned along
#' with the estimated group average potential outcome.
#' @param estimand Character, either '1' or '2.' If 1 is specified, then the
#' estimand with numerator depending on covariates is estimated. If estimand
#' is set equal to 2, the numerator considered is the product of Bernoulli.
#' @param verbose Whether printing of progress is wanted. Defaults to TRUE.
#' 
#' @export

CalcNumerator <- function(Ai_j, Xi_j, gamma_numer, alpha, re_alpha,
                          include_alpha = TRUE) {
  
  gamma_numer <- matrix(gamma_numer, nrow = length(gamma_numer), ncol = 1)
  
  lin_pred <- cbind(1, as.matrix(Xi_j)) %*% gamma_numer
  lin_pred <- lin_pred + re_alpha
  probs <- expit(lin_pred)
  
  if (include_alpha) {
    r <- (probs / alpha) ^ Ai_j * ((1 - probs) / (1 - alpha)) ^ (1 - Ai_j)
  } else {
    r <- probs ^ Ai_j * (1 - probs) ^ (1 - Ai_j)
  }
  
  return(list(prob = prod(r), re_alpha = re_alpha))
}

expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

  
Denominator <- function(A, X, phi_hat, alpha = NULL, integral_bound = 10,
                        include_alpha = TRUE) {
  
  integral_bound <- abs(integral_bound)
  if (is.null(alpha) & include_alpha) {
    stop('No alpha provided.')
  }
  
  X <- as.matrix(cbind(1, X))
  re_sd <- sqrt(phi_hat[[2]])
  phi_hat[[1]] <- matrix(phi_hat[[1]], nrow = length(phi_hat[[1]]), ncol = 1)
  
  # Creating the function that we will integrate over.
  f_int <- function(b) {
    r <- 1
    lin_pred <- X %*% phi_hat[[1]]  # Includes intercept.
    for (ii in 1:length(A)) {
      prob_trt <- expit(lin_pred[ii] + b)
      if (include_alpha) {
        success_weight <- prob_trt / alpha
        failure_weight <- (1 - prob_trt) / (1 - alpha)
      } else {
        success_weight <- prob_trt
        failure_weight <- 1 - prob_trt
      }
      r <- r * success_weight ^ A[ii] * failure_weight ^ (1 - A[ii])
    }
    if (re_sd > 0) { # If re_sd = 0, there is no random effect in the ps model.
      r <- r * dnorm(b, mean = 0, sd = re_sd)
    }
    return(r)
  }
  
  if (re_sd > 0) {
    ans <- integrate(f_int, lower = - integral_bound * re_sd,
                     upper = integral_bound * re_sd)
  } else {
    ans <- list(value = f_int(0))
  }
  
  return(ans)
}

GroupIPW <- function(dta, cov_cols, phi_hat, gamma_numer = NULL, alpha,
                     neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                     alpha_re_bound = 10, integral_bound = 10,
                     keep_re_alpha = FALSE, estimand = c('1', '2'),
                     verbose = TRUE) {
  
  estimand <- match.arg(estimand)
  integral_bound <- abs(integral_bound)
  alpha_re_bound <- abs(alpha_re_bound)
  phi_hat[[1]] <- matrix(phi_hat[[1]], ncol = 1)
  dta <- as.data.frame(na.omit(dta))
  
  # We only return the ksi's if we are estimating estimand 1.
  keep_re_alpha <- keep_re_alpha & (estimand == '1')

  # Specifyling neigh_ind will avoid re-running the following lines.
  if (is.null(neigh_ind)) {
    neigh_ind <- sapply(1 : max(dta$block), function(x) which(dta$block == x))
  }
  
  n_neigh <- length(neigh_ind)
  
  yhat_group <- array(NA, dim = c(n_neigh, 2, length(alpha)))
  dimnames(yhat_group) <- list(neigh = 1:n_neigh, trt = c(0, 1), alpha = alpha)
  
  # Names of treatment and outcome column.
  if (!is.null(trt_col)) {
    names(dta)[trt_col] <- 'V'
  }
  if (!is.null(out_col)) {
    names(dta)[out_col] <- 'Y'
  }
  if (is.null(gamma_numer)) {
    gamma_numer <- matrix(phi_hat[[1]], ncol = 1)
  }
  
  # If we want to return the ksis that make average propensity alpha.
  if (keep_re_alpha) {
    re_alphas <- matrix(NA, nrow = n_neigh, ncol = length(alpha))
    dimnames(re_alphas) <- list(neigh = 1 : n_neigh, alpha = alpha)
  }
  
  for (aa in 1 : length(alpha)) {
    if (verbose) {
      print(paste('alpha =', round(alpha[aa], 3)))
    }
    curr_alpha <- alpha[[aa]]
    
    for (nn in 1 : n_neigh) {
      
      # For estimand 1, we need to calculate numerator depending on covariates.
      if (estimand == '1') {
        
        # Calculating the random effect that gives alpha.
        Xi <- dta[neigh_ind[[nn]], cov_cols]
        lin_pred <- cbind(1, as.matrix(Xi)) %*% gamma_numer
        re_alpha <- FromAlphaToRE(alpha = curr_alpha, lin_pred = lin_pred,
                                  alpha_re_bound = alpha_re_bound)
        
        # Keeping the intercept that makes cluster average propensity alpha.
        if (keep_re_alpha) {
          re_alphas[nn, aa] <- re_alpha
        }
        
      }
      
      for (curr_it in c(0, 1)) {
        
        bern_prob <- curr_alpha ^ curr_it * (1 - curr_alpha) ^ (1 - curr_it)
        prob_ind <- list(prob = 1)  # For estimand 2.
        y_curr <- 0
        
        for (ind in neigh_ind[[nn]]) {
          if (dta$V[ind] == curr_it) {
            
            if (estimand == '1') {
              wh_others <- setdiff(neigh_ind[[nn]], ind)
              Ai_j <- dta$V[wh_others]
              Xi_j <- dta[wh_others, cov_cols]
              prob_ind <- CalcNumerator(Ai_j = Ai_j, Xi_j = Xi_j,
                                        gamma_numer = gamma_numer,
                                        alpha = curr_alpha, re_alpha = re_alpha)
            }
            
            y_curr <- y_curr + dta$Y[ind] * prob_ind$prob
          }
        }
        
        denom <- Denominator(A = dta$V[neigh_ind[[nn]]],
                             X = dta[neigh_ind[[nn]], cov_cols],
                             phi_hat = phi_hat, alpha = curr_alpha,
                             integral_bound = integral_bound)
        denom <- length(neigh_ind[[nn]]) * denom$value * bern_prob
        
        yhat_group[nn, curr_it + 1, aa] <- y_curr / denom
      }
    }
  }
  if (keep_re_alpha) {
    return(list(yhat_group = yhat_group, re_alpha = re_alphas))
  }
  return(list(yhat_group = yhat_group))
}




# Using control data
glmer_fit_col <- lme4::glmer(
  data = datTND[datTND$Y==0, ],
  formula = V ~ C + (1 | block),
  family = stats::binomial,
  nAGQ = 2
)

#gamma (fixed-effects coefficients) and b(variance-covariance parameters for random effects) components of the model.
#phi_hat = list()
#phi_hat$coefs <- lme4::getME(glmer_fit_col, 'gamma') 
#phi_hat$random_effects <- lme4::getME(glmer_fit_col, 'b')
#phi_hat[[1]]

phi_hat <- list(coefs = summary(glmer_fit_col)$coef[, 1], re_var = re_var)


dta <- datTND
estimand = '1'
alpha = 0.5
cov_cols = 'C'


resB <- GroupIPW(dta = datTND, cov_cols = 'C', phi_hat = phi_hat, gamma_numer = NULL, alpha = c(0.2, 0.5),
                     neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                     alpha_re_bound = 10, integral_bound = 10,
                     keep_re_alpha = FALSE, estimand = '1',
                     verbose = TRUE)
resB$yhat_group

resGLMM <- GroupIPW(dta = datTND, cov_cols = 'C', phi_hat = phi_hat, gamma_numer = NULL, alpha = c(0.2, 0.5),
                 neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                 alpha_re_bound = 10, integral_bound = 10,
                 keep_re_alpha = FALSE, estimand = '2',
                 verbose = TRUE)
head(resGLMM$yhat_group)

#############################
#####
#####Score function
#############################

Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = 'C',
          trt_name = NULL, integral_bound = 10)
CalcScore <- function(dta, neigh_ind = NULL, phi_hat, cov_cols,
                      trt_name = NULL, integral_bound = 10) {
  
  dta <- as.data.frame(dta)
  if (is.null(neigh_ind)) {
    n_neigh <- max(dta$neigh)
    neigh_ind <- sapply(1 : n_neigh, function(x) which(dta$neigh == x))
  }
  n_neigh <- length(neigh_ind)
  
  re_var <- phi_hat[[2]]
  num_gamma <- length(phi_hat$coefs) + (re_var > 0)
  phi_hat$coefs <- matrix(phi_hat$coefs, ncol = 1)
  
  if (is.null(trt_name)) {
    trt_name <- 'V'
  }
  trt_col <- which(names(dta) == trt_name)
  
  scores <- matrix(NA, nrow = num_gamma, ncol = n_neigh)
  
  for (nn in 1 : n_neigh) {
    
    Ai <- dta[neigh_ind[[nn]], trt_col]
    Xi <- dta[neigh_ind[[nn]], cov_cols]
    
    hess_function <- function(gamma) {
      
      phi_hat_hess <- NULL
      if (re_var == 0) {
        phi_hat_hess$coefs <- gamma
        phi_hat_hess$re_var <- 0
      } else {
        phi_hat_hess$coefs <- gamma[- num_gamma]
        phi_hat_hess$re_var <- gamma[num_gamma]
      }
      
      likelihood <- Denominator(A = Ai, X = Xi, phi_hat = phi_hat_hess,
                                include_alpha = FALSE,
                                integral_bound = integral_bound)
      return(log(likelihood$value))
    }
    
    hess_x <- as.numeric(phi_hat$coefs)
    if (re_var > 0) {
      hess_x <- c(hess_x, re_var)
    }
    scores[, nn] <- numDeriv::grad(hess_function, x = hess_x)
    
  }
  return(scores)  
}


Ypop.est <- Ypop(ygroup, ps = 'estimated', scores = Score.est,
                 dta = NULL, use = 'everything')
Ypop <- function(ygroup, ps = c('true', 'estimated'), scores = NULL,
                 dta = NULL, use = 'everything') {
  
  ps <- match.arg(ps)
  n_neigh <- dim(ygroup)[1]
  alpha <- as.numeric(dimnames(ygroup)[[3]])
  ypop <- apply(ygroup, c(2, 3), mean, na.rm = TRUE)
  
  ypop_var <- apply(ygroup, 3, function(x) cov(x, use = use))
  ypop_var <- array(ypop_var, dim = c(2, 2, length(alpha)))
  # In order to get 1 / N, instead of 1 / (N - 1) in the variance estimates.
  ypop_var <- ypop_var * (n_neigh - 1) / n_neigh
  ypop_var <- ypop_var / n_neigh  # Since we have n_neigh clusters.
  dimnames(ypop_var) <- list(a = c(0, 1), a = c(0, 1), alpha = alpha)
  
  if (ps == 'true') {
    return(list(ypop = ypop, ypop_var = ypop_var))
  }
  
  
  # Else, the propensity score is estimated.
  
  neigh_ind <- sapply(1 : n_neigh, function(x) which(dta$neigh == x))
  
  if (is.null(scores)) {
    stop('scores needs to be specified for estimated propensity score.')
  }
  
  num_gamma <- dim(scores)[1]
  
  
  var_est_ps <- array(NA, dim = c(2, 2, length(alpha)))
  dimnames(var_est_ps) <- list(it = c(0, 1), it = c(0, 1), alpha = alpha)
  
  # --- Calculating B11, the information matrix of the cluster ps.
  B11 <- matrix(0, nrow = num_gamma, ncol = num_gamma)
  for (nn in 1 : n_neigh) {
    scores_nn <- scores[, nn, drop = FALSE]
    B11 <- B11 + scores_nn %*% t(scores_nn)
  }
  B11 <- B11 / n_neigh
  B11_inv <- chol2inv(chol(B11))
  
  
  # ---- Calculating A21, B12 for each alpha.
  for (aa in 1 : length(alpha)) {
    
    A21 <- array(0, dim = c(2, num_gamma))
    B12 <- array(0, dim = c(num_gamma, 2))
    
    for (it in c(1, 2)) {
      for (nn in 1 : n_neigh) {
        A21[it, ] <- A21[it, ] - ygroup[nn, it, aa] * scores[, nn]
        B12[, it] <- B12[, it] + scores[, nn] * (ygroup[nn, it, aa] -
                                                   ypop[it, aa])
      }
    }
    A21 <- A21 / n_neigh
    B12 <- B12 / n_neigh
    
    chol_B11_inv <- chol(B11_inv)
    mat1 <- A21 %*% t(chol_B11_inv) %*% t(A21 %*% t(chol_B11_inv))
    mat <- A21 %*% B11_inv %*% B12
    
    var_est_ps[, , aa] <- mat1 + mat + t(mat)
    var_est_ps[, , aa] <- var_est_ps[, , aa] / n_neigh
    var_est_ps[, , aa] <- ypop_var[, , aa] + var_est_ps[, , aa]
  }
  
  return(list(ypop = ypop, ypop_var = var_est_ps))
  
}

ygroup = resGLMM$yhat_group
ypop <- apply(ygroup, c(2, 3), mean, na.rm = TRUE)
ps = 'estimated'
  
###################################
## Delta method
delta_method <- function(x, vec)  {
  
  if (length(vec) != dim(x)[1] | length(vec) != dim(x)[2]) {
    stop('Wrong dimensions.')
  }
  
  vec <- matrix(vec, ncol = 1)
  return(t(vec) %*% x %*% vec)
}

#' Direct effect estimates and asymptotic variance.
#' 
#' @param ypop A matrix with rows corresponding to the potential outcome under
#' control and treatment, and columns corresponding to the cluster-average
#' propensity of treatment.
#' @param ypop_var An 3-dimensional array, where the first two dimensions are
#' equal to 2 and include the variance covariance matrix of the population
#' average potential outcome for each alpha. Dimension 3 is alpha.
#' @param boots The results of BootVar() function including estimates of the
#' potential outcomes from the bootstrap samples.
#' @param alpha The values of alpha we consider. If ypop has column names,
#' alpha can be left null.
#' @param alpha_level Numeric. The alpha level of the confidence intervals
#' based on the quantiles of the bootstrap estimates.
#' 
#' @return A matrix with rows including the estimate and variance of the direct
#' effect and columns corresponding to alpha.
#' 
#' @export
#' 

ypop = Ypop.est$ypop
ypop_var = Ypop.est$ypop_var

DE(ypop = Ypop.est$ypop, ypop_var = Ypop.est$ypop_var, boots = NULL, alpha = NULL,
               alpha_level = 0.05)
  
  
DE <- function(ypop, ypop_var, boots = NULL, alpha = NULL,
               alpha_level = 0.05) {
  
  quants <- c(0, 1) + c(1, - 1) * alpha_level / 2
  norm_quant <- - qnorm(alpha_level / 2)
  
  if (is.null(alpha)) {
    if (is.null(colnames(ypop))) {
      stop('Specify alpha.')
    }
    alpha <- as.numeric(colnames(ypop))
  }
  
  dim_names <- c('mu_0','mu_1','DE.est', 'var', 'low_int', 'high_int')
  if (!is.null(boots)) {
    dim_names <- c(dim_names, 'boot_var', 'boot_var_LB', 'boot_var_UB',
                   'boot_low_quant', 'boot_high_quant')
  }
  
  de <- array(NA, dim = c(length(dim_names), length(alpha)))
  dimnames(de) <- list(stat = dim_names, alpha = alpha)
  
  de[1, ] <- ypop[1, ]
  de[2, ] <- ypop[2, ]
  de[3, ] <- ypop[2, ]/ypop[1, ]
  de[4, ] <- sapply(seq_along(ypop[1, ]), function(i) {
    delta_method(ypop_var[, , i], c(1, -1/(ypop[1, ][i])^2))
  })
  de[5, ] <- de[1, ] - norm_quant * sqrt(de[2, ])
  de[6, ] <- de[1, ] + norm_quant * sqrt(de[2, ])
  
  if (!is.null(boots)) {
    de[7, ] <- apply(boots[2, , ] - boots[1, , ], 1, var)
    de[8, ] <- de[1, ] - norm_quant * sqrt(de[5, ])
    de[9, ] <- de[1, ] + norm_quant * sqrt(de[5, ])
    de[10 : 11, ] <- apply(boots[2, , ] - boots[1, , ], 1, quantile,
                         probs = quants)
  }
  
  return(de)
}


IEvar(ygroupV = ygroup[, 1 , ], scores = Score.est)
IEvar <- function(ygroupV, scores){
  n_neigh <- dim(ygroupV)[1]
  alpha <- as.numeric(dimnames(ygroupV)[[2]])
  
  ie_var <- cov(ygroupV)
  ie_var <- ie_var * (n_neigh - 1) / (n_neigh ^ 2)
  
# Based on the estimated propensity score.
  if (is.null(scores)) {
    stop('Provide score matrix.')
  }
  
  var_est_ps <- array(0, dim = dim(ie_var))
  num_gamma <- dim(scores)[1]
  
  # --- Calculating B11, the information matrix of the cluster ps.
  B11 <- matrix(0, nrow = num_gamma, ncol = num_gamma)
  for (nn in 1 : n_neigh) {
    scores_nn <- scores[, nn, drop = FALSE]
    B11 <- B11 + scores_nn %*% t(scores_nn)
  }
  B11 <- B11 / n_neigh
  B11_inv <- chol2inv(chol(B11))
  
  # Calculating C21, and D12.
  ypop <- apply(ygroupV, 2, mean)
  
  C21 <- array(0, dim = c(length(alpha), num_gamma))
  D12 <- array(0, dim = c(num_gamma, length(alpha)))
  
  for (nn in 1 : n_neigh) {
    C21 <- C21 - t(ygroupV[nn, , drop = FALSE]) %*% t(scores[, nn, drop = FALSE])
    D12 <- D12 + scores[, nn, drop = FALSE] %*% (ygroupV[nn, , drop = FALSE] - ypop)
  }
  C21 <- C21 / n_neigh
  D12 <- D12 / n_neigh
  
  chol_B11_inv <- chol(B11_inv)
  mat1 <- C21 %*% t(chol_B11_inv) %*% t(C21 %*% t(chol_B11_inv))
  mat <- C21 %*% B11_inv %*% D12
  
  var_est_ps <- mat1 + mat + t(mat)
  var_est_ps <- ie_var + var_est_ps / n_neigh
  
  return(var_est_ps) 
}


#' Indirect effect estimates and asymptotic variance.
#' 
#' @param ygroupV An matrix including the group average potential outcome
#' estimates where rows correspond to group, and columns to values of alpha.
#' #' @param boots The results of BootVar() function including estimates of the
#' potential outcomes from the bootstrap samples.
#' @param ps String. Can take values 'true', or 'estimated' for known or
#' estimated propensity score. Defaults to 'true'.
#' @param scores A matrix with rows corresponding to the parameters of the
#' propensity score model and columns for groups. Includes the score of the
#' propensity score evaluated for the variables of each group. Can be left
#' NULL for ps = 'true'.
#' @param alpha_level Numeric. The alpha level of the confidence intervals
#' based on the quantiles of the bootstrap estimates.
#' 
#' @export
#' 
IE(ygroupV = ygroup[, 1 , ], scores = Score.est)
IE(ygroupV = ygroup[, 2 , ], scores = Score.est)
IE <- function(ygroupV, boots = NULL, scores = NULL, alpha_level = 0.05) {
  
  alpha <- as.numeric(dimnames(ygroupV)[[2]])
  ypop <- apply(ygroupV, 2, mean)
  names(ypop) <- alpha
  quants <- c(0, 1) + c(1, - 1) * alpha_level / 2
  norm_quant <- - qnorm(alpha_level / 2)
  
  ie_var <- IEvar(ygroupV = ygroupV, scores = scores)
  dim_names <- c('mu_Q1','mu_Q2', 'SE.est', 'var', 'LB', 'UB')
  if (!is.null(boots)) {
    dim_names <- c(dim_names, 'boot_var', 'boot_var_LB', 'boot_var_UB',
                   'boot_low_quant', 'boot_high_quant')
  }
  
  ie <- array(NA, dim = c(length(dim_names), length(alpha), length(alpha)))
  dimnames(ie) <- list(stat = dim_names, alpha1 = alpha, alpha2 = alpha)
  
  for (a1 in 1 : length(alpha)) {
    for (a2 in 1 : length(alpha)) {
      ie[1, a1, a2] <- ypop[a1]
      ie[2, a1, a2] <- ypop[a2]
      ie[3, a1, a2] <- ypop[a2]/ypop[a1]
      ie[4, a1, a2] <- delta_method(ie_var[c(a1, a2), c(a1, a2)], vec = c(1, -1/(ypop[a1])^2))
      ie_sd <- sqrt(ie[4, a1, a2])
      ie[c(5, 6), a1, a2] <- ie[1, a1, a2] + norm_quant * c(- 1, 1) * ie_sd
    }
  }
  
  if (!is.null(boots)) {
    ie_var_boots <- array(NA, dim = c(length(alpha), length(alpha)))
    for (a1 in 1 : length(alpha)) {
      for (a2 in 1 : length(alpha)) {
        ie[5, a1, a2] <- var(boots[1, a1, ] - boots[1, a2, ])
        ie_sd <- sqrt(ie[5, a1, a2])
        ie[c(6, 7), a1, a2] <- ie[1, a1, a2] + norm_quant * c(- 1, 1) * ie_sd
        ie[c(8, 9), a1, a2] <- quantile(boots[1, a2, ] - boots[1, a1, ],
                                        probs = quants)
      }
    }
  }
  
  return(ie)
}
