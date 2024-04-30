  # What are the observed percentage of treated?
  obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
  hist(obs_alpha$V, breaks = 100)
  print(paste(sum(obs_alpha$V %in% c(0, 1)), 'all treated/control'))

glm_form <- paste('V ~ (1 | block) +', paste(cov_names, collapse = ' + '))
PS_model <- function(datTND, glm_form, method = c('ORG_IPW', 'TND_IPW')){
  if (method == 'TND_IPW') {
    glmer_fit_col <- lme4::glmer(
      data = datTND[datTND$Y==0, ],
      formula = as.formula(glm_form),
      family = stats::binomial,
      nAGQ = 2
    )
    phi_hat <- list(coefs = summary(glmer_fit_col)$coef[, 1],
                    re_var = as.numeric(summary(glmer_fit_col)$varcor))
  }
  
  if (method == 'ORG_IPW') {
    glmer_fit <- lme4::glmer(
      data = datTND,
      formula = as.formula(glm_form),
      family = stats::binomial,
      nAGQ = 2
    )
    phi_hat <- list(coefs = summary(glmer_fit)$coef[, 1],
                    re_var = as.numeric(summary(glmer_fit)$varcor))
  }
  return(phi_hat)
}


Policy_coef <- function(datTND, gamma_form, method = c('ORG_IPW', 'TND_IPW')){
  if (method == 'TND_IPW') {
    gamma_glmod <- lme4::glmer(as.formula(gamma_form), data = datTND[datTND$Y==0, ], family = 'binomial',
                               control = glmerControl(optimizer = "bobyqa",
                                                      optCtrl = list(maxfun = 2e5)))
    gamma_numer <- summary(gamma_glmod)$coef[, 1]
    rm(list = c('gamma_form', 'gamma_glmod'))
  }
  
  if (method == 'ORG_IPW') {
    gamma_glmod <- lme4::glmer(as.formula(gamma_form), data = datTND, family = 'binomial',
                         control = glmerControl(optimizer = "bobyqa",
                                                optCtrl = list(maxfun = 2e5)))
    gamma_numer <- summary(gamma_glmod)$coef[, 1]
    rm(list = c('gamma_form', 'gamma_glmod'))
  }
  return(gamma_numer)
}

r <- 2
# Initialize a list to store results
resDE.org <- vector("list", length = r)
resSE0.org <- vector("list", length = r)
resSE1.org <- vector("list", length = r)

# Run the function r times and store the results
for (i in 1:r) {
  datTND <- datagen_int(nblocks=1000)
  # ----------- PS estimates ----------- #
  cov_names <- c('C')
  cov_cols <- which(names(datTND) %in% cov_names)
  cov_names <- names(datTND)[cov_cols]
  glm_form <- paste('V ~ (1 | block) +', paste(cov_names, collapse = ' + '))
  phi_hat <- PS_model(datTND, glm_form = glm_form, method = 'ORG_IPW')
  

  # ---------- Coefficients of counterfactual treatment allocation ----------- #
  gamma_numer <- Policy_coef(datTND, gamma_form= glm_form, method = 'ORG_IPW')
  
  # ----------- Calculating the IPW ----------- #
  # Type B: estimand = '1'
  obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
  alpha_range <- quantile(obs_alpha$V, probs = c(0.3, 0.8))
  alpha <- seq(alpha_range[1], alpha_range[2], length.out = 10)
  alpha <- unique(sort(round(c(alpha, 0.1, 0.3), 2)))
  alpha <- unique(sort(round(c(0.1, 0.3), 2)))
  
  resB <-GroupIPW(dta = datTND, cov_cols = cov_cols, phi_hat = phi_hat,
           alpha = alpha, trt_col = which(names(datTND) == 'V'), out_col = which(names(datTND) == 'Y'),
           estimand = '1',
           gamma_numer = gamma_numer)
  ygroup = resB$yhat_group

  
  Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = cov_cols,
                         trt_name = 'V', integral_bound = 10)
  Ypop.est <- Ypop(ygroup, ps = 'estimated', scores = Score.est,
                   dta = datTND, use = 'everything')
  resDE.org[[i]] <- DE(ypop = Ypop.est$ypop, ypop_var = Ypop.est$ypop_var, boots = NULL, alpha = alpha,
                       alpha_level = 0.05)

  #. indices corresponding to treatment = 0,  and [ygroup[, 2, ] is for treatment = 1]
  resSE0.org[[i]] <- IE(ygroupV = ygroup[, 1, ], scores = Score.est)
  resSE1.org[[i]] <- IE(ygroupV = ygroup[, 2, ], scores = Score.est)
}

# Save the list to a file
saveRDS(resDE.org, file = "resDE_org.rds")
saveRDS(resSE0.org, file = "resSE0_org.rds")
saveRDS(resSE1.org, file = "resSE0_org.rds")
# Read the list from the file
resDE.org <- readRDS(file = "resDE_org.rds")

# Print the loaded list
print(my_loaded_list)


# Calculate the average among the first row of the output for the r replicates
apply(sapply(resDE.org, function(x) x[1, ]), 1, mean)
apply(sapply(resDE.org, function(x) x[2, ]), 1, mean)
apply(sapply(resDE.org, function(x) x[3, ]), 1, mean)




# Sample list of 3-d arrays
list_of_arrays <- resSE0.org

# Function to compute average of 3-d arrays in the list
average_of_arrays <- function(list_of_arrays) {
  # Combine the arrays into a single 4-d array
  combined_array <- array(unlist(list_of_arrays), dim = c(dim(list_of_arrays[[1]]), length(list_of_arrays)))
  # Compute the average over the number of elements of the dimension-wise of the 3-d array
  dimensionwise_average <- apply(combined_array, c(1, 2, 3), mean)
  dimnames(dimensionwise_average) <- list(stat = dimnames(list_of_arrays[[1]])[[1]],
                                          alpha1 = dimnames(list_of_arrays[[1]])[[2]], alpha2 = dimnames(list_of_arrays[[1]])[[3]])
  return(dimensionwise_average)
}

# Compute the average of the arrays in the list
average_of_arrays(resSE0.org)



# create two vectors and take them as input in array 
vector1 <- c(3, 2, 1) 
vector2 <- c(2, 4, 6, 8, 0, 1) 
new.array <- array(c(vector1, vector2), dim = c(3, 3, 2)) 
print(new.array) 

new.array[1,,]
# using apply and calculate the sum of rows in matrices 
result <- apply(new.array, c(1,2), mean) 
print(result) 
