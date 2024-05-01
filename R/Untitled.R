#' TND Partial Interference DATA generation
#' 
#'
library(truncnorm)
datagen_int<-function(rangeN = 400:500, nblocks=1000, OR_C=3,OR_WI=1,OR_WC=5,OR_H=1.5,em=0){
  data.list <- list()
  N = sample(x = rangeN, size = nblocks, replace = T)
  vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = 1)
  for(i in 1:nblocks){
    # Step1. Data generation
    #generate data (WITH clustering for U2)
    C<-runif(n=N[i], 1, 2)
    U1<-rbinom(n=N[i],size=1,prob=0.5) #affects both
    U2<-rbinom(n=N[i],size=1,prob=plogis(0.4+0.2*vaxpos[i])) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    # Step2. Treatment model
    p_trt <- plogis(-0.25 + 0.3*C+0.5*vaxpos[i])
    V<-rbinom(prob= p_trt,size=1,n=N[i])
    g.V = (sum(V)-V)/(N[i]-1)
    # Step3. Outcome model
    #infection
    #Infection (with something) has some common risk factors U1 and C
    Infec<-rbinom(prob=plogis(0.2+0.5*C-5+0.5*U1),size=1,n=N[i]) #current prevalence around 0.007
    #Infected with COVID #vaccine more effective with high vaccination rates
    Infprob=plogis(-3.5 + C -0.5*V -0.2* g.V - 1.2*log(OR_C)*V*(g.V) +em*V*C +log(2)*U2-2*U1)
    #range(Infprob)
    # which(is.na(Infprob))
    Infec_COVID<-rbinom(prob = Infprob, size=1,n=N[i]) #0.018
    #symptoms based on infection
    W=W1=W2=rep(0,N[i])
    W1[Infec==1]<-rbinom(prob=plogis(-2+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
    W2[Infec_COVID==1]<-rbinom(prob=plogis(-1.5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
    W=(W1|W2)+0
    #hospitalization
    H=rep(0,N[i])
    Hprob<-plogis(0.5*C[W==1]-log(OR_H)*V[W==1]-0.5*U1[W==1]+g.V[W==1])
    H[W==1]<-rbinom(prob=Hprob,size=1,n=sum(W==1))
    #selection on outcome for testing (does not condition on infectious status, just being in the hospital)
    R <- as.vector(na.omit(sample((1:N[i])[H == 1], min(rangeN[1], sum(H == 1)))))
    data.list[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V, N=N[i]))[R,]
  }
  return(dplyr::bind_rows(data.list))
}

(datTND<-datagen_int(nblocks=1000))

# What are the observed percentage of treated?
  obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
  hist(obs_alpha$V, breaks = 100)
  print(paste(sum(obs_alpha$V %in% c(0, 1)), 'all treated/control'))
  # ----------- PS estimates ----------- #
  cov_names <- c('C')
  cov_cols <- which(names(datTND) %in% cov_names)
  cov_names <- names(datTND)[cov_cols]
  glm_form <- paste('V ~ (1 | block) +', paste(cov_names, collapse = ' + '))
  phi_hat <- PS_model(datTND, glm_form = glm_form, method = 'TND_IPW')
  # ---------- Coefficients of counterfactual treatment allocation ----------- #
  gamma_numer <- Policy_coef(datTND, gamma_form= glm_form, method = 'TND_IPW')
  # ----------- Calculating the IPW ----------- #
  # Type B: estimand = '1'
  obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
  alpha_range <- quantile(obs_alpha$V, probs = c(0.3, 0.8))
  alpha <- seq(alpha_range[1], alpha_range[2], length.out = 10)
  alpha <- unique(sort(round(c(alpha, 0.1, 0.3), 2)))
  
  resB <-GroupIPW(dta = datTND, cov_cols = cov_cols, phi_hat = phi_hat,
                  alpha = alpha, trt_col = which(names(datTND) == 'V'), out_col = which(names(datTND) == 'Y'),
                  estimand = '1',
                  gamma_numer = gamma_numer)
  ygroup = resB$yhat_group
  
  # ----------- Estimates and asymptotic variance of the population average potential----------- #
  Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = cov_cols,
                         trt_name = 'V', integral_bound = 10)
  Ypop.est <- Ypop(ygroup, ps = 'estimated', scores = Score.est,
                   dta = datTND, use = 'everything')
  de <- DE(ypop = Ypop.est$ypop, ypop_var = Ypop.est$ypop_var, boots = NULL, alpha = alpha,
     alpha_level = 0.05)
  se0 <- IE(ygroupV = ygroup[, 1, ], scores = Score.est)
  se1 <- IE(ygroupV = ygroup[, 2, ], scores = Score.est)
  
r <- 50
# Initialize a list to store results
resDE.org <- vector("list", length = r)
resSE0.org <- vector("list", length = r)
resSE1.org <- vector("list", length = r)

resDE.tnd <- vector("list", length = r)
resSE0.tnd <- vector("list", length = r)
resSE1.tnd <- vector("list", length = r)

# Run the function r times and store the results
for (i in 1:r) {
  tryCatch({
    datTND <- datagen_int(nblocks=1000)
    # ----------- PS estimates ----------- #
    cov_names <- c('C')
    cov_cols <- which(names(datTND) %in% cov_names)
    cov_names <- names(datTND)[cov_cols]
    glm_form <- paste('V ~ (1 | block) +', paste(cov_names, collapse = ' + '))
    phi_hat <- PS_model(datTND, glm_form = glm_form, method = 'TND_IPW')
    
    
    # ---------- Coefficients of counterfactual treatment allocation ----------- #
    gamma_numer <- Policy_coef(datTND, gamma_form= glm_form, method = 'TND_IPW')
    
    # ----------- Calculating the IPW ----------- #
    # Type B: estimand = '1'
    #obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
    #alpha_range <- quantile(obs_alpha$V, probs = c(0.3, 0.8))
    #alpha <- seq(alpha_range[1], alpha_range[2], length.out = 10)
    #alpha <- unique(sort(round(c(alpha, 0.1, 0.3), 2)))
    alpha <- as.numeric(dimnames(de.org)$alpha)
    
    resB <-GroupIPW(dta = datTND, cov_cols = cov_cols, phi_hat = phi_hat,
                    alpha = alpha, trt_col = which(names(datTND) == 'V'), out_col = which(names(datTND) == 'Y'),
                    estimand = '1',
                    gamma_numer = gamma_numer)
    ygroup = resB$yhat_group
    
    # ----------- Estimates and asymptotic variance of the population average potential----------- #
    Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = cov_cols,
                           trt_name = 'V', integral_bound = 10)
    Ypop.est <- Ypop(ygroup, ps = 'estimated', scores = Score.est,
                     dta = datTND, use = 'everything')
    
    resDE.tnd[[i]] <- DE(ypop = Ypop.est$ypop, ypop_var = Ypop.est$ypop_var, boots = NULL, alpha = alpha,
                         alpha_level = 0.05)
    
    #. indices corresponding to treatment = 0,  and [ygroup[, 2, ] is for treatment = 1]
    resSE0.tnd[[i]] <- IE(ygroupV = ygroup[, 1, ], scores = Score.est)
    resSE1.tnd[[i]] <- IE(ygroupV = ygroup[, 2, ], scores = Score.est)
  }, error = function(e) {
    # Print the error message
    cat("Error in iteration", i, ":", conditionMessage(e), "\n")
    # Optionally, you can print more details about the error using conditionCall(e) or traceback()
  })
}


# Save the list to a file
saveRDS(resDE.tnd, file = "resDE_tnd.rds")
saveRDS(resSE0.tnd, file = "resSE0_tnd.rds")
saveRDS(resSE1.tnd, file = "resSE1_tnd.rds")

# Save the list to a file
saveRDS(resDE.org, file = "resDE_org.rds")
saveRDS(resSE0.org, file = "resSE0_org.rds")
saveRDS(resSE1.org, file = "resSE1_org.rds")

# Read the list from the file
resDE.org <- readRDS(file = "resDE_org.rds")

# Print the loaded list
print(my_loaded_list)

# Sample list of 3-d arrays
list_of_mat <- resDE.org
list_of_arrays <- resSE0.org


# Function to compute average of 3-d arrays in the list
average_of_mat <- function(list_of_mat) {
  # Combine the mat into a single 3-d array
  combined_array <- array(unlist(list_of_mat), dim = c(dim(list_of_mat[[1]]), length(list_of_mat)))
  # Compute the average over the number of elements of the dimension-wise of the 3-d array
  dimensionwise_average <- apply(combined_array, c(1, 2), median)
  dimnames(dimensionwise_average) <- list(stat = dimnames(list_of_mat[[1]])[[1]],
                                          alpha = dimnames(list_of_mat[[1]])[[2]])
  return(dimensionwise_average)
}


# Function to compute average of 3-d arrays in the list
average_of_arrays <- function(list_of_arrays) {
  # Combine the arrays into a single 4-d array
  combined_array <- array(unlist(list_of_arrays), dim = c(dim(list_of_arrays[[1]]), length(list_of_arrays)))
  # Compute the average over the number of elements of the dimension-wise of the 3-d array
  dimensionwise_average <- apply(combined_array, c(1, 2, 3), median)
  dimnames(dimensionwise_average) <- list(stat = dimnames(list_of_arrays[[1]])[[1]],
                                          alpha1 = dimnames(list_of_arrays[[1]])[[2]], alpha2 = dimnames(list_of_arrays[[1]])[[3]])
  return(dimensionwise_average)
}

# Compute the average of the arrays in the list
de.org <- average_of_mat(resDE.org)
de.tnd <- average_of_mat(resDE.tnd)

se0.org <- average_of_arrays(resSE0.org)
se0.tnd <- average_of_arrays(resSE0.tnd)

se1.org <- average_of_arrays(resSE1.org)
se1.tnd <- average_of_arrays(resSE1.tnd)


# Specify plot_boot to be 1, 2, or 3 for the different calculations of CIs.
# 1 corresponds to asymptotic, 2 uses variance of bootstrap samples, 3 uses
# bootstrap quantiles.
plot_boot <- 1
index_low <- ifelse(plot_boot == 1, 5, ifelse(plot_boot == 2, 6, 8))
index_high <- index_low + 1

# Creating data frames
de_plot <- data.frame(alpha = as.numeric(dimnames(de.org)$alpha), de_org = de.org[3, ], de_tnd = de.tnd[3, ],
                      low_org = de.org[index_low, ], low_tnd = de.tnd[index_low, ],
                      high_org =  de.org[index_high, ], high_tnd =  de.tnd[index_high, ])

a1 <- which(alpha == 0.13)
se1_plot <- data.frame(alpha = alpha, se1_org = se1.org[3, a1, ], se1_tnd = se1.tnd[3, a1, ],
                       low_org = se1.org[index_low, a1, ], low_tnd = se1.tnd[index_low, a1, ],
                       high_org = se1.org[index_high, a1, ], high_tnd = se1.tnd[index_high, a1, ])

a1 <- which(alpha == 0.21)
se2_plot <- data.frame(alpha = alpha, se2_org = se1.org[3, a1, ], se2_tnd = se1.tnd[3, a1, ],
                       low_org = se1.org[index_low, a1, ], low_tnd = se1.tnd[index_low, a1, ],
                       high_org = se1.org[index_high, a1, ],high_tnd = se1.tnd[index_high, a1, ])


res_array <- array(NA, dim = c(length(alpha), 6, 3))
dimnames(res_array) <- list(alpha = alpha, quant = c('DE.ORG', 'DE.TND', 'SE1.ORG', 'SE1.TND', 'SE2.ORG','SE2.TND'),
                            stat = c('est',  'LB', 'HB'))

# Assign values to the res_array
# For DE
res_array[, 1:2, 'est'] <- as.matrix(de_plot[, c('de_org', 'de_tnd')])
res_array[, 1:2, 'LB'] <- as.matrix(de_plot[, c('low_org', 'low_tnd')])
res_array[, 1:2, 'HB'] <- as.matrix(de_plot[, c('high_org','high_tnd') ])

# For SE1
res_array[, 3:4, 'est'] <- as.matrix(se1_plot[, c('se1_org','se1_tnd') ])
res_array[, 3:4, 'LB'] <- as.matrix(se1_plot[, c('low_org', 'low_tnd')])
res_array[, 3:4, 'HB'] <- as.matrix(se1_plot[, c('high_org','high_tnd')])

# For SE2
res_array[, 5:6, 'est'] <- as.matrix(se2_plot[, c('se2_org','se2_tnd') ])
res_array[, 5:6, 'LB'] <- as.matrix(se2_plot[, c('low_org', 'low_tnd')])
res_array[, 5:6, 'HB'] <- as.matrix(se2_plot[, c('high_org','high_tnd')])

# Create res_df from res_array
res_df <- plyr::adply(res_array[, , 1], 1:2)

# Combine lower and upper bounds for DE and SE plots
res_df$LB.org <- c(de_plot$low_org, se1_plot$low_org, se2_plot$low_org)
res_df$UB.org <- c(de_plot$high_org, se1_plot$high_org, se2_plot$high_org)
res_df$LB.tnd <- c(de_plot$low_tnd, se1_plot$low_tnd, se2_plot$low_tnd)
res_df$UB.tnd <- c(de_plot$high_tnd, se1_plot$high_tnd, se2_plot$high_tnd)

# Define facet labels
f_names <- list('DE.ORG' = expression(DE(alpha) ~ 'Org'),
                'DE.TND' = expression(DE(alpha) ~ 'Tnd'),
                'SE1.ORG' = expression(SE(0.13,alpha) ~ 'Org'),
                'SE1.TND' = expression(SE(0.13,alpha) ~ 'Tnd'),
                'SE2.ORG' = expression(SE(0.21,alpha) ~ 'Org'),
                'SE2.TND' = expression(SE(0.21,alpha) ~ 'Tnd'))

# Define labeller function
f_labeller <- function(variable, value){
  return(f_names[value])
}

# Plotting
library(ggplot2)

# For DE plot
ggplot(data = res_df[res_df$quant %in% c("DE.ORG"), ], aes(x = alpha, y = V1, group = quant)) +  
  geom_line(aes(color = quant, linetype = "org")) +
  geom_ribbon(aes(ymin = LB.org, ymax = UB.org, group = quant), alpha = 0.3) +
  facet_wrap(~ quant, nrow = 1, labeller = f_labeller) +
  xlab(expression(alpha)) + ylab('') +
  theme(axis.title = element_text(size = 12),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0.1, 0.4, by = 0.1)) +
  geom_hline(yintercept = 0, linetype = 2)

# For SE plot
ggplot(data = res_df[res_df$quant %in% c("SE1.ORG", "SE1.TND", "SE2.ORG", "SE2.TND"), ], aes(x = alpha, y = est, group = quant)) +  
  geom_line(aes(color = quant, linetype = "org")) +
  geom_ribbon(aes(ymin = LB, ymax = UB, group = quant), alpha = 0.3) +
  facet_wrap(~ quant, nrow = 1, labeller = f_labeller) +
  xlab(expression(alpha)) + ylab('') +
  theme(axis.title = element_text(size = 12),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0.1, 0.4, by = 0.1)) +
  geom_hline(yintercept = 0, linetype = 2)
