##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##    Functions to calculate stationary occupancy probabilities from 
##    the model estimates
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#install.packages("expm")
library(doParallel)

get_steady_state <- function(mm = NULL, data_list = NULL,
                             param_covs, covs.m, Nsim = 1000, ncores = 4){
    
  #locate each element in data_list to a single variable
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  nsite <- nrow(covs.m)
   
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  n <- nrow(mm)

  O <- sample(1:nrow(mm),Nsim)
  
  steady.array <- foreach (i = 1:Nsim,
                           .packages = c('LaplacesDemon', 'expm'))%dopar%{
    o <- O[i]
    
    steady.m <- array(0, dim = c(nsite, 4))
    
    b <- matrix(mm[o,grep("b\\[", colnames(mm))], ncol = ncov_gam, nrow = nspec)
    d <- matrix(mm[o,grep("d\\[", colnames(mm))], ncol = ncov_eps, nrow = nspec)
    g <- matrix(mm[o,grep("g\\[", colnames(mm))], ncol = ncov_pi, nrow = nspec)
    h <- matrix(mm[o,grep("h\\[", colnames(mm))], ncol = ncov_tau, nrow = nspec)
    
    gam     <- invlogit(b %*% t(matrix(covs.m[,param_covs[["gam"]]], nrow = nsite)))
    eps     <- invlogit(d %*% t(matrix(covs.m[,param_covs[["eps"]]], nrow = nsite)))
    
    gam_one <- invlogit(b %*% t(matrix(covs.m[,param_covs[["gam"]]], nrow = nsite)) +
                          g %*% t(matrix(covs.m[,param_covs[["pi"]]], nrow = nsite)))
    eps_one <- invlogit(d %*% t(matrix(covs.m[,param_covs[["eps"]]], nrow = nsite)) +
                          h %*% t(matrix(covs.m[,param_covs[["tau"]]], nrow = nsite)))
    
    # Transition probabilities for each site
    tpm <- array(0, dim = c(nsite, 4, 4))
    # U to ...
    tpm[, 1, 1] <- (1 - gam[1, ])     * (1 - gam[2, ])     #-----------------|U
    tpm[, 2, 1] <- gam[1, ]           * (1 - gam_one[2, ]) #-----------------|A
    tpm[, 3, 1] <- (1 - gam[1, ])     * gam[2, ]           #-----------------|B
    tpm[, 4, 1] <- gam[1,  ]          * gam_one[2, ]       #-----------------|AB
    # A to ...
    tpm[, 1, 2] <- eps[1, ]           * (1 - gam_one[2, ]) #-----------------|U
    tpm[, 2, 2] <- (1 - eps[1, ])     * (1 - gam_one[2, ]) #-----------------|A
    tpm[, 3, 2] <- eps_one[1, ]       * gam_one[2, ]       #-----------------|B
    tpm[, 4, 2] <- (1 - eps_one[1, ]) * gam_one[2, ]       #-----------------|AB 
    # B to ...
    tpm[, 1, 3] <- (1 - gam_one[1, ]) * eps[2, ]           #-----------------|U
    tpm[, 2, 3] <- gam_one[1, ]       * eps_one[2, ]       #-----------------|A
    tpm[, 3, 3] <- (1 - gam_one[1, ]) * (1 - eps[2, ])     #-----------------|B
    tpm[, 4, 3] <- gam_one[1, ]       * (1 - eps_one[2, ]) #-----------------|AB
    # AB to ..
    tpm[, 1, 4] <- eps_one[1, ]       * eps_one[2, ]       #-----------------|U
    tpm[, 3, 4] <- (1 - eps_one[1, ]) * eps_one[2, ]       #-----------------|A
    tpm[, 2, 4] <- eps_one[1, ]       * (1 - eps_one[2, ]) #-----------------|B
    tpm[, 4, 4] <- (1 - eps_one[1, ]) * (1 - eps[2, ])     #-----------------|AB

    for(i in 1:nsite){
      tpm[i,,] <- apply(tpm[i,,], 2, function(x) x/ sum(x))
      steady.m[i,]  <-   (tpm[i,,]%^%20)%*%c(0.25,0.25,0.25,0.25)
    }
    steady.m
  }%>%sapply(.,I,simplify="array")
  stopCluster(cl)
  return(steady.array)
}
