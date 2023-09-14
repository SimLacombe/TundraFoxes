##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##    JAGS file for the dynamic occupancy model
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model{
  
  ## Detection model ##
  
  for(j in 1:nsite) {
    for(ti in 1:nseason) {
      for(day in 1:nsurvey) {
        y[j, ti, day] ~ dcat( rdm[j, (1:nout) , z[j, ti], bait[j,ti, day] ] )
      }
    }
  }
  
  ## Transition model ##
  
  for(j in 1:nsite) {
    
    # for first season
    z[j, 1] ~ dcat( fsm[j, ( 1:nout )] )
    
    # for following seasons
    for(t in 2:nseason){
      z[j, t] ~ dcat( tpm[j, ( 1:nout ) , z[ j, t-1]] )
    }
  }
  
  ## Probability matrices ##
  
  for( j in 1:nsite ) {
    
    ## fsm : First season occupancy probabilities ##
    fsm[j, 1] <- (1 - psinit[1, j]) * (1 - psinit[2, j]) #--------------------|0
    fsm[j, 2] <- psinit[1, j]       * (1 - psinit[2, j]) #--------------------|RF
    fsm[j, 3] <- (1 - psinit[1, j]) * psinit[2, j]       #--------------------|AF
    fsm[j, 4] <- psinit[1, j]       * psinit[2, j]       #--------------------|RF + AF

    ## tpm : array of transition probabilities ##
    # 0 to ...
    tpm[j, 1, 1] <- (1 - gam[1, j])     * (1 - gam[2, j])     #---------------|0
    tpm[j, 2, 1] <- gam[1, j]           * (1 - gam_one[2, j]) #---------------|RF
    tpm[j, 3, 1] <- (1 - gam[1, j])     * gam[2, j]           #---------------|AF
    tpm[j, 4, 1] <- gam[1,  j]          * gam_one[2, j]       #---------------|RF + AF
    # RF to ...
    tpm[j, 1, 2] <- eps[1, j]           * (1 - gam_one[2, j]) #---------------|0
    tpm[j, 2, 2] <- (1 - eps[1, j])     * (1 - gam_one[2, j]) #---------------|RF
    tpm[j, 3, 2] <- eps_one[1, j]       * gam_one[2, j]       #---------------|AF
    tpm[j, 4, 2] <- (1 - eps_one[1, j]) * gam_one[2, j]       #---------------|RF + AF 
    # AF to ...
    tpm[j, 1, 3] <- (1 - gam_one[1, j]) * eps[2, j]           #---------------|0
    tpm[j, 2, 3] <- gam_one[1, j]       * eps_one[2, j]       #---------------|RF
    tpm[j, 3, 3] <- (1 - gam_one[1, j]) * (1 - eps[2, j])     #---------------|AF
    tpm[j, 4, 3] <- gam_one[1, j]       * (1 - eps_one[2, j]) #---------------|RF + AF
    # AF+RF to ..
    tpm[j, 1, 4] <- eps_one[1, j]       * eps_one[2, j]       #---------------|0
    tpm[j, 3, 4] <- (1 - eps_one[1, j]) * eps_one[2, j]       #---------------|RF
    tpm[j, 2, 4] <- eps_one[1, j]       * (1 - eps_one[2, j]) #---------------|AF
    tpm[j, 4, 4] <- (1 - eps_one[1, j]) * (1 - eps[2, j])     #---------------|RF + AF
    
    ## rdm : array of detection probabilities ##
    
    ## Bait absent
    # TS = 0
    rdm[j, 1, 1, 1] <- 1 #----------------------------------------------------|OS = 0
    rdm[j, 2, 1, 1] <- 0 #----------------------------------------------------|OS = RF
    rdm[j, 3, 1, 1] <- 0 #----------------------------------------------------|OS = AF
    rdm[j, 4, 1, 1] <- 0 #----------------------------------------------------|OS = RF + AF
    # TS = RF
    rdm[j, 1, 2, 1] <- (1 - rho[1, j])#---------------------------------------|OS = 0
    rdm[j, 2, 2, 1] <-      rho[1, j] #---------------------------------------|OS = RF
    rdm[j, 3, 2, 1] <- 0 #----------------------------------------------------|OS = AF
    rdm[j, 4, 2, 1] <- 0 #----------------------------------------------------|OS = RF + AF
    # TS = AF
    rdm[j, 1, 3, 1] <- (1 - rho[2, j]) #--------------------------------------|OS = 0
    rdm[j, 2, 3, 1] <- 0 #----------------------------------------------------|OS = RF
    rdm[j, 3, 3, 1] <-      rho[2, j]  #--------------------------------------|OS = AF
    rdm[j, 4, 3, 1] <- 0 #----------------------------------------------------|OS = RF + AF
    # TS = RF + AF
    rdm[j, 1, 4, 1] <- (1 - rho[1, j]) * (1 - rho[2, j]) #--------------------|OS = 0
    rdm[j, 2, 4, 1] <-      rho[1, j]  * (1 - rho[2, j]) #--------------------|OS = RF
    rdm[j, 3, 4, 1] <- (1 - rho[1, j]) *      rho[2, j]  #--------------------|OS = AF
    rdm[j, 4, 4, 1] <-      rho[1, j]  *      rho[2, j]  #--------------------|OS = RF + AF
    
    ## Bait present
    # TS = 0
    rdm[j, 1, 1, 2] <- 1 #----------------------------------------------------|OS = 0
    rdm[j, 2, 1, 2] <- 0 #----------------------------------------------------|OS = RF
    rdm[j, 3, 1, 2] <- 0 #----------------------------------------------------|OS = AF
    rdm[j, 4, 1, 2] <- 0 #----------------------------------------------------|OS = RF + AF
    # TS = RF
    rdm[j, 1, 2, 2] <- (1 - rho_bait[1, j])#----------------------------------|OS = 0
    rdm[j, 2, 2, 2] <-      rho_bait[1, j] #----------------------------------|OS = RF
    rdm[j, 3, 2, 2] <- 0 #----------------------------------------------------|OS = AF
    rdm[j, 4, 2, 2] <- 0 #----------------------------------------------------|OS = RF + AF
    # TS = AF
    rdm[j, 1, 3, 2] <- (1 - rho_bait[2, j]) #---------------------------------|OS = 0
    rdm[j, 2, 3, 2] <- 0 #----------------------------------------------------|OS = RF
    rdm[j, 3, 3, 2] <-      rho_bait[2, j]  #---------------------------------|OS = AF
    rdm[j, 4, 3, 2] <- 0 #----------------------------------------------------|OS = RF + AF
    # TS = RF + AF
    rdm[j, 1, 4, 2] <- (1 - rho_bait[1, j]) * (1 - rho_bait[2, j]) #----------|OS = 0
    rdm[j, 2, 4, 2] <-      rho_bait[1, j]  * (1 - rho_bait[2, j]) #----------|OS = RF
    rdm[j, 3, 4, 2] <- (1 - rho_bait[1, j]) *      rho_bait[2, j]  #----------|OS = AF
    rdm[j, 4, 4, 2] <-      rho_bait[1, j]  *      rho_bait[2, j]  #----------|OS = RF + AF
    
    ## Linear predictors ##
    
    for( i in 1:nspec ) {
      
      # initial occupancy
      psinit[i ,j]  <- ilogit(yr_psi[i, year_cov[j,1]] + inprod( a[i, ], psi_cov[j, ] ))
      
      # colonization
      gam[i, j] <-     ilogit(yr_gam[i, year_cov[j,1]] + inprod( b[i, ], gam_cov[j, ] ))
      
      # extinction
      eps[i, j] <-     ilogit(yr_eps[i, year_cov[j,1]] + inprod( d[i, ], eps_cov[j, ] ))
      
      # detection probability
      ## without bait
      rho[i,j] <-      ilogit(yr_rho[i, year_cov[j,1]] + inprod(f[i, ], rho_cov[j, ] )) 
      ## with bait
      rho_bait[i,j] <- ilogit(yr_rho[i, year_cov[j,1]] + inprod(f[i, ], rho_cov[j, ] ) + bait_[i]) 
      
      # inx on colonization
      gam_one[i, j] <- ilogit(yr_gam[i, year_cov[j,1]] + inprod( b[i, ], gam_cov[j, ] ) + inprod( g[i, ], pi_cov[j, ] ))
      
      # inx on extinction
      eps_one[i, j] <- ilogit(yr_eps[i, year_cov[j,1]] + inprod( d[i, ], eps_cov[j, ] ) + inprod( h[i, ], tau_cov[j, ] ))
    }
  } # closes for loop for j (sites) all the way up at the top of the model
  
  ## Random effects ##
  
  for(i in 1:nspec){
    for ( yrp in 1:nyear ){
      yr_rho[i, yrp] ~ dnorm(0, 1 / (sigma_rho[i] * sigma_rho[i]))
      yr_psi[i, yrp] ~ dnorm(0, 1 / (sigma_psi[i] * sigma_psi[i]))
      yr_gam[i, yrp] ~ dnorm(0, 1 / (sigma_gam[i] * sigma_gam[i]))
      yr_eps[i, yrp] ~ dnorm(0, 1 / (sigma_eps[i] * sigma_eps[i]))
      # NB : dnorm uses precision instead of sd (prec = 1/sd^2)
    }
  }
  
  ## Priors ##
  
  for(i in 1:nspec){
    
    # Initial Occupancy
    for( psip in 1:ncov_psi ){
      a[i, psip] ~ dlogis(0, 1)
    }
    
    # Colonization
    for( gamp in 1:ncov_gam ){
      b[i, gamp] ~ dlogis(0, 1)
    }  
    
    # Extinction
    for( epsp in 1:ncov_eps ){
      d[i, epsp] ~ dlogis(0, 1)
    }
    
    # Detection
    for( rhop in 1:ncov_rho ){
      f[i, rhop] ~ dlogis(0, 1)
    }
    
    # Inxs on colonization
    for( pip in 1:ncov_pi ){
      g[i, pip] ~ dlogis(0, 1)
    }
    
    # Inxs on extinction
    for( taup in 1:ncov_tau ){
      h[i, taup] ~ dlogis(0, 1)
    }

    # bait on detection
    bait_[i] ~ dlogis(0,1)
    
    # RE variance
    
    sigma_rho[i] ~ dunif(0,100)
    sigma_psi[i] ~ dunif(0,100)
    sigma_gam[i] ~ dunif(0,100)
    sigma_eps[i] ~ dunif(0,100)
    
  }
}