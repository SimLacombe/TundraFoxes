library(tidyverse)
library(foreach)

library(LaplacesDemon)

library(runjags)
library(coda)
library(mcmcplots)

library(ggh4x)
library(ggpubr)

rm(list = ls())

years <- as.character(2006:2021)

palette1 <- c(rgb(30/255,136/255,229/255),rgb(255/255,193/255,7/255))
palette2 <- c(rgb(0/255,77/255,64/255), rgb(216/255,27/255,96/255))

quantiles <- c(0.05,0.25,0.50,0.75 ,0.95)

### LOAD MODEL ###

model.path <- "outputs/"
model.filename <- "M_final.rds"

Mod <- readRDS(paste0(model.path, "/", model.filename))

M.mat_ <- as.matrix(as.mcmc.list(Mod), chains = T)
M.mat <- M.mat_[, -grep("z", colnames(M.mat_))]

Nsample <- nrow(M.mat)

YEARMIN = 2006
YEARMAX = 2021
T=7
Nobs_min = 3 
Nweek_min = 3 
Npic_min = 35 

source("src/format_data.R")
source("src/get_steady_state.R")

psi_covs <- c("int","CLG","TFG", "rodents_fall", "supp_feeding")
gam_covs <- c("int","CLG","TFG", "rodents_fall", "supp_feeding")
eps_covs <- c("int","CLG","TFG", "rodents_fall", "supp_feeding")
pi_covs  <- c("int", "rodents_fall", "supp_feeding")
tau_covs <- c("int", "rodents_fall", "supp_feeding")
rho_covs <- c("int")

param_covs <- list(psi = psi_covs,
                   gam = gam_covs,
                   eps = eps_covs,
                   pi  = pi_covs,
                   tau = tau_covs,
                   rho = rho_covs)

data_list <- list(psi_cov = covs[,psi_covs]%>%as.matrix(),
                  gam_cov = covs[,gam_covs]%>%as.matrix(),
                  eps_cov = covs[,eps_covs]%>%as.matrix(),
                  pi_cov = covs[,pi_covs]%>%as.matrix(),
                  tau_cov = covs[,tau_covs]%>%as.matrix(),
                  rho_cov = covs[,rho_covs]%>%as.matrix(),
                  year_cov = covs[,"year"]%>%as.matrix(),
                  ncov_psi = length(psi_covs), ncov_gam = length(gam_covs),
                  ncov_eps = length(eps_covs), ncov_pi = length(pi_covs),
                  ncov_tau = length(tau_covs), ncov_rho = length(rho_covs),
                  nspec = 2, nseason = T, nsite = M,
                  nsurvey = K, nout = 4, nyear = length(years),
                  y = ob_state,
                  bait=bait)

### 1. Model fit and convergence -------------------------------------------------------

M.theta <- as.mcmc.list(Mod)[,c(1:(ncol(M.mat)-1))]
gelman.diag(M.theta)
denplot(M.theta, parms= c("b[2,1]","g[2,1]"))

source("src/calculate_GOF.r")

Chi2_open.df <- get_gof_opened(mm = M.mat_, data_list = data_list, e = .0001, Nsim = 1000)

Chi2_closed.l <- get_gof_closed(mm = M.mat_, data_list = data_list, e = .0001, Nsim = 1000)

Chi2_closed_all.df <- Chi2_closed.l[[1]]
Chi2_closed_species.df  <- Chi2_closed.l[[2]]
Chi2_closed_sites.df    <- Chi2_closed.l[[3]]

P_open <- nrow(Chi2_open.df[Chi2_open.df$chi2.obs<Chi2_open.df$chi2.sim,])/
  nrow(Chi2_open.df)

P_closed_all <- nrow(Chi2_closed_all.df[Chi2_closed_all.df$chi2.obs<Chi2_closed_all.df$chi2.sim,])/
  nrow(Chi2_closed_all.df)

ggplot(Chi2_open.df)+
  geom_point(aes(x=chi2.obs, y = chi2.sim))+
  geom_abline(intercept = 0,slope=1)+
  xlim(c(min(Chi2_open.df),max(Chi2_open.df)))+
  ylim(c(min(Chi2_open.df),max(Chi2_open.df)))+
  xlab("Observed Chi²")+
  ylab("Expected Chi²")+
  geom_text(x=min(Chi2_open.df)+30, y=max(Chi2_open.df)-30, label=paste0("p = ", P_open))+
  ggtitle("Transition model")+
  theme_bw() -> gof_tr


ggplot(Chi2_closed_all.df)+
  geom_point(aes(x=chi2.obs, y = chi2.sim))+
  geom_abline(intercept = 0,slope=1)+
  xlim(c(min(Chi2_closed_all.df),max(Chi2_closed_all.df)))+
  ylim(c(min(Chi2_closed_all.df),max(Chi2_closed_all.df)))+
  xlab("Observed Chi²")+
  ylab("Expected Chi²")+
  geom_text(x=min(Chi2_closed_all.df)+30, y=max(Chi2_closed_all.df)-30, label=paste0("p = ", P_closed_all))+
  ggtitle("Observation model")+
  theme_bw()-> gof_obs

ggarrange(gof_tr, gof_obs, 
          ncol = 2, common.legend = TRUE, legend = "bottom")

ggplot(Chi2_closed_species.df)+
  geom_boxplot(aes(y=chi2.obs-chi2.sim, x = species))+
  ylab("Chi² residual")+
  theme_bw() -> res_sp

ggplot(Chi2_closed_sites.df)+
  geom_boxplot(aes(y=value, x = factor(loc)))+
  ylab("Chi² residual")+
  xlab("site")+
  ylim(quantile(Chi2_closed_sites.df$value, c(0.05,0.95)))+
  theme_bw() -> res_loc

ggarrange(res_sp, res_loc, 
          ncol = 2, widths = c(0.35,0.65))

### 2. Get Model Estimates (Table 1) -------------------------------------------

param_covs <- list(`Initial occupancy` = rep(c("Intercept","CLG","TFG", "Rodents", "Feeding"),
                                           each = 2),
                   Colonization = rep(c("Intercept","CLG","TFG", "Rodents", "Feeding"),
                                      each = 2),
                   Extinction = rep(c("Intercept","CLG","TFG", "Rodents", "Feeding"),
                                    each = 2),
                   Detection  = rep(c("Intercept"), each = 2),
                   Colonization = rep(c("Competition", "Competition x Rodents",
                                      "Competition x Feeding"), each = 2),
                   Extinction = rep(c("Competition", "Competition x Rodents",
                                    "Competition x Feeding"), each = 2),
                   Detection = rep(c("Bait"), each = 2))


Model_summary <- foreach(i = 1:length(param_covs), .combine = rbind) %do% {
  data.frame(Parameter = names(param_covs[i]),
             Covariate = param_covs[[i]], 
             Species = c("Red Fox", "Arctic Fox"))
}

Model_summary <- apply(M.mat[,2:(nrow(Model_summary)+1)], 2,
                   FUN = function(x){quantile(x, c(0.025,0.5,0.975))}) %>%
  t %>%
  as.data.frame %>%
  set_names("lower", "median", "upper")%>%
  cbind(Model_summary, .)

Model_summary<- Model_summary%>%
  mutate(q0 = apply(M.mat[,2:(nrow(Model_summary)+1)], 2, FUN = function(x){ecdf(x)(0)}),
         CI = ifelse(q0>=0.95|q0<=0.05, "90", 
                     ifelse(q0>=0.85|q0<=0.15, "70", 
                            "")),
         effect = ifelse(q0>=0.85, "-",
                         ifelse(q0<=0.15, "+",
                                "")))

Model_summary$Covariate <- factor(Model_summary$Covariate,
                                  levels = rev(c("Intercept", "CLG", "TFG", "Rodents", "Feeding",
                                            "Competition", "Competition x Rodents", "Competition x Feeding",
                                            "Bait")))
Model_summary$Parameter <- factor(Model_summary$Parameter, 
                                  levels = c("Initial occupancy", "Colonization", "Extinction", "Detection"))
ggplot(Model_summary)+
  geom_point(aes(x = median, y = Covariate, color = Species), position=position_dodge(width = 1)) + 
  geom_errorbarh(aes(y = Covariate, xmin = lower, xmax = upper, color = Species), position=position_dodge(width = 1)) +
  geom_vline(xintercept = 0, color = "red") +
  scale_color_manual(values = palette1)+
  facet_wrap(~Parameter, scale = "free")+
  xlab("")+
  ylab("")+
  theme(panel.background = element_rect(color="darkgrey",fill = "white"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color="darkgrey", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 14, hjust = 0),
        strip.placement = "outside",
        axis.text.y = element_text(face = "bold", size = 11),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold", size = 13))+
  guides(colour = guide_legend(title.position="top", title.hjust = 0))

### 3. Plot average probabilities (Fig 3.) -------------------------------------

probs.mat <- invlogit(M.mat[, c("a[1,1]", "b[1,1]", "d[1,1]", "f[1,1]",
                                "a[2,1]", "b[2,1]", "d[2,1]", "f[2,1]")])

probsx.mat <- invlogit(M.mat[, c("b[1,1]", "d[1,1]", "f[1,1]",
                                 "b[2,1]", "d[2,1]", "f[2,1]")] + 
                       M.mat[, c("g[1,1]", "h[1,1]", "bait_[1]",
                                 "g[2,1]", "h[2,1]", "bait_[2]")])

probs.df <- data.frame(species = rep(c("Red Fox", "Arctic Fox"),  
                                   each = Nsample * 4),
                     param = rep(c("Initial occupancy", "Colonization", "Extinction", "Detection"), each = Nsample, 2),
                     competitor  =rep(c("", "Competitor absent", "Competitor absent", "Carrion absent"), each = Nsample, 2),
                     proba = c(probs.mat))
probsx.df <- data.frame(species = rep(c("Red Fox", "Arctic Fox"),  
                                     each = Nsample * 3),
                       param = rep(c("Colonization", "Extinction", "Detection"), each = Nsample, 2),
                       competitor  =rep(c("Competitor present", "Competitor present", "Carrion present"), each = Nsample, 2),
                       proba = c(probsx.mat))

probs.df <- rbind(probs.df, probsx.df)

probs.df$param <- factor(probs.df$param, levels = c("Detection", "Initial occupancy", "Colonization", "Extinction"))

ggplot(probs.df)+
  geom_violin(aes(x= competitor, y=proba, fill = species), show.legend = TRUE)+
  scale_fill_manual(values = palette1)+
  facet_wrap(param ~ ., nrow = 2, scale = "free_x")+
  ylab("probability")+
  xlab("")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 15, hjust = .5),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))+
  guides(fill = guide_legend(title.position="top"))


### 4. Plot stationary occupancy (Fig 4.) --------------------------------------

cov_seq <- list(CLG = seq(-3, 1.5, length.out = 100),
                TFG = seq(-1.5, 4.2, length.out = 100),
                supp_feeding = seq(-0.5, 7, length.out = 100),
                rodents_fall = seq(-1.6, 1.7, length.out = 100))

covs.m <- cbind(rep(1, 100),0,0,0,0)
colnames(covs.m) <- c("int", "CLG", "TFG", "rodents_fall", "supp_feeding")

casual_covNames <- c("CLG", "TFG", "Supplementary feeding", "Rodent abundance")

stat_occ <- foreach(i = seq_along(cov_seq), .combine = rbind) %do%{
  covs.m[, 2:5] <- 0
  covs.m[, names(cov_seq[i])] <- cov_seq[[i]]
  
  stat_data.i <- get_steady_state(mm = M.mat,
                                   data_list = data_list,
                                   param_covs = param_covs,
                                   covs.m = covs.m,
                                   Nsim = 1000,
                                   ncores = 4) %>%
    apply(c(1,3), function(x){c(x[3]+x[4], x[2]+x[4])}) %>%
    apply(c(1,2), quantile, quantiles)
  
  data.frame(species = as.factor(rep(c("Arctic fox", "Red fox"),
                                   each = 100)),
             covariate = casual_covNames[i],
             x = rep(cov_seq[[i]], 2),
             low = c(stat_data.i[1,1,], stat_data.i[1,2,]),
             inf = c(stat_data.i[2,1,], stat_data.i[2,2,]),
             med = c(stat_data.i[3,1,], stat_data.i[3,2,]),
             sup = c(stat_data.i[4,1,], stat_data.i[4,2,]),
             upp = c(stat_data.i[5,1,], stat_data.i[5,2,]))
}

covs.points <- covs %>%
  dplyr::select(CLG, TFG, supp_feeding, rodents_fall)%>%
  rename(`TFG` = TFG,
         `Supplementary feeding` = supp_feeding,
         `Rodent abundance` = rodents_fall)%>%
  pivot_longer(cols = c("CLG", "TFG", "Supplementary feeding", "Rodent abundance"),
               names_to = "covariate", values_to = "x")

covs.points$covariate <- factor(covs.points$covariate, levels = c("CLG", "TFG", "Supplementary feeding", "Rodent abundance"))
stat_occ$covariate <- factor(stat_occ$covariate, levels = c("CLG", "TFG", "Supplementary feeding", "Rodent abundance"))

ggplot(stat_occ)+
  geom_line(aes(x=x, y=med, color = species), size = 1.5)+
  geom_ribbon(aes(x=x, y=med, ymin = inf, ymax = sup, fill = species), alpha = .2, show.legend = TRUE)+
  geom_ribbon(aes(x=x, y=med, ymin = low, ymax = upp, fill = species, color = species), linetype = "dashed", size = 1, alpha = .05, show.legend = FALSE)+
  geom_point(data = covs.points, aes(x = x, y = 0.5), pch = 18, color = "red", size = 2,)+
  ylab("occupancy probability")+
  scale_color_manual(values = palette1, name = "")+
  scale_fill_manual(values = palette1, name = "")+
  xlab("")+
  theme_bw()+
  facet_wrap(~ covariate, scale = "free_x", strip.position = "bottom")+
  theme(strip.text = element_text(face = "bold", size = 11),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold", size = 13))

### 5. Plot colext (Fig 5.) -----------------------------------------------

cov_seq <- list(`Rodent abundance` = seq(-1.6,1.7, length.out = 100),
                `Supplementary feeding` = seq(-0.5, 7, length.out = 100))

id.l <- list(`Rodent abundance` = list(c("b[2,1]", "b[2,4]"), c("b[1,1]", "b[1,4]"),
                                       c("g[2,1]", "g[2,2]"), c("g[1,1]", "g[1,2]"),
                                       c("d[2,1]", "d[2,4]"), c("d[1,1]", "d[1,4]"),
                                       c("h[2,1]", "h[2,2]"), c("h[1,1]", "h[1,2]")),
            `Supplementary feeding` = list(c("b[2,1]", "b[2,5]"), c("b[1,1]", "b[1,5]"),
                                           c("g[2,1]", "g[2,3]"), c("g[1,1]", "g[1,3]"),
                                           c("d[2,1]", "d[2,5]"), c("d[1,1]", "d[1,5]"),
                                           c("h[2,1]", "h[2,3]"), c("h[1,1]", "h[1,3]")))
            
plots.colext <- foreach(i = seq_along(cov_seq)) %do% {
  
  col.mat <- cbind(cbind(1, cov_seq[[i]]) %*% t(M.mat[,id.l[[i]][[1]]])%>%
                     apply(1, function(x){quantile(invlogit(x), quantiles)}),
                   cbind(1, cov_seq[[i]]) %*% t(M.mat[,id.l[[i]][[2]]])%>%
                     apply(1, function(x){quantile(invlogit(x), quantiles)}))
  
  colx.mat <- cbind(cbind(1, cov_seq[[i]]) %*% t(M.mat[,id.l[[i]][[1]]] + M.mat[,id.l[[i]][[3]]])%>%
                      apply(1, function(x){quantile(invlogit(x), quantiles)}),
                    cbind(1, cov_seq[[i]]) %*% t(M.mat[,id.l[[i]][[2]]] + M.mat[,id.l[[i]][[4]]])%>%
                      apply(1, function(x){quantile(invlogit(x), quantiles)}))
  
  ext.mat <- cbind(cbind(1, cov_seq[[i]]) %*% t(M.mat[,id.l[[i]][[5]]])%>%
                     apply(1, function(x){quantile(invlogit(x), quantiles)}),
                   cbind(1, cov_seq[[i]]) %*% t(M.mat[,id.l[[i]][[6]]])%>%
                     apply(1, function(x){quantile(invlogit(x), quantiles)}))
  
  extx.mat <- cbind(cbind(1, cov_seq[[i]]) %*% t(M.mat[,id.l[[i]][[5]]] + M.mat[,id.l[[i]][[7]]])%>%
                      apply(1, function(x){quantile(invlogit(x), quantiles)}),
                    cbind(1, cov_seq[[i]]) %*% t(M.mat[,id.l[[i]][[6]]] + M.mat[,id.l[[i]][[8]]])%>%
                      apply(1, function(x){quantile(invlogit(x), quantiles)}))
  
  data.frame(Param = rep(c("Colonization", "Extinction"), each = 100*2*2),
             Species = rep(c("Arctic Fox","Red Fox"), each = 100,4),
             Competitor = rep(c("Absent", "Present"), each = 100*2, 2),
             x = rep(cov_seq[[i]],4*2),
             low = c(col.mat[1,], colx.mat[1,], ext.mat[1,], extx.mat[1,]),
             inf = c(col.mat[2,], colx.mat[2,], ext.mat[2,], extx.mat[2,]),
             med = c(col.mat[3,], colx.mat[3,], ext.mat[3,], extx.mat[3,]),
             sup = c(col.mat[4,], colx.mat[4,], ext.mat[4,], extx.mat[4,]),
             upp = c(col.mat[5,], colx.mat[5,], ext.mat[5,], extx.mat[5,])) %>%
  ggplot()+
    geom_line(aes(x=x, y=med, color = Competitor), size = 1.5, show.legend = TRUE)+
    geom_ribbon(aes(x=x,y=med, ymin = inf, ymax = sup, fill = Competitor), alpha = .2, show.legend = TRUE)+
    geom_ribbon(aes(x=x,y=med, ymin = low, ymax = upp, color = Competitor), alpha = 0, linetype = "dashed", show.legend = TRUE)+
    geom_point(data = covs.points %>% filter(covariate == names(id.l[i])), aes(x = x, y = 0.5), pch = 18, color = "red", size = 2,)+
    scale_color_manual(values = palette2)+
    scale_fill_manual(values = palette2)+
    facet_nested_wrap(Param ~ Species, ncol = 4, scale = "free_x", nest_line = element_line(), remove_labels = "none")+
    xlab(names(id.l[i]))+
    ylab("Probability")+
    theme_bw()+
    theme(strip.text = element_text(face = "bold", size = 15, hjust = .5),
          strip.placement = "outside",
          strip.background = element_blank(),
          axis.text.y = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(face = "bold", size = 13),
          legend.title = element_text(face = "bold", size = 13))
}

ggarrange(plotlist = plots.colext, nrow = 2, common.legend = TRUE, legend = "bottom")

### 6. Plot Random Effect estimates --------------------------------------------

RE_sigma.df <- M.mat[, grep("sigma", colnames(M.mat))] %>% 
  as.data.frame %>% 
pivot_longer(cols = colnames(M.mat)[grep("sigma", colnames(M.mat))], values_to = "sigma") %>%
  mutate(param = rep(c("Detection", "Initial occupancy", "Colonization", "Extinction"), each = 2, nrow(M.mat)),
         species = rep(c("Red fox", "Arctic fox"), 4 * nrow(M.mat)))

RE_sigma.df$param <- factor(RE_sigma.df$param, 
                                  levels = c("Initial occupancy", "Colonization", "Extinction", "Detection"))
ggplot(RE_sigma.df) + 
  geom_violin(aes(x=param, y = sigma, fill = species), show.legend = FALSE)+
  ylim(c(0,2))+
  scale_fill_manual(values = palette1)+
  theme(panel.background = element_rect(color="darkgrey",fill = "white"),
        strip.text = element_text(face = "bold", size = 15, hjust = .5),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold", size = 13))

