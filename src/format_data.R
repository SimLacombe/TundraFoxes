##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##    Format CT-data and covariates to run the model
##    Called by fit_models.r
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
################################################################################

## Main data path and filenames ##

data.path <- "data/main_dat"
data.filenames <- paste0(data.path,"/",list.files(data.path))

## Covariates paths ##

# 1. geographical covariates
cov.path <- "data/covs/camera_attributes_2019.txt"
# 2. productivity information
vegprod.path <- "data/covs/cam-coords-vegprod-hum.txt"
# 3. rodent abundance
rodents.path <- "data/covs/storskala_04-21.txt"
# 4. density of feeding stations
food_stations.path <- "data/covs/foodStations_kernels_2018-2021.rds"

## Years used in the analysis ##

allyears <- as.character(2005:2021)
years <- as.character(YEARMIN:YEARMAX)
years_id <- which(allyears %in% years)

## Sites to focus on ##

sites_used <- c("komag", "nyborg", "stjernevann")

################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OBSERVATION DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
################################################################################

### 1. LOAD DATA ---------------------------------------------------------------

dat.df <- foreach(f = data.filenames[years_id], .combine = rbind)%do%{
  readRDS(f)%>%
    ungroup%>%
    filter(site %in% sites_used,
           vis == 1)%>%
    dplyr::select(c("year", "site", "loc", "julian2",
                    "ArcticFox", "RedFox", "bait_corr", "newbait"))
}

dat.df$site.year <- paste0(dat.df$loc, ".", dat.df$year)
dat.df$site.year <- factor(dat.df$site.year, levels = unique(dat.df$site.year))

### 2. GET DAILY OBSERVATIONS --------------------------------------------------

daily_dat <- dat.df%>%
  group_by(site.year, year, loc, julian2)%>%
  summarise(RedFox  = as.numeric(any(RedFox>0)),
            ArcticFox  = as.numeric(any(ArcticFox>0)),
            Bait = as.numeric(any(bait_corr>0)),
            newbait = as.numeric(any(newbait>0)),
            Npic = n())

### 3. FILTER DATA -------------------------------------------------------------

# 1. Remove days with less than Npic_min pictures
daily_dat <- daily_dat%>%
  filter(Npic >= Npic_min)

# 2. Remove weeks with less than Nobs_min observations
daily_dat <- daily_dat%>%
  group_by(site.year, year, loc)%>%
  mutate(week = (julian2 - julian2[1])%/%7+1)%>%
  group_by(site.year, week)%>%
  mutate(Nobs = n())%>%
  filter(Nobs>=Nobs_min)

# 3. Remove weeks above T
daily_dat <- daily_dat%>%
  filter(week <= T)

# 4. Remove time series with less than Nweek_min weeks sampled
daily_dat<- daily_dat %>% 
  group_by(site.year) %>% 
  mutate(Nweek = length(unique(week)))%>%
  filter(Nweek>=Nweek_min)

### 4. GET DATA ARRAY TO RUN THE MODEL -----------------------------------------

## Get the occupancy state ## 
#1 = no animal, 2 = RF only, 3 = AF only, 4 = both

daily_dat$state <- paste0(daily_dat$RedFox,daily_dat$ArcticFox)
categories <- data.frame(state = c("00", "10", "01", "11"),
                         lvls = c(seq(1:4)))
daily_dat <- left_join(daily_dat, categories, by = "state")

## Get array parameters ##

# 1. Number of time series  (K)
M <- length(unique(daily_dat$site.year))

# 2. Number of weeks (T)
T

# 3. Number of days per week (K)
K <- 7

## Harmonize all time series to have a length T*K ##

daily_dat <- daily_dat%>%
  group_by(site.year,week)%>%
  mutate(obs = julian2-julian2[1]+1)

daily_dat_all <- data.frame(site.year = rep(unique(daily_dat$site.year), each = K*T),
                            week = rep(1:T, each = K, M),
                            obs = rep(1:K, M*T))

daily_dat <- merge(daily_dat,daily_dat_all, by = c("site.year", "week", "obs"),all=T)

## To a M*T*K array ## 

ob_state <- array(daily_dat$lvls, dim=c(K,T,M))
ob_state <- aperm(ob_state, c(3,2,1))

### 5. GET NAIVE OCCUPANCY TO INITIALIZE THE MODEL -----------------------------

naive_occ <- daily_dat%>%
  group_by(site.year, week)%>%
  summarise(RedFox = as.numeric(any(RedFox>0, na.rm=T)),
            ArcticFox = as.numeric(any(ArcticFox>0, na.rm=T)))

naive_occ$state <- paste0(naive_occ$RedFox, naive_occ$ArcticFox)
naive_occ <- left_join(naive_occ, categories, by = "state")

init <- matrix(naive_occ$lvls, nrow = M, ncol = T, byrow = T)

################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COVARIATE DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
################################################################################

### 1. BAIT --------------------------------------------------------------------

bait <- array(c(daily_dat$Bait), dim=c(K,T,M))
bait <- aperm(bait, c(3,2,1))
bait <- bait + 1 #1: nobait 2: bait
bait[is.na(bait)] <- 1

### 2. GEOGRAPHICAL COVARIATES -------------------------------------------------

all_site.years <- daily_dat%>%
  filter(!is.na(loc))%>%
  group_by(site.year,loc,year)%>%
  summarize()

## Get dataset ##

covs <- read.table(file = cov.path,
                   header = T, row.names = 1)

## Arrange covariate dataset according to the site.year sequences ##

covs <- covs[all_site.years$loc,]%>%
  mutate(site.year = all_site.years$site.year,
         loc = all_site.years$loc,
         year = factor(all_site.years$year, levels = years))

## Add an intercept column and select the covariates of interest ##

covs <- covs%>%
  mutate(int = 1)%>%
  dplyr::select(c("site.year", "loc", "year", "region",
                  "int", "long","lat", "altitude", "distroad", "distforest", "distcoast"))

rownames(covs) <- 1:nrow(covs)

### 3. GET PRODUCTIVITY INFORMATION --------------------------------------------

# 1. Get dataset 
prod <- read.table(file = vegprod.path,
                   header = T, row.names = 2)

# 2. Match the row names with all_site.years$loc
row.names(prod) <- str_to_lower(row.names(prod)) 

# 3. Add proportion of productive areas to the cov dataset
covs$propprod5 <- prod[all_site.years$loc,]$propprod5
covs$propprod10 <- prod[all_site.years$loc,]$propprod10

### 4. RUN PCA ON GEOGRAPHICAL COVS --------------------------------------------

covs_pca <- prcomp(covs[,c("altitude","distroad", "distforest", "distcoast", "propprod5")], scale. = TRUE)

covs <- cbind(covs, covs_pca$x[,c(1,2)])
##PC1: coast to land gradient (CLG)
##PC2: Tundra to Forest gradient (TFG)
covs <- covs %>% 
  rename(CLG = PC1,
         TFG = PC2)

### 5. GET LEMMING ABUNDANCE ---------------------------------------------------

# 1. Get dataset 
rodents_data <- read.table(file = rodents.path,
                           header = T)
# 2. Filter regions 
rodents_data <- rodents_data %>%
  filter(region %in% c("komagdalen", "stjernevann", "vestre_jakobselv"))

# 3. Add rodent abundance the previous fall and the spring to the cov dataset
covs <- covs%>%
  dplyr::rowwise()%>%
  dplyr::mutate(rodents_fall = sum(apply(as.matrix(rodents_data[rodents_data$season=="fall"&
                                                                  rodents_data$year+1==year,5:9]),2,mean)),
                rodents_spr = sum(apply(as.matrix(rodents_data[rodents_data$season=="spring"&
                                                                 rodents_data$year==year,5:9]),2,mean)))

### 6. GET FEEDING STATION PROXIMITY INDEX -------------------------------------

# 1. Get dataset 
# kernel.l : list of 4 spatial density kernel, each one corresponding to a period
#            from 2018 to 2021
kernel.l <- readRDS(file = food_stations.path)

# 2. Get periods corresponding to the kernel layers in the cov dataset 
year_to_period <- data.frame(year_ = 2006:2021,
                         period = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,2,3,4))

covs <- covs%>%
  rowwise()%>%
  mutate(period = year_to_period[year_to_period$year_ == year,]$period)

# 3. Extract feeding station density and add it to the cov dataset

covs <- covs%>%
  rowwise()%>%
  mutate(supp_feeding = 
           ifelse(is.na(period), 0, 
                  raster::extract(kernel.l[[period]], matrix(c(long,lat),1,2))))

### 7. SOME FINAL ADJUSTMENTS  -------------------------------------------------

# 1. scale
covs[,c("altitude", "distroad", "distforest", "distcoast", "propprod5", "propprod10",  "CLG", "TFG", "rodents_fall",
        "rodents_spr", "supp_feeding")] <- scale(covs[,c("altitude","distroad", "distforest","distcoast",
                                                        "propprod5", "propprod10", "CLG", "TFG",
                                                        "rodents_fall","rodents_spr", "supp_feeding")])

# 2. get an index for each year (from 1 to 17)
covs$year <- as.numeric(covs$year)

# 3. Get a treatment covariate starting in 2018
covs$treatment <- as.numeric(all_site.years$year>=2018)


### 8. Plot PCA if required

# library("factoextra")
# fviz_pca_biplot(covs_pca, fill.ind = covs$region, addEllipses = TRUE, 
#                 pointsize = 2.5, pointshape = 21, geom.ind = "point",
#                 col.var = "red",labelsize = 5, col.ind = "black",
#                 repel = TRUE, title = "",legend.title = list(fill = "Region", color = "Region")) +
#   theme(axis.title = element_text(face = "bold"))
