#libraries
library('reshape2')
library(xlsx)
library(tidyverse)
library(ggridges)

#for melt to work in the beginning: detach("package:reshape", unload=TRUE)

#Load data
cam_loc <- read.xlsx("Data2.xlsx", sheetName = "CamLoc")
cam_dt <- read.xlsx("Data2.xlsx", sheetName = "daily_occ")
psi_CoVar <- read.xlsx("Data2.xlsx", sheetName = "PsiCov")
det_CoVar <- read.xlsx("Data2.xlsx", sheetName = "DetCov")
sp_dt <- read.xlsx("Data2.xlsx", sheetName = "Species")
colnames(sp_dt)[1] <- "SPECIES"

#scaling code for any variable that needs it
psi_CoVar<- cbind(psi_CoVar[1:2], cbind(scale(psi_CoVar[,3:6], center = TRUE, scale = TRUE)))
det_CoVar<-cbind(det_CoVar[1:3], cbind(scale(det_CoVar[,4:5], center = TRUE, scale = TRUE)))

### process the data in long format
cam_dt_long <- melt(cam_dt,id.vars = "StudySite",variable.name="SPECIES",value.name= "Detections")
#merge with detetection covariates
cam_dt_long2 <- merge(cam_dt_long,det_CoVar,by="StudySite",all=TRUE)
#merge with psi covariates
cam_dt_long3 <- merge(cam_dt_long,psi_CoVar,by="StudySite",all=TRUE)
#add species specific feeding guild
cam_dt_long4 <- merge(cam_dt_long3,sp_dt[,c("SPECIES","Size")],by="SPECIES",all=TRUE)

#Define the group variables
#G <- cbind(as.numeric(sp_dt$Size=="L"),as.numeric(sp_dt$Size=="M"),
           #as.numeric(sp_dt$Size=="S")) 
G <- cbind(as.numeric(sp_dt$Diet=="Carnivore"),as.numeric(sp_dt$Diet=="Herbivore"))
#Define the covariates for occupancy
X = cam_dt_long4[,c(5:8)] #X = D[,5:8]
#Define the group covariates
XG = cbind(X*G[,1],X*G[,2]) 
#Define the covariates for detection
dX = cam_dt_long2[,c(4:7)]#dX = D[,9:11]


#Load the necessary libraries
library(R2jags); library(reshape); library(plyr)
#Define the necessary arguments to run the jags command
#Load all the data including the detection array, number of sampling occasions, individual
#species sampled, total number of sampled species, and covariate information

#data <- list(D = D$Detections, N = ceiling(D[,"num.nights"]), Species =
#               as.numeric(D$SppCode), n = nrow(D), nspp = max(as.numeric(D$SppCode)),X =
#               X, XG = XG, dX = dX)

data <- list(D = cam_dt_long4$Detections, 
             N = ceiling(cam_dt_long4[,"num.nights"]), 
             Species = as.numeric(cam_dt_long4$SPECIES), 
             n = nrow(cam_dt_long4), 
             nspp = max(as.numeric(cam_dt_long4$SPECIES)),
             X = X, 
             XG = XG,
             dX = dX)

#Specify the initial values
inits = function() {list(Z = as.numeric(data$D>0))}

#inits.fn = NULL - supposed to choose random

#Specify the parameters to be monitored
params = c("rho","pbeta","spbeta","sigpbeta","mbeta","sigbeta","sbeta","gbeta","psi.mean","sigma.occ","p.mean","sigma.p","alpha","Z","P")


#Specify the number of chains (nc), number of iterations (ni), burn-in period (nb), and thinning
#rate (nthin)
nc = 3
ni = 60000
nb = 10000
nthin = 50


# Run occupancy model. The model file must be in the correct folder (see script AllMammals_pcov.txt)
out3 <- jags(data = data, 
             inits = inits, 
             parameters.to.save = params, 
             model.file ="AllMammals_pcov.txt", 
             n.chains = nc, 
             n.iter = ni,
             n.burnin = nb, 
             n.thin = nthin)

#Write the model code to a text file called "AllMammals.txt"

#See a summary of the parameter estimates
out3.sum <- output$BUGSoutput$summary

write.table(x=out3.sum,file="results/dietgrp-030920.csv",sep=",", row.names=F) # save output

#Name this output alpha and P.
alpha <- out3$BUGSoutput$sims.list$alpha
p <- out3$BUGSoutput$sims.list$P

# we are log transforming alphas b/c they are on logit scale.
expit <- function(x)  1/(1+exp(-x))
logit.alpha <- expit(alpha)
logit.p <- expit(p)

sppnames <- c("SppCode", as.character(unique(sp_dt$SPECIES)))
sppnames2 <- as.character(unique(sp_dt$SPECIES))

# Columns represent species, so we take the mean of each column to get the mean psi and p value (when covariates are at their mean values) for each species.

psimeans <- colMeans(logit.alpha)
names(psimeans) <- sppnames2
psimeans <- as.data.frame(psimeans)
write.table(x = psimeans, file = "results/dietgrp-030920_alphaspsi()p().csv", sep = ",")

pmeans <- colMeans(logit.p)
names(pmeans) <- sppnames2
write.table(x = pmeans, file = "results/dietgrp-030920_detection.csv", sep=",")

# Get the quantiles and 95% confidence intervals for psi and p.

apply(logit.alpha, 2, function(x) sort(x)[])
psiCI <- apply(logit.alpha, 2, function(x) quantile(x,probs = c(0.025,0.1,0.5,0.9,0.975)))
colnames(psiCI) <- sppnames2
write.table(x=psiCI, file="results/dietgrp-030920_alphaCI.psi()p().csv",sep=",")

apply(logit.p, 2, function(x) sort(x)[])
pCI <- apply(logit.p, 2, function(x) quantile(x,probs = c(0.025,0.1,0.5,0.9,0.975)))
colnames(pCI) <- sppnames2
write.table(x = pCI, file="results/dietgrp-030920_pCI.psi()p().csv", sep = ",")


# Define the occupancy covariate effects where mbeta is the community-level hyperparameter, gbeta is the group-level hyperparameter, and sbeta is the species-specific parameter.
mbeta <- out3$BUGSoutput$sims.list$mbeta
gbeta <- out3$BUGSoutput$sims.list$gbeta
sbeta <- out3$BUGSoutput$sims.list$sbeta

# Calculate group-level estimates.
covs <- colnames(X) # define covariates
sizes <- c("Omnivore", "Carnivore", "Herbivore") # define groups
group <- data.frame(expand.grid(covs, sizes), matrix(NA, length(covs) * length(sizes), 4)) # create data frame where number of rows is equal to the number of covariates * the number of groups
colnames(group) <- c("Factor", "Group", "Mean", "SD", "LCI", "UCI")

# Create a loop estimating the reference group values.
for (a in 1:length(covs)){
  group[a,3:6] <- c(mean(mbeta[,a]),sd(mbeta[,a]),quantile(mbeta[,a],c(0.025,0.975)))
}

# Create a second loop estimating the other group values.
for (a in 1:length(covs)){
  for (b in 1:(length(sizes)-1)){
    sims <- mbeta[,a] + gbeta[,((b-1)*ncol(X)+a)]
    group[(ncol(X)*(b)+a),3:6] <- c(mean(sims),sd(sims),quantile(sims,c(0.025,0.975)))
  }
}

# Export table with group values.
write.table(x = group, file = "results/dietgrp-030920_group.csv", sep = ",", row.names=F)

# Species level estimates 1 = carnivore, 2 = herbivore, 3, omnivore. 

# Define the species
spec <- sp_dt[,1]

# Define the group levels
levels(sp_dt$Diet) <- levels(sp_dt$Diet)[c(1,2,3)]
gg <- as.numeric(sp_dt$Diet)

# Define the occupancy covariates and groups
covs <- colnames(X)
sizes <- c("Omnivore","Carnivore","Herbivore")

# Create a data frame where the number of rows is equal to the number of covariates * the number of species
species <- data.frame(expand.grid(covs,spec), matrix(NA,length(covs)*length(spec),4))
colnames(species) <- c("Factor","Species","Mean","SD","LCI","UCI")

# Re-define gbeta
gbeta <- cbind(gbeta, matrix(0,nrow(gbeta),length(covs)))

# Create a loop that will estimate species-specific values for each of the covariates
for (a in 1:length(covs)){
  for (b in 1:length(spec)){
    sims <- mbeta[,a] + gbeta[,((gg[b]-1)*ncol(X)+a)] + sbeta[,b,a]
    species[(ncol(X)*(b-1)+a),3:6] <- c(mean(sims),sd(sims),quantile(sims,c(0.025,0.975)))
  }
}

# export the table
write.table(x=species,file="results/dietgrp-030920_species.csv",sep=",", row.names=F) 



# calculate distributions for community hyperparameters
comm_param <- as.data.frame(mbeta)
names(comm_param) <- covs
head(comm_param)
comm_param_long <- comm_param %>% 
  pivot_longer(cols = everything(), names_to = "Factor", values_to = "Value")
write.csv(comm_param_long, "results/dietgrp-030920_communityposteriors.csv", row.names = F)

# and for group hyperparameters
head(gbeta)
group_param <- as.data.frame(gbeta)

group_param_carn <- group_param[,1:4]
names(group_param_carn) <- covs
group_param_carn_long <- group_param_carn %>% 
  pivot_longer(cols = everything(), names_to = "Factor", values_to = "Value") %>% 
  mutate(Group = "Carnivore")

group_param_herb <- group_param[,5:8]
names(group_param_herb) <- covs
group_param_herb_long <- group_param_herb %>% 
  pivot_longer(cols = everything(), names_to = "Factor", values_to = "Value") %>% 
  mutate(Group = "Herbivore")

group_param_omni <- comm_param # reference group
names(group_param_omni) <- covs
group_param_omni_long <- group_param_omni %>% 
  pivot_longer(cols = everything(), names_to = "Factor", values_to = "Value") %>% 
  mutate(Group = "Omnivore")

group_param_long <- rbind(group_param_herb_long, group_param_carn_long, group_param_omni_long)
head(group_param_long)

ggplot(group_param_long, aes(x = Value, y = Factor, fill = Group)) +
  geom_density_ridges(scale=.8) +
  theme_ridges() +
  theme(legend.position = "none") +
  geom_vline(xintercept=0) +
  facet_wrap(~Group)

ggplot(group_param_long, aes(x = Value, group = Group, fill = Group)) +
  geom_density() +
  theme_ridges() +
  geom_vline(xintercept=0) +
  facet_grid(Factor~Group)




# Species richness for each site - FIX FOR TARA DATA

# Define the z matrix
z = out3$BUGSoutput$sims.list$Z


# Sort the data frame based on species, study site, and diet category
d <- sort_df(merge(data.frame(ID = 1:nrow(cam_dt_long),cam_dt_long[,1:2]),data.frame(SppCode = spec, Group = sp_dt$Diet)),"ID")[,c(1,3,4)]

# Create a new data frame
dz <- data.frame(d,t(z))

# Melt the data frame for easy casting
m.dz <- melt(dz,id.vars = c("SppCode","StudySite","Group") )

# Aggregate the data by summing the values in the z matrix for each StudySite station during each iteration
z.all <- acast(m.dz,StudySite ~ variable, fun.aggregate = sum)

# Use the aggregated values to create probability distributions and estimate mean, sd, and 95% credible interval values for StudySite-station specific species richness
z.all <- t(apply(z.all,1,function(x) c(mean(x),sd(x),quantile(x,c(0.025,0.975)))))
names <- rownames(z.all)
rownames(z.all) <- NULL
z.all <- cbind(names,z.all)
colnames(z.all) = c("StudySite", "Mean","SD","LCI","UCI")

# Export estimates of species richness as a table
write.table(x=z.all,file="results/dietgrp-030920_spprich.csv",sep=",", row.names=F)

# To estimate group richness for each site:

# Aggregate the data by summing the group-specific values in teh z matrix for each StudySite station during each iteration
z.group <- acast(m.dz,StudySite + Group ~ variable, fun.aggregate = sum)

# Use the aggregated values to create probability distributions representing estimated StudySite-station specific group richness
z.group <- t(apply(z.group,1,function(x) c(mean(x),sd(x),quantile(x,c(0.025,0.975)))))
names <- rownames(z.group)
rownames(z.group) <- NULL
z.group <- cbind(names,z.group)
colnames(z.group) = c("StudySite", "Mean","SD","LCI","UCI")

# Export estimates of group richness as a table
write.table(x=z.group,file="results/dietgrp-030920_grouprich.csv",sep=",", row.names=F)