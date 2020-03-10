#libraries
library('reshape2')
library(xlsx)

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
           #as.numeric(sp_dt$Size=="S")) #G <- cbind(as.numeric(D$Diet=="Carnivore"),as.numeric(D$Diet=="Herbivore"))
#Define the covariates for occupancy
X = cam_dt_long4[,c(5:8)] #X = D[,5:8]
#Define the group covariates
#XG = cbind(X*G[,1],X*G[,2],X*G[,3]) #added in reviewed paper, b/c he had 4 categories
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

data <- list(D = cam_dt_long4$Detections, N = ceiling(cam_dt_long4[,"num.nights"]), Species =
               as.numeric(cam_dt_long4$SPECIES), n = nrow(cam_dt_long4), nspp = unique(as.numeric(cam_dt_long4$SPECIES)),X = X, dX = dX)

#Specify the initial values
inits = function() {list(Z = as.numeric(data$D>0))}

#inits.fn = NULL - supposed to choose random

#Specify the parameters to be monitored
params = c("rho","pbeta","spbeta","sigpbeta","mbeta","sigbeta","sbeta",
           "psi.mean","sigma.occ","p.mean","sigma.p","alpha","Z","P")

#Specify the number of chains (nc), number of iterations (ni), burn-in period (nb), and thinning
#rate (nthin)
nc = 3
ni = 60000
nb = 10000
nthin = 50


#Write the model code to a text file called "AllMammals.txt"

#Run the model and call the results "output"
output <- jags(data = data, inits = inits, parameters.to.save = params, model.file
               ="comm_covariance.txt", n.chains =nc, n.iter =ni, n.burnin =nb, n.thin =nthin)

#See a summary of the parameter estimates
output.sum <- output$BUGSoutput$summary

####################################################################################
###Occupancy pull-outs:
mbeta <- output$BUGSoutput$sims.list$mbeta
#gbeta <- output$BUGSoutput$sims.list$gbeta
sbeta <- output$BUGSoutput$sims.list$sbeta
###Detection pull-outs:
pbeta <- output$BUGSoutput$sims.list$pbeta
spbeta <- output$BUGSoutput$sims.list$spbeta
#TO ESTIMATE SPECIES-SPECIFIC COVARIATE VALUES:
#Begin by defining the species
spec <- sp_dt[,1]
#Define the group levels where 1 = carnivore, 2 =other, 3 = primate, 4 = ungulate
#levels(sp_dt$Size) <- levels(sp_dt$Size)[c(1,2,3,4)]
#gg <- as.numeric(sp_dt$Size)
#Define the occupancy covariates and groups
covs <- colnames(X)
#sizes <- c("L", "M","S","XL")
#Create a data frame where the number of rows is equal to the number of covariates * the #number of species
species <- data.frame(expand.grid(covs,spec), matrix(NA,length(covs)*length(spec),4))
#ATTEMPT explained: this is how you pull out all the posterior draws for each covariate for each species.
ATTEMPT<-data.frame(expand.grid(covs,spec), matrix(NA,length(covs)*length(spec),3000))
colnames(species) <- c("Covariate","Species","Mean","SD","LCI","UCI")
#Re-define gbeta
gbeta <- cbind(gbeta, matrix(0,nrow(gbeta),length(covs)))
#Create a loop that will estimate species-specific values for each of the covariates
for (a in 1:length(covs)){
  for (b in 1:length(spec)){
    sims <- mbeta[,a] + sbeta[,b,a]
   # ATTEMPT[(4*(b-1)+a),3:3002]<-c(sims)
   species[(4*(b-1)+a),3:6] <- c(mean(sims), sd(sims),
                                    quantile(sims,c(0.025,0.975)))
  }
}
#Export the results as a table
write.table(x=ATTEMPT,file="cov_posterior.csv",sep=",")
####################################################################################
#TO ESTIMATE SPECIES-SPECIFIC DETECTION VALUES:
pbeta <- output$BUGSoutput$sims.list$pbeta
spbeta <- output$BUGSoutput$sims.list$spbeta
#Begin by defining the species
spec <- sp_dt[,1]
#Define the group levels where 1 = carnivore, 2 =other, 3 = primate, 4 = ungulate
#levels(sp_dt$Order) <- levels(sp_dt$Order[c(1,2,3,4,5)])
#gg <- as.numeric(sp_dt$Order)
#Define the occupancy covariates and groups
covs <- colnames(dX)
#sizes <- c("carnivore", "ungulate","rodent","other","primate")
#Create a data frame where the number of rows is equal to the number of covariates * the #number of species
species_det <- data.frame(expand.grid(covs,spec), matrix(NA,length(covs)*length(spec),4))
ATTEMPTd<-data.frame(expand.grid(covs,spec), matrix(NA,length(covs)*length(spec),3000))
colnames(species_det) <- c("Covariate","Species","Mean","SD","LCI","UCI")
#Re-define gbeta
#gbeta <- cbind(gbeta, matrix(0,nrow(gbeta),length(covs)))
#Create a loop that will estimate species-specific values for each of the covariates
for (a in 1:length(covs)){
  for (b in 1:length(spec)){
    sims <- pbeta[,a] + spbeta[,b,a]
    ATTEMPTd[(4*(b-1)+a),3:3002]<-c(sims)
    #species_det[(4*(b-1)+a),3:6] <- c(mean(sims), sd(sims),
     #                                 quantile(sims,c(0.025,0.975)))
  }
}
#Export the results as a table
write.table(x=ATTEMPTd,file="det_posterior.csv",sep=",")

############################################################################
#TO ESTIMATE SPECIES RICHNESS FOR EACH SITE:
#Begin by increasing memory to avoid errors
memory.limit(size=10000)

#Define the z matrix
z = output$BUGSoutput$sims.list$Z
alpha = output$BUGSoutput$sims.list$alpha

#Sort the data frame based on species, study site, and diet category
#d <- sort_df(merge(data.frame(ID = 1:nrow(D),D[,1:2]),data.frame(SppCode = spec, Group = Spp$Diet)),"ID")[,c(1,3,4)]
#d <- Data_PES[,c(1:2,7,23)]
d <- sort_df(merge(data.frame(ID = 1:nrow(cam_dt),cam_dt[,1:2]),data.frame(SppCode = spec, Group = sp_dt$Size)),"ID")[,c(1,3,4)]
d<-cam_dt_long4[,c(1:2,9)]
#Create a new data frame
dz <- data.frame(d,t(z))
dalpha<-data.frame(d,t(alpha))

#Load reshape2 library
library(reshape2)

#Melt the data frame for easy casting
m.dz <- melt(dz,id.vars = c("SPECIES","StudySite","Size"))
m.alpha<-melt(dalpha,id.vars=c("SPECIES","StudySite","Diet"))

#Aggregate the data by summing the values in the z matrix for each camera station during each #iteration
z.all_test <- acast(m.dz, StudySite ~ variable, fun.aggregate = sum)

#Use the aggregated values to create probability distributions and estimate mean, sd, and 95% 
#credible interval values for camera-station specific species richness
z.all <- t(apply(z.all,1,function(x) c(mean(x),sd(x),quantile(x,c(0.025,0.975)))));
colnames(z.all) = c("Mean","SD","LCI","UCI")

#Export estimates of species richness as a table
write.table(x=z.all,file="Community_richness.csv",sep=",")

#To get your occupancy probabilities from the saved parameter "Z"
#From species richness estimator below: after you create the m.dz table:
library(dplyr)
mdz_tib<-as_tibble(m.dz)
m.alpha<-as_tibble(m.alpha)
comm_occ<-mdz_tib%>%group_by(SPECIES,StudySite)%>%summarise(mean(value))
comm_occ$sd<-mdz_tib%>%group_by(SPECIES,StudySite)%>%summarise(sd(value))
comm_occ_sd<-mdz_tib%>%group_by(SPECIES)%>%summarise(sd(value))
comm_occ_meana<-m.alpha%>%group_by(SPECIES)%>%summarise(mean(plogis(value)))
comm_occ_meanz<-mdz_tib%>%group_by(SPECIES)%>%summarise(mean(value))
comm_occ_meanq<-mdz_tib%>%group_by(SPECIES)%>%quantile(value,0.025)

write.table(x=comm_occ, file="community_psi_est.csv", sep=",")
#This should be the same value as your output summary mean "z".

###################################################################################
###################################################################################
###################################################################################
#TO ESTIMATE GROUP-LEVEL HYPER-PARAMETERS:
#Define the occupancy covariate effects where mbeta is the community-level hyper-parameter,
#gbeta is the group-level hyper-parameter, and sbeta is the species-specific parameter
mbeta <- output$BUGSoutput$sims.list$mbeta
gbeta <- output$BUGSoutput$sims.list$gbeta
sbeta <- output$BUGSoutput$sims.list$sbeta
#Define the covariates and the groups
covs <- colnames(X)
sizes <- c("S","M","L","XL")
#Create a data frame where the number of rows is equal to the number of covariates * the
#number of groups
group <- data.frame(expand.grid(covs,sizes), matrix(NA,length(covs)*length(sizes),4))
colnames(group) <- c("Factor","Group","Mean","SD","LCI","UCI")
#Create a loop estimating the reference group values
for (a in 1:length(covs)){
  group[a,3:6] <- c(mean(mbeta[,a]),sd(mbeta[,a]),quantile(mbeta[,a],
                                                           c(0.025,0.975)))
}
#Create a second loop estimating the other group values
for (a in 1:length(covs)){
  for (b in 1:(length(sizes)-1)){
    sims <- mbeta[,a] + gbeta[,((b-1)*4+a)]
    group[(4*(b)+a),3:6] <- c(mean(sims),sd(sims),quantile(sims,
                                                           c(0.025,0.975)))
  }
}
#Export the results as a table
write.table(x=group,file="Size_psi2.csv",sep=",")

####################################################################################
#TO ESTIMATE GROUP-LEVEL HYPER-PARAMETERS FOR DETECTION:
#Define the detection covariate effects where mbeta is the community-level hyper-parameter,
#gbeta is the group-level hyper-parameter, and spbeta is the species-specific parameter
pbeta <- output$BUGSoutput$sims.list$pbeta
gbeta <- output$BUGSoutput$sims.list$gbeta
spbeta <- output$BUGSoutput$sims.list$spbeta
#Define the covariates and the groups
covs <- colnames(dX)
sizes <- c("S","M","L","XL")
#Create a data frame where the number of rows is equal to the number of covariates * the
#number of groups
group <- data.frame(expand.grid(covs,sizes), matrix(NA,length(covs)*length(sizes),4))
colnames(group) <- c("Factor","Group","Mean","SD","LCI","UCI")
#Create a loop estimating the reference group values
for (a in 1:length(covs)){
  group[a,3:6] <- c(mean(pbeta[,a]),sd(pbeta[,a]),quantile(pbeta[,a],
                                                           c(0.025,0.975)))
}

#Create a second loop estimating the other group values
for (a in 1:length(covs)){
  for (b in 1:(length(sizes)-1)){
    sims <- pbeta[,a] + gbeta[,((b-1)*4+a)]
    group[(4*(b)+a),3:6] <- c(mean(sims),sd(sims),quantile(sims,
                                                           c(0.025,0.975)))
  }
}
#Export the results as a table
write.csv(group,file="Size_det2.csv")

#To estimate group richness for each site:
#Aggregate the data by summing the group-specific values in the z matrix for each camera
#station during each iteration
z.group <- acast(m.dz,StudySite + Size ~ variable, fun.aggregate = sum)
#Use the aggregated values to create probability distributions representing estimated camera-
#station specific group richness
z.group <- t(apply(z.group,1,function(x) c(mean(x),sd(x),quantile(x,c(0.025,0.975)))));
colnames(z.group) = c("Mean","SD","LCI","UCI")
#Export estimates of group richness as a table
write.table(x=z.group,file="Size_rich_slope.csv",sep=",")
