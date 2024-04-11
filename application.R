# Application for the semi supervised project
# Programmed by Chen, Jun 2023
# Please read the file Readme.md first before using it
#

# Set working directory
setwd(" Please put in your own directory!")

# Load data, two CSV files required
# The data are available upon request, with the consent of Meteo France.
stat_dat <- read.csv("stat_dat.csv", head=TRUE,sep=";",stringsAsFactors=F, na.strings=c(NA,"NA"," NA"))
fc_dat <- read.csv("fc_dat.csv", head=TRUE,sep=";",stringsAsFactors=F, na.strings=c(NA,"NA"," NA"))

# Length of the forecasting and observed data, number of stations
lfc <- dim(fc_dat)[1]
lstat <- dim(stat_dat)[1]
nstat <- dim(fc_dat)[2]-1 


# Produce the tilde y series

data <- stat_dat[,(1:nstat)+1]

size <- dim(data)
n <- size[1]
m <- size[2]


########################################
## Run test first to select a subset of
## stations satisfying the assumptions
########################################

# Functions for the tests

# Testing for identical distribution over time
testCfun <- function(data,k){
  s <- (1:n)/n
  results <- array(dim=c(m,2))
  
  for(i in 1:m){
    tempdata <- data[,i]
    sdata <- sort(tempdata, decreasing = T)
    threshold <- sdata[k+1]
    Cest <- cumsum((tempdata>threshold)/k)
    results[i,1] <- sqrt(k)*max(abs(Cest-s))
    results[i,2] <- k/n * (t(Cest-s)%*% (Cest-s))
    
  }
  return(results)
}

# Testing for no serial dependence
slidingEstimator <- function(b,n,Y){
  Fhat = ecdf(Y)    # define empirical CDF
  
  # compute the sliding blocks estimator
  slSum = 0;
  for(i in 1:(n+1-b)){
    slSum = slSum + b*(1-Fhat(max(Y[(i):(i+b-1)])))
  }
  theta = (n+1-b)/(slSum)
  return(theta)
  
}

# Function to perform all tests
pretest <- function(){
  set.seed(19810527)
  
  grid <- 5000
  simbm <- 2000
  increments <- rnorm(grid*simbm)
  dim(increments) <- c(grid,simbm)
  simbms <- apply(increments, 2, cumsum)/sqrt(grid)
  
  gstep <- (1:grid)/grid
  
  simbms <- simbms - gstep %*% t(simbms[grid,])
  
  t1c <- apply(abs(simbms),2,max)
  t2c <- apply(simbms^2, 2, sum)/grid
  
  k <- 200
  testing <- testCfun(data,k)
  pvalue <- testing
  
  for(i in 1:m){
    pvalue[i,1] <- sum(testing[i,1]<t1c)/simbm
    pvalue[i,2] <- sum(testing[i,2]<t2c)/simbm
  }
  
  colnames(pvalue) <- c("KStest", "ADtest")
  

  seqb <- seq(from=10, to=100, by=5)
  seqb <- 80
  lb <- length(seqb)

  testTheta <- numeric( nstat*lb )
  dim(testTheta) <- c(lb, nstat)
  Thetasd <- numeric(lb)
  for (j in 1:nstat){
    for (i in 1:lb){
        k<-n/seqb[i]
        testTheta[i,j]<-slidingEstimator(seqb[i],n,data[,j])
        if (j==1) Thetasd[i]<-sqrt(0.2726)/sqrt(k)
    }
  }
  return(list(pvalue,testTheta,Thetasd))
}

# Run all tests to select stations

restest <- pretest()
statselect <- (t(restest[[2]])+qnorm(0.9875)*restest[[3]]<1) | (restest[[1]][,2]<0.025)

statselect[is.na(statselect)] <- TRUE
statselect <- !statselect

########################################
## Determine the basic choice of k
# This is conducted for 3 stations only.
# After making the choice, this part
# should be skipped in applications.
########################################
# selectstat <- c(40,101,83)

# nselstat <- length(selectstat)

# kseq <- seq(from = 30, to =200, by=1)

# for(i in 1:nselstat){
#   z <- selectstat[i]+1
#   pairdat <- data.frame(fc_dat[,z],stat_dat[1:lfc,z])
#   colnames(pairdat) <- c("fc","stat")
#   pairdat <- na.omit(pairdat)
#   npair <- dim(pairdat)[1]

#   dat <- pairdat[,1]
#   sdat <- sort(dat)

#   lengthk <- length(kseq)
#   gammak <- numeric(lengthk)

#   for (j in 1:lengthk){
#     k <- kseq[j]
#     ml <- gp.fit(dat,sdat[npair-k],method="zhang") 
#     gammak[j] <- ml$approx.mean[2] 
#   }

#   MLE_data <- data.frame(gamma = gammak, k = kseq)
#   p <- ggplot(MLE_data, aes(x = k, y = gamma), s) +
#     geom_line() +
#     labs(x = "k",
#        y = expression(gamma)) +
#     theme_minimal()

#   filename = paste0("MLE_station_",selectstat[i])
#   pdf(filename, width=5,height=5)
#   print(p)
#   dev.off()
# }


# selectstat <- (1:nstat)[statselect]

# nselstat <- length(selectstat)

# kseq <- seq(from = 30, to =200, by=1)

# result <- numeric(length(kseq)*nselstat)
# dim(result) <- c(length(kseq),nselstat)

# for(i in 1:nselstat){
#   z <- selectstat[i]+1
#   pairdat <- data.frame(fc_dat[,z],stat_dat[1:lfc,z])
#   colnames(pairdat) <- c("fc","stat")
#   pairdat <- na.omit(pairdat)
#   npair <- dim(pairdat)[1]

#   dat <- pairdat[,1]
#   sdat <- sort(dat)

#   lengthk <- length(kseq)
#   gammak <- numeric(lengthk)

#   for (j in 1:lengthk){
#     k <- kseq[j]
#     ml <- gp.fit(dat,sdat[npair-k],method="zhang") 
#     gammak[j] <- ml$approx.mean[2] 
#   }
#   result[,i] <- gammak
# }



#   MLE_data <- data.frame(gamma = rowMeans(result), k = kseq)
#   p <- ggplot(MLE_data, aes(x = k, y = gamma), s) +
#     geom_line() +
#     labs(x = "k",
#        y = expression(gamma)) +
#     theme_minimal()

#   filename = "MLE_average.pdf"
#   pdf(filename, width=5,height=5)
#   print(p)
#   dev.off()

########################################
## Calculations for the estimators
# This part conducts the estimations for
# all selected stations (passing the
# pretests).
########################################

# Load the functions needed
source("semisupervise.R")
# Specify p and k first
p=1/910
k=136
k1=600

res <- numeric(nstat*13)
dim(res)<-c(nstat,13)
for (s in (1:nstat)){
    if (!statselect[s]) next
    z <- s+1
    pairdat <- data.frame(fc_dat[,z],stat_dat[1:lfc,z])
    colnames(pairdat) <- c("fc","stat")
    pairdat <- na.omit(pairdat)
    npair <- dim(pairdat)[1]

    moredat <- data.frame(rep(NA,lstat-lfc),stat_dat[(lfc+1):lstat,z])
    colnames(moredat) <- c("fc","stat")

    dat <- rbind(pairdat,moredat)

    res_temp <- imp_quantile_mle(dat, p, npair, k, k1)

    res[s,1:4]=array(unlist(res_temp$res_x_org), dim = 4)
    res[s,5:8]=array(unlist(res_temp$res_x_semi), dim = 4)
    res[s,9:12]=array(unlist(res_temp$res_y_org), dim = 4)
    res[s,13]=res_temp$R11
}
var_redu <- 1-res[,6]/res[,2]


#####################################
# Organizing the results and outputs
# This part focuses on the three chosen stations.
# The end product is Table 5 in the paper
#####################################


load("../data/StationDataForChen_long.RData")
res <- data.frame(res,var_redu, stat.loc)
colnames(res) <- c("MLE","var","quantile","quantile_var", "MLE_improved","var_improved","quantile_improved","quantile_improved_var","MLE_y", "var_y", "quantile_y","quantile_y_var", "R11", "var_redu","LON", "LAT")

selectstat <- c(40,101,83)

nselstat <- length(selectstat)
alpha <- 0.90

output <- rep(NA, nselstat*2*15)
dim(output) <- c(15, nselstat*2)
output <- data.frame(output)
filename<-"stations"

for(i in 1:nselstat){
    statres <- res[selectstat[i],]
    # First row, location
    output[1,c(2*i-1,2*i)] <- statres[c("LAT","LON")]
    # Row 2, 4, 6, 7, 8, 10, 12, 14, point estimates
    output[c(2,4,6,7,8,10,12,14),2*i-1] <- as.numeric(statres[c("MLE","MLE_improved", "R11", "var_redu", "quantile", "quantile_improved","MLE_y","quantile_y")])
    # Row 3, 5, 9, 11, 13, 15, confidence intervals for the point estimates in row 2,4,7,10,12,14
    output[3, 2*i-1] <- statres["MLE"] + sqrt(statres["var"]/k) *qnorm((1-alpha)/2)
    output[3, 2*i] <- statres["MLE"] - sqrt(statres["var"]/k) *qnorm((1-alpha)/2)
    output[5, 2*i-1] <- statres["MLE_improved"] + sqrt(statres["var_improved"]/k) *qnorm((1-alpha)/2)
    output[5, 2*i] <- statres["MLE_improved"] - sqrt(statres["var_improved"]/k) *qnorm((1-alpha)/2)
    output[9, 2*i-1] <- statres["quantile"] + sqrt(statres["quantile_var"]) *qnorm((1-alpha)/2)
    output[9, 2*i] <- statres["quantile"] - sqrt(statres["quantile_var"]) *qnorm((1-alpha)/2)
    output[11, 2*i-1] <- statres["quantile_improved"] + sqrt(statres["quantile_improved_var"]) *qnorm((1-alpha)/2)
    output[11, 2*i] <- statres["quantile_improved"] - sqrt(statres["quantile_improved_var"]) *qnorm((1-alpha)/2)
    output[13, 2*i-1] <- statres["MLE_y"] + sqrt(statres["var_y"]/k1) *qnorm((1-alpha)/2)
    output[13, 2*i] <- statres["MLE_y"] - sqrt(statres["var_y"]/k1) *qnorm((1-alpha)/2)
    output[15, 2*i-1] <- statres["quantile_y"] + sqrt(statres["quantile_y_var"]) *qnorm((1-alpha)/2)
    output[15, 2*i] <- statres["quantile_y"] - sqrt(statres["quantile_y_var"]) *qnorm((1-alpha)/2)
    
    filename <- paste(filename, selectstat[i],sep="_")   
}

filename <- paste0(filename, ".csv")
write.csv(output, file=filename, row.names = FALSE)


###############################################
# Making the heatmap for all selected stations
# The end product is Figure 4 in the paper
###############################################



res1 <- res[statselect,c("var_redu","LON","LAT")]
res1 <- as.data.frame(res1)

# Define the map boundaries for France
min_lon <- -5.3
max_lon <- 9.7
min_lat <- 41
max_lat <- 51.5

# Get the map data for France
france_map <- ne_countries(scale = "medium", country = "France", returnclass = "sf")
france_map_mainland <- france_map %>%
  st_crop(xmin = min_lon, xmax = max_lon, ymin = min_lat, ymax = max_lat)

# Convert the res1 data frame to an sf object
res_sf <- st_as_sf(res1, coords = c("LON", "LAT"), crs = 4326)

# Create the heatmap plot
heatmap <- ggplot() +
  geom_sf(data = france_map_mainland, fill = "lightgray", color = "black", size = 0.1, alpha=0.1) +
  geom_sf(data = res_sf, aes(fill = var_redu), size = 4, shape = 21) +
  scale_fill_gradient(low = "yellow", high = "red", limits = c(min(res1$var_redu), max(res1$var_redu)))+
  # scale_fill_viridis_c(low="white", high="red") +
  coord_sf() +
  # geom_sf(data = res_sf, aes(fill = var_redu), size = 4) +
  # scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Display the heatmap plot
pdf("heatmap.pdf", width=5,height=5)
print(heatmap)
dev.off()
