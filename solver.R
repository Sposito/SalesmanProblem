source("ParseDataSet.R")
source("helperFunctions.R")

library(GA)

#' loads data from a file inside this data object
dat <- parseData()

#' Precomputing the distance matrix for all depots and customers 
#' could take some memory, but saves a lot of unecessary computations.
mtx = getDistMtxByPos(dat$getAllPos())

depotsData <- dat$depotsData
customerData <- dat$customersData

binGA <- ga(type = "binary", fitness = binFitWithTour, nBits = dat$noOfDepots, maxiter = 500, run = 20, popSize = 30, pmutation = 0.3, pcrossover = 0.6, elitism = 0.3 )

actvDepots <- vector("list", dat$noOfDepots)
closestDepots <- assignDepotsByProximity(which(binGA@solution == 1))
for (i in 1:dat$noOfDepots){
    if(binGA@solution[i])
        actvDepots[[i]] <- which( c(closestDepots)   %in%    c(i) )
}


routes <- vector("list", dat$noOfDepots)
for (i in 1:dat$noOfDepots){
    # browser()
    if(!is.null(actvDepots[i])){
        print(1)
        routes[[i]] <- recursiveTruckSolver(actvDepots[[i]], i)
    }
    
}

dat$routes <- routes

summary(binGA)


### Plotting

sol <- getDepotListFromChromossome(binGA@solution)
dpListByCst <- assignDepotsByProximity(sol)

dat$customersData <- cbind(dat$customersData,matrix(dpListByCst,nrow=dat$noOfCustomers))
colnames(dat$customersData)[4] <- "id"

dpIds <- c(binGA@solution)
dpIds <- replace(dpIds,dpIds==0,'Closed')
dpIds <- replace(dpIds,dpIds==1,'Open')
colnames(dat$depotsData)[5] <- "id"
dat$depotsData <- cbind(dat$depotsData, matrix(dpIds, nrow=dat$noOfDepots))
servdCustomers <- c()

for (i in c(1:dat$noOfDepots)){
    # This could go really wrong if the ids here would be anything but intagers. BEWARE!
    # something like    sum(abs(dat$customersData[,4] - i) < 1e-6) would handle floats
    servdCustomers <- c(servdCustomers, sum(dat$customersData[,4] == i)) 
}




dat$depotsData <- cbind(dat$depotsData, matrix(servdCustomers, nrow=dat$noOfDepots))
colnames(dat$depotsData)[6] <- "customersServed"


dat$chromossomes <- l <- vector("list", dat$noOfDepots)
for (i in c(1:dat$noOfDepots)){
    chrm <- which( c(dat$customersData[,4])   %in%    c(i) )
    dat$chromossomes[[i]] <- chrm
}

# cD <- dat$customersData
# dD <- dat$depotsData

source("advcVisualizer.R")
plotSolution(dat)
