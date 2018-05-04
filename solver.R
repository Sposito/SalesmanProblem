source("ParseDataSet.R")
source("getDistMtxByPos.R")

library(GA)

require(compiler)
enableJIT(3)

#' loads data from a file inside this data object
dat <- parseData()

#' Precomputing the distance matrix for all depots and customers 
#' could take some memory, but saves a lot of unecessary computations.
mtx = getDistMtxByPos(dat$getAllPos())

########### HELPER FUNCTIONS ########### 

#' Returns the closest depot from a given list.
#' This functions recieves a list with an arbitrary number of  
#' depot ids and returns a list with the closest depot of this list 
#' for each customer.
assignDepotsByProximity <- function(depots){
    # Array initializer
    customerNearestDepot <- c()
    
    for (i in 1:dat$noOfCustomers){
        depotID <- -1
        loggedDist <- -1.0
        for(j in depots){
            if (depotID != -1){ # if not initalized this is skipped
                currentDist <- mtx[j, dat$noOfDepots + i]
                if (currentDist < loggedDist){
                    loggedDist <- currentDist
                    depotID <- j
                }
            }
            else{ # Initializes with the first id and distance
                depotID <- j
                loggedDist <- mtx[j, dat$noOfDepots + i]
            }
        }
        customerNearestDepot <- c(customerNearestDepot, depotID)
    }
    return(customerNearestDepot)
}

generateClusteredChromossomeMatrix <- function(){
    base <- c(1: (dat$noOfDepots * dat$noOfCustomers))
    
    
    return(c(tail(base, -dat$noOfCustomers * n), head(base, dat$noOfCustomers * n)))
}


chrSize <- dat$noOfDepots * dat$noOfCustomer


fits <- function(chromossome){
    totalDist <- 0
    depotsCost <- 0
    counter <- 0
    t <- Sys.time()
    for(i in 1:dat$noOfDepots){
        currtDist <-0
        demand <- 0
        if(counter <= dat$noOfCustomers){
            for (j in 1:dat$noOfCustomers) {
                if(counter >= dat$noOfCustomers){
                    break
                }
                
                # skips fillers genes
                if(chromossome[i * j] > dat$noOfCustomers){
                    next
                }
                currtDist <- currtDist + mtx[j, i + dat$noOfCustomers]
                demand <- demand + dat$customersData[j,3]
                counter <- counter + 1
            }
        }
        
        if(currtDist > 0){
            depotsCost <- depotsCost + dat$depotsData[i,4]
        }
        
        totalDist <- totalDist + currtDist
        
        # Set the fitness as 0 in case of configuration cant handle demand
        if(demand > dat$depotsData[i,3]){
            return(0)
        }
        
    }
    
    # Avoid 0 division errors when returning the fitness as the inverse of the cost
    if ((totalDist + depotsCost) == 0){
        return(0)
    }
    # print(Sys.time() - t)
    #print()
    return(1/(totalDist + depotsCost))
}

depotsData <- dat$depotsData

depotsData <- dat$customersData

getDepotListFromChromossome <- function(chromossome){
    dpts = c()
    for( i in 1:length(chromossome)){
        if(!chromossome[i])
            dpts <- c(dpts , i)
    }
    return(dpts)
}

getTourDist <- function(chromossome, dpt, chrm_map ){
    # Lets go slow to understand this this: 
    #
    # GA has no straight foward way of getting a no continuous permutation
    # chrmosome, a work around on this is is using chrmossome sequential values
    # as indexes for another no sequential chromossome
    # because of this you gona find stuff like chrm_map[chromossome[i]]]
    dst <- 0
    lst <- dpt
    for(i in chromossome){
        
        if (i>dat$noOfCustomers) {
            next
        }
        if(dst == 0){
            dst <- dst + mtx[dpt, dat$noOfDepots + chrm_map[i]   ]
            lst <- i + dat$noOfDepots
            next
        }
        d <-mtx[lst,   dat$noOfDepots + chrm_map[i]]
        dst <- dst + d
        lst <- chrm_map[i] + dat$noOfDepots
    }
    dst <- dst + mtx[lst, dpt]
    return(dst)
}

#tst <- ga(type="permutation",min = 1, max = 3, fitness = getTourDist, dpt = 1, chrm_map = c(2,3,4), monitor = TRUE)

binFit <- function(binChromossome){
    
    dpts = c()
    cost <- 0.0
    for( i in 1:dat$noOfDepots){
        if(binChromossome[i]){
            dpts <- c(dpts , i)
            cost <- cost +  dat$depotsData[i,4]
        }
    }
    if (length(dpts) == 0){
        print("len 0")
        return(0)
    }
  
    closestDepots <- assignDepotsByProximity(dpts)
    depositLoad <- rep(0, dat$noOfDepots)
    
    clDp <- -1
    
    for (i in 1:dat$noOfCustomers){
        clDp <- closestDepots[i]
        cost <- cost + mtx[dat$noOfDepots + i, clDp ]
        
        depositLoad[clDp] <- depositLoad[clDp] + dat$customersData[i,3]
        #if a depot is serving beyond its capacity this solution has no value
        if(depositLoad[clDp] > dat$depotsData[clDp, 3]){
            return(0)
        }
    }
    
    
    if(cost != 0){
        cost <- 1 /cost
    }
    return(cost)
    
}

binFitWithTour <-  cmpfun(function(binChromossome){
    dpts = c()
    cost <- 0.0
    for( i in 1:dat$noOfDepots){
        if(binChromossome[i]){
            dpts <- c(dpts , i)
            cost <- cost +  dat$depotsData[i,4]
        }
    }
    if (length(dpts) == 0){
        return(0)
    }
    closestDepots <- assignDepotsByProximity(dpts)
    
    
    for (i in dpts){
        chrm <- which( c(closestDepots)   %in%    c(i) )
        
        
        # ----- this part of code is a waste in tuzun instances since depot capacity is infnite
        depositLoad <- 0
        for(j in chrm){
            depositLoad <- depositLoad + dat$customersData[j,3]
            #if a depot is serving beyond its capacity this solution has no value
            if(depositLoad > dat$depotsData[i, 3]){
                return(0)
            }
        }
        # --------------------------------------------------------------------------------------
        
        ft <-  function(chromossome) 1  / getTourDist(chromossome, i, chrm_map = chrm)
        orderGA <- ga(type="permutation",min = 1, max = length(chrm), fitness = ft, monitor = FALSE)
        
        if(gen@fitnessValue == 0){
            return(0)
        }
        
        # fitness value must be reinverted because the sum of inverses is diferent from inverse of sum
        cost <- cost + (1/orderGA@fitnessValue)
    }
    
    if(cost != 0){
        cost <- 1 /cost
    }
    return(cost)
})

GA.fit <- ga(type = "binary", fitness = binFit, nBits = dat$noOfDepots, maxiter = 300, run = 60, popSize = 100, pmutation = 0.3, pcrossover = 0.6, elitism = 0.3 )
summary(GA.fit)

### Plotting
sol <- getDepotListFromChromossome(GA.fit@solution)
dpListByCst <- assignDepotsByProximity(sol)

dat$customersData <- cbind(dat$customersData,matrix(dpListByCst,nrow=dat$noOfCustomers))
colnames(dat$customersData)[4] <- "id"

dpIds <- c(GA.fit@solution)
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

cD <- dat$customersData
dD <- dat$depotsData

source("advcVisualizer.R")
plotSolution(dat)


#source("visualize.R")
#plotSolution(dat = dat, chosenDepots = sol)

findRouteInsideCluster <- function(depotId, nTrucks){
    # the cromossome is as big as the number o customer to be visited times each truck
    chrm <- c(1: (dat$depotsData[depotId,6] * nTrucks))
    a[]
    routeFit <- function(chromossome){
        cost <- nTrucks * dat$vehicleCost
        for (i in 1:nTrucks){
            # For each 'gene' (truck) it follows the order computing the costs
            # genes with 'nitrogenated bases' higher than the number of customers are 
            # just fillers and are ignored during the cost computation.
            initialPos <- (i-1) * dat$noOfCustomers + 1
            finalPos <- initialPos + (dat$noOfCustomers-1)
            lastValidPos = depotId
            for(j in initialPos:finalPos){
                if(chromossome[j] <= dat$noOfCustomers){
                    cost <- cost + mtx[chromossome[j] + dat$noOfDepots, lastValidPos]
                    lastValidPos <- chromossome[j] + dat$noOfDepots
                }
            }
            # finishes the route computing the cost of returning from the last customer to the depot
            cost <- cost + mtx[depotId, lastValidPos]
            
        }
        return(cost)
    }
    print(routeFit(c(1:2,203:300)))
    
    # geneticOrder <- ga(type="permutation",)
}