require(compiler)
#force compilation of all functions
enableJIT(3)

########### HELPER FUNCTIONS ########### 

#' Returns a a distance matrix.
#' 
#' Given a matrix where the rows are points and the columns are coordinates
#' it returns another matrix that computes the distances between the points
getDistMtxByPos <- function(pos){
    rows = nrow(pos)
    distances = matrix( rep(0,rows ** 2), nrow=rows, ncol=rows, byrow = FALSE)
    
    for (i in 1:rows){
        for (j in 1:rows){
            xI <- pos[i, 1];
            yI <- pos[i, 2];
            xJ <- pos[j, 1];
            yJ <- pos[j, 2];
            
            distances[i,j] = sqrt( (xJ - xI) ** 2 + (yJ - yI) ** 2);
        } 
    }
    return (distances);
}


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

#' Returns the total tour (Circular) distance for a given cromossome
#' 
#'Lets go slow to understand this this one: 
# GA has no straight foward way of getting a non continuous permutation
# chrmosome, a work around on this  is using chrmossome sequential values
# as indexes for another no sequential chromossome
# because of this you gona find stuff like chrm_map[chromossome[i]]]
getTourDist <- function(chromossome, dpt, chrm_map ){
    # Lets go slow to understand this this: 
    #
    # GA has no straight foward way of getting a non continuous permutation
    # chrmosome, a work around on this  is using chrmossome sequential values
    # as indexes for another no sequential chromossome
    # because of this you gona find stuff like chrm_map[chromossome[i]]]
    dst <- 0
    lst <- dpt
    for(i in chromossome){
        
        if (i > length(chrm_map)) {
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

getDepotListFromChromossome <- function(chromossome){
    dpts = c()
    for( i in 1:length(chromossome)){
        if(!chromossome[i])
            dpts <- c(dpts , i)
    }
    return(dpts)
}

generateClusteredChromossomeMatrix <- function(){
    base <- c(1: (dat$noOfDepots * dat$noOfCustomers))
    
    
    return(c(tail(base, -dat$noOfCustomers * n), head(base, dat$noOfCustomers * n)))
}


binFitWithTour <-  function(binChromossome){
    #
    dpts = c()
    cost <- 0.0
    openCost <- 0.0
    for( i in 1:dat$noOfDepots){
        if(binChromossome[i]){
            dpts <- c(dpts , i)
            openCost <- openCost +  dat$depotsData[i,4]
        }
    }
    if (length(dpts) == 0){
        return(0)
    }
    closestDepots <- assignDepotsByProximity(dpts)
    chrmList = vector("list", dat$noOfDepots)
    
    for (i in dpts){
        chrm <- which( c(closestDepots)   %in%    c(i) )
        chrmList[[i]] <- chrm
        
        if (length(chrm) < 3){
            if(length(chrm) == 1){
                cost <- cost + 2 * mtx[dat$noOfDepots + chrm[1], i]
            }
            else{
                cost <- cost +  mtx[chrm[1], i]
                cost <- cost +  mtx[dat$noOfDepots + chrm[2], dat$noOfDepots + chrm[1]]
                cost <- cost +  mtx[i, chrm[2]]
            }
            next
        }
        
        # ----- this part of code is a waste in tuzun instances since depot capacity is infinite
        # depositLoad <- rep(0, dat$noOfDepots)
        # for(j in chrm){
        #     depositLoad[i] <- depositLoad[i] + dat$customersData[j,3]
        #     #if a depot is serving beyond its capacity this solution has no value
        #     if(depositLoad[i] > dat$depotsData[i, 3]){
        #         return(0)
        #     }
        # }
        # --------------------------------------------------------------------------------------
        # browser()
        
        #this is a lazy evaluation fucntion, considering one truck for the given route problem.
        ft <-  function(chromossome) 1  / getTourDist(chromossome, i, chrm_map = chrm)
        orderGA <- ga(type="permutation", min = 1, max = length(chrm), fitness = ft, pmutation = 0.0, pcrossover = 0.5, monitor = FALSE, maxiter = 20, run = 5, popSize = 10)
        
        if(orderGA@fitnessValue == 0){
            return(0)
        }
        
        # fitness value must be reinverted because the sum of inverses is diferent from inverse of sum
        cost <- cost  + (1/orderGA@fitnessValue)
    }
    
    if(cost != 0){
        # if((cost + openCost) < bestChromossomeCost){
        #     print(binChromossome)
        #     print(bestChromossomesList)
        #     bestChromossomeCost <<- cost + openCost
        #     bestChromossomesList <<- chrmList
        # }
        cost <- 1 /(cost+openCost)
    }
    
    return(cost)
}

recursiveTruckSolver <- function(cstm, dpId, noTruck = 1, pop = 50, maxExplo = 1, bstFt = 0, bestNoTruck = -1){
    #print( str(noTruck) + " Trucks, " + str(maxExplo) + " explo, " + str(bstFt) + " bstft" )
    size = length(cstm)
    if(size < 3){ 
        result <- vector("list", 1)
        result[[1]] <- cstm
        return(result)
    }
    
    cost <- dat$vehicleCost * noTruck
    
    ft <- function(chromossome){
        
        for (i in 1:noTruck){
            iniPos = 1+(size*(i-1))
            finPos = iniPos + size -1
            cost <- cost + getTourDist(chromossome[iniPos:finPos], dpId, cstm)
        }
        if(cost == 0) return(0)
        return(1/cost)
    }
    
    seqGA <- ga(type = "permutation",fitness = ft, min = 1, max = (size * noTruck),maxiter = 100, run = 40, popSize = 50,  pmutation = 0, monitor = FALSE )
    
    # DEBUG REGION
    # 
    # cat('best no truck: ', bestNoTruck)
    # cat(' | no of trucks: ', noTruck)
    # cat(' | Max Explo ', maxExplo)
    # cat(' | best fitness: ', bstFt)
    # cat(' | current ft: ', seqGA@fitnessValue, '\n')
    
    if(seqGA@fitnessValue >= bstFt  || noTruck == bestNoTruck){ 
        if(maxExplo > 0 ){
            return(recursiveTruckSolver(cstm, dpId, noTruck + 1, pop = pop, maxExplo, bstFt = seqGA@fitnessValue, bestNoTruck = noTruck) )
        }
        else{
            # if(noTruck > 1)
            #     browser()
            result <- vector("list", noTruck)
            for (i in 1:noTruck){
                iniPos = 1+(size*(i-1))
                finPos = iniPos + size
                
                region <- seqGA@solution[1,]
                region <-region[iniPos:finPos]
                
                result[[i]] <- cstm[subset(region, region <= size)]
            }
            return(result)
        }
    }
    else{
        # in case of not finding the best solution bust still in exploration mode, go on without updating the best fitness
        if(maxExplo > 0){
            return(recursiveTruckSolver(cstm, dpId, noTruck + 1, pop = pop, maxExplo - 1, bstFt = bstFt, bestNoTruck = bestNoTruck) )
        }
        
        #
        else{
            return(recursiveTruckSolver(cstm, dpId, noTruck -1, pop = pop, maxExplo, bstFt = bstFt, bestNoTruck = bestNoTruck) )
        }
    }
    
}


# chrSize <- dat$noOfDepots * dat$noOfCustomer
# "fits <- function(chromossome){
#     totalDist <- 0
#     depotsCost <- 0
#     counter <- 0
#     t <- Sys.time()
#     for(i in 1:dat$noOfDepots){
#         currtDist <-0
#         demand <- 0
#         if(counter <= dat$noOfCustomers){
#             for (j in 1:dat$noOfCustomers) {
#                 if(counter >= dat$noOfCustomers){
#                     break
#                 }
#                 
#                 # skips fillers genes
#                 if(chromossome[i * j] > dat$noOfCustomers){
#                     next
#                 }
#                 currtDist <- currtDist + mtx[j, i + dat$noOfCustomers]
#                 demand <- demand + dat$customersData[j,3]
#                 counter <- counter + 1
#             }
#         }
#         
#         if(currtDist > 0){
#             depotsCost <- depotsCost + dat$depotsData[i,4]
#         }
#         
#         totalDist <- totalDist + currtDist
#         
#         # Set the fitness as 0 in case of configuration cant handle demand
#         if(demand > dat$depotsData[i,3]){
#             return(0)
#         }
#         
#     }
#     
#     # Avoid 0 division errors when returning the fitness as the inverse of the cost
#     if ((totalDist + depotsCost) == 0){
#         return(0)
#     }
#     # print(Sys.time() - t)
#     #print()
#     return(1/(totalDist + depotsCost))
# }


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