#' Parses .dat files in a custom data object
#' This function accepts a string to the path to the
#' desired file, if none is provied it toogles a file
#' chooser window.
parseData <- function(path = ""){
    if(path == ""){
        filePath <- file.choose()
    }
    
    file = file(filePath, "r")
    print (filePath)
    
    #    Data Variables
    noOfDepots      <- -1
    noOfCustomers   <- -1
    vehicleCapacity <- -1
    vehicleCost     <- -1
    
    valueIsDecimal  <-  FALSE
    
    #      Data Matrizes
    depotsData      <- matrix()
    customersData   <- matrix()
    
    #  Parse Control Variables
    depotsPosDone       <- FALSE
    customersPosDone    <- FALSE
    depotCapDone        <- FALSE
    customerDemandsDone <- FALSE
    depotCostDone       <- FALSE
    
    # Number of lines to be read each step
    linesToRead     <-  1
    posSep <- "\t"
    #posSep = ";"
    # Reads a given amount of lines until the end of the file
    while(TRUE){
        line = readLines(file, n = linesToRead)
        
        if ( length(line) == 0 ) { 
            break; # breaks the while loop when there is no more lines to be read.
        }
        
        # Skips blank lines, when skipping the program can ajust the amount of lines
        # needed to read a lines block
        if(line[1] == ""){
            if(!depotsPosDone && noOfDepots > 0){
                linesToRead <- noOfDepots
                next
            }
            
            if(length(line) > 1 && !customersPosDone && line[1] == ""){
                line <- line[2:noOfCustomers+1]
                next
            }
            
            if(ncol(depotsData) == 2 && vehicleCapacity > 0 && !depotCapDone){
                linesToRead <- noOfDepots
                next
            }
            
            
            if(ncol(customersData) == 2 && !customerDemandsDone && depotCapDone){
                linesToRead <- noOfCustomers;
                next
            }
            
            if(ncol(depotsData) == 3 &&  !depotCostDone && customerDemandsDone){
                linesToRead <- noOfDepots;
                next
            }
            next # default case
        }
        
        # ---- Section where progressive lines or chucks are parsed and stored in ----
        # ---- apropriate variables. -------------------------------------------------
        
        if(noOfCustomers < 0){
            noOfCustomers <- strtoi(line)
            next
        }
        
        if(noOfDepots < 0 ){
            noOfDepots <- strtoi(line)
            next
        }
        
        if(!depotsPosDone){
            depotsData <- as.matrix(read.table(textConnection(line), posSep, header = FALSE))
            # browser()
            depotsData <- depotsData[1:noOfDepots,1:2]
            colnames(depotsData)[1] <- "X"
            colnames(depotsData)[2] <- "Y"
            
            linesToRead <- noOfCustomers
            depotsPosDone <- TRUE
            next
        }
        
        if(!customersPosDone){
            customersData <- as.matrix(read.table(textConnection(line), posSep, header = FALSE))
            customersData <- customersData[1:noOfCustomers,1:2]
            colnames(customersData)[1] <- "X"
            colnames(customersData)[2] <- "Y"
            
            linesToRead <- 1
            customersPosDone <- TRUE;
            next
        }
        
        if(vehicleCapacity < 0){
            vehicleCapacity <- strtoi(line)
            next
        }
        
        if(ncol(depotsData) == 2 && !depotCapDone && length(line) == noOfDepots ){
            depotCapDone = TRUE
            depotsData <- cbind(depotsData, as.matrix(read.table(textConnection(line), "\n", header = FALSE)));
            colnames(depotsData)[3] <- "Capacity"
            
            linesToRead <- 1
            next;
        }
        
        if(ncol(customersData) == 2 && !customerDemandsDone && length(line) == noOfCustomers){
            customersData <- cbind(customersData, as.matrix(read.table(textConnection(line), "\n", header = FALSE)));
            colnames(customersData)[3] <- "Demand"
            
            customerDemandsDone = TRUE
            linesToRead <- 1
            next;
        }
        
        
        if(ncol(depotsData) == 3 && !depotCostDone && length(line) == noOfDepots ){
            depotsData <- cbind(depotsData, as.matrix(read.table( textConnection(line), "\n", header = FALSE)));
            colnames(depotsData)[4] <- "Cost"
            linesToRead <- 1
            next;
        }
        
        if(vehicleCost < 0){
            vehicleCost <- as.numeric(line)
            next;
        }
        
        if(line == "0" || line == "1"){
            valueIsDecimal <- line == "1"
        }
    }
    
    close(file);
    
    #' Concatanates positions of depots and customers in a single matrix.
    getAllPos <- function(){
        return(rbind(depotsData[,1:2], customersData[,1:2]) )
    }
    
    dat <- list( noOfDepots = noOfDepots, noOfCustomers = noOfCustomers, vehicleCapacity = vehicleCapacity,
                 vehicleCost = vehicleCost , valueIsDecimal = valueIsDecimal, depotsData = depotsData,
                 customersData = customersData, getAllPos = getAllPos);
    class(dat) <- "tsmData"
    
    return(dat)
}

#dat = parseData()
#cstm = dat$customersData

