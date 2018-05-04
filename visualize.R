source("ParseDataSet.R")

plotSolution <- function(dat = NULL, chosenDepots = NULL){
    if( is.null(dat) ){
        dat <- parseData()
    }
    #browser()
    cs <- dat$customersData
    dp <- dat$depotsData
    
    #plot(c(cs[1:dat$noOfCustomers, 1:2]), xlim = c(0, 100) , col="blue", axes = FALSE, ann=FALSE, type = "p", pch=16)
    #par(new=TRUE)
    #plot(c(dp[1:dat$noOfDepots, 1:2]), xlim = c(0, 100), col="black", pch=15, ann = FALSE)
    plot(c(dp[1:dat$noOfDepots, 1]), c(dp[1:dat$noOfDepots, 2]), xlim=c(0, 100), ylim=c(0,100), asp=1, col="black", pch=15, ann=FALSE)
    par(new=TRUE)
    plot(c(cs[1:dat$noOfCustomers, 1]), c(cs[1:dat$noOfCustomers, 2]), xlim=c(0, 100), ylim=c(0,100), asp=1, col="blue", pch=16, ann=FALSE)
    if(!is.null(chosenDepots)){
        par(new=TRUE)
        dp <- dp[-c(chosenDepots),]
        noCol <- dat$noOfDepots - length(chosenDepots)
        plot(c(dp[1:noCol, 1]), c(dp[1:noCol, 2]), xlim=c(0, 100), ylim=c(0,100), asp=1, col="red", pch=15, ann=FALSE)
    }
}
