source("ParseDataSet.R")
library("ggplot2")
plotSolution <- function(dat = NULL, chosenDepots = NULL){
    if( is.null(dat) ){
        dat <- parseData()
    }
    
    cs <- dat$customersData[,c(1:2, 4)]
    dp <- dat$depotsData[,c(1:2, 5)]
    dt <- rbind(cs,dp)
    
    
    # df <- data.frame(cs[,1:2])
    df <- data.frame(dt)
    df$id <- as.factor(df$id)
    df$X <- as.numeric(df$X)
    df$Y <- as.numeric(df$Y)

    plt <- ggplot(df, aes(x=X, y=Y, color=id)) + geom_point()  +  scale_x_discrete(name="X") + scale_y_discrete(name="Y") + coord_fixed()
    
    # plt <- plt + scale_x_continuous(name="X", limits=c(0, 100)) + scale_y_continuous(name="Y", limits=c(0, 100))
    plt
    #plot(c(cs[1:dat$noOfCustomers, 1:2]), xlim = c(0, 100) , col="blue", axes = FALSE, ann=FALSE, type = "p", pch=16)
    #par(new=TRUE)
    #plot(c(dp[1:dat$noOfDepots, 1:2]), xlim = c(0, 100), col="black", pch=15, ann = FALSE)
    # plot(c(dp[1:dat$noOfDepots, 1]), c(dp[1:dat$noOfDepots, 2]), xlim=c(0, 100), ylim=c(0,100), asp=1, col="black", pch=15, ann=FALSE)
    # par(new=TRUE)
    # plot(c(cs[1:dat$noOfCustomers, 1]), c(cs[1:dat$noOfCustomers, 2]), xlim=c(0, 100), ylim=c(0,100), asp=1, col="blue", pch=16, ann=FALSE)
    # if(!is.null(chosenDepots)){
    #     par(new=TRUE)
    #     dp <- dp[-c(chosenDepots),]
    #     noCol <- dat$noOfDepots - length(chosenDepots)
    #     plot(c(dp[1:noCol, 1]), c(dp[1:noCol, 2]), xlim=c(0, 100), ylim=c(0,100), asp=1, col="red", pch=15, ann=FALSE)
    # }
    
}
