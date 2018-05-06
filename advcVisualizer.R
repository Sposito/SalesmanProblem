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

    # ct <- 0
    # for(listRout in dat$routes){
    #     ct <- ct + 1
    #     if(is.null(i)){
    #         next
    #     }
    #     for (j in i){
    #         route < matrix(
    #         for (k in j){
    #         
    #         }
    #     }
    # 
    # }
    browser()
    plt <- ggplot(df, aes(x=X, y=Y, color=id)) + geom_point()  +  scale_x_discrete(name="X") + scale_y_discrete(name="Y") + coord_fixed()
    # plt = plt + geom_segment(aes(x = 50 ,y = 0, xend = 50 , yend = 50, colour = as.factor(9) ))
    # plt = plt + geom_segment(aes(x = 30 ,y = 0, xend = 30 , yend = 50, colour = as.factor(10) ))
    plt
    # for(i in 1:dat$noOfDepots){
    #     size = length(dat$bestChrm[i])
    #     if(size == 0)
    #         next
    # 
    #     # plt = plt +  geom_segment(aes(x = dat$depotsData[i,1], y = dat$depotsData[i,2], xend = dat$customersData[dat$bestChrm[i][1],1],
    #     #                                 yend = dat$customersData[dat$bestChrm[i][1],2], colour = as.factor(i) ), inherit.aes = TRUE)
    #     #plt = plt + geom_segment(aes(x = i *5 ,y = 0, xend = 90 , yend = 90, colour = as.factor(i) ))
    #     for (j in 1:size){}
    # }
    # plt <- plt + scale_x_continuous(name="X", limits=c(0, 100)) + scale_y_continuous(name="Y", limits=c(0, 100))
    
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
