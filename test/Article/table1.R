library(MDSV)
if(!require(timeDate)){install.packages("timeDate")}; library(timeDate)

index_set<-c("aex", "aord", "bell", "bsesn", "bvsp", "cac40", "dax", "dji", "ftmib", "ftse",
             "hsi", "ibex35", "kospi", "kse", "mxx", "n225", "nasdaq", "nifty50", "omxc20", "omxhpi",
             "omxspi", "oseax", "psi", "rut", "smsi", "sp500", "ssec", "ssmi", "sti", "stoxx50", "tsx")

start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")
#-------------------------------------------------------------------------------
# Descriptive statitics
#-------------------------------------------------------------------------------

des_stats<-NULL
for(index in index_set){
  donne<-get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  ech<-donne[,"r"]
  ech<-ech-mean(ech)
  q<-quantile(ech)
  #descriptive statistics
  des_stat <-    c("N" = length(ech),
                   "Min." = q[[1]],
                   "Q1" = q[[2]],
                   "Mean" = mean(ech),
                   "Median" = q[[3]],
                   "StDev" = sd(ech),
                   "Skewness" = timeDate::skewness(ech, method="moment"),
                   "Kurtosis" = timeDate::kurtosis(ech, method="moment"),
                   "Q3" = q[[4]],
                   "Max." = q[[5]]
  )
  des_stats<-rbind(des_stats,des_stat)
  ech<-donne[,"rv"]
  q<-quantile(ech)
  #descriptive statistics
  des_stat <-    c("N" = length(ech),
                   "Min." = q[[1]],
                   "Q1" = q[[2]],
                   "Mean" = mean(ech),
                   "Median" = q[[3]],
                   "StDev" = sd(ech),
                   "Skewness" = timeDate::skewness(ech, method="moment"),
                   "Kurtosis" = timeDate::kurtosis(ech, method="moment"),
                   "Q3" = q[[4]],
                   "Max." = q[[5]]
  )
  des_stats<-rbind(des_stats,des_stat)
}

rownames(des_stats)<-rep(index_set,1,each=2)
write.csv(round(des_stats,3), "stat.csv")


View(round(des_stats,2))

#-------------------------------------------------------------------------------
# Log return and realized variance graphs
#-------------------------------------------------------------------------------
pdf( file ="RV.pdf", width =6.5 , height =4)
index_set<-c("S&P 500","NASDAQ 100", "FTSE 100", "NIKKEI 225")
par(mfrow=c(2,2),mar=c(2.5, 2.5, 2, 2))
for(index in c("sp500","nasdaq","ftse","n225")){
  donne<-get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_<-rownames(donne)
  ech<-donne[,"rv"]
  
  # Traitement des dates
  t_ <- as.POSIXlt(c(as.character(start.date),nam_), format ="%Y -%m -%d")
  n.days <- rep (365,length (t_))
  n.days [((1900+ t_$ year) %% 4)==0] <- 366
  t_ <- 1900 + t_$ year + (t_$ yday +1)/n.days
  
  x.start <- t_[1] # date de début pour les figures
  i.start <- which.min(abs(t_-x.start ))
  x.end <- t_[length (t_)] # date de fin pour les figures
  i.end <- which.min(abs(t_-x.end))
  t_<- t_[(i.start+1):i.end]
  donne <- donne [i.start:(i.end-1),]
  
  index<-index_set[which(index == c("sp500","nasdaq","ftse","n225"))]

  t.at <- seq(from = ceiling(x.start ),to=floor(x.end),by =2)
  xlim <- range(t_)
  cex.axis <- 0.85
  cex.mtext <- 0.85
  
  plot(t_,donne[,"rv"],type="l",xlab="", ylab="", xaxt ="n",main="", xlim =xlim , xaxs ="i")
  axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
  mtext (index, side =3, line=0.5, font=2, las=0, cex=cex.mtext)
  
}
dev.off ()
