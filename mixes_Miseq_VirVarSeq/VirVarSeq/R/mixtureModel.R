# Project: IlluminaPlasmide
# 
# Author: BVERBIS2
###############################################################################
#install.packages("/VirVar/R/packages/rmgt_0.9.001.tar.gz",repos=NULL,lib="/VirVar/R/lib")


#load arguments from pipeline
arguments <- commandArgs(trailingOnly=TRUE)
stats <- arguments [1]
png <- arguments[2]
QScore <- arguments[3]

#load necessary functions
library(rmgt)

#Read stats file and extract necessary information 
x <- read.delim(stats,header=FALSE)[-c(1:6),]
x <- x[-c((length(x)-1):length(x))]

dat <- as.factor(unlist(lapply(x, function(x) strsplit(as.character(x),split=" ")[[1]][2])))
dat <- as.numeric(gsub("\\]","",gsub("\\[","",dat)))

names(dat) <- as.factor(unlist(lapply(x, function(x) strsplit(as.character(x),split=" ")[[1]][4])))
names(dat) <- gsub("\\]","",gsub("\\[","",names(dat)))


#Check upper and lower limit of quality score range
max <- max(as.numeric(names(dat)))
min <- min(as.numeric(names(dat)))
bins <- length(seq(min,max))


#Fill possible gaps (quality scores where zero counts are observed)
if (length(dat)==bins) {
	dat <- dat
} else {
	ID <- as.numeric(seq(min,max))[!(as.character(seq(min,max)) %in% names(dat))]
	tmp <- rep(0,times=bins)
	tmp[-(ID-1)] <- dat
	dat <- tmp
	names(dat) <- seq(min,max)
}


#Presence of point probability? If yes, low variability.  
Quantile0.01 <- quantile(as.numeric(rep(names(dat), dat)),0.01)
varPoint <- unname(ifelse(Quantile0.01 >= 3, 40, 0.8)) 


#Determine parameters of mixture distribution
TruncFit <- mgt(
		nx = dat,
		truncated=TRUE,
		nr = bins,
		ng = 3,
		breaks = (min-0.5) + seq (0,bins),
		avg.est = c(min,10,35),
		var.est = c(varPoint,40,40),
		p.est = c(0.05,0.10,0.85),
		iter = 2000
)	



png(png)
barplot(dat/sum(dat),space=0,names.arg=names(dat),las=2,col="light grey")
#cumulative probability in point 40
z3 <- pnorm(40, mean=TruncFit$vals$xbar[3],sd=sqrt(TruncFit$vals$var[3]))
z2 <- pnorm(40, mean=TruncFit$vals$xbar[2],sd=sqrt(TruncFit$vals$var[2]))
z1 <- pnorm(40, mean=TruncFit$vals$xbar[1],sd=sqrt(TruncFit$vals$var[1]))
normC <- z1*TruncFit$vals$pr[1]+z2*TruncFit$vals$pr[2]+z3*TruncFit$vals$pr[3]

x_continuous <- seq (min, max, length=200)
y3 <- dnorm(x_continuous, mean=TruncFit$vals$xbar[3],sd=sqrt(TruncFit$vals$var[3]))
y2 <- dnorm(x_continuous, mean=TruncFit$vals$xbar[2],sd=sqrt(TruncFit$vals$var[2]))
y1 <- dnorm(x_continuous, mean=TruncFit$vals$xbar[1],sd=sqrt(TruncFit$vals$var[1]))
lines(x_continuous-(min-0.5),(y1*TruncFit$vals$pr[1]+y2*TruncFit$vals$pr[2]+y3*TruncFit$vals$pr[3])/normC)
lines(x_continuous-(min-0.5),y3*TruncFit$vals$pr[3]/normC,type="l",lwd=2,col="green",lty=2)
lines(x_continuous-(min-0.5),y2*TruncFit$vals$pr[2]/normC,type="l",lwd=2,col="red",lty=2)
lines(x_continuous-(min-0.5),y1*TruncFit$vals$pr[1]/normC,type="l",lwd=2,col="grey",lty=2)
cutIntersect <- ifelse(TruncFit$vals$xbar[2] > TruncFit$vals$xbar[1],
		uniroot(function(y)TruncFit$vals$pr[2]/normC*dnorm(y,TruncFit$vals$xbar[2],sqrt(TruncFit$vals$var[2]))-TruncFit$vals$pr[3]/normC*dnorm(y,TruncFit$vals$xbar[3],sqrt(TruncFit$vals$var[3])), interval=c(min, max))$root,
		uniroot(function(y)TruncFit$vals$pr[1]/normC*dnorm(y,TruncFit$vals$xbar[1],sqrt(TruncFit$vals$var[1]))-TruncFit$vals$pr[3]/normC*dnorm(y,TruncFit$vals$xbar[3],sqrt(TruncFit$vals$var[3])), interval=c(min, max))$root)
#optimise(f=function(y)(TruncFit$vals$pr[2]/normC*dnorm(y,TruncFit$vals$xbar[2],sqrt(TruncFit$vals$var[2]))-TruncFit$vals$pr[3]/normC*dnorm(y,TruncFit$vals$xbar[3],sqrt(TruncFit$vals$var[3])))^2,c(2:40),maximum=FALSE)$minimum
abline(v=cutIntersect-(min-0.5),lty=2,col="grey")
dev.off()

write(round(cutIntersect), QScore)

q()
