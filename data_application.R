data<-read.csv2("data_COVID.csv",header=T,sep=",") 

head(data)
dim(data)
str(data)
data$data <- as.Date(data$data)
data$operatori.sanitari <- as.numeric(data$operatori.sanitari) # na here ..
data$all.inf <- sum(data$operatori.sanitari+data$altri)# tutti= oper.sanit + altri
# incomplete data first wave
datawave1.incomp<- data[1:75, 3]
# plot all the waves
#problem with plotting 724 so must cut the tail
plot(data$data[1:600],data$altri[1:600],type="l",xlab="time",ylab="number of infected", col="black",main="number of infected",pch=20)
wavebreaks <- c(data$data[75],data$data[164],data$data[475])
#abline(v=wavebreaks,lty=2,lwd=1,col="darkgreen")
window.start<-c(data$data[2],data$data[164])
window.end<-c(data$data[2+10],data$data[164+10])
abline(v=window.start,lty=1,lwd=1,col="green")
abline(v=window.end,lty=1,lwd=1,col="red")
abline(h=0,lty=3,lwd=0.5,col="darkgreen")
midwavebreaks <- c(data$data[35],data$data[350])
text(midwavebreaks[c(1,2)],c(180,180),labels=c("1","2"),col="darkgreen",lwd=2)

# waves 
Zf1 <- data$altri[1:75] # first wave
Zf2 <- data$altri[164:475] # second wave
Zf3 <- data$altri[476:593] # third wave

Zf1[5]<-1; Zf1[47]<-6;Zf1[48]<-6
Zf1<-Zf1[5:56]
#  mle
mle1<-mle(Zf1)
mle2<-mle(Zf2)
mle3<-mle(Zf3)
# add mles to plot
text(midwavebreaks[c(1,3,4)],c(300,300,300),
labels=c(paste("mle=",round(mle1,4)),paste("mle2=",round(mle2,4)),
paste("mle3=",round(mle3,4))),col="darkgreen",lwd=2)

# plot the first wave only
plot(data$data[5:56],Zf1,type="l",xlab="time",ylab="number of infected", col="black",pch=20,main="")

# simulate trajectories using mle1
Nsim<-100
Zsim<-matrix(NA,nrow=Nsim,ncol=75)
for(j in 1:Nsim){
	Zsim[j,]<-gw(p,time=75)}
plot(data$data[1:75],data$altri[1:75],type="l",xlab="time",ylab="number of infected", col="black",pch=20, lwd=2,main="sardinia - first wave")
for (i in 1:Nsim){
lines(data$data[1:75],Zsim[i,],type="l", lwd=1, col="gray",pch=20)
}
lines(data$data[1:75],apply(Zsim,2,max),type="l", lwd=2,col="orange",pch=20)

