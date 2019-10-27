#PA lung cancer data
library(stats)
library(sp)
library(spdep)
library(GISTools)
library(SpatialEpi)
library(MASS)

#2. Generate the number of cases by county assuming a population of 1000 per county with a common rate of 0.50 using the rpois function.  
#i. Display these with a chloropleth map.  
set.seed(2090)
data(pennLC)
PA.agg <- aggregate(pennLC$data[,2:3],by=list(pennLC$data$county), FUN=sum, na.rm=TRUE)
names(PA.agg)[1] <- "county"
#Generate number of cases
length(PA.agg$county)
PA.agg[,3] <-rep(1000,67)
PA.agg[,2] <-rpois(67,0.5)
#View(PA.agg)
#Calculate cancer rates by county
PA.agg$rate <-PA.agg$cases/PA.agg$population
#Plot choropleth map
PA.latlong <- pennLC$spatial.polygon
PA.utm <- spTransform(PA.latlong,CRS("+init=epsg:3724 +units=km"))
choropleth(PA.utm,PA.agg$rate)

#ii. Test for constant risk using Geary's c (MC approach) 
#Geary's C test for constant risk
PA.nb <- poly2nb(PA.utm,queen=FALSE)
PA.lw <- nb2listw(PA.nb)
geary.test(PA.agg$rate,PA.lw)
geary.mc(PA.agg$rate,PA.lw,nsim=999)

#class(PA.nb)
#unlist(PA.nb[1])
#as.numeric(PA.nb[1])

#3.i. Assume that all counties sharing a border are correlated with each other at a level 0.80 (this is a population parameter), while counties not sharing a border are correlated at 0.10. 
#View(PA.nb)
#Generate the correlation matrix. If the ith row in the PA.nb list contains the jth county, then assign 0.8 to item [i,j] in the correlation matrix. If not, assign 0.1.
cor.mat<-matrix(nrow=67,ncol=67)
for (i in 1:length(PA.nb)){
  for (j in 1:length(PA.agg$county)){
    cor.mat[i,j]<-ifelse(j%in%as.numeric(unlist(PA.nb[i])),0.8,0.1)
  }
}
#Assign 1 to diagonal elements
for (i in 1:nrow(cor.mat)){
  cor.mat[i,i]<-1
}

#Approximate to nearest positive definite matrix
cor.mat<-nearPD(cor.mat,corr=T)
cor.mat$mat
cov.mat<-cor.mat$mat*0.5#approximated process variance is 0.5
mean<-rep(0.5,67)#approximated process mean is 0.5
Sigma <- cov.mat
#generate z realization of multivariate normal distribution, which is the number of cases
z<-mvrnorm(n = 1,mean, Sigma)

#ii. Display a chloropleth map with the realized case counts.
PA.agg$newcases<-z
#Calculate new cancer rates by county
PA.agg$newrate <-PA.agg$newcases/PA.agg$population
#Plot choropleth map with realized case counts 
choropleth(PA.utm,PA.agg$newrate)

#iii. Test for constant risk using Geary's c (MC approach) 
#Geary's C test for constant risk
geary.test(PA.agg$newrate,PA.lw)
geary.mc(PA.agg$newrate,PA.lw,nsim=999)

#4. Repeat (3) with correlation of 0.30 across borders; 0.10 when no border. 
#Generate the correlation matrix. If the ith row in the PA.nb list contains the jth county, then assign 0.3 to item [i,j] in the correlation matrix. If not, assign 0.1.
cor.mat2<-matrix(nrow=67,ncol=67)
for (i in 1:length(PA.nb)){
  for (j in 1:length(PA.agg$county)){
    cor.mat2[i,j]<-ifelse(j%in%as.numeric(unlist(PA.nb[i])),0.3,0.1)
  }
}
#Assign 1 to diagonal elements
for (i in 1:nrow(cor.mat2)){
  cor.mat2[i,i]<-1
}

#Approximate to nearest positive definite matrix
cor.mat2<-nearPD(cor.mat2,corr=T)
cor.mat2$mat
cov.mat2<-cor.mat2$mat*0.5#approximated process variance is 0.5
#generate z realization of multivariate normal distribution, which is the number of cases
z2<-mvrnorm(n = 1,mean,cov.mat2)

#ii. Display a chloropleth map with the realized case counts.
#Calculate new cancer rates by county
PA.agg$newrate2 <-z2/PA.agg$population
#Plot choropleth map with realized case counts 
choropleth(PA.utm,PA.agg$newrate2)

#iii. Test for constant risk using Geary's c (MC approach) 
#Geary's C test for constant risk
geary.test(PA.agg$newrate2,PA.lw)
geary.mc(PA.agg$newrate2,PA.lw,nsim=999)
