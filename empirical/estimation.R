path = "/home/lorenzo/"
setwd(path)
################################################################################
### Packages
################################################################################
{
  # Package names
  packages <-
    c("survival",
      "foreach",
      "haven",
      "doParallel",
      "nprobust",
      "KernSmooth",
      "xtable")
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages],repos='http://cran.us.r-project.org')
  }
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
} #LOAD PACKAGES
################################################################################
### Data
################################################################################
{
#import dataframe
illinois <- as.data.frame(
  read_dta("illinois.dta"))
head(illinois)

wkpaid= as.numeric(illinois$wkpaid) #spells
#Descriptive statistics
hist(wkpaid, breaks = max(wkpaid)+1)
summary(wkpaid)
#Transforming NA in 0 for the variable "agrees to participate"
illinois$lagree[is.na(illinois$lagree)] = 0

#Creating variables zs, ws and delta
zs = vector(mode = "numeric", length = nrow(illinois)) #Treatment variable
ws = vector(mode = "numeric", length = nrow(illinois)) #Instrumental variable
deltas = vector(mode = "numeric", length = nrow(illinois)) #Censoring indicator
for (i in 1:nrow(illinois)){
  if ((illinois$hie[i]==1)&(illinois$lagree[i]==1))
  {
    zs[i] = 2
  }
  if ((illinois$jsie[i]==1)&(illinois$lagree[i]==1))
  {
    zs[i] = 1
  }
  if ((illinois$hie[i]==1))
  {
    ws[i] = 2
  }
  if ((illinois$jsie[i]==1))
  {
    ws[i] = 1
  }
  if (illinois$wkpaid[i]<26)
  {
    deltas[i] = 1
  }
}
table(deltas)/nrow(illinois)
summary(zs)
summary(ws)
(z_w = table(zs,ws))
ages = illinois$age
{###### descriptive statistics
# table = matrix(nrow = 49, ncol = 3)
# rownames(table) = colnames(illinois)
# colnames(table) = c("Control group", "JSIE", "HIE")
# table[,1] = round(colMeans(illinois[W==0,]),3)
# table[,2] = round(colMeans(illinois[W==1,]),3)
# table[,3] = round(colMeans(illinois[W==2,]),3)
# table = table[c(1,22,6,46,15,24,27, 2),]
# rownames(table) = c("Age", "Male","Black", "White", "Hispanic", "Native American", "Others", "Base period earnings" )
# #print(xtable(table, type = "latex", digits = 3), file = "desc_stat.tex")
# 
# #Average unemployment duration
# avg_dur = matrix(ncol = 8, nrow = 1)
# colnames(avg_dur) = c("All", "Control group", "JSIE", "Complier JSIE", "Non-complier JSIE", "HIE", "Complier HIE", "Non-complier HIE")
# avg_dur[1,1] = round(mean(wkpaid),2)
# avg_dur[1,2] = round(mean(wkpaid[W==0]),2)
# avg_dur[1,3] = round(mean(wkpaid[W==1]),2)
# avg_dur[1,6] = round(mean(wkpaid[W==2]),2)
# avg_dur[1,4] = round(mean(wkpaid[Z==1]),2)
# avg_dur[1,7] = round(mean(wkpaid[Z==2]),2)
# avg_dur[1,5] = round(mean(wkpaid[W==1&Z!=1]),2)
# avg_dur[1,8] = round(mean(wkpaid[W==2&Z!=2]),2)
# 
# #print(xtable(avg_dur, type = "latex", digits = 2), file = "avg_dur.tex")
# #T-tests
# t.test(wkpaid[Z==0],wkpaid[Z==1])
# t.test(wkpaid[Z==0],wkpaid[Z==2])
}
cutoff = 26
wkpaid = pmin(wkpaid, cutoff)

ys = wkpaid
hist(ys)
xs = (ages  - mean(ages))/sd(ages)
hist(xs)
illinois$male
data = as.data.frame(cbind(ys,deltas,xs,zs,ws,illinois$black,illinois$male))
names(data) = c("y","delta","x","z","w","black","male")

}
################################################################################
### Estimations Functions
################################################################################
{
#epanechnikov function
epanechnikov = function(x){as.numeric(abs(x)<1)*3*(1-x^2)/4}#function(x){dbeta(x+.5,2,2)}#
uniform = function(x){as.numeric(abs(x)<1)*0.5}
comp.time =function(str.code){
  ctime = Sys.time()
  eval(parse(text=str.code))
  print(Sys.time()-ctime)
}



# conditional cdf of t|z,x,w. bw is the bandwidth. When j is not NA, means x = xs[j].
Ft_zxw = function(t,z,x,w,j=NA){
  if(t<=0){return(0)}
  iw = matws[which(rangew==w),]
  iz = matzs[which(rangez==z),]
  if(is.na(j)){
    nww = sapply(1:n,function(i){kernel((x-xs[i])/h)})*iw*iz #further division by h is useless
  }else{
    nww = xkern[j,]*iw*iz
  }
  if(sum(nww)==0){return(0)}
  nww = nww/sum(nww)
  dens = c(matys %*% matrix(nww,ncol=1))
  not_zeros = which(dens*as.numeric(ys<=t)*deltas>0)
  1- prod(1-nww[not_zeros]/dens[not_zeros])
}

# probability of z|x,w. If j is not NA, means x = xs[j]
pzxw = function(z,x,w,j=NA){
  iw = matws[which(rangew==w),]
  iz = matzs[which(rangez==z),]
  if(is.na(j)){
    kx = sapply(xs, function(xi) kernel((x-xi)/h)) #further division by h is useless
  }else{
    kx = xkern[j,]
  }
  
  res = sum(iz*iw*kx)/sum(iw*kx)
  ifelse(is.finite(res),res,0)
}

# conditional cdf of t,z|x,w. When j is not NA, means x = xs[j].
Ftz_xw = function(t,z,x,w,j=NA){
  Ft_zxw(t,z,x,w,j)*pzxw(z,x,w,j)
}


# solve for the inversoe of the fucntion f
inverse = function (f, lower = -100, upper = 100) {
  function (y) as.numeric(uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1])
}

estimate.phi = function(x,u,j =NA,z){
  res = tryCatch({
    phi0 = inverse(function(t){Ftz_xw(t,0,x,0,j=j)},lower = 0,upper = 100)(u)
    if(z==0){
      res = phi0
      }else if(z == 1) {
      auxfun = function(t){Ftz_xw(t,0,x,1,j=j)}
      point = u-auxfun(phi0)
      phi1 = ifelse(point<=0,NA,inverse(function(t){Ftz_xw(t,1,x,1,j)},lower = 0,upper = 100)(point))
      res = phi1
      }else{
        auxfun = function(t){Ftz_xw(t,0,x,2,j=j)}
        point = u-auxfun(phi0)
        phi2 = ifelse(point<=0,NA,inverse(function(t){Ftz_xw(t,2,x,2,j)},lower = 0,upper = 100)(point))
        res = phi2
    }
  },error=function(cond){return(NA)})
  return(res)
}
}

################################################################################
### EParameter setting
################################################################################
kernel.type = "epanech"
kernel = epanechnikov
data = data[data$black == 1 & data$male==0,]
table(data$w, data$z)
(n = nrow(data))
set.seed(1)
data$u = runif(n)
lower = 0; upper = 1; # lower and upper bound for beta set
rangez = c(0,1,2); # range of z
rangew = c(0,1,2); # range of w
verbose = FALSE
num_boostrap = 5000 #number of bootstrap sample
samples = rbind(1:n,#corrspond the the estimation
                t(sapply(1:num_boostrap,function(i)sample.int(n,n,replace = TRUE))))
################################################################################
### Estimation
################################################################################

data$z1 = as.numeric(data$z==1); data$z2=as.numeric(data$z==2)
cox.res = coxph(Surv(y,delta)~x+z1+z2,data = data)
print(summary(cox.res))
ubar = 0.6
results = c()

cores=detectCores()
cl <- makeCluster(cores[1]-6) #not to overload your computer
registerDoParallel(cl)
finalMatrix <- foreach(smp=1:nrow(samples), .combine=rbind, .packages =packages) %dopar% {
  #write(paste("ubar = ",ubar,"; running ", smp, " out of ", nrow(samples), "...", sep = ""),file="simulations_state.txt",append=TRUE)
  #retrive data
  ys = data$y[samples[smp,]]; 
  deltas = data$delta[samples[smp,]]; 
  xs = data$x[samples[smp,]];
  zs = data$z[samples[smp,]];
  ws = data$w[samples[smp,]];
  us = data$u[samples[smp,]];
  z1s = as.numeric(zs==1);  z2s = as.numeric(zs==2);
  nprobust.h = kdbwselect(xs,eval = xs)$bws
  hs = nprobust.h[,2]
  min.ys = 0; max.ys = max(ys)
  xkern = t(sapply(1:n, function(i) kernel((xs - xs[i])/hs[i])))
  matys = sapply(1:n, function(i) as.numeric(ys[i] >= ys))
  matws = rbind(as.numeric(ws==0),as.numeric(ws==1),as.numeric(ws==2))
  matzs = rbind(as.numeric(zs==0),as.numeric(zs==1),as.numeric(zs==2))
  vs = pmin(ubar,us)
  Deltas = as.numeric(us<=ubar)
  tzs = sapply(1:n, function(j) estimate.phi(x = xs[j],u = vs[j],j = j,z = zs[j]) )
  sum(is.na(tzs))
  nonas = !is.na(tzs)
  est.beta <- coxph(Surv(tzs[nonas], Deltas[nonas]) ~ xs[nonas]+z1s[nonas]+z2s[nonas])$coefficients
  est.beta
  return(c(smp,smp == 1, est.beta))
}
stopCluster(cl)
finalMatrix = as.data.frame(finalMatrix)
names(finalMatrix) = c("id","isEst","x","z1","z2")
write.csv(finalMatrix,"empirical_application_results.csv",row.names = F)
head(finalMatrix)

est = finalMatrix[finalMatrix$isEst==1,c("x","z1","z2")]
sds = apply(finalMatrix[finalMatrix$isEst==0,c("x","z1","z2")], 2, sd)
CL95 = est -1.96*sds
CU95 = est +1.96*sds
res = cbind(as.numeric(est),
            as.numeric(sds),
            as.numeric(CL95),
            as.numeric(CU95),
            as.numeric(cox.res$coefficients),
            as.numeric(sqrt(diag(cox.res$var))),
            as.numeric(cox.res$coefficients-1.96*sqrt(diag(cox.res$var))),
            as.numeric( cox.res$coefficients+1.96*sqrt(diag(cox.res$var))))
colnames(res) = c("est","sd","CI025","CI975","cox.est","cox.sd","cox.CI025","cox.CI975")
rownames(res) = c("x","z1","z2")

write.csv(res,"estimation.csv")
res
library(xtable)
xtable(res,digits = c(0,rep(3,8)) )
