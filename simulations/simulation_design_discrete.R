path = "/home/lorenzo/"
setwd(path)
################################################################################
### Packages
################################################################################
{
  # Package names
  packages <-
    c("abind",
      "foreach",
      "survival",
      "dplyr",
      "xtable",
      "nprobust",
      "doParallel")
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages],
                     repos='http://cran.us.r-project.org',INSTALL_opts = '--no-lock')
  }
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
}

################################################################################
#### Define estimation functions
################################################################################

{
  
  Fucens = function(u,umax=0.7){ifelse(u>umax,1,ifelse(u<0,0,(1-(1-u/umax)/(1-u))))}
  invFucens = function(x,umax = 0.7){umax*x/(1+umax*(x-1))}
  
  
  # Generate the data given the design 
  generate.data = function(n,verbose = TRUE,beta0,inv.Lambda0, 
                           wdens,
                           zdens,
                           xdens,
                           cdens,
                           ubar
  ){
    us = runif(n)  
    ws = eval(parse(text=wdens))
    xs = eval(parse(text=xdens))
    zs = eval(parse(text=zdens))
    ts = inv.Lambda0(-log(1-us)*exp(-beta0[1]*zs - beta0[2]*xs))
    cs = eval(parse(text=cdens))
    ys = pmin(ts,cs)
    deltas = as.numeric(ts<=cs)
    if(verbose){
      print(table(zs,ws)/n)
      print(paste(  (sum(deltas==0)/n*100), "% censoring",sep = ""))
    }
    us = runif(n)
    vs = pmin(us,ubar)
    Deltas = as.numeric(us<=ubar)
    data = cbind(ys,deltas,zs,xs,ws,vs,Deltas)
    data
  }

  comp.time =function(str.code){
    ctime = Sys.time()
    eval(parse(text=str.code))
    print(Sys.time()-ctime)
  }

  

  # conditional cdf of t|z,x,w. bw is the bandwidth. When j is not NA, means x = xs[j].
  Ft_zxw = function(t,z,x,w){
    if(t<= 0){return(0)}
    index = which(c("000","001","011","100","101","111")==paste(x,z,w,sep=""))
    val = hatF.funs[[index]](t)
    return(ifelse(is.na(val),1,val))
  }
  # true value of phi
  phi = function(u,z,x){inv.Lambda0(-log(1-u)*exp( -sum(c(z,x)*beta0)))}
  
  # probability of z|x,w. If j is not NA, means x = xs[j]
  pzxw = function(z,x,w){
    index = which(c("000","001","011","100","101","111")==paste(x,z,w,sep=""))
    return(probs[index])
  }
  
  # conditional cdf of t,z|x,w. When j is not NA, means x = xs[j].
  Ftz_xw = function(t,z,x,w){ 
    Ft_zxw(t,z,x,w)*pzxw(z,x,w)
  }

  # solve for the inversoe of the fucntion f
  inverse = function (f, lower = -100, upper = 100) {
    function (y) as.numeric(uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1])
  }
  
  # create the values of phi in a point x, for each z, at points c(0,grid.u) 
  
  estimate.phi = function(x,u,z){
    res = tryCatch({
      phi0 = inverse(function(t){Ftz_xw(t,0,x,0)},lower = 0,upper = 10)(u)
      if(z==0){res = phi0}else{
        auxfun = function(t){Ftz_xw(t,0,x,1)}
        point = u-auxfun(phi0)
        phi1 = ifelse(point<=0,NA,inverse(function(t){Ftz_xw(t,1,x,1)},lower = 0,upper = 10)(point))
        res = phi1
      }
    },error=function(cond){return(NA)})
    return(res)
  }
}
################################################################################
### Simulations
################################################################################
#zdens, xdens, wdens, beta01, beta02, cdens, inv.Lambda0
designs = rbind(
  c("as.numeric(ws==1)*as.numeric(.65-ws+.5*us+.5*rnorm(n)+.3*xs>=0)", "as.numeric(runif(n)>0.5)","as.numeric(runif(n)>0.5)", .7, .7,"rexp(n,.43)", "function(t){t}"),#40%cens
  c("as.numeric(ws==1)*as.numeric(.65-ws+.5*us+.5*rnorm(n)+.3*xs>=0)", "as.numeric(runif(n)>0.5)","as.numeric(runif(n)>0.5)", .7, .7,"rexp(n,1.15)", "function(t){t}")#20%cens
)

allResults = c()
lower = 0; upper = 1; # lower and upper bound for beta set
rangez = c(0,1); # range of z
rangew = c(0,1); # range of w
verbose = FALSE; # boolean for printing intermidate results
ubar = .9
ctime = Sys.time()

N = 500;
for(design in c(1,2)){
  for(n in c(500,1000)){
    #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  finalMatrix <- foreach(id=1:N, .combine=rbind, .packages =packages) %dopar% {
    #write(paste("ubar = ",ubar,"; design = ",design, "; n = ",n,"; running ", id, " out of ", N, "...", sep = ""),file="simulations_state.txt",append=TRUE)
    set.seed(id)
    zdens = designs[design,1];xdens=designs[design,2];wdens=designs[design,3];
    beta0 = as.numeric(designs[design,4:5]); cdens = designs[design,6]
    inv.Lambda0 = eval(parse(text=designs[design,7]))
    data = generate.data(n=n,inv.Lambda0 = inv.Lambda0,
                         beta0 = beta0,cdens = cdens,
                         wdens = wdens,zdens = zdens,xdens = xdens,verbose = TRUE,ubar = ubar)
    est.beta = tryCatch({
      # Extract columns of data and assign them to variables
      ys = data[,1]; deltas = data[,2]; zs = data[,3]; xs = data[,4]; ws = data[,5]; vs =data[,6]; Deltas =data[,7] 
      
      # Define categories
      categories = c("000","001","011","100","101","111")
      nctg = length(categories);
      
      # Create matrix of categories
      categories.matrix = rbind(c(0,0,0),c(0,0,1),c(0,1,1),c(1,0,0),c(1,0,1),c(1,1,1))
      
      # Calculate the probability of each category
      probs = array(0,dim = nctg)
      for(ctg in 1:nctg){
        x = categories.matrix[ctg,1];z = categories.matrix[ctg,2]; w = categories.matrix[ctg,3]
        probs[ctg] = sum(xs == x & zs==z & ws == w)/sum(ws==w & xs == x)
      }                         
      
      # Display the probabilities of each category
      probs
      
      # Fit a Kaplan-Meier survival curve with the specified formula
      KM=survfit(Surv(ys, deltas) ~ xs+zs+ws)
      
      # Set minimum and maximum values for ys and create a grid based on it
      min.ys = 0; max.ys = max(ys)
      grid.y = seq(0,max.ys,by = 0.01)
      
      # Create a matrix of values for the survival function
      hatF.points = matrix(0,nrow =length(grid.y),ncol = length(categories))
      for(i in 2:length(grid.y)){
        hatF.points[i,] = (1-summary(KM,times = grid.y[i],extend = TRUE)$surv)
      }
      
      # Generate a list of functions to approximate each component of the survival function for each category
      hatF.funs = list()
      for(j in 1:nctg){
        indices = which(duplicated(hatF.points[,j])==FALSE)
        hatF.funs[[j]] = approxfun(x = grid.y[indices], y = hatF.points[indices,j])
      }

      tzs = sapply(1:n, function(j) estimate.phi(x = xs[j],u = vs[j],z = zs[j]) )
      sum(is.na(tzs))
      nonas = !is.na(tzs)
      est.beta <- coxph(Surv(tzs[nonas], Deltas[nonas]) ~ zs[nonas]+xs[nonas])$coefficients
      est.beta
  
    },error = function(cond){return(c(NA,NA))})
    est.beta
    (cox.beta = coxph(Surv(ys, deltas) ~ zs+xs)$coefficients)
    #bootstrap
    bdata = data[sample.int(n,n,replace = TRUE),]
    est.beta.boot = tryCatch({
      # Extract columns of data and assign them to variables
      ys = bdata[,1]; deltas = bdata[,2]; zs = bdata[,3]; xs = bdata[,4]; ws = bdata[,5]; vs =bdata[,6]; Deltas = bdata[,7];
      # Define categories
      categories = c("000","001","011","100","101","111")
      nctg = length(categories);
      
      # Create matrix of categories
      categories.matrix = rbind(c(0,0,0),c(0,0,1),c(0,1,1),c(1,0,0),c(1,0,1),c(1,1,1))
      
      # Calculate the probability of each category
      probs = array(0,dim = nctg)
      for(ctg in 1:nctg){
        x = categories.matrix[ctg,1];z = categories.matrix[ctg,2]; w = categories.matrix[ctg,3]
        probs[ctg] = sum(xs == x & zs==z & ws == w)/sum(ws==w & xs == x)
      }                         
      
      # Display the probabilities of each category
      probs
      
      # Fit a Kaplan-Meier survival curve with the specified formula
      KM=survfit(Surv(ys, deltas) ~ xs+zs+ws)
      
      # Set minimum and maximum values for ys and create a grid based on it
      min.ys = 0; max.ys = max(ys)
      grid.y = seq(0,max.ys,by = 0.01)
      
      # Create a matrix of values for the survival function
      hatF.points = matrix(0,nrow =length(grid.y),ncol = length(categories))
      for(i in 2:length(grid.y)){
        hatF.points[i,] = (1-summary(KM,times = grid.y[i],extend = TRUE)$surv)
      }
      
      # Generate a list of functions to approximate each component of the survival function for each category
      hatF.funs = list()
      for(j in 1:nctg){
        indices = which(duplicated(hatF.points[,j])==FALSE)
        hatF.funs[[j]] = approxfun(x = grid.y[indices], y = hatF.points[indices,j])
      }
      tzs = sapply(1:n, function(j) estimate.phi(x = xs[j],u = vs[j],z = zs[j]) )
      sum(is.na(tzs))
      nonas = !is.na(tzs)
      est.beta.boot <- coxph(Surv(tzs[nonas], Deltas[nonas]) ~ zs[nonas]+xs[nonas])$coefficients
      est.beta.boot
    },error = function(cond){return(c(NA,NA))})
    est.beta.boot
    (cox.beta.boot = coxph(Surv(ys, deltas) ~ zs+xs)$coefficients)
    return(cbind(design,id,1:2,sum(deltas==0)/n,n,est.beta,cox.beta,beta0,est.beta.boot,cox.beta.boot))
  }
  #stop cluster
  stopCluster(cl)
  Sys.time()-ctime
  results = as.data.frame(finalMatrix)
  names(results) = c("design","id","index","cens","n","est","cox","b0","estboot","coxboot")
  allResults = rbind(allResults,results)
  }
}    


stats = c()
for(design in unique(allResults$design)){
  for(n in unique(allResults$n)){
    current = allResults[allResults$n==n & allResults$design==design, ]
    bias = c(mean(current[current$index==1,c("est")]-current[current$index==1,c("b0")]),
             mean(current[current$index==2,c("est")]-current[current$index==2,c("b0")]))
    biascox = c(mean(current[current$index==1,c("cox")]-current[current$index==1,c("b0")]),
                mean(current[current$index==2,c("cox")]-current[current$index==2,c("b0")]))
    sds = c(sd(current[current$index==1,c("estboot")]-current[current$index==1,c("b0")]),
            sd(current[current$index==2,c("estboot")]-current[current$index==2,c("b0")]))
    sdscox = c(sd(current[current$index==1,c("coxboot")]-current[current$index==1,c("b0")]),
               sd(current[current$index==2,c("coxboot")]-current[current$index==2,c("b0")]))
    mses = bias^2 + sds^2
    msescox =  biascox^2 + sdscox^2
    rmse = sqrt(mean((current[current$index==1,c("est")]-current[current$index==1,c("b0")])^2 + 
                       (current[current$index==2,c("est")]-current[current$index==2,c("b0")])^2
    ))
    rmsecox = sqrt(mean((current[current$index==1,c("cox")]-current[current$index==1,c("b0")])^2 + 
                          (current[current$index==2,c("cox")]-current[current$index==2,c("b0")])^2))
    q_025 = c(quantile(current[current$index==1,c("estboot")]-current[current$index==1,c("est")],0.025),
              quantile(current[current$index==2,c("estboot")]-current[current$index==2,c("est")],0.025))
    q_975 = c(quantile(current[current$index==1,c("estboot")]-current[current$index==1,c("est")],0.975),
              quantile(current[current$index==2,c("estboot")]-current[current$index==2,c("est")],0.975))
    
    CP95 = c( sum((current[current$index==1,c("est")]-current[current$index==1,c("b0")]) > q_025[1]
                  & (current[current$index==1,c("est")]-current[current$index==1,c("b0")]) < q_975[1]),
              sum((current[current$index==2,c("est")]-current[current$index==2,c("b0")]) > q_025[2]
                  & (current[current$index==2,c("est")]-current[current$index==2,c("b0")]) < q_975[2]))
    CP95 = CP95/ sum(current$index==1)
    
    q_025n = c(-1.96*sd(current[current$index==1,c("estboot")]-current[current$index==1,c("est")]),
               -1.96*sd(current[current$index==2,c("estboot")]-current[current$index==2,c("est")]))
    q_975n = c(1.96*sd(current[current$index==1,c("estboot")]-current[current$index==1,c("est")]),
               1.96*sd(current[current$index==2,c("estboot")]-current[current$index==2,c("est")]))
    CP95n = c( sum((current[current$index==1,c("est")]-current[current$index==1,c("b0")]) > q_025n[1]
                   & (current[current$index==1,c("est")]-current[current$index==1,c("b0")]) < q_975n[1]),
               sum((current[current$index==2,c("est")]-current[current$index==2,c("b0")]) > q_025n[2]
                   & (current[current$index==2,c("est")]-current[current$index==2,c("b0")]) < q_975n[2]))
    CP95n = CP95n/ sum(current$index==1)
    
    cens = mean(current[current$index==1,c("cens")])
    stats = rbind(stats, 
                  cbind(design,cens,n,c(1,2),bias,sds,mses,rmse,CP95,CP95n, biascox,sdscox,msescox,rmsecox,nrow(current)/2))
  }
}
stats = as.data.frame(stats)
names(stats) = c("design","cens","n","index",
                 "bias",
                 "sd",
                 "mse",
                 "rmse",
                 "CP95", 
                 "CP95n", 
                 "biascox",
                 "sdcox",
                 "msecox",
                 "rmsecox","N")
stats = stats[order(stats$design,stats$n),]
library(xtable)

stats = stats[,c("design","cens","n",
                 "bias",
                 "sd",
                 "mse",
                 "rmse",
                 "CP95n", 
                 "biascox",
                 "sdcox",
                 "msecox",
                 "rmsecox","N")]
stats
write.csv(stats, "simulation_design_discrete.csv",row.names = FALSE)
xtable(stats,digits = c(0,0,1,0,3,3,3,3,3,3,3,3,3,0) )
