#sim_parallel_functions
#Updated May 25, 2017 to reflect errors in torus wrapping found in previous version
#Updated August 23, 2017 to scale by number of conspecfics rather than conspecific density &
#to update the RSF to include \beta_3 for a quadratic conspecific term &
#to record landscape stats for each landscape generated

#############################################################################################
#############################################################################################
sampfun<-function(x){
  ## A convenience function to use with apply when sampling values from rows within a matrix.
  sample(1:length(x),1,prob=x) 
}

#############################################################################################
#############################################################################################
#' Calculate distance weight on a matrix
#'
#' Convert distances in x,y (row, col on a matrix) to distances and then weight according to a function.
#'  
#' @param celldist matrix. Should represent a two-column matrix of change in row and column ranks, respectively.
#' @param min.dist numeric. The distance to use for celldist==0 (i.e., for the cell of origin).
#' @param dist.weight.fun function. A function to use for calculating the distance weight (e.g., dgamma)
#' @param ... Other parameters to pass to dist.weight.fun.
#' @author James Forester
#' @export
get.dist.weight <-function(celldist, min.dist=0.5, dist.weight.fun=NULL,...){
  ##calculate distance from central cell to others in the neighborhood
  celldist=sqrt(rowSums(celldist^2)) #sqrt((x-0)^2+ (y-0)^2)
  celldist[celldist==0]=min.dist
  
  if(is.function(dist.weight.fun)){
    ##calculate weight, assuming a uniform bearing distribution
    cell.dist.weight=dist.weight.fun(celldist,...)/(2*pi*celldist)
  }else{
    ##default to inverse distance weight
    cell.dist.weight=1/(2*pi*celldist^2)
  }
  
  return(cell.dist.weight)
}


#############################################################################################
#############################################################################################
#' Calculate coordinates of neighboring cells
#' Given a neighborhood size and starting coordinates, calculate the coordinates of all cells in the neighborhoods around each starting location.
#' @param loc matrix. Two columns representing matrix coordinates of current locations.
#' @param mapdim vector. The number of rows and columns in the map of interest.
#' @param rowcol.delta matrix. The relative coordinates of the neighborhood of interest (identical for all locations).
#' @param n.offset integer. Used to define the neighborhood if needed (i.e., the number of cells on any side of the starting cell to include in the neighborhood)
#' @param torus logical. If TRUE will wrap the map matrix into a torus to remove edge effects.
#' @param na.val numeric. Value to assign locations that fall off of map (if torus == FALSE).
#' @param ... Other parameters to pass through.
#' @author James Forester
#' @modified May 12, 2017 by Lauren White; switched order of row/columns in Cmatrix conversion
#' @export
#' get.neighbors(loc=loc,mapdim=mapdim,rowcol.delta=rowcol.delta,torus=TRUE, na.val=0)

get.neighbors <- function(loc, mapdim, rowcol.delta=NULL, n.offset, torus=TRUE, na.val=NA, ...){
  
  if(is.null(rowcol.delta)){
    rowcol.delta<- expand.grid(-n.offset:n.offset,-n.offset:n.offset)
  }
  
  neighborhood<-nrow(rowcol.delta) #e.g. an offset of 1, yields 9 neighbors
  outmat<-matrix(NA,nrow(loc),neighborhood) #rows= number of inds./col.= size of neighborhood
  
  for(i in 1:neighborhood){
    newloc=cbind(loc[,1]+rowcol.delta[i,1], loc[,2]+rowcol.delta[i,2])
    
    if(torus){
      ##Wrap map to torus
      newloc[newloc[,1]<=0,1]=mapdim[1]+newloc[newloc[,1]<=0,1]
      newloc[newloc[,2]<=0,2]=mapdim[2]+newloc[newloc[,2]<=0,2]
      newloc[newloc[,1]>mapdim[1],1] = newloc[newloc[,1]>mapdim[1],1]-mapdim[1]
      newloc[newloc[,2]>mapdim[2],2] = newloc[newloc[,2]>mapdim[2],2]-mapdim[2]
    }else{
      newloc[newloc[,1]<=0,1]=na.val
      newloc[newloc[,2]<=0,2]=na.val
      newloc[newloc[,1]>mapdim[1],1] = na.val
      newloc[newloc[,2]>mapdim[2],2] = na.val
    }
    
    outmat[,i]=Cmatrix(newloc[,2],newloc[,1],nrows=mapdim[1]) #row, col, nrows
  }
  
  return(outmat)
}

#############################################################################################
#############################################################################################

Cmatrix<-function (row, col, nrows) #y,x
{
  vector_loc <- row + (col - 1) * nrows
  return(vector_loc)
}

#############################################################################################
#############################################################################################
Rmatrix<-function (vector_loc, nrows) 
{
  col <- floor(vector_loc/nrows) + 1
  row <- vector_loc - (col - 1) * nrows + 1
  return(list(col = col, row = row))
}

#############################################################################################
#############################################################################################

#' Function to randomly generate individuals across a landscape, and randomly assign 'nI' number of infected
#' individuals at t=0
#' @param n.initial- the number of desired individuals to start with
#' @param dim- the length (1D) of the landscape
#' @param nI- the number of infected individuals to start with at t=0
#' @author Lauren White
make.inds <- function(n.initial,dim, nI){
  id<-1:n.initial
  xloc<-round(runif(n.initial, min=1, max=dim))
  yloc<-round(runif(n.initial, min=1, max=dim))
  vec<-Cmatrix(yloc,xloc,dim)
  I<-sample(1:n.initial, nI)
  status<-rep("S", times=n.initial)
  status[I]<-"I"
  inds <- data.frame(id = id, xloc=xloc, yloc=yloc, vec=vec, status=status, stringsAsFactors=FALSE) 
  inds
}

#############################################################################################
#############################################################################################
#    inds<-infection3(inds, consdens=as.numeric(Num[[1]]), idens=as.numeric(Num[[2]]), transProb=inf_prob, hab1=hab1)
#' Function that determines whether not a transmission event occurs
#' First checks to see if there are any individuals + infected infividuals present in the same cell
#' Then calculates final transmission probability based on # of infected individuals present
#' @param inds- dataframe with individual traits, including infection status
#' @param consdens- number of conspecifics by cell
#' @param idens- number of infected conspecifics by cell
#' @param transProb- per interaction transmissioin probability
#' @param lxy- convenience data frame for vector to xy notation
#' @param hab1- landscape matrix
#' @author Lauren White
#' @date May 25, 2017
#' 
infection3<-function(inds, consdens, idens, transProb, lxy, hab1){
  i_cells<-which(idens>0 & consdens>1) #Which cells have at least one infected, and at least one other susceptible individual present
  #lxy<-longxy(nrow(hab1))
  for (i in 1:length(i_cells)){
    loc<-lxy[i_cells[i],]
    s_inds<-which(inds$xloc==loc$lx & inds$yloc==loc$ly & inds$status=="S")
    i_count<-length(which(inds$xloc==loc$lx & inds$yloc==loc$ly & inds$status=="I"))
    finalProb <- 1 - (1 - transProb)^i_count
    transmit <- rbinom(length(s_inds), 1, finalProb)
    new_i<-s_inds[which(transmit == 1)]
    inds$status[new_i]<-"I"
  } 
  return(inds)
}

# idens<-Num[[2]]
# consdens<-Num[[1]]
#############################################################################################
#############################################################################################

#' Function that allows infected individuals to recover at rate 'gamma'
#' @param inds: the dataframe describing the current status and location of individuals in the simulation
#' @param gamma: the recovery rate
#' @author Lauren White
recover.inds<-function(inds, gamma){
  infected<-which(inds$status=="I" ) #which individuals are currently infected?
  if(length(infected>0)){
    rec.prob<-runif(length(infected), min=0, max=1)
    for (i in 1:length(infected)){
      if(rec.prob[i]<= gamma){
        inds$status[infected[i]]<-"R"
      }
    }
  }
  return(inds)
}

#############################################################################################
#############################################################################################

#'Convenience function to calculate xy coordinates that correspond to numeric length of matrix
#'assigns numbers vertically per default settings in R
#'column-x 
#'row-y
#'@param dim: dimension (1D) of landscape
longxy<-function(dim){
  ly <- rep(1:dim, dim)
  lx <- rep(1:dim, each=dim)
  lxy <- data.frame(lx,ly)
  lxy
}

repmat<-function (x, m = 1, n = 1) 
{
  if (length(x) * n * m > 1e+07) {
    return(kronecker(matrix(1, n, m), x))
  }
  else {
    mx = dim(x)[1]
    nx = dim(x)[2]
    return(matrix(t(matrix(x, mx, nx * n)), mx * m, nx * 
                    n, byrow = T))
  }
}

#############################################################################################
#############################################################################################

#'Modified/simplified sample.cell function
#' newloc<-sample.mod(loc=newloc, map.mat=linear.pred, n.offset=neighborhood.size)
sample.mod<-function(loc,map.mat,n.offset){
  # loc<-newloc
  # map.mat<-linear.pred
  # n.offset<-neighborhood.size
  
  rowcol.delta <- expand.grid(-n.offset:n.offset,-n.offset:n.offset)
  
  ##compute the distance weight to each neighboring cell
  cell.dist.weight <- get.dist.weight(rowcol.delta, dist.weight.fun=NULL)
  
  ##calculate the map locations of the neighborhood surrounding each current animal location
  mapdim<-dim(map.mat)
  neighbors <- get.neighbors(loc=loc,mapdim=mapdim,rowcol.delta=rowcol.delta,torus=TRUE, na.val=0)
  
  ##extract the map values for each neighboring cell
  neighbor.vals=list()
  neighbor.vals<-matrix(map.mat[as.vector(neighbors)],nrow=nrow(neighbors), ncol=ncol(neighbors))
  
  ##calculate the selection function based on map values within the neighborhood
  rsf.vals <-exp(neighbor.vals)
  
  ##compute the probability of moving to any cell within the available neighborhood (scaled by competing values)
  pmat<-rsf.vals *repmat(matrix(cell.dist.weight,1,length(cell.dist.weight)),nrow(rsf.vals),1)
  pmat<-pmat/rowSums(pmat, na.rm=TRUE)
  pmat[is.na(pmat)]=0
  
  ##sample based on pmat probability
  which.cell<-apply(pmat,1,sampfun) #sample(1:length(x),1,prob=x) 
  newloc<-cbind(loc[,1]+rowcol.delta[which.cell,1], loc[,2]+rowcol.delta[which.cell,2])
  
  ## Alternative way to get new locations using previously calculated neighborhoods
  # newloc.vec<-rep(NA, times=nrow(neighbors))
  # for(i in 1:nrow(neighbors))
  # {
  #   newloc.vec[i]<-neighbors[i, which.cell[i]]
  # }
  # test<-Rmatrix(newloc.vec-1, mapdim[1])
  
  
  ## Record new location, wrapping on torus if needed
  ##Wrap map to torus
  if(any(newloc<=0)){
    newloc[newloc[,1]<=0,1]=mapdim[1]+newloc[newloc[,1]<=0,1]
    newloc[newloc[,2]<=0,2]=mapdim[2]+newloc[newloc[,2]<=0,2]
  }
  if(any(newloc>mapdim[1])){ #assuming symmetrical landscape
    newloc[newloc[,1]>mapdim[1],1] = newloc[newloc[,1]>mapdim[1],1]-mapdim[1]
    newloc[newloc[,2]>mapdim[2],2] = newloc[newloc[,2]>mapdim[2],2]-mapdim[2]
  }
  return(newloc)
}
#############################################################################################
#############################################################################################

calc.dens<-function(hab1, newloc.vec, n.initial, inds){
  consdens<-hab1*0
  idens<-hab1*0
  for(i in 1:length(newloc.vec)) {
    consdens[newloc.vec[i]]=consdens[newloc.vec[i]]+1
    if(inds$status[i]=="I"){
      idens[newloc.vec[i]]=idens[newloc.vec[i]]+1
    }
  }
  dens<-consdens/max(consdens) #normalize to 1 by dividiing conspecific density by maximum number of individuals found in a single cell in the landscape
  Num<-list(consdens,idens,dens) #returns (as list): number of conspecifcs per cell, number of infected conspecifcs per cell, and normalized density
  return(Num)
}

#############################################################################################
#############################################################################################

#' Create neutral landscape maps
#' 
#' Use standard methods to generate fractal maps. Binary and continuous surfaces may be produced.
#' 
#' @param k integer. The extent of the map (2^k+1)^2 pixels
#' @param h numeric. Level of aggregation in the map.
#' @param p numeric (0,1). The proportion of map in habitat=1
#' @param binary logical. If TRUE, a 0/1 categorical landscape is produced.
#' @author Shannon Pittman, James Forester, modified by Lauren White
#' @export
#' @example examples/neutral.landscape_example.R
fracland_mod <- function(k, h, p, binary = TRUE) {
  ## Function for creating neutral landscapes Shannon Pittman University of Minnesota May, 2013 k = the extent of the map (2^k+1)^2 pixels h =
  ## how clumped the map should be (ranging from ?? to ??) -- weird behavior at higher values p = proportion of map in habitat 1 binary =
  ## plotflag == if TRUE will plot a filled contour version of the matrix
  
  ## function call: testmap=land(6,1,.5,FALSE,TRUE)
  A <- 2^k + 1  # Scalar-determines length of landscape matrix
  
  #Right now, as written (1-p) represents the amount of habitat listed as "1"
  B <- matrix(0, A, A)  # Creates landscape matrix
  
  B[1, 1] <- 0
  B[1, A] <- 0
  B[A, 1] <- 0
  B[A, A] <- 0
  
  
  iter <- 1
  for (iter in 1:k) {
    scalef <- (0.5 + (1 - h)/2)^(iter)
    
    d <- 2^(k - iter)
    
    # ALL SQUARE STEPS#
    for (i in seq(d + 1, A - d, 2 * d)) {
      for (j in seq(d + 1, A - d, 2 * d)) {
        B[i, j] <- mean(c(B[i - d, j - d], B[i - d, j + d], B[i + d, j - d], B[i + d, j + d])) + scalef * rnorm(n = 1)
      }
    }
    
    # OUTSIDE DIAMOND STEP#
    for (j in seq(d + 1, A - d, 2 * d)) {
      B[1, j] <- mean(c(B[1, j - d], B[1, j + d], B[1 + d, j])) + scalef * rnorm(n = 1)
      B[A, j] <- mean(c(B[A, j - d], B[A, j + d], B[A - d, j])) + scalef * rnorm(n = 1)
    }
    
    for (i in seq(d + 1, A - d, 2 * d)) {
      B[i, 1] <- mean(c(B[i - d, 1], B[i + d, 1], B[i, 1 + d])) + scalef * rnorm(n = 1)
      B[i, A] <- mean(c(B[i - d, A], B[i + d, A], B[i, A - d])) + scalef * rnorm(n = 1)
    }
    
    # INSIDE DIAMOND STEP#
    if (2 * d + 1 <= A - 2 * d) {
      for (i in seq(d + 1, A - d, 2 * d)) {
        for (j in seq(2 * d + 1, A - 2 * d, 2 * d)) {
          B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
        }
      }
      
      for (i in seq(2 * d + 1, A - 2 * d, 2 * d)) {
        for (j in seq(d + 1, A - d, 2 * d)) {
          B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
        }
      }
    }
    
    iter <- iter + 1
  }
  
  if (binary == T) {
    R <- sort(B)
    PosR <- (1 - p) * length(R)  #larger values become habitat, designated as 1
    pval <- R[PosR]
    T1 <- which(B > pval)
    T2 <- which(B <= pval)
    B[T1] <- 1  #habitat is 1
    B[T2] <- 0
  } 
  return(B)
}



##################################################

#'Implements above functions in a loop that can be run in parallel
#'@param maxtime: maximum time duration of simulation
#'@param nsim: number of simulations to be repeated per parameter set
#'@param density: desired conspecific density for simulation
#'@param inf_prob: infection probability
#'@param rec_rate: recovery rate
#'@param p: proportion of available habitat
#'@param H: Hurst's exponent - degree of clustering
#'@param beta1- RSF covariate selecting for resource availablity
#'@param beta2- RSF covariate selecting for conspecific density
#'@param beta3- RSF covariate selecting for conspecific density^2
#'@param percep-perceptual distance (# of cells beyond current cell)
IBM_loop_function<-function(maxtime, nsim, k, density, inf_prob, rec_rate, p, H, beta1, beta2, beta3, percep){
  #Set up initial conditions for all simulations
  dim<-2^k + 1
  n.initial <- round(density*dim^2) # inital population size
  neighborhood.size<-percep
  
  ##Initiate data structure across simulations
  mydata<-list(inds=NULL, landscape=NULL, N=NULL, movedat=NULL, infdat=NULL)
  filename<-paste("k", k, "dens", density, "inf", inf_prob, "rec", rec_rate,"p",p,"H",H, "betas", beta1, beta2, beta3, "per", percep, sep="")
  summary<-data.frame(k=rep(k, times=nsim), density=rep(density, times=nsim), inf_prob=rep(inf_prob, times=nsim), rec_rate=rep(rec_rate, times=nsim), p=rep(p, times=nsim), H=rep(H, times=nsim), beta1=rep(beta1, times=nsim), beta2= rep(beta2, times=nsim), beta3=rep(beta3, times=nsim), percep=rep(percep, times=nsim), duration=NaN, max_I=NaN, max_prevalence=NaN)
  infected<-matrix(NaN, nrow=maxtime, ncol=nsim)
  landscape0<-matrix(NaN, nrow=nsim, ncol=38)
  landscape1<-matrix(NaN, nrow=nsim, ncol=38)
 
  #Iterate for 'nsim' number of simulations per given parameter set
  for(count in 1:nsim){
    ##Generate landscape with specified p & H
    hab1<-fracland_mod(k=k, h=H, p=p, binary=TRUE)
    landscape_vals<-ClassStat(hab1) #SDMTools package
    landscape0[count,]<-as.numeric(landscape_vals[1,]) #patch_value=0
    landscape1[count,]<-as.numeric(landscape_vals[2,]) #patch_value=1
    lxy<-longxy(nrow(hab1))
    ##Initiate individuals
    inds <- make.inds(n.initial=n.initial, dim=nrow(hab1), nI=1)
    mydata$inds[[count]]<-inds
    
    ##Record initial locations
    startloc<-inds$vec
    newloc.vec<-startloc
    
    ##Convert vector notation to matrix coordinates
    newloc<-as.data.frame(Rmatrix(newloc.vec-1,nrow(hab1)))
    
    #Initiate longterm storage for movement and infection data
    movedat<-startloc
    infdat<-inds$status
    
    # Object for storing results for single simulation (disease dynamics through time)
    N <-data.frame(S=NaN, I=NaN, R=NaN)
    N[1,] <- c(sum(inds$status=="S"), sum(inds$status=="I"), sum(inds$status=="R"))
    
    ##Run the simulation for a max of 'T' timesteps
    for(t in 2:maxtime){
      ##Calculate conspecific density and infected conspecific density (individuals/cell)
      Num<-calc.dens(hab1=hab1, newloc.vec=newloc.vec, n.initial=n.initial, inds=inds)
      
      ##Calculate updated rsf map
      linear.pred<-(beta1*hab1+beta2*as.numeric(Num[[1]])+ beta3*as.numeric(Num[[1]])^2)
      
      ##Calculate the new locations based on the neighborhood around each of previous locations
      newloc<-sample.mod(loc=newloc, map.mat=linear.pred, n.offset=neighborhood.size)
      
      ##Update location lists
      inds$xloc<-newloc[,1] #col
      inds$yloc<-newloc[,2] #row
      newloc.vec<-Cmatrix(newloc[,2],newloc[,1],nrow(hab1)) #row, col, nrow
      movedat<-cbind(movedat,newloc.vec)
      inds$vec<-newloc.vec
      
      ##Recalculate conspecfic (and infected conspecific) density (individuals/cell)
      Num<-calc.dens(hab1=hab1, newloc.vec=newloc.vec, n.initial=n.initial, inds=inds)
      
      #Disease process(es)
      inds<-infection3(inds, consdens=as.numeric(Num[[1]]), idens=as.numeric(Num[[2]]), transProb=inf_prob, lxy=lxy, hab1=hab1)
      inds<-recover.inds(inds, gamma=rec_rate)
      
      #Update disease status lists
      infdat<-cbind(infdat,inds$status)
      N[t,] <- c(sum(inds$status=="S"), sum(inds$status=="I"), sum(inds$status=="R"))
      if(sum(inds$status=="I")==0){ #If number of infected in individuals--> 0, stop running
        summary$duration[count]<-t-1
        for (i in t:maxtime){
          N[i,]<-N[t,] #automatically fill remaining time slots with current N values
        }
        break
      }
    }
    
    mydata$N[[count]]<-N
    mydata$landscape[[count]]<-hab1
    mydata$movedat[[count]]<-movedat
    mydata$infdat[[count]]<-infdat
    summary$max_I[count]<-max(N$I)
    summary$max_prevalence[count]<-max(N$I)/n.initial
    infected[,count]<-N$I
    print(Num[[1]])
  }
  #Record output
  write.csv(summary,file=paste(filename,"summary.csv", sep="_"))
  write.csv(infected,file=paste(filename,"infected.csv", sep="_"))
  
  colnames(landscape0)<-colnames(landscape_vals)
  colnames(landscape1)<-colnames(landscape_vals)
  landscape0<-as.data.frame(landscape0)
  landscape1<-as.data.frame(landscape1)
  
  write.csv(landscape0,file=paste(filename,"landscape0.csv", sep="_"))#, col.names=colnames(landscape_vals))
  write.csv(landscape1,file=paste(filename,"landscape1.csv", sep="_"))#, col.names=colnames(landscape_vals))
  print(filename)
}

#c(500, 100, 8, 0.5, 0.2, 0.2, 0.5, 0.5, 0, -6, 3)
#c(500, 100, 8, 0.5, 0.2, 0.2, 0.25, 0.1, 3, -6, 3)
#(500, 100, 8, 0.5, 0.2, 0.2, 0.25, 0.5, 6, -6, 3)
#k8dens0.25inf0.2rec0.1p0.25H0.1betas60per3
library(SDMTools)
tic=Sys.time()
# 20A
# "k8dens0.5inf0.2rec0.1p0.25H0.1betas6-3per3_summary.csv"
IBM_loop_function(maxtime=10000, nsim=1, k=5, density=0.5, inf_prob=0.2, rec_rate=0.1, p=0.5, H=0.5, beta1=0, beta2=2, beta3=-0.5, percep=3)
print(difftime(Sys.time(),tic,units="mins"))

# tic=Sys.time()
# #11A
# IBM_loop_function(100, 100, 8, 0.25, 0.2, 0.1, 0.25, 0.1, 6, -3, 3)
# print(difftime(Sys.time(),tic,units="mins"))

# # 
# ##Make plots of some individual trajectories; note that most are within small area.
# hab1<-test$landscape[[1]]
# movedat<-test$movedat[[1]]
# image(hab1,asp=1)
# for(i in 1:10){
#   testdat=as.data.frame(Rmatrix(movedat[i,]-1,nrow(hab1)))
#   plot(testdat,type="b",asp=1)#, xlim=c(1,dim(hab1)[1]), ylim=c(1,dim(hab1)[2]))
#   points(testdat[1,],pch=10,col="green",cex=2,lwd=2)
#   points(testdat[nrow(testdat),],pch=17,col="red",cex=2,lwd=2)
# }
# 
# ##Make plots of infection curves
# plot(x=NULL, y=NULL,xlab="Time", ylab="NUmber of Infected Inidivduals", xlim=c(0,100), ylim=c(0,10))
# for(i in 1:length(test$N))
# {
#   plot(1:100, test$N[[i]]$I, type="l", add=TRUE)
# }

# test<-fracland_mod(9, 0.5, 0.5)
# library(SDMTools)
# landscape_vals<-ClassStat(test)
# landscape_vals[1,] #patch_value=0
# landscape_vals[2,] #patch_value=1
