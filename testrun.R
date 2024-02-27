
file <- list.files("./data", full.names = TRUE)
#task <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#file <- files[task]
dt <- readRDS(file)
dt <- dt[,5:18]

loglik_lca <- function(param, data, n.class){
  gammahat <- matrix(c(param[1:n.class-1],1-sum(param[1:n.class-1])),ncol = 1)
  
  oldrho <- param[n.class:length(param)]
  rhohat <- numeric(length(oldrho) * 2)
  
  # Insert values in between each existing value
  for (i in 1:length(oldrho)) {
    rhohat[2 * i - 1] <- oldrho[i]
    if (i <= length(oldrho)) {
      rhohat[2 * i] <- 1 - oldrho[i]
    }
  }
  
  rhohat = array(rhohat, dim=c(max(data),ncol(data),n.class))
  
  loglik = 0
  for(i in 1:nrow(data)){
    
    loglik = loglik + log(sum(gammahat*apply(rhohat, 3, function(x) prod(diag(x[unlist(data[i,]),])))))
  }
  

  return(-loglik)
}

# change later 
#nc <- 4
#ni <- 5
# n.param = ni*nc+(nc-1) = 
n.param <- 29
nc <- 2

startvalues <- runif(n.param)
while(sum(startvalues[1:nc-1])>=1){
  startvalues <- runif(n.param)
}

res <- optim(par = startvalues,method = "L-BFGS-B", fn = loglik_lca, data = dt, n.class = nc,lower = rep(0.00001,n.param), upper = rep(.9999,n.param))


all.sum <- list()
params <- c()

for (i in 1:length(res$par[nc:length(res$par)])) {
  params[2 * i - 1] <- res$par[nc:length(res$par)][i]
  if (i <= length(res$par[nc:length(res$par)])) {
    params[2 * i] <- 1 - res$par[nc:length(res$par)][i]
  }
}
params <- array(params, dim=c(max(dt),ncol(dt),nc))
# ordering classes based on class prevalence  (combat label switching)
dimnames(params)[[3]] <- matrix(c(res$par[1:nc-1],1-sum(res$par[1:nc-1])),ncol = 1)
ordered_params <- params[, , order(as.numeric(dimnames(params)[[3]]),decreasing = T)]

all.sum <- list(startvalues,ordered_params,res$value,res$convergence)

all.sum
