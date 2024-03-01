
file <- list.files("./data", full.names = TRUE)
task <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#file <- files[task]
dt <- readRDS(file)
dt <- dt[,5:18]

# Assuming adl.05 is your data frame
dt <- as.data.frame(lapply(dt, function(x) {
  if(is.numeric(x)) {
    # Apply the transformation logic to numeric columns
    ifelse(x == 1, 2, ifelse(x == 0, 1, x))
  } else {
    # Return non-numeric columns unchanged
    x
  }
}))


loglik_lca <- function(param, data, n.class){
  
  if(length(param[param<0 | param> 1]) >0 | (1-sum(param[1:n.class-1]))<=0){return(9999)}else{
    
    
  
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
  

  return(-loglik)}
}

loglik_lca(startvalues,data = dt, n.class = 2)

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

res <- optim(par = startvalues,method = "SANN", fn = loglik_lca, data = dt, n.class = nc)


all.sum <- list()
#params <- c()

#for (i in 1:length(res$par[nc:length(res$par)])) {
#  params[2 * i - 1] <- res$par[nc:length(res$par)][i]
#  if (i <= length(res$par[nc:length(res$par)])) {
#    params[2 * i] <- 1 - res$par[nc:length(res$par)][i]
#  }
#}
#params <- array(params, dim=c(max(dt),ncol(dt),nc))
# ordering classes based on class prevalence  (combat label switching)
#dimnames(params)[[3]] <- matrix(c(res$par[1:nc-1],1-sum(res$par[1:nc-1])),ncol = 1)
#ordered_params <- params[, , order(as.numeric(dimnames(params)[[3]]),decreasing = T)]

all.sum <- list(startvalues,res$value,res$convergence,res$par,res$counts,res$message)

path <- Sys.getenv("BOOT_OUTPUT_DIR")
loc <- paste0(path,"/test",task,"SA",".Rds")
saveRDS(all.sum,file = loc)

