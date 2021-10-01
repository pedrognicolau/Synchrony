# Important functions ####

## synchrony function smooth ####
gps <- read.csv2("Data/gps.csv")

sync_df5 <- function(res, round=0, dist_mat=gps$cumdist, method="pearson")
  # residuals as a smooth function of distance
  {
  distmat <- stats :: dist(dist_mat, upper = TRUE, diag = TRUE)
  m <- data.frame(t(combn(as.numeric(rownames(gps)),2)), dist=as.numeric(distmat))
  
  CR2 <- cor(res, method=method)
  CR2[upper.tri(CR2)] <- NA
  diag(CR2) <- NA
  nCR <- reshape2::melt(CR2, varnames = c('X2', 'X1'), na.rm = TRUE)
  colnames(nCR)[3] <- "corres"
  
  mat2 <- dplyr::left_join(m,nCR)
  mat3 <- arrange(mat2, dist)
  mat3$dist <- round(mat3$dist, round)
  mat4 <- aggregate(corres~dist, data=mat3, FUN=mean)
  
  return(mat4)
}



## functions for R squared ####

fitted_samples <- function(inlafit, ns=200)
{
  Lp <- nrow(inlafit$summary.fixed) # length of predictor
  
  # number of elements of the predictor (N size)
  npred <- nrow(inlafit$summary.linear.predictor)
  
  # sample from the posterior
  postsample <- inla.posterior.sample(n=ns, result=inlafit)
  
  # matrix npred x nsamples
  samp_pred <- matrix(data=NA,nrow =npred,ncol=ns )
  
  # loop over postsample
  for(i in 1:ns)
  {
    samp_pred[,i] <- postsample[[i]]$latent[1:npred]
  }
  
  return(samp_pred)
}
#inlafit <- spring_Ab
#y_obs <- stackeddata$S_t.x
residual_samples <- function(y_obs, inlafit, ns=200, mean = FALSE)
  # y_obs is the observed Y
  # y_pred is the marginal distribution matrix given by fitted_margins
{
  samppred <- fitted_samples(inlafit, ns=ns)
  residual_samples <- samppred-y_obs
  if(mean==TRUE) residual_samples <- apply(residual_samples,1,mean)
  return(residual_samples)
}



inla_R2 <-function(y_obs, inlafit, ns=200)
{
  
  y_pred <- fitted_samples(inlafit,ns=ns)
  y_res <- residual_samples(y_obs,inlafit,ns=ns)
  var_ypred <- apply(y_pred,2,var)
  var_res <- apply(y_res,2,var)
  
  r_sample <- var_ypred/(var_ypred+var_res)
  summary(r_sample)
  return(r_sample)
  
} 


# m <- data.frame(t(combn(as.numeric(rownames(gps)),2)), dist=as.numeric(distmat))

rescor_distance <- function(res, round=0, dist_mat=gps$cumdist, method="spearman", nstations=19){
  res <- matrix(res, ncol=nstations)
  
  distmat <- stats :: dist(dist_mat, upper = TRUE, diag = TRUE)
  m <- data.frame(t(combn(as.numeric(rownames(gps)),2)), dist=as.numeric(distmat))
  
  CR2 <- cor(res, method=method)
  CR2[upper.tri(CR2)] <- NA
  diag(CR2) <- NA
  nCR <- reshape2::melt(CR2, varnames = c('X2', 'X1'), na.rm = TRUE)
  colnames(nCR)[3] <- "corres"
  
  mat2 <- dplyr::left_join(m,nCR)
  mat3 <- dplyr::arrange(mat2, dist)
  mat3$dist <- round(mat3$dist, round)
  mat4 <- aggregate(corres~dist, data=mat3, FUN=mean)
  
  mat5 <- data.frame(res=mat4[,2])
  rownames(mat5) <- round(mat4[,1],3)
  
  return(mat5)
}


residual_correlation <- function(inlamodel, y_obs,round=0, FUNCTION=mean)
{
  res1 <- residual_samples(y_obs,inlamodel)
  
  matest <- apply(res1,2,FUN=rescor_distance, round=round)
  xa <- bind_cols(matest)
  xa2 <- as.data.frame(apply(xa,1,FUNCTION))
  colnames(xa2) <- "mean.res"
  xa2$distance <- as.numeric(rownames(xa2))
  # plot(xa2$distance,xa2$mean.res)
  # lines(smooth.spline(xa2$distance,xa2$mean.res, spar=1.2))
  
  return(xa2)
  
}