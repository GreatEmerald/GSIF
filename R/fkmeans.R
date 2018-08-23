# Purpose        : Fit/predict distribution of soil types (memberships);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Bas Kempen (bas.kempen@wur.nl); Dainius Masiliūnas (dainius.masiliunas@wur.nl)
# Dev Status     : Pre-Alpha
# Note           : if the regression model is difficult to fit, it might lead to artifacts;


# Fit a supervised fuzzy kmeans model and predict memberships:
setMethod("spfkm", signature(formulaString = "formula"), function(formulaString, observations, covariates, class.c = NULL, class.sd = NULL, fuzzy.e = 1.2){
  
  ## generate formula if missing:
  if(missing(formulaString)) {  
    formulaString <- as.formula(paste(names(observations)[1], "~", paste(names(covariates), collapse="+"), sep=""))
  }
  ## check the formula string:
  if(!plyr::is.formula(formulaString)){
      stop("'formulaString' object of class 'formula' required")
  }
  
  ## get regular data.frames from the input
  if (class(observations) == "SpatialPointsDataFrame")
  {
    obs.df = observations@data
  } else if (is.data.frame(observations)) {
    obs.df = observations
  } else {
    stop("'observations' must be a SpatialPointsDataFrame or a regular data.frame")
  }
  
  if (class(covariates) == "SpatialPixelsDataFrame")
  {
    cov.df = covariates@data
  } else if (is.data.frame(covariates)) {
    cov.df = covariates
  } else {
    stop("'covariates' must be a SpatialPixelsDataFrame or a regular data.frame")
  }
  
  ## selected variables:
  tv = all.vars(formulaString)[1]
  sel = names(cov.df) %in% all.vars(formulaString)[-1]
  if(all(sel==FALSE)|length(sel)==0){
      stop("None of the covariates in the 'formulaString' matches the column names in the 'covariates' object")
  }
 
  ## if available, use class centres:
  check_tc <- !is.null(class.c)&!is.null(class.sd)
  if(check_tc){
    if(!class(class.c)=="matrix"){ stop("Object of type 'matrix' with column names for covariates and row names correspodning to the class names required") }
    if(!class(class.sd)=="matrix"){ stop("Object of type 'matrix' with column names for covariates and row names correspodning to the class names required") }
    mout = list(NULL)
  }
  ## otherwise, estimate class centres using the multinomial logistic regression:
  else {
    message("Trying to estimate the class centres using the 'multinom' method...")
    ## multinomial logistic regression:
    rout <- spmultinom(formulaString=formulaString, observations, covariates, class.stats=TRUE, predict.probs=FALSE)
    mout = rout$model
    if(length(unique(rout$fit))<2){ stop("Predictions resulted in <2 classes. See ?multinom for more info") }
    class.c = rout$class.c
    class.sd = rout$class.sd
  }
  
  cl <- as.list(row.names(class.c))
  dsf <- NULL
  ## derive distances in feature space:
  for(c in unlist(cl)){
      dsf[[c]] <- data.frame(lapply(names(cov.df)[sel], FUN=function(x){rep(NA, length(cov.df[,1]))}))
      names(dsf[[c]]) <- names(cov.df)[sel]
      for(j in names(cov.df)[sel]){
         dsf[[c]][,j] <- ((cov.df[,j]-class.c[c,j])/class.sd[c,j])^2
      }
  }
  ## sum up distances per class:
  ds <- NULL
  ds <- lapply(dsf, FUN=function(x){sqrt(rowSums(x, na.rm=TRUE, dims=1))})
  names(ds) <- unlist(cl)
  ds <- data.frame(ds)
  ## total sum:
  tt <- rowSums(ds^(-2/(fuzzy.e-1)), na.rm=TRUE, dims=1)
  ## derive the fuzzy membership:
  mm <- cov.df[1]
  for(c in unlist(cl)){
    mm[,c] <- (ds[,c]^(-2/(fuzzy.e-1))/tt)
  }
  mm[,names(cov.df)[1]] <- NULL
  
  ## Derive the dominant class:
  maxm <- sapply(data.frame(t(as.matrix(mm))), FUN=function(x){max(x, na.rm=TRUE)})
  ## class having the highest membership
  cout <- NULL
  for(c in unlist(cl)){
       cout[which(mm[,c] == maxm)] <- c
  }
  cout <- as.factor(cout)
  
  ## construct a map: overlay observations and covariates:
  if (class(observations) == "SpatialPointsDataFrame")
  {
    pm <- covariates[1]
    pm@data[,tv] <- cout
    pm@data[,names(covariates)[1]] <- NULL
    
    ov <- over(observations, pm)
  } else {
    pm <- data.frame(cout)
    names(pm) <- tv
    ov <- cbind(obs.df, cov.df[complete.cases(cov.df),])
  }
  sel.c <- !is.na(ov[,tv]) & !is.na(obs.df[,tv])

  ## kappa statistics:
  if(requireNamespace("mda", quietly = TRUE)&requireNamespace("psych", quietly = TRUE)){
    cf <- mda::confusion(ov[sel.c,tv], as.character(obs.df[sel.c,tv]))
    ## remove missing classes:
    a <- attr(cf, "dimnames")[[1]] %in% attr(cf, "dimnames")[[2]] 
    b <- attr(cf, "dimnames")[[2]] %in% attr(cf, "dimnames")[[1]]
    c.kappa = psych::cohen.kappa(cf[a,b])
    message(paste("Estimated Cohen Kappa (weighted):", signif(c.kappa$weighted.kappa, 4)))  
  } else {
    cf <- NULL
  }
  
  ## create the output object:
  out <- new("SpatialMemberships", predicted = pm, model = mout, mu = mm, class.c = class.c, class.sd = class.sd, confusion = cf)
  return(out)

})

# end of script;
