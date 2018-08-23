# Purpose        : Fit/predict distribution of soil types (via multinom);
# Maintainer     : Tomislav Hengl (tom.hengl@isric.org)
# Contributions  : Bas Kempen (bas.kempen@wur.nl); Dainius Masiliūnas (dainius.masiliunas@wur.nl)
# Dev Status     : Pre-Alpha
# Note           : if the regression model is difficult to fit, it might lead to artifacts (see also "mlogit" package);


## fit a multinomial logistic regression and make predictions:
setMethod("spmultinom", signature(formulaString = "formula"), function(formulaString, observations, covariates, class.stats = TRUE, predict.probs = TRUE, ...){

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
  tv <- all.vars(formulaString)[1]
  sel <- names(cov.df) %in% all.vars(formulaString)[-1]
  if(all(sel==FALSE)|length(sel)==0){
      stop("None of the covariates in the 'formulaString' matches the column names in the 'covariates' object")
  }

  ## over observations and covariates:
  if (class(observations) == "SpatialPointsDataFrame")
  {
    ov <- over(observations, covariates[sel])
    ov <- cbind(obs.df[tv], ov)
  } else {
    ov <- cbind(obs.df, cov.df[complete.cases(cov.df),])
  }
    
  message("Fitting a multinomial logistic regression model...")
  mout <- nnet::multinom(formulaString, ov, ...)
  cout <- as.factor(paste(predict(mout, newdata=cov.df, na.action = na.pass)))

  ## predict probabilities if required:
  if(predict.probs == TRUE){
     if (class(covariates) != "SpatialPixelsDataFrame")
        stop("predict.probs not implemented yet for non-SpatialPixelsDataFrame covariates")
     probs <- predict(mout, newdata=covariates, type="probs", na.action = na.pass)
     mm <- covariates[1]
     mm@data <- data.frame(probs)
     pm <- covariates[1]
     pm@data[,tv] <- cout
     pm@data[,names(covariates)[1]] <- NULL    
     
     ## kappa statistics:
     if(requireNamespace("mda", quietly = TRUE)&requireNamespace("psych", quietly = TRUE)){
       cout.m <- as.factor(paste(predict(mout, newdata=ov, na.action = na.pass)))
       cf <- mda::confusion(cout.m, as.character(ov[,tv]))
       ## remove missing classes:
       a = attr(cf, "dimnames")[[1]] %in% attr(cf, "dimnames")[[2]] 
       b = attr(cf, "dimnames")[[2]] %in% attr(cf, "dimnames")[[1]]
       c.kappa = psych::cohen.kappa(cf[a,b])
       ac <- sum(diag(cf))/sum(cf)*100
       message(paste("Estimated Cohen Kappa (weighted):", signif(c.kappa$weighted.kappa, 4)))
       message(paste("Map purity:", signif(ac, 3)))
     } else {
       stop("Packages 'mda' and 'psych' required but not available")
     }
  }
  
  ## remove object class for consistency:
  class(mout) = "list"

  # estimate class centres using the results of multinom:
  if(class.stats == TRUE){
    ca <- aggregate(cov.df[,sel], by=list(cout), FUN="mean")
    class.c <- as.matrix(ca[-1]); attr(class.c, "dimnames")[[1]] <- ca[,1]
    ca <- aggregate(cov.df[,sel], by=list(cout), FUN="sd")
    class.sd <- as.matrix(ca[-1]); attr(class.sd, "dimnames")[[1]] <- ca[,1]
    # mask out classes that result in NA:
    for(c in row.names(class.c)){ 
        if(any(is.na(class.c[c,]))|any(is.na(class.sd[c,]))|any(class.sd[c,]==0)){ 
            class.c <- class.c[!(row.names(class.c) %in% c),]; class.sd <- class.sd[!(row.names(class.sd) %in% c),]
        } 
    }
    } else {
      class.c = NULL; class.sd = NULL
    }
    
  # create the output object:
  if(predict.probs == TRUE){
    if(all(!is.na(rowSums(mm@data, na.rm=TRUE)))){
      out <- new("SpatialMemberships", predicted = pm, model = mout, mu = mm, class.c = class.c, class.sd = class.sd, confusion = cf)
    } else {
      stop("Predicted probabilities contain missing values. Consider removing some classes or setting 'predict.probs == FALSE'") 
    }
  } else {  
    out <- list(model=mout, fit=cout, class.c=class.c, class.sd=class.sd)
  }
    
  return(out)
})

# end of script;
