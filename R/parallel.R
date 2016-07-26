#This function contains a wrapper for parallelizing code
#inspired by the boot package (see boot package and parallel
#package vignette for further details)
#Also the validation routines of the parallel arguments

#Run functions in parallel
#@inheritparam gammaImpute
#@param lapply.list the list or vector to be iterated over when using parallel apply
#@param FUN the function to be called inside parallel apply
#@param ... Additional arguments to be passed to FUN
#@return list of answers
#See boot::boot and parallel package vignette for further details
parallelRun <- function(parallel,ncpus,cl,lapply.list,FUN,...){
  if(!require("parallel")){
    stop("parallel package is not available")
  }
  
  #code from boot package
  have_mc <- have_snow <- FALSE
  if (parallel == "multicore") 
    have_mc <- .Platform$OS.type != "windows"
  else if (parallel == "snow") 
    have_snow <- TRUE
  
  if(!have_mc && !have_snow) 
    stop("Invalid parallel option")
  
  loadNamespace("parallel")
  
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  
  if(RNGkind()[1] != "L'Ecuyer-CMRG"){
    warning("The L'Ecuyer-CMRG random number generator has not been used so",
            " reproducibility cannot be guaranteed. Use the command: RNGkind(\"L'Ecuyer-CMRG\")",
            " for reproducibility and check the parallel package documentation for further details.")
  }
  
  if (have_mc) {
    parallel::mclapply(lapply.list, FUN, ..., 
                       mc.preschedule = TRUE, mc.set.seed =  TRUE, mc.cores = ncpus)
  }
  else if (have_snow) {
    list(...) #force evaluation of args
    if (is.null(cl)) {
      cl <- parallel::makePSOCKcluster(rep("localhost",ncpus))
      if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
        parallel::clusterSetRNGStream(cl)
      #This cluster call is needed as snow doesn't seem
      #to pass libraries which a package depends on onto
      #the cluster
      clusterCall(cl, function() library("survival"))
      res <- parallel::parLapply(cl, lapply.list, FUN, ...)
      parallel::stopCluster(cl)
      res
    }
    else{
      parallel::parLapply(cl, lapply.list, FUN)  
    } 
  }
  
}


#validate function arguments for parallelization
validate.parallel.arguments <- function(parallel, ncpus, cl){
  if(!is.null(cl) && !"cluster" %in% class(cl)){
    stop("cl is not a cluster")
  }
  
  if(!parallel %in% c("no", "multicore", "snow")){
    stop("Invalid argument for parallel")
  }
  
  if(!.internal.is.finite.number(ncpus) ||!.internal.is.wholenumber(ncpus) || ncpus < 1){
    stop("ncpus must be a positive integer")
  }
  
  if(parallel=="no" &&  ncpus != 1) stop("Cannot have ncpus > 1 if parallel == no")
}
