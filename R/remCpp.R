#' remCpp  
#'
#' A function that returns maximum likelihood estimates of REM (interval timing only) by using the optim function
#'
#' @param stats an array of dimensions n_dyads*statistics*M (if stats is an array for relevent::rem then use aperm(a = stats, perm = c(2,3,1)) )
#' @param event_binary a matrix of dimensions M*n_dyads (for each time point it is a vector with 1 when the dyad is observed 0 for the others)
#' @param interevent_time a vector of interevent times 
#' @param threads number of cores to use in the parallelization (i would change this parameter in 'n_threads')
#'
#' @return  object list (the function saves also the output of optim)
#' @export
remCpp <- function(pars_init = NULL,
                    stats, 
                    event_binary, 
                    interevent_time, 
                    threads, 
                    method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")){

    # some check before running the optimization
    if(length(method)>1){stop("only ONE method can be specified")}
    else{if(!method%in%c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                 "Brent")){stop("the method defined is not one of those supported by optim(). Read optim() documentation (?optim)")}
                 }
    if(is.null(pars_init)) pars_init <- rep(0,dim(stats)[2])
    else{ if(length(pars_init) != dim(stats)[2]) stop("length of initial parameters vector and number of columns of the array stats must be the same")}
    
    # optimization
    ott <- optim(par = pars_init,
                 fn = nllik,
                 method = method,
                 stats = stats,
                 event_binary = event_binary,
                 interevent_time = interevent_time,
                 hessian = TRUE,
                 threads = threads)

    # storing output
    out <- list()
    out$coef <- ott$par # MLE estimates
    out$hessian <- ott$hessian # hessian matrix
    out$llik <- -ott$value # loglikelihood value at MLE values (we take the '-value' because we minimized the '-loglik')
    out$AIC <- 2*length(out$coef) - 2*out$llik # AIC
    out$BIC <- length(out$coef)*log(length(interevent_time)) - 2*out$llik # BIC
    out$optim_output <- ott # optim output (if more information are needed)

    return(out)
}
