
#######################################################################################
#######################################################################################


# print.remstimate
#' @title Print out a quick overview of a \code{remstimate} object
#' @rdname print.remstimate
#' @description A function that prints out the estimates returned by a 'remstimate' object.
#' @param x is a \code{remstimate} object.
#' @param ... further arguments to be passed to the print method.
#' @method print remstimate
#'
#' @return no return value. Prints out the main characteristics of a 'remstimate' object.
#'
#' @export
#'
#' @examples
#'
#' # ------------------------------------ #
#' #       method 'print' for the         #
#' # ------------------------------------ #
#'
#' # loading data
#' data(ao_data)
#'
#' # processing event sequence with remify
#' ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")
#'
#' # specifying linear predictor (for sender rate and receiver choice model)
#' rate_model <- ~ 1 + remstats::indegreeSender()
#' choice_model <- ~ remstats::inertia() + remstats::reciprocity()
#'
#' # calculating statistics
#' ao_reh_stats <- remstats::remstats(reh = ao_reh,
#'                                    sender_effects = rate_model,
#'                                    receiver_effects = choice_model)
#'
#' # running estimation
#' ao_mle <- remstimate::remstimate(reh = ao_reh,
#'                                  stats = ao_reh_stats,
#'                                  method = "MLE",
#'                                  nsim = 100,
#'                                  ncores = 1)
#'
#' # print
#' ao_mle
#'
#' # ------------------------------------ #
#' #   for more examples check vignettes  #
#' # ------------------------------------ #
#'
print.remstimate<-function(x, ...){
  if (is.null(attr(x,"formula")))
    stop("invalid 'remstimate' object:  no 'formula' attribute")
  #if (!inherits(x, "remstimate"))
  #    warning("calling print.remstimate(<fake-remstimate-object>) ...")
  cat("Relational Event Model",paste("(",attr(x,"model")," oriented)",sep=""),"\n")
  if(attr(x,"approach") == "Frequentist"){
    if(attr(x,"model") == "tie"){
      cat("\nCoefficients:\n\n")
      print(x$coefficients)
      cat("\nNull deviance:",x$null.deviance,"\nResidual deviance:",x$residual.deviance,"\n")
      cat("AIC:",x$AIC,"AICC:",x$AICC,"BIC:",x$BIC)
      if(!is.null(x$WAIC)){
        cat(" WAIC:", x$WAIC, "\n\n")
      }else{
        cat("\n\n")
      }
    }
    else if(attr(x,"model") == "actor"){
      if(!is.null(x$sender_model)){  # printing sender model
        cat("\nCoefficients rate model **for sender**:\n\n")
        print(x$sender_model$coefficients)
        cat("\nNull deviance:",x$sender_model$null.deviance,"\nResidual deviance:",x$sender_model$residual.deviance,"\n")
        cat("AIC:",x$sender_model$AIC,"AICC:",x$sender_model$AICC,"BIC:",x$sender_model$BIC)
        if(!is.null(x$sender_model$WAIC)){
          cat(" WAIC:", x$sender_model$WAIC, "\n\n")
        }else{
          cat("\n\n")
        }
      }
      if(!is.null(x$sender_model) & !is.null(x$receiver_model)){ # if both models are estimated, separate the two print by a "-"
        cat(paste0(rep("-", getOption("width")),collapse = ""))
      }
      if(!is.null(x$receiver_model)){ # printing receiver model
        cat("\n\nCoefficients choice model **for receiver**:\n\n")
        print(x$receiver_model$coefficients)
        cat("\nNull deviance:",x$receiver_model$null.deviance,"\nResidual deviance:",x$receiver_model$residual.deviance,"\n")
        cat("AIC:",x$receiver_model$AIC,"AICC:",x$receiver_model$AICC,"BIC:",x$receiver_model$BIC)
        if(!is.null(x$receiver_model$WAIC)){
          cat(" WAIC:", x$receiver_model$WAIC, "\n\n")
        }else{
          cat("\n\n")
        }
      }
    }

  }
  else{ # Bayesian
    if(attr(x,"model") == "tie"){
      cat("\nPosterior Modes:\n\n")
      print(x$coefficients)
    }
    else if(attr(x,"model") == "actor"){
      if(!is.null(x$sender_model)){  # printing sender model
        cat("\nPosterior Modes rate model **for sender**:\n\n")
        print(x$sender_model$coefficients)
        cat("\n")
        # write some more info here
      }
      if(!is.null(x$sender_model) & !is.null(x$receiver_model)){ # if both models are estimated, separate the two print by a "-"
        cat(paste0(rep("-", getOption("width")),collapse = ""))
      }
      if(!is.null(x$receiver_model)){ # printing receiver model
        cat("\n\nPosterior Modes choice model **for sender**:\n\n")
        print(x$receiver_model$coefficients)
        # write some more info here
      }
    }
  }
}


#######################################################################################
#######################################################################################


# summary.remstimate
#' @title Generate the summary of a \code{remstimate} object
#' @rdname summary.remstimate
#' @description A function that returns the summary of a \code{remstimate} object.
#' @param object is a \code{remstimate} object.
#' @param ... further arguments to be passed to the 'summary' method.
#' @method summary remstimate
#'
#' @return no return value. Prints out the summary of a 'remstimate' object. The output can be save in a list, which contains the information printed out by the summary method.
#'
#' @export
#'
#' @examples
#'
#' # ------------------------------------ #
#' #       method 'summary' for the       #
#' #         tie-oriented model.         #
#' # ------------------------------------ #
#'
#' # loading data
#' data(tie_data)
#'
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
#'
#' # specifying linear predictor
#' tie_model <- ~ 1 +
#'                remstats::indegreeSender()+
#'                remstats::inertia()+
#'                remstats::reciprocity()
#'
#' # calculating statistics
#' tie_reh_stats <- remstats::remstats(reh = tie_reh,
#'                                     tie_effects = tie_model)
#'
#' # running estimation
#' tie_mle <- remstimate::remstimate(reh = tie_reh,
#'                                   stats = tie_reh_stats,
#'                                   method = "MLE",
#'                                   nsim = 100,
#'                                   ncores = 1)
#'
#' # summary
#' summary(tie_mle)
#'
#' # ------------------------------------ #
#' #      method 'summary' for the        #
#' #      actor-oriented model: "MLE"     #
#' # ------------------------------------ #
#'
#' # loading data
#' data(ao_data)
#'
#' # processing event sequence with remify
#' ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")
#'
#' # specifying linear predictor (for sender rate and receiver choice model)
#' rate_model <- ~ 1 + remstats::indegreeSender()
#' choice_model <- ~ remstats::inertia() + remstats::reciprocity()
#'
#' # calculating statistics
#' ao_reh_stats <- remstats::remstats(reh = ao_reh,
#'                                    sender_effects = rate_model,
#'                                    receiver_effects = choice_model)
#'
#' # running estimation
#' ao_mle <- remstimate::remstimate(reh = ao_reh,
#'                                  stats = ao_reh_stats,
#'                                  method = "MLE",
#'                                  nsim = 100,
#'                                  ncores = 1)
#'
#' # summary
#' summary(ao_mle)
#'
#' # ------------------------------------ #
#' #   for more examples check vignettes  #
#' # ------------------------------------ #
#'
summary.remstimate<-function (object, ...)
{
  if (inherits(object, "diagnostics")) return(summary.diagnostics(object, ...))
  # (codings are based on the structure of summary.lm() and summary.glm() from package 'stats')
  #if (!inherits(object, "remstimate"))
  #    warning("calling summary.remstimate(<fake-remstimate-object>) ...")
  summary_out <- list()
  if(attr(object, "approach") == "Frequentist"){ # Frequentist
    if(attr(object,"model") == "tie"){
      coefsTab <- cbind(object$coefficients,
                        object$se,
                        object$coefficients/object$se,
                        2*(1-stats::pnorm(abs(object$coefficients/object$se))),
                        (stats::dnorm(0,mean=object$coefficients,sd=object$se) / stats::dnorm(0,mean=0,sd=object$se*sqrt(object$df.null)))/((stats::dnorm(0,mean=object$coefficients,sd=object$se) / stats::dnorm(0,mean=0,sd=object$se*sqrt(object$df.null)))+1)
      )
      colnames(coefsTab) <- c("Estimate","Std. Err", "z value", "Pr(>|z|)", "Pr(=0)")
      rownames(coefsTab) <- attr(object, "statistics")
      summary_out$coefsTab <- coefsTab
      keep <- match(c("formula","aic",
                      "contrasts", "df.residual","null.deviance","df.null",
                      "iter", "na.action"), names(object), 0L)

      keep <- list(formula = attr(object,"formula"),
                   model = attr(object,"model"),
                   ordinal = attr(object,"ordinal"),
                   method = attr(object, "method"),
                   approach = attr(object, "approach"),
                   residual.deviance = object$residual.deviance,
                   null.deviance = object$null.deviance,
                   model.deviance = object$model.deviance,
                   df.residual = (object$df.null - object$df.model),
                   df.null = object$df.null,
                   df.model = object$df.model,
                   AIC = object$AIC,
                   AICC = object$AICC,
                   BIC = object$BIC,
                   WAIC = object$WAIC,
                   epsilon = attr(object, "epsilon"),
                   chiP = 1 - stats::pchisq(object$model.deviance, object$df.model))
      summary_out <- do.call(c, list(keep, summary_out))
    }
    else if(attr(object,"model") == "actor"){
      coefsTab <- list()
      if(!is.null(object$sender_model)){ # summary sender model
        coefsTab$sender_model <- cbind(object$sender_model$coefficients,
                                       object$sender_model$se,
                                       object$sender_model$coefficients/object$sender_model$se,
                                       2*(1-stats::pnorm(abs(object$sender_model$coefficients/object$sender_model$se))),
                                       (stats::dnorm(0,mean=object$sender_model$coefficients,sd=object$sender_model$se) / stats::dnorm(0,mean=0,sd=object$sender_model$se*sqrt(object$sender_model$df.null)))/((stats::dnorm(0,mean=object$sender_model$coefficients,sd=object$sender_model$se) / stats::dnorm(0,mean=0,sd=object$sender_model$se*sqrt(object$sender_model$df.null)))+1)
        )
        colnames(coefsTab$sender_model) <- c("Estimate","Std. Err", "z value", "Pr(>|z|)", "Pr(=0)")
        rownames(coefsTab$sender_model) <- names(object$sender_model$coefficients)
      }
      if(!is.null(object$receiver_model)){ # summary receiver model
        coefsTab$receiver_model <- cbind(object$receiver_model$coefficients,
                                         object$receiver_model$se,
                                         object$receiver_model$coefficients/object$receiver_model$se,
                                         2*(1-stats::pnorm(abs(object$receiver_model$coefficients/object$receiver_model$se))),
                                         (stats::dnorm(0,mean=object$receiver_model$coefficients,sd=object$receiver_model$se) / stats::dnorm(0,mean=0,sd=object$receiver_model$se*sqrt(object$receiver_model$df.null)))/((stats::dnorm(0,mean=object$receiver_model$coefficients,sd=object$receiver_model$se) / stats::dnorm(0,mean=0,sd=object$receiver_model$se*sqrt(object$receiver_model$df.null)))+1)
        )
        colnames(coefsTab$receiver_model) <- c("Estimate","Std. Err", "z value", "Pr(>|z|)", "Pr(=0)")
        rownames(coefsTab$receiver_model) <- names(object$receiver_model$coefficients)
      }
      summary_out$coefsTab <- coefsTab
      #keep <- match(c("formula",
      #  "contrasts", "sender_model", "receiver_model", "na.action"), names(object), 0L) # check if it works without this line
      keep <- list(formula = attr(object,"formula"),
                   model = attr(object,"model"),
                   ordinal = attr(object,"ordinal"),
                   method = attr(object, "method"),
                   approach = attr(object, "approach"),
                   epsilon = attr(object, "epsilon"))
      if(!is.null(object$sender_model)){
        keep$sender_model <- list(residual.deviance = object$sender_model$residual.deviance,
                                  null.deviance = object$sender_model$null.deviance,
                                  model.deviance = object$sender_model$model.deviance,
                                  df.residual = object$sender_model$df.residual,
                                  df.null = object$sender_model$df.null,
                                  df.model = object$sender_model$df.model,
                                  AIC = object$sender_model$AIC,
                                  AICC = object$sender_model$AICC,
                                  BIC = object$sender_model$BIC,
                                  WAIC = object$sender_model$WAIC,
                                  chiP = 1 - stats::pchisq(object$sender_model$model.deviance, object$sender_model$df.model)
        )
      }
      if(!is.null(object$receiver_model)){
        keep$receiver_model <- list(residual.deviance = object$receiver_model$residual.deviance,
                                    null.deviance = object$receiver_model$null.deviance,
                                    model.deviance = object$receiver_model$model.deviance,
                                    df.residual = object$receiver_model$df.residual,
                                    df.null = object$receiver_model$df.null,
                                    df.model = object$receiver_model$df.model,
                                    AIC = object$receiver_model$AIC,
                                    AICC = object$receiver_model$AICC,
                                    BIC = object$receiver_model$BIC,
                                    WAIC = object$receiver_model$WAIC,
                                    chiP = 1 - stats::pchisq(object$receiver_model$model.deviance, object$receiver_model$df.model)
        )
      }
      summary_out <- do.call(c, list(keep, summary_out))
    }
  }
  else if(attr(object, "approach") == "Bayesian"){ # Bayesian
    if(attr(object,"model") == "tie"){
      coefsTab <- cbind(object$coefficients,
                        object$sd,
                        apply(object$draws,2,stats::quantile,0.025),
                        apply(object$draws,2,stats::quantile,0.5),
                        apply(object$draws,2,stats::quantile,0.975),
                        (stats::dnorm(0,mean=object$coefficients,sd=object$sd) / stats::dnorm(0,mean=0,sd=object$sd*sqrt(object$df.null)))/((stats::dnorm(0,mean=object$coefficients,sd=object$sd) / stats::dnorm(0,mean=0,sd=object$sd*sqrt(object$df.null)))+1)
      )
      colnames(coefsTab) <- c("Post.Mode","Post.SD","Q2.5%","Q50%","Q97.5%","Pr(=0|x)")
      rownames(coefsTab) <- attr(object, "statistics")
      summary_out$coefsTab <- coefsTab

      keep <- list(formula = attr(object,"formula"),
                   model = attr(object,"model"),
                   ordinal = attr(object,"ordinal"),
                   method = attr(object, "method"),
                   approach = attr(object, "approach"),
                   statistics = attr(object, "statistics"),
                   prior = attr(object, "prior"),
                   nsim = attr(object, "nsim"),
                   seed = attr(object, "seed"),
                   loglik = object$loglik,
                   WAIC = object$WAIC
      )
      keep_HMC <- list()
      if(attr(object, "method") == "HMC"){
        keep_HMC <- list(
          nchains = attr(object, "nchains"),
          burnin = attr(object, "burnin"),
          thin = attr(object, "thin")
        )
      }
      summary_out <- do.call(c, list(keep, keep_HMC, summary_out))
    }
    if(attr(object,"model") == "actor"){
      coefsTab <- list()
      if(!is.null(object$sender_model)){ # summary sender model
        coefsTab$sender_model <- cbind(object$sender_model$coefficients,
                                       object$sender_model$sd,
                                       apply(object$sender_model$draws,2,stats::quantile,0.025),
                                       apply(object$sender_model$draws,2,stats::quantile,0.5),
                                       apply(object$sender_model$draws,2,stats::quantile,0.975),
                                       (stats::dnorm(0,mean=object$sender_model$coefficients,sd=object$sender_model$sd) / stats::dnorm(0,mean=0,sd=object$sender_model$sd*sqrt(object$sender_model$df.null)))/((stats::dnorm(0,mean=object$sender_model$coefficients,sd=object$sender_model$sd) / stats::dnorm(0,mean=0,sd=object$sender_model$sd*sqrt(object$sender_model$df.null)))+1)
        )
        colnames(coefsTab$sender_model) <- c("Post.Mode","Post.SD","Q2.5%","Q50%","Q97.5%","Pr(=0|x)")
        rownames(coefsTab$sender_model) <- names(object$sender_model$coefficients)
      }
      if(!is.null(object$receiver_model)){ # summary receiver model
        coefsTab$receiver_model <- cbind(object$receiver_model$coefficients,
                                         object$receiver_model$sd,
                                         apply(object$receiver_model$draws,2,stats::quantile,0.025),
                                         apply(object$receiver_model$draws,2,stats::quantile,0.5),
                                         apply(object$receiver_model$draws,2,stats::quantile,0.975),
                                         (stats::dnorm(0,mean=object$receiver_model$coefficients,sd=object$receiver_model$sd) / stats::dnorm(0,mean=0,sd=object$receiver_model$sd*sqrt(object$receiver_model$df.null)))/((stats::dnorm(0,mean=object$receiver_model$coefficients,sd=object$receiver_model$sd) / stats::dnorm(0,mean=0,sd=object$receiver_model$sd*sqrt(object$receiver_model$df.null)))+1)
        )
        colnames(coefsTab$receiver_model) <- c("Post.Mode","Post.SD","Q2.5%","Q50%","Q97.5%","Pr(=0|x)")
        rownames(coefsTab$receiver_model) <- names(object$receiver_model$coefficients)
      }
      summary_out$coefsTab <- coefsTab
      #keep <- match(c("formula",
      #  "contrasts", "sender_model", "receiver_model", "na.action"), names(object), 0L)
      keep <- list(formula = attr(object,"formula"),
                   model = attr(object,"model"),
                   ordinal = attr(object,"ordinal"),
                   method = attr(object, "method"),
                   approach = attr(object, "approach"),
                   epsilon = attr(object, "epsilon")) # add other attributes from the objects?
      if(!is.null(object$sender_model)){
        keep$sender_model <- list(
          #residual.deviance = object$sender_model$residual.deviance,
          #null.deviance = object$sender_model$null.deviance,
          #model.deviance = object$sender_model$model.deviance,
          #df.residual = object$sender_model$df.residual,
          df.null = object$sender_model$df.null,
          #df.model = object$sender_model$df.model,
          loglik = object$sender_model$loglik,
          WAIC = object$sender_model$WAIC
        )
      }
      if(!is.null(object$receiver_model)){
        keep$receiver_model <- list(
          #residual.deviance = object$receiver_model$residual.deviance,
          #null.deviance = object$receiver_model$null.deviance,
          #model.deviance = object$receiver_model$model.deviance,
          #df.residual = object$receiver_model$df.residual,
          df.null = object$receiver_model$df.null,
          #df.model = object$receiver_model$df.model,
          loglik = object$receiver_model$loglik,
          WAIC = object$receiver_model$WAIC
        )
      }
      summary_out <- do.call(c, list(keep, summary_out))
    }
  }
  else{
    if(length(summary_out)==0) stop("invalid 'remstimate' object")
  }
  class(summary_out) <- "summary.remstimate"

  cat("Relational Event Model",paste("(",summary_out$model," oriented)",sep=""),"\n\n")
  if(summary_out$model == "tie"){
    cat("Call:\n",deparse(summary_out$formula),"\n\n",sep="")
    second_line <- paste("(",summary_out$method," with ",sep="")
    if(summary_out$ordinal) second_line <- paste(second_line,"ordinal likelihood):\n\n",sep="")
    else{
      second_line <- paste(second_line,"interval likelihood):\n\n",sep="")
    }
    if(summary_out$approach == "Frequentist"){
      cat("\nCoefficients",second_line)
    }
    else{ # Bayesian
      cat("\nPosterior Modes",second_line)
    }
    stats::printCoefmat(summary_out$coefsTab, P.values = TRUE, signif.stars = FALSE, ...)
    if(summary_out$approach == "Frequentist"){
      cat("Null deviance:", summary_out$null.deviance, "on", summary_out$df.null, "degrees of freedom\n")
      cat("Residual deviance:", summary_out$residual.deviance, "on", summary_out$df.residual, "degrees of freedom\n")
      cat("Chi-square:", summary_out$model.deviance, "on", summary_out$df.model,
          "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(summary_out$model.deviance,
                                                                      summary_out$df.model), "\n")
      cat("AIC:", summary_out$AIC, "AICC:", summary_out$AICC, "BIC:", summary_out$BIC)
      if(!is.null(summary_out$WAIC)){
        cat(" WAIC:", summary_out$WAIC, "\n")
      }else{
        cat("\n")
      }
    }
    if(summary_out$approach == "Bayesian"){
      cat("Log posterior:",summary_out$loglik)
      if(!is.null(summary_out$WAIC)){
        cat(" WAIC:", summary_out$WAIC, "\n")
      }else{
        cat("\n")
      }
      # cat("Prior parameters:",paste(names(summary_out$prior.param),unlist(summary_out$prior.param),sep="="),"\n")
    }
  }
  else if(summary_out$model == "actor"){
    second_line <- paste("(",summary_out$method," with ",sep="")
    if(summary_out$ordinal){
      second_line <- paste(second_line,"ordinal likelihood):",sep="")
    }
    else{
      second_line <- paste(second_line,"interval likelihood):\n\n",sep="")
    }
    if(!is.null(summary_out$sender_model)){
      # sender rate summary
      cat("Call rate model **for sender**:\n\n\t",deparse(summary_out$formula$rate_model_formula),"\n\n",sep="")
      if(summary_out$approach == "Frequentist"){
        cat("\nCoefficients rate model",second_line)
      }
      else{ # Bayesian
        cat("\nPosterior Modes rate model",second_line)
      }
      stats::printCoefmat(summary_out$coefsTab$sender_model, P.values = TRUE, signif.stars = FALSE)
      if(summary_out$approach == "Frequentist"){
        cat("Null deviance:", summary_out$sender_model$null.deviance, "on", summary_out$sender_model$df.null, "degrees of freedom\n")
        cat("Residual deviance:", summary_out$sender_model$residual.deviance, "on", summary_out$sender_model$df.residual, "degrees of freedom\n")
        cat("Chi-square:", summary_out$sender_model$model.deviance, "on", summary_out$sender_model$df.model,
            "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(summary_out$sender_model$model.deviance,
                                                                        summary_out$sender_model$df.model), "\n")
        cat("AIC:", summary_out$sender_model$AIC, "AICC:", summary_out$sender_model$AICC, "BIC:", summary_out$sender_model$BIC)
        if(!is.null(summary_out$sender_model$WAIC)){
          cat(" WAIC:", summary_out$sender_model$WAIC, "\n")
        }else{
          cat("\n")
        }
      }
      if(summary_out$approach == "Bayesian"){
        cat("Log posterior:",summary_out$sender_model$loglik)
        if(!is.null(summary_out$sender_model$WAIC)){
          cat(" WAIC:", summary_out$sender_model$WAIC, "\n")
        }else{
          cat("\n")
        }
        # cat("Prior parameters:",paste(names(summary_out$sender_model$prior.param),unlist(summary_out$sender_model$prior.param),sep="="),"\n")
      }
    }
    if(!is.null(summary_out$sender_model) & !is.null(summary_out$receiver_model)){ # if both models are estimated, separate the two print by a "-"
      cat(paste0(rep("-", getOption("width")),collapse = ""),"\n\n")
    }
    if(!is.null(summary_out$receiver_model)){
      # receiver choice summary
      cat("Call choice model **for receiver**:\n\n\t",deparse(summary_out$formula$choice_model_formula),"\n\n",sep="")
      if(summary_out$approach == "Frequentist"){
        cat("\nCoefficients choice model",second_line)
      }
      else{ # Bayesian
        cat("\nPosterior Modes choice model",second_line)
      }
      stats::printCoefmat(summary_out$coefsTab$receiver_model, P.values = TRUE, signif.stars = FALSE)
      if(summary_out$approach == "Frequentist"){
        cat("Null deviance:", summary_out$receiver_model$null.deviance, "on", summary_out$receiver_model$df.null, "degrees of freedom\n")
        cat("Residual deviance:", summary_out$receiver_model$residual.deviance, "on", summary_out$receiver_model$df.residual, "degrees of freedom\n")
        cat("Chi-square:", summary_out$receiver_model$model.deviance, "on", summary_out$receiver_model$df.model,
            "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(summary_out$receiver_model$model.deviance,
                                                                        summary_out$receiver_model$df.model), "\n")
        cat("AIC:", summary_out$receiver_model$AIC, "AICC:", summary_out$receiver_model$AICC, "BIC:", summary_out$receiver_model$BIC)
        if(!is.null(summary_out$receiver_model$WAIC)){
          cat(" WAIC:", summary_out$receiver_model$WAIC, "\n")
        }else{
          cat("\n")
        }
      }
      if(summary_out$approach == "Bayesian"){
        cat("Log posterior:",summary_out$receiver_model$loglik)
        if(!is.null(summary_out$receiver_model$WAIC)){
          cat(" WAIC:", summary_out$receiver_model$WAIC, "\n")
        }else{
          cat("\n")
        }
        # cat("Prior parameters:",paste(names(summary_out$receiver_model$prior.param),unlist(summary_out$receiver_model$prior.param),sep="="),"\n")
      }
    }
  }

  invisible(summary_out)
}


#######################################################################################
#######################################################################################

#' AIC for remstimate objects
#' @description Returns the AIC (Akaike's Information Criterion) value of a 'remstimate' object.
#' @param object a \code{remstimate} object.
#' @param ... further arguments passed to \code{\link[stats]{AIC}}.
#' @return AIC value.
#' @importFrom stats AIC
#' @method AIC remstimate
#' @export
AIC.remstimate <- function(object,...) {
  if(attr(object, "approach") == "Frequentist"){
    if(attr(object,"model") == "tie"){
      return(object$AIC)
    }
    else if(attr(object,"model") == "actor"){
      AIC <- NULL
      if(!is.null(object$sender_model)){
        AIC <- c("sender model" = object$sender_model$AIC)
      }
      if(!is.null(object$receiver_model)){
        AIC <- c(AIC, "receiver model" = object$receiver_model$AIC)
      }
      return(AIC)
    }
  }
  else{
    stop("'approach' must be 'Frequentist'")
  }
}


#######################################################################################
#######################################################################################


# AICC
#' @title AICC
#' @description A function that returns the AICC (Akaike's Information Corrected Criterion) value in a 'remstimate' object.
#' @param object is a \code{remstimate} object.
#' @param ... further arguments to be passed to the 'AICC' method.
#'
#' @return AICC value of a 'remstimate' object.
#'
#' @export
#'
#' @examples
#'
#' # ------------------------------------ #
#' #       tie-oriented model: "MLE"      #
#' # ------------------------------------ #
#'
#' # loading data
#' data(tie_data)
#'
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
#'
#' # specifying linear predictor
#' tie_model <- ~ 1 +
#'                remstats::indegreeSender()+
#'                remstats::inertia()+
#'                remstats::reciprocity()
#'
#' # calculating statistics
#' tie_reh_stats <- remstats::remstats(reh = tie_reh,
#'                                     tie_effects = tie_model)
#'
#' # running estimation
#' tie_mle <- remstimate::remstimate(reh = tie_reh,
#'                                   stats = tie_reh_stats,
#'                                   method = "MLE",
#'                                   ncores = 1)
#'
#' # AICC
#' AICC(tie_mle)
#'
AICC <- function(object,...){
  UseMethod("AICC")
}


#' @describeIn AICC AICC (Akaike's Information Corrected Criterion) value of a 'remstimate' object
#' @method AICC remstimate
#' @export
AICC.remstimate <- function(object,...) {
  if(attr(object, "approach") == "Frequentist"){
    if(attr(object,"model") == "tie"){
      return(object$AICC)
    }
    else if(attr(object,"model") == "actor"){
      AICC <- NULL
      if(!is.null(object$sender_model)){
        AICC <- c("sender model" = object$sender_model$AICC)
      }
      if(!is.null(object$receiver_model)){
        AICC <- c(AICC, "receiver model" = object$receiver_model$AICC)
      }
      return(AICC)
    }
  }
  else{
    stop("'approach' must be 'Frequentist'")
  }
}


#######################################################################################
#######################################################################################


#' BIC for remstimate objects
#' @description Returns the BIC (Bayesian Information Criterion) value of a 'remstimate' object.
#' @param object a \code{remstimate} object.
#' @param ... further arguments passed to \code{\link[stats]{BIC}}.
#' @return BIC value.
#' @importFrom stats BIC
#' @method BIC remstimate
#' @export
BIC.remstimate <- function(object,...) {
  if(attr(object, "approach") == "Frequentist"){
    if(attr(object,"model") == "tie"){
      return(object$BIC)
    }
    else if(attr(object,"model") == "actor"){
      BIC <- NULL
      if(!is.null(object$sender_model)){
        BIC <- c("sender model" = object$sender_model$BIC)
      }
      if(!is.null(object$receiver_model)){
        BIC <- c(BIC, "receiver model" = object$receiver_model$BIC)
      }
      return(BIC)
    }
  }
  else{
    stop("'approach' must be 'Frequentist'")
  }
}


#######################################################################################
#######################################################################################


# WAIC
#' @title WAIC
#' @description A function that returns the WAIC (Watanabe-Akaike's Information Criterion) value in a 'remstimate' object.
#' @param object is a \code{remstimate} object.
#' @param ... further arguments to be passed to the 'WAIC' method.
#'
#' @return WAIC value of a 'remstimate' object.
#'
#' @export
#'
#' @examples
#'
#' # No examples available at the moment
#'
WAIC <- function(object,...){
  UseMethod("WAIC")
}


#' @describeIn WAIC WAIC (Watanabe-Akaike's Information Criterion) value of a 'remstimate' object
#' @method WAIC remstimate
#' @export
WAIC.remstimate <- function(object,...) {
  if(attr(object,"model") == "tie"){
    if(is.null(object$WAIC)){
      stop("WAIC was not computed. Re-run remstimate() with WAIC = TRUE.")
    }
    return(object$WAIC)
  }
  else if(attr(object,"model") == "actor"){
    WAIC <- NULL
    if(!is.null(object$sender_model)){
      if(is.null(object$sender_model$WAIC)){
        stop("WAIC was not computed. Re-run remstimate() with WAIC = TRUE.")
      }
      WAIC <- c("sender model" = object$sender_model$WAIC)
    }
    if(!is.null(object$receiver_model)){
      if(is.null(object$receiver_model$WAIC)){
        stop("WAIC was not computed. Re-run remstimate() with WAIC = TRUE.")
      }
      WAIC <- c(WAIC, "receiver model" = object$receiver_model$WAIC)
    }
    return(WAIC)
  }
}


#######################################################################################
#######################################################################################




