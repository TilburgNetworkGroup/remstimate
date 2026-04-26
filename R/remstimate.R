
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



# diagnostics
#' @title Compute the diagnostics of a \code{remstimate} object
#' @description A function that returns the diagnostics of a \code{remstimate} object. The output object of the method \code{diagnostics} contains the residuals of the model estimated in the \code{remstimate} object, and the event rates estimated from the model at each tiem point. For tie-oriented modeling frameworks the object contains: a list \code{residuals} with two objects, \code{standardized_residuals} containing standardized Schoenfeld's residuals (Schoenfeld, D., 1982, <doi:10.2307/2335876>; Grambsch, P. M., & Therneau, T. M., 1994, <doi:10.2307/2337123>; Winnett, A., & Sasieni, P., 2001, <jstor.org/stable/2673500>), and \code{smoothing_weights} (a matrix of weights used for the red smooth splines in the plot of the residuals), an array structure \code{rates} with the event rates estimated from the optimized model parameters, and \code{.reh.processed} which is a pseudo-hidden object containing a further processed \code{remify} object that helps speed up the plotting function \code{plot.remstimate} and that the user is not supposed to modify. As to the actor-oriented modeling frameworks, in the diagnostics output there are two main list objects named after \code{sender_model} and \code{receiver_model}. After selecting the model, the structure of diagnostics is the same as for the tie-oriented model. Each model's diagnostics (sender or receiver) is available only if the corresponding model is found in the \code{remstimate} object.
#' @param object is a \code{remstimate} object.
#' @param reh is a \code{remify} object, the same used for the 'remstimate' object.
#' @param stats is a \code{remstats} object, the same used for the 'remstimate' object.
#' @param ... further arguments to be passed to the 'diagnostics' method.
#' @export
#'
#' @return a object of class \code{"diagnostics" "remstimate"} with standardized Schoenfeld's residuals and estimated event rates given the optimized model parameters.
#'
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
#' # diagnostics
#' tie_diagnostics <- diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats)
#' names(tie_diagnostics)
#'
diagnostics <- function(object,reh,stats,...){
  UseMethod("diagnostics")
}

denormalize_reh <- function(reh) {
  if (is.null(reh$meta)) return(reh)  # already old-style
  attr(reh, "model")    <- reh$meta$model
  attr(reh, "directed") <- reh$meta$directed
  attr(reh, "ordinal")  <- reh$meta$ordinal
  attr(reh, "weighted") <- reh$meta$weighted
  attr(reh, "riskset")  <- reh$meta$riskset_source %||% reh$meta$riskset
  attr(reh, "origin")   <- reh$meta$origin
  attr(reh, "dyadID")        <- reh$ids$dyad
  attr(reh, "dyadIDactive")  <- reh$ids$dyad_active
  attr(reh, "actor1ID")      <- reh$ids$actor1
  attr(reh, "actor2ID")      <- reh$ids$actor2
  attr(reh, "indices_simultaneous_events") <- reh$indices_simultaneous_events
  reh
}

#' @describeIn diagnostics diagnostics of a 'remstimate' object
#' @method diagnostics remstimate
#' @export
diagnostics.remstimate <- function(object,reh,stats,...) {

  # if(attr(tie_hmc,"method")=="HMC"){
  #   message("Note. diagnostics() not yet supported for Bayesian fitted model.")
  # }else{
  reh <- denormalize_reh(reh)

  # processing dyadID actor1ID and actor2ID depending on pt/pe and on start and stop values

  # ... processing input arguments

  # ... remify object ('reh' input argument)
  if(!inherits(reh,"remify")){
    stop("'reh' must be a 'remify' object (see ?remify::remify).")
  }

  # ... model
  model <- ifelse(is.null(attr(reh, "model")),"",attr(reh, "model")) # attribute from reh object
  if(!(model %in% c("actor","tie"))){
    stop("attribute 'model' of input 'reh' must be either 'actor' or 'tie'")
  }
  model_object <- ifelse(is.null(attr(object, "model")),"",attr(object, "model"))
  model_stats <- ifelse(is.null(attr(stats, "model")),"",attr(stats, "model"))
  if((model_stats != model) | (model_object != model)){
    stop("attribute 'model' of input 'object' and 'stats' must be the same as the attribute of input 'reh'")
  }


  # ... active riskset? then overwrite two objects (this prevents from coding many ifelse() to switch between active-oriented riskset objects and full-oriented riskset objects)
  if((attr(reh,"riskset") == "active")){
    reh$D <- reh$activeD
    if(model == "tie"){
      attr(reh,"dyadID") <- attr(reh,"dyadIDactive")
      reh$omit_dyad <- list() # because "reh$omit_dyad$time" and "reh$omit_dyad$riskset" for riskset="active" are obsolete (will be removed from remify output in the future 3.x.x version)
      # check line 113 remify.cpp
    }
  }

  # ... type of likelihood
  ordinal <- attr(reh,"ordinal")

  # ... ncores
  if(!is.null(attr(object,"ncores"))){
    ncores <- attr(object,"ncores")
  }

  # ... omit dyad
  if(is.null(reh$omit_dyad)){
    reh$omit_dyad <- list()
  }

  # ... processing start and stop values and method from ("subset" and "method" attribute of remstats object)
  # we do this now because later the stats object will change dimensions and won't be a remstats object anymore
  if(all(inherits(stats,c("remstats","tomstats"),TRUE)) | all(inherits(stats,c("remstats","aomstats"),TRUE))){
    stats_attr_method <- attr(stats,"method")
    if(is.null(stats_attr_method)){
      stop("attribute 'method' not found inside object 'remstats'. Input argument 'stats' must be an object of class 'remstats' from the package 'remstats' (>=3.2.0)")
    }
    omit_dyad_receiver <- NULL
    if(stats_attr_method == "pe"){
      if(!is.null(attr(reh,"evenly_spaced_interevent_time"))){
        reh$intereventTime <- attr(reh,"evenly_spaced_interevent_time")
      }
      if(!is.null(reh$E)){ # reh$E is NULL only when there are no simultaneous events
        reh$M <- reh$E # overwriting dimension (we can do it because remstimate works only with reh$M so if the method is "pt", reh$M will remain so. For method "pe" we assign reh$E to reh$M
      }
    }
    if(!is.null(attr(stats,"subset"))){
      start_stop <- as.numeric(unlist(attr(stats,"subset")))
    }
    else{
      start_stop <- c(1,reh$M)
    }
  }else{
    stop("'stats' must be a 'remstats' object from the package 'remstats' (>= 3.2.0), suitable for tie-oriented modeling ('tomstats') or actor-oriented modeling ('aomstats')")
  }

  # ... stats
  model_formula <- variables_names <- where_is_baseline <- NULL
  if(model == "tie"){
    if(all(inherits(stats,c("remstats","tomstats"),TRUE))){
      if(!is.null(dimnames(stats)[[3]])){
        variables_names <- dimnames(stats)[[3]]
      }
      if(is.null(attr(stats,"formula"))){
        model_formula <- stats::as.formula(paste("~ ",paste(variables_names,collapse=" + ")))
      }
      else{
        model_formula <- attr(stats,"formula")
      }
      # is there a baseline term?
      if(any(tolower(variables_names) %in% c("baseline"))){
        where_is_baseline <- which(variables_names == "baseline")
      }
      stats <- aperm(stats, perm = c(2,3,1)) # stats reshaped in [D*U*M]
      # start and stop for tie-oriented model
      if(stats_attr_method == "pt"){
        attr(reh,"dyadID") <- attr(reh,"dyadID")[start_stop[1]:start_stop[2]] # this is already working for dyadIDactive because we reassing attribute dyadID in line 123
        attr(reh,"actor1ID") <- attr(reh,"actor1ID")[start_stop[1]:start_stop[2]]
        attr(reh,"actor2ID") <- attr(reh,"actor2ID")[start_stop[1]:start_stop[2]]
        if((length(reh$omit_dyad)>0) & !is.null(attr(reh,"indices_simultaneous_events"))){
          reh$omit_dyad$time <-  reh$omit_dyad$time[-attr(reh,"indices_simultaneous_events")]
        }
      }
      else if(stats_attr_method =="pe"){
        attr(reh,"dyadID") <- unlist(attr(reh,"dyadID"))[start_stop[1]:start_stop[2]] # this is already working for dyadIDactive because we reassing attribute dyadID in line 123
        attr(reh,"actor1ID") <- unlist(attr(reh,"actor1ID"))[start_stop[1]:start_stop[2]]
        attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID"))[start_stop[1]:start_stop[2]]
      }
      reh$intereventTime <- reh$intereventTime[start_stop[1]:start_stop[2]] # in line 168 we already re-assigned the intereventTime variable in case of method="pe", so line 224 is a valid processing for "pt" and "pe"
      reh$M <- diff(start_stop)+1
      if(length(reh$omit_dyad)>0){
        reh$omit_dyad$time <- reh$omit_dyad$time[start_stop[1]:start_stop[2]]
      }
    }else if(all(inherits(stats,c("remstats","aomstats"),TRUE))){
      stop("'remstats' object supplied cannot work for tie-oriented modeling")
    }

    # .. check on dimensions
    if(length(attr(reh,"dyadID")) != dim(stats)[3]){ # if the dimension of the processed intereventTime are different from the dimensions of the input stats object, then throw error
      stop("the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object")
    }
  }
  if(model == "actor"){
    model_formula <- list() # becomes a list
    if(all(inherits(stats,c("remstats","aomstats"),TRUE))){
      variables_rate <- variables_choice <- NULL
      if(!is.null(stats$sender_stats)){ # sender model is specified
        variables_rate <- dimnames(stats$sender_stats)[[3]]
        if(!is.null(attr(stats,"formula")$rate)){
          model_formula[["rate_model_formula"]] <- attr(stats,"formula")$rate
        }
        else{
          model_formula[["rate_model_formula"]] <- stats::as.formula(paste("~ ",paste(variables_rate,collapse=" + ")))
        }
        # is there a baseline term?
        if(any(tolower(variables_rate) %in% c("baseline"))){
          where_is_baseline <- which(variables_rate == "baseline")
        }
        stats$sender_stats <- aperm(stats$sender_stats, perm = c(2,3,1)) # stats reshaped in [N*U*M]
      }
      if(!is.null(stats$receiver_stats)){ # receiver model is specified
        variables_choice <- dimnames(stats$receiver_stats)[[3]]
        if(!is.null(attr(stats,"formula")$choice)){
          model_formula[["choice_model_formula"]] <- attr(stats,"formula")$choice
        }
        else{
          model_formula[["choice_model_formula"]] <- stats::as.formula(paste("~ ",paste(variables_choice,collapse=" + ")))
        }
        stats$receiver_stats <- aperm(stats$receiver_stats, perm = c(2,3,1)) # stats reshaped in [N*U*M]
      }

      # vector of variable names and list of model formulas
      variables_names <- list(sender_model = variables_rate, receiver_model = variables_choice)

      # start and stop for actor-oriented model
      if(stats_attr_method == "pt"){
        attr(reh,"actor1ID") <- attr(reh,"actor1ID")[start_stop[1]:start_stop[2]]
        attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID")[start_stop[1]:start_stop[2]]) # unlist here because of receiver-choice model
        if(is.null(reh$E)){ # this scenario can still happen
          reh$E <- reh$M
        }
        if(length(reh$omit_dyad)>0){

          # for the receiver model
          if(!is.null(stats$receiver_stats)){
            start_stop_time <- unique(reh$edgelist$time)[start_stop]
            lb_time <- min(which(reh$edgelist$time>=start_stop_time[1]))
            ub_time <- max(which(reh$edgelist$time<=start_stop_time[2]))
            omit_dyad_receiver <- list(time = reh$omit_dyad$time[lb_time:ub_time], riskset = reh$omit_dyad$riskset) # for the receiver model
          }

          # for the sender model (we process now the sender model because this will modify the reh$omit_dyad$time object)
          if(!is.null(stats$sender_stats)){
            if(!is.null(attr(reh,"indices_simultaneous_events"))){
              reh$omit_dyad$time <-  reh$omit_dyad$time[-attr(reh,"indices_simultaneous_events")][start_stop[1]:start_stop[2]] # for the sender model
            }
            else{
              reh$omit_dyad$time <-  reh$omit_dyad$time[start_stop[1]:start_stop[2]] # for the sender model
            }
          }
        }
      }else if(stats_attr_method == "pe"){
        attr(reh,"actor1ID") <- unlist(attr(reh,"actor1ID"))[start_stop[1]:start_stop[2]]
        attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID"))[start_stop[1]:start_stop[2]]
        if(length(reh$omit_dyad)>0){
          if(!is.null(stats$receiver_stats)){
            omit_dyad_receiver <- list(time = reh$omit_dyad$time[start_stop[1]:start_stop[2]]  , riskset = reh$omit_dyad$riskset)
          }
          reh$omit_dyad$time <- reh$omit_dyad$time[start_stop[1]:start_stop[2]]
        }
      }
      if(!ordinal){
        reh$intereventTime <- reh$intereventTime[start_stop[1]:start_stop[2]]
      }
      reh$M <- diff(start_stop)+1

      # .. check on dimensions
      no_correct_dimensions <- FALSE
      if(!is.null(stats$sender_stats)){
        if((length(attr(reh,"actor1ID")) != dim(stats$sender_stats)[3])){ # if the dimension of the processed intereventTime is different from the dimensions of the input stats object, then throw error
          no_correct_dimensions <- TRUE
        }
      }
      if(!is.null(stats$receiver_stats)){
        if((length(attr(reh,"actor2ID")) != dim(stats$receiver_stats)[3])){ # if the dimension of the edgelist is different from the dimensions of the input stats object, then throw error
          no_correct_dimensions <- TRUE
        }
      }
      if(no_correct_dimensions){
        stop("the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object")
      }

    }else if(all(inherits(stats,c("remstats","tomstats"),TRUE))){
      stop("'remstats' object supplied cannot work for actor-oriented modeling")
    }
  }

  # ... adjusting the intereventTime
  if(ordinal){
    reh$intereventTime <- c(1) # we can assign a vector of length 1 because the intereventTime will not be used from any ordinal likelihood
  }

  ## diagnostics per model type (tie-oriented and actor-oriented)

  if(attr(object,"model") == "tie"){ # tie-oriented modeling
    length_comparison <- length(object$coefficients) == dim(stats)[2] # check number of coefficients and dimension of 'stats'object
    name_comparison <- FALSE
    if(length_comparison){
      name_comparison <- all(names(object$coefficients) == dimnames(stats)[[2]]) # check that names are the same
    }
    if(!name_comparison | !length_comparison){
      stop("input 'object' not compatible with input 'stats'")
    }

    variables_names <- attr(object, "statistics")
    where_is_baseline <- attr(object,"where_is_baseline")
    select_vars <- if(is.null(where_is_baseline)) 1:length(variables_names) else c(1:length(variables_names))[-where_is_baseline]
    baseline_value <- if(is.null(where_is_baseline)) 0 else as.vector(object$coefficients)[where_is_baseline]
    stats <- if(length(select_vars)==1) array(stats[,select_vars,],dim=c(dim(stats)[1],1,dim(stats)[3])) else stats[,select_vars,]
    diagnostics <- list()
    diagnostics$residuals <- computeDiagnostics(pars = as.vector(object$coefficients)[select_vars],
                                                stats = stats,
                                                actor1 = list(),
                                                actor2 = list(),
                                                dyad = attr(reh,"dyadID"),
                                                omit_dyad = reh$omit_dyad,
                                                model = attr(reh,"model"),
                                                ncores = ncores,
                                                baseline = baseline_value,
                                                N = reh$N
    )
    colnames(diagnostics$residuals$smoothing_weights) <- variables_names[select_vars]
    # diagnostics$residuals$standardized_residuals <- do.call(rbind, diagnostics$residuals$standardized_residuals)
    # diagnostics$residuals$rates <- do.call(rbind, lapply(diagnostics$residuals$rates, t))
    # lambdas (event rates)
    diagnostics$rates <- diagnostics$residuals$rates
    diagnostics$residuals$rates <- NULL
    if(!is.null(attr(object,"where_is_baseline")) & (length(select_vars) < 1)){ # only the baseline is in the model
      diagnostics$residuals <- NULL
    }
  }else if(attr(object,"model") == "actor"){ # actor-oriented modeling
    compare_input_sender <- compare_input_receiver <- FALSE
    if(!is.null(stats[["sender_stats"]])){
      length_comparison <- length(object[["sender_model"]]$coefficients) == dim(stats[["sender_stats"]])[2] # check number of coefficients and dimension of 'stats'object
      name_comparison <- FALSE
      if(length_comparison){
        name_comparison <- all(names(object[["sender_model"]]$coefficients) == dimnames(stats[["sender_stats"]])[[2]]) # check that names are the same
      }
      if(!name_comparison | !length_comparison){
        compare_input_sender <- TRUE
      }
    }
    if(!is.null(stats[["receiver_stats"]])){
      length_comparison <- length(object[["receiver_model"]]$coefficients) == dim(stats[["receiver_stats"]])[2] # check number of coefficients and dimension of 'stats'object
      name_comparison <- FALSE
      if(length_comparison){
        name_comparison <- all(names(object[["receiver_model"]]$coefficients) == dimnames(stats[["receiver_stats"]])[[2]]) # check that names are the same
      }
      if(!name_comparison | !length_comparison){
        compare_input_receiver <- TRUE
      }
    }
    if(compare_input_sender | compare_input_receiver){
      stop("input 'object' not compatible with input 'stats'")
    }
    variables_names <- attr(object, "statistics")
    where_is_baseline <- attr(object,"where_is_baseline")
    senderRate <- c(TRUE,FALSE)
    which_model <- c("sender_model","receiver_model")
    which_stats <- c("sender_stats","receiver_stats")
    actor1ID_condition <- c(TRUE,FALSE)
    if(stats_attr_method == "pe"){
      actor1ID_condition <- as.logical(actor1ID_condition*FALSE) # actor1ID is needed unlist both for sender and receiver model
    }
    diagnostics <- list()
    # residuals
    for(i in 1:2){
      if(!is.null(stats[[which_stats[i]]])){
        actor1ID_ls <- if(actor1ID_condition[i]) attr(reh,"actor1ID") else unlist(attr(reh,"actor1ID"))
        omit_dyad_actor <- if(senderRate[i]) reh$omit_dyad else omit_dyad_receiver
        diagnostics[[which_model[i]]] <- list()
        baseline_value <- 0
        select_vars <- c(1:dim(stats[[which_stats[i]]])[2])
        if(senderRate[i]){ # only for sender model
          baseline_value <- if(is.null(where_is_baseline)) 0 else as.vector(object[[which_model[i]]]$coefficients)[where_is_baseline]
          select_vars <- if(is.null(where_is_baseline)) c(1:dim(stats[[which_stats[i]]])[2]) else c(1:dim(stats[[which_stats[i]]])[2])[-where_is_baseline]
        }
        stats[[which_stats[i]]] <- if(length(select_vars)==1) array(stats[[which_stats[i]]][,select_vars,],dim=c(dim(stats[[which_stats[i]]])[1],1,dim(stats[[which_stats[i]]])[3])) else stats[[which_stats[i]]][,select_vars,]
        diagnostics[[which_model[i]]] <- list()
        diagnostics[[which_model[i]]]$residuals <- computeDiagnostics(pars = as.vector(object[[which_model[i]]]$coefficients)[select_vars],
                                                                      stats = stats[[which_stats[i]]],
                                                                      actor1 = actor1ID_ls,
                                                                      actor2 = attr(reh,"actor2ID"),
                                                                      dyad = list(),
                                                                      omit_dyad = omit_dyad_actor,
                                                                      model = attr(reh,"model"),
                                                                      N = reh$N,
                                                                      senderRate = senderRate[i],
                                                                      ncores = ncores,
                                                                      baseline = baseline_value)
        colnames(diagnostics[[which_model[i]]]$residuals$smoothing_weights) <- if(senderRate[i]) variables_names[["sender_model"]][select_vars] else variables_names[["receiver_model"]][select_vars]
        # diagnostics$residuals$standardized_residuals <- do.call(rbind, diagnostics$residuals$standardized_residuals)
        # diagnostics$residuals$rates <- do.call(rbind, lapply(diagnostics$residuals$rates, t))
        # lambdas (event rates)
        diagnostics[[which_model[i]]]$rates <- diagnostics[[which_model[i]]]$residuals$rates
        diagnostics[[which_model[i]]]$residuals$rates <- NULL
        if(senderRate[i] & !is.null(where_is_baseline) & (length(select_vars) < 1)){ # only the baseline is in the model
          diagnostics[[which_model[i]]]$residuals <- NULL
        }
      }
    }
  }
  diagnostics$.reh.processed <- reh
  diagnostics$.reh.processed$stats.method <- stats_attr_method
  class(diagnostics) <- c("diagnostics","remstimate")
  return(diagnostics)
  #}
}

#' @method print diagnostics
#' @export
print.diagnostics <- function(x, ...) {
  reh <- x$.reh.processed
  cat("Diagnostics for Relational Event Model\n")
  cat("Model          :", reh$meta$model, "\n")
  if(reh$meta$model=="tie"){
    cat("Stats          :", paste(colnames(x$residuals$smoothing_weights), collapse = ", "), "\n")
  }else if(reh$meta$model=="actor"){
    cat("Stats sender   :", paste(colnames(x$sender_model$residuals$smoothing_weights), collapse = ", "), "\n")
    cat("Stats receiver :", paste(colnames(x$receiver_model$residuals$smoothing_weights), collapse = ", "), "\n")
  }
  invisible(x)
}

#######################################################################################
#######################################################################################


# plot.remstimate
#' @title Plot diagnostics of a \code{remstimate} object
#' @rdname plot.remstimate
#' @description A function that returns a plot of diagnostics given a 'remstimate' object and depending on the 'approach' attribute.
#' @param x is a \code{remstimate} object.
#' @param reh a \code{remify} object, the same used for the \code{remstimate} object.
#' @param diagnostics is a \code{'diagnostics' 'remstimate'} object.
#' @param which one or more numbers between 1 and 2. Plots described in order: (1) two plots: a Q-Q plot of the waiting times where theoretical quantiles (Exponential distribution with rate 1) are plotted against observed quantiles (these are calculated as the multiplication at each time point between the sum of the event rates and the corresponding waiting time, which should be distributed as an exponential with rate 1). Next to the q-q plot, a density plot of the rescaled waiting times (in red) vs. the theoretical distribution (exponential distribution with rate 1, in black). The observed density is truncated at the 99th percentile of the waiting times, (2) standardized Schoenfeld's residuals (per each variable in the model, excluding the baseline) with smoothed weighted spline (line in red). The Schoenfeld's residuals help understand the potential presence of time dependence of the effects of statistics specified in the model, (3) distributions of posterior draws with histograms (only for HMC method), (4) trace plots of posterior draws after thinning (only for HMC method).
#' @param effects [\emph{optional}] for tie-oriented modeling (model =  "tie"), the names of the statistics which the user wants to plot the diagnostics for (default value is set to all the statistics available inside the object 'diagnostics'). The user can specify this argument for the standardized Schoenfeld's residuals (\code{which = 2}), histograms of posterior distributions (\code{which = 3}) and trace plots (\code{which = 4}). Default value is \code{NULL}, selecting all the effects available in the 'remstimate' object.
#' @param sender_effects [\emph{optional}] for actor-oriented modeling (model =  "actor"), the names of the statistics as to the sender model which the user wants to plot the diagnostics for (default value is set to all the statistics available inside the object 'diagnostics'). The user can specify this argument for the standardized Schoenfeld's residuals (\code{which = 2}), histograms of posterior distributions (\code{which = 3}) and trace plots (\code{which = 4}). If the user wants to plot only the diagnostics of one or more effects of the sender model and at the same time wants to exclude the plots of the receiver model, then set argument \code{receiver_effects = NA} and specify the vector of effects to \code{sender_effects} (or leave it \code{sender_effects = NULL} for selecting all effects of the sender model). Default value is \code{NULL}, selecting all the effects available for the sender model in the 'remstimate' object.
#' @param receiver_effects [\emph{optional}] for actor-oriented modeling (model =  "actor"), the names of the statistics as to the receiver model which the user wants to plot the diagnostics for (default value is set to all the statistics available inside the object 'diagnostics'). The user can specify this argument for the standardized Schoenfeld's residuals (\code{which = 2}), histograms of posterior distributions (\code{which = 3}) and trace plots (\code{which = 4}). If the user wants to plot only the diagnostics of one or more effects of the receiver model and at the same time wants to exclude the plots of the sender model, then set argument \code{sender_effects = NA} and specify the vector of effects to \code{receiver_effects} (or leave it \code{receiver_effects = NULL} for selecting all effects of the receiver model). Default value is \code{NULL}, selecting all the effects available for the receiver model in the 'remstimate' object (\code{x}).
#' @param ... further arguments to be passed to the 'plot' method, for instance, the remstats object with statistics ('stats') when the object 'diagnostics' is not supplied.
#' @method plot remstimate
#'
#' @return no return value. The function plots the diagnostics of a 'remstimate' object.
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
#' # diagnostics
#' tie_diagnostics <- diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats)
#'
#' # plot
#' plot(x = tie_mle, reh  = tie_reh, diagnostics = tie_diagnostics)
#'
plot.remstimate <- function(x,
                            reh,
                            diagnostics = NULL,
                            which = c(1:4),
                            effects = NULL,
                            sender_effects = NULL,
                            receiver_effects = NULL,
                            ...)
{
  reh <- denormalize_reh(reh)
  # which plot
  selected <- which
  which <- rep(FALSE,4)
  which[selected] <- TRUE

  # hdi function for Highest density intervals
  hdi <- function(y) {
    # initializing output with NA's
    interval <- rep(NA,2)
    # sorting vector of numeric values
    y <- sort.int(as.numeric(y), method='quick') # method = "quick" is fast for large vectors (which is the case of posterior draws from HMC method)
    size <- length(y) # length of input vector y
    windows_size <- floor(size * 0.975) # elements to include in the credibility interval (take the largest integer), creating a bit more conservative intervals, (can be substituted with 'ceiling', selecting the smallest integers and being less conservative). For large vectors, the choice should not matter.
    if(windows_size < 2 | (size - windows_size) < 1){ # either the window is too narrow or we do not have enough data points
      return(interval)
    }
    omit <- size - windows_size # number of elements to be omitted
    lb <- y[1:omit]  # calculating candidates lower bounds
    ub <- y[-c(1:(size - omit))]  # calculating candidates upper bounds
    which_interval <- which.min(ub - lb)  # calculating intervals and selecting the smallest one
    if(length(which_interval) > 1) { # if there are multiple minima take the rightmost
      which_interval <- max(which_interval)
      interval <- c(lb[which_interval], ub[which_interval])
    } else { # otherwise return the only minimum found
      interval <- c(lb[which_interval], ub[which_interval])
    }
    return(interval)
  }

  if(attr(x,"model") != attr(reh,"model")){
    stop("'x' and 'reh' have different attribute 'model'")
  }
  # if diagnostics is NULL, then look for object 'stats' in '...' argument and compute diagnostics
  if(is.null(diagnostics)){
    additional_input_args <- list(...) # check for stats argument
    if(!any(names(additional_input_args) %in% "stats")){
      stop("'stats' must be provided if argument 'diagnostics' is NULL")
    }
    else{
      diagnostics <- diagnostics(object = x, reh = reh, stats = additional_input_args$stats)
    }
  }else if(!all(inherits(diagnostics,c("diagnostics","remstimate"),TRUE))){
    stop("'diagnostics' must be an object of class 'remstimate' 'diagnostics'")
  }

  # overwriting reh with one coming from diagnostics(), because there pt/pe methods and start/stop from remstats are already processed in there
  reh <- diagnostics$.reh.processed

  # saving current graphic parameters
  op <- par(no.readonly = TRUE)
  on.exit(expr = par(op))

  # tie-oriented modeling
  if(attr(x,"model") == "tie"){
    if(is.null(effects)){
      effects <- names(x$coefficients)
      effects_to_check <- effects
      which_effects <- 1:length(effects)
    }else{
      effects <- as.character(effects)
      effects_to_check <- effects
      available_effects <- names(x$coefficients)
      which_effects <- unlist(sapply(1:length(effects), function(y) which(available_effects == effects[y])))
      if(length(which_effects) == 0){
        par(op)
        stop("effects not found in object 'remstimate'")
      }
      else{
        effects <- effects[order(which_effects)]
        which_effects <- sort(which_effects)
      }
    }
    # (1) waiting times vs. theoretical distribution
    if(which[1L]){
      if(!attr(reh,"ordinal")){
        sum_rates <- lapply(diagnostics$rates,sum)
        observed <- sort(reh$intereventTime*unlist(sum_rates))
        # sum_rates <- rowSums(diagnostics$rates)
        # observed <- sort(reh$intereventTime * sum_rates)
        theoretical <- stats::qexp(p = c(1:reh$M)/reh$M, rate = 1)
        par(mfrow=c(1,2))
        plot(theoretical,observed, xlab = "Theoretical Quantiles",
             ylab = "Observed Quantiles",cex=0.8) # use bquote() / Observed quatiles : ~(t[m]-t[m-1])*Sigma[e](lambda[e](t[m])) / Theoretical Quantiles ~Exp(1)
        mtext(text = "Q-Q waiting times", side = 3, line = 2,cex=1.5)
        abline(a=0,b=1,lty=2,lwd=1.5)
        density_observed <- density(observed)
        if(max(density_observed$y,na.rm=TRUE)>0 & max(density_observed$y,na.rm=TRUE)!=Inf){
          density_observed$y[!is.na(density_observed$y)] <- density_observed$y[!is.na(density_observed$y)]/max(density_observed$y,na.rm=TRUE) # rescaling max density to 1 (to compare with dexp)
        }
        curve(dexp,from=min(observed),to=as.numeric(quantile(observed,probs=c(0.99))),col=1,lty=2,lwd=1.5,xlab="Waiting times",ylab="Density",ylim=c(0,1))
        lines(density_observed,col=2,lwd=1.5)
        mtext(text = "Density plot of waiting times", side = 3, line = 2,cex=1.5)
        par(op)
      }
    }

    # (2) standardized Schoenfeld's residuals
    if(which[2L] & !is.null(diagnostics$residuals)){
      # checking effects from 'remstimate' with effects from 'remstimate' 'diagnostics'
      available_effects <- colnames(diagnostics$residuals$smoothing_weights)
      if(!is.null(attr(x,"where_is_baseline")) & ("baseline" %in% tolower(effects))){
        effects_to_check <- effects_to_check[-which(effects_to_check == "baseline")]
      }
      if(length(effects_to_check)>0){ # model
        compare_remstimate_with_diagnostics <- prod(unlist(sapply(1:length(effects_to_check),function(y) effects_to_check[y] %in% available_effects)))
        if(!compare_remstimate_with_diagnostics){
          par(op)
          stop("one or more effects not found inside the object 'diagnostics'.")
        }
      }
      # the object diagnostics doesn't have residuals on the intercept, therefore we process the effects once again
      effects_diagnostics <- effects_to_check
      which_effects_diagnostics <- unlist(sapply(1:length(effects_diagnostics), function(y) which(available_effects == effects_diagnostics[y])))
      effects_diagnostics <- effects_diagnostics[order(which_effects_diagnostics)]
      which_effects_diagnostics <- sort(which_effects_diagnostics)

      P <- length(effects_diagnostics) # number of statistics
      n_repeats_per_time_point <- rep(1,dim(diagnostics$residuals$standardized_residuals)[1]) # for attr(residuals, "stats.method") = "pe"
      if(reh$stats.method == "pt"){
        n_repeats_per_time_point <- sapply(1:dim(diagnostics$residuals$standardized_residuals)[1], function(v) {
          dim(diagnostics$residuals$standardized_residuals[[v]])[1]
        })
      }
      if(!attr(reh,"ordinal")){
        t_p <- rep(cumsum(reh$intereventTime),n_repeats_per_time_point)
      }
      else{
        t_p <- rep(cumsum(rep(1,reh$M)),n_repeats_per_time_point)
      }
      for(p in 1:P){
        y_p <-  unlist(sapply(1:dim(diagnostics$residuals$standardized_residuals)[1], function(v) diagnostics$residuals$standardized_residuals[[v]][,which_effects_diagnostics[p]]))
        qrt_p <- quantile(y_p, probs=c(0.25,0.75))
        lb_p <- qrt_p[1] - diff(qrt_p)*1.5
        ub_p <- qrt_p[2] + diff(qrt_p)*1.5
        ylim_p <- c(min(y_p),max(y_p))
        #if(length(which(y_p<ub_p & y_p>lb_p)) >= floor(0.95*length(y_p))){ # reduce the ylim only if the bound lb_p and ub_p include more than the 90% of the observations
        ylim_p <- c(lb_p,ub_p)
        #}
        w_p <- rep(diagnostics$residuals$smoothing_weights[,which_effects_diagnostics[p]],n_repeats_per_time_point)
        w_p[w_p<0] <- w_p[w_p>1e04] <- 0.0
        if(all(w_p <= 0.0)){
          w_p <- rep(1/length(w_p),length(w_p)) # assigning equal weights
        }
        par(mfrow=c(1,1))
        plot(t_p,y_p,xlab = "Time", ylab = "Scaled Schoenfeld's residuals",ylim=ylim_p,col=grDevices::rgb(128,128,128,200,maxColorValue = 255)) # standardized Schoenfeld's residuals
        #lines(smooth.spline(x = t_p, y = y_p,w=w_p,cv=FALSE),lwd=3.5,col=2) # smoothed weighted spline of the residuals
        tryCatch(
          lines(smooth.spline(x = t_p, y = y_p, w=w_p, cv=FALSE), lwd=3.5, col=2),# smoothed weighted spline of the residuals
          error = function(e) NULL
        )
        abline(h=0,col="black",lwd=3,lty=2)
        mtext(text = effects_diagnostics[p], side = 3, line = 1,cex=1.5)
        par(op)
      }
    }

    # (3) histograms distribution of posterior draws (HMC method)
    if(which[3L] & (attr(x,"method") %in% c("HMC"))){
      P <- length(effects)
      for(p in 1:P){
        title_p <- bquote("Posterior distribution of " ~ beta[.(effects[p])])
        par(mfrow=c(1,1))
        hist(x$draws[,which_effects[p]],freq = FALSE, col = "lavender", main = title_p, xlab =  "Posterior draw")
        # posterior mean
        abline(v = x$coefficients[which_effects[p]], col = 2, lwd = 3.5, lty = 2)
        # hdi
        ci <- hdi(y = x$draws[,which_effects[p]])
        if(!all(is.na(ci))){
          abline(v = ci, col = 4, lwd = 3.5, lty = 2)
        }
        par(op)
      }
    }

    # (4) trace plots posterior draws (HMC method only)
    if(which[4L] & (attr(x,"method") == "HMC")){
      nchains <- attr(x,"nchains")
      P <- length(effects)
      for(p in 1:P){
        title_p <- bquote("Trace plot of " ~ beta[.(effects[p])])
        par(mfrow=c(1,1))
        if(nchains==1){
          plot(x$draws[,which_effects[p]], type= "l", main = title_p, ylab =  "Posterior draw", xlab = "Iteration",lwd=1.8)
          # posterior mean
          abline(h = x$coefficients[which_effects[p]], col=2, lwd=3.5, lty=2)
          # hdi
          ci <- hdi(y = x$draws[,which_effects[p]])
          if(!all(is.na(ci))){
            abline(h = ci, col = 4, lwd = 3.5, lty = 2)
          }
        }
        else{
          ndraws_per_chain <- length(x$draws[,which_effects[p]])/nchains
          seq_chains <- seq(1,length(x$draws[,which_effects[p]]),by=ndraws_per_chain)
          seq_chains <- cbind(seq_chains,seq_chains+ndraws_per_chain-1)
          chain_colors <- hcl.colors(n=nchains,palette="BuPu")
          # plotting first chain
          y_lim <- range(x$draws[,which_effects[p]])
          plot(x$draws[seq_chains[1,1]:seq_chains[1,2],which_effects[p]], type= "l", main = title_p, ylab = "Posterior draw", xlab = "Iterations \n (per each chain)", col = chain_colors[1],lwd=1.8,ylim=y_lim)
          for(chain in 2:nchains){
            lines(x$draws[seq_chains[chain,1]:seq_chains[chain,2],which_effects[p]], col=chain_colors[chain],lwd=1.8)
          }
          # posterior mean
          abline(h = x$coefficients[which_effects[p]], col=2, lwd=3.5, lty=2)
          # hdi
          ci <- hdi(y = x$draws[,which_effects[p]])
          if(!all(is.na(ci))){
            abline(h = ci, col = 4, lwd = 3.5, lty = 2)
          }
        }
        par(op)
      }
    }
  }else if(attr(x,"model") == "actor"){ # actor-oriented modeling
    effects_ls <- list(sender_model = sender_effects, receiver_model = receiver_effects)
    which_model <- c("sender_model","receiver_model")
    title_model <- c("Rate model (sender)","Choice model (receiver)")
    senderRate <- c(TRUE,FALSE)
    for(i in 1:2){
      if(!is.null(x[[which_model[i]]])){
        if(is.null(effects_ls[[which_model[i]]])){
          effects <- names(x[[which_model[i]]]$coefficients)
          effects_to_check <- effects # for checking later with statistics names inside diagnostics
          which_effects <- 1:length(effects)
        }
        else if(is.na(effects_ls[[which_model[i]]])){
          next
        }
        else{
          effects <- as.character(effects_ls[[which_model[i]]])
          effects_to_check <- effects # for checking later with statistics names inside diagnostics
          available_effects <- names(x[[which_model[i]]]$coefficients)
          which_effects <- unlist(sapply(1:length(effects), function(y) which(available_effects == effects[y])))
          if(length(which_effects) == 0){
            par(op)
            stop("effects not found in object 'remstimate'")
          }
          else{
            effects <- effects[order(which_effects)]
            which_effects <- sort(which_effects)
          }
        }
        # (1) waiting times vs. theoretical distribution
        if(which[1L]){
          if(!attr(reh,"ordinal") & i==1){
            sum_rates <- lapply(diagnostics[[which_model[i]]]$rates,sum)
            observed <- sort(reh$intereventTime*unlist(sum_rates))
            theoretical <- stats::qexp(p = c(1:reh$M)/reh$M, rate = 1)
            par(mfrow=c(1,2))
            plot(theoretical,observed, xlab = "Theoretical Quantiles",
                 ylab = "Observed Quantiles",cex=0.8) # use bquote() / Observed quatiles : ~(t[m]-t[m-1])*Sigma[e](lambda[e](t[m])) / Theoretical Quantiles ~Exp(1)
            mtext(text = "Q-Q waiting times", side = 3, line = 2,cex=1.5)
            mtext(text = title_model[i], side = 3, line = 1,cex=1)
            abline(a=0,b=1,lty=2,lwd=1.5)
            density_observed <- density(observed)
            if(max(density_observed$y,na.rm=TRUE)>0 & max(density_observed$y,na.rm=TRUE)!=Inf){
              density_observed$y[!is.na(density_observed$y)] <- density_observed$y[!is.na(density_observed$y)]/max(density_observed$y,na.rm=TRUE) # rescaling max density to 1 (to compare with dexp)
            }
            curve(dexp,from=min(observed),to=as.numeric(quantile(observed,probs=c(0.99))),col=1,lty=2,lwd=1.5,xlab="Waiting times",ylab="Density",ylim=c(0,1))
            lines(density_observed,col=2,lwd=1.5)
            mtext(text = "Density plot of waiting times", side = 3, line = 2,cex=1.5)
            mtext(text = title_model[i], side = 3, line = 1,cex=1)
            par(op)
          }
        }

        # (2) standardized Schoenfeld's residuals
        if(which[2L] & !is.null(diagnostics[[which_model[i]]]$residuals)){
          available_effects <- colnames(diagnostics[[which_model[i]]]$residuals$smoothing_weights)
          # checking effects from 'remstimate' with effects from 'remstimate' 'diagnostics'
          if(senderRate[i]){
            if(!is.null(attr(x,"where_is_baseline")) & ("baseline" %in% tolower(effects))){
              effects_to_check <- effects_to_check[-which(effects_to_check == "baseline")]
            }
          }
          if(length(effects_to_check)>0){
            compare_remstimate_with_diagnostics <- prod(unlist(sapply(1:length(effects_to_check),function(y) effects_to_check[y] %in% available_effects)))
            if(!compare_remstimate_with_diagnostics){
              par(op)
              stop("one or more effects not found inside the object 'diagnostics'.")
            }
          }
          # the object diagnostics doesn't have residuals on the intercept, therefore we process the effects once again
          effects_diagnostics <- effects_to_check
          which_effects_diagnostics <- unlist(sapply(1:length(effects_diagnostics), function(y) which(available_effects == effects_diagnostics[y])))
          effects_diagnostics <- effects_diagnostics[order(which_effects_diagnostics)]
          which_effects_diagnostics <- sort(which_effects_diagnostics)

          P <- length(effects_diagnostics) # number of statistics
          n_repeats_per_time_point <- rep(1,dim(diagnostics[[which_model[i]]]$residuals$standardized_residuals)[1]) # for attr(residuals, "stats.method") = "pe"
          if(i == 1){
            if(reh$stats.method == "pt"){
              n_repeats_per_time_point <- sapply(1:dim(diagnostics[[which_model[i]]]$residuals$standardized_residuals)[1], function(v) dim(diagnostics[[which_model[i]]]$residuals$standardized_residuals[[v]])[1])
            }
            if(!attr(reh,"ordinal")){
              t_p <- rep(cumsum(reh$intereventTime),n_repeats_per_time_point)
            }
            else{
              t_p <- rep(cumsum(rep(1,reh$M)),n_repeats_per_time_point)
            }
          }
          else if(i == 2){
            t_p <- cumsum(n_repeats_per_time_point)
          }
          for(p in 1:P){
            y_p <- unlist(sapply(1:dim(diagnostics[[which_model[i]]]$residuals$standardized_residuals)[1], function(v) diagnostics[[which_model[i]]]$residuals$standardized_residuals[[v]][,which_effects_diagnostics[p]]))
            qrt_p <- quantile(y_p, probs=c(0.25,0.75))
            lb_p <- qrt_p[1] - diff(qrt_p)*1.5
            ub_p <- qrt_p[2] + diff(qrt_p)*1.5
            #ylim_p <- c(min(y_p),max(y_p))
            #if(length(which(y_p<ub_p & y_p>lb_p)) >= floor(0.95*length(y_p))){ # reduce the ylim only if the bound lb_p and ub_p include more than the 95% of the observations
            ylim_p <- c(lb_p,ub_p)
            #}
            w_p <- rep(diagnostics[[which_model[i]]]$residuals$smoothing_weights[,which_effects_diagnostics[p]],n_repeats_per_time_point)
            w_p[w_p<0] <- w_p[w_p>1e04] <- 0.0
            if(all(w_p <= 0.0)){
              w_p <- rep(1/length(w_p),length(w_p)) # assigning equal weights
            }
            par(mfrow=c(1,1))
            plot(t_p,y_p,xlab = "Time", ylab = "Scaled Schoenfeld's residuals",ylim=ylim_p,col=grDevices::rgb(128,128,128,200,maxColorValue = 255)) # standardized Schoenfeld's residuals
            #lines(smooth.spline(x = t_p, y = y_p,w=w_p,cv=FALSE),lwd=3.5,col=2) # smoothed weighted spline of the residuals
            tryCatch(
              lines(smooth.spline(x = t_p, y = y_p, w=w_p, cv=FALSE), lwd=3.5, col=2),
              error = function(e) NULL
            )
            abline(h=0,col="black",lwd=3,lty=2)
            mtext(text = effects_diagnostics[p], side = 3, line = 2,cex=1.5)
            mtext(text = title_model[i], side = 3, line = 1, cex = 1)
            par(op)
          }
        }

        # (3) histograms distribution of posterior draws (HMC method)
        if(which[3L] & (attr(x,"method") %in% c("HMC"))){
          P <- length(effects)
          for(p in 1:P){
            title_p <- bquote("Posterior distribution of " ~ beta[.(effects[p])])
            par(mfrow=c(1,1))
            hist(x[[which_model[[i]]]]$draws[,which_effects[p]],freq = FALSE, main = "", col = "lavender", xlab =  "Posterior draw")
            # posterior mean
            abline(v = x[[which_model[[i]]]]$coefficients[which_effects[p]], col = 2, lwd = 3.5, lty = 2)
            # hdi
            ci <- hdi(y = x[[which_model[[i]]]]$draws[,which_effects[p]])
            if(!all(is.na(ci))){
              abline(v = ci, col = 4, lwd = 3.5, lty = 2)
            }
            # title and subtitle
            mtext(text = title_p, side = 3, line = 2,cex=1.5)
            mtext(text = title_model[i], side = 3, line = 1,cex=1)
            par(op)
          }
        }

        # (4) trace plots posterior draws (HMC method only)
        if(which[4L] & (attr(x,"method") == "HMC")){
          nchains <- attr(x,"nchains")
          P <- length(effects)
          for(p in 1:P){
            title_p <- bquote("Trace plot of " ~ beta[.(effects[p])])
            par(mfrow=c(1,1))
            if(nchains==1){
              plot(x[[which_model[i]]]$draws[,which_effects[p]], type= "l", main = "",  ylab =  "Posterior draw", xlab = "Iteration",lwd=1.8)
              # posterior mean
              abline(h = x[[which_model[i]]]$coefficients[which_effects[p]], col=2, lwd=3.5, lty=2)
              # hdi
              ci <- hdi(y = x[[which_model[i]]]$draws[,which_effects[p]])
              if(!all(is.na(ci))){
                abline(h = ci, col = 4, lwd = 3.5, lty = 2)
              }
              # title and subtitle
              mtext(text = title_p, side = 3, line = 2,cex=1.5)
              mtext(text = title_model[i], side = 3, line = 1,cex=1)
            }
            else{
              ndraws_per_chain <- length(x[[which_model[i]]]$draws[,which_effects[p]])/nchains
              seq_chains <- seq(1,length(x[[which_model[i]]]$draws[,which_effects[p]]),by=ndraws_per_chain)
              seq_chains <- cbind(seq_chains,seq_chains+ndraws_per_chain-1)
              chain_colors <- hcl.colors(n=nchains,palette="BuPu")
              # plotting first chain
              y_lim <- range(x[[which_model[i]]]$draws[,which_effects[p]])
              plot(x[[which_model[i]]]$draws[seq_chains[1,1]:seq_chains[1,2],which_effects[p]], type= "l", main = "", ylab = "Posterior draw", xlab = "Iterations (per each chain)", col = chain_colors[1],lwd=1.8,ylim=y_lim)
              for(chain in 2:nchains){
                lines(x[[which_model[i]]]$draws[seq_chains[chain,1]:seq_chains[chain,2],which_effects[p]], col=chain_colors[chain],lwd=1.8)
              }
              # posterior mean
              abline(h = x[[which_model[i]]]$coefficients[which_effects[p]], col=2, lwd=3.5, lty=2)
              # hdi
              ci <- hdi(y = x[[which_model[i]]]$draws[,which_effects[p]])
              if(!all(is.na(ci))){
                abline(h = ci, col = 4, lwd = 3.5, lty = 2)
              }
              # title and subtitle
              mtext(text = title_p, side = 3, line = 2,cex=1.5)
              mtext(text = title_model[i], side = 3, line = 1,cex=1)
            }
            par(op)
          }
        }
      }
    }
  }
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




