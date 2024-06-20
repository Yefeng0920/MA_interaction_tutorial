############################################################################
#' @title i2_ml
#' @description I2 (I-squared) for mulilevel meta-analytic models, based on Nakagawa & Santos (2012). Under multilevel models, we can have multiple I2 (see also Senior et al. 2016). Alternatively, the method proposed by Wolfgang Viechtbauer can also be used.
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param method Method used to calculate I2. Two options exist: a ratio-based calculation proposed by Nakagawa & Santos (\code{"ratio"}), or Wolfgang Viechtbauer's matrix method (\code{"matrix"}).
#' @param data Data frame used to fit the model.
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng_1 <- i2_ml(english_MA, data = english, boot = 10)
#' I2_eng_2 <- i2_ml(english_MA, data = english, method = "ratio")
#' I2_eng_3 <- i2_ml(english_MA, data = english, method = "matrix")
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' I2_fish_1 <- i2_ml(model, data = warm_dat, boot = 10)
#' I2_fish_2 <- i2_ml(model, method = c("matrix"),data = warm_dat)
#' I2_fish_2 <- i2_ml(model, method = c("ratio"),data = warm_dat)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' I2_lim_1 <- i2_ml(lim_MR, data=lim, boot = 10)
#' I2_lim_2 <- i2_ml(lim_MR, data=lim)
#' }
#' @references Senior, A. M., Grueber, C. E., Kamiya, T., Lagisz, M., O’Dwyer, K., Santos, E. S. A. & Nakagawa S. 2016. Heterogeneity in ecological and evolutionary meta-analyses: its magnitudes and implications. Ecology 97(12): 3293-3299.
#'  Nakagawa, S, and Santos, E.S.A. 2012. Methodological issues and advances in biological meta-analysis.Evolutionary Ecology 26(5): 1253-1274.
#' @export

i2_ml <- function(model, method = c("ratio", "matrix"), data, boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment i2_ml cannot take models with heterogeneous variance.")}
  
  ## evaluate choices
  method <- match.arg(method)
  
  if (method == "matrix") {
    # Wolfgang Viechtbauer's method
    I2s <- matrix_i2(model)
  } else {
    # Nakagawa & Santos (2012)
    I2s <- ratio_i2(model)
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot)
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    
    I2_each <- sapply(sim, function(ysim) {
      
      # The model
      tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                       mods = mods_formula,
                                       random = random_formula,
                                       data = data))
      pb$tick()
      Sys.sleep(1 / boot)
      
      if(method == "matrix"){
        I2 <- matrix_i2(tmp)
      } else {
        I2 <- ratio_i2(tmp)
      }
      
      return(I2) })
    
    # Summarise the bootstrapped distribution.
    I2s_each_95 <- data.frame(t(apply(I2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    I2s <-  round(I2s_each_95, digits = 3)
    colnames(I2s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(I2s)
}

#' @title matrix_i2
#' @description Calculated I2 (I-squared) for mulilevel meta-analytic models, based on a matrix method proposed by Wolfgang Viechtbauer.
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @examples
#' \dontrun{
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng <- i2_ml(english_MA, data = english, method = "matrix")
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence. The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
#' random=list(~1|Article, ~1|Datapoint), data=lim)
#' I2_lim <- i2_ml(lim_MR, data=lim, method = "matrix")
#' }
#' @export
matrix_i2 <- function(model){
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  W <- solve(model$V)
  X <- model.matrix(model)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2_total <- 100* (sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
  I2_each <- 100* (model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
  names(I2_each) <- paste0("I2_", model$s.names)
  names(I2_total) <- "I2_Total"
  I2s <- c(I2_total, I2_each)
  return(I2s)
}


#' @title ratio_i2
#' @description I2 (I-squared) for mulilevel meta-analytic models based on Nakagawa & Santos (2012). Under multilevel models, we can have a multiple I2 (see also Senior et al. 2016).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @examples
#' \dontrun{
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt,
#' sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng_1 <- i2_ml(english_MA, data = english, boot = 1000)
#' I2_eng_2 <- i2_ml(english_MA, data = english, method = "ratio")
#' }
#' @export
ratio_i2 <- function(model){
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / model$vi) * (model$k - 1) /
    (sum(1 / model$vi)^2 - sum((1 / model$vi)^2))
  
  # s^2_t = total variance
  I2_total <- 100 * (sum(model$sigma2) / (sum(model$sigma2) + sigma2_v))
  I2_each <- 100 * (model$sigma2 / (sum(model$sigma2) + sigma2_v))
  names(I2_each) <- paste0("I2_", model$s.names)
  names(I2_total) <- "I2_Total"
  
  I2s <- c(I2_total, I2_each)
  return(I2s)
}

############################################################################
#' @title r2_ml
#' @description R2 (R-squared) for mixed (mulitlevel) models, based on Nakagawa & Schielzeth (2013).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @param data Data frame used to fit the \code{rma.mv} or \code{rma} model object
#' @param boot The number of parametric bootstrap iterations, if desired. Defaults to \code{NULL}. A setting of 1000 is recommended as a minimum number of iterations.
#' @return A data frame containing all model results, including: mean effect size estimate, confidence and prediction intervals, with estimates converted back to r.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @references Nakagawa, S, and Schielzeth, H. 2013. A general and simple method for obtaining R2 from generalized linear mixed‐effects models. *Methods in Ecology and Evolution* 4(2): 133-142.
#' @examples
#' \dontrun{
#' data(lim)
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' R2 <- r2_ml(lim_MR,data=lim, boot = 10)
#' }
#' @export
r2_ml <- function(model, data, boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment r2_ml cannot take models with heterogeneous variance.")}
  
  R2 <- R2_calc(model)
  
  if(!is.null(boot)){
    
    if(any(class(model) %in% c("robust.rma")) == TRUE){stop("Sorry, bootstrapping currently doesn't work for robust.rma objects. Please use rma.mv instead.")}
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot)
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    # Parametric bootstrap
    R2 <- sapply(sim, function(ysim) {
      # The model
      tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                       mods = mods_formula,
                                       random = random_formula,
                                       data = data))
      R2s <- R2_calc(tmp)
      pb$tick()
      Sys.sleep(1 / boot)
      return(R2s)
    })
    
    # Summarise the bootstrapped distribution.
    R2 <- data.frame(t(apply(R2, 1, stats::quantile, probs=c(0.5, .025, .975))))
    R2 <-  round(R2, digits = 3)
    colnames(R2) = c("Est.", "2.5%", "97.5%")
  }
  
  return(R2)
  
}

#' @title R2_calc
#' @description Calculated R2 (R-squared) for mixed (mulitlevel) models, based on Nakagawa & Schielzeth (2013).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(lim)
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
#' random=list(~1|Article, ~1|Datapoint), data=lim)
#' R2 <- R2_calc(lim_MR)
#' }
#' @export

R2_calc <- function(model){
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  # fixed effect variance
  fix <- stats::var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  
  # conditional. Need to remove 'residual' variance; assume this is the sigma level with the largest k. Really the only way we can get that to work.
  R2c <- (fix + sum(model$sigma2) - model$sigma2[which(model$s.nlevels.f == max(model$s.nlevels.f))]) /
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_conditional = R2c)
  return(R2s)
}
############################################################################


#############
#' @title cv_ml
#' @description CV (I-squared) for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple CV. TODO - we need to cite original CV paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' CV_eng_1 <- cv_ml(english_MA, boot = 10)
#' CV_eng_2 <- cv_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' CV_fish_1 <- cv_ml(model, boot = 10)
#' CV_fish_2 <- cv_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' CV_lim_1 <- cv_ml(lim_MR, boot = 10)
#' CV_lim_2 <- cv_ml(lim_MR)
#' }
#' @references TODO
#' @export

cv_ml <- function(model,
                  boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment cv_ml cannot take models with heterogeneous variance.")}
  
  CVs <- ml_cv(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      CV_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        CV <- ml_cv(tmp)
      })
    } else{
      CV_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        CV <- ml_cv(tmp)
        return(CV) })
    }
    # Summarise the bootstrapped distribution.
    CVs_each_95 <- data.frame(t(apply(CV_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    CVs <-  round(CVs_each_95, digits = 3)
    colnames(CVs) = c("Est.", "2.5%", "97.5%")
  }
  
  return(CVs)
}


#' @title ml_cv
#' @description Calculated CV for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cv <- function(model){
  
  # total cv
  CV_total <- sqrt(sum(model$sigma2)) / abs(model$beta[[1]])
  # cv at different levels
  CV_each <-  sqrt(model$sigma2) / abs(model$beta[[1]])
  names(CV_each) <- paste0("CV_", model$s.names)
  names(CV_total) <- "CV_Total"
  
  CVs <- c(CV_total, CV_each)
  
  return(CVs)
}