############# Key functions #############

#' @title mod_results
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, this function creates a table of model results containing the mean effect size estimates for all levels of a given categorical moderator, and their corresponding confidence and prediction intervals. The function is capable of calculating marginal means from meta-regression models, including those with multiple moderator variables of mixed types (i.e. continuous and categorical variables).
#' @param model \code{rma.mv} model object
#' @param mod Moderator variable of interest that one wants marginal means for. Defaults to the intercept, i.e. \code{"1"}.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding variables in \code{by}. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param subset Used when one wishes to only plot a subset of levels within the main moderator of interest defined by \code{mod}. Default is \code{FALSE}, but use \code{TRUE} if you wish to subset levels of a moderator plotted (defined by \code{mod}) for plotting. Levels one wishes to plot are specified as a list, with the level names as a character string in the \code{at} argument. For subsetting to work, the \code{at} argument also needs to be specified so that the \code{mod_results} function knows what levels one wishes to plot.
#' @param N The name of the column in the data specifying the sample size so that each effect size estimate is scaled to the sample size, N. Defaults to \code{NULL}, so that precision is used for scaling each raw effect size estimate instead of sample size.
#' @param upper Logical, defaults to \code{TRUE}, indicating that the first letter of the character string for the moderator variable should be capitalized.
#' @param ... Additional arguments passed to \code{emmeans::emmeans()}.
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples \dontrun{
#' # Simple eklof data
#' data(eklof)
#' eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
#' m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment, data = eklof)
#' # Add the unit level predictor
#' eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
#' # fit a MLMR - accouting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type, random=list(~1|ExptID,
#' ~1|Datapoint), data = eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")
#'
#' # Fish example demonstrating marginalised means
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t",
#' control=list(optimizer="optim", optmethod="Nelder-Mead"), data = warm_dat)
#'   overall <- mod_results(model, group = "group_ID")
#' across_trait <- mod_results(model, group = "group_ID", mod = "trait.type")
#' across_trait_by_degree_diff <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")
#' across_trait_by_degree_diff_at_treat_end_days10 <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = 10),
#' by = "deg_dif",data = warm_dat)
#' across_trait_by_degree_diff_at_treat_end_days10And50 <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15),
#'  treat_end_days = c(10, 50)), by = "deg_dif")
#' across_trait_by_treat_end_days10And50 <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)),
#' by = "treat_end_days")
#' across_trait_by_treat_end_days10And50_ordinaryMM <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)),
#' by = "treat_end_days", weights = "prop")
#'
#' # Fish data example with a heteroscedastic error
#' model_het <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 + trait.type| es_ID), mods = ~ trait.type + deg_dif, method = "REML", test = "t", rho = 0, struc = "HCS", control=list(optimizer="optim", optmethod="Nelder-Mead"), data = warm_dat)
#' HetModel <- mod_results(model_het, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", weights = "prop")
#' orchard_plot(HetModel, xlab = "lnRR")
#' }
#' @export
#'
#'
# We will need to make sure people use "1" or"moderator_names"

mod_results <- function(model, mod = "1", group,  N = NULL,  weights = "prop", by = NULL, at = NULL, subset = FALSE, upper = TRUE, ...){
  
  if(any(grepl("-1|0", as.character(model$formula.mods)))){
    warning("It is recommended that you fit the model with an intercept. Unanticipated errors can occur otherwise.")
  }
  
  if(any(model$struct %in% c("GEN", "HCS"))){
    warning("We noticed you're fitting an ~inner|outer rma model ('random slope'). There are circumstances where the prediction intervals for such models are calculated incorrectly. Please check your results carefully.")
  }
  
  if(missing(model)){
    stop("Please specify the 'model' argument by providing rma.mv or rma model object. See ?mod_results")
  }
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class rma.mv, rma, or robust.rma")}
  
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  
  
  if(is.null(stats::formula(model))){ ##**NOTE** Not sure we need this bit of code anymore. Left here for now
    #model <- stats::update(model, "~1")
    model$formula.mods <- ~ 1
    #dat_tmp <- model$data$`1` <- "Intrcpt"
    #model$data <- dat_tmp
  }
  
  if(model$test == "t"){
    df_mod = as.numeric(model$ddf[[1]])
  } else{
    df_mod = 1.0e6 # almost identical to z value
  }
  
  # Extract the data from the model object
  data <- model$data 
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(is.character(data[[mod]]) | is.factor(data[[mod]]) | is.null(data[[mod]])) {
    grid <- emmeans::qdrg(formula = stats::formula(model), at = at, data = data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k-1) ## NOTE: Added data argument emmeans >vers 1.7.4. Object is unstable so feeding in the relevant arguments from model object directly. Note, we should think about df!
    mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights, ...)
    
    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if(is.null(by)){
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1]), upper = upper),
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
      
    } else{
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1]), upper = upper),
                              condition = mm_pi[,2], estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    }
    
    # Extract data
    data2 <- get_data_raw(model, mod, group, N, at = at, subset)
    
    mod_table$name <- factor(mod_table$name,
                             levels = mod_table$name,
                             labels = mod_table$name)
    
  } else{
    at2 <- list(mod = seq(min(data[,mod], na.rm = TRUE), max(data[,mod], na.rm = TRUE), length.out = 100))
    names(at2) <- mod
    grid <- emmeans::qdrg(formula =  stats::formula(model), data = data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k-1, at = c(at2, at))  # getting 100 points. Fixing this to make it more general
    mm <- emmeans::emmeans(grid, specs = mod, by = c(mod, by), weights = weights, df = df_mod)
    
    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if(is.null(by)){
      mod_table <- data.frame(moderator = mm_pi[,1],
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    } else{
      mod_table <- data.frame(moderator = mm_pi[,1],
                              condition = mm_pi[,2],
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    }
    
    # extract data
    data2 <- get_data_raw_cont(model, mod, group, N, by = by)
    
  }
  
  
  output <- list(mod_table = mod_table,
                 data = data2)
  
  class(output) <- c("orchard", "data.frame")
  
  return(output)
}




############# Key Sub-functions #############

#' @title pred_interval_esmeans
#' @description Function to get prediction intervals (credibility intervals) from \code{esmeans} objects (\pkg{metafor}).
#' @param model \code{rma.mv} object.
#' @param mm result from \code{emmeans::emmeans} object.
#' @param mod Moderator of interest.
#' @param ... other arguments passed to function.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export


pred_interval_esmeans <- function(model, mm, mod, ...){
  
  tmp <- summary(mm)
  tmp <- tmp[ , ]
  test.stat <- stats::qt(0.975, tmp$df[[1]])
  
  if(length(model$tau2) <= 1 | length(model$gamma2) <= 1){ # Note this should fix #46 but code is repetitive and needs to be cleaned up. Other issue is how this plays with different rma. objects. uni models will treat slots for gamma NULL and we need to deal with this. 
    sigmas <- sum(model$sigma2)
    taus   <- model$tau2
    gamma2 <- ifelse(is.null(model$gamma2), 0, model$gamma2)
    PI <- test.stat * base::sqrt(tmp$SE^2 + sigmas + taus + gamma2)
  } else {
    sigmas <- sum(model$sigma2)
    taus   <- model$tau2
    gammas <- model$gamma2
    w_tau <- model$g.levels.k
    w_gamma <- model$g.levels.k
    
    if(mod == "1"){
      tau <- weighted_var(taus, weights = w_tau)
      gamma <- weighted_var(gamma, weights = w_gamma)
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + tau + gamma)
      
    } else {
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + taus + gammas)
    }
  }
  
  tmp$lower.PI <- tmp$emmean - PI
  tmp$upper.PI <- tmp$emmean + PI
  
  # renaming "overall" to ""
  if(tmp[1,1] == "overall"){tmp[,1] <- "intrcpt"}
  
  return(tmp)
}

#' @title get_data_raw
#' @description Collects and builds the data used to fit the \code{rma.mv} or \code{rma} model in \pkg{metafor}.
#' @param model \code{rma.mv} object.
#' @param mod the moderator variable.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or whatever other grouping variable one wishes to present sample sizes.
#' @param N The name of the column in the data specifying the sample size, N. Defaults to \code{NULL}, so precision is plotted instead of sample size.
#' @param at List of moderators. If \code{at} is equal to \code{mod} then levels specified within \code{at} will be used to subset levels when \code{subset = TRUE}. Otherwise, it will marginalise over the moderators at the specified levels.
#' @param subset Whether or not to subset levels within the \code{mod} argument. Defaults to \code{FALSE}.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#' @examples \dontrun{
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#'  test <- get_data_raw(model, mod = "trait.type", group = "group_ID", at = list(trait.type = c("physiology", "morphology")))
#'  test2 <- get_data_raw(model, mod = "1", group = "group_ID")
#'
#'  data(english)
#'  # We need to calculate the effect sizes, in this case d
#'  english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"))
#'  model <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#'  test3 <-  get_data_raw(model, mod = "1", group = "StudyNo")}

get_data_raw <- function(model, mod, group, N = NULL, at = NULL, subset = TRUE){
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  
  # Extract the data from the model object
  data <- model$data 
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(at) & subset){
    # Find the at slot in list that pertains to the moderator and extract levels
    at_mod <- at[[mod]]
    position2 <- which(data[,mod] %in% at_mod)
    # Subset the data to only the levels in the moderator
    data <- data[position2,]
    yi <- model$yi[position2]
    vi <- model$vi[position2]
    type <- attr(model$yi, "measure")
  } else {
    # Extract effect sizes
    yi <- model$yi
    vi <- model$vi
    type <- attr(model$yi, "measure")
  }
  if(mod == "1"){
    moderator <- "Intrcpt"
  }else{
    # Get moderator
    moderator <- as.character(data[[mod]]) # Could default to base instead of tidy
    moderator <- firstup(moderator)
  }
  # Extract study grouping variable to calculate the
  stdy <- data[[group]] # Could default to base instead of tidy
  data_reorg <- data.frame(yi, vi, moderator, stdy, type)
  #names(data_reorg)[4] <- "stdy" # sometimes stdy gets replaced by group's names
  row.names(data_reorg) <- 1:nrow(data_reorg)
  
  if(is.null(N) == FALSE){
    data_reorg$N <- data[ ,N]
  }
  
  return(data_reorg)
}

#' @title get_data_raw_cont
#' @description Collects and builds the data used to fit the \code{rma.mv} or \code{rma} model in \pkg{metafor} when a continuous variable is fit within a model object.
#' @param model \code{rma.mv} object.
#' @param mod the moderator variable.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species or whatever other grouping variable one wishes to present sample sizes.
#' @param N  The name of the column in the data specifying the sample size, N. Defaults to \code{NULL} so that precision is plotted instead of sample size.
#' @param by Character name(s) of the 'condition' variables to use for grouping into separate tables.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export

#TODO what if there is no "by"

get_data_raw_cont <- function(model, mod, group, N = NULL, by){
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  
  # Extract the data from the model object
  data <- model$data 
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  # Extract effect sizes
  yi <- model$yi
  vi <- model$vi
  type <- attr(model$yi, "measure")
  # Get moderator
  moderator <- data[[mod]] # Could default to base instead of tidy
  #names(moderator) <  "moderator"
  if(is.null(by)){
    condition <- data[ , by]
  }else{
    condition <- data[[by]]
  }
  #names(condition) <  "condition"
  # Extract study grouping variable to calculate the
  stdy <- data[[group]] # Could default to base instead of tidy
  data_reorg <- data.frame(yi, vi, moderator, condition, stdy, type)
  # if(!is.na(names(data_reorg)[names(data_reorg) == by]) == TRUE) {  ## FAILING HERE
  #   names(data_reorg)[names(data_reorg) == by] <- "condition"
  # }
  #names(data_reorg)[5] <- "stdy" # sometimes stdy gets replaced by group's names
  row.names(data_reorg) <- 1:nrow(data_reorg)
  
  if(is.null(N) == FALSE){
    data_reorg$N <- data[ ,N]
  }
  
  return(data_reorg)
}

############# Helper-functions #############

#' @title firstup
#' @description Uppercase moderator names
#' @param x a character string
#' @param upper logical indicating if the first letter of the character string should be capitalized. Defaults to TRUE.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a character string with all combinations of the moderator level names with upper case first letters
#' @export

firstup <- function(x, upper = TRUE) {
  if(upper){
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  } else{ x }
}


#' @title print.orchard
#' @description Print method for class 'orchard'
#' @param x an R object of class orchard
#' @param ... Other arguments passed to print
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#'

print.orchard <- function(x, ...){
  return(print(x$mod_table))
}

#' @title weighted_var
#' @description Calculate weighted variance
#' @param x A vector of tau2s to be averaged
#' @param weights Weights, or sample sizes, used to average the variance
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a vector with a single weighted variance
#' @export
#'

weighted_var <- function(x, weights){
  weight_var <- sum(x * weights) / sum(weights)
  return(weight_var)
}


#' @title num_studies
#' @description Computes how many studies are in each level of categorical moderators of a \code{rma.mv} model object.
#' @param mod Character string describing the moderator of interest.
#' @param data Raw data from object of class "orchard"
#' @param group A character string specifying the column name of the study ID grouping variable.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a table with the number of studies in each level of all parameters within a \code{rma.mv} or \code{rma} object.
#' @export
#' @examples
#' \dontrun{data(fish)
#'warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list( ~1 | es_ID,~1 | group_ID), mods = ~experimental_design-1, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' num_studies(model$data, experimental_design, group_ID)
#' }

num_studies <- function(data, mod, group){
  
  # Summarize the number of studies within each level of moderator
  table <- data        %>%
    dplyr::group_by({{mod}}) %>%
    dplyr::summarise(stdy = length(unique({{group}})))
  
  table <- table[!is.na(table$moderator),]
  # Rename, and return
  colnames(table) <- c("Parameter", "Num_Studies")
  return(data.frame(table))
  
}

#' @title orchard_plot
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, or a results table of class \code{orchard}, it creates an orchard plot from mean effect size estimates for all levels of a given categorical moderator, and their corresponding confidence and prediction intervals.
#' @param object model object of class \code{rma.mv}, \code{rma}, or \code{orchard} table of model results.
#' @param mod the name of a moderator. Defaults to \code{"1"} for an intercept-only model. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding varaibles in 'by'. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param xlab The effect size measure label.
#' @param N The name of the column in the data specifying the sample size so that each effect size estimate is scaled to the sample size, N. Defaults to \code{NULL}, so that precision is used for scaling each raw effect size estimate instead of sample size.
#' @param alpha The level of transparency for effect sizes represented in the orchard plot.
#' @param angle The angle of y labels. The default is 90 degrees.
#' @param cb If \code{TRUE}, it uses 20 colour blind friendly colors.
#' @param k If \code{TRUE}, it displays k (number of effect sizes) on the plot.
#' @param g If \code{TRUE}, it displays g (number of grouping levels for each level of the moderator) on the plot.
#' @param transfm If set to \code{"tanh"}, a tanh transformation will be applied to effect sizes, converting Zr to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD).  \code{"invlogit"} can be used to convert lnRR to the inverse logit scale. \code{"percentr"} can convert to the percentage change scale when using response ratios and \code{"percent"} can convert to the percentage change scale of an log transformed effect size. Defaults to \code{"none"}.
#' @param condition.lab Label for the condition being marginalized over.
#' @param tree.order Order in which to plot the groups of the moderator when it is a categorical one. Should be a vector of equal length to number of groups in the categorical moderator, in the desired order (bottom to top, or left to right for flipped orchard plot)
#' @param trunk.size Size of the mean, or central point.
#' @param branch.size Size of the confidence intervals.
#' @param twig.size Size of the prediction intervals.
#' @param legend.pos Where to place the legend. To remove the legend, use \code{legend.pos = "none"}.
#' @param k.pos Where to put k (number of effect sizes) on the plot. Users can specify the exact position or they can use specify \code{"right"}, \code{"left"},  or \code{"none"}. Note that numeric values (0, 0.5, 1) can also be specified and this would give greater precision.
#' @param colour Colour of effect size shapes. By default, effect sizes are colored according to the \code{mod} argument. If \code{TRUE}, they are colored according to the grouping variable
#' @param fill If \code{TRUE}, effect sizes will be filled with colours. If \code{FALSE}, they will not be filled with colours.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param upper Logical, defaults to \code{TRUE}, indicating that the first letter of the character string for the moderator variable should be capitalized.
#' @param flip Logical, defaults to \code{TRUE}, indicating whether the plot should be flipped.
#' @return Orchard plot
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(eklof)
#' eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
#' m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
#' data=eklof)
#' # Add the unit level predictor
#' eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
#' # fit a MLMR - accounting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type-1,
#' random=list(~1|ExptID, ~1|Datapoint), data=eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")
#' orchard_plot(results, mod = "Grazer.type",
#' group = "ExptID", xlab = "log(Response ratio) (lnRR)")
#' # or
#' orchard_plot(eklof_MR, mod = "Grazer.type", group = "ExptID",
#' xlab = "log(Response ratio) (lnRR)")
#'
#' # Example 2
#' data(lim)
#' lim$vi<- 1/(lim$N - 3)
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article,
#' ~1|Datapoint), data=lim)
#' orchard_plot(lim_MR, mod = "Phylum", group = "Article",
#' xlab = "Correlation coefficient", transfm = "tanh", N = "N")
#' }
#' @export

orchard_plot <- function(object, mod = "1", group, xlab, N = NULL,
                         alpha = 0.5, angle = 90, cb = TRUE, k = TRUE, g = TRUE,
                         tree.order = NULL, trunk.size = 0.5, branch.size = 1.2, twig.size = 0.5,
                         transfm = c("none", "tanh", "invlogit", "percent", "percentr"), condition.lab = "Condition",
                         legend.pos = c("bottom.right", "bottom.left",
                                        "top.right", "top.left",
                                        "top.out", "bottom.out",
                                        "none"), # "none" - no legends
                         k.pos = c("right", "left", "none"),
                         colour = FALSE,
                         fill = TRUE,
                         weights = "prop", by = NULL, at = NULL, upper = TRUE, flip = TRUE)
{
  ## evaluate choices, if not specified it takes the first choice
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  
  if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){
    
    if(mod != "1"){
      results <-  mod_results(object, mod, group,  N,
                                       by = by, at = at, weights = weights, upper = upper)
    } else {
      results <-  mod_results(object, mod = "1", group,  N,
                                       by = by, at = at, weights = weights, upper = upper)
    }
  }
  
  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }
  
  mod_table <- results$mod_table
  
  data_trim <- results$data
  # making sure factor names match
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)
  
  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
  legend <- "Precision (1/SE)"
  
  #if tree.order isn't equal to NULL, and length of tree order does not match number of categories in categorical moderator, then stop function and throw an error
  if(!is.null(tree.order)&length(tree.order)!=nlevels(data_trim[,'moderator'])){
    stop("Length of 'tree.order' does not equal number of categories in moderator")
  }
  
  #if tree.order isn't equal to NULL but passes above check, then reorder mod table according to custom order if there is one
  if (!is.null(tree.order)){
    data_trim$moderator<-factor(data_trim$moderator, levels = tree.order, labels = tree.order)
    mod_table <- mod_table %>% dplyr::arrange(factor(name, levels = tree.order))
  }
  
  if(is.null(N) == FALSE){
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
    #latex2exp::TeX()
  }
  
  if(transfm == "tanh"){
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- Zr_to_r(mod_table[,cols])
    data_trim$yi <- Zr_to_r(data_trim$yi)
    label <- xlab
  }
  
  if(transfm == "invlogit"){
    
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- lapply(mod_table[,cols], function(x) metafor::transf.ilogit(x))
    data_trim$yi <- metafor::transf.ilogit(data_trim$yi)
    label <- xlab
  }
  
  if(transfm == "percentr"){
    
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- lapply(mod_table[,cols], function(x) (exp(x) - 1)*100)
    data_trim$yi <- (exp(data_trim$yi) - 1)*100
    label <- xlab
  } 
  
  
  
  if(transfm == "percent"){
    
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- lapply(mod_table[,cols], function(x) exp(x)*100)
    data_trim$yi <- (exp(data_trim$yi)*100)
    label <- xlab
  } else{
    label <- xlab
  }
  
  # Add in total effect sizes for each level
  mod_table$K <- as.vector(by(data_trim, data_trim[,"moderator"], function(x) length(x[,"yi"])))
  
  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[,2])
  
  # the number of groups in a moderator & data points
  group_no <- length(unique(mod_table[, "name"]))
  
  #data_no <- nrow(data)
  
  # colour blind friendly colours with grey
  cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  
  # setting fruit colour
  if(colour == TRUE){
    color <- as.factor(data_trim$stdy)
    color2 <- NULL
  }else{
    color <- data_trim$mod
    color2 <- mod_table$name
  }
  
  # whether we fill fruit or not
  if(fill == TRUE){
    fill <- color
  }else{
    fill <- NULL
  }
  
  # whether marginal
  if(names(mod_table)[2] == "condition"){
    
    # the number of levels in the condition
    condition_no <- length(unique(mod_table[, "condition"]))
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +
      
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
      # creating CI
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL),
                              size = branch.size, position = ggplot2::position_dodge2(width = 0.3)) +
      # drowning point estimate and PI
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR,  shape = as.factor(condition), fill = color2), size = trunk.size, position = ggplot2::position_dodge2(width = 0.3), linewidth = twig.size) +
      # this will only work for up to 5 different conditions
      # flipping things around (I guess we could do use the same geoms but the below is the original so we should not change)
      ggplot2::scale_shape_manual(values =  20 + (1:condition_no))  +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::labs(shape = condition.lab) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle))
    if(flip){
      plot <- plot + ggplot2::coord_flip()
    }
    
  } else {
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +
      
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
      # creating CI
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL), size = branch.size) +
      # drowning point estimate and PI
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name,  ymin = lowerPR, ymax = upperPR, fill = color2), size = trunk.size, linewidth = twig.size, shape = 21) +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle)) #+
    #ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
    if(flip){
      plot <- plot + ggplot2::coord_flip()
    }
  }
  
  # adding legend
  if(legend.pos == "bottom.right"){
    plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
  } else if ( legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
  } else if ( legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
  } else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
  } else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position="top")
  } else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position="bottom")
  } else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position="none")
  }
  
  # putting colors in
  if(cb == TRUE){
    plot <- plot +
      ggplot2::scale_fill_manual(values = cbpl) +
      ggplot2::scale_colour_manual(values = cbpl)
  }
  
  # putting k and g in
  if(k == TRUE && g == FALSE && k.pos == "right"){
    plot <- plot +
      ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                        label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = 3.5)
  } else if(k == TRUE && g == FALSE && k.pos == "left") {
    plot <- plot +  ggplot2::annotate('text', y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                      label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "left", size = 3.5)
  } else if (k == TRUE && g == TRUE && k.pos == "right"){
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                     label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, hjust = "right", size = 3.5)
  } else if (k == TRUE && g == TRUE && k.pos == "left"){
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate('text',  y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                     label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, hjust = "left", size = 3.5)
  } else if (k == TRUE && g == FALSE && k.pos%in%c('right','left','none')==FALSE) {
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, group_no,
                                                                 1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no]),
                                     parse = TRUE, size = 3.5)
  } else if (k == TRUE && g == TRUE && k.pos%in%c('right','left','none')==FALSE) {
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, group_no,
                                                                 1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no],
                                                                                          "~", "(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, size = 3.5)
  }
  return(plot)
}


#' @title Zr_to_r
#' @description Converts Zr back to r (Pearson's correlation coefficient)
#' @param df data frame of results of class 'orchard'
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals with estimates converted back to r
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export

Zr_to_r <- function(df){
  return(sapply(df, tanh))
}

#' @title bubble_plot
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, or a results table of class \code{orchard}, the \code{bubble_plot} function creates a bubble plot from slope estimates. In cases when a model includes interaction terms, this function creates panels of bubble plots.
#' @param object Model object of class \code{rma}, \code{rma.mv}, or \code{orchard} table of model results
#' @param mod The name of a continuous moderator, to be plotted on the x-axis of the bubble plot.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding variables in \code{by}. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights How to marginalize categorical variables; used when one wants marginalised means. The default is \code{weights = "prop"}, which weights means for moderator levels based on their proportional representation in the data. For example, if \code{"sex"} is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal when, for example, males and females are typically roughly equally prevalent in a population. In cases such as these, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param transfm If set to \code{"tanh"}, a tanh transformation will be applied to effect sizes, converting Zr to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD).  \code{"invlogit"} can be used to convert lnRR to the inverse logit scale. \code{"percentr"} can convert to the percentage change scale when using response ratios and \code{"percent"} can convert to the percentage change scale of an log transformed effect size. Defaults to \code{"none"}.
#' @param xlab Moderator label.
#' @param ylab Effect size measure label.
#' @param k.pos The position of effect size number, k.
#' @param N The vector of sample size which an effect size is based on. Defaults to precision (the inverse of sampling standard error).
#' @param alpha The level of transparency for pieces of fruit (effect size).
#' @param cb If \code{TRUE}, it uses a colourblind-friendly palette of 20 colours (do not make this \code{TRUE}, when colour = \code{TRUE}).
#' @param k If \code{TRUE}, it displays k (number of effect sizes) on the plot.
#' @param g If \code{TRUE}, it displays g (number of grouping levels for each level of the moderator) on the plot.
#' @param est.lwd Size of the point estimate.
#' @param ci.lwd Size of the confidence interval.
#' @param pi.lwd Size of the prediction interval.
#' @param est.col Colour of the point estimate.
#' @param ci.col Colour of the confidence interval.
#' @param pi.col Colour of the prediction interval.
#' @param condition.nrow Number of rows to plot condition variable.
#' @param legend.pos Where to place the legend, or not to include a legend ("none").
#' @param cond_levels Order of the levels of the condition variable in the order to plot. Defaults to NULL.
#'
#' @return Orchard plot
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(lim)
#' lim[, "year"] <- as.numeric(lim$year)
#' lim$vi<- 1/(lim$N - 3)
#' model<-metafor::rma.mv(yi=yi, V=vi, mods= ~Environment*year,
#' random=list(~1|Article,~1|Datapoint), data=na.omit(lim))
#' test <- orchaRd::mod_results(model, mod = "year", group = "Article", data = lim, weights = "prop", by = "Environment")
#' orchaRd::bubble_plot(test, mod = "year", data = lim, group = "Article",legend.pos = "top.left")
#' # Or just using model directly
#' orchaRd::bubble_plot(model, mod = "year", legend.pos = "top.left", data = lim, group = "Article", weights = "prop", by = "Environment")
#'
#' }
#' @export

# TODO - make poly works for bubble???
# TODO - write to https://github.com/rvlenth/emmeans/issues (missing combinations or interaction not allowed)

bubble_plot <- function(object, mod, group = NULL, xlab = "Moderator", ylab = "Effect size", N = "none",
                        alpha = 0.5, cb = TRUE, k = TRUE, g = FALSE, transfm = c("none", "tanh", "invlogit", "percent", "percentr"), 
                        est.lwd = 1, ci.lwd = 0.5, pi.lwd = 0.5,
                        est.col = "black", ci.col = "black", pi.col = "black",
                        legend.pos = c("top.left", "top.right",
                                       "bottom.right", "bottom.left",
                                       "top.out", "bottom.out",
                                       "none"),
                        k.pos = c("top.right", "top.left",
                                  "bottom.right", "bottom.left",
                                  "none"),
                        condition.nrow = 2,
                        #condition.lab = "Condition",
                        weights = "prop", by = NULL, at = NULL, cond_levels = NULL)
{
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  #facet <- match.arg(NULL, choices = facet)
  
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?bubble_plot")
  }
  
  if(is.numeric(by)){
    k = FALSE
    g = FALSE
  }
  
  
  if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){
    
    if(mod != "1"){
      results <-  mod_results(object, mod, group,
                                       by = by, at = at, weights = weights)
    } else {
      results <-  mod_results(object, mod = "1", group,
                                       by = by, at = at, weights = weights)
    }
  }
  
  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }
  
  mod_table <- results$mod_table
  
  data_trim <- results$data
  
  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
  legend <- "Precision (1/SE)"
  
  if(any(N != "none")){
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
  }
  
  if(transfm == "tanh"){
    cols <- which(colnames(mod_table) %in% c("condition", "moderator"))
    mod_table[,-cols] <- Zr_to_r(mod_table[,-cols])
    data_trim$yi <- Zr_to_r(data_trim$yi)
    label <- xlab
  }
  
  if(transfm == "invlogit"){
    
    cols <- which(colnames(mod_table) %in% c("condition", "moderator"))
    mod_table[,-cols] <- lapply(mod_table[,-cols], function(x) metafor::transf.ilogit(x))
    data_trim$yi <- metafor::transf.ilogit(data_trim$yi)
    label <- xlab
  }
  
  if(transfm == "percentr"){
    
    cols <- which(colnames(mod_table) %in% c("condition", "moderator"))
    mod_table[,-cols] <- lapply(mod_table[,-cols], function(x) (exp(x) - 1)*100)
    data_trim$yi <- (exp(data_trim$yi) - 1)*100
    label <- xlab
  } 
  
  if(transfm == "percent"){
    
    cols <- which(colnames(mod_table) %in% c("condition", "moderator"))
    mod_table[,-cols] <- lapply(mod_table[,-cols], function(x) exp(x)*100)
    data_trim$yi <- (exp(data_trim$yi)*100)
    label <- xlab
  } else{
    label <- xlab
  }
  
  if(is.null(data_trim$condition) == TRUE){
    
    # the number of effect sizes
    effect_num <- nrow(data_trim)
    
    # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
    group_num <- length(unique(data_trim$stdy))
    
    dat_text <- data.frame(K = effect_num, G = group_num)
    
  }else{
    
    if(!is.null(cond_levels)){
      data_trim$condition <- factor(data_trim$condition, levels = cond_levels, labels = cond_levels, ordered = TRUE)
    } else {
      # making sure factor names match
      data_trim$condition <- factor(data_trim$condition, levels = mod_table$condition, labels = mod_table$condition)
    }
    
    effect_num <- as.vector(by(data_trim, data_trim[,"condition"], function(x) base::length(x[,"yi"])))
    
    # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
    #group_num <- c(2,4)
    group_num <- as.vector(by(data_trim, data_trim[,"condition"], function(x) base::length(base::unique(x[,"stdy"]))))
    
    
    dat_text <- data.frame(K = effect_num, G = group_num, condition = as.vector(base::levels(data_trim$condition)))
  }
  # the number of groups in a moderator & data points
  #group_no <- length(unique(mod_table[, "name"]))
  
  #data_no <- nrow(data)
  
  # # colour blind friendly colours with grey
  # cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  
  if(is.null(data_trim$condition) == TRUE){
    plot <-ggplot2::ggplot() +
      # putting bubbles
      ggplot2::geom_point(data = data_trim, ggplot2::aes(x = moderator, y = yi, size = scale), shape = 21, alpha = alpha, fill = "grey90" ) +
      # prediction interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x, se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
      # confidence interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x, se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x, se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
      # main line
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
      #facet_grid(rows = vars(condition)) +
      ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
      ggplot2::guides(fill = "none", colour = "none") +
      # themes
      ggplot2::theme_bw() +
      #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
      ggplot2::theme(legend.direction="horizontal") +
      #theme(legend.background = element_rect(fill = "white", colour = "black")) +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))
  } else if(is.character(data_trim$condition) == TRUE || is.factor(data_trim$condition) == TRUE){
    
    plot <-ggplot2::ggplot() +
      # putting bubbles
      ggplot2::geom_point(data = data_trim, ggplot2::aes(x = moderator, y = yi, size = scale, fill = condition), shape = 21, alpha = alpha) +
      # prediction interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x,se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
      # confidence interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x,se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x,se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
      # main line
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
      ggplot2::facet_wrap(ggplot2::vars(condition), nrow = condition.nrow) +
      ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
      ggplot2::guides(fill = "none", colour = "none") +
      # themses
      ggplot2::theme_bw() +
      #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
      ggplot2::theme(legend.direction="horizontal") +
      #theme(legend.background = element_rect(fill = "white", colour = "black")) +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))
    
    # if(facet == "rows"){
    #   plot <- plot + facet_grid(rows = vars(condition))
    # } else{
    #   plot <- plot + facet_grid(cols = vars(condition))
    # }
    
    
  } else{
    plot <-ggplot2::ggplot() +
      # putting bubbles
      #geom_point(data = data, aes(x = moderator, y = yi, size = scale), shape = 21, alpha = alpha, fill = "grey90" ) +
      # prediction interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x,se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
      # confidence interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x,se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x,se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
      # main line
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
      ggplot2::facet_wrap(ggplot2::vars(condition), nrow = condition.nrow) +
      ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
      ggplot2::guides(fill = "none", colour = "none") +
      # themses
      ggplot2::theme_bw() # +
    #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
    # theme(legend.direction="horizontal") +
    # #theme(legend.background = element_rect(fill = "white", colour = "black")) +
    # theme(legend.background = element_blank()) +
    # theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))
  }
  
  # adding legend
  if(legend.pos == "bottom.right"){
    plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
  } else if ( legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
  } else if ( legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
  } else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
  } else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position="top")
  } else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position="bottom")
  } else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position="none")
  }
  
  # putting k and g in
  # c("top.right", "top.left", "bottom.right", "bottom.left","none")
  if(k == TRUE && g == FALSE && k.pos == "top.right"){
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = Inf),
                         label =  paste("italic(k)==", dat_text$K),
                         parse = TRUE,
                         hjust   = 2,
                         vjust   = 2.5
      )
    
  } else if(k == TRUE && g == FALSE && k.pos == "top.left") {
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = Inf),
                         label =  paste("italic(k)==", dat_text$K),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = 2.5
      )
  } else if(k == TRUE && g == FALSE && k.pos == "bottom.right") {
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = -Inf),
                         label =  paste("italic(k)==", rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = 2,
                         vjust   = -1.5
      )
  } else if (k == TRUE && g == FALSE && k.pos == "bottom.left"){
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = -Inf),
                         label =  paste("italic(k)==", dat_text$K),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = -1.5
      )
    # below get g ----
    
  } else if (k == TRUE && g == TRUE && k.pos == "top.right"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = Inf),
                         label =  paste("italic(k)==", dat_text$K,
                                        "~","(", dat_text$G, ")"),
                         parse = TRUE,
                         hjust   = 1.5,
                         vjust   = 2)
    
  } else if (k == TRUE && g == TRUE && k.pos == "top.left"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = Inf),
                         label =  paste("italic(k)==", dat_text$K,
                                        "~","(", dat_text$G, ")"),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = 2)
  } else if (k == TRUE && g == TRUE && k.pos == "bottom.right"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = -Inf),
                         label =  paste("italic(k)==", dat_text$K,
                                        "~","(", dat_text$G, ")"),
                         parse = TRUE,
                         hjust   = 1.5,
                         vjust   = -0.5)
  } else if (k == TRUE && g == TRUE && k.pos == "bottom.left"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = -Inf),
                         label =  paste("italic(k)==", dat_text$K,
                                        "~","(", dat_text$G, ")"),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = -0.5)
  }
  
  # # putting colors in
  # if(cb == TRUE){
  #   plot <- plot +
  #     ggplot2::scale_fill_manual(values=cbpl) +
  #     ggplot2::scale_colour_manual(values=cbpl)
  # }
  
  return(plot)
}


#' @title submerge
#' @description Merge two model results tables (orchard objects).
#' @param object1  object of class \code{orchard}.
#' @param object2  object of class \code{orchard}.
#' @param ... Other arguments passed to submerge.
#' @param mix If \code{TRUE}, it will add the number to the moderator name.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame.
#' @export
#'
submerge <- function(object1, object2, ..., mix = FALSE){
  orchard_list <- list(object1, object2, ...)
  
  len <- length(orchard_list)
  # merging tables
  tables <- lapply(orchard_list, function(x) x$mod_table)
  tables <- do.call("rbind", tables)
  
  # merging data
  ## checking moderator names are the same or not
  datas <- lapply(orchard_list, function(x) x$data)
  datas <- do.call("rbind", datas)
  
  # renaming
  if(mix == TRUE){
    names <- lapply(orchard_list, function(x) x$data$moderator)
    names <- as.vector(unlist(mapply(function(x, y) paste0(x, y), x = names, y = 1:len)))
    datas$moderator <- factor(names)
    tables$name <- levels(factor(names))
  }
  
  model_results <- list(mod_table = tables, data = datas)
  
  class(model_results) <- "orchard"
  
  return(model_results)
  
}
