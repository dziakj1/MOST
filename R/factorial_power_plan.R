#' @title  sample size, power and effect size calculations for a factorial or fractional
#' factorial experiment
#'
#' @description There are three ways to use this function:
#'
#' \enumerate{
#' \item Estimate power available from a given sample size and a given effect size.
#' \item Estimate sample size needed for a given power and a given effect size.
#' \item Estimate effect size detectable from a given power at a given sample size.
#' }
#'
#'That is, there are three main pieces of information: power, sample size, and effect size.
#'The user provides two of them, and this function calculates the third.
#'
#' @param alpha Two sided Type I error level for the test to be performed(default=0.05).
#' @param assignment One of three options: (default=unclustered)
#'                \enumerate{
#'                   \item \dQuote{independent} or equivalently \dQuote{unclustered}
#'                   \item \dQuote{within} or equivalently \dQuote{within_clusters}
#'                   \item \dQuote{between} or equivalently \dQuote{between_clusters}
#'                }
#'        Clusters in this context are preexisting units within which responses may
#'        be dependent (e.g., clinics or schools).  A within-cluster experiment involves
#'        randomizing individual members, while a between-cluster experiment involves
#'        randomizing clusters as whole units (see Dziak, Nahum-Shani, and Collins, 2012)
#'        <DOI:10.1037/a0026972>
#' @param change_score_icc The intraclass correlation of the change scores (posttest
#'        minus pretest). Relevant only if assignment is between clusters and there is a pretest.
#' @param cluster_size The mean number of members in each cluster. Relevant only
#'        if assignment is between clusters or within clusters.
#' @param cluster_size_sd Relevant only if assignment is between clusters. The standard deviation
#'                        of the number of members in each cluster (the default is 0
#'                        which means that the clusters are expected to be of equal size).
#' @param d_main Effect size measure: standardized mean difference raw_main/sigma_y.
#' @param effect_size_ratio Effect size measure: signal to noise ratio raw_coef^2/sigma_y^2.
#' @param icc Relevant only if assignment is between clusters or within clusters. The intraclass
#'            correlation of the variable of interest in the absence of treatment.
#' @param model_order The highest order term to be included in the regression model in the planned analysis
#'                    (1=main effects, 2=two-way interactions, 3=three-way interactions, etc.); must be >= 1 and
#'                    <= nfactors (default=1).
#' @param nclusters The total number of clusters available (for between clusters or within clusters assignment).
#' @param nfactors The number of factors (independent variables) in the planned experiment(default=1).
#' @param ntotal The total sample size available (for unclustered assignment. For clustered assignment, use
#'               \dQuote{cluster_size} and \dQuote{nclusters.}
#' @param power If specified: The desired power of the test. If returned in the output list: The expected power of the test.
#' @param pre_post_corr Relevant only if there is a pretest. The correlation between the pretest and the posttest.
#' @param pretest One of three options:
#'                \enumerate{
  #'                \item \dQuote{no} or \dQuote{none} for no pretest.
  #'                \item \dQuote{covariate} for pretest to be entered as a covariate in the model.
#'                  \item \dQuote{repeated} for pretest to be considered as a repeated measure.
#'                }
#'                The option \dQuote{yes} is also allowed and is interpreted as \dQuote{repeated.}
#'                The option \dQuote{covariate} is not allowed if assignment is between clusters.  This
#'                is because predicting power for covariate-adjusted cluster-level randomization is
#'                somewhat complicated, although it can be approximated in practice by using the
#'                formula for the repeated-measures cluster-level randomization (see simulations in Dziak,
#'                Nahum-Shani, and Collins, 2012).
#' @param raw_coef Effect size measure: unstandardized effect-coded regression coefficient.
#' @param raw_main Effect size measure: unstandardized mean difference.
#' @param sigma_y The assumed standard deviation of the response variable after treatment, within each
#'                treatment condition (i.e., adjusting for treatment but not adjusting for post-test).
#'                This statement must be used if the effect size argument used is either \dQuote{raw_main} or
#'                \dQuote{raw_coef}.
#' @param std_coef Effect size measure: standardized effect-coded regression coefficient raw_coef/sigma_y.
#' @importFrom stats pf qf
#'
#' @return A list with power, sample size and effect size.
#'
#' @export FactorialPowerPlan
#' @examples
#' FactorialPowerPlan(assignment="independent",
#'                    model_order=2,
#'                    nfactors=5,
#'                    ntotal=300,
#'                    raw_main=3,
#'                    sigma_y=10)
#'
#' FactorialPowerPlan(assignment="independent",
#'                    model_order=2,
#'                    nfactors=5,
#'                    ntotal=300,
#'                    pre_post_corr=.6,
#'                    pretest="covariate",
#'                    raw_main=3,
#'                    sigma_y=10)
#'

FactorialPowerPlan <- function( alpha = .05,
                                assignment = "unclustered",
                                change_score_icc = NULL,
                                cluster_size = NULL,
                                cluster_size_sd = NULL,
                                d_main = NULL,
                                effect_size_ratio = NULL,
                                icc = NULL,
                                model_order = 1,
                                nclusters = NULL,
                                nfactors = 1,
                                ntotal = NULL,
                                power = NULL,
                                pre_post_corr = NULL,
                                pretest = "none",
                                raw_coef = NULL,
                                raw_main = NULL,
                                sigma_y = NULL,
                                std_coef = NULL) {

  #Example 7
  # alpha <- .05;
  # assignment<-"within";
  # change_score_icc <- NULL;
  # cluster_size <- 10;
  # cluster_size_sd <- NULL;
  # d_main <- NULL;
  # effect_size_ratio <- NULL;
  # icc <- .1;
  # model_order<-2;
  # nclusters <- 50;
  # nfactors<-5;
  # ntotal<-NULL;
  # power <- .80;
  # pre_post_corr <- NULL;
  # pretest <- "none";
  # raw_coef <- NULL;
  # raw_main<- NULL;
  # sigma_y<- 10;
  # std_coef <- NULL;

  if(!is.null(change_score_icc)) has_specified_change_score_icc <- 1 else has_specified_change_score_icc <- 0;
  if(!is.null(cluster_size)) has_specified_cluster_size <- 1 else has_specified_cluster_size <- 0;
  if(!is.null(cluster_size_sd)) has_specified_cluster_size_sd <- 1 else has_specified_cluster_size_sd <- 0;
  if(!is.null(d_main)) has_specified_d_main <- 1 else has_specified_d_main <- 0;
  if(!is.null(effect_size_ratio)) has_specified_effect_size_ratio <- 1 else has_specified_effect_size_ratio <- 0;
  if(!is.null(icc)) has_specified_icc <- 1 else has_specified_icc <- 0;
  if(!is.null(nclusters)) has_specified_nclusters <- 1 else has_specified_nclusters <- 0;
  if(!is.null(ntotal)) has_specified_ntotal <- 1 else has_specified_ntotal <- 0;
  if(!is.null(power)) has_specified_power <- 1 else has_specified_power <- 0;
  if(!is.null(pre_post_corr)) has_specified_pre_post_corr <- 1 else has_specified_pre_post_corr <- 0;
  if(!is.null(pretest)) has_specified_pretest <- 1 else has_specified_pretest<- 0;
  if(!is.null(raw_coef)) has_specified_raw_coef <- 1 else has_specified_raw_coef<- 0;
  if(!is.null(raw_main)) has_specified_raw_main <- 1 else has_specified_raw_main<- 0;
  if(!is.null(sigma_y)) has_specified_sigma_y <- 1 else has_specified_sigma_y<- 0;
  if(!is.null(std_coef)) has_specified_std_coef <- 1 else has_specified_std_coef<- 0;

  num_effect_sizes_specified <- has_specified_d_main +
    has_specified_effect_size_ratio +
    has_specified_raw_coef +
    has_specified_raw_main +
    has_specified_std_coef;
  has_specified_some_effect_size <- 1*(num_effect_sizes_specified>0);
  if(num_effect_sizes_specified > 1) stop("Error: Only one effect size measure should be used.");
  num_sample_sizes_specified <- has_specified_ntotal + has_specified_nclusters;
  has_specified_some_sample_size <- 1*(num_sample_sizes_specified>0);
  if(num_sample_sizes_specified > 1) stop("Error: Only one sample size measure should be used.");
  to_find <- "????????????????";
  if(has_specified_some_effect_size == 0) to_find <- "effect_size";
  if(has_specified_power == 0) to_find <-"power";
  if(has_specified_some_sample_size == 0) to_find <- "sample_size";
  num_big_three_specified <- has_specified_some_effect_size +
    has_specified_power +
    has_specified_some_sample_size;
  if(num_big_three_specified != 2) {
    stop("Error: Exactly two of the following should be specified: \n
         an effect size, a sample size, and a power.  The third will be \n
         the result to be calculated.");
  }
  if( assignment == "unclustered" |
      assignment == "independent" |
      assignment == "" |
      assignment == " " ) {
    assignment <- "unclustered";
  } else if((assignment == "betweenclusters")|
            (assignment == "between_clusters")|
            (assignment == "between")) {
    assignment <- "between_clusters";
  } else if((assignment == "withinclusters")|
            (assignment == "within_clusters")|
            (assignment == "within")) {
    assignment <- "within_clusters";
  } else {
    stop("Error: Could not understand the specified assignment type.");
  }

  if(assignment == "unclustered"){
    if(has_specified_cluster_size == 1) {
      stop("Error:  If the subjects are not clustered then please do \n
           not specify a cluster size");
    }
    if(has_specified_cluster_size_sd == 1) {
      stop("Error:  If the subjects are not clustered then please do \n
           not specify a cluster size standard deviation");
    }
    if(has_specified_nclusters == 1) {
      stop("Error:  If the subjects are not clustered then please do \n
           not specify a number of clusters.");
    }
    }
  if(assignment == "within_clusters"){
    if(has_specified_cluster_size_sd == 1) {
      stop("Error:  If assignment is within clusters then please do \n
           not specify a cluster size standard deviation.");
    }
    }

  if(  (pretest == "no_pretest") |
       (pretest == "none") |
       (pretest == "no") |
       (pretest == "") |
       (pretest == " ") ) {
    pretest <- "none";
  } else if ((pretest == "covariate")|
             (pretest == "ancova")) {
    pretest <- "covariate";
  } else if((pretest == "yes") |
            (pretest == "repeated") |
            (pretest == "repeated_measure") |
            (pretest == "repeated_measures")) {
    pretest <- "repeated_measure";
  } else {
    stop("Macro error: Could not understand the specified pretest type.");
  }

  if(has_specified_sigma_y == 0) {
    if(to_find == "effect_size" |
       has_specified_d_main |
       has_specified_std_coef |
       has_specified_effect_size_ratio) {
      sigma_y <- 1.0;
    } else {
      stop("Macro error: Please specify a non-blank sigma_y.");
    }
  }

  if((model_order<1)|(model_order>10)) {
    stop("Please choose a different value for model_order.");
  }
  if((nfactors<1)|(nfactors>99)) {
    stop("Please choose a different value for nfactors.");
  }
  if(model_order>nfactors) {
    stop("Macro error: The specified model order is too high for the specified \n
         number of factors.");
  }
  if(pretest != "none") {
    if(has_specified_pre_post_corr == 1) {
      if((pre_post_corr<0)|(pre_post_corr>=1)) {
        stop("Macro error: The pretest-posttest correlation should be between 0 and 1.");
      }
    } else {
      stop("Macro error: The pretest-posttest correlation needs to be specified.");
    }
  }
  if(assignment=="unclustered") {
    if(pretest=="none") sigma_effective <- sigma_y;
    if(pretest=="covariate"){
      sigma_effective <- sigma_y*sqrt(1-pre_post_corr^2);
    }
  }

  if((assignment=="between_clusters"|assignment=="within_clusters")) {
    if(pretest=="none") {
      sigma_effective <- sigma_y;
    } else if(pretest=="covariate") {
      if(assignment=="between_clusters") {
        stop("Macro error: Not able to calculate power for this model setup.");
      } else {
        sigma_effective <- sigma_y*sqrt(1-pre_post_corr^2);
      }
    } else if(pretest=="repeated_measure") {
      sigma_effective <- sigma_y;
    }
  }
  if(has_specified_d_main) the_coef <- .5*d_main*sigma_y;
  if(has_specified_std_coef) the_coef <- std_coef*sigma_y;
  if(has_specified_effect_size_ratio) the_coef <- sqrt(effect_size_ratio);
  if(has_specified_raw_coef) the_coef <- raw_coef;
  if(has_specified_raw_main) the_coef <- .5*raw_main;
  if(has_specified_d_main | has_specified_std_coef |
     has_specified_effect_size_ratio | has_specified_raw_coef |
     has_specified_d_main) the_coef <- abs(the_coef);

  nregcoefs <- 1; #start with intercept
  for(k in 1:model_order){
    nregcoefs <- nregcoefs + choose(nfactors,k);
    # add in number of k-way effects */
  }
  if((pretest=="covariate")& (assignment!="betweenclusters")) {
    nregcoefs <- nregcoefs + 1;
    # add one for the covariate, if there is one and if
    #it is at a measurement level at which the extra df
    #appreciably affects the power. */
  }

  success <- 0;
  raw_regression_coefficient <- NA;
  print("------------------------------------------------------------");
  print("FactorialPowerPlan Macro");
  print("The Methodology Center");
  print("(c) 2012 Pennsylvania State University");
  print("------------------------------------------------------------");
  #if(MacroError == 1 ){
  #CALL HandleDataErrors(DataErrorCode);
  #}
  print("Assumptions:");
  print(paste("There are ", as.character(nfactors)," dichotomous factors.",sep=""));
  if(assignment == "unclustered") print("There is independent random assignment.");
  if(assignment == "within_clusters") print("There is random assignment of individuals for each cluster (within-clusters effects).");
  if(assignment == "between_clusters") print("There is random assignment of clusters (between-clusters effects).");
  if(pretest == "covariate") print("There is a pretest modeled as a covariate.");
  if(pretest == "repeated_measure") print("There is a pretest modeled as a repeated measure.");
  if(model_order == 1) print("Analysis will be based on main effects only.");
  if(model_order == 2) print("Analysis will be based on main effects and 2-way interactions.");
  if(model_order == 3) print("Analysis will be based on main effects, 2-way, and 3-way interactions.");
  if(has_specified_power) {
    print(sprintf("Desired power: %9.3f",  power));
  }
  print(sprintf( "Two-sided alpha:%9.2f", alpha));
  if(has_specified_cluster_size){
    print(sprintf( "Cluster size:%9.2f", cluster_size));
  }
  if(has_specified_cluster_size_sd){
    print(sprintf( "Cluster size standard deviation:%9.2f", cluster_size_sd));
  }
  if(has_specified_nclusters){
    print(paste("Number of clusters: ", as.character(nclusters), sep=""));
  }
  if(has_specified_ntotal){
    print(paste("Total number of participants: ", as.character(ntotal), sep=""));
  }
  if(has_specified_d_main){
    print(sprintf("Effect size as standardized mean difference:%9.2f", d_main));
  }
  if(has_specified_effect_size_ratio){
    print(sprintf("Effect size as variance ratio:%6.2f", effect_size_ratio));
  }
  if(has_specified_raw_coef){
    print(sprintf("Effect size as unstandardized regression coefficient:%9.2f", raw_coef));
  }
  if(has_specified_raw_main){
    print(sprintf("Effect size as unstandardized difference in means:%9.2f", raw_main));
  }
  if(has_specified_std_coef){
    print(sprintf("Effect size as standardized regression coefficient:%9.2f", std_coef));
  }
  if(has_specified_icc){
    print(sprintf( "Intraclass correlation of response variable:%9.3f", icc));
  }
  if(has_specified_change_score_icc){
    print(sprintf("Intraclass correlation of change scores:%9.3f", change_score_icc));
  }
  if(has_specified_pre_post_corr){
    print(sprintf("Pretest-posttest correlation:%9.3f", pre_post_corr));
  }
  if(has_specified_sigma_y){
    print(sprintf("Assumed standard deviation for the response variable is %9.2f", sigma_y));
  }
  if(to_find == "effect_size"){
    print("Attempting to calculate the estimated detectable effect size.");
  }
  if(to_find == "power"){
    print("Attempting to calculate the estimated power.");
  }
  if(to_find == "sample_size"){
    print("Attempting to calculate the estimated required sample size.");
  }
  print("------------------------------------------------------------");
  print("Results:");
  # One-level (subject) model
  if(assignment == "unclustered" & ((pretest == "none")|(pretest == "covariate"))){
    if(to_find == "power") {
      power <- CalculatePowerUnclustered(alpha,
                                         nfactors,
                                         nregcoefs,
                                         ntotal,
                                         sigma_effective,
                                         the_coef,
                                         0);
      success <- 1;
    }
    if(to_find == "sample_size") {
      ntotal <- CalculateNUnclustered( alpha,
                                       nfactors,
                                       nregcoefs,
                                       power,
                                       sigma_effective,
                                       the_coef );
      success <- 1;
    }
    if(to_find == "effect_size") {
      the_coef <- CalculateDetectableUnclustered( alpha,
                                                  nfactors,
                                                  nregcoefs,
                                                  ntotal,
                                                  power,
                                                  sigma_effective,
                                                  sigma_y);
      success <- 1;
    }
  }
  # Two-level (cluster, subject) models
  if(((assignment == "between_clusters") & (pretest == "none"))
     |((assignment == "within_clusters")& ((pretest == "none")|(pretest == "covariate")))) {
    if(has_specified_cluster_size == 0) {
      stop("Macro error: The cluster size must be specified.");
    }
    if(has_specified_cluster_size_sd == 0){
      cluster_size_sd <- 0;
    }
    cluster_size_cv <- cluster_size_sd/cluster_size;
    if(has_specified_icc == 0) {
      stop("Macro error: icc must be specified.");
    }
    if(pretest == "none") {
      sigma_effective <- sigma_y;
    } else {
      sigma_effective <- sigma_y*sqrt(1-pre_post_corr^2);
    }
    if(to_find == "power") {
      power = CalculatePowerTwoLevel( alpha,
                                      assignment,
                                      cluster_size,
                                      cluster_size_sd,
                                      icc,
                                      nclusters,
                                      nfactors,
                                      nregcoefs,
                                      sigma_effective,
                                      the_coef,
                                      0);
      success <- 1;
    }
    if(to_find == "sample_size") {
      nclusters <- CalculateNTwoLevel(alpha,
                                      assignment,
                                      cluster_size,
                                      cluster_size_sd,
                                      icc,
                                      nfactors,
                                      nregcoefs,
                                      power,
                                      sigma_effective,
                                      the_coef);
      ntotal <- nclusters*cluster_size;
      success <- 1;
    }
    if(to_find == "effect_size") {
      the_coef <- CalculateDetectableTwoLevel( alpha,
                                               assignment,
                                               cluster_size,
                                               cluster_size_sd,
                                               icc,
                                               nclusters,
                                               nfactors,
                                               nregcoefs,
                                               power,
                                               sigma_effective,
                                               sigma_y);
      success <- 1;
    }
  }
  #Repeated measures models: either two-level (subject, time) models or three-level (cluster, subject, time) models
  if(((assignment=="between_clusters")|(assignment=="within_clusters")|(assignment=="unclustered")) &
     (pretest=="repeated_measure")) {
    if(assignment=="unclustered") {
      cluster_size <- 1;
      cluster_size_sd <- 0;
      cluster_size_cv <- 0;
      icc <- 0;
      change_score_icc <- 0;
      nclusters <- ntotal;
    }
    if(assignment=="within_clusters") {
      if(has_specified_cluster_size == 0) {
        stop("Macro error: The cluster size must be specified.");
      }
      if(has_specified_cluster_size_sd == 0) {
        cluster_size_sd <- 0;
      }
      cluster_size_cv <- cluster_size_sd/cluster_size;
      if(has_specified_icc == 0) {
        stop("Macro error: icc must be specified.");
      }
      if(has_specified_change_score_icc == 1) {
        stop("Macro error: change_score_icc should not be specified for this model.");
      }
      change_score_icc <- 0;
      if(has_specified_pre_post_corr == 0){
        stop("Macro error: pre_post_corr must be specified.");
      }
    }
    if(assignment == "between_clusters"){
      if(has_specified_cluster_size == 0) {
        stop("Macro error: The cluster size must be specified.");
      }
      if(has_specified_cluster_size_sd == 0) {
        cluster_size_sd <- 0;
      }
      cluster_size_cv <- cluster_size_sd/cluster_size;
      if(has_specified_icc == 0) {
        stop("Macro error: icc must be specified.");
      }
      if(has_specified_change_score_icc == 0) {
        stop("Macro error: change_score_icc must be specified.");
      }
      if(has_specified_pre_post_corr == 0) {
        stop("Macro error: pre_post_corr must be specified.");
      }
    }
    if(to_find == "power") {
      power = CalculatePowerThreeLevel( alpha,
                                        assignment,
                                        change_score_icc,
                                        cluster_size,
                                        cluster_size_sd,
                                        icc,
                                        nclusters,
                                        nfactors,
                                        nregcoefs,
                                        pre_post_corr,
                                        sigma_y,
                                        the_coef,
                                        0);
      success <- 1;
    }
    if(to_find == "sample_size") {
      nclusters = CalculateNThreeLevel( alpha,
                                        assignment,
                                        change_score_icc,
                                        cluster_size,
                                        cluster_size_sd,
                                        icc,
                                        nfactors,
                                        nregcoefs,
                                        power,
                                        pre_post_corr,
                                        sigma_y,
                                        the_coef);
      success <- 1;
      if(assignment == "unclustered") {
        ntotal <- nclusters;
      } else {
        ntotal <- nclusters*cluster_size;
      }
    }
    if(to_find == "effect_size") {
      the_coef = CalculateDetectableThreeLevel( alpha,
                                                assignment,
                                                change_score_icc,
                                                cluster_size,
                                                cluster_size_sd,
                                                icc,
                                                nclusters,
                                                nfactors,
                                                nregcoefs,
                                                power,
                                                pre_post_corr,
                                                sigma_y);
      success <- 1;
    }
  }
  # Finish up
  if(success == 0) {
    stop("Macro error: I'm sorry, I can't do that kind of calculation or\n that kind of model yet.");
  }
  results <- list();
  results$power <- power;
  if(has_specified_sigma_y == 0) {
    results$std_regression_coefficient <- the_coef/sigma_y;
  } else {
    results$raw_regression_coefficient <- the_coef;
    results$std_regression_coefficient <- the_coef/sigma_y;
  }
  if(assignment == "unclustered"){
    results$ntotal <- ntotal;
  } else {
    results$nclusters <- nclusters;
    results$ntotal <- ntotal;
  }

  return(results);

}

CalculatePowerUnclustered <- function(alpha,
                                      nfactors,
                                      nregcoefs,
                                      ntotal,
                                      sigma_effective,
                                      the_coef,
                                      silent){
  num <- ntotal*(the_coef^2);
  den <- sigma_effective^2;
  lambda <- num/den;
  if(is.infinite(lambda)){
    stop("Macro error: Could not calculate noncentrality parameter.");
  }
  # Would this be more accurate if it was
  #lambda <- ((ntotal-nregcoefs)*(the_coef^2))/(sigma_effective^2);
  #analogously to Cohen (1988)?
  df1 <- 1;
  df2 <- ntotal-nregcoefs;
  if(df2 < 1){
    stop("Macro error: The sample size would be too small to fit the model.");
  }
  crit <- qf(1-alpha,df1,df2);
  power <- round(1-pf(crit,df1,df2,min(100,lambda)),digits=4);
  if(silent == 0){
    print(paste("The calculated power is",as.character(power), sep=" "));
    ncells <- 2^nfactors;
    if(ntotal < ncells){
      print(paste("However, a complete factorial requires at least",as.character(ncells),"subjects.",sep=" "));
    }
  }
  return(power);
}


CalculateNUnclustered <- function( alpha,
                                   nfactors,
                                   nregcoefs,
                                   power,
                                   sigma_effective,
                                   the_coef ){
  bottom <- nregcoefs+1;
  top <- 10000000;
  BinarySearchDone <- 0;
  while(BinarySearchDone == 0){
    middle <- ceiling(bottom+0.5*(top-bottom));
    fmiddle <- CalculatePowerUnclustered( alpha,
                                          nfactors,
                                          nregcoefs,
                                          middle,
                                          sigma_effective,
                                          the_coef,
                                          1 );
    if(middle==top){
      BinarySearchDone <- 1;
    } else {
      if(fmiddle < power) bottom <- middle;
      if(fmiddle == power){
        bottom <- middle;
        top <- middle;
      }
      if(fmiddle > power) top <- middle;
    }
  }
  ntotal <- middle;
  print(paste("The calculated sample size is",as.character(ntotal), "subjects.",sep=" "));
  ncells <- 2^nfactors;
  if(ntotal < ncells) {
    print(paste("However, a complete factorial requires at least",as.character(ncells),"subjects.",sep=" "));
  }
  return(ntotal);
}


CalculateDetectableUnclustered <- function(alpha,
                                           nfactors,
                                           nregcoefs,
                                           ntotal,
                                           power,
                                           sigma_effective,
                                           sigma_y){
  bottom <- 0;
  top <- 5*sigma_y;
  BinarySearchDone <- 0;
  while(BinarySearchDone == 0){
    middle <- bottom+.5*(top-bottom);
    fmiddle <- CalculatePowerUnclustered( alpha,
                                          nfactors,
                                          nregcoefs,
                                          ntotal,
                                          sigma_effective,
                                          middle,
                                          1);
    if(abs(middle-top)<1e-5) {
      BinarySearchDone <- 1;
    } else {
      if(fmiddle < power) bottom <- middle;
      if(fmiddle == power){
        bottom <- middle;
        top <- middle;
      }
      if(fmiddle > power) top <- middle;
    }
  }
  the_coef <- middle;
  if(fmiddle < (power-.01)) stop("Macro error: Could not find a detectable effect size");
  print("The detectable effect size is estimated as follows:");
  if((sigma_y < 0.999999) | (sigma_y > 1.000001)){
    calculated_raw_coef <- the_coef;
    print("  As an unstandardized regression coefficient for either");
    print(sprintf("%-60s %8.4f","  a main effect or an interaction:", calculated_raw_coef));
    calculated_raw_main <- 2*the_coef;
    print(sprintf("%-60s %8.4f","  As an unstandardized mean difference for a main effect:", calculated_raw_main));
    calculated_raw_ixn <- 4*the_coef;
    print("  As an unstandardized difference in differences for ");
    print(sprintf("%-60s %8.4f","  a 2-way interaction:", calculated_raw_ixn));
  }
  calculated_std_coef <- the_coef/sigma_y;
  print("  As a standardized regression coefficient for either");
  print(sprintf("%-60s %8.4f","  a main effect or an interaction:", calculated_std_coef));
  calculated_d_main <- 2*the_coef/sigma_y;
  print("  As a standardized mean difference (Cohen d) for a ");
  print(sprintf("%-60s %8.4f","  main effect:", calculated_d_main));
  calculated_d_ixn <- 4*the_coef/sigma_y;
  print("  As a standardized difference in differences for ");
  print(sprintf("%-60s %8.4f","  a 2-way interaction:", calculated_d_ixn));
  calculated_effect_size_ratio <- (the_coef^2)/(sigma_y^2);
  print("  As a standardized effect size ratio (Cohen f squared) ");
  print(sprintf("%-60s %8.4f","  for a main effect or interaction:", calculated_effect_size_ratio));
  ncells <- 2^nfactors;
  if(ntotal < ncells){
    print(paste("However, a complete factorial requires at least ",as.character(ncells)," subjects.",sep=""));
  }

  return(the_coef);
}


CalculatePowerTwoLevel <- function(alpha,
                                   assignment,
                                   cluster_size,
                                   cluster_size_sd,
                                   icc,
                                   nclusters,
                                   nfactors,
                                   nregcoefs,
                                   sigma_effective,
                                   the_coef,
                                   silent ){
  num <- cluster_size*nclusters*the_coef*the_coef;
  cluster_size_cv <- cluster_size_sd/cluster_size;
  if(assignment == "between_clusters"){
    nadj <- cluster_size*(1+cluster_size_cv^2);
    deff <- (1+(nadj-1)*icc);
    df2 <- nclusters - nregcoefs;
  } else {
    deff <- 1;
    df2 <- floor(nclusters*cluster_size) - nregcoefs;
  }
  if(df2 < 1){
    stop("Macro error: The sample size would be too small to fit the model.");
  }
  den <- (sigma_effective^2)*deff;
  lambda <- num/den;
  if(is.infinite(lambda)){
    stop("Macro error: Could not calculate noncentrality parameter.");
  }
  df1 <- 1;
  crit <- qf(1-alpha,df1,df2,0);
  power <- round(1-pf(crit,df1,df2,min(100,lambda)),digits=4);
  if(silent == 0){
    print(paste("The calculated power is",as.character(power),sep=" "));
    if(assignment == "between_clusters"){
      ncells <- 2^nfactors;
      if(nclusters < ncells){
        print(paste("However, a complete factorial requires at least ",as.character(ncells)," clusters.",sep=""));
      }
    }
  }
  return(power);
}


CalculateNTwoLevel <- function(alpha,
                               assignment,
                               cluster_size,
                               cluster_size_sd,
                               icc,
                               nfactors,
                               nregcoefs,
                               power,
                               sigma_effective,
                               the_coef){
  if(assignment == "within_clusters"){
    bottom <- ceiling((nregcoefs+1)/cluster_size);
  } else {
    bottom <- nregcoefs+1;
  }
  top <- 10000000;
  BinarySearchDone <- 0;
  while(BinarySearchDone == 0){
    middle <- ceiling(bottom+.5*(top-bottom));
    fmiddle <- CalculatePowerTwoLevel(alpha,
                                      assignment,
                                      cluster_size,
                                      cluster_size_sd,
                                      icc,
                                      middle,
                                      nfactors,
                                      nregcoefs,
                                      sigma_effective,
                                      the_coef,
                                      1);
    if(middle == top){
      BinarySearchDone <- 1;
    } else {
      if(fmiddle < power) bottom <- middle;
      if(fmiddle == power){
        bottom <- middle;
        top <- middle;
      }
      if(fmiddle > power) top <- middle;
    }
  }
  nclusters <- middle;
  ncells <- 2^nfactors;
  if(cluster_size != 1) {
    print(paste("The calculated sample size is ",as.character(nclusters), " clusters.",sep=""));
  } else {
    print(paste("The calculated sample size is ",as.character(nclusters), " subjects.",sep=""));
  }
  ncells <- 2^nfactors;
  if(assignment == "between_clusters") {
    if(nclusters < ncells) {
      print(paste("However, a complete factorial requires at least ",as.character(ncells)," clusters.",sep=""));
    }
  }
  if(assignment == "within_clusters"){
    ntotal <- nclusters * cluster_size;
    if(ntotal < ncells) {
      print(paste("However, a complete factorial requires at least ",as.character(ncells)," subjects.",sep=""));
    }
  }
  if(assignment == "unclustered"){
    if(nclusters < ncells) {
      print(paste("However, a complete factorial requires at least ",as.character(ncells)," subjects.",sep=""));
    }
  }
  return(nclusters);
}


CalculateDetectableTwoLevel <- function( alpha,
                                         assignment,
                                         cluster_size,
                                         cluster_size_sd,
                                         icc,
                                         nclusters,
                                         nfactors,
                                         nregcoefs,
                                         power,
                                         sigma_effective,
                                         sigma_y){
  bottom <- 0;
  top <- 5*sigma_y;
  BinarySearchDone <- 0;
  while(BinarySearchDone == 0){
    middle <- bottom+.5*(top-bottom);
    fmiddle <- CalculatePowerTwoLevel(alpha,
                                      assignment,
                                      cluster_size,
                                      cluster_size_sd,
                                      icc,
                                      nclusters,
                                      nfactors,
                                      nregcoefs,
                                      sigma_effective,
                                      middle,
                                      1);
    if(abs(middle-top) < 1e-5){
      BinarySearchDone <- 1;
    } else {
      if(fmiddle < power) bottom <- middle;
      if(fmiddle == power) {
        bottom <- middle;
        top <- middle;
      }
      if(fmiddle > power) top <- middle;
    }
  }
  the_coef <- middle;
  if(fmiddle < (power-.01)){
    stop("Macro error: Could not find a detectable effect size");
  }
  print("The detectable effect size is estimated as follows:");
  if((sigma_y < 0.999999)|(sigma_y > 1.000001)) {
    calculated_raw_coef <- the_coef;
    print("  As an unstandardized regression coefficient for either");
    print(sprintf("%-60s %8.4f","  a main effect or an interaction:", calculated_raw_coef));
    calculated_raw_main <- 2*the_coef;
    print(sprintf("%-60s %8.4f","  As an unstandardized mean difference for a main effect:", calculated_raw_main));
    calculated_raw_ixn <- 4*the_coef;
    print("  As an unstandardized difference in differences for ");
    print(sprintf("%-60s %8.4f","  a 2-way interaction:", calculated_raw_ixn));

  }
  calculated_std_coef <- the_coef/sigma_y;
  print("  As a standardized regression coefficient for either");
  print(sprintf("%-60s %8.4f","  a main effect or an interaction:", calculated_std_coef));
  calculated_d_main <- 2*the_coef/sigma_y;
  print("  As a standardized mean difference (Cohen d) for a ");
  print(sprintf("%-60s %8.4f","  main effect:", calculated_d_main));
  calculated_d_ixn <- 4*the_coef/sigma_y;
  print("  As a standardized difference in differences for ");
  print(sprintf("%-60s %8.4f","  a 2-way interaction:", calculated_d_ixn));
  calculated_effect_size_ratio <- (the_coef^2)/(sigma_y^2);
  print("  As a standardized effect size ratio (Cohen f squared) ");
  print(sprintf("%-60s %8.4f","  for a main effect or interaction:", calculated_effect_size_ratio));
  ncells <- 2^nfactors;
  if(assignment == "between_clusters") {
    if(nclusters < ncells) {
      print(paste("However, a complete factorial requires at least ",as.character(ncells)," clusters.",sep=""));
    }
  }
  return(the_coef);
}


CalculatePowerThreeLevel <- function( alpha,
                                      assignment,
                                      change_score_icc,
                                      cluster_size,
                                      cluster_size_sd,
                                      icc,
                                      nclusters,
                                      nfactors,
                                      nregcoefs,
                                      pre_post_corr,
                                      sigma_y,
                                      the_coef,
                                      silent){
  cluster_size_cv <- cluster_size_sd/cluster_size;
  nadj <- cluster_size*(1+cluster_size_cv^2);
  deff = (1+(nadj-1)*change_score_icc);
  if((assignment == "between_clusters")| assignment == "unclustered"){
    df2 <- nclusters - nregcoefs;
  }
  if(assignment == "within_clusters") {
    df2 <- floor(nclusters*cluster_size) - nregcoefs;
  }
  if(df2 < 1){
    stop("Macro error: The sample size would be too small to fit the model.");
  }
  num <- cluster_size*nclusters*the_coef*the_coef*(1-change_score_icc);
  den <- 2*(sigma_y^2)*(1-pre_post_corr)*(1-icc)*deff;
  lambda <- num/den;
  if(is.infinite(lambda)){
    stop("Macro error: Could not calculate noncentrality parameter.");
  }
  df1 <- 1;
  crit <- qf(1-alpha,df1,df2,0);
  power <- round(1-pf(crit,df1,df2,min(100,lambda)),digits=4);
  if(silent == 0) {
    print(paste("The calculated power is ",as.character(power),sep=""));
    ncells <- 2^nfactors;
    ntotal <- nclusters*cluster_size;
    if(assignment == "between_clusters") {
      if(nclusters < ncells) {
        print(paste("However, a complete factorial requires at least ",as.character(ncells)," clusters."));
      } else {
        if(ntotal < ncells) {
          print(paste("However, a complete factorial requires at least ",as.character(ncells)," subjects.",sep=""));
        }
      }
    }
  }
  return(power);
}

CalculateNThreeLevel <- function(alpha,
                                 assignment,
                                 change_score_icc,
                                 cluster_size,
                                 cluster_size_sd,
                                 icc,
                                 nfactors,
                                 nregcoefs,
                                 power,
                                 pre_post_corr,
                                 sigma_y,
                                 the_coef){
  if(assignment == "between_clusters"){
    bottom <- nregcoefs+1;
  } else {
    bottom <- ceiling((nregcoefs+1)/cluster_size);
  }
  top <- 10000000;
  BinarySearchDone <- 0;
  while(BinarySearchDone == 0){
    middle = ceiling(bottom+.5*(top-bottom));
    fmiddle = CalculatePowerThreeLevel(alpha,
                                       assignment,
                                       change_score_icc,
                                       cluster_size,
                                       cluster_size_sd,
                                       icc,
                                       middle,
                                       nfactors,
                                       nregcoefs,
                                       pre_post_corr,
                                       sigma_y,
                                       the_coef,
                                       1);
    if(middle == top){
      BinarySearchDone <- 1;
    } else {
      if(fmiddle < power) bottom <- middle;
      if(fmiddle == power){
        bottom <- middle;
        top <- middle;
      }
      if(fmiddle > power) top <- middle;
    }
  }
  nclusters <- middle;
  ncells <- 2^nfactors;
  if(cluster_size != 1){
    print(paste("The calculated sample size is ", as.character(nclusters)," clusters.", sep=""));
  } else {
    print(paste("The calculated sample size is ", as.character(nclusters)," subjects.", sep=""));
  }
  if(assignment == "unclustered"){
    if(nclusters < ncells){
      print(paste("However, a complete factorial requires at least ",as.character(ncells)," subjects.", sep=""));
    }
  }
  if(assignment == "between_clusters"){
    if(nclusters < ncells){
      print(paste("However, a complete factorial requires at least ",as.character(ncells)," clusters.",sep=""));
    }
  }
  if(assignment == "within_clusters"){
    if(nclusters*cluster_size < ncells){
      print(paste("However, a complete factorial requires at least ",as.character(ncells)," total subjects.",sep=""));
    }
  }
  return(nclusters);
}

CalculateDetectableThreeLevel <- function(alpha,
                                          assignment,
                                          change_score_icc,
                                          cluster_size,
                                          cluster_size_sd,
                                          icc,
                                          nclusters,
                                          nfactors,
                                          nregcoefs,
                                          power,
                                          pre_post_corr,
                                          sigma_y){
  bottom <- 0;
  top <- 5*sigma_y;
  BinarySearchDone <- 0;
  while(BinarySearchDone == 0){
    middle <- bottom+.5*(top-bottom);
    fmiddle <- CalculatePowerThreeLevel( alpha,
                                         assignment,
                                         change_score_icc,
                                         cluster_size,
                                         cluster_size_sd,
                                         icc,
                                         nclusters,
                                         nfactors,
                                         nregcoefs,
                                         pre_post_corr,
                                         sigma_y,
                                         middle,
                                         1);
    if(abs(middle-top) < 1e-5){
      BinarySearchDone <- 1;
    } else {
      if(fmiddle < power) bottom <- middle;
      if(fmiddle == power) {
        bottom <- middle;
        top <- middle;
      }
      if(fmiddle > power) top <- middle;
    }
  }
  the_coef <- middle;
  if(fmiddle < (power-.01) ) {
    stop("Macro error: Could not find a detectable effect size");
  }
  print("The detectable effect size is estimated as follows:");
  if((sigma_y < 0.999999)|(sigma_y > 1.000001)){
    calculated_raw_coef <- the_coef;
    print("  As an unstandardized regression coefficient for either");
    print(sprintf("%-60s %8.4f","  a main effect or an interaction:", calculated_raw_coef));
    calculated_raw_main <- 2*the_coef;
    print(sprintf("%-60s %8.4f","  As an unstandardized mean difference for a main effect:", calculated_raw_main));
    calculated_raw_ixn <- 4*the_coef;
    print("  As an unstandardized difference in differences for ");
    print(sprintf("%-60s %8.4f","  a 2-way interaction:", calculated_raw_ixn));
  }
  calculated_std_coef <- the_coef/sigma_y;
  print("  As a standardized regression coefficient for either");
  print(sprintf("%-60s %8.4f","  a main effect or an interaction:", calculated_std_coef));
  calculated_d_main <- 2*the_coef/sigma_y;
  paste("  As a standardized mean difference (Cohen d) for a ");
  print(sprintf("%-60s %8.4f","  main effect:", calculated_d_main));
  calculated_d_ixn <- 4*the_coef/sigma_y;
  print("  As a standardized difference in differences for ");
  print(sprintf("%-60s %8.4f","  a 2-way interaction:", calculated_d_ixn));
  calculated_effect_size_ratio <- (the_coef^2)/(sigma_y^2);
  print("  As a standardized effect size ratio (Cohen f squared) ");
  print(sprintf("%-60s %8.4f","  for a main effect or interaction:", calculated_effect_size_ratio));
  ncells <- 2^nfactors;
  ntotal <- nclusters*cluster_size;
  if(assignment == "between_clusters") {
    if(nclusters<ncells){
      print(paste("However, a complete factorial requires at least ",as.character(ncells)," clusters.",sep=""));
    }
  } else {
    if(ntotal < ncells) {
      print(paste("However, a complete factorial requires at least ",as.character(ncells)," subjects.", sep=""));
    }
  }
  return(the_coef);
}



