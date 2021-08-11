#' @title  The relative cost of reduced factorial designs
#'
#' @description Draw a graph of the relative cost of complete factorial,
#'              fractional factorial, and unbalanced reduced factorial designs,
#'              as presented by Collins, Dziak and Li
#'              (2009; https://www.ncbi.nlm.nih.gov/pubmed/19719358).  For
#'              purposes of illustration, a normally distributed response variable,
#'              dichotomous factors, and negligible interactions are assumed in this function.
#'
#' @param number_of_factors The number of factors to be tested.
#' @param desired_fract_resolution The desired resolution of the fractional factorial experiment to be compared.
#'                                 The default value is set to be 4.
#' @param min_target_d_per_factor The minimum Cohen's d (standardized difference, i.e., response difference
#'                                between levels on a given factor, divided by response standard deviation)
#'                                that is desired to be detected with 80% probability for any given factor.
#'                                The default value is set to be 0.2.
#' @param condition_costlier_than_subject The default value is set to be 1.
#' @param max_graph_ratio The default value is set to be 5.
#'
#' @importFrom graphics legend matplot mtext
#'
#' @export RelativeCosts1
#' @examples
#' RelativeCosts1(number_of_factors = 9,
#'               desired_fract_resolution = 4,
#'               min_target_d_per_factor = .2,
#'               condition_costlier_than_subject=1,
#'               max_graph_ratio = 4)

RelativeCosts1 <- function(number_of_factors,
                           desired_fract_resolution=4,
                           min_target_d_per_factor=.2,
                           condition_costlier_than_subject=1,
                           max_graph_ratio = 5){
  if(desired_fract_resolution < 3 ){
    stop('ERROR! The requested resolution is too low.');
  }
  if(desired_fract_resolution > 6 ){
    stop('ERROR! The requested resolution is too high.');
  }
  if(min_target_d_per_factor < 0.01 ){
    stop('ERROR! The requested target effect size is too low.');
  }
  if(number_of_factors < 2 ){
    stop('ERROR! The requested number of factors is too low.');
  }
  if(number_of_factors > 20 ){
    stop('ERROR! The requested number of factors is too high.');
  }
  if(number_of_factors > 20 ){
    stop('ERROR! The requested number of factors is too high.');
  }
  variance_adjustment <- 1;

  M <- matrix( c(1,		2,	2,		2,		2,		2,		     2,
                 2,		2,	4,		4,		4,		4,		     4,
                 3,		2,	4,		8,		8,		8,		     8,
                 4,		2,	8,		8,		16,		16,		     16,
                 5,		2,	8,		16,		16,		32,		     32,
                 6,		2,	8,		16,		32,		32,		     64,
                 7,		2,	8,		16,		64,		64,		    128,
                 8,		2,	16,		16,		64,		128,		  256,
                 9,		2,	16,		32,		128,	128,		  512,
                 10,	2,	16,		32,		128,	256,	   1024,
                 11,	2,	16,		32,		128,	256,		 2048,
                 12,	2,	16,		32,		256,	256,	   4096,
                 13,	2,	16,		32,		256,	512,	   8192,
                 14,	2,	16,		32,		256,	512,	  16384,
                 15,	2,	16,		32,		256,	512,	  32768,
                 16,	2,	32,		32,		256,	512,	  65536,
                 17,	2,	32,		64,		256,	512,	 131072,
                 18,	2,	32,		64,		512,	512,	 262144,
                 19,	2,	32,		64,		512,	1024,	 524288,
                 20,	2,	32,		64,		512,	1024,	1048576),
                 nrow = 20, ncol = 7, byrow = TRUE);


  fract_conditions <- M[number_of_factors,desired_fract_resolution];
  n_star <- ceiling(2*variance_adjustment*( 0.841621 + 1.959964)^2/(min_target_d_per_factor^2));
  separate_expers_conditions <- number_of_factors*2;
  single_factor_conditions<- number_of_factors+1;
  complete_factorial_conditions <- 2^number_of_factors;
  fract_factorial_conditions <- fract_conditions;
  separate_expers_subjects <- n_star * separate_expers_conditions;
  single_factor_subjects <- n_star * single_factor_conditions;
  complete_factorial_min <- 2*n_star;
  complete_factorial_subjects <- ceiling(complete_factorial_min/complete_factorial_conditions)*complete_factorial_conditions;
  fract_factorial_min <- 2*n_star;
  fract_factorial_subjects <- ceiling(fract_factorial_min/fract_factorial_conditions)*fract_factorial_conditions;
  separate_expers_sub_per_cell <- separate_expers_subjects/separate_expers_conditions;
  single_factor_sub_per_cell <- single_factor_subjects/single_factor_conditions;
  complete_factorial_sub_per_cell <- complete_factorial_subjects/complete_factorial_conditions;
  fract_factorial_sub_per_cell <- fract_factorial_subjects/fract_factorial_conditions;

  print(paste("Doing separate experiments on each of the",as.character(number_of_factors),"factors requires at least:",sep=" "));
  print(paste(as.character(separate_expers_subjects),"subjects total, i.e.",as.character(separate_expers_sub_per_cell),
              "subjects in each of", as.character(separate_expers_conditions),"cells."));
  cat("\n");

  print(paste("A comparative setup with",as.character(number_of_factors),"groups plus a control group requires at least:",sep=" "));
  print(paste(as.character(single_factor_subjects),"subjects total, i.e.",as.character(single_factor_sub_per_cell),
              "subjects in each of", as.character(single_factor_conditions),"cells."));
  cat("\n");

  print(paste("A 2^",as.character(number_of_factors)," complete factorial experiment requires at least:",sep=""));
  print(paste(as.character(complete_factorial_subjects),"subjects total, i.e.",as.character(complete_factorial_sub_per_cell),
              "subjects in each of", as.character(complete_factorial_conditions),"cells."));
  cat("\n");

  print(paste("A resolution",as.character(desired_fract_resolution),"fractional factorial experiment with",as.character(number_of_factors),"factors requires at least:",sep=" "));
  print(paste(as.character(fract_factorial_subjects),"subjects total, i.e.",as.character(fract_factorial_sub_per_cell),
              "subjects in each of", as.character(fract_factorial_conditions),"cells."));
  cat("\n");

  str1 = paste(as.character(number_of_factors),"separate experiments",sep=" ");
  str2 = paste("Single factor, ",as.character(number_of_factors),"+1 levels",sep="");
  str3 = paste("Complete 2^",as.character(number_of_factors)," factorial",sep="");
  str4 = paste("Fract. fact., ",as.character(number_of_factors),"factors, resol.",as.character(desired_fract_resolution),sep=" ");

  if(condition_costlier_than_subject == 1){
    separate_expers_cost <- function(cost_ratio_for_graph_macro) separate_expers_subjects + separate_expers_conditions * cost_ratio_for_graph_macro;
    single_factor_cost <- function(cost_ratio_for_graph_macro) single_factor_subjects + single_factor_conditions * cost_ratio_for_graph_macro;
    complete_factorial_cost <- function(cost_ratio_for_graph_macro) complete_factorial_subjects + complete_factorial_conditions * cost_ratio_for_graph_macro;
    fract_factorial_cost <- function(cost_ratio_for_graph_macro) fract_factorial_subjects + fract_factorial_conditions * cost_ratio_for_graph_macro;
    cost_ratio_for_graph_macro <- seq(0,200,10)
    x <- 200;
    ylim1 <- separate_expers_subjects + separate_expers_conditions * x;
    ylim2 <- single_factor_subjects + single_factor_conditions * x;
    ylim3 <- complete_factorial_subjects + complete_factorial_conditions * x;
    ylim4 <- fract_factorial_subjects + fract_factorial_conditions * x;
    ylist <- c(ylim1,ylim2,ylim3,ylim4)
    max=max(ylist);
    second_max <- max(ylist[ylist!=max]);
    ymax <- max_graph_ratio*10^floor(log10(second_max));
    if(max<ymax) ymax=max;

    matplot(cost_ratio_for_graph_macro,cbind(separate_expers_cost(cost_ratio_for_graph_macro),single_factor_cost(cost_ratio_for_graph_macro),
                                             complete_factorial_cost(cost_ratio_for_graph_macro),fract_factorial_cost(cost_ratio_for_graph_macro)),
            xlab="Overhead cost per condition", ylab="Cost", main="Costs of Different Approaches",
            type="l", ylim =range(0:ymax), lwd=2, lty=c(1,2,3,6),col=c("black","red","green","blue"));
    mtext("relative to per-subject cost fixed at 1 unit",side=3, adj=0);
    legend("topright", legend=c(str1,str2,str3,str4), xpd = TRUE,
           bty = "n", lty=c(1,2,3,6), lwd=2, col=c("black","red","darkgreen","blue"))

  }else{
    separate_expers_cost <- function(cost_ratio_for_graph_macro) separate_expers_conditions + separate_expers_subjects * cost_ratio_for_graph_macro;
    single_factor_cost <- function(cost_ratio_for_graph_macro) single_factor_conditions + single_factor_subjects * cost_ratio_for_graph_macro;
    complete_factorial_cost <- function(cost_ratio_for_graph_macro) complete_factorial_conditions + complete_factorial_subjects * cost_ratio_for_graph_macro;
    fract_factorial_cost <- function(cost_ratio_for_graph_macro) fract_factorial_conditions +  fract_factorial_subjects * cost_ratio_for_graph_macro;
    cost_ratio_for_graph_macro <- seq(0,200,10)
    x <- 200;
    ylim1 <- separate_expers_conditions + separate_expers_subjects * x;
    ylim2 <- single_factor_conditions + single_factor_subjects * x;
    ylim3 <- complete_factorial_conditions + complete_factorial_subjects * x;
    ylim4 <- fract_factorial_conditions + fract_factorial_subjects * x;
    ylist <- c(ylim1,ylim2,ylim3,ylim4)
    max=max(ylist);
    second_max <- max(ylist[ylist!=max]);
    ymax <- max_graph_ratio*10^floor(log10(second_max));
    if(max<ymax) ymax=max;

    matplot(cost_ratio_for_graph_macro,cbind(separate_expers_cost(cost_ratio_for_graph_macro),single_factor_cost(cost_ratio_for_graph_macro),
                                             complete_factorial_cost(cost_ratio_for_graph_macro),fract_factorial_cost(cost_ratio_for_graph_macro)),
            xlab="Cost per individual subject", ylab="Cost", main="Costs of Different Approaches",
            type="l", ylim =range(0:ymax), lwd=2, lty=c(1,2,3,6),col=c("black","red","green","blue"));
    mtext("relative to per-condition overhead cost fixed at 1 unit", side=3, adj=0);
    legend("topright", legend=c(str1,str2,str3,str4), xpd = TRUE,
           bty = "n", lty=c(1,2,3,6), lwd=2, col=c("black","red","darkgreen","blue"))
  }

}
