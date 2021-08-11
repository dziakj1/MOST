#' @title  Random Assignment Generator for a Factorial Experiment with Many Conditions
#'
#' @description This function provides a list of random numbers that can be used to assign
#' participants to conditions (cells) in an experiment with many conditions, such as a factorial
#' experiment.  The randomization is restricted as follows:  if the number of participants
#' available is a multiple of the number of conditions, then cell sizes will be
#' balanced; otherwise, they will be as near balanced as possible.
#'
#'
#' @param N The total number of participants to be randomized.
#' @param C The total number of conditions for the experiment you are planning.  Note that f
#'        or a complete factorial experiment having k factors, there will be 2^k conditions.
#'

#' @return A dataframe with 1 variable ranList with N observations, each observation of ranList
#'         provides a random number for each participant. This will be a number from 1 to C.
#'         For example, if the 4th number in the list is 7, the 4th subject is randomly assigned
#'         to experiment condition 7. Random numbers will be generated so that the experiment is
#'         approximately balanced.
#'

#' @export RandomAssignmentGenerator
#' @examples
#' result <- RandomAssignmentGenerator(35,17)
#' print(result)

RandomAssignmentGenerator <- function(N,
                                      C){

  numloop <- N %/% C;
  size1 <- N %% C;
  ranList <- c();
  for(i in 1:numloop) {
     tmp <- sample(1:C,size=C,replace=FALSE);
     ranList <- append(ranList,tmp);
  }
  tmp <- sample(1:C,size=size1,replace=FALSE);
  ranList <- append(ranList,tmp);
  #print(ranList);
  result <- data.frame(ranList);

  return(result);
}

