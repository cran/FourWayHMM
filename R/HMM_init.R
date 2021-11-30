#' Initialization for the ECM algorithm
#'
#' Runs the initialization of the ECM algorithm used for fitting the parsimonious hidden Markov models.
#' Parallel computing is implemented and highly recommended for a faster calculation.
#'
#' @param X An array of dimension \code{p} x \code{r} x \code{n} x \code{t}, where \code{p} is the number of
#'     variables in the rows of each data matrix, \code{r} is the number of variables in the columns of each
#'     data matrix, \code{n} is the number of data observations and \code{t} is the number of times.
#' @param k An integer or a vector indicating the number of states of the models.
#' @param nstartR An integer specifying the number of random starts to be considered.
#' @param ncores A positive integer indicating the number of cores used for running in parallel.
#' @param verbose A logical indicating whether the running output should be displayed.
#'
#' @return
#' \item{init}{A list of objects to be used by the \code{HMM.fit()} function.}
#' @export
#' @importFrom foreach %dopar%
#' @examples
#' data(simX)
#'
#' init <- HMM.init(X = simX, k = 2, nstartR = 1)
HMM.init <- function(X, k = 1:3, nstartR = 100, ncores = 1, verbose = FALSE) {
  results <- vector(mode = "list", length = length(k))

  for (g in 1:length(k)) {
    if (verbose == TRUE) {
      print(paste("Initializing Parsimonious Matrix Normal HMM with k =", k[g]))
    }

    results[[g]] <- r_Pars_init(X = X, k = k[g], nstartR = nstartR, ncores = ncores)
  }

  return(results)
}
