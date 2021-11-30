#' Fitting for parsimonious hidden Markov models for four-way data
#'
#' Fits, by using an ECM algorithm, parsimonious hidden Markov models to the given four-way data.
#' Parallel computing is implemented and highly recommended for a faster model fitting. The Bayesian
#' information criterion (BIC) is used to select the best fitting model.

#' @param X An array of dimension \code{p} x \code{r} x \code{n} x \code{t}, where \code{p} is the number of
#'     variables in the rows of each data matrix, \code{r} is the number of variables in the columns of each
#'     data matrix, \code{n} is the number of data observations and \code{t} is the number of times.
#' @param k An integer or a vector indicating the number of states of the models.
#' @param init.par The initial values for starting the algorithms, as produced by the \code{HMM.init()} function.
#' @param mod.row A character vector indicating the parsimonious structure of the row covariance matrix.
#'      Possible values are: "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV",
#'      "EVV", "VVV" or "all". When "all" is used, all of the 14 row parsimonious structures are considered.
#' @param mod.col A character vector indicating the parsimonious structure of the column covariance matrix.
#'      Possible values are: "II", "EI", "VI", "EE", "VE", "EV", "VV", or "all". When "all" is used, all of
#'      the 7 column parsimonious structures are considered.
#' @param ncores A positive integer indicating the number of cores used for running in parallel.
#' @param verbose A logical indicating whether the running output should be displayed.
#' @param ret.all A logical indicating whether to report the results of all the models or only those of the best
#'     model according to the BIC.
#'
#' @return A list with the following elements:
#' \item{all.models}{The results related to the all the fitted models (only when \code{ret.all = TRUE}).}
#' \item{BicWin}{The best fitting model according to the BIC.}
#' \item{Summary}{A quick table showing summary results for the best fitting model according to the BIC.}
#' \item{c.time}{Provides information on the computational times required to fit all the models for each state.}
#' @export
#' @importFrom foreach %dopar%
#' @examples
#' data(simX)
#'
#' init <- HMM.init(X = simX, k = 2, nstartR = 1)
#' res <- HMM.fit(X = simX, k = 2, init.par = init, mod.row = "VII", mod.col = "EE")
HMM.fit <- function(X, k = 1:3, init.par = NULL, mod.row = "all", mod.col = "all", ncores = 1, verbose = FALSE, ret.all = FALSE) {
  if (any(mod.row == "all")) {
    model.row <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV")
  } else {
    model.row <- mod.row
  }
  if (any(mod.col == "all")) {
    model.col <- c("II", "EI", "VI", "EE", "VE", "EV", "VV")
  } else {
    model.col <- mod.col
  }

  pt.mod <- length(model.row) * length(model.col)
  list.mod <- expand.grid(model.row, model.col)
  l <- NULL

  oper <- vector(mode = "list", length = length(k))
  time <- numeric(length(k))

  for (g in 1:length(k)) {
    ptm <- proc.time()

    if (verbose == TRUE) {
      if (length(k) == 1) {
        print(paste("Fitting Parsimonious Matrix Normal HMM with k =", k))
      } else {
        print(paste("Fitting Parsimonious Matrix Normal HMM with k =", k[g]))
      }
    }

    cluster <- snow::makeCluster(ncores, type = "SOCK")
    doSNOW::registerDoSNOW(cluster)

    pb <- utils::txtProgressBar(max = pt.mod, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    oper[[g]] <- foreach::foreach(l = 1:pt.mod, .combine = "comb", .export = c("HMM_Pars_MVN"), .packages = c("tensor"), .multicombine = TRUE, .init = list(list()), .options.snow = opts) %dopar% {
      res <- tryCatch(HMM_Pars_MVN(X = X, k = k[g], init.par = init.par[[g]], mod.row = as.character(list.mod[l, 1]), mod.col = as.character(list.mod[l, 2])), error = function(e) {
        NA
      })

      list(res)
    }

    snow::stopCluster(cluster)
    foreach::registerDoSEQ()
    close(pb)

    ptm2 <- proc.time() - ptm
    time[g] <- ptm2[3]
  }

  BICres <- array(NA, dim = c(pt.mod, 1, length(k)))

  for (t in 1:length(k)) {
    for (m in 1:pt.mod) {
      if (length(oper[[t]][[1]][[m]]) > 1) {
        if (oper[[t]][[1]][[m]][["mark"]] == 0 & oper[[t]][[1]][[m]][["check"]] == 1) {
          BICres[m, 1, t] <- oper[[t]][[1]][[m]][["BIC"]]
        } else {
          oper[[t]][[1]][[m]] <- NA
        }
      } else {
        BICres[m, 1, t] <- NA
      }
    }
  }

  df <- data.frame(matrix(NA, nrow = 1, ncol = 3), row.names = c("BIC"))
  colnames(df) <- c("Model", "Value", "G")

  complete.model <- cbind(rep(as.character(list.mod$Var1), length(k)), rep(as.character(list.mod$Var2), length(k)))

  df[1, 1] <- paste(complete.model[which.min(BICres), 1], complete.model[which.min(BICres), 2], sep = "-")
  df[1, 2] <- min(BICres, na.rm = TRUE)

  if (length(k) == 1) {
    df[1, 3] <- k
  } else {
    df[1, 3] <- k[which(BICres == min(BICres, na.rm = TRUE), arr.ind = TRUE)[1, 3]]
  }

  b.w1 <- substr(df[1, 1], 1, 3)
  b.w2 <- substr(df[1, 1], 5, 6)

  if (length(k) == 1) {
    for (b in 1:length(oper[[1]][[1]])) {
      tryCatch(if (b.w1 == oper[[1]][[1]][[b]][["name"]][1] & b.w2 == oper[[1]][[1]][[b]][["name"]][2]) {
        BicWin <- oper[[1]][[1]][[b]]
      }, error = function(e) {})
    }
  } else {
    slb <- which(BICres == min(BICres, na.rm = TRUE), arr.ind = TRUE)[1, 3]

    for (b in 1:length(oper[[slb]][[1]])) {
      tryCatch(if (b.w1 == oper[[slb]][[1]][[b]][["name"]][1] & b.w2 == oper[[slb]][[1]][[b]][["name"]][2]) {
        BicWin <- oper[[slb]][[1]][[b]]
      }, error = function(e) {})
    }
  }

  if (ret.all == FALSE) {
    return(list(BicWin = BicWin, Summary = df, c.time = time))
  } else {
    return(list(all.models = oper, BicWin = BicWin, Summary = df, c.time = time))
  }
}
