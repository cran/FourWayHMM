comb <- function(x, ...) {
  lapply(
    seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
  )
}

r_Pars_init <- function(X, k, nstartR = 100, ncores = 14) {
  dMVnorm <- function(X, M, U, V) {
    num <- dim(X)[3] # sample size
    p <- nrow(X) # rows of X
    r <- ncol(X) # columns of X

    if (is.na(num)) {
      X <- as.matrix(X)
      delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
    } else {
      delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
    }

    pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)

    return(pdf)
  }
  transit <- function(Y) {

    # Y is a matrix N \times T

    N <- dim(Y)[1]
    T <- dim(Y)[2]
    K <- max(Y)

    if (K == 1) {
      PI <- matrix(1, 1, 1)
    }

    if (K > 1) {
      PI <- matrix(0, nrow = K, ncol = K)

      for (i in 1:N) {
        for (t in 2:T) {
          PI[Y[i, t - 1], Y[i, t]] <- PI[Y[i, t - 1], Y[i, t]] + 1
        }
      }
      PI <- diag(1 / rowSums(PI)) %*% PI
    }

    return(PI)
  }
  tr <- function(x) {
    return(sum(diag(x)))
  }

  # Dimensions

  p <- dim(X)[1]
  r <- dim(X)[2]
  num <- dim(X)[3]
  t <- dim(X)[4]

  # Create some objects

  prior <- numeric(k)
  M <- array(0, dim = c(p, r, k))
  sigmaU <- array(0, dim = c(p, p, k))
  sigmaV <- array(0, dim = c(r, r, k))
  nu <- numeric(k)

  WR <- array(0, dim = c(p, p, k))
  WC <- array(0, dim = c(r, r, k))
  tempWR <- array(0, dim = c(p, p, num * t))
  tempWC <- array(0, dim = c(r, r, num * t))
  tempM <- array(0, dim = c(p, r, num * t))

  l <- (p + r) * p^2
  if (l > .Machine$integer.max) {
    l <- .Machine$integer.max
  }

  post <- dens <- array(0, c(num * t, k), dimnames = list(1:(num * t), paste("comp.", 1:k, sep = "")))
  post2 <- array(NA, c(k, k, num, t - 1))

  Xresh <- array(X, dim = c(p, r, num * t))

  eu <- matrix(0, nrow = num * t, ncol = k)
  classy <- numeric(num * t)
  rand.start <- matrix(0, nstartR, k)

  withr::with_seed(l, {
    for (i in 1:nstartR) {
      rand.start[i, ] <- sample(c(1:num * t), k)
    }
  })

  cluster <- snow::makeCluster(ncores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)

  pb <- utils::txtProgressBar(max = nstartR, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  oper0 <- foreach::foreach(l = 1:nstartR, .combine = "comb", .packages = c("tensor"), .multicombine = TRUE, .init = list(list(), list(), list(), list(), list(), list(), list(), list(), list()), .options.snow = opts) %dopar% {

    ### part 0 ###

    sec <- rand.start[l, ]

    for (j in 1:k) {
      M[, , j] <- Xresh[, , sec[j]]
    }

    for (j in 1:k) {
      for (i in 1:(num * t)) {
        eu[i, j] <- base::norm((Xresh[, , i] - M[, , j]), type = "F")
      }
    }

    for (i in 1:(num * t)) {
      classy[i] <- which.min(eu[i, ])
    }

    z <- mclust::unmap(classy)

    ### part 1 ###

    for (j in 1:k) {
      M[, , j] <- rowSums(Xresh * z[, j][slice.index(Xresh, 3)], dims = 2) / sum(z[, j])

      pt1 <- sweep(Xresh, 1:2, M[, , j]) * z[, j][slice.index(sweep(Xresh, 1:2, M[, , j]), 3)]

      pt2 <- aperm(sweep(Xresh, 1:2, M[, , j]), c(2, 1, 3))

      WR[, , j] <- tensor::tensor(pt1, pt2, c(2, 3), c(1, 3))
    }

    phi <- tr(rowSums(WR, dims = 2)) / ((num * t) * p * r)

    for (j in 1:k) {
      sigmaU[, , j] <- phi * diag(1, p, p)
      sigmaV[, , j] <- diag(1, r, r)

      dens[, j] <- dMVnorm(X = Xresh, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j])
    }

    if (k == 1) {
      prior <- 1
    } else {
      prior <- colMeans(z)
    }

    ### part 2 ###

    f <- array(dens, c(num, t, k))
    PI <- transit(matrix(classy, num, t))

    A <- B <- array(NA, c(num, t, k)) # forward and backward variables
    A[, 1, ] <- matrix(rep(prior, each = num), ncol = k) * f[, 1, ] # log(foo) + lscale
    B[, t, ] <- 1
    for (T in 2:t) {
      A[, T, ] <- A[, T - 1, ] %*% PI * f[, T, ]
      tp <- t - T + 1
      B[, tp, ] <- t(PI %*% t(B[, tp + 1, ] * f[, tp + 1, ]))
    }

    part <- rowSums(matrix(A[, t, ], nrow = num, ncol = k))
    part[part == 0] <- 5e-324
    llk <- sum(log(part)) # log-likelihood

    # return

    list(sigmaU, sigmaV, f, A, B, PI, classy, llk, M)
  }

  snow::stopCluster(cluster)
  foreach::registerDoSEQ()
  close(pb)

  df <- data.frame(llk = unlist(oper0[[8]]), pos = c(1:nstartR))
  df <- tidyr::drop_na(df)
  df <- df[!is.infinite(rowSums(df)), ]

  bestR <- utils::head(data.table::setorderv(df, cols = "llk", order = -1), n = 1)$pos

  res <- vector(mode = "list", length = 9)

  for (i in 1:9) {
    res[[i]] <- oper0[[i]][bestR]
  }

  return(res)
}

HMM_Pars_MVN <- function(X, k, init.par = NULL, mod.row = NULL, mod.col = NULL, tol = 0.001, tol2 = 0.001, maxit = 500, maxit2 = 100) {
  ptm <- proc.time()

  # Fuctions

  dMVnorm <- function(X, M, U, V) {
    tr <- function(x) {
      return(sum(diag(x)))
    }

    num <- dim(X)[3] # sample size
    p <- nrow(X) # rows of X
    r <- ncol(X) # columns of X

    if (is.na(num)) {
      X <- as.matrix(X)
      delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
    } else {
      delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
    }

    pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)

    return(pdf)
  }
  tr <- function(x) {
    return(sum(diag(x)))
  }

  # Dimensions

  p <- dim(X)[1] # rows of each matrix;
  r <- dim(X)[2] # columns of each matrix;
  num <- dim(X)[3] # sample size;
  t <- dim(X)[4] # times.

  # Create objects related to M and prior

  prior <- numeric(k)
  M <- array(0, dim = c(p, r, k))
  tempM <- array(0, dim = c(p, r, num * t))

  ## Objects related to the two covariance matrices

  sigmaU <- array(0, dim = c(p, p, k))
  sigmaV <- array(0, dim = c(r, r, k))

  WR <- array(0, dim = c(p, p, k))
  WC <- array(0, dim = c(r, r, k))
  tempWR <- array(0, dim = c(p, p, num * t))
  tempWC <- array(0, dim = c(r, r, num * t))

  phi.k <- numeric(k)
  temp.phi <- vector("list", k) # for VEI, VEE, VEV
  temp.phi2 <- numeric(k) # for EVI
  temp.numdeltaR <- array(0, dim = c(p, p, k)) # for VEI, VEE, VEV
  deltaU.k <- array(0, dim = c(p, p, k))
  deltaV.k <- array(0, dim = c(r, r, k))
  ftemp.r <- array(0, dim = c(p, p, k)) # MM object
  ftemp.c <- array(0, dim = c(r, r, k)) # MM object
  numphi <- numeric(k) # for EVE, EVV
  gammaU.k <- array(0, dim = c(p, p, k))
  gammaV.k <- array(0, dim = c(r, r, k))
  tempW_EEV <- vector("list", k) # for EEV, VEV
  tempW_EV <- vector("list", k) # for EV
  tempomegaR <- array(0, dim = c(p, p, k)) # for EEV, VEV
  tempomegaC <- array(0, dim = c(r, r, k)) # for EEV, VEV
  V_EVV.U.K <- array(0, dim = c(p, p, k)) # for EVV

  ## Other objects

  post <- dens <- array(0, c(num * t, k), dimnames = list(1:(num * t), paste("comp.", 1:k, sep = "")))
  post2 <- array(NA, c(k, k, num, t - 1)) # transtion probabilities
  Xresh <- array(X, dim = c(p, r, num * t)) # data reshape

  # Preliminary definition of convergence criterions for EM/MM algorithms

  check <- 0
  check2 <- 0
  loglik.old <- -Inf
  loglik.new <- NULL
  ll <- NULL
  mark <- 1
  MM.r.old <- -Inf
  MM.c.old <- -Inf
  m.iter <- 0
  m.iter2 <- 0

  ### Algorithm ###

  print(paste(paste(paste("Fitting Matrix Normal HMM with k =", k), paste("and", mod.row)), paste("-", mod.col), paste("parsimonious structure")))

  oper0 <- init.par
  sigmaU <- oper0[[1]][[1]]
  sigmaV <- oper0[[2]][[1]]
  f <- oper0[[3]][[1]]
  A <- oper0[[4]][[1]]
  B <- oper0[[5]][[1]]
  PI <- oper0[[6]][[1]]
  classy <- oper0[[7]][[1]]

  if (mod.row == "VEI" | mod.row == "VEE" | mod.row == "VEV") {
    for (j in 1:k) {
      temp.phi[[j]] <- eigen(sigmaU[, , j])$values

      phi.k[j] <- (prod(temp.phi[[j]]))^(1 / p)
    }
  }
  if (mod.row == "EVE" | mod.row == "VVE") {
    TempW2R <- array(0, dim = c(p, p, k))
    Tempz <- mclust::unmap(classy)

    for (j in 1:k) {
      deltaU.k[, , j] <- diag(eigen(sigmaU[, , j])$values, p, p)
      TempW2R[, , j] <- sigmaU[, , j] * sum(Tempz[, j])
    }

    gammaU <- eigen(rowSums(TempW2R, dims = 2) / ((det(rowSums(TempW2R, dims = 2)))^(1 / p)))$vectors
  } else {
    gammaU <- matrix(0, p, p)
  }
  if (mod.row == "EII" | mod.row == "VII") {
    for (j in 1:k) {
      sigmaU[, , j] <- diag(1, p, p)
    }
  }

  if (mod.col == "VE") {
    deltaV.k <- array(0, dim = c(r, r, k))
    TempW2C <- array(0, dim = c(r, r, k))
    Tempz <- mclust::unmap(classy)

    for (j in 1:k) {
      deltaV.k[, , j] <- diag(eigen(sigmaV[, , j])$values, r, r)
      TempW2C[, , j] <- sigmaV[, , j] * sum(Tempz[, j])
    }

    gammaV <- eigen(rowSums(TempW2C, dims = 2) / ((det(rowSums(TempW2C, dims = 2)))^(1 / r)))$vectors
  } else {
    deltaV.k <- array(0, dim = c(r, r, k))
    gammaV <- matrix(0, r, r)
  }
  if (mod.col == "II") {
    for (j in 1:k) {
      sigmaV[, , j] <- diag(1, r, r)
    }
  }

  ### Estimation ###

  while (check < 1) {
    m.iter <- m.iter + 1

    ### E - STEP ###

    for (j in 1:k) {
      numer <- (A[, , j] * B[, , j])
      numer[numer < 5e-324] <- 5e-324
      denom <- apply(A * B, c(1, 2), sum)
      denom[denom < 5e-324] <- 5e-324

      post[, j] <- numer / denom # posterior probabilities (z)

      for (n in 1:num) {
        for (T in 1:(t - 1)) {
          post2[, , n, T] <- diag(A[n, T, ], nrow = k, ncol = k) %*% PI %*% diag(B[n, T + 1, ] * f[n, T + 1, ], nrow = k, ncol = k)
          post2[, , n, T] <- post2[, , n, T] / sum(post2[, , n, T]) # posterior probabilities (zz)
          post2[, , n, T][is.nan(post2[, , n, T])] <- 0
        }
      }
    }

    ### M - STEP ###

    for (j in 1:k) {
      M[, , j] <- rowSums(Xresh * post[, j][slice.index(Xresh, 3)], dims = 2) / sum(post[, j])
    }

    # ROWS COVARIANCE MATRIX

    for (j in 1:k) {
      pt1 <- aperm(tensor::tensor(sweep(Xresh, 1:2, M[, , j]) * post[, j][slice.index(sweep(Xresh, 1:2, M[, , j]), 3)], solve(sigmaV[, , j]), 2, 1), c(1, 3, 2))

      pt2 <- aperm(sweep(Xresh, 1:2, M[, , j]), c(2, 1, 3))

      WR[, , j] <- tensor::tensor(pt1, pt2, c(2, 3), c(1, 3))
    } # Scatter Matrix

    if (mod.row == "EII") {
      phi <- tr(rowSums(WR, dims = 2)) / ((num * t) * p * r)

      for (j in 1:k) {
        sigmaU[, , j] <- phi * diag(1, p, p)
      }
    }

    if (mod.row == "VII") {
      for (j in 1:k) {
        phi.k[j] <- tr(WR[, , j]) / (p * r * sum(post[, j]))
        sigmaU[, , j] <- phi.k[j] * diag(1, p, p)
      }
    }

    if (mod.row == "EEI") {
      deltaU <- diag(diag(rowSums(WR, dims = 2)), p, p) / (det(diag(diag(rowSums(WR, dims = 2)), p, p)))^(1 / p)

      phi <- (det(diag(diag(rowSums(WR, dims = 2)), p, p)))^(1 / p) / ((num * t) * r)

      for (j in 1:k) {
        sigmaU[, , j] <- phi * deltaU
      }
    }

    if (mod.row == "VEI") {
      for (j in 1:k) {
        temp.numdeltaR[, , j] <- (1 / phi.k[j]) * WR[, , j]
      }

      deltaU <- diag(diag(rowSums(temp.numdeltaR, dims = 2)), p, p) / (det(diag(diag(rowSums(temp.numdeltaR, dims = 2)), p, p)))^(1 / p)

      for (j in 1:k) {
        phi.k[j] <- (tr(solve(deltaU) %*% WR[, , j])) / (p * r * sum(post[, j]))

        sigmaU[, , j] <- phi.k[j] * deltaU
      }
    }

    if (mod.row == "EVI") {
      for (j in 1:k) {
        deltaU.k[, , j] <- diag(diag(WR[, , j]), p, p) / (det(diag(diag(WR[, , j]), p, p)))^(1 / p)

        temp.phi2[j] <- det(diag(diag(WR[, , j]), p, p))^(1 / p)
      }

      phi <- sum(temp.phi2) / ((num * t) * r)

      for (j in 1:k) {
        sigmaU[, , j] <- phi * deltaU.k[, , j]
      }
    }

    if (mod.row == "VVI") {
      for (j in 1:k) {
        deltaU.k[, , j] <- diag(diag(WR[, , j]), p, p) / (det(diag(diag(WR[, , j]), p, p)))^(1 / p)

        phi.k[j] <- det(diag(diag(WR[, , j]), p, p))^(1 / p) / (r * sum(post[, j]))

        sigmaU[, , j] <- phi.k[j] * deltaU.k[, , j]
      }
    }

    if (mod.row == "EEE") {
      for (j in 1:k) {
        sigmaU[, , j] <- rowSums(WR, dims = 2) / ((num * t) * r)
      }
    }

    if (mod.row == "VEE") {
      for (j in 1:k) {
        temp.numdeltaR[, , j] <- (1 / phi.k[j]) * WR[, , j]
      }

      deltaU <- rowSums(temp.numdeltaR, dims = 2) / ((det(rowSums(temp.numdeltaR, dims = 2)))^(1 / p))

      for (j in 1:k) {
        phi.k[j] <- tr(solve(deltaU) %*% WR[, , j]) / (p * r * sum(post[, j]))

        sigmaU[, , j] <- phi.k[j] * deltaU
      }
    }

    if (mod.row == "EVE") {
      while (check2 < 1) {
        m.iter2 <- m.iter2 + 1

        for (j in 1:k) {
          ftemp.r[, , j] <- tcrossprod(solve(deltaU.k[, , j]), gammaU) %*% WR[, , j] - max(eigen(WR[, , j])$values) * tcrossprod(solve(deltaU.k[, , j]), gammaU)
        }

        f <- rowSums(ftemp.r, dims = 2)

        MM.r.new <- tr(f %*% gammaU)

        if ((abs(MM.r.new - MM.r.old)) < tol2 | m.iter2 == maxit2) {
          check2 <- 1
          res.svd <- svd(f)
          gammaU <- tcrossprod(res.svd$v, res.svd$u)
        } else {
          res.svd <- svd(f)
          gammaU <- tcrossprod(res.svd$v, res.svd$u)
        }

        MM.r.old <- MM.r.new
      }

      m.iter2 <- 0
      check2 <- 0
      MM.r.old <- -Inf

      for (j in 1:k) {
        deltaU.k[, , j] <- diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p) / (det(diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p)))^(1 / p)

        numphi[j] <- tr(gammaU %*% tcrossprod(solve(deltaU.k[, , j]), gammaU) %*% WR[, , j])
      }

      phi <- sum(numphi) / ((num * t) * p * r)

      for (j in 1:k) {
        sigmaU[, , j] <- phi * gammaU %*% tcrossprod(deltaU.k[, , j], gammaU)
      }
    }

    if (mod.row == "VVE") {
      while (check2 < 1) {
        m.iter2 <- m.iter2 + 1

        for (j in 1:k) {
          ftemp.r[, , j] <- tcrossprod(solve(deltaU.k[, , j]), gammaU) %*% WR[, , j] - max(eigen(WR[, , j])$values) * tcrossprod(solve(deltaU.k[, , j]), gammaU)
        }

        f <- rowSums(ftemp.r, dims = 2)

        MM.r.new <- tr(f %*% gammaU)

        if ((abs(MM.r.new - MM.r.old)) < tol2 | m.iter2 == maxit2) {
          check2 <- 1
          res.svd <- svd(f)
          gammaU <- tcrossprod(res.svd$v, res.svd$u)
        } else {
          res.svd <- svd(f)
          gammaU <- tcrossprod(res.svd$v, res.svd$u)
        }

        MM.r.old <- MM.r.new
      }

      m.iter2 <- 0
      check2 <- 0
      MM.r.old <- -Inf

      for (j in 1:k) {
        deltaU.k[, , j] <- diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p) / (det(diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p)))^(1 / p)
        phi.k[j] <- (det(diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p))^(1 / p)) / (r * sum(post[, j]))
        sigmaU[, , j] <- phi.k[j] * gammaU %*% tcrossprod(deltaU.k[, , j], gammaU)
      }
    }

    if (mod.row == "EEV") {
      for (j in 1:k) {
        tempW_EEV[[j]] <- eigen(WR[, , j])

        gammaU.k[, , j] <- tempW_EEV[[j]][["vectors"]]

        tempomegaR[, , j] <- diag(tempW_EEV[[j]][["values"]], p, p)
      }

      deltaU <- rowSums(tempomegaR, dims = 2) / ((det(rowSums(tempomegaR, dims = 2)))^(1 / p))

      phi <- ((det(rowSums(tempomegaR, dims = 2)))^(1 / p)) / ((num * t) * r)

      for (j in 1:k) {
        sigmaU[, , j] <- phi * gammaU.k[, , j] %*% tcrossprod(deltaU, gammaU.k[, , j])
      }
    }

    if (mod.row == "VEV") {
      for (j in 1:k) {
        tempW_EEV[[j]] <- eigen(WR[, , j])

        gammaU.k[, , j] <- tempW_EEV[[j]][["vectors"]]

        tempomegaR[, , j] <- diag(tempW_EEV[[j]][["values"]], p, p)

        temp.numdeltaR[, , j] <- (1 / phi.k[j]) * tempomegaR[, , j]
      }

      deltaU <- rowSums(temp.numdeltaR, dims = 2) / ((det(rowSums(temp.numdeltaR, dims = 2)))^(1 / p))

      for (j in 1:k) {
        phi.k[j] <- tr(tempomegaR[, , j] %*% solve(deltaU)) / (p * r * sum(post[, j]))

        sigmaU[, , j] <- phi.k[j] * gammaU.k[, , j] %*% tcrossprod(deltaU, gammaU.k[, , j])
      }
    }

    if (mod.row == "EVV") {
      for (j in 1:k) {
        V_EVV.U.K[, , j] <- WR[, , j] / ((det(WR[, , j]))^(1 / p))

        numphi[j] <- det(WR[, , j])^(1 / p)
      }

      phi <- sum(numphi) / ((num * t) * r)

      for (j in 1:k) {
        sigmaU[, , j] <- phi * V_EVV.U.K[, , j]
      }
    }

    if (mod.row == "VVV") {
      for (j in 1:k) {
        sigmaU[, , j] <- WR[, , j] / (r * sum(post[, j]))
      }
    }

    # COLUMNS COVARIANCE MATRIX

    for (j in 1:k) {
      pt1 <- aperm(tensor::tensor(aperm(sweep(Xresh, 1:2, M[, , j]), c(2, 1, 3)) * post[, j][slice.index(aperm(sweep(Xresh, 1:2, M[, , j]), c(2, 1, 3)), 3)], solve(sigmaU[, , j]), 2, 1), c(1, 3, 2))

      pt2 <- sweep(Xresh, 1:2, M[, , j])

      WC[, , j] <- tensor::tensor(pt1, pt2, c(2, 3), c(1, 3))
    } # Scatter Matrix

    if (mod.col == "II") {
      for (j in 1:k) {
        sigmaV[, , j] <- diag(1, r, r)
      }
    }

    if (mod.col == "EI") {
      deltaV <- diag(diag(rowSums(WC, dims = 2)), r, r) / (det(diag(diag(rowSums(WC, dims = 2)), r, r)))^(1 / r)

      for (j in 1:k) {
        sigmaV[, , j] <- deltaV
      }
    }

    if (mod.col == "VI") {
      for (j in 1:k) {
        sigmaV[, , j] <- diag(diag(WC[, , j]), r, r) / (det(diag(diag(WC[, , j]), r, r)))^(1 / r)
      }
    }

    if (mod.col == "EE") {
      for (j in 1:k) {
        sigmaV[, , j] <- rowSums(WC, dims = 2) / ((det(rowSums(WC, dims = 2)))^(1 / r))
      }
    }

    if (mod.col == "VE") {
      while (check2 < 1) {
        m.iter2 <- m.iter2 + 1

        for (j in 1:k) {
          ftemp.c[, , j] <- tcrossprod(solve(deltaV.k[, , j]), gammaV) %*% WC[, , j] - max(eigen(WC[, , j])$values) * tcrossprod(solve(deltaV.k[, , j]), gammaV)
        }

        f.C <- rowSums(ftemp.c, dims = 2)

        MM.c.new <- tr(f.C %*% gammaV)

        if ((abs(MM.c.new - MM.c.old)) < tol2 | m.iter2 == maxit2) {
          check2 <- 1
          res.svd.C <- svd(f.C)
          gammaV <- tcrossprod(res.svd.C$v, res.svd.C$u)
        } else {
          res.svd.C <- svd(f.C)
          gammaV <- tcrossprod(res.svd.C$v, res.svd.C$u)
        }

        MM.c.old <- MM.c.new
      }

      m.iter2 <- 0
      check2 <- 0
      MM.c.old <- -Inf

      for (j in 1:k) {
        deltaV.k[, , j] <- diag(diag(crossprod(gammaV, WC[, , j]) %*% gammaV), r, r) / (det(diag(diag(crossprod(gammaV, WC[, , j]) %*% gammaV), r, r)))^(1 / r)
      }

      for (j in 1:k) {
        sigmaV[, , j] <- gammaV %*% tcrossprod(deltaV.k[, , j], gammaV)
      }
    }

    if (mod.col == "EV") {
      for (j in 1:k) {
        tempW_EV[[j]] <- eigen(WC[, , j])

        gammaV.k[, , j] <- tempW_EV[[j]][["vectors"]]

        tempomegaC[, , j] <- diag(tempW_EV[[j]][["values"]], r, r)
      }

      deltaV <- rowSums(tempomegaC, dims = 2) / ((det(rowSums(tempomegaC, dims = 2)))^(1 / r))

      for (j in 1:k) {
        sigmaV[, , j] <- gammaV.k[, , j] %*% tcrossprod(deltaV, gammaV.k[, , j])
      }
    }

    if (mod.col == "VV") {
      for (j in 1:k) {
        sigmaV[, , j] <- WC[, , j] / ((det(WC[, , j]))^(1 / r))
      }
    }

    prior <- colMeans(matrix(post[1:num, ], nrow = num, ncol = k)) # initial prob.
    post2a <- as.matrix(apply(post2, c(1, 2), sum))
    PI <- diag(1 / rowSums(post2a), nrow = k, ncol = k) %*% post2a # Transition probability matrix

    for (j in 1:k) {
      dens[, j] <- dMVnorm(X = Xresh, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j])
    }

    f <- array(dens, c(num, t, k))
    A[, 1, ] <- matrix(rep(prior, each = num), ncol = k) * f[, 1, ]

    B[, t, ] <- 1
    for (T in 2:t) {
      A[, T, ] <- A[, T - 1, ] %*% PI * f[, T, ]
      tp <- t - T + 1
      B[, tp, ] <- t(PI %*% t(B[, tp + 1, ] * f[, tp + 1, ]))
    }

    part <- rowSums(matrix(A[, t, ], nrow = num, ncol = k))
    part[part < 5e-324] <- 5e-324

    loglik.new <- sum(log(part)) # log-likelihood

    ll <- c(ll, loglik.new)

    # stopping rule

    if ((loglik.new - loglik.old) < tol) {
      check <- 1
    }

    if (loglik.new < loglik.old) {
      mark <- 1 ## Bad situation
    } else {
      mark <- 0
    } ## Good situation

    if (m.iter == maxit) {
      check <- 1
    }

    loglik.old <- loglik.new
  }

  #### Output ####

  if (mark == 1) {
    return(NA)
  } else {

    # --------------------- #
    # Classification Matrix #
    # --------------------- #

    group <- apply(post, 1, which.max)
    groupT <- array(group, c(num, t), dimnames = list(1:num, paste("time", 1:t, sep = " ")))

    # -------------------- #
    # Information criteria #
    # -------------------- #

    # Number of parameters

    meanpar <- (p * r) * k

    if (mod.row == "EII") {
      rowpar <- 1
    }
    if (mod.row == "VII") {
      rowpar <- k
    }
    if (mod.row == "EEI") {
      rowpar <- p
    }
    if (mod.row == "VEI") {
      rowpar <- k + (p - 1)
    }
    if (mod.row == "EVI") {
      rowpar <- 1 + k * (p - 1)
    }
    if (mod.row == "VVI") {
      rowpar <- k * p
    }
    if (mod.row == "EEE") {
      rowpar <- p * (p + 1) / 2
    }
    if (mod.row == "VEE") {
      rowpar <- k - 1 + p * (p + 1) / 2
    }
    if (mod.row == "EVE") {
      rowpar <- 1 + k * (p - 1) + p * (p - 1) / 2
    }
    if (mod.row == "VVE") {
      rowpar <- k * p + p * (p - 1) / 2
    }
    if (mod.row == "EEV") {
      rowpar <- p + k * p * (p - 1) / 2
    }
    if (mod.row == "VEV") {
      rowpar <- k + (p - 1) + (k * p * (p - 1) / 2)
    }
    if (mod.row == "EVV") {
      rowpar <- 1 + k * (p * ((p + 1) / 2) - 1)
    }
    if (mod.row == "VVV") {
      rowpar <- k * p * (p + 1) / 2
    }

    if (mod.col == "II") {
      colpar <- 0
    }
    if (mod.col == "EI") {
      colpar <- r - 1
    }
    if (mod.col == "VI") {
      colpar <- k * (r - 1)
    }
    if (mod.col == "EE") {
      colpar <- r * ((r + 1) / 2) - 1
    }
    if (mod.col == "VE") {
      colpar <- k * (r - 1) + r * (r - 1) / 2
    }
    if (mod.col == "EV") {
      colpar <- (r - 1) + k * r * (r - 1) / 2
    }
    if (mod.col == "VV") {
      colpar <- k * (r * ((r + 1) / 2) - 1)
    }

    weights <- k - 1
    pipar <- k * (k - 1)

    npar <- meanpar + rowpar + colpar + weights + pipar

    name <- c(mod.row, mod.col)

    # to be minimized

    BIC <- -2 * loglik.new + npar * log(num * t)

    ptm2 <- proc.time() - ptm
    time <- ptm2[3]

    return(list(
      name = name, prior = prior, M = M, sigmaU = sigmaU, sigmaV = sigmaV, PI = round(PI, digits = 3),
      A = A, B = B, loglik = loglik.new, ll = ll, mark = mark, check = check, npar = npar, iter = m.iter, time = time, BIC = BIC, groupT = groupT
    ))
  }
}
