# Code taken from github::e3bo/spaero@tdar   spaero::get_stats()

get_tdar <- function(x){

    f1 <- function(n) return(rep(1, length.out = n))
    f2 <- function(n) return(seq(from = -n / 2 + 0.5, to = n / 2 - 0.5))
    ar_values <- vector()
    for (i in 3:length(x)) { ## TODO What if x is of length 2 or less, or a matrix? Why is lowest i 3 and not 2?
      current_data <- x[1:i]
      new_calculation <- tryCatch({
        M <- get_TDAR_M(y = current_data, f1, f2)
        theta <- get_TDAR_theta(M = M, x = current_data[-1])
        get_TDAR_values(y = current_data, theta = theta, f1, f2)
      }, error = function(e) {
        replace <- 5})
      ar_values[i] <- new_calculation[i]
    }
    return(ar_values)
}

get_TDAR_M <- function(y, ...) {
  ## use Eq. 12 in BMC Bioinformatics 2008, 9(Suppl 9):S14
  basis_functions <- list(...)
  
  # Get M
  N <- length(y)
  F <- matrix(nrow = N, ncol = length(basis_functions))
  for (i in 1:length(basis_functions)) {
    F[ , i] <- basis_functions[[i]](N)
  }
  U <- as.numeric(y) * F
  Phi <- U[-N, ]
  x <- y[-1]
  M <- solve(t(Phi) %*% Phi) %*% t(Phi)
  return(M)
}

get_TDAR_theta <- function(M, x) {
  theta <- M %*% x
  return(theta)
}

get_TDAR_values <- function(y, theta, ...) {
  basis_functions <- list(...)
  N <- length(y)
  bf_matrix <- matrix(nrow = N, ncol = length(basis_functions))
  for (i in 1:length(basis_functions)) {
    bf_matrix[ , i] <- basis_functions[[i]](N)*as.numeric(theta[i])
  }
  ar_values <- rowSums(x = bf_matrix)
  
  return(ar_values)
}   

  