#' Calculate partial correlation coefficient between x and y given z
#' 
#' @param x vector, independent variable
#' @param y vector, dependent variable
#' @param z vector or matrix, control variable
#' @param use the approach used to calculate the partial correlation coefficient.
#'   This can be done either "mat", which uses a variance-covariance matrix
#'   approach, or "rec", which uses a recursive approach.
#' @param method the correlation method used. This can be "pearson", which is
#'   the default, "spearman" or "kendall".
#'   
#' na.rm: If na.rm is T, then all the missing samples are deleted from the whole
#' dataset, which is (x, y, z). If not, the missing samples will be removed just
#' when the correlation coefficient is calculated. However, the number of
#' samples for the p-value is the number of samples after removing all the
#' missing samples from the whole dataset. Default is "T".
#' 

pcor_test <- function(x, y, z, use = "mat", method = "p", na.rm = TRUE){
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  
  if(use == "mat"){
    p.use <- "Var-Cov matrix"
    pcor.out <- pcor_mat(x, y, z, method = method, na.rm = na.rm)
  } else if(use == "rec"){
    p.use <- "Recursive formula"
    pcor.out <- pcor_rec(x, y, z, method = method, na.rm = na.rm)
  } else {
    stop("use should be either rec or mat!\n")
  }

  # sample number
  n <- dim(na.omit(data.frame(x, y, z)))[1]

  # given variables' number
  gn <- dim(z)[2]
  
  # xy correlation p-value
  cor.pvalue <- cor_pvalue(pcor.out$r, n, gn = 0, method)
  # partial correlation p-value
  pcor.pvalue <- cor_pvalue(pcor.out$partial.r, n, gn, method)
  
  return(data.frame(n = n, partial.r = pcor.out$partial.r, 
    partial.pvalue = pcor.pvalue, r = pcor.out$r, pvalue = cor.pvalue,
    gn = gn, Method = method, Use = p.use))
}

#' Calculate partial correlation coefficient between x and y given z using
#' variance-covariance matrix approach
#'
pcor_mat <- function(x, y, z, method = "p", na.rm = T) {

  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)

  if(dim(z)[2] == 0){
    stop("There should be given data\n")
  }

  data <- data.frame(x, y, z)

  if(na.rm == T){
    data <- na.omit(data)
  }

  xdata <- na.omit(data.frame(data[, c(1, 2)]))
  Sxx <- cov(xdata, xdata, m = method)

  xzdata <- na.omit(data)
  xdata <- data.frame(xzdata[, c(1, 2)])
  zdata <- data.frame(xzdata[, -c(1, 2)])
  Sxz <- cov(xdata, zdata, m = method)

  zdata <- na.omit(data.frame(data[, -c(1, 2)]))
  Szz <- cov(zdata, zdata, m = method)

  # is Szz positive definite?
  zz.ev <- eigen(Szz)$values
  if(min(zz.ev)[1] < 0){
    stop("\'Szz\' is not positive definite!\n")
  }

  # partial correlation
  Sxx.z <- Sxx - Sxz %*% solve(Szz) %*% t(Sxz)

  rxx.z <- cov2cor(Sxx.z)[1, 2]

  return(data.frame(r = cov2cor(Sxx)[1, 2], partial.r = rxx.z))
}

#' Calculate partial correlation coefficient between x and y given z using
#' recursive approach
#'
pcor_rec <- function(x, y, z, method = "p", na.rm = T) {
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)

  if(dim(z)[2] == 0){
    stop("There should be given data\n")
  }

  data <- data.frame(x, y, z)

  if(na.rm == T){
    data <- na.omit(data)
  }

  # recursive formula
  if(dim(z)[2] == 1){
    tdata <- na.omit(data.frame(data[, 1], data[, 2]))
    rxy <- cor(tdata[, 1], tdata[, 2], m = method)

    tdata <- na.omit(data.frame(data[, 1], data[, -c(1, 2)]))
    rxz <- cor(tdata[, 1], tdata[, 2], m = method)

    tdata <- na.omit(data.frame(data[, 2], data[, -c(1, 2)]))
    ryz <- cor(tdata[, 1], tdata[, 2], m = method)

    rxy.z <- (rxy - rxz*ryz)/( sqrt(1 - rxz^2) * sqrt(1 - ryz^2) )

  } else {
    x <- c(data[, 1])
    y <- c(data[, 2])
    z0 <- c(data[, 3])
    zc <- as.data.frame(data[, -c(1, 2, 3)])

    rxy.zc <- pcor_rec(x, y, zc, method = method, na.rm = na.rm)
    rxz0.zc <- pcor_rec(x, z0, zc, method = method, na.rm = na.rm)
    ryz0.zc <- pcor_rec(y, z0, zc, method = method, na.rm = na.rm)

    rxy.z <- (rxy.zc - rxz0.zc * ryz0.zc) / 
              (sqrt(1 - rxz0.zc^2) * sqrt(1 - ryz0.zc^2) )
  }
  return(data.frame(r = rxy, partial.r = rxy.z))
}

#' Calculate p-value from correlation coefficient and sample size
#' 
#' @param r numeric, correlation coefficient
#' @param n numeric, sample size
#' 
cor_pvalue <- function(r, n, gn, method) {
  
  if (method == "kendall"){
    T <- r / sqrt(2 * (2 * (n - gn) + 5) / 
      (9 * (n - gn) * (n - 1 - gn)))
    
  } else {
    T <- r * sqrt((n - 2 - gn) / (1 - r^2))
  }
#   p.value <- 2 * pnorm(-abs(T))
  p.value <- 2 * pt(abs(T), n - 2 - gn, lower.tail = FALSE)
  return(p.value)
}
