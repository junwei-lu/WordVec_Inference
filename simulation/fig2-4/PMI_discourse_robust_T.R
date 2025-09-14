library(pracma)
library(MASS)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(latex2exp)
library(reshape2)

###############################################################################
## Original data_gene() for reference
###############################################################################
data_gene <- function(tm,p,alpha,V){
  # V: d*p, p=[log(d)^2]
  z <- rep(0,p)
  to_seq <- rep(0,tm)
  r <- matrix(rnorm(tm*p, sd = 1/sqrt(p)), tm, p)
  for(i in 1:tm){
    z <- sqrt(alpha)*z + sqrt(1-alpha)*r[i,]  # p*1
    se <- exp(V %*% as.matrix(z))             # d*1
    to_seq[i] <- sample(1:d, 1, replace = TRUE, prob = se / sum(se))
  }
  return(to_seq)
}

###############################################################################
## The two alternative generators
###############################################################################
data_gene_unit <- function(tm,p,alpha,V){
  # V: d*p
  # p = floor(log(d)^2) + 1 in your code
  z <- rep(0,p)
  to_seq <- rep(0,tm)
  r <- matrix(rnorm(tm*p, sd = 1/sqrt(p)), tm, p)
  for(i in 1:tm){
    z <- sqrt(alpha)*z + sqrt(1-alpha)*r[i,]
    # Normalize z to have unit length
    z_norm <- sqrt(sum(z^2))
    if(z_norm > 1e-12) {
      z <- z / z_norm
    }
    se <- exp(V %*% as.matrix(z))
    to_seq[i] <- sample(1:d, 1, replace = TRUE, prob = se / sum(se))
  }
  return(to_seq)
}

data_gene_arma <- function(tm, p, alpha, V, d, theta1, theta2, theta3) {
  # AR(1) coefficient: phi = sqrt(alpha)
  # noise std. factor = sqrt(1 - alpha)
  # MA(3) coefficients: theta1, theta2, theta3
  z <- rep(0, p)  # Current z state (p-dimensional)
  e <- matrix(0, nrow = tm, ncol = p)  # Store noise for each time
  to_seq <- numeric(tm)
  r <- matrix(rnorm(tm * p, sd = 1 / sqrt(p)), nrow = tm, ncol = p)
  phi <- sqrt(alpha)
  for (i in seq_len(tm)) {
    # Current noise
    e[i, ] <- sqrt(1 - alpha) * r[i, ]
    # AR(1) part
    z <- phi * z
    # Add the new noise (MA(0) part)
    z <- z + e[i, ]
    # MA(3) parts
    if (i > 1) z <- z + theta1 * e[i - 1, ]
    if (i > 2) z <- z + theta2 * e[i - 2, ]
    if (i > 3) z <- z + theta3 * e[i - 3, ]
    
    se <- exp(V %*% z)
    probs <- se / sum(se)
    to_seq[i] <- sample(seq_len(d), size = 1, prob = probs)
  }
  return(to_seq)
}

###############################################################################
## Other helper functions
###############################################################################
coocur_cal <- function(wd_seq, d ,q = 1){ # Input one patient's time series data, output the co-occurrence matrix
  co_m <- matrix(0,d,d)
  l <- length(wd_seq)
  for (i in 1:(l-q)) {
    for (j in 1:min(q,l-i)) {
      a = wd_seq[i]
      b = wd_seq[i+j]
      if(a == b){
        co_m[a,a] <- co_m[a,a] + 2 
      } else {
        co_m[a,b] <- co_m[a,b] + 1
        co_m[b,a] <- co_m[b,a] + 1
      }
    }
  }
  return(co_m)
}

SPPMI_calc_lr <- function(co, q = 1, d, len, n, p){
  # Compute the PMI low-rank estimator
  r <- rowSums(co)
  c <- (2*q*len - 2*q^2)*n  # total co-occurrence
  SPPMI <- matrix(0,d,d)
  for (i in 1:d) {
    for(j in 1:i){
      SPPMI[i,j] = log(c*co[i,j] / (r[i]*r[j]))
      SPPMI[j,i] = SPPMI[i,j]
    }
  }
  re <- svd(SPPMI, nu = p, nv = p)
  L = re$d[1:p]
  PMI_lr = re$u %*% diag(L) %*% t(re$v)
  return(PMI_lr)
}

SPPMI_calc <- function(co, q = 1, d, len, n){
  # Compute the (symmetric) PMI
  r <- rowSums(co)
  c <- (2*q*len - 2*q^2)*n
  if (any(co <= 1)) {
    message("Warning: The matrix 'co' contains entries less than or equal to 1.")
  }
  SPPMI <- matrix(0,d,d)
  for (i in 1:d) {
    for(j in 1:i){
      SPPMI[i,j] = log(c*co[i,j] / (r[i]*r[j]))
      SPPMI[j,i] = SPPMI[i,j]
    }
  }
  return(SPPMI)
}

prob_est <- function(mc_steps, p, V){
  z <- matrix(rnorm(mc_steps * p), mc_steps, p)/sqrt(p)
  p_est <- rep(0, d)
  for(i in 1:mc_steps){
    tmp <- exp(V %*% as.matrix(z[i,]))
    p_est = p_est + tmp / sum(tmp)
  }
  return (p_est / mc_steps)
}

pwwu_est <- function(mc_steps, p, V, u, alpha){
  z <- matrix(rnorm(mc_steps * p), mc_steps, p)/sqrt(p)
  r <- matrix(rnorm(mc_steps * p), mc_steps, p)/sqrt(p)
  # zu = alpha^(u/2)*z + sqrt(1 - alpha^u)*r
  # But be mindful of operator precedence; in R: alpha^(u/2) is alpha^(u/2)
  # We'll do the same as your code
  zu <- alpha^(u/2)*z + sqrt(1 - alpha^u)*r
  pmatrix_est <- matrix(0,d,d)
  for(i in 1:mc_steps){
    tmp1 <- exp(V %*% as.matrix(z[i,]))
    tmp1 <- tmp1/sum(tmp1)
    tmp2 <- exp(V %*% as.matrix(zu[i,]))
    tmp2 <- tmp2/sum(tmp2)
    pmatrix_est <- pmatrix_est + as.matrix(tmp1) %*% t(tmp2)
  }
  pmatrix_est <- pmatrix_est / mc_steps
  # Symmetrize
  pmatrix_est <- (pmatrix_est + t(pmatrix_est)) / 2
  return(pmatrix_est)
}

###############################################################################
## Setting up the environment
###############################################################################
set.seed(111)
d <- 100
p <- floor(log(d)^2) + 1  # dimension of word vector
q <- 5                    # window size
mc_steps <- 1e7
alpha <- 1 - log(d)/(p^2)

# We'll create an initial V
ord <- 0.5
kpa <- d^(-ord)
U <- orth(randn(n = d, m = p + 2))[, 1:p]
L <- kpa * diag(p)
V_init <- U %*% L   # d x p
# Centering step
p_est_init <- prob_est(mc_steps, p, V_init)
mu_center <- t(V_init) %*% p_est_init
V <- V_init - as.matrix(rep(1, d)) %*% t(mu_center)

# For the "true" matrix VVT2 that we compare to:
t_len <- 1000
alpha_p <- 0
for(u in 1:q){
  alpha_p <- alpha_p + alpha^(u/2)
}
alpha_p <- alpha_p/(p*q)
V_tild <- sqrt(alpha_p)*V
VVT2 <- V_tild %*% t(V_tild)

###############################################################################
## Wrapper to run the experiment for each data generation function
###############################################################################
run_experiment <- function(data_gene_fn,
                           d, p, alpha, V,
                           Tlist,
                           n, q,
                           num_runs,
                           # For data_gene_arma, supply extra ARMA params, otherwise ignored
                           theta1 = 0.2, theta2 = 0.1, theta3 = 0.05) {
  
  # Storage for the repeated runs
  results_d1 <- matrix(0, nrow = num_runs, ncol = length(Tlist))
  results_d2 <- matrix(0, nrow = num_runs, ncol = length(Tlist))
  results_e1 <- matrix(0, nrow = num_runs, ncol = length(Tlist))
  results_e2 <- matrix(0, nrow = num_runs, ncol = length(Tlist))
  
  for(run in 1:num_runs) {
    cat("Run =", run, "\n")
    d1 <- rep(0, length(Tlist))
    d2 <- rep(0, length(Tlist))
    e1 <- rep(0, length(Tlist))
    e2 <- rep(0, length(Tlist))
    
    i <- 1
    for(t_len in Tlist) {
      total_seq <- matrix(0, nrow = n, ncol = t_len)
      
      # Generate data
      for(l in 1:n){
        if( identical(data_gene_fn, data_gene_arma) ){
          # For ARMA, we need extra arguments
          total_seq[l, ] <- data_gene_fn(tm = t_len, p = p, alpha = alpha, V = V, 
                                         d = d, theta1 = theta1, theta2 = theta2, theta3 = theta3)
        } else {
          # For data_gene or data_gene_unit
          total_seq[l, ] <- data_gene_fn(t_len, p, alpha, V)
        }
      }
      
      # Build co-occurrence
      co <- matrix(0, d, d)
      for (j in 1:n) {
        tmp2 <- total_seq[j,] 
        tmp3 <- coocur_cal(tmp2, d, q)
        co <-  co + tmp3
      }
      
      # SPPMI
      SP_est <- SPPMI_calc(co, q, d, t_len, n)
      PMI_lr <- SPPMI_calc_lr(co, q, d, t_len, n, p)
      
      # Evaluate distance
      d1[i] <- norm(VVT2 - SP_est, "F")  # Frobenius
      d2[i] <- norm(VVT2 - PMI_lr, "F")  # Frobenius
      e1[i] <- norm(VVT2 - SP_est, "M")  # Max norm
      e2[i] <- norm(VVT2 - PMI_lr, "M")  # Max norm
      
      i <- i + 1
    }
    results_d1[run, ] <- d1
    results_d2[run, ] <- d2
    results_e1[run, ] <- e1
    results_e2[run, ] <- e2
  }
  
  # Return in a convenient list
  list(
    d1 = results_d1,
    d2 = results_d2,
    e1 = results_e1,
    e2 = results_e2
  )
}

###############################################################################
## Running the experiment for each generator
###############################################################################
n=1000
Tlist <- c(100, 200, 400, 800, 1000, 2000, 4000, 8000)
num_runs <- 50

# 1) Original data_gene
res_original <- run_experiment(
  data_gene_fn = data_gene,
  d = d, p = p, alpha = alpha, V = V,
  Tlist = Tlist, n = n, q = q, num_runs = num_runs
)

# 2) data_gene_unit
res_unit <- run_experiment(
  data_gene_fn = data_gene_unit,
  d = d, p = p, alpha = alpha, V = V,
  Tlist = Tlist, n = n, q = q, num_runs = num_runs
)

# 3) data_gene_arma
# You can choose your own theta1, theta2, theta3
res_arma <- run_experiment(
  data_gene_fn = data_gene_arma,
  d = d, p = p, alpha = alpha, V = V,
  Tlist = Tlist, n = n, q = q, num_runs = num_runs,
  theta1 = 0.2, theta2 = 0.1, theta3 = 0.05
)

###############################################################################
## A small helper function for quantiles
###############################################################################
get_quantiles <- function(M) {
  qt20  <- apply(M, 2, quantile, probs = 0.2)
  qt80  <- apply(M, 2, quantile, probs = 0.8)
  med   <- apply(M, 2, median)
  list(median = med, qt20 = qt20, qt80 = qt80)
}

###############################################################################
## Combine results into a single data frame for plotting
## We'll illustrate the max-norm ("e" values) as in your original code
###############################################################################
make_plot_data <- function(res_list, Tlist, label){
  # res_list is something like res_original$e1 or $e2, etc.
  quants <- get_quantiles(res_list)
  data.frame(
    Length = Tlist,
    Median     = quants$median,
    Lower      = quants$qt20,
    Upper      = quants$qt80,
    Estimator  = label
  )
}

# We illustrate combining e1 & e2 (Max-norm) for each generator
# e1 = SPPMI_calc  =>  \widehat{PMI}
# e2 = SPPMI_calc_lr => \tilde{PMI}

###############################################################################
## Build a unified data frame
###############################################################################
build_generator_df <- function(results, generator_name, Tlist){
  # results$e1 => widehat{PMI}
  # results$e2 => tilde{PMI}
  # We'll also add a third "Embedding" if you want to replicate the same style 
  #   as your original code. Here you need to decide how to define "Embedding"
  #   in your numeric experiment. In your original code, "Embedding" was 
  #   something like: rep(norm(PMI_tr - VVT2,"M"), length(Tlist)), but PMI_tr 
  #   was never defined. You may omit or define it in a suitable way.
  # Here, let's just illustrate two lines (widehat{PMI} & tilde{PMI}) as an example.

  qd1 <- make_plot_data(results$e1, Tlist, "\\widehat{PMI}")
  qd2 <- make_plot_data(results$e2, Tlist, "\\tilde{PMI}")
  
  # Combine
  df <- rbind(qd1, qd2)
  df$Generator <- generator_name
  return(df)
}

df_orig <- build_generator_df(res_original, "AR(1) (original)", Tlist)
df_unit <- build_generator_df(res_unit,     "Unit-sphere", Tlist)
df_arma <- build_generator_df(res_arma,     "ARMA(1,3)", Tlist)

plot_data_all <- rbind(df_orig, df_unit, df_arma)
csvfile_name = paste0('/home/zhiweixu/pmi/PMI_discourse_robust_plot_data0107_d',d,'_ord',abs(ord),'.csv')
write.csv(plot_data_all, file = csvfile_name, row.names = FALSE)
