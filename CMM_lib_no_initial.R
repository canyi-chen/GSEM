#library(truncnorm)
########################PROGRESS BAR#############################################################
sapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...) {
    curVal <- get("counter", envir = env)
    assign("counter", curVal + 1, envir = env)
    setTxtProgressBar(get("pb", envir = env), curVal + 1)
    FUN(...)
  }

  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}

print_time_min = function() {
  tm = substr(Sys.time(), 6, 16)
  tm = sub("-", "", tm)
  tm = sub(" ", "", tm)
  tm = sub(":", "", tm)
  return(tm)
}


########################GENERATE COVARIATES################################################
CMM_sim_covariates = function(rho, p, p1, n, sigma) {
  Z = matrix(rnorm(n, sd = sigma * sqrt(rho)), nrow = n, ncol = p - 1)
  X = matrix(rnorm(n * (p - 1), sd = sigma * sqrt(1 - rho)), nrow = n, ncol = p - 1)
  W_m = cbind(1, X + Z)
  W_y = W_m
  W_x = W_m[, 1:p1]
  return(list(W_x, W_m, W_y))
}

########################GENERATE DATA FOR A GIVEN DISTRIBUTION & GIVEN COV#############################
CMM_sim_dat = function(n, W_x, W_m, W_y, theta_mx, theta_yx, theta_ym,
                       beta_x, beta_m, beta_y, sigma_x, sigma_m, sigma_y,
                       X_distri, M_distri, Y_distri)
{

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)

  Z_x = rnorm(n)
  Z_m = theta_mx * Z_x + rnorm(n)
  Z_y = theta_yx * Z_x + theta_ym * Z_m + rnorm(n)

  Z_m_star = Z_m / std_Z_m
  Z_y_star = Z_y / std_Z_y

  eta_x = W_x %*% beta_x
  eta_m = W_m %*% beta_m
  eta_y = W_y %*% beta_y

  #generate X
  if (X_distri == "Gaussian") { X = eta_x + sigma_x * Z_x }

  #generate M
  if (M_distri == "Bernoulli")
  {
    prob = exp(eta_m) / (1 + exp(eta_m))
    M = as.numeric(pnorm(Z_m / std_Z_m) > 1 - prob)
  }
  if (M_distri == "Gaussian") {
    M = eta_m + sigma_m * Z_m_star
  }

  #generate Y
  if (Y_distri == "Poisson")
  {

    phi_y = pnorm(Z_y / std_Z_y)
    Y = rep(NA, length = length(phi_y))
    for (i in 1:length(phi_y)) {
      k = 1
      F_y = ppois(k - 1, exp(eta_y[i]))
      while (phi_y[i] >= F_y)
      {
        k = k + 1
        F_y = ppois(k - 1, exp(eta_y[i]))
      }
      Y[i] = k - 1 }
  }
  if (Y_distri == "Gaussian") {
    Y = eta_y + sigma_y * Z_y_star
  }

  if (Y_distri == "Bernoulli")
  {
    prob = exp(eta_y) / (1 + exp(eta_y))
    Y = as.numeric(pnorm(Z_y_star) > 1 - prob)
  }

  theta_beta = list(theta = c(theta_mx, theta_yx, theta_ym),
                    beta_x = beta_x,
                    beta_m = beta_m,
                    beta_y = beta_y,
                    sigma_x = sigma_x,
                    sigma_m = sigma_m,
                    sigma_y = sigma_y)

  return(list(
    X = X,
    M = M,
    Y = Y,
    Z_x = Z_x,
    Z_m_star = Z_m_star,
    Z_y_star = Z_y_star,
    W_x = W_x,
    W_y = W_y,
    W_m = W_m,
    X_distri = X_distri,
    M_distri = M_distri,
    Y_distri = Y_distri,
    theta_beta = theta_beta
  ))
}


##############################GLM FITTING TO GET BETA EST###################################################
CMM_glm_fit = function(X, M, Y, W_x, W_m, W_y, X_distri, M_distri, Y_distri)
{

  ############X distribution############################
  if (X_distri == "Gaussian") {
    mod_x = lm(X ~ W_x + 0) #fit a linear model for X
    eta_x = mod_x$fitted.values
    sigma_x = sigma(mod_x)

    Z_x = (X - eta_x) / sigma_x
    pi_Z_x = dnorm(Z_x, log = TRUE) #take the log of dnorm(Z_x)

    X_res = list(distri = X_distri,
                 Z_x = Z_x, eta_x = eta_x,
                 pi_Z_x = pi_Z_x, beta_x = mod_x$coef, sigma_x = sigma_x)
  }


  ############M distribution############################
  if (M_distri == "Bernoulli") {
    mod_m = glm(M ~ W_m + 0, family = binomial(link = 'logit')) #fit a logistic model for M
    p = mod_m$fitted.values #pi
    lm = qnorm(ifelse(M == 0, 0, 1 - p))
    um = qnorm(ifelse(M == 0, 1 - p, 1))

    M_res = list(distri = M_distri, p = p,
                 lm = lm,
                 um = um, beta_m = mod_m$coef)
  }


  if (M_distri == "Gaussian") {
    mod_m = lm(M ~ W_m + 0) #fit a linear model for M
    eta_m = mod_m$fitted.values
    sigma_m = sigma(mod_m)

    Z_m_star = (M - eta_m) / sigma_m

    M_res = list(distri = M_distri,
                 Z_m_star = Z_m_star, eta_m = eta_m,
                 beta_m = mod_m$coef, sigma_m = sigma_m)
  }


  ############Y distribution############################

  if (Y_distri == "Gaussian") {
    mod_y = lm(Y ~ W_y + 0) #fit a linear model for Y
    eta_y = mod_y$fitted.values
    sigma_y = sigma(mod_y)
    Z_y_star = (Y - eta_y) / sigma_y

    Y_res = list(distri = Y_distri,
                 Z_y_star = Z_y_star, eta_y = eta_y,
                 beta_y = mod_y$coef, sigma_y = sigma_y)

  }

  if (Y_distri == "Bernoulli") {
    mod_y = glm(Y ~ W_y + 0, family = binomial(link = 'logit')) #fit a logistic model for M
    p = mod_y$fitted.values #pi
    ly = qnorm(ifelse(Y == 0, 0, 1 - p))
    uy = qnorm(ifelse(Y == 0, 1 - p, 1))

    Y_res = list(distri = Y_distri, p = p,
                 ly = ly,
                 uy = uy,
                 Z_y_star = 0.5 * (ly + uy),
                 beta_y = mod_y$coef)
  }

  glm_res = list(X_res = X_res,
                 M_res = M_res,
                 Y_res = Y_res)


  return(glm_res)
}


###########################COMPUTE LIKELIHOOD FUNCTION#####################################
neg_log_like_fun = function(theta, glm_res, batch = NULL)
{
  theta_mx = theta[1]
  theta_yx = theta[2]
  theta_ym = theta[3]

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)


  #1 X: CONTINUOUS(NORMAL), M: CONTINUOUS(NORMAL), Y: CONTINUOUS(NORMAL)
  if (glm_res$X_res$distri == "Gaussian" &
    glm_res$M_res$distri == "Gaussian" &
    glm_res$Y_res$distri == "Gaussian")
  {
    pi_Z_x_Z_m_star_Z_y_star = sapply(1:length(glm_res$M_res$Z_m_star), function(i)
      dmvnorm(c(glm_res$X_res$Z_x[i], glm_res$M_res$Z_m_star[i], glm_res$Y_res$Z_y_star[i]), mean = c(0, 0, 0),
              sigma = matrix(c(1, theta_mx / std_Z_m, alpha / std_Z_y,
                               theta_mx / std_Z_m, 1, (theta_mx * alpha + theta_ym) / (std_Z_m * std_Z_y),
                               alpha / std_Z_y, (theta_mx * alpha + theta_ym) / (std_Z_m * std_Z_y), 1), nrow = 3, ncol = 3, byrow = T),
              log = T))
    neg_log_like = -sum(pi_Z_x_Z_m_star_Z_y_star)
  }

  #2 X: CONTINUOUS(NORMAL), M: DISCRETE(BERNOULLI), Y: CONTINUOUS(NORMAL)
  if (glm_res$X_res$distri == "Gaussian" &
    glm_res$M_res$distri == "Bernoulli" &
    glm_res$Y_res$distri == "Gaussian")
  {

    pi_Z_x_Z_y_star = sapply(1:length(glm_res$M_res$lm), function(i)
      dmvnorm(c(glm_res$X_res$Z_x[i], glm_res$Y_res$Z_y_star[i]), mean = c(0, 0),
              sigma = matrix(c(1, alpha / std_Z_y, alpha / std_Z_y, 1), nrow = 2, ncol = 2, byrow = T), log = T))

    mean_Z_m_star = ((theta_mx - theta_ym * theta_yx) * glm_res$X_res$Z_x + theta_ym * std_Z_y * glm_res$Y_res$Z_y_star) / (std_Z_m * (theta_ym^2 + 1))
    std_Z_m_star = 1 / sqrt((theta_ym^2 + 1) * (theta_mx^2 + 1))

    prob_Z_m_star_Z_x_Z_y_star = sapply(1:length(glm_res$M_res$lm), function(i) {
      pnorm(glm_res$M_res$um[i], mean = mean_Z_m_star[i], sd = std_Z_m_star) - pnorm(glm_res$M_res$lm[i], mean = mean_Z_m_star[i], sd = std_Z_m_star)
    })

    neg_log_like_vec = pi_Z_x_Z_y_star + log(prob_Z_m_star_Z_x_Z_y_star)
    neg_log_like = -sum(neg_log_like_vec)

  }

  #3 X: CONTINUOUS(NORMAL), M: CONTINUOUS(NORMAL), Y: BERNOULLI
  if (glm_res$X_res$distri == "Gaussian" &
    glm_res$M_res$distri == "Gaussian" &
    glm_res$Y_res$distri == "Bernoulli")
  {
    pi_Z_x_Z_m_star = sapply(1:length(glm_res$Y_res$ly), function(i)
      dmvnorm(c(glm_res$X_res$Z_x[i], glm_res$M_res$Z_m_star[i]), mean = c(0, 0),
              sigma = matrix(c(1, theta_mx / std_Z_m, theta_mx / std_Z_m, 1), nrow = 2, ncol = 2, byrow = T), log = T))

    mean_Z_y_star = (theta_yx * glm_res$X_res$Z_x + std_Z_m * theta_ym * glm_res$M_res$Z_m_star) / std_Z_y
    std_Z_y_star = 1 / std_Z_y

    prob_Z_y_star_Z_x_Z_m_star = sapply(1:length(glm_res$Y_res$ly), function(i) {
      pnorm(glm_res$Y_res$uy[i], mean = mean_Z_y_star[i], sd = std_Z_y_star) - pnorm(glm_res$Y_res$ly[i], mean = mean_Z_y_star[i], sd = std_Z_y_star)
    })

    neg_log_like_vec = pi_Z_x_Z_m_star + log(prob_Z_y_star_Z_x_Z_m_star)
    neg_log_like = -sum(neg_log_like_vec)
  }


  return(neg_log_like)
}


################COMBINE ALL THETA AND BETA ESTIMATES#####################################
CMM_theta_beta = function(dat) {
  glm_res = CMM_glm_fit(dat$X, dat$M, dat$Y, dat$W_x, dat$W_m, dat$W_y, X_distri = dat$X_distri,
                        M_distri = dat$M_distri, Y_distri = dat$Y_distri)
  #theta0=comp_initial(dat=dat,glm_res=glm_res)
  res_opt = optim(par = c(0, 0, 0), neg_log_like_fun, glm_res = glm_res, method = "Nelder-Mead", control = list(trace = 0))
  theta = res_opt$par
  return(list(theta = theta,
              beta_x = glm_res$X_res$beta_x,
              beta_m = glm_res$M_res$beta_m,
              beta_y = glm_res$Y_res$beta_y,
              sigma_x = glm_res$X_res$sigma_x,
              sigma_m = glm_res$M_res$sigma_m,
              sigma_y = glm_res$Y_res$sigma_y))
}


neg_log_like_fun_rho = function(theta, glm_res, batch = NULL, rho)
{
  theta_mx = theta[1]
  theta_yx = theta[2]
  theta_ym = theta[3]

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)


  #2 X: CONTINUOUS(NORMAL), M: DISCRETE(BERNOULLI), Y: CONTINUOUS(NORMAL)
  if (glm_res$X_res$distri == "Gaussian" &
    glm_res$M_res$distri == "Bernoulli" &
    glm_res$Y_res$distri == "Gaussian")
  {

    pi_Z_x_Z_y_star = sapply(1:length(glm_res$M_res$lm), function(i)
      dmvnorm(c(glm_res$X_res$Z_x[i], glm_res$Y_res$Z_y_star[i]), mean = c(0, 0),
              sigma = matrix(c(1, alpha / std_Z_y, alpha / std_Z_y, 1), nrow = 2, ncol = 2, byrow = T), log = T))

    mean_Z_m_star = ((theta_mx - theta_ym * theta_yx - alpha * rho) * glm_res$X_res$Z_x + (theta_ym + rho) * std_Z_y * glm_res$Y_res$Z_y_star) / (std_Z_m * (theta_ym^2 + 1))
    std_Z_m_star = sqrt(1 - rho^2 - 2 * theta_ym * rho) / sqrt((theta_ym^2 + 1) * (theta_mx^2 + 1))

    prob_Z_m_star_Z_x_Z_y_star = sapply(1:length(glm_res$M_res$lm), function(i) {
      pnorm(glm_res$M_res$um[i], mean = mean_Z_m_star[i], sd = std_Z_m_star) - pnorm(glm_res$M_res$lm[i], mean = mean_Z_m_star[i], sd = std_Z_m_star)
    })

    neg_log_like_vec = pi_Z_x_Z_y_star + log(prob_Z_m_star_Z_x_Z_y_star)
    neg_log_like = -sum(neg_log_like_vec)

  }


  return(neg_log_like)
}

CMM_theta_beta_rho = function(dat, rho) {
  glm_res = CMM_glm_fit(dat$X, dat$M, dat$Y, dat$W_x, dat$W_m, dat$W_y, X_distri = dat$X_distri,
                        M_distri = dat$M_distri, Y_distri = dat$Y_distri)
  #theta0=comp_initial(dat=dat,glm_res=glm_res)
  res_opt = optim(par = c(0, 0, 0), neg_log_like_fun_rho, glm_res = glm_res, method = "Nelder-Mead", control = list(trace = 0), rho = rho)
  theta = res_opt$par
  return(list(theta = theta,
              beta_x = glm_res$X_res$beta_x,
              beta_m = glm_res$M_res$beta_m,
              beta_y = glm_res$Y_res$beta_y,
              sigma_x = glm_res$X_res$sigma_x,
              sigma_m = glm_res$M_res$sigma_m,
              sigma_y = glm_res$Y_res$sigma_y))
}


#############GIVEN X0, X1, AND PARAS, CALCULATE THE NDE/NIE FOR AN AVERAGE PERSON (MEAN COV)######################
effect_cal_NDE = function(theta_beta,
                          X_distri, M_distri, Y_distri,
                          mean_W_x, mean_W_m, mean_W_y,
                          x0, x1, tol = 1e-30)
{
  theta_mx = theta_beta$theta[1]
  theta_yx = theta_beta$theta[2]
  theta_ym = theta_beta$theta[3]
  beta_x = theta_beta$beta_x
  beta_m = theta_beta$beta_m
  beta_y = theta_beta$beta_y
  sigma_x = theta_beta$sigma_x
  sigma_m = theta_beta$sigma_m
  sigma_y = theta_beta$sigma_y

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)

  #1 X: CONTINUOUS, M: CONTINUOUS, Y: CONTINUOUS (Normal)
  if (X_distri == "Gaussian" &
    M_distri == "Gaussian" &
    Y_distri == "Gaussian")
  {

    Z_x0 = (x0 - sum(beta_x * mean_W_x)) / sigma_x
    Z_x1 = (x1 - sum(beta_x * mean_W_x)) / sigma_x

    NDE = sigma_y * theta_yx * (x1 - x0) / sigma_x / std_Z_y

  }

  #2 X: CONTINUOUS (NORMAL), M: BERNOULLI, Y: CONTINUOUS (NORMAL)
  if (X_distri == "Gaussian" &
    M_distri == "Bernoulli" &
    Y_distri == "Gaussian")
  {
    Z_x0 = (x0 - sum(beta_x * mean_W_x)) / sigma_x
    Z_x1 = (x1 - sum(beta_x * mean_W_x)) / sigma_x

    p = 1 / (1 + exp(-sum(beta_m * mean_W_m)))
    qm = qnorm(1 - p)
    qm_Z_x1 = std_Z_m * qm - theta_mx * Z_x1
    qm_Z_x0 = std_Z_m * qm - theta_mx * Z_x0

    p_Z_x0 = pnorm(qm_Z_x0)
    p_Z_x1 = pnorm(qm_Z_x1)
    c_Z_x1 = (theta_mx * Z_x1 * p_Z_x1 - dnorm(qm_Z_x1)) / std_Z_m

    ratio_p_Z_x0_x1 = p_Z_x0 / p_Z_x1 * c_Z_x1
    ratio_1_p_Z_x0_x1 = (1 - p_Z_x0) / (1 - p_Z_x1) * (theta_mx * Z_x1 / std_Z_m - c_Z_x1)


    NDE = sigma_y * theta_yx * (Z_x1 - Z_x0) / std_Z_y +
      sigma_y * theta_ym * std_Z_m / std_Z_y *
        (ratio_p_Z_x0_x1 + ratio_1_p_Z_x0_x1 - theta_mx * Z_x0 / std_Z_m)

  }

  #3 X: CONTINUOUS (NORMAL), M: CONTINUOUS (NORMAL), Y: BERNOULLI

  if (X_distri == "Gaussian" &
    M_distri == "Gaussian" &
    Y_distri == "Bernoulli")
  {
    Z_x0 = (x0 - sum(beta_x * mean_W_x)) / sigma_x
    Z_x1 = (x1 - sum(beta_x * mean_W_x)) / sigma_x

    p = 1 / (1 + exp(-sum(beta_y * mean_W_y)))
    qy = qnorm(1 - p)

    Z_m_star = rnorm(1000000, mean = theta_mx * Z_x0 / std_Z_m, sd = 1 / std_Z_m)
    Phi_z_x1 = pnorm(std_Z_y * qy -
                       theta_yx * Z_x1 -
                       std_Z_y * theta_ym * Z_m_star)
    Phi_z_x0 = pnorm(std_Z_y * qy -
                       theta_yx * Z_x0 -
                       std_Z_y * theta_ym * Z_m_star)
    NDE = mean(Phi_z_x0 - Phi_z_x1)
  }


  #4 X: BERNOULLI, M: CONTINUOUS (NORMAL), Y: CONTINUOUS (NORMAL)
  if (X_distri == "Bernoulli" &
    M_distri == "Gaussian" &
    Y_distri == "Gaussian")
  {
    p = 1 / (1 + exp(-sum(beta_x * mean_W_x)))
    qx = qnorm(1 - p)

    Z_m_star = rnorm(1000000, mean = 0, sd = 1)
    beta = qx * std_Z_m - theta_mx * Z_m_star
    Phi_beta = dnorm(beta) / (1 - pnorm(beta))
    #Phi_beta = zeta(-beta,1)
    NDE = sigma_y * theta_yx / std_Z_y / std_Z_m / pnorm(qx) * mean(Phi_beta)
  }
  names(NDE) = NULL
  return(NDE)

}

effect_cal_NIE = function(theta_beta, X_distri, M_distri, Y_distri, mean_W_x, mean_W_m, mean_W_y, x0, x1, tol = 1e-30)
{
  theta_mx = theta_beta$theta[1]
  theta_yx = theta_beta$theta[2]
  theta_ym = theta_beta$theta[3]
  beta_x = theta_beta$beta_x
  beta_m = theta_beta$beta_m
  beta_y = theta_beta$beta_y
  sigma_x = theta_beta$sigma_x
  sigma_m = theta_beta$sigma_m
  sigma_y = theta_beta$sigma_y

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)

  #1 X: CONTINUOUS(Normal), M: CONTINUOUS(Normal), Y: CONTINUOUS (Normal)
  if (X_distri == "Gaussian" &
    M_distri == "Gaussian" &
    Y_distri == "Gaussian")
  {

    Z_x0 = (x0 - sum(beta_x * mean_W_x)) / sigma_x
    Z_x1 = (x1 - sum(beta_x * mean_W_x)) / sigma_x

    NIE = sigma_y * theta_ym * theta_mx * (x1 - x0) / sigma_x / std_Z_y

  }


  #2 X: CONTINUOUS (NORMAL), M: CONTINUOUS (NORMAL), Y: BERNOULLI
  if (X_distri == "Gaussian" &
    M_distri == "Gaussian" &
    Y_distri == "Bernoulli")
  {
    Z_x0 = (x0 - sum(beta_x * mean_W_x)) / sigma_x
    Z_x1 = (x1 - sum(beta_x * mean_W_x)) / sigma_x

    p = 1 / (1 + exp(-sum(beta_y * mean_W_y)))
    qy = qnorm(1 - p)

    Z_m_star_0 = rnorm(1000000, mean = theta_mx * Z_x0 / std_Z_m, sd = 1 / std_Z_m)
    Z_m_star_1 = rnorm(1000000, mean = theta_mx * Z_x1 / std_Z_m, sd = 1 / std_Z_m)

    Phi_0 = mean(pnorm(std_Z_y * qy -
                         theta_yx * Z_x1 -
                         std_Z_y * theta_ym * Z_m_star_0))
    Phi_1 = mean(pnorm(std_Z_y * qy -
                         theta_yx * Z_x1 -
                         std_Z_y * theta_ym * Z_m_star_1))
    NIE = Phi_0 - Phi_1
  }


  #3. X: CONTINUOUS (NORMAL), M: BERNOULLI, Y: CONTINUOUS (NORMAL)
  if (X_distri == "Gaussian" &
    M_distri == "Bernoulli" &
    Y_distri == "Gaussian")
  {

    Z_x0 = (x0 - sum(beta_x * mean_W_x)) / sigma_x
    Z_x1 = (x1 - sum(beta_x * mean_W_x)) / sigma_x

    p = 1 / (1 + exp(-sum(beta_m * mean_W_m)))
    qm = qnorm(1 - p)
    qm_Z_x1 = std_Z_m * qm - theta_mx * Z_x1
    qm_Z_x0 = std_Z_m * qm - theta_mx * Z_x0

    p_Z_x0 = pnorm(qm_Z_x0)
    p_Z_x1 = pnorm(qm_Z_x1)
    c_Z_x1 = (theta_mx * Z_x1 * p_Z_x1 - dnorm(qm_Z_x1)) / std_Z_m


    NIE = sigma_y * theta_ym * std_Z_m / std_Z_y * (theta_mx * Z_x1 / std_Z_m -
      p_Z_x0 * c_Z_x1 / p_Z_x1 -
      (1 - p_Z_x0) * (theta_mx * Z_x1 / std_Z_m - c_Z_x1) / (1 - p_Z_x1))

  }

  #4 X: BERNOULLI, M: CONTINUOUS (NORMAL), Y: CONTINUOUS (NORMAL)
  if (X_distri == "Bernoulli" &
    M_distri == "Gaussian" &
    Y_distri == "Gaussian")
  {
    p = 1 / (1 + exp(-sum(beta_x * mean_W_x)))
    qx = qnorm(1 - p)

    Z_m_star = rnorm(1000000, mean = 0, sd = 1)
    beta = qx * std_Z_m - theta_mx * Z_m_star
    phi_qx = pnorm(qx)
    phi_int = ((theta_yx * theta_mx / std_Z_y / std_Z_m + std_Z_m * theta_ym / std_Z_y) * Z_m_star + theta_yx * dnorm(beta) /
      std_Z_y /
      std_Z_m /
      (1 - pnorm(beta))) * (phi_qx - pnorm(beta))

    NIE = sigma_y / phi_qx / (1 - phi_qx) * mean(phi_int)
  }
  names(NIE) = NULL
  return(NIE)

}


#############GIVEN X0, X1, AND PARAS, CALCULATE THE NDE/NIE FOR AN AVERAGE PERSON (MEAN COV)######################
effect_cal_TE = function(theta_beta,
                         X_distri, M_distri, Y_distri,
                         mean_W_x, mean_W_m, mean_W_y,
                         x0, x1, tol = 1e-30)
{
  theta_mx = theta_beta$theta[1]
  theta_yx = theta_beta$theta[2]
  theta_ym = theta_beta$theta[3]
  beta_x = theta_beta$beta_x
  beta_m = theta_beta$beta_m
  beta_y = theta_beta$beta_y
  sigma_x = theta_beta$sigma_x
  sigma_m = theta_beta$sigma_m
  sigma_y = theta_beta$sigma_y

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)


  #2 X: CONTINUOUS (NORMAL), M: BERNOULLI, Y: CONTINUOUS (NORMAL)
  if (X_distri == "Gaussian" &
    M_distri == "Bernoulli" &
    Y_distri == "Gaussian")
  {
    Z_x0 = (x0 - sum(beta_x * mean_W_x)) / sigma_x
    Z_x1 = (x1 - sum(beta_x * mean_W_x)) / sigma_x

    TE = sigma_y * alpha * (Z_x1 - Z_x0) / std_Z_y

  }

  names(TE) = NULL
  return(TE)

}


##############################################################################################################
###############BEGIN: COVERAGE OF EFFECTS/THETAS BASED ON PARAMETRIC BOOTSTARP#################################

########################PARAMETRIC BOOTSTRAP: FOR A GIVEN DATA, GET ONE CI FOR THE EFFECTS/THETA##############
CMM_boot_all = function(dat, theta_beta, mean_W_x, mean_W_m,
                        mean_W_y, x0, x1,
                        boot_num, prob)
{


  sim_effect_NDE = matrix(0, nrow = boot_num, ncol = 1)
  sim_effect_NIE = matrix(0, nrow = boot_num, ncol = 1)
  sim_theta = matrix(0, nrow = boot_num, ncol = length(theta_beta$theta))

  for (i in 1:boot_num) {
    sim_dat = CMM_sim_dat(n = nrow(dat$X), W_x = dat$W_x, W_m = dat$W_m, W_y = dat$W_y,
                          theta_mx = theta_beta$theta[1],
                          theta_yx = theta_beta$theta[2],
                          theta_ym = theta_beta$theta[3],
                          beta_x = theta_beta$beta_x,
                          beta_m = theta_beta$beta_m,
                          beta_y = theta_beta$beta_y,
                          sigma_x = theta_beta$sigma_x,
                          sigma_m = theta_beta$sigma_m,
                          sigma_y = theta_beta$sigma_y,
                          X_distri = dat$X_distri,
                          M_distri = dat$M_distri,
                          Y_distri = dat$Y_distri)

    sim_theta_beta = CMM_theta_beta(sim_dat)
    sim_theta[i,] = sim_theta_beta$theta

    sim_effect_NDE[i,] = effect_cal_NDE(theta_beta = sim_theta_beta,
                                        X_distri = sim_dat$X_distri,
                                        M_distri = sim_dat$M_distri,
                                        Y_distri = sim_dat$Y_distri,
                                        mean_W_x = mean_W_x,
                                        mean_W_m = mean_W_m,
                                        mean_W_y = mean_W_y,
                                        x0 = x0, x1 = x1, tol = 1e-30)

    sim_effect_NIE[i,] = effect_cal_NIE(theta_beta = sim_theta_beta,
                                        X_distri = sim_dat$X_distri,
                                        M_distri = sim_dat$M_distri,
                                        Y_distri = sim_dat$Y_distri,
                                        mean_W_x = mean_W_x,
                                        mean_W_m = mean_W_m,
                                        mean_W_y = mean_W_y,
                                        x0 = x0, x1 = x1, tol = 1e-30)

  }

  CI_sim_theta = apply(sim_theta, 2, quantile, prob = prob)
  CI_sim_effect_NDE = apply(sim_effect_NDE, 2, quantile, prob = prob)
  CI_sim_effect_NIE = apply(sim_effect_NIE, 2, quantile, prob = prob)

  return(list(CI_sim_theta = CI_sim_theta,
              CI_sim_effect_NDE = CI_sim_effect_NDE,
              CI_sim_effect_NIE = CI_sim_effect_NIE,
              sim_effect_NDE = sim_effect_NDE,
              sim_effect_NIE = sim_effect_NIE,
              sim_theta = sim_theta)
  )
}


CMM_boot_all_real = function(dat, theta_beta, mean_W_x, mean_W_m,
                             mean_W_y, x0, x1,
                             boot_num, prob)
{


  sim_effect_NDE = matrix(0, nrow = boot_num, ncol = 1)
  sim_effect_NIE = matrix(0, nrow = boot_num, ncol = 1)
  sim_effect_TE = matrix(0, nrow = boot_num, ncol = 1)
  sim_theta = matrix(0, nrow = boot_num, ncol = length(theta_beta$theta))

  for (i in 1:boot_num) {
    sim_dat = CMM_sim_dat(n = nrow(dat$X), W_x = dat$W_x, W_m = dat$W_m, W_y = dat$W_y,
                          theta_mx = theta_beta$theta[1],
                          theta_yx = theta_beta$theta[2],
                          theta_ym = theta_beta$theta[3],
                          beta_x = theta_beta$beta_x,
                          beta_m = theta_beta$beta_m,
                          beta_y = theta_beta$beta_y,
                          sigma_x = theta_beta$sigma_x,
                          sigma_m = theta_beta$sigma_m,
                          sigma_y = theta_beta$sigma_y,
                          X_distri = dat$X_distri,
                          M_distri = dat$M_distri,
                          Y_distri = dat$Y_distri)

    sim_theta_beta = CMM_theta_beta(sim_dat)
    sim_theta[i,] = sim_theta_beta$theta

    sim_effect_NDE[i,] = effect_cal_NDE(theta_beta = sim_theta_beta,
                                        X_distri = sim_dat$X_distri,
                                        M_distri = sim_dat$M_distri,
                                        Y_distri = sim_dat$Y_distri,
                                        mean_W_x = mean_W_x,
                                        mean_W_m = mean_W_m,
                                        mean_W_y = mean_W_y,
                                        x0 = x0, x1 = x1, tol = 1e-30)

    sim_effect_NIE[i,] = effect_cal_NIE(theta_beta = sim_theta_beta,
                                        X_distri = sim_dat$X_distri,
                                        M_distri = sim_dat$M_distri,
                                        Y_distri = sim_dat$Y_distri,
                                        mean_W_x = mean_W_x,
                                        mean_W_m = mean_W_m,
                                        mean_W_y = mean_W_y,
                                        x0 = x0, x1 = x1, tol = 1e-30)

    sim_effect_TE[i,] = effect_cal_TE(theta_beta = sim_theta_beta,
                                      X_distri = sim_dat$X_distri,
                                      M_distri = sim_dat$M_distri,
                                      Y_distri = sim_dat$Y_distri,
                                      mean_W_x = mean_W_x,
                                      mean_W_m = mean_W_m,
                                      mean_W_y = mean_W_y,
                                      x0 = x0, x1 = x1, tol = 1e-30)

  }

  CI_sim_theta = apply(sim_theta, 2, quantile, prob = prob)
  CI_sim_effect_NDE = apply(sim_effect_NDE, 2, quantile, prob = prob)
  CI_sim_effect_NIE = apply(sim_effect_NIE, 2, quantile, prob = prob)
  CI_sim_effect_TE = apply(sim_effect_TE, 2, quantile, prob = prob)

  return(list(CI_sim_theta = CI_sim_theta,
              CI_sim_effect_NDE = CI_sim_effect_NDE,
              CI_sim_effect_NIE = CI_sim_effect_NIE,
              CI_sim_effect_TE = CI_sim_effect_TE,
              sim_effect_NDE = sim_effect_NDE,
              sim_effect_NIE = sim_effect_NIE,
              sim_theta = sim_theta)
  )
}

#################END:COVERAGE OF EFFECTS/THETAS BASED ON PARAMETRIC BOOTSTRP################################
##############################################################################################################


##############################################################################################################
############BEGIN: COVERAGE OF EFFECTS/THETAS BASED ON NONPARAMETRIC BOOTSTRP#################################

########################NONPARAMETRIC BOOTSTRAP: RESAMPLE THE DATA FROM A GIVEN DATA#######################
CMM_sim_dat_boot_np = function(n, X, M, Y, W_x, W_m, W_y, X_distri, M_distri, Y_distri) {
  idx = sample(1:nrow(X), n, replace = T)
  X = X[idx]
  M = M[idx]
  Y = Y[idx]
  W_x = W_x[idx,]
  W_y = W_y[idx,]
  W_m = W_m[idx,]
  return(list(X = X, M = M, Y = Y, W_x = W_x, W_m = W_m, W_y = W_y, X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri))
}


########################NONPARAMETRIC BOOTSTRAP: FOR A GIVEN DATA, RESAMPLE AND GET ONE CI FOR THE EFFECTS/THETA#######################
CMM_boot_all_np = function(dat, mean_W_x, mean_W_m,
                           mean_W_y, x0, x1,
                           boot_num, prob)
{

  sim_effect_NDE = matrix(0, nrow = boot_num, ncol = 1)
  sim_effect_NIE = matrix(0, nrow = boot_num, ncol = 1)
  sim_theta = matrix(0, nrow = boot_num, ncol = 3)

  for (i in 1:boot_num) {
    sim_dat = CMM_sim_dat_boot_np(n = nrow(dat$X), X = dat$X, M = dat$M, Y = dat$Y,
                                  W_x = dat$W_x, W_m = dat$W_m, W_y = dat$W_y,
                                  X_distri = dat$X_distri,
                                  M_distri = dat$M_distri,
                                  Y_distri = dat$Y_distri)

    sim_theta_beta = CMM_theta_beta(sim_dat)
    sim_theta[i,] = sim_theta_beta$theta

    sim_effect_NDE[i,] = effect_cal_NDE(theta_beta = sim_theta_beta, X_distri = sim_dat$X_distri,
                                        M_distri = sim_dat$M_distri, Y_distri = sim_dat$Y_distri,
                                        mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                                        x0 = x0, x1 = x1, tol = 1e-30)

    sim_effect_NIE[i,] = effect_cal_NIE(theta_beta = sim_theta_beta, X_distri = sim_dat$X_distri,
                                        M_distri = sim_dat$M_distri, Y_distri = sim_dat$Y_distri,
                                        mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                                        x0 = x0, x1 = x1, tol = 1e-30)

  }

  CI_sim_theta = apply(sim_theta, 2, quantile, prob = prob)
  CI_sim_effect_NDE = apply(sim_effect_NDE, 2, quantile, prob = prob)
  CI_sim_effect_NIE = apply(sim_effect_NIE, 2, quantile, prob = prob)

  return(list(CI_sim_theta = CI_sim_theta,
              CI_sim_effect_NDE = CI_sim_effect_NDE,
              CI_sim_effect_NIE = CI_sim_effect_NIE)
  )
}


##############END:COVERAGE OF EFFECTS/THETAS BASED ON NONPARAMETRIC BOOTSTRAP################################
##############################################################################################################

sim_boot_nboot = function(seed, n_rep, n_list, theta_mx, theta_yx, theta_ym,
                          beta_x, beta_m, beta_y, sigma_x, sigma_m, sigma_y, rho_w, sigma_w,
                          X_distri, M_distri, Y_distri, save_path,
                          boot_num, prob = c(0.025, 0.975), x0, x1)
{
  set.seed(seed)

  coverage_boot_all = list()
  coverage_nboot_all = list()
  bias_all = list()
  MSE_all = list()


  for (n in n_list) {
    coverage_boot_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 5)
    coverage_nboot_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 5)
    bias_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 5)
    MSE_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 5)


    W = CMM_sim_covariates(n = n, rho = rho_w, p = length(beta_m), p1 = length(beta_x),
                           sigma = sigma_w)
    mean_W_x = apply(W[[1]], 2, mean, na.rm = FALSE)
    mean_W_m = apply(W[[2]], 2, mean, na.rm = FALSE)
    mean_W_y = apply(W[[3]], 2, mean, na.rm = FALSE)
    cat("sample size = ", n, "\n")

    theta_true_beta_true = list(theta = c(theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym),
                                beta_x = beta_x,
                                beta_m = beta_m,
                                beta_y = beta_y,
                                sigma_x = sigma_x,
                                sigma_m = sigma_m,
                                sigma_y = sigma_y)

    NDE_true_all = effect_cal_NDE(theta_beta = theta_true_beta_true, mean_W_x = mean_W_x,
                                  mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                                  X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                                  x0 = x0, x1 = x1, tol = 1e-30)


    NIE_true_all = effect_cal_NIE(theta_beta = theta_true_beta_true, mean_W_x = mean_W_x,
                                  mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                                  X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                                  x0 = x0, x1 = x1, tol = 1e-30)

    pb = txtProgressBar(style = 3)
    for (rp in 1:n_rep) {
      dat = CMM_sim_dat(n = n, W_x = W[[1]], W_m = W[[2]], W_y = W[[3]],
                        theta_mx = theta_mx,
                        theta_yx = theta_yx,
                        theta_ym = theta_ym,
                        beta_x = beta_x,
                        beta_m = beta_m,
                        beta_y = beta_y,
                        sigma_x = sigma_x,
                        sigma_m = sigma_m,
                        sigma_y = sigma_y,
                        X_distri = X_distri,
                        M_distri = M_distri,
                        Y_distri = Y_distri)

      #####CMM###############
      theta_beta_est = CMM_theta_beta(dat)
      NDE_est = effect_cal_NDE(theta_beta = theta_beta_est, mean_W_x = mean_W_x,
                               mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                               X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                               x0 = x0, x1 = x1, tol = 1e-30)


      NIE_est = effect_cal_NIE(theta_beta = theta_beta_est, mean_W_x = mean_W_x,
                               mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                               X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                               x0 = x0, x1 = x1, tol = 1e-30)

      # coverage boot
      CI_boot_all = CMM_boot_all(dat = dat, theta_beta = theta_beta_est,
                                 mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                                 x0 = x0, x1 = x1, boot_num = boot_num,
                                 prob = prob)

      coverage_boot = list(NDE = CI_boot_all$CI_sim_effect_NDE[1,] < NDE_true_all & CI_boot_all$CI_sim_effect_NDE[2,] > NDE_true_all,
                           NIE = CI_boot_all$CI_sim_effect_NIE[1,] < NIE_true_all & CI_boot_all$CI_sim_effect_NIE[2,] > NIE_true_all,
                           theta = CI_boot_all$CI_sim_theta[1,] < theta_true_beta_true$theta & CI_boot_all$CI_sim_theta[2,] > theta_true_beta_true$theta)
      coverage_boot_all[[paste(n)]][rp,] = unlist(coverage_boot)

      #coverage nboot
      CI_nboot_all = CMM_boot_all_np(dat = dat, mean_W_x = mean_W_x, mean_W_m = mean_W_m,
                                     mean_W_y = mean_W_y, x0 = x0, x1 = x1,
                                     boot_num = boot_num, prob = prob)

      coverage_nboot = list(NDE = CI_nboot_all$CI_sim_effect_NDE[1,] < NDE_true_all & CI_nboot_all$CI_sim_effect_NDE[2,] > NDE_true_all,
                            NIE = CI_nboot_all$CI_sim_effect_NIE[1,] < NIE_true_all & CI_nboot_all$CI_sim_effect_NIE[2,] > NIE_true_all,
                            theta = CI_nboot_all$CI_sim_theta[1,] < theta_true_beta_true$theta & CI_nboot_all$CI_sim_theta[2,] > theta_true_beta_true$theta)
      coverage_nboot_all[[paste(n)]][rp,] = unlist(coverage_nboot)


      # get bias
      bias_all[[paste(n)]][rp,] = c(NDE_est, NIE_est, theta_beta_est$theta) - c(NDE_true_all, NIE_true_all, unlist(theta_true_beta_true$theta))
      # get MSE
      MSE_all[[paste(n)]][rp,] = (c(NDE_est, NIE_est, theta_beta_est$theta) - c(NDE_true_all, NIE_true_all, unlist(theta_true_beta_true$theta)))^2

      setTxtProgressBar(pb, rp / n_rep)
    }
    close(pb)

    accuracy_all = list(coverage_boot = coverage_boot_all, coverage_nboot = coverage_nboot_all, bias = bias_all, MSE = MSE_all)
    save(accuracy_all, file = file.path(save_path, sprintf("accuracy_n%d_seed%d.RData", n, seed)))

  }
  return(list(accuracy_all = accuracy_all))

}


###################COMPARISON W/ PACKAGE########################################################
#######################GENERATE DATA VIA SEM####################
GEN_SEM_data = function(n, beta_xm, beta_xy, beta_my, W, beta_x, beta_m, beta_y,
                        sigma_x, sigma_m, sigma_y, X_distri, M_distri, Y_distri)
{
  W_x = W[[1]]
  W_m = W[[2]]
  W_y = W[[3]]

  if (X_distri == "Gaussian" &
    M_distri == "Gaussian" &
    Y_distri == "Gaussian")
  {
    X = matrix(W_x %*% beta_x + rnorm(n, mean = 0, sd = sigma_x), nrow = n, ncol = 1)
    M = matrix(W_m %*% beta_m + beta_xm * X + rnorm(n, mean = 0, sd = sigma_m), nrow = n, ncol = 1)
    Y = matrix(W_y %*% beta_y +
                 beta_xy * X +
                 beta_my * M +
                 rnorm(n, mean = 0, sd = sigma_y), nrow = n, ncol = 1)
  }

  if (X_distri == "Gaussian" &
    M_distri == "Bernoulli" &
    Y_distri == "Gaussian")
  {
    X = matrix(W_x %*% beta_x + rnorm(n, mean = 0, sd = sigma_x), nrow = n, ncol = 1)
    pm = 1 / (1 + exp(-W_m %*% beta_m - beta_xm * X))
    M = matrix(as.numeric(runif(n) < pm), nrow = n, ncol = 1)
    Y = matrix(W_y %*% beta_y +
                 beta_xy * X +
                 beta_my * M +
                 rnorm(n, mean = 0, sd = sigma_y), nrow = n, ncol = 1)
  }

  if (X_distri == "Gaussian" &
    M_distri == "Gaussian" &
    Y_distri == "Bernoulli")
  {
    X = matrix(W_x %*% beta_x + rnorm(n, mean = 0, sd = sigma_x), nrow = n, ncol = 1)
    M = matrix(W_m %*% beta_m + beta_xm * X + rnorm(n, mean = 0, sd = sigma_m), nrow = n, ncol = 1)
    py = 1 / (1 + exp(-W_y %*% beta_y - beta_xy * X - beta_my * M))
    Y = matrix(as.numeric(runif(n) < py), nrow = n, ncol = 1)

  }

  return(list(X = X, M = M, Y = Y, W_x = W_x, W_m = W_m, W_y = W_y,
              X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri))

}

#######################GENERATE DATA VIA CMM####################
GEN_CMM_data = function(n, theta_mx, theta_yx, theta_ym, W,
                        beta_x, beta_m, beta_y, sigma_x, sigma_m, sigma_y,
                        X_distri, M_distri, Y_distri)
{

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)

  Z_x = rnorm(n)
  Z_m = theta_mx * Z_x + rnorm(n)
  Z_y = theta_yx * Z_x + theta_ym * Z_m + rnorm(n)

  Z_m_star = Z_m / std_Z_m
  Z_y_star = Z_y / std_Z_y

  W_x = W[[1]]
  W_m = W[[2]]
  W_y = W[[3]]

  eta_x = W_x %*% beta_x
  eta_m = W_m %*% beta_m
  eta_y = W_y %*% beta_y

  #generate X
  if (X_distri == "Gaussian") { X = eta_x + sigma_x * Z_x }


  #generate M
  if (M_distri == "Bernoulli")
  {
    prob = exp(eta_m) / (1 + exp(eta_m))
    M = as.numeric(pnorm(Z_m / std_Z_m) > 1 - prob)
  }
  if (M_distri == "Gaussian") {
    M = eta_m + sigma_m * Z_m_star
  }

  #generate Y
  if (Y_distri == "Gaussian") {
    Y = eta_y + sigma_y * Z_y_star
  }
  if (Y_distri == "Bernoulli")
  {
    prob = exp(eta_y) / (1 + exp(eta_y))
    Y = as.numeric(pnorm(Z_y_star) > 1 - prob)
  }

  return(list(X = X, M = M, Y = Y, W_x = W_x, W_m = W_m, W_y = W_y,
              X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri))
}


########################SUMMARIZE PACKAGE RESULTS###############################
pkg_summary = function(dat, x0, x1)
{
  if (dat$X_distri == "Gaussian" &
    dat$M_distri == "Gaussian" &
    dat$Y_distri == "Gaussian")
  {
    med.fit <- lm(dat$M ~ dat$X + dat$W_m + 0)
    out.fit <- lm(dat$Y ~ dat$M + dat$X + dat$W_m + 0)
    med.out <- mediate(med.fit, out.fit, treat = "dat$X", mediator = "dat$M",
                       robustSE = TRUE, sims = 500, control.value = x0, treat.value = x1)
    res = summary(med.out)
  }

  if (dat$X_distri == "Gaussian" &
    dat$M_distri == "Bernoulli" &
    dat$Y_distri == "Gaussian")
  {
    med.fit <- glm(dat$M ~ dat$X + dat$W_m + 0, family = "binomial")
    out.fit <- lm(dat$Y ~ dat$M + dat$X + dat$W_m + 0)
    med.out <- mediate(med.fit, out.fit, treat = "dat$X", mediator = "dat$M",
                       robustSE = TRUE, sims = 500, control.value = x0, treat.value = x1)
    res = summary(med.out)
  }

  if (dat$X_distri == "Gaussian" &
    dat$M_distri == "Gaussian" &
    dat$Y_distri == "Bernoulli")
  {
    med.fit <- lm(dat$M ~ dat$X + dat$W_m + 0)
    out.fit <- glm(dat$Y ~ dat$M + dat$X + dat$W_m + 0, family = "binomial")
    med.out <- mediate(med.fit, out.fit, treat = "dat$X", mediator = "dat$M",
                       robustSE = TRUE, sims = 500, control.value = x0, treat.value = x1)
    res = summary(med.out)
  }

  return(res)
}


#############find true value for SEM #############
true_val_SEM = function(n = 100000, beta_xm, beta_xy, beta_my, W, beta_x, beta_m, beta_y,
                        sigma_x, sigma_m, sigma_y, X_distri, M_distri, Y_distri, x0, x1)
{
  dat = GEN_SEM_data(n = n, beta_xm = beta_xm, beta_xy = beta_xy, beta_my = beta_my, W = W,
                     beta_x = beta_x,
                     beta_m = beta_m,
                     beta_y = beta_y,
                     sigma_x = sigma_x,
                     sigma_m = sigma_m,
                     sigma_y = sigma_y,
                     X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri)
  pkg_res = pkg_summary(dat = dat, x0 = x0, x1 = x1)

  return(list(NDE = pkg_res$z0, NIE = pkg_res$d1))
}

##############################################################################################################
############BEGIN: COVERAGE OF EFFECTS/THETAS BASED ON NONPARAMETRIC BOOTSTRP#################################


########################Type 1/Power given data are generated by SEM############
Power_SEM_comparison = function(seed, n_rep, n_list, beta_xm, beta_xy, beta_my,
                                beta_x, beta_m, beta_y,
                                sigma_x, sigma_m, sigma_y, rho_w, sigma_w,
                                X_distri, M_distri, Y_distri, save_path,
                                boot_num, prob = c(0.025, 0.975), x0, x1)
{
  set.seed(seed)

  tab = matrix(NA, nrow = 2, ncol = length(n_list))
  coverage_all = list()
  coverage_pkg_all = list()

  for (n in n_list) {
    coverage_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    coverage_pkg_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)

    W = CMM_sim_covariates(n = n, rho = rho_w, p = length(beta_m), p1 = length(beta_x),
                           sigma = sigma_w)

    mean_W_x = apply(W[[1]], 2, mean, na.rm = FALSE)
    mean_W_m = apply(W[[2]], 2, mean, na.rm = FALSE)
    mean_W_y = apply(W[[3]], 2, mean, na.rm = FALSE)
    cat("sample size = ", n, "\n")


    pb = txtProgressBar(style = 3)
    for (rp in 1:n_rep) {
      dat = GEN_SEM_data(n = n, beta_xm = beta_xm, beta_xy = beta_xy, W = W, beta_my = beta_my,
                         beta_x = beta_x,
                         beta_m = beta_m,
                         beta_y = beta_y,
                         sigma_x = sigma_x,
                         sigma_m = sigma_m,
                         sigma_y = sigma_y,
                         X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri)

      ################CMM########################################################
      theta_beta_est = CMM_theta_beta(dat)

      CI_all = CMM_boot_all(dat = dat, theta_beta = theta_beta_est,
                            mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                            x0 = x0, x1 = x1, boot_num = boot_num,
                            prob = prob)

      coverage = list(NDE = CI_all$CI_sim_effect_NDE[1,] <= 0 & CI_all$CI_sim_effect_NDE[2,] >= 0,
                      NIE = CI_all$CI_sim_effect_NIE[1,] <= 0 & CI_all$CI_sim_effect_NIE[2,] >= 0)
      coverage_all[[paste(n)]][rp,] = unlist(coverage) #1 substracts coverage of 0 is power

      ################PACKAGE####################################################
      pkg_res = pkg_summary(dat = dat, x0 = x0, x1 = x1)

      coverage_pkg = c(NDE = (pkg_res$z0.ci[1] <= 0) & (pkg_res$z0.ci[2] >= 0),
                       NIE = (pkg_res$d1.ci[1] <= 0) & (pkg_res$d1.ci[2] >= 0))
      coverage_pkg_all[[paste(n)]][rp,] = coverage_pkg #1 substracts coverage of 0 is power

      setTxtProgressBar(pb, rp / n_rep)
    }
    close(pb)
    save(coverage_all, file = file.path(save_path, sprintf("SEM_all_n%d_seed%d.RData", n, seed)))
    save(coverage_pkg_all, file = file.path(save_path, sprintf("SEM_pkg_all_n%d_seed%d.RData", n, seed)))
  }
  return(list(coverage_all = coverage_all, coverage_pkg_all = coverage_pkg_all))

}


########################Type 1/Power given data are generated by CMM############
Power_CMM_comparison = function(seed, n_rep, n_list,
                                theta_mx, theta_yx, theta_ym, beta_x, beta_m, beta_y,
                                sigma_x, sigma_m, sigma_y, rho_w, sigma_w,
                                X_distri, M_distri, Y_distri, save_path,
                                boot_num, prob = c(0.025, 0.975), x0, x1)
{
  set.seed(seed)

  tab = matrix(NA, nrow = 2, ncol = length(n_list))
  coverage_all = list()
  coverage_pkg_all = list()

  for (n in n_list) {
    coverage_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    coverage_pkg_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)

    W = CMM_sim_covariates(n = n, rho = rho_w, p = length(beta_m), p1 = length(beta_x),
                           sigma = sigma_w)

    mean_W_x = apply(W[[1]], 2, mean, na.rm = FALSE)
    mean_W_m = apply(W[[2]], 2, mean, na.rm = FALSE)
    mean_W_y = apply(W[[3]], 2, mean, na.rm = FALSE)
    cat("sample size = ", n, "\n")

    pb = txtProgressBar(style = 3)
    for (rp in 1:n_rep) {
      dat = GEN_CMM_data(n = n, theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym, W = W,
                         beta_x = beta_x, beta_m = beta_m, beta_y = beta_y,
                         sigma_x = sigma_x,
                         sigma_m = sigma_m,
                         sigma_y = sigma_y,
                         X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri)

      ################CMM########################################################
      theta_beta_est = CMM_theta_beta(dat)

      CI_all = CMM_boot_all(dat = dat, theta_beta = theta_beta_est,
                            mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                            x0 = x0, x1 = x1, boot_num = boot_num,
                            prob = prob)

      coverage = list(NDE = CI_all$CI_sim_effect_NDE[1,] <= 0 & CI_all$CI_sim_effect_NDE[2,] >= 0,
                      NIE = CI_all$CI_sim_effect_NIE[1,] <= 0 & CI_all$CI_sim_effect_NIE[2,] >= 0)
      coverage_all[[paste(n)]][rp,] = unlist(coverage) #1 substracts coverage of 0 is power

      ################PACKAGE####################################################
      pkg_res = pkg_summary(dat = dat, x0 = x0, x1 = x1)

      coverage_pkg = c(NDE = (pkg_res$z0.ci[1] <= 0) & (pkg_res$z0.ci[2] >= 0),
                       NIE = (pkg_res$d1.ci[1] <= 0) & (pkg_res$d1.ci[2] >= 0))
      coverage_pkg_all[[paste(n)]][rp,] = unlist(coverage_pkg) #1 substracts coverage of 0 is power

      setTxtProgressBar(pb, rp / n_rep)
    }
    close(pb)
    save(coverage_all, file = file.path(save_path, sprintf("CMM_all_n%d_seed%d.RData", n, seed)))
    save(coverage_pkg_all, file = file.path(save_path, sprintf("CMM_pkg_all_n%d_seed%d.RData", n, seed)))
  }
  return(list(coverage_all = coverage_all, coverage_pkg_all = coverage_pkg_all))

}


########################################
Coverage_CMM_comparison = function(seed, n_rep, n_list,
                                   theta_mx, theta_yx, theta_ym, beta_x, beta_m, beta_y,
                                   sigma_x, sigma_m, sigma_y, rho_w, sigma_w,
                                   X_distri, M_distri, Y_distri, save_path,
                                   boot_num, prob = c(0.025, 0.975), x0, x1)
{
  set.seed(seed)

  coverage_all = list()
  coverage_pkg_all = list()
  bias_all = list()
  bias_pkg_all = list()
  MSE_all = list()
  MSE_pkg_all = list()


  for (n in n_list) {
    coverage_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    coverage_pkg_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    bias_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    bias_pkg_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    MSE_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    MSE_pkg_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)

    W = CMM_sim_covariates(n = n, rho = rho_w, p = length(beta_m), p1 = length(beta_x),
                           sigma = sigma_w)

    mean_W_x = apply(W[[1]], 2, mean, na.rm = FALSE)
    mean_W_m = apply(W[[2]], 2, mean, na.rm = FALSE)
    mean_W_y = apply(W[[3]], 2, mean, na.rm = FALSE)
    cat("sample size = ", n, "\n")

    theta_true_beta_true = list(theta = c(theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym),
                                beta_x = beta_x,
                                beta_m = beta_m,
                                beta_y = beta_y,
                                sigma_x = sigma_x,
                                sigma_m = sigma_m,
                                sigma_y = sigma_y)

    NDE_true_all = effect_cal_NDE(theta_beta = theta_true_beta_true, mean_W_x = mean_W_x,
                                  mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                                  X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                                  x0 = x0, x1 = x1, tol = 1e-30)


    NIE_true_all = effect_cal_NIE(theta_beta = theta_true_beta_true, mean_W_x = mean_W_x,
                                  mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                                  X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                                  x0 = x0, x1 = x1, tol = 1e-30)

    pb = txtProgressBar(style = 3)
    for (rp in 1:n_rep) {
      dat = GEN_CMM_data(n = n, theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym, W = W,
                         beta_x = beta_x, beta_m = beta_m, beta_y = beta_y,
                         sigma_x = sigma_x,
                         sigma_m = sigma_m,
                         sigma_y = sigma_y,
                         X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri)

      ################CMM########################################################
      theta_beta_est = CMM_theta_beta(dat)

      NDE_est = effect_cal_NDE(theta_beta = theta_beta_est, mean_W_x = mean_W_x,
                               mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                               X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                               x0 = x0, x1 = x1, tol = 1e-30)


      NIE_est = effect_cal_NIE(theta_beta = theta_beta_est, mean_W_x = mean_W_x,
                               mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                               X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                               x0 = x0, x1 = x1, tol = 1e-30)

      CI_all = CMM_boot_all(dat = dat, theta_beta = theta_beta_est,
                            mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                            x0 = x0, x1 = x1, boot_num = boot_num,
                            prob = prob)

      coverage = list(NDE = CI_all$CI_sim_effect_NDE[1,] <= NDE_true_all & CI_all$CI_sim_effect_NDE[2,] >= NDE_true_all,
                      NIE = CI_all$CI_sim_effect_NIE[1,] <= NIE_true_all & CI_all$CI_sim_effect_NIE[2,] >= NIE_true_all)
      coverage_all[[paste(n)]][rp,] = unlist(coverage) #1 substracts coverage of 0 is power

      # get bias
      bias_all[[paste(n)]][rp,] = c(NDE_est, NIE_est) - c(NDE_true_all, NIE_true_all)
      # get MSE
      MSE_all[[paste(n)]][rp,] = (c(NDE_est, NIE_est) - c(NDE_true_all, NIE_true_all))^2

      ################PACKAGE####################################################
      pkg_res = pkg_summary(dat = dat, x0 = x0, x1 = x1)

      coverage_pkg = c(NDE = (pkg_res$z0.ci[1] <= NDE_true_all) & (pkg_res$z0.ci[2] >= NDE_true_all),
                       NIE = (pkg_res$d1.ci[1] <= NIE_true_all) & (pkg_res$d1.ci[2] >= NIE_true_all))
      coverage_pkg_all[[paste(n)]][rp,] = unlist(coverage_pkg) #1 substracts coverage of 0 is power

      # get bias
      bias_pkg_all[[paste(n)]][rp,] = c(pkg_res$z0, pkg_res$d1) - c(NDE_true_all, NIE_true_all)
      # get MSE
      MSE_pkg_all[[paste(n)]][rp,] = (c(pkg_res$z0, pkg_res$d1) - c(NDE_true_all, NIE_true_all))^2


      setTxtProgressBar(pb, rp / n_rep)
    }
    close(pb)

    accuracy_all = list(coverage = coverage_all, bias = bias_all, MSE = MSE_all)
    accuracy_pkg_all = list(coverage = coverage_pkg_all, bias = bias_pkg_all, MSE = MSE_pkg_all)


    save(accuracy_all, file = file.path(save_path, sprintf("accuracy_n%d_seed%d.RData", n, seed)))
    save(accuracy_pkg_all, file = file.path(save_path, sprintf("accuracy_pkg_n%d_seed%d.RData", n, seed)))

  }
  return(list(accuracy_all = accuracy_all, accuracy_pkg_all = accuracy_pkg_all))
}


####################################################################
Coverage_SEM_comparison = function(seed, n_rep, n_list, beta_xm, beta_xy, beta_my,
                                   beta_x, beta_m, beta_y,
                                   sigma_x, sigma_m, sigma_y, rho_w, sigma_w,
                                   X_distri, M_distri, Y_distri, save_path,
                                   boot_num, prob = c(0.025, 0.975), x0, x1, NDE_true_all,
                                   NIE_true_all)
{
  set.seed(seed)

  coverage_all = list()
  coverage_pkg_all = list()
  bias_all = list()
  bias_pkg_all = list()
  MSE_all = list()
  MSE_pkg_all = list()

  for (n in n_list) {
    coverage_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    coverage_pkg_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    bias_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    bias_pkg_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    MSE_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)
    MSE_pkg_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 2)

    W = CMM_sim_covariates(n = n, rho = rho_w, p = length(beta_m), p1 = length(beta_x),
                           sigma = sigma_w)

    mean_W_x = apply(W[[1]], 2, mean, na.rm = FALSE)
    mean_W_m = apply(W[[2]], 2, mean, na.rm = FALSE)
    mean_W_y = apply(W[[3]], 2, mean, na.rm = FALSE)
    cat("sample size = ", n, "\n")

    pb = txtProgressBar(style = 3)
    for (rp in 1:n_rep) {
      dat = GEN_SEM_data(n = n, beta_xm = beta_xm, beta_xy = beta_xy, W = W, beta_my = beta_my,
                         beta_x = beta_x,
                         beta_m = beta_m,
                         beta_y = beta_y,
                         sigma_x = sigma_x,
                         sigma_m = sigma_m,
                         sigma_y = sigma_y,
                         X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri)


      ################CMM########################################################
      theta_beta_est = CMM_theta_beta(dat)

      NDE_est = effect_cal_NDE(theta_beta = theta_beta_est, mean_W_x = mean_W_x,
                               mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                               X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                               x0 = x0, x1 = x1, tol = 1e-30)


      NIE_est = effect_cal_NIE(theta_beta = theta_beta_est, mean_W_x = mean_W_x,
                               mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                               X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                               x0 = x0, x1 = x1, tol = 1e-30)

      CI_all = CMM_boot_all(dat = dat, theta_beta = theta_beta_est,
                            mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                            x0 = x0, x1 = x1, boot_num = boot_num,
                            prob = prob)

      coverage = list(NDE = CI_all$CI_sim_effect_NDE[1,] <= NDE_true_all & CI_all$CI_sim_effect_NDE[2,] >= NDE_true_all,
                      NIE = CI_all$CI_sim_effect_NIE[1,] <= NIE_true_all & CI_all$CI_sim_effect_NIE[2,] >= NIE_true_all)
      coverage_all[[paste(n)]][rp,] = unlist(coverage) #1 substracts coverage of 0 is power

      # get bias
      bias_all[[paste(n)]][rp,] = c(NDE_est, NIE_est) - c(NDE_true_all, NIE_true_all)
      # get MSE
      MSE_all[[paste(n)]][rp,] = (c(NDE_est, NIE_est) - c(NDE_true_all, NIE_true_all))^2

      ################PACKAGE####################################################
      pkg_res = pkg_summary(dat = dat, x0 = x0, x1 = x1)

      coverage_pkg = c(NDE = (pkg_res$z0.ci[1] <= NDE_true_all) & (pkg_res$z0.ci[2] >= NDE_true_all),
                       NIE = (pkg_res$d1.ci[1] <= NIE_true_all) & (pkg_res$d1.ci[2] >= NIE_true_all))
      coverage_pkg_all[[paste(n)]][rp,] = unlist(coverage_pkg) #1 substracts coverage of 0 is power

      # get bias
      bias_pkg_all[[paste(n)]][rp,] = c(pkg_res$z0, pkg_res$d1) - c(NDE_true_all, NIE_true_all)
      # get MSE
      MSE_pkg_all[[paste(n)]][rp,] = (c(pkg_res$z0, pkg_res$d1) - c(NDE_true_all, NIE_true_all))^2


      setTxtProgressBar(pb, rp / n_rep)
    }
    close(pb)

    accuracy_all = list(coverage = coverage_all, bias = bias_all, MSE = MSE_all)
    accuracy_pkg_all = list(coverage = coverage_pkg_all, bias = bias_pkg_all, MSE = MSE_pkg_all)

    save(accuracy_all, file = file.path(save_path, sprintf("accuracy_n%d_seed%d.RData", n, seed)))
    save(accuracy_pkg_all, file = file.path(save_path, sprintf("accuracy_pkg_n%d_seed%d.RData", n, seed)))

  }
  return(list(accuracy_all = accuracy_all, accuracy_pkg_all = accuracy_pkg_all))
}


####################Odds ratio comparison##############################
OR_cal_NDE = function(theta_beta,
                      X_distri, M_distri, Y_distri,
                      mean_W_x, mean_W_m, mean_W_y,
                      x0, x1, tol = 1e-30)
{
  theta_mx = theta_beta$theta[1]
  theta_yx = theta_beta$theta[2]
  theta_ym = theta_beta$theta[3]
  beta_x = theta_beta$beta_x
  beta_m = theta_beta$beta_m
  beta_y = theta_beta$beta_y
  sigma_x = theta_beta$sigma_x
  sigma_m = theta_beta$sigma_m
  sigma_y = theta_beta$sigma_y

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)


  #3 X: CONTINUOUS (NORMAL), M: CONTINUOUS (NORMAL), Y: BERNOULLI

  if (X_distri == "Gaussian" &
    M_distri == "Gaussian" &
    Y_distri == "Bernoulli")
  {
    Z_x0 = (x0 - sum(beta_x * mean_W_x)) / sigma_x
    Z_x1 = (x1 - sum(beta_x * mean_W_x)) / sigma_x

    p = 1 / (1 + exp(-sum(beta_y * mean_W_y)))
    qy = qnorm(1 - p)

    Z_m_star = rnorm(1000000, mean = theta_mx * Z_x0 / std_Z_m, sd = 1 / std_Z_m)
    Phi_z_x1 = mean(pnorm(std_Z_y * qy -
                            theta_yx * Z_x1 -
                            std_Z_y * theta_ym * Z_m_star))
    Phi_z_x0 = mean(pnorm(std_Z_y * qy -
                            theta_yx * Z_x0 -
                            std_Z_y * theta_ym * Z_m_star))
    OR_NDE = Phi_z_x0 * (1 - Phi_z_x1) / (Phi_z_x1 * (1 - Phi_z_x0))
  }

  names(OR_NDE) = NULL
  return(OR_NDE)

}

OR_cal_NIE = function(theta_beta, X_distri, M_distri, Y_distri, mean_W_x, mean_W_m, mean_W_y, x0, x1, tol = 1e-30)
{
  theta_mx = theta_beta$theta[1]
  theta_yx = theta_beta$theta[2]
  theta_ym = theta_beta$theta[3]
  beta_x = theta_beta$beta_x
  beta_m = theta_beta$beta_m
  beta_y = theta_beta$beta_y
  sigma_x = theta_beta$sigma_x
  sigma_m = theta_beta$sigma_m
  sigma_y = theta_beta$sigma_y

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)


  #2 X: CONTINUOUS (NORMAL), M: CONTINUOUS (NORMAL), Y: BERNOULLI
  if (X_distri == "Gaussian" &
    M_distri == "Gaussian" &
    Y_distri == "Bernoulli")
  {
    Z_x0 = (x0 - sum(beta_x * mean_W_x)) / sigma_x
    Z_x1 = (x1 - sum(beta_x * mean_W_x)) / sigma_x

    p = 1 / (1 + exp(-sum(beta_y * mean_W_y)))
    qy = qnorm(1 - p)

    Z_m_star_0 = rnorm(1000000, mean = theta_mx * Z_x0 / std_Z_m, sd = 1 / std_Z_m)
    Z_m_star_1 = rnorm(1000000, mean = theta_mx * Z_x1 / std_Z_m, sd = 1 / std_Z_m)

    Phi_0 = mean(pnorm(std_Z_y * qy -
                         theta_yx * Z_x1 -
                         std_Z_y * theta_ym * Z_m_star_0))
    Phi_1 = mean(pnorm(std_Z_y * qy -
                         theta_yx * Z_x1 -
                         std_Z_y * theta_ym * Z_m_star_1))
    OR_NIE = Phi_0 * (1 - Phi_1) / (Phi_1 * (1 - Phi_0))
  }

  names(OR_NIE) = NULL
  return(OR_NIE)

}

OR_VW = function(dat, x0, x1)
{

  med.fit <- lm(dat$M ~ dat$X + dat$W_m + 0)
  out.fit <- glm(dat$Y ~ dat$M + dat$X + dat$W_m + 0, family = "binomial")
  OR_NDE_VW = exp(coef(out.fit)[2] * (x1 - x0))
  OR_NIE_VW = exp(coef(out.fit)[1] * coef(med.fit)[1] * (x1 - x0))

  return(c(OR_NDE_VW = OR_NDE_VW, OR_NIE_VW = OR_NIE_VW))
}


OR_CMM_comparison = function(seed, n_rep, n_list,
                             theta_mx, theta_yx, theta_ym, beta_x, beta_m, beta_y,
                             sigma_x, sigma_m, sigma_y, rho_w, sigma_w,
                             X_distri, M_distri, Y_distri, save_path,
                             x0, x1)
{
  set.seed(seed)

  OR = list()

  for (n in n_list) {
    OR[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 6)

    W = CMM_sim_covariates(n = n, rho = rho_w, p = length(beta_m), p1 = length(beta_x),
                           sigma = sigma_w)

    mean_W_x = apply(W[[1]], 2, mean, na.rm = FALSE)
    mean_W_m = apply(W[[2]], 2, mean, na.rm = FALSE)
    mean_W_y = apply(W[[3]], 2, mean, na.rm = FALSE)
    cat("sample size = ", n, "\n")

    pb = txtProgressBar(style = 3)
    for (rp in 1:n_rep) {
      dat = GEN_CMM_data(n = n, theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym, W = W,
                         beta_x = beta_x, beta_m = beta_m, beta_y = beta_y,
                         sigma_x = sigma_x,
                         sigma_m = sigma_m,
                         sigma_y = sigma_y,
                         X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri)

      ################CMM########################################################
      theta_true_beta_true = list(theta = c(theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym),
                                  beta_x = beta_x,
                                  beta_m = beta_m,
                                  beta_y = beta_y,
                                  sigma_x = sigma_x,
                                  sigma_m = sigma_m,
                                  sigma_y = sigma_y)

      OR_NDE_true = OR_cal_NDE(theta_beta = theta_true_beta_true,
                               X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                               mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                               x0 = x0, x1 = x1, tol = 1e-30)


      OR_NIE_true = OR_cal_NIE(theta_beta = theta_true_beta_true,
                               X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                               mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                               x0 = x0, x1 = x1, tol = 1e-30)

      theta_beta = CMM_theta_beta(dat)


      OR_NDE_CMM = OR_cal_NDE(theta_beta = theta_beta,
                              X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                              mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                              x0 = x0, x1 = x1, tol = 1e-30)
      OR_NIE_CMM = OR_cal_NIE(theta_beta = theta_beta,
                              X_distri = X_distri, M_distri = M_distri, Y_distri = Y_distri,
                              mean_W_x = mean_W_x, mean_W_m = mean_W_m, mean_W_y = mean_W_y,
                              x0 = x0, x1 = x1, tol = 1e-30)

      temp = OR_VW(dat = dat, x0 = x0, x1 = x1)
      OR_NDE_VW = temp[1]
      OR_NIE_VW = temp[2]

      OR[[paste(n)]][rp,] = c(OR_NDE_true, OR_NDE_CMM, OR_NDE_VW, OR_NIE_true, OR_NIE_CMM, OR_NIE_VW)

    }
    close(pb)
    save(OR, file = file.path(save_path, sprintf("OR_n%d_seed%d.RData", n, seed)))
  }
  return(OR)

}


#####FULL MLE NEG LIKELI FUN FOR CCC####################

MLE_neg_log_like_fun = function(theta, X, M, Y, X_distri, M_distri, Y_distri)
{
  theta_mx = theta[1]
  theta_yx = theta[2]
  theta_ym = theta[3]

  alpha = theta_yx + theta_mx * theta_ym
  std_Z_m = sqrt(theta_mx^2 + 1)
  std_Z_y = sqrt(alpha^2 + theta_ym^2 + 1)


  #1 X: CONTINUOUS(NORMAL), M: CONTINUOUS(NORMAL), Y: CONTINUOUS(NORMAL)
  if (X_distri == "Gaussian" &
    M_distri == "Gaussian" &
    Y_distri == "Gaussian")
  {

    mu_x = theta[4]
    mu_m = theta[5]
    mu_y = theta[6]
    sigma_x = exp(theta[7])
    sigma_m = exp(theta[8])
    sigma_y = exp(theta[9])

    Z_x = (X - mu_x) / sigma_x
    Z_m_star = (M - mu_m) / sigma_m
    Z_y_star = (Y - mu_y) / sigma_y
    D = diag(c(sigma_x, sigma_m, sigma_y))

    Sigma = D %*%
      matrix(c(1, theta_mx / std_Z_m, alpha / std_Z_y,
               theta_mx / std_Z_m, 1, (theta_mx * alpha + theta_ym) / (std_Z_m * std_Z_y),
               alpha / std_Z_y, (theta_mx * alpha + theta_ym) / (std_Z_m * std_Z_y), 1), nrow = 3, ncol = 3, byrow = T) %*%
      D

    pi_Z_x_Z_m_star_Z_y_star = dmvnorm(cbind(X, M, Y), mean = c(mu_x, mu_m, mu_y),
                                       sigma = Sigma, log = T)
    MLE_neg_log_like = -sum(pi_Z_x_Z_m_star_Z_y_star)
  }

  #2 X: CONTINUOUS(NORMAL), M: DISCRETE(BERNOULLI), Y: CONTINUOUS(NORMAL)
  if (X_distri == "Gaussian" &
    M_distri == "Bernoulli" &
    Y_distri == "Gaussian")
  {
    mu_x = theta[4]
    odds_M = theta[5]
    mu_y = theta[6]
    sigma_x = exp(theta[7])
    sigma_y = exp(theta[8])
    Z_x = (X - mu_x) / sigma_x
    p = 1 / (1 + exp(-odds_M))
    lm = qnorm(ifelse(M == 0, 0, 1 - p))
    um = qnorm(ifelse(M == 0, 1 - p, 1))
    Z_y_star = (Y - mu_y) / sigma_y


    pi_Z_x_Z_y_star = sapply(1:length(lm), function(i)
      -log(sigma_x) - log(sigma_y) + dmvnorm(c(Z_x[i], Z_y_star[i]), mean = c(0, 0),
                                             sigma = matrix(c(1, alpha / std_Z_y, alpha / std_Z_y, 1), nrow = 2, ncol = 2, byrow = T), log = T))

    mean_Z_m_star = ((theta_mx - theta_ym * theta_yx) * Z_x + theta_ym * std_Z_y * Z_y_star) / (std_Z_m * (theta_ym^2 + 1))
    std_Z_m_star = 1 / sqrt((theta_ym^2 + 1) * (theta_mx^2 + 1))

    prob_Z_m_star_Z_x_Z_y_star = sapply(1:length(lm), function(i) {
      pnorm(um[i], mean = mean_Z_m_star[i], sd = std_Z_m_star) - pnorm(lm[i], mean = mean_Z_m_star[i], sd = std_Z_m_star)
    })

    neg_log_like_vec = pi_Z_x_Z_y_star + log(prob_Z_m_star_Z_x_Z_y_star)
    MLE_neg_log_like = -sum(neg_log_like_vec)

  }

  return(MLE_neg_log_like)
}


MLE_theta = function(dat, para0 = c(0, 0, 0, mean(dat$X),
                                    mean(dat$M), mean(dat$Y),
                                    log(sd(dat$X)),
                                    log(sd(dat$Y))), trace = 0, method = "BFGS") {

  res_opt = optim(par = para0, MLE_neg_log_like_fun, X = dat$X, M = dat$M, Y = dat$Y,
                  X_distri = dat$X_distri, M_distri = dat$M_distri, Y_distri = dat$Y_distri,
                  method = method, control = list(trace = trace))
  theta = res_opt$par
  theta[7:9] = exp(theta[7:9])
  loglik = -res_opt$value
  return(list(theta = theta, loglik = loglik))
}


MLE_theta_CDC = function(dat, para0 = c(0, 0, 0, mean(dat$X),
                                        log(mean(dat$M) / (1 - mean(dat$M))), mean(dat$Y),
                                        log(sd(dat$X)),
                                        log(sd(dat$Y))), trace = 0, method = "BFGS") {

  res_opt = optim(par = para0, MLE_neg_log_like_fun, X = dat$X, M = dat$M, Y = dat$Y,
                  X_distri = dat$X_distri, M_distri = dat$M_distri, Y_distri = dat$Y_distri,
                  method = method, control = list(trace = trace))
  theta = res_opt$par
  theta[7:8] = exp(theta[7:8])
  loglik = -res_opt$value
  return(list(theta = theta, loglik = loglik))
}


MLE_CMM_comparison = function(seed, n_rep, n_list,
                              theta_mx, theta_yx, theta_ym, beta_x, beta_m, beta_y,
                              sigma_x, sigma_m, sigma_y, rho_w, sigma_w,
                              X_distri, M_distri, Y_distri, save_path)
{
  set.seed(seed)

  est_all = list()
  est_mle_all = list()
  bias_all = list()
  bias_mle_all = list()
  MSE_all = list()
  MSE_mle_all = list()

  for (n in n_list) {
    est_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    est_mle_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    bias_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    bias_mle_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    MSE_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    MSE_mle_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)

    W = CMM_sim_covariates(n = n, rho = rho_w, p = 4, p1 = 3,
                           sigma = sigma_w)

    cat("sample size = ", n, "\n")

    theta_true_beta_true = list(theta = c(theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym),
                                beta_x = beta_x,
                                beta_m = beta_m,
                                beta_y = beta_y,
                                sigma_x = sigma_x,
                                sigma_m = sigma_m,
                                sigma_y = sigma_y)


    pb = txtProgressBar(style = 3)
    for (rp in 1:n_rep) {
      dat = CMM_sim_dat(n = n, W_x = W[[1]][, 1], W_m = W[[2]][, 1], W_y = W[[3]][, 1],
                        theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym,
                        beta_x = matrix(beta_x),
                        beta_m = matrix(beta_m),
                        beta_y = matrix(beta_y),
                        sigma_x = sigma_x,
                        sigma_m = sigma_m,
                        sigma_y = sigma_y,
                        X_distri = X_distri,
                        M_distri = M_distri,
                        Y_distri = Y_distri)

      ################CMM########################################################
      theta_beta_est = CMM_theta_beta(dat)

      est_all[[paste(n)]][rp,] = theta_beta_est$theta

      # get bias
      bias_all[[paste(n)]][rp,] = theta_beta_est$theta - theta_true_beta_true$theta
      # get MSE
      MSE_all[[paste(n)]][rp,] = (theta_beta_est$theta - theta_true_beta_true$theta)^2

      ################MLE####################################################
      theta_mle = MLE_theta(dat = dat, trace = 0)

      est_mle_all[[paste(n)]][rp,] = theta_mle$theta[1:3]

      # get bias
      bias_mle_all[[paste(n)]][rp,] = theta_mle$theta[1:3] - theta_true_beta_true$theta
      # get MSE
      MSE_mle_all[[paste(n)]][rp,] = (theta_mle$theta[1:3] - theta_true_beta_true$theta)^2


      setTxtProgressBar(pb, rp / n_rep)
    }
    close(pb)

    accuracy_all = list(est = est_all, bias = bias_all, MSE = MSE_all)
    accuracy_mle_all = list(est = est_mle_all, bias = bias_mle_all, MSE = MSE_mle_all)


    save(accuracy_all, file = file.path(save_path, sprintf("accuracy_n%d_seed%d.RData", n, seed)))
    save(accuracy_mle_all, file = file.path(save_path, sprintf("accuracy_mle_n%d_seed%d.RData", n, seed)))

  }
  return(list(accuracy_all = accuracy_all, accuracy_mle_all = accuracy_mle_all))
}


MLE_CMM_CDC_comparison = function(seed, n_rep, n_list,
                                  theta_mx, theta_yx, theta_ym, beta_x, beta_m, beta_y,
                                  sigma_x, sigma_m, sigma_y, rho_w, sigma_w,
                                  X_distri, M_distri, Y_distri, save_path)
{
  set.seed(seed)

  est_all = list()
  est_mle_all = list()
  bias_all = list()
  bias_mle_all = list()
  MSE_all = list()
  MSE_mle_all = list()

  for (n in n_list) {
    est_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    est_mle_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    bias_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    bias_mle_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    MSE_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)
    MSE_mle_all[[paste(n)]] = matrix(NA, nrow = n_rep, ncol = 3)

    W = CMM_sim_covariates(n = n, rho = rho_w, p = 4, p1 = 3,
                           sigma = sigma_w)

    cat("sample size = ", n, "\n")

    theta_true_beta_true = list(theta = c(theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym),
                                beta_x = beta_x,
                                beta_m = beta_m,
                                beta_y = beta_y,
                                sigma_x = sigma_x,
                                sigma_m = sigma_m,
                                sigma_y = sigma_y)


    pb = txtProgressBar(style = 3)
    for (rp in 1:n_rep) {
      dat = CMM_sim_dat(n = n, W_x = W[[1]][, 1], W_m = W[[2]][, 1], W_y = W[[3]][, 1],
                        theta_mx = theta_mx, theta_yx = theta_yx, theta_ym = theta_ym,
                        beta_x = matrix(beta_x),
                        beta_m = matrix(beta_m),
                        beta_y = matrix(beta_y),
                        sigma_x = sigma_x,
                        sigma_m = sigma_m,
                        sigma_y = sigma_y,
                        X_distri = X_distri,
                        M_distri = M_distri,
                        Y_distri = Y_distri)

      ################CMM########################################################
      theta_beta_est = CMM_theta_beta(dat)

      est_all[[paste(n)]][rp,] = theta_beta_est$theta

      # get bias
      bias_all[[paste(n)]][rp,] = theta_beta_est$theta - theta_true_beta_true$theta
      # get MSE
      MSE_all[[paste(n)]][rp,] = (theta_beta_est$theta - theta_true_beta_true$theta)^2

      ################MLE####################################################
      theta_mle = MLE_theta_CDC(dat = dat, trace = 0)

      est_mle_all[[paste(n)]][rp,] = theta_mle$theta[1:3]

      # get bias
      bias_mle_all[[paste(n)]][rp,] = theta_mle$theta[1:3] - theta_true_beta_true$theta
      # get MSE
      MSE_mle_all[[paste(n)]][rp,] = (theta_mle$theta[1:3] - theta_true_beta_true$theta)^2


      setTxtProgressBar(pb, rp / n_rep)
    }
    close(pb)

    accuracy_all = list(est = est_all, bias = bias_all, MSE = MSE_all)
    accuracy_mle_all = list(est = est_mle_all, bias = bias_mle_all, MSE = MSE_mle_all)


    save(accuracy_all, file = file.path(save_path, sprintf("accuracy_n%d_seed%d.RData", n, seed)))
    save(accuracy_mle_all, file = file.path(save_path, sprintf("accuracy_mle_n%d_seed%d.RData", n, seed)))

  }
  return(list(accuracy_all = accuracy_all, accuracy_mle_all = accuracy_mle_all))
}








