# packages
library(adaptr)
library(ggplot2)
library(reshape)

# simulate theta (freq or bayesian)
get_theta <- function(bayesian, theta_mean, theta_sd){
  if (bayesian==TRUE){
    bayesian_theta_sim <- rnorm(1, mean=theta_mean, sd=theta_sd)
    return(bayesian_theta_sim)
  }
  else{
    freq_theta_sim <- theta_mean
    return(freq_theta_sim)
  }
}

# calculate gain of fixed sampling trial
gain_fixed <- function(mu_x, mu_y, sigma, a, m, r_x, n, theta_sd=3){
  theta <- rnorm(1, mean=mu_x-mu_y, sd=theta_sd)
  theta_hat_sd <- sqrt((sigma^2)/(n*r_x*(1-r_x)))
  total_gain = 0
  for (i in 1:m){
    set.seed(i)
    X = rnorm(r_x*n, mean=mu_x, sd = sigma)
    Y = rnorm((1-r_x)*n, mean=mu_x-theta, sd = sigma)
    X_bar = mean(X)
    Y_bar = mean(Y)
    theta_hat <- X_bar - Y_bar
    z_stat <- theta_hat/theta_hat_sd
    p_value <- pnorm(z_stat, lower.tail = FALSE)
    comp_1 <- a*ifelse(p_value<0.05, 1, 0)
    comp_2 <- (1-a)*ifelse(mu_x > mu_y, r_x, 1-r_x)
    trial_gain <- comp_1 + comp_2
    total_gain = total_gain + trial_gain
  }
  mean_gain <- total_gain/m
  return(mean_gain)
}

# calculate optimal gain for fixed sampling trial
gain_fixed_opt <- function(mu_x, mu_y, sigma, a, m, n, cap, ret=F, theta_sd=3){
  opt_gain <- 0
  opt_r_x <- 0
  if (cap==T){
    r_xs <- seq(0.4, 0.6, 0.1)
  }
  else{
    r_xs <- seq(0.1, 0.9, 0.1)
  }
  for (r_x in r_xs){
    theta_hat_sd <- sqrt((sigma^2)/(n*r_x*(1-r_x)))
    total_gain <- 0
    for (i in 1:m){
      set.seed(i)
      theta <- rnorm(1, mean=mu_x-mu_y, sd=theta_sd)
      theta_hat_sd <- sqrt((sigma^2)/(n*r_x*(1-r_x)))
      theta_hat <- rnorm(1, mean=theta, sd=theta_hat_sd)
      z_stat <- theta_hat/theta_hat_sd
      p_value <- pnorm(z_stat, lower.tail = FALSE)
      comp_1 <- a*ifelse(p_value<0.05, 1, 0)
      comp_2 <- (1-a)*ifelse(mu_x > mu_y, r_x, 1-r_x)
      trial_gain <- comp_1 + comp_2
      total_gain <- total_gain + trial_gain
    }
    if (total_gain > opt_gain){
      opt_gain <- total_gain
      opt_r_x <- r_x
    }
  }
  mean_gain <- opt_gain/m
  if (ret==T){
    return(c(mean_gain, opt_r_x))
  }
  return(mean_gain)
}

#calculate gain of RAR trial
gain_adaptive <- function(mu_x, mu_y, sigma, a, m, n, theta_sd=3){
  total_gain <- 0
  matrix <- matrix(NA, nrow = m, ncol = 5)
  for (i in 1:m){
    set.seed(i)
    theta <- rnorm(1, mean=mu_x-mu_y, sd=theta_sd)
    trial_setup <- setup_trial_norm(arms=c("Treatment X", "Control Y"), true_ys=c(mu_x, mu_x-theta), sds=c(sigma, sigma), 
                                    data_looks=c(0.5*n, n), highest_is_best=TRUE, control="Control Y",
                                    soften_power=1, inferiority=c(0, 0.05), superiority=c(1, 0.95))
    trial_sim <- run_trial(trial_spec=trial_setup, seed=i, sparse=T)
    prob_x <- trial_sim$trial_res$ns[1]/n
    status_x <- trial_sim$trial_res$final_status[1]
    comp_1 <- a*ifelse(mu_x > mu_y & status_x=="superior" | mu_y >= mu_x & status_x=="inferior", 1, 0)
    comp_2 <- (1-a)*ifelse(mu_x > mu_y, prob_x, 1-prob_x)
    trial_gain <- comp_1 + comp_2
    total_gain <- total_gain + trial_gain
  }
  expected_gain <- total_gain/m
  return(expected_gain)
}

# function to compare trial designs
loop_compare <- function(mu_x, mu_y, sigma, m, n, r_x, gran, cap){
  weights <- seq(from = 0, to = 1, by = gran)
  fixed_gains <- rep(0, length(weights))
  fixed_gains_opt <- rep(0, length(weights))
  adaptive_gains <- rep(0, length(weights))
  
  for (i in 1:length(weights)){
    fixed_gains[i] <- gain_fixed(mu_x, mu_y, sigma, weights[i], m, r_x, n)
    fixed_gains_opt[i] <- gain_fixed_opt(mu_x, mu_y, sigma, weights[i], m, n, cap)
    adaptive_gains[i] <- gain_adaptive(mu_x, mu_y, sigma, weights[i], m, n)
  }
  results <- data.frame(weights, fixed_gains, fixed_gains_opt, adaptive_gains) #fixed_gains_opt_new
  mresults <- melt(results, id=c("weights"))
  ylabel <- paste0("Expected gain over ", m, " simulations")
  plot <- ggplot(results) + geom_line(aes(x=weights, y=fixed_gains, col="fixed", linetype="fixed")) + 
    geom_line(aes(x=weights, y=fixed_gains_opt, col="per weight optimal fixed", linetype="per weight optimal fixed")) + 
    geom_line(aes(x=weights, y=adaptive_gains, col="adaptive", linetype="adaptive")) + xlab("weighting on statistical power") + 
    ylab(ylabel) + labs(col="Trial design", linetype="Trial design") + scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  return(plot)
}

# function to compare all fixed ratios
loop_compare_fast_plot <- function(mu_x, mu_y, sigma, m, n, r_x, gran, cap, opt, theta_sd=3){
  weights <- seq(from = 0, to = 1, by = gran)
  comp1_fixed_vec <- rep(0, m)
  comp2_fixed_vec <- rep(0, m)
  comp1_rar_vec <- rep(0, m)
  comp2_rar_vec <- rep(0, m)
  fixed_gains <- rep(0, length(weights))
  adaptive_gains <- rep(0, length(weights))
  if (opt==T){
    fixed_gains_opt <- rep(0, length(weights))
  }
    for (j in 1:m){
      set.seed(j)
      theta <- rnorm(1, mean=mu_x-mu_y, sd=theta_sd)
      theta_hat_sd <- sqrt((sigma^2)/(n*r_x*(1-r_x)))
      theta_hat <- rnorm(1, mean=theta, sd=theta_hat_sd)
      z_stat <- theta_hat/theta_hat_sd
      p_value <- pnorm(z_stat, lower.tail = FALSE)
      comp1_fixed <- ifelse(p_value<0.05, 1, 0)
      comp1_fixed_vec[j] <- comp1_fixed
      comp2_fixed <- ifelse(mu_x > mu_y, r_x, 1-r_x)
      comp2_fixed_vec[j] <- comp2_fixed
      trial_setup <- setup_trial_norm(arms=c("Treatment X", "Control Y"), true_ys=c(theta, 0), sds=c(sigma, sigma), 
                                      data_looks=c(0.5*n, n), highest_is_best=TRUE, control="Control Y",
                                      soften_power=1, inferiority=c(0, 0.05), superiority=c(1, 0.95))
      trial_sim <- run_trial(trial_spec=trial_setup, seed=j, sparse=T)
      prob_x <- trial_sim$trial_res$ns[1]/n
      status_x <- trial_sim$trial_res$final_status[1]
      comp1_rar <- ifelse(theta>0 & status_x=="superior" | theta<=0 & status_x=="inferior", 1, 0)
      comp1_rar_vec[j] <- comp1_rar
      comp2_rar <- ifelse(theta>0, prob_x, 1-prob_x)
      comp2_rar_vec[j] <- comp2_rar
    }
    for (i in 1:length(weights)){
      a1 <- weights[i]
      fixed_gains[i] <- a1*mean(comp1_fixed_vec) + (1-a1)*mean(comp2_fixed_vec)
      adaptive_gains[i] <- a1*mean(comp1_rar_vec) + (1-a1)*mean(comp2_rar_vec)
      if (opt==T){
        fixed_gains_opt[i] <- gain_fixed_opt(mu_x, mu_y, sigma, a1, m, n, cap)
      }
    }
  if (opt==T){
    results <- data.frame(weights, fixed_gains, fixed_gains_opt, adaptive_gains)
    ylabel <- paste0("Expected gain over ", m, " simulations")
    plot <- ggplot(results) + geom_line(aes(x=weights, y=fixed_gains, col="fixed", linetype="fixed")) + 
      geom_line(aes(x=weights, y=fixed_gains_opt, col="optimal fixed", linetype="optimal fixed")) + 
      geom_line(aes(x=weights, y=adaptive_gains, col="adaptive", linetype="adaptive")) + xlab("w_1") + 
      ylab(ylabel) + labs(col="Trial design", linetype="Trial design") + scale_x_continuous(breaks = seq(0, 1, by = 0.1))
    }
  else{
    results <- data.frame(weights, fixed_gains, adaptive_gains)
    ylabel <- paste0("Expected gain over ", m, " simulations")
    plot <- ggplot(results) + geom_line(aes(x=weights, y=fixed_gains, col="fixed", linetype="fixed")) + 
      geom_line(aes(x=weights, y=adaptive_gains, col="adaptive", linetype="adaptive")) + xlab("w_1") + 
      ylab(ylabel) + labs(col="Trial design", linetype="Trial design") + scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  }
  #return(results)
  return(plot)
}

loop_compare_fast_plot_best_trial <- function(theta_mean, sigma, m, n, a, cap, opt, theta_sd=3){
  comp1_fixed_vec <- rep(0, m)
  comp2_fixed_vec <- rep(0.5, m)
  comp1_rar_vec <- rep(0, m)
  comp2_rar_vec <- rep(0, m)
  fixed_gain_opt=0
  for (j in 1:m){
    set.seed(j)
    theta <- rnorm(1, mean=theta_mean, sd=theta_sd)
    theta_hat_sd <- sqrt((sigma^2)/(n*0.5*0.5))
    theta_hat <- rnorm(1, mean=theta, sd=theta_hat_sd)
    z_stat <- theta_hat/theta_hat_sd
    p_value <- pnorm(z_stat, lower.tail = FALSE)
    comp1_fixed_vec[j] <- ifelse(p_value<0.05, 1, 0)
    trial_setup <- setup_trial_norm(arms=c("Treatment X", "Control Y"), true_ys=c(theta, 0), sds=c(sigma, sigma), 
                                    data_looks=c(0.5*n, n), highest_is_best=TRUE, control="Control Y",
                                    soften_power=1, inferiority=c(0, 0.05), superiority=c(1, 0.95))
    trial_sim <- run_trial(trial_spec=trial_setup, seed=j, sparse=T)
    prob_x <- trial_sim$trial_res$ns[1]/n
    status_x <- trial_sim$trial_res$final_status[1]
    comp1_rar_vec[j] <- ifelse(theta>0 & status_x=="superior" | theta <= 0 & status_x=="inferior", 1, 0)
    comp2_rar_vec[j] <- ifelse(theta>0, prob_x, 1-prob_x)
  }
  fixed_gain <- a*mean(comp1_fixed_vec) + (1-a)*mean(comp2_fixed_vec)
  adaptive_gain <- a*mean(comp1_rar_vec) + (1-a)*mean(comp2_rar_vec)
  if (opt==T){
    fixed_gain_opt <- gain_fixed_opt(theta, 0, sigma, a, m, n, cap)
  }
  if(fixed_gain > adaptive_gain & fixed_gain > fixed_gain_opt){
    best="fixed"
  }
  else if(adaptive_gain>fixed_gain & adaptive_gain>fixed_gain_opt){
    best="adaptive"
  }
  else if(fixed_gain_opt>fixed_gain & fixed_gain_opt>adaptive_gain){
    best="optimal fixed"
  }
  else{
    best="inconclusive"
  }
  return(best)
}

get_comps <- function(thetas, sigmas, m, n, a, theta_sd=3){
  param_space <- matrix(NA, nrow=length(thetas), ncol=length(sigmas))
  comp1_fixed_vec <- rep(0, m)
  comp2_fixed_vec <- rep(0.5, m)
  comp1_rar_vec <- rep(0, m)
  comp2_rar_vec <- rep(0, m)
  theta <- rep(thetas, each=length(sigmas))
  sigma <- rep(sigmas, times=length(thetas))
  for (i in 1:length(thetas)){
    the <- theta[i]
    for (j in 1:length(sigmas)){
      sig <- sigma[j]
      print(c(i, j))
      for (k in 1:m){
        set.seed(k)
        theta <- rnorm(1, mean=the, sd=theta_sd)
        theta_hat_sd <- sqrt((sig^2)/(n*0.5*0.5))
        theta_hat <- rnorm(1, mean=theta, sd=theta_hat_sd)
        z_stat <- theta_hat/theta_hat_sd
        p_value <- pnorm(z_stat, lower.tail = FALSE)
        comp1_fixed_vec[k] <- ifelse(p_value<0.05, 1, 0)
        trial_setup <- setup_trial_norm(arms=c("Treatment X", "Control Y"), true_ys=c(the, 0), sds=c(sig, sig), 
                                        data_looks=c(0.5*n, n), highest_is_best=TRUE, control="Control Y",
                                        soften_power=1, inferiority=c(0, 0.05), superiority=c(1, 0.95))
        trial_sim <- run_trial(trial_spec=trial_setup, seed=k, sparse=T)
        prob_x <- trial_sim$trial_res$ns[1]/n
        status_x <- trial_sim$trial_res$final_status[1]
        comp1_rar_vec[k] <- ifelse(theta>0 & status_x=="superior" | theta <= 0 & status_x=="inferior", 1, 0)
        comp2_rar_vec[k] <- ifelse(theta>0, prob_x, 1-prob_x)
      }
      fixed_gain <- a*mean(comp1_fixed_vec) + (1-a)*mean(comp2_fixed_vec)
      adaptive_gain <- a*mean(comp1_rar_vec) + (1-a)*mean(comp2_rar_vec)
      max_gain <- max(fixed_gain, adaptive_gain)#, fixed_gain_opt)
      if(max_gain==fixed_gain){
        best="fixed"
      }
      else if(max_gain==adaptive_gain){
        best="adaptive"
      }
      param_space[i, j] <- best
    }
  }
  results <- c(t(param_space))
  df <- data.frame(theta, sigma, results)
  return(df)
}

make_best_plot <- function(thetas, sigmas, m, n, a, cap, opt){
  param_space <- matrix(NA, nrow=length(thetas), ncol=length(sigmas))
  mu_theta <- rep(thetas, each=length(sigmas))
  sigma_p <- rep(sigmas, times=length(thetas))
  best <- rep(0, length(mu_theta))
  for (i in 1:length(best)){
    #print(loop_compare_fast_plot_best_trial(theta[i], sigma[i], m, n, a))
    best[i] <- loop_compare_fast_plot_best_trial(mu_theta[i], sigma_p[i], m, n, a, cap, opt)
    print(i)
  }
  df <- data.frame(mu_theta, sigma_p, best)
  #return(df)
  plot <- ggplot(df) + geom_point(aes(x=mu_theta, y=sigma_p, col=best, shape = best)) + 
    labs(col="Best trial design", shape="Best trial design")
  return(plot)
}

fixed_compare <- function(thetas, sigmas, m, n, theta_sd=3){
  theta_plot <- rep(thetas, each=length(sigmas))
  sigma_plot <- rep(sigmas, times=length(thetas))
  comp1_plot <- rep(0, length(theta_plot))
  comp2_plot <- rep(0.5, length(theta_plot))
  for (i in 1:length(theta_plot)){
    theta <- theta_plot[i]
    sigma <- sigma_plot[i]
    comp1_fixed_vec <- rep(0, m)
    comp2_fixed_vec <- rep(0, m)
    for (k in 1:m){
      set.seed(k)
      theta_sim <- rnorm(1, mean=theta, sd=theta_sd)
      theta_hat_sd <- sqrt((sigma^2)/(n*0.5*0.5))
      theta_hat <- rnorm(1, mean=theta_sim, sd=theta_hat_sd)
      z_stat <- theta_hat/theta_hat_sd
      p_value <- pnorm(z_stat, lower.tail = FALSE)
      comp1_fixed_vec[k] <- ifelse(p_value<0.05, 1, 0)
    }
    comp1_plot[i] <- mean(comp1_fixed_vec)
  }
  res_df <- data.frame(theta_plot, sigma_plot, comp1_plot, comp2_plot)
  return(res_df)
}

rar_compare <- function(thetas, sigmas, m, n, theta_sd=3){
  theta_plot <- rep(thetas, each=length(sigmas))
  sigma_plot <- rep(sigmas, times=length(thetas))
  comp1_plot <- rep(0, length(theta_plot))
  comp2_plot <- rep(0, length(theta_plot))
  for (i in 1:length(theta_plot)){
    theta <- theta_plot[i]
    sigma <- sigma_plot[i]
    comp1_rar_vec <- rep(0, m)
    comp2_rar_vec <- rep(0, m)
    for (k in 1:m){
      set.seed(k)
      theta_sim <- rnorm(1, mean=theta, sd=theta_sd)
      trial_setup <- setup_trial_norm(arms=c("Treatment X", "Control Y"), true_ys=c(theta, 0), sds=c(sigma, sigma), 
                                      data_looks=c(0.5*n, n), highest_is_best=TRUE, control="Control Y",
                                      soften_power=1, inferiority=c(0, 0.05), superiority=c(1, 0.95))
      trial_sim <- run_trial(trial_spec=trial_setup, seed=k, sparse=T)
      prob_x <- trial_sim$trial_res$ns[1]/n
      status_x <- trial_sim$trial_res$final_status[1]
      comp1_rar_vec[k] <- ifelse(theta>0 & status_x=="superior" | theta <= 0 & status_x=="inferior", 1, 0)
      comp2_rar_vec[k] <- ifelse(theta>0, prob_x, 1-prob_x)
    }
    comp1_plot[i] <- mean(comp1_rar_vec)
    comp2_plot[i] <- mean(comp2_rar_vec)
  }
  res_df <- data.frame(theta_plot, sigma_plot, comp1_plot, comp2_plot)
  return(res_df)
}
