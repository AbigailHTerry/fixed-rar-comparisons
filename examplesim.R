source("functions.R")
library(ggpmisc)

mu_x <- 20
mu_y <- 10
sigma <- 28
r_x <- 0.5
m <- 500
a <- 0.5
n <- 200
theta_sd <- 3

theta <- rnorm(1, mean=mu_x-mu_y, sd=theta_sd)
trial_setup <- setup_trial_norm(arms=c("Treatment X", "Control Y"), true_ys=c(mu_x, mu_x-theta), sds=c(sigma, sigma), 
                                data_looks=c(0.5*n, n), highest_is_best=TRUE, control="Control Y",
                                soften_power=1, inferiority=c(0, 0.05), superiority=c(1, 0.95))
# print(trial_setup)
#trial_sim <- run_trial(trial_spec=trial_setup, seed=2, sparse=F)
#rar <- rar_sim(10, 10, 1, 100, 200, 3)
#plot_history(trial_sim, x_value = "total n", y_value = "n all")

trial_sim <- run_trial(trial_spec=trial_setup, seed=4, sparse=F)
alloc_x <- trial_sim$trial_res$final_alloc[1]
final_alloc <- ifelse(mu_x >= mu_y, alloc_x, 1-alloc_x)
final_allocs[i] <- final_alloc

# maximum expected gain for optimised
for (i in seq(0, 1, 0.1)){
  print(gain_fixed_opt(20, 10, 28, i, 500, 200, F, T))
}

# maximum expected gain for optimised (capped)
for (i in seq(0, 1, 0.1)){
  print(gain_fixed_opt(20, 10, 28, i, 500, 200, T, T))
}

#example 1
set.seed(18)
theta1 <- rnorm(1, mean=mu_x-mu_y, sd=theta_sd)
trial_setup1 <- setup_trial_norm(arms=c("Treatment X", "Control Y"), true_ys=c(mu_x, mu_x-theta1), sds=c(sigma, sigma), 
                                data_looks=c(0.5*n, n), highest_is_best=TRUE, control="Control Y",
                                soften_power=1, inferiority=c(0, 0.05), superiority=c(1, 0.95))
trial_sim1 <- run_trial(trial_spec=trial_setup1, seed=18, sparse=F)
print(trial_sim1$trial_res$final_alloc[1])
print(trial_sim1$trial_res$ns_all[1])
print(trial_sim1$all_looks[[1]]$sum_ys[1]-trial_sim1$all_looks[[1]]$sum_ys[2])
plot_rules <- theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
                    legend.position = "bottom", plot.margin = unit(c(t=0.2, r=1, b=0.2, l=1), "cm"))

pdf(file = "plots/rarex1alln.pdf", width = 5, height = 4)
plot_history(trial_sim1, x_value = "total n", y_value = "n all") + plot_rules
dev.off()
pdf(file = "plots/rarex1prob.pdf", width = 5, height = 4)
plot_history(trial_sim1, x_value = "total n", y_value = "prob") + plot_rules
dev.off()
pdf(file = "plots/rarex1sumys.pdf", width = 5, height = 4)
plot_history(trial_sim1, x_value = "total n", y_value = "sum ys all") + plot_rules
dev.off()

# example 2
set.seed(3)
theta2 <- rnorm(1, mean=mu_x-mu_y, sd=theta_sd)
trial_setup2 <- setup_trial_norm(arms=c("Treatment X", "Control Y"), true_ys=c(mu_x, mu_x-theta2), sds=c(sigma, sigma), 
                                data_looks=c(0.5*n, n), highest_is_best=TRUE, control="Control Y",
                                soften_power=1, inferiority=c(0, 0.05), superiority=c(1, 0.95))
trial_sim2 <- run_trial(trial_spec=trial_setup2, seed=3, sparse=F)
print(trial_sim2$trial_res$final_alloc[1])
print(trial_sim2$trial_res$ns_all[1])
print(trial_sim2$all_looks[[1]]$sum_ys[1]-trial_sim2$all_looks[[1]]$sum_ys[2])
pdf(file = "plots/rarex2alln.pdf", width = 5, height = 4)
plot_history(trial_sim2, x_value = "total n", y_value = "n all") + plot_rules
dev.off()
pdf(file = "plots/rarex2prob.pdf", width = 5, height = 4)
plot_history(trial_sim2, x_value = "total n", y_value = "prob") + plot_rules
dev.off()
pdf(file = "plots/rarex2sumys.pdf", width = 5, height = 4)
plot_history(trial_sim2, x_value = "total n", y_value = "sum ys all") + plot_rules
dev.off()

# overall behaviour of example
rejections_rar <- rep(0, m)
final_allocs <- rep(0, m)
overall_allocs <- rep(0, m)
diff_in_sum_ys_interim <- rep(0, m)

for (i in 1:m){
  set.seed(i)
  theta <- rnorm(1, mean=mu_x-mu_y, sd=theta_sd)
  trial_setup <- setup_trial_norm(arms=c("Treatment X", "Control Y"), true_ys=c(mu_x, mu_x-theta), sds=c(sigma, sigma), 
                                  data_looks=c(0.5*n, n), highest_is_best=TRUE, control="Control Y",
                                  soften_power=1, inferiority=c(0, 0.05), superiority=c(1, 0.95))
  trial_sim <- run_trial(trial_spec=trial_setup, seed=i, sparse=F)
  status_x <- trial_sim$trial_res$final_status[1]
  rejections_rar[i] <- ifelse(mu_x > mu_y & status_x=="superior" | mu_y >= mu_x & status_x=="inferior", 1, 0)
  alloc_x <- trial_sim$trial_res$final_alloc[1]
  final_allocs[i] <- ifelse(mu_x > mu_y, alloc_x, 1-alloc_x)
  overall_allocs[i] <- trial_sim$trial_res$ns_all[1]/n
  sum_ys <- trial_sim$all_looks[[1]]$sum_ys
  diff_in_sum_ys_interim[i] <- sum_ys[1] - sum_ys[2]
}
comp_1 <- mean(rejections_rar)
df <- data.frame(final_allocs, overall_allocs, diff_in_sum_ys_interim)

pdf(file = "plots/interimprobhist.pdf", width = 5, height = 4)
ggplot(df) + geom_histogram(aes(final_allocs), fill="grey", colour="black", bins=30) + plot_rules + 
  xlab("r_x calculated at interim analysis") + scale_x_continuous(breaks = seq(0, 1, by = 0.1))
dev.off()
pdf(file = "plots/overallprobhist.pdf", width = 5, height = 4)
ggplot(df) + geom_histogram(aes(overall_allocs), fill="grey", colour="black", bins=30) + plot_rules + 
  xlab("overall R_x") + scale_x_continuous(breaks = seq(0, 1, by = 0.1))
dev.off()
pdf(file = "plots/interimdiffsumyshist.pdf", width = 5, height = 4)
ggplot(df) + geom_histogram(aes(diff_in_sum_ys_interim), fill="grey", colour="black", bins=30) + plot_rules + 
  xlab("diff in sum ys at interim") + scale_x_continuous(limits = symmetric_limits)
dev.off()

# fixed allocation components
r_xs <- seq(0.1, 0.9, 0.1)
comp_1s <- rep(0, length(r_xs))
for (i in 1:length(r_xs)){
  rejections <- rep(0, m)
  r_x <- r_xs[i]
  theta_hat_sd <- sqrt((sigma^2)/(n*r_x*(1-r_x)))
  for (j in 1:m){
    set.seed(j)
    theta <- rnorm(1, mu_x-mu_y, theta_sd)
    X = rnorm(r_x*n, mean=mu_x, sd = sigma)
    Y = rnorm((1-r_x)*n, mean=mu_x-theta, sd = sigma)
    X_bar = mean(X)
    Y_bar = mean(Y)
    theta_hat <- X_bar - Y_bar
    z_stat <- theta_hat/theta_hat_sd
    p_value <- pnorm(z_stat, lower.tail = FALSE)
    rejections[j] <- ifelse(p_value < 0.05, 1, 0)
  }
  comp_1s[i] <- mean(rejections)
}

table_df <- data.frame(r_xs, comp_1s)
table_df

# fixed behaviour
fixed1 <- fixed_compare(seq(-40, 40, 5), 28, 500, 200, 3)
fixed1
fixed2 <- fixed_compare(10, c(1, seq(5, 60, 5)), 500, 200, 3)
fixed2

# rar behaviour
rar1 <- rar_compare(seq(-40, 40, 5), 28, 500, 200, 3)
rar1
pdf(file = "plots/thetacomp1plot.pdf", width = 5, height = 4)
ggplot(rar1) + geom_line(aes(x=theta_plot, y=comp1_plot))
dev.off()
pdf(file = "plots/thetacomp2plot.pdf", width = 5, height = 4)
ggplot(rar1) + geom_line(aes(x=theta_plot, y=comp2_plot))
dev.off()
rar2 <- rar_compare(10, c(1, seq(5, 60, 5)), 500, 200, 3)
rar2
pdf(file = "plots/sigmacomp1plot.pdf", width = 5, height = 4)
ggplot(rar2) + geom_line(aes(x=sigma_plot, y=comp1_plot))
dev.off()
pdf(file = "plots/sigmacomp2plot.pdf", width = 5, height = 4)
ggplot(rar2) + geom_line(aes(x=sigma_plot, y=comp2_plot))
dev.off()

# example of gain comparison for one case
pdf(file = "plots/gainplot1.pdf", width = 5, height = 4)
loop_compare_fast_plot(20, 10, 28, 100, 200, 0.5, 0.05, cap=F, opt=F) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

pdf(file = "plots/gainplot2.pdf", width = 5, height = 4)
loop_compare_fast_plot(30, 10, 28, 100, 200, 0.5, 0.05, cap=F, opt=F) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

pdf(file = "plots/gainplot3.pdf", width = 5, height = 4)
loop_compare_fast_plot(15, 10, 28, 100, 200, 0.5, 0.05, cap=F, opt=F) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

pdf(file = "plots/gainplot4.pdf", width = 5, height = 4)
loop_compare_fast_plot(11, 10, 70, 100, 200, 0.5, 0.05, cap=F, opt=F) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

# example of gain comparison for one case with opt
pdf(file = "plots/gainplot1opt.pdf", width = 5, height = 4)
loop_compare_fast_plot(20, 10, 28, 100, 200, 0.5, 0.05, cap=F, opt=T) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

pdf(file = "plots/gainplot2opt.pdf", width = 5, height = 4)
loop_compare_fast_plot(30, 10, 28, 100, 200, 0.5, 0.05, cap=F, opt=T) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

pdf(file = "plots/gainplot3opt.pdf", width = 5, height = 4)
loop_compare_fast_plot(15, 10, 28, 100, 200, 0.5, 0.05, cap=F, opt=T) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

pdf(file = "plots/gainplot4opt.pdf", width = 5, height = 4)
loop_compare_fast_plot(11, 10, 70, 100, 200, 0.5, 0.05, cap=F, opt=T) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

# gain plots with capped opt
pdf(file = "plots/gainplot1cap.pdf", width = 5, height = 4)
loop_compare_fast_plot(20, 10, 28, 100, 200, 0.5, 0.1, cap=T, opt=T) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

pdf(file = "plots/gainplot2cap.pdf", width = 5, height = 4)
loop_compare_fast_plot(30, 10, 28, 100, 200, 0.5, 0.05, cap=T, opt=T) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

pdf(file = "plots/gainplot3cap.pdf", width = 5, height = 4)
loop_compare_fast_plot(20, 10, 60, 100, 200, 0.5, 0.05, cap=T, opt=T) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

pdf(file = "plots/gainplot4cap.pdf", width = 5, height = 4)
loop_compare_fast_plot(11, 10, 70, 100, 200, 0.5, 0.05, cap=T, opt=T) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) 
dev.off()

# gain plots fixed and opt only
pdf(file = "plots/gainplot1fixed.pdf", width = 5, height = 4)
loop_compare_fast_plot(20, 10, 28, 100, 200, 0.5, 0.05, cap=F, opt=T) + plot_rules + 
  scale_color_manual(values = c("optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_linetype_manual(values=c("optimal fixed"=2, "fixed"=3))
dev.off()

pdf(file = "plots/gainplot2fixed.pdf", width = 5, height = 4)
loop_compare_fast_plot(30, 10, 28, 100, 200, 0.5, 0.05, cap=F, opt=T) + plot_rules + 
  scale_color_manual(values = c("optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_linetype_manual(values=c("optimal fixed"=2, "fixed"=3)) 
dev.off()

pdf(file = "plots/gainplot3fixed.pdf", width = 5, height = 4)
loop_compare_fast_plot(15, 10, 28, 100, 200, 0.5, 0.05, cap=F, opt=T) + plot_rules + 
  scale_color_manual(values = c("optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_linetype_manual(values=c("optimal fixed"=2, "fixed"=3)) 
dev.off()

pdf(file = "plots/gainplot4fixed.pdf", width = 5, height = 4)
loop_compare_fast_plot(11, 10, 70, 100, 200, 0.5, 0.05, cap=F, opt=T) + plot_rules + 
  scale_color_manual(values = c("optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_linetype_manual(values=c("optimal fixed"=2, "fixed"=3)) 
dev.off()

# wider characterisation of behaviour
plot <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 500, 200, 0.3, cap=F, opt=F) + plot_rules
plot2 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 500, 200, 0.5, cap=F, opt=F) + plot_rules
plot3 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 500, 200, 0.7, cap=F, opt=F) + plot_rules
plot4 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 500, 200, 0.8, cap=F, opt=F) + plot_rules
plot5 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 500, 200, 0.9, cap=F, opt=F) + plot_rules
plot6 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 500, 200, 1.0, cap=F, opt=F) + plot_rules

pdf(file = "plots/bigplot0_3.pdf", width = 5, height = 4)
plot + scale_color_manual(values = c("adaptive" = "hotpink", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_5.pdf", width = 5, height = 4)
plot2 + scale_color_manual(values = c("adaptive" = "hotpink", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_7.pdf", width = 5, height = 4)
plot3 + scale_color_manual(values = c("adaptive" = "hotpink", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_8.pdf", width = 5, height = 4)
plot4 + scale_color_manual(values = c("adaptive" = "hotpink", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_9.pdf", width = 5, height = 4)
plot5 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot1_0.pdf", width = 5, height = 4)
plot6 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "fixed"=17))
dev.off()

plot7 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.1, cap=F, opt=T) + plot_rules
plot8 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.3, cap=F, opt=T) + plot_rules
plot9 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.5, cap=F, opt=T) + plot_rules
plot10 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.7, cap=F, opt=T) + plot_rules
plot11 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.9, cap=F, opt=T) + plot_rules
plot12 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.8, cap=F, opt=T) + plot_rules
plot13 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 1.0, cap=F, opt=T) + plot_rules

pdf(file = "plots/bigplot0_1opt.pdf", width = 5, height = 4)
plot7 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_3opt.pdf", width = 5, height = 4)
plot8 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_5opt.pdf", width = 5, height = 4)
plot9 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_7opt.pdf", width = 5, height = 4)
plot10 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_9opt.pdf", width = 5, height = 4)
plot11 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_8opt.pdf", width = 5, height = 4)
plot12 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot1_0opt.pdf", width = 5, height = 4)
plot13 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3", "inconclusive"="black")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17, "inconclusive"=18))
dev.off()

# capping optimisation function
plot14 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.1, cap=T, opt=T) + plot_rules
plot15 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.3, cap=T, opt=T) + plot_rules
plot16 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.5, cap=T, opt=T) + plot_rules
plot17 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.7, cap=T, opt=T) + plot_rules
plot18 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.9, cap=T, opt=T) + plot_rules
plot19 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 0.8, cap=T, opt=T) + plot_rules
plot20 <- make_best_plot(seq(0, 50, 5), c(1, 10, 20, 30, 40, 50, 60, 70, 80), 300, 200, 1.0, cap=T, opt=T) + plot_rules

pdf(file = "plots/bigplot0_1cap.pdf", width = 5, height = 4)
plot14 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3"), aesthetics = c("color", "fill")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17)) + geom_encircle(aes(fill=best),alpha=0.3)
dev.off()
pdf(file = "plots/bigplot0_3cap.pdf", width = 5, height = 4)
plot15 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17)) + 
  geom_encircle(data=plot15$data[plot15$data$best=="RAR"], aes(x=theta, y=sigma))
dev.off()
pdf(file = "plots/bigplot0_5cap.pdf", width = 5, height = 4)
plot16 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_7cap.pdf", width = 5, height = 4)
plot17 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_9cap.pdf", width = 5, height = 4)
plot18 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot0_8cap.pdf", width = 5, height = 4)
plot19 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17))
dev.off()
pdf(file = "plots/bigplot1_0cap.pdf", width = 5, height = 4)
plot20 + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3", "inconclusive"="black")) + 
  scale_shape_manual(values = c("adaptive" = 15, "optimal fixed" = 16, "fixed"=17, "inconclusive"=18))
dev.off()




# capping optimisation function
pdf(file = "plots/gainplot1cap.pdf", width = 5, height = 4)
loop_compare_fast_plot(20, 10, 28, 100, 200, 0.5, 0.05, T, T) + plot_rules + 
  scale_color_manual(values = c("adaptive" = "hotpink", "optimal fixed"="dodgerblue2", "fixed"="green3"))
dev.off()

#
test <- loop_compare_fast_plot(30, 0, 20, 100, 200, 0.5, 0.05, F, T) + plot_rules
test

loop_compare_fast_plot_best_trial(15, 40, 100, 200, 0.9, F, T)
loop_compare_fast_plot(30, 0, 20, 100, 200, 0.5, 0.05, F, T)

