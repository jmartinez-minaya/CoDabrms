### --- Simulation for Compositional data using splines --- ####



### --- 1. Loading libraries --- ####
library(compositions)
library(dplyr)
library(MASS)
library(brms)
library(ggplot2)
library(brms)
library(pbapply)
library(gridExtra)
library(tidybayes)
library(mgcv)
library(viridis)
library(patchwork)

### --- 2. Functions to simulate and fit --- ####

### ----- 2.1. Simulation with splines --- ####
#Function to simulate. It simulates a one and two dimensional spline
simulation_spline <- function(beta = matrix(c(-1, 1, 2 ,1), 
                                            byrow = TRUE, 
                                            ncol = 2), 
                       sigma = c(10, 9), 
                       rho   = 0.9,
                       sim1,
                       sim2,
                       n,
                       x = NULL)
{
  dim(beta)
  D1 <- dim(beta)[1]
  D2 <- dim(beta)[2]
  x <- cbind(1,x)
  
  #Simulating
  set.seed(314)
  
  mu <- apply(beta, 1, function(y)x %*% y)
  
  #We add the spline term and remove the fixed term
  #mu <- mu
  mu <- mu + sim1[, c("sim11", "sim12")] + sim2[, c("sim21", "sim22")]
 # mu <- mu + sim2$sim2
  
  #mu <- mu + sim1$sim1 + sim2$data$f
  
  
  #Correlation terms
  Sigma <- rho*matrix(sigma, ncol = 1) %*% matrix(sigma, nrow = 1)
  diag(Sigma) <- diag(Sigma)/rho
  
  mu %>%
    apply(., 1, function(x) 
      MASS::mvrnorm( n  = 1,
               mu = x,
               Sigma = Sigma)) %>%
    t(.)-> y
  
  z <- ilrInv(y)
  y_ilr <- ilr(z)
  colnames(x) <- c("intercept", paste0("x", 1:(D2 -1)))
  
  #We add covariates corresponding to splines
  x <- cbind(x, x_s1 = sim1$x_s1, x_s2 = sim2$x_s2, x_s3 = sim2$x_s3)
  data <- data.frame(y_ilr, x[,-1])
  colnames(data)[1:D1] <- paste0("yilr", c(1:D1))
  data
}


### ----- 2.2. Fitting using splines --- ####
#Function to fit. Note that you can include also fixed effects
fit_spline <- function(data, covariate = TRUE, formula1 = "x1 + x2")
{
  if(covariate == FALSE)
  {
    bf_ilr_1 <- bf(yilr1 ~ 1) + gaussian()
    bf_ilr_2 <- bf(yilr2 ~ 1) + gaussian()
  }else{
    bf_ilr_1 <- bf(as.formula(paste0("yilr1 ~ ", formula1))) + gaussian()
    bf_ilr_2 <- bf(as.formula(paste0("yilr2 ~ ", formula1))) + gaussian()
    #bf_ilr_3 <- bf(as.formula(paste0("yilr3 ~", formula1))) + gaussian()
    
  }
  #Fitting
  
  bprior <- c(#prior(normal(0, 10), class = Intercept, resp = yilr1),
              prior(normal(0, 10), class = b, resp = yilr1),
              #prior(normal(0, 10), class = Intercept, resp = yilr2),
              prior(normal(0, 10), class = b, resp = yilr2))
  
  model1 <- brm(
    bf_ilr_1 + bf_ilr_2 +  set_rescor(TRUE),
    data   = data,            
    prior  = bprior,
    chains = 2, 
    cores  = 2,
    iter   = 5000,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE)
  )
  
  # get_prior(bf_ilr_1 + bf_ilr_2 +  set_rescor(TRUE),
  #           data   = data)
  summary(model1)
  model1
}


### --- 3. Simulating --- ####
n <- 100
m <- 2
betas <- c(3, 0,
           -4, 0)

betas <- c(0, 0,
           0, 0)

sigma <- c(0.03, 0.05)
#sigma <- c(0.5, 0.1)

rho <- rep(0.5, m-1)
set.seed(10)
x <- matrix(runif(n, 0, 1), ncol = 1)

### ----- 3.1. One dimensional spline: one per category --- ####
set.seed(10)
x_s1 <- runif(n, 0, 1)
x_s1 <- x_s1[order(x_s1)]
#sim1 <- 5*(sin(x_s/(4*pi))) + rnorm(n, sd = 0.1)

#Sim coordinate 1
sim11 <- 0.2 * (x_s1)^11 * (10 * (1 - (x_s1)))^6 + 10 * (10 * x_s1)^3 * 
  (1 - x_s1)^10
sim11 <- sim11/8
sim11 <- (sim11 - mean(sim11))
sim1 <- data.frame(x_s1 = x_s1, sim11 = sim11)
#plot(sim1$x_s1, sim1$sim11)

#Sim coordinate 2
sim12 <- sin(2*pi*x_s1)/2
sim12 <- (sim12 - mean(sim12))
sim1 <- cbind(sim1, sim12 = sim12)
#plot(sim1$x_s1, sim1$sim12)

### ----- 3.2. Two dimensional spline: one different per category --- ####
set.seed(10)
sim21 <- mgcv::gamSim(eg=2, n = n, dist="normal", scale=2, verbose = TRUE)
#contour(sim21$truth$x, sim21$truth$z, sim21$truth$f)
sim21_center11 <- sim21$truth$f - mean(sim21$truth$f)
#contour(sim21$truth$x, sim21$truth$z, sim21_center11)

sim21_center <- sim21$data$f - mean(sim21$data$f)
sim2 <- data.frame(x_s2 = sim21$data$x, x_s3 = sim21$data$z, sim21 = sim21_center)


#Plotting
plot_data_ilr1_sim2 <- data.frame(xs2 = sim21$pr$x, 
                                  xs3 = sim21$pr$z,
                                  z   = as.numeric(sim21_center11))
p3_data <- ggplot(plot_data_ilr1_sim2, 
                  aes(xs2, xs3, z = z)) +
  geom_raster(aes(fill = z)) +
  geom_contour(colour = "white") +
  theme_bw() +
  xlab(expression(xs[2])) +
  ylab(expression(xs[3])) +
  theme(legend.position="bottom")
p3_data






set.seed(10)
source("functions/gamSim2.R")
sim22 <- gamSim2(eg=2, n = n, dist="normal", scale=300, verbose = TRUE)
contour(sim22$truth$x, sim22$truth$z, sim22$truth$f)
sim22_center11 <- sim22$truth$f - mean(sim22$truth$f)
contour(sim22$truth$x, sim22$truth$z, sim22_center11)

sim22_center <- sim22$data$f - mean(sim22$data$f)
sim2 <- cbind(sim2, sim22 = sim22_center)

#Plotting
plot_data_ilr2_sim2 <- data.frame(xs2 = sim22$pr$x, 
                                  xs3 = sim22$pr$z,
                                  z   = as.numeric(sim22_center11))
p4_data <- ggplot(plot_data_ilr2_sim2, 
                  aes(xs2, xs3, z = z)) +
  geom_raster(aes(fill = z)) +
  geom_contour(colour = "white") +
  theme_bw() +
  xlab(expression(xs[2])) +
  ylab(expression(xs[3])) +
  theme(legend.position="bottom")

### ----- 3.3. Simulating --- ####
data_spline <- simulation_spline(beta = matrix(betas, byrow = TRUE, ncol = dim(x)[2] + 1 ),
                         sigma = sigma, 
                         rho = rho, 
                         n = n, 
                         x = x,
                         sim1 = sim1,
                         sim2 = sim2)


data_spline$x1
data_p <- data_spline[,c("yilr1", "yilr2")] %>% 
  compositions::ilrInv(.) %>%
  .[,1:3] %>% as.data.frame(.)*100
colnames(data_p) <- paste0("y", 1:3)




### ----- 3.4. Ggtern --- ####
library(ggtern)
var = ""
p_data <- ggtern(data = data_p, aes(x = y1, y = y2, z = y3)) +
  geom_point(size = 3) + 
  labs(x = "y1", y = "y2", z = "y3") +
  #theme_clockwise() +
  theme_showarrows() +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold")) +
  #guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
  scale_color_viridis(option = "E", direction = -1, name = var)

p_data


pdf("figures/data_spline.pdf", width = 6, height = 5)
p_data
dev.off()

### --- 4. Fitting --- ####
mod_spline_1 <- fit_spline(data = data_spline, 
                    formula1 = "-1 + t2(x_s2, x_s3, m = 2, k = 10, bs = 'ps') + s(x_s1, bs = 'ps', k = 10)")

summary(mod_spline_1)
loo1 <- loo::loo(mod_spline_1, moment_match = TRUE, reloo = TRUE)
#saveRDS(mod_spline_1, "rds/simulation_mod1_spline.rds")

pdf("fitting_splines.pdf")
print(plot(conditional_smooths(mod_spline_1), ask = FALSE))
dev.off()


### --- 5. Posterior distributions --- ####
### ----- 5.1. Posterior distributions: one dimension --- ####

source("functions/posterior_smooths_plot.R")
p1 <- posterior_smooths_plot(mod = mod_spline_1,
                             smooth = 's(x_s1,bs="ps",k=10)',
                             resp = "yilr1",
                             x = "x_s1",
                             lab.x = expression(xs[1]),
                             lab.y = expression(paste(s^1, "(x", s[1], ")")))


p2 <- posterior_smooths_plot(mod = mod_spline_1,
                             smooth = 's(x_s1,bs="ps",k=10)',
                             resp = "yilr2",
                             x = "x_s1",
                             lab.x = expression(xs[1]),
                             lab.y = expression(paste(s^2, "(x", s[1], ")")))



summary(mod_spline_1)
#plot(conditional_effects(mod_spline_1))
p1 <- p1 + geom_point(data = sim1, aes(x = x_s1, y = sim11), size = 0.5)
p2 <- p2 + geom_point(data = sim1, aes(x = x_s1, y = sim12), size = 0.5)

pdf("spline_ilr1.pdf", width = 8, height = 5)
p1 / p2 +   plot_annotation(tag_levels = 'a') 

dev.off()

### ----- 5.2. Posterior distributions: two dimensions --- ####
p3 <- posterior_smooths_plot(mod = mod_spline_1,
                             smooth = 't2(x_s2,x_s3,m=2,k=10,bs="ps")',
                             resp = "yilr1",
                             x = c("x_s2", "x_s3"),
                             lab.x  = expression(xs[2]),
                             lab.y  = expression(xs[3]),
                             lab.z  = expression(paste(t^1, "(", xs[2],",", xs[3], ")")))

p4 <- posterior_smooths_plot(mod = mod_spline_1,
                             smooth = 't2(x_s2,x_s3,m=2,k=10,bs="ps")',
                             resp   = "yilr2",
                             x      = c("x_s2", "x_s3"),
                             lab.x  = expression(xs[2]),
                             lab.y  = expression(xs[3]),
                             lab.z  = expression(paste(t^2, "(", xs[2],",", xs[3], ")")))



# Calcular rango común con margen
common_min <- min(p3[[2]][1], p4[[2]][1], sim2$sim21, sim2$sim22) - 0.05
common_max <- max(p3[[2]][2], p4[[2]][2], sim2$sim21, sim2$sim22) + 0.05
common_limits <- c(common_min, common_max)

# Escala común para todos los gráficos
shared_scale <- scale_fill_viridis(name = "", limits = common_limits)

# Aplicar a los cuatro gráficos
p3_data <- p3_data +  theme(legend.position="bottom",
                            panel.background = element_blank(),
                            legend.key.height = unit(0.4, "cm"),
                            legend.key.width = unit(2, "cm")) + shared_scale
p4_data <- p4_data + theme(legend.position="bottom",
                           panel.background = element_blank(),
                           legend.key.height = unit(0.4, "cm"),
                           legend.key.width = unit(2, "cm")) + shared_scale
p3[[1]] <- p3[[1]] + shared_scale
p4[[1]] <- p4[[1]] + shared_scale





pdf("spline_ilr2.pdf", width = 8, height = 7)
(p3_data | p3[[1]]) /
  (p4_data | p4[[1]]) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
dev.off()





### ----- 5.3. Posterior distributions: fixed effect intercept --- ####
# #Intercepts
# data_plot <- mod_spline_1 %>%
#   gather_draws(b_yilr1_Intercept, b_yilr2_Intercept,)
# colnames(data_plot) <- c("chain", "iteration", "draw", "variable_p",  "value_p")
# 
# values_plot <- data.frame(values_int = betas[seq(1,length(betas), 2)], variable_p = c("b_yilr1_Intercept",  "b_yilr2_Intercept"))
# 
# #pdf("figures/sim_spline_intercept.pdf", width = 4, height = 3)
# ggplot(data = data_plot) +
#   geom_histogram(aes(x = value_p, y = ..density.., group = variable_p),
#                  color = "#e9ecef", alpha=0.6, position = 'identity', bins = 30, fill = "cornflowerblue") +
#   theme_bw() +
#   geom_vline(data = values_plot, aes(xintercept = values_int),
#              color = "blue4", size=1) +
#   facet_wrap( ~ variable_p, scales = "free")
# #dev.off()

### ----- 5.3. Fitted values --- ####
# Data_inverse
z_comp <- ilrInv(data_spline[, c("yilr1", "yilr2")])
z_comp <- z_comp %>% as.data.frame() 
colnames(z_comp) <- paste0("Cat", 1:3)  
z_comp %>%
  tidyr::pivot_longer(., cols = 1:3,  names_to = "Category", values_to = "prop_data") -> z_comp

# Prediction
y_pred <- predict(mod_spline_1)
z_comp_pred <- ilrInv(y_pred[,"Estimate",])
z_comp_pred <- z_comp_pred %>% as.data.frame() 
colnames(z_comp_pred) <- paste0("Cat", 1:3)  
z_comp_pred %>%
  tidyr::pivot_longer(., cols = 1:3,  names_to = "Category", values_to = "prop_pred") -> z_comp_pred

#Join
z_comp_total <- cbind(z_comp, prop_pred = z_comp_pred$prop_pred)

z_comp_total$Category <- factor(z_comp_total$Category, 
                                levels = c(paste0("Cat",1:3)),
                                ordered = TRUE,
                                labels = c(paste0("Category-", 1:3)))


# Plotting
pdf("figures/sim_spline_fitted.pdf", width = 8, height = 3)
ggplot(z_comp_total, aes(x = prop_data, y = prop_pred)) +
  geom_smooth(method = lm, se = TRUE, level = 0.99) +
  geom_point() + 
  scale_y_continuous(breaks= scales::pretty_breaks()) +
  scale_x_continuous(breaks= scales::pretty_breaks()) +
  facet_wrap(~Category, scales = "free", labeller = label_parsed) +
  theme_bw() +
  xlab("Data proportions") +
  ylab("Predicted proportions")
dev.off()


### ----- 5.4. Posterior distributions: variances and correlation parameters --- ####
#Variances
#correlation parameters
b <- summary(mod_spline_1)
b$spec_pars
data_plot <- mod_spline_1 %>%
  gather_draws(sigma_yilr1, sigma_yilr2, rescor__yilr1__yilr2)  
colnames(data_plot) <- c("chain", "iteration", "draw", "variable_p",  "value_p")

values_plot <- data.frame(values_int = c(sigma, rho), variable_p = c(paste0("sigma_yilr", 1:2), 
                                                                     "rescor__yilr1__yilr2"))

data_plot$variable_p <- factor(data_plot$variable_p, 
                               levels = c(paste0("sigma_yilr", 1:2), "rescor__yilr1__yilr2"),
                               ordered = TRUE,
                               labels = c(paste(expression(sigma[1])), paste(expression(sigma[2])), paste(expression(rho[12]))))

values_plot$variable_p <- factor(values_plot$variable_p, 
                                 levels = c(paste0("sigma_yilr", 1:2), "rescor__yilr1__yilr2"),
                                 ordered = TRUE,
                                 labels = c(paste(expression(sigma[1])), paste(expression(sigma[2])), paste(expression(rho[12]))))



pdf("sim_splines_variances.pdf", width = 8, height = 3)
ggplot(data = data_plot) +
  geom_histogram(aes(x = value_p, y = ..density.., group = variable_p),
                 color = "#e9ecef", alpha=0.6, position = 'identity', bins = 15, fill = "cornflowerblue") +
  theme_bw() +
  facet_wrap( ~ variable_p, scales="free",
              labeller = label_parsed) +
  geom_vline(data = values_plot, aes(xintercept = values_int),
             color = "blue4", size=1) +
  xlab("x")
dev.off()

### --- 6. Fitting more models to compute R-squared --- ####
# mod_spline_1 <- fit_spline(data = data_spline, 
#                            formula1 = "t2(x_s2, x_s3, m = 2, k = 40, bs = 'ps') + s(x_s1, bs = 'ps', k = 40)")

mod_spline_2 <- fit_spline(data = data_spline, 
                           formula1 = "s(x_s1, bs = 'ps', k = 10)")

mod_spline_3 <- fit_spline(data = data_spline, 
                           formula1 = "t2(x_s2, x_s3, m = 2, k = 10, bs = 'ps')")

mod_spline_4 <- fit_spline(data = data_spline, 
                           formula1 = "s(x_s1, bs = 'ps', k = 5)")

mod_spline_5 <- fit_spline(data = data_spline, 
                           formula1 = "t2(x_s2, x_s3, m = 2, k = 5, bs = 'ps')")

mod_spline_6 <- fit_spline(data = data_spline, 
                           formula1 = "t2(x_s2, x_s3, m = 2, k = 5, bs = 'ps') + s(x_s1, bs = 'ps', k = 5)")    



### --- 7. Computing R-squared and WAIC --- ####
rm(var)
source("functions/bayes_R2.R")
m1 <- 6

# List of models
list_mod <- list(mod_spline_1, mod_spline_6, mod_spline_3, mod_spline_5,
                 mod_spline_2, mod_spline_4)

# --- BM-CoDa-R²: Based on variance ---
list_mod %>%
  sapply(., bayes_R2.CoDa, resp_names = c("yilr1", "yilr2")) %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols = 1:m1, names_to = "fit", values_to = "R_squared") -> r_squared

r_squared$fit <- as.factor(r_squared$fit)
levels(r_squared$fit) <- paste0("M", 1:m1)

r_squared %>%
  group_by(fit) %>%
  summarise(
    m = mean(R_squared),
    sd = sd(R_squared)
  ) -> r_squared_mean

# --- BR-CoDa-R²: Based on residuals ---
list_mod %>%
  sapply(., bayes_R2.CoDa, type = "res") %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols = 1:m1, names_to = "fit", values_to = "R_squared") -> r_squared_res

r_squared_res$fit <- as.factor(r_squared_res$fit)
levels(r_squared_res$fit) <- paste0("M", 1:m1)

r_squared_res %>%
  group_by(fit) %>%
  summarise(
    m = mean(R_squared),
    sd = sd(R_squared)
  ) -> r_squared_res_mean

# --- Combine all R² values ---
r_squared_all <- cbind(r_squared_res, r_sq = "r_sq_res") %>%
  rbind(., cbind(r_squared, r_sq = "r_sq_var"))

r_squared_all$r_sq <- factor(r_squared_all$r_sq,
                             levels = c("r_sq_res", "r_sq_var"),
                             ordered = TRUE,
                             labels = c(expression(paste("BR-CoDa-", R^2)),
                                        expression(paste("BM-CoDa-", R^2))))

# --- Combine mean and SD summaries ---
r_squared_mean_data <- rbind(
  cbind(r_squared_res_mean, r_sq = 1),
  cbind(r_squared_mean, r_sq = 2)
)

r_squared_mean_data$r_sq <- factor(r_squared_mean_data$r_sq,
                                   levels = c(1, 2),
                                   ordered = TRUE,
                                   labels = c(expression(paste("BR-CoDa-", R^2)),
                                              expression(paste("BM-CoDa-", R^2))))

r_squared_mean_data$fit <- as.factor(r_squared_mean_data$fit)
levels(r_squared_mean_data$fit) <- paste0("M", 1:m1)

# --- Compute WAIC for each model ---
waic_res <- list_mod %>%
  pblapply(., function(x) {
    waic(x)
  })

waic_res1 <- waic_res %>%
  sapply(., function(x) x$estimates["waic", "Estimate"]) %>%
  round(., 3)

# --- Format BR and BM R² as mean (SD) strings ---
BR_values <- paste0(
  round(r_squared_mean_data$m[1:m1], 3),
  " (", round(r_squared_mean_data$sd[1:m1], 3), ")"
)

BM_values <- paste0(
  round(r_squared_mean_data$m[(m1 + 1):(2 * m1)], 3),
  " (", round(r_squared_mean_data$sd[(m1 + 1):(2 * m1)], 3), ")"
)

# --- Final table combining all results ---
total_models <- data.frame(
  Model = paste0("M", 1:m1),
  WAIC  = waic_res1,
  BR    = BR_values,
  BM    = BM_values
)

# Display table (optional)
xtable::xtable(total_models, digits = 3)


### --- 8. Plotting R-squared --- ####
p1 <- r_squared %>%
  ggplot( aes(x = R_squared, fill = fit, y = ..density..)) +
  geom_histogram( color="#e9ecef", alpha = 0.9, position = 'identity', bins = 150)+ 
  # scale_fill_brewer(palette="Dark2") +
  geom_vline(data =  r_squared_mean_data, 
             aes(xintercept = m, color = fit),
             linetype = "dashed") +  #scale_color_brewer(palette="Dark2")+
  #scale_color_brewer(palette="Dark2") + 
  theme_bw() +
  labs(fill="") +
  theme(legend.position="bottom") +
  guides(fill = guide_legend( ncol = 7),
         col = "none")+
  xlim(c(0,1.01)) +
  xlab(expression(CoDa-R^2)) +
  facet_wrap(~r_sq, ncol = 1, labeller = label_parsed)

p1

pdf("r-squared_sim_splines.pdf", width = 10, height = 6)
p1 
dev.off()




### --- 9. Computing probabilities to compare --- ####
lapply(paste0("M", 1:6), function(x){ r_squared_all %>% 
    dplyr::filter(fit == x) %>%
    dplyr::filter(r_sq == 'paste("BR-CoDa-", R^2)') %>%
    dplyr::select(R_squared) %>% pull(.)
}) -> m_sim_br

lapply(paste0("M", 1:6), function(x){ r_squared_all %>% 
    dplyr::filter(fit == x) %>%
    dplyr::filter(r_sq == 'paste("BM-CoDa-", R^2)') %>%
    dplyr::select(R_squared) %>% pull(.)
}) -> m_sim_bm


all_mod <- expand.grid(paste0("M", 1:6), paste0("M", 1:6))
all_mod <- expand.grid( 1:6, 1:6)

# probs_br
probs_br <- apply(all_mod, 1, function(x, m_sim){
  res <- (m_sim[[x[1]]] > m_sim[[x[2]]]) %>% table(.)/5000
  if(length(names(res)) == 2){
    res <- res[names(res)==TRUE]
  }else if(length(names(res)) == 1 & names(res)==TRUE){
    res <- res
  }else{
    res <- 0
  }
  
}, m_sim = m_sim_br) 

probs_bm <- apply(all_mod, 1, function(x, m_sim){
  res <- (m_sim[[x[1]]] > m_sim[[x[2]]]) %>% table(.)/5000
  if(length(names(res)) == 2){
    res <- res[names(res)==TRUE]
  }else if(length(names(res)) == 1 & names(res)==TRUE){
    res <- res
  }else{
    res <- 0
  }
  
}, m_sim = m_sim_bm)


cbind(all_mod, probs_br = round(probs_br,3), 
      probs_bm = round(probs_bm,3))






