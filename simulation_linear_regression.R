# Simulation examples: Linear regression

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
library(dplyr)


### --- 2. Loading the functions for computing R2 --- ####
source("functions/bayes_R2.R")

### --- 3. Functions to simulate and fit --- ####
simulation <- function(beta = matrix(c(-1, 1, 2 ,1), byrow = TRUE, ncol = 2), 
                       sigma = c(10, 9), 
                       rho   = 0.9,
                       n,
                       x = NULL)
{
  dim(beta)
  D1 <- dim(beta)[1]
  D2 <- dim(beta)[2]
  if(is.null(x)){
    set.seed(314)
    x <- cbind(1, matrix(rnorm((D2-1)*n), ncol = D2-1))
  }else{
    x <- cbind(1,x)
  }
  #Simulating
  set.seed(314)
  
  mu <- apply(beta, 1, function(y)x %*% y)
  Sigma <- matrix(sigma, ncol = 1) %*% matrix(sigma, nrow = 1)
  
  for(i in 1:length(rho))
  {
    for (j in 1:length(rho))
    {
      if(i == j)
      {
        Sigma[i,j] <- Sigma[i,j]
      }else{
        Sigma[i,j] <- rho[i + j - 2] * Sigma[i,j]
      }
    }
  }

  
  # Sigma <- rho*matrix(sigma, ncol = 1) %*% matrix(sigma, nrow = 1)
  # diag(Sigma) <- diag(Sigma)/rho
  
  mu %>%
    apply(., 1, function(x) 
      mvrnorm( n  = 1,
               mu = x,
               Sigma = Sigma)) %>%
    t(.)-> y
  
  z <- ilrInv(y)
  yilr <- ilr(z)
  colnames(x) <- c("intercept", paste0("x", 1:(D2 -1)))
  data <- data.frame(yilr, x[,-1])
  colnames(data)[1:D1] <- paste0("yilr", c(1:D1))
  data
}


fit <- function(data, covariate = TRUE, formula1 = "x1 + x2")
{
  if(covariate == FALSE)
  {
    bf_ilr_1 <- bf(yilr1 ~ 1) + gaussian()
    bf_ilr_2 <- bf(yilr2 ~ 1) + gaussian()
  }else{
    bf_ilr_1 <- bf(as.formula(paste0("yilr1 ~", formula1))) + gaussian()
    bf_ilr_2 <- bf(as.formula(paste0("yilr2 ~", formula1))) + gaussian()
    bf_ilr_3 <- bf(as.formula(paste0("yilr3 ~", formula1))) + gaussian()
    
  }
  #Fitting
  
  bprior <-  c(prior(normal(0, 10), class = Intercept, resp = yilr1),
               prior(normal(0, 10), class = b, resp = yilr1),
               prior(normal(0, 10), class = Intercept, resp = yilr2),
               prior(normal(0, 10), class = b, resp = yilr2),
               prior(normal(0, 10), class = Intercept, resp = yilr3),
               prior(normal(0, 10), class = b, resp = yilr3))
  
  formula.1 <- bf_ilr_1 + bf_ilr_2 + bf_ilr_3+  set_rescor(TRUE)
  model1 <- brm(
    bf_ilr_1 + bf_ilr_2 + bf_ilr_3+  set_rescor(TRUE),
    data   = data,            
    prior  = bprior,
    chains = 2, 
    cores  = 2,
    iter   = 5000,
  )
  
  
  summary(model1)
  model1
}



### --- 4. Simulating the data --- ####
#We fixed the variance-covariance function, and we modify the covariate $x$. 
#Then, we see how the R-squared works.

n <- 100
m <- 3
# betas <- c(-1, 1, 1, 2,
#            2 ,1, -1, 1,
#            1, 2, -2, 1)

betas <- c(1 ,-0.5, +1, -0.5,
           -0.5, 1, -1, 0,
           -2, 1, -0.5, -0.5)

sigma <- c(0.1, 0.05, 0.08)
#rho <- rep(0.9, m)
rho <- c(0.5,0.2,0.8)
#sigma <- c(0.5, 0.9, 1)


set.seed(100)
x <- matrix(rnorm(m*n, 10, 1), ncol = m)
data_test3 <- simulation(beta = matrix(betas, byrow = TRUE, ncol = m +1 ),
                         sigma = sigma, 
                         rho = rho, 
                         n = n, 
                         x = x)

data_test3[,c("x1", "x2", "x3")] %>% cor(.)



### ----- 3.2. Plotting in the Euclidean space --- ####
data_test3$yilr1
data_p <- data_test3[,c("yilr1", "yilr2", "yilr3")] %>% compositions::ilrInv(.) %>%
  .[,1:4] %>% as.data.frame(.) *100
colnames(data_p) <- paste0("y", 1:4)

summary(data_p)
### ----- 3.3. Ggtern in three dimensional space --- ####
library(ggtern)
# var = "Elevation (m)"
p_data <- ggtern(data = data_p, aes(x = y1, y = y2, z = y3,
                                   col = y4)) +
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
  scale_color_viridis(option = "E", direction = -1)
p_data


# pdf("figures/data_linear_regression.pdf", width = 6, height = 5)
# p_data
# dev.off()
### --- 4. Fitting the models --- ####
#We fit model with all the covariates
mod_test3_1 <- fit(data = data_test3, 
                   formula1 = "x1 + x2 + x3")
#mod_test3_2 <- fit(data = data_test3, formula1 = "x1 + x2")
mod_test3_3 <- fit(data = data_test3, 
                   formula1 = "x1 + x2")
#mod_test3_4 <- fit(data = data_test3, formula1 = "x1 + x3")
mod_test3_5 <- fit(data = data_test3, 
                   formula1 = "x1")
mod_test3_6 <- fit(data = data_test3, 
                   formula1 = "x2")
# mod_test3_7 <- fit(data = data_test3, formula1 = "x3")

summary(mod_test3_5)
summary(mod_test3_6)


#Covariate with noise
set.seed(20)
data_test3$x3 <- rnorm(n)
mod_test3_8 <- fit(data = data_test3, formula1 = "x1 + x2 + x3")
#mod_test3_9 <- fit(data = data_test3, formula1 = "x1 + x3")
mod_test3_10 <- fit(data = data_test3, formula1 = "x3")

### --- 5. Computing the R-squared --- ####
list_mod <- list(mod_test3_1, mod_test3_8, mod_test3_3, mod_test3_5, mod_test3_6, mod_test3_10)
m1 <- 6
list_mod %>%
  sapply(., bayes_R2.CoDa) %>%
  as.data.frame() %>%
  tidyr::pivot_longer(., cols = 1:m1,  names_to = "fit", values_to = "R_squared")-> r_squared


####################3

# summary(mod_test3_1)
# summary(mod_test3_8)
# summary(mod_test3_3)
# summary(mod_test3_5)
# summary(mod_test3_10)

r_squared$fit <- as.factor(r_squared$fit)
levels(r_squared$fit) <- paste0("M", 1:m1)

r_squared %>%
  group_by(fit) %>%
  summarise(m = mean(R_squared)) -> r_squared_mean


r_squared %>%
  group_by(fit) %>%
  summarise(sd = sd(R_squared)) -> r_squared_sd


# Based on residuals
list_mod %>%
  sapply(., bayes_R2.CoDa, type = "res") %>%
  as.data.frame() %>%
  tidyr::pivot_longer(., cols = 1:m1,  names_to = "fit", values_to = "R_squared")-> r_squared_res

r_squared_res$fit <- as.factor(r_squared_res$fit)
levels(r_squared_res$fit) <- paste0("M", 1:m1)


r_squared_res %>%
  group_by(fit) %>%
  summarise(m = mean(R_squared)) -> r_squared_res_mean

r_squared_res %>%
  group_by(fit) %>%
  summarise(sd = sd(R_squared)) -> r_squared_res_sd

#Plotting
# p1 <- r_squared %>%
#   ggplot( aes(x = R_squared, fill = fit, y = ..density..)) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 100) +
#   theme_bw() +
#   labs(fill="") +
#   #theme(legend.position="bottom") +
#   xlim(c(0,1.01)) +
#   xlab("CoDa-R-squared")
# 
# 
# 
# p2 <- r_squared_res %>%
#   ggplot( aes(x = R_squared, fill = fit, y = ..density..)) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 100) +
#   theme_bw() +
#   labs(fill="") +
#   theme(legend.position="bottom") +
#   xlim(c(0, 1.01)) +
#   xlab("CoDa-R-squared")
# 
# 
# g_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
# 
# mylegend<-g_legend(p2)


#Join both r-squared
r_squared <- cbind(r_squared_res, r_sq = "r_sq_res") %>% rbind(., cbind(r_squared, r_sq = "r_sq_var"))
# r_squared$fit <- factor(r_squared_res$fit, labels = c(paste0("M", 1:m1)),
#                             levels = c(paste0("V", 1:m1)))
r_squared$r_sq <- factor(r_squared$r_sq, 
                               levels = c("r_sq_res", "r_sq_var"),
                               ordered = TRUE,
                               labels = c(expression(paste("BR-CoDa-", R^2)), expression(paste("BM-CoDa-", R^2))))



r_squared_mean_data <- rbind(cbind(r_squared_res_mean, r_sq = 1),
                               cbind(r_squared_mean, r_sq = 2))


r_squared_mean_data$r_sq <- factor(r_squared_mean_data$r_sq, 
                                     levels = c(1, 2),
                                     ordered = TRUE,
                                     labels = c(expression(paste("BR-CoDa-", R^2)), 
                                                expression(paste("BM-CoDa-", R^2))))



r_squared_mean_data
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

# pdf("figures/r-squared_sim_linear_regression.pdf", width = 10, height = 6)
# p1 
# dev.off()





### --- 6. Computing waic --- ####
# loo_res <- list_mod %>%
#   pbsapply(., function(x){
#     loo1 <- loo::loo(x)
#     loo1$estimates[3,1]
#   })
   

waic_res <- list_mod %>%
  pblapply(., function(x){
    waic1 <- waic(x)
    waic1
  })
waic_res1 <- waic_res %>% sapply(., function(x) x$estimates[3]) %>% round(., 3) 

### --- 6. Plotting posterior distributions --- ####
## Checking recovering of the posterior distributions
#We use the first model of test1

fit1 <- mod_test3_1
fit1

### ----- 6.1. Intercepts --- ####
data_plot <- fit1 %>%
  gather_draws(b_yilr1_Intercept, b_yilr2_Intercept, b_yilr3_Intercept) 
colnames(data_plot) <- c("chain", "iteration", "draw", "variable_p",  "value_p")

values_plot <- data.frame(values_int = betas[seq(1,length(betas), 4)], variable_p = c("b_yilr1_Intercept",  "b_yilr2_Intercept", "b_yilr3_Intercept"))

ggplot(data = data_plot) +
  geom_histogram(aes(x = value_p, y = ..density.., group = variable_p),
                 color = "#e9ecef", alpha=0.6, position = 'identity', bins = 30, fill = "cornflowerblue") +
  theme_bw() +
  geom_vline(data = values_plot, aes(xintercept = values_int), 
             color = "blue4", size=1) +
  facet_wrap( ~ variable_p)

fit1 <- mod_test3_1
fit1

### ----- 6.2. Fixed effects --- ####
data_plot <- fit1 %>%
  gather_draws(b_yilr1_x1, b_yilr1_x2, b_yilr1_x3, 
               b_yilr2_x1, b_yilr2_x2, b_yilr2_x3,
               b_yilr3_x1, b_yilr3_x2, b_yilr3_x3) 
colnames(data_plot) <- c("chain", "iteration", "draw", "variable_p",  "value_p")

values_plot <- data.frame(values_int = betas[-seq(1,12, 4)], variable_p = as.vector(outer(paste0("b_yilr", 1:3), paste0("_x", 1:3), paste, sep="")) %>% 
                            matrix(., ncol = 3) %>% t() %>% as.character())

ggplot(data = data_plot) +
  geom_histogram(aes(x = value_p, y = ..density.., group = variable_p),
                 color = "#e9ecef", alpha=0.6, position = 'identity', bins = 15, fill = "cornflowerblue") +
  theme_bw() +
  facet_wrap( ~ variable_p, scales="free") +
  geom_vline(data = values_plot, aes(xintercept = values_int),
             color = "blue4", size=1) 

### ----- 6.3. Correlation parameters --- ####

get_variables(fit1)
#correlation parameters
data_plot <- fit1 %>%
  gather_draws(sigma_yilr1, sigma_yilr2, sigma_yilr3, rescor__yilr1__yilr2,
               rescor__yilr1__yilr3, rescor__yilr2__yilr3)  
colnames(data_plot) <- c("chain", "iteration", "draw", "variable_p",  "value_p")

values_plot <- data.frame(values_int = c(sigma, rho), variable_p = c(paste0("sigma_yilr", 1:3), 
                                                                     "rescor__yilr1__yilr2", 
                                                                     "rescor__yilr1__yilr3",
                                                                   "rescor__yilr2__yilr3"))

ggplot(data = data_plot) +
  geom_histogram(aes(x = value_p, y = ..density.., group = variable_p),
                 color = "#e9ecef", alpha=0.6, position = 'identity', bins = 15, fill = "cornflowerblue") +
  theme_bw() +
  facet_wrap( ~ variable_p, scales="free") +
  geom_vline(data = values_plot, aes(xintercept = values_int),
             color = "blue4", size=1) 



### --- 7. Posterior distributions in tables --- ####

### ----- 7.1. Fixed effects --- ####
library(xtable)
table_1 <- summary(mod_test3_1)
table_1$fixed <- as.data.frame(table_1$fixed)

table_1_fixed <- rbind(table_1$fixed[stringr::str_detect(rownames(table_1$fixed), "^yilr1"),],
      table_1$fixed[stringr::str_detect(rownames(table_1$fixed), "^yilr2"),],
      table_1$fixed[stringr::str_detect(rownames(table_1$fixed), "^yilr3"),])
      
table_1_fixed <- cbind(real = betas, table_1_fixed)
xtable::xtable(table_1_fixed[, c("real", "Estimate", "l-95% CI", "u-95% CI")], digits = 3)


### ----- 7.2. Hyperparameters --- ####
table_1$spec_pars <- as.data.frame(table_1$spec_pars)
table_1$rescor_pars <- as.data.frame(table_1$rescor_pars)

table_1_hyper <- cbind(real = c(sigma, rho), rbind(table_1$spec_pars, table_1$rescor_pars))
xtable::xtable(table_1_hyper[, c("real", "Estimate", "l-95% CI", "u-95% CI")], digits = 3)


### --- 8. Plotting fitted against data --- ####
# Data_inverse
z_comp <- ilrInv(data_test3[, c("yilr1", "yilr2", "yilr3")])
z_comp <- z_comp %>% as.data.frame() 
colnames(z_comp) <- paste0("Cat", 1:4)  
z_comp %>%
  tidyr::pivot_longer(., cols = 1:4,  names_to = "Category", values_to = "prop_data") -> z_comp

# Prediction
y_pred <- predict(mod_test3_1)
z_comp_pred <- ilrInv(y_pred[,"Estimate",])
z_comp_pred <- z_comp_pred %>% as.data.frame() 
colnames(z_comp_pred) <- paste0("Cat", 1:4)  
z_comp_pred %>%
  tidyr::pivot_longer(., cols = 1:4,  names_to = "Category", values_to = "prop_pred") -> z_comp_pred

#Join
z_comp_total <- cbind(z_comp, prop_pred = z_comp_pred$prop_pred)

z_comp_total$Category <- factor(z_comp_total$Category, 
                         levels = c(paste0("Cat",1:4)),
                         ordered = TRUE,
                         labels = c(paste0("Category-", 1:4)))


# Plotting
# pdf("figures/sim_linear_regression_fitted.pdf", width = 8, height = 7)
# ggplot(z_comp_total, aes(x = prop_data, y = prop_pred)) +
#   geom_smooth(method = lm, se = TRUE, level = 0.99) +
#   geom_point() + 
#   scale_y_continuous(breaks= scales::pretty_breaks()) +
#   scale_x_continuous(breaks= scales::pretty_breaks()) +
#   facet_wrap(~Category, scales = "free", labeller = label_parsed) +
#   theme_bw() +
#   xlab("Data proportions") +
#   ylab("Predicted proportions")
# dev.off()

### --- 9. Plotting prior distributions --- ####
prior_summary(mod_test3_1)

### ----- 9.1. Standard deviation --- ####
# student_t(3, 0, 2.5)     sigma            yilr1                       default
# student_t(3, 0, 2.5)     sigma            yilr2                       default
# student_t(3, 0, 3.1)     sigma            yilr3                       default
x <- seq(-10, 10, length.out = 100)
y1 <- dstudent_t(x, df = 3, mu = 0, sigma = 2.5, log = FALSE)
plot(x,y1)

#lkj correlation cholesky


### --- 10. Table for models --- ####
r_squared_mean_data
#sapply(list_mod, function(x){x[[1]][[1]][[1]][[1]]}) %>% unlist(.)
total_models <- data.frame(Model = r_squared_mean_data$fit[1:6],
                           waic   = waic_res1, 
                           BR    = r_squared_mean_data[1:6, c("m")],
                           BM    = r_squared_mean_data[7:12, c("m")])

xtable::xtable(total_models, digits = 3)

### --- 11. Computing probabilities to compare --- ####
lapply(paste0("M", 1:6), function(x){ r_squared %>% 
  dplyr::filter(fit == x) %>%
  dplyr::filter(r_sq == 'paste("BR-CoDa-", R^2)') %>%
  dplyr::select(R_squared) %>% pull(.)
}) -> m_sim_br

lapply(paste0("M", 1:6), function(x){ r_squared %>% 
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



(m_sim_bm[[all_mod[1,1]]] > m_sim_bm[[all_mod[1,2]]]) %>% 
  table(.)/10000
hist(m_sim_bm[[all_mod[7,1]]])
hist(m_sim_bm[[all_mod[7,2]]])

(m_sim_bm[[all_mod[7,1]]] > m_sim_bm[[all_mod[7,2]]]) %>% 
  table(.)


