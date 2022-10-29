## Load ----

library(tidyverse)
library(spBayesSurv)
library(LVBart)
library(gridExtra)
source("fit_icbart.R")

new_data <- readRDS("simulated_psa.rds")
fitted_new_psa <- readRDS("CacheFits/fitted_new_psa.rds")

## RMST Plots ----

mean_survs <- apply(fitted_new_psa$surv, c(2,3), mean)
rmst_survs <- apply(mean_survs, 1, function(x) cumsum(x) * 90 / 101)
grid_rmst  <- seq(0, 90, length = 101)

rmst_data  <- data.frame(time = grid_rmst, X = rmst_survs)

## Delta X1 ----

translate_1 <- function(x) {
  y <- "$X_3 = 1, X_4 = 1$"
  y <- ifelse(x == "Delta2", "$X_3 = 1, X_4 = 0$", y)
  y <- ifelse(x == "Delta3", "$X_3 = 0, X_4 = 1$", y)
  y <- ifelse(x == "Delta4", "$X_3 = 0, X_4 = 0$", y)
  return(y)
}

rmst_data_1 <- rmst_data %>%
  mutate(Delta1 = -X.1 + X.5, 
         Delta2 = -X.2 + X.6, 
         Delta3 = -X.3 + X.7, 
         Delta4 = -X.4 + X.8) %>%
  pivot_longer(cols = starts_with("D"), 
               values_to = "Delta", names_to = "Comp") %>%
  mutate(Covariate = translate_1(Comp))



p1 <- ggplot(rmst_data_1, aes(x = time, y = Delta, color = Covariate, linetype = Covariate)) + 
  geom_line(size = 1.4) + scale_color_viridis_d() + theme_bw() + 
  xlab("Time since surgery") + ylab("Difference in RMST") + 
  ggtitle("Difference for $X_2$") +
  theme(legend.position = "bottom")  + 
  guides(col=guide_legend(title = "", nrow = 4), linetype = guide_legend(title = "", nrow = 4))

## Delta X2 ----

translate_2 <- function(x) {
  y <- "$X_2 = 1, X_4 = 1$"
  y <- ifelse(x == "Delta2", "$X_2 = 1, X_4 = 0$", y)
  y <- ifelse(x == "Delta3", "$X_2 = 0, X_4 = 1$", y)
  y <- ifelse(x == "Delta4", "$X_2 = 0, X_4 = 0$", y)
  return(y)
}

rmst_data_2 <- rmst_data %>%
  mutate(Delta1 = -X.1 + X.3, 
         Delta2 = -X.2 + X.4, 
         Delta3 = -X.5 + X.7, 
         Delta4 = -X.6 + X.8) %>%
  pivot_longer(cols = starts_with("D"), values_to = "Delta", names_to = "Comp") %>%
  mutate(Covariate = translate_2(Comp))

p2 <- ggplot(rmst_data_2, aes(x = time, y = Delta, color = Covariate, linetype = Covariate)) + 
  geom_line(size = 1.4) + scale_color_viridis_d() + theme_bw()+ 
  xlab("Time since surgery") + ylab("Difference in RMST") + 
  ggtitle("Difference for $X_3$") +
  theme(legend.position = "bottom",
  legend.text = element_text(size=12))  + 
  guides(col=guide_legend(title = "", nrow = 4), linetype = guide_legend(title = "", nrow = 4))

## Delta X3 ----

translate_3 <- function(x) {
  y <- "$X_2 = 1, X_3 = 1$"
  y <- ifelse(x == "Delta2", "$X_2 = 1, X_3 = 0$", y)
  y <- ifelse(x == "Delta3", "$X_2 = 0, X_3 = 1$", y)
  y <- ifelse(x == "Delta4", "$X_2 = 0, X_3 = 0$", y)
  return(y)
}

rmst_data_3 <- rmst_data %>%
  mutate(Delta1 = -X.1 + X.2, 
         Delta2 = -X.3 + X.4, 
         Delta3 = -X.5 + X.6, 
         Delta4 = -X.7 + X.8) %>%
  pivot_longer(cols = starts_with("D"), values_to = "Delta", names_to = "Comp") %>%
  mutate(Covariate = translate_3(Comp))

p3 <- ggplot(rmst_data_3, aes(x = time, y = Delta, color = Covariate, linetype = Covariate)) + 
  geom_line(size = 1.4) + scale_color_viridis_d() + theme_bw() + 
  xlab("Time since surgery") + ylab("Difference in RMST") + 
  ggtitle("Difference for $X_4$") +
  theme(legend.position = "bottom") + 
  guides(col=guide_legend(title = "", nrow = 4), linetype = guide_legend(title = "", nrow = 4))


## Plot results ----

gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
