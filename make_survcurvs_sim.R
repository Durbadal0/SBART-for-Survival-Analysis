## Load ----

library(ProjectTemplate)
load.project()
new_data <- readRDS("CacheData/simulated_psa.rds") %>% arrange(cluster)

# Defining covariates for new subjects to be predicted ----

cov1 = c(1, 1, 1, 1, 0.2142857, 0.5517241) 
cov2 = c(1, 1, 1, 0, 0.2142857, 0.5517241) 
cov3 = c(1, 1, 0, 1, 0.2142857, 0.5517241) 
cov4 = c(1, 1, 0, 0, 0.2142857, 0.5517241) 
cov5 = c(1, 0, 1, 1, 0.2142857, 0.5517241) 
cov6 = c(1, 0, 1, 0, 0.2142857, 0.5517241) 
cov7 = c(1, 0, 0, 1, 0.2142857, 0.5517241) 
cov8 = c(1, 0, 0, 0, 0.2142857, 0.5517241) 
cov_new = rbind(cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8)

## Fit IC-BART ----

set.seed(digest::digest2int("fit the stuff"))

# fitted_new_psa <- fit_icbart(
#   X = as.matrix(new_data[,c("X.1", "X.2", "X.3", "X.4", "X.5", "X.6")]),
#   left = new_data$left,
#   right = new_data$right,
#   cluster = new_data$cluster,
#   status = new_data$status,
#   cov_new = cov_new,
#   num_burn = 4000,
#   num_thin = 1,
#   num_save = 1000,
#   do_loglik = TRUE
# )
# saveRDS(object = fitted_new_psa, file = "CacheFits/fitted_new_psa.rds")

fitted_new_psa <- readRDS("CacheFits/fitted_new_psa.rds")

## Make some plots ----

grid_size <- length(fitted_new_psa$surv[1,1,])
num_pred <- nrow(cov_new)

mean_surv <- do.call(
  c, 
  lapply(1:length(fitted_new_psa$surv[1,,1]), 
         FUN = function(x) colMeans(fitted_new_psa$surv[,x,])))

grid      = rep(seq(0,90, length=grid_size), num_pred)
covariate = rep(1:num_pred, each = grid_size)
psa.data  = data.frame(cbind(grid, mean_surv, covariate), row.names = NULL)

## Plot! ----

lcols=c("red","green","black","blue","pink","orange","yellow", "purple")
# linetype = rep(c('solid', 'dashed'),4)
linetype= rep("solid",num_pred)
LegendTitle = "Covariates"
cov_names= c("$X_2 = 1, X_3 = 1, X_4 = 1$",
             "$X_2 = 1, X_3 = 1, X_4 = 0$",
             "$X_2 = 1, X_3 = 0, X_4 = 1$",
             "$X_2 = 1, X_3 = 0, X_4 = 0$",
             "$X_2 = 0, X_3 = 1, X_4 = 1$",
             "$X_2 = 0, X_3 = 1, X_4 = 0$",
             "$X_2 = 0, X_3 = 0, X_4 = 1$",
             "$X_2 = 0, X_3 = 0, X_4 = 0$")
cohort = factor(rep(cov_names, each = grid_size))
psa.data1 <- psa.data %>% mutate(Covariates = cohort)
# psa.data1 = data.frame(cbind(psa.data, cohort))
g1=ggplot(data=psa.data1, mapping = aes(x=grid, y=mean_surv))
g1=g1+
  geom_line(aes(col=Covariates), size=1.5)+
  # geom_point(aes(col = Covariates, shape = Covariates)) + 
  ylim(0,1)+
  scale_color_manual(values = lcols) +
  # scale_color_viridis_d() + 
  # scale_linetype_manual(name = LegendTitle) +
  xlab("Time since surgery")+
  ylab("Survival Probability")+ggtitle("Estimated survival curves for time to recurrence of PSA")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + 
  theme(legend.position = "bottom", legend.text = element_text(size=12),
        legend.title = element_text(size=15) , legend.title.align = 0.5,
        axis.text=element_text(size=15),
        axis.title=element_text(size=14))+
  guides(col=guide_legend(title = LegendTitle, nrow = 4))
g1

