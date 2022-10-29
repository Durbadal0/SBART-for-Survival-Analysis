# posterior predictions

sigma1 <- read_csv("sigma1.csv")
sigma1=unlist(sigma1)


rho <- read_csv("rho.csv")
rho=unlist(rho)

var_imp <- read_csv("var_imp.csv")
frailty <- read_csv("frailty.csv")
alpha <- read_csv("alpha.csv")
frailty=frailty[,-1]
prediction_ij <- readRDS("C:/Users/durbadal.ghosh/OneDrive - Florida State University/Desktop/Survival Bart My functions/DURBART/prediction_ij.RDS")

cutoff <- readRDS("C:/Users/durbadal.ghosh/OneDrive - Florida State University/Desktop/Survival Bart My functions/DURBART/cutoff.RDS")

cum_hazard_ij <- readRDS("C:/Users/durbadal.ghosh/OneDrive - Florida State University/Desktop/Survival Bart My functions/DURBART/cum_hazard_ij.RDS")
bart_tree <- readRDS("C:/Users/durbadal.ghosh/OneDrive - Florida State University/Desktop/Survival Bart My functions/DURBART/bart_tree.RDS")

plot(mcmc(sigma1[1:1000]))
plot(mcmc(rho[1:1000]))

plot(mcmc(alpha[1:1000,2]))
plot(mcmc(frailty[1:1000,1:12]))



cum_hazard_ij[[1]]=rep(1,5113)


cu.haz=unlist(cum_hazard_ij)


str(cu.haz)


cum.haz=matrix(cu.haz,ncol=5113,byrow = T)


cum.haz_ij=colMeans(cum.haz)



m.residual_ij=cen-cum.haz_ij


Index=1:N
m.residual_ij
died=cen
dat_mres=data.frame(cbind(Index,m.residual_ij,died))

ggplot(data = dat_mres,aes(Index,m.residual_ij,color=(died)))+
  geom_point()
abline(0,0)
