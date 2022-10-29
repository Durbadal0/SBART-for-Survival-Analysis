# data import and cleaning


library(readxl)
 Surv_rust <- read_excel("C:/Users/durbadal.ghosh/OneDrive - Florida State University/Desktop/Dr Rust/Surv_rust.xlsx")
 View(Surv_rust)
Surv_rust$county[which(Surv_rust$county==86)]=13
 
library(readr)
W_mat <- read_csv("C:/Users/durbadal.ghosh/OneDrive - Florida State University/Desktop/Dr Rust/W.mat.csv", 
                        col_types = cols(...1 = col_skip()))

 View(W_mat)
 
 
 Surv_rust_A=Surv_rust[which(Surv_rust$Race==2 & Surv_rust$stage_at_diagnosis=="L"),]
 
 Surv_rust_A[6,3]=2*34-1
 Surv_rust_A[1,3]=13*2-1
 Surv_rust_A[16,3]=21*2-1
 
 
 Surv_rust_A[573,3]=30*2-1

 
 
 age=Surv_rust_A$Age
 HR_p=Surv_rust_A$HR_p
 Tgrade=Surv_rust_A$Tgrade
 Smuprop=Surv_rust_A$SMU_Yes/(Surv_rust_A$SMU_Yes+Surv_rust_A$SMU_No)
 Adjuvant=Surv_rust_A$adjuvant
 Stagediag=as.numeric((as.factor(Surv_rust_A$stage_at_diagnosis)))
 Trtdelay=as.numeric((as.factor(Surv_rust_A$TX_Delay)))
 Biodelay=as.numeric((as.factor(Surv_rust_A$BX_Delay)))
 
 cen=Surv_rust_A$cen
 
 
 XDA=model.matrix(~age+HR_p+as.factor(Surv_rust_A$Tgrade)+Smuprop+Surv_rust_A$TX_Delay+Surv_rust_A$BX_Delay)
 X_Cov=XDA[,-c(1,7,9)]
 Yt=Surv_rust_A$as.numeric.date_diff.
 
 county=Surv_rust_A$county
county =(county+1)/2

Adj.M=W_mat

D=diag(rowSums(Adj.M),67)

 rou=.6
 s_W=(1/det(as.matrix(D-rou*Adj.M)))^(1/67)
 VCov_W=(as.matrix(s_W*(D-rou*Adj.M)))
 
 det(VCov_W)
 #start_time <- Sys.time() 

#end_time <- Sys.time()
#end_time-start_time


 f=function(i){
  qij=rpois(1,H(Yt[i])*county[i])
  Cij=runif(qij,0,H(Yt[i])*county[i])
  G_s_ij=H_inv(Yt[i]/county[i])
  Uij=runif(qij,0,1)
  return(list(Uij,G_s_ij))
}

 
