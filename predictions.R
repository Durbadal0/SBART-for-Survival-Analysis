

Y.time=sapply(x.time,function(x){mean_survival_t_given_xij(5,x)})


plot(x.time,Y.time,type='l')




mean_hazard_t_given_x_county(t=.98,x=c(.99,33,5,8,0,98,0),ci=4)
fun_z=function(t){mean_hazard_t_given_x_county(t,x=c(.99,33,5,8,0,98,0),ci=4)}
y=(z[1]+z[2])/2
k=
fun_z(x.time[106])
Z[2]=y



z=unique(sort(runif(10,z[4],z[5])))
k=sapply(z,fun_z)
k



k=sapply(x.time,fun_z)



tree.p=extract(bart_save[[5]],"trees",sampleNums=1000)
barrt=extract(bart_save[[5]],"bart",sampleNums=1000)





getPredictionsForTree <- function(tree, x) {
  predictions <- rep(NA_real_, nrow(x))
  
  getPredictionsForTreeRecursive <- function(tree, indices) {
    if (tree$var[1] == -1) {
      # Assigns in the calling environment by using <<-
      predictions[indices] <<- tree$value[1]
      return(1)
    }
    
    goesLeft <- x[indices, tree$var[1]] <= tree$value[1]
    headOfLeftBranch <- tree[-1,]
    n_nodes.left <- getPredictionsForTreeRecursive(
      headOfLeftBranch, indices[goesLeft])
    
    headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
    n_nodes.right <- getPredictionsForTreeRecursive(
      headOfRightBranch, indices[!goesLeft])
    
    return(1 + n_nodes.left + n_nodes.right)
  }
  
  getPredictionsForTreeRecursive(tree, seq_len(nrow(x)))
  
  return(predictions)
}


getPredictionsForTree(tree[which(tree$tree==1),],matrix(c(0,0,0,
                             0,0,0,0,0),nrow = 1))






all_tree_bart=function(b_tree,testy){
  fun_z=function(i){
    getPredictionsForTree(b_tree[which(b_tree$tree==i),],testy)
  }
  res=matrix(NA,nrow=10,ncol = dim(testy)[1])
  for(j in 1:10){
    res[j,]=c(getPredictionsForTree(b_tree[which(b_tree$tree==j),],testy))
  }
  
  return(colSums(res))
}


all_tree_bart(tree.p,testy =matrix(rep(0,16),ncol=8))



bart_save[[5]]$fit$plotTree(treeNum = 10,sampleNum = 1000)

bart_res$fit[sampleNums=1000]

