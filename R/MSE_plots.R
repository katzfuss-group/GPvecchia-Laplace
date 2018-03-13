

# header = c("Mod", "Domain", "Sample", "C_Smoothness", "C_Range","Neighbors", "Neighbors_LR","Seed_off", "MSE_Laplace", "MSE_VL",  "MSE_LowRank")
# mse_df = as.data.frame(mse_results)
# names(mse_df)=header
# write.csv(mse_df,file = "mse_results_mult_nbrs.csv")
setwd("/Users/danielzilber/Desktop/Katzfuss/code_CPP/")
library(reshape)
data_df = read.csv("1D_sim.csv" )

#  plot1:  mse vs neighbors, per model type
#header = c("Mod", "Domain", "Sample", "C_Smoothness", "C_Range","Neighbors","Seed_off", "MSE_Laplace", "MSE_VL",  "MSE_LowRank")
#mse_df = mse_df[,-1]
#names(mse_df)=header
head(data_df)
#average samples
agg_mse = aggregate( cbind(MSE_Laplace, MSE_VL, MSE_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean )
agg_ls = aggregate( cbind(LS_Laplace, LS_VL, LS_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean ); agg_ls[,4:6] = -agg_ls[,4:6]
agg_iters = aggregate( cbind(Iter_Laplace, Iter_VL, Iter_LowRank) ~ Sample+Neighbors+C_Smoothness,    data = data_df, mean )

agg_t = aggregate( cbind(Time_Laplace, Time_VL, Time_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean )
head(agg_mse)


sub_fun = function(mod, nbr) intersect(which(agg_mse$Mod == mod),which(agg_mse$Neighbors == nbr))
agg_mse[sub_fun(2,1),]

########################################################################
############### Plot set 1, effect of smoothing on MSE   #################
########################################################################
par(mfrow=c(2,2))
for (modl in c(1,2,3,4)){
  sub1 =agg_mse[sub_fun(modl,1),]
  plot(sub1$C_Smoothness, sub1$MSE_VL, ylim = c(.00, .4), type=  "l", main = paste("Model",modl))#, ylim = c(.04, .4))
  points(sub1$C_Smoothness,  sub1$MSE_LowRank, col=2, type = "l")
  points(sub1$C_Smoothness,  sub1$MSE_Laplace, col=3, type = "l")
  
  for( i in 2:6){
    sub1 =agg_mse[sub_fun(modl,i),]
    #points(sub1$C_Smoothness, sub1$MSE_Laplace, type = "l")
    points(sub1$C_Smoothness, sub1$MSE_VL, type = "l")
    points(sub1$C_Smoothness, sub1$MSE_LowRank, type = "l", col=2)
    points(sub1$C_Smoothness, sub1$MSE_Laplace, col=3, type = "l")
  }
}



########################################################################
############### Plot set 1, effect of neighbts on MSE   #################
########################################################################

agg_mse[sub_fun_sm(2,1),]

sm_vals = unique(agg_mse$C_Smoothness)


plot_by_nbr = function(agg_df, y_lab, dimen){
  sub_fun_mod = function(mod ) which(agg_df$Mod == mod)
  sub_fun_jsm = function(sm) which(agg_df$C_Smoothness == sm)
  sub_fun_sm = function(mod, smooth) intersect(which(agg_df$Mod == mod),which(agg_df$C_Smoothness ==smooth))
  
  par(mfrow=c(4,2))
  for (i in c(.5, 1.5)){
    
    sub_m =agg_df[sub_fun_jsm(i),]  
    
    sub1 =agg_df[sub_fun_sm(1,i),]
    ymax = max(sub1[,4:6])
    ymin = min(sub1[,4:6])
    
    # 5 = VL, 6 = LowRank, 4 = Laplace
    plot(sub1$Neighbors, sub1[,5], type=  "l", main = paste("Smoothness",i, "Gaussian D =",dimen), ylim = c(ymin, ymax), ylab = y_lab, xlab = "Conditioning Set Size")
    points(sub1$Neighbors, sub1[,5], type=  "p", pch=1)
    points(sub1$Neighbors,  sub1[,6], col=2, type = "l")
    points(sub1$Neighbors,  sub1[,6], col=2, type = "p", pch=2)
    points(sub1$Neighbors,  sub1[,4], col=3, type = "l")
    points(sub1$Neighbors,  sub1[,4], col=3, type = "p", pch=3)
    
    for(modl in c(2,3,4)){
      sub1 =agg_df[sub_fun_sm(modl,i),]
      ymax = max(sub1[,4:6])
      ymin = min(sub1[,4:6])
      mod_id = ": Gamma"
      if (modl==2) mod_id = ": Poisson"
      if (modl==3) mod_id = ": Logistic"
      plot(sub1$Neighbors, sub1[,5], type = "l",main = paste("Smoothness",i, mod_id, "D =", dimen), ylim = c(ymin, ymax), ylab = y_lab, xlab = "Conditioning Set Size")
      points(sub1$Neighbors, sub1[,5], type = "p", pch=1)
      points(sub1$Neighbors, sub1[,6], type = "l", col=2)
      points(sub1$Neighbors, sub1[,6], type = "p", col=2, pch=2)
      points(sub1$Neighbors, sub1[,4], col=3, type = "l")
      points(sub1$Neighbors, sub1[,4], col=3, type = "p", pch=3)
    }
  }
  
legend("topright", legend = c("VL", "Laplace", "Low Rank"), col=c(1,3,2), pch = c(1,3,2), lwd=1)
}



agg_df = agg_ls

plot_by_nbr(agg_mse, "MSE", 2)
plot_by_nbr(agg_t, "Time", 1)


pdf("LogScore_nbr_vs_smooth_2D_extended.pdf")
plot_by_nbr(agg_ls, "LS", 2)
dev.off()
#  Plot 3:  % increase of MSE vs neighbors



########################################################################
############### Plot set 2, effect of sample on time   #################
########################################################################
data_df = read.csv("results_time_vs_n.csv" )

agg_t = aggregate( cbind(Time_Laplace, Time_VL, Time_LowRank) ~ Mod+C_Smoothness+Sample,    data = data_df, mean )
update_funs(agg_t)
mod=1;smooth=.5

#pdf("Runtime_compare.pdf")
par(mfrow=c(1,1))
for(mod in 1:1){
agg_t_sub = agg_t[intersect(which(agg_t$Mod == mod),which(agg_t$C_Smoothness ==smooth)),]
x_vals = (agg_t_sub$Sample)
agg_t_sub[,4:6]=(agg_t_sub[,4:6] )
plot(x_vals, agg_t_sub[,4], type = "l", main = "Run Time, Logistic 1D", ylab = "Time (mins)", xlab = "Sample Size")
points(x_vals, agg_t_sub[,5], type = "l", col=3)
points(x_vals, agg_t_sub[,6], type = "l", col=2)
}
dev.off()

