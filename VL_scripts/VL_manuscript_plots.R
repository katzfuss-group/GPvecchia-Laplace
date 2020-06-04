require(ggplot2)
require(directlabels)
require(reshape)
require(scales) # transparent plots

#######################################################################################
#####  Contains all plots that do not involve HMC/MCMC, LGCP, or satellite data  ######
#######################################################################################

map_mod = function(x) c("Gaussian", "Logistic", "Poisson","Gamma")[x]

create_mod_factor = function(agg_df){
  agg_df[,1] = map_mod(agg_df[,1])
  agg_df$Mod<- factor(agg_df$Mod, levels =c("Gaussian", "Logistic", "Poisson","Gamma") )
  return(agg_df)
}




###### Simple Example Plot #####
create_simple_example = function(){
  # library(GPVecchia) or source("../GPvecchia/VL_scripts/importer.R")
  source(file.path("VL_scripts/data_generating_functions.R"))

  pdf("VL_scripts/plots/simple_example.pdf", width = 7, height = 7)
  par(mfrow=c(1,1))
  covparms= c(1,.2,1)
  covfun <- function(locs) Matern(rdist(locs), range = covparms[2], smoothness = covparms[3])
  ls=logistic_sample(100, covfun, seed = 13)

  m=3
  vecchia.approx=vecchia_specify(ls$locs, m)#, cond.yz = "z"  )
  vecchia.approx_LL=vecchia_specify(ls$locs, m, conditioning = "firstm")#, cond.yz = "z"  )
  # for quick evaluation, compare posteriors
  ll_pred = calculate_posterior_VL(ls$z, vecchia.approx_LL, ls$type, covparms, return_all = T)
  vl_pred = calculate_posterior_VL(ls$z, vecchia.approx, ls$type, covparms, return_all = T)
  lap_pred = calculate_posterior_laplace(ls$z, ls$type, covfun(ls$locs), return_all = TRUE)

  plot(ls$locs, ls$y, type="p", ylim=c(-3, 2), ylab = "y", xlab = "s",main = "Sample Realization and Posterior Estimates")
  points(ls$locs,vl_pred$mean, col=3, type="l")
  points(ls$locs,vl_pred$mean- 1.96*sqrt(diag(solve(vl_pred$W))), col=3, type="l", lty = 3)
  points(ls$locs,vl_pred$mean+ 1.96*sqrt(diag(solve(vl_pred$W))), col=3, type="l", lty = 3)

  points(ls$locs,ll_pred$mean, col=2, type="l")
  #points(ls$locs, ls$score(ll_pred$mean,ls$z), type="l", col=4)
  points(ls$locs,ll_pred$mean+ 1.96*sqrt(diag(solve(ll_pred$W))), col=2, type="l", lty=3)
  points(ls$locs,ll_pred$mean- 1.96*sqrt(diag(solve(ll_pred$W))), col=2, type="l", lty = 3)

  points(ls$locs,lap_pred$mean, col=1, type="l")
  points(ls$locs,lap_pred$mean - 1.96*lap_pred$sd, col=1, type="l", lty = 3)
  points(ls$locs,lap_pred$mean + 1.96*lap_pred$sd, col=1, type="l", lty = 3)
  legend("bottomleft", legend = c("Truth", "Laplace", "VL (m=3)", "LowRank (m=3)"), lty = c(NA, 1, 1, 1), pch = c(1, NA, NA, NA), col=c(1,1,3,2))
  dev.off()
}


####### MSE Plot #######
create_MSE_plot=function(data_df, dim, colScale, shapeScale, exclude_gauss=FALSE, exclude_iw = FALSE){

  # aggregate data
  agg_mse_2d = aggregate( cbind(MSE_Laplace, MSE_VL, MSE_VL_z, MSE_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean )
  agg_mse_1d = aggregate( cbind(MSE_Laplace, MSE_VL, MSE_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean )
  colnames(agg_mse_1d)[5] = c("MSE_VL-IW")
  colnames(agg_mse_2d)[5:6] = c("MSE_VL-IW", "MSE_VL-RF")
  if(exclude_iw) agg_mse_2d = agg_mse_2d[,-5]
  agg_mse = agg_mse_1d
  if(dim==2) agg_mse = agg_mse_2d

  if(ncol(agg_mse)>6)  agg_mse[,7] = sqrt(agg_mse[,7]/agg_mse[,4])
  agg_mse[,5] = sqrt(agg_mse[,5]/agg_mse[,4])
  agg_mse[,6] = sqrt(agg_mse[,6]/agg_mse[,4])
  agg_mse[,4] = agg_mse[,4]/agg_mse[,4]


  agg_mse = create_mod_factor(agg_mse)
  MSE_melted = melt(agg_mse, id = c("Mod", "Neighbors", "C_Smoothness"))
  colnames(MSE_melted)[4:5] <- c("Algorithm", "RRMSE")
  colnames(MSE_melted)[2] <- c("m")
  MSE_melted$Algorithm <- substr(MSE_melted$Algorithm,5, 99)

  p1 = ggplot(MSE_melted, aes(x = m, y = RRMSE, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
    geom_line() + geom_point() +theme_bw() + theme(legend.position = "top")+ scale_shape(solid=FALSE)+
    facet_grid(C_Smoothness~ Mod, scales = "free_y",labeller = label_bquote(nu == .(C_Smoothness))) +colScale+shapeScale

  if(exclude_gauss){
    MSE_melted = MSE_melted[which(MSE_melted$Mod!="Gaussian"),]
    p1 = ggplot(MSE_melted, aes(x = m, y = RRMSE, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
      geom_line() + geom_point() +theme_bw() + theme(legend.position = "top")+ scale_shape(solid=FALSE)+
      facet_grid(C_Smoothness~ Mod, scales = "free_y",labeller = label_bquote(nu == .(C_Smoothness))) +colScale+shapeScale
    #facet_wrap(.~C_Smoothness, scales = "free_y",labeller = label_bquote(nu == .(C_Smoothness))) +colScale+shapeScale
  }

  p1
}

####### LS Plot #######
create_LS_plot = function(data_df, dim, colScale, shapeScale, exclude_gauss=FALSE){

  #  data_df$d1 = data_df$LS_VL - data_df$LS_Laplace
  # data_df$d1 = data_df$LS_LowRank - data_df$LS_Laplace

  # aggregate and negate logscores (for consistence appearance wrt MSE)
  agg_ls_1d = aggregate( cbind(LS_Laplace, LS_VL, LS_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean )
  agg_ls_1d[,4:6] = -agg_ls_1d[,4:6]
  agg_ls_2d = aggregate( cbind(LS_Laplace, LS_VL, LS_VL_z, LS_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean )
  agg_ls_2d[,4:7] = -agg_ls_2d[,4:7]

  colnames(agg_ls_1d)[5] ="LS_VL-IW"
  colnames(agg_ls_2d)[5:6] = c("LS_VL-IW", "LS_VL-RF")
  agg_ls_2d = agg_ls_2d[,-5]

  agg_ls = agg_ls_1d
  if(dim==2) agg_ls = agg_ls_2d


  if(ncol(agg_ls)>6) agg_ls[,7] = (agg_ls[,7]-agg_ls[,4])
  agg_ls[,6] = (agg_ls[,6]-agg_ls[,4])
  agg_ls[,5] = (agg_ls[,5]-agg_ls[,4])
  agg_ls[,4] = agg_ls[,4]-agg_ls[,4]


  agg_ls = create_mod_factor(agg_ls)
  LS_melted = melt(agg_ls, id = c("Mod", "Neighbors", "C_Smoothness"))
  colnames(LS_melted)[4:5] <- c("Algorithm", "dLS")
  colnames(LS_melted)[2] <- c("m")
  LS_melted$Algorithm <- substr(LS_melted$Algorithm,4,99)


  p1 = ggplot(LS_melted, aes(x = m, y = dLS, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
    geom_line() + geom_point() +theme_bw() + theme(legend.position = "bottom")+ scale_shape(solid=FALSE)+
    facet_grid(C_Smoothness~Mod, scales = "free_y", labeller = label_bquote(nu == .(C_Smoothness))) +colScale+shapeScale

  if(exclude_gauss){
    colnames(LS_melted)[2] <- c("m")
    LS_melted = LS_melted[which(LS_melted$Mod!="Gaussian"),]
    LS_melted = LS_melted[which(LS_melted$Algorithm!="VL_I"),]
    ggplot(LS_melted, aes(x = m, y = LS, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
      geom_line() + geom_point() +theme_bw() + theme(legend.position = "bottom")+ scale_shape(solid=FALSE)+
      facet_wrap(.~C_Smoothness, scales = "free_y", labeller = label_bquote(nu == .(C_Smoothness))) +colScale+shapeScale
    #ggsave("LGCP_LS_pois_v4.pdf", device= "pdf",width = 7, height = 4)

    pois_subset = LS_melted[which(LS_melted$Mod=="Poisson" & LS_melted$C_Smoothness==.5) ,]
    pois_subset[which(pois_subset$Algorithm=="VL_ZY"),][,4] = "LOCA"
    pois_subset = pois_subset[which(pois_subset$Algorithm!="VL_I"),]


    ptit = expression("2D Poisson Samples, n=2500,"~nu == .5)
    ggplot(pois_subset, aes(x = Neighbors, y = LS, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
      geom_line() + geom_point() +theme_bw() + theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5))+
      scale_shape(solid=FALSE)+colScale
    #ggsave("2D_poisson_alone.pdf", device= "pdf",width = 4, height = 4)
  }
  p1

}



####### Time Complexity Plot #######

create_time_plot = function(data_df, log_scale= FALSE){

  colScale <- scale_colour_manual(name = "Algorithm",values = c(1,3,2))
  shapeScale <- scale_shape_manual(name = "Algorithm",values = c(1,2,0))
  agg_time = aggregate( cbind(Time_Laplace, Time_VL_z, Time_LowRank) ~ Mod+Sample,    data = data_df, mean )
  colnames(agg_time)[c(2:5)] = c("n", "Laplace", "VL", "LowRank")

  agg_time = create_mod_factor(agg_time)
  time_melted = melt(agg_time, id = c("Mod", "n"))
  colnames(time_melted)[3:4] <- c("Algorithm", "Time")

  p1 = ggplot(time_melted, aes(x = n, y = Time, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
    geom_line() + geom_point() +theme_bw() + theme(legend.position = "top")+ scale_shape(solid=FALSE)+
    facet_wrap(.~ Mod, scales = "free_y", nrow=1) +colScale+shapeScale+ylab("Time (s)")

  if(log_scale) p1 = p1+scale_y_log10()+scale_x_log10()

  p1


}



####### MSE with varying sample size Plot #######
create_MSE_vs_n_plot = function(data_df){

  colScale <- scale_colour_manual(name = "Algorithm",values = c(1,2, 3))
  shapeScale <- scale_shape_manual(name = "Algorithm",values = c(1,0,2))

  agg_tmse = aggregate( cbind(MSE_Laplace, MSE_VL_z, MSE_LowRank) ~ Mod+Sample,    data = data_df, mean )
  colnames(agg_tmse)[4] = c("MSE_VL-RF")

  agg_tmse = create_mod_factor(agg_tmse)
  MSE_melted = melt(agg_tmse, id = c("Mod", "Sample"))
  colnames(MSE_melted)[3:4] <- c("Algorithm", "RMSE")
  colnames(MSE_melted)[2] <- c("n")
  MSE_melted$Algorithm <- substr(MSE_melted$Algorithm,5, 99)
  MSE_melted$RMSE[which(MSE_melted$RMSE>.9)]<-NA
  #MSE_melted<- MSE_melted[which(MSE_melted$n<160000),]
  p1 = ggplot(MSE_melted, aes(x = n, y = RMSE, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
    geom_line() + geom_point() +theme_bw() + theme(legend.position = "top")+ scale_shape(solid=FALSE)+
    facet_wrap(.~ Mod, scales = "free_y", nrow = 1) +colScale+shapeScale+scale_x_log10()

  p1
}




####### Parameter Estimation Contour Plot  #######
create_llh_contour_plot = function(data_df){


  sub_df = data_df[,c(1,2,3,7,10)]
  colnames(sub_df)[4:5]<- c("VL", "LowRank")
  param_est_melted = melt(sub_df, id = c("Smooth", "Range"))
  colnames(param_est_melted) <- c("Smoothness", "Range","Algorithm", "LLH")
  true_pt = data.frame(list("Smoothness" = rep(.5,3),"Range"=rep(.05,3),"Algorithm" = unique(param_est_melted$Algorithm), "LLH"=rep(0,3)))

  #custom_breaks = c(-650, -900,-2000,-5000,-10000, -15000)
  custom_breaks = c(-935, -955, -975, -1000, -1035, -1070, -1090)
  cplt = ggplot(param_est_melted, aes(Smoothness, Range, z =LLH)) + #geom_raster(aes(fill=LLH)) +
    geom_contour(color = "black", breaks = custom_breaks)+
    facet_grid(.~Algorithm)+theme_bw() + theme(legend.position = "none",text = element_text(size=20))
  cplt+ geom_dl(aes(label=..level..), method=list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust = .7, vjust = -.5, fill = "transparent"),
                stat="contour", breaks = custom_breaks)+geom_point(data=true_pt, color="red")

}


###### LGCP plots #######
# involved, see separate file

#### Generate 1D Plots ####

### 1D plots
data_df = read.csv("VL_scripts/saved_data/1D_nonpara_VLZY_fdg.csv")
colScale_1D <- scale_colour_manual(name = "Algorithm",values = c(1,2,3))
shapeScale_1D <- scale_shape_manual(name = "Algorithm",values = c(1,0,2))

# 1D simulation MSE plot
create_MSE_plot(data_df, dim=1, colScale=colScale_1D, shapeScale=shapeScale_1D )
#ggsave("VL_scripts/plots/1D_MSE_GG_nug.pdf", device= "pdf",width = 8, height = 5)

# 1D simulation LS plot
create_LS_plot(data_df, dim=1, colScale=colScale_1D, shapeScale=shapeScale_1D )
#ggsave("VL_scripts/plots/1D_LS_GG_nug.pdf", device= "pdf",width = 8, height = 5)




#### Generate 2D Plots ####

### 2D Simulation MSE + LS
data_df = read.csv("VL_scripts/saved_data/2D_nonpara_VLZY_fdg.csv")
colScale_2D <- scale_colour_manual(name = "Algorithm",values = c(1,2,4,3))
shapeScale_2D <- scale_shape_manual(name = "Algorithm",values = c(NA,0,4,2))
create_MSE_plot(data_df, dim=2, colScale=colScale_2D, shapeScale=shapeScale_2D, exclude_iw = TRUE, exclude_gauss = TRUE )
#ggsave("VL_scripts/plots/2D_MSE_pres.pdf", device= "pdf",width = 6, height = 3.5)
create_LS_plot(data_df, dim=2, colScale=colScale_2D, shapeScale=shapeScale_2D )
#ggsave("VL_scripts/plots/2D_LS.pdf", device= "pdf",width = 8, height = 5)


# 2D Simulation time
data_df = read.csv("VL_scripts/saved_data/time_nonparatest_VLZY.csv")
data_df2 = read.csv("VL_scripts/saved_data/gauss_time.csv")
data_df =rbind(data_df, data_df2)
create_time_plot(data_df)
#ggsave("VL_scripts/plots/time_log_2D_vlzy_GG.pdf", device= "pdf", width = 8, height = 3)


#2D simulation MSE vs sample
data_df = read.csv("VL_scripts/saved_data/2D_time_combined_large.csv")
create_MSE_vs_n_plot(data_df)
#ggsave("VL_scripts/plots/mse_vs_sample_2D_horiz_log.pdf", device= "pdf",width = 8, height = 3)


#2D Param estimation, pois example
data_df = read.csv("VL_scripts/saved_data/param_est_pois_2D_VLZY_LR_yinit_2.csv")
create_llh_contour_plot(data_df)
#ggsave("VL_scripts/plots/2D_pois_contour.pdf", device= "pdf",width = 7, height = 4)


#### 3D, 4D plots ####
data_df = read.csv("VL_scripts/saved_data/3D_nonpara.csv")
create_MSE_plot(data_df, dim=2, colScale=colScale_2D, shapeScale=shapeScale_2D, exclude_gauss = TRUE )
#ggsave("VL_scripts/plots/3d_MSE.pdf", device= "pdf",width = 7, height = 4)
#ggsave("VL_scripts/plots/4d_MSE.pdf", device= "pdf",width = 7, height = 4)



##### Pseudodata Plots ####
# library(GPVecchia) or source("../GPvecchia/VL_scripts/importer.R")
source(file.path("VL_scripts/data_generating_functions.R"))

dimen=1 # number of spatial dimensions
samp_size=30  # number of observed locs
domn = 1# domain over which to generate points

sig2=1; range=.05; smooth=1.5
covparms=c(sig2,range,smooth)
covfun <- function(locs) sig2*Matern(rdist(locs),range=range,smoothness=smooth)
m=4



pdf("VL_scripts/plots/pseudo_plots_paper.pdf", width = 10, height = 5)
# setup multiple plots
m_lay <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m_lay,heights = c(0.4,0.11))
# create logistic pseudo-data plot
ls = logistic_sample(samp_size, covfun, seed = 16, dom = domn, dimen=dimen)
vecchia.approx=vecchia_specify(ls$locs, m)#, cond.yz = "z"  )
pred1 = calculate_posterior_VL(ls$z, vecchia.approx, ls$type, covparms, return_all = TRUE)
pred = pred1
U = array(pred$t+.5*sqrt(pred$D))
L =array(pred$t-.5*sqrt(pred$D))
pred$sd = sqrt(diag(solve(pred$W)))
df = data.frame("loc" =ls$locs,"t" = array(pred$t), "U" = U, "L" = L,
                "z" = ls$z, "y" = pred$mean, "U_y" =(pred$mean+.5*pred$sd),
                "L_y" =(pred$mean-.5*pred$sd))
par(mar = c(2.5,3.1,2.,2.), mgp = c(1.5,.5,0))

plot(df$loc, df$t , main="Logistic",
     ylim=range(c(min(df$L), max(df$U))),
     pch=20, xlab="Location", ylab="")
arrows(df$loc, df$L, df$loc, df$U, length=0.01, angle=90, code=3, col=3)
points(df$loc, df$z, col=1)
points(df$loc, df$U_y, type ="l", lty = 2, col=4)
points(df$loc, df$L_y, type ="l", lty = 2, col=4)
points(df$loc, df$y, type = "l", lwd = 2, col=2)


# create Poisson pseudo-data plot
ls = pois_sample(samp_size, covfun, seed = 16, dom = domn, dimen=dimen)
vecchia.approx=vecchia_specify(ls$locs, m)#, cond.yz = "z"  )
pred2 = calculate_posterior_VL(ls$z,vecchia.approx, ls$type, covparms, return_all = TRUE)
pred = pred2
U = array(pred$t+.5*sqrt(pred$D))
L =array(pred$t-.5*sqrt(pred$D))
pred$sd = sqrt(diag(solve(pred$W)))
df = data.frame("loc" =ls$locs,"t" = array(pred$t), "U" = U, "L" = L,
                "z" = ls$z, "y" = pred$mean, "U_y" =(pred$mean+.5*pred$sd),
                "L_y" =(pred$mean-.5*pred$sd))

par(mar = c(2.5,2.1,2.,2.))

plot(df$loc, df$t, main="Poisson",
     ylim=range(c(min(df$L), max(df$z))),
     pch=20, xlab="Location", ylab="")
arrows(df$loc, df$L, df$loc, df$U, length=0.01, angle=90, code=3, col=3)
points(df$loc, df$z, col=1)
points(df$loc, df$U_y, type ="l", lty = 2, col=4)
points(df$loc, df$L_y, type ="l", lty = 2, col=4)
points(df$loc, df$y, type = "l", lwd = 2, col=2)
par(mar = c(2.1,2.1,2.,2.), mgp = c(0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom", legend = c("Data", "Pseudo data          ", "Pseudo Var  ", "Latent Post  ", "Post CI"),
       col=c(1,1,3,2,4), pch = c(1,20, NA, NA, NA), lty=c(NA,NA, 1, 1, 2), horiz = TRUE)
dev.off()



##### Parameter estimation variance table and plots ####
VL_df= read.csv("VL_scripts/saved_data/101_seed_param_est_VL.csv")
lap_df = read.csv("VL_scripts/saved_data/101_seed_param_est_lap.csv")
LR_df = read.csv("VL_scripts/saved_data/101_seed_param_est_LR.csv")
my_df = rbind(VL_df, LR_df)
# cast the data frame and keep complete cases
library(reshape2)
#  Version 1, by seed and m
check_cases = reshape(my_df, 
                      idvar = c("seed", "m"), 
                      timevar =c("method"), 
                      v.names = c("Range","Smoothness"), 
                      direction = "wide")

incomplete_case = merge(check_cases, unique(lap_df), by = "seed", all.x=F, all.y=F)
incomplete_case$m = incomplete_case$m.x


# relative to Truth: compute MSE
true_range = 0.05; true_smooth = 0.5
incomplete_case$Range.VL=(incomplete_case$Range.VL-true_range)^2
incomplete_case$Smoothness.VL=(incomplete_case$Smoothness.VL-true_smooth)^2
incomplete_case$Range.LR=(incomplete_case$Range.LR-true_range)^2
incomplete_case$Smoothness.LR=(incomplete_case$Smoothness.LR-true_smooth)^2
incomplete_case$Range =(incomplete_case$Range -true_range)^2
incomplete_case$Smoothness =  (incomplete_case$Smoothness-true_smooth)^2


require(dplyr)
## Aggregate MSE by method and m, take square root
MSE_cases = incomplete_case[,c(11,3,4,5,6, 9,10)]
fitable = MSE_cases %>% 
  na.omit() %>% 
  filter(.data$Range.LR<100) %>% 
  group_by(m) %>% 
  summarise(Range_VL_MSE = sqrt(mean(Range.VL)), 
            Range_VL_low = quantile(Range.VL, c(.05)), 
            Range_VL_hi = quantile(Range.VL, c(.95)),
            Smoothness_VL_MSE = sqrt(mean(Smoothness.VL)), 
            Smoothness_VL_low = quantile(Smoothness.VL, c(.05)), 
            Smoothness_VL_hi = quantile(Smoothness.VL, c(.95)),
            Range_LR_MSE = sqrt(mean(Range.LR)), 
            Range_LR_low = quantile(Range.LR, c(.05)),
            Range_LR_hi = quantile(Range.LR, c(.95)),
            Smoothness_LR_MSE = sqrt(mean(Smoothness.LR)), 
            Smoothness_LR_low = quantile(Smoothness.LR, c(.05)),
            Smoothness_LR_hi = quantile(Smoothness.LR, c(.95)),
            Range_Lap_MSE = sqrt(mean(Range)), 
            Range_Lap_low = quantile(Range, c(.05)), 
            Range_Lap_hi = quantile(Range, c(.95)),
            Smoothness_Lap_MSE = sqrt(mean(Smoothness)), 
            Smoothness_Lap_low = quantile(Smoothness, c(.05)),
            Smoothness_Lap_hi = quantile(Smoothness, c(.95)), 
            n = sum(!is.na(Range.VL))
)
fitable$n
# table 2
fitable[,c(1,2,5,8,11,14, 17,20)]

# Scatter plot for sample parameter estimates, m=20
scatter_df = my_df[my_df$m==20,]
scatter_m20 = rbind(lap_df, scatter_df)
require(ggplot2)
true_pt = data.frame(list("Range" = rep(.05,3),
                          "Smoothness"=rep(.5,3),
                          "method" = unique(scatter_dats$variable)))
ggplot(scatter_m20, aes(y = Range, x = Smoothness))+geom_point()+
  facet_grid(.~ method)+theme_bw()+coord_trans(x="log2", y="log2")+
  geom_point(data=true_pt[1:2], color="red", pch = 16, size = 2)+#geom_point(data=truncated_pt, pch=1)+
  scale_x_continuous(name="Smoothness", breaks = c(1, 5 ,10, 20, 40)) +
  scale_y_continuous(name="Range", breaks = c(1, 2,3,5)) 
ggsave("param_est_scatter_m20.pdf", device = "pdf", width = 8, height =3)


