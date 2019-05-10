require(ggplot2)
require(directlabels)
library(reshape)
library(scales) # transparent plots



map_mod = function(x) c("Gaussian", "Logistic", "Poisson","Gamma")[x]

create_mod_factor = function(agg_df){
  agg_df[,1] = map_mod(agg_df[,1])
  agg_df$Mod<- factor(agg_df$Mod, levels =c("Gaussian", "Logistic", "Poisson","Gamma") )
  return(agg_df)
}




####### MSE Plot #######
create_MSE_plot=function(data_df, dim, colScale, shapeScale, exclude_gauss=FALSE){

  # aggregate data
  agg_mse_2d = aggregate( cbind(MSE_Laplace, MSE_VL, MSE_VL_z, MSE_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean )
  agg_mse_1d = aggregate( cbind(MSE_Laplace, MSE_VL, MSE_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean )
  colnames(agg_mse_1d)[5] = c("MSE_VL-IW")
  colnames(agg_mse_2d)[5:6] = c("MSE_VL-IW", "MSE_RF-full")

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
    p1 = ggplot(MSE_melted, aes(x = m, y = RMSE, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
      geom_line() + geom_point() +theme_bw() + theme(legend.position = "bottom")+ scale_shape(solid=FALSE)+
      facet_wrap(.~C_Smoothness, scales = "free_y",labeller = label_bquote(nu == .(C_Smoothness))) +colScale+shapeScale
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
  colnames(agg_ls_2d)[5:6] = c("LS_VL-IW", "LS_RF-full")

  agg_ls = agg_ls_1d
  if(dim==2) agg_ls = agg_ls_2d


  if(ncol(agg_mse)>6) agg_ls[,7] = (agg_ls[,7]-agg_ls[,4])
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


    ptit = expression("2D Ppoisson Samples, n=2500,"~nu == .5)
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
  colnames(agg_time)[c(2:5)] = c("n", "Laplace", "RF-full", "LowRank")

  agg_time = create_mod_factor(agg_time)
  time_melted = melt(agg_time, id = c("Mod", "n"))
  colnames(time_melted)[3:4] <- c("Algorithm", "Time")

  p1 = ggplot(time_melted, aes(x = n, y = Time, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
    geom_line() + geom_point() +theme_bw() + theme(legend.position = "top")+ scale_shape(solid=FALSE)+
    facet_wrap(.~ Mod, scales = "free_y", nrow=1) +colScale+shapeScale

  if(log_scale) p1 = p1+scale_y_log10()+scale_x_log10()

  p1


}



####### MSE with varying sample size Plot #######
create_MSE_vs_n_plot = function(data_df){

  colScale <- scale_colour_manual(name = "Algorithm",values = c(1,2, 3))
  shapeScale <- scale_shape_manual(name = "Algorithm",values = c(1,0,2))

  agg_tmse = aggregate( cbind(MSE_Laplace, MSE_VL_z, MSE_LowRank) ~ Mod+Sample,    data = data_df, mean )
  colnames(agg_tmse)[4] = c("MSE_RF-full")

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

#### Generate 1D Plots ####

### 1D plots
data_df = read.csv("../CPP_analysis/1D_nonpara_VLZY_fdg.csv")
colScale_1D <- scale_colour_manual(name = "Algorithm",values = c(1,2,3))
shapeScale_1D <- scale_shape_manual(name = "Algorithm",values = c(1,0,2))

# 1D simulation MSE plot
create_MSE_plot(data_df, dim=1, colScale=colScale_1D, shapeScale=shapeScale_1D )
#ggsave("1D_MSE_GG_nug.pdf", device= "pdf",width = 8, height = 5)

# 1D simulation LS plot
create_LS_plot(data_df, dim=1, colScale=colScale_1D, shapeScale=shapeScale_1D )
#ggsave("4D_LS_GG_nug.pdf", device= "pdf",width = 8, height = 5)




#### Generate 2D Plots ####

### 2D Simulation MSE + LS
data_df = read.csv("../CPP_analysis/2D_nonpara_VLZY_fdg.csv")
colScale_2D <- scale_colour_manual(name = "Algorithm",values = c(1,2,3,4))
shapeScale_2D <- scale_shape_manual(name = "Algorithm",values = c(NA,0,2,4))
create_MSE_plot(data_df, dim=2, colScale=colScale_2D, shapeScale=shapeScale_2D )
#ggsave("2D_MSE.pdf", device= "pdf",width = 8, height = 5)
create_LS_plot(data_df, dim=2, colScale=colScale_2D, shapeScale=shapeScale_2D )
#ggsave("2D_LS.pdf", device= "pdf",width = 8, height = 5)


# 2D Simulation time
data_df = read.csv("../CPP_analysis/time_nonparatest_VLZY.csv")
data_df2 = read.csv("../CPP_analysis/gauss_time.csv")
data_df =rbind(data_df, data_df2)
time_plot(data_df)
#ggsave("time_log_2D_vlzy_GG.pdf", device= "pdf", width = 8, height = 3)


#2D simulation MSE vs sample
data_df = read.csv("../CPP_analysis/2D_time_combined_large.csv")
create_MSE_vs_n_plot(data_df)
#ggsave("mse_vs_sample_2D_horiz_log.pdf", device= "pdf",width = 8, height = 3)


#2D Param estimation, pois
data_df = read.csv("../CPP_analysis/param_est_pois_2D_VLZY_LR_yinit_2.csv")
create_llh_contour_plot(data_df)


#Generate simulated grid search plot



#ggsave("LGCP_MSE_pois_v4.pdf", device= "pdf",width = 7, height = 4)
