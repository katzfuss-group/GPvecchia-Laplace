library(fields)


dimen=1 # number of spatial dimensions
samp_size=10^2  # number of observed locs
domn = 1# domain over which to generate points

sig2=1; range=.05; smooth=1.5
covparms=c(sig2,range,smooth)
covfun <- function(locs) sig2*Matern(rdist(locs),range=range,smoothness=smooth)
par(mfrow=c(1,1))
# current data generation is over a uniform lattice
ls = pois_sample(samp_size, covfun, seed = 16, dom = domn, dimen=dimen)
pred = estimate_laplace_vec_cpp(ls, covparms, m=4,return_all = TRUE)
ylim_D = c(min(pred$t-sqrt(diag(pred$D))), max(pred$t+sqrt(diag(pred$D))))


################  Scaling effect of likelihood on pseudodata  #############
#pdf("gamma_pseudodata__400_seed14.pdf")
plot(ls$locs, pred$t,  col=2, main = "Pseudodata")
points(ls$locs, ls$z, pch = ".")
points(ls$locs, pred$mean, type = "l", col=3)
legend("topright", legend = c("Pseudodata", "Observation", "Posterior"), col = c(2,1,3), pch=c("o",".","-"))
points(ls$locs, pred$mean+(pred$sd), col=3, type = "l", lty = 3)
points(ls$locs, pred$mean-(pred$sd), col=3, type = "l", lty = 3)

pseudo_offset = function(y, z)(ls$score(y_o = y, z = z)/diag(ls$hess(y_o = y, z = z)))
y = seq(-3, 3, length.out=100)
par(mfrow=c(1,1))
plot(y, y+pseudo_offset(y,1), col=1, main = "Scaling of Pseudodata",type = "l", ylab = "pseudo = y+pseudo offset at y|z")
points(y, y+pseudo_offset(y,1), col=2, type = "l")
points(y,y+pseudo_offset(y,3), col=3, type = "l")
points(y, y+pseudo_offset(y,5), col=4, type = "l")
points(y, y+pseudo_offset(y,40), col=5, type = "l")
legend("topright", legend = c("z=0", "z=1","z=3","z=5","z=7"), col = c(1,2,3,4,5), lwd=1)
#text(x = -2, y = -2,labels = "Note that slope = -1")
#legend("topright", legend = c("z=0", "z=1"), col = c(1,2), lty=1)
dev.off()



################  Scaling effect of likelihood on pseudodata  #############


gs = gauss_sample(samp_size, covfun, seed = 14, dom = domn, dimen=dimen)
pred_g = estimate_laplace_vec_cpp(ls, covparms, m=4,return_all = TRUE)
ls = logistic_sample(samp_size, covfun, seed = 14, dom = domn, dimen=dimen)
pred_l = estimate_laplace_vec_cpp(ls, covparms, m=4,return_all = TRUE)
ps = pois_sample(samp_size, covfun, seed = 14, dom = domn, dimen=dimen)
pred_p = estimate_laplace_vec_cpp(ls, covparms, m=4,return_all = TRUE)

sig_nois = function(y, z, ls){  
  sig_var =  (1/diag(ls$hess(y_o = y, z = z))) 
  return(sig_var)#/(1+sig_var))
}
y = seq(-3, 3, length.out=100)
plot(y, sig_nois(y,0,ls), col=1, main = "Signal",type = "l", ylab = "D")
plot(y, sig_nois(y,1,ls), col=1, main = "Signal",type = "l", ylab = "pseudo = y+pseudo offset at y|z")
sig_nois(c(0,1),1,ls)

plot(y, sig_nois(y,1,gs), col=1, main = "Signal",type = "l", ylab = "pseudo = y+pseudo offset at y|z")
plot(y, sig_nois(y,0,gs), col=1, main = "Signal" , type = "l", ylab =  "D/(D+1)")

plot(y, sig_nois(y,0,ps), col=1, main = "Signal" , type = "l", ylab =  "D/(D+1)")

