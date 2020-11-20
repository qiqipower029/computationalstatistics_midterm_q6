require(PoisBinNonNor)
require(moments)
require(boot)
library(tidyverse)


# ----------- Functions to calculate estimates and CI--------------# 


# Function to cumpute point estimates for mean, variance, skewness and kurtosis for non-nomarl distribution;
non.normal.est = function(data){
  c(mean(data), var(data), skewness(data), kurtosis(data)-3)
}

# Function to compute 95% CI
# 1.Poisson 95% CI (variance stabilization)
poi.ci.f = function(est, n.obs){
  (c(-1.96, 1.96)*sqrt(1/(4*n.obs)) + sqrt(est))^2
}
# 2.Binary 95% CI
bin.ci.f = function(est, n.obs){
  c(-1.96, 1.96)*sqrt(est*(1 - est)/n.obs) + est 
}
# 3. Mean Variance Skewness (Kurtosis-3) 95% CI (bootstrap percentile)
nonnor_f = function(data, indices){
  d = data[indices]
  c(mean(d), var(d), skewness(d), kurtosis(d)-3)
}

nonnor.ci.f = function(data){
  re = boot(data = data, statistic = nonnor_f, R = 1000)
  mean.ci = boot.ci(re, type = "perc", index = 1)$percent[4:5]
  var.ci = boot.ci(re, type = "perc", index = 2)$percent[4:5]
  skew.ci = boot.ci(re, type = "perc", index = 3)$percent[4:5]
  kur.ci = boot.ci(re, type = "perc", index = 4)$percent[4:5]
  c(mean.ci, var.ci, skew.ci, kur.ci)
}

# 4.Correlation 95% CI (Fisher Transformation)
cor.ci.f = function(rho, n.obs){
  x = 1/2*log((1 + rho)/(1 - rho)) + c(-1.96, 1.96)*sqrt(1/n.obs)
  (exp(2*x) - 1)/(exp(2*x) + 1)
}


# ----------- Generate random correlation matrix within boundaries--------------# 
gen.corr = function(para.list){
  no.pois = para.list[[1]]
  no.bin = para.list[[2]]
  no.nonn = para.list[[3]]
  pois.list = para.list[[4]] %>% unlist()
  bin.list = para.list[[5]] %>% unlist()
  nonn.list = para.list[[6]] 
  
  skew = c(nonn.list[[1]][3], nonn.list[[2]][3])
  kurt = c(nonn.list[[1]][4], nonn.list[[2]][4])
  # Compute the coefficients of Fleishman third order polynomials
  coef.mat = fleishman.coef(no.nonn, skewness.vec = skew, kurtosis.vec = kurt)
  
  # Calculates the approximate upper and lower correlation bounds
  re = correlation.limits(no.pois, no.bin, no.nonn,
                          pois.list, bin.list, coef.mat)
  L = re[lower.tri(re)] %>% as.vector()
  U = re[upper.tri(re)] %>% as.vector()
  
  # Generate one correlation matrix within bounds
  stat = F
  
  while (stat == F) {
    rho = numeric(length(L))
    for (i in 1:length(L)) {
      rho[i] = runif(1, max(L[i], -0.35), min(U[i],0.35))
    }
    cor.mat = matrix(0, nrow = 6, ncol = 6)
    cor.mat[upper.tri(cor.mat)] = rho
    cor.mat = cor.mat + t(cor.mat) + diag(1,6)
    diag(cor.mat) = 1
    stat = is.positive.definite(cor.mat)
  }
  return(cor.mat)
}

# ----------- Create a simmulation function--------------# 
PoisBinNonNor.sim = function(n.obs, n.sim, cor.mat, para.list) {
  no.pois = para.list[[1]]
  no.bin = para.list[[2]]
  no.nonn = para.list[[3]]
  pois.list = para.list[[4]] %>% unlist()
  bin.list = para.list[[5]] %>% unlist()
  nonn.list = para.list[[6]] 
  
  TV = c(unlist(pois.list), unlist(bin.list), 
         unlist(nonn.list), cor.mat[upper.tri(cor.mat)])
  
  skew = c(nonn.list[[1]][3], nonn.list[[2]][3])
  kurt = c(nonn.list[[1]][4], nonn.list[[2]][4])
  mean.vec = c(nonn.list[[1]][1], nonn.list[[2]][1])
  var.vec = c(nonn.list[[1]][2], nonn.list[[2]][2])
  # Compute the coefficients of Fleishman third order polynomials
  coef.mat = fleishman.coef(no.nonn, skewness.vec = skew, kurtosis.vec = kurt)
  
  cmat.star = overall.corr.mat(no.pois, no.bin, no.nonn, 
                               pois.list, bin.list, coef.mat, corr.mat = cor.mat)
  
  est_m = matrix(nrow = n.sim, ncol = 27)
  cover_m = matrix(nrow = n.sim, ncol = 27)
  for (i in 1:n.sim) {
    # Generate data
    print(i)
    start_time <- Sys.time()
    data = gen.PoisBinNonNor(n.obs, no.pois, no.bin, no.nonn, 
                             pois.list, bin.list, mean.vec = mean.vec, variance.vec = var.vec,
                             coef.mat = coef.mat, final.corr.mat = cmat.star)
    
    # Separate data into (pois, binary) (non-normal) and (ordinal)
    data.poibi = data[, 1:4]
    data.nonor = data[, 5:6]
    
    # Caculate estimates and 95% CI for non-ordinal parameters
    est.poibi = apply(data.poibi, 2, mean) # Poisson(lambda) Binary(p) Non-nomarl(mean)
    ci.poi = sapply(est.poibi[1:2], poi.ci.f, n.obs) # CI for Poisson
    ci.bin = sapply(est.poibi[3:4], bin.ci.f, n.obs) # CI for Binary

    
    # Calculate estimate and 95% CI for non-normal paramters
    est.nonor = sapply(c(list(data.nonor[,1]), list(data.nonor[,2])), 
                       non.normal.est, simplify = F)
    ci.nonnor = sapply(c(list(data.nonor[,1]), list(data.nonor[,2])), 
                       nonnor.ci.f, simplify = F)
    
    # Caculate estimates and 95% for all correlation
    est.cor = cor(data)
    est.cor.vec = est.cor[upper.tri(est.cor)] # only upper traingle
    ci.cor = sapply(est.cor.vec, cor.ci.f, n.obs)
    
    # Organize Results
    est_m[i, ] = c(est.poibi, unlist(est.nonor), est.cor.vec)
    est.ci = c(as.vector(ci.poi), as.vector(ci.bin), unlist(ci.nonnor), as.vector(ci.cor))
    
    # Test whether CI cover TV 
    cover_vec = numeric(27)
    for (j in 1:27) {
      L = est.ci[2*j - 1]
      U = est.ci[2*j]
      if ((TV[j] >= L) && (TV[j] <= U)) {
        cover_vec[j] = 1
      }else cover_vec[j] = 0
    }
    cover_m[i, ] = cover_vec
    end_time <- Sys.time()
    t = end_time - start_time
    print(t)
  }
  
  # Summarize Result
  AE = apply(est_m, 2, mean)
  sd = apply(est_m, 2, sd)
  RB = (AE - TV)*100/TV
  SB = abs(AE - TV)/sd*100
  RMSE = sqrt(apply((est_m - TV)^2, 2, mean))
  CR = apply(cover_m,2,mean)
  
  reuslt_table = data.frame(
    TV = TV,
    AE = round(AE, 4),
    RB = round(RB, 4),
    SB = round(SB, 4),
    RMSE = round(RMSE, 4),
    CR = round(CR, 4)
  )
  
  name_dist = c("Poi1", "Poi2", "Bin1", "Bin2", 
               "Non-nor1", "Non-nor2")
  names_cor = NULL
  for (i in 1:5) {
    for (j in (i+1):6) {
      names_cor = c(names_cor, paste(name_dist[i], name_dist[j]))
    }
  }
  names = c("lambda1", "lambda2", "p1", "p2", "mean1","var1", "skew1", "kur1",
            "mean2", "var2", "skew2", "kur2", names_cor)
  row.names(reuslt_table) = names
  
  return(reuslt_table)
}




# ----------- Setting 1. CorA --------------# 
# sample 100
# sample 1000
# sample 10000
seed = 1
set.seed(seed)
para.list.1 = list(no.pois  = 2, no.bin = 2, 
                   no.nonn = 2, pois.list = list(0.4, 4), 
                   bin.list = list(0.2, 0.5), 
                   nonn.list = list(c(0.6667, 0.0317, -0.4677, -0.3750),
                                    c(32, 64, 0.5, 0.375)))

cor.mat.1 = gen.corr(para.list.1)


# Simulation
re1 = PoisBinNonNor.sim(100, 1000, cor.mat.1, para.list.1)
write.csv(re1, "result_s1.csv")

re2 = PoisBinNonNor.sim(1000, 1000, cor.mat.1, para.list.1)
write.csv(re2, "result_s2.csv")

re3 = PoisBinNonNor.sim(10000, 1000, cor.mat.1, para.list.1)
write.csv(re3, "result_s3.csv")



# ----------- Setting 2. Cor B --------------# 
# sample 100
# sample 1000
# sample 10000

# Simulation
set.seed(100)
cor.mat.2 = gen.corr(para.list.1)

re4 = PoisBinNonNor.sim(100, 1000, cor.mat.2, para.list.1)
write.csv(re4, "result_s4.csv")

re5 = PoisBinNonNor.sim(1000, 1000, cor.mat.2, para.list.1)
write.csv(re5, "result_s5.csv")

re6 = PoisBinNonNor.sim(10000, 1000, cor.mat.2, para.list.1)
write.csv(re6, "result_s6.csv")

# ----------- Setting 3. Cor C --------------# 
set.seed(2)

para.list.2 = list(no.pois  = 2, no.bin = 2, 
                   no.nonn = 2, pois.list = list(0.1, 8), 
                   bin.list = list(0.1, 0.7), 
                   nonn.list = list(c(0, 3, 1, 2),
                                    c(2, 2, 0, -.9582)))

cor.mat.1 = gen.corr(para.list.2)

re7 = PoisBinNonNor.sim(100, 1000, cor.mat.1, para.list.2)
write.csv(re7, "result_s7.csv")

re8 = PoisBinNonNor.sim(1000, 1000, cor.mat.1, para.list.2)
write.csv(re8, "result_s8.csv")

re9 = PoisBinNonNor.sim(10000, 1000, cor.mat.1, para.list.2)
write.csv(re9, "result_s9.csv")

# ----------- Setting 4. Cor D --------------# 

set.seed(3)
cor.mat.2 = gen.corr(para.list.2)

re10 = PoisBinNonNor.sim(100, 1000, cor.mat.2, para.list.2)
write.csv(re10, "result_s10.csv")

re11 = PoisBinNonNor.sim(1000, 1000, cor.mat.2, para.list.2)
write.csv(re11, "result_s11.csv")

re12 = PoisBinNonNor.sim(10000, 1000, cor.mat.2, para.list.2)
write.csv(re12, "result_s12.csv")