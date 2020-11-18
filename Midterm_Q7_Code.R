# ----------- Title: Simulation Code for PoisBinOrdNonNor--------------# 
# Jun Lu #

# Library packages
require(PoisBinOrdNonNor)
require(moments)

# ----------- Functions to calculate estimates and CI--------------# 
# Function to compute percentage estimates for ordinal variable
ord.est = function(data, n.obs){
    per_table = data.frame(cumsum(table(data))/n.obs)
    per_table[1:3,1]
}

# Function to compute 95% CI
# Poisson 95% CI (variance stabilization)
poi.ci.f = function(est, n.obs){
    (c(-1.96, 1.96)*sqrt(1/(4*n.obs)) + sqrt(est))^2
}
# Binary 95% CI
bin.ci.f = function(est, n.obs){
    c(-1.96, 1.96)*sqrt(est*(1 - est)/n.obs) + est 
}
# Ordinal 95% CI (cumulative)
ord.ci.f = function(est, n.obs){
    p0 = est[1]
    p1 = est[2] - p0
    p2 = est[3] - est[2]
    t0.ci = c(-1.96, 1.96)*sqrt(p0*(1 - p0)/n.obs) + p0 
    t1.ci = c(-1.96, 1.96)*sqrt((p0*(1 - p0) + p1*(1 - p1) - 2*p0*p1)/n.obs) + p0 + p1
    t2.ci = c(-1.96, 1.96)*sqrt((p0*(1 - p0) + p1*(1 - p1) + 
                                     p2*(1 - p2) - 2*p0*p1 - 2*p1*p2 - 2*p0*p2)/n.obs) + p0 + p1 + p2
    list(t0.ci, t1.ci, t2.ci)
}
# Correlation 95% CI (Fisher Transformation)
cor.ci.f = function(rho, n.obs){
    x = 1/2*log((1 + rho)/(1 - rho)) + c(-1.96, 1.96)*sqrt(1/n.obs)
    (exp(2*x) - 1)/(exp(2*x) + 1)
}


# ----------- Generate random correlation matrix within boundaries--------------# 
gen.corr = function(seed, no.pois, no.bin, no.ord, no.nonn,
                    pois.list, bin.list, ord.list, 
                    is.ord.list.cum, nonn.list){
    set.seed(seed)
    
    # Calculates the approximate upper and lower correlation bounds
    re = lower.upper.cors(no.pois, no.bin, no.ord, no.nonn,
                          pois.list, bin.list, ord.list, 
                          is.ord.list.cum, nonn.list)
    L = re$min[upper.tri(re$min)]
    U = re$max[upper.tri(re$max)]
    
    # Generate one correlation matrix within bounds
    stat = F
    while (stat == F) {
        rho = numeric(length(L))
        for (i in 1:length(L)) {
            rho[i] = runif(1, L[i], U[i])
        }
        cor.mat = matrix(0, nrow = 8, ncol = 8)
        cor.mat[upper.tri(cor.mat)] = rho
        cor.mat = cor.mat + t(cor.mat) + diag(1,8)
        diag(cor.mat) = 1
        stat = is.positive.definite(cor.mat)
    }
    return(cor.mat)
}

# ----------- Create a simmulation function--------------# 
PoisBinOrdNonNor.sim = function(seed, n.obs, cor.mat, no.pois, no.bin, no.ord, 
                                no.nonn, pois.list, bin.list, 
                                ord.list, is.ord.list.cum, nonn.list) {
    set.seed(seed)
    
    TV = c(unlist(pois.list), 1 - unlist(bin.list), unlist(ord.list), cor.mat[upper.tri(cor.mat)])
    
    
    
    cmat.star = find.cor.mat.star(cor.mat, no.pois, no.bin, no.ord, 
                                  no.nonn, pois.list, bin.list, ord.list, 
                                  is.ord.list.cum, nonn.list)
    
    est_m = matrix(nrow = n.sim, ncol = 38)
    cover_m = matrix(nrow = n.sim, ncol = 38)
    for (i in 1:n.sim) {
        # Generate data
        data = genPBONN(n.obs, cmat.star, no.pois, no.bin, no.ord, 
                        no.nonn, pois.list, bin.list, 
                        ord.list, is.ord.list.cum, nonn.list)
        
        # Separate data into non-ordinal and ordinal
        data.ord = data[, 5:6]
        data.no.ord = data[, -(5:6)]
        
        # Caculate estimates and 95% CI for non-ordinal parameters
        est.non = apply(data.no.ord, 2, mean) # Poisson(lambda) Binary(p) Non-nomarl(mean)
        ci.poi = sapply(est.non[1:2], poi.ci.f, n.obs) # CI for Poisson
        ci.bin = sapply(est.non[3:4], bin.ci.f, n.obs) # CI for Binary
        
        # Calculate estimates and 95% CI for ordinal parameters 
        est.ord = sapply(c(list(data.ord[,1]), list(data.ord[,2])), ord.est, n.obs, simplify = F)
        # (p1, p2, p3)
        ci.ord = sapply(est.ord, ord.ci.f, 100, simplify = T) # CI for ordinal
        
        # Caculate estimates and 95% for all correlation
        est.cor = cor(data)
        est.cor.vec = est.cor[upper.tri(est.cor)] # only upper traingle
        ci.cor = sapply(est.cor.vec, cor.ci.f, n.obs)
        
        # Organize Results
        est_m[i, ] = c(est.non[1:4], unlist(est.ord), est.cor.vec)
        est.ci = c(as.vector(ci.poi), as.vector(ci.bin), unlist(ci.ord), as.vector(ci.cor))
        
        # Test whether CI cover TV 
        cover_vec = numeric(38)
        for (j in 1:38) {
            L = est.ci[2*j - 1]
            U = est.ci[2*j]
            if ((TV[j] >= L) && (TV[j] <= U)) {
                cover_vec[j] = 1
            }else cover_vec[j] = 0
        }
        cover_m[i, ] = cover_vec
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
    return(reuslt_table)
}

# haven't calculated skewness and kurtosis (No CR) add later


# ----------- Setting 1. --------------# 

# Input variables
seed = 123
set.seed(seed)
n.obs = 100
n.sim = 1000
no.pois = 2
no.bin = 2
no.ord = 2
no.nonn = 2
pois.list = list(1, 2)
bin.list = list(.3, .6) # !!! percentage for 0
ord.list = list(c(.1, .2, .3), c(.4, .5, .7))
is.ord.list.cum = T
nonn.list = list(c(-1, 1, 0, 1), c(0, 3, 0, 2))
cor.mat = diag(0.8, 8) + 0.2

# Validates The Target Correlation Matrix
check = validate.cor.mat(cor.mat, no.pois, no.bin, no.ord, 
                         no.nonn, pois.list, bin.list, 
                         ord.list, is.ord.list.cum, nonn.list)
# Simulation
PoisBinOrdNonNor.sim(seed, n.obs, cor.mat, no.pois, no.bin, no.ord, 
                                no.nonn, pois.list, bin.list, 
                                ord.list, is.ord.list.cum, nonn.list) 

# At least for this setting everything works well.

# Q1 package problem upper bounds is not symmetric (wrong code .make.sym.matrix)
# Q2 continuous variables (estimates for mean var skew kurtosis? or just v1 v2)
# Q3 nearPD will produce matrix where diagonal >1







