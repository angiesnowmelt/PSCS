setwd("J:/PSCS/SIM2")
library(psych)
library(nlme)
library(lme4)
library(sampling)
library(plyr)


#~~~~~~~~~~~~~~~~~ Global variables ~~~~~~~~~~~~~~~#

#~~~ coefficients
b0 = 0
b1 = 0.8
b2 = -0.25
b3 = 0.6
b4 = -0.4
b5 = -0.8
b6 = -0.5
b7 = 0.7

b4_var = 0 # variance for the random slope; 
          # when b4 is not a constant b4 only, it has random slope

#~~ coefficients for the treatment effect formula
a0 = -3.85
a1 = 0.3
a2 = -0.36
a3 = -0.73
a4 = -0.2
a5 = -0.71
a6 = -0.19
a7 = 0.26
te = -0.4 # treatment effect

x1_var = 1
x2_var = 1
x3_var = 1
x4_var = 1
x5_var = 1 # on the same scale with x1
x6_var = 1 # on the same scale with x2
x7_var = 1
x8_var = 1 # on the same scale with x3
x9_var = 1 # on the same scale with x4
x10_var = 1

x1x5_cor = 0.2
x2x6_cor = 0.9
x1x5_cov = x1x5_cor*sqrt(x1_var)*sqrt(x5_var)
x2x6_cov = x2x6_cor*sqrt(x2_var)*sqrt(x6_var)
x3x8_cor = 0.2
x4x9_cor = 0.9
x3x8_cov = x3x8_cor*sqrt(x3_var)*sqrt(x8_var)
x4x9_cov = x4x9_cor*sqrt(x4_var)*sqrt(x9_var)


#~~~~~~~~~~~~~~ Specify grouping effects ~~~~~~~~~~~~~~~~~~~~#

# true ICC = 0.25; clustering and stratification each takes 0.125.
# true individual level variance 0.5 (arbitrary)

# Explained variance in Y through path tracing
# This is based on all X vars are continuous
# When generating Y_true, 1356 are dummy
# Thus, the explained variance is deflated(lower var and low cor between vars).
# Consequently, the empirical ICC is inflated, over 1/3.
# Since now we don't know what the yvar_ex is, to maintain ICC around 0.25, try larger yvar_ind
# so the calculation of all the other random effects are based on yvar_ind_prior. yvar_ind is flexible

yvar_ex = sum(b1*x1_var*b1, b2*x2_var*b2, b3*x3_var*b3, b4*x4_var*b4, 
    b5*x5_var*b5, b6*x6_var*b6, b7*x7_var*b7,
    2*b1*x1x5_cov*b5, 2*b2*x2x6_cov*b6) #2.5715
yvar_ind_prior = 0.5143 # predefined, assuming yvar_ind = yvar_cl = yvar_st

yvar_st = yvar_ind_prior # 0.5143, keep this value based on prior = 0.5143
yvar_cl = yvar_st
yvar_ind = 2 # set a random value to calculate the actual explained var, 
             # knowing yvar_ex, yvar_cl, and yvar_st, adjust this value to get ICC = 0.25

ysd_st = sqrt(yvar_st) # 0.7171471
ysd_cl = sqrt(yvar_cl)
ysd_ind = sqrt(yvar_ind)

# yvar_tot = yvar_ex + yvar_st + yvar_cl + yvar_ind # the tot var varies depending on yvar_ind

# Define the random effects for the treatment effect model
# ovar_ex = 


# icc = (yvar_st+yvar_cl)/yvar_tot


#~~~~~~~~~~ Function--Generate correlated random variable ~~~~~~~~~~~~~#

fun.sample.cor = function(x, rho) {
  y = (rho*(x-mean(x))/sqrt(var(x)) + sqrt(1-rho^2)*rnorm(length(x)))
  return(y)
}

#~~~~~~~~~~~~~~~~~~~ Generate Population ~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 7 predictors, X1, X2, X5, and X6 are binary
# Y_true-log odds of the treatment
# propensity scores are probabilities, and later transformed to binary
# correlations: (X1, X5) = 0.2, (X2, X6) = 0.9, (X3, X8) = 0.2, (X4, X9) = 0.9
# random slope: slope(X4), maybe define later

# 50 counties (K)
# 30 schools per county? (J)
# number of students per school (50, 11) (I)



K = 50 
J = 30 
num_col = 31 # consistent with the following columnnames  
columnnames = c("STUID", "SCHID",  "STRID", "School",  "Student", "n_j", "Nh","Y_true", "PS_true", "Treat", "Outcome",
  "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X11", "X33", "X55", "X66", "X8", "X9", "X10", "X88", "X99",
  "ELL", "St_eff", "Cl_eff", "Ind_eff")

#~~~ run these all together, index needs to be rerun each time

F.GENERATE = function(scenario, seed){


index =1 
N = 76000 # change this number when the final population size is decided
mydata = as.data.frame(matrix(NA, nrow = N, ncol = num_col+2)) # ncol+PRI+CLP(added later)


for (k in 1: K){
  # set.seed(999 + k)
  set.seed(seed + k)
  st_eff = ysd_st*rnorm(1) # for PS model
  st_eff2 = ysd_st*rnorm(1) # for TE model
  cl_size = round(rnorm(n = J, mean = 50, sd = 11), 0)
  st_size = matrix (NA, nrow = 1, ncol = K)
  st_size[, k] = sum(cl_size)
  Nh = sum(cl_size)
  mydata.h = as.data.frame(matrix(NA, nrow = Nh, ncol = num_col))
  index.h = 1 # within-strata index
  
  for (j in 1: J){
    cl_eff = ysd_cl*rnorm(1) # for PS model
    cl_eff2 = ysd_cl*rnorm(1) # for TE model
    cl_eff_slope = rnorm(1)*b4_var # much smaller variance
    ind_eff = ysd_ind*rnorm(cl_size[j]) # for PS model
    ind_eff2 = ysd_ind*rnorm(cl_size[j]) # for TE model
    
    #~~~ Stratify schools within each strata (PRI) later outside the cl loop based on cl_eff
    # select 1/3 (10) private and 2/3(20) public schools for each county
    
      X1 = rnorm(cl_size[j], 0, sqrt(x1_var)) 
      X2 = rnorm(cl_size[j], 0, sqrt(x2_var)) 
      X3 = rnorm(cl_size[j], 0, sqrt(x3_var))
      X4 = rnorm(cl_size[j], 0, sqrt(x4_var)) # ran slope is only rel. to the coeff not the var itself
      X5 = fun.sample.cor(X1, x1x5_cor)
      X6 = fun.sample.cor(X2, x2x6_cor)
      X7 = rnorm(cl_size[j], 0, sqrt(x7_var))
      X8 = fun.sample.cor(X3, x3x8_cor)
      X9 = fun.sample.cor(X4, x4x9_cor)
      X10 = rnorm(cl_size[j], 0, sqrt(x10_var))
    
    #~~~ Dichotomize X1, X2, X5, X6, X8, and X9 (will attenuate correlation)
      X11 = ifelse(X1 > mean(X1), 1, 0)  
      X33 = ifelse(X3 > mean(X3), 1, 0) 
      X55 = ifelse(X5 > mean(X5), 1, 0) 
      X66 = ifelse(X6 > mean(X6), 1, 0) # proportion of ELL about 25.14%
      X88 = ifelse(X8 > mean(X8), 1, 0)
      X99 = ifelse(X9 > mean(X9), 1, 0)
    
    #~~~ pop 1
    if (scenario == "A"){
        Y_true = (b0 + b1*X11 + b2*X2 + b3*X33 + (b4+cl_eff_slope)*X4 + b5*X55 + b6*X66 + b7*X7 
                  + st_eff + cl_eff + ind_eff)
    } else
      
    #~~~ pop 2
    if (scenario == "B"){
        Y_true = (b0 + b1*X11 + b2*X2 + b3*X33 + (b4+cl_eff_slope)*X4 + b5*X55 + b6*X66 + b7*X7 
                  + b2*X2*X2 + st_eff + cl_eff + ind_eff)
    } else 
 
    #~~~ pop 3
    if (scenario == "C"){
      Y_true = (b0 + b1*X11 + b2*X2 + b3*X33 + (b4+cl_eff_slope)*X4 + b5*X55 + b6*X66 + b7*X7 
                + b2*X2*X2 + b4*X4*X4 + b7*X7*X7 + st_eff + cl_eff + ind_eff)
    } else
        
    #~~~ pop 4
    if (scenario == "D"){
       Y_true = (b0 + b1*X11 + b2*X2 + b3*X33 + (b4+cl_eff_slope)*X4 + b5*X55 + b6*X66 + b7*X7 
                 + b1*0.5*X11*X33 + b2*0.7*X2*X4 + b4*0.5*X4*X55 + b5*0.5*X55*X66 
                 + st_eff + cl_eff + ind_eff)
    } else
       
    #~~~ pop 5
    if (scenario == "E"){
      Y_true = (b0 + b1*X11 + b2*X2 + b3*X33 + (b4+cl_eff_slope)*X4 + b5*X55 + b6*X66 + b7*X7 
                + b2*X2*X2 
                + b1*0.5*X11*X33 + b2*0.7*X2*X4 + b4*0.5*X4*X55 + b5*0.5*X55*X66 
                + st_eff + cl_eff + ind_eff)
    } else
    
    #~~~ pop 6
    if (scenario == "F"){
      Y_true = (b0 + b1*X11 + b2*X2 + b3*X33 + (b4+cl_eff_slope)*X4 + b5*X55 + b6*X66 + b7*X7 
                + b1*0.5*X11*X33 + b2*0.7*X2*X4 + b3*0.5*X33*X55 + b4*0.7*X4*X66 + b5*0.5*X55*X7 
                + b1*0.5*X11*X66 + b2*0.7*X2*X33 + b3*0.5*X33*X4 + b4*0.5*X4*X55 + b5*0.5*X55*X66
                + st_eff + cl_eff + ind_eff) 
    } else
      
    #~~~ pop7
    {
      Y_true = (b0 + b1*X11 + b2*X2 + b3*X33 + (b4+cl_eff_slope)*X4 + b5*X55 + b6*X66 + b7*X7 
                + b2*X2*X2 + b4*X4*X4 + b7*X7*X7 
                + b1*0.5*X11*X33 + b2*0.7*X2*X4 + b3*0.5*X33*X55 + b4*0.7*X4*X66 + b5*0.5*X55*X7 
                + b1*0.5*X11*X66 + b2*0.7*X2*X33 + b4*0.5*X4*X55 + b5*0.5*X55*X66
                + st_eff + cl_eff + ind_eff)      
    }

      # PS_true = 1/(1+exp(Y_true)^(-1))
      PS_true = exp(Y_true)/(1+exp(Y_true))
    
    #~~~ Assign treatment based on PS_true
      prob.exposure = runif(cl_size[j])
      treat = ifelse(PS_true > prob.exposure, 1, 0)
      # treat = ifelse(PS_true > 0.5, 1, 0) # this is wrong!!!
    
    #~~~ Treatment effect estimation model
    #~~~ the same random effects for the outcome?
    outcome = (te*treat + a0 + a1*X11 + a2*X2 + a3*X33 + a4*X4 + a5*X88 + a6*X99 +a7*X10 
               + st_eff2 + cl_eff2 + ind_eff2)

    #~~~ Stratify students within each school (ELL)
      ELL = ifelse(X6 < mean(X6) - 0.6745, 1, 0) # z = -0.6745, proportion 25%
    

#~~~~~~~~~~~ Save all variables in each school to a matrix

        temp = as.data.frame(matrix(NA, nrow = cl_size[j], ncol = num_col))
        colnames(temp) = columnnames        

        temp[, "Student"] = 1: cl_size[j]
        temp[, "School"] = rep(j, nrow(temp))
        temp[, "STRID"] = rep(k, nrow(temp))
        temp[, "STUID"] = k*100000 + j*1000 + temp[, "Student"]
        temp[, "SCHID"] = k*100 + j
        temp[, "n_j"] = rep(cl_size[j], nrow(temp))
        temp[, "Nh"] = rep(Nh, nrow(temp))
        temp[, "Y_true"] = Y_true
        temp[, "PS_true"] = PS_true
        temp[, "Treat"] = treat
        temp[, "Outcome"] = outcome
        temp[, "X1"] = X1
        temp[, "X2"] = X2
        temp[, "X3"] = X3
        temp[, "X4"] = X4
        temp[, "X5"] = X5
        temp[, "X6"] = X6
        temp[, "X7"] = X7
        temp[, "X8"] = X8
        temp[, "X9"] = X9
        temp[, "X10"] = X10
        temp[, "ELL"] = ELL
        temp[, "St_eff"] = st_eff
        temp[, "Cl_eff"] = cl_eff
        temp[, "Ind_eff"] = ind_eff
        temp[, "X11"] = X11
        temp[, "X33"] = X33
        temp[, "X55"] = X55
        temp[, "X66"] = X66
        temp[, "X88"] = X88
        temp[, "X99"] = X99        
        mydata.h[index.h:(index.h + cl_size[j]-1), ] = temp
        index.h = index.h + cl_size[j]

  } # end of j loop
    

    colnames(mydata.h) = colnames(temp)

    # calculate CLP (the prop for sampling proportional to size) for each strata : 
    # clustersize/total clustersize for either pri or pub
    # aggregate by school first (so you can directly add up clustersize)

    aggsch = ddply(mydata.h, .(SCHID), summarize, 
                      n_j = mean(n_j),
                      cl_eff = mean(Cl_eff),
                      PRI = NA,
                      CLP = NA)
    # for each strata, set the top 10 cl_eff to be private schools.
    aggsch$PRI [ which (aggsch$cl_eff %in%  aggsch[order(aggsch$cl_eff, decreasing = TRUE)[1:10], ]$cl_eff )   ] = 1
    aggsch$PRI [ which (aggsch$cl_eff %in%  aggsch[order(aggsch$cl_eff, decreasing = TRUE)[11:30], ]$cl_eff)   ] = 0

    aggsch$CLP[which(aggsch$PRI==1)] = aggsch$n_j[which(aggsch$PRI==1)]/
                                           sum(aggsch$PRI[which(aggsch$PRI==1)]*aggsch$n_j[which(aggsch$PRI==1)])
                                           # if private, CLP = n_j/sum (n_j) given private
    aggsch$CLP[which(aggsch$PRI==0)] = aggsch$n_j[which(aggsch$PRI==0)]/
                                           sum(abs(aggsch$PRI[which(aggsch$PRI==0)]-1)*aggsch$n_j[which(aggsch$PRI==0)])

    mydata.h.new = merge(mydata.h, aggsch[,c("SCHID", "PRI", "CLP")], by = "SCHID")

    mydata[index:(index + st_size[k] -1), ] = mydata.h.new
    index = index + st_size[k]

} # end of k loop

colnames(mydata) = colnames(mydata.h.new)
return(mydata)

} # end of F.GENERATE

pop1 = F.GENERATE("A", seed = 999)
pop2 = F.GENERATE("B", seed = 999+100)
pop3 = F.GENERATE("C", seed = 999+200)
pop4 = F.GENERATE("D", seed = 999+300)
pop5 = F.GENERATE("E", seed = 999+400)
pop6 = F.GENERATE("F", seed = 999+500)
pop7 = F.GENERATE("G", seed = 999+600)


write.table(pop1[1:74882, ], "pop1.csv", sep = ",", col.names = TRUE, row.names = FALSE)
write.csv(pop2[1:75638, ], "pop2.csv", row.names = FALSE)
write.csv(pop3[1:75164, ], "pop3.csv", row.names = FALSE)
write.csv(pop4[1:74894, ], "pop4.csv", row.names = FALSE)
write.csv(pop5[1:74821, ], "pop5.csv", row.names = FALSE)
write.csv(pop6[1:74769, ], "pop6.csv", row.names = FALSE)
write.csv(pop7[1:74614, ], "pop7.csv", row.names = FALSE)




######################################################
#~~~~~~~~~~~~~~~~~~~~~~ Check ~~~~~~~~~~~~~~~~~~~~~~~#
######################################################

pop1 = read.csv("pop1.csv")[1:74882, ]
pop7 = read.csv("pop7.csv")


#~~~ Check descriptives, correlation, and coefficients ~~~#

# For both: ELL = 0.2497, PRI = 0.3353
# pop1: TREAT = 0.4837; pop7: TREAT = 0.5184

descrip1 = round (describe(pop1), 4)
descrip7 = round (describe(pop7), 4)

# run a correlation matrix
# pretty close the the generated corrs

corr1 = cor(pop1[, 12:27]) # deflated to about 0.16 and 0.71
corr7 = cor(pop7[, 12:27])

# do the coefficients match?
# how close should they be?
# perfectly match if no random effects
# a little (how much) when there is random effects

fit_check_ps = lm(Y_true ~ X11+X2+X33+X4+X55+X66+X7, data = pop1)
summary(fit_check_ps)

fit_check_te = lm(Outcome ~ Treat+X11+X2+X33+X4+X88+X99+X10, data = pop1)
summary(fit_check_te) # of course this is bad! unbalanced treat groups!
# #~~~ coefficients
# b0 = 0
# b1 = 0.8
# b2 = -0.25
# b3 = 0.6
# b4 = -0.4
# b5 = -0.8
# b6 = -0.5
# b7 = 0.7
# 
# #~~ coefficients for the treatment effect formula
# a0 = -3.85
# a1 = 0.3
# a2 = -0.36
# a3 = -0.73
# a4 = -0.2
# a5 = -0.71
# a6 = -0.19
# a7 = 0.26
# te = -0.4 # treatment effect


# null = lme(fixed = Y_true~1, random = ~1|STRID/SCHID,
#            data = pop1)
# summary(null)

###~~~ the empirical cl and st var are
###~~~0.5300 and 0.5543
###~~~ In order to achieve ICC = 0.25, fix the above values and change yvar_ind only


fit_3lev_null = lmer(Y_true ~ (1|SCHID) + (1|STRID), data = pop1)
summary(fit_3lev_null)

fit_3lev_null_pop7 = lmer(Y_true ~ (1|SCHID) + (1|STRID), data = pop7)
summary(fit_3lev_null_pop7)

# when yvar_ind = 2
# pop1
(0.5300+0.5543)/(0.5300+0.5543+3.2266) # 0.2515252
# pop7
(0.4586+0.4818)/(0.4586+0.4818+5.3196) # 0.1502236

# when yvar_ind = 1.86
# pop1
(0.5297+0.5542)/(0.5297+0.5542+3.0871) # 0.2598657
# pop7
(0.4592+0.4816)/(0.4592+0.4816+5.1787) # 0.1537381



fit_3lev = lmer(Y_true ~ X11+X2+X33+X4+X55+X66+X7 + (1|SCHID) + (1|STRID), data = pop1) 
summary(fit_3lev)
# pop1: the coefficients are close to the simulated parameters

#~~~ Check random effects ~~#

# true strata and cluster level effects are ysd_cl = 0.7171471, var = 0.5143
# true individual effects is sd = 0.7171471
# true ysd_ind + ysd_ex = 
# true PRI_pp = 0.3353

# the above fit_3lev produced random effects as well
# pretty good?

# Random effects:
#   Groups   Name        Variance Std.Dev.
# SCHID    (Intercept) 0.5314   0.7290  
# STRID    (Intercept) 0.5617   0.7495  
# Residual             1.9900   1.4107  
# Number of obs: 74882, groups:  SCHID, 1500; STRID, 50 


# The following check the random effects one by one
# Similar but different results

#~~~ Individual level: 
sd(pop1$Ind_eff) # 1.411919
# true ysd_ind = 1.414214

#~~~ school-level
stat_sch = (ddply(pop1, .(SCHID), summarize, 
                     strid = mean(STRID),
                     cl_size = mean(n_j),
                     cl_eff = mean(Cl_eff),
                     PRI = mean(PRI),
                     ELL = mean(ELL))) 
head(stat_sch)


# the cl_eff column
mean(stat_sch$PRI) # check prob. of private schools across all schools
                   # 0.3333
mean(stat_sch$ELL) # check prob. of ELL students across all schools
                   # 0.2497
                   # ELL: each school should have about 0.25. Since school sizes are different,
                   # ELL prop for each school fluctuates, the mean is around 0.25 
names(stat_sch)

CL_EFF = rbind(mean = apply(stat_sch[,3:5], 2, mean), 
      sd = apply(stat_sch[,3:5], 2, sd)) 
CL_EFF

# true ysd_cl = 0.7171471

# cl_size     cl_eff       PRI
# mean 49.92133 -0.0218220 0.3333333
# sd   11.34235  0.7254926 0.4715617

#~~~ strata-level

stat_str = ddply(pop1, .(STRID), summarize, 
                 st_eff = mean(St_eff), 
                 ave.cl_eff = sd(Cl_eff),
                 st_size = length(unique(SCHID)),
                 pri_size = mean(PRI))
STR_EFF = sd(stat_str[, 2])
names(stat_str)
mean(stat_str[,"st_eff"])
mean(stat_str[, "pri_size"]) # check the prob. of private schools across all counties 0.3354155
STR_EFF # 0.787755
# CL_EFF.MEAN = mean(stat_str[,3]) # averaged cluster effect across strata # 0.7077405


# STR_EFF = 0.787755 (across 50 counties), CL_EFF = 0.7254926 (across all 1500 schools), true ysd_cl = 0.7171471.
# not perfectly match; only 50 counties and 30 schools per county.


mydata = pop1
names(mydata)
data1 = subset(mydata, Treat==1)
descrip1 = describe(data1)
data2 = subset(mydata, Treat==0)
descrip2 = describe(data2)




