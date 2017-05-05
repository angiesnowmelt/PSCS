setwd("J:\\PSCS\\SIM2")


library(pps)
library(plyr)
library(psych)
# library(sampling)

F.sample = function(mydata, nh_c, nh_d, nh_a, nh_b){
  
  # mydata = read.csv("pop7.csv") 
  
  #~~~~~~~~~~~~~~~~~~~~~ A function for sampling schools from counties ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # First, I need to aggregate by school first
  
  stat_sch = (ddply(mydata, .(SCHID), summarize, 
                    strid = mean(STRID),
                    cl_size = mean(n_j),
                    CLP = mean(CLP),
                    cl_eff = mean(Cl_eff),
                    PRI = mean(PRI),
                    ELL = mean(ELL))) 
  agg_sch = stat_sch # 50*30 = 1500 schools in total
  agg_sch[, "SCH_P"] = NA
  agg_sch[, "SCHWT"] = NA
  # head(agg_sch)
  
  
  
  samp.schools = function(whichcounty, nh_c, nh_d){
    
  #   nh_c = 5
  #   nh_d = 5
  #   whichcounty = unique(agg_sch$strid)[1]
    
    single_county = agg_sch[agg_sch$strid == whichcounty, ] # get the units in the first selected schools
    table(single_county$PRI)
    
    # if select by probability, then replace nh_c and nh_d with the following values in the function.
  #   nh_c = e * as.numeric(table(single_county$PRI))[1] # nh_c = 0.3
  #   nh_d = f * as.numeric(table(single_county$PRI))[2] # nh_d = 0.3
    
    single_county = single_county[order(single_county$PRI),] # order by strata (private schools)
    
    num_pub = ifelse(as.numeric(table(single_county$PRI))[1] >= nh_c, nh_c, as.numeric(table(single_county$PRI))[1])
    num_pri = ifelse(as.numeric(table(single_county$PRI))[2] >= nh_d, nh_d, as.numeric(table(single_county$PRI))[2])
    nh_pri = c(num_pub, num_pri)
  
    # add selection weights for private and public schools 
    single_county$SCH_P[which(single_county$PRI == 1)] <- single_county$CLP[which(single_county$PRI == 1)]*num_pri
    single_county$SCH_P[which(single_county$PRI == 0)] <- single_county$CLP[which(single_county$PRI == 0)]*num_pub
    single_county$SCHWT = 1/single_county$SCH_P
  
    samp_sch = single_county[ppssstrat(single_county$cl_size, single_county$PRI, nh_pri), ]
  }
  
  whichcounty= as.list(unique(mydata$STRID))
  
  samp_sch_all = ldply(whichcounty, samp.schools, nh_c, nh_d)
  
  sum(samp_sch_all$SCHWT) 
  
  # get the dataset with selected schools (with all students)
  mydata_selsch = subset(mydata, SCHID %in% samp_sch_all$SCHID)
  mydata_selsch[, "STU_P"] = NA
  mydata_selsch[, "STUWT"] = NA
  
  
  #~~~~~~~~~~~~~~~~~~~~~ A function for sampling students from schools ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  samp.students = function(whichschool, nh_a, nh_b){
    # nh_a = 10
    # nh_b = 10
    # whichschool = unique(mydata_selsch$SCHID)[1]
    single_sch = mydata_selsch[mydata_selsch$SCHID == whichschool, ] # get the units in the first selected schools
    table(single_sch$ELL)
    single_sch = single_sch[order(single_sch$ELL),] # order by strata (ELL students)
    
    num_nonell = ifelse(as.numeric(table(single_sch$ELL))[1] >= nh_a, nh_a, as.numeric(table(single_sch$ELL))[1])
    num_ell = ifelse(as.numeric(table(single_sch$ELL))[2] >= nh_b, nh_b, as.numeric(table(single_sch$ELL))[2])
    nh_ell = c(num_nonell, num_ell)
    
    # add selection weights for ELL and non-ELL students
    
    single_sch$STU_P[which(single_sch$ELL==1)] <- num_ell/as.numeric(table(single_sch$ELL))[2]
    single_sch$STU_P[which(single_sch$ELL==0)] <- num_nonell/as.numeric(table(single_sch$ELL))[1]
    single_sch$STUWT = 1/single_sch$STU_P
    
    #   single_sch$STU_P = ifelse(single_sch$ELL==1, num_ell/as.numeric(table(single_sch$ELL))[2], 
    #                             num_nonell/as.numeric(table(single_sch$ELL))[1]) # this didn't work
    
    samp_stu = single_sch[stratsrs(single_sch$ELL, nh_ell),]
  }
  
  whichschool= as.list(unique(mydata_selsch$SCHID))
  samp_stu_all = ldply(whichschool, samp.students, nh_a, nh_b)
  
  #~~~~~~~ merge school weight and selected students ~~~~~~~~~#
  
  samp_sch_all_sub = samp_sch_all[ ,c("SCHID", "SCH_P", "SCHWT")]
  sample = merge(samp_stu_all, samp_sch_all_sub,  by = "SCHID")
  sample[, "SAMPWT"] = NA
  sample$SAMPWT = sample$SCHWT*sample$STUWT
  return(sample)

} # end of F.sample




#################################################
####### Save Samples for Each Population ########
#################################################

F.save = function(nrep, allsamp, directory){
  r = 1
  for (r in 1:nrep) {
    mydat = paste("samp",r, sep = "")
    mydat = as.data.frame(allsamp[, r])
    myfilename = file.path(directory, paste("pop-samp-", r, ".csv", sep = ""))
    write.table(mydat, file = myfilename, sep = ",", col.names = TRUE, row.names=FALSE)
  }  
}



#########################################################


directory = "J:\\PSCS\\SIM2"
nrep = 100

directory_samp1 = "J:\\PSCS\\SIM2\\POP1_SAMP"
directory_samp2 = "J:\\PSCS\\SIM2\\POP2_SAMP"
directory_samp3 = "J:\\PSCS\\SIM2\\POP3_SAMP"
directory_samp4 = "J:\\PSCS\\SIM2\\POP4_SAMP"
directory_samp5 = "J:\\PSCS\\SIM2\\POP5_SAMP"
directory_samp6 = "J:\\PSCS\\SIM2\\POP6_SAMP"
directory_samp7 = "J:\\PSCS\\SIM2\\POP7_SAMP"

seed_samp1 = 999999
seed_samp2 = 999999*100
seed_samp3 = 999999*200
seed_samp4 = 999999*300
seed_samp5 = 999999*400
seed_samp6 = 999999*500
seed_samp7 = 999999*600



##############################
#~~~ Sample data for pop1
##############################


filename = file.path(directory, "pop1.csv") #~~~ change here
mydata = read.csv(filename)
set.seed(seed_samp1) #~~~ change here
sample = replicate(nrep, F.sample(mydata = mydata, nh_c = 5, nh_d = 5, nh_a = 10, nh_b = 10))

#~~~ save the fitted data for each samp
F.save(nrep = nrep, allsamp = sample, directory = directory_samp1) #~~~ change here



##############################
#~~~ Sample data for pop2
##############################


filename = file.path(directory, "pop2.csv") #~~~ change here
mydata = read.csv(filename)
set.seed(seed_samp2) #~~~ change here
sample = replicate(nrep, F.sample(mydata = mydata, nh_c = 5, nh_d = 5, nh_a = 10, nh_b = 10))

#~~~ save the fitted data for each samp
F.save(nrep = nrep, allsamp = sample, directory = directory_samp2) #~~~ change here



##############################
#~~~ Sample data for pop3
##############################


filename = file.path(directory, "pop3.csv") #~~~ change here
mydata = read.csv(filename)
set.seed(seed_samp3) #~~~ change here
sample = replicate(nrep, F.sample(mydata = mydata, nh_c = 5, nh_d = 5, nh_a = 10, nh_b = 10))

#~~~ save the fitted data for each samp
F.save(nrep = nrep, allsamp = sample, directory = directory_samp3) #~~~ change here



##############################
#~~~ Sample data for pop4
##############################


filename = file.path(directory, "pop4.csv") #~~~ change here
mydata = read.csv(filename)
set.seed(seed_samp4) #~~~ change here
sample = replicate(nrep, F.sample(mydata = mydata, nh_c = 5, nh_d = 5, nh_a = 10, nh_b = 10))

#~~~ save the fitted data for each samp
F.save(nrep = nrep, allsamp = sample, directory = directory_samp4) #~~~ change here


##############################
#~~~ Sample data for pop5
##############################


filename = file.path(directory, "pop5.csv") #~~~ change here
mydata = read.csv(filename)
set.seed(seed_samp5) #~~~ change here
sample = replicate(nrep, F.sample(mydata = mydata, nh_c = 5, nh_d = 5, nh_a = 10, nh_b = 10))

#~~~ save the fitted data for each samp
F.save(nrep = nrep, allsamp = sample, directory = directory_samp5) #~~~ change here


##############################
#~~~ Sample data for pop6
##############################


filename = file.path(directory, "pop6.csv") #~~~ change here
mydata = read.csv(filename)
set.seed(seed_samp6) #~~~ change here
sample = replicate(nrep, F.sample(mydata = mydata, nh_c = 5, nh_d = 5, nh_a = 10, nh_b = 10))

#~~~ save the fitted data for each samp
F.save(nrep = nrep, allsamp = sample, directory = directory_samp6) #~~~ change here



##############################
#~~~ Sample data for pop7
##############################


filename = file.path(directory, "pop7.csv") #~~~ change here
mydata = read.csv(filename)
set.seed(seed_samp7) #~~~ change here
sample = replicate(nrep, F.sample(mydata = mydata, nh_c = 5, nh_d = 5, nh_a = 10, nh_b = 10))

#~~~ save the fitted data for each samp
F.save(nrep = nrep, allsamp = sample, directory = directory_samp7) #~~~ change here

















  
  
  





  
  
########################################################
#~~~~~~~~~~~~~~~~~~ Check sample ~~~~~~~~~~~~~~~~~~~~~~#
########################################################

sum(samp_sch_all$SCHWT) # how to check this one? add up to be the total # of schools
sum(sample$SAMPWT)
  e.table(sample, "sample.csv", sep = ",", col.names = TRUE, row.names = FALSE)


descrip_pop = describe(mydata)
descrip_pop = describe(pop)
  ta = pop
descrip_samp = describe(sample)
  ite.table (descrip_pop, "descrip_pop.csv", sep = ",", col.names = NA)
write.table (descrip_samp, "descrip_samp.csv", sep = ",", col.names = NA)

sample = read.csv("sample.csv")

#~~~ Individual level
sd(mydata$Ind_eff)
sd(sample1$Ind_eff)

#~~~ school-level
stat_sch = (ddply(mydata, .(SCHID), summarize, 
                  strid = mean(STRID),
                  cl_size = mean(n_j),
                  cl_eff = mean(Cl_eff),
                  PRI = mean(PRI),
                  ELL = mean(ELL)))
head(stat_sch)

stat_sch_samp = (ddply(sample1, .(SCHID), summarize, 
                  strid = mean(STRID),
                  cl_size = mean(n_j),
                  cl_eff = mean(Cl_eff),
                  PRI = mean(PRI),
                  ELL = mean(ELL),
                  SCHWT = mean(SCHWT))) 
stat_sch_samp

# this only gives you useful info on ELL rate... most 0.5 (if there were at least 10 in that school)

stat_sch_popsamp = cbind(subset(stat_sch, SCHID %in% stat_sch_samp$SCHID), stat_sch_samp)

# the cl_eff column
cbind (rbind(mean = apply(stat_sch[,3:5], 2, mean), 
      sd = apply(stat_sch[,3:5], 2, sd)),
      rbind(mean = apply(stat_sch_samp[,3:7], 2, mean), 
      sd = apply(stat_sch_samp[, 3:7], 2, sd))) 

#~~~ strata-level
  
stat_str = ddply(mydata, .(STRID), summarize, 
                 st_eff = mean(St_eff), 
                 cl_eff = sd(Cl_eff),
                 st_size = length(unique(SCHID)),
                 pri_size = mean(PRI),
                 nh = mean(Nh))

stat_str_samp = ddply(sample1, .(STRID), summarize, 
                 st_eff = mean(St_eff), 
                 cl_eff = sd(Cl_eff),
                 st_size = length(unique(SCHID)),
               pri_size = mean(PRI),
                 nh = sum(unique(n_j)) # not guaranteed
                 #, SCHWT = sum(SCHWT)
                 )
strat_sum = cbind(stat_str, stat_str_samp)
strat_sum


# STR_EFF, of course exactly the same because the whole strata were selected
cbind(sd(stat_str[, 2]), sd(stat_str_samp[,2]))

# CL_EFF
CL_EFF.MEAN = cbind(mean(stat_str[,3]), mean(stat_str_samp[,3]))
CL_EFF.MEAN

# Check the weights

  m(samp_sch_all$SCHWT) # how to check this one? twice the total number of schools
sum(sample1$SAMPWT)

