#extra bits of old code I didn't want to delete from encounter_risk_v2.R

#compute Morisita Horn Index of overlap

both = left_join(sealeff,fisheff,by=c('GRIDCODE'))
both$Pf[is.na(both$Pf)] <- 0
# both$Ps[is.na(both$Ps)] <- 0
# both$n[is.na(both$n)] <- 0
both$PfPs<-both$Pf*both$Ps  #plot Pf*Ps like Cronin et al. Figure 1

Cmhf = 2*sum(both$Pf*both$Ps)/(sum(both$Pf^2)+ sum(both$Ps^2))
#0.1918607
#at smaller scale, Cmhf = 0.04037908


regdat = both %>% mutate(region = case_when(GRIDCODE %in% c(869, 913, 917, 919, 920, 962, 963, 967, 968, 969, 1013:1019, 1064, 1067, 1021) ~ 1,
                                            GRIDCODE %in% c(964, 914:916, 861:867, 811:816, 761:765, 712:714) ~ 2, 
                                            GRIDCODE %in% c(766, 715, 716) ~ 3,
                                            TRUE ~ 4))

#compute observed Morisita Horn Index of overlap for each region
mho1=regdat[regdat$region==1,]
Cmhf1 = 2*sum(mho1$Pf*mho1$Ps)/(sum(mho1$Pf^2)+ sum(mho1$Ps^2))  #0.1328
mho2=regdat[regdat$region==2,]
Cmhf2 = 2*sum(mho2$Pf*mho2$Ps)/(sum(mho2$Pf^2)+ sum(mho2$Ps^2))  #0.2174
mho3=regdat[regdat$region==3,]
Cmhf3 = 2*sum(mho3$Pf*mho3$Ps)/(sum(mho3$Pf^2)+ sum(mho3$Ps^2))  #0.2001
mho4=regdat[regdat$region==4,]
Cmhf4 = 2*sum(mho4$Pf*mho4$Ps)/(sum(mho4$Pf^2)+ sum(mho4$Ps^2))  #0.0363

#randomly shuffle seal data and compute statistic again:
ranseals=sample_n(seals,98852,replace=T)      #not sure if this is getting rid of temporal autocorrelation      

#compute Cmhf for each seal
#seal1<-ranseals[ranseals$id==142351,]
#s1cm=seal1 %>%
#  group_by(GRIDCODE)%>%
#  summarise(wtsum=sum(wtpos))%>%
#  mutate(Ps=wtsum/sum(wtsum))  #Ps is proportion of seal effort in each grid cell irrespective of time
#s1f = full_join(s1cm,fisheff,by=c('GRIDCODE'))
#s1f$Pf[is.na(s1f$Pf)] <- 0
#s1f$Ps[is.na(s1f$Ps)] <- 0
#s1fmh = 2*sum(s1f$Pf*s1f$Ps)/(sum(s1f$Pf^2)+ sum(s1f$Ps^2))
#s1fmh
#results recorded for each seal in deploymentsummary2.xls

#now permute the seal data within each region
#1000x, compute overlap index for each permutation, and compare Cdata to Crandom to examine likelihood
#that the actual seal effort is/is not distributed randomly
#this method cannot capture whether observed spatial pattern is significantly different from random because
#the seal tracks are always associated with the same gridcode

## function to compute Cmhf
cmf = function(d, i) {
  
  ## get data
  use_data = d[i,]
  
  new_data=use_data %>%
    group_by(GRIDCODE)%>%
    summarise(wtsum=sum(wtpos))%>%
    mutate(Ps=wtsum/sum(wtsum))  
  
  
  #join with fish effort
  both = left_join(new_data,fisheff,by=c('GRIDCODE'))
  both$Pf[is.na(both$Pf)] <- 0
  #both$Ps[is.na(both$Ps)] <- 0
  
  ## calculate the overlap statistic
  cmhf=with(both,2*sum(Pf*Ps)/(sum(Pf^2)+sum(Ps^2)))
  
  ## output the indices
  return(cmhf)
  
}

## bootstrap resampled seal data to get 1000 cmhf statistics
#boot = boot(data = seals, statistic = cmf, R = 1000)
boot1 = boot(data = sealreg, statistic = cmf, R = 1000,strata=(sealreg$region))  #remove factor


#calculate pvalue by counting # of statistics that are >= to observed value
#overall cmhf: pval = sum(boot$t>=0.1915665)/1000
#p=0.18
#stratified bootsrap:
pval1 = sum(boot1$t>=0.1918607)/1000  #0.506
