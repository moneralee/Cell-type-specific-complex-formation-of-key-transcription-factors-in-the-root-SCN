library(viridis)
library(fmsb)
library(ggplot2)

time<-1000
zoomplot<-0.2

##############################################################################################################
##################### Predicting binding rates from experimentally determined relative binding rates #########
##############################################################################################################

########################## Heterodimers and Trimeric complexes - assuming donor represents a dimeric complex that has been formed (NLS seq) ############
#Protein complexes, and relative binding rates
names2complexes<-c("WOX5-PLT3","BRAVO-PLT3","BRAVO-WOX5","BRAVO-PLT3_prion","WOX5-PLT3+BRAVO","BRAVO-PLT3+WOX5","WOX5-PLT3dPrD")
target2Complex<-c(
  0.84166 #  0.96 #WOX5PLT3 
  ,0.65 #BRAVOPLT3
  ,0.61 #BRAVOWOX5
  ,0.38292011 #BRAVO-PLTprion
  ,1  #WOX5PLT3,BRAVO 
  ,0.82 #BRAVOPLT3,WOX5
  ,0.5093 #WOX5-PLTprion
)
# Here we only have data of how BRAVO-PLT+WOX5 and WOX5-PLT3+BRAVO leads to the formation of WOX5-BRAVO-PLT3 trimeric complex 
# the BRAVO-PLT3 + WOX5 was not measured

target2Complex<-target2Complex*0.8 #Relative binding rates 

#Range of association and dissociation rates considered
ass<-seq(0.01,0.5,by=0.0005/3)
diss<-seq(0.01,0.5,by=0.0005/3)

#Matrices to store ass and diss parameters that can produce reported binding rates
asssuccess=matrix(ncol=length(ass)*length(diss),nrow=length(target2Complex),"NA")
disssuccess=matrix(ncol=length(ass)*length(diss),nrow=length(target2Complex),"NA")
cont2<-rep(0,length(target2Complex))
pdf(paste("1-proteinComplexes.pdf"),width=10,height=10)
coloresassdiss<-magma(length(ass)+length(diss))
par(mfrow=c(3,3))
for(c in 1:length(target2Complex)){
  plot(-100,-100,ylim=c(0,0.5),xlim=c(0,0.5),type="l",main=names2complexes[c],xlab="Association rate",ylab="Dissociation rate")
  for(i in 1:length(ass)){
    for(j in 1:length(diss)){
      #Assume equal levels of donor and acceptor, try different rates and select the one that produces reported binding
      target<-1
      donor<-1
      targetdonor<-0
      for(t in 1:time){
        target[t+1]<-target[t]-ass[i]*target[t]*donor[t]+diss[j]*targetdonor[t]
        donor[t+1]<-donor[t]-ass[i]*target[t]*donor[t]+diss[j]*targetdonor[t]
        targetdonor[t+1]<-targetdonor[t]+ass[i]*target[t]*donor[t]-diss[j]*targetdonor[t]
      }
      #Evaluation: is the protein complex levels similar (0.00001 deviation allowed) and is it not changing (steady state)?
      if((abs(targetdonor[t+1]-target2Complex[c])<0.00001)&((targetdonor[t+1]-targetdonor[t])==0)){ 
        cont2[c]<-cont2[c]+1
        asssuccess[c,cont2[c]]<-as.double(ass[i])
        disssuccess[c,cont2[c]]<-as.double(diss[j])
        points(ass[i],diss[j],pch=20,col=coloresassdiss[i+j])
      }
    }
  }
}
dev.off() #We end pdf saving ass diss plots for all predicted complexes (di and tri)

# Saving successful binding rates
outputRates=matrix(ncol=max(cont2)+1,nrow=length(names2complexes)*2+6+2,'')

# WOX5PLT3
outputRates[1,1]<-"WOX5PLT3 ass"
outputRates[1,2:(1+cont2[1])]<-asssuccess[1,1:cont2[1]]
outputRates[2,1]<-"WOX5PLT3 diss"
outputRates[2,2:(1+cont2[1])]<-disssuccess[1,1:cont2[1]]

# BRAVOPLT3
outputRates[4,1]<-"BRAVOPLT3 ass"
outputRates[4,2:(1+cont2[2])]<-asssuccess[2,1:cont2[2]]
outputRates[5,1]<-"BRAVOPLT3 diss"
outputRates[5,2:(1+cont2[2])]<-disssuccess[2,1:cont2[2]]

# BRAVOWOX5
outputRates[7,1]<-"BRAVOWOX5 ass"
outputRates[7,2:(1+cont2[3])]<-asssuccess[3,1:cont2[3]]
outputRates[8,1]<-"BRAVOWOX5 diss"
outputRates[8,2:(1+cont2[3])]<-disssuccess[3,1:cont2[3]]

# BRAVOPLT3dPrD
outputRates[10,1]<-"BRAVOPLT3dPrD ass"
outputRates[10,2:(1+cont2[4])]<-asssuccess[4,1:cont2[4]]
outputRates[11,1]<-"BRAVOPLT3dPrD diss"
outputRates[11,2:(1+cont2[4])]<-disssuccess[4,1:cont2[4]]

# WOX5PLT3+BRAVO
outputRates[13,1]<-"WOX5PLT3+BRAVO ass"
outputRates[13,2:(1+cont2[5])]<-asssuccess[5,1:cont2[5]]
outputRates[14,1]<-"WOX5PLT3+BRAVO diss"
outputRates[14,2:(1+cont2[5])]<-disssuccess[5,1:cont2[5]]

# BRAVOPLT3+WOX5
outputRates[16,1]<-"BRAVOPLT3+WOX5 ass"
outputRates[16,2:(1+cont2[6])]<-asssuccess[6,1:cont2[6]]
outputRates[17,1]<-"BRAVOPLT3+WOX5 diss"
outputRates[17,2:(1+cont2[6])]<-disssuccess[6,1:cont2[6]]

# WOX5PLT3dPrD
outputRates[19,1]<-"WOX5PLT3dPrD ass"
outputRates[19,2:(1+cont2[7])]<-asssuccess[7,1:cont2[7]]
outputRates[20,1]<-"WOX5PLT3dPrD diss"
outputRates[20,2:(1+cont2[7])]<-disssuccess[7,1:cont2[7]]

write.table(outputRates,"BindingRates.csv",sep=",",col.names = F,row.names = F)
saveRDS(cont2,"bindingRates.rds")
