library(viridis)
library(fmsb)
library(ggplot2)
zoomplot<-0.2

rates<-read.table("BindingRates.csv",sep=",",header=F)
ratesOptions<-readRDS("bindingRates.rds")

ass<-0
diss<-0
helper<-0

#Here I am reading the binding rates, selecting one set at random per simulated complex
for(i in 1:length(ratesOptions)){
  helper[i]<-1+3*(i-1) # Because in table, there is a space in between ass/diss per complex
  sampleassdiss<-sample(1:ratesOptions[i],1) # Random selection of random rates to be used in the simulation
  ass[i]<-as.double(rates[helper[i],sampleassdiss+1]) 
  diss[i]<-as.double(rates[helper[i]+1,sampleassdiss+1]) 
}

#############################Protein complex simulation#######################
time<-2000
heatmapdata=matrix(nrow=7,ncol=4*4+4,0)
heatmapcont<-0

###### Root stem cell niche protein levels (normalized)
PLT3values<-c(0.438535048,0.421883943,0.46759333,0.156499291)
WOX5values<-c(0.289363605,0.456050971,0.210445826,0.02819338)
BRAVOvalues<-c(0.607472928,0.285295553,0.110637236,0.018408146)

# BRAVO+WOX5
abw<-ass[3]
dbw<-diss[3]
# WOX5PLT3+BRAVO
awpb1<-ass[5]
dwpb1<-diss[5]
# BRAVOPLT3+WOX5
awpb2<-ass[6]
dwpb2<-diss[6]

# WOX5+PLT3
awp<-ass[1]
dwp<-diss[1]
# BRAVO+PLT3
abp<-ass[2]
dbp<-diss[2]


BRAVO=matrix(ncol=4,nrow=time+1,0)
PLT3=matrix(ncol=4,nrow=time+1,0)
WOX5=matrix(ncol=4,nrow=time+1,0)
WOX5PLT3=matrix(ncol=4,nrow=time+1,0)
BRAVOPLT3=matrix(ncol=4,nrow=time+1,0)
BRAVOWOX5=matrix(ncol=4,nrow=time+1,0)
WOX5PLT3BRAVO=matrix(ncol=4,nrow=time+1,0) 

for(cell in 1:4){
  BRAVO[1,cell]<-BRAVOvalues[cell]
  PLT3[1,cell]<-PLT3values[cell]
  WOX5[1,cell]<-WOX5values[cell]
  WOX5PLT3[1,cell]<-0
  BRAVOPLT3[1,cell]<-0
  BRAVOWOX5[1,cell]<-0
  WOX5PLT3BRAVO[1,cell]<-0
  
  for(t in 1:time){
    BRAVO[t+1,cell]<-BRAVO[t,cell]-abw*BRAVO[t,cell]*WOX5[t,cell]+dbw*BRAVOWOX5[t,cell]-abp*BRAVO[t,cell]*PLT3[t,cell]+dbp*BRAVOPLT3[t,cell]-awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]+dwpb1*WOX5PLT3BRAVO[t,cell]
    PLT3[t+1,cell]<-PLT3[t,cell]-abp*BRAVO[t,cell]*PLT3[t,cell]+dbp*BRAVOPLT3[t,cell]-awp*WOX5[t,cell]*PLT3[t,cell]+dwp*WOX5PLT3[t,cell]
    WOX5[t+1,cell]<-WOX5[t,cell]-awp*WOX5[t,cell]*PLT3[t,cell]+dwp*WOX5PLT3[t,cell]-abw*BRAVO[t,cell]*WOX5[t,cell]+dbw*BRAVOWOX5[t,cell]-awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]+dwpb2*WOX5PLT3BRAVO[t,cell]
    WOX5PLT3[t+1,cell]<-WOX5PLT3[t,cell]+awp*WOX5[t,cell]*PLT3[t,cell]-dwp*WOX5PLT3[t,cell]-awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]+dwpb1*WOX5PLT3BRAVO[t,cell]
    BRAVOPLT3[t+1,cell]<-BRAVOPLT3[t,cell]+abp*BRAVO[t,cell]*PLT3[t,cell]-dbp*BRAVOPLT3[t,cell]-awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]+dwpb2*WOX5PLT3BRAVO[t,cell]
    BRAVOWOX5[t+1,cell]<-BRAVOWOX5[t,cell]+abw*BRAVO[t,cell]*WOX5[t,cell]-dbw*BRAVOWOX5[t,cell]
    WOX5PLT3BRAVO[t+1,cell]<-WOX5PLT3BRAVO[t,cell]+awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]-dwpb1*WOX5PLT3BRAVO[t,cell]+awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]-dwpb2*WOX5PLT3BRAVO[t,cell]
  }
}

###heatmap
heatmapdata[,heatmapcont*4+1]<-c(WOX5PLT3[t,1],BRAVOPLT3[t,1],BRAVOWOX5[t,1],WOX5PLT3BRAVO[t,1],
                                 WOX5[t,1],PLT3[t,1],BRAVO[t,1])
heatmapdata[,heatmapcont*4+2]<-c(WOX5PLT3[t,2],BRAVOPLT3[t,2],BRAVOWOX5[t,2],WOX5PLT3BRAVO[t,2],
                                 WOX5[t,2],PLT3[t,2],BRAVO[t,2])
heatmapdata[,heatmapcont*4+3]<-c(WOX5PLT3[t,3],BRAVOPLT3[t,3],BRAVOWOX5[t,3],WOX5PLT3BRAVO[t,3],
                                 WOX5[t,3],PLT3[t,3],BRAVO[t,3])
heatmapdata[,heatmapcont*4+4]<-c(WOX5PLT3[t,4],BRAVOPLT3[t,4],BRAVOWOX5[t,4],WOX5PLT3BRAVO[t,4],
                                 WOX5[t,4],PLT3[t,4],BRAVO[t,4])
heatmapcont<-heatmapcont+1




########### Control1-Ass-diss-same ###############
abw<-0.1
dbw<-0.1
# WOX5PLT3+BRAVO
awpb1<-0.1
dwpb1<-0.1
# BRAVOPLT3+WOX5
awpb2<-0.1
dwpb2<-0.1
# WOX5+PLT3
awp<-0.1
dwp<-0.1
# BRAVO+PLT3
abp<-0.1
dbp<-0.1


BRAVO=matrix(ncol=4,nrow=time+1,0)
PLT3=matrix(ncol=4,nrow=time+1,0)
WOX5=matrix(ncol=4,nrow=time+1,0)
WOX5PLT3=matrix(ncol=4,nrow=time+1,0)
BRAVOPLT3=matrix(ncol=4,nrow=time+1,0)
BRAVOWOX5=matrix(ncol=4,nrow=time+1,0)
WOX5PLT3BRAVO=matrix(ncol=4,nrow=time+1,0) 

for(cell in 1:4){
  BRAVO[1,cell]<-BRAVOvalues[cell]
  PLT3[1,cell]<-PLT3values[cell]
  WOX5[1,cell]<-WOX5values[cell]
  WOX5PLT3[1,cell]<-0
  BRAVOPLT3[1,cell]<-0
  BRAVOWOX5[1,cell]<-0
  WOX5PLT3BRAVO[1,cell]<-0
  
  for(t in 1:time){
    BRAVO[t+1,cell]<-BRAVO[t,cell]-abw*BRAVO[t,cell]*WOX5[t,cell]+dbw*BRAVOWOX5[t,cell]-abp*BRAVO[t,cell]*PLT3[t,cell]+dbp*BRAVOPLT3[t,cell]-awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]+dwpb1*WOX5PLT3BRAVO[t,cell]
    PLT3[t+1,cell]<-PLT3[t,cell]-abp*BRAVO[t,cell]*PLT3[t,cell]+dbp*BRAVOPLT3[t,cell]-awp*WOX5[t,cell]*PLT3[t,cell]+dwp*WOX5PLT3[t,cell]
    WOX5[t+1,cell]<-WOX5[t,cell]-awp*WOX5[t,cell]*PLT3[t,cell]+dwp*WOX5PLT3[t,cell]-abw*BRAVO[t,cell]*WOX5[t,cell]+dbw*BRAVOWOX5[t,cell]-awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]+dwpb2*WOX5PLT3BRAVO[t,cell]
    WOX5PLT3[t+1,cell]<-WOX5PLT3[t,cell]+awp*WOX5[t,cell]*PLT3[t,cell]-dwp*WOX5PLT3[t,cell]-awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]+dwpb1*WOX5PLT3BRAVO[t,cell]
    BRAVOPLT3[t+1,cell]<-BRAVOPLT3[t,cell]+abp*BRAVO[t,cell]*PLT3[t,cell]-dbp*BRAVOPLT3[t,cell]-awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]+dwpb2*WOX5PLT3BRAVO[t,cell]
    BRAVOWOX5[t+1,cell]<-BRAVOWOX5[t,cell]+abw*BRAVO[t,cell]*WOX5[t,cell]-dbw*BRAVOWOX5[t,cell]
    WOX5PLT3BRAVO[t+1,cell]<-WOX5PLT3BRAVO[t,cell]+awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]-dwpb1*WOX5PLT3BRAVO[t,cell]+awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]-dwpb2*WOX5PLT3BRAVO[t,cell]
  }
}

###heatmap
heatmapdata[,1+heatmapcont*4]<-c(WOX5PLT3[t,1],BRAVOPLT3[t,1],BRAVOWOX5[t,1],WOX5PLT3BRAVO[t,1],
                                 WOX5[t,1],PLT3[t,1],BRAVO[t,1])
heatmapdata[,2+heatmapcont*4]<-c(WOX5PLT3[t,2],BRAVOPLT3[t,2],BRAVOWOX5[t,2],WOX5PLT3BRAVO[t,2],
                                 WOX5[t,2],PLT3[t,2],BRAVO[t,2])
heatmapdata[,3+heatmapcont*4]<-c(WOX5PLT3[t,3],BRAVOPLT3[t,3],BRAVOWOX5[t,3],WOX5PLT3BRAVO[t,3],
                                 WOX5[t,3],PLT3[t,3],BRAVO[t,3])
heatmapdata[,4+heatmapcont*4]<-c(WOX5PLT3[t,4],BRAVOPLT3[t,4],BRAVOWOX5[t,4],WOX5PLT3BRAVO[t,4],
                                 WOX5[t,4],PLT3[t,4],BRAVO[t,4])
heatmapcont<-heatmapcont+1

########### Control2-Ass>diss ###############
abw<-0.1
dbw<-0.05
# WOX5PLT3+BRAVO
awpb1<-0.1
dwpb1<-0.05
# BRAVOPLT3+WOX5
awpb2<-0.1
dwpb2<-0.05
# WOX5+PLT3
awp<-0.1
dwp<-0.05
# BRAVO+PLT3
abp<-0.1
dbp<-0.05



BRAVO=matrix(ncol=4,nrow=time+1,0)
PLT3=matrix(ncol=4,nrow=time+1,0)
WOX5=matrix(ncol=4,nrow=time+1,0)
WOX5PLT3=matrix(ncol=4,nrow=time+1,0)
BRAVOPLT3=matrix(ncol=4,nrow=time+1,0)
BRAVOWOX5=matrix(ncol=4,nrow=time+1,0)
WOX5PLT3BRAVO=matrix(ncol=4,nrow=time+1,0) 

for(cell in 1:4){
  BRAVO[1,cell]<-BRAVOvalues[cell]
  PLT3[1,cell]<-PLT3values[cell]
  WOX5[1,cell]<-WOX5values[cell]
  WOX5PLT3[1,cell]<-0
  BRAVOPLT3[1,cell]<-0
  BRAVOWOX5[1,cell]<-0
  WOX5PLT3BRAVO[1,cell]<-0
  
  for(t in 1:time){
    BRAVO[t+1,cell]<-BRAVO[t,cell]-abw*BRAVO[t,cell]*WOX5[t,cell]+dbw*BRAVOWOX5[t,cell]-abp*BRAVO[t,cell]*PLT3[t,cell]+dbp*BRAVOPLT3[t,cell]-awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]+dwpb1*WOX5PLT3BRAVO[t,cell]
    PLT3[t+1,cell]<-PLT3[t,cell]-abp*BRAVO[t,cell]*PLT3[t,cell]+dbp*BRAVOPLT3[t,cell]-awp*WOX5[t,cell]*PLT3[t,cell]+dwp*WOX5PLT3[t,cell]
    WOX5[t+1,cell]<-WOX5[t,cell]-awp*WOX5[t,cell]*PLT3[t,cell]+dwp*WOX5PLT3[t,cell]-abw*BRAVO[t,cell]*WOX5[t,cell]+dbw*BRAVOWOX5[t,cell]-awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]+dwpb2*WOX5PLT3BRAVO[t,cell]
    WOX5PLT3[t+1,cell]<-WOX5PLT3[t,cell]+awp*WOX5[t,cell]*PLT3[t,cell]-dwp*WOX5PLT3[t,cell]-awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]+dwpb1*WOX5PLT3BRAVO[t,cell]
    BRAVOPLT3[t+1,cell]<-BRAVOPLT3[t,cell]+abp*BRAVO[t,cell]*PLT3[t,cell]-dbp*BRAVOPLT3[t,cell]-awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]+dwpb2*WOX5PLT3BRAVO[t,cell]
    BRAVOWOX5[t+1,cell]<-BRAVOWOX5[t,cell]+abw*BRAVO[t,cell]*WOX5[t,cell]-dbw*BRAVOWOX5[t,cell]
    WOX5PLT3BRAVO[t+1,cell]<-WOX5PLT3BRAVO[t,cell]+awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]-dwpb1*WOX5PLT3BRAVO[t,cell]+awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]-dwpb2*WOX5PLT3BRAVO[t,cell]
  }
}
###heatmap
heatmapdata[,1+heatmapcont*4]<-c(WOX5PLT3[t,1],BRAVOPLT3[t,1],BRAVOWOX5[t,1],WOX5PLT3BRAVO[t,1],
                                 WOX5[t,1],PLT3[t,1],BRAVO[t,1])
heatmapdata[,2+heatmapcont*4]<-c(WOX5PLT3[t,2],BRAVOPLT3[t,2],BRAVOWOX5[t,2],WOX5PLT3BRAVO[t,2],
                                 WOX5[t,2],PLT3[t,2],BRAVO[t,2])
heatmapdata[,3+heatmapcont*4]<-c(WOX5PLT3[t,3],BRAVOPLT3[t,3],BRAVOWOX5[t,3],WOX5PLT3BRAVO[t,3],
                                 WOX5[t,3],PLT3[t,3],BRAVO[t,3])
heatmapdata[,4+heatmapcont*4]<-c(WOX5PLT3[t,4],BRAVOPLT3[t,4],BRAVOWOX5[t,4],WOX5PLT3BRAVO[t,4],
                                 WOX5[t,4],PLT3[t,4],BRAVO[t,4])
heatmapcont<-heatmapcont+1

###########Control3-Ass<diss###############
abw<-0.05
dbw<-0.1
# WOX5PLT3+BRAVO
awpb1<-0.05
dwpb1<-0.1
# BRAVOPLT3+WOX5
awpb2<-0.05
dwpb2<-0.1
# WOX5+PLT3
awp<-0.05
dwp<-0.1
# BRAVO+PLT3
abp<-0.05
dbp<-0.1


BRAVO=matrix(ncol=4,nrow=time+1,0)
PLT3=matrix(ncol=4,nrow=time+1,0)
WOX5=matrix(ncol=4,nrow=time+1,0)
WOX5PLT3=matrix(ncol=4,nrow=time+1,0)
BRAVOPLT3=matrix(ncol=4,nrow=time+1,0)
BRAVOWOX5=matrix(ncol=4,nrow=time+1,0)
WOX5PLT3BRAVO=matrix(ncol=4,nrow=time+1,0) #FORMED BY ASS 1 AND 5

for(cell in 1:4){
  BRAVO[1,cell]<-BRAVOvalues[cell]
  PLT3[1,cell]<-PLT3values[cell]
  WOX5[1,cell]<-WOX5values[cell]
  WOX5PLT3[1,cell]<-0
  BRAVOPLT3[1,cell]<-0
  BRAVOWOX5[1,cell]<-0
  WOX5PLT3BRAVO[1,cell]<-0
  
  for(t in 1:time){
    BRAVO[t+1,cell]<-BRAVO[t,cell]-abw*BRAVO[t,cell]*WOX5[t,cell]+dbw*BRAVOWOX5[t,cell]-abp*BRAVO[t,cell]*PLT3[t,cell]+dbp*BRAVOPLT3[t,cell]-awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]+dwpb1*WOX5PLT3BRAVO[t,cell]
    PLT3[t+1,cell]<-PLT3[t,cell]-abp*BRAVO[t,cell]*PLT3[t,cell]+dbp*BRAVOPLT3[t,cell]-awp*WOX5[t,cell]*PLT3[t,cell]+dwp*WOX5PLT3[t,cell]
    WOX5[t+1,cell]<-WOX5[t,cell]-awp*WOX5[t,cell]*PLT3[t,cell]+dwp*WOX5PLT3[t,cell]-abw*BRAVO[t,cell]*WOX5[t,cell]+dbw*BRAVOWOX5[t,cell]-awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]+dwpb2*WOX5PLT3BRAVO[t,cell]
    WOX5PLT3[t+1,cell]<-WOX5PLT3[t,cell]+awp*WOX5[t,cell]*PLT3[t,cell]-dwp*WOX5PLT3[t,cell]-awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]+dwpb1*WOX5PLT3BRAVO[t,cell]
    BRAVOPLT3[t+1,cell]<-BRAVOPLT3[t,cell]+abp*BRAVO[t,cell]*PLT3[t,cell]-dbp*BRAVOPLT3[t,cell]-awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]+dwpb2*WOX5PLT3BRAVO[t,cell]
    BRAVOWOX5[t+1,cell]<-BRAVOWOX5[t,cell]+abw*BRAVO[t,cell]*WOX5[t,cell]-dbw*BRAVOWOX5[t,cell]
    WOX5PLT3BRAVO[t+1,cell]<-WOX5PLT3BRAVO[t,cell]+awpb1*WOX5PLT3[t,cell]*BRAVO[t,cell]-dwpb1*WOX5PLT3BRAVO[t,cell]+awpb2*WOX5[t,cell]*BRAVOPLT3[t,cell]-dwpb2*WOX5PLT3BRAVO[t,cell]
  }
}

###heatmap
heatmapdata[,1+heatmapcont*4]<-c(WOX5PLT3[t,1],BRAVOPLT3[t,1],BRAVOWOX5[t,1],WOX5PLT3BRAVO[t,1],
                                 WOX5[t,1],PLT3[t,1],BRAVO[t,1])
heatmapdata[,2+heatmapcont*4]<-c(WOX5PLT3[t,2],BRAVOPLT3[t,2],BRAVOWOX5[t,2],WOX5PLT3BRAVO[t,2],
                                 WOX5[t,2],PLT3[t,2],BRAVO[t,2])
heatmapdata[,3+heatmapcont*4]<-c(WOX5PLT3[t,3],BRAVOPLT3[t,3],BRAVOWOX5[t,3],WOX5PLT3BRAVO[t,3],
                                 WOX5[t,3],PLT3[t,3],BRAVO[t,3])
heatmapdata[,4+heatmapcont*4]<-c(WOX5PLT3[t,4],BRAVOPLT3[t,4],BRAVOWOX5[t,4],WOX5PLT3BRAVO[t,4],
                                 WOX5[t,4],PLT3[t,4],BRAVO[t,4])
heatmapcont<-heatmapcont+1
#########Control4-protein-levels-same########
# BRAVO+WOX5
abw<-ass[3]
dbw<-diss[3]
# WOX5PLT3+BRAVO
awpb1<-ass[5]
dwpb1<-diss[5]
# BRAVOPLT3+WOX5
awpb2<-ass[6]
dwpb2<-diss[6]
# WOX5+PLT3
awp<-ass[1]
dwp<-diss[1]
# BRAVO+PLT3
abp<-ass[2]
dbp<-diss[2]


BRAVO=matrix(ncol=1,nrow=time+1,0)
PLT3=matrix(ncol=1,nrow=time+1,0)
WOX5=matrix(ncol=1,nrow=time+1,0)
WOX5PLT3=matrix(ncol=1,nrow=time+1,0)
BRAVOPLT3=matrix(ncol=1,nrow=time+1,0)
BRAVOWOX5=matrix(ncol=1,nrow=time+1,0)
WOX5PLT3BRAVO=matrix(ncol=1,nrow=time+1,0) 

BRAVO[1]<-1
PLT3[1]<-1
WOX5[1]<-1
WOX5PLT3[1]<-0
BRAVOPLT3[1]<-0
BRAVOWOX5[1]<-0
WOX5PLT3BRAVO[1]<-0

for(t in 1:time){
  BRAVO[t+1]<-BRAVO[t]-abw*BRAVO[t]*WOX5[t]+dbw*BRAVOWOX5[t]-abp*BRAVO[t]*PLT3[t]+dbp*BRAVOPLT3[t]-awpb1*WOX5PLT3[t]*BRAVO[t]+dwpb1*WOX5PLT3BRAVO[t]
  PLT3[t+1]<-PLT3[t]-abp*BRAVO[t]*PLT3[t]+dbp*BRAVOPLT3[t]-awp*WOX5[t]*PLT3[t]+dwp*WOX5PLT3[t]
  WOX5[t+1]<-WOX5[t]-awp*WOX5[t]*PLT3[t]+dwp*WOX5PLT3[t]-abw*BRAVO[t]*WOX5[t]+dbw*BRAVOWOX5[t]-awpb2*WOX5[t]*BRAVOPLT3[t]+dwpb2*WOX5PLT3BRAVO[t]
  WOX5PLT3[t+1]<-WOX5PLT3[t]+awp*WOX5[t]*PLT3[t]-dwp*WOX5PLT3[t]-awpb1*WOX5PLT3[t]*BRAVO[t]+dwpb1*WOX5PLT3BRAVO[t]
  BRAVOPLT3[t+1]<-BRAVOPLT3[t]+abp*BRAVO[t]*PLT3[t]-dbp*BRAVOPLT3[t]-awpb2*WOX5[t]*BRAVOPLT3[t]+dwpb2*WOX5PLT3BRAVO[t]
  BRAVOWOX5[t+1]<-BRAVOWOX5[t]+abw*BRAVO[t]*WOX5[t]-dbw*BRAVOWOX5[t]
  WOX5PLT3BRAVO[t+1]<-WOX5PLT3BRAVO[t]+awpb1*WOX5PLT3[t]*BRAVO[t]-dwpb1*WOX5PLT3BRAVO[t]+awpb2*WOX5[t]*BRAVOPLT3[t]-dwpb2*WOX5PLT3BRAVO[t]
}

###heatmap
heatmapdata[,1+heatmapcont*4]<-c(WOX5PLT3[t],BRAVOPLT3[t],BRAVOWOX5[t],WOX5PLT3BRAVO[t],
                                 WOX5[t],PLT3[t],BRAVO[t])

###########Control5-Ass==diss-and protein levels same###############
#All same ass diss
abw<-0.1
dbw<-0.1
# WOX5PLT3+BRAVO
awpb1<-0.1
dwpb1<-0.1
# BRAVOPLT3+WOX5
awpb2<-0.1
dwpb2<-0.1
# WOX5+PLT3
awp<-0.1
dwp<-0.1
# BRAVO+PLT3
abp<-0.1
dbp<-0.1

BRAVO=matrix(ncol=1,nrow=time+1,0)
PLT3=matrix(ncol=1,nrow=time+1,0)
WOX5=matrix(ncol=1,nrow=time+1,0)
WOX5PLT3=matrix(ncol=1,nrow=time+1,0)
BRAVOPLT3=matrix(ncol=1,nrow=time+1,0)
BRAVOWOX5=matrix(ncol=1,nrow=time+1,0)
WOX5PLT3BRAVO=matrix(ncol=1,nrow=time+1,0) 

BRAVO[1]<-1
PLT3[1]<-1
WOX5[1]<-1
WOX5PLT3[1]<-0
BRAVOPLT3[1]<-0
BRAVOWOX5[1]<-0
WOX5PLT3BRAVO[1]<-0

for(t in 1:time){
  BRAVO[t+1]<-BRAVO[t]-abw*BRAVO[t]*WOX5[t]+dbw*BRAVOWOX5[t]-abp*BRAVO[t]*PLT3[t]+dbp*BRAVOPLT3[t]-awpb1*WOX5PLT3[t]*BRAVO[t]+dwpb1*WOX5PLT3BRAVO[t]
  PLT3[t+1]<-PLT3[t]-abp*BRAVO[t]*PLT3[t]+dbp*BRAVOPLT3[t]-awp*WOX5[t]*PLT3[t]+dwp*WOX5PLT3[t]
  WOX5[t+1]<-WOX5[t]-awp*WOX5[t]*PLT3[t]+dwp*WOX5PLT3[t]-abw*BRAVO[t]*WOX5[t]+dbw*BRAVOWOX5[t]-awpb2*WOX5[t]*BRAVOPLT3[t]+dwpb2*WOX5PLT3BRAVO[t]
  WOX5PLT3[t+1]<-WOX5PLT3[t]+awp*WOX5[t]*PLT3[t]-dwp*WOX5PLT3[t]-awpb1*WOX5PLT3[t]*BRAVO[t]+dwpb1*WOX5PLT3BRAVO[t]
  BRAVOPLT3[t+1]<-BRAVOPLT3[t]+abp*BRAVO[t]*PLT3[t]-dbp*BRAVOPLT3[t]-awpb2*WOX5[t]*BRAVOPLT3[t]+dwpb2*WOX5PLT3BRAVO[t]
  BRAVOWOX5[t+1]<-BRAVOWOX5[t]+abw*BRAVO[t]*WOX5[t]-dbw*BRAVOWOX5[t]
  WOX5PLT3BRAVO[t+1]<-WOX5PLT3BRAVO[t]+awpb1*WOX5PLT3[t]*BRAVO[t]-dwpb1*WOX5PLT3BRAVO[t]+awpb2*WOX5[t]*BRAVOPLT3[t]-dwpb2*WOX5PLT3BRAVO[t]
  
}

###heatmap
heatmapdata[,2+heatmapcont*4]<-c(WOX5PLT3[t],BRAVOPLT3[t],BRAVOWOX5[t],WOX5PLT3BRAVO[t],
                                 WOX5[t],PLT3[t],BRAVO[t])

###########Control6-Ass>diss-and protein levels same###############
#All same ass diss
abw<-0.1
dbw<-0.05
# WOX5PLT3+BRAVO
awpb1<-0.1
dwpb1<-0.05
# BRAVOPLT3+WOX5
awpb2<-0.1
dwpb2<-0.05
# WOX5+PLT3
awp<-0.1
dwp<-0.05
# BRAVO+PLT3
abp<-0.1
dbp<-0.05
BRAVO=matrix(ncol=1,nrow=time+1,0)
PLT3=matrix(ncol=1,nrow=time+1,0)
WOX5=matrix(ncol=1,nrow=time+1,0)
WOX5PLT3=matrix(ncol=1,nrow=time+1,0)
BRAVOPLT3=matrix(ncol=1,nrow=time+1,0)
BRAVOWOX5=matrix(ncol=1,nrow=time+1,0)
WOX5PLT3BRAVO=matrix(ncol=1,nrow=time+1,0) 

BRAVO[1]<-1
PLT3[1]<-1
WOX5[1]<-1
WOX5PLT3[1]<-0
BRAVOPLT3[1]<-0
BRAVOWOX5[1]<-0
WOX5PLT3BRAVO[1]<-0

for(t in 1:time){
  BRAVO[t+1]<-BRAVO[t]-abw*BRAVO[t]*WOX5[t]+dbw*BRAVOWOX5[t]-abp*BRAVO[t]*PLT3[t]+dbp*BRAVOPLT3[t]-awpb1*WOX5PLT3[t]*BRAVO[t]+dwpb1*WOX5PLT3BRAVO[t]
  PLT3[t+1]<-PLT3[t]-abp*BRAVO[t]*PLT3[t]+dbp*BRAVOPLT3[t]-awp*WOX5[t]*PLT3[t]+dwp*WOX5PLT3[t]
  WOX5[t+1]<-WOX5[t]-awp*WOX5[t]*PLT3[t]+dwp*WOX5PLT3[t]-abw*BRAVO[t]*WOX5[t]+dbw*BRAVOWOX5[t]-awpb2*WOX5[t]*BRAVOPLT3[t]+dwpb2*WOX5PLT3BRAVO[t]
  WOX5PLT3[t+1]<-WOX5PLT3[t]+awp*WOX5[t]*PLT3[t]-dwp*WOX5PLT3[t]-awpb1*WOX5PLT3[t]*BRAVO[t]+dwpb1*WOX5PLT3BRAVO[t]
  BRAVOPLT3[t+1]<-BRAVOPLT3[t]+abp*BRAVO[t]*PLT3[t]-dbp*BRAVOPLT3[t]-awpb2*WOX5[t]*BRAVOPLT3[t]+dwpb2*WOX5PLT3BRAVO[t]
  BRAVOWOX5[t+1]<-BRAVOWOX5[t]+abw*BRAVO[t]*WOX5[t]-dbw*BRAVOWOX5[t]
  WOX5PLT3BRAVO[t+1]<-WOX5PLT3BRAVO[t]+awpb1*WOX5PLT3[t]*BRAVO[t]-dwpb1*WOX5PLT3BRAVO[t]+awpb2*WOX5[t]*BRAVOPLT3[t]-dwpb2*WOX5PLT3BRAVO[t]
}

###heatmap
heatmapdata[,3+heatmapcont*4]<-c(WOX5PLT3[t],BRAVOPLT3[t],BRAVOWOX5[t],WOX5PLT3BRAVO[t],
                                 WOX5[t],PLT3[t],BRAVO[t])

###########Control7-Ass<diss-and protein levels same###############
#All same ass diss
abw<-0.05
dbw<-0.1
# WOX5PLT3+BRAVO
awpb1<-0.05
dwpb1<-0.1
# BRAVOPLT3+WOX5
awpb2<-0.05
dwpb2<-0.1
# WOX5+PLT3
awp<-0.05
dwp<-0.1
# BRAVO+PLT3
abp<-0.05
dbp<-0.1

BRAVO=matrix(ncol=1,nrow=time+1,0)
PLT3=matrix(ncol=1,nrow=time+1,0)
WOX5=matrix(ncol=1,nrow=time+1,0)
WOX5PLT3=matrix(ncol=1,nrow=time+1,0)
BRAVOPLT3=matrix(ncol=1,nrow=time+1,0)
BRAVOWOX5=matrix(ncol=1,nrow=time+1,0)
WOX5PLT3BRAVO=matrix(ncol=1,nrow=time+1,0) 

BRAVO[1]<-1
PLT3[1]<-1
WOX5[1]<-1
WOX5PLT3[1]<-0
BRAVOPLT3[1]<-0
BRAVOWOX5[1]<-0
WOX5PLT3BRAVO[1]<-0

for(t in 1:time){
  BRAVO[t+1]<-BRAVO[t]-abw*BRAVO[t]*WOX5[t]+dbw*BRAVOWOX5[t]-abp*BRAVO[t]*PLT3[t]+dbp*BRAVOPLT3[t]-awpb1*WOX5PLT3[t]*BRAVO[t]+dwpb1*WOX5PLT3BRAVO[t]
  PLT3[t+1]<-PLT3[t]-abp*BRAVO[t]*PLT3[t]+dbp*BRAVOPLT3[t]-awp*WOX5[t]*PLT3[t]+dwp*WOX5PLT3[t]
  WOX5[t+1]<-WOX5[t]-awp*WOX5[t]*PLT3[t]+dwp*WOX5PLT3[t]-abw*BRAVO[t]*WOX5[t]+dbw*BRAVOWOX5[t]-awpb2*WOX5[t]*BRAVOPLT3[t]+dwpb2*WOX5PLT3BRAVO[t]
  WOX5PLT3[t+1]<-WOX5PLT3[t]+awp*WOX5[t]*PLT3[t]-dwp*WOX5PLT3[t]-awpb1*WOX5PLT3[t]*BRAVO[t]+dwpb1*WOX5PLT3BRAVO[t]
  BRAVOPLT3[t+1]<-BRAVOPLT3[t]+abp*BRAVO[t]*PLT3[t]-dbp*BRAVOPLT3[t]-awpb2*WOX5[t]*BRAVOPLT3[t]+dwpb2*WOX5PLT3BRAVO[t]
  BRAVOWOX5[t+1]<-BRAVOWOX5[t]+abw*BRAVO[t]*WOX5[t]-dbw*BRAVOWOX5[t]
  WOX5PLT3BRAVO[t+1]<-WOX5PLT3BRAVO[t]+awpb1*WOX5PLT3[t]*BRAVO[t]-dwpb1*WOX5PLT3BRAVO[t]+awpb2*WOX5[t]*BRAVOPLT3[t]-dwpb2*WOX5PLT3BRAVO[t]
}

###heatmap
heatmapdata[,4+heatmapcont*4]<-c(WOX5PLT3[t],BRAVOPLT3[t],BRAVOWOX5[t],WOX5PLT3BRAVO[t],
                                 WOX5[t],PLT3[t],BRAVO[t])

###########Plots###############
###Figure Controls, reordering results per cell type
heatmapdata2=matrix(nrow=7,ncol=20,0)
heatmapdata2[,1]<-heatmapdata[,1]
heatmapdata2[,2]<-heatmapdata[,5]
heatmapdata2[,3]<-heatmapdata[,9]
heatmapdata2[,4]<-heatmapdata[,13]
heatmapdata2[,5]<-heatmapdata[,2]
heatmapdata2[,6]<-heatmapdata[,6]
heatmapdata2[,7]<-heatmapdata[,10]
heatmapdata2[,8]<-heatmapdata[,14]
heatmapdata2[,9]<-heatmapdata[,3]
heatmapdata2[,10]<-heatmapdata[,7]
heatmapdata2[,11]<-heatmapdata[,11]
heatmapdata2[,12]<-heatmapdata[,15]
heatmapdata2[,13]<-heatmapdata[,4]
heatmapdata2[,14]<-heatmapdata[,8]
heatmapdata2[,15]<-heatmapdata[,12]
heatmapdata2[,16]<-heatmapdata[,16]
heatmapdata2[,17]<-heatmapdata[,17]
heatmapdata2[,18]<-heatmapdata[,18]
heatmapdata2[,19]<-heatmapdata[,19]
heatmapdata2[,20]<-heatmapdata[,20]
heatmapcont<-4*4+4
for(i in 1:(4*4+4)){
  heatmapdata[,i]<-heatmapdata2[,heatmapcont]
  heatmapcont<-heatmapcont-1
}
x<-c("WOX5-PLT3","BRAVO-PLT3","BRAVO-WOX5","WOX5-PLT3-BRAVO","WOX5","PLT3","BRAVO")
y<-c("Control 7","Control 6","Control 5","Control 4",
     "CC: Control 3","CC: Control 2","CC: Control 1","CC: Wt",
     "CSC: Control 3","CSC: Control 2","CSC: Control 1","CSC: Wt",
     "QC: Control 3","QC: Control 2","QC: Control 1","QC: Wt",
     "SI: Control 3","SI: Control 2","SI: Control 1","SI: Wt")

data<-expand.grid(X=x,Y=y)
data$Z<-c(heatmapdata)
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()+#scale_fill_viridis(discrete=FALSE)+
  scale_fill_gradientn(colours = colorspace::diverge_hcl(7))+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

# ggsave(paste("heatmap-controls.pdf"), width = 1.5*3, height = 2.3*3)
