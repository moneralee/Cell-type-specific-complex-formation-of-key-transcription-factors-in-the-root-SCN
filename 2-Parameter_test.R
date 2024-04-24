library(viridis)
library(fmsb)
library(ggplot2)
zoomplot<-0.2

rates<-read.table("BindingRates.csv",sep=",",header=F)
ratesOptions<-readRDS("bindingRates.rds")

pdf(paste("2-ParameterTestWt.pdf"))
par(mfrow=c(2,2))
for(tests in 1:100){
  
  ass<-0
  diss<-0
  helper<-0
  
  for(i in 1:length(ratesOptions)){
    helper[i]<-1+3*(i-1) # Because in table, there is a space in between ass/diss per complex
    sampleassdiss<-sample(1:ratesOptions[i],1) # Random selection of random rates to be used in the simulation
    ass[i]<-as.double(rates[helper[i],sampleassdiss+1]) 
    diss[i]<-as.double(rates[helper[i]+1,sampleassdiss+1]) 
  }
  
  #############################Protein complex simulation#######################
  time<-2000
  heatmapdata=matrix(nrow=8,ncol=14+8+2,0)
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
  
  #Check conservation of matter
  cellid<-2
  BRAVOvalues[cellid]-sum(BRAVO[t,cellid],BRAVOPLT3[t,cellid],BRAVOWOX5[t,cellid],WOX5PLT3BRAVO[t,cellid])
  PLT3values[cellid]-sum(PLT3[t,cellid],BRAVOPLT3[t,cellid],WOX5PLT3[t,cellid],WOX5PLT3BRAVO[t,cellid])
  WOX5values[cellid]-sum(WOX5[t,cellid],BRAVOWOX5[t,cellid],WOX5PLT3[t,cellid],WOX5PLT3BRAVO[t,cellid])
  
  ########## Radar plot protein complexes
  data<-as.data.frame(matrix(c(WOX5PLT3[t,],BRAVOPLT3[t,],BRAVOWOX5[t,],WOX5PLT3BRAVO[t,]),ncol=4))
  rownames(data)<-c("Stele initials","QC","CSC","CC")
  colnames(data)<-c("WOX5-PLT3","BRAVO-PLT3","BRAVO-WOX5","WOX5-PLT3-BRAVO")
  data <- rbind(rep(zoomplot,4) , rep(0,4) , data)
  radarchart(data,pfcol=adjustcolor(viridis(4),alpha.f = 0.2),pcol=viridis(4),caxislabels = seq(0, zoomplot, length.out = 5),calcex = 1,vlcex = 0.8)
}


dev.off()
