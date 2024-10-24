library(viridis)
library(fmsb)
library(ggplot2)

rates<-read.table("BindingRates.csv",sep=",",header=F)
ratesOptions<-readRDS("bindingRates.rds")

replicates<-100
cellsComplexes=matrix(ncol=5,nrow=4*replicates)
for(tests in 1:replicates){
  
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
  PLT3values<-c(0.711044,0.68113,0.790596,0.294725)
  WOX5values<-c(0.524895,0.764438,0.327142,0.053709)
  BRAVOvalues<-c(1.001838,0.486001,0.220858,0.033101)
  
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
  
  ########## Saving end state
  cellsComplexes[tests,1]<-WOX5PLT3[t,4]
  cellsComplexes[tests,3]<-BRAVOPLT3[t,4]
  cellsComplexes[tests,4]<-BRAVOWOX5[t,4]
  cellsComplexes[tests,5]<-WOX5PLT3BRAVO[t,4]
  cellsComplexes[tests,2]<-"CC"
  cellsComplexes[tests,2]<-"CC"
  cellsComplexes[tests,2]<-"CC"
  cellsComplexes[tests,2]<-"CC"
  
  cellsComplexes[tests+replicates,1]<-WOX5PLT3[t,3]
  cellsComplexes[tests+replicates,3]<-BRAVOPLT3[t,3]
  cellsComplexes[tests+replicates,4]<-BRAVOWOX5[t,3]
  cellsComplexes[tests+replicates,5]<-WOX5PLT3BRAVO[t,3]
  cellsComplexes[tests+replicates,2]<-"CSC"
  cellsComplexes[tests+replicates,2]<-"CSC"
  cellsComplexes[tests+replicates,2]<-"CSC"
  cellsComplexes[tests+replicates,2]<-"CSC"
  
  cellsComplexes[tests+2*replicates,1]<-WOX5PLT3[t,2]
  cellsComplexes[tests+2*replicates,3]<-BRAVOPLT3[t,2]
  cellsComplexes[tests+2*replicates,4]<-BRAVOWOX5[t,2]
  cellsComplexes[tests+2*replicates,5]<-WOX5PLT3BRAVO[t,2]
  cellsComplexes[tests+2*replicates,2]<-"QC"
  cellsComplexes[tests+2*replicates,2]<-"QC"
  cellsComplexes[tests+2*replicates,2]<-"QC"
  cellsComplexes[tests+2*replicates,2]<-"QC"
  
  cellsComplexes[tests+3*replicates,1]<-WOX5PLT3[t,1]
  cellsComplexes[tests+3*replicates,3]<-BRAVOPLT3[t,1]
  cellsComplexes[tests+3*replicates,4]<-BRAVOWOX5[t,1]
  cellsComplexes[tests+3*replicates,5]<-WOX5PLT3BRAVO[t,1]
  cellsComplexes[tests+3*replicates,2]<-"SI"
  cellsComplexes[tests+3*replicates,2]<-"SI"
  cellsComplexes[tests+3*replicates,2]<-"SI"
  cellsComplexes[tests+3*replicates,2]<-"SI"
  
}
write.table(cellsComplexes,file="3-test.csv",sep=",",col.names = c("WOX5PLT3","CellType","BRAVOPLT3","BRAVOWOX5","WOX5PLT3BRAVO"),row.names = FALSE)
a<-read.csv(paste("3-test.csv"))
ggplot(a, aes(x=CellType, y=WOX5PLT3)) + geom_violin(trim=FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.2))+ stat_summary(fun.y=median, geom="point", size=2, color="red")
ggsave(paste("2-Test-parameters-WOX5-PLT3.pdf"),height = 3,width = 7)

ggplot(a, aes(x=CellType, y=BRAVOPLT3)) + geom_violin(trim=FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.2))+ stat_summary(fun.y=median, geom="point", size=2, color="red")
ggsave(paste("2-Test-parameters-BRAVO-PLT3.pdf"),height = 3,width = 7)# dev.off()

ggplot(a, aes(x=CellType, y=BRAVOWOX5)) + geom_violin(trim=FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.2))+ stat_summary(fun.y=median, geom="point", size=2, color="red")
ggsave(paste("2-Test-parameters-BRAVO-WOX5.pdf"),height = 3,width = 7)# dev.off()

ggplot(a, aes(x=CellType, y=WOX5PLT3BRAVO)) + geom_violin(trim=FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.2))+ stat_summary(fun.y=median, geom="point", size=2, color="red")
ggsave(paste("2-Test-parameters-WOX5-PLT3-BRAVO.pdf"),height = 3,width = 7)# dev.off()

