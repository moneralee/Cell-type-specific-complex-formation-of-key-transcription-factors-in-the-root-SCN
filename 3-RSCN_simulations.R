library(viridis)
library(fmsb)
library(ggplot2)
zoomplot<-0.2
rates<-read.table("BindingRates.csv",sep=",",header=F)
ratesOptions<-readRDS("bindingRates.rds")

ass<-0
diss<-0
helper<-0

for(i in 1:length(ratesOptions)){
  helper[i]<-1+3*(i-1) # Because in table, there is a space in between ass/diss per complex
  sampleassdiss<-sample(1:ratesOptions[i],1) # Random selection of random binding rates to be used in the simulation
  ass[i]<-as.double(rates[helper[i],sampleassdiss+1]) 
  diss[i]<-as.double(rates[helper[i]+1,sampleassdiss+1]) 
}

#############################Protein complex simulation#######################
time<-2000
heatmapdata=matrix(nrow=7,ncol=4,0)
heatmapcont<-0

###### Root stem cell niche protein levels (normalized)
PLT3values<-c(0.711044,0.68113,0.790596,0.294725)
WOX5values<-c(0.524895,0.764438,0.327142,0.053709)
BRAVOvalues<-c(1.001838,0.486001,0.220858,0.033101)

for(prion in 0:1){
  # BRAVO+WOX5
  abw<-ass[3]
  dbw<-diss[3]
  # WOX5PLT3+BRAVO
  awpb1<-ass[5]
  dwpb1<-diss[5]
  # BRAVOPLT3+WOX5
  awpb2<-ass[6]
  dwpb2<-diss[6]
  
  if (prion==0){
    # WOX5+PLT3
    awp<-ass[1]
    dwp<-diss[1]
    # BRAVO+PLT3
    abp<-ass[2]
    dbp<-diss[2]
  }else{
    # WOX5+PLT3dPrD
    awp<-0
    dwp<-0
    # BRAVO+PLT3dPrD
    abp<-ass[4]
    dbp<-diss[4]
  }  
  
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
  
  ########## Radar plot protein complexes
  data<-as.data.frame(matrix(c(WOX5PLT3[t,],BRAVOPLT3[t,],BRAVOWOX5[t,],WOX5PLT3BRAVO[t,]),ncol=4))
  rownames(data)<-c("Stele initials","QC","CSC","CC")
  colnames(data)<-c("WOX5-PLT3","BRAVO-PLT3","BRAVO-WOX5","WOX5-PLT3-BRAVO")
  data <- rbind(rep(zoomplot,4) , rep(0,4) , data)
  pdf(paste("3-Protein-Complexes-prion-simulation-",prion,"-zoom-",zoomplot,"-cells-radarplot.pdf"))
  radarchart(data,pfcol=adjustcolor(viridis(4),alpha.f = 0.2),pcol=viridis(4),caxislabels = seq(0, zoomplot, length.out = 5),calcex = 1,vlcex = 0.8)
  
  legend(1,1.25,
         legend=c("SI","QC","CSC","CC"),
         pch=c(15,16),
         col=c(viridis(4)),
         lty=c(1,2),cex = 0.7)
  
  dev.off()
  
  ########## Free protein barplots
  pdf(paste("3-Free_protein-prion-simulation-",prion,"-zoom-",zoomplot,"-cells-barplot.pdf"))
  par(mfrow=c(2,2))
  cellcolors<-viridis(4)
  for(i in 1:4){
    barplot(c(WOX5[time+1,i],PLT3[time+1,i],BRAVO[time+1,i]),col=cellcolors[i],names.arg = c("WOX5","PLT3","BRAVO"),ylim=c(0,0.3))
  }
  dev.off()
  
  ###heatmap
  heatmapdata[,1]<-c(WOX5PLT3[t,1],
                     BRAVOPLT3[t,1],
                     BRAVOWOX5[t,1],
                     WOX5PLT3BRAVO[t,1],
                     WOX5[t,1],
                     PLT3[t,1],
                     BRAVO[t,1])
  heatmapdata[,2]<-c(WOX5PLT3[t,2],
                     BRAVOPLT3[t,2],
                     BRAVOWOX5[t,2],
                     WOX5PLT3BRAVO[t,2],
                     WOX5[t,2],
                     PLT3[t,2],
                     BRAVO[t,2])
  heatmapdata[,3]<-c(WOX5PLT3[t,3],
                     BRAVOPLT3[t,3],
                     BRAVOWOX5[t,3],
                     WOX5PLT3BRAVO[t,3],
                     WOX5[t,3],
                     PLT3[t,3],
                     BRAVO[t,3])
  heatmapdata[,4]<-c(WOX5PLT3[t,4],
                     BRAVOPLT3[t,4],
                     BRAVOWOX5[t,4],
                     WOX5PLT3BRAVO[t,4],
                     WOX5[t,4],
                     PLT3[t,4],
                     BRAVO[t,4]) #max value in other plots for comparison
  
  x<-c("WOX5-PLT3","BRAVO-PLT3","BRAVO-WOX5","WOX5-PLT3-BRAVO","WOX5","PLT3","BRAVO")
  y<-c("CC: Wt","CSC: Wt","QC: Wt","SI: Wt")
  data<-expand.grid(X=x,Y=y)
  heatmapdata2=matrix(nrow=7,ncol=4,0)
  heatmapcont<-4
  for(i in 1:(4)){
    heatmapdata2[,i]<-heatmapdata[,heatmapcont]
    heatmapcont<-heatmapcont-1
  }
  data<-expand.grid(X=x,Y=y)
  data$Z<-c(heatmapdata2)
  ggplot(data, aes(X, Y, fill= Z)) + geom_tile()+scale_fill_gradientn(colours = colorspace::diverge_hcl(7),limits=range(0,0.43), oob = scales::squish)+
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  
  ggsave(paste("3-heatmap-prion",prion,".pdf"), width = 1.5*3, height = 1*3)
}