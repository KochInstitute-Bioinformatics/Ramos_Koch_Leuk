pdf(file='invitro.poolB_hEGFRv3.CAR_ET.1.10_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolB_hEGFRv3.CAR_ET.1.10_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Dph6","Mcat","Pcyt2","Mthfd2","Elovl6","Taz","Alad","Coq6","Dnajc11","Vkorc1l1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='2,3,4_vs_0 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(565.7569502298023,297.98765860642584,293.29786806899295,340.52800964048504),c(469.37586841893307,167.338414893098,180.7827376459741,161.0605451002294),c(286.2518129782816,102.01379303643407,87.23083145155394,170.26400482024252),c(537.8064365046503,232.6630367497619,223.76604589746444,202.4761138402884))
targetgene="Dph6"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(417.33008424106373,211.18644874209156,270.54199899467454,248.49341244035395),c(437.57011142134627,181.65614023154487,168.14058816024163,168.42331287623992),c(212.03837998391228,48.32232301725824,51.83281289150306,69.02594790009832),c(880.9230877513446,358.8379912948251,343.86646601192274,322.1210902004588))
targetgene="Mcat"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(329.62329979317275,191.4995764017271,166.8763732116684,191.4319621762727),c(578.2864908652153,405.37059864477743,395.6992789034258,387.4656542125519),c(313.23851588532494,175.39213539597435,174.46166290310788,144.49431760420583),c(320.94900243019447,161.07441005752747,184.57538249169383,173.02504273624646))
targetgene="Pcyt2"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(270.8308398885425,133.33381721428663,150.4415788802162,167.5029669042386),c(327.69567815695535,168.2332727267509,198.48174692599954,193.2726541202753),c(302.6365968861293,190.60471856807416,115.04356032016534,150.0163934362137),c(145.53543353441253,41.163460348034796,42.98330825149034,67.1852559560957))
targetgene="Mthfd2"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(333.4785430656075,205.81730174017397,142.85628918877674,235.6085688323356),c(369.1395433356291,238.9270415853324,222.5018309488912,240.21029869234215),c(224.56792061932526,97.5395038681694,82.17397165726095,91.11425122812979),c(248.6631910720426,93.96007253355769,116.30777526873858,106.76013275215207))
targetgene="Elovl6"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(403.836732787542,304.2516634419963,225.0302608460377,258.6172181323684),c(645.7532481328237,475.16950966970603,482.9301103549798,422.4388011486017),c(179.26881216821675,110.96237137296336,80.90975670868771,69.02594790009832),c(200.47265016660796,76.95777369415201,92.28769124584691,39.57487679605637))
targetgene="Taz"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(184.08786625876022,91.2754990325989,80.90975670868771,100.31771094814289),c(207.21932589336882,112.75208704026923,115.04356032016534,102.15840289214552),c(283.3603805239555,149.44125822003937,108.72248557729911,96.63632706013765),c(128.18683880845606,40.26860251438187,60.68231753151578,40.49522276805768))
targetgene="Alad"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(213.96600162012967,99.32921953547527,67.00339227438201,42.3359147120603),c(316.12994833965104,165.54869922579212,125.1572799087513,118.72463038816912),c(223.6041098012166,212.97616440939743,219.9734010517447,186.83023231626612),c(645.7532481328237,362.41742262943677,312.26109229759163,363.53665894051784))
targetgene="Coq6"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(310.3470834309989,180.76128239789193,156.76265362308243,158.2995071842255),c(392.2710029702377,198.65843907095055,188.36802733741357,181.30815648425826),c(255.40986679880342,88.59092553164011,99.87298093728639,125.16705219217829),c(124.3315955360213,68.90405319127564,51.83281289150306,43.256260684061616))
targetgene="Dnajc11"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(266.01178579799904,101.11893520278113,109.98670052587235,105.83978678015076),c(266.01178579799904,156.60012088926283,140.32785929163023,177.626772596253),c(218.78505571067313,140.49267988351008,122.62885001160481,113.20255455616125),c(130.11446044467345,71.58862669223443,29.076943817184645,48.77833651606948))
targetgene="Vkorc1l1"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Ralgapb","Cnot6l","Rasa3","Dcp2","Dazap1","Edc4","Sms","Stk10","Dcp1a","Zfp36l2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='2,3,4_vs_0 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(307.4556509766728,532.4404110234935,458.9100263320881,518.1547822367381),c(338.29759715615097,619.2416208878278,624.5221845951833,632.2776827649006),c(462.62919269217224,658.6153655685567,690.261361920992,638.7201045689098),c(1002.3632508330398,1681.4378694338561,1754.7303486196647,1656.6227496023598))
targetgene="Ralgapb"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(821.1668170286057,1279.6467021236904,1225.0242851674748,1263.6350195578),c(692.9799782201496,893.0681179856246,819.211286675463,826.4706828571773),c(731.5324109444973,881.4349661481364,1065.733201647246,1086.9285929335483),c(847.1897091175405,1234.903810441044,1241.459079498927,1244.3077541457724))
targetgene="Cnot6l"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(209.1469475295862,283.66993326797893,340.073821166203,290.82932715241424),c(903.0907365678446,1212.5323645997207,1219.9674253731819,1090.6099768215536),c(385.5243272434769,569.1295822032637,510.74283922359115,594.5434979128469),c(389.37957051591167,476.9592253370119,520.8565588121771,493.3054409927027))
targetgene="Rasa3"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(586.9607882281936,988.817906186488,817.9470717268898,981.0888061533975),c(176.37737971389066,306.04137910930217,370.41497993196094,213.52026550430415),c(554.191220412498,728.4142765934853,716.8098758410301,874.3286734012454),c(482.8692198724548,465.3260734995238,420.9835778748907,431.64226086861487))
targetgene="Dcp2"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(976.3403587441052,1166.8946150834213,1024.0141083443289,1229.5822185937513),c(263.12035334367295,338.25626112080766,460.17424128066136,410.4743035125847),c(246.7355694358252,369.57628529866025,371.67919488053417,336.84662575247984),c(379.74146233482475,483.2232301725824,391.90663405770607,474.89852155267647))
targetgene="Dazap1"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(816.3477629380623,1165.1048994161154,1286.9708176475638,1114.5389720935875),c(246.7355694358252,401.79116731016575,379.26448457197364,309.23624659244047),c(242.88032616339044,324.8333936160137,288.24100827469994,377.3418485205375),c(553.2274095943893,701.5685415838974,635.9001191323424,664.4897917849465))
targetgene="Edc4"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(351.7909486096726,544.0735628609816,475.3448206635403,521.8361661247433),c(608.1646262265848,1003.135631524935,938.0474918413481,839.3555264651956))
targetgene="Sms"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(792.252492485345,970.9207495134294,965.8602207099595,1042.7519862774852),c(985.0146561070834,1443.4056856821767,1231.3453599103411,1376.8375741139612),c(252.51843434447736,308.725952610261,265.4851392003815,351.5721613045008),c(819.2391953923883,983.4487591845705,1045.505762470074,1120.9813938975967))
targetgene="Stk10"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1077.5404946455178,1583.8983655656868,1653.593152733805,1652.0210197423532),c(678.5228159485193,670.2485174060448,493.04382994356575,529.1989339007538),c(1083.32335955417,1378.9759216591658,1359.0310697162388,1415.4921049380162),c(473.23111169136786,603.134179882075,726.9235954296162,572.4551945848154))
targetgene="Dcp1a"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(309.3832726128902,466.22093133317674,452.5889515892219,505.26993862871973),c(494.4349496897591,612.9776160522573,596.7094557265718,686.578095112978),c(266.97559661610774,318.5693887804432,308.4684474518719,431.64226086861487),c(665.9932753131063,1202.6889284295385,1122.6228743330419,1072.2030573815273))
targetgene="Zfp36l2"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("invitro.poolB_hEGFRv3.CAR_ET.1.10_v_input_summary.Rnw");
library(tools);

texi2dvi("invitro.poolB_hEGFRv3.CAR_ET.1.10_v_input_summary.tex",pdf=TRUE);

