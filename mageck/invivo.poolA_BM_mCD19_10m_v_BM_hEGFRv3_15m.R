pdf(file='invivo.poolA_BM_mCD19_10m_v_BM_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolA_BM_mCD19_10m_v_BM_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Cer1","Exosc6","Ptgir","Map3k7","Gclc","Atp5a1","Pabpc1l2a_Pabpc1l2b","Atp6ap2","Stub1","Akt2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='6,7,8_vs_1,2,3,4,5 neg.'


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
targetmat=list(c(1.6847076630845468,4.126502478160202,14.046380768674057,0.0,19.594483431256595,0.0,11.593119892779795,0.0),c(30.184345630264797,35.95952159539605,10.534785576505543,7.807815600634233,12.246552144535372,0.0,17.003242509410367,0.8931477942843442),c(8.423538315422734,2.063251239080101,3.5115951921685142,0.0,17.14517300234952,0.0,0.0,2.2328694857108604),c(16.14511510456024,153.5648422229618,31.604356729516628,0.0,24.493104289070743,0.0,0.0,23.66841654853512))
targetgene="Cer1"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(15.021976662503876,2.063251239080101,7.0231903843370285,23.4234468019027,56.33413986486271,0.0,0.0,2.2328694857108604),c(63.878498891955736,10.905756549423392,14.046380768674057,109.30941840887925,51.43551900704856,0.0,2.318623978555959,0.0),c(66.12477577606846,5.89500354022886,35.115951921685145,23.4234468019027,48.986208578141486,24.541522192521775,1.545749319037306,60.28747611419323),c(26.25336108306752,14.737508850572151,21.069571153011086,0.0,4.898620857814149,0.0,0.0,2.679443382853033))
targetgene="Exosc6"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(48.15456070316663,217.82038081145637,143.97540287890908,15.615631201268465,122.46552144535372,0.0,0.0,138.8844820112155),c(57.70123746064573,19.158761505743797,38.62754711385366,70.2703404057081,75.92862329611931,107.71001406717889,42.50810627352591,112.0900481826852),c(94.20323682747758,100.21506018389063,35.115951921685145,124.92504961014772,48.986208578141486,2.726835799169086,33.23361035930208,65.19978898275713),c(120.4565979105451,220.47313240455938,196.6493307614368,382.5829644310774,200.8434551703801,489.46702595085094,210.99478204859227,410.8479853707983))
targetgene="Ptgir"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(15.162368967760923,133.22708000917225,38.62754711385366,78.07815600634233,24.493104289070743,0.0,0.0,4.912312868563893),c(0.8423538315422734,77.22454637699806,7.0231903843370285,0.0,2.4493104289070744,0.0,0.0,0.0),c(122.98365940517192,22.106263275858225,52.673927882527714,0.0,63.68207115158393,4.090253698753629,0.0,404.1493769136658),c(27.23610721986684,6.484503894251747,63.208713459033255,117.1172340095135,58.78345029376979,6.817089497922715,13.1388692118171,4.019165074279549))
targetgene="Map3k7"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(9.125499841707962,10.611006372411948,14.046380768674057,7.807815600634233,36.739656433606115,0.0,0.0,6.25203455999041),c(87.88558309091053,103.4573121310165,35.115951921685145,70.2703404057081,48.986208578141486,0.0,0.0,3.126017279995205),c(20.778061178042744,13.263757965514936,59.69711826686474,31.23126240253693,29.391725146884895,0.0,44.82673025208187,156.7474378969024),c(40.57337621928617,140.3010842574469,24.5811663451796,7.807815600634233,29.391725146884895,0.0,0.0,1.3397216914265164))
targetgene="Gclc"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(5.194515294510686,2.063251239080101,10.534785576505543,39.03907800317116,12.246552144535372,324.49346010112123,24.731989104596895,12.057495222838646),c(23.305122672669565,30.064518055167188,7.0231903843370285,7.807815600634233,9.797241715628298,0.0,2.318623978555959,0.0),c(27.376499525123887,0.294750177011443,7.0231903843370285,46.8468936038054,12.246552144535372,0.0,0.0,4.019165074279549),c(74.68870639674824,24.46426469194977,49.1623326903592,39.03907800317116,4.898620857814149,0.0,0.0,0.0))
targetgene="Atp5a1"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(19.514530430729334,1.473750885057215,24.5811663451796,0.0,9.797241715628298,0.0,0.0,0.0),c(11.792953641591827,17.979760797698024,49.1623326903592,0.0,17.14517300234952,144.52229735596154,19.321866487966325,123.2543956112395),c(30.184345630264797,0.294750177011443,0.0,15.615631201268465,36.739656433606115,0.0,0.0,15.183512502833851),c(42.53886849288481,25.053765045972657,45.650737498190686,187.3875744152216,56.33413986486271,0.0,0.0,18.309529782829056))
targetgene="Pabpc1l2a_Pabpc1l2b"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(39.309845471972764,13.558508142526378,28.092761537348114,0.0,26.94241471797782,0.0,0.772874659518653,0.4465738971421721),c(65.56320655504028,58.36053504826572,73.7434990355388,0.0,19.594483431256595,2.726835799169086,3.091498638074612,20.98897316568209),c(8.56393062067978,1.473750885057215,10.534785576505543,54.65470920443963,0.0,0.0,0.0,0.0),c(0.0,0.0,0.0,0.0,4.898620857814149,5.453671598338172,1.545749319037306,0.0))
targetgene="Atp6ap2"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(244.84418036828748,244.6426469194977,147.4869980710776,203.00320561649005,169.00241959458813,27.26835799169086,4.637247957111918,100.03255295984656),c(32.57101481963457,215.46237939536485,70.23190384337029,140.5406808114162,85.7258650117476,5.453671598338172,1.545749319037306,1.3397216914265164),c(157.80095110891924,22.40101345286967,115.88264134156097,7.807815600634233,75.92862329611931,103.61976036842526,29.369237061708812,50.01627647992328),c(128.59935161545374,13.263757965514936,7.0231903843370285,46.8468936038054,58.78345029376979,1.363417899584543,0.0,2.2328694857108604))
targetgene="Stub1"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(86.62205234359712,138.5325831953782,112.37104614939246,148.34849641205042,137.16138401879616,289.0445947119231,435.90130796852026,519.8120162734883),c(165.1013509822856,66.90829018159756,35.115951921685145,93.6937872076108,110.21896930081834,0.0,0.0,16.076660297118195),c(26.534145693581614,240.81089461834895,203.67252114577383,70.2703404057081,46.53689814923441,184.0614164439133,136.79881473480157,186.66788900542795),c(72.44242951263551,71.62429301378066,73.7434990355388,132.73286521078197,117.56690058753958,117.25393936427069,110.52107631116738,91.99422281128746))
targetgene="Akt2"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetgenelist=c("Cd19","Trappc4","Tas2r135","Gpx6","Pklr","Ccl17","Dut","Casp1","Dio1","Gfap")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='6,7,8_vs_1,2,3,4,5 pos.'


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
targetmat=list(c(35.94043014580367,40.97027460459058,24.5811663451796,148.34849641205042,19.594483431256595,151.33938685388426,83.47046322801452,1734.4930165001965),c(179.98293533953242,77.51929655400951,49.1623326903592,15.615631201268465,46.53689814923441,2336.8982798879065,522.4632698346094,998.5392340098969),c(133.65347460470738,153.5648422229618,98.3246653807184,109.30941840887925,51.43551900704856,13958.672455946551,4283.271363052375,2594.1477684988777),c(122.5624824894008,188.05061293330064,214.20730672227936,117.1172340095135,195.94483431256594,6447.603247135304,1986.287874962938,5541.088915740072))
targetgene="Cd19"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(15.302761273017968,1.179000708045772,7.0231903843370285,0.0,4.898620857814149,0.0,0.0,2.2328694857108604),c(32.149837903863435,42.73877566665924,42.13914230602217,31.23126240253693,41.638277291420266,5.453671598338172,0.772874659518653,6.25203455999041),c(0.0,16.80076008965225,14.046380768674057,0.0,0.0,269.9567441177395,273.59762946960313,20.542399268539917),c(36.3616070615748,4.421252655171645,3.5115951921685142,23.4234468019027,14.695862573442447,0.0,0.0,3.126017279995205))
targetgene="Trappc4"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(267.16655690415774,395.84948772636795,354.67111440901994,554.3549076450305,274.3227680375923,2736.379724466178,488.4567848157887,1441.9871138720737),c(189.52961209701152,402.3339916206197,259.85804422047005,304.5048084247351,244.93104289070743,2051.9439388747373,440.5385559256322,290.71960703955403),c(261.691256999133,155.6280934620419,101.83626057288691,257.6579148209297,176.35035088130937,1614.286793108099,456.7689237755239,762.3016424216878),c(139.83073603601738,392.31248560223065,126.41742691806651,312.3126240253693,166.55310916568106,339.4910569965512,844.7520028538877,1007.9172858498824))
targetgene="Tas2r135"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(15.86433049404615,38.317523011487594,115.88264134156097,93.6937872076108,24.493104289070743,69.53431287881169,9.274495914223836,674.3265846846799),c(43.80239924019822,50.40228026895676,59.69711826686474,15.615631201268465,93.07379629846882,6731.194170248888,692.4956949287131,168.80493311974107),c(142.49818983590126,136.76408213330956,217.71890191444788,226.42665241839276,129.81345273207495,16.361014795014515,218.7235286437788,478.2806438392663),c(288.50618730322867,300.6451805516719,231.76528268312194,437.237673635517,176.35035088130937,3150.8587659398786,4154.974169572279,1084.727996158336))
targetgene="Gpx6"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(337.78388644845165,151.79634116089315,263.36963941263855,132.73286521078197,210.64069688600839,118.61735726385524,675.4924524193027,509.9873905363605),c(58.262806681673915,93.73055628963888,98.3246653807184,148.34849641205042,117.56690058753958,1355.2373921870358,280.55350140527105,149.1556816454855),c(61.21104509207187,86.65655204136425,59.69711826686474,163.96412761331888,56.33413986486271,159.51989425139152,1113.712384366379,255.44026916532243),c(83.11224471217098,68.38204106665478,126.41742691806651,85.88597160697655,107.76965887191128,458.10841426040645,305.28549050986794,85.74218825129705))
targetgene="Pklr"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(70.05576032326574,67.49779053562045,63.208713459033255,23.4234468019027,19.594483431256595,20.451268493768143,101.24658039694354,28.580729417099015),c(9.687069062736144,0.0,0.0,0.0,2.4493104289070744,32.72202959002903,2.318623978555959,112.0900481826852),c(27.376499525123887,39.79127389654481,42.13914230602217,0.0,41.638277291420266,699.4333824868705,63.37572208052954,549.7324673820139),c(167.48802017165536,292.68692577236294,294.9739961421552,242.04228361966122,213.09000731491548,203.1492670380969,390.30170305691973,212.56917503967392))
targetgene="Ccl17"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(21.90119962009911,267.9279109034017,80.76668941987583,31.23126240253693,4.898620857814149,20.451268493768143,6.955871935667877,25.00813823996164),c(31.167091767064118,1.768501062068658,7.0231903843370285,0.0,29.391725146884895,0.0,0.0,0.0),c(0.0,0.0,0.0,0.0,0.0,1.363417899584543,0.0,2.679443382853033),c(0.0,10.905756549423392,31.604356729516628,0.0,0.0,268.59332621815497,193.21866487966324,16.523234194260368))
targetgene="Dut"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(183.77352758147265,65.72928947355179,98.3246653807184,359.15951762917473,80.82724415393345,132.25153625970066,1068.1127794547783,265.26489490245024),c(73.42517564943483,127.33207646894338,59.69711826686474,78.07815600634233,80.82724415393345,259.0494009210632,262.777384236342,24.114990445677293),c(139.69034373076033,392.90198595625355,98.3246653807184,312.3126240253693,102.87103801409712,511.2817123442036,60.28422344245493,1642.052219791767),c(114.13894417397805,172.7236037287056,77.25509422770732,46.8468936038054,88.17517544065468,145.8857152555461,558.0155041724674,377.8015169822776))
targetgene="Casp1"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(96.44951371159031,322.7514438275301,150.9985932632461,406.0064112329801,144.5093153055174,1243.4371244211031,1238.9180792084007,129.5064301712299),c(14.320015136218649,72.50854354481498,77.25509422770732,0.0,19.594483431256595,6.817089497922715,602.8422344245494,542.5872850277391),c(253.40811098896725,221.94688328961658,133.44061730240355,0.0,169.00241959458813,357.21548969115025,152.25630792517464,101.37227465127307),c(13.61805360993342,54.8235329241284,31.604356729516628,117.1172340095135,26.94241471797782,12.270761096260888,24.731989104596895,21.435547062824263))
targetgene="Dio1"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(446.8687076331761,482.80078994474366,277.4160201813126,577.7783544469332,244.93104289070743,1308.8811836011612,316.10573574312906,1016.8487637927259),c(62.19379122887119,36.549021949418936,28.092761537348114,15.615631201268465,112.66827972972543,114.52710356510161,831.6131336420706,138.43790811407337),c(85.07773698576962,58.950035402288606,112.37104614939246,15.615631201268465,63.68207115158393,115.89052146468615,240.3640191103011,191.13362797684965),c(121.72012865785851,4.421252655171645,70.23190384337029,7.807815600634233,71.03000243830516,452.65474266206826,222.58790194137205,276.42924233100456))
targetgene="Gfap"
collabel=c("In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1823_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1825_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1828_Pool.A..1.8._BM_mCD19_10m_1360..1361..1362..1363..1406..1407..1408..1409")

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
Sweave("invivo.poolA_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolA_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.tex",pdf=TRUE);

