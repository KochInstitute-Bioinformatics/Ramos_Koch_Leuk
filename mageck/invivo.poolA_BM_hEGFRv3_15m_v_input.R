pdf(file='invivo.poolA_BM_hEGFRv3_15m_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolA_BM_hEGFRv3_15m_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Cox15","Cox17","Myh9","Cox5a","Cxcr4","Nme6","Aldoa","Hscb","H2-K1","Ugp2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4,5_vs_0 neg.'


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
targetmat=list(c(40.07697480916438,40.34012169011508,0.6432275763189099,0.0,0.0,5.331068417906823),c(29.51777186543245,7.940968836636827,1.6080689407972748,11.658987731086892,0.0,2.6655342089534115),c(61.43536258171307,79.88614649656648,4.502593034232369,3.8863292436956307,0.0,42.648547343254585),c(137.74960203868477,6.670413822774935,18.975213501407843,38.86329243695631,0.0,18.658739462673882))
targetgene="Cox15"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(44.6366306257759,5.558678185645779,0.6432275763189099,11.658987731086892,0.0,7.996602626860234),c(32.63753637153507,0.7940968836636828,0.6432275763189099,15.545316974782523,0.0,0.0),c(32.15757260136544,3.494026288120204,3.8593654579134595,3.8863292436956307,43.24354712509756,7.996602626860234),c(92.39302575765443,6.3527750693094625,7.718730915826919,31.090633949565046,0.0,18.658739462673882))
targetgene="Cox17"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(33.83744579695915,0.6352775069309462,2.2512965171161845,19.431646218478154,8.648709425019511,0.0),c(38.39710161357067,6.511594446042198,6.753889551348554,0.0,0.0,26.655342089534116),c(54.47588791425339,8.893885097033246,5.145820610551279,27.204304705869415,0.0,31.986410507440937),c(31.437626946110985,1.5881937673273656,0.9648413644783649,3.8863292436956307,0.0,0.0))
targetgene="Myh9"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(40.55693857933402,14.929021412877235,0.0,0.0,0.0,15.993205253720468),c(42.956757430182186,6.511594446042198,7.075503339508009,7.772658487391261,0.0,15.993205253720468),c(9.359293518307851,1.9058325207928386,36.98558563833732,15.545316974782523,8.648709425019511,15.993205253720468),c(39.357029153909934,0.0,1.9296827289567298,3.8863292436956307,0.0,0.0))
targetgene="Cox5a"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(81.11387715866805,84.80954717528132,1.6080689407972748,11.658987731086892,0.0,37.317478925347764),c(76.07425757188689,17.94658957079923,5.78904818687019,50.522280168043196,0.0,15.993205253720468),c(17.998641381361253,47.64581301982096,0.0,0.0,0.0,15.993205253720468),c(102.23228304613191,8.73506572030051,12.221323950059288,46.63595092434757,492.9764372261121,10.662136835813646))
targetgene="Cxcr4"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(127.43038098003765,32.71679160694373,11.256482585580924,38.86329243695631,34.594837700078045,74.63495785069553),c(70.79465610002092,4.446942548516623,9.326799856624193,27.204304705869415,0.0,10.662136835813646),c(51.836087178320405,9.687981980696929,3.2161378815945496,15.545316974782523,0.0,5.331068417906823),c(275.9791678475392,102.91495612281328,46.63399928312097,93.27190184869514,216.2177356254878,125.28010782081034))
targetgene="Nme6"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(18.23862326644607,0.15881937673273655,0.0,0.0,0.0,0.0),c(29.51777186543245,3.1763875346547312,5.78904818687019,38.86329243695631,25.946128275058534,2.6655342089534115),c(6.719492782374867,0.0,0.0,0.0,0.0,0.0),c(18.478605151530886,0.0,0.6432275763189099,7.772658487391261,0.0,10.662136835813646))
targetgene="Aldoa"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(40.55693857933402,0.47645813019820965,0.0,0.0,0.0,0.0),c(38.637083498655485,2.064651897525575,772.8379329471703,34.97696319326067,0.0,7.996602626860234),c(39.59701103899475,2.699929404456521,0.0,0.0,8.648709425019511,7.996602626860234),c(17.038713841021984,1.5881937673273656,0.0,0.0,0.0,0.0))
targetgene="Hscb"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(131.75005491156435,19.217144584661124,165.95271469027875,54.40860941173883,17.297418850039023,37.317478925347764),c(68.39483724917275,13.81728577574808,108.38384660973632,38.86329243695631,0.0,15.993205253720468),c(54.955851684423024,2.541110027723785,34.41267533306168,27.204304705869415,8.648709425019511,0.0),c(66.2350002834094,0.9529162603964193,20.90489623036457,7.772658487391261,0.0,2.6655342089534115))
targetgene="H2-K1"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(59.275525615949725,7.305691329705882,1.2864551526378198,15.545316974782523,8.648709425019511,5.331068417906823),c(113.27144976003348,45.58116112229539,14.151006679016017,77.72658487391261,17.297418850039023,23.989807880580702),c(111.35159467935495,8.2586075901023,30.231696086988766,15.545316974782523,34.594837700078045,13.327671044767058),c(75.59429380171726,18.264228324264703,13.186165314537654,34.97696319326067,0.0,34.65194471639435))
targetgene="Ugp2"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetgenelist=c("Sipa1","Ncor1","Nxf2","Itga2","Pck1","Dll1","Pak2","Atp2b1","Ikbkg","Git1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4,5_vs_0 pos.'


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
targetmat=list(c(103.432192471556,295.5628600996227,371.4639253241705,909.4010430247776,847.5735236519121,290.5432287759219),c(294.93773676923973,1673.6385920095777,371.78553911232996,730.6298978147786,1193.5219006526925,1372.750117611007),c(5.759565242035601,132.9318183253005,50.81497852919389,11.658987731086892,0.0,285.212160358015),c(46.55648570645444,45.10470299209718,905.9860412451847,237.06608386543348,17.297418850039023,125.28010782081034))
targetgene="Sipa1"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(158.1480622708942,461.8467475387979,396.87141458876744,217.63443764695532,2006.5005866045267,373.1747892534776),c(278.858950468557,316.84465658180943,1089.6275142842335,501.33647243673636,1383.7935080031218,573.0898549249835),c(43.43672120035182,73.21573267379155,54.67434398710734,73.84025563021699,665.9506257265024,271.88448931324797),c(27.837898669838737,13.182008268817134,0.32161378815945496,11.658987731086892,0.0,23.989807880580702))
targetgene="Ncor1"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(21.838351542718318,11.117356371291558,4.180979246072915,7.772658487391261,8.648709425019511,21.324273671627292),c(82.31378658409213,211.07095167780687,319.3624916423388,101.04456033608639,138.37935080031218,194.58399725359905),c(113.03146787494866,51.61629743813938,83.61958492145828,62.18126789913009,69.18967540015609,82.63156047755575),c(139.6694571193633,388.4721954882736,79.76021946354483,167.1121574789121,294.0561204506634,295.8742971938287))
targetgene="Nxf2"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(51.59610529323559,138.96695464114447,13.829392890856564,69.95392638652135,86.48709425019511,45.314081552207995),c(89.2732612515518,24.775822770306902,67.53889551348554,69.95392638652135,86.48709425019511,58.64175259697505),c(86.87344240070364,57.01615624705242,16.08068940797275,104.93088957978203,69.18967540015609,101.29029994022964),c(90.4731706769759,516.1629743813938,487.5665028497337,101.04456033608639,432.4354712509756,90.62816310441599))
targetgene="Itga2"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(28.557844325093185,27.952210304961632,40.20172351993187,19.431646218478154,17.297418850039023,31.986410507440937),c(127.67036286512248,106.25016303420075,128.00228768746308,112.70354806717329,103.78451310023414,93.2936973133694),c(120.7108881976628,130.3907082975767,236.06452050903994,205.97544991586844,242.16386390054632,77.30049205964893),c(119.99094254240835,523.1510269576341,701.1180581876118,120.47620655456456,501.6251466511317,223.90487355208657))
targetgene="Pck1"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(74.87434814646281,47.486993643088226,172.38499045346785,170.99848672260777,311.35353930070244,39.983013134301174),c(37.19719218814659,20.80533835198849,16.08068940797275,66.06759714282572,60.54096597513658,47.979615761161405),c(24.958116048820937,82.26843714755753,355.06162212803827,81.61291411760824,51.89225655011707,55.97621838802164),c(55.19583356950784,70.99226139953323,35.05590290938059,27.204304705869415,518.9225655011707,50.64514997011482))
targetgene="Dll1"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(8.159384092883768,4.288123171783887,15.437461831653838,11.658987731086892,121.08193195027316,0.0),c(133.9098918773277,598.1137727754858,506.863330139301,411.9508998317369,259.4612827505853,511.782568119055),c(44.396648740691084,69.0864288787404,467.6264479838475,136.02152352934706,8.648709425019511,50.64514997011482),c(81.8338228139225,97.67391669063298,142.796521942798,178.77114520999902,129.73064137529266,82.63156047755575))
targetgene="Pak2"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(36.95721030306177,21.2817964821867,20.261668654045664,34.97696319326067,43.24354712509756,2.6655342089534115),c(70.31469232985128,20.646518975255752,293.3117748014229,101.04456033608639,8.648709425019511,50.64514997011482),c(86.393478630534,396.09552557144497,93.26799856624194,128.2488650419558,164.3254790753707,42.648547343254585),c(149.26873252275598,519.8158200462467,372.1071529004894,614.0400205039097,1141.6296441025754,407.82673396987195))
targetgene="Atp2b1"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(50.39619586781151,26.046377784168794,243.14002384854794,120.47620655456456,138.37935080031218,77.30049205964893),c(79.67398584815913,401.8130231338235,300.38727814093096,194.31646218478153,259.4612827505853,98.62476573127623),c(99.35250042511412,448.8235586467135,58.533709445020804,124.36253579826018,51.89225655011707,170.59418937301834),c(53.2759784888293,48.75754865695012,129.2887428401009,38.86329243695631,233.5151544755268,50.64514997011482))
targetgene="Ikbkg"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
targetmat=list(c(84.71360543494029,148.49611724510868,77.8305367345881,62.18126789913009,43.24354712509756,146.60438149243762),c(29.51777186543245,23.346448379712275,10.934868797421469,101.04456033608639,25.946128275058534,15.993205253720468),c(131.0301092563099,309.53896525210354,421.314062488886,132.13519428565144,259.4612827505853,130.61117623871718),c(82.31378658409213,58.921988767845264,400.08755247036197,73.84025563021699,181.62289792540975,255.8912840595275))
targetgene="Git1"
collabel=c("In.vivo.input_Pool.A..1.8._1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1813_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1814_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1810_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1811_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409","In.vivo_AR1812_Pool.A..1.8._BM_hEGFRv3_15m_1360..1361..1362..1363..1406..1407..1408..1409")

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
Sweave("invivo.poolA_BM_hEGFRv3_15m_v_input_summary.Rnw");
library(tools);

texi2dvi("invivo.poolA_BM_hEGFRv3_15m_v_input_summary.tex",pdf=TRUE);

