pdf(file='invivo.poolC_BM_hEGFRv3_15m_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolC_BM_hEGFRv3_15m_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Mars2","Trim28","Kif18b","Prmt1","Ndufb5","Pfdn5","Mto1","Mcm4","Ccnh","Srp54a_Srp54b_Srp54c")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4_vs_0 neg.'


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
targetmat=list(c(109.87044807703386,7.7748673840782265,0.0,1.8551900669845274,12.538036736567395),c(108.77174359626352,6.803008961068448,2.010267062371475,61.2212722104894,0.0),c(57.1326330000576,0.0,6.030801187114425,0.0,23.147144744432115),c(216.4447827117567,10.690442653107562,27.138605342014912,212.4192626697284,5.786786186108029))
targetgene="Mars2"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(476.8377446543269,63.17079749563559,23.11807121727196,543.5706896264666,180.35483613370022),c(522.9833328466812,133.14460395233962,54.277210684029825,215.20204777020518,159.13662011797078),c(752.6125693276819,56.36778853456714,209.0677744866334,209.63647756925158,533.34879348629),c(363.6711831349821,231.30230467632725,5.025667655928688,26.900255971275648,135.02501100918732))
targetgene="Trim28"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(90.09376742316776,0.0,26.133471810829175,36.17620630619828,4.822321821756691),c(337.3022755964939,112.73557706913428,101.51848664975948,45.45215664112092,21.21821601572944),c(164.80567211555078,6.803008961068448,543.7772403714839,15.769115569368482,57.867861861080286),c(126.35101528858894,25.268318998254237,1.0051335311857374,10.2035453684149,2.8933930930540144))
targetgene="Kif18b"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(277.97223363489564,29.15575269029335,338.73000000959354,46.37975167461318,81.97947096986374),c(293.3540963656804,89.4109749168996,621.1725222727857,47.30734670810545,52.081075674972254),c(181.28623932710587,107.87628495408539,38.195074185058026,24.117470870798854,12.538036736567395),c(195.56939757712027,24.296460575244456,28.143738873200647,618.7058873393398,63.65464804718832))
targetgene="Prmt1"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(323.01911734647956,23.32460215223468,5.025667655928688,458.2319465451783,87.76625715597177),c(113.16656151934487,41.78991218942047,19.097537092529013,11.131140401907164,92.58857897772846),c(76.9093136539237,5.83115053805867,4.02053412474295,38.03139637318281,18.324822922675423),c(53.83651955774659,5.83115053805867,1.0051335311857374,0.9275950334922637,4.822321821756691))
targetgene="Ndufb5"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(424.0999295773507,87.46725807088005,79.40554896367325,499.9737230523301,68.476969868945),c(432.8895654235134,23.32460215223468,116.59548961754554,36.17620630619828,72.33482732635035),c(305.4398456541541,107.87628495408539,77.39528190130179,211.49166763623612,35.68518148099951),c(329.61134423110155,380.9685018198331,9.046201780671637,38.03139637318281,70.40589859764768))
targetgene="Pfdn5"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(437.2843833465947,114.67929391515383,19.097537092529013,317.2375014543542,4.822321821756691),c(342.79579800034566,394.57451974197,46.23614243454392,317.2375014543542,39.54303893840486),c(441.67920126967607,192.4279677559361,32.1642729979436,93.68709838271863,15.43142982962141),c(217.54348719252704,13.606017922136896,29.148872404386385,207.78128750226708,8.680179279162044))
targetgene="Mto1"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(104.37692567318216,0.0,17.087270030157537,68.64203247842751,1.9289287287026762),c(212.04996478867534,56.36778853456714,9.046201780671637,16.696710602860747,42.43643203145888),c(104.37692567318216,0.0,103.52875371213096,17.62430563635301,16.395894193972747),c(82.40283605777539,0.9718584230097783,53.27207715284408,60.29367717699714,4.822321821756691))
targetgene="Mcm4"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(294.4528008464507,2.915575269029335,28.143738873200647,20.4070907368298,10.60910800786472),c(158.21344523092876,135.08832079835918,56.287477746401294,245.81268387544986,53.0455400393236),c(405.4219534042549,61.227080649616035,106.54415430568817,8.348355301430374,126.34483173002529),c(136.239355615522,9.718584230097782,11.056468843043112,679.9271595498293,24.11160910878345))
targetgene="Ccnh"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(117.56137944242623,529.6628405403292,93.47741840027358,38.958991406675075,148.52751211010607),c(59.33004196159828,0.9718584230097783,5.025667655928688,4.6379751674613185,0.0),c(239.5175768079338,85.5235412248605,116.59548961754554,73.28000764588883,72.33482732635035),c(116.46267496165589,53.45221326553781,373.9096736010943,13.913925502383956,248.83180600264524),c(125.25231080781859,81.63610753282138,47.24127596572966,42.66937154064413,12.538036736567395),c(119.7587884039669,13.606017922136896,10.051335311857375,24.117470870798854,3.8578574574053524),c(45.04688371158388,0.9718584230097783,0.0,0.0,18.324822922675423),c(384.5465682696185,639.4828423404341,152.7802967402321,275.4957249472023,264.2632358322666))
targetgene="Srp54a_Srp54b_Srp54c"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetgenelist=c("Fibp","Tmem94","Ubr4","Fubp1","Gm46977","Gbp4","9530053A07Rik","4930555G01Rik_Gm21560_Gm3159_Gm3594_Gm3642","Gcc1","Eomes")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4_vs_0 pos.'


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
targetmat=list(c(114.2652660001152,330.43186382332465,31.15913946675786,11.131140401907164,34.720717116648174),c(343.894502481116,4800.980609668305,904.6201780671637,2273.535427089538,1925.0708712452708),c(381.25045482730746,124.39787814525162,1028.2516024030094,77.91798281335015,73.2992916907017),c(689.9864139237726,1714.358258189249,4995.513649993115,4943.153933480273,1184.3622394234433))
targetgene="Fibp"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(165.90437659632113,353.7564659755593,448.2895549088389,724.451721157458,537.2066509436953),c(360.37506969267105,561.7341684996519,1107.6571513666827,1487.862433721591,2961.8700629229593),c(202.1616244617423,949.5056792805534,101.51848664975948,583.4572760666339,164.9234063040788),c(293.3540963656804,197.287259870985,1227.2680415777854,1336.664443262352,234.36484053737516))
targetgene="Tmem94"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(829.5218829816056,537.4377079244074,1532.8286350582496,2264.2594767546157,1085.0224098952554),c(759.2047962123039,1542.3393173165182,482.46409496915396,3465.495045127097,2284.81607914832),c(552.6483538274803,1488.8871040509803,651.3265282083579,1261.5292455494787,1681.0613870643822),c(291.1566874041397,604.4959391120822,1439.351216657976,291.2648405165708,721.4193445348009))
targetgene="Ubr4"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(162.6082631540101,35.9587616513618,28.143738873200647,11.131140401907164,9.644643643513382),c(281.26834707720667,964.0835556257001,274.4014540137063,1420.1479962766557,1215.225099082686),c(535.0690821351549,87.46725807088005,2763.112077229592,401.6486495021502,864.1600704587989),c(350.486729365738,463.57646777566424,1229.2783086401569,2081.5232551566396,822.6881027916913))
targetgene="Fubp1"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(339.49968455803463,493.7040788889674,530.7105044660693,155.8359656267003,1042.5859778637964),c(325.2165263080202,587.9743459209159,951.8614540328933,930.3778185927405,1626.086918296356),c(369.16470553883374,825.1078011353018,2056.503204806019,576.0365157986957,2191.26303580624),c(319.72300390416854,556.874876384603,332.6991988224791,1316.2573525255223,329.8468126081576))
targetgene="Gm46977"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(592.2017151352125,421.78655558624376,465.37682493899644,641.8957631766465,876.6981071953663),c(382.34915930807784,409.1523960871167,338.73000000959354,636.3301929756929,719.4904158060982),c(369.16470553883374,1176.9205502648415,1576.0493768992362,1425.7135664776092,377.1055664613732),c(431.79086094274305,1509.2961309341856,813.1530267292616,1916.4113391950168,670.30273322418))
targetgene="Gbp4"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(448.2714281542981,1594.8196721590461,234.19611276627683,721.6689360569811,2719.7895074707735),c(150.52251386553638,850.376120133556,206.05237389307618,822.7767947076379,60.7612549541343),c(312.03207253877616,190.48425090991654,161.82649852090373,509.24967338725276,557.4604025950734),c(275.774824673355,95.24212545495827,206.05237389307618,67.71443744493524,648.1200528440992))
targetgene="9530053A07Rik"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(97.78469878856014,837.7419606344289,42.21560830980097,270.85774977974097,964.4643643513381))
targetgene="4930555G01Rik_Gm21560_Gm3159_Gm3594_Gm3642"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(260.39296194257025,457.7453172376056,346.7710682590794,372.89320346389,1000.1495458323376),c(182.3849438078762,1153.595948112607,780.988753731318,276.4233199806946,1380.1485053867648),c(161.50955867323978,10.690442653107562,17.087270030157537,489.7701776839152,203.50198087813234),c(382.34915930807784,111.7637186461245,317.622195854693,94.6146934162109,271.97895074707736))
targetgene="Gcc1"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(309.83466357723546,218.66814517720013,138.70842730363177,476.7838472150235,52.081075674972254),c(268.0838933079626,128.28531183729075,523.6745697477692,699.4066552531668,630.7596942857751),c(176.8914214040245,161.3284982196232,313.6016617299501,6.493165234445846,168.78126376148415),c(151.62121834630673,306.1354032480802,1230.2834421713426,307.0339560859393,1472.7370843644933))
targetgene="Eomes"
collabel=c("In.vivo.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
Sweave("invivo.poolC_BM_hEGFRv3_15m_v_input_summary.Rnw");
library(tools);

texi2dvi("invivo.poolC_BM_hEGFRv3_15m_v_input_summary.tex",pdf=TRUE);

