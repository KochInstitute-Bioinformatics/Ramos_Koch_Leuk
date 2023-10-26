pdf(file='invivo.poolB_BM_mCD19_15m_v_BM_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolB_BM_mCD19_15m_v_BM_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Rps15a","Capzb","Polr2c","Hrh3","Rasal1","Gtf2h2","Tada3","Ppp1r8","Rpl27a","Haus8")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='9,10,11,12_vs_1,2,3,4,5 neg.'


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
targetmat=list(c(102.42798827322086,43.118437848342424,16.628621617989058,0.0,15.597748757652067,48.196559804117776,0.0,0.0,0.0),c(15.800062020869175,131.20987001162266,11.877586869992186,0.7541131285626209,5.793449538556482,0.0,11.891074646238506,1.300710512813602,93.1605199184138),c(205.40080627129927,160.88277347714862,428.3849664443848,6.787018157063589,17.82599858017379,0.0,0.0,0.0,11.943656399796641),c(2.724148624287789,286.0653349723363,123.52690344791873,36.197430171005806,80.21699361078205,13.770445658319364,0.0,0.0,0.0))
targetgene="Rps15a"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(20.158699819729637,1572.2002445562275,209.8373680365286,12.819923185564555,26.738997870260686,897.3740420671453,4.756429858495403,5.202842051254408,0.0),c(56.11746166032845,537.8213753126582,193.20874641853956,455.48432965182303,251.34657998045043,9.180297105546243,204.52648391530232,6.50355256406801,57.32955071902387),c(39.77256991460172,3.245473816541903,400.67059708106973,137.248589398397,61.49969510159958,0.0,0.0,0.0,0.0),c(445.1258852086247,142.80084792784373,94.22885583527133,35.443317042443184,168.4556865826423,0.0,0.0,1.300710512813602,4.777462559918656))
targetgene="Capzb"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(75.73133175520053,101.07332742944783,0.0,127.44511872708293,15.597748757652067,0.0,0.0,1.300710512813602,14.332387679755968),c(35.413932115741254,3.245473816541903,3.9591956233307286,15.082262571252418,8.021699361078205,0.0,0.0,0.0,0.0),c(2.724148624287789,1.39091734994653,55.4287387266302,0.0,11.586899077112964,2.295074276386561,0.0,2.601421025627204,52.55208815910522),c(224.4698466413138,40.33660314844936,365.82967559575934,13.574036314127177,7.57604939657386,39.01626269857153,0.0,0.0,2.388731279959328))
targetgene="Polr2c"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(546.4642140321305,189.6283987093769,654.8509560989025,708.8663408488636,40.999796734399716,369.5069584982363,154.58397040110057,420.12949563879346,54.94081943906455),c(0.5448297248575578,311.10184727137386,1345.3346728077815,206.62699722615812,251.79222994495478,0.0,0.0,1.300710512813602,90.77178863845447),c(143.83504736239524,229.50136274117742,899.5292456207416,198.3317528119693,323.09622426564994,220.32713053310982,356.7322393871552,303.0655494855693,905.3291551045853),c(28.87597541745056,45.900272548235485,7.918391246661457,26.393959499691732,37.43459701836496,32.13103986941185,7.134644787743104,59.83268358942569,0.0))
targetgene="Hrh3"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(109.51077469636911,365.3476239192885,699.1939470802066,242.82442739716393,491.9975608127966,4.590148552773122,11.891074646238506,189.9037348707859,47.774625599186564),c(406.98780446859564,459.46636459900367,1209.9301824898707,164.39666202665137,120.32549041617308,32.13103986941185,2501.8821055685817,3.9021315384408064,31.053506639471266),c(17.43455119544185,4.6363911664884325,19.795978116653643,272.9889525396688,11.586899077112964,0.0,0.0,0.0,0.0),c(70.28303450662496,11.590977916221082,1148.1667307659113,141.77326816977273,26.293347905756338,778.0301796950441,225.93041827853162,0.0,9.554925119837312))
targetgene="Rasal1"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(19.61387009487208,138.1644567613553,2.375517373998437,44.492674585194635,529.877807795666,0.0,0.0,0.0,16.721118959715298),c(16.344891745726734,9.736421449625709,2.375517373998437,8.29524441418883,32.97809737332151,654.0961687701698,0.0,1.300710512813602,16.721118959715298),c(3.813808074002904,70.00950661397533,6.334712997329166,52.03380587082084,8.46734932558255,34.426114145798415,21.40393436322931,318.6740756393325,0.0),c(90.98656405121214,289.3108087888782,11.08574774532604,255.6443505827285,136.81453910283383,20.655668487479048,0.0,0.0,0.0))
targetgene="Gtf2h2"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(45.220867163177296,1.854556466595373,10.293908620659893,5.278791899938346,9.35864925459124,153.76997651789958,49.94251351420173,226.32362922956676,11.943656399796641),c(0.0,89.48234951322675,3.9591956233307286,12.819923185564555,1.3369498935130342,0.0,0.0,16.909236666576827,0.0),c(17.979380920299405,127.03711796178305,571.7078480089572,392.1388268525629,175.14043605020748,34.426114145798415,0.0,1.300710512813602,2.388731279959328),c(0.0,4.6363911664884325,37.216438859308845,8.29524441418883,0.8912999290086895,224.91727908588297,0.0,0.0,0.0))
targetgene="Tada3"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(384.1049560245782,130.28259177832496,16.628621617989058,99.54293297026597,1.782599858017379,0.0,0.0,0.0,16.721118959715298),c(1.0896594497151155,0.0,5.54287387266302,1.5082262571252418,10.695599148104273,0.0,0.0,3.9021315384408064,11.943656399796641),c(15.800062020869175,86.7005148133337,11.877586869992186,293.35000701085954,117.65159062914701,263.9335417844545,0.0,2.601421025627204,0.0),c(537.2021087095519,98.75513184620361,52.26138222796562,65.60784218494803,468.8237626585707,50.49163408050434,223.55220334928393,74.14049923037531,965.0474371035685))
targetgene="Ppp1r8"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(32.68978349145347,0.0,11.877586869992186,54.29614525650871,85.11914322032985,130.81923375403397,0.0,0.0,0.0),c(1.6344891745726733,0.0,2.375517373998437,42.98444832806939,84.22784329132115,41.311336974958095,0.0,0.0,11.943656399796641),c(1.0896594497151155,74.64589778046377,38.00827798397499,417.02456009512935,1.3369498935130342,18.360594211092486,2.3782149292477013,1.300710512813602,14.332387679755968),c(141.65572846296502,44.04571608164011,7.918391246661457,81.44421788476306,27.18464783476503,0.0,0.0,2.601421025627204,0.0))
targetgene="Rpl27a"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(9.262105322578481,5.100030283137276,53.845060477297906,153.83907822677466,12.032549041617308,29.83596559302529,0.0,28.615631281899244,169.59992087711228),c(111.14526387094179,52.39122018131929,14.253104243990622,20.361054471190766,299.922426111424,1888.8461294661395,152.20575547185288,477.36075820259197,1086.8727323814942),c(232.64229251417717,0.9272782332976865,234.38438090117913,15.082262571252418,192.9664346303813,0.0,14.269289575486209,0.0,0.0),c(51.21399413661043,186.382924892835,163.118859681226,61.083163413572294,18.271648544678133,22.95074276386561,49.94251351420173,234.12789230644836,90.77178863845447))
targetgene="Haus8"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetgenelist=c("Ints10","Sephs1","Tmprss13","E130309D02Rik","Dcp2","Arid1b","Kmt2d","Gnpnat1","Grik2","Mgat4d")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='9,10,11,12_vs_1,2,3,4,5 pos.'


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
targetmat=list(c(233.18712223903472,50.07302459807507,11.08574774532604,404.95875003812745,48.13019616646923,215.7369819803367,4214.196854626926,2475.252105884285,169.59992087711228),c(390.6429127228689,190.55567694267458,753.0390075575045,963.0024651744669,415.3457669180493,100.98326816100868,1902.571943398161,303.0655494855693,83.60559479857649),c(204.85597654644172,223.0104151080936,63.34712997329166,64.09961592782278,254.0204797674765,3635.397653796312,1086.8442226661996,7773.046024574086,40.60843175930858),c(286.03560555021784,466.8845904653852,954.9579843473717,1137.956711000995,282.9877274602589,771.1449568658844,4221.331499414669,6481.440485350179,1352.0219044569797))
targetgene="Ints10"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(153.09715268497374,17.618286432656046,46.718508355302596,18.85282821406552,17.82599858017379,798.6858481825232,856.1573745291724,721.8943346115492,845.6108731056022),c(270.23554352934866,57.02761134780772,1339.7917989351185,92.00180168463976,49.02149609547792,169.8354964526055,756.272347500769,4640.935109718932,22454.074031617685),c(33.234613216311025,175.25558609326276,227.25782877918382,120.65810057001934,209.90113328154638,73.44237684436995,1177.2163899776122,584.0190202533073,499.24483751149955),c(187.42142535099987,170.1555558101255,254.9721981424989,25.63984637112911,132.3580394577904,80.32759967352963,1610.0515071006937,1065.2819099943401,7457.619056033022))
targetgene="Sephs1"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(28.87597541745056,108.49155329582932,115.60851220125727,96.52648045601548,416.237066847058,371.80203277462283,49.94251351420173,122.26678820447859,243.65059055585147),c(290.93907307393584,439.5298825831034,342.86634098044107,104.8217248702043,177.81433583723356,1920.9771693355515,1717.0711789168404,20508.302655532065,57.32955071902387),c(13.075913396581386,38.482046681853994,56.22057785129635,49.01735335657036,64.61924485312998,78.03252539714306,768.1634221470075,518.9834946126272,5202.656727751417),c(87.7175857020668,83.4550409967918,35.63276060997656,351.4167179101813,135.92323917382515,291.4744331010932,2171.3102304031513,339.48544384435013,7732.323153228345))
targetgene="Tmprss13"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(692.4785802939559,121.0098094453481,627.1365867355873,484.14062853720264,164.44483690210322,401.63799836764815,6024.018415784427,5288.688945100106,18713.320847201376),c(476.18117952550546,234.13775390766585,521.0301440303239,815.1962919761933,117.65159062914701,20.655668487479048,901.3434581848788,13.00710512813602,64.49574455890186),c(300.2011783965143,444.6299128662407,58.59609522529478,711.8827933631142,79.32569368177336,367.21188422184974,21751.153742899474,1.300710512813602,2.388731279959328),c(1486.2954894114175,311.5654863880227,308.0254194951307,3330.9176888610964,457.2368635814577,6706.207035601531,28160.442977222032,2031.7098210148465,2388.7312799593283))
targetgene="E130309D02Rik"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(669.5957318499385,316.2018775545111,1568.6333059636347,204.36465784047027,169.346986511651,4121.953400390264,30581.46577519619,6862.548665604564,21990.660163305576),c(171.07653360527314,40.33660314844936,257.3477155164974,61.837276542134916,198.3142342044334,2.295074276386561,23.782149292477012,83.24547282007053,0.0),c(423.8775259391799,410.3206182342263,670.6877385922254,574.6342039647171,224.1619321456854,3745.561219062867,61.833588160440236,349.89112794685894,599.5715512697914),c(750.7753608537146,255.0015141568638,1337.41628156112,337.84268159605415,159.09703732805107,619.6700546243715,187.8789794105684,1207.0593558910227,7.166193839877984))
targetgene="Dcp2"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(627.6438430359066,342.62930720349516,227.25782877918382,82.95244414188831,37.8802469828693,3713.4301791934554,6475.87925234149,534.5920207663904,1987.424424926161),c(213.57325214416264,114.51886181226429,45.1348301059703,288.0712151109212,174.2491361211988,240.98279902058889,425.70047233533853,102.75613051227457,54.94081943906455),c(243.5388870113283,127.96439619508074,471.1442791763567,125.9368924699577,86.90174307834722,814.751368117229,351.9758095286598,4996.0290797170455,2278.849641081199),c(277.3183299524969,533.6486232628187,288.229441378477,177.2165852122159,223.71628218118107,2106.878185722863,1070.1967181614655,4752.796213820902,757.227815747107))
targetgene="Arid1b"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(81.72445872863366,45.900272548235485,18.212299867321352,66.36195531351063,23.619448118730272,3695.069584982363,1809.8215611575006,97.55328846102016,7.166193839877984),c(79.54513982920344,181.74653372634657,718.1980860721942,441.91029333769586,143.49928857039902,571.4734948202537,3182.0515753334244,1692.2243771704964,2567.8861259562777),c(748.0512122294268,69.5458674973265,581.2099175049509,333.31800282467844,327.107073946189,94.098045331849,33.29500900946782,123.5674987172922,138.54641423764102),c(127.49015561666852,366.73854126923504,193.20874641853956,175.70835895509066,76.20614393024294,206.55668487479048,658.7655354016133,5920.834254327517,511.18849391129623))
targetgene="Kmt2d"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(0.0,10.663699682923395,94.22885583527133,112.36285615583051,6.239099503060826,114.75371381932804,3241.506948564617,2685.967208960088,9.554925119837312),c(58.84161028461624,232.28319744107048,24.547012864650515,1.5082262571252418,139.93408885436426,114.75371381932804,19.02571943398161,58.53197307661209,0.0),c(1044.4385825519382,384.8204668185399,346.0336974791057,89.73946229895189,479.51936180667497,1023.6031272684061,347.2193796701644,788.2305707650429,458.63640575219097),c(38.682910464886604,84.38231923008948,517.8627875316593,23.377506985441247,23.619448118730272,13.770445658319364,42.80786872645862,9.104973589695215,85.9943260785358))
targetgene="Gnpnat1"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(56.66229138518601,71.40042396392187,23.75517373998437,31.67275139963008,134.58628928031212,810.161219564456,19.02571943398161,304.3662599983829,2761.3733596329835),c(203.22148737186905,58.418528697754255,52.26138222796562,55.050258385071324,205.89028360100727,231.80250191504265,2157.040940827665,128.77034076854662,487.30118111170293),c(547.5538734818456,245.72873182388693,152.03311193589997,92.75591481320237,471.0520124810924,358.0315871163035,111.77610167464196,1238.2764081985492,1039.0981067823077),c(189.05591452557255,472.4482598651713,862.3128067614326,203.61054471190764,293.6833266083632,5471.457074905561,83.23752252366954,1127.716014609393,487.30118111170293))
targetgene="Grik2"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(129.66947451609875,6.954586749732649,155.9923075592307,260.1690293541042,51.69539588250399,50.49163408050434,33.29500900946782,71.53907820474811,57.32955071902387),c(663.6026048765053,437.2116869998592,133.82081206857862,364.99075422430855,812.8655352559248,433.76903823706,10613.973229232492,117.06394615322418,441.9152867924757),c(37.048421290313925,146.04632174438564,62.55529084862551,126.69100559852032,67.29314464015606,105.5734167137818,2297.3556216532793,650.355256406801,23.887312799593282),c(13.075913396581386,2.3181955832442163,23.75517373998437,217.18458102603483,7.57604939657386,794.09569962975,149.8275405426052,2003.0941897329471,2854.533879551397))
targetgene="Mgat4d"
collabel=c("In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1836_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1840_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1841_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1842_Pool.B..9.16._BM_mCD19_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
Sweave("invivo.poolB_BM_mCD19_15m_v_BM_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolB_BM_mCD19_15m_v_BM_hEGFRv3_15m_summary.tex",pdf=TRUE);

