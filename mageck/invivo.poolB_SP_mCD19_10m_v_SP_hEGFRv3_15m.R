pdf(file='invivo.poolB_SP_mCD19_10m_v_SP_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolB_SP_mCD19_10m_v_SP_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Capzb","Rasa3","Fkrp","Jmjd1c","Irgm1","Tsc22d4","Slc26a6","Cnot8","Vim","Rpl7a")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='19,20,21,22,23_vs_13,14,15,16,17,18 neg.'


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
targetmat=list(c(190.15536931922625,498.0728732454788,349.98346545317816,130.9079042941825,281.14128844271784,128.4566248592635,13.352404695943893,2.4992114544338024,1.7351044316743645,5.879097716155285,0.0),c(610.5962995347883,840.3540773556228,1051.193676876242,449.61617761799386,899.1874266721636,702.8275185566687,271.4988954841925,197.4377049002704,23.42390982760392,149.91699176195976,464.41113007274623),c(275.3548529752432,49.11658534315969,192.08683983131803,42.53126004494536,173.09938833869816,197.32987758929258,5.563501956643289,0.0,6.940417726697458,47.03278172924228,33.65298044005407),c(269.79836491072035,432.0724616906079,187.11371776448777,71.80602345250516,224.21598623737412,386.9715781296517,85.67793013230664,27.491325998771828,181.3184131099711,48.5025561582811,346.62569853255695))
targetgene="Capzb"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(257.4506136562252,218.722294106258,299.0089642681682,316.49885721758045,141.73238508269245,176.50773141509774,93.46683287160725,989.6877359557858,7.8079699425346405,63.20030044866932,400.47046723664346),c(510.579514373377,828.8423776658198,899.5134538379197,829.0833938631557,852.7177922188218,724.290346151608,322.68311348531074,219.93060799017462,628.9753564819571,424.76480999221934,185.0913924202974),c(212.9987091400424,599.3758305157456,230.0068955908986,571.134063460695,462.3728628107508,408.11406501421874,36.719112913845706,94.9700352684845,136.2056978864376,457.09984743107344,0.0),c(176.57284293928151,283.1878123691551,162.86974768869035,147.47852509091445,312.50829169872355,271.64892239564955,142.4256500900682,3253.973313672811,0.8675522158371822,174.90315705561974,434.12344767669754))
targetgene="Rasa3"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(88.28642146964076,99.76806397829313,123.08477115404844,21.54180703575155,238.15687657337668,251.14711693182693,0.0,112.46451544952112,43.37761079185911,10.288421003271749,0.0),c(11.112976129045691,0.0,39.16333627628814,18.780036902962888,5.808704306667724,46.76974371434532,0.0,0.0,0.0,83.77714245521281,33.65298044005407),c(463.04067204357045,135.83805633967603,333.82081873597986,298.82352836773305,354.33096270673116,269.72687813341616,124.62244382880966,387.3777754372394,46.84781965520784,201.3590967783185,141.34251784822712),c(51.86055526887989,72.1399847227658,103.81392314508126,128.69848818795157,206.78987331737096,42.284973769134126,77.88902739300605,7.497634363301407,8.675522158371823,251.33142736563843,0.0))
targetgene="Fkrp"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(74.08650752697127,204.14080783250748,506.6368105583307,197.74274150766806,465.8580853947514,350.7730778575899,154.66535439468342,42.486594725374644,242.91462043441103,44.09323287116464,67.30596088010815),c(311.7807191760041,334.6067376502754,291.5492811679228,115.44199155056599,196.33420556536905,304.96435627436125,50.071517609789595,449.85806179808446,91.09298266290413,44.09323287116464,23.557086308037853),c(441.43210734820383,405.9792757270543,714.886297106847,783.238009658864,741.1906695308015,200.53328469301488,349.38792287719855,354.8880265295999,120.58975800136832,293.95488580776424,737.0002716371843),c(443.2842700363781,296.23440535093187,171.57271130564328,364.55365752810314,206.78987331737096,169.78057649728095,32.26831134853107,164.94795599263097,50.31802851855657,77.89804473905753,673.0596088010815))
targetgene="Jmjd1c"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(613.0658497856873,735.9813335014085,660.8035946300681,631.3406523554877,678.4566630187901,602.2405354997891,258.14649078824857,377.3809296195042,195.199248563366,534.997892170131,333.16450635653536),c(253.7462882798766,367.6069434277108,330.712617444211,235.85516934015158,475.1520122854198,138.70752759117482,52.29691839244691,157.45032162932955,82.41746050453231,51.442105016358745,121.15072958419466),c(24.695502508990423,30.697865839474808,22.37904930073608,44.188322124618566,12.779149474668992,61.18507568109559,8.901603130629262,1524.5189872046194,0.0,0.0,0.0),c(254.36367584260137,232.53633373402167,256.11578644175734,461.21561217570627,772.5576727868072,249.22507266959354,459.54526161873565,262.41720271554925,402.54422814845253,130.8099241844551,925.4569621014871))
targetgene="Irgm1"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(238.92898677448235,644.6551826289709,320.7663733105505,401.0090232809134,362.46314873606593,405.2309986208687,149.10185243804014,282.4108943510197,108.44402697964777,465.91849400530634,111.05483545217845),c(133.97310111127305,181.88485509888824,101.32736211166615,55.2354026557732,47.63137531467533,151.2008152956917,67.87472387104812,27.491325998771828,76.34459499367203,124.93082646829981,0.0),c(204.9726708246205,325.397377898433,211.3576878402852,371.1819058467959,513.4894607094268,275.49301092011626,172.46856065594196,54.982651997543655,8.675522158371823,160.20541276523153,97.59364327615681),c(399.4497530829201,321.5601446684986,491.71744435784,270.6534730132887,513.4894607094268,487.55856118653134,301.54180605006627,82.47397799631548,170.9077865199249,511.4815013055098,134.6119217602163))
targetgene="Tsc22d4"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(498.84915068160655,416.72352877087053,338.17230054445633,189.4574311093021,512.3277198480932,206.29941747971498,145.76375126405415,87.47240090518308,76.34459499367203,116.11217989406688,50.47947066008111),c(493.2926626170837,237.1410136099429,262.3321890252952,247.45460389786393,276.4943249973836,375.4393125562515,430.61505144419056,167.44716744706477,541.3525826824017,827.4830035488563,181.726094376292),c(191.39014444467577,161.16379565724273,72.73191022739226,170.6773942063392,166.12894317069689,186.1179527262646,158.0034555686694,187.44085908253518,251.59014259278285,183.72180362985267,811.0368286053032),c(150.02517774211682,178.81506851494075,137.38249709618538,256.8446223493454,166.12894317069689,121.08878852070227,364.96572835579974,0.0,7.8079699425346405,70.54917259386342,0.0))
targetgene="Slc26a6"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(1026.7155168112768,952.4012876697059,620.3969778370724,937.8971370950289,1696.1416575469752,1257.0169475006235,643.1408261879642,407.3714670727098,226.43112833350457,507.07217801839334,1016.320009289633),c(698.882721004429,959.3083074835877,1119.5741052951578,1147.2393131604094,880.5995728908268,850.504586038266,708.790149276355,247.42193398894645,491.9021063796823,1233.1407459635711,693.2513970651139),c(805.6907693558126,429.7701217526473,1219.6581868901164,957.7818820511073,535.5625370747641,520.8739950652431,206.96227278713033,259.91799126111545,244.6497248660854,671.6869140707413,649.5025224930437),c(736.5433623306394,668.446028654564,775.1854021671637,1153.8675614791023,678.4566630187901,613.7728010731893,893.4984142369121,399.87383270940836,694.9093248855829,483.5557871537722,3.3652980440054074))
targetgene="Cnot8"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(110.51237372773214,38.37233229934351,66.51550764385446,95.55724659448764,89.45404632268294,39.7222480861563,200.28607043915838,682.284727060428,38.172297496836016,1521.21653405518,861.5162992653843),c(245.72024996445472,277.0482392012601,261.0889085085876,169.020332126666,58.08704306667723,171.3822800491421,429.50235105286185,57.48186345197745,81.54990828869514,151.3867661909986,450.9499378967246),c(437.7277819718553,12910.754925437117,3696.8946164299296,16998.695167314203,1502.1309337042733,2053.704294196355,134.6367473507676,107.4660925406535,519.6637772864722,831.8923268359729,2486.9552545199963),c(305.6068435487565,330.00205777435417,192.7084800896718,380.01957027171966,296.2439196400539,271.64892239564955,199.17337004782974,64.97949781527886,309.7161410538741,368.91338168874415,413.9316594126651))
targetgene="Vim"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(35.808478638036114,94.39593745638503,82.67815436105275,145.26910898468353,118.49756785602156,16.65771693935587,318.2323119199961,0.0,0.8675522158371822,10.288421003271749,0.0),c(80.87777071694363,11.511699689803052,36.05513498451924,8.837664424923712,17.42611292000317,48.05110655583423,153.55265400335477,44.98580617980844,172.64289095159927,4.409323287116464,47.114172616075706),c(206.82483351279478,227.93165385810045,286.57615910109257,397.1425450950093,659.8688092374534,355.89852922354555,1037.036764718309,182.44243617366757,271.54384355703803,108.76330774887278,373.5480828846002),c(45.06929207890752,99.00061733230625,25.48725059250498,16.01826677017423,112.68886354935384,43.56633661062304,0.0,0.0,14.748387669232098,29.395488580776426,0.0))
targetgene="Rpl7a"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetgenelist=c("Sephs1","Hdac7","Ewsr1","Pdcd5","Vim","Aldh1a1","Kmt2b","Bpifb9a","Gm9141_Tdpoz2","B3gnt2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='19,20,21,22,23_vs_13,14,15,16,17,18 pos.'


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
targetmat=list(c(121.00796229405307,68.30275149283145,77.7050322942225,73.46308553217835,80.16011943201458,86.17165109012939,124.62244382880966,9.99684581773521,592.5381634167954,574.6818017541791,710.077887285141),c(308.07639379965553,212.58272093836305,190.2219190562567,122.07023986925877,147.54108938936017,93.53948742869063,923.5413248027859,484.8470221601577,604.6838944385161,97.0051123165622,551.9088792168868),c(201.88573301099672,282.42036572316823,218.1957306821768,222.04631867620827,317.1552551440577,290.549024307611,1008.1065545437639,297.40616307762247,33.834536417650106,1012.6745816077479,1457.1740530543414),c(250.04196290352803,178.81506851494075,225.03377352406835,203.2662817732454,185.87853781336716,113.40061147176878,484.0246702279661,724.7713217858027,360.9017217882678,668.7473652126637,2813.3891647885207))
targetgene="Sephs1"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(68.53001946244842,134.30316304770227,124.94969192910978,197.74274150766806,74.35141512534686,193.16544835445362,291.5275025281083,339.89275780299715,463.2728832570553,338.0481186789289,656.2331185810544),c(113.59931154135595,202.60591454053372,164.1130282053979,81.74839593054433,74.35141512534686,135.50412048745252,221.42737787440288,832.2374143264562,199.5370096425519,867.1669131329046,588.9271577009463),c(98.78201003596169,100.53551062428,159.76154639692146,98.87137075383403,59.24878392801078,169.78057649728095,406.1356428349601,339.89275780299715,1066.221673263897,623.1843579124602,814.4021266493086),c(152.49472799301586,462.77032753008274,211.97932809863897,148.03087911747218,332.25788634139377,146.71604535048053,389.4451369650302,532.3320397944,697.5119815330945,617.305260196305,922.0916640574817))
targetgene="Hdac7"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(318.57198236597645,317.72291143856427,196.43832163979448,467.2915064678413,144.05586680535956,321.94241392408935,371.64193070377166,414.8691014360112,1006.3605703711314,114.64240546502806,390.37457310462725),c(627.2657637283568,568.6779646762708,387.90352121275873,583.285852044965,580.8704306667723,642.923805717062,745.5092621902007,14825.322347701316,709.657712554815,1316.9178884187838,2140.329555987439),c(289.55476691791273,275.5133459092864,345.0103433863479,281.1481995178856,347.3605175387299,466.4160743019643,1055.9526713708963,949.7003526848449,896.1814389598093,429.1741332793358,642.7719264050328),c(301.28513060968316,443.584161380411,387.2818809544049,150.79264925026084,231.1864314053754,326.42718386930056,235.89248296167543,664.7902468793915,275.01405242038675,809.8457104003905,1295.639746942082))
targetgene="Ewsr1"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(301.28513060968316,248.65271329974593,248.656103341512,173.99151836568558,255.58298949337984,220.39440873609303,362.7403275731424,2306.7721724423996,919.6053487874132,568.8027040380239,821.1327227373195),c(212.9987091400424,283.95525901514196,241.81806049962043,118.75611570991238,211.43683676270513,229.68428933688764,133.52404695943892,842.2342601441915,552.6307614882851,589.3795460445673,1413.4251784822711))
targetgene="Pdcd5"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(110.51237372773214,38.37233229934351,66.51550764385446,95.55724659448764,89.45404632268294,39.7222480861563,200.28607043915838,682.284727060428,38.172297496836016,1521.21653405518,861.5162992653843),c(245.72024996445472,277.0482392012601,261.0889085085876,169.020332126666,58.08704306667723,171.3822800491421,429.50235105286185,57.48186345197745,81.54990828869514,151.3867661909986,450.9499378967246),c(437.7277819718553,12910.754925437117,3696.8946164299296,16998.695167314203,1502.1309337042733,2053.704294196355,134.6367473507676,107.4660925406535,519.6637772864722,831.8923268359729,2486.9552545199963),c(305.6068435487565,330.00205777435417,192.7084800896718,380.01957027171966,296.2439196400539,271.64892239564955,199.17337004782974,64.97949781527886,309.7161410538741,368.91338168874415,413.9316594126651))
targetgene="Vim"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(590.2225099648712,480.42160038778076,593.0448064695062,616.4270936384289,536.7242779360977,380.8851046325793,1023.6843600223651,404.872255618276,805.9560085127423,1668.1939769590622,625.9454361850057),c(504.4056387461294,477.35181380383324,304.6037265933522,325.8888756690619,446.10849075208114,400.7462286756575,329.35931583328266,1464.5379122982083,1932.0387846694048,392.4297725533653,743.730867725195),c(547.0053805741379,448.95628790231905,443.22950420624517,261.26345456180724,365.9483713200666,389.2139631022573,626.4503203180343,264.91641416998306,758.2406366416973,636.4123277738096,346.62569853255695),c(330.9197336204717,535.6777588988354,418.36389387209397,245.24518779163301,333.4196272027273,518.9519508030097,662.0567328405514,212.4329736268732,222.96091947015583,768.6920263873035,844.6898090453573))
targetgene="Aldh1a1"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(188.30320663105198,306.9786583947481,316.414891502074,261.26345456180724,454.24067678141597,384.7291931570461,487.3627714019521,349.8896036207323,1113.937045134942,314.53172781430777,827.8633188253302),c(174.10329268838248,171.14060205507204,135.51757632112404,149.13558717058766,59.24878392801078,270.0472188437884,441.7420553574771,357.38723798403373,449.3920478036604,293.95488580776424,67.30596088010815),c(443.2842700363781,322.32759131448546,173.43763208070462,309.31825487232993,300.89088308538805,284.4625508105387,163.5669575253127,474.85017634242246,631.5780131294687,865.6971387038657,1389.8680921742332),c(363.0238868821592,273.2110059713258,346.87526416140923,246.3498958447485,152.18805283469436,226.80122294353757,212.52577474377364,532.3320397944,1604.10404708295,655.5193953513143,706.7125892411356))
targetgene="Kmt2b"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(280.2939534770413,318.49035808455113,369.8759537204991,600.9611808948124,379.8892616560691,148.95843032308613,447.3055573141204,1034.6735421355943,287.1597834421073,859.8180409877104,467.7764281167516),c(1067.463095951111,1009.9597861187211,1127.6554286537569,1355.4767811726745,1296.502801248236,1007.1511934102855,853.4412001490805,2206.8037142650473,2062.171617044982,959.7627021623503,4317.677390458938),c(196.32924494647386,239.44335354790348,165.3563087221055,122.62259389581651,138.2471624986918,233.5283778613544,457.3198608360783,1002.1837932279548,743.4922489724652,402.718193556637,679.7902048890923))
targetgene="Bpifb9a"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(48.77361745525609,48.34913869717282,42.27153756805704,85.61487411644846,153.3497936960279,74.31904480635694,131.2986461767816,687.2831499692957,23.42390982760392,820.1341314036623,464.41113007274623))
targetgene="Gm9141_Tdpoz2"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(81.4951582796684,107.44253043816182,162.24810743033657,116.54669960368146,160.32023886402916,220.07406802572078,110.15733874153712,4.998422908867605,510.12070291226314,52.91187944539757,100.95894132016222),c(343.2674848749669,290.09483218303694,316.414891502074,195.53332540143714,283.4647701653849,211.42486884567063,405.0229424436314,954.6987755937125,2908.902579702072,471.79759172146163,763.9226559892275),c(414.26705458831435,598.6083838697588,317.65817201878156,418.6843521307609,322.96395945072544,407.4733835934743,475.12306709733684,1946.885723003932,779.0618898217897,2413.3696124817448,599.0230518329626),c(182.12933100380437,72.90743136875267,131.78773477100137,17.12297482328969,128.95323560802348,91.9377838768295,142.4256500900682,444.85963888921685,337.47781196066387,271.90826937218196,33.65298044005407))
targetgene="B3gnt2"
collabel=c("In.vivo_AR1830_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1834_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._SP_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1844_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1845_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1846_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1848_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1843_Pool.B..9.16._SP_mCD19_10m_1410..1411..1412..1413..1414..1415..1416..1417")

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
Sweave("invivo.poolB_SP_mCD19_10m_v_SP_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolB_SP_mCD19_10m_v_SP_hEGFRv3_15m_summary.tex",pdf=TRUE);

