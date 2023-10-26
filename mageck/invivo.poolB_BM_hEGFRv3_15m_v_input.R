pdf(file='invivo.poolB_BM_hEGFRv3_15m_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolB_BM_hEGFRv3_15m_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Hmbs","Hsd17b10","Cdipt","Ppcs","Sec61a1","Coa6","Slc25a3","Atic","Timm10","Mtif2")
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
targetmat=list(c(426.90305402700096,50.0840346903402,69.27410038862683,43.249683074325354,42.01051472282916,31.03585640318409),c(219.9459953810512,13.912231858427832,103.50365587477185,2.790302133827442,137.8470014342832,5.431274870557216),c(213.88449157133718,0.0,17.114777743072512,0.0,13.128285850884112,6.2071712806368184),c(221.67785361239805,309.77902938099305,0.8149894163367862,8.370906401482326,0.0,5.431274870557216))
targetgene="Hmbs"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(338.5782842283111,38.95424920359793,5.704925914357504,29.29817240518814,148.34963011499048,386.39641221964195),c(361.95837035149367,7.41985699116151,37.48951315149217,318.0944432563284,42.01051472282916,11.638446151194035),c(290.95218286627244,5.564892743371133,413.1996340827506,2.790302133827442,11.8154572657957,1.5517928201592046),c(561.1220669563826,95.53065876120445,332.5156818654088,15.346661736050931,36.75920038247551,5.431274870557216))
targetgene="Hsd17b10"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(162.79467374660481,12.057267610637455,15.484798910398938,1.395151066913721,15.753943021060934,306.4790819814429),c(142.0123749704425,17.62216035400859,4.889936498020718,0.0,17.066771606149345,35.691234863661705),c(296.14775756031304,50.0840346903402,22.819703657430015,227.40962390693653,22.31808594650299,24.05278871246767),c(310.0026234110879,4.637410619475944,35.859534318818596,30.693323472101863,22.31808594650299,242.8555763549155))
targetgene="Cdipt"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(159.3309572839111,1.8549642477903776,47.2693861475336,1.395151066913721,7.876971510530467,1.5517928201592046),c(307.40483606406764,20.404606725694155,11.409851828715007,18.13696386987837,82.7082008605699,16.29382461167165),c(204.35927129892946,48.22907044254982,4.889936498020718,4.185453200741163,127.34437275357588,1.5517928201592046),c(226.0074991907652,39.88173132749312,34.229555486145024,27.90302133827442,81.39537227548149,11.638446151194035))
targetgene="Ppcs"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(104.77742299648502,4.637410619475944,0.0,16.741812802964652,213.991059369411,32.5876492233433),c(103.04556476513815,78.83598053109105,4.889936498020718,0.0,10.50262868070729,107.84960100106473),c(79.66547864195554,0.0,6.51991533069429,9.766057468396047,14.441114435972523,228.11354456340308),c(109.9729976905256,4.637410619475944,41.5644602331761,0.0,39.384857552652335,6.983067690716421))
targetgene="Sec61a1"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(435.56234518373526,10.202303362847077,22.00471424109323,193.9259983010072,65.64142925442056,8.534860510875625),c(158.46502816823767,8.3473391150567,8.964883579704649,44.64483414123907,410.9153471326727,20.17330666206966),c(107.3752103435053,11.129785486742266,9.779872996041435,18.13696386987837,149.66245870007887,22.500995892308467),c(48.49203047771208,4.637410619475944,107.57860295645578,12.55635960222349,7.876971510530467,3.8794820503980114))
targetgene="Coa6"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(275.3654587841507,11.129785486742266,54.60429089456468,131.14420028988977,53.82597198862486,139.6613538143284),c(157.59909905256424,9.274821238951889,42.379449649512885,27.90302133827442,7.876971510530467,16.29382461167165),c(125.55972177264734,22.259570973484532,15.484798910398938,13.95151066913721,27.569400286856634,10.862549741114432),c(247.65572708260095,600.0809341601872,101.87367704209828,15.346661736050931,31.507886042121868,19.397410251990056))
targetgene="Slc25a3"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(268.43802585876324,2.7824463716855665,53.78930147822789,768.7282378694603,31.507886042121868,524.5059732138111),c(244.19201061990722,75.1260520355103,441.72426365453816,15.346661736050931,22.31808594650299,19.397410251990056),c(252.85130177664152,61.21382017708246,25.264671906440373,50.22543840889396,30.195057457033457,24.05278871246767),c(187.90661810113428,76.05353415940549,63.569174474269325,1.395151066913721,22.31808594650299,23.27689230238807))
targetgene="Atic"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(144.6101623174628,10.202303362847077,5.704925914357504,0.0,368.90483240984355,82.24501946843785),c(134.21901292938162,80.69094477888143,145.88310552428473,1.395151066913721,6.564142925442056,3.1035856403184092),c(89.19069891436328,0.0,0.0,0.0,9.189800095618878,532.2649373146072),c(166.25839020929854,127.06505097364087,44.824417898523244,62.781798011117445,171.98054464658185,13.966135381432842))
targetgene="Timm10"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(155.00131170554394,194.77124601798965,14.669809494062152,1420.2637861181681,15.753943021060934,38.79482050398011),c(90.0566280300367,12.984749734532643,19.55974599208287,0.0,6.564142925442056,0.0),c(471.06543892634585,5.564892743371133,17.114777743072512,432.4968307432535,68.26708642459738,42.67430255437813),c(544.6694137585874,34.31683858412199,79.05397338466827,37.66907880667047,340.0226035378985,20.17330666206966))
targetgene="Mtif2"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetgenelist=c("Cnot8","Cnot6l","Coro1a","Ap2s1","Edc3","Rasal1","Zfp36l2","Dcp1a","Slc44a4","E130309D02Rik")
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
targetmat=list(c(333.38270953427053,2896.5266729246746,3969.813446976486,2437.3289138982705,1596.399559467508,2215.9601471873443),c(240.7282941572135,1649.0632162856457,1533.8100815458317,1799.7448763187,1110.652982984796,1807.8386354854733),c(166.25839020929854,1038.7799787626116,210.26726941489085,3927.350253362125,1375.8443571726548,1246.0896345878414),c(223.4097118437449,1921.7429607108313,931.5329028729467,867.7839636203345,3830.8338112879837,894.6085608217815))
targetgene="Cnot8"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(643.3853329453584,470.23343681486074,1632.4238009225828,2947.9542043886927,647.2244924485867,2884.7828526759613),c(459.80836042259125,1468.2042021260838,3484.0797548397613,2243.4029155972635,1147.4121833672714,1158.4133402488462),c(491.8477377025082,1048.9822821254586,1781.5668641122147,2268.5156348017103,1164.4789549734207,1338.421307387314),c(596.6251606989932,980.3486049572145,2992.641136788679,739.4300654642722,1591.1482451271543,1322.9033791857219))
targetgene="Cnot6l"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(269.3039549744367,336.6760109739535,1306.4280343878684,361.34412633065375,204.80125927379214,285.52987890929364),c(347.2375753850454,137.26735433648795,1041.5564740784127,665.4870589178449,1184.1713837497468,1125.825691025503),c(232.0690030004792,256.9125483189673,992.6571090982056,671.0676631854998,833.6461515311411,277.77091480849765),c(155.00131170554394,271.75226230129033,138.54820077725367,570.6167863677119,945.236581263656,717.7041793236322))
targetgene="Coro1a"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(268.43802585876324,492.49300778834527,489.80863921840853,212.0629621708856,501.50051950377303,137.3336645840896),c(160.19688639958454,1184.3946722141561,1062.7461989031692,1471.8843755939756,456.8643476107671,616.0617496032042),c(189.63847633248113,1183.4671900902608,280.3563592198545,542.7137650294375,779.8201795425163,1170.0517864000403),c(162.79467374660481,32.46187433633161,180.92765042676655,22.322417070619537,102.40062963689607,27.15637435278608))
targetgene="Ap2s1"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(252.85130177664152,423.8593306201013,142.62314785893759,150.67631522668188,1442.798615012164,282.42629326897526),c(897.9684929533468,1918.0330322152504,2189.8765616969445,3060.961440808704,1643.6613885306908,3700.2499796696234),c(248.52165619827437,601.0084162840824,215.97219532924836,48.83028734198024,410.9153471326727,464.76194963768177),c(90.92255714571014,3007.8245277920973,523.2232052882167,991.9524085756557,655.1014639591172,20.17330666206966))
targetgene="Edc3"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(193.9681219108483,186.42390690293294,642.2116600733875,1231.9183920848157,422.7308043984684,856.5896367278809),c(208.6889168772966,692.829146549706,807.6545115897552,2131.7908302441656,286.19663154927366,209.49203072149263),c(131.62122558236135,29.67942796464604,8.149894163367863,34.878776672843024,475.24394780200487,20.17330666206966),c(146.34202054880964,119.64519398247936,20.374735408419657,2022.9690470248954,246.8117739966213,45.77788819469654))
targetgene="Rasal1"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(287.48846640357874,2203.6975263749687,342.2955548614502,87.89451721556442,3602.4016374826,808.4840593029456),c(298.74554490733334,366.3554389385996,234.71695190499443,172.99873229730142,76.14405793512785,2444.0736917507475),c(219.08006626537775,4120.803076466324,552.5628242763411,256.7077963121247,224.4936880501183,55.86454152573137),c(439.8919907621024,608.4282732752439,1285.238309563112,248.33688991064236,1318.0798994287647,1186.345611011712))
targetgene="Zfp36l2"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(781.0680623374337,1756.6511426574875,1771.7869911161733,3353.9431648605855,691.8606643415927,505.1085629618211),c(407.85261348218546,1380.093400356041,303.1760628772845,689.2046270553782,422.7308043984684,363.11951991725385),c(743.8331103634763,1067.5319246033623,969.8374054407756,2367.571360552585,1097.5246971339118,728.5667290647466),c(266.7061676274164,178.07656778787626,1036.666537580392,2123.4199238426836,1059.4526681663478,1580.5009873321499))
targetgene="Dcp1a"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(174.0517522503594,66.7787129204536,35.859534318818596,270.65930698126186,195.61145917817328,382.51693016924395),c(276.23138789982414,242.0728343366443,317.0308829550099,842.6712444158875,957.0520385294517,765.8097567485675),c(277.09731701549754,1189.032082833632,324.3657877020409,2537.7797907160584,1492.6861012455236,1128.9292766658214),c(624.3348924005429,295.8667975225652,1448.2361928304692,385.061694468187,1837.9600191237755,1008.665333103483))
targetgene="Slc44a4"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(317.7959854521488,1178.829779470785,212.71223766390122,1104.959644995667,842.83595162676,286.30577531937325),c(253.71723089231497,810.619376284395,411.56965525007706,918.0094020292285,1419.1677004805724,204.83665226101502),c(219.08006626537775,511.04265026624904,781.574850266978,103.24117895161535,1239.3101843234601,138.10956099416921),c(629.5304670945835,2530.171233986075,547.6728877783204,542.7137650294375,5798.763860335512,796.069716741672))
targetgene="E130309D02Rik"
collabel=c("In.vivo.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1830_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1831_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1832_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1833_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417","In.vivo_AR1835_Pool.B..9.16._BM_hEGFRv3_15m_1410..1411..1412..1413..1414..1415..1416..1417")

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
Sweave("invivo.poolB_BM_hEGFRv3_15m_v_input_summary.Rnw");
library(tools);

texi2dvi("invivo.poolB_BM_hEGFRv3_15m_v_input_summary.tex",pdf=TRUE);

