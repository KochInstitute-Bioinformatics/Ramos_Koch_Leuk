pdf(file='invitro.poolF_hEGFRv3.CAR_ET.1.10_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolF_hEGFRv3.CAR_ET.1.10_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Dph3","Znhit1","Vrk1","Aifm1","Fastkd5","Mtg1","Slc25a32","Slc25a19","Commd4","Ubxn7")
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
targetmat=list(c(444.2075613884075,266.46958310212415,235.12073242087027,229.72854235949697),c(1465.697523230948,854.270134062692,790.1598384635804,827.8149198816357),c(337.37283143423355,204.7505252512645,150.32309121990068,197.47601301345924),c(314.881309338618,133.23479155106207,162.85001548822572,203.70018534339633))
targetgene="Dph3"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(1313.879749085543,1044.3256455399423,863.3941649553269,865.1599538612584),c(2020.4884015894654,1610.573509631956,1520.5758842628413,1465.5096667761015),c(790.9518603624808,476.11844627806,524.2036001514484,441.35040157735875),c(1244.5308892907283,997.3016014630969,958.7915113064178,884.9641385474218))
targetgene="Znhit1"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(762.8374577429614,570.1665344317508,561.7843729564237,475.8662663161009),c(864.0493071732315,632.8652598675448,577.2021259020545,565.2680143279248),c(223.04092744818774,110.70243709757362,110.81509929672164,123.35177890239001),c(470.447670499959,252.7542369130442,272.70150522584544,228.59687466314477))
targetgene="Vrk1"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(507.9335406593183,381.0906905394349,331.481688331063,357.6069920472958),c(434.83609384856766,257.6525748377156,271.7378956667435,303.8527764705662),c(421.71603929279195,271.3679210267955,239.9387802163799,290.2727641143398),c(352.36717949797725,210.62853076087018,238.975170657278,221.80686848503157))
targetgene="Aifm1"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(772.2089252828013,561.3495261673423,512.6402854422254,627.509737627296),c(206.17228587647605,151.84847566481338,117.56036621043513,111.46926809069188),c(601.6482160577165,306.6359540844296,370.0260706951401,405.13703529408826),c(290.51549373503445,177.31983287310464,186.9402544657739,172.01348984553468))
targetgene="Fastkd5"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(532.2993562629018,333.08697887765516,380.62577584526133,387.030352152453),c(734.7230551234419,501.5898034863513,422.0609868866442,406.2687029904405),c(376.7329951015608,218.46587144034442,223.55741771164713,261.41523785735865),c(253.02962357567515,109.72276951263935,128.16007136055634,144.28763128490573))
targetgene="Mtg1"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(708.4829460118905,564.2885289221452,419.1701582093384,497.3679525467927),c(959.6382760795976,698.5029880581416,717.8891215309359,635.4314115017613),c(172.4350027330527,96.00742332355942,107.92427061941586,95.06008649358496),c(376.7329951015608,255.69323966784705,253.4293140438069,326.4861303976102))
targetgene="Slc25a32"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(226.78951446412367,119.51944536198214,149.35948166079874,173.710991390063),c(704.7343589959545,388.9280312189091,450.0056641006001,393.25452448239014),c(631.6369121852039,266.46958310212415,306.4278397944129,265.94190864276743),c(309.2584288147141,156.7468135894848,182.12220667026426,162.39431442654097))
targetgene="Slc25a19"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(682.242836900339,594.6582240551079,561.7843729564237,511.5137987511952),c(479.8191380397988,296.8392782350868,291.973696407884,262.5469055537108),c(571.6595199302291,317.4122975187067,326.6636405355534,375.14784134075495),c(581.0309874700689,262.550912762387,310.28227803082063,368.3578351626417))
targetgene="Commd4"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(462.9504964680871,263.5305803473213,315.1003258263302,235.386880841258),c(1223.9136607030805,1064.8986648235623,1063.824953248528,1053.5826253039),c(406.7216912290482,193.97418181698743,247.64765668919534,245.5718901084278),c(890.289416284783,416.35872359706894,463.4961979280271,432.29706000654113))
targetgene="Ubxn7"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetgenelist=c("Cyrib","Pum2","Gm46735","Dnajc12","Chd6","Arl15","Kcns3","4930467E23Rik_Gm15319_Gm20778_Gm21119","Btbd35f1_Btbd35f10_Btbd35f11_Btbd35f12_Btbd35f13_Btbd35f14_Btbd35f15_Btbd35f16_Btbd35f17_Btbd35f18_Btbd35f19_Btbd35f2_Btbd35f20_Btbd35f21_Btbd35f22_Btbd35f23_Btbd35f26_Btbd35f27_Btbd35f3_Btbd35f4_Btbd35f5_Btbd35f6_Btbd35f7_Btbd35f8_Btbd35f9_Gm14367_Gm5925","G3bp1")
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
targetmat=list(c(517.3050081991581,694.5843177184045,651.4000619529029,649.0114238579878),c(603.5225095656845,917.9485270834202,889.411623051079,849.3166061123275),c(1177.0563230038815,1274.547527999498,1458.904872480318,1586.598110285787),c(1002.7470267628609,1547.8747841961622,1585.1377247226706,1309.9053585276736))
targetgene="Cyrib"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(1799.3217676492457,2212.089406781604,2280.863826394262,2460.245571869687),c(787.203273346545,951.2572249711858,847.0128024505941,742.3740088070444),c(2339.118297944019,3080.074887033376,2912.9916971651264,2920.834324285033),c(1283.8910529580555,1115.8413792401448,1219.92970182304,1340.4603863291832))
targetgene="Pum2"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(566.0366394063252,539.7968392987882,507.82223764671573,482.6562724942141),c(350.4928859900093,524.1221579398397,539.6213530970793,510.38213105484306),c(444.2075613884075,752.384705229527,642.7275759209856,691.4489624711954),c(1593.1494817727696,1639.9635371799845,1482.9951114578662,1469.4705037133342))
targetgene="Gm46735"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(1128.3246917967144,1455.7860312123398,1327.8539724424559,1532.8438947090574),c(1047.7300709540918,1253.9745087158783,1298.945685669398,1162.7885580018874),c(674.7456628684671,604.4548999044507,678.3811296077569,583.9405313177361),c(1774.9559520456621,2614.732784189593,2298.2087984580967,2135.456943016605))
targetgene="Dnajc12"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(1205.170725623401,1270.628857659761,1447.3415577710948,1497.7621961221391),c(1570.657959677154,1951.4978291890855,1989.85373954548,1933.454259217737),c(659.7513148047234,799.4087493063723,820.0317347957401,811.9715721327049),c(734.7230551234419,838.5954527037436,722.7071693264455,815.9324090699375))
targetgene="Chd6"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(189.3036443047644,235.12022038422717,265.956238312132,255.19106552742153),c(678.494249884403,1016.8949531617825,921.2107385014425,860.6332830758495),c(1823.6875832528292,2431.5349458068827,2428.296088936857,2247.492044955473),c(749.7174031871856,604.4548999044507,631.1642612117624,617.3247283601261))
targetgene="Arl15"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(558.5394653744534,641.6822681319533,693.7988825533877,631.4705745645286),c(993.375559223021,1417.5789953999029,1084.0607539896682,1219.9377766676737),c(721.6030005676662,868.9651478367063,830.6314399458614,852.145775353208),c(905.2837643485267,1049.2239834646136,924.1015671787484,1028.1201021359755))
targetgene="Kcns3"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(29.988696127487426,94.04808815369087,93.47012723288695,72.99256641471702))
targetgene="4930467E23Rik_Gm15319_Gm20778_Gm21119"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(260.526797607547,612.2922405839249,377.7349471679555,462.28625395987444))
targetgene="Btbd35f1_Btbd35f10_Btbd35f11_Btbd35f12_Btbd35f13_Btbd35f14_Btbd35f15_Btbd35f16_Btbd35f17_Btbd35f18_Btbd35f19_Btbd35f2_Btbd35f20_Btbd35f21_Btbd35f22_Btbd35f23_Btbd35f26_Btbd35f27_Btbd35f3_Btbd35f4_Btbd35f5_Btbd35f6_Btbd35f7_Btbd35f8_Btbd35f9_Gm14367_Gm5925"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(311.13272232268207,369.3346795202235,411.461281736523,378.5428444298115),c(875.2950682210393,1068.8173351632993,1034.91666647547,999.2625758789943),c(4593.89338802948,6263.994538069785,5845.255585512291,5497.641668878997),c(369.23582106968894,374.2330174448949,487.5864369055752,371.75283825169834))
targetgene="G3bp1"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
Sweave("invitro.poolF_hEGFRv3.CAR_ET.1.10_v_input_summary.Rnw");
library(tools);

texi2dvi("invitro.poolF_hEGFRv3.CAR_ET.1.10_v_input_summary.tex",pdf=TRUE);

