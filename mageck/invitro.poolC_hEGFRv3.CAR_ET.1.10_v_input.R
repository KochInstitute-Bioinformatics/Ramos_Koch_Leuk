pdf(file='invitro.poolC_hEGFRv3.CAR_ET.1.10_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolC_hEGFRv3.CAR_ET.1.10_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Mto1","Sepsecs","Ctc1","Tmem208","Polg","Tubd1","Mrm2","Pot1b","Dcps","Trim28")
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
targetmat=list(c(426.1814714395743,118.14685428164505,236.824235280007,116.90574589236789),c(410.82358057688697,97.56065997499478,83.72573974545702,88.65352396837898),c(362.83017163098896,119.04190620802115,84.9218217418207,81.83402212465752),c(158.37824952146343,34.907025128667854,37.07854188727382,61.37551659349314))
targetgene="Mto1"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(326.35518083210644,133.36273703003874,209.31434936364255,145.1579678163568),c(232.28809929814636,111.88149079701236,136.35334758545858,113.00888769595562),c(262.0440128446031,132.46768510366263,156.686741523641,138.33846597263533),c(765.0149385976142,422.46450924951864,416.2365347345578,482.2362018060175))
targetgene="Sepsecs"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(431.9406805130821,234.5036047105379,247.58897324728005,221.1467026463959),c(436.74002140767186,211.23225462475932,312.1774010509183,204.5850553116438),c(499.13145303733927,331.1692127591566,367.1971728836472,290.3159356327136),c(331.15452172669626,190.64606031810905,208.1182673672789,212.37877170446833))
targetgene="Ctc1"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(908.0352972563903,626.5363484632692,647.0803600327464,604.9872349930038),c(293.7196627488958,209.44215077200712,147.11808555273163,180.22969158406715),c(332.11438990561425,117.25180235526895,145.92200355636794,116.90574589236789),c(417.54265782931265,183.48564490710027,141.13767557091327,178.28126248586102))
targetgene="Tmem208"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(502.0110575740932,304.31765496787364,289.45184312000856,226.01777539191124),c(851.4030747002307,548.6668308685486,534.6486523745613,548.482791145026),c(253.4051992343415,153.0538794103129,162.66715150545934,177.30704793675795),c(308.1176854326652,136.94294473554314,154.29457753091364,137.36425142353227))
targetgene="Polg"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(227.48875840355657,89.50519263760988,95.68655970909374,71.1176620845238),c(172.77627220523283,91.29529649036209,212.90259535273356,84.75666577196672),c(274.5222991705366,147.6835678520563,75.35316577091132,158.7969715037997),c(478.9742212800621,196.01637187636564,196.15744740364215,208.48191350805607))
targetgene="Tubd1"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(205.41179028844348,154.8439832630651,135.15726558909492,147.10639691456294),c(335.9538626212861,166.4796583059544,185.39270943636913,166.59068789662425),c(761.1754658819424,407.248626501125,374.37366486182924,377.9952450519895),c(539.4459165518936,179.01038527521976,202.13785738546053,166.59068789662425))
targetgene="Mrm2"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(785.1721703548914,548.6668308685486,541.8251443527433,522.1789983192432),c(404.10450332446123,217.49761810939202,258.3537112145531,193.86869527151006),c(474.1748803854723,189.75100839173297,242.80464526182536,185.1007643295825),c(553.843939235663,681.1345159722113,598.0409981818359,624.4715259750651))
targetgene="Pot1b"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(137.2611495852683,74.28930988921621,39.47070588000117,84.75666577196672),c(310.0374217905011,176.32522949609148,122.00036362909452,180.22969158406715),c(361.87030345207097,234.5036047105379,229.64774330182496,226.01777539191124),c(202.53218575168958,102.03591960687527,110.0395436654578,66.24658933900847))
targetgene="Dcps"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(274.5222991705366,148.57861977843243,143.52983956364062,141.26110961994453),c(359.95056709423505,182.59059298072418,110.0395436654578,236.73413543204498),c(490.4926394270776,383.9772764153464,388.7266488181933,342.9235212842791),c(287.0005854964701,171.84996986421098,135.15726558909492,140.28689507084147))
targetgene="Trim28"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetgenelist=c("Kctd5","Ddx6","Hexim1","F8a","Arid4a","Dip2b","Slc39a5","Ddx6_Gm13410","Pitpnm2","Elf4")
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
targetmat=list(c(503.9307939319291,778.695175947206,736.7865097600218,720.9187663362686),c(517.3689484367806,733.0475277020249,784.6297896145686,852.4377304651825),c(603.7570845393969,1015.8839364368722,1107.57192863276,924.5296070988094),c(846.6037338056409,1175.2031793318179,1194.885914367308,1280.1179175214284))
targetgene="Kctd5"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(399.3051624298714,640.8571792852867,615.9822281272909,619.6004532295498),c(568.2419619194325,827.0279799715154,929.3557111745729,726.764053630887),c(213.09073571978715,326.6939531272761,355.2363529200105,257.1926409632093))
targetgene="Ddx6"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(564.4024892037605,732.1524757756489,583.6880142254718,806.6496466573384),c(397.3854260720355,566.5678693960706,687.7471479091113,427.68018705624587),c(354.1913580207273,391.1376918263552,454.5111586181953,443.2676198418949),c(506.810398468683,617.5858291995082,478.4327985454687,681.9501843721459))
targetgene="Hexim1"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(235.16770383490024,311.4780703788824,343.2755329563738,371.17574320826805),c(502.0110575740932,603.2649983774907,673.3941639527471,583.5545149127364),c(325.3953126531885,357.12571862406344,448.5307486363769,463.7261253730593),c(185.2545585311663,274.78094139746236,312.1774010509183,246.47628092307562))
targetgene="F8a"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(393.54595335636367,512.8647538135046,517.9035044254699,544.5859329486137),c(337.873598979122,466.3220536419475,369.58933687637455,574.7865839708088),c(358.03083073639914,436.78534007153627,418.6286987272851,319.54237210580555),c(413.7031851136408,516.4449615190091,571.727194261835,551.4054347923352))
targetgene="Arid4a"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(412.7433169347229,598.7897387456102,486.8053725200144,525.1016419665524),c(342.67293987371175,434.10018429240796,528.6682423927429,371.17574320826805),c(322.5157081164346,362.49603018232006,392.3148948072843,402.35060877956613),c(383.9472715671841,473.4824690529563,588.4723422109265,403.3248233286692))
targetgene="Dip2b"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(371.4689852412506,569.2530251751989,539.4329803600159,493.9267763952543),c(371.4689852412506,413.5139899857577,537.0408163672886,551.4054347923352),c(333.0742580845322,293.5770318513604,399.4913867854664,341.94930673517604),c(366.6696443466608,361.60097825594397,361.2167629018289,349.7430231280006))
targetgene="Slc39a5"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(649.830757127459,978.2917555290761,948.4930231163917,1032.6674220492496))
targetgene="Ddx6_Gm13410"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(2963.113068319743,2967.992187863144,3256.931276098278,3234.392303022178),c(423.30186690282045,528.0806365618984,534.6486523745613,538.7406456539953),c(601.837348181561,577.3084925125838,638.7077860582007,650.7753188008479),c(541.3656529097295,705.3009179843659,522.6878324109246,673.1822534302185))
targetgene="Pitpnm2"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(714.1419251149624,870.8855243639442,846.8260534254796,1027.7963493037344),c(407.94397604013307,597.894686819234,642.2960320472918,604.0130204439007),c(443.4590986600976,547.7717789421725,531.0604063854703,561.1475802833659),c(349.3920171261375,398.298107237364,325.3343030109187,392.60846328853546))
targetgene="Elf4"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
Sweave("invitro.poolC_hEGFRv3.CAR_ET.1.10_v_input_summary.Rnw");
library(tools);

texi2dvi("invitro.poolC_hEGFRv3.CAR_ET.1.10_v_input_summary.tex",pdf=TRUE);

