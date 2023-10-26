pdf(file='invitro.poolD_hEGFRv3.CAR_ET.1.10_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolD_hEGFRv3.CAR_ET.1.10_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Mtx2","Scd2","Gpr89","Yme1l1","Mtrr","Gtpbp3","Bloc1s1","Gdi2","Aagab","Thoc7")
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
targetmat=list(c(253.37010101462917,189.40077434235198,179.32943329136324,171.8027701780829),c(1018.1379426801092,542.5415729548663,591.9612361074127,502.03444219703266),c(360.4934893112555,138.07927419797272,134.9323405833073,132.70992426929757),c(802.959658362538,497.32977520862744,422.2076463413164,411.50364114510876))
targetgene="Mtx2"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(318.5756417169234,146.63285755536927,167.1419960773871,201.6367841611033),c(202.13717617711225,80.6480716554531,114.03959107363391,97.73211477196332),c(182.57551396642396,139.3012146776008,133.19127812416784,146.0837926065136),c(338.13730392761175,100.19911932950234,90.53524787525134,124.47985144639539))
targetgene="Scd2"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(212.38376114461562,111.19658364615503,118.39224722148253,96.70335566910055),c(340.0003193762487,149.07673851462542,209.7980263263036,171.8027701780829),c(284.1098559171393,206.50794105714505,235.04343198382563,189.29167492675),c(310.192072198057,229.72481017007854,264.64149378919626,195.46422954392665))
targetgene="Gpr89"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(362.35650475989246,178.40331002569928,216.7622761628614,236.6145936584375),c(328.8222266844268,212.61764345528545,219.3738698515706,209.86685698400547),c(156.49329768550626,120.97210748317966,134.06180935373757,123.45109234353262),c(382.84967469489925,212.61764345528545,239.39608813167425,248.95970289279077))
targetgene="Yme1l1"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(430.3565686351422,294.4876555903666,258.54777518220817,280.8512350815367),c(336.27428847897477,227.28092921082236,234.1729007542559,257.189775715693),c(294.3564408846427,249.27585784412776,208.92749509673388,220.15444801263317),c(210.52074569597866,107.5307622072708,125.35649705804033,142.99751529792528))
targetgene="Mtrr"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(255.23311646326616,205.28600057751697,175.84730837308436,186.2053976181617),c(543.0690032776795,287.15601271259817,330.8018672364953,318.91532188745924),c(148.10972816663985,65.98478589991618,55.71399869246237,108.01970580059104),c(66.13704842661276,17.10716671479308,20.022218280103665,27.77649577729484))
targetgene="Gtpbp3"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(267.3427168794065,168.62778618867466,187.16421435749078,147.11255170937636),c(279.4523172955469,118.5282265239235,172.36518345480545,160.4864200465924),c(110.84941919390026,86.7577740535935,81.8299355795541,61.72554617176631),c(540.2744801047239,318.92646518292815,346.47142936875036,402.2448092193438))
targetgene="Bloc1s1"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(547.7265418992719,394.686774919869,326.4492110886467,446.48145064244295),c(216.1097920418896,109.97464316652696,127.9680907467495,129.62364696070924),c(150.9042513395953,94.08941693136195,134.06180935373757,140.93999709219975),c(346.5208734464781,172.29360762755888,195.86952665318802,217.06817070404486))
targetgene="Gdi2"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(272.00025550099895,118.5282265239235,96.62896648223942,175.917806589534),c(371.67158200307733,250.49779832375583,255.936181493499,291.1388261101644),c(54.0274480104724,50.09955966475117,47.87921762633485,36.00656860019701),c(336.27428847897477,230.94675064970662,212.4096200150128,189.29167492675))
targetgene="Aagab"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(340.9318271005672,301.8192984681351,242.00768182038343,254.10349840710464),c(331.6167498573823,254.16361976264008,261.1593688709174,341.54802215044026),c(162.0823440314172,114.86240508503927,103.59321631879722,112.13474221204213),c(313.918103095331,213.83958393491352,195.86952665318802,273.6499213614973))
targetgene="Thoc7"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetgenelist=c("Lztr1","Rbp7","Gm17455","Ptbp3","E2f7","1700019D03Rik","Cdk17","B020004C17Rik","Zmym4","Myef2")
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
targetmat=list(c(352.10991979238906,608.5263588547825,519.7071440531255,561.7024701630734),c(402.4113369055875,562.0926206289156,550.175737088066,600.7953160718588),c(364.21952020852945,432.56692978833934,478.7921762633485,549.3573609287201),c(232.8769310796224,320.14840566255623,226.33811968812836,300.3976580359294))
targetgene="Lztr1"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.9315077243184896,10.997464316652696,6.964249836557796,3.0862773085883157),c(3.7260308972739584,15.885226235165005,16.540093361824766,27.77649577729484),c(4.657538621592448,15.885226235165005,7.834781066127521,18.517663851529893),c(0.9315077243184896,8.55358335739654,4.352656147848623,0.0))
targetgene="Rbp7"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.9315077243184896,15.885226235165005,12.187437213976143,6.1725546171766315),c(0.9315077243184896,6.109702398140387,17.41062459139449,8.230072822902175),c(5.589046345910938,0.0,7.834781066127521,7.201313720039403),c(7.452061794547917,0.0,4.352656147848623,14.402627440078806))
targetgene="Gm17455"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(488.11004754288854,585.309489741849,490.9796134773246,564.7887474716617),c(593.3704203908778,679.398906673211,713.8356082471742,748.9366268840979),c(635.2882679852099,716.0571210620533,721.6703893133016,622.399257231977),c(530.0278951372206,657.4039780399056,654.6394846364328,624.4567754377025))
targetgene="Ptbp3"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(485.3155243699331,505.883358566024,574.5506115160182,495.86188757985605),c(362.35650475989246,549.8732158326347,589.3496424187035,408.4173638365204),c(490.90457071584405,683.0647281120952,594.5728297961218,640.9169210835068))
targetgene="E2f7"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.9315077243184896,13.44134527590885,5.223187377418347,13.373868337216035),c(2.794523172955469,6.109702398140387,3.482124918278898,3.0862773085883157),c(0.0,6.109702398140387,4.352656147848623,3.0862773085883157),c(8.383569518866407,7.331642877768464,12.187437213976143,4.115036411451087))
targetgene="1700019D03Rik"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(293.4249331603242,301.8192984681351,356.9178041235871,285.99503059585055),c(446.1921999485565,577.9778468640806,537.11776864452,519.5233469456998),c(416.38395277036483,609.7482993344106,544.0820184810779,532.8972152829158),c(530.0278951372206,711.169359143541,662.4742657025604,758.1954588098629))
targetgene="Cdk17"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.9315077243184896,7.331642877768464,5.223187377418347,9.258831925764946),c(3.7260308972739584,1.2219404796280773,8.705312295697246,1.0287591028627718),c(4.657538621592448,1.2219404796280773,3.482124918278898,9.258831925764946),c(3.7260308972739584,15.885226235165005,7.834781066127521,17.48890474866712))
targetgene="B020004C17Rik"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.9315077243184896,23.21686911293347,2.6115936887091733,12.345109234353263),c(2.794523172955469,6.109702398140387,1.741062459139449,4.115036411451087),c(1.8630154486369792,2.4438809592561546,0.0,4.115036411451087),c(2.794523172955469,3.665821438884232,6.964249836557796,0.0))
targetgene="Zmym4"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(50.30141711319844,41.54597630735463,80.95940434998438,74.07065540611957),c(272.00025550099895,311.5948223051597,349.0830230574595,350.8068540762052),c(484.3840166456146,596.3069540585017,548.4346746289265,578.1626158088777))
targetgene="Myef2"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
Sweave("invitro.poolD_hEGFRv3.CAR_ET.1.10_v_input_summary.Rnw");
library(tools);

texi2dvi("invitro.poolD_hEGFRv3.CAR_ET.1.10_v_input_summary.tex",pdf=TRUE);

