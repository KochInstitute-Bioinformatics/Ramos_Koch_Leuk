pdf(file='invitro.poolD_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolD_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Fbxl14","Mki67","POTENTIALLY_ACTIVE_29","Amigo2","Nucks1","Fbxw21","Stxbp3","Slc6a20a","Mbd2","Ccdc71l")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='11,12,13_vs_5,6,7 neg.'


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
targetmat=list(c(178.59611250082722,156.27353951547443,225.61968978425463,96.92678414971199,61.69810538366181,19.052573517007907),c(164.43910358307872,149.98064530679758,217.52052143302498,4.690005684663483,16.36888510178783,10.717072603316948),c(422.5322661604937,362.8902327003634,351.7353112534021,153.7279641084142,164.9479960257081,216.72302375596496),c(408.3752572427452,431.0632532943623,421.156754263942,203.23357966875096,256.86558159728594,188.14416348045307))
targetgene="Fbxl14"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(1367.7848615917012,1586.8581562880054,1604.792357593647,2451.830749593521,1867.3120466116422,2154.1315932667067),c(2026.6302766099968,2353.5424340451314,2171.7341421797228,1432.536180793324,1444.2393239808184,1284.8579265532208),c(1193.5447518347967,1146.3555616806277,1165.1232185268943,886.9321861441388,804.5936600032633,801.3988735591452),c(203.64312827838225,174.10340644005873,252.2312429382949,51.59006253129832,62.95725039149165,41.6775045684548))
targetgene="Mki67"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(761.211479500477,681.7302059399891,765.9499212162901,815.0187656459653,868.8100554025848,591.8205648720581),c(676.269425993986,889.3957148263241,765.9499212162901,245.9647425734627,329.8959920514162,227.44009635928188))
targetgene="POTENTIALLY_ACTIVE_29"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(229.77914474191795,235.98353282538082,187.4378961284577,100.05345460615432,86.88100554025847,95.26286758503954),c(252.64815914751168,363.93904840180954,325.1237580993618,194.89579178490476,201.46320125277327,163.1376607393802),c(547.7673450482689,639.7775778821435,698.8425263061015,434.60719344548284,433.14588269346257,439.39997673599487),c(733.9864623509607,718.4387554906039,591.2392896397647,326.7370626982227,439.4416077326117,403.67640139160505))
targetgene="Amigo2"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(512.919323096888,570.5557415866986,514.8757023281709,296.51258161928024,307.23138191047923,354.85418175427225),c(693.6934369696764,670.1932332240815,630.5781073457373,380.93268394322297,348.7831671688637,319.13060640988243),c(743.7874685247865,711.0970455804809,686.1152617541692,519.0272957694256,555.2829484529564,342.94632330614235),c(369.17123254744166,353.4508913873482,340.16507075164543,194.37468004216439,167.4662860413678,161.9468748945672))
targetgene="Nucks1"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(384.4172421511708,370.2319426104864,374.8757922569154,249.6125247726454,220.35037637022077,153.61137398087627),c(271.1611708091828,282.13142368901083,253.38826698847058,219.90915543644334,101.99074563421647,65.49322146471468),c(567.3693573959206,539.0912705433144,600.49548204117,502.8728317444735,397.88982247422723,357.23575344389826),c(614.1963868930887,663.9003390154047,562.3136883853731,368.42600211745366,396.6306774663974,294.1241036688096))
targetgene="Fbxw21"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(168.79510632700135,131.10196268076712,149.25610247266076,76.08231444009651,41.551785258384484,65.49322146471468),c(360.4592270595964,426.8679904885778,448.92533146815794,145.390176224568,260.6430166207754,209.57830868708697),c(174.24010975690462,224.44656010947332,217.52052143302498,110.47568946096206,141.0242408769413,103.5983684987305),c(401.84125312686126,342.9627343728868,292.7270846944432,174.05132207528928,182.57602613532578,180.99944841157512))
targetgene="Stxbp3"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(621.8193916949533,579.9950828997138,516.0327263783465,362.172661204569,357.5971822236726,238.15716896259883),c(351.7472215717512,375.47602111771704,395.7022251600774,200.10690921230864,235.46011646417875,182.19023425638812),c(154.63809740925285,170.95695933572034,230.2477859849573,72.95564398365418,85.62186053242864,89.30893836097457),c(798.2375028238192,777.1724347715875,705.7846706071555,441.90275784384824,493.5848430692945,618.0178534579439))
targetgene="Slc6a20a"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(191.66412073259508,224.44656010947332,201.32218473056568,64.0967443570676,70.51212043847065,33.342003654763836),c(291.85218384281524,287.37550219624154,357.5204315042804,105.2645720335582,163.68885101787828,196.47966439414404),c(125.23507888777519,155.22472381402827,148.0990784224851,68.7867500417311,56.66152535234248,90.49972420578756),c(333.23420991008004,409.03812356399345,468.59474032114423,190.72689784298166,135.98766084562197,185.7625917908271))
targetgene="Mbd2"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(570.6363594538626,521.2614036187301,622.4789389945076,356.44043203442476,333.67342707490576,357.23575344389826),c(512.919323096888,613.5571853459901,524.1318947295762,165.7135341914431,201.46320125277327,244.1110981866638),c(274.4281728671248,268.4968195702111,267.2725555905786,145.390176224568,103.2498906420463,71.44715068877966),c(343.03521608390594,357.6461541931327,402.64436946113136,776.4564966831767,672.3834341811308,537.0444160106604))
targetgene="Ccdc71l"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetgenelist=c("Hmbox1","Ints14","Lsm12","Carmil3","Wdhd1","Gfy","Prss50","Nek11","Gm32856","Slc30a7")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='11,12,13_vs_5,6,7 pos.'


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
targetmat=list(c(308.187194132525,328.27931455264087,230.2477859849573,873.3832808328888,805.8528050110931,930.0037447989484),c(521.6313285847332,528.6031135288531,559.9996402850218,1881.2133912927973,1304.474228111707,2051.7240106127892),c(645.7774067865276,660.7538919110663,599.3384579909944,2290.8072210867417,2370.9700497435756,2840.024239878991),c(491.13930937727486,499.2362738883612,490.57819727448185,1811.3844177655856,1469.422224137415,1749.2644060302885))
targetgene="Hmbox1"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(401.84125312686126,467.7718028449771,520.6608225790492,2345.0028423317417,1862.275466580323,1863.579847132336),c(129.5910816316978,105.93038584605983,150.41312652283642,707.1486348987053,1140.7853770938286,783.5370858869502),c(623.9973930669146,793.9534859947257,724.2970554099661,3187.119418600207,2559.8418009180505,2514.9397042450437),c(356.1032243156738,403.79404505676274,401.4873454109557,980.211188094668,822.2216901128809,902.6156703682497))
targetgene="Ints14"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(629.4423964968179,640.8263935835897,697.6855022559258,1910.395648886259,1751.4707058912977,1779.0340521506134),c(270.07217012320217,264.30155676442655,251.07421888811925,778.0198319113979,489.807408045805,791.8725868006411),c(254.82616051947298,202.42143037910444,197.8511125800387,366.3415551464921,514.9903082024017,646.5967137334559),c(361.54822774557704,345.0603657757791,377.18984035726675,885.3688509159176,888.9563755278621,915.7143146611926))
targetgene="Lsm12"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(559.746352594056,734.1709910122959,573.8839288871297,3897.394723955355,3268.7404403262462,2936.4778933088437),c(828.7295220312775,869.4682164988476,710.4127668078581,4635.288951675743,4473.742212819397,4349.940691101868),c(362.63722843155773,443.64904171171594,385.2890087084964,275.1470001669244,275.75275671473344,220.29538129040392),c(300.56418933066044,317.7911575381795,343.6361429021724,210.008032324376,211.53636131541194,180.99944841157512))
targetgene="Carmil3"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(105.63306654012342,124.8090684720903,96.03299616458018,33.87226327812516,23.923755148766826,36.91436118920282),c(138.30308711954302,188.78682626030465,216.36349738284932,76.08231444009651,137.2468058534518,69.06557899915366),c(153.54909672327219,199.27498327476604,171.23955942599838,1988.0412985545768,1764.062155969596,953.8194616952084),c(141.570089177485,169.90814363427418,159.66931892424174,1899.9734140314513,1480.7545292078835,1488.4823060162428))
targetgene="Wdhd1"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(287.4961810988926,184.5915634545201,259.1733872393489,1059.420172991207,868.8100554025848,1194.3582023474332),c(732.89746166498,695.3648100587889,733.5532478113714,1492.985142951209,1261.6632978454927,1380.1207941382602),c(438.86727645020346,447.8443045175005,509.0905820772925,1043.7868207089953,857.4777503321162,648.9782854230818),c(370.2602332334223,388.06180953507067,439.6691390667526,821.27210655885,758.0052947135595,445.3539059600598))
targetgene="Gfy"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(657.7564143323149,789.7582231889412,807.602787022614,787.3998432807249,654.7554040715131,602.537637475375),c(824.373519287355,799.1975645019564,856.197797129992,4037.052671009779,3047.130918948196,3378.2594417344644),c(610.9293848351468,674.3884960298661,639.8342997471426,1007.8301104599086,1033.758051428293,887.1354543856806),c(323.4332037362542,532.7983763346376,440.8261631169283,2938.5491173130426,2474.2199403856216,3755.738554540184))
targetgene="Prss50"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(446.4902812520681,478.2599598594385,448.92533146815794,1884.86117349198,1218.8523675792783,1236.0357069158879),c(673.002423936044,742.561516623865,703.4706225068041,10719.789659912503,11216.463729748151,10552.744156732755),c(503.1183169230621,469.8694342478694,553.0574959839678,2724.893302789484,2153.137963389014,1996.9478617513912),c(815.6615137995097,810.7345372178639,838.842436377357,741.0208981768304,829.77656015986,576.3403488894892))
targetgene="Nek11"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(129.5910816316978,156.27353951547443,165.45443917512006,97.44789589245238,79.32613549327948,107.17072603316947),c(163.35010289709808,167.81051223138192,214.049449282498,522.6750779686082,803.3345149954334,922.8590297300705),c(251.55915846153104,229.690638616704,208.26432903161967,1116.221352949909,1482.0136742157133,1231.272563536636),c(413.82026067264843,370.2319426104864,459.3385479197389,1117.7846881781302,1179.8188723365536,1498.0085927747466))
targetgene="Gm32856"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(1147.806723023609,1233.407264900657,1261.1562146914746,1035.449032825149,985.9105411307592,1026.457398228801),c(601.1283786613209,463.57654003919254,546.1153516829138,2467.985213618473,1265.4407328689822,2076.730513353862),c(556.4793505361141,516.0173251114994,590.082265589589,1445.0428626190933,833.5539951833495,1645.666037531558),c(1743.4900982550269,1759.9127470266178,1589.7510449413635,3831.7346443700662,4114.885885587894,4962.004615335747))
targetgene="Slc30a7"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
Sweave("invitro.poolD_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2_summary.Rnw");
library(tools);

texi2dvi("invitro.poolD_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2_summary.tex",pdf=TRUE);

