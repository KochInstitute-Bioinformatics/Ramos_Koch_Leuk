pdf(file='invitro.poolE_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolE_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Usp18","Klhl8","Znfx1","Zcchc10","Ehf","Meis2","Tnfrsf23","Scgb2b4-ps","Gm6829_Golga7","Vmn1r233")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='10,11,12_vs_4,5,6 neg.'


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
targetmat=list(c(339.88751361604756,310.572221007204,409.876730468762,37.00805696216725,6.6726606450993255,31.887000941759677),c(228.28265839883792,269.96737232428336,269.347565736615,11.775290851598673,22.877693640340546,104.10167954515659),c(427.14221860404786,417.02277025702307,412.21888321429776,1195.192021437265,1445.107648281511,887.2089085560192))
targetgene="Usp18"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(357.13553669507087,365.443638146286,374.74443928572526,122.79946173810043,124.8740777868588,143.49150423791855),c(422.0692706396292,380.8076349452289,357.1782936942069,143.82676683024093,142.0323480171142,218.5197417479413),c(292.2018027505125,283.136512437663,282.22940583706185,1497.9852147640881,957.0501839542461,930.3501451242823),c(617.885062066188,623.3392986999713,632.3812412946614,347.37108012216083,304.0826779695264,148.18076908229497))
targetgene="Klhl8"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(758.9130154770256,832.9481121712645,837.3196065290424,587.0823581725624,563.3632058933858,426.7231008382545),c(244.5160918849775,217.29081187076463,278.71617671875816,50.46553222113717,35.26977769552501,55.33332516364179),c(237.41396473479142,321.5465044350204,292.76909319197284,98.40778783121748,104.85609585156082,201.6383883081862),c(614.8412932875368,671.6261457823634,681.5664489509128,277.5604272162544,342.212167370094,246.65533081419986))
targetgene="Znfx1"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(799.4965991923746,755.030699833768,805.7005444643094,352.4176333442745,576.7085271835846,503.6270442860278),c(778.1902177418164,982.1983667895674,922.8081817410985,291.85899467890994,354.60425142527845,283.23159660033593),c(289.15803397186136,200.82938672904004,234.2152745535783,142.9856746265553,159.19061824736963,227.89827143669416),c(859.3573851725143,762.7126982332395,916.952799877259,643.4355358194989,709.2085028505569,673.3784316524543))
targetgene="Zcchc10"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(404.8212475606059,347.88478466177975,510.5892985268007,206.06758990297675,179.2086001826676,118.16947407828586),c(367.28143262390813,468.6019023677601,419.24534145090513,228.77707940248848,232.58988534346219,359.197687079234),c(217.12217287711695,194.2448166723502,244.75496190848932,59.71754646167898,86.74458838629123,114.41806220278472),c(681.8042064178626,674.9184308107084,745.9756494531468,466.8061730455188,405.12582488103044,376.07904051898913))
targetgene="Ehf"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(282.0559068216753,332.5207878628368,323.21707888393803,286.81244145679625,204.9460055280507,322.6214212930979),c(289.15803397186136,150.34768296108462,258.80787838170403,17.662936277398007,25.737405345383113,47.830501412639514),c(110.59026562432592,93.28140913643936,129.98947737723594,98.40778783121748,91.51077456136218,94.72314985640374),c(480.9154670268852,601.3907318443385,584.3671100111778,166.53625632975266,72.44602986107839,149.11862205117023))
targetgene="Meis2"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(356.12094710218713,383.00249163079224,377.08659203126103,179.15263938503693,130.59350119694395,150.99432798892082),c(459.609085576327,423.6073403137129,368.8890574218858,293.5411790862812,172.53593953756828,187.57059377505692),c(519.4698715564667,523.4733195068421,495.3653056808181,146.35004344129777,167.7697533624973,174.44065221080294),c(440.3318833115363,398.36648842973517,379.4287447767968,313.72739197473607,276.4387981541149,461.42366068664))
targetgene="Tnfrsf23"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(383.5148661100477,384.09991997357383,379.4287447767968,230.45926380985972,153.4711948372845,181.00562299292991),c(219.1513520628844,295.2082242082611,318.5327733928665,100.08997223858871,74.35250433110677,157.55929877104782),c(812.686263899863,852.7018223413339,950.9140146875278,361.66964758481635,257.3740534538311,419.2202770872522))
targetgene="Scgb2b4-ps"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(478.8862878411178,571.7601665892342,490.6810001897465,136.25693699707034,183.02154912272437,148.18076908229497),c(917.1889919668865,944.8858031349916,920.4660289955626,517.271705266656,548.1114101331589,394.83609989649483))
targetgene="Gm6829_Golga7"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(761.9567842556768,751.7384148054231,673.3689143415376,541.6633791735389,451.83444939672574,663.9999019637015),c(703.1105878684208,666.1390040684553,696.7904417968954,283.4480726420538,423.23733234630004,219.4575947168166),c(522.5136403351179,505.91446602233583,453.206556261174,182.5170081997794,256.4208162188169,169.7513873664265),c(406.85042674637333,459.82247562550697,491.8520765625144,6787.614083742948,7783.182023890856,6812.563965910067))
targetgene="Vmn1r233"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetgenelist=c("Kpna6","Prss55","Ddx17","4930426I24Rik_Zbtb1","Mrgprb13","Vmn2r55","Dnajc22","Ankrd13c","Tbc1d31","Rhox4a_Rhox4a2_Rhox4c_Rhox4f")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='10,11,12_vs_4,5,6 pos.'


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
targetmat=list(c(604.6953973586996,545.4218863624749,571.4852699107311,2630.095320924932,1953.183094544074,2533.140868932144),c(449.46318964748974,369.83335151741255,439.1536397879593,1111.9238932723888,1094.3163457962894,1326.1240979896525),c(766.0151426272117,758.322984862113,782.2790170089515,2948.028173918096,2754.8556091910073,3037.7057661870467),c(124.79451992469806,159.12710970333774,124.13409551339649,434.84466930546523,1492.7695100322205,692.13549102996))
targetgene="Kpna6"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(265.82247333553573,267.77251563872005,255.29464926340034,318.77394519684975,282.15822156420006,228.83612440556945),c(498.1634901059085,456.53019059716206,528.155444118319,669.5093941337531,770.215685891465,1001.626970758804),c(813.7008534927467,785.7586934316539,850.2014466294892,4813.5706816928005,6650.736188693999,6681.264550267528),c(148.1300805610237,122.91197439154364,195.56975425223789,2041.3307783449984,3078.9562690958314,1540.8924278620925))
targetgene="Prss55"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(289.15803397186136,280.94165575209973,254.12357289063243,85.79140477593319,136.31292460702906,109.7287973584083),c(304.37687786511725,311.66964934998566,282.22940583706185,796.5143168902816,703.4890794404718,822.4970537036246),c(661.5124145601881,671.6261457823634,567.9720407924274,958.8451122016062,954.1904722492035,663.9999019637015),c(527.5865882995365,489.45304088061124,553.9191243192126,19483.05980617369,22009.294519242616,18625.75996186315))
targetgene="Ddx17"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(212.04922491269832,248.01880546865056,176.8325322879516,2175.0644387310117,1948.416908369003,1822.2483185246779),c(363.22307425237324,310.572221007204,350.15183545759953,1976.5666786612057,2202.9312501177915,2321.1860979663293))
targetgene="4930426I24Rik_Zbtb1"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(605.7099869515832,617.852156986063,569.1431171651952,1971.520125439092,1170.5753245974245,1390.835952842047),c(447.4340104617223,423.6073403137129,544.5505133370696,1232.2000783994324,1298.309114089326,1304.5534797055209),c(397.71912041041986,403.8536301436434,535.1819023549264,963.0505732200343,1028.5429765803103,843.1298190188809),c(375.3981493669779,500.42732430842767,473.11485459822813,843.6154802966763,1076.2048383310198,1105.7286503039604))
targetgene="Mrgprb13"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(386.55863488869886,429.0944820276211,381.77089752233263,11259.701330739388,9416.077407470162,6854.767349509455))
targetgene="Vmn2r55"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(366.2668430310244,421.41248362814963,494.1942293080502,3911.0787471381304,3404.01016623567,3184.9486823004663),c(444.39024168307117,443.3610504837824,460.2330144977813,2890.833904067474,3755.754705955906,2988.937411805532),c(265.82247333553573,242.53166375474234,229.53096906250673,115.22963190492986,91.51077456136218,84.40676719877561),c(322.63949053702424,409.34077185755154,319.70384976563435,1016.0393820522282,993.2731988847853,523.3219566324088))
targetgene="Dnajc22"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(194.801201833675,254.6033755253404,262.32110750000766,322.13831401159223,373.6689961255622,333.8756569196013),c(361.19389506660576,362.1513531179411,316.1906206473307,296.9055479010237,302.176203499498,265.4123901917055),c(358.1501262879546,433.48419539874766,295.11124593750867,14239.691008397536,15365.230991193717,16814.765878964976),c(861.3865643582817,905.3783827948527,860.7411339844002,852.8674945372181,1557.5896420131853,1978.8697643268504))
targetgene="Ankrd13c"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(312.49359460818704,229.36252364136269,319.70384976563435,3098.583678377822,3915.8985614382896,2439.355572044615),c(316.55195297972193,255.70080386812202,297.45339868304444,232.9825404209166,299.31649179445543,396.71180583424535),c(1008.5020553264218,903.1835261092893,1068.021651964317,5643.7286867305065,4589.837286593322,4519.513457009996),c(304.37687786511725,325.93621780614694,249.4392673995609,253.16875330937145,207.80571723309328,205.38980018368733))
targetgene="Tbc1d31"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(280.02672763590783,288.6236541515712,330.2435371205454,2234.7819851926906,1549.963744133072,753.0959340068536),c(408.8796059321408,354.46935471846956,446.18009802456663,3453.5245883331536,3775.772687891204,3174.6322996428385))
targetgene="Rhox4a_Rhox4a2_Rhox4c_Rhox4f"
collabel=c("In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._hEGFRv3.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep1_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep2_1431..1432..1433..1434..1435..1436..1437..1438","In.vitro_Pool.E..33.40._mCD19.CAR_ET.1.2_rep3_1431..1432..1433..1434..1435..1436..1437..1438")

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
Sweave("invitro.poolE_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2_summary.Rnw");
library(tools);

texi2dvi("invitro.poolE_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2_summary.tex",pdf=TRUE);

