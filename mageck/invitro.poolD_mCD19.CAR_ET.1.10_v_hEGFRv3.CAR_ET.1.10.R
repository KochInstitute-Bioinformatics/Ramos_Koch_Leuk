pdf(file='invitro.poolD_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolD_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Mbd2","Dennd10","Chic2","Usp22","E2f3","Glmp","Zbed6","Smg9","Dym","Gpr108")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='8,9,10_vs_2,3,4 neg.'


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
targetmat=list(c(181.3526226407379,240.5628134071484,229.503388608258,111.36853380888675,160.12742800923178,160.59135377290906),c(378.4750385545835,392.5462962961608,379.94340889214004,189.32650747510746,255.17080463406614,197.48396207209086),c(166.8969788070559,194.8745271398993,155.9305319730748,87.06994461422055,69.21637210721632,36.16922382272726),c(420.5278206162038,416.78906043796644,549.0511689192775,179.20209531066322,233.476120839267,88.97629060390906))
targetgene="Mbd2"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(1169.593001088817,1192.5575129757474,1194.7353435683478,955.7445083235372,894.6474364902886,1077.119485440818),c(437.61176332873714,589.2856514469682,535.8739408652149,342.20513115821564,472.11764258205756,386.2873104267272),c(285.17042835536324,248.95453945623498,274.52558445963876,180.21453652710764,179.7559514426215,246.67410647099993),c(186.60922039844044,151.05106888355832,129.5760758649495,74.92065001688745,91.94413608272019,133.8261281440909))
targetgene="Dennd10"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(511.2041319365728,494.1794228906537,486.4593356624799,334.10560142666026,276.86548842886526,248.84425990036357),c(148.49888665509698,226.57660332533746,196.56031847310135,161.990594631108,106.40725861258628,213.39842055409085),c(187.9233698378661,188.34762910172086,184.48119275687725,94.15703312933152,102.2749378897674,120.80520756790906),c(115.64515066945606,131.47037476902298,122.98746183791816,82.00773853199843,47.521688312417176,47.74337544599999))
targetgene="Chic2"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(202.3790136715481,213.52280724898057,226.20908159474234,152.87862368310817,143.59814511795625,131.65597471472725),c(515.1465802548497,387.8842262688905,504.02897306789674,385.7401034653259,485.54768493121895,379.05346566218174),c(356.1344980843476,367.3711181489011,387.6301252570099,213.62509666977368,195.2521541531923,225.69595665381814),c(190.55166871671736,163.17245095446114,152.63622495955914,163.0030358475524,108.47341897399572,74.50860107481816))
targetgene="Usp22"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(254.9449912485736,232.17108735806184,213.0318535406797,143.76665273510835,100.20877752835796,130.20920576181814),c(214.20635862637883,147.32141286174206,211.93375120284114,144.77909395155277,193.18599379178286,150.46397110254543),c(304.8826699467478,248.02212545078092,258.05404939206045,189.32650747510746,190.0867532496687,149.74058662609087),c(576.9116039078547,561.3132312833463,498.53846137870397,579.1163758062111,428.7282749924593,391.351001761909))
targetgene="E2f3"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(349.56375088721944,314.22351983801946,308.56675693263395,227.7992736999956,269.63392716393224,313.9488627812727),c(381.1033374334347,376.6952582034417,315.1553709596653,219.69974396844023,232.44304065856227,220.6322653186363),c(482.2928442692088,491.3821808742915,479.87072163544855,281.45865817155016,265.5016064411133,290.07717505827264),c(366.6476935997527,343.1283540070954,383.2377159056557,378.6530149502149,390.5043083063846,404.3719223380908))
targetgene="Glmp"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(494.1201892240395,515.6249450160972,543.5606572300848,296.6452764182165,349.18110107819575,265.4821028588181),c(749.0651804726131,765.5118984777862,755.4944084329259,649.9872609573208,592.98802372451,604.7494223159998),c(785.8613647765309,706.7698161341802,790.6336832437596,751.2313826017634,777.9093760706551,819.5946118229998),c(483.60699370863443,445.6938946070424,432.6523211083907,341.1926899417712,303.72557312718806,255.35472018845448))
targetgene="Zbed6"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(76.22066748668695,102.565540599947,62.591833256797635,88.08238583066498,29.959325240436915,71.61506316899998),c(1030.2931605096994,970.6429796776803,1017.9408671763405,565.9546399924336,554.7640570384352,681.4281768201816),c(78.84896636553822,66.20139438723852,47.21840052705787,63.783796635998776,99.17569734765323,49.19014439890908),c(111.70270235117914,64.3365663763304,77.9652659865374,29.360795276888325,73.3486928300352,31.105532487545446))
targetgene="Smg9"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(528.2880746491061,413.0594044161502,461.2029818921931,370.55348521865955,389.4712281256799,369.64946746827263),c(491.4918903451882,422.3835444706908,395.3168416218798,318.91898317999386,352.2803416203099,354.4583934627272),c(465.20890155667547,548.2594352069895,556.7378852841474,468.7602832137688,429.761355173164,420.2863808200908),c(425.78441837390636,505.3683909561025,429.35801409487505,353.3419845391043,395.6697092099082,322.6294764987272))
targetgene="Dym"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(296.99777331019396,261.07592152713784,293.1933242028942,212.61265545332924,268.6008469832275,163.48489167872722),c(576.9116039078547,587.4208234360601,644.5860723112318,561.9048751266558,552.6978966770258,586.6648104046362),c(632.1058803637314,703.040160112364,661.0576073788102,453.5736649671024,483.4815245698095,437.6476082549999),c(605.8228915752187,593.0153074687845,616.0354115274293,429.2750757724362,325.4202569219872,457.90237359572717))
targetgene="Gpr108"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetgenelist=c("Ints14","Lsm12","Caprin1","Cwc27","Zfp608","Lyrm9","Gm7233","Pkib","Otub1","Prss50")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='8,9,10_vs_2,3,4 pos.'


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
targetmat=list(c(461.2664532383986,413.99181842160425,479.87072163544855,952.7071846742039,871.9196725147847,932.4425901499088),c(126.15834618486116,106.29519662176325,136.16468989198083,380.6778973831038,333.68489836762495,452.83868226054534),c(629.4775814848801,639.6360077414877,647.8803793247474,1651.2916240208572,1457.6761349743615,1670.294756133545),c(391.61653294883985,446.6263086124965,383.2377159056557,536.5938447155453,445.2575578837348,584.4946569752726))
targetgene="Ints14"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(612.3936387723469,670.4056699214717,616.0354115274293,995.2297157648698,1043.4109825117685,1080.7364078230908),c(253.63084180914794,229.37384534169965,259.152151729899,377.6405737337705,390.5043083063846,401.47838443227266),c(202.3790136715481,193.00969912899117,187.7754997703929,287.5333054702167,248.9723235498378,292.97071296409086),c(349.56375088721944,297.4400677398463,319.5477803110195,503.1832845728792,455.58835969078206,542.538357340909))
targetgene="Lsm12"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(211.57805974752756,248.02212545078092,234.99390029745078,250.07298046177297,257.2369649954756,275.60948552918177),c(436.2976138893115,486.7201108470212,441.43713981109914,635.8130839270989,651.873594024679,822.488149728818),c(270.7147845216812,311.42627782165727,299.7819382299255,456.6109886164357,347.1149407167863,486.837752653909),c(381.1033374334347,390.68146828525266,362.3737714867232,539.6311683648785,626.0465895070611,598.238962027909))
targetgene="Caprin1"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(333.7939576141118,289.04834169075974,334.9212130407593,388.7774271146592,479.34920384699063,475.9869855070908),c(164.26867992820462,170.63176299809365,210.83564886500255,293.60795276888325,209.7152766830584,251.73779780618176),c(488.86359146633697,524.9490850706378,423.8675024056823,556.8426690444337,610.5503867964902,755.2133934185453),c(851.5688367478127,873.6719231104576,687.4120634869355,963.8440380550926,1025.848619439788,1066.2687182939997))
targetgene="Cwc27"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(1876.6053994998097,1699.790731942758,1839.3214158795797,2387.336388375954,2637.453701339153,2221.5137271919084),c(388.98823406998855,405.6000923725177,408.49406967594246,438.387046720436,416.3313128240026,299.48117325218175),c(538.8012701645112,543.5973651797191,533.6777361895378,633.78820149421,760.3470129986748,820.3179962994543),c(190.55166871671736,199.53659716716962,214.12995587851825,319.9314243964383,208.6821965023537,347.94793317463626))
targetgene="Zfp608"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(406.0721767825218,370.16836016526327,344.8041340813063,502.1708433564348,395.6697092099082,468.0297562660908),c(377.16088911515783,340.33111199073323,327.2344966758894,567.9795224253224,542.3670948699786,609.8131136511817),c(813.4585030044693,892.3202032195388,763.1811247977957,828.1769150515396,1045.4771428731779,978.7391966429998),c(588.7389488626853,550.1242632178976,654.4689933517787,719.8457048919862,788.2401778777023,609.8131136511817))
targetgene="Lyrm9"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(758.2642265485925,722.6208542268993,737.924771027509,1015.4785400937583,761.3800931793795,1083.6299457289088),c(406.0721767825218,456.88286267249117,406.29786500026535,534.5689622826563,575.4256606525297,578.7075811636362),c(523.0314768914035,509.0980469779187,580.8961367165956,558.8675514773225,521.7054912558842,554.8358934406363),c(546.686166801065,587.4208234360601,540.2663502165691,880.8238583066498,891.5481959481743,821.0413807759089))
targetgene="Gm7233"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(639.9907770002852,605.1366895396873,602.8581834733667,986.11774481687,672.5351976387735,964.2715071139089),c(759.5783759880181,720.7560262159911,733.5323616761548,804.8907670733179,909.1105590201546,1001.8874998895452),c(722.7821916841003,609.7987595669576,742.3171803788632,1206.8299300017545,1366.7650790723462,1110.3951713577271),c(638.6766275608596,655.4870458342067,658.8614027031331,701.6217629959865,725.2222868547143,820.3179962994543))
targetgene="Pkib"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(562.4559600741726,517.4897730270053,522.6967128111522,725.9203521906527,628.1127498684705,694.4490973963635),c(525.6597757702548,553.8539192397138,640.1936629598775,823.1147089693175,673.5682778194782,749.4263176069089),c(533.5446724068087,630.311867686947,552.3454759327932,763.3806771990965,769.6447346250172,618.4937273686362),c(427.098567813332,382.2897422361661,408.49406967594246,507.2330494386569,463.85300113641983,572.1971208755454))
targetgene="Otub1"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(738.551984957208,767.3767264886943,708.276007905868,652.0121433902096,801.6702202268636,810.1906136290908),c(772.7198703822745,856.8884710122844,744.5133850545403,1168.3571637768664,1290.3171457001968,1459.7898734852724),c(733.2953871995054,693.7160200578234,685.2158588112584,697.5719981302088,647.7412733018601,548.3254331525453),c(428.41271725275766,356.1821500834523,418.37699071648944,705.6715278617642,665.3036363738404,804.4035378174543))
targetgene="Prss50"
collabel=c("In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
Sweave("invitro.poolD_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10_summary.Rnw");
library(tools);

texi2dvi("invitro.poolD_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10_summary.tex",pdf=TRUE);

