pdf(file='invivo.poolF_BM_mCD19_10m_v_BM_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolF_BM_mCD19_10m_v_BM_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Gins2","Dok4","Kmt5b","Ndrg3","2310022A10Rik","Sprtn","Ddx23","Polr1c","Gm1673","Brd9")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='5,6,7,8,9_vs_1,2,3,4 neg.'


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
targetmat=list(c(25.042175354981307,0.0,24.803767623941333,12.70236252750086,0.0,0.0,0.0,0.0,0.0),c(10.173383737961156,24.55101923328975,4.960753524788267,62.534707827696536,0.7438254159563782,0.0,23.321009653660404,0.0,11.649973676268349),c(28.17244727435397,22.96708250856138,7.4411302871824,1.9542096196155168,105.6232090658057,0.0,0.0,0.0,365.0325085230749),c(3.9128398992158293,0.7919683623641854,0.8267922541313778,19.542096196155168,0.0,0.0,11.046794046470719,0.0,148.53716437242144))
targetgene="Gins2"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(76.69166202463025,5.543778536549298,41.33961270656889,23.4505154353862,113.80528864132586,20.00522161806558,20.86616653222247,65.80572970204317,28.154103050981842),c(57.127462528551106,91.07636167188133,38.859235944174756,6.839733668654309,63.96898577224852,123.06242389294886,115.37762670758306,2.437249248223821,137.85802183584212),c(101.73383737961156,4.751810174185112,119.88487684904977,4.885524049038792,0.0,0.0,3.682264682156906,177.91919512033894,0.970831139689029),c(13.303655657333818,13.463462160191153,120.71166910318117,111.38994831808445,2.9753016638255128,0.0,0.0,0.0,11.649973676268349))
targetgene="Dok4"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(86.08247778274824,69.69321588804831,367.92255308846313,39.084192392310335,188.93165565292006,26.673628824087437,20.86616653222247,82.86647443960992,64.07485521947592),c(3.1302719193726634,43.5582599300302,28.110936640466846,29.31314429423275,45.37335037333907,25.461191150265282,23.321009653660404,0.0,0.970831139689029),c(237.9006658723224,253.42987595653935,299.29879599555875,178.8101801948198,194.1384335646147,96.99501390577251,368.2264682156906,36.558738723357315,249.50360290008047),c(222.2493062754591,302.53191442311885,42.993197214831646,313.65064394829045,0.0,41.829099746864394,4.909686242875875,0.0,303.8701467226661))
targetgene="Kmt5b"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(83.73477384321875,62.56550062677065,48.78074299375129,38.10708758250258,63.96898577224852,27.886066497909596,116.60504826830203,26.809741730462033,83.4914780132565),c(22.694471415451808,13.463462160191153,36.378859181780626,259.90987940886373,0.0,32.12959835628714,0.0,0.0,147.5663332327324),c(36.78069505262879,59.39762717731391,89.2935634461888,22.473410625578442,0.0,45.46641276833086,0.0,0.0,3.883324558756116),c(39.12839899215829,68.10927916331995,89.2935634461888,32.244458723656024,57.27455702864112,10.305720227488328,23.321009653660404,14.623495489342925,95.14145168952484))
targetgene="Ndrg3"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(209.72821859796844,200.36799567813893,241.42333820636233,202.260695630206,53.55542994885923,116.394016686927,141.1534794826814,199.85443835435333,191.25373451873872),c(54.77975858902161,38.0144813934809,177.76033463824623,42.01550682173361,21.570937062734966,89.72038786283956,273.71500804033,31.684240226909672,15.533298235024464),c(10.955951717804322,33.26267121929579,0.8267922541313778,43.96971644134913,0.0,5.4559695321997035,0.0,0.0,0.0),c(9.39081575811799,60.189595539678095,91.77394020858294,72.30575592577412,32.72831830208064,0.0,0.0,0.0,8.73748025720126))
targetgene="2310022A10Rik"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(4.695407879058995,0.7919683623641854,15.709052828496178,0.9771048098077584,8.18207957552016,0.0,9.81937248575175,0.0,3.883324558756116),c(52.43205464949211,129.8828114277264,40.51282045243751,126.04652046520083,0.0,0.0,0.0,0.0,216.49534415065347),c(181.55577132361447,129.09084306536224,353.8670847682297,96.73337617096809,211.2464181316114,491.64347673488436,136.2437932398055,243.7249248223821,761.1316135161987),c(15.651359596863317,7.127715261277669,4.960753524788267,13.679467337308617,7.438254159563781,0.0,0.0,0.0,1.941662279378058))
targetgene="Sprtn"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(9.39081575811799,25.342987595653934,14.8822605743648,59.60339339827326,0.0,0.0,3.682264682156906,0.0,1.941662279378058),c(2.3477039395294974,1.5839367247283709,0.0,11.7252577176931,0.0,4.243531858377547,0.0,58.4939819573717,1.941662279378058),c(104.86410929898422,133.05068487718316,48.78074299375129,42.992611631541365,426.2119633430047,118.21267319766024,8.59195092503278,43.87048646802878,159.21630690900076),c(0.7825679798431658,26.926924320382305,125.67242262796942,130.93204451423964,151.74038485510115,0.0,0.0,34.121489475133494,0.970831139689029))
targetgene="Ddx23"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(11.738519697647488,0.0,0.0,24.42762024519396,0.0,0.0,0.0,0.0,5.8249868381341745),c(0.0,3.1678734494567418,15.709052828496178,4.885524049038792,0.0,0.0,0.0,0.0,4.854155698445145),c(3.9128398992158293,11.879525435462782,4.133961270656889,49.832345300195676,0.0,0.0,27.00327433581731,0.0,15.533298235024464),c(4.695407879058995,0.0,2.4803767623941333,0.9771048098077584,0.0,0.0,0.0,12.186246241119106,0.0))
targetgene="Polr1c"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(81.38706990368925,85.53258313533203,58.70225004332782,3.9084192392310335,5.9506033276510255,10.911939064399407,0.0,0.0,2.9124934190670873),c(32.86785515341297,9.503620348370225,3.307169016525511,170.99334171635772,133.88857487214807,43.64775625759763,189.02292035072117,0.0,47.57072584476242),c(48.519214750276284,60.98156390204228,14.8822605743648,83.05390883365946,37.93509621377529,0.0,0.0,0.0,7.766649117512232),c(0.7825679798431658,0.0,0.8267922541313778,0.9771048098077584,9.669730407432917,0.0,0.0,0.0,0.0))
targetgene="Gm1673"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(35.99812707278563,98.99604529552317,49.607535247882666,192.4896475321284,56.53073161268474,13.336814412043719,17.18390185006556,43.87048646802878,15.533298235024464),c(165.90441172675116,76.0289627869618,62.009419059853336,94.77916655135256,0.0,0.0,13.501637167908655,0.0,98.05394510859193),c(52.43205464949211,88.70045658478877,32.244897911123736,697.6528342027394,43.88569954142631,30.917160682464985,39.277489943007,99.92721917717667,0.0),c(0.7825679798431658,0.0,0.0,0.9771048098077584,0.0,0.0,39.277489943007,0.0,6.795817977823203))
targetgene="Brd9"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetgenelist=c("Csn3","Nrd1","Lhfpl5","Rpap3","Shbg","Ndufb4_Ndufb4b_Ndufb4c","Klf12","Cisd2","L3mbtl3","Elk3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='5,6,7,8,9_vs_1,2,3,4 pos.'


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
targetmat=list(c(40.693534951844626,28.510861045110676,99.21507049576533,211.0546389184758,310.9190238697661,358.88155145135823,87.14693081104677,31.684240226909672,304.8409778623551),c(18.78163151623598,116.41934926753525,161.22448955561867,77.19127997481291,57.27455702864112,31.523379519376064,87.14693081104677,21.93524323401439,110.6747499245493),c(107.99438121835689,86.32455149769622,205.87127127871307,38.10708758250258,688.7823351756062,43.64775625759763,309.3102333011801,182.79369361678658,729.0941859064608),c(91.5604536416504,82.36470968587528,72.75771836356125,170.99334171635772,33.47214371803702,278.86066497909593,371.9087328978475,38.99598797158114,299.98682216390995))
targetgene="Csn3"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(0.0,8.71165198600604,0.0,0.0,0.0,27.279847660998517,15.956480289346594,87.74097293605756,168.92461830589104),c(49.30178273011945,38.0144813934809,0.0,16.610781766731893,168.84836942209785,92.75148204739496,28.23069589653628,7.311747744671463,52.42488154320757),c(72.77882212541442,42.766291567666016,18.189429590890313,88.91653769250601,22.314762478691346,0.0,7.364529364313812,0.0,8.73748025720126),c(27.389879294510806,102.95588710734411,81.02564090487502,152.42835033001032,152.48421027105752,129.73083109897073,49.09686242875875,255.9111710635012,522.3071531526977))
targetgene="Nrd1"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(118.9503329361612,44.350228292394384,17.362637336758933,24.42762024519396,199.34521147630934,277.0420084683627,119.05989138973996,90.17822218428138,60.1915306607198),c(72.77882212541442,108.49966564389341,92.60073246271432,77.19127997481291,317.6134526133735,11.518157901310484,158.33738133274696,155.98395188632455,13.591635955646407),c(57.127462528551106,264.5174330296379,62.83621131398471,50.80945011000344,167.3607185901851,34.55447370393146,27.00327433581731,146.23495489342926,8.73748025720126),c(1.5651359596863317,11.087557073098596,0.8267922541313778,0.9771048098077584,80.33314492328884,56.37835183273027,2.4548431214379374,24.37249248223821,59.22069952103077))
targetgene="Lhfpl5"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(111.90722111757272,22.175114146197192,75.23809512595538,37.12998277269482,1119.4572510143491,186.10918293170099,381.72810538359926,31.684240226909672,5.8249868381341745),c(16.433927576706484,205.11980585232402,62.009419059853336,53.74076453942671,165.12924234231596,26.06740998717636,93.28403861464162,497.19884663765947,216.49534415065347),c(75.12652606494392,129.09084306536224,28.110936640466846,85.98522326308273,18.595635398909454,9.69950139057725,30.685539017974218,116.9879639147434,179.60376084247036),c(68.08341424635543,25.342987595653934,56.22187328093369,22.473410625578442,158.43481359870856,50.31616346361949,228.30041029372816,46.3077357162526,155.33298235024463))
targetgene="Rpap3"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(2.3477039395294974,45.14219665475857,2.4803767623941333,3.9084192392310335,24.546238726560482,99.41988925341681,204.97940064006778,804.2922519138609,16.504129374713493),c(177.64293142439865,28.510861045110676,42.993197214831646,28.336039484424994,98.92878032219829,184.89674525787882,256.53110619026444,2.437249248223821,24.270778492225727),c(70.43111818588493,68.10927916331995,105.00261627468498,53.74076453942671,52.067779116946475,10.305720227488328,60.14365647522946,2.437249248223821,12.620804815957378),c(174.512659505026,167.89729282120732,251.34484525593885,99.66469060039135,20.82711164677859,148.52361504321414,106.78567578255027,41.43323721980496,300.957653303599))
targetgene="Shbg"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(46.17151081074678,194.82421714158963,90.12035570032018,139.72598780250945,242.48708560177928,149.7360527170363,413.64106596229243,299.78165753153,88.34563371170164),c(27.389879294510806,1.5839367247283709,8.267922541313778,0.0,275.2154039038599,119.42511087148239,13.501637167908655,14.623495489342925,346.58671686898333),c(4.695407879058995,32.4707028569316,30.59131340286098,42.01550682173361,22.314762478691346,6.062188369110782,36.82264682156906,12.186246241119106,0.0))
targetgene="Ndufb4_Ndufb4b_Ndufb4c"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(88.43018172227774,122.75509616644874,29.7645211487296,72.30575592577412,159.17863901466492,1188.7951391826243,18.41132341078453,492.32434814121183,366.00333966276395),c(35.21555909294246,191.65634369213288,117.40450008665564,161.22229361828013,72.89489076372506,55.77213299581919,38.05006838228803,34.121489475133494,85.43314029263455),c(75.12652606494392,41.18235484293764,46.30036623135716,140.70309261231722,261.8265464166451,410.4101525887999,336.31350763699743,46.3077357162526,152.42048893117754),c(107.21181323851371,108.49966564389341,100.86865500402808,212.03174372828357,98.18495490624193,85.47685600446202,200.0697143971919,80.4292251913861,44.65823242569533))
targetgene="Klf12"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(323.20057567522747,182.15272334376266,137.2475141858087,350.78062672098525,167.3607185901851,298.25966776025047,79.78240144673296,68.24297895026699,111.64558106423834),c(68.08341424635543,7.127715261277669,9.921507049576533,13.679467337308617,111.57381239345673,137.00545714190366,61.371078035948436,321.7169007655444,64.07485521947592),c(76.69166202463025,58.60565881494972,15.709052828496178,45.92392606096465,81.07697033924522,121.8499862191267,268.80532179745416,165.73294887921983,145.62467095335435),c(156.51359596863315,187.69650188031196,151.30298250604213,182.71859943405082,127.93797154449705,468.00094209535234,123.96957763261584,368.02463648179696,432.99068830130693))
targetgene="Cisd2"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(15.651359596863317,33.26267121929579,100.04186274989671,24.42762024519396,5.2067779116946475,56.98457066964134,42.9597546251639,12.186246241119106,45.629063565384364),c(41.47610293168779,68.10927916331995,133.9403451692832,77.19127997481291,60.99368410842301,441.93353210817594,355.95225260850094,297.34440828330617,141.74134639459822),c(18.78163151623598,19.799209059104637,7.4411302871824,51.786554919811195,8.92590499147654,0.0,12.274215607189687,0.0,1.941662279378058),c(99.38613344008206,118.79525435462782,156.26373603083042,120.18389160635428,73.63871617968144,1227.5931447449332,93.28403861464162,485.0126003965404,609.6819557247102))
targetgene="L3mbtl3"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(29.7375832340403,8.71165198600604,19.843014099153066,2.931314429423275,14.132682903171185,52.7410388112638,82.23724456817091,29.24699097868585,176.69126742340327),c(11.738519697647488,13.463462160191153,0.0,8.793943288269826,61.73750952437939,223.08853198327677,137.47121480052448,90.17822218428138,10.679142536579318),c(171.3823875856533,81.5727413235111,33.07169016525511,121.16099641616204,112.3176378094131,75.77735461388477,19.6387449715035,338.7776455031111,29.12493419067087),c(189.38145112204614,102.95588710734411,105.82940852881636,523.7281780569585,356.29237424310514,400.1044323613116,497.1057320911823,677.5552910062222,52.42488154320757))
targetgene="Elk3"
collabel=c("In.vivo_AR1911_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._BM_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._BM_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
Sweave("invivo.poolF_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolF_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.tex",pdf=TRUE);

