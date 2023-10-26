pdf(file='invitro.poolC_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolC_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Ecpas","Arid4b","Rnf183","Abi3bp","Ptbp1","Zfp869","Kctd5","Znrf1","Jade2","Zbtb25")
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
targetmat=list(c(460.99603489882475,559.5440715338403,429.5997578796704,356.48023555123723,239.40150858052732,315.1200129741543),c(936.8066871762506,994.2977386248818,845.620213067742,490.83752797784337,540.4399727284292,365.88402722227033),c(694.5441206320668,544.4484580931792,785.1305919869839,256.79579278311013,238.50821936940596,346.72779543052843),c(615.2423452524959,546.461206551934,620.9444776249259,348.895549688445,197.4169156578229,259.5669407781028))
targetgene="Ecpas"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(279.73483403123396,304.93139150135545,297.51017715230046,42.257535521271286,91.11549953437981,77.5827387565547),c(460.1245868177306,597.786292250182,450.5859529485049,546.0973821210442,535.080237461701,242.32633216553506),c(845.3046386613611,847.3671011357799,832.0409103761432,434.494147282815,662.820594652057,537.3323017583604),c(555.9838757380911,577.6588076626338,722.1720067804804,255.71226623128265,173.2981069575459,363.01059245350905))
targetgene="Arid4b"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(372.1083306272177,416.638930962248,434.53768613116085,331.55912485920544,201.88336171342976,178.1529556631997),c(331.1502708157909,322.0397534007714,307.3860336552814,125.68908001198638,142.03298456829793,99.61240531705789),c(410.45204619536196,483.05963010115715,288.8688027121922,450.74704556022704,284.9592583477172,257.6513175989286),c(570.7984931166923,548.4739550106889,375.2825471132753,218.8723634691487,309.07806704799424,330.4449984075478))
targetgene="Rnf183"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(370.36543446502935,460.91939705485413,411.08252693658113,260.0463724385925,401.08685579349543,210.71854970916092),c(277.99193786904556,260.6509254087494,211.09643275121735,123.52202690833144,135.77996009044833,78.5405503461418),c(598.6848317117062,534.3847157994051,561.6893386070403,470.2505234931215,294.7854396700523,329.4871868179607),c(240.51967038199552,218.38320777489812,278.99294620921125,106.18560207909195,106.3014161234431,132.17799936301913))
targetgene="Abi3bp"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(1677.5375561063095,1549.8163132412126,1638.1576974319617,1136.6193528670149,1199.6874105360007,1307.4128197863847),c(393.02308457347823,313.98875956575216,324.66878253549805,128.9396596674688,150.96587667951164,234.6638394488383),c(337.25040738345024,509.22536006496983,528.3583229094797,357.56376210306473,258.1605820140761,158.9967238714578),c(275.377593625763,214.3577108573885,261.7101973289946,627.3618735081044,696.7655846746691,830.4226481720115))
targetgene="Ptbp1"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(399.9946692222317,549.4803292400662,524.6548767208618,325.05796554824065,167.0450824796963,237.5372742175996),c(372.9797787083119,452.8684032198348,422.1928655024347,375.98371348413167,268.88005254753256,580.4338232897796),c(701.5157052808203,775.9145308499836,812.2891973701813,1121.4499811414303,1258.6444984700113,972.1787634309015),c(502.82554279134575,685.3408502060167,555.5169282926772,201.5359386399092,275.13307702538214,214.5497960675093))
targetgene="Zfp869"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(734.6307323623994,879.571076475857,756.7375045409136,802.8931749041544,1642.7588592522006,1797.8123536549774),c(862.7336002832449,728.6149420692453,827.1029821246527,937.2504673307606,1894.6664167884271,847.6632567845792),c(954.2356487981344,905.7368064396696,1211.0269036780364,883.0741397393871,2034.0195337233608,1688.6218324420488),c(1284.5144715328313,1211.6745721704026,1431.9991929322346,737.8815817945062,485.9493308500256,412.816795112038))
targetgene="Kctd5"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(515.8972640077585,418.6516794210029,522.1859125951166,265.4640051977299,192.95046960221606,110.14833280251594),c(467.96761954757824,537.4038384875373,528.3583229094797,352.1461293439274,560.0923353730993,330.4449984075478),c(333.7646150590735,306.94413996011025,328.37222872411587,85.59859759437003,114.34101902353544,168.57483976732874),c(370.36543446502935,424.68992479726734,457.9928453257406,442.0788331456073,870.9569808433364,560.3197799084506))
targetgene="Znrf1"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(269.27745705810366,242.53618927995598,229.61366369430658,173.36424829239502,50.917485033918126,167.61702817774164),c(223.962156841206,371.3520906402646,304.91706952953615,127.85613311564133,100.04839164559351,185.81544837989645),c(279.73483403123396,170.07724476478242,262.9446793918672,115.93734104553917,68.78326925634553,156.1232891026965),c(289.32076292326997,235.49156967431412,143.19991929322347,164.69603587777527,84.8624750565302,75.66711557738051))
targetgene="Jade2"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(291.06365908545837,267.69554501439126,303.68258746666356,487.58694832236097,521.6808992948804,227.00134673214154),c(473.1963080341434,588.7289241857852,312.32396190677184,349.97907624027243,219.74914593585717,250.94663647181892),c(462.73893106101315,586.7161757270304,576.5031233615117,370.5660807249943,298.3585965145378,336.1918679450704),c(449.66720984460034,447.83653207294776,438.2411323197787,217.78883691732122,399.30027737125266,254.7778828301673))
targetgene="Zbtb25"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetgenelist=c("Ddx6","Ints13","Irf2","Rlim","Gatad2b","Zfp341","Hexim1","Rcor2","Fam234b","Grk3")
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
targetmat=list(c(849.6618790668321,961.0873890554273,1053.0131996303414,10830.931412067379,12958.053296526603,11411.367278340651),c(874.9338734185634,1145.2538730314934,918.4546547772263,10720.411703780977,7868.091371557032,11651.777987327012),c(369.4939863839352,283.7975326844298,460.4618094514858,3784.7582455333486,6909.592048023802,5427.918278190068))
targetgene="Ddx6"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(756.416934389754,833.277861924496,872.7788184509395,4407.786012834144,4012.6551363571966,4144.45074814336),c(903.6916600946716,852.3989722826669,1006.1028812411821,3108.637677193008,2913.0161174667896,2853.3207253799565),c(423.5237674117747,488.0915012480442,475.2755942059572,1316.4847604703746,1448.9151004388632,1778.6561218632355),c(849.6618790668321,719.5575740048487,801.1788588043279,2352.336144017435,1642.7588592522006,3202.921955579246))
targetgene="Ints13"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(276.24904170685716,286.8166553725621,253.0688228888863,847.3177635290806,668.1803299187852,830.4226481720115),c(310.23551686953044,415.6325567328706,371.57910092465744,1122.5335076932577,1167.5289989356313,1013.3646617831466),c(172.54672005664898,321.033379171394,251.8343408260137,730.296895931714,933.4872256218323,1342.8518486011073),c(247.49125503074902,285.8102811431846,308.620515718154,1463.8443715189103,1678.4904276970553,920.4569375931984))
targetgene="Irf2"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(1606.078813456586,1486.4147367904357,1433.2336749951073,4646.161854236187,6107.418336436811,5125.249815880546),c(1277.5428868840777,1090.9096646451133,1111.0338565853544,3045.793137187015,3030.0370041236893,3526.662272859684),c(143.7889333805408,150.9561344066116,96.28960090406406,119.18792070102157,427.8855321271365,216.4654192466835),c(1194.75531918013,1406.91117266962,1282.6268633246482,6284.454000599319,5660.773730876125,7298.524312653664))
targetgene="Rlim"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(583.8702143331051,628.9838933608817,655.5099753853591,988.1762152666515,1343.5069735265413,1823.673266573829),c(599.5562797928004,757.7997947211903,567.8617489214034,2887.598260620204,2484.237296128532,2769.033305496292),c(386.9229480058189,523.3145992762536,460.4618094514858,1565.6958673906925,1991.141651589535,2024.813700387119),c(442.69562519584684,443.8110351554381,465.39973770297627,1385.8304597873325,1445.3419435943777,1528.6672969810038))
targetgene="Gatad2b"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(352.06502476205145,387.45407831030315,334.544639038479,2473.691117822111,2723.638804709059,2073.662091456061),c(579.5129739276341,565.5823169101049,738.2202735978244,1490.932535314597,1716.008574564153,1827.5045129321775),c(684.9581917400308,828.245990777609,830.8064283132707,1424.8374156531215,1383.704988027003,1810.2639043196098),c(587.3560066574819,657.1623717834492,695.0134013972829,1541.858283250488,1305.988826659444,915.6678796452629))
targetgene="Zfp341"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(721.5590111459866,810.1312546488157,848.0891771934872,1435.672681171396,1777.6455301315275,2042.0543089996868),c(610.0136567659307,627.9775191315043,661.6823856997222,1034.7678569952327,1037.1087741119113,1025.816212447779),c(359.9080574918991,326.0652503182811,407.3790807479633,821.3131262852214,461.8305221497486,703.0337067569278),c(711.9730822539506,697.4173409585457,827.1029821246527,1804.0717087927355,2175.1592290805374,1654.1406152169134))
targetgene="Hexim1"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(521.9974005754178,485.07237855991195,444.4135426341418,1100.8629766567083,784.3079273645634,1203.0113565213915),c(480.1678926828969,561.5568199925951,427.13079375392516,913.4128831905563,1028.1758820006976,815.0976627386179),c(511.5400236022876,601.8117891676916,561.6893386070403,982.7585825075142,1083.5598130902226,990.3771836330563),c(428.75245589833986,376.38396178715163,449.35147088563224,760.6356393828831,827.1858094983892,618.7462868732634))
targetgene="Rcor2"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(498.4683023858748,524.320973505631,481.4480045203203,2315.496241255301,1460.527860183441,2784.3582909296856),c(819.1611962285355,868.5009599527054,807.3512691186909,2196.3083205542794,2719.172358653452,3361.918679450704),c(720.6875630648924,717.5448255460939,804.8823049929457,1466.0114246225653,1168.4222881467529,1487.4813986287586),c(270.1489051391979,226.4342016099174,207.3929865625995,235.12526174656074,282.2793907143531,171.44827453609003))
targetgene="Fam234b"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(500.2111985480632,437.77278977917365,540.7031435382058,1112.7817687268105,1155.0229499799323,959.7272127662693),c(269.27745705810366,301.9122688132232,292.57224890081,796.3920155931896,788.7743734201703,1159.9098349899723),c(461.867482979919,508.21898583559243,409.84804487370855,406.3224569353008,457.36407609414175,365.88402722227033),c(328.5359265725084,421.6708021091351,391.3308139306193,1073.7748128610217,1213.0867487028213,933.8662998474177))
targetgene="Grk3"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
Sweave("invitro.poolC_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2_summary.Rnw");
library(tools);

texi2dvi("invitro.poolC_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2_summary.tex",pdf=TRUE);

