pdf(file='invitro.poolB_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolB_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Kat6a","Rab1a","Kmt2a","Kdsr","Fdft1","Ube2m","Xrn1","Trib3","Aff1","Ube2f")
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
targetmat=list(c(205.32166330542046,293.68957485130227,262.039875546315,110.59627948167856,128.82736325822617,198.39535386578345),c(411.36375349963185,377.74555661908875,330.68694153450457,161.43489182406307,118.72325633601237,228.15465694565097),c(174.34330708740964,194.4427530050001,223.65656983248857,88.29864248940466,138.93147018044002,84.31802539295796),c(214.68721285970278,202.54453438020846,215.5370243930253,125.75867263642482,105.25111377306061,111.59738654950318))
targetgene="Kat6a"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(556.8899850354036,505.34861327862006,540.3188419715566,282.7340370620331,348.59168881637675,215.75494732903948),c(431.53570638577844,529.6539574042451,516.6983461476634,319.3021617293623,283.75700273217143,188.47558617249427),c(382.5466779479939,372.68194325958353,376.4516521932976,198.89492197108322,177.6638800489263,297.59303079867516))
targetgene="Rab1a"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(432.2561332745694,373.6946659314846,433.2884702695406,302.3559576152341,235.76249485165576,166.15610886259364),c(293.2137437379162,226.84987850583346,284.92223087571153,153.40774250684444,164.19173748597456,136.39680578272612),c(484.1268692675177,401.03817807281274,461.33780906041375,280.9502261026512,318.2793680497353,334.79215964850954),c(208.20337086058424,241.02799591244806,239.89566071141513,103.46103564415091,134.72142562951757,101.67761885621401))
targetgene="Kmt2a"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(342.9231990644917,300.77863355460954,360.21256131437104,172.13775758035453,201.2401295340919,203.35523771242802),c(388.3100930583215,371.6692205876825,343.2353299409478,239.03066855717626,215.55428100722813,181.0357604025274),c(234.85916574584937,300.77863355460954,270.89756148027493,140.92106579117106,150.71959492302278,143.836631552693),c(404.15948461172235,278.49873477278663,327.7343795565179,255.08496719161346,230.71044139054885,270.3136696421299))
targetgene="Kdsr"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(222.61190863640323,239.00255056864597,233.99053675544184,155.19155346622637,133.8794167193331,121.51715424279236),c(304.7405739585714,306.8549695860158,298.20875977665145,189.97586717417366,175.13785331837286,238.07442463894012),c(252.14941107683214,229.88804652153658,221.4421483489986,206.03016580861086,149.03557710265383,198.39535386578345),c(125.35427864962512,140.76845139424486,184.5351236241655,83.83911509094987,47.994507880515634,59.51860615973503))
targetgene="Fdft1"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(190.19269864081053,239.00255056864597,210.37004093154866,122.19105071766099,109.46115832398303,203.35523771242802),c(195.23568686234717,181.27735827028656,184.5351236241655,159.65108086468115,111.98718505453648,111.59738654950318),c(150.56921975730833,143.806619409948,133.60342950389582,50.8386123423845,51.36254352125358,39.67907077315669),c(222.61190863640323,196.4681983488022,196.34537153611208,137.35344387240724,118.72325633601237,166.15610886259364))
targetgene="Ube2m"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(447.38509793917933,406.10179143231795,434.02661076403723,346.059326120091,323.3314215108422,250.4741342555516),c(625.3305394705437,582.3155363430993,640.7059492231026,463.7908494392972,459.7368649607287,411.6703592715006),c(246.38599596650454,279.51145744468766,333.6395035124912,258.6525891103773,187.76798697114012,208.31512155907262),c(263.6762412974873,261.2824493504689,295.2561977986648,198.89492197108322,241.65655722294716,151.27645732265987))
targetgene="Xrn1"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(497.8149801545457,390.9109513538023,479.79132142283026,361.22171927483726,424.3724907329804,597.6660035206726),c(349.4070410636102,341.28754043065123,413.3586769181307,208.70588224768375,244.18258395350063,342.2319854184764),c(408.48204594446804,381.7964473066929,426.6452058190706,253.30115623223153,274.4949047201421,324.8723919552204),c(335.7189301765822,347.3638764620575,259.82545406282503,211.3815986867566,292.17709183401627,166.15610886259364))
targetgene="Trib3"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(825.6092145544275,713.9694836902348,724.1158251012254,620.7662138649055,566.6719965541583,818.3808346963567),c(319.1491117343904,393.9491193695054,397.11958603920414,268.4635493869778,261.8647710673748,334.79215964850954),c(218.28934730365754,227.86260117773452,219.96586736000526,119.51533427858813,125.45932761748824,148.79651539933758),c(198.83782130630192,176.21374491078134,209.631900437052,105.24484660353282,126.30133652767273,138.8767477060484))
targetgene="Aff1"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(348.6866141748193,326.0967003521356,438.4554537310172,194.43539457262844,302.2811987562301,255.43401810219618),c(120.31129042808848,204.56997972401052,190.4402475801388,110.59627948167856,122.93330088693479,173.5959346325605),c(252.14941107683214,282.5496254603908,305.59016472161807,182.84062333664602,149.8775860128383,163.67616693927133))
targetgene="Ube2f"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetgenelist=c("Irf8","Cnot8","Ints10","E130309D02Rik","Edc4","Hcfc2","Dgkd","Cnot9","Zfp36l2","Dcp1a")
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
targetmat=list(c(201.7195288614657,157.9847368165626,169.0341732397356,1165.7204619560796,1241.9631425221153,1192.852065118023),c(330.6759419550456,291.66412950750015,338.0683464794712,1555.4831565810275,1519.8260828829953,1822.7573136418853),c(197.39696752872,200.51908903640637,191.17838807463545,613.6309700273778,660.1349855846362,528.2276296676484),c(376.0628359488754,382.80916997859396,253.18218961235505,2263.6561074556466,1943.356564705791,2447.702678319103))
targetgene="Irf8"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(312.66526973527186,294.7022975232033,335.8539249959812,1146.0985414028787,1364.0544344988655,1450.7660251435414),c(283.12776729484295,331.1603137116408,321.8292556005446,819.6611358359887,986.8344427362163,622.4654227538955),c(226.93446996914892,220.77354247442722,175.67743769020555,812.525891998461,814.2226161483969,602.6258873673172),c(337.15978395416414,248.11705461575536,283.4459498867182,735.8220207450388,767.9121260882501,1031.6558401020739))
targetgene="Cnot8"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(259.3536799647416,284.5750708041929,299.68504076564477,650.199094694707,735.9157875012398,677.024145066986),c(348.6866141748193,429.3944128860419,342.49718944645116,684.0915029229634,670.2390925068499,634.865132370507),c(353.009175507565,326.0967003521356,380.14235466578094,704.6053289558554,786.4363221123089,1041.575607795363),c(1168.5324136189192,962.0865383059902,997.2278080649903,1898.8667662620455,1952.6186627178204,2348.5050013862115))
targetgene="Ints10"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(376.7832628376663,394.9618420414065,338.8064869739678,559.2247357662295,685.3952528901707,848.1401377762243),c(293.93417062670716,258.2442813347658,246.5389251618851,519.9808946598274,503.52132829032195,726.6229835334318),c(282.407340406052,296.72774286700536,341.02090845745784,514.6294617816817,658.4509677642671,597.6660035206726),c(613.8037092498885,756.5038359100786,680.5655359259223,1909.569632018337,1597.2909026199677,2006.273015967735))
targetgene="E130309D02Rik"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(937.9958092058156,1030.951679995261,893.8881388354577,2301.1161376026666,1961.8807607298497,2589.059367948474),c(323.47167306713607,303.8168015703127,248.0152061508784,621.6581193445965,660.1349855846362,543.1072812075822),c(261.5149606311145,230.90076919343764,302.6376027436314,672.4967316869809,462.2628916912822,495.9883846644586),c(564.814680812104,509.39950396622424,532.9374370265899,909.7435892847752,808.3285537771055,1006.8564208688509))
targetgene="Edc4"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(1200.9516236145118,1132.2239471853652,1063.6604525696898,1391.3725483178916,1506.3539403200434,1773.1584751754394),c(559.0512657017764,555.9847468736722,612.6566104322294,587.7657111163401,852.9550260168832,843.1802539295796),c(404.8799115005133,422.3053541827346,445.0987181814872,668.0372042885261,578.4601212967411,726.6229835334318),c(240.62258085617697,327.10942302403663,239.15752021691847,415.62795353598557,462.2628916912822,657.1846096804077))
targetgene="Hcfc2"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(597.9543176964877,754.4783905662765,752.1651638920986,1105.0708893370947,1071.0353337546649,1076.2947947218752),c(636.8573696911989,607.6336031406254,742.5693374636419,887.4459522925014,924.5257833825644,1051.4953754886521),c(447.38509793917933,464.8397064025784,411.14425543464074,617.1985919461416,737.5998053216088,669.5843192970191),c(773.0180516726882,860.814271115886,880.6016099345177,1277.2086469174492,1211.6508217554738,835.7404281596127))
targetgene="Dgkd"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(524.4707750398109,635.9898379538546,507.8406602137034,1052.4484660353282,1190.6005990008616,885.3392666260586),c(494.93327259938195,469.9033197620836,535.1518585100799,1153.2337852404062,1286.5896147618928,1490.4450959166982),c(172.1820264210368,194.4427530050001,197.8216525251054,492.3318247894078,479.1030698949719,565.4267585174828),c(95.09634932040527,88.10687245539067,84.14801637261947,167.67823018189975,209.66021863593673,171.1159927092382))
targetgene="Cnot9"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(375.34240906008444,362.5547165405731,405.23913147866745,548.521870009938,513.6254352125358,756.3822866132994),c(493.4924188218,478.00510113729194,550.6528088945098,842.8506783079536,849.5869903761452,748.9424608433325),c(256.4719724095778,247.1043319438543,346.18789191893444,500.35897410662636,468.15695406257356,543.1072812075822),c(968.2537385350354,899.2977326481255,859.9336760886112,1321.803920901997,1326.1640335405637,1775.6384170987617))
targetgene="Zfp36l2"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
targetmat=list(c(1275.1555931599796,1324.6412548465632,1324.9621876215083,2088.8426334362193,2000.613170598336,2321.225640229666),c(539.5997397044208,394.9618420414065,424.43078433558065,430.79034669073184,607.9304331531981,734.0628093033987),c(1110.1778356268524,1088.6768722936204,1135.260080535866,1714.2423319660177,1518.9840739728106,1368.9279416739057),c(485.5677230450996,582.3155363430993,459.12338757692373,1031.9346400024363,1010.4106922213819,840.7003120062574))
targetgene="Dcp1a"
collabel=c("In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x")

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
Sweave("invitro.poolB_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10_summary.Rnw");
library(tools);

texi2dvi("invitro.poolB_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10_summary.tex",pdf=TRUE);

