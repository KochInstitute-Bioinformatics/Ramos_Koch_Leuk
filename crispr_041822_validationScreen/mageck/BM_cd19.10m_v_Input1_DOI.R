pdf(file='BM_cd19.10m_v_Input1_DOI.pdf',width=4.5,height=4.5);
gstable=read.table('BM_cd19.10m_v_Input1_DOI.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Kdsr","H2-K1","Rps10","B2m","Rpl13a","Tmx2","Trmt112","H2-T23","Rfk","Elob")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2473_BM,AR2474_BM,AR2475_BM,AR2476_BM_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 neg.'


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
targetmat=list(c(539.8498910218691,555.3064076300514,592.6032706272189,6.860884766068447,402.36435923335625,19.490312867818293,75.4212359164916),c(898.8720136689333,877.2478525444063,871.1160266830677,252.13751515301544,391.48964682164393,41.76495614532491,621.4709839518908),c(625.4358493546044,712.4444938382485,613.0874346210039,264.1440634936352,36.24904137237444,80.7455818809615,541.5244738804097),c(1264.4776715723617,1283.5072949363305,1220.7842997699593,445.9575097944491,382.4273864785503,5382.110681927537,438.9515930339811),c(653.9645021321829,642.1794959402741,661.9619311675787,145.7938012789545,38.06149344099316,4797.401295892988,301.6849436659664),c(657.9146232860014,703.5016759239608,712.6332842048363,37.73486621337646,77.93543895060505,139.21652048441638,270.00802458103993))
targetgene="Kdsr"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(745.6950933708582,734.588614387913,721.6175666582508,879.9084712482784,317.1791120082763,2720.290810265496,184.02781563623952),c(631.5804822605444,623.442163167481,615.60303370796,8.57610595758556,5.437356205856166,8.352991229064983,78.43808535315127),c(563.1117155943562,536.1432263851493,608.0562364470918,10.29132714910267,402.36435923335625,768.4751930739784,122.1824021847164),c(665.8148655936386,650.6964653824529,608.774979043365,3.4304423830342237,25.37432896066211,55.686608193766546,9.050548309978993),c(443.29137392852664,425.84847210893514,403.5739678073783,6.860884766068447,67.06072653889271,44.54928655501324,3.016849436659664),c(478.84246431289364,467.5816223756108,484.073138589972,138.93291651288607,41.686397578230604,219.96210236537786,55.811714578203784))
targetgene="H2-K1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(337.9548098266985,352.17668643408933,310.13743029186776,108.05893506557804,54.37356205856166,590.2780468539254,58.82856401486345),c(506.49331238962355,538.6983172178029,521.0883822980397,1111.4633321030885,97.87241170541098,30.627634506571603,288.10912120099795),c(514.3935546972607,517.8317420844651,496.6511340247523,68.60884766068448,128.68409687192926,64.03959942283153,72.40438647983194),c(481.0369760650151,483.33801584364136,451.7297217576799,34.30442383034224,50.74865792132422,36.196295325948256,355.98823352584037),c(348.04956388645707,340.25292921503916,300.4344052421801,99.48282910799249,76.12298688198632,19.490312867818293,96.53918197310925),c(362.53334145045847,305.75920297421544,337.0902776521112,24.013096681239567,14.499616548949776,13.921652048441636,52.79486514154412))
targetgene="Rps10"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(646.0642598245458,666.0270103783745,678.4930108818613,63.46318408613314,67.06072653889271,44.54928655501324,78.43808535315127),c(718.0442452941284,708.1860091171591,756.476582577499,214.40264893963896,500.2367709387673,367.5316140788592,199.11206281953784),c(506.49331238962355,483.7638643157503,473.65137094401126,116.63504102316361,14.499616548949776,16.705982458129967,10.558973028308824),c(570.5730555515689,637.0693142749669,627.1029152483305,75.46973242675293,1038.5350353185277,27.843304096883273,159.8930201429622),c(892.288478412569,945.8094565539449,842.0069515340048,823.3061719282136,1946.5735216965074,378.66893571761256,766.2797569115547),c(555.2114732867191,527.6262569429706,531.869521242137,125.21114698074916,369.7402219982193,55.686608193766546,158.38459542463235))
targetgene="B2m"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(291.43116068172446,318.53465713748346,307.2624599067751,307.024593281563,115.99693239159821,11.13732163875331,7.54212359164916),c(376.57821666403555,363.24874670892166,336.01216375770144,0.0,14.499616548949776,66.82392983251987,0.0),c(283.5309183740874,267.8586889565202,285.3408107204438,41.16530859641068,0.0,108.58888597784477,88.99705838146009),c(406.4235764928868,373.4691100395361,383.44917511172986,48.026193362479134,155.8708779012101,36.196295325948256,377.106179582458),c(291.43116068172446,276.8015068708078,247.606824416103,12.006548340619783,23.561876892043387,2.7843304096883275,13.575822464968489),c(227.79031987020332,209.09159980548714,217.41963537263035,12.006548340619783,0.0,0.0,0.0))
targetgene="Rpl13a"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(442.8524715781024,389.6513519796756,424.4175030992999,63.46318408613314,10.874712411712332,5.568660819376655,699.9090693050421),c(359.46102499748844,333.4393536612962,373.0274074657691,475.11627005023996,28.999233097899552,2.7843304096883275,75.4212359164916),c(617.096704696543,608.5374666436683,582.9002455775312,46.31097217096202,21.749424823424665,217.17777195568954,16.59267190162815),c(389.74528717676407,352.17668643408933,425.49561699370963,116.63504102316361,18.12452068618722,698.8669328317702,57.32013929653362),c(671.9594984995786,775.0442192382619,681.7273525650905,56.60229932006469,690.5442381437331,467.767508827639,250.39850324275213),c(416.51833055264535,400.297563782399,390.2772297763249,60.032741703098914,195.74482341082197,55.686608193766546,9.050548309978993))
targetgene="Tmx2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(342.3438333309414,383.68947337015055,349.6682730868915,53.17185693703047,10.874712411712332,144.78518130379302,12.067397746638656),c(479.28136666331795,439.475623216421,447.41726618004094,65.17840527765026,309.9293037338015,11.13732163875331,36.20219323991597),c(397.64552948440115,384.96701878647735,401.41774001855885,205.8265429820534,63.43582240165527,8.352991229064983,25.643220211607144),c(596.907196577026,595.7620124804002,556.3067695154244,355.05078664404215,146.80861755811648,161.491163761923,149.33404711465337),c(238.7628786308104,227.40308410617135,245.09122532914694,193.81999464143362,1.812452068618722,36.196295325948256,36.20219323991597),c(375.700411963187,331.73595977286044,362.2462685216717,82.33061719282136,616.2337033303655,1378.243552795722,113.1318538747374))
targetgene="Trmt112"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(590.762563671086,644.7345867729277,636.087197701745,53.17185693703047,219.30670030286535,75.17692106158484,291.1259706376576),c(976.1188273436072,959.0107591893219,933.6466325588325,710.1015732880843,482.11225025258005,3480.4130121104095,360.51350768082983),c(802.3134965755909,720.5356148083182,796.3667966706593,214.40264893963896,250.11838546938364,353.6099620304176,28.66006964826681),c(785.6352072594681,749.4933109117258,738.50801767067,120.06548340619783,90.6226034309361,30.627634506571603,126.7076763397059),c(443.73027627895095,502.92704556065235,526.8383230682249,29.158760255790902,351.61570131203206,398.1592485854308,203.63733697452733),c(469.1866126035594,454.8061682123427,442.74543930426546,32.58920263882513,94.24750756817355,289.57036260758605,259.44905155273113))
targetgene="H2-T23"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(402.034552988644,410.9437755851224,388.83974458377855,246.9918515784641,105.12221997988587,247.80540646226115,46.76116626822479),c(506.9322147400478,529.3296508314063,563.8535667762926,22.297875489722454,150.43352169535393,1035.7709124040578,242.85637965110297),c(528.4384299108377,550.6220744368532,478.6825691179233,217.8330913226732,253.7432896066211,103.02022515846812,1131.318538747374),c(492.8873395264707,521.6643783334455,496.29176272661573,185.24388868384807,357.05305751788825,41.76495614532491,113.1318538747374),c(744.3783863195854,677.0990706532068,733.1174481986213,219.54831251419031,105.12221997988587,80.7455818809615,306.2102178209559),c(475.3312455094994,465.0265315429572,455.3234347390457,315.6006992391486,201.18217961667813,111.37321638753309,98.04760669143909))
targetgene="Rfk"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(354.6330991428213,327.4774750517711,333.4965646707454,5.145663574551335,38.06149344099316,25.058973687194946,76.92966063482143),c(565.7451296969018,578.3022251239339,541.9319175899612,27.44353906427379,7.249808274474888,27.843304096883273,25.643220211607144),c(765.0067967895267,758.4361288260135,847.3975210060535,677.5123706492592,302.6794954593266,765.6908626642901,343.9208357792017),c(694.343518371217,635.7917688586401,683.88358035391,533.4337905618218,77.93543895060505,1016.2805995362395,131.2329504946954),c(522.2937970048978,559.5648923511408,504.91667388189364,13.721769532136895,132.3090010091667,128.07919884566306,45.252741549894964),c(461.7252726463466,468.4333193198286,399.62088352787595,20.58265429820534,262.8055499497147,19.490312867818293,90.50548309978993))
targetgene="Elob"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetgenelist=c("Arl15","Cnot8","Ints10","Cd19","Tmem94","Setd5","Ewsr1","Zc3hav1","Tsc2","Thrap3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2473_BM,AR2474_BM,AR2475_BM,AR2476_BM_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 pos.'


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
targetmat=list(c(685.5654713627314,700.5207366191983,613.4468059191405,1653.4732286224958,3369.348395562204,3566.7272548107476,1363.6159453701682),c(941.4455416600888,855.9554289389596,874.7097396644335,1958.7826007125418,8212.22032291143,4098.534363061218,6362.535461915231),c(732.0891205077054,737.1437052205666,733.8361907948945,2171.470028460664,7367.617658935105,2795.467731327081,2802.653126656828),c(772.0292343963152,804.4277638137785,679.2117534781345,801.0082964384912,5529.791261355721,7286.592682154353,4715.335669499055),c(967.7796826855458,977.3222434900061,927.8966917886472,12124.898602834464,3697.402219982193,19326.037373646683,2695.5549716554096),c(766.3235038407996,821.0358542260269,767.2577215215963,3073.6763751986646,3382.035560042535,7345.063620757808,4992.885817671744))
targetgene="Arl15"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(733.8447299094025,812.5188847838482,784.1481725340155,1614.0231412176022,2883.6112411723866,4591.360845576052,2034.8649450269434),c(548.1890356799305,620.8870723348274,601.94692437877,744.4059971184265,1843.2637537852402,2906.840947714614,2170.623169676628),c(564.8673249960533,528.4779538871885,495.93239142847915,2718.625588554622,1946.5735216965074,14097.064864252003,1770.8906193192229),c(769.8347226441938,725.2199480015165,784.1481725340155,2543.673027019877,1834.2014934421466,6665.687000793856,2775.501481726891),c(571.4508602524176,548.4928320763084,538.697575906732,1005.1196182290275,2174.9424823424665,1470.1264563154368,3175.2340320842964),c(683.8098619610342,619.1836784463917,657.2901042918031,1653.4732286224958,2194.8794550972725,665.4549679155102,1837.2613069257354))
targetgene="Cnot8"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(679.4208384567913,610.240860532104,691.4303776147782,1159.4895254655676,6252.959636734591,6576.58842768383,4620.304912244275),c(651.7699903800615,683.9126462069498,711.19579901229,3368.6944201396077,11422.072936435186,462.19884800826236,8263.15060701082),c(495.0818512785922,486.7448036205128,568.1660223539316,751.266881884495,5176.36310797507,197.68745908787125,2012.238574251996),c(706.632784183097,837.6439446382753,739.5861315650798,3130.278674518729,688.7317860751143,16719.904110178406,2881.0912120099792),c(402.47345533906827,420.312441971519,376.6211204471349,2041.1132179053632,210.24443995977174,478.9048304663923,4111.965782167123),c(693.0268113199442,689.022827872257,665.1962728508079,3214.3245129030674,7614.111140267251,359.17862284979424,3722.7922048380256))
targetgene="Ints10"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(524.0494064065949,471.4142586245912,557.7442547079708,1008.5500606120618,781.1668415746692,8756.719138469789,1282.1610105803572),c(619.2912164486645,647.2896776055813,643.6339949626132,590.0360898818865,2903.5482139271926,7957.61631088924,413.308372822374),c(603.051829482966,690.3003732885838,608.774979043365,2294.9659542498957,1245.154571141062,48968.018915188615,1781.4495923475317),c(680.7375455080643,677.9507675974247,670.22747102472,3119.9873473696266,1455.3990111008338,12632.507068755942,2344.092012284559),c(653.9645021321829,587.6708915103304,610.9312068321844,852.4649321840046,7537.988153385265,16054.449142262896,657.6731771918068),c(456.0195420908309,422.8675328041726,477.6044552235136,1430.4944737252713,1116.4704742691326,31913.99515584761,724.0438647983194))
targetgene="Cd19"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(639.0418222177573,693.7071610654554,675.9774117949053,10078.63972135455,3838.773481334453,2628.407906745781,2043.9154933369225),c(378.77272841615695,413.9247148898849,451.7297217576799,396.21609524045283,868.1645408683678,735.0632281577184,333.3618627508929),c(664.0592561919415,627.2747994164614,757.9140677700453,3541.931760482836,1036.722583249909,426.0025526823141,3034.950533279622),c(617.5356070469674,568.9335587375373,635.0090838073353,1665.4797769631157,2789.363733604213,6582.157088503206,4734.945190837343),c(545.994523927809,462.0455922381946,518.2134119129471,2440.7597555288503,1071.1591725536648,2102.169459314687,547.5581727537291),c(859.809704481172,780.1544009035691,850.9912339874193,2349.853032378443,1788.8901917266785,7640.20264418477,1092.0994960707983))
targetgene="Tmem94"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(669.7649867474571,608.9633151157772,619.9154892855989,1384.1835015543093,636.1706760851714,2909.625278124302,1627.5902710778887),c(652.6477950809101,649.418919966126,600.150067888087,1329.2964234257618,1181.7187487394067,3302.2158658903563,9687.103541114182),c(484.10929251798507,557.4356499905961,550.5568287452392,1337.8725293833472,1317.6526538858109,913.2603743777714,2359.176259467857),c(705.7549794822484,680.5058584300783,725.9300222358897,8462.90135894543,11030.583289613542,272.8643801494561,13424.979993135505),c(396.7677247835526,382.41192795382375,428.72995867693885,1258.97235457356,5546.103329973289,245.02107605257282,1908.1572686872375),c(660.9869397389714,611.9442544205398,632.8528560185158,998.2587334629591,1424.5873259343155,868.7110878227581,683.3163974034139))
targetgene="Setd5"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(700.488151277157,785.6904310409853,752.8828695961332,2797.5257633644096,944.2875277503541,1381.0278832054105,1976.03638101208),c(568.3785437994476,549.3445290205264,490.1824506582939,300.1637085154946,1145.4697073670322,484.473491285769,1226.3492960021536),c(514.3935546972607,548.9186805484173,561.6973389874731,1591.7252657278798,826.4781432901373,2372.249509054455,2828.296346868435),c(687.3210807644285,655.8066470477601,648.6651931365252,2126.8742774812185,1681.955519678174,3605.707880546384,1333.4474510035716),c(545.5556215773848,563.3975286001212,582.1815029812582,1862.7302139875835,1792.515095863916,1219.5367194434875,570.1845435286765),c(477.5257572616208,489.2998944531665,524.3227239812688,1066.8675811236435,4822.93495459442,83.52991229064982,4890.312936825316))
targetgene="Ewsr1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(756.2287497810411,745.6606746627454,771.9295483973718,2339.5617052293405,1167.219132190457,13158.745516187035,728.5691389533089),c(613.5854858931488,558.7131954069229,622.7904596706916,445.9575097944491,2492.1215943507427,874.2797486421348,862.8189388846639),c(880.8770173015376,888.3199128192387,936.5216029439251,2389.3031197833366,3014.1077901129347,570.7877339861071,1401.326563328414),c(456.0195420908309,456.50956210077845,455.3234347390457,180.09822510929675,1261.4666397586304,4229.39789231657,1398.3097138917542),c(671.0816937987299,637.4951627470758,677.4148969874516,3646.56025316538,3606.779616551257,612.5526901314321,266.9911751443803),c(597.3460989274503,620.0353753906095,603.7437808694528,799.2930752469741,6869.193340064956,7161.297813718378,1316.8547791019434))
targetgene="Zc3hav1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(734.7225346102512,750.3450078559437,790.2574846023374,3509.342557844011,1381.0884762874662,2664.6042020717296,1117.7427162824056),c(720.2387570462498,738.4212506368935,743.8985871427187,1879.8824259027547,532.8609081739042,932.7506872455897,1849.328704672374),c(764.1289920886782,759.2878257702313,829.0695848010879,4277.7616516436765,7023.251765897548,130.86352925535138,2013.7469989703259),c(494.2040465777436,551.89961985318,479.0419404160599,552.30122366851,398.73945509611883,2380.60250028352,224.75528303114498),c(940.5677369592403,926.6462753090428,930.7716621737399,4264.0398821115405,975.0992129168724,3614.060871775449,2561.3051717240546),c(750.084116875101,822.3133996423537,769.0545780122792,8735.62152839665,1169.0315842590758,6451.293559247854,2452.698592004307))
targetgene="Tsc2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
targetmat=list(c(612.2687788418759,575.3212858191713,623.8685735651013,668.9362646916736,1758.0785065601603,2299.8569184025587,1345.5148487502101),c(649.57547862794,560.4165892953586,650.1026783290715,379.0638833252817,364.3028657923631,3753.2773922598653,4114.982631603782),c(572.7675673036904,645.5862837171456,670.9462136209931,557.4468872430614,1091.0961453084706,1999.149234156219,1077.0152488875),c(635.9695057647873,666.8787073225924,596.5563549067213,5739.130106816257,821.0407870842811,1308.6352925535139,963.8833950127627),c(825.5753211480779,755.029341049142,808.2260495091664,2077.1328629272225,3171.7911200827634,8545.110027333478,615.4372850785714),c(578.0343955087818,597.4654063688359,572.8378492297071,1282.9854512547997,1919.3867406672266,2633.976567565158,2071.0671382668593))
targetgene="Thrap3"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_BM","AR2474_BM","AR2475_BM","AR2476_BM")

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
Sweave("BM_cd19.10m_v_Input1_DOI_summary.Rnw");
library(tools);

texi2dvi("BM_cd19.10m_v_Input1_DOI_summary.tex",pdf=TRUE);

