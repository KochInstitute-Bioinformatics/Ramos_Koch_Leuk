pdf(file='BM_cd19.15m_v_egfr.15m.pdf',width=4.5,height=4.5);
gstable=read.table('BM_cd19.15m_v_egfr.15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Pak2","Fitm2","Ptpn2","Zbtb7a","Kat6a","Tecr","Btbd35f27","Sys1","Elob","Map3k7")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2480_BM,AR2481_BM,AR2482_BM,AR2483_BM_vs_AR2466_BM,AR2467_BM,AR2468_BM,AR2469_BM neg.'


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
targetmat=list(c(3361.956418564977,2751.711491694909,3410.226885924859,3417.482874779378,932.1007978821056,1672.0376364741076,1397.1426567388507,2337.1557990034007),c(501.83981791853967,754.0956437347448,978.6390535312535,1577.4945648215996,463.0969743596773,504.9349131554071,253.0413409315607,344.17345031977055),c(1547.1851212949812,1175.8338288384111,769.502486876824,1585.9348781700298,178.38684471507977,256.9415634031312,212.67287798138955,651.3695877952682),c(2825.858162891851,1975.0537228311614,336.21955213343676,849.9395541869186,567.0575196240947,639.1581179182368,2183.835385938528,392.5283978853582),c(1312.9314793514561,1214.015886749031,444.29012246204144,1359.7344804321012,385.12656541136425,67.75076049933311,81.72152255766358,541.3857854892259),c(3958.238416239404,2190.2616856001105,1210.7906490519597,2640.130015388959,67.33808045536124,200.69564902632638,521.8362283802614,1049.5868030412837))
targetgene="Pak2"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(772.2037050628452,518.9288779670626,159.1038952060013,587.4458090507402,36.62246480905611,414.1744604110175,484.42155540205397,283.492731806092),c(624.9850130904323,613.5162487001894,236.15420923658058,813.6462067886689,62.612601125160445,25.566324716729476,18.70733648910371,71.11021700821706),c(699.9832146612841,702.029201129354,467.30515132831835,1050.819011879557,120.4997229201201,8.948213650855315,597.6501709939974,563.1929187050791),c(353.6952222230298,565.7886763119144,327.2136712727197,1215.4051221739453,216.19067935668608,24.288008480893,41.35305960749241,108.08752985248992),c(3991.570950270894,207.39799637814042,1907.2454356140788,2442.6266830356926,274.07780115164576,31.957905895911843,368.2391498381467,1409.8785692162503),c(240.73496800520354,653.4338546976558,79.05162088851638,61.61428744354028,250.45040450064178,2.5566324716729474,30.522496376958685,35.08104039072042))
targetgene="Fitm2"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(111.1084467716324,208.2657704215636,464.3031910414127,140.10920158394092,987.6251800119647,203.25228149799932,95.50587576016105,141.2722977896579),c(615.7259758594629,610.0451525264967,901.5887395006741,726.7109792998381,118.13698325501971,115.04846122528264,59.075799439274874,348.9141314536517),c(917.5705895890643,946.7414813746911,402.26267844536187,968.9479723997841,93.32821677146558,46.019384490113055,375.1313264393955,942.4474094155701),c(134.2560398490558,363.5973241943131,31.020256298025416,118.16438687802246,34.25972514395572,15.339794830037684,17.722739831782462,42.666130204930234),c(836.0910619565338,1384.0995992599746,1283.8383493666647,924.2143116531041,205.5583508637343,635.3231692107274,941.274404399113,422.86875714219747),c(624.9850130904323,471.20130557878764,200.13068579371236,750.3438566754426,2.3627396651003942,2.5566324716729474,147.68949859818719,78.69530682242687))
targetgene="Ptpn2"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(637.0217614906925,1437.0338159087887,2332.5231429257174,3773.664098083131,481.9988916804804,442.2974175994199,409.5922094456391,175.40520195360207),c(1005.5314432832732,1364.1407962612416,961.6279452387879,1427.2569872195427,226.82300784963786,548.3976651738473,172.30441503121838,762.3015263280869),c(50.924704770331516,334.09300671792494,25.016335724214045,108.88004219474927,250.45040450064178,2.5566324716729474,73.84474929909359,280.64832312576334),c(925.9037230969367,568.391998442184,1188.7762736146515,842.3432721733315,41.347944139256896,86.92550403688021,887.1215882464444,327.10699823779845),c(350.917511053739,372.27506462854495,345.22543299415383,238.01683642573093,163.0290368919272,2256.2281562513763,257.96432421816695,325.210725784246),c(1051.82662943812,1401.4550801284383,1283.8383493666647,968.1039410649411,1121.1199710901371,940.8407495756446,498.20590860455144,1173.7926487489697))
targetgene="Zbtb7a"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(1252.7477373501554,886.8650723784916,1135.7416418793177,319.8878759055036,90.96547710636519,801.504279869469,80.73692590034233,361.23990240174265),c(833.313350787243,667.3182393924267,343.2241261362167,348.5849412901662,1161.2865453968438,29.401273424238894,823.1228055205632,694.0357180001985),c(399.99040837787663,412.1926706260113,764.4992197319812,530.895709616258,1033.6986034814224,448.68899877860224,1517.263448932043,818.2415637078843),c(1530.5188542792364,295.0431747638818,109.07122375757324,55.706068099639154,17.720547488252958,24.288008480893,35.445479663564925,77.74717059565064),c(840.7205805720185,96.32291881997317,54.03528516430234,261.64971380133545,49.61753296710828,497.26501574038826,82.70611921498482,36.97731284427287),c(134.2560398490558,424.34150723393583,661.4319165482193,380.6581320142009,35.441094976505916,44.74106825427658,127.99756545176223,775.5754335029541))
targetgene="Kat6a"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(439.8042684710449,550.1687435302972,672.4391042668735,1483.8070866540247,79.15177878086321,196.86070031881695,732.5399130470084,968.9952237653044),c(339.8066663765758,252.52224663614587,1091.7128910047009,106.34794819022021,378.03834641606306,793.8343824544502,1410.9270099413482,546.126466623107),c(112.96025421782628,404.38270423520265,887.5795914951143,752.8759506799717,77.97040894831301,171.29437560208748,271.7486774206644,868.4927837270243),c(2401.7942577134536,1194.9248577937212,533.3482776402434,1188.3961194589688,1907.9122795685682,54.96759814096837,259.93351753280945,322.3663171039173),c(314.8072658529585,496.3667528380599,1129.7377213055063,357.0252546385964,76.78903911576282,127.83162358364737,73.84474929909359,329.95140691812713),c(824.9802172793705,184.83587124913768,1191.778233901557,704.7661645939196,85.0586279436142,80.53392285769785,22.645723118388702,190.57538158202172))
targetgene="Tecr"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(358.3247408385145,420.0026370168199,486.317566478721,264.1818078058645,379.2197162486133,608.4785282581615,1405.0194299974207,1043.8979856806263),c(152.77411431099455,569.2597724856072,199.13003236474378,386.566351358102,406.3912223972678,3171.5025811102914,654.7567771186299,511.99356245916283),c(468.50728388705,324.54749224026995,902.5893929296427,711.5184152726638,220.91615868688686,489.5951183253694,453.8990590250953,415.28366732798764),c(169.4403813267394,69.42192347385453,78.05096745954782,48.953817420895014,73.24492961811222,67.75076049933311,286.5176272804831,220.91574083886098),c(104.62712070995384,9.545514477654999,218.14244751514647,111.41213619927831,22.446026818453745,1.2783162358364737,0.0,0.0),c(86.10904624801512,227.3567993768736,28.01829601111973,26.164971380133544,17.720547488252958,242.88008480893,92.5520857881973,56.888173606573645))
targetgene="Btbd35f27"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(115.73796538711709,87.64517838574135,246.1607435262662,45.57769208152295,7.088218995301183,2.5566324716729474,53.16821949534739,2.8444086803286823),c(46.29518615484683,451.24250258005446,731.4776565760187,258.27358846196336,735.9934056787728,86.92550403688021,300.30198048298064,107.1393936257137),c(99.07169837137222,85.90963029889498,343.2241261362167,33.7612533937207,108.68602459461813,57.52423061264132,12.799756545176223,46.458675112035145),c(265.7343685288208,98.92624095024271,81.0529277464535,167.96223563376049,90.96547710636519,9054.313898429744,665.5873403491636,755.6645727406533),c(49.07289732413764,47.72757238827499,126.08233205003879,131.66888823551074,0.0,386.05150322261505,119.136195535871,267.37441595089615),c(246.29039034378516,250.7866985492995,111.07253061551036,209.31977104106835,259.9013631610434,0.0,10.830563230533727,30.340359256839278))
targetgene="Sys1"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(13.88855584645405,34.71096173692727,172.11238978259263,21.100783371075437,33.07835531140552,379.6599220434327,455.86825233973775,23.703405669405686),c(169.4403813267394,45.99202430142863,63.04116602501939,397.53875871106123,3.5441094976505916,1.2783162358364737,19.691933146424958,45.51053888525892),c(181.47712972699958,327.1508143705395,728.475696289113,252.36536911806223,109.86739442716834,493.43006703287887,451.92986571045276,488.29015678975713),c(196.29158929655057,168.34816442409723,123.0803717631331,152.76967160658617,272.89643131909554,371.99002462841383,67.9371693551661,194.36792648912663),c(159.25544037267312,135.37275077401634,79.05162088851638,143.48532692331298,4.7254793302007885,15.339794830037684,61.04499275391737,6.636953587433592),c(143.51507708002518,348.84516545611905,337.22020556240534,22.788846040761474,72.06355978556202,8.948213650855315,70.89095932712985,52.14749247269251))
targetgene="Elob"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(312.02955468366764,321.0763960665772,158.10324177703276,197.5033323532661,656.8416268979096,29.401273424238894,108.30563230533727,340.38090541266564),c(259.2530424671423,573.598642702723,231.15094209173776,357.0252546385964,1180.1884627176469,648.1063315690922,241.22618104370574,597.3258228690233),c(215.73556748158626,257.72889089668496,284.18557382707155,236.3287737560449,33.07835531140552,8.948213650855315,22.645723118388702,91.96921399729406),c(198.14339674274444,716.7813598675481,253.16531752904615,1011.993570476778,86.23999777616439,29.401273424238894,236.3031977570995,402.9578963798967),c(71.29458667846413,83.30630816862543,86.05619489129631,218.60411572434154,283.5287598120473,393.7214006376339,163.44304511532715,55.94003737979742),c(771.2778013397483,602.2351861356881,679.4436782696534,539.3360229646881,70.88218995301183,63.915811791823685,182.15038160443086,205.74556121044137))
targetgene="Map3k7"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetgenelist=c("Cd19","Maea","Ints10","Arl15","Ints13","Psmd12","Psma7","Unc5c","Tsc2","Kifc3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2480_BM,AR2481_BM,AR2482_BM,AR2483_BM_vs_AR2466_BM,AR2467_BM,AR2468_BM,AR2469_BM pos.'


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
targetmat=list(c(128.7006175104742,400.9116080615099,614.4012053866969,155.30176561111523,1253.4333923357592,4288.7509712313695,521.8362283802614,399.16535147279177),c(94.44217975588754,279.42324198226447,276.18034639532306,409.3551973988635,3084.5566327885645,455.0805799577846,345.593426719758,570.778008519289),c(624.0591093673353,504.17671922886854,641.4188479688481,260.8056824664924,1193.183530875699,2056.8108234608862,1703.352217165759,1207.9255529129139),c(423.13800145530007,113.6783996884368,230.1502886627692,201.72348902748118,10394.873156609185,9134.847821287442,763.0624094239671,708.2577614018419),c(1146.2688091940076,563.1853541816449,213.13918037030368,236.3287737560449,1222.717776689454,16241.007776302398,654.7567771186299,636.1994081668486),c(253.69762012856066,513.7222337065235,542.3541585009605,745.2796686663845,3397.619638414367,13964.32656027764,1980.0238778730295,1697.1638459294472))
targetgene="Cd19"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(329.6217254225095,373.1428386719681,727.4750428601444,107.19197952506322,1156.561066066643,1977.5552168390248,2589.489208754882,486.3938843362047),c(124.07109889498952,359.2584539771972,448.2927361779157,353.64912929922434,465.45971402477767,336.1971700249926,307.19415708422935,475.01624961488994),c(278.69702065217797,234.29899172425905,530.3463173533378,524.9874902723569,1247.5265431730081,773.3813226810665,1656.091577614339,1225.940141221662),c(345.36208871515737,971.9069286339635,521.3404364926207,766.3804520374599,576.5084782844962,2511.891403418671,759.1240227946821,3911.061935451938),c(387.02775625451955,399.1760599746636,180.11761721434112,300.47515520411423,1520.4229744921038,1870.176653028761,251.07214761691822,979.424722259843),c(202.77291535822914,674.2604317398121,237.15486266554916,636.3996264716352,114.59287375736912,276.1163069406783,511.00566514972763,446.5721628116031))
targetgene="Maea"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(1082.381452300319,646.4916623502703,553.3613462196147,502.19864423159544,5130.689182765506,273.55967446900536,635.0648439722049,1489.5220122654532),c(232.40183449733112,966.7002843734243,302.19733554850563,1250.854438237352,3181.428959057681,2676.794197841576,1591.1081982311366,1408.930432989474),c(20.36988190813261,577.0697388764158,1789.1683309957884,440.5843567880551,2328.4799399564386,511.3264943345895,8079.60016997816,862.8039663663669),c(3397.1407600426605,341.0351990653104,804.5253568907237,1749.6769571295754,5722.5554688731545,640.4364341540734,3409.6582243034813,641.8882255275059),c(369.43558551567776,543.2265511829117,490.3201801945953,552.8405243221765,3343.2766261170577,1035.4361510275437,71.8755559844511,28.444086803286822),c(47.22108987794377,633.4750516989226,1080.7057032860466,281.90646583756785,3567.7368943015954,2244.723310128848,7818.682055788029,1499.9515107599918))
targetgene="Ints10"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(4504.521612866597,813.9720527309444,760.496606016107,1271.9552216084273,6504.622298021385,4414.025962343344,4135.305960749241,811.6046101204507),c(1425.8917335692825,1699.969351066013,2127.3891899871624,4693.658253062021,4841.253573790707,2657.6194543040287,4529.14462367774,2090.6403800415815),c(2869.375637877407,404.38270423520265,444.29012246204144,1369.8628564502174,11250.184915375527,1798.5909438219185,2126.7287798138955,676.9692659182264),c(1706.4405616676543,970.1713805471171,839.5482269046233,1709.1634530571105,16947.931617765127,948.5106469906635,1617.6923079788103,951.9287716833323),c(3020.2979447422076,3812.131372758037,2236.4604137447354,2443.470714370536,718.2728581905199,4811.582311688487,5963.701953394799,5889.822240733924),c(1806.4381637621234,2357.7420759807846,939.6135698014796,2980.2746433306947,1571.2218772917622,20377.63911546923,1966.239524670532,3474.9192711348737))
targetgene="Arl15"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(724.9826151849014,78.96743795150952,174.11369664052975,793.3894547524364,354.4109497650591,579.0772548339226,1072.2257598228389,2476.531824339506),c(195.36568557345365,370.53951654169856,244.15943666832908,823.7745828067851,178.38684471507977,1236.1318000538702,2618.0425118171984,556.5559651176455),c(305.5482286219891,418.26708892997357,578.3776819438287,849.0955228520756,5937.564778397291,1090.4037491685121,3851.742123440722,140.32416156288167),c(316.65907329915234,400.9116080615099,1154.7540570297203,719.1146972862509,430.0186190482718,262.0548283464771,1593.0773915457792,2321.037483148205),c(1673.1080276361647,420.8704110602431,590.3855230914514,1201.0565894816139,3098.733070779167,3083.2987608375747,3945.27880588624,621.9773647652052),c(765.7223790011667,219.54683298606497,219.14310094411502,1868.6853753424407,3036.1204696540067,256.9415634031312,677.4025002370186,289.1815491667494))
targetgene="Ints13"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(124.99700261808646,165.7448422938277,103.06730318376187,62.458318778383294,8.26958882785138,159.7895294795592,66.95257269784486,29.39222303006305),c(143.51507708002518,118.88504394897589,377.2463427211478,196.6593010184231,51.980272632208674,148.28468335703096,153.59707854211467,18.96272453552455),c(99.99760209446916,97.19069286339635,351.2293535679652,119.85244954770849,464.27834419222745,934.4491683964623,84.67531252962732,1149.1411068527877),c(92.59037230969366,212.6046406386795,224.14636808895784,121.54051221739452,43.710683804357295,5.113264943345895,78.76773258569983,312.88495483615503),c(417.5825791167184,236.03453981110542,218.14244751514647,165.43014162923143,94.50958660401577,268.4464095256595,200.85771809353457,261.6855985902388),c(148.14459569550988,544.962099269758,245.16009009729765,356.18122330375337,1824.0350214575044,2240.8883614213382,474.5755888288415,1010.7132177434585))
targetgene="Psmd12"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(38.88795637007134,38.182057910619996,233.1522489496749,97.90763484179003,37.80383464160631,804.060912341142,24.614916433031198,120.41330080058088),c(549.9868115195804,62.479731126469076,164.10716235084413,67.5225067874414,841.1353207757403,143.17141841368505,91.56748913087606,27.495950576510594),c(80.5536239094335,85.90963029889498,155.10128149012706,125.76066889160961,1683.4520113840308,393.7214006376339,199.87312143621332,89.12480531696538),c(14.814459569550987,20.82657704215636,4.002613715874247,17.724658031703367,33.07835531140552,147.00636712119447,113.22861559194351,17.066452081972095),c(38.88795637007134,71.1574715607009,484.3162596207839,60.77025610869726,307.1561564630513,320.8573751949549,107.32103564801602,336.58836050556073),c(15.740363292647924,5.20664426053909,2.0013068579371236,124.91663755676659,22.446026818453745,72.864025442679,14.768949859818719,151.7017962841964))
targetgene="Psma7"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(324.06630308392783,251.65447259272267,239.15616952348626,124.91663755676659,2259.960489668527,3330.013794354014,1553.6935252529292,321.4181808771411),c(372.21329668496855,370.53951654169856,348.2273932810595,1100.6168606352949,577.6898481170464,226.26197374305585,554.3279180718625,316.67749974325994),c(1544.4074101256904,438.22589192870674,450.29404303585284,275.1542151588237,414.6608112251192,282.5078881198607,1899.2869519726871,1152.9336517598927),c(124.07109889498952,118.88504394897589,478.3123390469725,193.283175679051,187.83780337548134,357.9285460342126,1212.0384851624563,1525.5511888829499),c(307.400036068183,906.8238753772248,744.48615115261,156.98982828080125,147.67122906877464,478.09027220284116,220.54965123995953,648.5251791149395),c(287.9560578831473,400.0438340180867,388.253530439802,344.3647846159511,342.59725143955717,48.576016961786,188.05796154835835,249.35982764214782))
targetgene="Unc5c"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(864.7940773725388,168.34816442409723,493.322140481501,913.2419043001449,523.3468358197373,210.92217891301817,1172.6546188696063,830.5673346559753),c(634.2440503214017,248.18337641902994,442.28881560410434,1241.5700935540788,1196.7276403733497,464.02879360863994,1713.1981837389712,955.7213165904373),c(302.7705174526983,268.14217941776315,172.11238978259263,395.0066647065322,3346.8207356147086,329.80558884581023,280.6100473365557,732.9093032980238),c(190.73616695796895,524.1355222276018,273.1783861084174,156.98982828080125,138.22027040837307,363.04181097755855,715.8017698725472,302.4554563416166),c(1108.3067565470333,848.6830144678717,566.369840796206,826.3066768113141,2280.0437768218803,4081.6637410258604,972.7814974333929,1133.0227909975918),c(345.36208871515737,171.81926059778996,375.24503586321066,1093.8646099565508,3271.2130663314956,3099.9168719034487,4071.30717802336,588.7925968280373))
targetgene="Tsc2"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(171.2921887729333,223.88570320318087,326.2130178437511,169.65029830344653,64.97534079026084,30.67958966007537,15.753546517139966,103.3468487186088),c(811.0916614329166,262.0677611138009,289.1888409719144,604.3264357476005,1058.5073699649765,2141.179695026093,776.8467626264646,1340.6646246615855),c(667.5765843528914,606.574056352804,300.19602869056854,684.5094125576873,404.0284827321674,62.63749555598721,640.9724239161324,1089.4085245658853),c(165.73676643435167,589.2185754843404,87.05684832026488,178.09061165187669,135.85753074327266,915.2744248589152,213.6574746387108,895.9887343035349),c(160.18134409577004,441.6969881023995,428.27966759854445,455.77692081522946,138.22027040837307,515.1614430420989,297.34819051101687,609.6515938171142),c(127.77471378737727,123.22391416609179,113.07383747344748,420.3276047518227,647.3906682375081,278.67293941235124,338.7012501185093,546.126466623107))
targetgene="Kifc3"
collabel=c("AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
Sweave("BM_cd19.15m_v_egfr.15m_summary.Rnw");
library(tools);

texi2dvi("BM_cd19.15m_v_egfr.15m_summary.tex",pdf=TRUE);

