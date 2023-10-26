pdf(file='cd19_1to10_v_Input1_DOI.pdf',width=4.5,height=4.5);
gstable=read.table('cd19_1to10_v_Input1_DOI.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Kdsr","Sys1","Fitm2","Elob","Ptpn2","Azin1","Gtpbp4","Rnaseh2b","Zfp592","Ndufaf8")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 neg.'


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
targetmat=list(c(773.6592155313073,796.718260303653,857.0040057557965,286.07610119305826,254.88049990612438,258.35901209376067),c(1288.1740434212336,1258.6193375962616,1259.78029712071,547.4104855261763,558.1324815462577,532.323602099225),c(896.3125057984657,1022.1699766012357,886.6275523465063,330.92046300169983,399.9949451081514,329.4510892470774),c(1812.1237397932491,1841.4945065607437,1765.4594345375624,746.890577709444,794.4085654008402,849.6370196371995),c(937.1969358875186,921.358233541341,957.3082950892523,360.30125177287874,370.2278794256843,379.7357291847892),c(942.857856976772,1009.3393911208855,1030.5875945504818,669.5727125221309,749.7579668771396,629.4249757720478))
targetgene="Kdsr"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(789.3839963347891,728.2884710751183,832.5775726020534,177.83108993081999,156.27709483295217,197.67065354824643),c(1025.8846996191562,1049.0531080838743,1056.573161735315,771.6322945693842,798.1294486111486,749.0677397617759),c(783.7230752455357,843.7637403982704,827.3804591650867,259.7880270293718,271.6244743525121,260.0929651950611),c(732.7747854422545,788.1645366500861,745.266066861014,269.0661708518494,279.06624077312887,260.0929651950611),c(858.5730318701093,771.0570893429525,823.74247975921,306.17874614175963,275.3453575628205,329.4510892470774),c(575.5269774074359,587.1520307912657,589.872375095712,153.08937307087982,150.6957700174896,161.25763842093787))
targetgene="Sys1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(745.9836013171793,711.1810237679847,701.6103139904943,340.1986068241774,381.39052905660947,273.96459000546434),c(1088.7838228330836,1002.0076279892568,1144.4043788200509,542.7714136149375,632.5501457524255,511.51616488362015),c(822.0915404060314,877.9786350125378,853.8857376936165,412.87740010025163,431.62245239577265,494.17663387061606),c(942.2288657446328,974.513516245649,928.7241711859359,340.1986068241774,331.1586057174463,341.58876095618024),c(867.3789091200591,899.9739244074239,748.9040462668906,225.76816634695408,228.83431743396568,201.13855975084724),c(490.6131610686339,488.17322851427815,484.3709723252895,184.01651914580503,173.0210692793399,218.47809076385133))
targetgene="Fitm2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(508.22491556853356,469.84382068520637,482.2921269505028,77.31786518731305,109.76605470409736,84.96370196371996),c(810.7696982275245,829.7111943959821,783.724706294567,265.97345624435684,234.41564224942826,239.2855279794562),c(1096.331717618755,1088.155844785894,1225.47934843673,459.2681192126395,470.69172610401074,449.0938532368055),c(995.0641292443318,912.1935296268051,989.0106870547487,539.6786990074451,517.2027662328655,650.2324129876528),c(748.4995662457363,802.8280629133435,730.1944378938108,199.48009218326763,215.81122619788633,223.67995006775254),c(661.6987762105165,672.0782870659649,577.9190141906887,245.87081129565547,264.1827079318953,234.08366867555497))
targetgene="Elob"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(811.3986894596637,788.1645366500861,732.7929946122941,256.6953124218793,279.06624077312887,279.1664493093656),c(1235.9677711536738,1258.6193375962616,1240.5509774039333,763.9005080506528,768.3623829286815,735.1961149513727),c(796.9318911204605,903.6398059732383,850.2477582877399,239.68538208067042,258.60138311643277,225.41390316905293),c(579.9299160324108,574.9324255718846,587.2738183772286,81.95693709855182,79.99898902163028,72.82603025461711),c(1364.2819825100858,1289.7793309056835,1441.1595560708452,1000.4931755238307,1036.2659740708852,1010.8946580581375),c(710.7600923173799,670.2453462830578,741.1083761114407,149.9966584633873,199.0672517514986,164.7255446235387))
targetgene="Ptpn2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(661.6987762105165,669.0233857611196,679.2627262115379,340.1986068241774,277.2057991679747,275.69854310676476),c(652.8928989605666,685.5198528072842,731.2338605812041,208.7582360057452,284.6475655885915,220.2120438651517),c(376.136756819286,359.2563934498067,402.2565800212168,89.68872361728313,87.44075544224705,98.83532677412322),c(614.5244338000709,760.0594446455094,693.8146438350444,184.01651914580503,124.6495875453309,157.78973221833706),c(705.0991712281264,662.30260289046,652.2377363393115,301.53967423052086,359.06522979475915,365.86410437438593),c(586.8488195859428,689.7967146340676,610.1411174998817,184.01651914580503,208.36945977726955,242.75343418205702))
targetgene="Azin1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1146.0220249577576,1107.7072131369039,1096.5909351999578,641.7382810546982,705.107368353439,683.1775219123605),c(594.3967143716142,630.5316293200689,667.3093653065147,281.4370292818195,299.531098429825,279.1664493093656),c(765.4823295134967,712.4029842899228,747.3449122358006,358.7548944691325,396.274061897843,416.14874431209773),c(678.6815394782769,690.4076948950367,552.4531583495523,310.8178180529984,359.06522979475915,423.0845567172994),c(429.6010115511243,369.6430578862807,374.7118788052937,126.80129890719338,102.32428828348058,123.11067019232891),c(506.3379418721157,430.13010372221754,471.3781887328729,129.8940135146859,148.8353284123354,173.3953101300407))
targetgene="Gtpbp4"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(895.0545233341873,812.6037470888484,861.6814078490664,408.23832818901286,377.6696458463011,469.90129045241036),c(703.2121975317086,726.4555302922112,598.1877565948586,182.4701618420588,171.16062767418572,156.05577911703665),c(756.0474610314076,722.7896487263968,703.6891593652811,484.0098360725796,500.4587917864778,393.60735399519245),c(832.1554001202597,815.0476681327247,959.9068518077356,304.63238883801336,336.73993053290883,313.8455113353737),c(879.9587337628446,811.9927668278794,745.266066861014,208.7582360057452,245.57829188035342,185.5329818391436),c(1081.8649192795517,979.4013583334016,1003.5626046782553,312.3641753567447,299.531098429825,265.2948244989623))
targetgene="Rnaseh2b"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(698.1802676745943,739.8970960335305,711.484829520731,154.6357303746261,180.46283569995668,178.59716943394196),c(670.5046534604663,755.171602557757,728.6353038627208,244.3244539919092,182.32327730511088,301.70783962627087),c(864.2339529593628,947.6303847630105,890.7852430960795,388.13568324031144,344.1816969535256,372.79991677958753),c(652.2639077284274,585.9300702693275,639.7646640905915,375.7648248103414,334.87948892775466,310.3776051327729),c(708.873118620962,680.6320107195318,717.7213656450909,296.9006023192821,258.60138311643277,317.31341753797454),c(470.4854416401771,505.28067582141176,506.19884876054925,230.40723825819285,264.1827079318953,246.22134038465782))
targetgene="Zfp592"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(680.5685131746947,742.3410170774067,639.2449527468949,126.80129890719338,111.62649630925155,123.11067019232891),c(705.7281624602657,683.686912024377,650.6786023082215,445.3509034789231,450.2268684473146,452.5617594394063),c(821.4625491738922,761.2814051674475,778.5275928576004,153.08937307087982,126.5100291504851,147.3860136105346),c(718.9369783351905,778.9998327355502,661.0728291821547,344.83767873541615,355.34434658445076,332.91899544967816),c(667.9886885319092,622.5888859274711,706.2877160837643,561.3277012598927,548.8302735204868,617.2873040629449),c(749.7575487100148,694.6845567218201,762.4165412030039,259.7880270293718,200.92769335665278,195.93670044694602))
targetgene="Ndufaf8"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetgenelist=c("Jak2","Nprl3","Jak1","Tsc1","Irf1","Stat1","Nprl2","Tsc2","Zfp36l2","Mark2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 pos.'


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
targetmat=list(c(537.1585122469402,612.8132017519662,598.7074679385552,1516.9765149750817,1516.2599082006668,1447.85083958584),c(857.94404063797,902.4178454513002,821.1439230407267,2726.2279265046577,2816.7085902034473,2765.6551965741496),c(718.3079871030511,801.6061023914054,710.9651181770342,2447.8836118303307,2465.085126829305,2509.0301375816894),c(1234.080797457256,1352.7102977854966,1424.0090817288553,4300.419661718352,4478.082943606141,4353.9562373653225),c(567.9790826217646,637.8633924516976,620.5353443738151,2336.5458859606,2212.065068528335,2330.4329681477475),c(993.8061467800533,930.5229374558769,1002.523181990862,4586.49576291141,5084.5869068864085,4912.289135984054))
targetgene="Jak2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(953.5507079231397,969.0146938969276,1192.2178224401437,3168.4861153760885,3380.4223965651677,3013.6104900601076),c(1017.7078136013457,1069.8264369568224,1073.2039247336081,2256.1353061657946,2390.6674626231375,2354.708311565953),c(1172.4396567076071,1242.1228705500969,1230.1567505300002,3661.7740952711456,3622.2798052352127,3519.9247956398267),c(1097.5897000830334,1092.4327066126775,1072.6842133899113,2602.519342204957,2658.5710537653413,2561.0487306207015),c(1274.9652275463088,1277.5597256863023,1132.9707292587243,2948.9033782441193,3136.7045462899687,3093.3723327199264),c(1059.2212349225379,1138.2562261853568,1006.6808727404352,2192.7346567121976,2238.1112510004937,2158.771611119007))
targetgene="Nprl3"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(964.8725501016467,909.7496085829289,976.5376148060288,3372.6052794705947,3276.237666676533,3155.794644366741),c(1001.3540415657245,1088.155844785894,1137.6481313519942,4178.257434722396,4359.014680876273,4315.809269136714),c(882.4746986914017,928.0790164120007,902.2188926574061,2945.810663636627,2943.218619353933,2767.38914967545),c(940.9708832803542,852.9284443128064,973.4193467438489,2557.6749803963153,2340.4355392839743,2427.5343418205703),c(1517.1268519199293,1573.885152256296,1467.664834599375,4348.356738134486,4411.10704582059,4579.370140534375),c(779.3201366205608,694.073576460851,850.2477582877399,1544.8109464425145,1774.8612913170996,1678.4666020587943))
targetgene="Jak1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(936.5679446553793,929.3009769339387,1010.8385634900086,1673.1586026534542,1836.255864287188,1792.907506744621),c(971.1624624230394,939.6876413704127,977.0573261497254,2738.598784934628,2600.897364005561,2585.3240740389074),c(854.7990844772737,934.7997992826603,922.4876350615759,2103.0459330949147,2234.3903677901853,2372.047842578957),c(661.0697849783772,650.0829976710788,664.1910972443346,1854.0824071917666,1929.2779445448978,2087.6795339656906),c(970.5334711909001,950.6852860678558,906.896294750676,1999.4399937439152,2106.019897034546,1978.4404885837646),c(682.4554868711125,673.300247587903,691.7357984602578,1329.8672812217842,1248.3563170584632,1146.1429999595691))
targetgene="Tsc1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(971.7914536551787,955.5731281556083,964.5842539010056,2466.439899475286,2321.831123232432,2401.525045301064),c(1226.5329026715847,1083.2680026981416,1131.4115952276343,2848.3901535006125,2775.7788748900552,2817.6737896131617),c(1000.7250503335853,1043.5542857351527,1044.100089486595,2202.0128005346755,2392.5279042282914,2281.882281311336),c(1023.9977259227384,1103.4303513101206,1066.967388609248,2669.0127062660463,2492.991750906618,2479.5529348595824),c(1098.2186913151727,1237.8460087233136,1156.8774510687708,1425.7414340540524,1609.2819884583766,1492.9336202196507),c(885.619654852098,835.8209970056727,905.3371607195861,2232.9399466096006,2212.065068528335,2496.8924658725864))
targetgene="Irf1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(948.5187780660256,928.6899966729696,910.5342741565527,2941.171591725388,3064.147323688955,2890.499819867779),c(1025.8846996191562,993.4539043356899,1014.9962542395818,3819.502540253264,3962.74061897843,3884.0549469129123),c(764.2243470492182,848.6515824860229,813.8679642289735,1464.400366647709,1454.8653352305785,1425.3094492689347),c(800.705838513296,898.1409836245167,958.8674291203423,3290.648342372043,3151.588079131202,3108.97791063163),c(799.4478560490176,895.6970625806405,829.4593045398733,2331.9068140493614,2368.342163361287,2441.4059666309736),c(911.4082953698083,961.6829307652988,926.1256144674526,1302.0328497543517,1307.8904484233974,1354.217372115618))
targetgene="Stat1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1110.169524725819,1036.222522603524,1080.4798835453614,2064.387000501258,2074.392389746925,2160.5055642203074),c(597.5416705323105,534.6077283479266,565.9656532856654,963.3806002339205,961.8483098647175,924.197002993117),c(920.8431638518974,975.7354767675872,928.7241711859359,4334.439522400769,4132.040805047462,4149.349771411875),c(685.6004430318089,700.1833790705416,712.0045408644277,1091.72825644486,1118.1254046976696,1080.2527821101537),c(981.2263221372677,968.4037136359585,1050.8563369546516,2724.6815692009113,2595.3160391900988,2665.085916698726),c(947.260795601747,878.5896152735069,964.064542557309,2438.6054680078532,2429.7367363313756,2583.590120937607))
targetgene="Nprl2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1052.931322601145,1076.547219827482,1142.8452447889608,1776.7645420044537,1756.2568752655577,1813.714943960226),c(1032.174611940549,1059.4397725203482,1075.8024814520913,2463.3471848677937,2319.9706816272783,2314.827390236044),c(1095.0737351544765,1089.3778053078322,1198.9740699082004,4156.608432469949,4442.734553108212,4633.122686674688),c(708.2441273888228,791.8304182159004,692.7752211476511,774.7250091768767,786.9667989802234,802.8202859020886),c(1347.9282104744645,1329.4930478686724,1346.052380174356,4687.008987654916,4569.244582258697,4643.5264052824905),c(1074.9460157260196,1179.8028839312528,1112.1822755108578,2797.3603624769858,2757.174458838513,2828.0775082209643))
targetgene="Tsc2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(863.6049617272234,714.846905333799,788.402108387837,890.7018069578462,1010.2197915987266,847.9030665358991),c(1127.7812792257187,1114.4279960075637,1139.2072653830842,2370.565746643018,2394.388345833446,2242.0013599814265),c(1155.4568934398467,1240.2899297671897,1112.7019868545544,3098.900036707507,3004.613192324021,2923.4449287924867),c(889.3936022449337,1000.7856674673186,946.3943568716224,1577.284449821186,1598.1193388274514,1583.0991814872718),c(1099.4766737794514,1176.1370023654385,1164.153409880524,1843.2579060655428,1860.4416051541925,1895.2107397213451),c(790.6419787990677,759.4484643845403,745.7857782047107,1315.9500654880678,1378.5872294192566,1300.4648259753055))
targetgene="Zfp36l2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1181.245533957557,1103.4303513101206,1221.8413690308535,1919.0294139491095,1944.1614773861313,1947.2293327603572),c(629.6202233714135,620.1449648835949,661.5925405258514,1473.6785104701864,1439.981802389345,1492.9336202196507),c(737.8067152993686,794.2743392597768,841.4126654448966,1260.2812025532025,1190.6826272986832,1191.2257805933798),c(1159.8598320648216,1187.1346470628816,1257.1817404022265,2352.0094589980627,2496.7126341169264,2470.88316935308),c(945.3738219053291,860.8711877054042,913.6525422187327,1275.744775590665,1462.3071016511954,1404.50201205333),c(498.1610558543052,508.3355771262571,600.2666019696452,773.1786518731304,703.2469267482847,764.6733176734796))
targetgene="Mark2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
Sweave("cd19_1to10_v_Input1_DOI_summary.Rnw");
library(tools);

texi2dvi("cd19_1to10_v_Input1_DOI_summary.tex",pdf=TRUE);

