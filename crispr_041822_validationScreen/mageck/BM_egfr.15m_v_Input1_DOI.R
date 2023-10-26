pdf(file='BM_egfr.15m_v_Input1_DOI.pdf',width=4.5,height=4.5);
gstable=read.table('BM_egfr.15m_v_Input1_DOI.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Cxcr4","H2-K1","Rasgrp2","Ndufaf8","Hdac7","Stk17b","Supt4a","Lsm4","Psma7","Cnot9")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2466_BM,AR2467_BM,AR2468_BM,AR2469_BM_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 neg.'


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
targetmat=list(c(832.7873206053545,909.990383198323,825.9652912663685,73.52110073353785,413.11825652953837,192.5703491492203,459.1110055320198),c(603.7586250614595,546.9396744729661,667.6384481483165,152.42179420367603,49.10012065310087,40.11882273942089,115.20128552094961),c(786.0069913027717,871.6998787624456,815.0601260516047,17.931975788667767,37.24836739200756,10.029705684855223,254.12048276680062),c(530.6643605261738,587.5937902937743,595.7451367324662,68.14150799693752,84.65538043638081,248.73670098440954,210.0729324205552),c(740.2012521939927,719.956027849894,746.8018697073425,139.86941115160857,612.904954359397,343.01593442204864,103.34232965849893),c(617.4028877747128,603.193625434317,622.8061022653986,166.76737483461022,18.62418369600378,138.4099384510021,272.755984836366))
targetgene="Cxcr4"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(827.9143696363354,815.4459278010946,811.0211759720626,64.55511283920396,704.3327652306883,116.3445859443206,237.17911724901393),c(701.2176444418403,692.0654135077117,691.8721486255693,227.73609251608065,79.57605761019796,381.1288160244985,171.10779172964575),c(625.1996093251432,595.1573467255525,683.3903534585307,86.07348378560528,159.15211522039593,1484.3964413585732,315.1093986308328),c(739.2266620001889,722.3196392348248,684.1981434744392,159.59458451914313,71.11051956655989,14.041587958797313,1.6941365517786708),c(492.1680478709234,472.7222769861418,453.5740939325828,17.931975788667767,35.55525978327994,20.059411369710446,15.247228966008038),c(531.6389507199776,519.0490601307837,544.0465757143268,30.484358840735204,30.47593695709709,14.041587958797313,15.247228966008038))
targetgene="H2-K1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(639.3311671352984,636.2841848233469,668.8501331721791,213.39051188514642,316.6111228320642,686.0318688440973,144.00160690118702),c(656.873790623767,557.8122868436474,570.7036462393049,511.0613099770313,108.35888695856744,116.3445859443206,216.84947862766987),c(555.5164104681709,522.3581160696867,482.250639497332,39.45034673506909,99.89334891492936,116.3445859443206,870.7861876142368),c(877.1311744234278,837.6638748194433,852.6223617913466,136.28301599387504,186.24183696003777,58.1722929721603,249.03807311146463),c(370.3442736454473,329.48742705934086,343.71465176903644,98.62586683767272,57.56565869673895,429.27140331180357,37.27100413913076),c(542.3594428518195,550.7214526888553,546.0660507540979,77.1074958912714,82.9622728276532,272.8079946280621,181.27261104031777))
targetgene="Rasgrp2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(527.2532948478605,574.3575665381624,496.79085978368374,23.311568525268097,82.9622728276532,6.0178234109131346,33.882731035573414),c(546.7450987239366,528.9762279474927,505.67654995867645,437.5402092434935,76.18984239274273,443.3129912706009,708.1490786434844),c(636.407396553887,589.0119571247327,605.0347219154132,177.52656030781088,30.47593695709709,22.065352506681492,47.43582344980278),c(556.9782957588766,602.7209031573308,513.7544501177607,139.86941115160857,204.86602065604157,22.065352506681492,52.5182331051388),c(517.5073929098223,481.7040002488785,548.8933158097774,283.32521746095074,126.98307065457122,62.184175246102384,54.21236965691747),c(580.85575550707,537.4852289332432,592.5139766688326,1335.9321962557485,611.2118467506695,443.3129912706009,1294.3203255589044))
targetgene="Ndufaf8"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(1244.0643823905616,1133.588020212768,1113.1346419218148,735.2110073353784,895.6539250169089,778.3051611447653,443.8637765660118),c(6.82213135662666,7.563556431778269,7.2701101431758595,0.0,0.0,0.0,0.0),c(695.8573983759194,642.9022967011529,740.339549580075,390.91707219295733,120.21064021966075,122.36240935523372,67.76546207114683),c(935.1192909547543,904.3177158744893,924.919568215151,147.0422014670757,174.39008369894447,212.62976051893074,355.76867587352086),c(992.6201123891791,964.8261673287154,1012.9686799491698,157.80138694027636,550.2599728364753,60.17823410913134,218.54361517944855),c(826.4524843456297,920.8629955690043,892.6079675788138,145.2490038882089,431.7424402255421,320.95058191536714,470.9699613944705))
targetgene="Hdac7"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(879.5676499079373,865.0817668846396,875.6443772447368,283.32521746095074,123.59685543711599,180.53470232739403,310.02698897549675),c(442.951243083831,504.39466954421334,493.96359472800424,331.7415520903537,204.86602065604157,70.20793979398657,103.34232965849893),c(449.2860793435558,510.5400591450332,480.2311644575609,37.65714915620231,209.94534348222442,54.16041069821821,108.42473931383493),c(534.0754262044871,552.6123417967998,537.1803605791051,68.14150799693752,201.47980543858634,110.32676253340746,120.28369517628563),c(760.1803511669707,800.7915372145243,786.3835804868554,112.97144746860693,96.50713369747413,12.035646821826269,130.44851448695766),c(441.9766528900272,369.19609832617675,410.76122308943604,181.11295546554445,198.0935902211311,246.7307598474385,386.26313380553694))
targetgene="Stk17b"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(434.66722643649865,429.2318275034168,435.39881857464314,68.14150799693752,242.1143880480491,288.85552372383046,55.906506208696136),c(543.8213281425252,612.1753486970537,553.3361608972738,313.8095763016859,182.85562174258254,395.1704039832958,181.27261104031777),c(455.13362050637863,399.4503240532899,373.6028823576483,208.0109191485461,40.63458260946279,322.9565230523382,52.5182331051388),c(750.4344492289326,757.7738100087854,824.7536062425058,728.0382170199114,485.92188370482586,176.52282005345194,177.88433793676043),c(479.01108025457194,509.5946145910609,452.76630391667436,295.87760051301814,226.87641956950057,180.53470232739403,164.33124552253108),c(383.01394616489677,350.2872072467311,360.27434709515927,168.56057241347702,32.169044565824706,417.2357564899773,20.32963862134405))
targetgene="Supt4a"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(479.49837535147384,506.75828092914406,550.10500083364,272.5660319877501,108.35888695856744,190.56440801224926,121.9778317280643),c(516.0455076191166,512.9036705299638,527.486880388204,378.3646891408899,128.67617826329882,140.41587958797314,123.67196827984297),c(472.1889488979453,499.6674467743519,491.13632967232473,114.7646450474737,11.851753261093315,98.29111571158118,16.941365517786707),c(252.41886019518643,309.6330914259229,314.2303161883788,28.691161261868427,22.01039891345901,44.130705013362984,22.023775173122722),c(754.3328100041479,807.4096490923303,735.8967044925787,1228.340341523742,912.5850011041852,2310.8441897906437,1026.6467503778745),c(297.2500091101616,319.0875369656457,295.24725081453073,209.80411672741286,8.46553804363808,56.16635183518925,20.32963862134405))
targetgene="Lsm4"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(518.4819831036261,568.6848992143287,576.3581763506639,75.31429831240462,74.49673478401512,467.38428491425344,196.5198400063258),c(442.951243083831,492.10389034257366,487.905169608691,1065.1593618468653,121.90374782838836,328.9743464632513,135.53092414229366),c(602.7840348676556,639.5932407622499,642.5969576551552,156.00818936140956,167.617653264034,310.9208762305119,252.42634621502197),c(275.3217297495759,296.39686767031094,281.91871555204165,28.691161261868427,40.63458260946279,8.023764547884179,35.57686758735209),c(491.19345767711957,460.43149778450214,429.3403934553299,75.31429831240462,138.83482391566454,970.8755102939856,121.9778317280643),c(258.26640135800926,245.81558403279377,284.34208559976696,30.484358840735204,10.158645652365697,4.011882273942089,250.73220966324328))
targetgene="Psma7"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(564.2877222124052,487.3766675727122,531.121935459792,139.86941115160857,487.61499131355345,306.9089939565698,106.73060276205626),c(728.0188747714451,729.883195666603,687.4293035380729,487.74974145176327,345.3939521804337,92.27329230066806,147.38988000474436),c(503.37583509966714,444.8316626439595,411.5690131053445,68.14150799693752,113.43820978475028,22.065352506681492,120.28369517628563),c(496.5537037430405,523.303560623659,533.141410499563,109.38505231087338,11.851753261093315,310.9208762305119,271.0618482845873),c(374.2426344206625,423.0864379025969,396.22100280308433,41.243544313935864,23.70350652218663,66.19605752004448,404.89863587510234),c(577.4446898287566,566.3212878293979,634.1151624881167,30.484358840735204,336.92841413679565,208.61787824498865,350.68626621818487))
targetgene="Cnot9"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetgenelist=c("Arl15","Pak2","Sipa1","Cnot8","Zfp36l2","Nosip","Pou2af1","Setd5","Nprl3","Zbtb7a")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2466_BM,AR2467_BM,AR2468_BM,AR2469_BM_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 pos.'


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
targetmat=list(c(761.1549413607745,777.6281456422033,689.448778577844,8723.906221186868,1588.134936986504,1524.515264097994,2553.063783530457),c(1045.2479828545847,950.1717767421451,983.0804493605579,2761.5242714548363,3316.7978054974,4264.630857200441,9421.093364441189),c(812.8082216323764,818.2822614630115,824.7536062425058,5557.119296908141,788.9881456670691,890.6378648151439,2749.583623536783),c(857.1520754504497,892.9723812268219,763.3615650334652,3304.863137851469,1892.894306557475,1682.9846139187066,3430.6265173518086),c(1074.485688668699,1084.8976256831954,1042.8569105377817,5849.410502263426,7437.8217251404185,4483.278441130285,4904.525317399252),c(850.8172391907249,911.4085500292814,862.3158419822478,3498.528476369081,4600.173372912934,1883.578727615811,5981.996164330487))
targetgene="Arl15"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(1184.1270854716274,1087.7339593451125,1200.7798586478793,6511.100408865266,5368.844227275271,6836.24739479732,6859.558898151839),c(891.7500273304848,902.8995490435309,912.3988229685704,971.913087745793,1471.3105119842985,1961.8104319576817,3166.341215274336),c(924.3987988229125,918.4993841840736,927.7468332708305,2996.433154286384,2294.16080982592,1542.5687343307334,3183.2825807921226),c(924.3987988229125,953.4808326810481,881.7028023640501,5472.839010701402,3853.5129174640547,673.996222022271,1705.9955076411215),c(517.5073929098223,562.5395096135088,481.8467444893778,2542.754166833089,2368.657544609935,890.6378648151439,2729.2539849154387),c(777.2356795585374,725.1559728967416,759.3226149539231,7665.9196496554705,4273.403604428503,2427.188775734964,5299.259133963682))
targetgene="Pak2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(656.3864955268651,618.7934605748596,599.3801918040542,281.53201988208394,1818.3975717734597,708.0972213507788,2017.716633168397),c(898.0848635902097,862.7181554997088,838.8899315209034,1466.8356195130234,2278.9228413473716,2974.8107061280593,4870.642586363679),c(494.6045233554329,504.8673918211995,527.8907753961582,1708.9172926600381,707.7189804481436,1570.651910248328,1229.9431365913151),c(742.6377276785022,704.8289149863375,670.465713203996,2114.1799454839297,3958.485589205167,2009.9530192449868,1519.6404869454677),c(799.651254016025,726.1014174507138,839.6977215368117,5474.632208280269,3460.7119522392477,2595.687831240532,4611.439693941542),c(449.7733744404577,444.8316626439595,449.5351438530406,1542.149917825428,4049.913400076458,2565.5987141859664,857.2330952000075))
targetgene="Sipa1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(814.757402019984,901.9541044895586,881.2989073560958,1450.6968413032223,1632.155734813422,1345.986502907571,1916.0684400616767),c(608.6315760304785,689.2290798457948,676.5241383233091,1707.1240950811714,927.8229695827337,3426.1474619465444,1018.1760676189812),c(627.1487897127508,586.648345739802,557.3751109768159,4952.811712830037,7121.210602308354,6577.480988128055,1973.6690828221515),c(854.7155999659402,805.0460377073996,881.2989073560958,4500.925922955609,1796.3871728600009,3560.5455181236043,2417.532859388163),c(634.4582161662794,608.8662927581507,605.4386169233674,2736.419505350701,1921.6771359058444,1129.3448601146981,2344.68498766168),c(759.2057609731669,687.3381907378503,738.7239695482582,7429.2175692450555,4601.866480521661,3783.2049843273903,3799.948285639559))
targetgene="Cnot8"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(669.0561680463146,553.0850640737859,612.7087270665432,2571.445328094958,1537.3417087246755,7175.251446945427,138.919197245851),c(873.7201087451144,862.2454332227227,885.3378574356379,6017.971074676902,2146.8604478666175,4936.621138085741,1536.5818524632546),c(895.1610930087983,959.626222281868,864.739212029973,2512.269807992354,1996.1738706898595,2932.6859422516673,7123.84420022931),c(689.0352670192927,774.3190897033003,735.4928094846244,792.5933298591153,1571.203860899228,1121.321095566814,1328.2030565944779),c(851.7918293845287,909.990383198323,904.7248178174402,4233.739483704459,5059.005534878117,7941.520961268366,2690.2888442245294),c(612.5299368056938,587.5937902937743,579.5893364142977,2501.5106225191535,656.9257521863151,461.3664615033403,775.9145407146312))
targetgene="Zfp36l2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(517.9946880067242,617.3752937439012,551.7205808654569,1144.0600553170036,1083.5888695856743,655.9427517895316,3108.740572513861),c(459.51927637849576,561.5940650595365,512.5427650938981,667.0694993384409,1791.3078500338179,1470.3548533997757,911.4454648569249),c(492.1680478709234,541.2670071491324,476.596109385973,1603.1186355068983,372.48367392007555,1693.0143196035617,1802.5612910925058),c(547.7196889177404,554.9759531817305,634.1151624881167,2341.9160380000103,965.0713369747413,1295.8379744832948,2212.5423366229443),c(809.397155954063,775.7372565342588,801.3276957811614,1301.8614422572798,1991.0945478636768,633.8773992828501,1645.0065917770894),c(394.2217333936406,427.3409383954722,419.24301825647456,277.9456247243504,2038.50156090805,351.03969896993283,2198.989244208715))
targetgene="Nosip"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(923.4242086291086,949.2263321881728,894.2235476106307,1687.3989217136368,1175.0166804569656,2477.3373041592404,3056.2223394087223),c(709.5016610891727,747.3739199150903,743.5707096437087,785.4205395436481,1212.2650478489732,1843.45990487639,4186.211419445096),c(666.1323974649032,739.3376412063258,714.8941640789595,1794.9907764456434,961.685121757286,970.8755102939856,892.8099627873595),c(514.096327231509,600.8300140493863,508.90771002231014,2465.646670941818,2112.998295692065,2489.3729509810664,804.7148620948686),c(595.474608414127,559.2304536746058,537.5842555870594,15498.60667414555,2553.206273961245,2338.927365708238,3205.306355965245),c(541.3848526580157,574.8302888151485,583.2243914858856,936.0491361684574,5744.714116412802,1674.9608493708224,2742.807077329668))
targetgene="Pou2af1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(743.612317872306,675.9928560901828,696.7188887210199,1534.9771275099608,1425.596606548653,2112.25601723051,3784.7010566735507),c(724.6078090931317,720.9014724038663,674.5046632835381,1744.7812442373736,7078.882912090164,3570.5752238084597,3579.7105339083314),c(537.4864918828005,618.7934605748596,618.7671521858565,832.0436765941844,702.6396576219607,411.21793307906415,1905.9036207510046),c(783.5705158182622,755.4101986238546,815.8679160675131,7622.882907762668,1295.2273206766265,6152.221467090194,7186.527252645122),c(440.5147675993215,424.50460473355537,481.8467444893778,1104.6097085819345,1263.0582761108017,631.871458145879,886.0334165802449),c(733.8664159342679,679.3019120290858,711.2591090073715,3946.8278710857753,2324.636746783017,2581.6462432817343,2712.312619397652))
targetgene="Setd5"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(738.739366903287,749.7375313000209,926.5351482469679,3057.4018719678543,1904.7460598185683,1715.079672110243,1709.383780744679),c(788.4434667872812,827.7367070027343,834.0431914254527,973.7062853246597,416.5044717469936,784.3229845556784,432.00482070356105),c(908.3180606251497,961.0443891128264,956.0194838276254,1375.3825429908177,1031.1025337151184,884.6200414042307,1626.371089707524),c(850.329944093823,845.2274312512216,833.6392964174985,1651.5349701363014,2365.27132939248,3139.297879359685,2170.188922828477),c(987.7471614201601,988.4622811780225,880.4911173401874,2711.314739246566,2199.3467837371736,1181.4993296759453,2363.3204897312457),c(820.6049431828069,880.6816020251822,782.3446304073133,405.26265282389153,2756.3791870085593,2240.6362499966567,4604.663147734427))
targetgene="Nprl3"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
targetmat=list(c(502.88854000276524,448.14071858286246,455.99746398030805,1233.7199342603424,2803.7862000529326,4675.848790279505,7574.4845230024375),c(896.622978299504,873.5907678703901,849.3912017277129,1947.4125706493194,2661.5651609198126,1927.709432629174,2864.7849090577324),c(323.5639443428645,334.2146498292023,339.2718066815401,98.62586683767272,651.8464293601322,50.148528424276115,218.54361517944855),c(771.8754334926165,720.4287501268801,797.2887457016192,1793.1975788667767,1108.9854837165885,2383.0580707216013,1690.7482786751134),c(621.301248549928,649.5204085789588,645.0203277028804,679.6218823905084,726.3431641441474,692.0496922550104,477.74650760158517),c(762.1295315545783,739.810363483312,753.6680848425641,2037.0724495926584,2734.3687880951,2573.62247873385,1943.1746248901354))
targetgene="Zbtb7a"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_BM","AR2467_BM","AR2468_BM","AR2469_BM")

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
Sweave("BM_egfr.15m_v_Input1_DOI_summary.Rnw");
library(tools);

texi2dvi("BM_egfr.15m_v_Input1_DOI_summary.tex",pdf=TRUE);

