pdf(file='SP_egfr.15m_v_Input1_DOI.pdf',width=4.5,height=4.5);
gstable=read.table('SP_egfr.15m_v_Input1_DOI.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("H2-K1","Cnot9","Rasgrp2","B4galt1","Sec63","Gclc","Xrcc3","Ddx3x","Tmx2","Cxcr4")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2466_SP,AR2467_SP,AR2468_SP,AR2469_SP_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 neg.'


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
targetmat=list(c(877.7403953603928,847.8477458946574,852.5245428565798,233.7364843227504,174.33651119860434,175.3979388839782,356.04213800798465),c(743.4187339161891,719.5646956462483,727.2781583233672,240.2291644428268,76.91316670526663,224.3069795343183,154.54576320126807),c(662.8257370496668,618.8059780181876,718.3623140006639,111.9987320713179,189.7191445396577,477.28477600159454,240.62188447792371),c(783.7152323494502,751.0210757837892,719.2114420313976,147.7084727317381,59.821351881874044,6.7460745724607,1.9562754835603553),c(521.7879925332529,491.50593964907677,476.7853892569418,66.5499712307831,87.16825559930217,64.08770843837665,142.80811029990593),c(563.6343562908703,539.6735217346862,571.8877286991101,68.1731412608022,150.40797044585474,69.14726436772217,58.68826450681066))
targetgene="H2-K1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(598.248015201492,506.7426237781981,558.3016802073718,108.75239201127971,377.7291075969761,141.6675660216747,269.96601673132903),c(771.832931529386,758.8851708181745,722.6079541543321,387.93763717456494,531.5554410075093,150.10015923725058,262.1409147970876),c(533.6702933533171,462.5070892097812,432.6307316587923,102.25971189120331,104.26007042269477,101.1911185869105,271.9222922148894),c(526.4375885063215,544.097075191528,560.4245002842059,87.65118162103141,76.91316670526663,139.98104737855954,140.8518348163456),c(396.7655230351864,439.89781598592367,416.497299074853,108.75239201127971,80.33152966994514,111.31023044560156,66.51336644105209),c(612.1968031206977,588.8241156995939,666.5655041259115,55.1877810206494,136.73451858714066,252.97779646727625,174.10851803687163))
targetgene="Cnot9"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(677.8077685184434,661.5669947676573,703.0780094474583,202.8962537523875,126.47942969310512,165.27882702528717,58.68826450681066),c(696.4061524107178,579.9770087859106,599.9089537133204,209.38893387246392,124.77024821076586,156.8462338097113,230.84050706012192),c(588.9488232553548,543.1140633122297,506.9294343479862,58.4341210806876,76.91316670526663,236.11261003612452,254.3158128628462),c(929.9191946137181,870.948525058164,896.2546364393626,160.6938329718909,174.33651119860434,111.31023044560156,121.28907998074203),c(392.6325488369032,342.5796399354065,361.30397707716605,40.5792507504775,37.601992611463686,82.63941351264359,21.519030319163907),c(575.000035336149,572.6044196911744,574.0105487759442,228.8669742326931,271.7598556919421,312.00594897630737,191.7149973889148))
targetgene="Rasgrp2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(734.119541970052,814.9168479381692,755.2993833375775,150.9548127917763,215.35686677474655,537.9994471537409,551.6696863640202),c(521.7879925332529,490.52292776977856,561.6981923303064,110.3755620412988,126.47942969310512,79.26637622641323,95.8574986944574),c(292.40792452853583,258.5321242554144,244.97343686665667,29.2170605403438,46.14790002315998,77.57985758329805,19.562754835603553),c(765.1168484571758,700.3959639999343,742.5624628765728,35.7097406604202,109.38761486971254,121.4293423042926,144.76438578346628),c(757.3675218353949,834.0855795844832,868.2334114251523,376.5754469644312,319.6169371974413,86.01245079887393,195.62754835603553),c(692.78979998722,703.3449996378288,661.0461719261428,144.46213267169992,83.74989263462366,89.38548808510429,76.29474385885386))
targetgene="B4galt1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(714.4879145282067,645.8388046988869,760.3941515219793,321.3876659437818,536.6829854545272,421.6296607787938,197.58382383959588),c(690.206691113293,633.5511562076599,609.6739260667573,228.8669742326931,266.6323112449243,234.42609139300933,197.58382383959588),c(552.2686772455914,550.486652406966,515.4207146553227,58.4341210806876,297.39757792703097,32.043854219188326,39.125509671207105),c(483.5579811991334,434.49125064978386,442.3957040122292,99.0133718311651,133.31615562246216,99.50459994379533,72.38219289173314),c(963.499609974769,905.3539408335994,971.4024671592902,506.42904936595926,360.6372927735835,80.95289486952841,436.24943283395925),c(521.7879925332529,546.0630989501243,524.7611229933927,149.33164276175722,123.0610667284266,75.89333894018287,248.44698641216513))
targetgene="Sec63"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(635.4447829860407,643.3812750006415,606.2774139438227,292.170605403438,182.88241861030065,210.8148303893969,185.84617093823374),c(513.5220441366865,532.7924385795992,531.9787112546288,245.0986745328841,58.11217039953478,231.05305410677897,152.58948771770773),c(548.1357030473083,663.5330185262536,678.4532965561826,267.82305495315154,259.7955853155673,330.55765405057434,224.97168060944085),c(520.7547489836821,628.6360968111692,562.54732036104,116.8682421613752,222.1935927041036,254.66431511039144,158.45831416838877),c(709.3216967803528,604.5523057683644,647.8846874497714,152.57798282179542,164.0814223045688,381.15321334402955,311.04780188609647),c(667.4753330227354,577.0279731480161,592.6913654520844,287.3010953133807,244.41295197451396,256.3508337535066,314.9603528532172))
targetgene="Gclc"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(361.6352423497793,291.9545281515516,349.41618464689503,123.36092228145161,51.27544447017775,87.6989694419891,74.3384683752935),c(533.1536715785317,571.6214078118762,553.6314760383367,180.1718733321201,264.92312976258506,266.46994561219765,494.9376973407699),c(790.4313154216604,777.0708905851903,813.0400894274653,780.7447844391871,230.73950011579987,467.1656641429035,279.7473941491308),c(294.9910334024628,278.6838677810265,330.73536797075485,24.3475504502865,15.382633341053324,26.9842982898428,105.63887611225918),c(575.5166571109344,568.1808662343327,512.8733305631217,386.3144671445458,167.49978526924733,224.3069795343183,240.62188447792371),c(875.1572864864659,743.6484866890531,875.8755637017551,545.3851300864176,464.89736319627826,293.45424390204045,326.69800575457936))
targetgene="Xrcc3"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(797.6640202686559,772.1558311886996,801.5768610125612,241.8523344728459,218.77522973942507,431.7487726374848,156.50203868482842),c(708.8050750055673,702.8534936981797,748.5063590917083,360.3437466642402,690.5093188650604,462.10610821355795,436.24943283395925),c(582.2327401831446,568.1808662343327,546.4138877771007,225.62063417265492,447.8055483728857,553.1781149417774,659.2648379598397),c(625.1123474903327,610.4503770441534,574.0105487759442,178.54870330210102,119.64270376374809,118.05630501806226,215.1903031916391),c(612.7134248954832,693.5148808448473,729.8255424155682,707.7021330883276,196.5558704690147,229.3665354636638,240.62188447792371),c(656.626275752242,608.9758592252061,635.9968950195004,220.7511240825976,189.7191445396577,165.27882702528717,142.80811029990593))
targetgene="Ddx3x"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(521.2713707584676,449.7279347789052,501.41010214821756,118.4914121913943,133.31615562246216,327.18461676434396,91.9449477273367),c(423.1132335492417,384.84915074522706,440.6974479507619,394.4303172946413,170.91814823392585,193.94964395824513,54.775713539689946),c(726.3702153482709,702.3619877585306,688.6428329249864,97.390201801146,165.79060378690806,227.68001682054864,113.46397804650061),c(458.7601360094343,406.47541208978646,502.683794194318,284.05475525334253,338.4179335031732,59.02815250903113,242.57815996148406),c(790.9479371964458,894.5408101613197,805.3979371508625,384.6912971145267,143.5712445164977,143.3540846647899,295.39759801761363),c(490.2740642713436,462.01558327013214,461.0765206883694,79.5353314709359,100.84170745801624,72.52030165395253,84.11984579309528))
targetgene="Tmx2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(882.9066131082468,946.1489338244727,868.2334114251523,400.9229974147177,839.2081078285759,387.89928791649027,350.1733115573036),c(640.0943789591093,568.6723721739818,701.8043174013578,194.780403602292,249.5404964215317,327.18461676434396,150.63321223414735),c(833.3109227288485,906.3369527128975,856.770183010248,149.33164276175722,158.95387785755102,111.31023044560156,156.50203868482842),c(562.6011127412994,610.9418829838024,626.2319226660634,301.9096255835526,343.54547795019096,285.02165068646457,500.80652379145096),c(784.748475899021,748.5635460855439,785.0188644132551,467.4729686455008,323.0353001621198,537.9994471537409,223.0154051258805),c(654.5597886531004,627.161578992222,654.6777116956405,126.60726234148981,234.1578630804784,165.27882702528717,199.54009932315623))
targetgene="Cxcr4"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetgenelist=c("Sipa1","Arl15","Setd5","Pak2","Pou2af1","Cnot8","Pkn1","Nosip","Ncor1","Ubr4")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2466_SP,AR2467_SP,AR2468_SP,AR2469_SP_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 pos.'


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
targetmat=list(c(695.8895306359324,643.3812750006415,630.0529988043648,722.3106633584996,982.7793523450736,919.1526604977704,1424.1685520319386),c(952.1339309294902,896.998339859565,881.8194599168907,6830.299486320373,4025.1223909089535,4831.875912524976,6798.057305372235),c(524.3711014071799,524.9283435452139,554.9051680844372,4939.306401348122,2490.2774197682993,4379.88891617011,6207.262109337007),c(787.331584772948,732.8353560167734,704.7762655089256,4536.760233903385,5912.058747411495,5496.364257912355,5391.495232692339),c(847.7763324228397,754.9531233009818,882.6685879476242,13034.055341053374,5389.0492138156815,7467.904551713995,8069.6363696864655),c(476.84189812692324,462.5070892097812,472.5397491032736,3160.312048447188,4972.0089321249025,5297.355058024765,2263.410734479331))
targetgene="Sipa1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(806.9632122147932,808.5272707227313,724.7307742311663,2439.6245551187076,8422.846344967866,2015.3897785226343,2322.098998986142),c(1108.1537069146807,987.9269386946443,1033.3888134028464,3923.2019625561647,4642.1369060334255,4776.220797302176,2629.2342499051174),c(861.7251203420454,850.7967815325518,866.9597193790519,3923.2019625561647,5624.916258378499,3619.269008125166,3171.1225588513357),c(908.7377018475167,928.4547199971059,802.4259890432947,3022.3425958955645,2404.8183456513366,3231.3697202086755,3355.012454306009),c(1139.1510134018047,1128.006131494631,1096.224287677136,6151.814413772389,5847.109851082603,5319.279800385262,4035.796322585013),c(902.0216187753066,947.6234516434199,906.4441728081664,5095.130724229955,4795.963239443959,6337.937060826828,5031.540543717234))
targetgene="Arl15"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(788.3648283225187,702.8534936981797,732.3729265077691,2709.070780101878,2071.527956595181,2099.715710678393,1418.2997255812577),c(768.2165791058883,749.546557964842,709.0219056625938,1980.267436623302,4021.704027944275,4718.87916343626,3536.9460742771225),c(569.833817588295,643.3812750006415,650.4320715419723,2067.9186182443336,1580.992871163814,963.002145218765,1705.87222166463),c(830.7278138549215,785.4264915592246,857.6193110409818,4785.1052484963075,7934.020441018837,6084.9592643595515,4288.1558599642985),c(467.0260844060006,441.3723338048709,506.5048703326194,1990.0064568034168,1430.5849007179593,1524.6128533761182,3184.8164872362586),c(778.0323928268108,706.2940352757233,747.6572310609747,2041.9478977640279,3067.9807607989687,2813.113096716112,2793.5613905241876))
targetgene="Setd5"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(1255.3909127285194,1130.9551671325255,1262.2288176855636,9161.1716494278,6722.210770040303,6324.444911681907,7623.605559434704),c(945.4178478572801,938.7763447297366,959.0901107136524,1197.899482154096,1437.4216266473163,964.6886638618802,1529.807428144198),c(980.0315067679018,954.9960407381561,975.2235432975916,2574.347667610293,3047.470583010898,1661.2208634684475,2822.9055227775925),c(980.0315067679018,991.3674802721878,926.8232455457738,5104.86974441007,2844.077986612526,2008.6437039501734,2357.311957690228),c(548.6523248220936,584.8920681824013,506.5048703326194,3309.6436912089453,3201.296916421431,1951.3020700842576,2026.701400968528),c(824.0117307827113,753.9701114216837,798.1803488896265,4082.2726254980366,4045.6325686970245,2689.9972357687043,4063.184179354858))
targetgene="Pak2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(978.998263218331,986.9439268153461,939.9847300221453,1905.6016152424236,1430.5849007179593,2069.3583751023198,2302.536244150538),c(752.2013040875408,777.0708905851903,781.6223522903205,1873.1382146420415,2201.4257492529646,1479.0768500120084,1383.0867668771712),c(706.2219661316403,768.715289611156,751.4783071992761,3056.4291665259657,1943.3393454197367,880.3627317061214,856.8486617994356),c(545.0359723985958,624.7040492939765,534.9506593621966,1483.5774074374574,1507.4980674232258,1863.6031006422684,727.7344798844522),c(631.3118087877575,581.4515266048578,565.094704453241,5035.073433119249,1743.3651119860435,2276.8001682054864,1694.1345687632677),c(573.9667917865783,597.6712226132773,613.0704381896919,1181.667781853905,1117.804689449875,1519.5532974467728,1349.8300836566452))
targetgene="Pou2af1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(863.791607441187,937.7933328504384,926.398681530407,1652.387090559444,1069.9476079443757,1011.911185869105,1496.5507449236718),c(645.2605967069633,716.6156600083539,711.1447257394279,2572.7244975802737,1435.712445164977,2000.2111107345977,1118.9895765965232),c(664.8922241488084,609.9588711045043,585.8983412062153,1925.0796556026528,3334.613072043893,2464.0037375912707,2312.3176215683397),c(906.1545929735897,837.0346152223777,926.398681530407,2504.5513563194713,2080.0738640068776,2108.1483038939687,1919.1062493727086),c(672.6415507705894,633.0596502680108,636.4214590348671,2306.5246126571415,1789.5130120092035,2062.612300529859,2308.4050706012194),c(804.8967251156515,714.6496362497576,776.5275841059187,2754.5195409424127,2797.930086589366,3701.908421637809,2425.7815996148406))
targetgene="Cnot8"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(748.584951664043,816.3913657571164,740.8642068151055,1314.767724315471,1461.3501674000659,1398.12395514248,1003.5693230664623),c(645.7772184817487,677.2951848364278,630.9021268350984,1647.5175804693865,1004.9987116154839,426.68921670813927,991.8316701651001),c(960.3998793260566,983.5033852378026,928.946065622608,2697.7085898917444,1712.599845303937,1726.9950905499393,1762.6042106878801),c(653.0099233287443,657.1434413108155,677.1796045100821,1470.5920471973047,1037.4731597799298,1180.5630501806224,618.1830528050723),c(982.6146156418288,1177.6482313991878,1086.0347513083323,2324.3794829873514,1830.5333675853458,2175.6090496185757,2120.602624179425),c(941.2848736589968,804.1037172658896,890.3107402242271,2696.0854198617253,2917.572790353114,2679.878123910013,1750.866557786518))
targetgene="Pkn1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(549.1689465968791,641.9067571816943,579.9544449910798,1081.0312399927207,1714.309026786276,1052.3876333038693,1770.4293126221216),c(487.1743336226312,583.9090563031032,538.7717355004979,2077.657638424448,1119.5138709322143,1337.4092839903337,897.930446954203),c(521.7879925332529,562.7743008981929,500.9855381328507,686.6009226980793,680.2542299710249,497.5229997189766,958.5749869445741),c(580.6828748587884,577.0279731480161,666.5655041259115,938.1922773510398,1040.8915227446084,1354.2744704214856,2249.7168060944086),c(858.1087679185476,806.561246964135,842.3350064877761,1395.926225816426,1396.401271071174,1231.1586094740778,823.5919785789096),c(417.94701580138775,444.3213694427654,440.6974479507619,1637.778560289272,711.0194966531315,5477.812552838089,1623.7086513550948))
targetgene="Nosip"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(431.8958037205935,452.1854644771506,433.0552956741591,327.8803460638582,311.071029785745,426.68921670813927,1475.031714604508),c(945.9344696320654,875.3720785150057,850.8262867951125,2649.0134889911715,2081.783045489217,2396.5429918666637,2509.9014454079356),c(640.6110007338947,720.0562015858974,658.498787833942,999.8727384917656,681.9634114533641,1195.741717968659,717.9531024666504),c(303.77360357381457,312.1062716771637,321.8195236480516,264.5767148931133,174.33651119860434,269.84298289842803,465.59356508736454),c(771.832931529386,636.5001918455544,745.1098469687737,1485.2005774674765,3568.7709351243716,2696.743310341165,2670.316035059885),c(942.3181172085676,966.792183289734,958.2409826829187,2077.657638424448,2562.0630420265484,2727.100645917238,2509.9014454079356))
targetgene="Ncor1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
targetmat=list(c(611.1635595711269,618.3144720785385,622.410846527762,938.1922773510398,974.2334449333773,1374.5126941388676,410.8178515476746),c(757.8841436101802,735.2928857150188,783.3206083517878,808.3386749495119,1121.2230524145534,1226.0990535447322,991.8316701651001),c(830.2111920801361,911.7435180490373,780.34866024422,1628.0395401091573,1664.7427637984376,1006.8516299397595,2549.026955079143),c(726.3702153482709,647.8048284574832,794.3592727513252,1309.8982142254138,712.7286781354708,1050.701114660754,835.3296314802717),c(810.0629428635056,763.8002302146652,763.3660996295471,1425.14328635677,1131.478141308589,1151.8922332476645,1551.3264584633619),c(849.8428195219813,781.494444042032,761.243279552713,1071.2922198126062,1495.533797046851,1334.0362467041034,1445.6875823511025))
targetgene="Ubr4"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP")

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
Sweave("SP_egfr.15m_v_Input1_DOI_summary.Rnw");
library(tools);

texi2dvi("SP_egfr.15m_v_Input1_DOI_summary.tex",pdf=TRUE);

