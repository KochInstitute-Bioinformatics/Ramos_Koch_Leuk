pdf(file='cd19_1to10_v_cd19_1to2.pdf',width=4.5,height=4.5);
gstable=read.table('cd19_1to10_v_cd19_1to2.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Jak2","Jak1","Tsc1","Irf1","Stat1","Ifngr1","Tsc2","Nprl3","Etv5","Cbfa2t3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3_vs_amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3 neg.'


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
targetmat=list(c(5102.881895553063,5135.699508500145,5312.133800657815,1573.7654676753477,1561.7392278813384,1500.4961388581357),c(9464.362005283041,8896.222798680315,9490.075324849484,2828.285952611252,2901.194099401652,2866.2171754236247),c(8588.965273851902,8397.970545665832,8530.725787223962,2539.5216466157754,2539.0238980892923,2600.260973566135),c(13138.149047057626,12812.770222947307,12686.383432935943,4461.408527630115,4612.400394491266,4512.2704247578195),c(8309.347722066495,8097.310900418235,8422.16330232292,2424.0159242175846,2278.4146527005046,2415.169833084233),c(15709.523127238663,16334.701862826827,16012.966522903653,4758.194064347687,5237.0960856438005,5090.90486393425))
targetgene="Jak2"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(13197.394746148791,13492.10158048591,13212.054412456777,3498.8608409785256,3374.5064788945238,3270.5424823015655),c(12252.785749424234,12124.897398214169,12396.121631200527,4334.673082220988,4489.760749602424,4472.736394751976),c(8148.775266585768,8112.116110222094,8269.604441961983,3056.0889051187946,3031.4987220960456,2868.0141767875266),c(8190.3026257618185,8681.54725652436,8308.458383926567,2653.4231228695467,2410.6355198462866,2515.8019094627425),c(10497.009003460827,10811.219743694839,10608.268919330214,4511.140158107113,4543.415594241292,4745.880602065074),c(4410.759242618891,4305.46889719144,4446.490828946877,1602.6418982748955,1828.0972066242905,1739.4973202570964))
targetgene="Jak1"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(3185.9789959865807,2897.835103532236,2863.1926938901056,1735.7943282616986,1891.3332735200993,1858.0994102746256),c(5416.275032801656,5196.059210008185,5223.569668238544,2841.1199217666067,2678.909743040627,2679.329033577821),c(5245.736011118676,4894.830133614286,4821.317092605211,2181.7747564102683,2301.4095861171622,2458.29786581788),c(4435.121960002174,4515.588990176976,4637.903631272398,1923.4911271587584,1987.1454960895064,2163.589642137959),c(4365.909694708757,4443.271234596588,4595.05001881146,2074.290264734174,2169.1887189713802,2050.378556212135),c(2696.5098558315344,2791.9209103200146,2879.7627573750015,1379.6516842006108,1285.8000268814453,1187.817901539195))
targetgene="Tsc1"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(5906.851569201398,5851.474459402094,5852.660699165107,2558.772600348807,2391.473075332405,2488.846889004213),c(4425.155393799922,4500.214349226815,4640.760538769794,2955.021398020378,2859.0367214711127,2920.1272163406834),c(5501.544543643146,5902.1538314229965,5760.096896249482,2284.4465096531044,2464.2903644851544,2364.8537948949784),c(4646.634642738857,4914.190792588563,4716.754278200523,2768.9288452677374,2567.7675648601144,2569.7119503798017),c(2587.985023851456,2505.4970437299744,2497.508534223439,1479.1149451546082,1657.5514504507455,1547.2181743195868),c(5061.90823449936,5263.252085271853,5209.8565122510445,2316.5314325414906,2278.4146527005046,2587.681964018821))
targetgene="Irf1"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(7741.807146660474,7724.3334995902505,7576.518683093754,3051.2761666855367,3156.054611436275,2995.6012736245657),c(10782.717234592054,11254.806606664304,10901.387628563027,3962.4879767157076,4081.600681456749,4025.2830551403886),c(2343.80415189628,2412.679766882705,2372.375985837501,1519.221098765091,1498.5031609855296,1477.1351211274105),c(9093.384263310325,9075.024178619227,9135.247413672922,3413.8357953243017,3246.118100651518,3222.0234454762126),c(6589.561354055665,6414.072431948734,6576.601059005212,2419.2031857843267,2439.3791866171086,2530.1779203739584),c(3647.7632300242603,3691.052690331294,3570.5629902453143,1350.775253601063,1347.119849325866,1403.45806520743))
targetgene="Stat1"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(3086.867032086407,2912.640313336095,3092.3166751812514,1487.1361758767048,1532.995561110516,1611.9102234200573),c(1742.4879910270715,1715.69604380873,1738.713902915105,1074.8449167609408,1161.2441375412159,1344.1570201986654),c(2793.9607253646654,2935.986990334488,2921.4736068369807,888.7523640083004,804.8226695830209,896.7036805870775),c(6091.786742065408,6223.8824290837765,6070.3570504666695,2152.8983258107205,1983.3130071867301,2230.078692602331),c(5820.474662115213,6203.952338963198,5965.222854562503,1832.0490969268576,1745.6986952146003,1802.392367993665),c(21528.89039310918,21630.41152343791,21633.64633328022,5011.664955165939,5580.103842442279,5534.764200818034))
targetgene="Ifngr1"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(3621.185720151588,3346.5468468184226,3467.142938839585,1843.2788199377926,1808.934762110409,1879.6634266414492),c(5530.336846005208,5296.848522903687,5409.2686555692735,2555.5641080599685,2389.556830881017,2398.9968208091154),c(11127.117466692096,11452.968645577494,10797.967577157297,4312.213636199118,4575.9917499148905,4801.587644346035),c(1192.6657555361653,1301.719600446982,1162.7613514401048,803.7273183540767,810.5714029371854,832.011631486607),c(12692.42205856802,12752.979952585569,12794.945917836985,4862.470063734943,4706.296372609285,4812.369652529446),c(7627.191635334576,8198.669644460038,8198.753136026566,2902.0812752545403,2839.874276957231,2930.9092245240954))
targetgene="Tsc2"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(4174.330144376579,4356.717700358644,4163.08560520521,3287.100349915176,3481.81616817226,3123.188370461605),c(3409.119339292558,3278.784540408453,3488.284054320314,2340.59512470778,2462.3741200337663,2440.32785217886),c(4478.310413545267,4480.853690252538,4477.345429918752,3798.8548699849375,3730.9279468527184,3647.912768720977),c(3880.316441410142,3866.437483392392,3983.6718143687517,2699.9462610577066,2738.3133210336596,2654.1710144831936),c(4833.2309099699105,4659.085639045147,4705.32664821094,3059.2973974076335,3230.788145040413,3205.850433201095),c(3099.0483907780485,3140.982203003304,3211.164027072918,2274.821032786588,2305.2420750199385,2237.2666980579393))
targetgene="Nprl3"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1521.5624402104838,1606.365263718695,1463.3080201661464,1158.2657162707453,1073.0968927773613,1108.7498415275088),c(2538.1521928401958,2656.396297500075,2587.786811141147,1737.3985744061179,1755.279917471541,1947.9494784697235),c(2625.0827980487275,2824.947916805546,2645.496342588543,1601.037652130476,1590.4828946521604,1446.586097941077),c(1727.5381417236933,1690.6410733714304,1600.4395800411467,1211.2058390365826,1212.9827377286958,1317.201999740136),c(1191.5583592914707,1204.9163055755967,1233.0412758760422,652.9281807786612,666.8530690830745,745.755566019313),c(2079.1364494142526,2057.9241627363926,1964.9809767088552,968.9646712292662,967.7034479510133,975.7717405987638))
targetgene="Etv5"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(2188.2149795166783,2158.1440444855916,2046.688531134376,1719.7518668175055,1657.5514504507455,1646.0532493341946),c(1844.3684455389816,1867.164728725133,1715.8586429359384,1830.4448507824382,1843.4271622353956,1804.1893693575669),c(2719.7651769701224,2648.9936925981456,2634.06871259896,1342.7540228789665,1444.8483163466615,1421.4280788464496),c(2957.8553695794776,3024.2488180113396,2863.1926938901056,2026.1628804015945,1866.4220956520533,1908.4154484638805),c(1606.8319510519739,1606.365263718695,1602.7251060390633,944.9009790629764,917.8810922149215,1049.448796518744),c(2668.27125159182,2774.8379759309464,2598.643059631251,1576.9739599641864,1615.3940725202065,1529.2481606805673))
targetgene="Cbfa2t3"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetgenelist=c("H2-T23","Brd4","Eif2ak4","Stub1","Rnf31","Keap1","Rbck1","Cbx4","Phf6","Zmiz1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3_vs_amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3 pos.'


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
targetmat=list(c(517.1540462724133,452.1283301639999,454.81967358541687,691.4300882447246,663.0205801802982,632.5444800934896),c(809.5066548718075,796.0647425305691,814.2186367578129,1175.9124238593577,1188.07155986065,1223.7579288172342),c(1190.4509630467758,1188.9722334791331,1231.3271313776047,1719.7518668175055,1670.9651616104625,1681.9932766122338),c(733.096313987875,702.6780345369973,733.6538453312504,1270.5629463800974,1402.6909384161222,1290.2469792816066),c(569.201669773063,550.0704873279898,613.6637304406253,787.6848569098836,967.7034479510133,718.8005455607836),c(270.20468370550077,315.46485505145586,310.8315357166668,421.9167359822798,452.2336905276023,436.6713314281761))
targetgene="H2-T23"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(970.6328084748828,1009.0319912476169,995.3465720927088,624.0517501791135,525.0509796803518,566.0554296291172),c(522.6910274958867,506.22428906271523,442.8206620963544,636.8857193344679,716.6754248191663,569.649432356921),c(167.21683294889596,148.05209803858938,132.56050787916672,332.07895189479814,417.7412904026156,472.61135870621524),c(97.45086953313142,104.7753309196171,93.1351844151042,319.2449827394436,343.007756798478,296.50522504382326),c(942.947902357516,1003.3376797845942,937.6370406453129,952.9222097850729,933.2110478260266,1036.8697869714304),c(240.30498509874454,227.77245852090675,283.4052237416668,206.9477526300916,208.8706452013078,195.87314866531355))
targetgene="Brd4"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1076.9428479655717,1122.348789361768,1070.7689300239588,1182.329408437035,1224.4802044370247,1318.999001104038),c(509.95597068189795,460.6697973585339,442.8206620963544,826.1867643759471,760.7490472010936,790.680600116862),c(930.7665436658746,976.0049847620854,1020.4873580697922,1381.25593034503,1506.1681387910821,1604.7222179644496),c(468.98230962819497,568.2922840096624,456.5338180838544,585.5498427130499,565.2921131595028,575.0404364486269),c(780.1606543873987,849.5912702829822,843.3590932312504,1224.0398081919373,1408.4396717702866,1326.1870065596459),c(525.4595181076234,530.1403972074105,488.53118205468775,583.9455965686305,607.449491090042,537.3034078066858))
targetgene="Eif2ak4"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(791.7883149566928,792.0787245064532,778.7929837901046,673.7833806561122,718.5916692705545,734.9735578359013),c(262.452909992638,281.8684174196221,265.6923972578126,487.6908279034717,517.3860018747991,546.2884146261956),c(808.3992586271129,808.0227966029167,733.6538453312504,755.5999340214972,837.3988252566194,871.5456614924501),c(146.17630429969714,132.67745708842818,174.27135734114592,316.036490450605,249.11177868045885,282.12921413260756),c(394.2330631113044,409.4209941913299,445.6775695937502,832.6037489536244,732.0053804302714,765.5225810222346),c(366.5481569939375,349.0612926832896,308.54600971875016,551.8606736802443,611.2819799928183,587.6194459959406))
targetgene="Stub1"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(917.4777887295385,915.0758521077429,866.7857347098962,1304.2521154129029,1071.1806483259732,1092.5768292523912),c(896.4372600803396,909.3815406447202,921.0669771604171,1155.0572239819066,1180.4065820550975,1257.9009547313713),c(796.7715980578188,869.5213604035615,785.6495617838546,1142.2232548265522,1153.5791597356633,1171.6448892640774),c(1376.4935321554813,1502.7287950916823,1481.59222814948,1524.0338371983491,1761.0286508257054,1739.4973202570964),c(640.0750294335223,650.859800223491,680.5153658796878,851.8547026866562,950.45724788852,893.1096778592737),c(387.58868564313633,428.2122220193047,390.253564144271,802.1230722096574,908.2998699579808,882.3276696758619))
targetgene="Rnf31"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(562.0035941825477,612.1384822749369,581.0949849703128,903.1905793080742,927.4623144718622,941.6287146846265),c(246.39566444456523,228.91132081351128,253.69338576875012,495.71205862556826,500.1398018123059,519.3333941676661),c(496.11351762321453,493.1273726977631,535.3844650119794,896.773594730397,914.0486033121452,993.7417542377834),c(871.5208445747094,853.0078571607958,894.7834281843755,898.3778408748162,831.650091902455,938.0347119568227),c(425.7938560851026,397.46294011898226,411.3946796250002,514.9630123586,530.7997130345162,465.4233532506074),c(649.487897513427,636.0545904196321,574.2384069765628,779.663626187787,845.063803062172,891.3126764953718))
targetgene="Keap1"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(660.0081618380265,615.5550691527505,605.0930079484378,1140.6190086821327,1111.4217818051243,995.5387556016854),c(1466.19262797575,1599.5320899630676,1422.1685522036464,1352.3794997454825,1176.5740931523212,1265.0889601869792),c(983.9215634112189,1011.309715832826,991.3469015963547,968.9646712292662,944.7085145343556,981.1627446904697),c(881.4874107769615,880.9099833296068,876.4992202010421,806.9358106429154,799.0739362288565,823.0266246670973),c(1197.6490386372914,1254.456815303894,1151.3337214505214,1034.738763150458,1092.2593372912427,1051.245797882646),c(845.4970328243845,928.172768472695,918.7814511625004,1538.4720524981228,1561.7392278813384,1635.2712411507828))
targetgene="Rbck1"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(438.5289128990914,383.2271614614256,457.67658108281273,442.7719358597309,354.50522350680683,391.7462973306271),c(694.8911435459087,670.2204591977682,703.9420073583336,903.1905793080742,846.9800475135602,787.086597389058),c(611.2827270714607,553.4870742058034,607.3785339463544,834.2079950980437,814.4038918399617,934.4407092290187),c(689.3541623224353,695.2754296350679,737.0821343281253,675.3876268005315,557.6271353539503,612.7774650905681),c(396.44785560069374,351.90844841480094,424.536454113021,771.6423954656905,691.7642469511204,796.0716042085679),c(444.06589412256477,458.9615039196271,450.24862158958354,903.1905793080742,917.8810922149215,853.5756478534306))
targetgene="Cbx4"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(648.3805012687324,644.0266264678638,615.3778749390628,798.9145799208187,779.9114917149751,763.7255796583327),c(404.19962931355644,451.5588990176976,444.53480659479186,603.1965503016623,603.6170021872657,575.0404364486269),c(702.089219136424,658.8318362717228,733.6538453312504,968.9646712292662,946.6247589857437,970.3807365070579),c(593.0106890339986,593.9166855932643,606.2357709473961,1132.5977779600362,1036.6882482009864,1223.7579288172342),c(629.5547651089229,629.2214166640049,664.5166838942712,830.9995028092051,797.1576917774684,803.2596096641757),c(450.7102715907328,438.4619826527455,454.81967358541687,543.8394429581477,521.2184907775754,610.9804637266661))
targetgene="Phf6"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(867.6449577182781,763.0377360450376,813.0758737588545,1097.3043627828113,1122.919248513453,1121.3288510748225),c(796.2178999354715,805.7450720177077,858.7863937171879,1140.6190086821327,1124.8354929648413,1090.779827888489),c(525.4595181076234,476.61386945499737,475.3894075666669,510.1502739253421,515.469757423411,469.01735597841133),c(953.4681666821153,999.9210929067806,982.7761791041671,474.8568587481172,469.4798905900956,465.4233532506074),c(169.43162543828532,149.7603914774962,122.84702238802089,332.07895189479814,283.6041788054455,307.287233227235),c(122.36728503876161,142.927217721869,146.8450453661459,277.5345829845415,243.36304532629444,256.9711950379801))
targetgene="Zmiz1"
collabel=c("amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
Sweave("cd19_1to10_v_cd19_1to2_summary.Rnw");
library(tools);

texi2dvi("cd19_1to10_v_cd19_1to2_summary.tex",pdf=TRUE);

