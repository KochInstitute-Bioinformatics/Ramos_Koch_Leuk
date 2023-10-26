pdf(file='cd19_1to10_v_Input1_ACT.pdf',width=4.5,height=4.5);
gstable=read.table('cd19_1to10_v_Input1_ACT.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Ptpn2","Kdsr","Usp22","Stub1","Fitm2","Chic2","Socs1","Kat6a","Rcor1","Pcyt2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3_vs_Input2_ACT_rep1,Input2_ACT_rep2,Input2_ACT_rep3 neg.'


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
targetmat=list(c(870.8335985998186,845.8829457634338,815.3586420842145,274.56664695354397,297.706570744492,297.2037798691153),c(1506.3635469145236,1501.980729573062,1430.9381746537708,817.0838770786187,819.685424783168,782.6981531956825),c(918.6296525639824,878.4840778781979,860.836813355844,256.3724715530079,275.8747555565626,239.97820734773282),c(479.0110023661246,441.27960969627173,370.3222517832687,87.66284511167368,85.3425502800877,77.53142083542137),c(1655.529253791694,1780.2546786955131,1619.8891481513745,1070.14831674062,1105.4837326978802,1076.2099606440634),c(975.8798710485303,862.7656748942939,944.2134606871647,160.43954671381786,212.3640204644043,175.36868998488168))
targetgene="Ptpn2"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(734.2734444164937,704.9994819817745,691.9178914897916,305.9929499181062,271.90533461330267,275.05194534470917),c(1491.657068771704,1425.7173669474532,1387.6256305855522,585.52164470816,595.413141488984,566.7177665827229),c(952.2444597475701,865.6764902616835,789.9125224441361,353.9594123377012,426.71275140043855,350.73737996976337),c(1731.1625699547662,1617.2490181216924,1679.443896245175,798.8897016780827,847.4713713859873,904.533243079916),c(869.257904513088,870.333794849507,901.4423234197989,385.3857153022635,394.95738385435936,404.27098007041144),c(1077.7747553237805,1022.2783570272469,971.2838007298013,716.1889044029189,799.8383200668685,670.0929943632848))
targetgene="Kdsr"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(851.400038196807,827.8358904856179,880.3274581865423,711.226856566409,652.9697451662524,731.0105393054015),c(1560.9876085878534,1508.3845233813195,1535.4296872183481,967.5993281194169,922.8903693079252,956.2208569701969),c(1273.6860534406273,1432.7033238291883,1313.452898868728,686.4166173838598,595.413141488984,651.633132259613),c(1239.5460148947961,1186.4483437480235,1203.5473182956234,641.7581868552714,553.7342215847551,546.411918268684),c(993.7377373648112,1114.2601226367601,1044.3737188449202,491.24273581447324,506.10117026563637,494.7243043784031),c(801.5030587836691,819.685607456927,819.6898964910364,334.1112209916619,295.72186027286205,311.9716695520527))
targetgene="Usp22"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1095.1073902778178,1040.9075753785407,1115.298009756628,694.6866971113762,744.26642686123,755.0083600401748),c(1449.6385597922192,1460.647151356129,1403.3264278102815,502.8208474329962,535.8718273400856,561.1798079516213),c(1351.420295052674,1374.4870164813951,1303.7075764533788,779.0415103320433,867.3184761022867,895.3033120280802),c(660.215822340152,659.5907622504959,669.1788058539769,325.8411412641455,258.0123613118931,289.81983502764655),c(1372.4295495424162,1462.9758036500407,1328.6122892926046,858.4342757162007,758.1594001626396,786.3901256164168),c(1069.3710535278835,954.747440503807,1023.8002604125163,568.9814852531272,633.122640449953,603.6374907900664))
targetgene="Stub1"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(646.034575559576,622.3323255479083,629.6561093917275,363.8835080107209,406.86564668413905,291.66582123801373),c(978.506027859748,983.8555941777034,968.0353599246849,580.5595968716501,674.8015603541819,544.5659320583168),c(851.400038196807,813.8639767221476,717.9054179307227,441.6222574493749,460.45282941814764,526.106069954645),c(961.1733929057106,1035.6681077172393,1016.7619720014309,363.8835080107209,353.2784639501305,363.6592834423336),c(619.7730074473981,664.8302299117972,702.2046207059935,241.4863280434784,244.11938801048345,214.13440040259238),c(399.7010666673474,380.152486981089,372.4878789866797,196.82789751488994,184.57807386158504,232.59426250626413))
targetgene="Fitm2"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1066.2196653544222,1026.9356616150703,959.914257911894,693.032681165873,684.7251127123316,642.4032012077771),c(1056.2402694717946,1084.5698058893856,1020.0104128065473,463.1244647409175,504.11645979400646,437.49873185702063),c(951.193997023083,979.780452663358,980.4877163442978,745.9611914219778,672.8168498825519,686.7068702565894),c(915.4782643905211,876.7375886577641,944.2134606871647,380.4236674657537,383.0491210245797,404.27098007041144),c(779.4433415694397,760.304973962178,883.034492190806,517.7069909425256,482.28464460607705,494.7243043784031),c(1259.5048066600514,1173.0585930580312,1268.5161343979512,555.749357689101,537.8565378117155,568.5637527930901))
targetgene="Chic2"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(595.6123647841945,526.8575814975276,565.2287000902522,248.1023918254915,261.98178225515295,188.2905934574519),c(702.7595626818802,656.6799468831063,681.0897554727369,370.499571792734,373.12556866642996,378.427173125271),c(427.0130975040124,535.0078645262187,539.7825804501739,193.5198656238834,194.50162621973476,156.90882788120993),c(1044.6851795024363,915.1603515073076,937.1751722760793,871.6664032802269,867.3184761022867,828.8478084548618),c(601.9151411311171,593.2241718740117,598.254514942269,41.35039863758192,31.755367546079146,66.45550357321832),c(1046.7861049514106,958.2404189446745,949.627528695692,302.6849180270996,321.52309640405133,310.1256833416855))
targetgene="Socs1"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(640.7822619371404,672.9805129404883,565.2287000902522,226.60018453394892,184.57807386158504,251.05412460993588),c(1013.6965291300663,893.6203177886242,995.1056999673216,621.9099955092321,647.0156137513626,551.9498768997855),c(1248.4749480529367,1351.782656615756,1436.3522426622983,810.4678132966056,760.1441106342695,828.8478084548618),c(660.7410537023956,659.5907622504959,617.2037529721146,296.0688542450865,208.3945995211444,177.21467619524884),c(487.93993552426514,476.79155717842553,454.7817127162949,218.33010480643253,206.40988904951445,262.1300418721389),c(716.9408094624563,637.4685654583345,693.5421118923498,1091.6505240321626,1105.4837326978802,1052.21213990929))
targetgene="Kat6a"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(528.3827504170191,546.6511259957773,552.7763436706394,375.46161962924384,325.49251734731126,334.1235040764588),c(1324.108264216009,1304.6274476640435,1310.7458648644642,605.3698360541993,545.7953796982354,651.633132259613),c(1432.305924838182,1445.5109114457027,1383.2943761787305,643.4122028007747,643.0461928081027,701.4747599395267),c(538.3621462996466,614.7642055926952,575.5154293064542,287.79877451757017,375.1102791380599,313.81765576241986),c(689.6287786257913,673.5626760139662,625.3248549849055,279.5286947900538,279.8441764998225,267.6680005032405),c(568.3003339475295,557.1300613183801,557.1075980774613,244.79435993448496,238.1652565955936,302.7417385002168))
targetgene="Rcor1"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1187.0228786704404,1203.9132359523614,1187.3051142700415,1205.7776242718887,1222.581650524047,1270.0385127326167),c(1285.2411434099856,1314.5242199131683,1292.3380336354714,671.5304738743304,615.2602462052835,647.9411598388787),c(693.3053981614962,744.0044079047959,707.0772819136681,506.1288793240027,506.10117026563637,485.4943733265672),c(1471.6982770064487,1609.6808981664794,1470.4608711160204,1088.342492141156,928.844500722815,1037.4442502263528),c(679.1241513809201,561.7873659062035,700.0389935025826,213.36805696992272,222.287572822554,166.1387589330458),c(986.3844982934014,972.7944957816228,937.1751722760793,327.4951572096488,341.3702011203508,334.1235040764588))
targetgene="Pcyt2"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetgenelist=c("Jak2","Jak1","Nprl3","Tsc1","Irf1","Nprl2","Stat1","Tsc2","Ifngr1","Dcp2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3_vs_Input2_ACT_rep1,Input2_ACT_rep2,Input2_ACT_rep3 pos.'


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
targetmat=list(c(496.34363732016203,547.8154521427331,512.1708336066846,1622.5896425387145,1617.5390343784065,1541.3984856565917),c(879.2373003957156,845.8829457634338,933.3853246701101,2916.030111922277,3004.851654047739,2944.348005535645),c(787.8470433653366,802.220715252589,772.5875048168487,2618.307241731687,2629.7413749096795,2671.142046401303),c(1728.011181781305,1638.789051840376,1580.366451689125,4599.818344444613,4777.198105213281,4635.2713742319775),c(795.2002824367464,731.1968202882814,756.3453007912667,2499.2180936554514,2359.8207507680067,2481.005466733484),c(1222.7386113030022,1121.828242591973,1214.3754543126781,4905.811294362719,5424.213718964644,5229.678933970208))
targetgene="Jak2"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1116.6418761298037,1043.8183907459304,1069.8198384849986,3607.4087771426466,3495.075140540336,3359.6949028682598),c(1219.0619917672975,1234.7678788466917,1254.43955757578,4469.151084749854,4650.176635028965,4594.6596776039),c(1095.1073902778178,1107.856328828503,1082.8136017054642,3150.900376183742,3139.8119661185756,2946.193991746012),c(908.1250253191113,954.165277430329,897.6524758138298,2735.7423738624198,2496.765773310473,2584.380694514046),c(1840.9359246636698,1982.8474282658333,1833.7448344882034,4651.092838755214,4705.748528234603,4875.249581579711),c(736.8996012277115,774.8590507991263,794.7851836518107,1652.3619295577735,1893.4137899349691,1786.914651635426))
targetgene="Jak1"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1385.5603335985052,1444.9287483722248,1443.9319378742364,3389.0786723362144,3606.218926951613,3208.324033618151),c(1374.5304749913905,1355.2756350566235,1329.69510289431,2413.209264489281,2550.3529560444813,2506.8492736786243),c(1564.138996761315,1652.7609656038462,1725.463474317657,3916.7097589517593,3864.231288263506,3747.3520070453665),c(1348.2689068792126,1383.8016256570422,1312.3700852670224,2783.708836282015,2836.151263959194,2726.5216327123185),c(1534.2008091134321,1620.159833489082,1414.1545638273362,3154.208408074749,3346.22185516809,3293.239399295041),c(1324.108264216009,1287.1625554597056,1302.0833560508206,2345.3946107236466,2387.6066973708257,2298.2528319071334))
targetgene="Nprl3"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1156.0342282980705,1220.7959650832215,1130.9988069813573,1789.6452530345455,1958.9092354987574,1908.7497415196597),c(1073.572904425832,1132.8893409880538,1032.4627692261602,2929.2622394863033,2774.6252393386653,2752.365439657459),c(1132.3988169971103,1118.3352641511055,1080.1065677012004,2249.4616858844565,2383.637276427566,2525.3091357822964),c(765.2620947888636,708.4924604226421,761.7593687997941,1983.1651186584288,2058.1447590802545,2222.5673972820796),c(1190.1742668439017,1150.3542331923918,1097.431585328488,2138.642617535737,2246.6922538850995,2106.2702660289474),c(735.3239071409807,638.0507285318124,713.5741635239009,1422.453713132818,1331.7407264636943,1220.196885052703))
targetgene="Tsc1"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1205.405976348965,1142.2039501637007,1099.5972125318988,2638.1554330777267,2476.9186685941736,2556.690901358538),c(1427.57884257799,1346.5431889544545,1274.4716092073313,3046.697371617036,2961.1880236718803,2999.7275918466603),c(1135.0249738083282,1089.227110477209,1185.1394870666304,2355.318706396666,2552.3376665161113,2429.317852843203),c(1183.3462591347354,1149.772070118914,1106.6355009429844,2854.831521938656,2659.5120319841285,2639.7602808250613),c(1358.7735341240839,1496.1590988382827,1376.7974945684978,1525.0027017540212,1716.7745579599039,1589.394127126138),c(888.6914649160997,1021.1140308802911,887.3657465976279,2388.399025306732,2359.8207507680067,2658.220142928733))
targetgene="Irf1"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1240.5964776192834,1265.6225217410222,1245.2356419612838,2208.1112872468743,2212.9521758673905,2300.0988181175007),c(563.5732516873375,485.5240032805945,517.0434948143592,1030.4519340485415,1026.0953138326825,983.9106501257046),c(1239.5460148947961,1256.3079125653753,1127.208959375388,4636.206695245685,4408.041957490112,4417.4450014086515),c(680.1746141054073,702.0886666143848,730.8991811511883,1167.7352575253135,1192.8109934495978,1150.0494090587504),c(1227.9909249254379,1404.1773332287696,1336.1919845045427,2914.3760959767737,2768.6711079237757,2837.2808053343488),c(1221.1629172162716,1235.3500419201698,1195.967623083685,2608.3831460586675,2592.0318759487104,2750.5194534470916))
targetgene="Nprl2"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1142.9034442419816,1148.02558089848,1186.222300668336,3145.9383283472325,3268.818146774522,3077.259012682082),c(1128.196966099162,1141.039624016745,1093.100330921666,4085.4193853930938,4227.433304571787,4135.009111222474),c(677.0232259319458,808.6245090608462,768.2562504100268,1566.3531003916032,1552.0435888146183,1517.4006649218184),c(1028.403007272886,1025.1891723946367,1051.9534140568585,3519.745932030973,3362.0995389411296,3309.853275188346),c(902.8727116966757,795.2347583708538,812.1102012790981,2494.2560458189414,2526.536430384922,2599.1485841969834),c(922.3062720996874,972.2123327081449,910.1048322334426,1392.681426113759,1395.2514615558525,1441.715230296764))
targetgene="Stat1"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1271.5851279916533,1300.552306149698,1223.5793699271744,1900.464321383265,1873.5666852186696,1930.9015760440657),c(1487.980449235999,1478.1120435604669,1340.5232389113646,2634.84740118672,2474.9339581225436,2464.3915908401796),c(1386.0855649607488,1405.3416593757256,1319.9497804789607,4445.994861512808,4739.488606252313,4932.475154101093),c(618.1973133606674,646.2010115605035,576.0568361073069,828.6619886971416,839.5325294994674,854.6916154000023),c(1827.2799092453372,1666.7328793673166,1811.0057488523887,5013.322330820432,4874.448918323149,4943.551071363297),c(1436.5077757361303,1380.8908102896526,1410.364716221367,2992.1148454154277,2941.340918955581,3010.803509108863))
targetgene="Tsc2"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1433.8816189249126,1443.1822591517912,1416.320191030747,1533.2727814815375,1587.7683773039573,1655.8496306993566),c(808.3310664928354,812.1174875017139,849.4672705379367,1108.1906834871954,1202.7345458077477,1380.7976853546472),c(390.7721335092069,344.05837642545725,339.462064134663,916.3248338088154,833.5783980845775,921.1471189732206),c(724.8192798961096,706.7459712022084,729.8163675494828,2219.6893988653974,2054.175338136995,2290.8688870656647),c(719.5669662736741,711.4032757900318,738.4788763631266,1888.886209764742,1808.0712396548813,1851.5241689982772),c(1271.5851279916533,1266.2046848145,1233.3246923425236,5167.145813752237,5779.476893386404,5685.6375279309))
targetgene="Ifngr1"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(865.5812849773831,927.9679391238221,885.7415261950697,1554.7749887730802,1613.5696134351467,1797.990568897629),c(848.2486500233457,872.6624471434186,856.5055589490221,1561.3910525550932,1419.067987215412,1495.2488303974121),c(1528.9484954909965,1431.5389976822325,1446.6389718785001,2198.1871915738548,2405.4690916154955,2163.4958385503296),c(1052.5636499360896,1107.2741657550248,1045.9979392474784,1498.5384466259688,1409.144434857262,1528.4765821840215),c(244.23258344325433,281.1847644898407,259.3338576084587,400.271858811793,345.3396220636107,433.80675943628626),c(914.4278016660339,927.3857760503441,900.9009166189462,1291.7864534380592,1278.1535437296857,1390.0276164064833))
targetgene="Dcp2"
collabel=c("Input2_ACT_rep1","Input2_ACT_rep2","Input2_ACT_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
Sweave("cd19_1to10_v_Input1_ACT_summary.Rnw");
library(tools);

texi2dvi("cd19_1to10_v_Input1_ACT_summary.tex",pdf=TRUE);

