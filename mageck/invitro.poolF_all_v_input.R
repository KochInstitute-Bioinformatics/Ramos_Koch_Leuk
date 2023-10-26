pdf(file='invitro.poolF_all_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolF_all_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Dph3","Ncln","Commd4","Ilf3","Trit1","Ubxn7","Mrpl49","Slc25a19","Fitm2","Tada1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4,5,6,7,8,9,10,11,12,13_vs_0 neg.'


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
targetmat=list(c(231.79209411623071,133.4388831815999,139.64030012327143,122.3689580211026,120.07642034694052,95.14417785580517,158.25441320306135,105.9167330799772,158.48765863777552,195.46322932188724,137.98427995953767,98.62330414058404,219.65496765465497,152.32097984400724),c(764.8161080122043,411.06764961156364,447.6703739246054,411.23994089059073,432.68916987087186,479.0789190857013,462.72214295242935,446.2393508451499,453.6241253136331,406.70920386852316,336.1421556909038,297.2117260835288,406.36169016111165,420.03664138802),c(176.04462844270688,70.83011808691558,107.29714237413135,78.2358911938197,103.21840074158187,92.90549131802152,104.92955658029068,72.92627523539414,59.276880199168005,84.3769840861563,105.82253049528452,99.29421097147237,76.87923867912923,99.23942626200471),c(164.30831987985974,70.1977063182824,69.82015006163572,84.75554879330467,106.47170277068618,71.63796920907683,79.12720660153067,78.13529489506516,51.7892742792731,81.34184077370463,47.723886301794984,28.178086897309726,76.87923867912923,71.5447026540034))
targetgene="Dph3"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(145.72583132201848,115.73135365987099,123.72541297686917,127.88559137451297,96.41604195345471,142.15659514926182,86.00783326253334,114.59843251276223,61.148781679141734,88.01915606109829,31.124273675083685,45.621664500406226,65.89649029639648,50.77365994800241),c(268.9570712319133,166.95670691915817,169.92992404706928,179.03982792431816,178.9316116007365,195.88507205606945,178.89629318606936,187.52470774815637,153.4959213578456,135.36739173534426,127.6095220678431,323.3770924881735,142.7757289755257,136.16572440600646),c(261.1328655233485,162.52982453872593,135.0198490162514,165.4990006023109,186.91698930853795,88.4281182424542,151.37378654205867,154.5342499035733,120.42566187830974,134.76036307285395,130.72194943535146,124.78867054522878,109.82748382732748,120.01046896800571),c(89.97836564849463,49.96052972202081,80.6012026446824,47.14213956550674,60.33396490338883,47.01241729345667,25.802349978760002,38.19947750425408,32.4462923195446,55.23960828662031,13.48718525920293,10.063602463324901,32.94824514819824,18.463149072000878))
targetgene="Ncln"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(356.0013597396961,278.8935899672301,311.6237579956829,292.3815677307492,267.3622758463897,235.06208646728334,301.02741641886666,309.0684998071466,239.603389436637,261.62935353333353,159.77127153209625,131.49773885411204,208.67221927192222,193.8630652560092),c(250.374582674072,136.60094202476577,155.5551872696737,151.95817328030364,137.23019468221773,122.00841630920897,129.01174989380002,123.28013194554724,154.74385567782807,122.61978982304727,146.28408627289332,84.53426069192918,153.75847735825846,143.0894053080068),c(298.2978426390311,164.42705984462546,166.33623985272038,170.01260970964665,196.08538593601372,151.11134130039645,120.41096656754668,166.68862910947234,154.74385567782807,173.00316880974492,128.64699785701256,63.736148934391046,131.79298059279296,140.78151167400668),c(303.1879712068841,212.49035426074676,137.5867662979292,161.48690361801246,192.53632917699085,111.93432688918254,163.41488319881336,131.96183137833225,124.1694648382572,146.29390766017028,116.1972883869791,57.697987456396106,175.72397412372396,143.0894053080068))
targetgene="Commd4"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(885.1132707813873,676.6805924374971,697.174733703686,582.255574846312,615.7613476904684,723.0957517041193,608.935459498736,638.9730782529773,475.46297591332655,481.3737293548341,483.46371775296655,689.0213153223116,538.1546707539046,588.5128766700279),c(75.30797994493572,36.67988258072414,52.87849600256234,48.145163811581355,37.26509596974016,34.69964133564659,48.16438662701867,27.781438184912055,64.8925846390892,40.06389172436198,29.049322096744774,11.405416125101556,54.91374191366374,32.31051087600154),c(290.47363693046634,270.03982520636566,197.1392472328538,233.20313721234717,197.85991431552515,189.1690124427185,201.258329834328,217.04248581962543,207.1570971170924,125.04790447300861,86.11049050106486,112.7123475892389,120.81023221006022,124.62625623600593),c(204.4073741362541,149.88158916606244,127.31909717121806,134.40524897399794,114.45708047848764,128.72447592255992,84.28767659728267,111.12575273964822,88.6033367187564,69.8082961863883,51.87378945847281,28.178086897309726,43.93099353093099,80.77627719000384))
targetgene="Ilf3"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(506.6173196295676,314.94106077932105,353.207817958863,320.4662466208384,325.62595764034853,310.05808548303565,295.8669464231147,312.54117958026063,276.4174518761203,319.9041051324055,250.03166518983895,179.13212384718327,285.5514579510514,323.10510876001536),c(139.8576770405949,42.37158849842271,71.36030043064238,77.7343790707824,81.92406018744464,98.50220766248064,73.96673660577868,92.02601398752118,79.24382931888776,63.130980898994636,47.723886301794984,36.22896886796965,43.93099353093099,43.84997904600208),c(375.561874011108,211.85794249211358,259.25864544945614,198.09728859973578,191.35331025731654,265.28435472736265,178.89629318606936,217.04248581962543,177.83064059750401,199.71242995931956,129.68447364618203,148.27040962632023,131.79298059279296,200.78674615800955),c(123.23123990989481,104.34794182447385,67.25323277995793,109.32964282213266,91.0924568149204,97.38286439358882,108.369869910792,41.672157277368086,74.87605919894906,81.94886943619497,46.68641051262553,40.92531668418793,54.91374191366374,41.542085412001974))
targetgene="Trit1"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(241.57235125193665,151.77882447196197,138.10014975426475,163.99446423319898,123.03396764612624,110.81498362029072,228.78083647833867,98.9713735337492,87.3554023987739,86.80509873611763,79.88563576604813,62.39433527261439,120.81023221006022,120.01046896800571),c(638.6507909615977,624.1904156409436,558.0478170367502,553.6693838331855,550.6953071083824,598.8486488571266,536.688879558208,534.792685059557,500.42166231297625,464.98395546759514,406.6905093544268,371.01147748124475,516.1891739884392,572.3576212320272),c(212.23157984481884,136.60094202476577,101.64992435444023,128.8886156205876,128.35755278466056,102.97958073804794,120.41096656754668,131.96183137833225,104.20251571853744,72.84343949883997,115.15981259780963,140.2195276556603,87.86198706186198,55.38944721600264),c(464.56221394603205,244.11094269240547,218.1879689426116,241.22733118094408,225.95661365778955,203.72047493831224,182.33660651657067,289.96876105501957,155.36782283781932,170.57505415978358,175.33340836963808,283.1226826348739,186.7067225064567,207.71042706000986))
targetgene="Ubxn7"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(487.0568053581557,327.58929615198457,367.06917127992307,384.6597983696135,386.255677273656,362.6672191209515,352.6321163763867,359.4223565172997,198.421556877215,296.83701595777285,311.24273675083685,242.86827278157432,296.5342063337842,304.6419596880145),c(1.9560514271411875,1.2648235372663497,1.0267669126711134,1.0030242460746115,0.2957547299185727,3.3580298066754763,1.7201566652506668,0.0,1.8719014799737266,0.0,0.0,0.0,0.0,0.0),c(255.26471124192497,168.22153045642452,173.0102247850826,238.21825844272024,188.39576295813083,163.4241172582065,225.34052314783736,218.77882570618243,186.5661808373814,158.4344809099769,134.8718525920293,164.37217356764006,120.81023221006022,117.7025753340056),c(113.45098277418887,112.56929481670512,80.08781918834684,93.78276700797618,123.6254771059634,68.27993940240135,67.086109944776,161.47960944980133,71.1322562390016,59.48880892405264,38.386604199269875,51.659825978401166,10.982748382732748,101.54731989600482))
targetgene="Mrpl49"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(118.34111134204184,84.74317699684542,62.63278167293792,77.7343790707824,90.79670208500183,80.59271536021143,110.09002657604267,102.4440533068632,94.84300831866881,103.80190128584695,76.77320839853975,76.48337872126926,120.81023221006022,96.9315326280046),c(367.7376683025432,232.72753085700833,203.81323216521602,234.20616145842178,205.54953729340804,204.83981820720405,184.05676318182134,177.10666842881434,155.36782283781932,158.4344809099769,131.75942522452092,112.04144075835057,120.81023221006022,138.47361804000658),c(329.5946654732901,156.83811862102735,139.64030012327143,159.48085512586323,139.00472306172918,101.86023746915612,158.25441320306135,138.90719092456027,119.17772755832726,135.9744203978346,64.32349892850628,108.68690660390894,98.84473544459473,115.39468170000549),c(161.37424273914797,93.59694175770987,82.14135301368907,94.78579125405079,84.88160748663037,101.86023746915612,120.41096656754668,81.60797466817917,103.5785485585462,98.94567198592429,101.6726273386067,140.89043448654863,87.86198706186198,99.23942626200471))
targetgene="Slc25a19"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(436.1994682524848,395.2573553957343,226.40210424398052,208.6290431835192,216.49246230039523,226.10734031614874,233.94130647409068,152.7979100170163,169.09510035762662,145.0798503351896,232.3945767739582,57.02708062550778,197.68947088918947,242.32883157001152),c(480.2106253631615,459.7633557963181,282.36090098455617,319.9647344978011,307.8806738452342,315.6548018274948,287.26616309686136,256.97830321043654,365.0207885948767,227.63574843387488,377.641187257682,297.8826329144171,362.43069663018065,417.7287477540199),c(45.9672085378179,63.87358863195066,34.3966915744823,37.11189710476062,36.08207705006587,32.46095479786294,43.00391663126667,50.3538567101531,29.326456519588383,23.67411783712299,46.68641051262553,10.734509294213229,54.91374191366374,41.542085412001974),c(174.08857701556568,132.17405964433354,92.40902214040021,102.80998522264768,112.38679736905763,58.20584998237492,123.851279898048,114.59843251276223,56.78101155920304,84.3769840861563,46.68641051262553,81.17972653748754,120.81023221006022,85.39206445800406))
targetgene="Fitm2"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(82.15415993992987,39.20952965525684,60.065864391260135,55.66784565714094,42.88443583819304,66.0412528646177,43.00391663126667,32.99045784458306,44.3016683593782,40.670920386852316,91.29786944691215,2.0127204926649807,43.93099353093099,48.4657663140023),c(390.2322597146669,287.74735472809454,298.27578813095846,283.3543495160777,301.3740697870256,268.6423845340381,280.3855364358587,302.1231402609186,393.09931079448256,303.51433124516655,251.06914097900838,270.3754528479957,285.5514579510514,302.3340660540144),c(541.8262453181089,457.86612049041855,506.19608794685894,423.2762318434861,421.74624486388467,472.36285947235035,469.60276961343203,416.7215727736808,347.5497081151219,388.49834399381314,342.3670104259205,369.6696638194681,373.4134450129134,394.64981141401876),c(331.55071690043127,235.8895897001742,214.5942847482627,240.22430693486945,234.82925555534672,220.51062397168963,189.21723317757335,241.35124423142346,294.5124995158663,250.70283760850754,243.8068104548222,314.6553036866253,241.62046442012044,216.9420015960103))
targetgene="Tada1"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetgenelist=c("Cyrib","Arl15","Gm46567","Defa20_Defa32","Smyd5","4930467E23Rik_Gm15319_Gm20778_Gm21119","Pum2","Pdk2","Gm46735","Fbxw5")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4,5,6,7,8,9,10,11,12,13_vs_0 pos.'


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
targetmat=list(c(269.93509694548385,309.24935486162246,363.9888705419097,339.0221951732187,339.2306752166029,433.18584506113643,338.87086305438135,385.4674548156548,565.9382141120567,367.8593694691418,743.8701408345,856.0771162135051,779.7751351740251,454.6550458980216),c(314.9242797697312,421.18623790969446,481.0402985864166,462.8956895634332,443.92784960777766,498.1077546568623,478.2035529396854,460.1300699376059,363.7728542748942,406.1021752060328,245.8817620331611,253.60278207578753,406.36169016111165,330.0287896620157),c(614.2001481223328,727.2735339281511,667.9118766925593,759.2893542784809,829.2962626916778,606.6840517393694,808.4736326678134,800.4526877027786,861.0746807879142,887.4759045608669,1643.3616500444186,1777.903101854066,1295.9643091624641,1130.8678806600537),c(523.2437567602676,727.9059456967842,811.1458610101796,824.9874423963679,684.6721997614958,760.0340795775495,870.3992726168374,774.4075894044236,704.4589236301124,668.945586064347,455.45187144539125,531.3582100635548,483.2409288402409,576.9734085000274))
targetgene="Cyrib"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(98.78059707062997,109.40723597353924,123.2120295205336,138.41734595829638,133.38538319327628,109.6956403513989,110.09002657604267,125.01647183210424,185.31824651739893,231.88494907130723,130.72194943535146,222.7410678549245,197.68947088918947,122.31836260200582),c(354.04530831255494,507.1942384438062,532.8920276763079,479.4455896236643,449.84294420614907,536.1654257991844,497.1252762574427,595.5645810890521,544.0993635123632,654.3768981645791,641.1600377067239,1058.0200723108915,626.0166578157666,761.6048992200363),c(951.6190193041878,1269.882831415415,1274.2177386248518,1263.8105500540105,1174.7377872365707,1365.598788048027,1142.1840257264428,1302.2549149177526,2088.4180844906878,1941.277662644085,2548.0405382001845,3406.8648872509234,2767.6525924486523,2204.038420470105),c(391.2102854282375,389.5656494780357,316.75759255903847,328.49044058943525,322.6684103411628,334.68363739865583,342.3111763848827,371.57673572319874,339.43813503523575,419.45680578082016,481.3887661746277,372.35329114302135,472.25818045750816,447.7313649960213))
targetgene="Arl15"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(275.80325122690743,244.11094269240547,308.030073801334,289.3724949925254,207.32406567291946,327.96757778530485,330.27007972812805,276.07804196256353,281.4091891560502,235.52712104624922,307.092833594159,270.3754528479957,351.4479482474479,318.48932149201516),c(123.23123990989481,155.573295083761,214.08090129192715,193.58367949240002,209.39434878234948,223.8686537783651,213.29942649108267,284.75974139534856,194.05378675727633,214.28111785908757,113.08486101947072,292.5153782673105,230.6377160373877,175.39991618400833),c(59.65956852780622,56.28464740835256,105.24360854878913,93.78276700797618,85.76887167638608,89.54746151134603,130.73190655905069,112.86209262620523,129.16120211818713,97.12458599845328,172.22098100212972,108.01599977302061,153.75847735825846,163.8604480140078),c(243.52840267907786,232.09511908837516,248.990976322745,246.24245241131712,182.77642308967793,223.8686537783651,211.57926982583203,246.56026389109448,257.07446991639176,240.99037900866222,343.40448621509,383.0878004372346,252.6032128028532,353.10772600201676))
targetgene="Gm46567"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(44.989182824247315,138.4981773306653,103.70345817978246,91.27520639278964,100.55660817231472,123.1277595781008,91.16830325828533,100.7077134203062,106.07441719851117,132.3322484228926,142.1341831162155,191.87935363406146,186.7067225064567,154.62887347800736),c(243.52840267907786,235.25717793154104,272.09323185784507,241.22733118094408,233.35048190575387,226.10734031614874,235.66146313934135,220.51516559273944,167.22319887765292,176.03831212219657,129.68447364618203,151.62494378076187,164.7412257409912,198.47885252400943))
targetgene="Defa20_Defa32"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(216.14368269910122,204.90141303714864,200.21954797086713,274.82864342444356,253.4618035402168,185.81098263604304,221.90020981733602,251.7692835507655,206.53312995710115,277.4120987580822,146.28408627289332,129.48501836144706,175.72397412372396,138.47361804000658),c(421.5290825489259,524.269356196902,550.3470651917168,504.0196836524923,457.82832191395056,627.9515738483141,514.3268429099494,512.220266534316,655.1655179908043,657.4120414770307,676.4342145384854,791.6700604482256,724.8613932603613,611.5918130100291),c(363.8255654482609,387.6684141721362,385.5509757080031,384.6597983696135,365.84860090927447,307.819398945252,359.5127430373894,375.04941549631275,377.50013179470153,449.8082389053368,904.6788881557658,994.2839233765003,713.8786448776286,683.1365156640325),c(156.484114171295,266.2453545945666,296.2222543056162,217.1547492751534,267.3622758463897,280.9551604918482,208.1389564953307,246.56026389109448,298.880269635805,214.8881465215779,409.80293672193517,536.7254647106614,373.4134450129134,422.34453502202007))
targetgene="Smyd5"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(15.6484114171295,36.047470812090964,49.28481180821345,48.646675934618656,38.15236015949588,53.72847690680762,36.123289970264004,45.14483705048209,59.276880199168005,71.02235351136896,39.42407998843933,26.165366404644747,43.93099353093099,41.542085412001974))
targetgene="4930467E23Rik_Gm15319_Gm20778_Gm21119"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(938.90468502777,1212.9657722384293,1159.2198444056871,1187.0791952293027,1285.9415656859542,1105.9111496651235,1190.3484123534615,1257.1100778672705,2192.620600209225,2423.2584206614097,2017.8904099345923,4843.276412182831,2361.2909022875406,2007.8674615800956),c(410.7707996996494,402.84629661933235,498.49533610182556,440.82915614979174,388.0302056531674,519.375276765807,497.1252762574427,428.87595197957984,460.48776407353677,510.5111051543701,760.4697534612113,607.1706819539357,669.9476513466976,546.970791258026),c(1220.576090536101,1422.2940676560102,1614.0775867189902,1516.0711479417753,1526.6859158396724,1756.2495888912742,1598.0255420178694,1517.561060850821,1594.8600609376151,1738.5300893723138,1559.3261111216925,1533.6930154107151,1834.118979916369,1663.9913101140792),c(669.9476137958567,634.3090039390744,584.7437567661991,634.9143477652291,700.6429551770988,570.865067134831,540.1291928887093,663.2818366647754,617.7274883913298,702.9391911638056,925.4284039391549,729.9466320064996,680.9303997294304,662.3654729580314))
targetgene="Pum2"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(118.34111134204184,185.2966482095202,212.54075092292047,183.0519249086166,187.80425349829366,248.49420569398524,221.90020981733602,215.30614593306842,241.47529091661073,277.4120987580822,354.816719895954,426.02583761408755,527.1719223711718,288.4867042500137),c(957.4871735856113,1054.8628300801356,979.5356346882422,960.3957156164405,970.9627783226742,1168.5943727230658,990.810239184384,1165.0840638797492,1001.4672917859438,1014.9519236838369,893.2666544749018,1182.8087428561203,933.5336125322835,1045.4758162020496),c(258.19878838263674,321.2651784656528,211.0006005539138,303.41483443757,293.38869207922414,228.3460268539324,318.22898307137336,305.5958200340326,250.21083115648813,309.58461787006985,421.2151704027992,284.46449629665057,373.4134450129134,325.41300239401545),c(182.89080843770103,115.09894189123781,138.6135332106003,172.52017032483317,125.40000548547484,152.23068456928826,159.97456986831202,125.01647183210424,149.75211839789813,167.53991084733192,293.6056483349561,415.2913283198743,296.5342063337842,242.32883157001152))
targetgene="Pdk2"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(295.3637654983193,250.43506037873723,282.87428444089176,264.29688884066013,252.27878462054252,329.0869210541967,326.8297663976267,239.6149043448665,303.87200691573497,389.7124013187938,570.6116840432009,676.2740855354334,571.1029159021028,459.2708331660218),c(182.89080843770103,264.3481192886671,274.66014913952284,280.8467889008912,266.77076638655257,288.790563374091,393.91587634240267,260.4509829835505,235.23561931669832,231.88494907130723,306.05535780498957,331.42797445883343,219.65496765465497,270.0235551780128),c(231.79209411623071,307.352119555723,394.2784944657076,334.5085860658829,361.4122799604959,379.45736815432883,318.22898307137336,375.04941549631275,358.78111699496424,379.39291405645815,369.3413809443264,348.87155206192995,373.4134450129134,350.79983236801667),c(831.3218565350047,796.2064167091671,859.4039059057219,771.8271573544135,768.0750335985333,948.0837487513761,813.6341026635654,656.3364771185474,891.6490716274851,962.1404300471779,939.9530649875272,1011.0565941487085,1043.361096359611,930.0811345020442))
targetgene="Gm46735"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(532.045988182403,583.7160624484204,643.7828542447882,645.9476144720498,593.2839882166569,695.1121699818236,605.4951461682347,729.2627523539414,861.0746807879142,889.2969905483379,1920.3676857526634,2027.4804429445237,1702.3259993235758,1449.3572021520688),c(564.3208367302326,579.2891800679881,660.211124847526,699.1078995140042,655.3924814995571,753.3180199641986,703.5440760875227,727.5264124673845,880.4176627476427,773.3545160126843,1346.643574341954,1707.4578846107918,1175.154076952404,1179.3336469740561),c(363.8255654482609,350.35611982277885,382.4706749699898,385.16131049265084,345.7372792748115,441.02124794337925,383.5949363508987,328.1682385592737,323.21498887546346,370.8945127815935,467.90158091542474,420.65858296698093,307.5169547165169,348.49193873401657),c(528.1338853281206,581.1864153738877,610.9263130393125,657.9839054249452,542.7099294005809,511.5398738835642,705.2642327527734,590.3555614293812,548.4671336323019,529.3289936915704,600.6984819291151,673.5904582118801,691.9131481121631,602.3602384740286))
targetgene="Fbxw5"
collabel=c("In.vitro.input.Pool.F..41.47._1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._no.CAR_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._hEGFRv3.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.10_rep3_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep1_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep2_1439..1440..1441..1442..1443..1444..1445..1446","In.vitro_Pool.F..41.48._mCD19.CAR_ET.1.2_rep3_1439..1440..1441..1442..1443..1444..1445..1446")

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
Sweave("invitro.poolF_all_v_input_summary.Rnw");
library(tools);

texi2dvi("invitro.poolF_all_v_input_summary.tex",pdf=TRUE);

