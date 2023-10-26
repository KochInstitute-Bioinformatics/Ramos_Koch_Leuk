pdf(file='invivo.poolF_SP_mCD19_10m_v_SP_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolF_SP_mCD19_10m_v_SP_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Fitm2","Pa2g4","Aktip","Elf1","Vps16","Spinkl","Mup-ps19_Mup1_Mup10_Mup11_Mup12_Mup13_Mup14_Mup15_Mup16_Mup17_Mup18_Mup19_Mup2_Mup22_Mup3_Mup7_Mup8_Mup9","Kctd10","Tlnrd1","Capn8")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='19,20,21,22,23,24_vs_15,16,17,18 neg.'


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
targetmat=list(c(280.5395329035659,355.3957738225747,294.2130352230873,197.18592500019938,102.06241730268127,51.87349058063619,58.123066524487484,83.87946304773025,81.00850074937618,54.397011447841095),c(212.52994916936808,350.5273385647312,363.9398239175762,205.4884902633657,144.39512447936156,85.75005585778635,214.2250166188253,138.40111402875493,129.6136011990019,135.21542845606214),c(34.004791867098895,73.72201961877283,57.822215014942,85.10129394745448,32.474405505398586,23.290138628040737,31.552521827578918,0.0,25.92272023980038,94.80621995195162),c(30.361421309909726,72.33103811653183,82.48168906543198,110.00898973695334,23.77590403073825,49.756205250814304,0.0,4.193973152386513,38.88408035970057,20.204604252055265))
targetgene="Fitm2"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(178.52515730226918,223.25253110967998,238.94180028233387,166.0513052633258,79.4463134685644,125.97847712440218,29.891862784022134,137.70211850335716,112.33178770580164,97.91462060611397),c(174.88178674508,229.51194786976447,201.52742586090082,197.18592500019938,135.11672290639055,99.51241050162861,68.08702078582819,153.77901558750546,46.44487376297568,115.0108242040069),c(167.5950456307017,139.09815022409967,139.45357738897778,166.0513052633258,72.48751228883613,52.93213324554713,74.72965696005534,62.90959728579769,51.84544047960076,60.61381275616579),c(217.38777657895363,166.9177802689196,291.6620551488987,290.58978421082014,227.32083853779008,142.91675976297725,308.88258210156204,148.18705138432344,211.7022152917031,124.33602616649394))
targetgene="Pa2g4"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(183.38298471185473,172.4817062778836,111.39279657290298,201.33720763178255,145.55492467598293,98.45376783671766,127.87074635387246,186.6318052811998,166.33745487205243,376.11647915364415),c(70.43849743899057,167.61327102004012,133.50129054920433,135.954506184348,63.2091107158651,177.85196770503836,29.891862784022134,36.34776732068311,50.765327136275744,85.48101798946458),c(197.95646694061142,119.62440919272572,158.16076459969432,97.55514184220391,73.06741238714682,49.756205250814304,142.81667774588354,139.10010955415268,64.80680059950095,66.83061406449049),c(162.73721822111614,129.3612797084127,143.70521084595882,184.73207710544995,76.54681297701094,60.34263189992373,28.231203740465347,91.56841382710553,83.16872743602622,94.80621995195162))
targetgene="Aktip"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(34.004791867098895,39.64297281386841,67.17580862030027,65.38270144743453,8.118601376349647,94.2191971770739,136.17404157165637,68.5015614889797,31.32328695642546,37.30080784994818),c(180.95407100706197,145.35756698418416,203.22807924369323,239.73657197392663,38.27340648850547,74.10498654376599,9.96395426134071,115.3342616906291,206.30164857507802,87.03521831654575),c(241.6769136268814,227.425475616403,210.8810194662591,250.1147785528845,124.67852113679814,102.68833849636144,91.33624739562319,186.6318052811998,213.86244197835313,94.80621995195162),c(365.5515125713131,480.5841090242644,464.27837350232846,330.02696921086005,334.02245662695685,335.58972477676883,425.128715150537,334.81885666552324,466.60896431640685,433.6218912556476))
targetgene="Elf1"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(199.1709237930078,290.7151339683683,123.29737025244987,136.99232684224378,150.77402556077914,123.86119179458028,101.30020165696389,39.84274494767187,79.92838740605117,97.91462060611397),c(85.01197966774723,221.861549607439,136.05227062339296,113.1224517106407,56.83020963444752,113.27476514547087,107.94283783119104,90.86941830170778,76.68804737607613,35.74660752286701),c(155.4504771067378,299.7565137329348,112.24312326429919,241.8122132897182,60.88951032262234,136.5649037735116,106.28217878763425,291.4811340908626,52.92555382292578,74.60161569989636),c(123.87459894443168,41.03395431610941,94.38626274497886,117.27373434222385,128.73782182497297,15.879639973664139,0.0,11.183928406364034,31.32328695642546,18.65040392497409))
targetgene="Vps16"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(205.24320805498974,281.6737542038019,237.24114689954146,215.86669684232353,301.5480511215583,185.26246635941496,373.64828480027666,327.129905886148,300.2715094443544,271.9850572392055),c(174.88178674508,209.34271608727002,122.44704356105366,138.03014750013958,111.34081887565229,34.935207942061105,31.552521827578918,153.77901558750546,43.20453373300063,49.73441046659757),c(335.1900912614034,522.3135540914943,448.97249305719674,365.3128715793168,249.93694237190695,208.55260498745568,385.27289810517414,260.0263354479638,231.14425547155338,359.02027555575125),c(359.4792283093311,337.3130142934417,307.81826228542656,387.10710539512826,265.59424502629554,562.1392550677106,280.6513783610967,407.5143913068895,437.4459040466314,341.9240719578583))
targetgene="Spinkl"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(17.002395933549447,0.0,0.0,23.869875131603084,5.219100884796201,6.351855989465656,0.0,0.0,0.0,0.0),c(2.428913704792778,15.300796524650965,27.21045412467859,10.378206578957862,0.579900098310689,2.117285329821885,0.0,11.882923931761786,0.0,3.1084006541623483))
targetgene="Mup-ps19_Mup1_Mup10_Mup11_Mup12_Mup13_Mup14_Mup15_Mup16_Mup17_Mup18_Mup19_Mup2_Mup22_Mup3_Mup7_Mup8_Mup9"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(115.37340097765696,75.11300112101382,39.96535449562168,143.2192507896185,37.693506390194784,25.407423957862623,33.213180871135705,36.34776732068311,19.442040179850284,101.02302126027632),c(245.32028418407057,223.94802186080048,298.4646686800683,323.8000452634853,187.88763185266325,230.7841009505855,68.08702078582819,90.86941830170778,210.6221019483781,184.94983892265972),c(306.04312680389,255.9405964123434,205.77905931788186,306.1570940792569,175.70972978813876,185.26246635941496,83.03295217783926,267.71528622733905,112.33178770580164,239.34685037050082),c(179.73961415466556,199.60584557158305,223.63591983720218,211.7154142107404,294.0093498435193,191.61432234888062,229.17094801083636,406.116400256094,156.6164347821273,195.82924121222794))
targetgene="Kctd10"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(352.19248719495283,424.249358183504,430.2653058464802,188.8833597370331,254.57614315839248,231.84274361549643,333.79246775491384,330.62488351313675,422.3243172400812,212.92544481012087),c(134.80471061599917,209.34271608727002,174.31697173622223,308.2327353950485,89.3046151398461,129.154405119135,195.95776713970065,160.06997531608525,200.90108185845295,88.58941864362693),c(178.52515730226918,390.1703113785996,215.98297961463632,243.88785460550977,73.06741238714682,56.10806124027996,44.8377941760332,104.84932880966282,74.5278206894261,170.96203597892915),c(115.37340097765696,108.49655717479776,114.7941033384878,110.00898973695334,34.21410580033065,58.22534657010184,36.53449895824927,47.53169572704714,60.48634722620089,76.15581602697753))
targetgene="Tlnrd1"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(415.34424351956505,523.0090448426148,409.00713856157506,501.26737776366474,276.03244679588795,314.41687147855,574.5880290706476,255.13336677017952,282.98969595115415,295.2980621454231),c(234.39017251250309,177.3501415357271,212.5816728490515,209.63977289494883,584.5392990971745,200.08346366816815,302.23994592733493,500.4807961847905,439.6061307332815,161.63683401644212),c(78.93969540576529,83.4588901344598,129.24965709222332,74.7230873684966,67.26841140403992,148.20997308753198,59.78372556804427,46.832700201649395,36.72385367305054,118.11922485816923),c(134.80471061599917,123.79735369944872,142.0045574631664,99.63078315799548,40.59300688174823,53.990775910458076,34.87383991469249,66.40457491278644,111.25167436247663,57.50541210200345))
targetgene="Capn8"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetgenelist=c("Arl15","Slc23a2","Calm1","Rpn2","Setd5","Adgrb3","Gal3st2_Gal3st2b","Gm10354","4930555G01Rik_Gm16434_Gm16440_Gm21560_Gm3278_Gm6676_Gm9602_Gm9603_LOC108168187","Zfp931")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='19,20,21,22,23,24_vs_15,16,17,18 pos.'


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
targetmat=list(c(121.4456852396389,68.15809360980884,68.87646200309268,30.096799078977803,309.08675239959723,362.0557913995424,259.0628107948585,488.59787225302875,406.12261709020595,66.83061406449049),c(632.7320200985187,992.4653018489512,346.08296339825586,300.967990789778,440.72407471612365,790.8060706884742,539.7141891559552,653.5608162468982,778.7617205373364,446.055493872297),c(1363.8350452411448,1110.0032387883155,792.504476381264,790.8193413165891,2756.2651672707048,1261.9020565738435,1482.9685258962093,2837.222837589476,1497.037093848472,3027.582237154127),c(195.52755323581863,468.06527550409544,382.6470111282927,130.76540289486906,288.2103488604124,268.89523688737944,815.3835903863816,610.2230936722376,572.4600719622584,1073.9524260130913))
targetgene="Arl15"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(239.24799992208864,284.45571720828383,335.87904310150134,320.68658328979797,646.5886096164182,599.1917483395936,370.3269667131631,383.7485434433659,328.3544563708048,595.2587252720897),c(335.1900912614034,404.0801264010096,338.43002317569,180.5807944738668,565.4025958529218,378.99407403811745,209.24303948815492,492.09284988001747,387.7606902536807,833.0513753155094),c(514.9297054160689,621.7687315017256,420.91171224112196,500.229557105769,790.9837340957798,1009.9451023250392,498.1977130670356,732.5473106168442,977.5025757091394,1162.5418446567182),c(253.8214821508453,308.1024027463808,215.98297961463632,227.2827240791772,426.8064723566671,335.58972477676883,712.4227296858609,359.98269557984236,476.329984406332,391.6584824244559))
targetgene="Slc23a2"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(1102.7268219759212,994.5517741023127,777.1985959361323,991.1187282904759,2017.472442022887,1446.1058802683476,1378.3470061521316,1776.8466255610858,1866.4358572656274,1445.406304185492),c(527.0742739400329,577.9528141811342,385.1979912024813,615.4276501322013,1243.3058107781171,869.1456278918839,830.3295217783926,692.0055701437745,558.4185984990332,756.8955592885318),c(64.36621317700862,130.0567704595332,70.5771153858851,121.42501697380699,234.85953981582904,449.92313258715063,254.08083366418813,236.26048758444023,141.49484797557707,357.46607522867004),c(156.66493395913417,122.40637219720772,93.53593605358266,161.90002263174264,317.2053537759469,339.82429543641257,49.81977130670356,67.80256596358196,143.65507466222712,97.91462060611397))
targetgene="Calm1"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(75.29632484857612,132.8387334640152,54.42090824935718,18.68077184212415,147.294624970915,200.08346366816815,119.56745113608854,116.03325721602685,89.64940749597632,206.70864350179616),c(2.428913704792778,36.86100980938642,24.659474050489973,11.41602723685365,56.83020963444752,158.7963997366414,18.267249479124636,36.34776732068311,65.88691394282597,63.72221341032814),c(72.86741114378334,11.823342769048473,9.353593605358267,15.567309868436794,52.7709089462727,10.586426649109427,28.231203740465347,24.464843388921324,19.442040179850284,35.74660752286701),c(78.93969540576529,62.594167600844855,76.52940222565854,55.00449486847667,219.20223716144045,131.2716904489569,190.9757900090303,82.48147199693474,100.45054092922648,183.39563859557856))
targetgene="Rpn2"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(608.4428830505909,413.81699691669655,387.7489712766699,514.75904631631,712.117320725526,571.667039051909,528.0895758510577,916.383133796453,537.8964449758579,842.3765772779964),c(576.8670048882848,445.114080717119,351.18492354663306,459.7545514478333,1134.8644923940183,768.5745747253443,287.2940145353238,1114.1988674840168,499.0123646161573,503.56090597430045),c(712.8861723566804,909.701902465612,708.3221339330396,425.50646973727237,1632.4187767445896,1221.6736353072279,1009.6806984825254,1140.0617019237336,911.6156617663133,786.4253655030741),c(445.70566482947476,193.34642881149855,220.23461307161736,368.4263335530041,315.46565348101484,607.660889658881,461.6632141087863,180.34084555262004,335.91524977407994,256.4430539683937))
targetgene="Setd5"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(272.03833493679116,135.6206964684972,225.3365732199946,226.2449034212814,305.60735180973313,341.94158076623444,371.98762575671986,297.07309829404466,450.4072641665316,228.4674480809326),c(299.9708425419081,147.44403923754567,240.64245366512628,171.24040855280472,451.16227648571606,845.8554892638432,254.08083366418813,199.9127202637571,478.490211092982,354.3576745745077),c(85.01197966774723,105.71459417031576,67.17580862030027,60.193598157955606,240.07864070062524,174.67603971030553,58.123066524487484,120.22723036841336,180.37892833527764,198.9376418663903),c(274.4672486415839,289.3241524661273,229.5882066769756,343.51863776350524,342.1410580033065,138.68218910333349,149.45931392011067,188.02979633199533,162.01700149875236,200.49184219347146))
targetgene="Adgrb3"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(102.01437560129668,80.67692712997781,81.63136237403577,145.29489210541007,156.57302654388604,31.759279947328277,149.45931392011067,89.47142725091227,169.57779490202748,164.74523467060448),c(136.01916746839558,324.09869002215225,134.35161724060055,140.10578881593113,233.119839520897,236.0773142751402,137.83470061521317,265.6182996511458,280.8294692645041,368.3454775182383),c(126.30351264922446,98.75968665911077,161.56207136527914,78.87437000007975,178.0293301813815,195.84889300852439,225.84962992372277,125.81919457159538,180.37892833527764,172.51623630601034))
targetgene="Gal3st2_Gal3st2b"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(114.15894412526056,63.28965835196536,68.02613531169648,60.193598157955606,155.41322634726464,207.49396232254475,370.3269667131631,202.7087023653481,173.89824827532755,149.2032313997927),c(178.52515730226918,102.23714041471327,137.75292400618537,120.38719631591121,490.0155830725322,143.97540242788818,347.0777401033681,80.3844854207415,250.58629565140367,102.5772215873575),c(75.29632484857612,79.28594562773682,79.93070899124336,57.08013618426824,129.31772192328364,70.92905854903316,132.85272348454282,99.2573646064808,187.93972173855275,222.2506467726079))
targetgene="Gm10354"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(3.643370557189167,0.0,0.0,0.0,0.0,14.820997308753196,9.96395426134071,13.28091498255729,71.28748065945105,0.0))
targetgene="4930555G01Rik_Gm16434_Gm16440_Gm21560_Gm3278_Gm6676_Gm9602_Gm9603_LOC108168187"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
targetmat=list(c(60.72284261981945,97.36870515686978,80.78103568263957,167.08912592122158,111.92071897396298,184.20382369450402,205.92172140104137,181.73883660341556,220.34312203830322,60.61381275616579),c(208.8865786121789,166.2222895177991,188.77252548995773,119.34937565801542,202.96503440874116,172.55875438048363,164.40524531212174,212.49463972091664,141.49484797557707,132.10702780189982),c(155.4504771067378,182.91406754469108,147.10651761154364,209.63977289494883,310.24655259621863,293.2440181803311,338.77444488558416,327.82890141154576,248.42606896475365,281.31025920169253))
targetgene="Zfp931"
collabel=c("In.vivo_AR1910_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1912_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1913_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1914_Pool.F..41.48._SP_hEGFRv3_15m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1925_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1926_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1927_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1928_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1929_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446","In.vivo_AR1924_Pool.F..41.48._SP_mCD19_10m_1439..1440..1441..1442..1443..1444..1445..1446")

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
Sweave("invivo.poolF_SP_mCD19_10m_v_SP_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolF_SP_mCD19_10m_v_SP_hEGFRv3_15m_summary.tex",pdf=TRUE);

