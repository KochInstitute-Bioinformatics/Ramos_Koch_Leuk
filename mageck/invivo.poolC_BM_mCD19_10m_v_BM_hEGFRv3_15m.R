pdf(file='invivo.poolC_BM_mCD19_10m_v_BM_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolC_BM_mCD19_10m_v_BM_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Med15","4930555G01Rik_Gm21560_Gm3159_Gm3594_Gm3642","Srpx","Vars2","Mrps18b","Gm20738_Gm20747_Gm20781_Gm20798_Gm20800_Gm20801_Gm20803_Gm20804_Gm20805_Gm20806_Gm20809_Gm20810_Gm20813_Gm20815_Gm20816_Gm20818_Gm20823_Gm20840_Gm20841_Gm20842_Gm20844_Gm20846_Gm20847_Gm20848_Gm20852_Gm20853_Gm20854_Gm20859_Gm20860_Gm20861_Gm20862_Gm20863_Gm20865_Gm20866_Gm20867_Gm20868_Gm20879_Gm20880_Gm20881_Gm20882_Gm20884_Gm20889_Gm20891_Gm20892_Gm20893_Gm20895_Gm20898_Gm20902_Gm20907_Gm20909_Gm20912_Gm20913_Gm20917_Gm20918_Gm20919_Gm20921_Gm20924_Gm20925_Gm20926_Gm20927_Gm20930_Gm20932_Gm20933_Gm20934_Gm20935_Gm20936_Gm21076_Gm21111_Gm21114_Gm21118_Gm21127_Gm21151_Gm21155_Gm21160_Gm21163_Gm21180_Gm21184_Gm21201_Gm21245_Gm21247_Gm21249_Gm21268_Gm21275_Gm21281_Gm21282_Gm21285_Gm21287_Gm21292_Gm21302_Gm21308_Gm21316_Gm21330_Gm21333_Gm21340_Gm21344_Gm21380_Gm21387_Gm21394_Gm21396_Gm21412_Gm21427_Gm21435_Gm21443_Gm21462_Gm21469_Gm21470_Gm21476_Gm21506_Gm21524_Gm21530_Gm21573_Gm21617_Gm21634_Gm21642_Gm21654_Gm21661_Gm21672_Gm21683_Gm21724_Gm21745_Gm21751_Gm21753_Gm21763_Gm21783_Gm21805_Gm21810_Gm21843_Gm21908_Gm21918_Gm21943_Gm28093_Gm28371_Gm28459_Gm29317_Gm31186","Ptar1","Actl6a","Plekho2","Dmxl1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='5,6,7,8_vs_1,2,3,4 neg.'


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
targetmat=list(c(265.1228516417483,228.9739209833356,217.90042279614707,164.90191249947168,214.80720781389655,0.0,34.891501952635096,2.3642566059414314),c(0.0,8.428487888957138,89.72370350429584,68.82460773787473,50.12168182324253,0.0,0.0,0.0),c(1.3666126373285994,3.277745290149998,2.5635343858370243,5.54292142855367,0.0,143.51864128069442,0.0,0.0),c(379.46277563157446,40.73769146329283,236.27241922797907,186.14977797559408,18.616624677204367,0.0,24.553279151854326,137.12688314460303))
targetgene="Med15"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(392.67336445908427,19.66647174089999,124.7586734440685,461.91011904613913,0.0,0.0,0.0,16.54979624159002))
targetgene="4930555G01Rik_Gm21560_Gm3159_Gm3594_Gm3642"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(478.3144230650098,405.50391732427124,137.14908964228079,290.5414648800215,164.68552599065404,4.484957540021701,10.33822280078077,9.457026423765726),c(334.3645585997307,54.7851712782214,55.9705007574417,33.71943869036816,256.3366013245832,672.7436310032551,45.229724753415866,309.71761537832754),c(95.66288461300196,263.1561218663284,30.335156899071453,175.98775535657902,167.54962209483932,524.740032182539,12.922778500975962,4.728513211882863),c(225.4910851592189,85.68962687106423,160.64815484578685,179.22112618990198,4.296144156277931,58.30444802028211,0.0,0.0))
targetgene="Srpx"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(106.14024816585456,248.6403927242356,166.20247934843374,143.65404702334928,78.7626428650954,0.0,2.5845557001951924,30.735335877238608),c(183.12609340203232,54.7851712782214,628.9204359920166,62.81977619027492,20.048672729297014,62.78940556030381,16.799612051268753,44.92087551288719),c(3.6443003662095985,117.06233179107137,23.92632093447889,243.88854285636148,0.0,161.45847144078124,10.33822280078077,0.0),c(9.110750915523996,4.682493271642855,183.29270858734722,29.100337499906768,0.0,0.0,0.0,85.11323781389153))
targetgene="Vars2"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(18.677039376824194,28.563208957021413,138.85811256617214,405.09517440346406,10.024336364648507,31.394702780151906,0.0,0.0),c(20.04365201415279,27.626710302692842,132.02202087060675,20.78595535707626,50.12168182324253,71.75932064034721,126.64322930956443,26.006822665355745),c(256.4676382720005,35.11869953732141,86.73291338748598,207.8595535707626,27.208912989760233,4.484957540021701,112.42817295849088,7.092769817824294),c(5.921988095090597,29.0314582841857,13.244927660157959,26.79078690467607,42.961441562779314,0.0,0.0,0.0))
targetgene="Mrps18b"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(44.64267948606758,1.8729973086571419,242.68125519257163,295.62247618952904,2.8640961041852875,0.0,0.0,0.0))
targetgene="Gm20738_Gm20747_Gm20781_Gm20798_Gm20800_Gm20801_Gm20803_Gm20804_Gm20805_Gm20806_Gm20809_Gm20810_Gm20813_Gm20815_Gm20816_Gm20818_Gm20823_Gm20840_Gm20841_Gm20842_Gm20844_Gm20846_Gm20847_Gm20848_Gm20852_Gm20853_Gm20854_Gm20859_Gm20860_Gm20861_Gm20862_Gm20863_Gm20865_Gm20866_Gm20867_Gm20868_Gm20879_Gm20880_Gm20881_Gm20882_Gm20884_Gm20889_Gm20891_Gm20892_Gm20893_Gm20895_Gm20898_Gm20902_Gm20907_Gm20909_Gm20912_Gm20913_Gm20917_Gm20918_Gm20919_Gm20921_Gm20924_Gm20925_Gm20926_Gm20927_Gm20930_Gm20932_Gm20933_Gm20934_Gm20935_Gm20936_Gm21076_Gm21111_Gm21114_Gm21118_Gm21127_Gm21151_Gm21155_Gm21160_Gm21163_Gm21180_Gm21184_Gm21201_Gm21245_Gm21247_Gm21249_Gm21268_Gm21275_Gm21281_Gm21282_Gm21285_Gm21287_Gm21292_Gm21302_Gm21308_Gm21316_Gm21330_Gm21333_Gm21340_Gm21344_Gm21380_Gm21387_Gm21394_Gm21396_Gm21412_Gm21427_Gm21435_Gm21443_Gm21462_Gm21469_Gm21470_Gm21476_Gm21506_Gm21524_Gm21530_Gm21573_Gm21617_Gm21634_Gm21642_Gm21654_Gm21661_Gm21672_Gm21683_Gm21724_Gm21745_Gm21751_Gm21753_Gm21763_Gm21783_Gm21805_Gm21810_Gm21843_Gm21908_Gm21918_Gm21943_Gm28093_Gm28371_Gm28459_Gm29317_Gm31186"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(251.0011877226861,28.09495962985713,130.74025367768823,3.6952809523691132,77.33059481300276,0.0,0.0,0.0),c(27.33225274657199,37.45994617314284,614.3937411389402,389.85214047494145,1225.833132591303,148.00359882071612,23.261001301756732,0.0),c(2.7332252746571988,31.372704920007127,5.981580233619723,2.771460714276835,0.0,125.57881112060763,0.0,0.0),c(38.265153845200786,42.61068877194998,8.972370350429586,2.771460714276835,4.296144156277931,0.0,43.937446903318275,0.0))
targetgene="Ptar1"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(8.199675823971596,0.46824932716428547,37.17124859463685,0.0,0.0,0.0,0.0,0.0),c(256.4676382720005,132.98280891465708,11.963160467239446,71.59606845215157,44.39348961487195,35.879660320173606,232.61001301756733,435.02321549322335),c(25.51010256346719,29.0314582841857,86.30565765651315,254.51247559442268,267.7929857413244,0.0,0.0,0.0),c(90.19643406368756,63.68190849434282,23.92632093447889,52.65775357125986,20.048672729297014,0.0,3.8768335502927886,7.092769817824294))
targetgene="Actl6a"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(371.2630998076028,244.426148779757,356.3312796313464,300.7034874990366,204.78287144924806,228.73283454110674,18.091889901366347,330.9959248318004),c(180.84840567315132,196.1964680818356,183.29270858734722,381.53775833211097,12.888432468833793,0.0,0.0,0.0),c(201.3475952330803,22.4759677038857,45.28910748312076,43.41955119033708,24.344816885574943,4.484957540021701,0.0,54.37790193665292))
targetgene="Plekho2"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(17.310426739495593,165.76026181615705,43.152828828256574,51.27202321412145,199.05467924087748,0.0,60.73705895458702,7.092769817824294),c(143.03878937372673,5.618991925971425,108.09569993612786,3.233370833322974,0.0,0.0,6.461389250487981,0.0),c(331.1757957792973,383.0279496203855,18.371996431832006,528.8870863078293,60.14601818789104,4.484957540021701,232.61001301756733,373.5525437387462),c(5.4664505493143976,0.9364986543285709,38.02576005658253,50.34820297602917,0.0,0.0,10.33822280078077,0.0))
targetgene="Dmxl1"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetgenelist=c("Plagl2","Fibp","Zfp42","Frk","Dennd6a","Ocm","Ccdc125","Slc25a26","Vmn2r51","Sirt7")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='5,6,7,8_vs_1,2,3,4 pos.'


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
targetmat=list(c(130.28373809199314,40.26944213612855,201.2374492882064,216.63584583263926,1313.1880637689544,215.27796192104165,664.2308149501645,113.4843170851887),c(111.60669871516896,116.59408246390709,518.6884574010246,366.29472440358836,1008.1618286732212,340.85677304164926,220.97951236668897,347.54572107339044),c(29.609940475452987,55.253420605385685,179.02015127761885,107.62505773775042,2523.2686677872384,834.2021024440363,81.41350455614857,3250.8528331694683),c(891.486977084023,199.9424626991499,38.8802715185282,133.03011428528808,1459.256965082404,3000.436594274518,1210.8643455414476,352.2742342852733))
targetgene="Plagl2"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(154.88276556390792,14.51572914209285,5.127068771674049,16.62876428566101,63.010114292076324,0.0,112.42817295849088,186.7762718693731),c(2250.355476134427,421.4243944478569,1047.2037966144244,921.9725976160937,48.689633771149886,7848.675695037977,439.3744690331827,226.9686341703774),c(58.308805859353576,479.019061689064,35.88948140171834,35.10516904750658,2125.1593093054835,5346.0693877058675,15051.160120086703,92.20600763171582),c(803.5682307492165,2327.199156006499,2276.8457903542503,567.2256261886589,292.1378026268993,9848.966757887654,158.95017556200435,333.3601814377418))
targetgene="Fibp"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(415.90577929367043,294.52882678633557,411.0200131958695,120.55854107104231,1661.1757404274667,816.2622722839495,109.84361725829568,3721.339897751813),c(88.37428388058277,20.60297039522856,64.94287110787128,168.13528333279464,482.6001935552209,1318.57751676638,21.968723451659137,1052.094189643937),c(113.42884889827376,37.928195500307126,117.06807028655744,653.1409083312408,2139.4797898264096,421.58600876203985,673.2767599008477,182.04775865749022),c(224.1244725218903,159.20477123585707,85.45114619456747,166.2876428566101,157.5252857301908,49.33453294023871,827.0578240624616,163.13370580995877))
targetgene="Zfp42"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(334.3645585997307,296.87007342215696,465.7087467603927,329.3419148798972,832.019918265826,224.24787700108504,1378.8604660541353,716.3697516002537),c(801.2905430203355,378.81370567590693,783.5870106041838,1115.512937496426,1705.5692300423386,3260.5641315957764,2485.0503057376777,650.1705666338936),c(390.85121427597943,471.0588231272712,306.7696148384972,183.37831726131725,1376.1981780610306,179.39830160086802,191.25712181444425,3662.233482603277),c(79.71907051083497,77.72938830927139,158.08462045994983,93.76775416636625,345.12358055432713,565.1046500427343,246.82506936864087,248.24694362385029))
targetgene="Frk"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(75.61923259884917,95.99111206867852,1060.8759800055552,258.2077565467918,1145.638441674115,112.12393850054252,324.36174037449666,768.3833969309652),c(403.15072801193685,714.5484732526996,1754.3120313744703,286.3842738086063,1003.8656845169432,892.5065504643185,2899.871495619006,1255.4202577549001),c(91.56304670101616,127.83206631584993,251.22636981202837,157.51135059473344,253.47250522039795,22.424787700108503,1115.2357846342256,763.6548837190824),c(225.4910851592189,263.1561218663284,142.27615841395485,543.6682101173058,935.1273780164963,130.06376866062934,519.4956957392337,680.9059025111322))
targetgene="Dennd6a"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(154.42722801813173,41.20594079045712,52.552454909659,91.9201136901817,54.417825979520465,71.75932064034721,737.8906524057275,879.5034574102125),c(225.4910851592189,67.4279031116571,431.9555440135386,119.1728107139039,1363.3097455921968,161.45847144078124,6.461389250487981,1922.1406206303836),c(210.9138836943805,171.84750306929277,397.7750855357116,182.45449702322497,22.9127688334823,94.18410834045572,12.922778500975962,37.8281056950629),c(313.8653690398017,122.68132371704279,8.11785888848391,12.471573214245756,953.7440026937007,242.18770716117183,50.39883615380625,470.48706458234483))
targetgene="Ocm"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(183.58163094780852,199.00596404482133,621.2298328345055,122.40618154722688,65.87421039626162,1435.1864128069442,90.45944950683173,1572.2306429510518),c(59.219880950905974,359.6154832621712,70.92445134149101,66.0531470235979,34.36915325022345,130.06376866062934,1104.8975618334448,1539.1310504678718),c(18.677039376824194,23.412466358214274,17.090229238913494,226.3359583326082,267.7929857413244,1668.4042048880726,24.553279151854326,52.01364533071149),c(87.91874633480657,87.0943748525571,244.8175338474358,167.6733732137485,204.78287144924806,816.2622722839495,323.06946252439906,252.97545683573316))
targetgene="Ccdc125"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(112.51777380672135,176.9982456680999,131.16750940866106,244.81236309445376,54.417825979520465,130.06376866062934,50.39883615380625,392.4665965862776),c(4.099837911985798,189.17272817437131,31.189668361017127,0.9238202380922783,949.4478585374228,910.4463806244053,45.229724753415866,1536.7667938619304),c(84.27444596859696,14.047479814928565,42.2983173663109,302.55112797522116,27.208912989760233,771.4126968837326,0.0,479.94409100611057),c(0.0,47.761431370757116,6.40883596459256,8.314382142830505,34.36915325022345,58.30444802028211,20.67644560156154,44.92087551288719))
targetgene="Slc25a26"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(201.8031327788565,54.31692195105711,39.307527249501035,25.405056547537654,650.1498156500603,94.18410834045572,441.95902473337793,68.56344157230151),c(106.14024816585456,73.51514436479282,194.82861332361384,35.10516904750658,665.9023442230794,1372.3970072466404,2024.9993911029333,158.4051925980759),c(50.56466758115818,14.51572914209285,93.99626081402423,54.50539404744442,309.32237925201105,542.6798623426258,429.03624623240194,167.86221902184164),c(120.71744963069295,1617.333176025442,54.68873356452318,202.7785422612551,5.728192208370575,0.0,120.18184005907645,567.4215854259435))
targetgene="Vmn2r51"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(237.3350613494001,178.40299364959276,191.41056747583113,69.74842797596702,38.66529740650138,0.0,445.8358582836707,761.2906271131409),c(194.51453204643732,196.1964680818356,81.60584461581193,331.18955535608177,671.63053643145,112.12393850054252,748.2288752065082,451.5730117348134),c(141.21663919062195,94.58636408718566,143.1306698759005,347.8183196417428,167.54962209483932,403.64617860195307,3288.8471284983825,1593.5089524045247),c(54.66450549314398,33.71395155582855,29.907901168098615,46.191011904613916,526.9936831700929,49.33453294023871,684.907260551726,383.00957016251186))
targetgene="Sirt7"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1866_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1868_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1869_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399")

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
Sweave("invivo.poolC_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolC_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.tex",pdf=TRUE);

