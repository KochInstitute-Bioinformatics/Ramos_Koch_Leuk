pdf(file='invivo.poolC_BM_mCD19_15m_v_BM_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolC_BM_mCD19_15m_v_BM_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("LTO1","Med15","Defa23_Defa31","Gm20738_Gm20747_Gm20776_Gm20798_Gm20800_Gm20801_Gm20803_Gm20804_Gm20805_Gm20806_Gm20808_Gm20809_Gm20810_Gm20813_Gm20815_Gm20816_Gm20818_Gm20823_Gm20825_Gm20826_Gm20840_Gm20841_Gm20842_Gm20844_Gm20846_Gm20847_Gm20848_Gm20852_Gm20853_Gm20854_Gm20859_Gm20860_Gm20861_Gm20862_Gm20863_Gm20865_Gm20866_Gm20867_Gm20868_Gm20879_Gm20880_Gm20881_Gm20882_Gm20884_Gm20889_Gm20891_Gm20892_Gm20893_Gm20895_Gm20898_Gm20902_Gm20907_Gm20909_Gm20912_Gm20913_Gm20914_Gm20917_Gm20919_Gm20921_Gm20924_Gm20925_Gm20927_Gm20930_Gm20932_Gm20933_Gm20934_Gm20935_Gm20936_Gm21065_Gm21076_Gm21111_Gm21114_Gm21118_Gm21127_Gm21151_Gm21155_Gm21160_Gm21163_Gm21180_Gm21184_Gm21201_Gm21242_Gm21244_Gm21245_Gm21247_Gm21249_Gm21257_Gm21268_Gm21275_Gm21281_Gm21282_Gm21285_Gm21287_Gm21292_Gm21302_Gm21308_Gm21316_Gm21330_Gm21333_Gm21340_Gm21344_Gm21380_Gm21387_Gm21394_Gm21396_Gm21412_Gm21427_Gm21435_Gm21443_Gm21462_Gm21469_Gm21470_Gm21476_Gm21506_Gm21524_Gm21530_Gm21573_Gm21617_Gm21634_Gm21642_Gm21654_Gm21661_Gm21672_Gm21683_Gm21724_Gm21745_Gm21751_Gm21753_Gm21763_Gm21783_Gm21805_Gm21806_Gm21810_Gm21838_Gm21843_Gm21908_Gm21918_Gm21943_Gm28093_Gm28371_Gm28459_Gm29317_Gm31186_Gm38185_Ssty2","Fbrsl1","Ndufb5","Spin1","Sepsecs","Pacsin3","Zc3hc1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='9,10,11,12,13_vs_1,2,3,4,5 neg.'


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
targetmat=list(c(35.8853413385924,64.37499767425199,99.12506247958072,128.5076545020525,6.001105757844376,1.1566748773828543,0.0,0.0,9.132452575144063,103.11670454537459),c(219.9724819716313,0.9980619794457671,17.04412839015415,146.66634481212512,127.52349735419298,15.036773405977106,0.0,293.94655907404587,45.66226287572031,9.497591208126607),c(0.46604339400769346,27.446704434758598,28.705900446575413,92.19027388190722,7.501382197305469,3.4700246321485633,0.0,11.876628649456398,3.652981030057625,10.854389952144693),c(5.126477334084628,58.38662579757738,1.7941187779109633,0.0,1.500276439461094,34.70024632148563,11.534884310003712,2.9691571623640995,9.132452575144063,16.28158492821704))
targetgene="LTO1"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(271.2372553124776,244.02615397449006,228.75014418364782,166.22185745374182,225.0414659191641,277.60197057188503,32.682172211677184,65.32145757201019,34.703319785547436,29.849572368397908),c(0.0,8.982557815011905,94.19123584032558,69.37550913335443,52.50967538113829,3.4700246321485633,15.379845746671617,2.9691571623640995,0.0,9.497591208126607),c(1.3981301820230803,3.493216928060185,2.691178166866445,5.587289326176196,0.0,12.723423651211398,0.0,2.9691571623640995,14.6119241202305,24.422377392325558),c(388.21414720840863,43.41569610589087,248.03692104619068,187.63979987075055,19.50359371299422,5.7833743869142715,144.1860538750464,20.784100136548698,10.958943090172875,23.065578648307472))
targetgene="Med15"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(0.0,12.974805732794973,59.65444936553953,50.285603935585755,72.01326909413251,15.036773405977106,0.0,0.0,5.4794715450864375,18.995182416253215),c(27.03051685244622,216.08041855000857,99.57359217405846,439.0678195486793,1414.7606824118116,84.43726604894837,24.992249338341377,47.50651459782559,9.132452575144063,16.28158492821704),c(3.262303758053854,0.9980619794457671,28.705900446575413,18.62429775392065,1.500276439461094,0.0,0.0,0.0,9.132452575144063,4.07039623205426),c(54.527077098900136,30.440890373095897,165.0589275678086,80.0844803418588,91.51686280712673,347.0024632148563,24.992249338341377,47.50651459782559,279.4530487994083,115.32789324153737))
targetgene="Defa23_Defa31"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(62.449814797030925,82.83914429399867,82.5294637839043,23.74597963624883,97.5179685649711,1.1566748773828543,0.0,0.0,9.132452575144063,12.211188696162779))
targetgene="Gm20738_Gm20747_Gm20776_Gm20798_Gm20800_Gm20801_Gm20803_Gm20804_Gm20805_Gm20806_Gm20808_Gm20809_Gm20810_Gm20813_Gm20815_Gm20816_Gm20818_Gm20823_Gm20825_Gm20826_Gm20840_Gm20841_Gm20842_Gm20844_Gm20846_Gm20847_Gm20848_Gm20852_Gm20853_Gm20854_Gm20859_Gm20860_Gm20861_Gm20862_Gm20863_Gm20865_Gm20866_Gm20867_Gm20868_Gm20879_Gm20880_Gm20881_Gm20882_Gm20884_Gm20889_Gm20891_Gm20892_Gm20893_Gm20895_Gm20898_Gm20902_Gm20907_Gm20909_Gm20912_Gm20913_Gm20914_Gm20917_Gm20919_Gm20921_Gm20924_Gm20925_Gm20927_Gm20930_Gm20932_Gm20933_Gm20934_Gm20935_Gm20936_Gm21065_Gm21076_Gm21111_Gm21114_Gm21118_Gm21127_Gm21151_Gm21155_Gm21160_Gm21163_Gm21180_Gm21184_Gm21201_Gm21242_Gm21244_Gm21245_Gm21247_Gm21249_Gm21257_Gm21268_Gm21275_Gm21281_Gm21282_Gm21285_Gm21287_Gm21292_Gm21302_Gm21308_Gm21316_Gm21330_Gm21333_Gm21340_Gm21344_Gm21380_Gm21387_Gm21394_Gm21396_Gm21412_Gm21427_Gm21435_Gm21443_Gm21462_Gm21469_Gm21470_Gm21476_Gm21506_Gm21524_Gm21530_Gm21573_Gm21617_Gm21634_Gm21642_Gm21654_Gm21661_Gm21672_Gm21683_Gm21724_Gm21745_Gm21751_Gm21753_Gm21763_Gm21783_Gm21805_Gm21806_Gm21810_Gm21838_Gm21843_Gm21908_Gm21918_Gm21943_Gm28093_Gm28371_Gm28459_Gm29317_Gm31186_Gm38185_Ssty2"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(186.8834009970851,41.419572146999336,17.04412839015415,12.105793540048424,181.53344917479237,209.35815280629663,46.139537240014846,608.6772182846404,1137.90359086295,282.214138755762),c(60.58564122100015,18.963177609469575,138.14714589914416,100.1056004273235,87.01603348874345,0.0,0.0,0.0,43.8357723606915,56.98554724875964),c(311.3169871971392,218.57557349862302,316.661964301285,47.957566716345674,235.54340099539175,555.2039411437701,1049.6744722103379,136.58122946874857,20.091395665316938,28.49277362437982),c(116.97689189593106,499.0309897228836,300.9634249945641,35.85177317629725,304.55611721060205,864.0361334049921,201.86047542506498,68.29061473437429,69.40663957109487,226.58539025102047))
targetgene="Fbrsl1"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(11.185041456184642,2.495154948614418,221.57366907200398,42.370277390169484,1.500276439461094,9.253399019062835,15.379845746671617,65.32145757201019,20.091395665316938,28.49277362437982),c(20.03986594233082,9.481588804734788,5.38235633373289,44.698314609409564,0.0,0.0,0.0,0.0,3.652981030057625,6.783993720090433),c(2.7962603640461605,1.9961239588915343,18.389717473587375,8.846541433112309,7.501382197305469,2.3133497547657087,0.0,0.0,5.4794715450864375,0.0),c(2.7962603640461605,0.49903098972288357,0.4485296944777408,2.3280372192400813,63.01161045736595,0.0,0.0,0.0,0.0,5.427194976072347))
targetgene="Ndufb5"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(110.91832777383104,217.07848052945434,235.0295599063362,269.1211025441534,300.05528789221876,1088.4310596172659,46.139537240014846,35.62988594836919,208.2199187132846,118.04149072957354),c(257.25595349224676,480.5668431031369,97.7794733961475,14.899438203136521,120.02211515688751,124.92088675734827,17.302326465005567,0.0,27.397357725432187,21.708779904289386),c(214.379961243539,667.7034642492182,135.0074380378,298.4543715065784,336.061922439285,13.880098528594253,11.534884310003712,3809.4286393131397,553.4266260537302,37.99036483250643),c(147.26971250643112,179.65115630023809,197.8015952646837,198.34877107925493,520.5959244929995,212.8281774384452,38.44961436667904,5.938314324728199,23.744376695374562,39.347163576524515))
targetgene="Spin1"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(871.5011467943867,193.62402401247883,130.07361139854484,162.49699790295767,70.51299265467141,23.133497547657086,1.922480718333952,2.9691571623640995,14.6119241202305,70.5535346889405),c(263.3145176143468,231.55037923141796,139.04420528809965,42.370277390169484,1.500276439461094,25.446847302422796,7.689922873335808,2.9691571623640995,7.30596206011525,18.995182416253215),c(64.78003176706939,114.27809664654033,136.80155681571094,140.6134480421009,10.501935076227657,160.77780795621675,34.604652930011135,32.660728786005095,27.397357725432187,37.99036483250643),c(488.41347692006275,431.1627751205714,366.0002306938365,726.3476124029054,118.52183871742642,135.33096065379397,174.94574536838965,35.62988594836919,43.8357723606915,177.74063546636935))
targetgene="Sepsecs"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(17.243605578284658,164.68022660855158,18.838247168065116,904.6752633966956,531.0978595692272,17.350123160742815,5.767442155001856,14.845785811820498,14.6119241202305,31.206371112415994),c(116.97689189593106,340.83816598072946,478.1326543132717,236.52858147479228,399.073532896651,33.543571444102774,9.61240359166976,32.660728786005095,135.16029811213212,274.0733462916535),c(59.65355443298476,71.36143153037234,135.45596773227774,530.7924859867386,457.58431403563367,284.5420198361822,401.79847013179597,118.76628649456399,93.15101626646944,563.0714787675059),c(30.758864004507767,40.92054115727645,78.0441668391269,26.07401685548891,303.055840771141,104.1007389644569,3.844961436667904,59.38314324728199,29.223848240461,89.54871710519372))
targetgene="Pacsin3"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(137.01675783826187,177.65503234134655,394.2576014459342,1134.2197332137675,31.505805228682974,133.01761089902826,9.61240359166976,160.3344867676614,25.570867210403375,36.63356608848834),c(34.487211156569316,132.24321227656415,269.56634638112223,101.96803020271557,33.00608166814406,85.59394092633123,0.0,11.876628649456398,16.438414635259313,23.065578648307472),c(175.23231614689274,82.83914429399867,218.88249090513753,117.79868329354812,376.56938630473456,8.09672414167998,3.844961436667904,17.814942974184596,14.6119241202305,25.779176136343647),c(81.09155055733866,56.889532828408726,134.55890834332226,6.984111657720244,13.502487955149846,409.4629065935304,0.0,23.753257298912796,25.570867210403375,35.27676734447025))
targetgene="Zc3hc1"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetgenelist=c("Ccnh","Ralgapa2","Epc2","Slc43a2","Gpr141","Lurap1l","Nmu","Rab2b","Fubp1","Cc2d2b")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='9,10,11,12,13_vs_1,2,3,4,5 pos.'


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
targetmat=list(c(1.3981301820230803,13.97286771224074,9.867653278510298,5.121681882328179,0.0,241.74504937301657,128.80620812837478,745.258447753389,409.133875366454,103.11670454537459),c(64.78003176706939,27.94573542448148,118.86036903660131,25.608409411640896,16.50304083407203,3.4700246321485633,5.767442155001856,121.73544365692808,21.91788618034575,14.924786184198954),c(29.360733822484686,52.89728491062566,4.036767250299667,60.994575144090135,511.594265856233,1530.2808627775164,67.28682514168833,6526.207442876291,1053.8850271716249,33.919968600452165),c(4.660433940076935,5.4893408869517195,328.772266052184,11.640186096200408,1.500276439461094,46.26699509531417,1199.627968240386,255.34751596331256,16.438414635259313,16.28158492821704))
targetgene="Ccnh"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(25.63238667042314,67.36918361258928,126.03684414824517,227.68204004167995,10.501935076227657,50.89369460484559,84.58915160669389,148.45785811820497,65.75365854103725,77.33752840903094),c(120.23919565398491,134.23933623545568,352.09581016502653,369.6923104153249,30.005528789221877,684.7515274106498,3537.364521734472,2523.7835880094844,365.29810300576247,147.89106309797145),c(316.90950792523154,203.6046438069365,306.345781328297,1251.5528090634677,1377.2537714252842,288.0120444683307,190.32559111506126,5154.4568338640765,1441.101016357733,1527.7553857643657),c(514.0458635904859,315.8866164945853,935.6329426805673,121.52354284433225,679.6252270758755,1723.445567300453,2008.9923506589798,442.40441719225083,1501.375203353684,1947.0061976659542))
targetgene="Ralgapa2"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(707.4538721036787,638.759666845291,468.71353072923915,383.1949262869174,306.05639365006317,618.821059399827,3124.0311672926723,970.9143920930605,542.4676829635573,2104.3948519720525),c(78.7613335873002,268.47867247091136,143.9780319273548,493.54389047889725,162.02985546179815,145.74103455023965,167.25582249505382,1003.5751208790656,447.49017618205903,776.0888815783455),c(732.6202153800941,409.20541157276455,333.25756299696144,210.45456461930337,4412.313008455078,1167.0849512793,4014.139739881292,2773.192789648069,1196.3512873438722,351.4108747006844),c(728.8918682280325,19.961239588915344,216.63984243274882,145.7351299244291,106.51962720173766,154.9944335693025,438.3256037801411,519.6025034137174,2379.9171410825425,645.8362021526092))
targetgene="Epc2"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(306.1905098630546,225.0629763650205,194.6618874033395,302.6448385012106,1135.709264672048,1633.2249268645903,855.5039196586087,638.3687899082814,261.18814364912015,971.46790071695),c(96.47098255959254,63.8759666845291,43.05885066986312,222.56035815935178,22.50414659191641,1892.3200993983496,540.2170818518405,14281.645950971319,2646.5847562767494,237.43978020316516),c(957.71917468581,370.78002536410247,653.9562945485461,1159.8281426254086,580.6069820714433,197.79140403246808,173.0232646500557,608.6772182846404,1554.3434282895194,522.3675164469634),c(45.20620921874627,42.916665116167984,28.705900446575413,223.0259656031998,874.6611642058177,89.06396555847978,126.88372741004083,62.35230040964609,31.05033875548981,314.7773086121961))
targetgene="Slc43a2"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(61.983771403023226,89.82557815011904,404.57378441892223,101.03681531501954,157.52902614341485,145.74103455023965,930.4806676736328,3099.80007750812,557.0796070837878,131.6094781697544),c(440.877050731278,281.9525091934292,161.91921970646445,644.4007022856546,150.02764394610938,1666.768498308693,2839.5040209792473,1549.90003875406,400.00142279130995,937.5479321164979),c(349.0665021117624,120.26646852321494,105.40447820226909,59.597752812546084,211.53897796401424,160.77780795621675,105.73643950836737,864.024734247953,157.07818429247786,435.5323968298058),c(91.34450522550792,282.95057117287496,100.0221218685362,595.0463132377648,43.50801674437172,1080.334335475586,499.84498676682756,190.02605839130237,571.6915312040182,495.2315415666016))
targetgene="Gpr141"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(85.75198449741559,226.56006933418914,168.64716512363054,448.3799684256397,4.500829318383282,40.483620708399904,434.48064234347316,825.4256911372197,670.3220190155741,301.2093211720152),c(253.52760634018523,187.6356521358042,305.8972516338192,188.10540731459858,22.50414659191641,143.42768479547394,194.17055255172914,635.3996327459173,652.057113865286,1014.8854605255289),c(43.342035642715494,146.2160799888049,293.3384201884425,80.0844803418588,78.01437485197688,28.916871934571358,797.8294981085901,1291.5833656283833,272.14708673929306,43.41755980857877),c(53.59499031088475,124.25871644099801,11.213242361943522,23.280372192400815,3.000552878922188,86.75061580371407,1336.1240992420967,641.3379470706456,129.68082656704567,29.849572368397908))
targetgene="Lurap1l"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(218.57435178960824,856.3371783644682,618.9709783792823,490.75024581580914,909.1675223134229,672.0281037594384,699.7829814735586,917.4695631705067,3335.1716804426114,1458.5586498194432),c(206.4572235454082,188.13468312552712,137.69861620466642,523.3427668851702,195.0359371299422,3030.4881787430786,822.8217474469315,1451.9178523960447,263.014634164149,120.75508821760971),c(192.00987833116972,320.87692639181415,591.1621373216624,787.3421875469955,207.03814864563097,205.88812817414808,1449.5504616238,1413.3188092853113,1221.9221545542755,407.039623205426),c(478.1605222518935,390.24223396329495,347.1619835257714,293.33268962425024,202.53731932724767,283.3853449587993,449.86048809014477,2416.8939301643773,1103.2002710774027,1432.7794736830995))
targetgene="Nmu"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(966.5739991719562,1365.3487878818094,419.8237940311654,71.7035463525945,1279.735802860313,532.070443596113,251.8449741017477,605.7080611222763,396.3484417612523,639.0522084325188),c(545.736814383009,343.83235191906675,61.00003844897275,478.6444522757607,364.5671747890458,159.6211330788339,836.2791124752691,2654.426503153505,425.57229000171327,111.2574970094831),c(184.08714063303893,259.99514564562236,151.15450703899864,566.1786517191878,132.02432667257625,732.1751973833468,3095.1939565176626,4581.409501527805,5466.6861114812355,1534.539379484456),c(175.23231614689274,631.2742019994478,280.33105904858803,39.11102528323337,469.5865255513224,207.04480305153092,467.16281455515036,671.0295186942865,1936.0799459305413,843.9288187792499))
targetgene="Rab2b"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(17.243605578284658,13.97286771224074,5.38235633373289,4.656074438480163,0.0,34.70024632148563,169.17830321338778,11.876628649456398,235.6172764387168,162.8158492821704),c(462.3150468556319,136.2354601943472,686.6989622454212,586.6653792485005,2197.9049838105025,238.275024740868,2230.0776332673845,2289.220172182721,1488.5897697484822,339.19968600452165),c(41.94390546069241,1371.836190748207,194.21335770886176,417.18426968782256,843.1553589771348,380.54603465895906,1270.7597548187423,3076.046820209207,1609.1381437403838,350.05407595666634),c(222.30269894166977,610.3149004310866,1006.5006344080504,397.16314960235786,475.5876313091668,77.49721678465124,73.05426729669018,380.05211678260474,708.6783198311792,299.85252242799714))
targetgene="Fubp1"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(5.126477334084628,140.22770811213027,56.51474150419534,12.105793540048424,1.500276439461094,89.06396555847978,430.63568090680525,148.45785811820497,213.69939025837107,326.9884973083589),c(84.35385431539251,209.09398469388822,87.46329042315946,143.87270014903703,414.07629729126194,630.3878081736556,993.9225313786532,100.95134352037938,299.54444446472525,188.59502541851404),c(162.18310111467733,311.8943685768022,226.5074957112591,127.11083217050844,174.0320669774869,448.78985242454746,1687.93807069721,3996.485540542078,1411.877168117272,271.3597488036173),c(45.20620921874627,340.3391349910066,10.316182972988038,51.21681882328179,84.01548060982125,0.0,1.922480718333952,11.876628649456398,894.9803523641181,279.50054126772585))
targetgene="Cc2d2b"
collabel=c("In.vivo_AR1850_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1851_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1852_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1853_Pool.C..17.24._BM_hEGFRv3_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1865_Pool.C..17.24._BM_mCD19_10m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1856_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1859_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1860_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1861_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399","In.vivo_AR1862_Pool.C..17.24._BM_mCD19_15m_1392..1393..1394..1395..1396..1397..1398..1399")

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
Sweave("invivo.poolC_BM_mCD19_15m_v_BM_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolC_BM_mCD19_15m_v_BM_hEGFRv3_15m_summary.tex",pdf=TRUE);

