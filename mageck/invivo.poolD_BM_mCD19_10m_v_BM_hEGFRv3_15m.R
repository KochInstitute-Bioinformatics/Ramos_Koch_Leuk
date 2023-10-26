pdf(file='invivo.poolD_BM_mCD19_10m_v_BM_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolD_BM_mCD19_10m_v_BM_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Nsmce1","Dhps","Nr1h5","Tssk4","Ang","Sirt4","Klhl29","Tacstd2","Zpld2","Vmn2r112")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='6,7,8,9_vs_1,2,3,4,5 neg.'


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
targetmat=list(c(378.6463428403889,91.95723197948848,323.02333035659507,842.2006823926387,121.92806589059623,563.5432409345161,83.66881565769575,154.01250563992753,21.33226479332245),c(76.80090916102229,249.36150293536966,68.00491165402002,133.9670791647865,189.99833580497815,0.0,2.497576586796888,3.1431123599985207,118.3940696029396),c(148.24361535732208,101.07011082430266,136.00982330804004,113.48803521602933,81.53471890843551,163.00837547692615,0.0,4.714668539997781,164.25843890858286),c(16.074608894167454,40.59373303599041,56.67075971168335,41.81138139537923,2.2440748323422617,0.0,0.0,260.8783258798772,0.0))
targetgene="Nsmce1"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(33.04225161578866,6.627548250773945,0.0,2.5598804935946466,0.0,0.0,0.0,0.0,0.0),c(112.52226225917218,41.422176567337154,14.734397525037672,268.7874518274379,14.212473938167658,0.0,0.0,1.5715561799992603,10.666132396661226),c(22.325845686343687,59.6479342569655,10.200736748103003,43.517968391108994,8.976299329369047,26.779947399780724,0.0,0.0,4.26645295866449),c(16.074608894167454,41.422176567337154,18.134643107738672,0.8532934978648822,65.0781701379256,119.9275905294528,58.69304978972687,97.43648315995414,97.06180480961716))
targetgene="Dhps"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.8930338274537475,14.911983564241375,105.40761306373103,0.8532934978648822,5.236174608798611,0.0,0.0,0.0,0.0),c(0.0,0.0,0.0,0.0,61.338045417355154,0.0,0.0,17.287117979991862,0.0),c(0.0,49.70661188080459,5.667075971168335,13.652695965838115,0.0,0.0,0.0,0.0,0.0),c(0.0,0.0,0.0,0.0,2.992099776456349,0.0,0.0,0.0,0.0))
targetgene="Nr1h5"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.0,0.0,0.0,1.7065869957297644,0.0,0.0,0.0,0.0,0.0),c(0.0,89.47190138544825,0.0,97.27545875659656,79.29064407609324,0.0,0.0,0.0,0.0),c(0.0,7.455991782120687,0.0,0.0,0.0,0.0,0.0,56.576022479973375,0.0),c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0))
targetgene="Tssk4"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(207.1838479692694,106.86921554372987,160.94495758118072,376.30243255841305,285.7455286515813,513.4763827523174,197.30855035695416,237.30498317988832,293.3186409081837),c(17.86067654907495,84.50124019736779,48.73685335204768,40.958087897514346,39.64532203804662,171.15879425077244,2.497576586796888,334.7414663398425,231.45507300754858),c(210.7559832790844,391.02534679566276,233.4835300121354,151.03294912208415,309.6823268632321,115.2702083729692,162.34247814179773,50.28979775997633,811.6926753859193),c(1660.1498852365166,1603.0382331559479,1840.6662754354752,448.83237987692803,507.16091210935116,679.9777948466062,91.1615454180864,47.14668539997781,551.4390449073853))
targetgene="Ang"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(108.9501269493572,90.300344916795,124.67567136570337,39.25150090178458,97.24324273483134,587.994497256055,171.08399619558682,136.72538765993565,290.1188011891853),c(82.15911212574477,110.18298966911684,2.266830388467334,20.479043948757173,90.51101823780456,84.99722435582578,4.995153173593776,4.714668539997781,2.133226479332245),c(23.218879513797436,28.167080065789264,26.06854946737434,52.904196867622694,34.409147429248016,0.0,148.60580691441484,1.5715561799992603,0.0),c(57.15416495703984,0.8284435313467431,0.0,0.0,0.7480249441140873,0.0,0.0,0.0,6.399679437996735))
targetgene="Sirt4"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(135.74114177296963,579.0820284113735,487.36853352047683,559.7605345993627,745.0328443376309,516.96941936968,88.66396883128952,122.58138203994231,148.25924031359102),c(48.22382668250236,153.26205329914748,172.27910952351738,174.07187356443598,13.464448994053571,0.0,136.1179239804304,6.286224719997041,17.06581183465796),c(73.2287738512073,140.83540032894632,32.86904063277634,76.7964148078394,224.4074832342262,11.64345539120901,18.73182440097666,383.45970791981955,25.59871775198694),c(5.358202964722485,0.0,0.0,0.0,29.92099776456349,8.150418773846306,78.67366248410197,0.0,0.0))
targetgene="Klhl29"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(81.26607829829102,3.3137741253869724,0.0,59.730544850541754,78.54261913197917,0.0,0.0,0.0,13.865972115659593),c(0.0,0.0,0.0,2.5598804935946466,0.0,0.0,0.0,0.0,0.0),c(0.0,49.70661188080459,0.0,0.0,0.0,0.0,0.0,0.0,0.0),c(0.0,0.8284435313467431,0.0,0.0,0.0,0.0,0.0,0.0,0.0))
targetgene="Tacstd2"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.0,0.0,81.60589398482402,53.75749036548758,92.00706812603273,0.0,0.0,0.0,0.0),c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.133226479332245),c(0.0,1.6568870626934862,4.533660776934668,0.0,2.2440748323422617,0.0,0.0,0.0,0.0))
targetgene="Zpld2"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.0,0.0,2.266830388467334,5.119760987189293,23.936798211650792,0.0,2.497576586796888,0.0,0.0),c(0.0,0.0,0.0,0.0,5.984199552912698,0.0,0.0,0.0,0.0),c(0.8930338274537475,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),c(0.0,1.6568870626934862,9.067321553869336,3.4131739914595287,0.0,0.0,0.0,0.0,6.399679437996735))
targetgene="Vmn2r112"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetgenelist=c("Rab4b","Nucks1","Fastk","1700109H08Rik","Gm9726_Lysmd1","Kdm1b","Cfap221","Vmn1r104_Vmn1r111_Vmn1r114_Vmn1r115_Vmn1r116_Vmn1r120_Vmn1r121_Vmn1r123_Vmn1r126_Vmn1r127_Vmn1r128_Vmn1r130_Vmn1r131_Vmn1r135_Vmn1r137_Vmn1r138_Vmn1r142_Vmn1r143_Vmn1r152_Vmn1r159_Vmn1r160_Vmn1r163_Vmn1r165_Vmn1r166_Vmn1r178_Vmn1r180_Vmn1r183_Vmn1r242_Vmn1r244_Vmn1r246_Vmn1r247_Vmn1r248_Vmn1r250_Vmn1r252_Vmn1r254_Vmn1r256_Vmn1r257_Vmn1r258_Vmn1r259_Vmn1r93_Vmn1r94_Vmn1r95","Gm6346_Gm6348_Gm6468_LOC108168681_LOC108168682_LOC620551_LOC620639_Pramel40_Pramel44","Best1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='6,7,8,9_vs_1,2,3,4,5 pos.'


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
targetmat=list(c(16.074608894167454,122.60964263931798,6.800491165402002,65.70359933559592,26.928897988107142,123.4206271468155,23.726977574570437,0.0,0.0),c(118.77349905134841,37.27995891060344,674.3820405690319,95.5688717608668,201.96673491080355,68.69638680813316,560.7059437359013,110.00893259994822,175.99118454491023),c(20.539778031436192,52.191942474844815,13.600982330804005,22.185630944486938,175.7858618668105,2.328691078241802,103.64942835207086,55.00446629997411,21.33226479332245),c(26.791014823612425,5.799104719427202,222.14937806979873,11.092815472243469,3.7401247205704364,79.17549666022127,961.5669859168019,6.286224719997041,842.6244593362368))
targetgene="Rab4b"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(553.6809730213234,595.6508990383082,278.82013778148206,962.5150655915871,294.72182798095037,317.866332180006,88.66396883128952,265.59299441987497,357.31543528815104),c(350.069260361869,323.0929772252298,290.15428972381875,251.72158187014026,270.78502976929957,2892.234319176318,719.3020569975038,160.29873035992455,1318.3339642273274),c(368.8229707383977,554.2287224709711,287.8874593353514,391.66171551998093,418.89396870388885,1283.1087841112328,183.57187912957127,1060.8004214995008,435.178201783778),c(116.98743139644093,191.37045574109766,116.7417650060677,223.56289644059913,145.86486410224703,908.1895205143028,524.4910832273465,746.4891854996487,141.8595608755943))
targetgene="Nucks1"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(80.37304447083727,121.78119910797123,90.67321553869336,338.7575186523582,160.07733804041467,18.629528625934416,24.97576586796888,111.58048877994749,190.92376990023592),c(828.7353918770776,637.0730756056454,760.5215953307905,553.7874801143086,363.5401228394464,1391.3929192494768,363.3973933789472,985.3657248595363,2069.2296849522777),c(97.34068719245847,104.38388494968963,121.27542578300236,134.82037266265138,334.367150018997,555.3928221606698,319.68980311000166,295.45256183986095,623.9687452046817),c(72.33574002375354,132.55096501547888,38.53611660394468,40.958087897514346,50.865696199757934,372.5905725186883,349.6607221515643,14.144005619993344,156.79214623092003))
targetgene="Fastk"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(33.04225161578866,2.4853305940402293,43.069777380879344,43.517968391108994,11.220374161711309,48.90251264307784,96.15669859168018,34.574235959983724,44.79775606597715),c(249.15643785959554,27.33863653444252,5.667075971168335,52.05090336975781,296.96590281329264,53.55989479956145,810.4636024155901,350.4570281398351,215.45587441255674),c(46.43775902759487,183.08602042763022,49.870268546281345,30.71856592313576,139.8806645493343,763.810673663311,550.7156373887138,23.573342699988906,36.264850148648165),c(131.27597263570087,138.3500697349061,387.6279964279141,45.22455538683876,210.9430342401726,628.7465911252865,49.95153173593776,1596.7010788792486,151.4590800325894))
targetgene="1700109H08Rik"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(157.17395363185955,36.4515153792567,60.07100529438435,255.13475586159976,48.62162136741567,32.601675095385225,104.89821664546929,473.0384101797774,446.91094742010534),c(228.61665982815936,34.79462831656321,123.5422561714697,112.63474171816445,86.02286857312004,528.6128747608891,494.5201641857838,183.87207305991345,44.79775606597715),c(253.62160699686427,207.11088283668576,129.20933214263803,226.12277693419378,445.0748417478819,89.65460651230937,268.4894830806655,419.6055000598025,57.59711494197062),c(218.79328772616813,291.6121230340536,89.53980034445969,120.31438319894839,115.94386633768353,331.83847864945676,571.9450383764873,535.9006573797478,106.66132396661226))
targetgene="Gm9726_Lysmd1"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(392.04185025219516,154.91894036184095,189.28033743702238,160.41917759859786,237.12390728416565,96.64067974703478,287.2213074816421,465.18062927978104,241.0545921645437),c(36.61438692560365,66.27548250773944,44.203192575113015,151.88624261994903,281.2573789868968,386.5627189881391,42.4588019755471,278.1654438598691,113.06100340460898),c(97.34068719245847,4.9706611880804585,4.533660776934668,62.2904253441364,18.700623602852183,0.0,31.2197073349611,48.71824157997707,260.2536304785339),c(0.8930338274537475,24.02486240905555,0.0,6.826347982919057,0.7480249441140873,0.0,0.0,284.4516685798661,214.38926117289063))
targetgene="Kdm1b"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(51.79596199231735,55.505716600231786,71.40515723672102,187.72456953027407,124.17214072293848,742.8524539591349,564.4523086160967,1.5715561799992603,94.9285783302849),c(69.6566385413923,196.3411169291781,145.07714486190937,48.63772937829828,41.889396870388886,31.437329556264327,364.6461816723456,245.16276407988462,351.98236908982045),c(58.04719878449359,130.8940779527854,153.01105122154505,369.476084575494,107.71559195242857,151.36492008571713,27.473342454765767,73.86314045996524,43.73114282631102),c(109.84316077681095,107.6976590750766,194.94741340819073,54.61078386335246,273.02910460164185,589.1588427951759,369.6413348459394,246.73432025988387,280.51928203219023))
targetgene="Cfap221"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(82.15911212574477,0.0,0.0,0.8532934978648822,0.0,1.164345539120901,0.0,121.00982585994305,212.25603469355838))
targetgene="Vmn1r104_Vmn1r111_Vmn1r114_Vmn1r115_Vmn1r116_Vmn1r120_Vmn1r121_Vmn1r123_Vmn1r126_Vmn1r127_Vmn1r128_Vmn1r130_Vmn1r131_Vmn1r135_Vmn1r137_Vmn1r138_Vmn1r142_Vmn1r143_Vmn1r152_Vmn1r159_Vmn1r160_Vmn1r163_Vmn1r165_Vmn1r166_Vmn1r178_Vmn1r180_Vmn1r183_Vmn1r242_Vmn1r244_Vmn1r246_Vmn1r247_Vmn1r248_Vmn1r250_Vmn1r252_Vmn1r254_Vmn1r256_Vmn1r257_Vmn1r258_Vmn1r259_Vmn1r93_Vmn1r94_Vmn1r95"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(117.88046522389467,128.4087473587452,62.33783568285168,89.59581727581264,54.60582092032837,299.23680355407157,57.44426149632842,656.9104832396908,7.4662926776628575),c(172.35552869857327,65.4470389763927,92.9400459271607,123.72755719040792,279.01330415455453,322.5237143364896,1090.1921801368417,847.0687810196014,212.25603469355838),c(151.81575066713708,14.911983564241375,15.867812719271338,21.332337446622056,81.53471890843551,27.944292938901626,43.70759026894554,583.0473427797256,53.33066198330613),c(189.32317142019446,71.24614369581991,116.7417650060677,25.598804935946465,63.58212024969742,93.14764312967208,22.478189281171993,6.286224719997041,12.79935887599347),c(110.73619460426468,29.82396712848275,95.20687631562802,77.64970830570428,12.716424049939484,3.493036617362703,27.473342454765767,47.14668539997781,31.998397189983677),c(66.97753705903106,62.133264851005734,31.735625438542677,20.479043948757173,11.220374161711309,436.6295771703379,61.19062637652375,886.3576855195828,34.13162366931592))
targetgene="Gm6346_Gm6348_Gm6468_LOC108168681_LOC108168682_LOC620551_LOC620639_Pramel40_Pramel44"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(0.0,0.8284435313467431,0.0,41.81138139537923,14.212473938167658,220.0613068938503,23.726977574570437,3.1431123599985207,174.9245713052441),c(112.52226225917218,133.37940854682563,45.33660776934668,104.10180673951562,91.25904318191864,0.0,0.0,45.57512921997855,67.19663409896572),c(116.09439756898718,61.30482131965899,47.60343815781401,325.958116184385,98.73929262305951,59.38162249516595,231.02583427871213,411.7477191598062,1443.1277132682637),c(23.218879513797436,19.882644752321834,24.935134273140672,163.83235159005739,41.889396870388886,33.76602063450613,814.2099672957855,1.5715561799992603,25.59871775198694))
targetgene="Best1"
collabel=c("In.vivo_AR1870_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1871_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1872_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1873_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1874_Pool.D..25.32._BM_hEGFRv3_15m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1883_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1884_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1885_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430","In.vivo_AR1886_Pool.D..25.32._BM_mCD19_10m_1400..1401..1402..1403..1427..1428..1429..1430")

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
Sweave("invivo.poolD_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolD_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.tex",pdf=TRUE);

