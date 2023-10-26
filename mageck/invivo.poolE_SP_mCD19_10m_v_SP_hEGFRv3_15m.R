pdf(file='invivo.poolE_SP_mCD19_10m_v_SP_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolE_SP_mCD19_10m_v_SP_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Odr4","Atp13a1","Igfbp4","Chid1","Gprasp1","Nhlrc2","Vmn1r124","Slc28a2","Sybu","Tmem145")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='23,24,25,26,27_vs_17,18,19,20,21,22 neg.'


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
targetmat=list(c(191.35518566894737,148.29433931680794,122.9339570660172,137.1768567424801,250.00312625886727,225.4317901931959,71.81077019182278,66.69900194691267,41.538433284633214,118.1714504704808,95.43487994888902),c(439.9196536512914,350.1609792245046,307.0260133759324,348.8549374054451,339.82460874708306,317.9500393476289,205.50103384681202,177.86400519176712,223.9817481034144,231.37771394640362,144.8565142081351),c(270.26454058397724,121.119983944618,100.69464825005429,196.30481223492845,206.5894097228963,121.18587565298971,94.72910110410666,156.39765973758833,180.81435665075634,160.87205862367978,51.12582854404769),c(92.71849202516007,83.85229657704323,32.74120464572317,159.054200274686,59.88098832547719,59.941400860618565,21.390442184798278,22.23300064897089,12.21718626018624,5.958224393469621,37.492274265634975))
targetgene="Odr4"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(171.62784694018993,104.03896056781291,137.14240436510462,134.81173852278218,154.19354493810377,88.60902735917527,63.40738219065204,149.49776298445943,65.15832672099327,48.65883254666857,180.6445941889685),c(485.2925327274336,563.6737714345684,353.975665320743,457.65037551155,464.0776595224482,222.82564232969074,229.1833091228387,335.02832012414757,249.23059970779929,324.72322944409433,385.14790836515925),c(136.11863722842648,173.13946422852445,137.14240436510462,121.21230875951906,70.3601612824357,79.48750983690722,67.22710400936602,84.33207142713096,56.1990567968567,132.07397405524327,143.15231992333352),c(199.24612116045037,358.70149091290716,261.3118785875642,272.5798748201868,254.49420038327804,259.3117124187629,193.27792402692728,203.9302818146985,160.45237955044595,144.98346024109412,83.50551995527789))
targetgene="Atp13a1"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(205.1643227790776,250.00406942414742,230.4239496765046,210.49552155311605,128.74412489977595,264.5240081457732,146.67731783861677,93.53193376463615,147.42071420624728,202.5796293779671,182.3487884737701),c(538.5563472950787,385.8758462850971,321.852219253241,386.1055493656876,497.01220310146067,411.7713624338145,457.60267388193455,350.36142401998956,294.02694932848215,471.6927644830117,603.2847768197628),c(424.1377826682854,211.18356174959038,258.2230856964582,445.8247844130604,410.18477002951875,341.4053701191753,142.8575960199028,149.49776298445943,248.41612062378687,206.55177897361352,120.99779422091287),c(126.25496786404776,128.8840854795294,174.20791905837612,126.53382475383941,161.6786684787884,139.4289106975258,58.82371600819526,129.56472791986485,161.26685863445837,116.18537567265761,153.37748563214308))
targetgene="Igfbp4"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(88.77302427940857,163.82254238663077,200.1537793436662,125.94254519891493,191.619162641527,155.0657978785567,198.62553457312686,86.63203701150726,139.27592336612312,88.38032850313272,81.8013256704763),c(244.61900023659254,402.9568696619022,213.7444680645324,371.3235604925755,279.94362042160583,509.50190731525777,291.06280258600515,441.5933922002494,227.23966443946406,316.7789302528015,231.77042273301618),c(76.9366210421541,173.9158743820156,104.40119971938144,114.11695410042526,200.60131089034857,229.34101198845363,126.81476438130407,49.83258766148648,59.45697313290636,252.23149932354727,47.71743997444451),c(195.30065341469887,229.81740543337776,206.3313651258781,293.2746592425437,270.9614721727843,269.73630387278354,165.77592693218665,79.73214025837837,210.95008275921575,136.04612365088968,76.68874281607154))
targetgene="Chid1"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(193.32791954182312,39.59691782804819,74.74878796476423,33.11165507577106,115.27090252654358,117.27665385773197,67.99104837310881,45.23265649273388,49.68322412475737,82.42210410966308,68.16777139206359),c(1025.821613895388,679.3588843047485,744.3990867565363,883.371655057178,1086.839938107411,874.3626082059794,740.262088466769,532.8253603805092,560.3616098005422,886.7823972280619,487.39956545325464),c(435.9741859055399,412.2737915037959,371.89066408915755,250.7025312879809,303.8960157517967,192.85494189938146,346.0667967754864,350.36142401998956,289.1400748244077,332.66752863538716,167.01103991055578),c(323.52835515162235,243.79278819621828,274.28480873020925,290.90954102284576,191.619162641527,269.73630387278354,168.83170438715783,282.1291116834927,219.0948735993399,260.1757985148401,342.5430512451195))
targetgene="Gprasp1"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(114.41856462679328,105.59178087479519,48.185169101252974,132.44662030308425,112.27685311026973,45.607587611340215,6.87549927368516,65.16569155732847,36.651558780558716,49.651869945580174,35.78807998083338),c(0.0,0.7764101534911411,5.559827203990728,36.65933240531796,25.449420038327805,13.030739317525775,46.600606188310536,20.699690259386692,28.50676794043456,20.853785377143673,13.633554278412717),c(13.809137110130223,22.515894451243092,15.443964455529798,1.7738386647734496,17.964296497643154,3.909221795257732,0.7639443637427956,0.7666551947920996,19.547498016297983,9.930373989116035,0.0),c(9.863669364378731,0.7764101534911411,17.914998768414566,4.730236439395866,17.964296497643154,67.75984445113403,26.738052730997847,22.99965584376299,17.918539848273152,13.902523584762449,34.08388569603179))
targetgene="Nhlrc2"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(295.9100809313619,270.1907334149171,170.50136758904898,235.9205424148688,109.28280369399586,276.25167353154643,95.49304546784946,68.99896753128897,61.0859313009312,311.8137432582435,110.77262851210332),c(205.1643227790776,271.7435537218994,183.474297731694,195.71353268000396,169.16379201947305,347.9207397779382,276.54785967489204,202.39697142511432,179.18539848273153,286.9878082854534,238.58719987222256),c(256.455403473847,307.4584207824919,444.168417741037,214.04319888266295,151.1994955218299,281.46396925855674,362.10962841408514,243.7963519438877,603.5290012532002,454.8111287015144,391.9646855043656))
targetgene="Vmn1r124"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(311.6919519143679,241.46355773574487,188.41636635746355,134.81173852278218,197.6072614740747,121.18587565298971,94.72910110410666,151.03107337404364,60.27145221691878,81.42906671075148,122.70198850571445),c(228.83712925358657,165.37536269361306,174.82567763659733,119.43847009474561,113.77387781840666,62.547548724123715,48.128494915796125,36.79944935002078,117.2849880977879,58.589206535784605,136.33554278412717),c(128.2277017369235,153.72921039124594,190.26964209212713,131.85534074815976,226.0507309286764,101.63976667670104,165.77592693218665,117.29824480319125,113.21259267772582,71.49869272163545,281.1920569922623),c(203.19158890620187,144.41228854935224,164.94154038505826,245.38101529366054,110.7798284021328,217.61334660268042,100.84065601404902,252.9962142813929,251.67403695983654,210.52392856925994,202.79911989138915))
targetgene="Slc28a2"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(266.31907283822574,138.2010073214231,135.90688720866223,179.15770514211843,221.5596568042656,109.4582102672165,170.35959311464342,122.66483116673595,107.5112390896389,68.51958052490065,80.09713138567471),c(49.318346821893655,36.49127721408363,110.57878550159336,108.79543810610492,112.27685311026973,119.88280172123713,16.04283163859871,14.566448701049893,21.176456184322817,63.55439353034262,15.337748563214307),c(98.63669364378731,55.125120897871014,37.06551469327152,34.29421418562003,31.43751887087552,24.75840470329897,32.08566327719742,39.86607012918918,63.52936855296845,84.4081789074863,25.562914272023846),c(90.74575815228432,41.149738135030475,27.79913601995364,77.4576216951073,127.24710019163902,37.78914402082474,46.600606188310536,39.86607012918918,37.466037864571135,29.791121967348104,11.929359993611127))
targetgene="Sybu"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(0.0,3.1056406139645643,17.914998768414566,85.73553546405007,25.449420038327805,27.364552566804125,2.291833091228387,6.899896753128897,0.0,0.0,0.0),c(370.8739681006403,246.89842881018285,211.89119232986883,346.4898191857472,417.6698935702034,280.16089532680417,151.26098402107354,340.39490648769225,231.31205985952613,188.67710579320465,403.8940454979767),c(126.25496786404776,82.29947627006095,117.37412986202648,185.07050069136326,110.7798284021328,58.638326928865986,119.93926510761891,125.73145194590434,68.41624305704295,66.53350572707743,120.99779422091287),c(59.182016186272385,51.24307013041531,98.84137251539072,34.29421418562003,98.80363073703735,56.03217906536083,39.72510691462537,10.733172727089395,18.733018932285567,4.965186994558017,98.8432685184922))
targetgene="Tmem145"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetgenelist=c("C2cd5","Grk2","Rerg","Get4","Gm7452_Zfp59","Zic4","Bop1","Tmem161a","Nosip","Noc2l")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='23,24,25,26,27_vs_17,18,19,20,21,22 pos.'


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
targetmat=list(c(260.4008712195985,172.36305407503332,234.74825972405296,131.85534074815976,269.4644474646473,205.88568121690724,248.28191821640857,210.0635233730353,282.62424215230834,211.51696596817155,173.82781704976213),c(266.31907283822574,490.69121700640113,557.2182375555151,426.90383865547693,238.02692859377183,509.50190731525777,642.4772099076911,462.2930824596361,452.0358916268909,373.3820619907629,691.9028796294454),c(311.6919519143679,421.5907133456896,406.4851444695443,518.5521696687719,285.93171925415356,220.21949446618558,830.4075233884188,627.8906045347296,891.8545969935955,724.9173012054705,603.2847768197628),c(59.182016186272385,64.4420427397647,101.93016540649667,157.28036160991255,119.76197665095438,61.24447479237114,126.05082001756128,387.16087337001034,324.977154520954,44.686682951022156,340.8388569603179))
targetgene="C2cd5"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(266.31907283822574,428.57840472710984,536.2144458959946,545.7510291952981,443.1193136085312,379.19451414,584.4174382632386,949.1191311526194,713.4836775948764,664.3420198718627,371.51435408674655),c(491.2107343460608,606.3763298765812,560.3070304466211,488.9881919225476,302.3989910436598,555.109494926598,656.9921528188042,1132.3497227079313,949.6826119584771,804.3602931183988,632.2560796613898),c(976.5032670734944,385.8758462850971,651.7353000233575,629.7127259945746,823.3635894753113,796.1781723008248,961.0420095884369,1380.7460058205716,1168.7774855578168,1410.113106454477,1121.359839399446),c(767.3934765486653,886.660395286883,735.7504666614396,1008.131641146244,547.9110431781163,602.0201564696907,579.0698277170391,1025.7846506318294,935.0219884462535,680.2306182544484,783.9293710087312))
targetgene="Grk2"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(201.21885503332612,184.00920637740043,229.80619109828342,102.88264255686009,115.27090252654358,204.58260728515467,336.8994644105729,292.86228441058205,253.30299512786138,242.30112533443125,305.05077697948457),c(27.618274220260446,63.66563258627357,12.35517156442384,43.75468706441176,23.952395330190875,22.152256839793818,65.69921528188043,41.399380518773384,76.5610338971671,52.630982142314984,69.87196567686517),c(122.30950011829627,122.67280425160028,134.67137005221986,119.43847009474561,121.2590013590913,132.9135410387629,223.0717542128963,209.29686817824322,171.04060764260734,138.03219844871288,262.4459198594448),c(128.2277017369235,107.9210113352686,114.28533697092051,177.38386647734498,217.0685826798548,118.57972778948455,220.77992112166794,113.46496882923074,99.36644824951475,228.39860174966879,122.70198850571445))
targetgene="Rerg"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(57.20928231339664,24.068714758225372,45.714134788368206,17.147107092810014,32.934543579012455,32.576848293814436,73.33865891930839,96.59855454380455,140.90488153414796,103.27588948680676,228.362034163413),c(220.94619376208357,227.48817497290432,164.32378180683705,182.11410291674085,158.68461906251454,207.18875514865982,253.62952876260815,146.43114220529102,319.27580093286707,177.75369440517701,219.84106273940506),c(276.18274220260446,175.46869468899789,126.02274995712315,104.65648122163354,164.67271789506228,96.42747094969073,126.81476438130407,389.46083895438665,285.882158488358,289.96692048218824,129.51876564492082),c(17.754604855881716,13.975382762840539,35.82999753682913,34.88549374054451,40.4196671196971,13.030739317525775,74.86654764679398,59.03244999899167,5.701353588086912,44.686682951022156,30.675497126428613))
targetgene="Get4"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(126.25496786404776,149.0707494702991,241.54360408448605,144.8634909564984,188.62511322525313,177.21805471835054,401.0707909649677,413.9938051877338,248.41612062378687,297.911219673481,405.5982397827783))
targetgene="Gm7452_Zfp59"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(65.10021780489963,13.198972609349399,9.26637867331788,24.242461751903814,64.37206244988798,6.5153696587628875,45.836661824567734,22.99965584376299,30.13572610845939,53.62401954122659,69.87196567686517),c(118.36403237254477,113.3558824097066,218.06877811208076,120.62102920459458,223.05668151240252,129.00431924350517,89.38149055790709,107.33172727089395,198.7328964990295,202.5796293779671,253.92494843543685),c(41.42741133039067,46.58460920946846,113.04981981447813,70.95354659093799,80.8393342393942,96.42747094969073,90.14543492164988,55.96582921982328,48.86874504074496,78.44995451401667,144.8565142081351),c(84.82755653365709,121.119983944618,165.55929896327945,128.30766341861286,248.50610155073034,229.34101198845363,288.00702513103397,348.06145843561325,236.19893436360064,347.5630896190612,100.54746280329378))
targetgene="Zic4"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(57.20928231339664,12.422562455858257,32.74120464572317,92.23961056821939,56.88693890920333,53.42603120185568,20.62649782105548,103.49845129693345,84.70582473729127,5.958224393469621,56.238411398452456),c(185.43698405032015,96.27485903290149,112.43206123625694,101.70008344701112,161.6786684787884,142.03505856103095,180.29086984329976,80.49879545317046,111.58363450970099,118.1714504704808,88.61810280968265),c(82.85482266078134,57.45435135834444,50.03844483591655,41.98084839963831,34.43156828714938,24.75840470329897,117.64743201639052,117.29824480319125,43.167391452658045,192.6492553888511,28.971302841627022),c(45.37287907614216,90.83998795846351,54.362754883464895,162.0105980493084,44.91074124410789,65.15369658762887,83.26993564796473,161.76424610113304,176.74196123069427,145.9764976400057,131.2229599297224))
targetgene="Bop1"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(86.80029040653284,93.16921841893692,66.71792644788873,26.016300416677264,61.37801303361412,71.66906624639176,90.90937928539267,154.09769415321202,68.41624305704295,245.28023753116605,248.81236558103208),c(380.737637465019,136.64818701444082,220.53981242496553,241.24205840918916,209.58345913917015,185.036498308866,158.1364832947587,259.1294558397297,173.4840448946446,132.07397405524327,143.15231992333352),c(55.23654844052089,107.14460118177746,92.6637867331788,92.83089012314387,53.89288949292947,123.79202351649487,164.24803820470106,82.03210584275466,61.0859313009312,62.56135613143102,66.463577107262),c(181.49151630456865,118.01434333065345,328.0298050354529,148.4111682860453,395.2145229481494,254.0994166917526,359.05385095911396,308.96204350121616,319.27580093286707,388.27762297443695,431.1611540548022))
targetgene="Tmem161a"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(55.23654844052089,85.40511688402552,223.0108467378503,204.58272600387122,109.28280369399586,87.3059534274227,237.58669712400945,235.3631448011746,316.8323636808298,172.788507410619,235.17881130261938),c(171.62784694018993,34.162046753610205,46.94965194481059,13.599429763263116,82.33635894753114,37.78914402082474,166.53987129592946,112.69831363443865,100.99540641753958,4.965186994558017,17.041942848015896),c(218.97345988920782,194.87894852627642,224.86412247251388,272.5798748201868,166.16974260319918,230.6440859202062,126.05082001756128,194.73041947719332,370.5879832256493,177.75369440517701,195.9823427521828),c(106.5276291352903,49.69024982343303,90.19275242029403,57.35411682767487,13.473222373232367,79.48750983690722,15.278887274855911,300.52883635850304,58.64249404889395,57.596169136873,46.01324568964292))
targetgene="Nosip"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(124.282233991172,53.57230059088873,50.03844483591655,47.893643948883145,85.33040836380499,24.75840470329897,66.46315964562322,54.432518830239076,75.74655481315469,14.895560983674052,25.562914272023846),c(25.6455403473847,97.82767933988377,32.12344606750198,18.32966620265898,19.461321205780084,15.636887181030929,97.02093419533504,100.43183051776505,60.27145221691878,75.47084231728186,139.74393135373035),c(37.481943584639176,24.068714758225372,16.06172303375099,18.32966620265898,52.395864784792536,39.09221795257732,25.210164003512254,59.799105193783774,39.909475116608384,10.923411388027638,6.8167771392063585),c(69.04568555065111,55.125120897871014,85.25068379452449,27.19885952652623,73.35421069870955,62.547548724123715,29.029885822226234,117.29824480319125,17.918539848273152,61.568318732519415,139.74393135373035))
targetgene="Noc2l"
collabel=c("In.vivo_AR1895_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1890_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._SP_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1905_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._SP_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
Sweave("invivo.poolE_SP_mCD19_10m_v_SP_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolE_SP_mCD19_10m_v_SP_hEGFRv3_15m_summary.tex",pdf=TRUE);

