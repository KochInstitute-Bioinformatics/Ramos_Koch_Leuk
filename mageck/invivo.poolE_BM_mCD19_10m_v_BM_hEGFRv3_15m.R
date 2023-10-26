pdf(file='invivo.poolE_BM_mCD19_10m_v_BM_hEGFRv3_15m.pdf',width=4.5,height=4.5);
gstable=read.table('invivo.poolE_BM_mCD19_10m_v_BM_hEGFRv3_15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("D830030K20Rik_Gm10408_Gm3043_Gm8050_Gm8271","Tmem169","Terf2","Rrp15","Psmb4","Polg2","Med31","Tshz3","Gm20737_Gm20772_Gm20773_Gm20777_Gm20793_Gm20795_Gm20807_Gm20827_Gm20831_Gm20833_Gm20836_Gm20851_Gm21198_Gm21310_Gm21425_Gm21440_Gm21660_Gm21794_Gm28646_Gm29109","Eif3j1_Eif3j2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='7,8,9,10,11_vs_1,2,3,4,5,6 neg.'


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
targetmat=list(c(99.766450509106,107.93686525363704,5.6385519165511635,153.8095148174042,110.74444482687852,10.484452669672798,68.79839016153122,0.0,0.0,0.0,40.81843103464889))
targetgene="D830030K20Rik_Gm10408_Gm3043_Gm8050_Gm8271"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(10.121234109619449,128.6776746553163,11.74698315948159,19.513147402207995,104.23006571941508,15.726679004509197,22.359476802497646,32.336134704369016,30.46628523648061,4.008795529147334,87.26699048787003),c(32.773519974005836,3.3862545961925345,92.09634797033567,6.886993200779292,231.26045831495222,53.29596773750339,0.8599798770191402,0.0,90.12942715792181,0.0,92.8971189064423),c(45.78653525780227,16.931272980962675,115.59031428929886,32.139301603636696,0.0,27.084836063321397,0.0,104.03799861405683,0.0,0.0,0.0),c(104.58608579940098,47.40756434669548,89.2770720120601,51.65244900584469,221.48888965375704,35.82187995471539,474.7088921145654,231.97661853134295,15.233142618240304,88.19350164124135,56.30128418572261))
targetgene="Tmem169"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(6.265525877383468,2.1164091226203343,14.566259117757173,82.6439184093515,57.00081719030512,3.4948175565575994,9.459778647210543,30.93021580417906,0.0,0.0,8.44519262785839),c(11.085161167678445,3.3862545961925345,181.37341998239577,21.80881180246776,0.0,0.0,0.8599798770191402,0.0,12.694285515200255,0.0,2.81506420928613),c(0.9639270580589951,8.888918315005403,0.0,4.591328800519528,13.028758214926885,7.863339502254599,0.0,0.0,130.75114080656263,0.0,0.0),c(7.229452935442463,0.8465636490481336,42.28913937413373,4.591328800519528,6.5143791074634425,100.47600475103098,0.0,7.029594500949786,0.0,10.021988822868336,0.0))
targetgene="Terf2"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(32.29155644497634,84.65636490481336,81.28912346361261,8.034825400909174,4.885784330597582,0.8737043891393999,15.479637786344524,0.0,0.0,22.04837541031034,0.0),c(96.39270580589951,81.6933921331449,71.89153693602734,10.330489801168937,0.0,96.98118719447338,6.019859139133982,0.0,90.12942715792181,154.33862787217237,60.5238804996518),c(99.766450509106,52.06366441646022,192.6505238154981,2.295664400259764,14.657352991792745,11.358157058812198,5.159879262114842,0.0,0.0,0.0,178.75657728966928))
targetgene="Rrp15"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(15.422832928943922,11.428609262149804,42.28913937413373,3.443496600389646,60.25800674403684,0.0,9.459778647210543,0.0,31.735713788000634,26.05717093945767,22.52051367428904),c(0.0,0.0,0.0,1.147832200129882,26.05751642985377,0.0,0.0,0.0,0.0,0.0,0.0),c(9.157307051560453,0.4232818245240668,7.518069222068219,20.660979602337875,166.11666724031778,415.8832892303543,23.219456679516785,0.0,0.0,10.021988822868336,0.0),c(20.242468219238898,18.201118454534875,0.9397586527585273,42.46979140480563,0.0,79.50709941168539,23.219456679516785,0.0,0.0,0.0,1.407532104643065))
targetgene="Psmb4"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(23.61621292244538,3.8095364207166016,8.92770720120601,22.95664400259764,30.94330076045135,0.8737043891393999,0.0,1.4059189001899572,0.0,0.0,226.61266884753348),c(127.23837166378736,44.86787339955108,66.72286434585544,259.4100772293533,39.08627464478066,27.084836063321397,18.919557294421086,4.217756700569872,151.061997630883,4.008795529147334,8.44519262785839),c(4.3376717612654785,15.238145682866406,32.42167352016919,50.504616805714804,175.88823590151296,34.074471176436596,34.39919508076561,2.8118378003799145,2.538857103040051,80.17591058294668,0.0),c(29.88173879982885,2.9629727716684675,37.590346110341095,9.182657601039056,21.17173209925619,20.968905339345596,0.0,8.435513401139744,8.885999860640178,0.0,35.188302616076626))
targetgene="Polg2"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(0.0,18.201118454534875,7.518069222068219,2.295664400259764,4.885784330597582,26.211131674181996,8.599798770191402,28.118378003799144,0.0,8.017591058294668,32.3732384067905),c(131.57604342505283,0.0,125.92765946964266,257.11441282909357,52.11503285970754,13.105565837090998,1.7199597540382805,101.22616081367693,17.771999721280356,12.026386587442001,12.667788941787586),c(16.868723516032414,107.0903016045889,6.108431242930427,321.39301603636693,8.142973884329303,55.91708090492159,3.439919508076561,5.623675600759829,7.616571309120152,220.48375410310337,0.0),c(40.484936438477796,0.8465636490481336,28.662638909135083,0.0,19.54313732239033,7.863339502254599,0.0,0.0,0.0,112.24627481612535,115.41763258073134))
targetgene="Med31"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(69.40274818024766,154.49786595128438,49.33732926982268,28.69580500324705,325.71895537317215,104.84452669672798,70.5183499155695,36.55389140493889,63.47142757600127,639.4028868989998,84.45192627858391),c(86.27147169628006,92.69871957077063,157.40957433705333,99.86140141129972,26.05751642985377,145.03492859714038,903.8388507471165,0.0,0.0,0.0,112.60256837144522),c(44.34064467071378,62.22242820503782,213.3252141761857,40.174127004545866,57.00081719030512,55.91708090492159,23.219456679516785,50.61308040683846,24.119142478880484,308.67725574434473,42.225963139291956),c(197.1230833730645,17.35455480548674,3.7590346110341093,228.41860782584652,86.31552317389061,281.33281330288673,63.63851089941638,170.11618692298484,57.124284818401144,1038.2780420491595,61.93141260429486))
targetgene="Tshz3"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(90.60914345754554,157.46083872295284,119.81922822671224,14.921818601688466,17.914542545524466,46.30633262438819,0.0,0.0,44.42999930320089,0.0,14.075321046430652))
targetgene="Gm20737_Gm20772_Gm20773_Gm20777_Gm20793_Gm20795_Gm20807_Gm20827_Gm20831_Gm20833_Gm20836_Gm20851_Gm21198_Gm21310_Gm21425_Gm21440_Gm21660_Gm21794_Gm28646_Gm29109"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(10.603197638648947,9.735481964053537,1.8795173055170546,0.0,0.0,31.453358009018395,0.8599798770191402,0.0,0.0,26.05717093945767,0.0),c(127.23837166378736,29.629727716684677,0.0,0.0,19.54313732239033,131.05565837090998,2.579939631057421,97.00840411310705,0.0,0.0,1.407532104643065))
targetgene="Eif3j1_Eif3j2"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetgenelist=c("Dnlz","Gm17175","Gm40814","Fam168a","4930578I06Rik","Atp2a2","Armc6","Tmem182","Gen1","C77080")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='7,8,9,10,11_vs_1,2,3,4,5,6 pos.'


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
targetmat=list(c(67.95685759315916,12.698454735722004,9.397586527585274,41.32195920467575,115.63022915747611,13.979270226230398,475.56887199158456,196.82864602659401,27.92742813344056,282.62008480488703,302.619402498259),c(66.51096700607067,7.619072841433203,0.0,66.57426760753316,0.0,22.716314117624396,121.25716265969878,52.018999307028416,182.79771141888367,208.45736751566136,49.26362366250728),c(4.8196352902949755,3.3862545961925345,29.60239756189361,3.443496600389646,0.0,8.737043891393999,0.8599798770191402,94.19656631272713,41.891142200160836,20.04397764573667,14.075321046430652),c(106.03197638648946,5.079381894288802,23.024086992583918,43.61762360493552,9.771568661195165,123.19231886865538,354.31170933188577,8.435513401139744,192.95313983104387,10.021988822868336,374.4035398350553))
targetgene="Dnlz"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(85.78950816725057,95.6616923424391,40.87950139499594,123.96587761402725,105.85866049628095,198.33089633464377,70.5183499155695,790.1264219067559,481.1134210260896,1074.3572018114855,23.928045778932105),c(48.678316431979255,22.010654875251475,350.99985680530995,268.5927348303924,71.65817018209786,329.38655470555375,110.9374041354691,29.524296903989104,224.6888536190445,32.07036423317867,360.32821878862467),c(65.54703994801167,38.94192785621415,74.24093356792366,32.139301603636696,255.6893799679401,46.30633262438819,19.779537171440225,582.0504246786423,431.60570751680865,394.8663596210124,88.6745225925131))
targetgene="Gm17175"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(51.088134077126746,19.047682103583007,179.02402335049945,104.45273021181926,91.2013075044882,64.65412479631559,705.183499155695,186.98721372526433,101.55428412160204,72.15831952465201,147.79087098752183),c(44.82260819974327,1.2698454735722005,0.0,0.0,6.5143791074634425,0.8737043891393999,142.7566595851773,0.0,600.4397048689721,194.4265831636457,14.075321046430652),c(53.97991525130373,134.18033837412918,77.5300888525785,233.00993662636603,620.4946099858929,124.93972764693417,202.09527109949795,126.53270101709616,30.46628523648061,807.7722991231878,157.64359572002328))
targetgene="Gm40814"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(116.63517402513841,305.6094773063762,187.01197189894694,172.1748300194823,30.94330076045135,110.96045742070378,203.81523085353624,123.72086321671624,81.24342729728163,76.16711505379935,38.00336682536276),c(122.41873637349238,61.799146380513754,74.71081289430292,25.252308402857402,37.45767986791479,47.180037013527595,1217.7315058591025,298.05480684027094,12.694285515200255,106.23308152240435,356.10562247469545),c(108.44179403163696,507.5149076043561,109.01200371998917,10.330489801168937,96.08709183508577,348.6080512666205,67.07843040749295,77.32553951044765,210.7251395523242,206.4529697510877,242.0955219986072),c(57.35365995451021,75.34416476528389,135.7951253236072,121.67021321376748,138.43055603359815,82.12821257910359,208.11513023863193,486.4479394657252,629.6365615539327,595.3061360783792,8.44519262785839))
targetgene="Fam168a"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(140.73335047661328,105.8204561310167,109.9517623727477,120.52238101363761,190.5455888933057,136.29788470574638,85.99798770191403,111.06759311500662,92.66828426096185,142.31224128473036,291.3591456611145),c(16.386759987002918,80.00026483504863,52.62648455447753,58.53944220662398,136.80196125673228,16.6003833936486,78.25816880874176,35.14797250474893,512.8491348140902,1980.3449913987831,608.0538692058042),c(178.8084692699436,521.4832078136503,553.0479671483934,343.20182783883473,281.7468963977939,125.81343203607358,143.61663946219642,216.51151062925342,168.83399735216338,563.2357718452005,508.1190897761465),c(61.20936818674619,336.5090504966331,70.01201963051028,55.09594560623434,45.6006537522441,50.67485457008519,167.69607601873236,309.3021580417906,60.93257047296122,809.7766968877614,333.58510880040643))
targetgene="4930578I06Rik"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(75.18631052860162,1504.7668861830575,143.31319454567543,275.47972803117165,166.11666724031778,196.58348755636496,1117.1138602478632,125.1267821169062,510.3102777110502,60.13193293721001,18.297917360359847),c(350.86944913347423,179.47149359820432,725.4936799295831,204.31413162311898,342.0049031418307,511.99077203568834,226.1747076560339,2004.840351670879,2116.1373953838824,990.1724956993916,775.5501896583289),c(442.44251964907875,359.7895508454568,389.5299615684096,51.65244900584469,78.17254928956132,519.8541115379429,461.80919395927833,125.1267821169062,2041.2411108442009,3682.078693521826,2088.777643290309),c(160.97581869585218,69.41821922194696,33.36143217292772,165.287836818703,92.82990228135405,15.726679004509197,61.05857126835896,94.19656631272713,71.08799888512142,30.065966468605005,128.0854215225189))
targetgene="Atp2a2"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(157.12011046361621,115.97921991959431,42.759018700512996,125.11370981415713,86.31552317389061,92.61266524877638,539.207382891001,205.26415942773374,312.27942367392626,759.6667527734198,130.90048573180505),c(95.91074227687001,190.053539211306,108.07224506723064,25.252308402857402,154.71650380225677,19.221496561066797,91.15786696402887,278.3719422376115,48.23828495776097,1172.5726922755953,522.1944108225772),c(203.87057277947747,59.68273725789342,209.09630023877233,66.57426760753316,52.11503285970754,138.91899787316459,16.339617663363665,12.653270101709616,429.0668504137686,10.021988822868336,346.252897742194),c(145.55298576690828,291.217895272558,158.8192123161911,142.33119281610536,399.00572033213587,79.50709941168539,448.9094958039912,59.04859380797821,883.5222718579377,372.8179842107021,88.6745225925131))
targetgene="Armc6"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(320.02378327558637,115.55593809507025,183.7228166142921,410.92392764649776,368.06241957168453,435.10478579142114,397.3107031828428,291.0252123393212,340.20685180736683,280.6156870403134,216.75994411503203),c(253.9947797985452,97.77810146505944,63.90358838757986,83.79175060948138,48.85784330597582,350.35546004489936,63.63851089941638,70.29594500949786,53.31599916384107,252.55411833628204,209.7222835918167),c(85.30754463822107,121.05860181388312,228.36135262032215,226.12294342558675,47.22924852910996,33.2007667872972,0.8599798770191402,175.73986252374465,106.63199832768214,6.013193293721001,9.852724732501455),c(0.9639270580589951,0.0,109.9517623727477,0.0,0.0,0.0,85.99798770191403,0.0,50.77714206080102,0.0,61.93141260429486))
targetgene="Tmem182"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(79.0420187608376,24.550345822395876,7.987948548447482,19.513147402207995,55.37222241343926,160.7616076016496,73.95826942364606,108.25575531462671,143.44542632176288,362.7959953878337,271.6536961961116),c(181.7002504441206,98.2013832895835,202.04811034308338,160.69650801818347,501.6071912746851,226.28943678710456,186.61563331315344,80.13737731082756,219.6111394129644,134.2946502264357,188.60930202217074),c(79.0420187608376,173.5455480548674,204.86738630135895,101.00923361142961,109.11585005001267,375.6928873299419,749.9024527606903,184.1753759248844,260.2328530616052,817.7942879460561,691.0982633797449),c(344.12195972706127,187.09056643963754,338.7829943194491,493.56784605584926,348.5192822492942,331.13396348383253,270.8936612610292,778.8790707052364,335.1291376012867,1713.7600887104852,1183.7345000048178))
targetgene="Gen1"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
targetmat=list(c(296.8895338821705,193.43979380749855,293.20469966066054,137.73986401558582,149.83071947165917,122.31861447951599,506.5281475642736,459.735480362116,413.8337077955283,761.6711505379935,121.0477609993036),c(4.8196352902949755,96.08497416696316,34.30119082568625,141.1833606159755,148.20212469479333,83.87562135738239,13.759678032306244,26.712459103609188,7.616571309120152,10.021988822868336,150.60593519680796),c(146.51691282496725,72.38119199361543,70.01201963051028,167.58350121896277,0.0,441.22071651539693,105.77752487335425,337.42053604558976,802.2788445606561,78.17151281837302,280.09888882396996),c(422.20005142983985,137.14331114579764,78.46984750533703,215.79245362441782,112.37303960374439,244.63722895903197,94.59778647210543,463.9532370626859,1186.9156956712238,294.64647139232903,567.2354381711552))
targetgene="C77080"
collabel=c("In.vivo_AR1890_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1891_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1892_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1893_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1894_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1895_Pool.E..33.40._BM_hEGFRv3_15m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1903_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1904_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1906_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1907_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438","In.vivo_AR1909_Pool.E..33.40._BM_mCD19_10m_1431..1432..1433..1434..1435..1436..1437..1438")

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
Sweave("invivo.poolE_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.Rnw");
library(tools);

texi2dvi("invivo.poolE_BM_mCD19_10m_v_BM_hEGFRv3_15m_summary.tex",pdf=TRUE);

