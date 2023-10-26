pdf(file='cd19_1to10_v_NOCAR_T.pdf',width=4.5,height=4.5);
gstable=read.table('cd19_1to10_v_NOCAR_T.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Ptpn2","Stub1","Chic2","Socs1","Usp22","Fitm2","Kdsr","Ccdc71l","Kat6a","Atg12")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3_vs_NoCAR_T_rep1,NoCAR_T_rep2,NoCAR_T_rep3 neg.'


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
targetmat=list(c(503.92274697945,560.0262738875184,538.8022285242245,151.9240345966206,164.2878853338961,165.31164347678344),c(787.3169254787842,813.9449966358018,871.2754729856334,452.11128367909987,452.3393109526606,435.3548871686719),c(522.8822166677857,576.0843749308881,574.9831404214955,141.85677929202527,152.24010707607707,133.48145125454565),c(263.4368419852966,290.04945009586527,231.75340863927622,48.50586646759574,47.09586046238355,43.1247765591609),c(1001.8582930046887,1003.6313152106064,996.4418709005167,592.1376529157442,610.0556808732008,598.6129698569239),c(753.3894534049202,708.5637085386882,755.8876999078503,88.7748876859771,117.19202487151256,97.54413745524488))
targetgene="Ptpn2"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(646.617703054819,728.6363348429003,764.6884622612405,384.38611163000394,410.7197133347403,419.9531812546859),c(921.0310801228363,966.496956547814,1039.467820183758,278.22232841790765,295.71819360101296,312.1412398567836),c(991.8796247476698,1020.6930475691868,1079.560182015869,431.0615680422187,478.62537260608394,497.9884912188818),c(311.3344496189869,411.48883923634867,378.43278119578014,180.29539045502565,142.3828339560433,161.2045218997205),c(864.152671057829,968.5042191782352,910.3899723340345,474.991409371362,418.38648131698875,437.4084479572034),c(553.8160882645441,558.0190112570972,620.9426771558667,314.830529525527,349.3855694767524,335.75718892489556))
targetgene="Stub1"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(651.6070371833283,634.2949912131032,647.3449642160374,383.47090660231345,377.86213626796103,357.319577204476),c(696.511044339913,720.6072843212154,785.2235744191511,256.257407753336,278.1941524987307,243.34695344097935),c(618.6774319351663,677.4511377671594,652.2342766345876,412.757467488409,371.2906208546052,381.96230666685364),c(487.95687776821984,520.8846525943047,518.2671163663139,210.4971563688117,211.38374579627964,224.8649063441961),c(420.10193362049193,429.55420291013957,418.52514302789126,286.459173667122,266.14637424091165,275.1771456632172),c(641.6283689263096,644.3313043652093,628.7655770255469,307.50888930400316,296.8134461699056,316.2483614338466))
targetgene="Chic2"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(276.40911071942105,259.9405106395471,271.8457704713873,137.28075415357284,144.57333909382857,104.73160021510503),c(405.1339312349637,365.32179873666075,352.03049413560944,205.00592620266877,205.90748295181643,210.4899808244758),c(254.45604055397968,220.79888934633343,200.4618091605554,107.07898823978682,107.33475175147879,87.27633351258753),c(495.9398123738349,446.61593526871985,507.5106290455036,482.31304959288593,478.62537260608394,461.0243970253153),c(365.21925820688847,387.4016876712941,324.6503445917287,22.88012569226214,17.52404110228225,36.96409419356648),c(512.9035484107669,562.0335365179396,523.156428784864,167.48252006735888,177.43091616060778,172.4991062366436))
targetgene="Socs1"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(368.2128586839941,333.20559664992135,359.85339400528966,393.53816190690884,360.3380951656788,406.60503612923134),c(737.4235841936902,657.3785114629472,723.6182379454194,535.3949411989341,509.2924445350779,531.8722442296511),c(629.653967017887,697.5237640713715,618.0090897047367,379.81008649155154,328.5757706677922,362.45347917580466),c(444.0507374373371,464.6812989425108,428.3037678649915,355.0995507439084,305.5754667210467,303.92699670265773),c(440.0592701345296,396.43436950818955,423.41445544644137,271.8158932240742,279.28940506762336,275.1771456632172),c(306.3451154904775,250.9078288026516,296.29233256413795,184.8714155934781,163.19263276500345,173.52588663090933))
targetgene="Usp22"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(278.40484437082483,244.88604091138797,265.0007330854171,201.34510609190684,224.52677662299135,162.23130229398623),c(386.174461546628,407.4743139755062,442.9717051206419,321.23696471936046,372.3858734234978,302.90021630839203),c(374.2000596382054,433.568728170982,399.9457558374007,244.35974239335965,254.09859598309262,292.63241236573464),c(348.2555221699565,349.263697693291,367.6762938749699,201.34510609190684,194.95495726289005,202.2757376703499),c(274.4133770680173,260.9441419547577,243.48775844379654,133.6199340428109,134.7160659737948,119.10652573482534),c(187.59896323195366,196.71173778127886,190.68318432345512,108.90939829516779,101.85848890701558,129.37432967748268))
targetgene="Fitm2"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(288.38351262784363,219.7952580311228,272.8236329550973,169.31293012273983,150.04960193829177,152.99027874559462),c(489.9526114196236,504.82655155093505,517.2892538826039,323.9825798024319,328.5757706677922,315.2215810395808),c(256.45177420538346,271.98408642207437,290.4251576618778,195.85387592576393,235.47930231191773,195.08827491048976),c(513.9014152364688,531.9245970616214,511.4220789803437,442.04402837450453,467.6728469171576,503.12239319021046),c(312.3323164446888,318.1511269217622,307.0488198849483,213.24277145188316,217.9552612096355,224.8649063441961),c(504.92061380515185,491.77934445319715,508.48849152921366,396.28377698998025,441.38678526373417,372.72128311846205))
targetgene="Kdsr"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(437.06566965742394,482.7466626163017,428.3037678649915,333.1346300793368,337.3377912189333,360.3999183872732),c(609.6966305038494,584.113425452573,614.0976397698965,480.48263953750495,536.6737587573939,548.3007305379028),c(796.2977269101011,806.9195774293275,729.4854128476796,587.5616277772917,658.246793904477,622.2289189250358),c(329.2960524816208,394.4271068777683,385.27781858175035,304.7632742209317,296.8134461699056,330.6232869535669),c(318.3195173989001,300.08576324797133,311.9381323034984,266.32466305793133,199.3359675384606,176.60622781370654),c(436.0678028317221,447.61956658393046,487.95337937130313,351.43873063314646,453.43456352155323,364.5070399643362))
targetgene="Ccdc71l"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(293.37284675635306,312.1293390304986,336.384694396249,125.38308879359653,101.85848890701558,139.64213362014004),c(355.2405899498697,384.39079372566226,369.6320188423899,344.1170904116226,357.05233745900085,307.00733788545494),c(662.583572266049,689.4947135496866,705.038850754929,448.45046356833797,419.4817338858814,461.0243970253153),c(241.48377181985524,228.82793986801826,215.1297464162058,163.82169995659692,115.00151973372726,98.57091784951062),c(161.65442576370475,182.66089936833038,213.17402144878574,120.8070636551441,113.90626716483463,145.80281598573447),c(408.12753171206936,438.586884747035,441.01598015322185,604.0353182757206,610.0556808732008,585.2648247314693))
targetgene="Kat6a"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(538.8480858790158,569.0589557244139,565.2045155843952,448.45046356833797,408.529208196955,426.1138636202803),c(338.2768539129377,341.2346471716062,353.9862191030295,237.95330719952625,263.9558691031264,246.42729462377656),c(495.9398123738349,516.8701273334623,497.7320042084034,365.16680604850376,334.0520335122554,324.46260458797246),c(805.2785283414181,799.8941582228533,832.1609736372324,800.8043992291749,847.7254883229039,869.682993943078),c(487.95687776821984,498.8047636596714,516.3113913988939,308.42409433169365,378.9573888368537,329.59650655930113),c(466.00380760277847,451.6340918447729,462.52895479484243,369.7428311869562,415.1007236103108,365.5338203586019))
targetgene="Atg12"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetgenelist=c("Jak2","Jak1","Stat1","Tsc1","Irf1","Nprl3","Tsc2","Ifngr1","Nprl2","Cnot8")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3_vs_NoCAR_T_rep1,NoCAR_T_rep2,NoCAR_T_rep3 pos.'


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
targetmat=list(c(246.47310594836463,243.88240959617735,217.08547138362584,897.8161321643664,892.6308436475022,857.3616292118893),c(514.8992820621706,495.7938697140396,481.1083419853329,1613.5064638183262,1658.212389303458,1637.7147288538483),c(454.0294056943559,456.6522484208259,433.19308028354163,1448.7695588340387,1451.2096537827488,1485.7512305025195),c(936.9969493340664,913.3044968416518,915.2792847525845,2545.1851820072407,2636.272933324586,2578.245570001262),c(401.1424639321562,409.48157660592744,430.2594928324116,1382.8747968403238,1302.2553044133497,1379.9928498931488),c(791.3083927815917,773.7997440273775,784.245711935441,2714.49811212998,2993.325270783587,2908.868856954829))
targetgene="Jak2"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(592.7328944669174,599.1678951807321,600.4075649979561,1996.062165392949,1928.7397738199402,1868.7403175636389),c(679.547308302981,630.2804659522609,656.1457265694277,2472.883984819692,2566.176768915457,2555.656401327416),c(649.6113035319246,566.0480617787821,599.4297025142461,1743.4655777503751,1732.6895639881575,1638.741509248114),c(371.20645916109976,384.39079372566226,397.99003086998067,1513.7491158000632,1377.8277316669419,1437.4925519720298),c(1073.7047044552241,1173.2450074811989,1084.4494944344192,2573.5565378656456,2596.843840844451,2711.727021255808),c(368.2128586839941,427.54694027971834,404.8350682559509,914.2898226627951,1044.8709507235792,993.9234216492321))
targetgene="Jak1"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(674.5579741744716,653.3639862021048,692.3266384666986,1740.7199626673037,1803.880980966179,1711.6429172409812),c(719.4619813310563,717.5963903755836,670.813663825078,2260.5564183954993,2332.8879717413247,2299.9880831552478),c(540.8438195304195,518.8773899638835,517.2892538826039,866.6991612228899,856.487508874045,844.0134840864347),c(564.7926233472647,628.2732033218397,543.6915409427746,1947.5562989253533,1855.3578517041333,1841.017246918464),c(446.0464710887409,443.6050413230881,480.1304795016229,1380.1291817572524,1394.2565202003316,1445.7067951261558),c(473.98674220839354,482.7466626163017,552.4923032961649,770.6026333153889,769.9625559315264,801.9154879215396))
targetgene="Stat1"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(742.4129183221995,669.4220872454745,729.4854128476796,990.2518399611055,1081.0142854970363,1061.6909276707706),c(928.0161479027495,1013.6676283627125,1012.0876706398772,1620.8281040398501,1531.1630913119116,1530.929567850212),c(739.4193178450939,738.6726479950064,770.5556371635007,1244.6788376590605,1315.3983352400614,1404.6355793555263),c(542.8395531818234,593.1461072894684,564.2266531006851,1097.3308282008923,1135.7769139416685,1236.2435946959456),c(736.4257173679882,659.3857740933685,753.9319749404302,1183.3601008037979,1239.8259079864692,1171.5564298572044),c(400.14459710645434,437.5832534318244,412.6579681256311,787.0763238138177,734.9144737269619,678.7018406096513))
targetgene="Tsc1"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(611.6923641552531,604.1860517567851,573.0274154540754,1459.7520191663245,1366.8752059780156,1422.0908460580438),c(847.188935020897,837.0285168856458,828.2495237023923,1685.8076610058745,1634.11683278782,1668.5181406818203),c(563.7947565215628,652.3603548868942,565.2045155843952,1303.2519594312514,1408.494803595936,1351.242998853708),c(614.6859646323587,612.2151022784699,637.5663393789371,1579.643877793778,1467.6384423161385,1468.295963800002),c(795.2998600843993,767.777956136114,819.448761349002,843.8190355306277,947.3934720921342,884.0579194627984),c(568.7840906500722,673.4366125063169,641.4777893137772,1321.5560599850612,1302.2553044133497,1478.5637677426594))
targetgene="Irf1"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1253.3207330815626,1367.9494826320565,1305.446415752885,1875.255101737805,1990.0739176779282,1784.5443252338484),c(876.1270729662516,960.4751686565504,922.1243221385548,1335.2841354004186,1407.3995510270433,1394.367775412869),c(1146.5489827314614,1151.1651185465655,1229.173142023503,2167.20550557107,2132.456751633971,2084.364200359443),c(1232.3655297418231,1289.6662400456294,1177.3464303868716,1540.2900616030872,1565.1159209475836,1516.5546423304916),c(1365.0818175601735,1243.4991995459413,1256.5532915673837,1745.295987805756,1846.595831152992,1831.7762233700723),c(1001.8582930046887,944.4170676131806,1014.0433956072973,1297.7607292651087,1317.5888403778467,1278.3415908608408))
targetgene="Nprl3"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(695.5131775142112,605.1896830719957,620.9426771558667,1051.570576816368,1033.9184250346527,1074.0122924019595),c(940.9884166368739,886.2064513309655,971.0174463240561,1457.9216091109436,1365.779953409123,1370.751826344757),c(1095.6577746206656,1322.7860734475794,1162.6784931312213,2460.0711144320253,2615.463134515626,2743.5572134780455),c(381.1851274181186,384.39079372566226,372.56560629352003,458.5177188729333,463.291836641587,475.3993225450356),c(1520.7490423696668,1664.0207206191856,1590.0043985125028,2773.986438929862,2689.9403092003254,2749.7178958436402),c(959.9478863252097,1033.7402546669246,1046.312857569728,1655.6058950920885,1623.1643070988935,1674.6788230474149))
targetgene="Tsc2"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(709.4833130740375,725.6254408972685,666.9022138902379,848.3950606690802,876.2020551141126,921.0220136563648),c(543.8374200075252,513.8592333878305,558.359478198425,613.1873685526253,663.7230567489403,768.0317349107703),c(163.6501594151085,155.562853857644,176.01524706780472,507.02358534052905,460.00607893490906,512.363416738602),c(400.14459710645434,356.2891168997653,406.79079322337094,1228.2051471606317,1133.5864088038832,1274.234469283778),c(408.12753171206936,434.57235948619257,485.9976544038831,1045.1641416225345,997.7750902611956,1029.8607354485328),c(913.0481455172212,891.2246079070185,892.788447627254,2859.100506505077,3189.3754806153697,3162.483614338466))
targetgene="Ifngr1"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1008.8433607846018,1065.8564567536641,970.039583840346,1221.7987119667982,1221.2066143152942,1279.3683712551067),c(404.13606440926185,428.55057159492895,472.3075796319427,570.1727322511725,566.2455781174953,547.2739501436371),c(1386.037020899913,1426.1600989142717,1427.6792262166384,2565.319692616431,2432.555955510555,2457.085483477905),c(501.9270133280462,549.9899607354123,505.55490407808355,646.1347495494829,658.246793904477,639.6841856275533),c(1152.5361836856728,1166.2195882747246,1130.4090311687903,1612.5912587906357,1527.8773336052336,1578.1614659864356),c(945.9777507653833,895.2391331678609,995.4640084168068,1443.2783286678957,1430.3998549737887,1529.902787455946))
targetgene="Nprl2"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
targetmat=list(c(1021.8156295187263,905.275446319967,942.6594342964653,1141.2606695300356,1100.728831737104,1138.699457240701),c(436.0678028317221,387.4016876712941,398.9678933536907,684.5733607124832,674.6755824378666,659.1930131186023),c(351.24912264706217,298.0785006175501,296.29233256413795,574.748757389625,536.6737587573939,596.5594090683924),c(581.7563593841967,602.1787891263639,482.08620446904297,1107.3980835054876,973.6795337455576,1076.065853190491),c(333.2875197844283,341.2346471716062,336.384694396249,431.9767730699092,348.29031690785973,394.2836713980425),c(391.1637956751374,417.5106271276123,387.2335435491704,781.5850936476747,738.2002314336398,695.1303269179031))
targetgene="Cnot8"
collabel=c("NoCAR_T_rep1","NoCAR_T_rep2","NoCAR_T_rep3","amCD19_1to10rep1","amCD19_1to10rep2","amCD19_1to10rep3")

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
Sweave("cd19_1to10_v_NOCAR_T_summary.Rnw");
library(tools);

texi2dvi("cd19_1to10_v_NOCAR_T_summary.tex",pdf=TRUE);

