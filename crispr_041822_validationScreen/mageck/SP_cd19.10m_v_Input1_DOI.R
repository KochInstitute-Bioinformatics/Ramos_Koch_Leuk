pdf(file='SP_cd19.10m_v_Input1_DOI.pdf',width=4.5,height=4.5);
gstable=read.table('SP_cd19.10m_v_Input1_DOI.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("H2-K1","H2-T23","Sec63","Sptssa","Elob","B4galt1","Gclc","Xrcc3","Cnot9","Rps10")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2473_SP,AR2474_SP,AR2475_SP,AR2476_SP_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 neg.'


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
targetmat=list(c(608.3052574305932,588.3373357749003,585.3834793721384,105.93779778294775,150.81302945200017,422.82153362296776,105.39561167783715),c(515.2155770704082,499.31933888374147,499.38341641657024,6.485987619364149,243.2029934406129,86.48622278651614,69.92593467087272),c(459.3617688542973,429.40098883513014,493.26137803668234,90.80382667109808,154.8890572750272,48.047901548064516,322.26735109184824),c(543.1424811784638,521.1475067037958,493.84442931095737,38.915925716184894,19.021463174126147,0.0,11.147612773617391),c(361.6176044761031,341.065122188348,327.383290505434,34.591933969942126,32.60822258421625,57.657481857677425,55.738063868086954),c(390.6186202806222,374.4895041628061,392.6850332242382,15.13397111184968,17.662787233117136,28.828740928838712,18.241548175010276))
targetgene="H2-K1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(481.91811448003443,516.3725949931588,516.0003777334089,69.18386793988425,423.9068935948113,19.219160619225807,289.8379321140522),c(796.2748043117359,768.0786551681597,757.3836052832746,181.60765334219616,158.96508509805423,403.60237300374195,139.85186934174544),c(654.4920603785312,577.0821867426848,646.0208118967423,56.21189270115595,236.40961373556783,192.19160619225806,87.15406350282687),c(640.886645556658,600.2746150514924,599.0851843176018,84.31783905173393,44.83630605329735,67.26706216729032,63.845418612535966),c(361.9756417082576,402.797909304439,427.37658404360303,192.41763270780308,364.12515219041484,365.1640517652903,290.8513514571083),c(382.741801173222,364.25755049715565,359.1595849534235,36.75392984306351,58.42306546338745,115.31496371535485,35.469677006964424))
targetgene="H2-T23"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(495.16549206975304,448.15957055548927,522.1224161132967,250.79152128208042,220.1055024434597,9.609580309612904,97.28825693338814),c(478.3377421584888,439.6329425007806,418.6308149294774,510.23102605664633,211.95344679740566,76.87664247690323,118.57006313756679),c(382.741801173222,381.9929368509498,353.9121234849482,140.5297317528899,52.988361699351415,0.0,56.75148321114308),c(335.12284929666583,301.50156801449964,303.76971389729493,25.943950477456596,44.83630605329735,19.219160619225807,32.42941897779605),c(667.7394379682497,628.241955070937,667.0106577706438,252.9535171552018,339.66898525225264,999.396352199742,351.6565120404759),c(361.6176044761031,378.92335075125465,360.3256875019736,131.88174826040435,129.07421439585602,96.09580309612903,47.63070912363794))
targetgene="Sec63"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(415.3231892992867,404.1621697931924,405.8036868954266,32.42993809682074,38.042926348252294,0.0,34.4562576639083),c(438.59560938933294,419.51010029166804,416.5901354695148,164.31168635722508,187.49727985924346,57.657481857677425,240.18038430430198),c(650.9116880569856,596.5228987074206,578.3868640808379,239.9815419164735,322.0061980191355,0.0,172.2812883195415),c(511.2771675167081,485.33566887401923,528.5359801303222,298.35543049075085,206.51874303336962,38.438321238451614,199.6436105820569),c(537.0558482318363,467.94134764241346,512.5020700877586,30.26794222369936,187.49727985924346,182.58202588264515,40.53677372224506),c(482.27615171218895,422.57968639136317,436.9969300691412,125.39576064104021,199.72536332832456,76.87664247690323,163.16051423203635))
targetgene="Sptssa"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(289.2940835808825,262.2790789628396,270.53579126361774,47.563909208670424,17.662787233117136,0.0,100.32851496255653),c(461.50999224722466,463.16643593177656,439.6206608033788,307.0034139832364,6.793379705045053,38.438321238451614,35.469677006964424),c(624.0588956453937,607.4369826174478,687.4174523702701,588.0628774890162,248.63769720464893,153.75328495380646,265.5158678807051),c(566.41490126851,509.21022742720356,554.7732874726988,369.70129430375647,453.7977642970095,297.896989598,169.24103029037312),c(426.06430626392347,448.15957055548927,409.59352017821436,32.42993809682074,233.6922618535498,67.26706216729032,168.227610947317),c(376.6551682265945,375.1716344071828,324.17650849692126,337.2713562069357,107.33539933971184,19.219160619225807,107.42245036394941))
targetgene="Elob"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(508.7709068916262,565.485972588281,518.6241084676465,518.8790095491319,73.36850081448657,259.45866835954837,199.6436105820569),c(361.6176044761031,340.3829919439713,385.6884179329378,77.83185143236979,70.65114893246854,57.657481857677425,32.42941897779605),c(202.64907339947956,179.40025427107105,168.21029262834853,15.13397111184968,104.61804745769382,0.0,128.70425656812807),c(530.2531408208997,486.0177991183959,509.8783393535209,155.66370286473958,29.89087070219823,96.09580309612903,46.617289780581814),c(524.8825823385813,578.7875123536265,596.1699279462266,105.93779778294775,279.88724384785615,67.26706216729032,157.0799981736996),c(480.1279283192616,488.064189851526,453.90541702311725,118.90977302167606,168.4758166851173,9.609580309612904,96.27483759033201))
targetgene="B4galt1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(440.38579555010574,446.45424494454755,416.2986098323773,123.23376476791883,88.31393616558569,672.6706216729033,81.07354744449012),c(355.8890087616302,369.71459245216926,365.2816233333115,142.69172762601127,218.7468265024507,153.75328495380646,88.167482845883),c(379.8775033159855,460.4379149542698,465.85796814575554,84.31783905173393,92.38996398861272,240.23950774032258,136.81161131257707),c(360.90153001179397,436.2222912788971,386.2714692072128,121.07176889479744,370.91853189545986,76.87664247690323,217.8851587570672),c(491.5851197482075,419.51010029166804,444.8681222718542,116.74777714855468,161.68243698007225,201.80118650187097,146.94580474313832),c(462.5841039436883,400.41045344912055,406.96978944397665,131.88174826040435,171.19316856713533,124.92454402496774,163.16051423203635))
targetgene="Gclc"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(250.62606250819027,202.59268257987873,239.92559936417823,30.26794222369936,14.945435351099116,0.0,21.281806204178658),c(369.49442358350336,396.6587371050487,380.1494308273249,60.53588444739872,66.57512110944151,259.45866835954837,166.20077226120475),c(547.796965196473,539.2239581797782,558.2715951183491,170.79767397658924,107.33539933971184,461.2598548614194,168.227610947317),c(204.43925956025234,193.3839242807933,227.09847133012738,19.457962858092447,4.076027823027031,105.70538340574194,1.0134193430561265),c(398.85347662017705,394.2712812497303,352.1629696621231,131.88174826040435,135.86759410090104,9.609580309612904,74.99303138615336),c(606.5150712698204,516.0315298709705,601.4173894147019,324.2993809682074,404.88543042068517,499.698176099871,208.76438466956205))
targetgene="Xrcc3"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(414.6071148349776,351.6381409761868,383.35621283583765,1701.4907521465284,133.15024221888302,182.58202588264515,83.10038613060237),c(534.9076248389089,526.6045486588093,496.1766344080575,235.65755017023073,104.61804745769382,903.300549103613,310.1063189751747),c(369.85246081565793,320.9422799792355,297.06462424313196,105.93779778294775,92.38996398861272,345.94489114606455,91.20774087505139),c(364.8399395654941,377.55909026250123,384.8138410215252,25.943950477456596,66.57512110944151,19.219160619225807,83.10038613060237),c(274.9725942947002,305.2532843585715,285.98665003190627,21.619958731213828,47.55365793531537,240.23950774032258,62.83199926947984),c(424.27412010315066,408.5960163816409,457.695250305905,38.915925716184894,258.14842879171204,19.219160619225807,149.98606277230672))
targetgene="Cnot9"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(275.68866875900926,282.0608560497638,251.586624849679,144.85372349913266,296.1913551399643,249.8490880499355,145.9323854000822),c(413.1749659063594,431.4473795682602,422.7121738494027,183.76964921531754,110.05275122172985,9.609580309612904,89.18090218893913),c(419.6196360851414,414.7351885810312,402.8884305240514,73.50785968612702,187.49727985924346,28.828740928838712,169.24103029037312),c(392.40880644139503,387.108913683775,366.4477258818615,92.96582254421946,141.3022978649371,76.87664247690323,162.14709488898023),c(283.9235250985641,272.51103262849006,243.71543264696598,19.457962858092447,99.18334369365778,0.0,139.85186934174544),c(295.7387537596645,244.88475773123386,273.4510476349929,58.37388857427734,33.96689852522526,115.31496371535485,31.41599963473992))
targetgene="Rps10"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetgenelist=c("Cd19","Arl15","Zfp36l2","Setd5","Sipa1","Sephs1","Tmem94","Ncor1","Nosip","Thrap3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2473_SP,AR2474_SP,AR2475_SP,AR2476_SP_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 pos.'


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
targetmat=list(c(427.4964551925417,377.55909026250123,452.44778883742964,50067.500429744985,4646.671718250816,798950.1165215264,1816.0474627565786),c(505.1905345700807,518.418985726289,522.1224161132967,21323.7652965962,5196.935474359466,405130.2962729704,1437.0286284535873),c(491.943156980362,552.8665630673121,493.84442931095737,47862.26463916117,6214.583754175214,784237.849067509,3933.080470400827),c(555.3157470717188,542.97567452385,543.6953132614732,114293.9118325619,18683.152864814903,3341789.210149745,8553.259255393707),c(533.4754759102907,470.66986861992024,495.5935831337825,36976.61541799501,7843.636207445018,442540.39241829346,3164.908608364283),c(372.00068420858526,338.67766633302955,387.4375717557629,47475.26737787245,9679.20740374819,2183911.6594838668,3084.848480262849))
targetgene="Cd19"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(559.2541566254189,561.0521259998325,497.6342625937451,2010.656162002886,4561.075133967249,1604.7999117053548,2226.4822966943098),c(767.9898629715259,685.5408955985795,709.5734007927215,1615.010917221673,4173.85249077968,5938.720631340774,4149.952209814838),c(597.206103233802,590.3837265080303,595.2953510348141,4047.256274483229,8641.178984817307,3680.469258581742,3933.080470400827),c(629.7874913598666,644.2720158137894,550.9834541899111,2408.4634026572203,3047.5101356832106,1258.8550205592903,2172.7710715123353),c(789.4720969007993,782.7444554222586,752.7191950890743,6185.470193000277,6459.145423556836,1585.580751086129,3437.518411646381),c(625.1330073418575,657.5735555791349,622.4072352886033,5024.478409134093,4827.375618405014,2661.8537457627745,3999.9661470425312))
targetgene="Arl15"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(491.5851197482075,399.04619296036714,442.2443915376165,1375.0293753051994,751.3477953779828,1249.2454402496774,1230.2910824701376),c(641.9607572531216,622.1027828715468,639.024196605442,2704.6568372748497,1860.0273632413355,1268.4646008689033,2079.5364919511717),c(657.7143954679221,692.3621980423464,624.1563891114284,1900.3943724736955,835.5857037205415,1921.9160619225806,2967.2918364683383),c(506.26464626654433,558.664670144514,530.8681852274223,1697.1667604002855,1171.178661149767,1335.7316630361936,1102.6002452450657),c(625.8490818061665,656.5503602125699,653.0174271880428,1450.6992308644478,1472.8047200537674,1220.4166993208387,1224.2105664118008),c(450.0528008182788,423.9439468801166,418.3392892923399,862.6363533754318,1679.323463087137,4564.550647066129,697.232508022615))
targetgene="Zfp36l2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(546.3648162678547,487.72312472933766,502.8817240622205,492.9350590716753,1388.5668117112089,855.2526475555484,1554.585272248098),c(532.401364213827,520.1243113372307,486.84781401965694,7947.4968295942035,3123.5959883797154,1095.492155295871,3456.7733791644473),c(394.9150670664769,446.45424494454755,446.61727609467926,665.8947289213859,1358.6759410090106,932.1292900324516,851.2722481671462),c(575.7238693045285,545.0220652569801,588.8817870177886,4326.153742115887,5729.536443234998,4458.845263660387,4758.003815648513),c(323.66565786772,306.2764797251365,347.7900851050603,884.2563121066456,2721.427909841048,566.9652382671613,1464.3909507161027),c(539.2040716247636,490.1105805846561,513.3766469991712,1148.0198086274543,1186.1240965008662,3305.695626506839,1154.284631740928))
targetgene="Setd5"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(482.27615171218895,446.45424494454755,432.62404551207834,307.0034139832364,205.1600670923606,884.0813884843872,486.4412846669407),c(659.8626188608495,622.4438479937351,605.4987483346272,1275.5775651416159,2141.2732830302007,1585.580751086129,1090.439213128392),c(363.40779063687586,364.25755049715565,381.02400773873745,1498.2631400731184,1508.1302945200018,787.9855853882581,1204.9555988937343),c(545.6487418035457,508.52809718282685,483.9325576482817,2544.6691426638677,1148.0811701526138,1402.998725203484,2395.723326984683),c(587.5390979656289,523.8760276813025,606.0817996089022,2146.861902009533,1940.1892437608672,365.1640517652903,1556.6121109342102),c(330.46836527865656,320.9422799792355,324.46803413405877,4395.337610055772,1547.5318968092631,432.4311139325807,526.9780583891858))
targetgene="Sipa1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(453.27513590766983,476.1269105749338,496.468160045195,633.4647908245652,706.5114893246855,2162.1555696629034,606.0247671475636),c(490.15297081958926,435.19909591233204,456.8206733944925,1558.7990245205172,1857.3100113593175,1575.9711707765161,1332.6464361188064),c(407.44637019188644,425.30820736886994,436.4138787948661,1537.1790657893032,1235.0364303771905,528.5269170287097,1200.9019215215098),c(491.2270825160529,456.68619861019795,517.7495315562339,4401.823597675136,3008.1085333939495,2104.498087805226,2146.422168592876),c(511.6352047488627,468.62347788679017,561.4783771268618,949.116188300287,922.5409639451182,970.5676112709033,1633.6319810064758),c(453.63317313982435,445.77211470017085,461.19355795155525,2689.5228661630003,1775.7894548987767,1825.8202588264517,968.828891961657))
targetgene="Sephs1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(521.3022100170357,555.5950840448189,548.3597234556735,1055.0539860832348,832.8683518385235,1441.4370464419355,1734.9739153120886),c(308.98613134938313,331.5152987670743,366.4477258818615,1245.3096229179166,917.1062601810821,269.0682486691613,369.89806021548617),c(541.7103322498456,502.3889249834366,614.8275687230278,1405.297317528899,1084.2234009251904,1835.4298391360646,1891.040494142732),c(503.75838564146244,455.66300324363294,515.1258008219962,425.9131870049124,1025.800335461803,1268.4646008689033,655.6823149573138),c(445.3983168002695,370.0556575743576,420.37996875230255,1297.1975238728296,646.729747920289,1787.381937588,1469.4580474313834),c(701.3949377907782,624.8313038490535,690.3327087416452,1954.4442693017302,1004.0615204056588,1739.3340360399357,1048.8890200630908))
targetgene="Tmem94"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(299.3191260812101,313.77991241328016,297.35614988026947,1876.6124178693603,115.4874549857659,547.7460776479355,692.1654113073344),c(655.5661720749948,607.4369826174478,584.2173768235883,2049.572087719071,2437.464638170165,1076.272994676645,1875.83920399689),c(443.9661678716513,499.6604040059298,452.15626320029213,778.3185143236979,800.2601292543072,605.403559505613,602.9845091183953),c(210.5258925068798,216.576352589601,220.97643295023948,350.243331445664,169.83449262612632,0.0,203.69728795428142),c(534.9076248389089,441.67933323391065,511.62749317634604,4151.032076393055,1589.6508509805424,2825.2166110261937,1785.6448824648949),c(653.0599114499129,670.8750953444805,657.9733630193806,1602.0389419829446,2937.457384461481,3930.3183466316777,3512.5114430325343))
targetgene="Ncor1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(380.59357778029465,445.4310495779825,398.2240203298511,3565.1311947771605,3043.4341078601838,720.7185232209678,2441.3271974222084),c(337.6291099217477,405.18536515975745,369.94603352751176,9616.557643643911,614.1215253360728,663.0610413632903,402.3274791932822),c(361.6176044761031,390.51956490565846,344.00025182227256,2689.5228661630003,690.2073780325774,67.26706216729032,375.9785762738229),c(402.43384894172266,400.41045344912055,457.695250305905,1122.0758581499977,1874.9727985924346,1364.5604039650323,1281.9754689659999),c(594.69984260872,559.687865511079,578.3868640808379,1212.8796848210957,1035.311067048866,1835.4298391360646,1421.8273383077453),c(289.652120813037,308.3228704582666,302.60361134874483,1532.8550740430605,1312.4809590147042,816.8143263170969,1643.766174437037))
targetgene="Nosip"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(499.46193885560774,460.77898007645814,506.08850607073316,784.804501943062,904.878176712001,1691.286134491871,1171.5127605728821),c(529.8951035887451,448.841700799866,527.369877581772,1524.207090550575,1036.669742989875,874.4718081747742,888.7687638602229),c(467.23858796169753,517.0547252375355,544.2783645357482,471.31510034046147,457.87379212003657,470.86943517103225,684.0580565628853),c(518.7959493919539,534.107981346953,483.9325576482817,2146.861902009533,1211.9389393800375,1393.389144893871,840.1246353935288),c(673.4680336827226,604.708461639941,655.6411579222805,2062.5440629577993,3052.9448394472465,1864.2585800649033,1109.6941806464586),c(471.53503474755223,478.51436643025227,464.69186559720544,1236.661639425431,1601.8789344496236,1095.492155295871,1622.4843682328585))
targetgene="Thrap3"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
Sweave("SP_cd19.10m_v_Input1_DOI_summary.Rnw");
library(tools);

texi2dvi("SP_cd19.10m_v_Input1_DOI_summary.tex",pdf=TRUE);

