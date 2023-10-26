pdf(file='SP_cd19.10m_v_egfr.15m.pdf',width=4.5,height=4.5);
gstable=read.table('SP_cd19.10m_v_egfr.15m.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Sipa1","Pak2","Ubr4","Fitm2","Irf1","Ptpn2","H2-T23","Tlnrd1","Capzb","Ifngr1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2473_SP,AR2474_SP,AR2475_SP,AR2476_SP_vs_AR2466_SP,AR2467_SP,AR2468_SP,AR2469_SP neg.'


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
targetmat=list(c(301.2367681873479,416.58413051731105,378.6922227412089,573.7371530607197,182.91212643826498,126.22940792478109,527.9996080444457,291.8842315888952),c(2848.5490349041793,1706.1836997709,1990.7398498230523,2738.649185282968,759.9870042153264,1317.4671979434106,946.9558187753646,654.3071524784401),c(2059.91794515528,1055.5879620238647,1804.5205549704947,2500.6428388209665,892.6627015613917,927.9115416987219,470.6083463004842,723.0215653316592),c(1892.0376788396343,2506.0252303641373,2264.510007180917,2172.00493658701,1516.1096677312526,706.3831105724505,837.9124214618378,1437.529840575309),c(5435.800558526749,2284.3300235149245,3076.7874537579323,3250.914500515754,1279.0967715013883,1193.7456590502477,218.08679462705365,934.0295410844647),c(1317.9954778893625,2107.553453347579,2182.5179295965822,911.8322611143581,2618.734880626709,952.1542756710309,258.2606778478267,316.2079175546365))
targetgene="Sipa1"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(3820.629931796385,2849.4354527384075,2605.6804317055658,3071.227589941792,1017.6097175086572,851.8395144063041,390.2605798589381,651.2666917327225),c(499.5791796005904,609.2995717653192,397.453121849489,616.2945792492895,425.0774769339961,738.9854079834867,659.9995100555572,1005.7844146834015),c(1073.6213805508623,1291.7730516736792,684.425393394662,1137.2289998167837,191.92892140353155,234.06777628436228,149.21728053429987,142.29356289958642),c(2128.9654740431665,1205.5582490100967,827.5641051096877,949.6610843930869,596.3965812740612,511.60528245010613,361.56494898695735,331.41022128322476),c(1380.2736412000052,1356.9775242763887,803.9392691955572,816.4721024325626,319.45216448373037,274.193680790253,57.391261743961486,108.24040254754864),c(1702.4954426768088,1714.8776294512613,1108.2827436187674,1636.8847072899928,807.6472061745925,469.80746525647004,286.95630871980745,335.05877417808597))
targetgene="Pak2"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(391.26933036469,412.9616598171605,566.3012138240097,165.50110184443838,91.45606321913249,212.3329113436715,177.9129114062806,138.0369178555817),c(337.11440574673986,475.26815585974964,505.1545796933191,399.5669458815726,256.33459972686427,361.1331405530161,269.738930196619,205.53514641051373),c(678.9673673975503,705.6572923893234,414.8243247275261,1026.8949319204914,425.0774769339961,387.0477872130705,143.47815435990373,282.76284935174226),c(546.2878020835724,302.1140563925543,432.8903757206847,336.5189070836914,358.09557147773,169.69913780616267,149.21728053429987,217.08889724424083),c(594.3502976820032,479.6151206999303,474.58126262797373,624.9636845839982,162.3023093747985,294.2566330431983,195.13028992946906,403.1650948821615),c(446.778128098089,633.9323725263429,549.6248590610941,582.4062583954284,288.5374388885307,205.6452605926897,355.8258228125612,209.18369930537492))
targetgene="Ubr4"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(305.29838753369415,252.12396073047694,199.42140903986598,364.8905245427379,51.52454265866619,46.813555256872455,103.30427113913068,66.28204425664495),c(444.07038186719143,241.98104277005547,162.59445893842732,376.7120318173407,142.98060587779867,275.02963713412566,212.3476684526575,108.84849469669217),c(541.5492461795018,281.1037263316812,380.7767670865734,328.6379022339562,87.59172251973253,74.4001146046723,183.65203758067676,116.14560048641455),c(318.1601821304573,213.0012771688512,410.6552360367972,449.21727643490414,128.81135664666547,138.76875308287194,74.60864026714994,122.8346141269934),c(450.8397474444352,139.8273690258105,480.14004754894563,695.1046277466412,182.91212643826498,99.47880492085397,252.52155167343054,117.36178478470163),c(196.3116017400694,97.80670890406434,96.58388800188631,241.94684888686942,50.236429092199536,64.36863847819963,183.65203758067676,95.47046741553447))
targetgene="Fitm2"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(224.06600060676885,681.7489857683299,459.9894522104226,366.46672551268495,453.41597539626247,119.54175717379931,252.52155167343054,179.9952761464854),c(613.3045212982857,557.1359936831516,245.28138463788395,235.6420450070813,180.33589930533168,189.762090059108,154.95640670869602,437.8263473833428),c(256.55895537753895,452.08434337878623,472.4967182826093,409.81225218622836,274.36818965739747,50.15738063236335,309.91281341739204,144.72593149616054),c(364.19186805571496,594.8096889647171,557.9630364425518,425.57426188569866,78.57492755446594,40.125904505890674,63.130387918357634,88.78145377495564),c(340.49908853536175,363.6960582951133,607.2972526161773,498.86760698823565,398.0270920381963,396.2433069956704,361.56494898695735,473.09569203366766),c(226.096810279942,328.1958454336381,248.75562521349136,286.86857653035986,90.16794965266584,336.05445023683444,74.60864026714994,133.17218066243345))
targetgene="Irf1"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(109.66372235134911,127.51096864529868,227.21533364472535,243.52304985681644,185.48835357119827,68.54842019756325,45.91300939516919,39.52598969432956),c(374.3459164215806,252.84845487050706,403.70675488558237,329.4260027189297,105.6253124502657,193.9418717784716,200.8694161038652,209.79179145451843),c(452.193620559884,609.2995717653192,377.99737462608743,469.70788904421556,78.57492755446594,152.9800109287082,45.91300939516919,272.42528281630223),c(136.06424810259983,146.34781628608144,87.55086250530701,208.84662851798174,9.016794965266584,44.30568622525429,120.52164966231912,13.378027281157697),c(536.1337537177068,508.5948863011345,570.4703025147386,421.6337594608311,374.8410478417965,650.374035532978,476.34747247488036,370.3281188284108),c(268.0668768588534,231.83812480963397,106.3117616135871,278.19947119565114,36.067179861066336,50.993336976236066,154.95640670869602,158.71205092646179))
targetgene="Ptpn2"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(193.6038555091719,210.82779474876088,183.43990239207184,215.15143239776987,41.21963412693295,260.8183792882894,11.478252348792298,173.91435465505006),c(253.85120914664145,277.48125563153064,385.64070389242374,245.887351311737,108.201539583199,97.80689223310853,241.04329932463824,83.91671658180738),c(254.52814570436584,160.83769908668356,258.48349882519216,194.66081978845847,33.490952728133024,145.4564038338537,114.78252348792297,52.29592482634373),c(178.7112512392356,271.6853025112898,150.7820409813621,173.38210669417353,50.236429092199536,27.58655934779984,40.17388322077304,38.3098053960425),c(216.61969847180072,302.1140563925543,190.3883835432867,284.5042750754393,114.64210741553228,224.03630015788963,218.08679462705365,174.5224468041936),c(185.48061681647937,113.74557998472666,45.85997559801796,33.100220368887676,21.897930629933132,35.94612278652706,68.86951409275379,21.28322522002361))
targetgene="H2-T23"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(647.1513491845046,641.901808066674,368.26950101438666,460.25068322453336,266.6395082585975,292.5847203554529,723.1298979739147,305.87035101919645),c(230.83536618401263,173.15409946719538,168.84809197452068,169.44160426930594,88.87983608619918,80.25180901178135,206.60854227826135,142.90165504872996),c(163.8186469692993,192.71544124800823,121.59842014625974,212.78713094284933,130.09947021313212,183.07443930812622,68.86951409275379,222.56172658653261),c(511.0871010819048,338.3387633940596,293.2259045812664,330.2141032039032,92.74417678559914,147.12831652159915,195.13028992946906,199.4542249190784),c(387.8846475760681,157.21522838653303,310.59710745930346,301.8424857448566,114.64210741553228,171.3710504939081,68.86951409275379,175.13053895333712),c(511.0871010819048,495.5539917805926,308.512563113939,457.0982812846393,115.93022098199893,285.8970696044711,86.08689261594223,252.35824189456565))
targetgene="Tlnrd1"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(176.68044156606246,321.6753981733671,496.8164023118612,423.99806091575164,279.52064392326406,184.74635199587166,80.34776644154609,164.18488026875357),c(884.756080945761,603.5036186450784,284.8877271998085,360.16192163289685,303.99480168613053,617.7717381219419,258.2606778478267,322.2888390460718),c(879.340588483966,643.3507963467342,403.70675488558237,654.9115030129918,553.8888335806615,308.4678908890346,872.3471785082146,429.9211494444769),c(1163.6539427282044,979.5160773207035,908.8613345789014,1072.6047600489553,69.55813258919936,125.39345158090836,91.82601879033838,196.41376417336073),c(166.5263932001968,111.57209756463635,191.08323165840818,182.83931251385573,133.9638109125321,109.51028104732664,57.391261743961486,110.06467899497925),c(934.1724496596405,700.5858334091126,682.3408490492975,656.4877039829389,28.338498462266404,188.92613371523527,235.3041731502421,99.11902031039567))
targetgene="Capzb"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(481.9788290997566,488.3090503802915,578.1136317810749,528.8154254172292,1339.6381091253209,591.0211351180147,2066.0854227826135,546.0667499308914),c(282.95948112878966,265.1648552510189,332.13739902806947,136.34138390041826,261.4870539927309,394.571394307925,223.8259208014498,314.3836411072059),c(270.0976865320265,197.06240608818888,119.51387580089529,158.40819747967672,81.15115468739926,34.27421009878162,45.91300939516919,53.51210912463079),c(206.46565010593505,226.76666582942323,148.69749663599762,274.25896877078355,170.03099077359843,70.22033288530868,80.34776644154609,51.07974052805666),c(222.03519093359571,207.92981818864047,155.64597778721247,247.46355228168403,119.79456168139889,114.52601911056297,86.08689261594223,109.4565868458357),c(468.44009794526903,241.25654863002535,257.78865071007067,446.85297497998357,112.06588028259897,131.2451459880174,154.95640670869602,89.9976380732427))
targetgene="Ifngr1"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetgenelist=c("Cd19","Sephs1","Rpn2","Etv5","Uggt1","B3gnt2","Tsc1","Ptpn7","Ints13","Ccdc134")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2473_SP,AR2474_SP,AR2475_SP,AR2476_SP_vs_AR2466_SP,AR2467_SP,AR2468_SP,AR2469_SP pos.'


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
targetmat=list(c(304.62145097596976,205.75633576855014,347.42405756074214,192.2965183335379,29830.13397223479,2858.9706960447106,477156.6892654702,1089.701131265209),c(189.54223616282562,212.27678302882111,186.91414296767925,218.30383433766394,12704.664106060616,3197.5330153131636,241955.82038636724,862.274667485528),c(398.71563249965817,758.5453646115211,343.2549688700132,406.65985024633426,28516.2581344388,3823.6643168738324,468370.0870924697,2360.0056308260464),c(488.07125811927597,202.85835920842973,471.80187016748783,347.5523138733206,68096.12369125971,11495.235684593807,1995815.5619033072,5132.297738771407),c(377.7305992102025,399.19627115658847,323.1043735314902,345.9761129033735,22030.606327279198,4825.9759731772265,264298.2385832915,1899.0717817752495),c(488.74819467700036,489.75803866035176,350.20345002122804,779.4313796388074,28285.68580604127,5955.352993749275,1304296.7708979663,1851.0325019929105))
targetgene="Cd19"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(292.43659293693094,339.0632575340897,257.78865071007067,161.5605994195708,377.4172749747298,434.6972988138157,1291.3033892391336,363.63910518783194),c(251.8203994734683,445.5638961185153,523.2206306864776,330.2141032039032,928.7298814224581,1142.7523220740115,941.2166926009684,799.6411761237442),c(156.37234483433113,194.88892366809856,300.8692338476027,135.55328341544475,915.8487457577916,759.8843165803047,315.65193959178816,720.5891967350851),c(388.5615841337925,544.0950991626097,512.1030608445338,471.2840900141626,2622.5992213261093,1850.8073453342074,1256.8686321927566,1287.9391718860002),c(124.55632662128541,257.91991385071776,184.13475050719333,96.93635965174248,565.4818556788614,567.6143574895785,579.651743614011,980.2445444193731),c(181.4189974701331,365.8695407152036,185.5244467374363,231.7015425822137,1602.4132766845185,1092.5949414416482,1090.4339731352682,581.3360945812163))
targetgene="Sephs1"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(382.4691551142731,249.22598417035653,507.2391240386835,358.58572066294977,2389.4506657956445,440.54899322092473,298.43456106859975,499.85174659598306),c(65.66284609926458,55.786048782318176,86.16116627506405,65.41234025280183,681.4120766608604,120.37771351767203,332.86931811497664,140.4692864521558),c(467.08622482982025,506.42140388104417,400.23251430997493,311.2996915645388,966.0851748499911,494.88615557265166,1124.868730181645,837.9509815197868),c(308.683070322316,357.90010517487246,263.3474356310425,290.8090789552274,1188.9288218487222,763.2281419557955,1767.650861714014,419.58358290903686),c(159.08009106522866,99.98019132415465,150.7820409813621,150.52719262994157,571.9224235111947,234.90373262823502,1061.7383422632874,226.2102794813938),c(162.46477385385055,163.735675646804,282.10833473932263,232.48964306718722,703.3100072907935,1375.1481856706284,1021.5644590425145,875.0446026175422))
targetgene="Rpn2"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(261.2975112816096,381.0839176558358,288.36196777541596,208.84662851798174,480.46636029206223,839.3001692482134,430.4344630797112,598.3626747572353),c(452.193620559884,636.8303490864633,434.97492006604915,367.2548259976585,999.5761275781241,596.8728295251238,763.3037811946878,638.4967566007083),c(386.5307744606194,336.88977511399935,462.073996555787,289.23287798528037,1990.1354601909816,1142.7523220740115,625.5647530091802,549.7153028257527),c(144.18748679529236,222.41970098924259,369.65919724462964,293.17338041014796,296.2661202873306,583.4975280231603,1543.824940912564,465.7985862439453),c(186.15755337420373,292.6956325721629,291.1413602359019,182.0512120288822,656.9379188979939,763.2281419557955,384.521453684542,515.0540503245713),c(261.2975112816096,220.24621856915226,337.0013358339199,198.60132221332603,1561.1936425575855,581.8256153354148,912.5210617289877,551.5395792731832))
targetgene="Etv5"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(104.92516644727847,44.91863668186658,55.58784920971874,308.14728962464477,303.99480168613053,189.762090059108,321.39106576618434,58.98493846692258),c(253.85120914664145,244.15452519014576,233.4689666808187,111.12216838126577,390.29841063939637,390.3916125885614,1107.6513516584566,412.89456926845804),c(117.78696104404163,110.12310928457613,116.73448334040935,199.38942269829957,1102.6252128954566,783.2910942087409,1463.4771744710179,341.1396956695213),c(282.95948112878966,69.55143744289019,99.36328046237224,170.22970475427945,225.41987413166459,212.3329113436715,143.47815435990373,431.74542589190753),c(112.37146858224662,69.55143744289019,134.80053433356795,78.8100484973516,189.35269427059825,94.46306685761763,200.8694161038652,147.15830009273466),c(383.8230282297219,264.44036111098876,340.47557640952726,316.81639495935343,533.279016517195,576.8098772721785,1130.6078563560413,330.1940369849377))
targetgene="Uggt1"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(215.94276191407633,197.06240608818888,315.4610442651539,290.8090789552274,329.7570730154636,1047.4532988725211,527.9996080444457,217.08889724424083),c(345.23764443943236,489.75803866035176,289.75166400565894,226.9729396723726,1642.3447972449849,377.85226743047053,1182.2599919256068,333.2344977306554),c(439.3318259631208,386.1553766360466,275.85470170322924,369.619127452579,373.5529342753299,654.5538172523417,878.0863046826107,598.9707669063787),c(284.31335424423844,252.12396073047694,242.501992177398,195.44892027343198,276.94441679033076,181.40252662038077,172.17378523188447,130.13171991671578),c(298.5290219564504,449.91086095869593,319.63013295588274,218.30383433766394,753.5464363829931,844.3159073114497,510.78222952125725,413.5026614176016),c(66.33978265698897,136.20489832565997,34.74240575607421,56.74323491809315,28.338498462266404,238.24755800372589,614.086500660388,259.655347684288))
targetgene="B3gnt2"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(310.71387999548915,510.7683687212248,397.453121849489,406.65985024633426,811.5115468739925,616.9357817780692,1974.259403992275,670.1175483561719),c(310.71387999548915,197.06240608818888,587.1466572776542,290.8090789552274,1657.8021600425848,511.60528245010613,1262.6077583671527,539.9858284394562),c(161.78783729612616,368.04302313529394,228.60502987496832,543.7893346317261,479.17824672559556,517.4569768572152,1796.3464925859946,621.4701764246894),c(145.5413599107411,157.93972252656314,132.7159899882035,61.47183782793425,95.32040391853245,113.69006276669025,5.739126174396149,128.3074434692852),c(270.7746230897509,359.34909345493264,334.9167914885554,277.41137071067766,445.6872939974625,669.6010314420507,315.65193959178816,385.5304225569991),c(200.37322108641567,115.91906240481698,201.50595338523044,172.59400620920002,443.1110668645292,529.1603656714333,384.521453684542,236.54784601683383))
targetgene="Tsc1"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(283.63641768651405,253.57294901053714,195.94716846425857,214.36333191279635,465.0089974944624,606.0683493077238,177.9129114062806,435.3939787867687),c(108.30984923590036,136.20489832565997,260.5680431705566,172.59400620920002,280.80875748973074,203.13739156107155,315.65193959178816,349.65298575753076),c(51.447178387052666,102.87816788427507,49.33421617362538,160.77249893459728,92.74417678559914,138.76875308287194,189.3911637550729,113.1051397406969),c(218.65050814497386,296.3181032723134,353.6776905968355,258.49695907131326,1020.1859446415906,356.9533588336525,63.130387918357634,802.0735447203183),c(199.69628452869128,156.49073424650294,107.00660972870858,210.4228294879288,693.0050987590603,433.86134246994294,895.3036832057992,384.92233040785555),c(133.3565018717023,26.806283181113926,43.080583137532024,35.46452182380822,69.55813258919936,51.82929332010879,28.695630871980743,46.82309548405194))
targetgene="Ptpn7"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(134.71037498715106,59.408519482468705,129.93659752771757,171.8059057242265,450.83974826332917,371.16461667948874,91.82601879033838,273.64146711458926),c(245.72797045394893,204.30734748848994,204.28534584571636,168.65350378433243,136.5400380454654,227.3801255333805,45.91300939516919,135.60454925900757),c(161.78783729612616,137.65388660572017,187.60899108280074,171.8059057242265,127.52324308019882,479.83894138294266,384.521453684542,274.2495592637328),c(448.1320012135377,345.58370479436064,300.1743857324812,342.82371096347947,253.758372593931,391.2275689324341,430.4344630797112,688.3603128304779),c(304.62145097596976,299.94057397246394,414.8243247275261,280.5637726505717,318.16405091726375,572.6300955528149,1021.5644590425145,577.6875416863551),c(267.389940301129,224.5931834093329,107.70145784383006,365.67862502771146,513.9573130201952,241.59138337921678,1360.1729033318873,372.1523952758414))
targetgene="Ints13"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
targetmat=list(c(373.66897986385624,454.98231993890664,522.5257825713561,378.2882327872877,507.516745187862,483.18276675843356,780.5211597178762,894.5035513901352),c(85.97094283099591,149.97028698623197,109.7860021891945,144.22238875015344,369.68859357592993,455.5962074106337,510.78222952125725,100.33520460868273),c(192.9269189514475,231.11363066960385,120.20872391601678,159.98439844962377,164.8785365077318,281.71728788510745,91.82601879033838,196.41376417336073),c(242.34328766532704,174.60308774725559,146.6129522906332,117.42697226105389,251.18214546099767,248.27903413019857,350.08669663816505,300.39752167690466),c(175.3265684506137,296.3181032723134,252.92471390422025,162.3486999045443,258.9108268597976,536.6839727662878,407.4779583821266,241.4125832099821),c(600.4427267015226,574.5238530438742,394.67372938900303,533.5440283270704,1290.689793599588,583.4975280231603,843.6515476362339,887.2064456004127))
targetgene="Ccdc134"
collabel=c("AR2466_SP","AR2467_SP","AR2468_SP","AR2469_SP","AR2473_SP","AR2474_SP","AR2475_SP","AR2476_SP")

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
Sweave("SP_cd19.10m_v_egfr.15m_summary.Rnw");
library(tools);

texi2dvi("SP_cd19.10m_v_egfr.15m_summary.tex",pdf=TRUE);

