pdf(file='BM_cd19.15m_v_Input1_DOI.pdf',width=4.5,height=4.5);
gstable=read.table('BM_cd19.15m_v_Input1_DOI.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("H2-K1","H2-T23","Rfk","Cxcr4","Gclc","Paf1","Ndufs2","Kdsr","Tap1","Rpia")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2480_BM,AR2481_BM,AR2482_BM,AR2483_BM_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 neg.'


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
targetmat=list(c(744.7075980063084,735.5043477037707,726.9954292380494,1293.7386285964458,2.2889574912971837,98.7638153559017,516.2395960613356),c(630.7441044915114,624.2193420512002,620.1908218549694,45.9810174945425,2.2889574912971837,345.6733537456559,8.38051292307363),c(562.3660083826331,536.8115789907521,612.5877820073603,282.1562437165108,18.31165993037747,25.542366040319404,43.578667199982874),c(664.9331525459505,651.50761929934,613.3118810404659,10.450231248759659,384.54485853792687,15.325419624191642,10.056615507688356),c(442.70434019209625,426.3793320021859,406.58160708880945,238.2652724717202,13.733744947783102,34.0564880537592,8.38051292307363),c(478.20835163324455,468.16450653840013,487.68069879663966,4.180092499503863,233.47366411231275,108.98076177202945,868.221138830428))
targetgene="H2-K1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(589.9802395035263,645.5383086513095,640.8276442984798,336.497446210061,13.733744947783102,83.43839573171005,236.33046443067636),c(974.8261906804179,960.2062556689227,940.6046440042093,202.73448622593736,354.7884111510635,126.00900579890906,321.8116962460274),c(801.2510236348039,721.4338297476986,802.3017286810345,91.962034989085,629.4633101067255,168.57961586610807,53.63528270767123),c(784.5948207364873,750.4276243238472,744.0117565160316,45.9810174945425,0.0,207.74457712793114,611.7774433843749),c(443.1426613209993,503.5539910945816,530.7645912664245,227.81504122296056,1046.053573522813,645.3704486187369,968.7872939073116),c(468.5652867973771,455.37312657833456,446.0450043930661,41.800924995038635,4.5779149825943675,146.44289863116458,447.5193900921318))
targetgene="H2-T23"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(401.50215407520807,411.4560553821094,391.7375769101441,18.810416247767385,991.1185937316806,44.273434469886965,98.89005249226884),c(506.2609038830407,529.9895096787171,568.0556914713643,125.4027749851159,224.31783414712402,74.92427371827024,65.36800079997431),c(527.7386391992909,551.3084762788264,482.2499560483475,163.02360748065067,45.779149825943676,444.43716910155763,53.63528270767123),c(492.2346277581426,522.3146817026777,499.99038235943533,167.20369998015454,93.84725714318454,35.75931245644716,30.16984652306507),c(743.3926346195992,677.9431378834756,738.5810137677394,252.89559621998373,444.05775331165364,98.7638153559017,742.5134449843235),c(474.70178260202,465.606230546387,458.7167374724146,50.161109994046356,1496.9781993083582,56.193205288702686,264.8242083691267))
targetgene="Rfk"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(749.0908092953391,820.7802141042079,740.3912613505034,629.1039211753314,865.2259317103354,68.1129761075184,83.8051292307363),c(543.0798787108982,493.3208871265291,598.4678508618006,33.440739996030906,217.45096167323246,156.65984504729235,38.5503594461387),c(707.0119809206448,786.2434882120308,730.6159244035775,20.900462497519317,25.17853240426902,17.0282440268796,429.08226166136984),c(477.3317093754384,529.9895096787171,534.0230369153998,58.521294993054084,6477.74970037103,144.74007422847663,68.72020596920376),c(665.8097948037566,649.3757226393292,669.4295561061521,342.76758495931676,654.6418425109946,69.81580051020637,217.89333599991437),c(555.352870320184,544.0600276347892,558.2803545244383,127.49282123486783,0.0,80.03274692633413,10.056615507688356))
targetgene="Cxcr4"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(539.1349885507707,558.1305455908614,517.0067096374175,12.540277498511589,11.444787456485919,194.12198190642746,112.29887316918663),c(435.6912021296472,462.1951958903695,453.64804424067523,71.06157249156567,270.09698397306767,63.00450289945453,522.9440063997945),c(465.0587177661526,575.612098202951,578.5551274513958,35.530786245782835,153.36015191691132,11.91977081881572,93.86174473842465),c(441.8276979342901,545.3391656307958,479.7156094324778,12.540277498511589,20.600617421674652,304.8055680811449,98.89005249226884),c(601.814909983909,524.4465783626887,552.4875622595932,769.1370199087108,18.31165993037747,119.19770818815721,129.0598990153339),c(566.3108985427607,500.5693357705663,505.42112510772756,133.76295998412363,27.467489895566203,825.8698353036607,254.76759286143835))
targetgene="Gclc"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(548.339732257735,582.4341675149859,500.35243187598815,163.02360748065067,32.045404878160575,34.0564880537592,186.04738689223458),c(287.10033943150796,333.8550169577116,314.9830794009477,169.29374622990647,18.31165993037747,6.811297610751841,491.0980572921147),c(405.4470442353357,447.27191927029304,394.99602255911947,179.74397747866612,11.444787456485919,170.28244026879602,88.83343698458047),c(297.18172539627847,350.05743157379464,322.9481687651096,20.900462497519317,0.0,119.19770818815721,1144.7780652918577),c(394.05069488385595,438.3179532982471,394.63397304256665,114.95254373635625,109.86995958226481,679.4269366724961,16.76102584614726),c(502.7543348518162,492.89450779452693,516.2826106043119,169.29374622990647,48.06810731724086,20.433892832255523,2554.3803389528425))
targetgene="Paf1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(277.01895346673746,244.31535723725253,221.93635364687464,0.0,7075.167605599595,32.353663651071244,209.51282307684073),c(358.10836231380455,362.8488115338602,353.72237767209873,33.440739996030906,27.467489895566203,56.193205288702686,38.5503594461387),c(296.3050831384723,290.7907044254908,263.5720480504482,6.2701387492557945,359.3663261336578,0.0,5.028307753844178),c(462.42879099273415,477.5448518424482,452.19984617446397,186.0141162279219,66.37976724761833,153.25419624191642,83.8051292307363),c(510.20579404316834,489.0570938065072,496.3698871939072,181.83402372841806,25.17853240426902,45.976258872574924,35.19815427690924),c(554.0379069334749,603.3267547830931,588.6925139148747,158.8435149811468,176.24972682988314,338.8620561349041,102.24225766149829))
targetgene="Ndufs2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(539.1349885507707,555.9986489308504,597.0196527955893,2675.2591996824726,22.889574912971838,68.1129761075184,154.2014377845548),c(897.6816719934783,878.341423924503,877.6080281240197,376.2083249553477,2893.24226899964,1776.0458520035424,492.7741598767294),c(624.6076086868684,713.332622439657,617.6564752390997,129.58286748461975,157.93806689950569,83.43839573171005,95.53784732303937),c(1262.803172369732,1285.1073066545885,1229.8822077299071,568.4925799325254,5401.939679461354,2162.5869914137093,2033.1124351376625),c(653.0984820655677,642.9800326592964,666.8952094902824,56.431248743302156,469.2362857159227,282.6688508462014,83.8051292307363),c(657.0433722256953,704.3786564676111,717.944191324229,194.37430122692965,29.75644738686339,105.57511296665353,940.2935499688613))
targetgene="Kdsr"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(740.7627078461808,729.5350370557401,723.3749340725212,217.36480997420088,309.0092613251198,56.193205288702686,125.70769384610445),c(629.8674622337053,557.7041662588592,614.0359800735715,744.0564649116876,0.0,22.136717234943482,1540.3382752609332),c(706.5736597917416,752.5595209838582,701.6519630793524,353.21781620807644,1252.0597477395595,245.20671398706628,274.88082387681504),c(639.5105270695726,677.9431378834756,742.5635584498203,440.99975869765757,1023.1639986098411,357.5931245644716,177.66687396916095),c(575.077321120822,581.1550295189794,524.2476999684739,1065.9235873734851,9.155829965188735,39.16496126182309,73.74851372304794),c(603.5681944995213,567.9372702269117,580.3653750341599,20.900462497519317,38.91227735205212,151.55137183922847,249.73928510759416))
targetgene="Tap1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(427.801421809392,408.4714000580941,410.2021022543376,6.2701387492557945,0.0,5.108473208063881,20.113231015376712),c(765.3086910647525,788.8017642040439,749.080449747771,252.89559621998373,141.9153644604254,490.41342797413256,184.37128430761985),c(604.8831578862305,613.5598587511455,545.9706709616427,246.62545747072792,352.4994536597663,1234.5476919487712,734.13293206125),c(355.0401144114831,333.42863762570937,380.5140418970069,148.39328373238715,0.0,37.46213685913512,261.47200319989724),c(451.47076277015753,544.0600276347892,524.9717990015795,22.99050874727125,2.2889574912971837,136.2259522150368,30.16984652306507),c(600.4999465971998,631.0414113632352,640.1035452653741,330.2273074608052,112.15891707356201,187.31068429567563,169.2863610460873))
targetgene="Rpia"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetgenelist=c("Arl15","Ints10","Cnot8","Ints13","Cd19","Nosip","Ints14","Sipa1","Maea","Setd5")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='AR2480_BM,AR2481_BM,AR2482_BM,AR2483_BM_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 pos.'


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
targetmat=list(c(684.6576033465884,701.3940011435958,618.0185247556525,11507.794651134136,7903.770217449175,7151.862491289433,1434.7438124302055),c(940.1988214970756,857.0224573243937,881.2285232895479,8565.009531483416,4758.742624406845,7832.992252364617,3695.806199075471),c(731.1196430103133,738.0626236957838,739.305112800845,19903.510436387645,3220.5631902551377,3678.100709805994,1196.7372454149142),c(771.0068657404923,805.4305581521292,684.2735862848174,29983.80349894121,1698.4064585425103,2797.7404936163184,1682.8069949531848),c(966.4980892312595,978.5405669450167,934.8118517393642,1270.7481198491744,8615.6359972426,10314.007407080975,10411.949255626678),c(765.3086910647525,822.0593521002145,772.9757178402566,2779.761512170069,36488.27136876841,3400.5403321678564,6142.915972612971))
targetgene="Arl15"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(678.5211075419455,611.0015827591325,696.583269847613,9077.07086267264,489.83690313759735,1098.3217397337344,2633.1571604297346),c(650.9068764210524,684.7652071955106,716.4959932580177,5628.494550581952,4793.076986776303,2751.7642347437436,2490.688440737483),c(494.42623340265794,487.3515764784985,572.400285669998,4119.481158261057,915.5829965188735,13973.377048457402,1525.2533519994006),c(705.6970175339355,838.6881460482997,745.09790506569,10124.184033798358,1146.7677031398891,5896.880906508406,1134.7214497841694),c(401.9404752041111,420.8364006861575,379.42789334734846,5914.830886797967,1854.0555679507188,124.30618139622109,50.28307753844178),c(692.1090625379405,689.8817591795369,670.1536551392577,6311.939674250833,4019.4093547178545,13522.128581745092,2651.5942888604964))
targetgene="Ints10"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(732.8729275259257,813.5317654601707,789.9920451182389,808.8478986539975,524.1712655070551,1197.085555089636,174.3146687999315),c(547.4630899999289,621.6610660591871,606.4329402259625,815.1180374032533,6116.094416746075,1576.815396889051,787.7682147689212),c(564.1192928982454,529.1367510147127,499.6283328428825,2244.7096722335746,5445.429871796,1311.1747900697294,3171.1860900910615),c(768.815260095977,726.1240023997226,789.9920451182389,2765.1311884218057,12568.665584712835,6096.111361622898,556.466058092089),c(570.6941098317914,549.1765796188155,542.7122253126673,2482.974944705295,3392.2350021024263,1253.2787603783388,3647.1992241216435),c(682.9043188309762,619.9555487311783,662.1885657750958,3927.1969032838797,498.9927331027861,3953.9582630414434,3905.3190221523114))
targetgene="Cnot8"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(482.5915629222752,486.49881781449415,442.7865587440908,627.0138749255794,1036.8977435576242,1854.3757745271887,4377.9799510136645),c(517.2189321056173,540.6489929787717,565.5213448554946,315.59698371254166,2213.4218940843766,4527.810086747286,983.8722171688441),c(636.0039580383482,586.6979608350078,605.3467916763041,10504.572451253209,1952.4807400764978,6661.4490633153,248.06318252297945),c(673.2612539951087,546.6183036268023,635.0348520336347,760.7768349097031,469.2362857159227,2755.1698835491197,4103.099127136849),c(747.7758459086299,776.0103842439784,813.5252636941717,5482.191313099317,5520.965469008807,6823.2173815706565,1099.5232955072602),c(510.20579404316834,598.6365821310691,513.0241649553366,5371.418861862464,460.0804557507339,1171.5431890493167,511.21128830749143))
targetgene="Ints13"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(523.3554279102603,472.00192052641984,561.9008496899664,2217.5390709867993,7679.452383302051,902.4969334246189,705.6391881227996),c(618.4711128822255,648.0965846433226,648.4306841460889,5457.110758102293,814.8688669017974,597.691365343474,1009.0137559380651),c(602.253231112812,691.1608971755434,613.3118810404659,2110.946712249451,3682.9326034971687,2945.8862166501713,2135.354692799161),c(679.8360709286547,678.79589654748,675.222348370997,18390.316951567245,16356.890232809676,1319.6889120831693,1252.0486307072003),c(653.0984820655677,588.4034781630165,615.4841781397828,2163.197868493249,29081.20492693072,1132.3782277874936,1124.664834276481),c(455.41565293028515,423.39467667817064,481.16380749868904,6010.973014286556,25004.571634930435,3424.379873805488,3000.2236264603594))
targetgene="Cd19"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(465.9353600239587,556.8514075948548,494.55963961114315,5442.48043435403,4490.934597925075,1263.4957067944665,953.7023706457791),c(413.33682455559085,506.5386464185969,459.4408365055202,1281.1983510979342,450.9246257855452,194.12198190642746,781.0638044304623),c(442.70434019209625,488.2043351425029,427.2184295323198,1394.0608485845385,59.51289477372678,1100.0245641364222,1962.716126583844),c(492.6729488870457,500.5693357705663,568.4177409879171,2936.514980901464,812.5799094105002,1857.7814233325646,1153.1585782149314),c(728.0513951079919,699.6884838155871,718.3062408407818,2039.8851397578853,1915.8574202157429,1323.0945608885452,1726.3856621531677),c(354.60179328258005,385.44691612997605,375.8073981818203,1626.0559823070028,656.9308000022917,1861.1870721379405,77.1007188922774))
targetgene="Nosip"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(787.6630686388088,770.4674529279499,810.9909170783021,1711.747878546832,3385.368129628535,2859.042172113085,1932.546280060779),c(735.502854299344,711.6271051116483,692.600725165532,5434.120249355023,1091.8327233487566,1215.8166235192036,6287.060794889837),c(860.4243760367177,860.4334919804112,864.2121960115657,620.7437361763238,4742.719921967765,5571.641445595006,1240.3159126148971),c(489.16637985582116,429.7903666582034,466.6818268365765,432.6395736986498,320.4540487816057,2191.535006259405,1005.6615507688356),c(449.71747825454526,472.00192052641984,442.7865587440908,131.6729137343717,613.4406076676453,3807.515364410279,169.2863610460873),c(754.7889839710789,864.2709059684308,809.5427190120909,599.8432736788044,6624.24297981405,6033.106858723443,564.8465710151627))
targetgene="Ints14"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(590.4185606324294,558.1305455908614,537.2814825643751,332.31735371055714,1014.0081686446524,1019.9918172100881,160.9058481230137),c(807.8258405683498,778.1422809039893,751.9768458801934,3185.230484621944,7178.170692707969,1227.7363943380194,2869.487624860411),c(444.89594583661153,455.37312657833456,473.19871813452716,1379.4305248362748,3483.7933017543137,1772.6402031981665,936.9413447996318),c(668.001400448272,635.7315840152592,601.0021974776703,2677.349245932224,2776.505436943484,2319.2468364610017,3922.0800479984587),c(719.2849725299307,654.9186539553575,752.700944913299,1845.5108385309557,979.6738062751947,1515.5137183922845,2274.471207322183),c(404.5704019775295,401.2229514140569,402.96111192328135,3950.187412031151,4754.164709424251,1195.382730686948,1225.2309893533647))
targetgene="Sipa1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(653.5368031944707,710.3479671156417,674.8602988544442,2046.155278507141,3541.0172390367434,4478.428179069336,859.8406259073544),c(621.97768191345,667.7100339154232,658.9301201261204,823.478222402261,601.9958202111593,531.2812136386436,839.7273948919777),c(695.177310440262,746.5902103358276,654.9475754440394,2207.08883973804,1384.8192822347962,2864.150645321149,2167.2006419068407),c(884.0937169974833,917.9947018007063,978.2577937257018,1019.9425698789427,4497.801470398966,1312.8776144724172,6913.923161535745),c(710.9568710807723,761.0871076239018,726.2713302049438,2689.889523430736,3348.7448097677798,434.22022268542986,1731.413969907012),c(592.1718451480416,613.1334794191433,539.0917301471392,202.73448622593736,494.4148181201917,883.7658649950514,789.4443173535359))
targetgene="Maea"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
targetmat=list(c(668.8780427060781,609.7224447631258,624.5354160536032,758.6867886599512,1391.6861547086878,294.5886216650171,1040.8597050457447),c(651.7835186788585,650.2284813033335,604.6226926431984,3822.694590796283,8844.531746372319,1794.77692043311,2100.1565385222516),c(483.4682051800813,558.1305455908614,554.6598593589101,1306.2789060949574,1192.5468529658326,233.28694316825056,1273.8379643071917),c(704.8203752761294,681.3541725394931,731.3400234366832,2294.870782227621,6427.392635562492,14325.86169981381,4253.948359752175),c(396.2423005283713,382.88864013796297,431.9250732475064,1793.2596822871574,6617.376107340158,9857.650467160602,6615.576901474324),c(660.1116201280167,612.7071000871412,637.5691986495044,1111.9046048680277,2298.1133212623727,1254.9815847810266,2998.547523875745))
targetgene="Setd5"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","AR2480_BM","AR2481_BM","AR2482_BM","AR2483_BM")

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
Sweave("BM_cd19.15m_v_Input1_DOI_summary.Rnw");
library(tools);

texi2dvi("BM_cd19.15m_v_Input1_DOI_summary.tex",pdf=TRUE);

