pdf(file='cd19_1to2_v_Input1_DOI.pdf',width=4.5,height=4.5);
gstable=read.table('cd19_1to2_v_Input1_DOI.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Zfp592","Elob","Stub1","Kdsr","Ube2l3","Rps10","Sys1","Ptpn2","Baz1b","Eif2b2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 neg.'


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
targetmat=list(c(1226.6558924790045,1304.6988523727823,1257.5626429753695,222.44713827403507,219.27473254486915,222.64372808077445),c(1178.0316949392964,1331.6331804564484,1287.87642472715,277.5836084445224,307.37618758521836,274.8566377634339),c(1518.4010777172543,1671.0057143106403,1574.479452198527,496.2282315343859,544.2712111381574,557.5944694412316),c(1145.983928379034,1033.2008252894288,1130.7959192861067,437.28924617972706,409.1823134096219,397.01514343607124),c(1245.44251425571,1200.193659408158,1268.5858363396533,512.3889210671149,478.68457238589735,520.158798348004),c(826.6113581750409,890.9875730076722,894.71586140103,234.80531262259257,203.61225164880707,192.10410166261514))
targetgene="Zfp592"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(892.9170820928249,828.499931853567,852.4602868379422,154.95249375499023,196.75991625677992,198.0149970983879),c(1424.467968833727,1463.0727015047385,1385.2479661116563,423.98044303512665,351.42691510539294,442.3320084436625),c(1926.1812798116262,1918.8015326803677,2166.0574960817544,589.3898535465886,521.7563948500681,581.2380511843227),c(1748.2609206322388,1608.5180731565351,1748.0947476859958,916.406159385341,985.7573913959072,934.9066280913936),c(1315.0635243693832,1415.6682840774863,1290.6322230682208,208.1877063333918,178.16072019270618,297.5150702672296),c(1162.56035935848,1185.110435681305,1021.4825850902929,139.74243301830407,151.73028368060142,205.8961910127516))
targetgene="Elob"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1532.7673178994407,1609.5954462798818,1625.9210212318512,1359.3991783413253,1361.6569329013973,1342.758413159715),c(1680.8501013158252,1732.4159823413988,1894.1520597627555,450.5980493243274,484.5580027219207,458.09439627238993),c(1974.8054773513343,2024.3840987683386,2000.7095956174983,1387.9180422226118,1389.066274469506,1264.9316232553733),c(1148.1941191762935,1137.7060182540529,1176.725891637289,250.9660021553216,228.08487804890407,300.470517985116),c(2201.35003407043,2220.466007217427,2333.2425954400574,676.847702782534,703.8327352667898,768.4164066504605),c(1748.2609206322388,1799.2131159888904,1789.43172280206,629.3162629803898,600.0687993303785,531.9805892195495))
targetgene="Stub1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1359.2673403145727,1404.89455284402,1514.7704881419902,434.4373597915984,408.203408353618,403.9111881111395),c(2263.235376393695,2219.3886340940803,2226.685059585315,783.3181279393372,793.8920004191468,745.7579741466649),c(1574.7609430473708,1802.4452353589304,1567.130656622338,852.7140300504677,805.6388610911933,852.1540919905748),c(3183.7798434522633,3247.202593766776,3120.482321539321,1394.572443794912,1252.0195666289626,1364.431696424215),c(1646.5921439583035,1624.6786700067346,1692.0601814175536,519.0433226394151,454.21194598580036,449.22805311873077),c(1656.538002545971,1779.820399768651,1821.5827034478875,1109.3838049820467,1145.3189155245398,1134.8919236683723))
targetgene="Kdsr"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1669.7991473295278,1502.9355070685642,1506.5030931187773,589.3898535465886,701.8749251547821,736.8916309930057),c(2031.1653426814507,2051.3184268520044,1820.6641040008637,1024.7778421342298,1166.854826756625,1085.6344617035993),c(1522.8214593117732,1524.482969535497,1369.6317755122543,585.5873383624171,600.0687993303785,631.4806623883912),c(1738.3150620445713,1703.3269080110394,1704.001974228861,515.2408074552436,467.91661676985467,486.6637242119583),c(1296.2769025926777,1275.609778042423,1171.214294955147,523.7964666196295,523.7142049620759,499.47066432279934),c(1276.3851854173427,1407.0492990907132,1234.5976567997784,605.5505430793177,652.9296723545881,523.1142460658904))
targetgene="Ube2l3"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(850.9234569448951,890.9875730076722,792.7513227814054,203.43456235317737,266.2621752330554,217.71798188429713),c(1275.2800900187128,1362.8770010335008,1331.9691981842848,702.514680275692,558.9547869782156,652.1687964135959),c(1295.1718071940481,1310.0857179895156,1269.504435786677,331.7694498189668,260.3887448970321,279.7823839599113),c(1211.1845568981882,1222.8184949984375,1154.6795049087214,263.32417650387913,241.7895488329584,263.0348468918884),c(876.340651113379,860.8211255539662,767.949137711767,163.5081529193762,146.83575840058202,178.31201231247866),c(912.8087992681601,773.5539025628883,861.6462813081787,173.01444087980505,181.0974353607178,144.81693817643293))
targetgene="Rps10"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1386.894725280316,1284.2287630291962,1471.596314131879,157.8043801431189,221.23254265687692,167.47537068022856),c(1802.4105951650959,1849.8496527861826,1867.5126757990697,1270.0400715132942,1218.7367947248308,1230.4513998800323),c(1376.9488666926484,1487.8522833417112,1462.4103196616425,308.0037299178947,331.84881398531536,375.34186017157106),c(1287.4361394036398,1389.811329117167,1317.2716070319066,415.42478387074067,460.0853763218237,460.06469475098083),c(1508.4552191295866,1359.644881663461,1455.980123532477,578.9329367901169,598.1109892183707,615.7182745596639),c(1011.1622897462064,1035.355571536122,1042.6103723718368,249.06474456323582,285.840276353133,291.6041748314568))
targetgene="Sys1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1425.5730642323567,1389.811329117167,1295.225220303339,360.28831370025335,336.7433392653348,422.6290236577533),c(2171.512458307427,2219.3886340940803,2192.69688004544,1227.2617756913644,1207.968839108788,1297.4415481521237),c(1400.1558700638727,1593.4348494296821,1502.8286953306826,446.7955341401559,527.6298251860914,487.6488734512538),c(1018.8979575366146,1013.8081090691892,1038.0173751367186,117.87797070931772,157.6037140166247,115.26246099756908),c(2396.9519196278925,2274.334663384759,2547.276266596567,2033.3949947357307,2097.793535016315,2157.476834057062),c(1248.7578004515992,1181.8783163112653,1309.9228114557175,343.17699537148144,260.3887448970321,288.64872711357043))
targetgene="Ptpn2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1346.0061955310157,1284.2287630291962,1280.5276291509608,547.5621865207016,646.0773369625608,571.3865587913681),c(1320.589001362532,1293.925121139316,1330.1319992902374,206.286448741306,215.35911232085363,168.46051991952402),c(1749.3660160308686,1948.9679801340737,1882.210266951448,803.2813326562377,877.0989301794766,826.5402117688927),c(1588.0220878309276,1624.6786700067346,1871.1870735871644,404.9678671142689,368.06830105745894,458.09439627238993),c(1235.4966556680424,1098.9205858135738,1228.167460670613,616.0074598357894,632.3726661785065,589.1192450986864),c(931.5954210448656,762.7801713294218,949.8318282224486,399.26409433801166,358.2792504974201,378.29730788945744))
targetgene="Baz1b"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(815.5604041887435,845.7379018271132,797.3443200165236,318.46064667436644,333.8066240973231,349.72797994988906),c(836.5572167627085,691.6735451885436,802.8559166986655,35.17326545358674,37.19839212814745,39.40596957181849),c(1027.7387207256525,1024.5818403026556,942.4830326462595,264.274805299922,317.16523814525715,299.4853687458205),c(1350.4265771255348,1447.9894777778854,1330.1319992902374,585.5873383624171,595.1742740503591,577.2974542271409),c(1090.7291584475472,1093.5337201968407,1068.3311568884988,549.4634441127874,544.2712111381574,564.4905141162998),c(922.7546578558278,812.3393350033673,955.3434249045905,440.14113256785566,398.4143577935792,404.896337350435))
targetgene="Eif2b2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetgenelist=c("Jak2","Jak1","Tsc1","Nprl3","Stat1","Irf1","Ifngr1","Tsc2","Nprl2","Cbfa2t3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3_vs_Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 pos.'


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
targetmat=list(c(943.7514704297927,1080.605242716681,1058.2265629712388,8760.994984331226,8828.744700098994,9158.932477729912),c(1507.350123730957,1591.2801031829888,1451.3871262973587,16249.09801076103,15293.433689948619,16362.343715458332),c(1262.018945235156,1413.513537830793,1256.644043528346,14746.15388421723,14436.891765945224,14708.27814268125),c(2168.197172111538,2385.3040950894633,2516.9624848447866,22556.520072505573,22026.342665143307,21873.268560077147),c(997.9011449626496,1124.7775407738932,1096.8077397462318,14266.086342215573,13920.029896375176,14521.099787215113),c(1746.0507298349794,1640.8392668569343,1771.9783333086107,26971.24020132873,28080.870436527304,27608.807431255325))
targetgene="Jak2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1695.2163414980116,1604.2085806631485,1726.0483609574285,22658.23735368216,23194.176396955936,22779.60586022897),c(1759.3118746185362,1918.8015326803677,2010.8141895347583,21036.464627633,20843.82535749062,21372.81274651505),c(1550.4488442775166,1636.5297743635476,1594.6886400330472,13990.403991363137,13945.481427831277,14258.064940323224),c(1653.2227163500818,1504.0128801919109,1720.5367642752865,14061.701151066352,14924.386483835156,14325.055088595316),c(2665.490101494918,2775.313165740947,2594.1248383947727,18022.02071538101,18585.491393289667,18290.28077675955),c(1369.2131989022403,1223.8958681217841,1502.8286953306826,7572.708989277621,7401.501128445338,7666.431380197287))
targetgene="Jak1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1645.4870485596737,1638.684520610241,1786.6759244609889,5469.918092430759,4981.647830003746,4936.582838109561),c(1706.267295484309,1656.999863707134,1726.966960404452,9299.050882891499,8932.508636035405,9006.234345639115),c(1501.8246467378083,1648.380878720361,1630.5140184669694,9006.25721371029,8414.667861409353,8312.68928117511),c(1161.4552639598503,1146.3250032408262,1173.970093296218,7614.536656303508,7762.71709411077,7996.4563753612665),c(1705.1622000856794,1676.3925799273734,1602.95603505626,7495.7080567981475,7638.396151998277,7922.570182414107),c(1199.0285075132613,1187.2651819279984,1222.655863988471,4629.56223672885,4799.571489587024,4965.152166049129))
targetgene="Tsc1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1675.3246243226765,1708.7137736277728,2107.267131472241,7166.790493367309,7489.602583485686,7177.797357506737),c(1788.0443549829092,1886.4803389799686,1896.9078581038264,5853.021497236042,5636.5353124703415,6014.336105898797),c(2059.8978230458238,2190.2995597637214,2174.324891104967,7688.685702394852,7703.003885694533,7719.629439119241),c(1928.3914706088856,1926.343144543794,1895.9892586568026,6662.006602668537,6646.765330266346,6868.460496367962),c(2240.0283730224705,2252.7872009178263,2002.5467945115454,8298.038760658343,8009.401168223747,8112.703985598131),c(1860.9806512924717,2007.1461287947923,1779.3271288847998,5320.669371452027,5399.640288917402,5536.538724840498))
targetgene="Nprl3"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1666.4838611336386,1637.6071474868943,1609.3862311854255,13291.691826271615,13278.847084692634,13063.078913057829),c(1802.4105951650959,1751.8086985616383,1794.0247200371782,18512.54517413914,19348.058431916692,18795.662336518122),c(1342.6909093351267,1496.4712683284845,1438.5267340390276,4024.0116936495315,4147.62072228844,4090.339641554759),c(1406.7864424556512,1583.7384913195624,1694.8159797586245,15612.176717412298,15600.809877533837,15750.56603785585),c(1404.5762516583916,1579.4289988261758,1466.084717449737,11313.433301706373,11026.386550827705,11339.067744290769),c(1601.2832326144844,1695.785296147613,1636.9442145961348,6262.742508330525,6345.262573017151,6156.197596357343))
targetgene="Stat1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1707.3723908829388,1685.0115649141467,1704.9205736758845,10141.307996185496,10059.228355495872,10090.88365810342),c(2154.936027327981,1910.1825476935946,1999.7909961704745,7597.425337974736,7736.2866575986645,8001.382121557744),c(1758.2067792199066,1840.1532946760628,1845.4662890705022,9445.447717482104,10146.350905480218,9931.289481337553),c(1799.0953089692066,1945.7358607640338,1885.8846647395426,7977.6768563918895,8447.950633313485,8132.4069703840405),c(1929.4965660075152,2182.7579479002948,2044.8023690746331,4443.238992704444,4307.182246417073,4306.087324960465),c(1555.9743212706653,1473.8464327382048,1600.200236715189,8690.648453424054,9048.019432643863,8982.590763896023))
targetgene="Irf1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(2066.5283954376023,1996.372397561326,2142.1739104591393,5299.755537939083,5007.0993614598465,5331.627683067041),c(1291.8565209981589,1361.7996279101542,1388.0037644527272,2991.6288211469587,2949.440933739691,2997.8091351760913),c(971.378855395536,1068.754138359868,956.2620243516142,4796.872904832397,5047.234468756006,5037.068060517698),c(1245.44251425571,1344.561657936608,1246.539449611086,10458.81801406382,10699.43226212241,10466.22551827499),c(1136.0380697913663,1243.2885843420238,1230.923259011684,9993.009904002805,10665.170585162274,10284.958058244625),c(1852.1398881034338,1954.3548457508068,1891.3962614216844,36962.34884773945,37184.68745736339,37299.72049820479))
targetgene="Ifngr1"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1849.9296973061744,1898.3314433367816,2020.0001840049947,6217.112326120467,5753.025014134803,5977.885584044864),c(1813.4615491513932,1868.1649958830756,1901.5008553389446,9494.880414876334,9105.774830948092,9326.407848410141),c(1923.9710890143665,1920.956278927061,2119.208924283548,19103.836285277815,19688.71739140604,18617.350324205643),c(1244.3374188570804,1396.2755678572469,1224.4930628825184,2047.654426676374,2237.77695802487,2004.7787019662655),c(2368.2194392635197,2344.3639164022907,2379.17256779124,21791.26389169105,21923.5576342629,22060.446915543285),c(1888.6080362582152,2080.407501182364,1965.8028166305996,13094.91166549074,14094.274996343865,14135.906434650587))
targetgene="Tsc2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1950.49337858148,1827.2248171959031,1909.7682503621575,5983.257642293917,5943.91150005556,5766.07849759634),c(1049.8406286982472,942.701482928311,1000.3547978087491,1757.712643883294,1577.0160452222508,1615.644752444558),c(1617.8596635939302,1720.5648779845858,1641.537211831253,10225.913959033312,10408.697460489257,10561.78499448665),c(1204.55398450641,1234.6695993552507,1258.4812424223933,2813.8612362869394,2673.3897079465964,2837.229809170931),c(1723.9488218623849,1707.636400504426,1857.4080818818095,4176.112301016393,4197.544880144638,4417.409189000852),c(1664.2736703363792,1549.2625513724697,1704.001974228861,5455.658660490116,5651.2188883104,5519.791187772475))
targetgene="Nprl2"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
targetmat=list(c(1512.8756007241057,1391.9660753638602,1483.5381069431862,3756.885001961481,3710.050162254706,3528.8045751563454),c(2310.7544785347736,2596.4692272654047,2290.98702087697,3166.5445196188493,3209.8296786367227,2958.403165604273),c(1684.1653875117142,1583.7384913195624,1607.5490322913784,4669.48864616265,4553.86632053005,4541.53799315208),c(1775.8883055979823,1821.83795157917,1783.0015266728944,5078.259028461091,5198.964752436607,4936.582838109561),c(1485.2482157583622,1414.5909109541396,1546.9214687878177,2758.724766116452,2761.4911629869457,2763.3436162237713),c(1838.878743319877,2062.092158085471,1980.500407782978,4581.080168130662,4770.204337906908,4480.458740315762))
targetgene="Cbfa2t3"
collabel=c("Input1_DOI_rep1","Input1_DOI_rep2","Input1_DOI_rep3","amCD19_1to2rep1","amCD19_1to2rep2","amCD19_1to2rep3")

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
Sweave("cd19_1to2_v_Input1_DOI_summary.Rnw");
library(tools);

texi2dvi("cd19_1to2_v_Input1_DOI_summary.tex",pdf=TRUE);

