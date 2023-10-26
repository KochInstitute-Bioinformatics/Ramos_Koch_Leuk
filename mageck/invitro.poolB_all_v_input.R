pdf(file='invitro.poolB_all_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolB_all_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Dph6","Mcat","Alg3","Dnajc11","Vps39","Ippk","Tpk1","Tecr","Alad","Ube2m")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4,5,6,7,8,9,10,12,13_vs_0 neg.'


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
targetmat=list(c(489.4739455134216,257.4027521656195,256.33376774275195,253.61992972573088,293.04175585089445,345.77241796410254,328.417996014739,265.21512826450964,379.21554870248633,290.47107897941595,270.0951089779086,253.4114313942133,210.08547990971354),c(406.0882648467399,181.03050701757854,143.9471908945784,156.32607737404965,138.60083047001766,179.61175944588393,189.37349380590146,207.08578508324726,161.28661311903215,199.36058526226995,368.0707857640126,164.76957267607287,372.99074611088975),c(247.65547158004466,71.42219222177903,87.75390247049165,75.43006530635961,146.52087792544722,122.64239081106612,95.53978065883317,96.27672464396584,106.56436938221768,65.85213902328374,63.55179034774319,95.9417764949285,372.1005533994079),c(465.2920981200839,186.68771036187786,200.14047931866517,193.4945153510964,174.24104401945075,223.13002715303642,209.84630394708,188.92036533910277,172.80708548467732,187.63349201154819,103.27165931508269,306.5965466250976,256.3755009067691))
targetgene="Dph6"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(361.0599972867318,235.4810892064596,181.66597353540377,233.94252138494141,213.84128129659865,199.39279022186233,199.60989887649072,219.8015789041484,353.29448587978476,318.4356859619063,262.15113518444065,364.9958888394019,466.4609808164826),c(378.5709902267349,106.07256270561243,156.26352808341935,145.39418385138882,144.93686843436132,102.86136003508771,98.09888192648049,141.69027400432708,93.123818288965,147.0397076821267,164.17545839833656,203.35485235338106,83.67811487929269),c(183.44849746669976,40.307573828132725,41.56763801233815,44.82076344290934,59.400355915721846,59.34309232793522,66.53663295883024,74.47822095099244,31.681299005524174,47.81045709909641,42.36786023182879,8.342763173472044,40.05867201668267),c(762.1451212934708,367.7182173794564,308.6782007953259,297.3475038163741,277.2016609400353,307.0015976431849,346.3317048882702,259.7655023412663,317.7730294190455,443.8253753350082,370.71877702850196,217.95468790695713,97.03100555152024))
targetgene="Mcat"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(200.12563360003608,118.09411981224851,130.86108263143493,102.75979911301165,163.9449823273923,121.85114958002698,116.86562455589414,63.57896910450574,90.24370019755371,55.929213964980704,66.19978161223248,80.29909554466842,67.65464607261961),c(443.61182114674665,239.01684129664667,284.045526417644,207.70597693055547,254.2335233192895,269.81325978434546,243.96765418237757,247.04970852036516,239.04980158713695,175.00431466461706,135.0475544889543,280.5254117079975,144.21121926005762),c(130.08166184002346,107.48686354168726,47.72580660675862,42.63438473837717,62.568374897893676,37.97957908987854,68.24270047059512,61.76242713009129,86.40354274233866,39.69170023321212,45.01585149631809,8.342763173472044,22.254817787045926),c(273.50503258671597,145.67298611570774,135.47970907725028,126.80996486286544,103.7526216661275,109.1912898834008,105.77618572942244,116.25868636252478,84.48346401473113,85.6979891398898,68.84777287672179,63.61356919772433,36.49790117075532))
targetgene="Alg3"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(268.5018917467151,161.2302953125309,155.4937570091168,135.5554796809941,136.22481623338876,186.73293052523616,133.92629967354293,139.87373202991265,135.36555029633055,168.68972599115148,198.59934483669747,114.7129936352406,92.58004199411106),c(339.37972031339456,213.55942624729968,170.88917849516795,162.88521348764613,156.02493487196273,226.29499207719297,133.07326591766048,174.38802954378718,208.32854194541653,141.6272031048705,150.93550207589007,95.9417764949285,239.46183938861418),c(220.97205376670652,38.89327299205789,76.20733635595327,86.36195882902042,107.71264539384228,69.62922833144398,92.98067939118586,94.46018266955139,87.36358210614242,110.95634383375204,121.80759816650777,65.69925999109235,80.11734403336534),c(107.5675280600194,43.136175500282384,59.27237272129699,44.82076344290934,37.22422304051902,18.98978954493927,34.97438399118,18.165419744144497,22.080905367486547,54.125045772561975,39.719868967339494,41.71381586736022,16.02346880667307))
targetgene="Dnajc11"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(391.91269913340403,284.2744680510413,309.4479718696284,370.59119041820156,363.5301782042177,189.8978954493927,332.6831647941512,219.8015789041484,228.48936858529555,231.83561272580712,251.55917012648345,203.35485235338106,210.08547990971354),c(155.9312228466948,94.75815601701377,72.35848098444048,83.08239077222218,79.2004745542958,72.0029520245614,89.5685443676561,67.21205305333464,60.48247991963706,56.831298061190076,13.239956322446497,56.313651420936296,71.21541691854696),c(253.4924692267124,127.99422566477233,114.6958900710812,124.62358615833328,144.14486368881836,119.47742588690957,128.80809713824829,108.99251846486699,112.32460556504026,124.48760527689254,71.49576414121108,68.82779618114436,73.88599505299248),c(210.9657720867047,135.06572984514648,97.76092643642491,103.85298846527773,84.7445077730965,174.07307082860999,156.10517732648634,63.57896910450574,57.602361828225774,115.46676431479888,108.56764184406128,32.328207297204166,61.42329709224676))
targetgene="Vps39"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(205.9626312467038,113.14406688598659,97.76092643642491,79.80282271542394,127.51276403241623,124.2248732731444,116.86562455589414,72.66167897657799,60.48247991963706,103.73967106407711,108.56764184406128,100.11315808166452,139.76025570264844),c(256.82789645337965,144.96583569767031,143.17741982027584,169.4443496012426,175.82505351053666,122.64239081106612,146.7218060117795,161.67223572288603,121.92499920307789,132.60636214277685,103.27165931508269,113.67014823855659,149.55237552894863),c(196.79020637336882,99.70820894327568,75.43756528165072,74.33687595409353,71.28042709886621,62.50805725209176,46.063822817651705,43.597007385946796,50.882086281599435,49.614625291515146,45.01585149631809,30.242516503836157,32.937130324827976),c(185.1162110800334,120.21557106636075,86.21436032188653,123.5303968060672,91.87255048298313,72.0029520245614,69.09573422647756,119.89177031135368,74.88307037669351,55.929213964980704,105.91965057957198,75.08486856124838,86.3486930137382))
targetgene="Ippk"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(240.1507603200433,131.52997775495942,142.4076487459733,170.5375389535087,140.18483996110356,122.64239081106612,136.48540094119025,107.17597649045254,130.56535347731176,123.58552118068317,185.35938851425098,63.61356919772433,72.10560963002881),c(146.7587979733598,91.22240392682669,77.7468785045584,95.10747364714908,49.89629896920635,36.397096627800266,76.77303802941951,99.90980859279473,46.08188946258062,103.73967106407711,60.90379908325389,39.62812507399221,67.65464607261961),c(64.20697411334491,24.750264631309566,19.244276857563957,17.491029636257302,22.176132875202825,9.494894772469635,46.063822817651705,29.064671590631196,8.640354274233866,16.237513731768594,2.6479912644892996,29.199671107152152,5.341156268891023),c(502.8156544200907,286.39591930515354,317.145682612654,270.0177700097221,219.38531451539936,323.6176634950067,284.9132744647346,228.88428877622067,268.8110218650536,282.3523221135316,278.03908277137646,324.3249183687257,125.51717231893903))
targetgene="Tpk1"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(198.45791998670245,132.94427859103425,68.50962561292769,77.61644401089178,71.28042709886621,97.32267141781375,80.18517305294927,78.11130489982133,53.76220437301072,84.79590504368042,60.90379908325389,216.91184251027312,126.40736503042086),c(375.23556300006766,262.3528050918814,202.44979254157283,183.65581118070168,178.9930724927085,179.61175944588393,169.75371742060537,199.81961718558946,104.64429065461016,91.11049371714599,105.91965057957198,79.25625014798442,131.7485212993119),c(85.88725108668216,61.522086369255206,40.02809586373303,43.72757409064326,48.312289478120434,28.484684317408906,24.73797892059073,18.165419744144497,26.88110218650536,46.908373002887046,42.36786023182879,27.11398031378414,8.011734403336535),c(319.3671569533909,166.18034823879282,121.62382973980421,116.97126069247071,133.05679725121695,92.57522403157894,138.19146845295512,101.72635056720918,112.32460556504026,101.03341877544902,121.80759816650777,102.19884887503254,82.78792216781085))
targetgene="Tecr"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(159.26665007336206,65.76498887747971,78.51664957886095,69.96411854502921,86.32851726418242,56.178127403778674,76.77303802941951,38.14738146270344,82.56338528712361,45.10420481046831,37.071877702850195,95.9417764949285,118.39563062708433),c(179.27921343336567,84.85805016448994,96.99115536212234,99.48023105621341,87.91252675526833,115.52121973171388,86.15640934412635,98.09326661838028,101.76417256319887,110.95634383375204,60.90379908325389,195.01208917990903,90.79965657114738),c(245.1539011600442,134.35857942710908,128.55176940852724,94.014284294883,83.1604982820106,98.90515387989203,129.66113089413074,83.5609308230647,114.24468429264778,126.29177346931128,45.01585149631809,63.61356919772433,58.75271895780125),c(110.90295528668666,57.986334279068124,34.63969834361512,52.473088908771906,34.84820880389015,30.858408010526315,38.386519014709755,52.679717258019046,46.08188946258062,50.51670938772451,5.295982528978599,29.199671107152152,72.10560963002881))
targetgene="Alad"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(406.0882648467399,240.4311421327215,203.2195636158754,257.99268713479523,225.72135247974302,212.0526499184885,190.2265275617839,207.08578508324726,131.52539284111552,117.27093250721761,217.13528368812257,53.18511523088428,99.70158368596576),c(260.9971804867137,207.90222290300036,208.6079611359933,195.68089405562856,198.0011863857395,169.32562344237516,201.31596638825562,207.08578508324726,171.84704612087356,119.97718479584572,119.15960690201848,153.2982733125488,124.62697960745719),c(276.00660300671643,161.9374457305683,160.8821545292347,155.23288802178357,143.35285894327538,99.69639511093116,110.89438826471707,134.42410610666929,54.722243736814484,55.02712986877134,42.36786023182879,21.899753330364113,10.682312537782046),c(424.43311459340987,275.0815126165549,237.8592619594905,212.07873433961979,210.6732623144268,212.0526499184885,305.38608460591314,234.33391469946403,147.84606202577947,127.19385756552064,177.41541472078308,31.285361900520165,154.0033390863578))
targetgene="Ube2m"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetgenelist=c("Ralgapb","Zfp36l2","Edc4","Cnot6l","Dcp2","E130309D02Rik","Dcp1a","Zfp36l1","Sms","Edc3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4,5,6,7,8,9,10,12,13_vs_0 pos.'


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
targetmat=list(c(266.00032132671464,376.2040223959054,458.0137892100222,396.82773487258754,445.89867174068536,444.67757184399454,500.7308147029917,390.5565244991067,624.9856258362496,526.8171121862699,574.614104394178,848.8761529007804,738.8599505299248),c(292.6837391400528,490.76239011796685,532.6815834173703,540.0355400194442,544.1072601880121,515.0980414064777,614.1843042353561,425.07082201298124,757.4710580411689,837.1340412822918,627.573929683964,870.7759062311445,1056.6587485289406),c(400.25126720007216,627.9495712172255,566.551510686683,596.8813863372804,549.6512934068128,666.2251165349527,613.3312704794737,721.1671638425365,795.8726325933194,771.2819022590081,924.1489513067655,1355.699015689207,1372.6771611049928),c(867.2110789334897,1388.8434210254854,1446.3998486145072,1517.346820945321,1425.6085419773244,1491.4897205087718,1531.195591808978,1649.4201127683205,2175.449198379327,1915.1245362524846,2332.880304015073,3507.08906904831,2747.1347076329494))
targetgene="Ralgapb"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(267.6680349400483,424.2902508224497,401.05072971163287,391.3617881112571,434.81060530308395,424.105299836977,390.68946019415705,445.0527837315402,590.4242087393142,550.2712986877134,807.6373356692363,860.3474522643045,1339.7400307801647),c(427.76854182007713,663.3070921190964,527.2931858972524,515.9853742695905,590.8355401750466,725.5682088628879,769.43644780596,650.322026840373,907.2371987945559,910.2028530752506,799.6933618757685,1952.2065825924583,1632.6134328576893),c(230.97833544670831,249.6240975672079,274.03850245171077,266.7382019529239,371.4502256596473,295.1329791775978,391.5424939500395,517.7144627081182,538.582083093911,501.5587574924076,579.9100869231567,1213.8720417401823,1363.7752339901745),c(576.1950534067706,1011.2250977935051,1034.5723238626383,970.7521448122802,922.685528557546,993.007744954116,973.3115154618629,950.0514526187573,1422.7783371571766,1420.7824515297518,1895.9617453743385,3496.66061508147,3323.9795846731795))
targetgene="Zfp36l2"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(706.2767152467941,1010.5179473754678,1002.2419387419309,1112.866760606871,959.1177468525221,1148.0910262377868,1201.0715282824742,1008.1807958000196,2476.9015586137084,2101.8559441678235,2764.5028801268286,4029.554612786997,3620.4137575966315),c(213.46734250670517,360.6467131990822,345.6272123618487,327.9568056798244,266.1135945024339,270.6045010153846,350.5968736676824,310.62867762487093,669.1474365712227,707.2339314281431,579.9100869231567,1528.811351538752,1725.1934748518004),c(210.1319152800379,262.3528050918814,279.4268999718287,249.24717231666656,324.7219456726128,390.0819269022942,311.35732089709023,212.53541100649062,723.8696803080372,495.24416881894206,529.5982528978599,1162.772617302666,1240.038447094199),c(478.633807026753,620.8780670368515,603.5005222532058,549.8742441898389,571.8274262820156,582.3535460448043,675.6027346588917,708.4513700216354,979.2401510798381,866.0007323609916,1075.0844533826555,1165.901153492718,1216.8934365956713))
targetgene="Edc4"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(710.4459992801281,953.9459139324745,1100.7726362526585,1059.300482345833,1087.4225156304813,1173.410745631039,1142.2121991265858,1113.5402303160577,1269.1720389485745,1434.3137129728923,1403.4353701793289,1555.9253318525361,1993.1414810078334),c(599.5430439934414,712.100470963678,768.2315321539531,708.3867002684208,711.2202614975763,719.2382790145748,747.2575701530166,942.7852847210994,1323.894282685389,982.3695807719998,834.1172483141294,1625.7959734303645,1266.744228438654),c(632.8973162601142,973.0389752194847,758.22450818802,921.5586239603066,935.3576044862334,793.6149547322536,868.388363488323,886.4724835142515,1250.931291036303,1118.584279299614,1128.0442786724416,1415.1412033001955,1469.708166656513),c(732.9601330601322,1063.554228728274,1062.2840825375304,1073.511943925292,1070.7904159740792,1069.7581443649121,980.9888192648049,930.0694909001983,1533.1828639946093,1796.0494355528483,1456.3951954691147,1856.2648060975298,1813.3225532885021))
targetgene="Cnot6l"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(507.8187952600916,700.078913857042,850.5970371043269,707.2935109161547,844.2770587487931,876.695283991363,836.8261145206727,806.5446366400157,1462.1399510731308,1498.3616838037574,1331.9396060381177,2943.9525548389474,3128.1371881471755),c(152.59579562002753,214.9737270833745,263.26170741147496,320.3044802139618,183.74510096596626,284.0556019430499,192.7856288294312,250.68279246919406,426.2574775288707,437.5107866615426,323.05493426769453,606.9360208700912,972.0904409381661),c(479.4676638334198,693.7145600947052,626.5936544822824,619.8383627348682,752.40450826581,699.4572482385964,783.9380216559615,904.637903258396,1423.7383765209804,1305.315687214953,1827.1139724976167,2658.21291614753,2639.421389543647),c(417.76226014007534,417.92589706011296,400.28095863733034,364.0320543046051,371.4502256596473,422.52281737489875,320.74069221179707,405.0888602944223,365.77499760923365,412.25243196768037,471.34244507909534,685.1494256213916,582.1860333091215))
targetgene="Dcp2"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(348.55214518672955,454.6977187980586,402.590271860238,426.34384738377173,363.5301782042177,447.05129553711197,416.28047287063026,457.76857755244134,601.9446811049593,734.2964543144242,905.6130124553405,1371.341696639467,1801.7500480392382),c(262.6648941000474,207.19507248496294,314.0665983154438,278.7632848278507,264.52958501134793,260.3183650118758,401.7788990206288,457.76857755244134,559.7029490975938,539.446289533201,775.8614404953647,1048.0596236674255,938.2631179018563),c(228.47676502670788,327.4106435513237,301.7502611266029,320.3044802139618,365.90619244084655,301.4629090259109,400.92586526474634,425.07082201298124,553.9427129147712,705.4297632357244,638.1658947419212,752.9343764058519,1094.0468424111777),c(617.0540369334446,795.5442202920932,655.8449553057796,816.6124461427628,730.2283753906072,921.7960341605938,743.8454351294869,955.5010785420005,2055.4442779038563,1711.2535305091678,2142.2249329718434,5440.524434500457,4921.875501783077))
targetgene="E130309D02Rik"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(932.2519098535015,1260.1420449426757,1362.4948015155283,1429.8916727640344,1421.6485182496097,1516.8094399020242,1416.0360347648489,1505.913296789579,2248.4121900284126,2143.3518125934543,2478.5198235619846,4141.139070232186,3667.593971305169),c(587.0351918934392,478.03368259329335,576.5585346526162,426.34384738377173,455.4027286872008,481.86590970283396,522.0566586000526,599.4588515567684,463.69901271721744,651.3047174631624,783.8054142888327,745.6344586290639,766.4559245858618),c(937.2550506935023,1288.4280616641722,1186.2172255002424,1175.1785536860375,1218.1032986450693,1197.1479825622132,1124.2984902530545,1220.7162068065102,1845.1956572308322,1627.3597095616967,1461.6911779980933,3232.8207297204167,2470.284774362098),c(409.42369207340715,595.4206519875045,518.8257040799243,628.5838775529968,492.62695172771987,626.6630549829958,621.0085742824156,675.7536144821753,1110.7655439209536,1082.5009154512395,897.6690386618726,2355.7877511091683,2565.5353944906547))
targetgene="Dcp1a"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(364.39542451339906,485.8123371917049,481.1069214390989,389.17540940672495,539.3552317147544,597.3871294345479,474.2867682706361,588.5595997102818,697.9486174853356,548.4671304952947,585.2060694521352,881.2043601979846,932.0317689214835),c(218.47048334670606,183.8591086897282,338.69927269312564,294.0679357595759,313.63387923501136,286.42932563616733,267.85259934708586,465.03474545009914,474.2594457190589,465.475393644033,606.3899995680496,617.3644748369312,765.5657318743799),c(394.41426955340444,372.66827030571835,351.0156098819666,391.3617881112571,372.24223040519024,386.1257207470985,383.0121563912151,454.1354936036124,378.25550933868254,405.03575919800545,360.12681197054474,817.5907910002603,551.919481118739),c(528.665215426762,632.8996241434875,791.3246643830299,624.2111201439325,582.9154927196171,625.8718137519568,682.4270047059512,730.2498737146088,964.8395606227817,961.6216465591845,738.7895627925146,1390.1129137797793,1387.810437200184))
targetgene="Zfp36l1"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(304.3577344333882,439.84756001927286,468.0208131759555,411.0391964520466,449.0666907228572,503.22942294089063,526.3218273794648,608.5415614288406,660.5070822969889,674.7589039646059,680.53375497375,1108.5446566750977,633.817210575068),c(526.1636450067616,732.6078330867632,862.9133742931679,811.1464993814324,722.3083279351777,918.6310692364372,786.4971229236088,690.2859502774909,1411.2578647915313,1471.2991609174765,1080.3804359116343,2197.2752508131994,2781.8522233807407))
targetgene="Sms"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
targetmat=list(c(268.5018917467151,374.7897215598306,431.8415726837352,399.0141135771197,401.5464059902797,460.5023964647773,435.0472155000439,452.318951629198,689.3082632111017,460.96497316298615,585.2060694521352,678.8923532412875,622.2447053258041),c(1174.0703837868784,1424.9080923453937,1631.9146775214235,1572.006288558625,1428.7765609594962,1746.2693969033737,1536.3137943442728,1449.600495582731,2260.8927017578617,2368.872836645796,2687.711133456639,3095.165137358128,3457.508491395455),c(264.332607713381,360.6467131990822,291.7432371606696,308.279397339035,365.1141876953036,386.1257207470985,447.8427218382805,339.6933492155021,395.5362178871503,255.28979922725065,638.1658947419212,514.1227805652147,381.00248051422625),c(150.09422520002707,151.33018946000706,138.5587933744605,214.26511304415195,167.9050060551071,181.98548313900133,134.77933342942538,203.45270113441836,222.729132402473,237.2481173030633,195.95135357220818,514.1227805652147,472.6923297968555))
targetgene="Edc3"
collabel=c("In.vitro.input_Pool.B..9.16._1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._no.CAR_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep3_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._hEGFRv3.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep2_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.10_rep3_1410..1411..1412..1413..1414..1415..1416..1417.x","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep1_1410..1411..1412..1413..1414..1415..1416..1417","In.vitro_Pool.B..9.16._mCD19.CAR_ET.1.2_rep2_1410..1411..1412..1413..1414..1415..1416..1417")

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
Sweave("invitro.poolB_all_v_input_summary.Rnw");
library(tools);

texi2dvi("invitro.poolB_all_v_input_summary.tex",pdf=TRUE);

