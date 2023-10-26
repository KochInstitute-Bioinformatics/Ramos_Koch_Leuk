pdf(file='invitro.poolD_all_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolD_all_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Mtx2","Bloc1s1","Rabif","Cbx3","Scd2","Gtpbp3","Mbd2","Srrd","Gdi2","Eif3l")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4,5,6,7,8,9,10,11,12,13_vs_0 neg.'


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
targetmat=list(c(281.9217458923781,194.32282497664616,212.62908814971465,200.5822089411275,191.097870448011,222.91898145084696,145.7512096611396,219.33661154539024,227.58732187552314,296.2055343764982,299.7579027159541,471.4798511948517,362.9458671500532,532.7140086157151),c(1132.8693686042986,622.3495125497905,609.0794525062794,662.1160295144015,558.4177292133495,625.7573438696363,598.1197789057137,665.1234544700752,669.8796643883322,699.8262625378803,714.1068011537034,525.4249142606242,581.2393959432011,653.8415558183603),c(401.1166016924644,153.65060579548768,155.01346426398553,150.9235067275474,147.61452268139772,159.55114916025087,161.94578851237733,179.02609915326445,128.8230123823716,186.6203366767681,121.42092261912066,207.68849280322414,215.66348627756784,184.16331156320544),c(893.4431799971688,528.0932585744072,558.3228314640894,472.2445210507128,457.71945017487667,531.8371638675029,521.465439009855,477.79813217725547,601.1740577844007,577.2210413490734,587.3737131699961,1198.659301321465,1158.533728113032,871.3767426312742))
targetgene="Mtx2"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(297.4689009967372,194.96841575729948,189.30847848168145,209.34550933175927,163.6347034375184,184.44565470298505,136.03446235039698,136.34438014983718,65.48503129437222,74.86513506219185,62.98710360866884,64.19462504826927,51.28582905381186,32.13587987008954),c(310.9431020871817,183.9933724861932,133.06465516466014,192.79260859389925,178.5105855682019,174.26153879913923,109.04349759833407,161.2420495685031,121.30833666006657,166.00529948572975,148.7406302084228,188.26827009954604,274.83944287812,201.46724687786906),c(123.34076382791541,70.36939509121073,97.39784037825639,91.52780407993197,68.6579175262315,99.5780221709367,53.98192950412578,69.95059503339472,18.24992675416931,62.93011353053807,42.49732291669223,8.091759459865875,14.465233835690526,21.011921453520085),c(601.156664035218,358.9484740432401,358.03994843274535,387.53261727460557,447.42076254594195,437.9169838653694,306.61735958343445,391.24909086475014,299.51350378901395,239.78543259049852,375.64597935290453,279.9748773113593,253.79910275347922,244.72708516452806))
targetgene="Bloc1s1"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(3.109431020871817,0.6455907806533095,0.0,0.9737000434035316,0.0,1.1315684337606444,2.1592771801650312,1.185603305650758,0.0,0.0,4.553284598217024,0.0,0.0,0.0),c(4.14590802782909,1.291181561306619,2.743601137415673,0.9737000434035316,0.0,1.1315684337606444,5.398192950412578,1.185603305650758,2.14705020637286,0.0,0.0,2.6972531532886252,0.0,0.0),c(4.14590802782909,0.6455907806533095,0.0,0.0,0.0,2.263136867521289,0.0,2.371206611301516,0.0,0.0,0.0,0.0,0.0,0.0),c(9.328293062615451,4.519135464573167,0.0,3.8948001736141262,16.020180756120684,1.1315684337606444,0.0,3.556809916952274,0.0,0.0,6.071046130956033,0.0,0.0,1.2359953796188285))
targetgene="Rabif"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(457.08636006815715,324.7321626686147,286.7063188599378,350.53201562527136,361.59836563815253,259.1291713311876,354.1214575470651,353.3097850839259,244.76372352650603,220.25539735688326,264.84938746295694,175.32145496376063,191.99310363734696,243.49108978490924),c(278.81231487150626,166.56242140855386,193.42388018780494,166.5027074220039,186.52067594626223,196.89290747435214,148.99012543138716,180.21170245891523,215.77854574047242,177.94032101738352,146.46398790931428,120.29749063667268,151.22744464585548,90.22766271217448),c(421.8461418316099,295.03498675856247,260.6421080544889,270.68861206618175,252.88999622161936,221.7874130170863,247.23723712889608,230.00704129624705,162.10229058115092,261.4854717389599,162.40048400307387,810.7942978785608,578.609353427621,430.12639210735233))
targetgene="Cbx3"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(354.47513637938715,166.56242140855386,164.61606824494038,186.95040833347807,224.28253058568956,145.97232795512312,163.02542710245984,197.9957520436766,140.6317885174223,162.75029361346054,176.81921856409446,141.33606523232396,191.99310363734696,176.74733928549247),c(224.9155105097281,72.30616743317067,90.5388375347172,127.55470568586263,108.7083694165332,123.34095927991024,123.07879926940679,125.67395039898035,209.33739512135384,143.22025837984526,110.03771112357809,162.91409045863296,139.39225332574506,131.0155102395958),c(203.1494933636254,194.96841575729948,156.38526483269337,148.97610664074034,162.4904048120812,166.34055976281473,147.91048684130465,149.3860165119955,149.21998934291375,101.99018399776861,105.48442652536107,37.22209351538303,84.16136049856306,30.899884490470715),c(376.2411535254899,94.90184475603651,112.48764663404259,101.26480451396728,138.4601336779002,135.78821205127733,128.47699221981935,106.70429750856822,190.0139432639981,175.77031710253738,179.85474162957246,228.18761676821768,269.5793578469598,201.46724687786906))
targetgene="Scd2"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(283.9946999062926,210.46259449297892,230.46249554291654,196.68740876751338,207.1180512041317,170.8668334978573,197.57386198510036,173.09808262501068,197.52861898630312,124.77522510365309,111.5554726563171,212.54354847914365,114.40684942773416,113.71157492493222),c(604.2660950560899,357.0117017012802,322.37313364634156,370.006016493342,354.7325738855294,358.70719350212426,369.2363978082203,350.9385784726244,333.8663070909797,388.43070075745914,392.34135621303363,291.84279118582924,323.4952294163517,281.8069465530929),c(164.7998441062063,79.40766602035707,74.07723071022318,62.31680277782602,120.15135567090512,88.26233783333026,83.1321714363537,90.10585122945761,67.63208150074509,87.88515855126869,99.41338039440504,97.64056414904823,85.4763817563531,64.27175974017908),c(73.58986749396634,33.5707205939721,19.20520796190971,22.395100998281226,30.896062886804174,30.5523477115374,22.67241039173283,39.124909086475014,33.27927819877933,17.360031318769124,22.007542224715618,80.91759459865875,65.75106288950239,28.427893731233056))
targetgene="Gtpbp3"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(257.0462977254036,298.90853144248234,189.30847848168145,251.21461119811113,239.15841271637305,199.1560443418734,231.04265827765835,206.29497518323188,118.0877613505073,168.1753034005759,168.4715301340299,66.35242757090018,73.64119043624268,34.6078706293272),c(555.5516757290981,340.22634140429415,395.0785637878569,409.9277182728868,395.9273244012683,303.2603402478527,295.8209736826093,366.3514214460842,200.7491942958624,267.99548348349833,207.1744492188746,108.96902739286045,170.9527635127062,203.93923763710671),c(183.45643023143722,191.74046185403293,174.21867222589523,203.5033090713381,162.4904048120812,130.1303698824741,159.7865113322123,151.75722312329702,92.32315887403297,72.6951311473457,37.9440383184752,71.2074832468197,59.17595660055215,93.93564885103098),c(636.3968822717652,524.2197138904874,438.9761819865077,435.2439194013786,572.1493127185959,346.2599407307572,421.0590501321811,480.169338788557,190.0139432639981,245.21044237761387,93.342334263449,197.43893082072736,142.02229584132516,192.81527922053726))
targetgene="Mbd2"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(134.7420109044454,61.331124162064405,82.30803412247019,67.18530299484368,45.77194501748767,64.49940072435673,73.41542412561107,27.268876029967434,39.720428817897904,94.39517029580712,86.51240736612347,28.051432794201702,42.08068024928153,102.58761650836277),c(309.90662508022444,205.94345902840575,159.12886597010905,165.52900737860037,164.77900206295558,207.07702337819794,163.02542710245984,231.1926446018978,157.8081901684052,137.79524859272993,149.49951097479232,134.86265766443125,216.9785075353579,108.76759340645691),c(97.4288386539836,79.40766602035707,116.6030483401661,83.73820373270371,56.07063264642239,52.05214795298964,69.096869765281,98.40507436901291,67.63208150074509,58.5901057008458,72.09367280510288,49.629458020510704,130.18710452121474,143.37546403578412),c(544.150428652568,345.3910676495206,336.09113933341996,407.9803181860797,376.47424776883605,370.0228778397307,346.5639874164875,359.2378016121797,329.57220667823395,317.9055735249596,353.6384371281889,228.72706739887542,224.86863508209817,184.16331156320544))
targetgene="Srrd"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(609.4484800908762,427.3810967924909,443.0915836926312,365.13751627632433,496.62560343974116,448.1010997692152,428.6165202627587,394.8059007817024,404.7189639012841,486.08087692553545,315.6943988097137,170.4663992878411,188.04803986397684,174.27534852625482),c(240.4626656140872,78.11648445905045,123.46205118370528,143.13390638031913,144.18162680508615,121.07782241238895,117.68060631899421,101.96188428596518,104.1319350090837,107.41519378488395,179.09586086320297,220.0958573083518,124.92701949005453,164.3873854893042),c(167.90927512707813,98.77538943995636,105.62864379050342,149.94980668414385,156.76891168489524,89.39390626709091,126.31771503965433,115.00352064812353,83.73495804854153,96.56517421065325,91.82457273070999,148.88837406153212,176.2128485438664,113.71157492493222),c(385.56944658810534,204.65227746709914,193.42388018780494,219.0825097657946,241.44700996724742,238.76093952349598,213.7684408363381,170.72687601370916,138.48473831104945,227.85041105884474,210.96885305072215,90.08825531984007,132.81714703679484,132.25150561921467))
targetgene="Gdi2"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(297.4689009967372,222.72881932539178,296.3089228408927,230.76691028663697,227.71542646200115,220.65584458332566,316.33410689417707,271.5031569940236,205.0432947086081,233.2754208459601,239.8063221727633,230.34541929084858,326.1252719319318,281.8069465530929),c(204.18597037058265,147.8402887696079,176.9622733633109,189.87150846368866,138.4601336779002,161.81428602777214,145.7512096611396,168.35566940240764,142.7788387237952,58.5901057008458,88.03016889886247,36.682642884725304,97.31157307646353,39.551852147802514),c(324.41730317762625,192.38605263468625,253.78310521094974,178.18710794284627,192.2421690734482,183.3140862692244,195.41458480493532,156.49963634590006,162.10229058115092,159.49528774119133,116.10875725453413,66.8918782015579,127.55706200563463,88.99166733255565),c(288.1406079341217,168.4991937505138,167.35966938235606,150.9235067275474,180.79918281907626,119.9462539786283,177.06072877353256,173.09808262501068,132.04358769193087,88.97016050869176,96.37785732892702,57.72121748037658,52.60085031160191,29.663889110851883))
targetgene="Eif3l"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetgenelist=c("Slc43a1","Gm17455","Rbp7","B020004C17Rik","Psmd9","Add1","Atf7ip","Snx31","Sp3","C1qtnf6")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2,3,4,5,6,7,8,9,10,11,12,13_vs_0 pos.'


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
targetmat=list(c(11.401247076529996,3.8735446839198575,13.718005687078366,15.579200694456505,8.010090378060342,15.841958072649021,20.513133211567798,1.185603305650758,5.36762551593215,10.850019574230704,7.588807663695041,3.776154414604075,3.9450637733701432,3.7079861388564854),c(0.0,21.304495761559217,10.974404549662692,2.921100130210595,3.4328958763115747,6.7894106025638665,6.477831540495094,13.041636362158338,10.7352510318643,3.255005872269211,3.0355230654780163,19.959673334335825,0.0,7.415972277712971),c(0.0,3.8735446839198575,2.743601137415673,4.868500217017658,0.0,9.052547470085155,3.238915770247547,0.0,6.441150619118579,10.850019574230704,8.347688430064546,6.4734075678927,14.465233835690526,27.191898351614228),c(1.0364770069572724,0.0,9.602603980954855,4.868500217017658,5.721493127185958,0.0,4.3185543603300625,10.670429750856822,0.0,0.0,9.106569196434048,5.3945063065772505,0.0,9.887963036950628))
targetgene="Slc43a1"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(1.0364770069572724,8.392680148493024,17.833407393201874,13.631800607649442,6.865791752623149,9.052547470085155,15.114940261155219,7.113619833904548,10.7352510318643,11.935021531653772,0.0,0.0,0.0,0.0),c(1.0364770069572724,20.658904980905906,6.859002843539183,19.474000868070632,9.154389003497533,14.710389638888378,12.955663080990188,8.299223139555306,0.0,6.510011744538422,12.90097302828157,18.880772073020378,0.0,6.179976898094143),c(6.218862041743634,9.683861709799643,0.0,8.763300390631784,8.010090378060342,2.263136867521289,0.0,1.185603305650758,2.14705020637286,4.340007829692281,12.142092261912065,0.0,0.0,0.0),c(8.29181605565818,5.164726245226476,0.0,4.868500217017658,16.020180756120684,14.710389638888378,12.955663080990188,8.299223139555306,16.10287654779645,15.190027403922983,3.0355230654780163,1.07890126131545,5.260085031160191,25.9559029719954))
targetgene="Gm17455"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(1.0364770069572724,9.038270929146334,12.346205118370529,7.7896003472282525,3.4328958763115747,9.052547470085155,8.637108720660125,3.556809916952274,4.29410041274572,0.0,8.347688430064546,0.0,1.3150212577900477,2.471990759237657),c(4.14590802782909,37.44426527789195,17.833407393201874,18.5003008246671,30.896062886804174,35.07862144657998,33.46879629255798,22.526462807364403,13.955826341423588,27.125048935576757,18.213138392868096,26.433080902228525,35.50557396033129,7.415972277712971),c(5.182385034786362,9.683861709799643,17.833407393201874,8.763300390631784,20.597375257869448,0.0,5.398192950412578,10.670429750856822,9.66172592867787,16.275029361346053,9.106569196434048,2.1578025226309,0.0,0.0),c(1.0364770069572724,1.291181561306619,9.602603980954855,4.868500217017658,0.0,11.315684337606445,3.238915770247547,5.92801652825379,1.07352510318643,5.425009787115352,14.418734561020578,0.0,0.0,3.7079861388564854))
targetgene="Rbp7"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(1.0364770069572724,6.455907806533095,8.230803412247019,5.84220026042119,10.298687628934724,7.9209790363245105,14.035301671072704,7.113619833904548,7.514675722305009,11.935021531653772,0.7588807663695041,7.0128581985504255,0.0,0.0),c(4.14590802782909,6.455907806533095,1.3718005687078365,9.737000434035316,1.1442986254371916,14.710389638888378,9.71674731074264,15.412842973459854,25.764602476474316,14.105025446499914,23.525303757454626,14.025716397100851,0.0,0.0),c(5.182385034786362,1.9367723419599288,1.3718005687078365,3.8948001736141262,10.298687628934724,4.526273735042578,1.0796385900825156,5.92801652825379,2.14705020637286,0.0,0.0,0.0,0.0,0.0),c(4.14590802782909,27.114812787439,17.833407393201874,8.763300390631784,19.453076632432257,7.9209790363245105,25.911326161980377,9.484826445206064,15.029351444610018,16.275029361346053,20.48978069197661,10.789012613154501,9.205148804530335,17.3039353146636))
targetgene="B020004C17Rik"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(25.91192517393181,7.747089367839715,20.577008530617547,3.8948001736141262,14.875882130683491,9.052547470085155,10.796385900825156,9.484826445206064,15.029351444610018,0.0,6.071046130956033,21.038574595651276,9.205148804530335,0.0),c(6.218862041743634,7.101498587186405,8.230803412247019,2.921100130210595,1.1442986254371916,5.6578421688032225,7.5574701305776095,2.371206611301516,0.0,10.850019574230704,0.0,0.0,0.0,0.0),c(1.0364770069572724,1.291181561306619,5.487202274831346,5.84220026042119,2.2885972508743833,16.973526506409666,5.398192950412578,8.299223139555306,9.66172592867787,9.765017616807633,15.177615327390082,0.0,0.0,16.06793993504477),c(2.072954013914545,11.620634051759572,2.743601137415673,16.552900737860035,4.5771945017487665,2.263136867521289,15.114940261155219,4.742413222603032,0.0,0.0,12.90097302828157,12.407364505127676,34.19055270254124,17.3039353146636))
targetgene="Psmd9"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(1.0364770069572724,7.747089367839715,1.3718005687078365,0.0,2.2885972508743833,0.0,0.0,1.185603305650758,1.07352510318643,3.255005872269211,0.0,4.3156050452618,0.0,0.0),c(1.0364770069572724,3.2279539032665476,5.487202274831346,5.84220026042119,9.154389003497533,3.3947053012819333,7.5574701305776095,0.0,3.2205753095592895,0.0,9.106569196434048,1.07890126131545,17.09527635127062,23.483912212757744),c(27.984879187846357,18.076541858292668,19.20520796190971,9.737000434035316,14.875882130683491,14.710389638888378,10.796385900825156,3.556809916952274,5.36762551593215,10.850019574230704,22.766422991085122,3.776154414604075,0.0,2.471990759237657),c(1.0364770069572724,1.291181561306619,23.32060966803322,1.9474000868070631,3.4328958763115747,6.7894106025638665,11.876024490907671,1.185603305650758,11.80877613505073,29.2950528504229,11.383211495542561,9.71011135183905,52.60085031160191,0.0))
targetgene="Add1"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(17.620109118273632,3.2279539032665476,1.3718005687078365,5.84220026042119,0.0,4.526273735042578,0.0,8.299223139555306,2.14705020637286,3.255005872269211,0.0,10.249561982496775,3.9450637733701432,0.0),c(1.0364770069572724,14.848587955026119,6.859002843539183,14.605500651052973,3.4328958763115747,16.973526506409666,18.353856031402767,15.412842973459854,35.426328405152184,4.340007829692281,24.28418452382413,1.07890126131545,2.6300425155800955,16.06793993504477),c(10.364770069572725,8.392680148493024,5.487202274831346,8.763300390631784,9.154389003497533,9.052547470085155,14.035301671072704,9.484826445206064,0.0,2.1700039148461405,0.0,0.0,1.3150212577900477,6.179976898094143),c(4.14590802782909,12.91181561306619,16.461606824494037,7.7896003472282525,2.2885972508743833,7.9209790363245105,4.3185543603300625,8.299223139555306,1.07352510318643,0.0,0.0,0.0,3.9450637733701432,0.0))
targetgene="Atf7ip"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(7.255339048700907,1.9367723419599288,12.346205118370529,10.710700477438847,10.298687628934724,10.1841159038458,1.0796385900825156,5.92801652825379,0.0,0.0,9.865449962803552,0.0,3.9450637733701432,0.0),c(4.14590802782909,14.848587955026119,12.346205118370529,1.9474000868070631,14.875882130683491,10.1841159038458,11.876024490907671,13.041636362158338,16.10287654779645,14.105025446499914,5.312165364586528,18.341321442362652,26.300425155800955,27.191898351614228),c(1.0364770069572724,10.975043271106262,4.115401706123509,11.68440052084238,10.298687628934724,6.7894106025638665,24.83168757189786,20.155256196062886,0.0,20.615037191038336,30.355230654780165,4.3156050452618,23.670382640220858,0.0),c(2.072954013914545,0.6455907806533095,10.974404549662692,4.868500217017658,8.010090378060342,2.263136867521289,0.0,2.371206611301516,0.0,0.0,0.0,0.0,0.0,0.0))
targetgene="Snx31"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(15.547155104359087,21.950086542212524,27.43601137415673,5.84220026042119,16.020180756120684,10.1841159038458,4.3185543603300625,8.299223139555306,22.54402716691503,15.190027403922983,0.0,1.07890126131545,10.520170062320382,46.96782442551549),c(4.14590802782909,5.810317025879786,5.487202274831346,3.8948001736141262,10.298687628934724,1.1315684337606444,7.5574701305776095,1.185603305650758,0.0,5.425009787115352,0.0,0.0,7.8901275467402865,6.179976898094143),c(1.0364770069572724,12.91181561306619,0.0,0.0,2.2885972508743833,4.526273735042578,2.1592771801650312,2.371206611301516,13.955826341423588,5.425009787115352,4.553284598217024,1.618351891973175,32.875531444751196,7.415972277712971),c(4.14590802782909,14.848587955026119,9.602603980954855,13.631800607649442,8.010090378060342,4.526273735042578,29.15024193222792,5.92801652825379,24.69107737328789,13.020023489076843,36.42627678573619,6.4734075678927,1.3150212577900477,12.359953796188286))
targetgene="Sp3"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
targetmat=list(c(1.0364770069572724,6.455907806533095,12.346205118370529,6.815900303824721,0.0,3.3947053012819333,7.5574701305776095,1.185603305650758,3.2205753095592895,4.340007829692281,0.7588807663695041,8.091759459865875,0.0,14.831944555425942),c(15.547155104359087,1.9367723419599288,5.487202274831346,17.526600781263568,4.5771945017487665,11.315684337606445,11.876024490907671,23.71206611301516,1.07352510318643,5.425009787115352,12.142092261912065,16.722969550389475,2.6300425155800955,1.2359953796188285),c(1.0364770069572724,0.6455907806533095,4.115401706123509,2.921100130210595,1.1442986254371916,9.052547470085155,1.0796385900825156,8.299223139555306,0.0,2.1700039148461405,3.7944038318475206,0.0,0.0,0.0),c(1.0364770069572724,1.291181561306619,4.115401706123509,0.0,0.0,7.9209790363245105,10.796385900825156,0.0,6.441150619118579,0.0,27.319707589302148,7.0128581985504255,23.670382640220858,40.78784752742134))
targetgene="C1qtnf6"
collabel=c("In.vitro.input.Pool.D..25.32._1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._no.CAR_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._hEGFRv3.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.10_rep3_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep1_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep2_1400..1401..1402..1403..1427..1428..1429..1430","In.vitro_Pool.D..25.32._mCD19.CAR_ET.1.2_rep3_1400..1401..1402..1403..1427..1428..1429..1430")

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
Sweave("invitro.poolD_all_v_input_summary.Rnw");
library(tools);

texi2dvi("invitro.poolD_all_v_input_summary.tex",pdf=TRUE);

