pdf(file='invitro.poolC_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolC_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Senp3","Ptbp1","Kctd5","Foxf2","Ecpas","Csrnp1","Commd3","Fhod1","Slc25a22","Dot1l")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='8,9,10_vs_2,3,4 neg.'


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
targetmat=list(c(874.7718010501738,1102.9756778300402,1092.714176613208,730.5567873729357,677.9940697230943,598.5302520635298),c(340.53938290109755,328.23783847896163,365.2171934379271,230.75811187762307,200.74135010669264,178.84012036132495),c(402.70133374812326,283.5878383917499,447.46449705397504,297.75240242273946,213.53369104486424,214.78788324802343),c(572.971025198672,586.4837849293211,673.6445819981068,555.0955502309642,359.1695724948177,442.1574835063914))
targetgene="Senp3"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(1563.958647397633,1508.4459488922869,1577.3857872077763,1327.1249936556387,1286.1222773992513,1276.145582477796),c(298.197184498051,366.8540547706042,343.6762329670574,155.25660983471414,195.82121897662665,211.1931069593536),c(528.8270311189001,473.0486495726212,562.0232413763275,288.1817894877228,173.18861577832305,214.78788324802343),c(180.17956767253838,345.13243310655525,234.0131614789935,316.8936282927727,301.11202516003897,277.6964682997458))
targetgene="Ptbp1"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(783.781119375542,743.362163614119,724.5595794747079,672.0697083256118,688.8183582092395,710.8670110844625),c(737.8353296190446,791.6324339786721,856.7427460004992,466.8332309413665,587.4636569298799,534.7229729396399),c(1022.5190465416553,1117.4567589394062,929.19870394797,804.9948879786206,712.4349876335563,890.6058255179549),c(1182.8788617702144,1205.5500023547156,1286.5828208510354,892.1938058309943,930.8888098084865,942.7300817036678))
targetgene="Kctd5"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(197.29662660142952,267.90000052327014,181.13989486867698,57.42367761009975,104.30677995739912,43.13731546403818),c(318.9178347803929,226.87027071339995,213.45133557498153,221.18749894260645,208.61355991479823,151.8792981963011),c(272.07114718553294,278.7608113552946,319.19786879561457,82.94531210347742,122.0192520256367,108.74198273226291),c(198.1975244397922,257.03918969124567,174.28595290067298,140.36898971357718,200.74135010669264,211.1931069593536))
targetgene="Foxf2"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(480.17854784731475,453.7405414267999,476.83853405970643,453.0090122574536,378.8500970150817,409.8044969083627),c(1024.3208422183807,903.8608125762584,878.2837064713689,626.3434465249769,605.1761289981175,587.7459231975203),c(563.0611489766824,567.1756767834999,590.4181438152011,479.5940481880553,411.3229624735173,358.57893479481737),c(631.529384692247,685.4378391766552,607.0634314517823,395.5853346473538,388.69035927521367,387.33714510417616))
targetgene="Ecpas"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(183.78315902598914,171.3594597941638,156.66153069723413,145.6859968996975,203.69342878473225,162.66362706231064),c(361.2600331834394,434.4324332809786,323.1144070630454,223.31430181705457,205.66148123675865,177.04273221699003),c(587.3853906124751,442.8797305947754,521.8787241351612,530.6373171748106,556.9588439234708,672.2231659812617),c(802.6999739811585,1005.22838034182,1004.592065596014,555.0955502309642,593.3678142859592,640.7688734554005))
targetgene="Csrnp1"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(255.8549860950045,289.6216221873191,286.88642808931,147.81279977414565,145.63588144995347,186.92836701083212),c(230.62984662084912,201.52837877200952,205.61825904011982,185.0318500769881,89.54638656720113,149.18321597979872),c(191.89123957125338,222.04324367694463,168.4111454995267,80.81850922902927,69.86586204693714,128.51325231994707),c(196.39572876306684,207.56216256757867,153.724126996661,215.8704917564861,99.38664882733312,164.46101520664556))
targetgene="Commd3"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(259.4585774484553,259.4527032094733,338.7805601327688,284.9915851760506,361.13762494684414,361.2750170113198),c(509.9081765132836,470.6351360543935,403.4034415453779,356.23948147006325,227.31005820904903,381.0462865990039),c(655.8536263280397,464.60135225882436,562.0232413763275,312.6400225438764,272.5752646056562,507.7621507746161),c(448.64712350462054,438.0527035583201,501.31689823114925,311.57662110665234,282.4155268657882,272.304303866741))
targetgene="Fhod1"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(136.93647143112918,168.94594627593614,163.51547266523812,91.45252360126997,114.14704221753111,112.33675902093276),c(142.34185846130532,89.30000017442339,132.1831665257913,59.55048048454789,107.25885863543871,68.30074948472712),c(193.69303524797874,234.11081126808293,231.07575777842035,131.8617782157846,206.64550746277183,113.23545309310023),c(232.4316422975745,225.66351395428612,172.32768376695756,143.55919402524938,110.21093731347831,119.52631159827246))
targetgene="Slc25a22"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(227.92715310576105,189.46081118087125,270.24114045272887,176.52463857919554,160.39627484015148,193.21922551600434),c(200.90021795488028,173.77297331239146,248.7001799818592,197.7926673236769,123.0032782516499,116.83022938177007),c(149.54904116820686,114.6418921158138,155.68239613037642,47.853064675083125,89.54638656720113,70.99683170122951),c(257.65678177172987,213.59594636314782,234.99229604585122,155.25660983471414,77.73807185504273,70.09813762906205))
targetgene="Dot1l"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetgenelist=c("Ints13","Ddx6","Rlim","Gatad2b","Ubr4","Zfp341","Hexim1","Gskip","Ccnc","Mettl5")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='8,9,10_vs_2,3,4 pos.'


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
targetmat=list(c(680.1778679638323,564.7621632652722,706.9351572712691,2117.2322615131225,1951.324006184174,2085.8689415006797),c(832.4296026471274,794.0459474968998,822.4730361604793,2064.0621896519187,1688.58900383865,1512.5021234578387),c(324.32322181056907,340.30540607009993,450.40190075454814,794.3608736063799,675.0419910450547,689.2983533524434),c(673.8715830952935,813.354055642721,648.1870832598063,1454.7331661225269,1249.7133070367631,1393.8745059317337))
targetgene="Ints13"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(645.0428522676874,621.4797309436221,622.7295845215057,2859.4864646955225,2849.739950534225,2816.5072221728265),c(832.4296026471274,937.6500018314455,730.4343868758542,3442.2304522943127,2263.260319830358,3344.0406425351266),c(328.82771100238256,358.4067574568074,258.49152565043636,1079.3524587824304,1225.112651386433,1110.7858731989832))
targetgene="Ddx6"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(1326.1216180698825,1472.243246118872,1551.9282884694758,2626.6015499434516,2719.8484887004824,2851.5562909873574),c(1162.1582114878724,1059.5324345019424,916.4699545788197,1721.6469268657684,1742.7104462693758,1750.6560525822163),c(72.97272490737804,131.5364867434074,118.47528258978332,160.5736170208345,109.22691108746511,118.627617526105),c(1103.5998519942975,997.987839787137,1023.1956223663104,2049.174569530782,2013.3176584230057,2049.022484541814))
targetgene="Rlim"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(672.9706852569309,632.3405417756467,494.46295626314526,851.7845512164796,850.1986592754041,888.80843737362),c(533.3315203107136,640.7878390894434,567.8980487774737,1099.5570860896878,1013.5470127935952,1051.4720644359306),c(509.00727867492094,395.81621698933606,377.94594280707736,749.698013242969,734.0835646058466,900.491460311797),c(366.6654202136156,307.7229735740265,402.4243069785202,688.0207298839729,680.9461484011339,880.7201907241129))
targetgene="Gatad2b"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(617.1150192784439,611.8256768707116,520.8995895683036,828.3897195975501,1018.4671439236612,938.2366113428304),c(491.89021974602974,401.8500007849052,423.9652674493899,694.4011385073173,632.728863326487,671.3244719090942),c(372.07080724379176,438.0527035583201,419.0695946151013,580.6171847243419,623.8726272923683,612.0106631460417),c(311.7106520734914,310.1364870922542,275.13681328701745,460.4528223180221,564.8310537315764,457.43528273323824))
targetgene="Ubr4"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(324.32322181056907,279.96756811440844,357.3841169030654,717.7959701262469,401.4827002133853,553.5955484551566),c(639.6374652375113,613.0324336298254,556.1484339751812,1014.4849711117622,1130.646133689166,1048.7759822194282),c(645.9437501060501,592.5175687248903,653.0827560940949,903.891221640459,918.0964688703149,764.7886554145102),c(672.9706852569309,704.7459473224765,665.8115054632451,683.7671241350766,869.8791837956682,870.8345559302708))
targetgene="Zfp341"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(736.934431780682,588.8972984475488,810.7234213581867,1195.263215439854,1152.2947106614563,945.4261639201701),c(570.2683316835839,693.8851364904519,429.84007485053615,846.4675440303592,926.9527049044336,823.2037701053953),c(393.69235536449634,458.5675684632552,445.50622792025956,471.08683669026277,399.5146477613589,424.18360206304214),c(621.6195084702574,482.7027036455318,685.3941968003994,966.6319064366791,912.1923115142357,1257.2730069622794))
targetgene="Hexim1"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(254.05319041827912,269.10675728238397,214.43047014183924,161.63701845805855,150.55601258001948,174.34665000048764),c(463.0614889184236,439.25946031743393,539.5031463386001,786.9170635458114,806.9015053308234,764.7886554145102),c(592.7907776426513,452.53378466768606,524.8161278357344,772.0294434246744,660.2815976548567,983.1713149512035),c(699.9976204078116,724.0540554682977,656.020159794668,1146.3467493275468,1337.2916411519377,1286.9299113438058))
targetgene="Gskip"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(308.10706072004064,288.41486542820525,309.4065231270374,585.9341919104622,704.5627778254507,733.3343628886491),c(440.5390429593563,477.8756766090765,544.3988191728887,702.9083500051099,619.9365223883156,833.9880989714048),c(227.02625526739834,260.65945996858716,201.70172077268896,550.841944482068,455.6041426441112,603.9224164965345),c(109.9095362802484,127.91621646606593,125.32922455778731,197.7926673236769,164.33237974420427,157.27146262930586))
targetgene="Ccnc"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(890.9879621407023,1008.8486506191615,995.7798544942946,1127.2055234575137,1007.642855437516,1231.2108788694231),c(405.40402726321133,526.1459469736296,394.5912304436585,613.5826292782881,655.3614665247907,701.8800703627879),c(605.403347379729,481.495946886418,487.60901429514126,647.6114752694583,802.9654004267705,868.1384737137683),c(706.3039052763504,529.7662172509712,558.1067031088967,649.7382781439064,791.1570857146122,876.2267203632756))
targetgene="Mettl5"
collabel=c("In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
Sweave("invitro.poolC_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10_summary.Rnw");
library(tools);

texi2dvi("invitro.poolC_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10_summary.tex",pdf=TRUE);

