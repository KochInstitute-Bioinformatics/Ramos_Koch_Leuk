pdf(file='invitro.poolC_all_v_input.pdf',width=4.5,height=4.5);
gstable=read.table('invitro.poolC_all_v_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Mto1","Sepsecs","Tmem208","Ctc1","Rabl3","Commd3","Manf","Dcps","Tubd1","Mrm2")
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
targetmat=list(c(433.62004571128284,186.30844600770746,120.78212204405128,244.08212410365573,119.43118183549575,142.5969108742914,200.10634552196106,176.93988020199347,181.63973875906032,152.66539138997211,268.91131878185746,320.9576549921902,259.1217970897258,181.741657369732),c(417.99409811808346,111.34410086851156,99.73675229395143,86.29166003664596,90.56864622525094,121.377132470379,94.63157226654317,96.29381235482639,117.84959240915224,131.57346231635756,148.9914063521102,152.58642614382813,116.12816324422495,119.91202135734895),c(369.1630118893354,166.4649428826262,121.69713812014257,87.52439803716948,83.60182728484702,124.77229701500498,92.66008117765685,107.12686624474435,170.82784954721149,135.59097261609367,129.00475428048568,70.50545208025163,179.3920133698102,99.3021426865546),c(161.14258455486862,80.47642934060737,35.6856269675606,38.21487801622893,62.701370463635264,66.20570862020672,45.34429504438527,39.72119759636588,37.8416122414709,44.192613297097196,64.50237714024284,44.19744757269505,34.66512335648506,74.0082006814888))
targetgene="Mto1"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(332.05138635548684,128.98277031302825,136.33739533760334,215.7291500916149,148.29371744574055,158.72394246126484,136.03288513315582,143.23704587780423,188.12687228616963,148.647881090236,92.66538687753196,97.86577676811046,60.66396587384886,63.70326134609163),c(236.34245734714065,91.50059774343029,114.3770095114122,140.53213205968058,115.45014244097922,97.61098065799709,128.14692077761055,120.36726544353297,71.35846879820227,42.18385814722914,55.41753528950441,54.72064937571768,45.93128844734271,54.335134677548744),c(266.61773080896444,140.00693871585116,135.42237926151205,161.4886780685803,141.32689850533663,148.53844882738687,254.32235046633477,167.31049896651083,111.36245888204293,89.38960416912842,117.19445987452572,98.91809694841272,75.39664330035501,65.57688667980021),c(778.3675144862442,444.27398663376397,431.8875879150924,428.9928241821828,492.6536250714199,529.6456689616538,478.0865890549316,421.28542905236543,416.2577346561799,435.89986752136775,393.37365213697393,295.7019706649359,388.2493815926327,456.22776875803856))
targetgene="Sepsecs"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(923.8841514479134,696.727443058409,640.5112532639083,666.911258283221,618.0563659986905,606.0368712157385,606.2335098325422,586.1885827100056,474.64193640016356,556.425176513451,466.96087112795516,487.2242434799478,298.1200608657715,363.483314739464),c(298.84624771993816,153.23594079923868,214.11376180536362,151.6267740643922,184.12307199638926,162.96789814204732,172.50547027755266,209.43904187174738,132.9862373057406,244.0637507089686,152.6253430924056,91.55185568629689,188.05829420893147,89.9340160180117),c(337.9111167029366,103.62718298653552,119.86710596795997,150.39403606386867,119.43118183549575,106.94768315571855,122.23244751095159,114.34890217135633,65.95252419227785,124.54281929181937,87.2144817670889,57.877609916624465,38.13163569213357,32.788443339900105),c(424.8304501901082,238.1220375009752,187.57829559871598,145.46308406177462,182.13255229913102,203.70987267755916,284.88046234407267,222.679441070536,147.04169328114406,223.97619921028806,144.448985426741,207.30707551954583,82.32966797165201,323.2003700647296))
targetgene="Tmem208"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(439.47977605873257,280.0138774317023,239.73421193591994,255.17676610836736,225.92398563881278,272.46195470623536,290.7949356107316,296.1034729910911,209.75065070986727,269.17319008231925,187.1477421252116,176.7897902907802,338.85158080964146,292.2855520585381),c(444.3628846816074,214.97128385504706,215.9437939575462,321.7446181366371,209.00456821211756,220.6856954006891,234.6074395774716,178.1435528564288,381.65968917826365,294.28262945566996,256.19254019082365,869.2164689296693,839.7626133108506,860.9308408390913),c(507.8432967789799,270.0921258691617,338.5559481537801,378.4505661607187,296.5874348914811,237.661518123819,290.7949356107316,312.9548901531857,322.1942985130951,265.15567978258315,273.45373970722665,228.3534791255911,180.2586414537223,125.53289735847468),c(336.93449497836167,249.14620590379815,194.89842420744637,214.4964120910914,216.9666470011506,179.09492972902075,184.33441681087055,223.88311372497134,168.66547170484174,223.97619921028806,183.51380538491622,152.58642614382813,116.99479132813708,179.8680320360234))
targetgene="Ctc1"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(428.736937088408,313.0863826401711,264.439645990385,411.7344921748536,349.33620686882506,236.81272698766253,314.4528286773674,293.6961276822205,324.35667635546486,241.0506179841665,260.7349611161929,250.4522029119386,275.5877306840562,227.64547804559214),c(370.13963361391035,267.88729218859714,224.1789386423679,225.59105409580303,170.18943411558143,250.39338516616647,225.7357296774832,231.1051496515833,190.28925012853938,150.65663624010406,214.4022676774269,107.33665839083083,103.12874198554306,184.55209537029486),c(283.2203001267388,102.52476614625323,134.50736318542073,159.02320206753328,121.421701532754,116.28438565344001,187.29165344420002,77.0350498838611,101.63175859137898,182.79671863799294,136.27262776107642,320.9576549921902,231.3896984045378,169.56309270062624),c(947.3230728377125,572.15434010651,636.851188959543,626.230904265945,576.2554523562669,591.607421901078,637.7773672547232,670.4456685204786,495.18452590267634,630.749117058569,628.6710560710992,410.40487031788257,779.9652755209139,859.0572155053827))
targetgene="Rabl3"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(461.9420757239567,285.5259616331138,259.8645656099285,295.8571201256433,291.61113564833545,331.8773342371901,241.5076583885737,239.53085823263064,150.28526004469873,148.647881090236,188.9647104953593,146.27250506201457,128.26095641899474,176.12078136860626),c(355.4903077452859,248.04378906351585,234.24411547937217,205.8672460874268,209.00456821211756,219.8369042645326,215.8782742330516,174.53253489312283,188.12687228616963,91.39835931899647,150.8083747222579,225.1965185846843,140.3937495937645,174.2471560348977),c(436.5499108850077,195.12778072996582,194.89842420744637,226.82379209632654,171.18469396421057,151.93361337201287,157.71928711090527,243.14187619593662,82.1703580100511,71.31080782031593,129.91323846555952,71.55777226055389,96.19571731424604,72.13457534778023),c(267.5943525335394,156.54319132008555,199.47350458790285,212.03093609004438,156.25579623477358,128.16746155963096,207.00656433306318,148.05173649554555,219.48135100053122,101.44213506833674,166.25260586851323,82.08097406357652,97.06234539815817,102.11258068711746))
targetgene="Commd3"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(242.2021876945904,184.10361232714288,172.02302230516395,192.30712808166814,164.21787502380664,155.32877791663887,134.0613940442695,149.2554091499809,130.82385946337084,142.62161564063186,119.91991242974724,70.50545208025163,123.92781599943409,133.9642113601633),c(459.0122105502318,358.28547309174513,306.53038549058465,357.494020151819,329.4310098962424,293.6817331101478,337.12497619956,332.213652624151,232.4556180547498,246.07250585883665,272.54525552215284,227.30115894528882,174.19224486633743,231.3927287130093),c(378.929229135085,306.47188159847735,228.7540190228244,223.125578094756,204.0282689689719,272.46195470623536,190.2488900775295,226.290459033842,179.47736091669054,173.7573204635867,168.97805842373475,137.85394361959646,105.72862623727944,48.71425867642301),c(800.8298141514683,497.189994967314,491.36363286102676,602.8088822559982,617.0611061500614,569.5388523610092,517.5164108326579,522.3939320249332,382.74087809944854,575.5083504371976,435.1639246503707,474.5964013163206,301.58657320142004,290.41192672482947))
targetgene="Manf"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(139.65690661421948,42.99425677100942,75.9463343155777,40.68035401727595,86.58760683073442,53.473841577859275,33.51534851106737,30.091816360883243,73.52084664057203,59.2582769211076,138.0895961312241,71.55777226055389,13.866049342594025,144.26915069556046),c(315.4488170377125,186.30844600770746,180.25816698998563,125.73927605339841,184.12307199638926,119.67955019806601,151.80481384424633,193.79129736408808,136.22980406929523,170.7441877387846,93.5738710626058,169.42354902866435,74.53001521644288,87.12357801744885),c(368.18639016476044,253.55587326492733,239.73421193591994,236.68569610051463,230.90028488195844,216.4417397199066,269.10853363298213,239.53085823263064,257.3229632420021,289.2607415809998,157.1677640177748,173.6328297498734,168.99247636286466,319.45311939731243),c(206.06718388531684,98.11509878512405,104.31183267440792,113.41189604816327,67.67766970678092,162.11910700589084,82.80262573322527,128.79297402458027,110.28126996085805,117.51217626728118,107.20113383871345,11.575521983324894,34.66512335648506,132.0905860264547))
targetgene="Dcps"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(231.45934872426582,73.86192829891361,91.50160760912975,98.6190400418811,72.65396894992658,85.72790475180615,88.71709899988421,77.0350498838611,105.95651427611853,64.28016479577774,67.22782969546437,46.302087933299575,56.330825454288224,79.62907668261454),c(175.79191042349305,85.98851354201884,93.33163976131235,219.42736409318545,86.58760683073442,136.65537292119595,110.4035009776337,84.25708581047309,72.43965771938716,102.44651264327076,88.12296595216274,47.354408113601835,48.531172699079086,88.06039068430313),c(279.31381322843896,155.44077447980328,150.97765255506408,77.66249403298137,162.22735532654838,202.86108154140265,139.97586731092844,169.7178442753815,238.9427515818591,175.76607561345475,223.48710952816532,352.5272604012581,800.7643495348049,301.65367872708094),c(487.3342405629057,208.35678281335333,200.38852066399417,202.16903208585626,212.98560760663406,195.22196131599418,185.3201623553137,199.80966063626474,190.28925012853938,298.30013975540606,255.28405600574982,152.58642614382813,262.58830942537435,192.98340937198347))
targetgene="Tubd1"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(208.99704905904173,124.57310295189907,158.29778116379447,139.29939405915707,150.28423714299882,157.02636018895186,156.73354156646212,172.12518958425215,211.91302855223705,251.0943937335068,146.26595379688868,287.28340922251783,335.385068473993,206.0987867079435),c(341.81760360123644,264.58004166775027,170.19299015298134,191.07439008114463,170.18943411558143,240.2078915322885,170.53397918866634,231.1051496515833,202.1823282615731,219.95868891055193,169.8865426088086,122.06914091506252,328.452043802696,162.06859136579195),c(774.4610275879444,339.54438680694614,416.3323146215404,385.8469941638598,386.1608212681029,293.6817331101478,423.87058411055796,311.7512174987504,326.51905419783463,248.0812610087047,307.06765455495884,444.079116087555,243.52249157930757,389.71406941138406),c(548.8614092111283,227.0978690981523,183.0032152182595,208.33272208847384,170.18943411558143,223.23206880915856,179.40568908865475,172.12518958425215,273.54079705977534,187.81860651266308,288.897970853482,169.42354902866435,173.3256167824253,188.299346037712))
targetgene="Mrm2"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetgenelist=c("Ddx6","Kctd5","Ppp6r1","Hexim1","Ccny","Rlim","Mark3","Ints13","Blmh","Ankrd28")
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
targetmat=list(c(406.2746374231839,496.0875781270317,655.1515104813691,634.8600702696097,632.9852637281274,827.571357752584,941.3869949432159,1026.7327742333364,2907.31700906615,2908.6774570089424,2847.1894360214233,10518.992522301422,12571.306985229308,11161.186112901994),c(578.1600609483771,636.0945168428829,845.4748543083589,957.8374264067702,742.4638470773318,852.1863007011225,1121.7784295763138,895.5324548998854,3499.808537875466,2310.068422348262,3380.469652659769,10411.65586391059,7633.260163098011,11396.326092282421),c(216.81002285564142,179.6939449660137,333.9808677733236,366.1231861554836,262.7486000380906,359.8874417303545,277.9802435329706,448.969900104378,1097.4067550026562,1250.4500807928637,1122.8864527512696,3675.754389795805,6703.368229060298,5308.9173830632535))
targetgene="Ddx6"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(512.7264054018547,693.4201925375621,796.0639861994289,759.3666083224845,736.492287985557,715.5309277799265,861.5416058433201,737.8513371688572,683.3113981888459,703.064302453819,718.6109903934097,779.769253603977,1593.7290463144006,1758.3973756854998),c(526.3991095459041,817.9932954894612,749.3981663187727,808.676128343425,870.8523675504898,840.3032247949316,713.6797741768464,806.4606784716709,474.64193640016356,599.6134122356142,540.5480901189364,910.2569559614575,1838.1181659776205,829.0792101660454),c(614.2950647576507,899.5721416703508,1038.5432463636228,1141.5153884847737,944.5015963490455,929.4262940913636,887.1709899988422,1180.8028740010586,818.4600133369563,727.1693642522356,900.3078274081781,857.6409469463443,1973.312147067912,1651.600731664111),c(861.3803610751158,1192.8150211854406,1201.4161079078738,1231.5052625229903,1307.7714410986785,1251.1181346946757,1186.8376355095622,1396.2602791449826,907.1175048741167,950.1411858875897,952.9999101424611,716.6300427858412,471.4456776481968,403.76625941419843))
targetgene="Kctd5"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(496.1238360840804,611.8413463566725,645.0863336443648,565.8267422402928,616.0658463014322,654.4179659766588,559.9034692437137,609.0583631442769,1020.6423415985295,945.1192980129196,730.4212847993696,1599.5266740594398,1313.808175210784,2052.5565530777462),c(380.8824725842349,490.57549392562026,511.49398653503533,504.18984221411716,605.1179879665118,502.4843526046459,568.7751791437022,656.0015966672547,610.8717404694588,758.3050690751904,682.271622990456,1156.4998781521872,1113.6170878270825,1378.0514329426587),c(119.14785039814528,217.17611753561167,136.33739533760334,150.39403606386867,172.1799538128397,156.17756905279535,278.9659890774137,210.64271452618271,220.5625399217161,198.86675983693738,275.27070807737437,446.18375644815956,285.98726769100176,347.5574994029411),c(341.81760360123644,362.6951404528743,366.006430436519,382.14878016228926,468.7673887043208,409.1173276274313,533.2883395437485,430.91481028784807,581.679639597467,474.0662153688608,655.9255816233144,690.3220382782846,668.1702526962496,591.1287927850561))
targetgene="Ppp6r1"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(574.2535740500773,751.8482850725236,748.4831502426814,601.5761442554747,824.0751546649207,702.7990607375791,793.5251632767422,826.9231135970715,1215.2563474118083,1176.1261402477458,955.7253626976826,1394.3242389004986,1724.589886985132,1997.2846057333434),c(404.321393974034,489.47307708533793,579.2051761657914,708.8243503010204,436.9190735481886,594.1537953095475,615.1052197325306,645.1685427773367,860.6263812631668,946.1236755878535,832.17151352764,1004.9657721886612,1006.155205421979,1003.3263662009432),c(360.37341636816075,498.2924118075963,399.862025251897,468.4404401989352,452.8432311262547,350.550739232633,319.3815563995832,397.2119759636588,478.9666920849031,407.77729542321504,428.8045353548538,797.6586966691153,448.0467193825694,687.6204974710479),c(515.6562705755796,684.6008578153038,631.3610925029953,493.0952002094055,696.6818940403919,693.4623582398576,683.1216622991085,806.4606784716709,982.8007293570586,931.0580119638431,1270.9693749183061,1752.113100203268,2110.239384326028,1617.8754756573564))
targetgene="Hexim1"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(336.93449497836167,346.1588878486399,472.1482952631095,517.7499602198758,422.9854356673808,393.8390871766144,433.7280395549895,482.67273442856725,520.0518710899287,533.3244922899684,497.8493334204658,498.7997654632727,530.3763873542215,468.4063334271443),c(475.61477986800617,659.245270488811,559.9898385678741,671.842210285315,549.3834364432804,679.0329089251971,534.2740850881917,636.7428341962894,776.2936454107459,819.5721011461661,863.0599758201506,958.6636842553617,500.91103250120915,842.1945875020054),c(253.92164838948995,270.0921258691617,253.45945307728942,208.33272208847384,239.85762351962063,281.79865720395685,303.60962768849265,204.62435125400606,305.97646469532185,276.2038331068575,361.5767056593894,357.7888613027694,173.3256167824253,271.6756733877437),c(181.6516407709428,263.47762482746793,249.79938877292423,304.4862861293079,328.4357500476133,404.02458081049235,265.1655514552095,259.9932933580312,225.96848452764053,350.52777365197545,352.491863808651,370.4167034663966,298.1200608657715,181.741657369732))
targetgene="Ccny"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(1237.379725036476,1387.9428019154066,1346.90366400639,1503.9403606386868,1577.4868600771729,1564.322063936423,1455.9461691425442,1397.4639517994178,2670.536635326661,2776.099617117651,2882.620319239303,4512.348933136104,5925.136209707209,5012.884580337299),c(914.1179342021638,1058.320166671001,1180.3707381577738,1082.3439644596451,931.5632183168668,1244.327805605424,1068.5481701763833,1083.3053889917967,1750.4448633983254,1778.7526852081621,1769.7271925238458,2958.0720268296614,2939.6024606299334,3449.344239357491),c(111.33487660154559,87.09093038230112,74.1163021633951,134.368442057063,120.42644168412487,140.05053746582192,147.8618316664737,93.88646704595573,163.25952709891732,111.48591081767701,119.91991242974724,115.75521983324893,415.11485219390863,211.71966270906924),c(976.6217245749614,1133.2845118101968,1120.8946932118395,1019.4743264329459,1040.0465418174422,1163.6926476705567,1378.0722711315348,1250.6158879583077,2083.4510511232693,2054.9565183150194,2071.3439419683614,6103.457045753125,5491.822167751146,7138.5125214296795))
targetgene="Rlim"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(712.9338589397217,779.408706079581,850.0499346888154,1065.0856324523158,947.4873758949329,1132.2873756327663,903.9286642543759,956.9197602760872,1095.2443771602864,975.2506252609403,1032.9465184289593,995.4948905659409,835.42947289129,885.2879701773028),c(784.2272448336939,935.9518973996664,876.5854008954631,840.7273163570364,970.3783524134029,950.646072495276,1079.391371165258,964.1417962026992,1034.6977975739328,1065.6446070050029,1221.0027447392447,1153.3429176112802,1191.613615379174,1019.2521815374661),c(327.16827773261207,450.8884876754577,438.29270044773153,430.2255621827063,545.4023970487639,488.0549032899855,514.5591741993285,501.93149689953253,558.9746722525845,444.939265695774,546.9074794144533,938.6696008296186,1204.6130366378559,1201.9306515740523),c(730.513049982071,733.1071987877247,711.8825071990295,637.3255462706567,601.1369485719953,765.6096048131599,819.1545474322643,718.5925746978919,781.6995900166703,778.392620573871,786.7473042739479,1044.9539390401471,642.1714101788857,1008.0104295352146))
targetgene="Mark3"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(681.681963753323,750.7458682322414,690.8371374489296,576.9213842450044,718.5776107102328,736.7507061838389,816.1973107989348,850.9965666857781,2152.6471420791017,1991.6807310941758,2108.591793556389,4280.838493469606,3892.8933529332726,4053.588409478507),c(670.9391247829984,787.125623961557,845.4748543083589,811.1416043444721,836.0182728484702,880.1964081942868,834.9264761433548,980.9932133647937,2098.5876960198575,1723.5119185867904,1528.9788834792773,3019.1065972871925,2826.074181637445,2790.7649345589257),c(278.337191503864,396.8700625016254,329.4057873928671,347.6321161476309,457.81953036940035,412.5124921720573,478.0865890549316,463.413971957602,807.6481241251075,689.0030164047427,696.8073699516374,1278.5690190672497,1405.6707521054693,1739.661122348414),c(641.6404730457496,770.5893713573226,684.4320249162905,830.8654123528482,658.8620197924848,827.571357752584,704.808064276858,781.183552728529,1479.0664441809197,1275.5595201662145,1409.05897104953,2284.587111436213,1593.7290463144006,3132.7015579607414))
targetgene="Ints13"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(680.7053420287481,809.1739607672029,817.1093559495287,754.4356563203904,795.2126190546758,787.6781743532288,875.3420434655243,869.0516565023081,897.3868045834528,955.1630737622598,901.216311593252,1510.0794587337475,1498.3999570840667,1334.9580502673614),c(725.6299413591962,772.7942050378872,689.007105296747,938.113618398394,848.9566508806489,824.1761932079581,786.62494446564,834.1451495236836,1151.4662010619002,1044.5526779313882,1241.897880995943,1575.3233099124877,1268.7435148473533,1411.776688949413),c(142.58677178794434,281.11629427198466,261.6945977621111,177.51427207538597,153.27001668888622,182.49009427364675,210.9495465108358,252.77125743141926,248.67345187252306,246.07250585883665,212.5852993072792,271.4986065179839,295.52017661403517,435.61789008724423),c(451.19923675363214,551.2084201411463,483.1284881762051,465.9749641978882,567.2981137186048,461.7423780691341,578.6326345881338,505.5425148628385,579.5172617550973,567.4733298377254,618.6777300352869,1084.9421058916332,987.0893875759122,1018.3153688706118))
targetgene="Blmh"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
targetmat=list(c(139.65690661421948,154.33835763952098,178.42813483780301,204.63450808690328,173.17521366146883,157.87515132510833,85.75986236655474,182.95824347417013,211.91302855223705,190.83173923746514,129.91323846555952,396.7247079739532,456.71300022169066,550.8458481103218),c(82.03622486429676,51.81359149326776,128.10225065278166,138.06665605863355,142.32215835396576,70.44966430098921,89.70284454432738,99.90483031813237,91.90105830071504,102.44651264327076,86.30599758201507,77.87169334236746,171.59236061460106,109.60708202195177),c(1276.4445940194744,1484.9554838602483,1390.8244356587722,1466.9582206229813,1694.9275222154104,1512.5458046308768,1482.5612988425096,1665.8829537384963,1986.1440482166297,1897.2692390503773,2006.8415648281186,3089.6120493674443,2573.018781135104,2991.2428452657437),c(520.5393791984544,582.0760916690506,627.7010281986301,679.2386382884561,565.3075940213465,689.2184025590751,776.7674890212085,752.2954090220811,864.9511369479063,928.044879239041,757.6758103515849,1003.913452008359,1471.534486482791,1219.730092244284))
targetgene="Ankrd28"
collabel=c("In.vitro.input_Pool.C..17.24._1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._no.CAR_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._hEGFRv3.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.10_rep3_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep1_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep2_1392..1393..1394..1395..1396..1397..1398..1399","In.vitro_Pool.C..17.24._mCD19.CAR_ET.1.2_rep3_1392..1393..1394..1395..1396..1397..1398..1399")

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
Sweave("invitro.poolC_all_v_input_summary.Rnw");
library(tools);

texi2dvi("invitro.poolC_all_v_input_summary.tex",pdf=TRUE);

