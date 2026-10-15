args<-commandArgs(trailingOnly=TRUE);out<-normalizePath(args[1]);tab<-file.path(out,'tables');fig<-file.path(out,'figures')
pal<-hcl.colors(256,'viridis');colprob<-function(p)pal[pmax(1,pmin(256,1+floor((log10(p)+5)/5*255)))]
symbols<-function(x,y,status,p) {
 if(status=='BH q < 0.05')points(x,y,pch=21,bg=colprob(p),col='black',lwd=1.2,cex=1.7)
 else if(status=='Finite nominal only')points(x,y,pch=16,col=colprob(p),cex=1.5)
 else if(status=='Tested not nominal')points(x,y,pch=1,col='#999999',cex=1.35)
 else if(status=='Not selected by model')points(x,y,pch=4,col='#999999',cex=1.1)
 else points(x,y,pch=15,col='#dddddd',cex=1.5)
}
panel<-function(d,xs,ys,ylabels,title,xcol='subtype',paircol='interaction_name',rotate=FALSE) {
 plot(NA,xlim=c(.5,length(xs)+.5),ylim=c(length(ys)+.6,.4),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
 abline(h=seq_along(ys),col='#eeeeee',lwd=.6)
 axis(1,at=seq_along(xs),labels=xs,las=if(rotate)2 else 1,tick=FALSE,cex.axis=.83)
 axis(2,at=seq_along(ys),labels=ylabels,las=1,tick=FALSE,cex.axis=.86)
 title(main=title,adj=0,cex.main=1.05,line=1)
 for(i in seq_len(nrow(d)))symbols(match(d[[xcol]][i],xs),match(d[[paircol]][i],ys),d$status[i],d$probability[i])
}
key<-function() {
 par(mar=c(.1,.2,.1,.2));plot.new();plot.window(c(0,1),c(0,1))
 legend(.01,.98,legend=c('BH q < 0.05','Finite P < 0.05 only','Tested, not nominal','Not selected','Ineligible'),pch=c(21,16,1,4,15),pt.bg=c('#21918c',NA,NA,NA,NA),col=c('black','#21918c','#999999','#999999','#dddddd'),ncol=3,bty='n',cex=.8,pt.cex=1.15,xjust=0,yjust=1)
 x0<-.38;x1<-.88;y0<-.22;y1<-.32
 for(i in 1:256)rect(x0+(x1-x0)*(i-1)/256,y0,x0+(x1-x0)*i/256,y1,col=pal[i],border=NA)
 text(.01,.27,'log10 modeled probability',adj=0,cex=.8)
 text(seq(x0,x1,length.out=6),.1,-5:0,cex=.75)
}
draw5<-function() {
 a<-read.csv(file.path(tab,'S17_Figure5_subtype_data.csv'));b<-read.csv(file.path(tab,'S17_Figure5_donor_data.csv'))
 b$donor<-sub('^donor__ATRT-(.*?)__.*$','\\1',b$run)
 layout(matrix(1:4,4,1),heights=c(2.6,2.4,3.5,1.4));par(oma=c(2,0,.4,0),family='sans')
 par(mar=c(2.3,17,2.2,1.2))
 panel(a[a$panel=='A',],c('MYC','SHH','TYR'),c('SPP1_CD44','SPP1_ITGAV_ITGB1','TGFB1_TGFBR1_TGFBR2','GAS6_AXL','GAS6_MERTK'),c('SPP1–CD44','SPP1–(ITGAV + ITGB1)','TGFB1–(TGFBR1 + TGFBR2)','GAS6–AXL','GAS6–MERTK'),'A  Source myeloid → tumor')
 panel(a[a$panel=='B',],c('MYC','SHH','TYR'),c('APP_CD74','MIF_CD74_CXCR4','MIF_CD74_CD44','TGFB1_TGFBR1_TGFBR2'),c('APP–CD74','MIF–(CD74 + CXCR4)','MIF–(CD74 + CD44)','TGFB1–(TGFBR1 + TGFBR2)'),'B  Tumor → source myeloid')
 par(mar=c(4.5,17,2.2,1.2))
 panel(b,sort(unique(b$donor)),c('APP_CD74','SPP1_CD44','SPP1_ITGAV_ITGB1','NCAM1_FGFR1','PTN_PTPRZ1'),c('Tumor → myeloid: APP–CD74','Myeloid → tumor: SPP1–CD44','Myeloid → tumor:\nSPP1–(ITGAV + ITGB1)','Astrocyte → tumor: NCAM1–FGFR1','Astrocyte → tumor: PTN–PTPRZ1'),'C  Same-donor candidate support',xcol='donor',rotate=TRUE)
 key();mtext('All models: 10,000 permutations. Model support does not establish active signaling or subtype differences.',side=1,outer=TRUE,line=.35,cex=.69)
}
draw9<-function() {
 d<-read.csv(file.path(tab,'S18_Figure_S9_source_data.csv'));d$lib<-sub('SCPCL','',d$library_id)
 layout(matrix(1:4,4,1),heights=c(2.5,2.5,3.2,1.4));par(oma=c(2,0,.4,0),family='sans')
 pairs<-c('APP_CD74','SPP1_ITGAV_ITGB1','SPP1_CD44');labs<-c('A  Tumor → myeloid: APP–CD74','B  Myeloid → tumor: SPP1–(ITGAV + ITGB1)','C  Myeloid → tumor: SPP1–CD44')
 for(k in 1:3) {
  z<-d[d$interaction_name==pairs[k],];z$pop<-if(k==1)z$target else z$source
  par(mar=c(if(k==3)4.8 else 2.3,12,2.2,1.2))
  panel(z,sort(unique(z$lib)),c('Microglia','Macrophage_like','BAM_like','Myeloid_unresolved'),c('Microglia','Macrophage-like','BAM-like','Unresolved myeloid'),labs[k],xcol='lib',paircol='pop',rotate=(k==3))
  abline(v=c(2.5,4.5,6.5),col='#dddddd',lty=2)
 }
 key();mtext('Primary conservative models: ≥20 nuclei per group; 10,000 permutations. Adjacent specimens share one patient.',side=1,outer=TRUE,line=.35,cex=.7)
}
for(name in c('Figure_5_fresh_models','Figure_S9_external_CellChat'))for(ext in c('png','pdf','svg')) {
 p<-file.path(fig,paste0(name,'.',ext));w<-9;h<-if(grepl('Figure_5',name))8.5 else 9.5
 if(ext=='png')png(p,width=w,height=h,units='in',res=230,type='cairo')
 else if(ext=='pdf')cairo_pdf(p,width=w,height=h) else svg(p,width=w,height=h)
 if(grepl('Figure_5',name))draw5() else draw9()
 dev.off()
}
