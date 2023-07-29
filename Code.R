library(ggplot2)
library(ggnewscale)
library(ggpattern)
library(reshape2)

library(stringi)

library(foreach)
library(doParallel)
library(dplyr)
library(doRNG)

library(diptest)
library(evt0)
library(Epi)

library(shades)
library(ggpubr)
library(tidyr)
library(multiplex)
library(gridExtra)

#parallel setup
cores=detectCores()
mcl <- makeCluster(cores[1]-8) #set number of cores to something suitable to your machine 
registerDoParallel(mcl)

library(RColorBrewer)
library(scales)


#### WEEKLY DATA #### -------------------------------------------------------

data2 <- read.csv("C:/Users/yasse/Downloads/new CL plots/weekly.csv")
data2 <- data2[,-1]

datelength<-max(data2$terminus)
dates<-seq(min(base::as.Date(data2$date)), max(base::as.Date(data2$date)), by='weeks')



#full graph
library(leidenAlg)
library(igraph)

#create list of graphs
try<-foreach(i=1:(datelength-1), .packages = c("igraph"),
             .combine = c) %dopar% 
  {
    s<-data2[data2$onset==i,][,c(2,3,4)]
    g<-graph_from_data_frame(s[,c(1,2)], directed = F)
    E(g)$weight<-s$count
    return(list(g))
  }


clusterExport(mcl, list('strength'))
deg <- parSapply(mcl, try, strength)
deg1<-deg

clusterExport(mcl, list('E'))

tweights<-parSapply(mcl, try, function(x){E(x)$weight})


#largest connected subgraph
clusterExport(mcl, list('V','components','induced_subgraph'))

#create largest connected subgraph from each graph in the list
tryc<-parLapply(mcl,try, function(x){comp<-igraph::components(x);
biggest_cluster_id <- which.max(comp$csize);
vert_ids <- V(x)[comp$membership == biggest_cluster_id];
x<-induced_subgraph(x, vert_ids); return(x)})

vl<-sapply(tryc,vcount)
vd<-matrix(vl,nrow=length(vl),ncol=length(vl), byrow = T)
vd[upper.tri(vd, diag = T)]<-NA;vd[1,1]<-vl[1]
minvp<-apply(vd, 1, function(x){min(x,na.rm = T)})


el<-sapply(tryc,ecount)
ed<-matrix(el,nrow=length(el),ncol=length(el), byrow = T)
ed[upper.tri(ed, diag = T)]<-NA;ed[1,1]<-el[1]
minep<-apply(ed, 1, function(x){min(x,na.rm = T)})



degc <- parSapply(mcl, tryc, strength)
degce<-vector('list', length = length(degc))
for(i in 1:length(degc)){
  degce[[i]]<-list(degc[[i]],minvp[i])
}


trycd<-vector('list', length(tryc))
for (i in 1:length(tryc)) {
  trycd[[i]]<-list(tryc[[i]],degc[[i]])
}


tweightsc<-parSapply(mcl, tryc, function(x){E(x)$weight})
tweightsce<-vector('list', length = length(tweightsc))
for(i in 1:length(tweightsc)){
  tweightsce[[i]]<-list(tweightsc[[i]],minep[i])
}





#### Figure 1 ####  --------------------------------

#loss from using largest connected subgraph
#number of vertices
loss<-foreach(i=1:length(tryc), .packages = c('igraph'), .combine = rbind) %dopar% {
  return(c(length(V(try[[i]])),length(V(tryc[[i]]))))
}

#number of edges
losse<-foreach(i=1:length(tryc), .packages = c('igraph'), .combine = rbind)%dopar% {
  return(c(length(E(try[[i]])),length(E(tryc[[i]]))))
}

#degree comparison - no loss
lossd<-foreach(i=1:length(tryc), .packages = c('igraph'), .combine = rbind)%dopar% {
  return(c(sum(strength(try[[i]])),sum(strength(tryc[[i]]))))
}

#weight comparison - no loss
lossw<-foreach(i=1:length(tryc), .packages = c('igraph'), .combine = rbind)%dopar% {
  return(c(sum(E(try[[i]])$weight),sum(E(tryc[[i]])$weight)))
}



par(mfrow=(c(1,3)))
plot(dates,loss[,1], ylim=c(min(loss,2), max(loss[,1])), type='l', ylab = '',
     xlab = 'Date', main = '(a) Number of Nodes')
lines(dates,loss[,2], lty=2)
legend('bottomright', legend = c('Full Graph', "Largest Connected"), lty=c(1,2))

plot(dates,losse[,1], ylim=c(min(losse,2), max(losse[,1])), type='l', ylab = '',
     xlab = 'Date', main = '(b) Number of Edges')
lines(dates,losse[,2], lty=2)
legend('bottomright', legend = c('Full Graph', "Largest Connected"), lty=c(1,2))

plot(dates,lossw[,1], ylim=c(min(lossw,2), max(lossw[,1])), type='l', ylab = '',
     xlab = 'Date', main = '(c) Sum of Weights')
lines(dates,lossw[,2], lty=2)
legend('bottomright', legend = c('Full Graph', "Largest Connected"), lty=c(1,2))
par(mfrow=c(1,1))






#clustering-------------------
library(tidyr)
library(multiplex)

lmaccuracy<-function(testgr, d, itr=5, gma=1){
  
  
  
  ec<-data.frame("key"=V(testgr)$name)
  #cl<-leiden(testgr, partition_type=type, resolution_parameter = gma, n_iterations = itr )
  cl<-find_partition(testgr,E(testgr)$weight, resolution = gma, niter = itr)
  mmbr<-data.frame('indicator'=as.character(ec$key),"cluster"=cl)
  
  ec<-ec %>%
    separate(key, c("company", "kpi"), "-")
  
  clstr<-data.frame('indicator'=as.character(ec[,1]),"cluster"=cl+1)
  
  acc1<-c()
  
  
  for(i in sort(unique(clstr$cluster))) {
    subst<-NA
    subst<-clstr[clstr$cluster==i,]
    tbl<-as.matrix(sort(table(subst$indicator), decreasing = T))
    order<-cbind(row.names(tbl),as.numeric(tbl))
    
    result<-c(i,nrow(subst),apply(order, 1, c))
    
    outlier<-NA
    outlier<-c(order[-1,1])
    out<-c()
    out<-list(result, outlier)
    acc1<-c(acc1,out)
  }
  
  acc2<-c()
  
  for(i in unique(ec$company)) {
    subst<-NA
    subst<-ec[ec$company==i,]
    out<-c(i,0)
    try(out<-c(i,nrow(subst)))
    acc2<-rbind(acc2,out)
  }
  
  output<-list(list(clustering=acc1,real=acc2,
                    ncluster=length(unique(as.numeric(as.character(cl)))),
                    nreal=nrow(acc2), graph=testgr, clusters=mmbr, optg=gma))
  return(output)
}

clusterExport(mcl, list('lmaccuracy'))
registerDoRNG(seed = 2385186, once = TRUE)
partitionhuge<-foreach(i=1:length(trycd), .packages = c("igraph", "leidenAlg", "tidyr","multiplex"),
                       .combine = c) %dopar% {
                         x<-trycd[[i]]
                         out<-lmaccuracy(testgr=x[[1]],d=x[[2]],gma=1,itr=20)
                         return(out)
                       }

#cluster the clusters----------------------------
clhuge<-lapply(partitionhuge, function(x){n<-length(x$clustering);n1<-seq(1,n-1,2);
nmax<-max(sapply(x$clustering, length));
acc<-t(data.frame(sapply(x$clustering[n1], function(x){out<-rep(NA,nmax);
out[1:length(x)]<-x;return(out)})));
acc<-as.data.frame(acc);
colnames(acc)<-c('cl.nbr', 'cl.size', 'prefix', 'count',
                 paste(c('prefix','count'),c(apply(data.frame(v1=2:(ncol(acc)/2-1)), 1,
                                                   function(x){return(rep(x,2))})),sep=''));
acc$cl.size<-as.numeric(as.character(acc$cl.size));
acc$count<-as.numeric(as.character(acc$count));
acc$accuracy<-acc$count/acc$cl.size;
try(acc$count2<-as.numeric(as.character(acc$count2)));
try(acc$accuracy2<-acc$count2/acc$cl.size)
return(acc)})


clhuge<-lapply(clhuge, function(x){
  x$name<-NA;
  for (i in 1:nrow(x)) {
    x$name[i]<-paste(x$prefix[i], round(x$accuracy[i],3)*100, "%",
                     if(is.na(x$prefix2[i])){NULL} else {
                       paste(",", x$prefix2[i], round(x$accuracy2[i],3)*100, "%")}, "Size:", x$cl.size[i])
    x$name1[i]<-paste(x$prefix[i], 
                      if(is.na(x$prefix2[i])){NULL} else {
                        paste(",", x$prefix2[i])})
    
  };
  return(x)}
)


#25% cluster dominance threshold
thresh<-0.25

clcl<-vector("list", length(partitionhuge))

for (i in 1:length(partitionhuge)) {
  #print(i)
  tst<-NA
  tstm<-NA
  tst<-partitionhuge[[i]]$graph
  tstm<-as.numeric(partitionhuge[[i]]$clusters$cluster+1)
  
  testcol<-NA
  #testing<-contract(tst, (tstm))
  testcol<-contract(tst, (tstm),  vertex.attr.comb=toString)
  V(testcol)$clname<-V(testcol)$name
  V(testcol)$name<-make.unique(c(clhuge[[i]][order(clhuge[[i]]$cl.size, decreasing = T),]$prefix))
  
  V(testcol)$name[which(clhuge[[i]][order(clhuge[[i]]$cl.size, decreasing = T),]$accuracy<thresh)]<-
    'no_dominance'
  
  #V(testcol)$name[which(clhuge[[i]][order(clhuge[[i]]$cl.size, decreasing = T),]$accuracy<thresh)]<-
  #  clhuge[[i]][order(clhuge[[i]]$cl.size, decreasing = T),]$name1[
  #    which(clhuge[[i]][order(clhuge[[i]]$cl.size, decreasing = T),]$accuracy<thresh)]
  
  V(testcol)$name<-make.unique(V(testcol)$name)
  
  V(testcol)$namelong<-make.unique(c(clhuge[[i]][order(clhuge[[i]]$cl.size, decreasing = T),]$name))
  V(testcol)$size<-(c(clhuge[[i]][order(clhuge[[i]]$cl.size, decreasing = T),]$cl.size))
  
  x<-partitionhuge[[i]]$clusters
  
  sumw<-c()
  for (j in sort(unique(x$cluster))) {
    gr<-induced_subgraph(tst, V(tst)[which(x$cluster==j)])
    sumw<-c(sumw,sum(E(gr)$weight))
    
  }
  
  V(testcol)$sumw<-sumw
  
  testcol<-igraph::simplify(testcol, edge.attr.comb=list(weight="sum"), remove.loops = F)
  d<-strength(testcol)
  testcol<-igraph::simplify(testcol, edge.attr.comb=list(weight="sum"), remove.loops = T)
  d1<-strength(testcol)
  b<-igraph::betweenness(testcol, weights = NA, normalized = T, directed = F)
  
  V(testcol)$d<-d
  V(testcol)$d1<-d1
  V(testcol)$b<-b
  clcl[[i]]<-list(testcol,d,d1)
}




#### NODE Analysis #### ----------------------------------------------------------


#### Figure 3 #### ------------------------------------------------------------------
len<-10

cldom<-foreach(i=1:length(partitionhuge), .packages = c('igraph', 'tidyr'), .combine=c) %dopar%{
  x<-partitionhuge[[i]]$clusters
  y<-partitionhuge[[i]]$graph
  clusterdom<-c()
  for (j in sort(unique(x$cluster))) {
    gr<-induced_subgraph(y, V(y)[which(x$cluster==j)])
    btw<-igraph::betweenness(gr, weights = NA)
    eigen<-eigen_centrality(gr, weights = NULL)
    str<-igraph::strength(gr)
    clusterdom<-rbind(clusterdom,c(j,
                                   V(gr)$name[which(btw==max(btw))[1]],
                                   V(gr)$name[which(eigen$value==max(eigen$value))[1]],
                                   V(gr)$name[which(str==max(str))[1]],
                                   length(str)))
  }
  gr<-y
  btw<-igraph::betweenness(gr, weights = NA)
  eigen<-eigen_centrality(gr, weights = NULL)
  str<-igraph::strength(gr)
  clusterdom<-rbind(cbind(rep('full',len), V(gr)$name[order(btw, decreasing=T)][1:len],
                          V(gr)$name[order(eigen$value, decreasing=T)][1:len],
                          V(gr)$name[order(str, decreasing=T)][1:len],
                          rep(length(str),len)),
                    clusterdom)
  return(list(clusterdom))
}



a<-sort(table(sapply(cldom[1:112], function(x){x[1:5,2]})), decreasing = T)[1:10]
a1<-sort(table(sapply(cldom[113:length(cldom)], function(x){x[1:5,2]})), decreasing = T)[1:10]

b<-sort(table(sapply(cldom[1:112], function(x){x[1:5,4]})), decreasing = T)[1:10]
b1<-sort(table(sapply(cldom[113:length(cldom)], function(x){x[1:5,4]})), decreasing = T)[1:10]


c<-sort(table(sapply(cldom[1:112], function(x){x[11:15,2]})), decreasing = T)[1:10]
c1<-sort(table(sapply(cldom[113:length(cldom)], function(x){x[11:15,2]})), decreasing = T)[1:10]

sort(table(sapply(cldom[1:112], function(x){x[11:15,4]})), decreasing = T)[1:10]
sort(table(sapply(cldom[113:length(cldom)], function(x){x[11:15,4]})), decreasing = T)[1:10]

a<-a/length(1:112)
b<-b/length(1:112)
c<-c/length(1:112)

a1<-a1/(length(cldom)-length(1:112))
b1<-b1/(length(cldom)-length(1:112))
c1<-c1/(length(cldom)-length(1:112))

a<-data.frame(a)
b<-data.frame(b)
c<-data.frame(c)

a1<-data.frame(a1)
b1<-data.frame(b1)
c1<-data.frame(c1)

f<-cbind(b,a,c)
f1<-cbind(b1,a1,c1)


#nested pie plot

#25

col.c<-c("#1e2b00",
  "#f148f9",
  "#5aff96",
  "#5f61ff",
  "#ffbe01",
  "#001e6f",
  "#fdff8d",
  "#015aba",
  "#c02100",
  "#06ffec",
  "#e60040",
  "#00a25b",
  "#ae005c",
  "#b6fff2",
  "#600053",
  "#2f5a00",
  "#ff74be",
  "#008b97",
  "#ffac5c",
  "#9ba9ff",
  "#773500",
  "#ffaae6",
  "#0c0006",
  "#005675",
  "#2e001e")




potential<-na.omit(unique(c(apply(f[,seq(1,ncol(f)-1,2)],2,as.character), apply(f1[,seq(1,ncol(f)-1,2)],2,as.character))))
potential<-sort(potential)

dfcol<-data.frame("browser"=potential,'cl.hue'=col.c[1:length(potential)])
dfcol$eh<-1

gcol<-ggplot(data=dfcol, aes(eh, fill=browser))+
  geom_bar()+
  scale_fill_manual(values=brightness(dfcol$cl.hue, 0.8)) +
  guides(fill=guide_legend(nrow=3)) + 
  theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
        legend.title=element_text(size=11), 
        legend.text=element_text(size=10))+
  labs(fill='Nodes')
#gcol  

f.legend<-get_legend(gcol)

plots<-list()



for (i in seq(1,ncol(f)-1,2)) {
  col.ch<-c(i,i+1)
  
  l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
         rep("(c) Cluster Dominance",2))
  
  #l.n<-c(rep("(a) Cluster Size",2),rep("(b) Sum of Weights",2))
  
  #l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
  #       rep("(c) Top Prefixes",2))
  
  #l.n<-c(rep("(a) Strength",2), rep("(b) Betweenness Centrality",2))
  
  # l.n<-c(rep("(a) Edge Weights",2),rep("(b) Edge Betweenness",2),
  #        rep("(c) Cluster Dominance",2))
  
  thr<-0.1
  
  
  eh<-merge(na.omit(f[,c(col.ch)]),na.omit(f1[,c(col.ch)]), by='Var1', all=T)
  eh[is.na(eh)]<-0
  df<-rbind(as.matrix(eh[,c(1,2)]),as.matrix(eh[,c(1,3)]))
  df<-data.frame(df)
  df$period<-c(rep(2019,nrow(df)/2),rep(2020,nrow(df)/2))
  colnames(df)<-c("browser","cc","year")
  #df<-na.omit(df)
  df$cc<-as.numeric(df$cc)
  
  # Data for plot
  pdat = df %>% 
    group_by(year) %>% 
    arrange(browser) %>% 
    mutate(ccc = cc/sum(cc)) %>%
    mutate(cc_cum = cumsum(ccc)-0.5*ccc) %>% 
    ungroup
  
  pdat1<-merge(pdat,dfcol)
  if(nrow(pdat1)==0){pdat$cl.hue<-dfcol$cl.hue[1:nrow(pdat)]}
  else(pdat<-pdat1)

  pdat$browser<-gsub(", ",",\n ", pdat$browser)
  
  
  d<-pdat[pdat$year==2020,]
  m<-which(d$cc>0.1)
  m<-unlist(lapply(m,function(x){rep(x,2)}))
  
  pdat$m<-NA
  pdat$m[1:length(m)]<-m
  
  pdat$ang<- -pdat$cc_cum*360+90
  #pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]<-pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]-180
  pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]<-pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]-180
  
  
  cl.hue<-unique(pdat$cl.hue)
  
  
  g<-NA
  
  g<-ggplot(data=pdat, aes(x=cc_cum, y=year, fill=browser)) +
    geom_tile(aes(width=ccc), colour="white", size=0.4) +
    geom_text(aes(label=replace(round(100*cc,1),which(cc<thr),""))
              , size=2.5, colour="white", position=position_nudge(y=-0.2)) +
    geom_text(data=pdat %>% filter(year==unique(year)[2]), size=3.5, 
              aes(label=replace(browser,which(cc<thr),""),angle=ang, colour=browser), 
              position=position_nudge(y=0.5)) +
    geom_text(data=pdat %>% filter(year==unique(year)[1]), size=3.5, 
              aes(label=replace(browser,c(m,which(cc<thr)),""),angle=(ang), colour=browser),
              position=position_nudge(y=0.5)) +
    scale_y_continuous(#breaks=min(pdat$year):max(pdat$year),
      #labels=c("Before","After")) +
      breaks=NULL) +
    coord_polar() +
    theme_void() +
    theme(axis.text.y=element_text(angle=0, colour="grey40", size=9),
          axis.ticks.y=element_line(),
          axis.ticks.length=unit(0.1,"cm")) +
    #guides(fill=FALSE, colour=FALSE) +
    guides(colour=FALSE) +
    #scale_fill_manual(values=brightness(col.c, 0.8)) +
    scale_fill_manual(values=brightness(cl.hue, 0.8)) +
    #scale_colour_manual(values=brightness(col.c, 1.5)) +
    scale_colour_manual(values=rep("black",(nrow(pdat)/2))) +
    #scale_fill_manual(values=sample(hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,70),
    #                  size=(nrow(pdat)/2))) +
    #scale_colour_manual(values=hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,40)) +
    theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
          legend.title=element_text(size=11), 
          legend.text=element_text(size=10)) +
    ggtitle(l.n[i])
  
  #print(g)
  #print(g1)
  #assign(paste("g",i,sep = ""),g)
  
  plots<-c(plots,list(g))
  
}


ggarrange(plotlist=plots, ncol=ncol(f)/2, nrow=1, legend = "bottom", 
          legend.grob = f.legend)








#### Figure 4 #### ------------------------------------------------
potential<-unique(c(apply(f[,c(1,3,5)],2,as.character), apply(f1[,c(1,3,5)],2,as.character)))
potential

#topn <- potential[-5]
topn<-potential[c(1:6,8,16,17,20,19,21)]

clusterExport(mcl, list('topn'))

topo<-parSapply(mcl, trycd, function(x){
  d<-c()
  for (i in topn) {
    try(di<-x[[2]][which(V(x[[1]])$name==i)])
    if(length(di)==0){di<-0}
    d<-c(d,di)
  }
  return(c(d,max(x[[2]])))
})

rownames(topo)<-c(topn,'max')

topodf<-melt(topo)
colnames(topodf)[1]<-'name'
topodf<-topodf[order(topodf$name),]
topodf$dates<-rep(dates,length(topn)+1)

ggplot(topodf, aes(x=dates,y=value, col=name, linetype=name)) +
  geom_line(lwd = 1.5) +
  scale_linetype_manual(values = c(rep("solid", length(topn)), rep("dashed", 1))) +
  scale_color_manual(values = c(hue_pal()(length(topn)), 'black')) +
  xlab("Date") + ylab("Weighted Degree") + labs(col="Node", linetype="Node") +
  theme(legend.position="bottom", plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=11), 
        legend.text=element_text(size=7)) +
  guides(col=guide_legend(nrow=2), linetype=guide_legend(nrow=2))



#### Figure 5 #### ------------------------------------------

topo<-parSapply(mcl, trycd, function(x){
  d<-c()
  bet<-igraph::betweenness(x[[1]], weights=NA, directed=F, normalized = T)
  for (i in topn) {
    try(di<-bet[which(V(x[[1]])$name==i)])
    if(length(di)==0){di<-0}
    d<-c(d,di)
  }
  return(c(d,max(bet)))
})

rownames(topo)<-c(topn,'max')

topodf<-melt(topo)
colnames(topodf)[1]<-'name'
topodf<-topodf[order(topodf$name),]
topodf$dates<-rep(dates,length(topn)+1)


ggplot(topodf, aes(x=dates,y=value, col=name, linetype=name)) +
  geom_line(lwd = 1.5) +
  scale_linetype_manual(values = c(rep("solid", length(topn)), rep("dashed", 1))) +
  scale_color_manual(values = c(hue_pal()(length(topn)), 'black')) +
  xlab("Date") + ylab("Normalized Betweenness Centrality") + labs(col="Node", linetype="Node")+
  theme(legend.position="bottom", plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=11), 
        legend.text=element_text(size=7)) +
  guides(col=guide_legend(nrow=2), linetype=guide_legend(nrow=2))



#### EDGE Analysis #### ----------------------
edges<-parLapply(mcl, tryc, function(x){y<-data.frame(igraph::as_edgelist(x),
                                                      'weight'=as.numeric(E(x)$weight),
                                                      'bet'=as.numeric(edge_betweenness(x, directed = F, 
                                                                                        weights = NA)));
y$e<-apply(y[,c(1:2)],1,function(z){toString(sort(as.character(z), decreasing = F))});
y1<-data.frame(y$e[order(y$weight, decreasing = T)][1:10],y$e[order(y$bet, decreasing = T)][1:10]);
return(y1)})

len<-10

ecldom<-foreach(i=1:length(partitionhuge), .packages = c('igraph', 'tidyr'), .combine=c) %dopar%{
  x<-partitionhuge[[i]]$clusters
  y<-partitionhuge[[i]]$graph
  clusterdom<-c()
  for (j in sort(unique(x$cluster))) {
    gr<-induced_subgraph(y, V(y)[which(x$cluster==j)])
    y1<-data.frame(as_edgelist(gr))
    y2<-apply(y1[,c(1:2)],1,function(z){toString(sort(as.character(z), decreasing = F))})
    btw<-igraph::edge_betweenness(gr, weights = NA, directed = F)
    eigen<-strength(gr, weights = NULL)
    str<-E(gr)$weight
    clusterdom<-rbind(clusterdom,c(j,
                                   y2[which(btw==max(btw))[1]],
                                   y2[which(eigen==max(eigen))[1]],
                                   y2[which(str==max(str))[1]],
                                   length(str)))
  }
  gr<-y
  y1<-data.frame(as_edgelist(gr))
  y2<-apply(y1[,c(1:2)],1,function(z){toString(sort(as.character(z), decreasing = F))})
  btw<-igraph::edge_betweenness(gr, weights = NA, directed=F)
  eigen<-strength(gr, weights = NULL)
  str<-E(gr)$weight
  clusterdom<-rbind(cbind(rep('full',len), y2[order(btw, decreasing=T)][1:len],
                          y2[order(eigen, decreasing=T)][1:len],
                          y2[order(str, decreasing=T)][1:len],
                          rep(length(str),len)),
                    clusterdom)
  return(list(clusterdom))
}



a<-sort(table(sapply(ecldom[1:112], function(x){x[1:5,2]})), decreasing = T)[1:10]
a1<-sort(table(sapply(ecldom[113:length(ecldom)], function(x){x[1:5,2]})), decreasing = T)[1:10]

b<-sort(table(sapply(ecldom[1:112], function(x){x[1:5,4]})), decreasing = T)[1:10]
b1<-sort(table(sapply(ecldom[113:length(ecldom)], function(x){x[1:5,4]})), decreasing = T)[1:10]


c<-sort(table(sapply(ecldom[1:112], function(x){x[11:15,2]})), decreasing = T)[1:10]
c1<-sort(table(sapply(ecldom[113:length(ecldom)], function(x){x[11:15,2]})), decreasing = T)[1:10]

sort(table(sapply(ecldom[1:112], function(x){x[11:15,4]})), decreasing = T)[1:10]
sort(table(sapply(ecldom[113:length(ecldom)], function(x){x[11:15,4]})), decreasing = T)[1:10]

a<-a/length(1:112)
b<-b/length(1:112)
c<-c/length(1:112)

a1<-a1/(length(ecldom)-length(1:112))
b1<-b1/(length(ecldom)-length(1:112))
c1<-c1/(length(ecldom)-length(1:112))

a<-data.frame(a)
b<-data.frame(b)
c<-data.frame(c)

a1<-data.frame(a1)
b1<-data.frame(b1)
c1<-data.frame(c1)

f<-cbind(b,a,c)
f1<-cbind(b1,a1,c1)

ra<-c(rep("(a) Edge Weights",2),rep("(b) Edge Betweenness",2),rep("(c) Cluster Dominance",2))


f2<-c()
for (i in seq(1,ncol(f)-1,2)) {
  t<-na.omit(f[,c(i,i+1)][1:6,])
  f2<-c(f2, list(list(t,ra[i])))
}


f3<-c()
for (i in seq(1,ncol(f)-1,2)) {
  t<-na.omit(f1[,c(i,i+1)][1:6,])
  f3<-c(f3, list(list(t,ra[i])))
}


#### Figure 6 #### ----------------------------------
par(mfrow=c(1,3))
sapply(f2,function(x){y<-separate(x[[1]]["Var1"],Var1, c("l1", "l2"), ",");
barplot(x[[1]][,2], col=1:nrow(x[[1]]), legend.text = paste(y$l1,"\n", y$l2, sep = ""),
        beside = F, xlab = paste(x[[2]]), cex.lab=2, ylim = c(0,1),
        args.legend = list(x='topright', cex=1, inset=c(-0.01,0.1),
                           y.intersp=1.5))})
par(mfrow=c(1,1))


#### Figure 7 #### -------------------------------------------------
par(mfrow=c(1,3))
sapply(f3,function(x){y<-separate(x[[1]]["Var1"],Var1, c("l1", "l2"), ",");
barplot(x[[1]][,2], col=1:nrow(x[[1]]), legend.text = paste(y$l1,"\n", y$l2, sep = ""),
        beside = F, xlab = paste(x[[2]]), cex.lab=2, ylim = c(0,1),
        args.legend = list(x='topright', cex=1, inset=c(-0.03,0.1),
                           y.intersp=1.5))})
par(mfrow=c(1,1))


#### Cluster Analysis #### ------------------------------------------

#### Figure 8 #### --------------------------------------------------

plot(dates,sapply(partitionhuge, function(x){x$nclust}), type='l',
     xlab='Dates', ylab= "Number of Clusters")




#### Figure 9 #### -------------------------------------------------

cld<-sapply(clcl, function(x){
  x<-x[[1]]
  y<-V(x)$name
  return(c(y[order(V(x)$size, decreasing=T)][1:10],
           y[order(V(x)$sumw, decreasing=T)][1:10]))
})

a<-sort(table(as.character(cld[1:5,1:112])), decreasing = T)[1:10]
a1<-sort(table(as.character(cld[1:5,113:ncol(cld)])), decreasing = T)[1:10]

b<-sort(table(as.character(cld[11:15,1:112])), decreasing = T)[1:10]
b1<-sort(table(as.character(cld[11:15,113:ncol(cld)])), decreasing = T)[1:10]

a<-a/length(1:112)
b<-b/length(1:112)

a1<-a1/(ncol(cld)-length(1:112))
b1<-b1/(ncol(cld)-length(1:112))

a<-data.frame(a)
b<-data.frame(b)

a1<-data.frame(a1)
b1<-data.frame(b1)

f<-cbind(a,b)
f1<-cbind(a1,b1)

col.c<-c("#c97252",
         "#ac47d8",
         "#5cbc40",
         "#736cd5",
         "#9fae3f",
         "#c953aa",
         "#58bb7e",
         "#d44365",
         "#65b9ad",
         "#d94c28",
         "#78a8d5",
         "#d79b39",
         "#6d699d",
         "#8b6c2a",
         "#cf99bc",
         "#4f7e45",
         "#9e5c6b",
         "#4f7b82",
         "#bcac84",
         "#787053")


potential<-na.omit(unique(c(apply(f[,seq(1,ncol(f)-1,2)],2,as.character), apply(f1[,seq(1,ncol(f)-1,2)],2,as.character))))
potential<-sort(potential)

dfcol<-data.frame("browser"=potential,'cl.hue'=col.c[1:length(potential)])
dfcol$eh<-1

gcol<-ggplot(data=dfcol, aes(eh, fill=browser))+
  geom_bar()+
  scale_fill_manual(values=brightness(dfcol$cl.hue, 0.8)) +
  guides(fill=guide_legend(nrow=2)) + 
  theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
        legend.title=element_text(size=11), 
        legend.text=element_text(size=10))+
  labs(fill='Clusters')
gcol  

f.legend<-get_legend(gcol)

plots<-list()

for (i in seq(1,ncol(f)-1,2)) {
  col.ch<-c(i,i+1)
  
  # l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
  #        rep("(c) Cluster Dominance",2))
  
  l.n<-c(rep("(a) Cluster Size",2),rep("(b) Sum of Weights",2))
  
  #l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
  #       rep("(c) Top Prefixes",2))
  
  #l.n<-c(rep("(a) Strength",2), rep("(b) Betweenness Centrality",2))
  
  # l.n<-c(rep("(a) Edge Weights",2),rep("(b) Edge Betweenness",2),
  #        rep("(c) Cluster Dominance",2))
  
  thr<-0.1
  
  
  eh<-merge(na.omit(f[,c(col.ch)]),na.omit(f1[,c(col.ch)]), by='Var1', all=T)
  eh[is.na(eh)]<-0
  df<-rbind(as.matrix(eh[,c(1,2)]),as.matrix(eh[,c(1,3)]))
  df<-data.frame(df)
  df$period<-c(rep(2019,nrow(df)/2),rep(2020,nrow(df)/2))
  colnames(df)<-c("browser","cc","year")
  #df<-na.omit(df)
  df$cc<-as.numeric(df$cc)
  
  # Data for plot
  pdat = df %>% 
    group_by(year) %>% 
    arrange(browser) %>% 
    mutate(ccc = cc/sum(cc)) %>%
    # Get cumulative value of cc
    mutate(cc_cum = cumsum(ccc)-0.5*ccc) %>% 
    ungroup
  
  pdat1<-merge(pdat,dfcol)
  if(nrow(pdat1)==0){pdat$cl.hue<-dfcol$cl.hue[1:nrow(pdat)]}
  else(pdat<-pdat1)
  
  #854adb
  
  
  
  
  #fix the colors
  
  #only for links
  pdat$browser<-gsub(", ",",\n ", pdat$browser)
  
  
  d<-pdat[pdat$year==2020,]
  m<-which(d$cc>0.1)
  m<-unlist(lapply(m,function(x){rep(x,2)}))
  
  pdat$m<-NA
  pdat$m[1:length(m)]<-m
  
  pdat$ang<- -pdat$cc_cum*360+90
  #pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]<-pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]-180
  pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]<-pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]-180
  
  
  cl.hue<-unique(pdat$cl.hue)
  
  
  g<-NA
  
  g<-ggplot(data=pdat, aes(x=cc_cum, y=year, fill=browser)) +
    geom_tile(aes(width=ccc), colour="white", size=0.4) +
    geom_text(aes(label=replace(round(100*cc,1),which(cc<thr),""))
              , size=2.5, colour="white", position=position_nudge(y=-0.2)) +
    geom_text(data=pdat %>% filter(year==unique(year)[2]), size=3.5, 
              aes(label=replace(browser,which(cc<thr),""),angle=ang, colour=browser), 
              position=position_nudge(y=0.4)) +
    geom_text(data=pdat %>% filter(year==unique(year)[1]), size=3.5, 
              aes(label=replace(browser,c(m,which(cc<thr)),""),angle=(ang), colour=browser),
              position=position_nudge(y=0.4)) +
    scale_y_continuous(#breaks=min(pdat$year):max(pdat$year),
      #labels=c("Before","After")) +
      breaks=NULL) +
    coord_polar() +
    theme_void() +
    theme(axis.text.y=element_text(angle=0, colour="grey40", size=9),
          axis.ticks.y=element_line(),
          axis.ticks.length=unit(0.1,"cm")) +
    #guides(fill=FALSE, colour=FALSE) +
    guides(colour=FALSE) +
    #scale_fill_manual(values=brightness(col.c, 0.8)) +
    scale_fill_manual(values=brightness(cl.hue, 0.8)) +
    #scale_colour_manual(values=brightness(col.c, 1.5)) +
    scale_colour_manual(values=rep("black",(nrow(pdat)/2))) +
    #scale_fill_manual(values=sample(hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,70),
    #                  size=(nrow(pdat)/2))) +
    #scale_colour_manual(values=hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,40)) +
    theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
          legend.title=element_text(size=11), 
          legend.text=element_text(size=10)) +
    ggtitle(l.n[i])
  
  #print(g)
  #print(g1)
  #assign(paste("g",i,sep = ""),g)
  
  plots<-c(plots,list(g))
  
}




ggarrange(plotlist=plots, ncol=ncol(f)/2, nrow=1, legend = "bottom", 
          legend.grob = f.legend)




#### Figure 10 #### ----------------------------------------------------


nod<-vector("list", length(partitionhuge))

for(i in 1:length(clcl)){
  testcol<-clcl[[i]][[1]]
  testgr<-partitionhuge[[i]]$graph
  m<-which(V(testcol)$name=='no_dominance'|V(testcol)$name=='no_dominance.1')
  dom<-c()
  if(length(m)!=0){
    for (j in m) {
      gr<-induced_subgraph(testgr, V(testgr)[which(partitionhuge[[i]]$clusters$cluster==(j-1))])
      d<-igraph::strength(gr)
      b<-igraph::betweenness(gr, directed = F, weights=NA, normalized = T)
      dom<-rbind(dom, cbind(V(gr)$name[which(d==max(d))],
                            V(gr)$name[which(b==max(b))]))
    }
  } else {
    dom<-cbind(NA,NA)
  }
  nod[[i]]<-dom
}

a<-sort(table((na.omit(unlist(sapply(nod[1:112], function(x){return(as.character(x[,1]))}))))),
        decreasing = T)[1:10]
a1<-sort(table((na.omit(unlist(sapply(nod[113:length(nod)], function(x){return(as.character(x[,1]))}))))),
         decreasing = T)[1:10]

b<-sort(table((na.omit(unlist(sapply(nod[1:112], function(x){return(as.character(x[,2]))}))))),
        decreasing = T)[1:10]
b1<-sort(table((na.omit(unlist(sapply(nod[113:length(nod)], function(x){return(as.character(x[,2]))}))))),
         decreasing = T)[1:10]

l.nd1<-length(which(is.na(unlist(sapply(nod[1:112], function(x){return(as.character(x[,1]))})))))
l.nd2<-length(which(is.na(unlist(sapply(nod[113:length(nod)], function(x){return(as.character(x[,1]))})))))

a<-a/(length(1:112)-l.nd1)
b<-b/(length(1:112)-l.nd1)

a1<-a1/(ncol(cld)-length(1:112)-l.nd2)
b1<-b1/(ncol(cld)-length(1:112)-l.nd2)

a<-data.frame(a)
b<-data.frame(b)

a1<-data.frame(a1)
b1<-data.frame(b1)

f<-cbind(a,b)
f1<-cbind(a1,b1)


nod<-vector("list", length(partitionhuge))

for(i in 1:length(clcl)){
  testcol<-clcl[[i]][[1]]
  testgr<-partitionhuge[[i]]$graph
  #m<-which(V(testcol)$name=='no_dominance')
  m<-which(V(testcol)$name=='no_dominance'|V(testcol)$name=='no_dominance.1')
  dom<-c()
  if(length(m)!=0){
    for (j in m) {
      dom<-c(dom, as.character(clhuge[[i]][j,c(3)]))
      #dom<-c(dom, as.character(clhuge[[i]][j,c(3,5)]))
    }
  } else {
    dom<-c(NA,NA)
  }
  nod[[i]]<-dom
}

c<-sort(table(na.omit(unlist(nod[1:112]))),
        decreasing = T)[1:10]
c1<-sort(table(na.omit(unlist(nod[113:length(nod)]))),
         decreasing = T)[1:10]



c<-c/(length(1:112)-l.nd1)

c1<-c1/(ncol(cld)-length(1:112)-l.nd2)

c<-data.frame(c)
c1<-data.frame(c1)

f<-cbind(a,b,c)
f1<-cbind(a1,b1,c1)




potential<-na.omit(unique(c(apply(f[,seq(1,ncol(f)-3,2)],2,as.character), apply(f1[,seq(1,ncol(f)-3,2)],2,as.character))))
potential<-sort(potential)

dfcol<-data.frame("browser"=potential,'cl.hue'=col.c[1:length(potential)])
dfcol$eh<-1

potential2<-na.omit(unique(c(f[,5], f1[,5])))
potential2<-sort(potential2)

dfcol2<-data.frame("browser"=potential2,'cl.hue'=col.c[1:length(potential2)] )
dfcol2$eh<-1

gcol<-ggplot()+
  geom_bar(data=dfcol, aes(eh, fill=browser))+
  scale_fill_manual(values=brightness(dfcol$cl.hue, 0.8)) +
  guides(fill=guide_legend(nrow=4)) + 
  labs(fill='Nodes') +
  new_scale_fill() +
  geom_bar(data=dfcol2, aes(eh, fill=browser))+
  scale_fill_manual(limits = potential2,
                    values = brightness(dfcol2$cl.hue, 0.8))+
  guides(fill=guide_legend(nrow=4)) +
  theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10))+
  labs(fill='Prefixes')
gcol  

f.legend<-get_legend(gcol)

plots<-list()





for (i in seq(1,ncol(f)-1,2)) {
  col.ch<-c(i,i+1)
  
  #l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
  #       rep("(c) Cluster Dominance",2))
  
  #l.n<-c(rep("(a) Cluster Size",2),rep("(b) Sum of Weights",2))
  
  l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
         rep("(c) Top Prefixes",2))
  
  #l.n<-c(rep("(a) Strength",2), rep("(b) Betweenness Centrality",2))
  
  #l.n<-c(rep("(a) Edge Weights",2),rep("(b) Edge Betweenness",2),
  #       rep("(c) Cluster Dominance",2))
  
  thr<-0.1
  
  
  eh<-merge(na.omit(f[,c(col.ch)]),na.omit(f1[,c(col.ch)]), by='Var1', all=T)
  eh[is.na(eh)]<-0
  df<-rbind(as.matrix(eh[,c(1,2)]),as.matrix(eh[,c(1,3)]))
  df<-data.frame(df)
  df$period<-c(rep(2019,nrow(df)/2),rep(2020,nrow(df)/2))
  colnames(df)<-c("browser","cc","year")
  #df<-na.omit(df)
  df$cc<-as.numeric(df$cc)
  
  # Data for plot
  pdat = df %>% 
    group_by(year) %>% 
    arrange(browser) %>% 
    mutate(ccc = cc/sum(cc)) %>%
    # Get cumulative value of cc
    mutate(cc_cum = cumsum(ccc)-0.5*ccc) %>% 
    ungroup
  
  pdat1<-merge(pdat,dfcol)
  if(nrow(pdat1)==0){pdat$cl.hue<-dfcol2$cl.hue[1:nrow(pdat)]}
  else(pdat<-pdat1)
  
  #854adb
  
  
  
  
  #fix the colors
  
  #only for links
  #pdat$browser<-gsub(", ",",\n ", pdat$browser)
  
  
  d<-pdat[pdat$year==2020,]
  m<-which(d$cc>0.1)
  m<-unlist(lapply(m,function(x){rep(x,2)}))
  
  pdat$m<-NA
  pdat$m[1:length(m)]<-m
  
  pdat$ang<- -pdat$cc_cum*360+90
  #pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]<-pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]-180
  pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]<-pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]-180
  
  
  cl.hue<-unique(pdat$cl.hue)
  
  
  g<-NA
  
  g<-ggplot(data=pdat, aes(x=cc_cum, y=year, fill=browser)) +
    geom_tile(aes(width=ccc), colour="white", size=0.4) +
    geom_text(aes(label=replace(round(100*cc,1),which(cc<thr),""))
              , size=2.5, colour="white", position=position_nudge(y=-0.2)) +
    geom_text(data=pdat %>% filter(year==unique(year)[2]), size=3.5, 
              aes(label=replace(browser,which(cc<thr),""),angle=ang, colour=browser), 
              position=position_nudge(y=0.4)) +
    geom_text(data=pdat %>% filter(year==unique(year)[1]), size=3.5, 
              aes(label=replace(browser,c(m,which(cc<thr)),""),angle=(ang), colour=browser),
              position=position_nudge(y=0.4)) +
    scale_y_continuous(#breaks=min(pdat$year):max(pdat$year),
      #labels=c("Before","After")) +
      breaks=NULL) +
    coord_polar() +
    theme_void() +
    theme(axis.text.y=element_text(angle=0, colour="grey40", size=9),
          axis.ticks.y=element_line(),
          axis.ticks.length=unit(0.1,"cm")) +
    #guides(fill=FALSE, colour=FALSE) +
    guides(colour=FALSE) +
    #scale_fill_manual(values=brightness(col.c, 0.8)) +
    scale_fill_manual(values=brightness(cl.hue, 0.8)) +
    #scale_colour_manual(values=brightness(col.c, 1.5)) +
    scale_colour_manual(values=rep("black",(nrow(pdat)/2))) +
    #scale_fill_manual(values=sample(hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,70),
    #                  size=(nrow(pdat)/2))) +
    #scale_colour_manual(values=hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,40)) +
    theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
          legend.title=element_text(size=11), 
          legend.text=element_text(size=10)) +
    ggtitle(l.n[i])
  
  #print(g)
  #print(g1)
  #assign(paste("g",i,sep = ""),g)
  
  plots<-c(plots,list(g))
  
}



ggarrange(plotlist=plots, ncol=ncol(f)/2, nrow=1, legend = "bottom", 
          legend.grob = f.legend)








#### Figure 11 #### ------------------------------------------------------

cld<-sapply(clcl, function(x){
  x<-x[[1]]
  y<-V(x)$name
  return(c(y[order(V(x)$d1, decreasing=T)][1:10],
           y[order(V(x)$b, decreasing=T)][1:10]))
})

a<-sort(table(as.character(cld[1:5,1:112])), decreasing = T)[1:10]
a1<-sort(table(as.character(cld[1:5,113:ncol(cld)])), decreasing = T)[1:10]

b<-sort(table(as.character(cld[11:15,1:112])), decreasing = T)[1:10]
b1<-sort(table(as.character(cld[11:15,113:ncol(cld)])), decreasing = T)[1:10]

a<-a/length(1:112)
b<-b/length(1:112)

a1<-a1/(ncol(cld)-length(1:112))
b1<-b1/(ncol(cld)-length(1:112))

a<-data.frame(a)
b<-data.frame(b)

a1<-data.frame(a1)
b1<-data.frame(b1)

f<-cbind(a,b)
f1<-cbind(a1,b1)


potential<-na.omit(unique(c(apply(f[,seq(1,ncol(f)-1,2)],2,as.character), apply(f1[,seq(1,ncol(f)-1,2)],2,as.character))))
potential<-sort(potential)

dfcol<-data.frame("browser"=potential,'cl.hue'=col.c[1:length(potential)])
dfcol$eh<-1

gcol<-ggplot(data=dfcol, aes(eh, fill=browser))+
  geom_bar()+
  scale_fill_manual(values=brightness(dfcol$cl.hue, 0.8)) +
  guides(fill=guide_legend(nrow=2)) + 
  theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
        legend.title=element_text(size=11), 
        legend.text=element_text(size=10))+
  labs(fill='Clusters')
gcol  

f.legend<-get_legend(gcol)

plots<-list()



for (i in seq(1,ncol(f)-1,2)) {
  col.ch<-c(i,i+1)
  
  # l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
  #        rep("(c) Cluster Dominance",2))
  
  #l.n<-c(rep("(a) Cluster Size",2),rep("(b) Sum of Weights",2))
  
  #l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
  #       rep("(c) Top Prefixes",2))
  
  l.n<-c(rep("(a) Strength",2), rep("(b) Betweenness Centrality",2))
  
  # l.n<-c(rep("(a) Edge Weights",2),rep("(b) Edge Betweenness",2),
  #        rep("(c) Cluster Dominance",2))
  
  thr<-0.1
  
  
  eh<-merge(na.omit(f[,c(col.ch)]),na.omit(f1[,c(col.ch)]), by='Var1', all=T)
  eh[is.na(eh)]<-0
  df<-rbind(as.matrix(eh[,c(1,2)]),as.matrix(eh[,c(1,3)]))
  df<-data.frame(df)
  df$period<-c(rep(2019,nrow(df)/2),rep(2020,nrow(df)/2))
  colnames(df)<-c("browser","cc","year")
  #df<-na.omit(df)
  df$cc<-as.numeric(df$cc)
  
  # Data for plot
  pdat = df %>% 
    group_by(year) %>% 
    arrange(browser) %>% 
    mutate(ccc = cc/sum(cc)) %>%
    # Get cumulative value of cc
    mutate(cc_cum = cumsum(ccc)-0.5*ccc) %>% 
    ungroup
  
  pdat1<-merge(pdat,dfcol)
  if(nrow(pdat1)==0){pdat$cl.hue<-dfcol$cl.hue[1:nrow(pdat)]}
  else(pdat<-pdat1)
  
  #854adb
  
  
  
  
  #fix the colors
  
  #only for links
  pdat$browser<-gsub(", ",",\n ", pdat$browser)
  
  
  d<-pdat[pdat$year==2020,]
  m<-which(d$cc>0.1)
  m<-unlist(lapply(m,function(x){rep(x,2)}))
  
  pdat$m<-NA
  pdat$m[1:length(m)]<-m
  
  pdat$ang<- -pdat$cc_cum*360+90
  #pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]<-pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]-180
  pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]<-pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]-180
  
  
  cl.hue<-unique(pdat$cl.hue)
  
  
  g<-NA
  
  g<-ggplot(data=pdat, aes(x=cc_cum, y=year, fill=browser)) +
    geom_tile(aes(width=ccc), colour="white", size=0.4) +
    geom_text(aes(label=replace(round(100*cc,1),which(cc<thr),""))
              , size=2.5, colour="white", position=position_nudge(y=-0.2)) +
    geom_text(data=pdat %>% filter(year==unique(year)[2]), size=3.5, 
              aes(label=replace(browser,which(cc<thr),""),angle=ang, colour=browser), 
              position=position_nudge(y=0.4)) +
    geom_text(data=pdat %>% filter(year==unique(year)[1]), size=3.5, 
              aes(label=replace(browser,c(m,which(cc<thr)),""),angle=(ang), colour=browser),
              position=position_nudge(y=0.4)) +
    scale_y_continuous(#breaks=min(pdat$year):max(pdat$year),
      #labels=c("Before","After")) +
      breaks=NULL) +
    coord_polar() +
    theme_void() +
    theme(axis.text.y=element_text(angle=0, colour="grey40", size=9),
          axis.ticks.y=element_line(),
          axis.ticks.length=unit(0.1,"cm")) +
    #guides(fill=FALSE, colour=FALSE) +
    guides(colour=FALSE) +
    #scale_fill_manual(values=brightness(col.c, 0.8)) +
    scale_fill_manual(values=brightness(cl.hue, 0.8)) +
    #scale_colour_manual(values=brightness(col.c, 1.5)) +
    scale_colour_manual(values=rep("black",(nrow(pdat)/2))) +
    #scale_fill_manual(values=sample(hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,70),
    #                  size=(nrow(pdat)/2))) +
    #scale_colour_manual(values=hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,40)) +
    theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
          legend.title=element_text(size=11), 
          legend.text=element_text(size=10)) +
    ggtitle(l.n[i])
  
  #print(g)
  #print(g1)
  #assign(paste("g",i,sep = ""),g)
  
  plots<-c(plots,list(g))
  
}


ggarrange(plotlist=plots, ncol=ncol(f)/2, nrow=1, legend = "bottom", 
          legend.grob = f.legend)





#### Company Analysis ####---------------------------------------------
#### Figure 12 #### ---------------------------------------------------

names<-data.frame(prefix=unique(unlist(sapply(clcl, function(x){return(V(x[[1]])$name)}))))
companies<-read.csv("C:/Users/yasse/Downloads/companies.csv")
m<-as.numeric(unique(na.omit(match(companies$company,names$prefix))))
#names$prefix[m]
ncomp<-data.frame(prefix=names[m,])

ncomp<-ncomp[-which(ncomp$prefix=='usa'|ncomp$prefix=='icra'|ncomp$prefix=='asia'|
                      ncomp$prefix=='australia'|ncomp$prefix=='aerospace_defense'|
                      ncomp$prefix=='africa'|ncomp$prefix=='fd'|ncomp$prefix=='alignable'|
                      ncomp$prefix=='royal_bank_of_canada'),]
ncomp<-data.frame('prefix'=ncomp)

ncomp1<-data.frame(prefix=ncomp$prefix, 'tag'=1)

clusterExport(mcl, list('ncomp1',"%>%","separate"))

topcomp<-parSapply(mcl,partitionhuge, function(x,k=10){
  names<-x$clusters[,1]
  ec<-data.frame("key"=names)
  ec<-ec %>%
    separate(key, c("prefix", "kpi"), "-")
  #mer<-as.numeric(unique(na.omit(match(ncomp$prefix,ec$prefix))))
  ec$order<-1:nrow(ec)
  
  merged<-merge(ec, ncomp1, by='prefix',all.x = T)
  merged<-merged[order(merged$order),]
  mer<-which(merged$tag==1)
  #mer<-which(!is.na(merged$tag))
  
  
  #de<-strength(x$graph,vids = V(x$graph))
  de<-igraph::betweenness(x$graph, directed = F, weights=NA)
  de<-de[mer]
  
  se<-igraph::strength(x$graph)
  se<-se[mer]
  
  return(cbind(names[mer][order(se, decreasing=T)][1:k],names[mer][order(de, decreasing=T)][1:k]))
})


a<-sort(table(topcomp[1:3,1:112]),decreasing = T)[1:10]
a1<-sort(table(topcomp[1:3,113:ncol(topcomp)]),decreasing = T)[1:10]

b<-sort(table(topcomp[11:13,1:112]),decreasing = T)[1:10]
b1<-sort(table(topcomp[11:13,113:ncol(topcomp)]),decreasing = T)[1:10]

a<-a/length(1:112)
b<-b/length(1:112)

a1<-a1/(ncol(topcomp)-length(1:112))
b1<-b1/(ncol(topcomp)-length(1:112))

a<-data.frame(a)
b<-data.frame(b)

a1<-data.frame(a1)
b1<-data.frame(b1)

f<-cbind(a,b)
f1<-cbind(a1,b1)


col.c<-c("#5dcf77",
         "#9149d6",
         "#74d243",
         "#d243cb",
         "#bbc340",
         "#5661de",
         "#629732",
         "#db3d92",
         "#63cfa2",
         "#da4630",
         "#61c2e3",
         "#dca236",
         "#6787df",
         "#457536",
         "#e073cd",
         "#adc27e",
         "#9e3f91",
         "#539f88",
         "#d04569",
         "#8dd0c8",
         "#7b5db7",
         "#856f2c",
         "#cf99db",
         "#a35229",
         "#4670a8",
         "#e28262",
         "#41657c",
         "#d9af80",
         "#7d638b",
         "#808662",
         "#9dafde",
         "#916354",
         "#5d92a6",
         "#a65b73",
         "#386a5b",
         "#d9a5b4")





potential<-na.omit(unique(c(apply(f[,seq(1,ncol(f)-1,2)],2,as.character), apply(f1[,seq(1,ncol(f)-1,2)],2,as.character))))
potential<-sort(potential)

dfcol<-data.frame("browser"=potential,'cl.hue'=col.c[1:length(potential)])
dfcol$eh<-1

gcol<-ggplot(data=dfcol, aes(eh, fill=browser))+
  geom_bar()+
  scale_fill_manual(values=brightness(dfcol$cl.hue, 0.8)) +
  guides(fill=guide_legend(nrow=4)) + 
  theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
        legend.title=element_text(size=11), 
        legend.text=element_text(size=10))+
  labs(fill='Nodes')
gcol  

f.legend<-get_legend(gcol)

plots<-list()



for (i in seq(1,ncol(f)-1,2)) {
  col.ch<-c(i,i+1)
  
  # l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
  #        rep("(c) Cluster Dominance",2))
  
  #l.n<-c(rep("(a) Cluster Size",2),rep("(b) Sum of Weights",2))
  
  #l.n<-c(rep("(a) Strength",2),rep("(b) Betweenness Centrality",2),
  #       rep("(c) Top Prefixes",2))
  
  l.n<-c(rep("(a) Strength",2), rep("(b) Betweenness Centrality",2))
  
  # l.n<-c(rep("(a) Edge Weights",2),rep("(b) Edge Betweenness",2),
  #        rep("(c) Cluster Dominance",2))
  
  thr<-0.1
  
  
  eh<-merge(na.omit(f[,c(col.ch)]),na.omit(f1[,c(col.ch)]), by='Var1', all=T)
  eh[is.na(eh)]<-0
  df<-rbind(as.matrix(eh[,c(1,2)]),as.matrix(eh[,c(1,3)]))
  df<-data.frame(df)
  df$period<-c(rep(2019,nrow(df)/2),rep(2020,nrow(df)/2))
  colnames(df)<-c("browser","cc","year")
  #df<-na.omit(df)
  df$cc<-as.numeric(df$cc)
  
  # Data for plot
  pdat = df %>% 
    group_by(year) %>% 
    arrange(browser) %>% 
    mutate(ccc = cc/sum(cc)) %>%
    # Get cumulative value of cc
    mutate(cc_cum = cumsum(ccc)-0.5*ccc) %>% 
    ungroup
  
  pdat1<-merge(pdat,dfcol)
  if(nrow(pdat1)==0){pdat$cl.hue<-dfcol$cl.hue[1:nrow(pdat)]}
  else(pdat<-pdat1)
  
  #854adb
  
  
  
  
  #fix the colors
  
  #only for links
  pdat$browser<-gsub(", ",",\n ", pdat$browser)
  
  
  d<-pdat[pdat$year==2020,]
  m<-which(d$cc>0.1)
  m<-unlist(lapply(m,function(x){rep(x,2)}))
  
  pdat$m<-NA
  pdat$m[1:length(m)]<-m
  
  pdat$ang<- -pdat$cc_cum*360+90
  #pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]<-pdat$ang[pdat$cc_cum>=0.25 & pdat$cc_cum<=0.75]-180
  pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]<-pdat$ang[pdat$cc_cum>=0.5 & pdat$cc_cum<=1]-180
  
  
  cl.hue<-unique(pdat$cl.hue)
  
  
  g<-NA
  
  g<-ggplot(data=pdat, aes(x=cc_cum, y=year, fill=browser)) +
    geom_tile(aes(width=ccc), colour="white", size=0.4) +
    geom_text(aes(label=replace(round(100*cc,1),which(cc<thr),""))
              , size=2.5, colour="white", position=position_nudge(y=-0.2)) +
    geom_text(data=pdat %>% filter(year==unique(year)[2]), size=3.5, 
              aes(label=replace(browser,which(cc<thr),""),angle=ang, colour=browser), 
              position=position_nudge(y=0.6)) +
    geom_text(data=pdat %>% filter(year==unique(year)[1]), size=3.5, 
              aes(label=replace(browser,c(m,which(cc<thr)),""),angle=(ang), colour=browser),
              position=position_nudge(y=0.6)) +
    scale_y_continuous(#breaks=min(pdat$year):max(pdat$year),
      #labels=c("Before","After")) +
      breaks=NULL) +
    coord_polar() +
    theme_void() +
    theme(axis.text.y=element_text(angle=0, colour="grey40", size=9),
          axis.ticks.y=element_line(),
          axis.ticks.length=unit(0.1,"cm")) +
    #guides(fill=FALSE, colour=FALSE) +
    guides(colour=FALSE) +
    #scale_fill_manual(values=brightness(col.c, 0.8)) +
    scale_fill_manual(values=brightness(cl.hue, 0.8)) +
    #scale_colour_manual(values=brightness(col.c, 1.5)) +
    scale_colour_manual(values=rep("black",(nrow(pdat)/2))) +
    #scale_fill_manual(values=sample(hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,70),
    #                  size=(nrow(pdat)/2))) +
    #scale_colour_manual(values=hcl(seq(0,2000,length=(nrow(pdat)/2))[1:(nrow(pdat)/2)],100,40)) +
    theme(legend.position="bottom", plot.title = element_text(size = 15, hjust=0.5, face = "bold"),
          legend.title=element_text(size=11), 
          legend.text=element_text(size=10)) +
    ggtitle(l.n[i])
  
  #print(g)
  #print(g1)
  #assign(paste("g",i,sep = ""),g)
  
  plots<-c(plots,list(g))
  
}

ggarrange(plotlist=plots, ncol=ncol(f)/2, nrow=1, legend = "bottom", 
          legend.grob = f.legend)





#### Figure 13 #### ---------------------------------------------------
cldcomp<-sapply(clcl, function(x){
  x<-x[[1]]
  y<-data.frame("prefix"=V(x)$name)
  y<-unique(merge(y,ncomp1))
  return(y[,1])
})

cldcomp1<-do.call("c",cldcomp[1:112])
cldcomp2<-do.call("c",cldcomp[113:length(cldcomp)])

a<-sort(table(as.character(cldcomp1)),decreasing = T)[1:10]
b<-sort(table(as.character(cldcomp2)),decreasing = T)[1:10]


a<-a/length(1:112)
b<-b/(length(cldcomp)-length(1:112))

a<-data.frame(a)
b<-data.frame(b)

f<-cbind(a,b)


ra<-c(rep("(a) Before Covid-19",2),rep("(b) After Covid-19",2))


f2<-c()
for (i in seq(1,ncol(f)-1,2)) {
  t<-na.omit(f[,c(i,i+1)][1:10,])
  f2<-c(f2, list(list(t,ra[i])))
}


par(mfrow=c(1,2))
sapply(f2,function(x){barplot(x[[1]][,2], col=1:nrow(x[[1]]), legend.text = x[[1]][,1],
                              args.legend = list(x='topright', inset=c(0.07,0.01)),
                              beside = F, xlab = paste(x[[2]]), cex.lab=1.5, ylim=c(0,1))})

par(mfrow=c(1,1))






#### Figure 14 #### ----------------------------------------------------

compmain<-'amazon'
k<-1
l<-1
amnodes<-c("amazon-revenue","amazon-demand","amazon-employment")

ifull<-1:length(partitionhuge)
ibefore<-1:112
iafter<-113:length(partitionhuge)

whereamnodes<-foreach(i=iafter, .packages = c('igraph','tidyr'), .combine = c)%dopar%{
  y<-clcl[[i]]
  namescl<-V(y[[1]])$name
  mer1<-which(namescl==compmain)[1]
  
  wheref<-'nowhere'
  
  #if(is.na(mer1)){
  x<-clhuge[[i]]
  z<-partitionhuge[[i]]$clusters
  clwhere<-c()
  for(j in amnodes){
    try(clwhere<-c(clwhere,z$cluster[z$indicator==j]))
  }
  
  if(length(clwhere)!=0){
    #wheref<-unique(x$prefix[clwhere+1])
    wheref<-unique(V(y[[1]])$name[clwhere+1])
  }
  #}
  return(wheref)
}

length(whereamnodes)

whereafnodes<-whereamnodes[whereamnodes!='nowhere']
length(whereafnodes)
whereafnodest<-sort(table(whereafnodes),decreasing = T)
sum(whereafnodest)
whereafnodest

barplot(prop.table(whereafnodest)[1:5], col=1:length(whereafnodest), legend.text = names(whereafnodest)[1:5],
        args.legend = list(x='topright', inset=c(0.1,0.1)),
        beside = F, xlab = "Clusters housing main Amazon nodes", cex.lab=1.3, ylim=c(0,1),
        xaxt='n')








#### Figure 15 #### --------------------------------------------------------

#astrazeneca
compmain<-'astrazeneca'
fbnodes<-c("astrazeneca-usage","astrazeneca-production","astrazeneca-risk")

wherefbnodes<-foreach(i=iafter, .packages = c('igraph','tidyr'), .combine = c)%dopar%{
  y<-clcl[[i]]
  namescl<-V(y[[1]])$name
  mer1<-which(namescl==compmain)[1]
  
  wheref<-'nowhere'
  
  #if(is.na(mer1)){
  x<-clhuge[[i]]
  z<-partitionhuge[[i]]$clusters
  clwhere<-c()
  for(j in fbnodes){
    try(clwhere<-c(clwhere,z$cluster[z$indicator==j]))
  }
  
  if(length(clwhere)!=0){
    #wheref<-unique(x$prefix[clwhere+1])
    wheref<-unique(V(y[[1]])$name[clwhere+1])
  }
  #}
  return(wheref)
}

length(wherefbnodes)

wherefbfnodes<-wherefbnodes[wherefbnodes!='nowhere']
wherefbfnodest<-sort(table(wherefbfnodes),decreasing = T)
sum(wherefbfnodest)
wherefbfnodest

par(mfrow=c(1,2))

barplot(prop.table(wherefbfnodest)[1:5], col=1:length(wherefbfnodest), legend.text = names(wherefbfnodest)[1:5],
        args.legend = list(x='topright', inset=c(0.1,0.1)),
        beside = F, xlab = "(a) Clusters housing main Astrazeneca nodes", cex.lab=1.5, ylim=c(0,1),
        xaxt='n')


#pfizer
compmain<-'pfizer'
fbnodes<-c("pfizer-production","pfizer-risk")

wherefbnodes<-foreach(i=iafter, .packages = c('igraph','tidyr'), .combine = c)%dopar%{
  y<-clcl[[i]]
  namescl<-V(y[[1]])$name
  mer1<-which(namescl==compmain)[1]
  
  wheref<-'nowhere'
  
  #if(is.na(mer1)){
  x<-clhuge[[i]]
  z<-partitionhuge[[i]]$clusters
  clwhere<-c()
  for(j in fbnodes){
    try(clwhere<-c(clwhere,z$cluster[z$indicator==j]))
  }
  
  if(length(clwhere)!=0){
    #wheref<-unique(x$prefix[clwhere+1])
    wheref<-unique(V(y[[1]])$name[clwhere+1])
  }
  if(length(wheref)!=1){wheref<-'split'}
  #}
  return(wheref)
}

length(wherefbnodes)

wherefbfnodes<-wherefbnodes[wherefbnodes!='nowhere']
length(wherefbfnodes)
wherefbfnodest<-sort(table(wherefbfnodes),decreasing = T)
sum(wherefbfnodest)
wherefbfnodest

barplot(prop.table(wherefbfnodest)[1:5], col=1:length(wherefbfnodest), legend.text = names(wherefbfnodest)[1:5],
        args.legend = list(x='topright', inset=c(0.1,0.1)),
        beside = F, xlab = "(b) Clusters housing main Pfizer nodes", cex.lab=1.5, ylim=c(0,1),
        xaxt='n')

par(mfrow=c(1,1))
