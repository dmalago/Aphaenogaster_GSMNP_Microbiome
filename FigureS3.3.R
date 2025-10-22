library('ape')
library('phytools')
library('stringr')
library('ggtree')
library('rstatix')
library('caper')
library('spdep')
library('raster')
library('dplyr')
library('betareg')
library('ggplot2')
library('adegenet')
library('geosphere')
library('adespatial')
library('vegan')
library('ggpubr')
library('qiime2R')
library('ggpattern')
library("gridExtra")
library("grid")
library("cowplot")

Wtree<-read.tree(file='Wolbachia-rooted_tree_out/tree.nwk')
Wtree$edge.length <- sqrt(Wtree$edge.length+1e-2)

fixme<-c()
for (k in 1:length(Wtree$tip.label)){
  if (grepl('_',Wtree$tip.label[k])){
    temp<-Wtree$tip.label[k]
  }else{
    temp<-tolower(Wtree$tip.label[k])
  }
  fixme<-c(fixme,temp)
}

Wtree$tip.label<-fixme
cols<-c('black','grey',rep('yellow',5),rep('green',5),rep('blue',13),rep('purple',4),rep('orange',2),rep('red',32))
ggtree(Wtree,layout='rectangular')+geom_tiplab(size=2,offset=0.1)+xlim(-1, 3)+ geom_tippoint(size=3,color='black',fill=cols,shape=21)

