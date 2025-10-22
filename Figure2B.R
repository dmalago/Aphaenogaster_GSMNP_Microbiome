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

meta<-read.csv('metadata.csv',header=TRUE)

host_treeCO1<-read.tree("Full_Trees/COI_TREE_FULL/iqt.contree") ## tree
host_treeCO1<-root(host_treeCO1, outgroup = "A_umphreyi", resolve.root = FALSE)

#Wtree<-read.tree(file='WolbachiaA-rooted_tree_out/tree.nwk')


host_treeCAD<-read.tree("Full_Trees/CAD_TREE_FULL/iqt.contree") ## tree
host_treeCAD<-root(host_treeCAD, outgroup = "A_umphreyi", resolve.root = FALSE)

for (k in 1:2){
  if (k==1){host_tree<-host_treeCO1}else{host_tree<-host_treeCAD}

  #Relabel the Aphaenogaster tree tips
  tip_names <- host_tree$tip.label
  simplified_tip_names <- gsub("([0-9]+)$", "", tip_names)
  simplified_tip_names <- gsub("(A)(carolinensis|rudis|fulva|tennesseeinsis|picea)", "\\2", simplified_tip_names) # Remove the leading "A" from specific taxa, including picea
  simplified_tip_names <- gsub("_$", "", simplified_tip_names)
  simplified_tip_names <- gsub("picea", "PICEA", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("rudis", "RUDIS", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("carolinensis", "CARO", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("fulva", "FULVA", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("OCO_4ASP_tennesseeinsis", "OCO_4_ASP_2_TENN", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_53_PK7A_PICEA", "PK_7A_PICEA", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_49_AG7E_RUDIS", "AG_7E_RUDIS", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("FR_6ASP_PICEA", "FR_6_ASP_PICEA", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_51_CCD_CARO", "CC_1D_CARO", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_54_TCI_CARO", "TC_1I_CARO", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_50_BMC_PICEA", "BM_1C_PICEA", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_68_TCI_CARO", "TC_1I_CARO", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_63_AG7E_RUDIS", "AG_7E_RUDIS", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_65_CCD_CARO", "CC_1D_CARO", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_64_BMC_PICEA", "BM_1C_PICEA", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_67_PK7A_PICEA", "PK_7A_PICEA", simplified_tip_names, ignore.case = FALSE)
  simplified_tip_names <- gsub("DM_66_CCB_RUDIS", "CC_1B_RUDIS", simplified_tip_names, ignore.case = FALSE)
  host_tree$tip.label <- simplified_tip_names



  #Normalized Aphaenogaster tree
  host_tree2<-host_tree
  # Apply sqrt transformation, adding a small constant
  host_tree2$edge.length <- sqrt(host_tree2$edge.length+1e-2)


if (k==1){
  host_treeCO11<-host_tree
  host_treeCO12<-host_tree2
}else{
  host_treeCAD1<-host_tree
  host_treeCAD2<-host_tree2
}
}

####Wolbachia
#Make phyloseq object from qza files
data<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza","taxonomy.qza",metadata = "metadata.tsv")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")

#Rarefy your data, here I am rarefying to 22533 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=22533, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified)

#Remove samples that are contols
data_rarified<-subset_samples(data_rarified,Type != "Negative")

#Rename to keep in the outgroup
tax_table(data_rarified)[which(rownames(tax_table(data_rarified))=='73b413c294026326492469ffdf8f1373'),6]<-'Wolbachia'

Wdata<-subset_taxa(data_rarified,Genus =='Wolbachia')



microbiome<-sample_names(Wdata)
host1<-host_treeCO12$tip.label
host2<-host_treeCAD2$tip.label

good<-c(intersect(intersect(microbiome,host1),host2),'A_umphreyi')

host_treeCO12good<-drop.tip(host_treeCO12,host1[which(!(host1 %in% good))])
host_treeCAD2good<-drop.tip(host_treeCAD2,host2[which(!(host2 %in% good))])
Wdatagood<-prune_samples(good,Wdata)
Wv<-colSums(otu_table(Wdatagood))
#Wv<-log(Wv+1)/max(log(Wv+1))
Wv<-sqrt(Wv)

cophylo2<-cophylo(host_treeCO12good,host_treeCAD2good,plot=FALSE)

col_list<-cophylo2$assoc[,1]
col_list[which(cophylo2$assoc[,1] %in% 'A_umphreyi')]<-'black'

f <- colorRamp(c("blue",'lavender', "pink","red"))

tt<-sample_names(Wdatagood)


for (k in 1:length(col_list)){
  #if (col_list[k]!='black'){
    if (col_list[k] %in% tt){
      hit<-which(tt==col_list[k])
 col_list[k]<-rgb(f((Wv[hit]/sqrt(22533)))/255)
  }else(col_list[k]<-'lightgrey')
}


cop1<-cophylo2[[1]][[1]]$tip.label
cop2<-cophylo2[[1]][[2]]$tip.label



tip_colors1<-rep('black',length(cop1))
tip_colors2<-rep('black',length(cop2))


tip_colors1[which(cop1 %in% red)] <- "red"
tip_colors1[which(cop1 %in% blue)] <- "blue"
tip_colors1[which(cop1 %in% purple)] <- "purple"
tip_colors2[which(cop2 %in% allgreen)] <- "green"
tip_colors2[which(cop2 %in% allyellow)] <- "yellow"

plot(cophylo2,link.type='curved',link.col=col_list,fsize=0.5,link.lty='solid',link.lwd=2)+tiplabels.cophylo(pch=16, col=tip_colors2,which=c("right"))+tiplabels.cophylo(pch=16, col=tip_colors1,which=c("left"))

