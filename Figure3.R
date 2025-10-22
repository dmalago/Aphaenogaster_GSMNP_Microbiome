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
library('pheatmap')

which_tree<-'CAD'

meta<-read.csv('metadata.csv',header=TRUE)

host_treeCO1<-read.tree("Full_Trees/COI_TREE_FULL/iqt.contree") ## tree
host_treeCO1<-root(host_treeCO1, outgroup = "A_umphreyi", resolve.root = FALSE)

Wtree<-read.tree(file='WolbachiaA-rooted_tree_out/tree.nwk')


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


Wtree$edge.length <- sqrt(Wtree$edge.length+1e-2)

Wtree<-root(Wtree, outgroup = "73b413C294026326492469ffdf8f1373", resolve.root = FALSE)

#Normalized Aphaenogaster tree
Wtree2<-Wtree
# Apply sqrt transformation, adding a small constant
Wtree2$edge.length <- sqrt(Wtree2$edge.length+1e-2)

if (which_tree=='CO1'){
  host_tree<-host_treeCO12
}else{
  host_tree<-host_treeCAD2
}

#Drop tree tips that aren't in the microbiome dataset
host_tree_preserve <- drop.tip(host_tree, 'A_umphreyi')

tips_to_drop <- setdiff(host_tree$tip.label, colnames(otu_table(Wdata)))
host_tree <- drop.tip(host_tree, setdiff(host_tree$tip.label, colnames(otu_table(data_rarified))))
tree_tips<-host_tree$tip.label

#Drop microbiome samples that aren't in the Aphaenogaster tree
phyloseq_samples <- sample_names(Wdata)
matched_samples <- intersect(tree_tips, phyloseq_samples)
Wdata <- prune_samples(matched_samples, Wdata)
#Wdata<- prune_taxa(taxa_sums(Wdata) > 0, Wdata)

Wtable<-otu_table(Wdata)

Wtable2<-Wtable
Wtable2[which(Wtable2<0.005*colSums(otu_table(data_rarified))[1])]<-0


# Find tip labels for the Aphaenogaster tree
tree_tips <- host_tree$tip.label

# Find sample names from the phyloseq object
phyloseq_samples <- sample_names(data_rarified)

# Find the intersection of tree tips and phyloseq sample names
matched_tips <- intersect(tree_tips, phyloseq_samples)
if (length(matched_tips)!=length(tree_tips)){print('WARNING!')}

rowme<-c()
colme<-c()

for (k in 1:46){
  for (j in 1:33){
    if (Wtable2[j,k]>0){
      rowme<-c(rowme,rownames(Wtable2[j,k]))
      colme<-c(colme,colnames(Wtable2[j,k]))
    }
  }
}

dimmer<-cbind(colme,rowme)

fixme<-c()
for (k in 1:length(Wtree2$tip.label)){
  if (grepl('_',Wtree2$tip.label[k])){
    temp<-Wtree2$tip.label[k]
  }else{
    temp<-tolower(Wtree2$tip.label[k])
  }
  fixme<-c(fixme,temp)
}

Wtree2$tip.label<-fixme

rorder<-c()
for (k in 1:length(Wtree2$tip.label)){
  rorder<-c(rorder,which(rownames(Wtable)==Wtree2$tip.label[k]))
}

corder<-c()
for (k in 1:length(host_tree$tip.label)){
  corder<-c(corder,which(colnames(Wtable)==host_tree$tip.label[k]))
}

Wtable<-as(otu_table(Wdata), "matrix")
ord.mat <- as(Wtable[rorder,corder],"matrix")

Wdist<-as.dist(cophenetic.phylo(Wtree2))
Adist<-as.dist(cophenetic.phylo(host_tree))

hc<-hclust(Adist)

callback_rev <- function(hc, Adist){
  hc<-reorder(hc,1:1:46)
  hc
}

scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
pheatmap(log10(ord.mat+1)-1 ,cluster_rows = as.hclust(chronos(Wtree2)),cluster_cols=as.hclust(chronos(host_tree)),fontsize=5)
