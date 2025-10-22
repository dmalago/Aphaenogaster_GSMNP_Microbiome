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


cophylo1<-cophylo(host_treeCO11,host_treeCAD1)
cophylo2<-cophylo(host_treeCO12,host_treeCAD2)

col_list<-cophylo2$assoc[,1]
col_list[which(cophylo2$assoc[,1] %in% 'A_umphreyi')]<-'black'

f <- colorRamp(c("purple","lavender","peachpuff","darkorange"))

col_long<-c()
for (k in 1:length(col_list)){
  col_long<-c(col_long,meta$Long[which(grepl(paste0(str_split(col_list[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
}


col_long[which(cophylo2$assoc[,1] %in% 'A_umphreyi')]<-col_long[1]

long<-(col_long-min(col_long))/max(col_long-min(col_long))

col_list[which(cophylo2$assoc[,1] %in% 'A_umphreyi')]<-'black'
for (k in 1:length(col_list)){
  if (k!=which(cophylo2$assoc[,1] %in% 'A_umphreyi')){
    col_list[k]<-rgb(f(long[k])/255)
  }
}


plot.cophylo(cophylo2,fsize=c(0.3,0.3),
     link.lty="solid",link.col=col_list, link.type='curved',tip.color=col_list,tip.color=col_list,link.lwd=2)


#CAD Structure

#ggtree(host_treeCAD2)+geom_nodelab(size=2)+geom_tiplab(size=2)


lightyellow<-clade.members(x=MRCA(host_treeCAD2,'CC_1D_CARO', 'PK_6B_PICEA'),phy = host_treeCAD2,tip.labels=TRUE)
darkyellow<-clade.members(x=MRCA(host_treeCAD2,'CL_5C_CARO', 'BM_7I_PICEA'),phy = host_treeCAD2,tip.labels=TRUE)

lightgreen<-clade.members(x=MRCA(host_treeCAD2,'AG_5A_PICEA', 'OCO_4C_CARO'),phy = host_treeCAD2,tip.labels=TRUE)
darkgreen<-clade.members(x=MRCA(host_treeCAD2,'CL_4A_PICEA', 'TM_2E_RUDIS'),phy = host_treeCAD2,tip.labels=TRUE)

allgreen<-union(lightgreen,darkgreen)
allyellow<-union(lightyellow,darkyellow)

#COI Structure

#ggtree(host_treeCO12)+geom_nodelab(size=2)+geom_tiplab(size=2)

purple<-clade.members(x=MRCA(host_treeCO12,'EM_9D_PICEA', 'RC_5B_PICEA'),phy = host_treeCO12,tip.labels=TRUE)
blue<-clade.members(x=MRCA(host_treeCO12,'BMM_5B_PICEA', 'TC_6A_PICEA'),phy = host_treeCO12,tip.labels=TRUE)
red<-setdiff(host_treeCO12$tip.label,c('A_umphreyi',purple,blue))

#Make phyloseq object from qza files
data<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza","taxonomy.qza",metadata = "metadata.tsv")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")

host_tree<-read.tree("Full_Trees/COI_TREE_FULL/iqt.contree") ## tree
host_tree<-root(host_tree, outgroup = "A_umphreyi", resolve.root = TRUE)



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

purple<-clade.members(x=MRCA(host_tree,'EM_9D_PICEA', 'RC_5B_PICEA'),phy = host_tree,tip.labels=TRUE)
blue<-clade.members(x=MRCA(host_tree,'BMM_5B_PICEA', 'TC_6A_PICEA'),phy = host_tree,tip.labels=TRUE)
red<-setdiff(host_tree$tip.label,c('A_umphreyi',purple,blue))


#Drop tree tips that aren't in the microbiome dataset
tips_to_drop <- setdiff(host_tree$tip.label, colnames(otu_table(data)))
host_tree <- drop.tip(host_tree, setdiff(host_tree$tip.label, colnames(otu_table(data))))
tree_tips<-host_tree$tip.label

#Drop microbiome samples that aren't in the Aphaenogaster tree
phyloseq_samples <- sample_names(data)
matched_samples <- intersect(tree_tips, phyloseq_samples)
data <- prune_samples(matched_samples, data)
data<- prune_taxa(taxa_sums(data) > 0, data)

# Find tip labels for the Aphaenogaster tree
tree_tips <- host_tree$tip.label

# Find sample names from the phyloseq object
phyloseq_samples <- sample_names(data)
tt<-phyloseq_samples[which(phyloseq_samples %in% union(allgreen,allyellow))]
ss<-phyloseq_samples[which(!(phyloseq_samples %in% union(allgreen,allyellow)))]
host_tree<-drop.tip(host_tree,ss)
data<-prune_samples(tt,data)

# Find the intersection of tree tips and phyloseq sample names
matched_tips <- intersect(tree_tips, phyloseq_samples)
if (length(matched_tips)!=length(tree_tips)){print('WARNING!')}


#Rarefy your data, here I am rarefying to 22533 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=22533, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified)

#Remove samples that are contols
data_rarified<-subset_samples(data_rarified,Type != "Negative")

#Find OTU table in biom file
OTU_biom<-otu_table(data_rarified)


#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-phy_tree(data)

#Rename the tree tips to match the ASV table names
cut_names<-c()
for (k in 1:length(taxa_names(tree_file))){
  #Depending on the Qiime2 output, you may need to split the names with a space or with an underscore
  #cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],' ')[[1]][1])
  cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],'_')[[1]][1])
}
taxa_names(tree_file)<-cut_names

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)
colnames(TAX)<-c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
TAX<-gsub('p__', '', TAX)
TAX<-gsub('c__', '', TAX)
TAX<-gsub('o__', '', TAX)
TAX<-gsub('f__', '', TAX)
TAX<-gsub('g__', '', TAX)
TAX<-gsub('s__', '', TAX)



metadata<-sample_data(data_rarified)

temp_phylo<-phyloseq(OTU_biom,TAX,tree_file)
Mphylo<-temp_phylo

M2phylo<-phyloseq(otu_table(Mphylo),tax_table(Mphylo),sample_data(metadata),phy_tree(Mphylo))


Group<-sample_names(M2phylo)
Group[intersect(which(Group %in% red),which(Group %in% allgreen))]<-'R-G'
Group[intersect(which(Group %in% red),which(Group %in% allyellow))]<-'R-Y'
Group[intersect(which(Group %in% purple),which(Group %in% allyellow))]<-'P-Y'
Group[intersect(which(Group %in% purple),which(Group %in% allgreen))]<-'P-G'
Group[intersect(which(Group %in% blue),which(Group %in% allgreen))]<-'B-G'
Group[intersect(which(Group %in% blue),which(Group %in% allyellow))]<-'B-Y'

GY<-sample_names(M2phylo)
GY[which(GY %in% allgreen)]<-'Green'
GY[which(GY %in% allyellow)]<-'Yellow'

RPB<-sample_names(M2phylo)
RPB[which(RPB %in% red)]<-'Red'
RPB[which(RPB %in% purple)]<-'Purple'
RPB[which(RPB %in% blue)]<-'Blue'



rMphylo<-M2phylo
sample_data(rMphylo)$Group<-Group
sample_data(rMphylo)$GY<-GY
sample_data(rMphylo)$RPB<-RPB
#rMGphylo<-rMphylo
rMGphylo<-tax_glom(rMphylo,taxrank = 'Genus')

Wreads<-otu_table(rMGphylo)[which(tax_table(rMGphylo)[,6]=='Wolbachia'),]
#Wreads<-otu_table(rMGphylo)[which(rownames(tax_table(rMGphylo))==tstrains[1]),]

diversity_df<-data.frame(sample_data(rMGphylo)$Group,factor(sample_data(rMGphylo)$GY),factor(sample_data(rMGphylo)$RPB),t(Wreads)/22533,t(Wreads))
colnames(diversity_df)<-c('Group','GY','RPB','Wolbachia','W')
diversity_df$Group<-factor(diversity_df$Group,levels=c('R-G','R-Y','P-G','P-Y','B-G','B-Y'))



#Define violin plot
p_ants <- ggplot(diversity_df, aes(x=Group, y=Wolbachia))+geom_violin_pattern(aes(fill=Group,pattern_fill = Group),pattern='stripe',scale='width',pattern_spacing=0.05,pattern_density=0.5,width=0.75) + theme(axis.title.x = element_blank())+ylab('Wolbachia Relative Abundance')
#Choose the size of font for the axes titles and labels
p_ants<-p_ants+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_ants<-p_ants+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_ants<-p_ants+scale_fill_manual(values=c('red', 'red','purple','purple','blue','blue'))+scale_pattern_fill_manual(values = c('green','yellow', 'green','yellow','green','yellow'))
#Add boxplots inside the violins
p_ants<-p_ants+geom_boxplot(aes(fill=Group),width=0.1) + theme(axis.title.x = element_blank())

p_ants




fit_betareg <- betareg(Wolbachia ~ GY*RPB, data = diversity_df)
AIC(fit_betareg)
summary(fit_betareg)
fit_betareg <- betareg(Wolbachia ~ RPB, data = diversity_df)
AIC(fit_betareg)
summary(fit_betareg)
fit_betareg <- betareg(Wolbachia ~ GY, data = diversity_df)
AIC(fit_betareg)
summary(fit_betareg)
Ronly<-RPB
Ponly<-RPB
Ronly[Ronly=='Purple']<-'Blue'
Ponly[Ponly=='Red']<-'Blue'
fit_betareg <- betareg(Wolbachia ~ GY*Ronly, data = diversity_df)
AIC(fit_betareg)
fit_betareg <- betareg(Wolbachia ~ GY*Ponly, data = diversity_df)
AIC(fit_betareg)

fit_betareg <- betareg(Wolbachia ~ GYYellow, data = diversity_df)
