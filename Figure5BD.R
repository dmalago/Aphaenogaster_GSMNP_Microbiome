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

tree_choice<-'CO1'
set.seed(1)

#Read in metadata
meta<-read.csv('metadata.csv',header=TRUE)

data_path<-"./Soil_tifs"
files2 <- list.files(path = data_path, pattern = "*.tif", full.names = TRUE)
bioclim_layers2 <- stack(files2)

data_path2<-"./Soil_tifs2"
files3 <- list.files(path = data_path2, pattern = "*.tif", full.names = TRUE)
bioclim_layers3 <- stack(files3)

data_path3<-"./Soil_tifs3"
files4 <- list.files(path = data_path3, pattern = "*.tif", full.names = TRUE)
bioclim_layers4 <- stack(files4)

#Read in host (Aphaenogaster) phylogeny
if (tree_choice == 'CO1'){
  host_tree<-read.tree("Full_Trees/COI_TREE_FULL/iqt.contree") ## tree
  host_tree<-root(host_tree, outgroup = "A_umphreyi", resolve.root = TRUE)
}else if (tree_choice == 'CAD'){
  host_tree<-read.tree("Full_Trees/CAD_TREE_FULL/iqt.contree") ## tree
  host_tree<-root(host_tree, outgroup = "A_umphreyi", resolve.root = TRUE)
}else{
  print('WARNING')
}

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

host_tree<-drop.tip(host_tree, 'A_umphreyi')

#Normalized Aphaenogaster tree
host_tree2<-host_tree
host_tree
# Apply sqrt transformation, adding a small constant
host_tree2$edge.length <- sqrt(host_tree2$edge.length+1e-2)

tipper<-host_tree$tip.label

long<-c()
lat<-c()
elevation<-c()
rownamer<-c()
for (k in 1:length(tipper)){
  long<-c(long,meta$Long[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
  lat<-c(lat,meta$Lat[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
  elevation<-c(elevation,meta$Elevation[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
  rownamer<-c(rownamer,tipper[k])
}

lonlat<-data.frame(long,lat)
colnames(lonlat)=c('longitude','latitude')

values2 <- raster::extract(bioclim_layers2, lonlat)
bioclim_values2 <- as.data.frame(values2)

values3 <- raster::extract(bioclim_layers3, lonlat)
bioclim_values3 <- as.data.frame(values3)

values4 <- raster::extract(bioclim_layers4, lonlat)
bioclim_values4 <- as.data.frame(values4)

bioclim_values<-bind_cols(bioclim_values2, bioclim_values3, bioclim_values4,as.data.frame(elevation),)
rownames(bioclim_values)<-tipper
colnames(bioclim_values) <- c(" min_temp_soil", "soil_moisture_summer", "soil_temp_range","elevation")

# Update specific rows and columns by index
bioclim_values[which(grepl('OCO_',rownames(bioclim_values))==TRUE), 1] <- -5.193484
bioclim_values[which(grepl('OCO_',rownames(bioclim_values))==TRUE), 2] <-0.071459
bioclim_values[which(grepl('OCO_',rownames(bioclim_values))==TRUE), 3] <-32.10958
bioclim_values[which(grepl('PK_',rownames(bioclim_values))==TRUE),2]<- 0.160193

bioclim_values<-scale(bioclim_values)
df_bioclim_values<-data.frame(bioclim_values)

e1<-c()
e2<-c()
e3<-c()
e4<-c()
rownamer2<-c()
for (k in 1:length(tipper)){
  e1<-c(e1,df_bioclim_values$X.min_temp_soil[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
  e2<-c(e2,df_bioclim_values$soil_temp_range[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
  e3<-c(e3,df_bioclim_values$soil_moisture_summer[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
  e4<-c(e4,df_bioclim_values$elevation[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
  rownamer2<-c(rownamer2,tipper[k])
}
edf<-data.frame(e1,e1+e2,e3,e4)


sedf<-scale(edf)


colnames(sedf)<-c('minT','maxT','moisture','elevation')

pca<-prcomp(sedf)

environment<-pca$x[,1]
names(environment)<-tipper

if (tree_choice=='CAD'){
  yellow<-clade.members(x=MRCA(host_tree2,'CL_5C_CARO', 'PK_6B_PICEA'),phy = host_tree2,tip.labels=TRUE)
  green<-clade.members(x=MRCA(host_tree2,'TM_2E_RUDIS', 'AG_5A_PICEA'),phy = host_tree2,tip.labels=TRUE)
}else{
  purple<-clade.members(x=MRCA(host_tree2,'EM_9D_PICEA', 'RC_5B_PICEA'),phy = host_tree2,tip.labels=TRUE)
  blue<-clade.members(x=MRCA(host_tree2,'BMM_5B_PICEA', 'TC_6A_PICEA'),phy = host_tree2,tip.labels=TRUE)
  red<-setdiff(host_tree2$tip.label,c('A_umphreyi',purple,blue))
}

if (tree_choice=='CO1'){
a<-environment[which(names(environment) %in% red)]
b<-environment[which(names(environment) %in% purple)]
c<-environment[which(names(environment) %in% blue)]

diversity_df<-data.frame(c(rep('Red',length(a)),rep('Purple',length(b)),rep('Blue',length(c))),c(a,b,c))
colnames(diversity_df)<-c('Group','Environment')
diversity_df$Group<-factor(diversity_df$Group,levels=c('Red','Purple','Blue'))

kw_env<-kruskal.test(Environment ~ Group, data = diversity_df)
if (kw_env$p.value<0.05){
  pw_env<-pairwise.wilcox.test(diversity_df$Environment, diversity_df$Group,p.adjust.method = "BH")
  print(pw_env)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Do you want to show p-values for non-significant pairwise comparisons? (yes or no)
p_value_ns<-'no'
pshow<-'num'

##########################env Violin Plot#########################################################################################################################################

#Define violin plot
p_env <- ggplot(diversity_df, aes(x=Group, y=Environment)) + geom_violin(aes(fill=Group),scale='width') +labs(x ="Group", y = "Environment")
#Choose the size of font for the axes titles and labels
p_env<-p_env+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 13))
#Choose the size of font for the legend title and lables
p_env<-p_env+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_env<-p_env+scale_fill_manual(values=c("red","purple", "blue"))
#Add boxplots inside the violins
p_env<-p_env+geom_boxplot(aes(fill=Group),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)

#If there are significant differences in env between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_env$p.value<0.01){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-5.5
  y_step<-0.5
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_env$p.value))){
    for (j in 1:k){
      if (rownames(pw_env$p.value)[k]!=colnames(pw_env$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_env$p.value[k,j]<0.1 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_env$p.value)[k])
          group2<-c(group2,colnames(pw_env$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_env$p.value[k,j])),5)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_env<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_env<-p_env+stat_pvalue_manual(stat.test_env,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_env

}else{
  a<-environment[which(names(environment) %in% green)]
  b<-environment[which(names(environment) %in% yellow)]

  diversity_df<-data.frame(c(rep('Green',length(a)),rep('Yellow',length(b))),c(a,b))
  colnames(diversity_df)<-c('Group','Environment')
  diversity_df$Group<-factor(diversity_df$Group,levels=c('Green','Yellow'))

  kw_env<-kruskal.test(Environment ~ Group, data = diversity_df)
  if (kw_env$p.value<0.05){
    pw_env<-pairwise.wilcox.test(diversity_df$Environment, diversity_df$Group,p.adjust.method = "BH")
    print(pw_env)
  }


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  #Do you want to show p-values for non-significant pairwise comparisons? (yes or no)
  p_value_ns<-'no'
  pshow<-'num'

  ##########################env Violin Plot#########################################################################################################################################

  #Define violin plot
  p_env <- ggplot(diversity_df, aes(x=Group, y=Environment)) + geom_violin(aes(fill=Group),scale='width') +labs(x ="Group", y = "Environment")
  #Choose the size of font for the axes titles and labels
  p_env<-p_env+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 13))
  #Choose the size of font for the legend title and lables
  p_env<-p_env+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
  #Choose the violin colors for each group
  p_env<-p_env+scale_fill_manual(values=c("green", "yellow"))
  #Add boxplots inside the violins
  p_env<-p_env+geom_boxplot(aes(fill=Group),width=0.1)
  #Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)

  #If there are significant differences in env between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
  if (kw_env$p.value<0.01){
    #Names of the groups you're comparing
    group1<-c()
    group2<-c()
    #p-value for the pairwise comparisons
    p.adj<-c()
    #locations of the p-value brackets
    ypos<-c()
    new_y<-3.5
    y_step<-0.5
    #For each pairwise comparison...
    for (k in 1:length(rownames(pw_env$p.value))){
      for (j in 1:k){
        if (rownames(pw_env$p.value)[k]!=colnames(pw_env$p.value)[j]){
          #If there is a significant difference or you want to also show non-significant differences...
          if (pw_env$p.value[k,j]<0.1 || p_value_ns=='yes'){
            #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
            group1<-c(group1,rownames(pw_env$p.value)[k])
            group2<-c(group2,colnames(pw_env$p.value)[j])
            p.adj<-round(c(p.adj,as.numeric(pw_env$p.value[k,j])),5)
            #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
            ypos<-c(ypos,new_y)
            new_y<-new_y+y_step
          }
        }
      }
    }
    pstar<-c()
    for (k in 1:length(p.adj)){
      if (p.adj[k]<=0.001){
        pstar<-c(pstar,'***')
      }
      else if (p.adj[k]<=0.01){
        pstar<-c(pstar,'**')
      }
      else if (p.adj[k]<=0.05){
        pstar<-c(pstar,'*')
      }
      else if (p.adj[k]<=0.1){
        pstar<-c(pstar,'.')
      }
      else{
        pstar<-c(pstar,'ns')
      }
    }
    if (pshow=='star'){
      pdisplay<-"{pstar}"
    }
    else{
      pdisplay<-"p = {p.adj}"
    }
    #Create your tibble (what is needed for the stat_pvalue_manual function)
    stat.test_env<-as_tibble(data.frame(group1,group2,p.adj))
    #Add the pairwise comparisons to your plot
    p_env<-p_env+stat_pvalue_manual(stat.test_env,label=pdisplay,y.position=ypos,size=5)
  }
  #Make your plot
  p_env
}



