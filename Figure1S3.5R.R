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

set.seed(1)

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



#Make phyloseq object from qza files
data<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza","taxonomy.qza",metadata = "metadata.tsv")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")
#data<-subset_taxa(data, Genus != "Wolbachia")
#data<-subset_taxa(data, Genus != "Spiroplasma")
#data<-subset_taxa(data, Genus != "Entomoplasma")
#data<-subset_taxa(data, Genus != "Candidatus_Sulcia")

#Rarefy your data, here I am rarefying to 22533 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=22533, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified)

#Remove samples that are contols
data_rarified<-subset_samples(data_rarified,Type != "Negative")


lonlat<-data.frame(sample_data(data_rarified)$Long,sample_data(data_rarified)$Lat)
colnames(lonlat)=c('longitude','latitude')

values2 <- raster::extract(bioclim_layers2, lonlat)
bioclim_values2 <- as.data.frame(values2)

values3 <- raster::extract(bioclim_layers3, lonlat)
bioclim_values3 <- as.data.frame(values3)

values4 <- raster::extract(bioclim_layers4, lonlat)
bioclim_values4 <- as.data.frame(values4)

elevation<-sample_data(data_rarified)$Elevation

tipper<-sample_names(data_rarified)


bioclim_values<-bind_cols(bioclim_values2, bioclim_values3, bioclim_values4,as.data.frame(elevation),)
rownames(bioclim_values)<-tipper
colnames(bioclim_values) <- c(" min_temp_soil", "soil_moisture_summer", "soil_temp_range","elevation")

# Update specific rows and columns by index
bioclim_values[which(grepl('OCO_',rownames(bioclim_values))==TRUE), 1] <- -5.193484
bioclim_values[which(grepl('OCO_',rownames(bioclim_values))==TRUE), 2] <-0.071459
bioclim_values[which(grepl('OCO_',rownames(bioclim_values))==TRUE), 3] <-32.10958
bioclim_values[which(grepl('PK_',rownames(bioclim_values))==TRUE),2]<- 0.160193

#bioclim_values<-scale(bioclim_values)
df_bioclim_values<-data.frame(bioclim_values)


# e1<-c()
# e2<-c()
# e3<-c()
# e4<-c()
# rownamer2<-c()
# for (k in 1:length(tipper)){
#   e1<-c(e1,df_bioclim_values$X.min_temp_soil[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
#   e2<-c(e2,df_bioclim_values$soil_temp_range[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
#   e3<-c(e3,df_bioclim_values$soil_moisture_summer[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
#   e4<-c(e4,df_bioclim_values$elevation[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
#   rownamer2<-c(rownamer2,tipper[k])
# }
# edf<-data.frame(e1,e1+e2,e3,e4)
#colnames(edf)<-c('min_soil_temp','max_soil_temp','soil_moisture','elevation')
sedf<-scale(df_bioclim_values)
colnames(sedf)<-c('minT','maxT','moisture','elevation')

pca<-prcomp(sedf)

environment<-pca$x[,1]
site<-sample_data(data_rarified)$Site

dwr<-tax_glom(data_rarified,taxrank='Genus')
wolbachia<-otu_table(dwr)[which(tax_table(dwr)[,6]=='Wolbachia')]
df<-data.frame(as.vector(wolbachia),environment,site)
colnames(df)<-c('wolbachia','environment','site')
df$wolbachia<-df$wolbachia/22533

model <- betareg(wolbachia ~ environment, data = df)
predicted_data <- data.frame(environment = seq(min(df$environment), max(df$environment), length.out = 100))
predicted_data$wolbachia = predict(model, newdata = predicted_data, type = "response")


p <- ggplot(df, aes(x = environment, y = wolbachia)) +
  geom_point(color='black',size = 3, shape=16) +scale_fill_manual(values=c('black','black','black'))+# Add points
  geom_line(data = predicted_data, aes(x = environment, y = wolbachia), color = "black", size = 1,linetype="dashed") +labs(x = "Environment", y = "Wolbachia Relative Abundance")+
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 22),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




# Display the plot
print(p)


#Test for spatial autocorrelation using Moran's I


set.seed(10)

latlong<-jitter(as.matrix(lonlat),amount = 0.001)
nb <- chooseCN(coordinates(latlong), type = 1, plot.nb = FALSE)
distnb <- nbdists(nb, latlong)
lw <- nb2listw(nb, style = 'W', zero.policy = FALSE)

moran.test(as.vector(wolbachia), lw,alternative = 'greater', randomisation = TRUE)
ff<-moran.plot(as.vector(wolbachia), lw,labels = FALSE,xlim=c(-5000,25000))

if (require(ggplot2, quietly=TRUE)) {
  xname <- attr(ff, "xname")
  ggplot(ff, aes(x=x, y=wx)) + geom_point(shape=1) + xlim(-100,19000)+
    geom_smooth(formula=y ~ x, method="lm") +
    geom_hline(yintercept=mean(ff$wx), lty=2) +
    geom_vline(xintercept=mean(ff$x), lty=2) + theme_minimal() +
    geom_point(data=ff[ff$is_inf,], aes(x=x, y=wx), shape=9) +
    geom_text(data=ff[ff$is_inf,], aes(x=x, y=wx, label='', vjust=1.5)) +theme_bw()+theme(panel.grid = element_blank())+
    xlab('Wolbachia') + ylab(paste0("Spatially lagged Wolbachia"))+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 13))
}

#cd<-as.matrix(vegdist(sign(t(otu_table(data_rarified))),metric='jaccard',diag=TRUE,upper=TRUE))

cd<-matrix(nrow=length(df$wolbachia),ncol=length(df$wolbachia))
for (k in 1:length(df$wolbachia)){
  for (j in 1:length(df$wolbachia)){
    cd[k,j]<-abs(df$wolbachia[k]-df$wolbachia[j])
  }
}

min_soil_temp<-sedf[,1]
max_soil_temp<-sedf[,2]
soil_moisture<-sedf[,3]
elevation<-sedf[,4]


#Spatial distance matrix
sd<-distm(data.frame(sample_data(data_rarified)$Long,sample_data(data_rarified)$Lat), fun=distHaversine)
dbmem = dbmem(as.dist(sd))


Y<-cbind(dbmem,environment)
rda1<-dbrda(cd~environment+MEM1+MEM2+MEM3+MEM4+MEM5+MEM6,data=Y)
vif.cca(rda1)
anova.cca(rda1)
vp<-varpart(cd,dbmem,environment)
par(mar=c(3,3,3,3))
plot(vp, bg = c('tan', 'blue'), Xnames=c('',''),cex=1.5)


colorme<-rep('black',length(as.vector(wolbachia)))
colorme[which(as.vector(wolbachia)>113)]<-'red4'
colorme[which(as.vector(wolbachia)>1130)]<-'red'
llat<-sample_data(data_rarified)$Lat
llong<-sample_data(data_rarified)$Long

lld<-data.frame(llat,llong)

shapefile_path <- "GRSM_BOUNDARY_POLYGON.shp"
shape <- st_read(shapefile_path)
shape <- st_transform(shape, crs = st_crs("+proj=longlat +datum=WGS84"))

# Simplify the geometry with a small tolerance
shape <- st_simplify(shape, dTolerance = 0.001)  # Adjust dTolerance as needed
# Filter for only the largest polygon or relevant geometries
shape <- shape[st_area(shape) == max(st_area(shape)), ]

latlong2<-jitter(as.matrix(lld),amount = 0.015)

ggplot(shape)+geom_sf()+theme_classic2()+geom_point(data = latlong2, mapping = aes(x = llong, y = llat), colour = colorme,size=1.5)+ylab('Latitude')+xlab('Longitude')+theme_bw()


