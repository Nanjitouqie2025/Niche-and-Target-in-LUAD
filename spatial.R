
devtools::install_github("dmcable/RCTD", build_vignettes = TRUE)
pth <-  "D:/R432/Rbase432/R-4.3.2/library/RCTD-master.zip"
devtools::install_local(pth)

#### single 

if(T){
  library(tidyverse)
  library(Seurat)
  library(SeuratObject)
  library(DropletUtils)
  library(org.Hs.eg.db) 
  library(clusterProfiler)
  library(stringr)
}


sc = Read10X("D:/Spatial/sp_and_luad_ctl/sc/P16T1/")
sc <- CreateSeuratObject(counts = sc, 
                         project = 'LUAD', 
                         min.cells = 3, 
                         min.features =50)

sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

range(sc$nCount_RNA)
range(sc$nFeature_RNA)
range(sc$percent.mt)

sc <- subset(sc, 
             subset = percent.mt < 20)

# p1 <- VlnPlot(sc,pt.size = 0, 
#               features = c("nFeature_RNA")) 
# p2 <- VlnPlot(sc,pt.size = 0, 
#               features = c("nCount_RNA"))
# p3 <- VlnPlot(sc,pt.size = 0, 
#               features = c( "percent.mt")
#              )
# p4 = FeatureScatter(sc,feature1 = "nCount_RNA",pt.size=0.8,raster=T,raster.dpi = c(800, 800),
#                    feature2 = "nFeature_RNA")
# p4
# 
# p5 = FeatureScatter(sc,feature1 = "nCount_RNA",pt.size=0.8,raster=T,raster.dpi = c(800, 800),
#                     feature2 = "percent.mt")
# p5


sc <- NormalizeData(sc, 
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)

sc<- FindVariableFeatures(sc, 
                          selection.method="vst", nfeatures=3000) 

all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)

sc <- RunPCA(sc, 
             features = VariableFeatures(object = sc))

ElbowPlot(sc)

pc.num = 1:16
sc <- RunUMAP(sc, dims=pc.num)
sc <- RunTSNE(sc,  dims = pc.num )

sc <- FindNeighbors(sc, dims = pc.num) %>%
  FindClusters(resolution = 1)

cols=Getcolor::getrdmcol(length(Idents(sc) %>% unique()))

DimPlot(sc,reduction = "tsne",raster = T,label=T,
        group.by = 'seurat_clusters',
        label.box = T,pt.size =2,cols = cols ) 

DimPlot(sc,reduction = "umap",raster = T,label=T,
        group.by = 'seurat_clusters',
        label.box = T,pt.size =2,cols = cols )


mar = c("TRAC","CD3D",#t
        'NKG7',
        'CD14','LYZ',
        "MS4A1","CD79A",#b
        "EPCAM", "KRT8", #epi
        "MZB1", "JCHAIN", #pla
        "MKI67", "TOP2A"#cycle,
        #'RACGAP1'
)


DotPlot(sc,cols = c('lightgrey','black'),  features = mar)+
  coord_flip() 


#### gptcelltype  ####
markers <- FindAllMarkers(sc,
                          logfc.threshold = 0.25,
                          min.pct = 0.1,
                          test.use = "wilcox")    

top20 <- markers %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC) %>%
  ungroup()

library(GPTCelltype)
#GPTCelltype::gptcelltype()
Sys.setenv(OPENAI_API_KEY =
             'sk-')
base_url<-"https://..."


res2 <- gptcelltype(markers,
                    topgenenumber = 20,
                    model = 'gpt-4o',base_url = base_url)


mar = c("TRAC","CD3D",'CD8A',
        'CD14','LYZ',
        "MS4A1","CD79A",#b
        "EPCAM", "KRT8", #epi
        "MZB1", "JCHAIN", #pla
        "MKI67", "TOP2A", #cycle,
        "VWF", "PECAM1", 
        "COL1A1","COL1A2" #fib
        #'RACGAP1'
)


DotPlot(sc,cols = c('lightgrey','black'),  features = mar)+
  coord_flip()+theme_bw()


new.cluster.ids = c( "unknow",              #0
                     "unknow" ,  #1
                     "Epithelial cells", #2
                     "myeloid cells", #3            
                     "T cells",  #4
                     "Epithelial cells", #5
                     "Epithelial cells", #6
                     "CD8+ cells",#7
                     "Myeloid cells"  , #8                  
                     "Cycling Epithelial cells", #9
                     "Myeloid cells", #10
                     "B cells", #11
                     "Cycling Epithelial cells", #12
                     "unknow",
                     "unknow",
                     "T cells", #15
                     "Epithelial cells", #16
                     "unknow",
                     "unknow", #18
                     "unknow",
                     "endothelial",
                     "Fibroblast", #21
                     "unknow", #22
                     "unknow"
                     
)


res2 <- c(res2
          0                            1 
          "- Monocyte"                 "- Monocyte" 
          2                            3 
          "- T cell"                   "- T cell" 
          4                            5 
          "- Macrophage" "- Natural killer (NK) cell" 
          6                            7 
          "- Macrophage"           "- Dendritic cell" 
          8                            9 
          "- B cell"                   "- T cell" 
          10                           11 
          "- T cell"                   "- T cell" 
          12                           13 
          "- Macrophage"                   "- T cell" 
          14                           15 
          "- T cell"                   "- T cell" 
          16                           17 
          "- B cell"             "- Myeloid cell" 
          18                           19 
          "- Natural killer (NK) cell"             "- Myeloid cell" 
          20 
          "- Proliferating cell")

new.cluster.ids <- res2

names(new.cluster.ids) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids)

sc$cell_type <- Idents(sc)

sc2 <- subset(sc,subset=cell_type!='unknow')

mar = c("TRAC","CD3D",'CD8A',
        'CD14','LYZ',
        "MS4A1","CD79A",#b
        "EPCAM", "KRT8", #epi
        "MZB1", "JCHAIN", #pla
        "MKI67", "TOP2A", #cycle,
        "VWF", "PECAM1", 
        "COL1A1","COL1A2" #fib
        #'RACGAP1'
)

DotPlot(sc2,cols = c('lightgrey','black'),  features = mar)+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle=60,hjust=1) )


DimPlot(sc2,reduction = "umap",raster = T,label=T,
        
        label.box = T,pt.size =2,cols = cols )











#### RCTD  ####
pth1 <- "D:/Spatial/sp_and_luad_ctl/P16_T1/"
pth2 <- "D:/Spatial/sp_and_luad_ctl/P16_T1/spatial"

st <- Load10X_Spatial(data.dir =pth1, 
                      filename = paste0("P16_T1-filtered_feature_bc_matrix.",'h5'),  
                      slice =pth2)

ct <- as.matrix(st@assays$Spatial@counts)

coord <- st@images[[1]]@coordinates

coords <- GetTissueCoordinates(st)
head(coords)

coord <- setNames(as.data.frame(coords),c('x','y'))

head(coord)

library(viridis)


nUMI <- colSums(ct)
head(nUMI)

library(spacexr)


obj <- SpatialRNA(coord,ct,nUMI)

barcodes <- colnames(obj@counts)

plot_puck_continuous(
  obj,barcodes,obj@nUMI,
  ylimit = c(0,round(quantile(obj@nUMI,0.9))),
  size=2
  
)

sc$celltype <- Idents(sc)



sc_ct <- as.matrix(sc@assays$RNA@counts)

cell_types <- sc$cell_type

numi <- sc$nCount_RNA

ref <- Reference(sc_ct,cell_types,numi )

rctd <- create.RCTD(obj,ref,max_cores = 1)

rctd <- run.RCTD(rctd,doublet_mode = 'full')

barcodes <- colnames(rctd@spatialRNA@counts)

weight <- rctd@results$weights

norm_weight <- normalize_weights(weight)


idx <- unique(cell_types)

i=1

pf <- function(celltype){
  
  p <- plot_puck_continuous( rctd@spatialRNA,
                             small_point = T,
                             barcodes,norm_weight[,celltype],
                             size=1,alpha=0.8 )
  
  pp <- p +labs(title = celltype)+
    theme( panel.background = element_rect(fill = 'black', color = NA),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.text = element_blank(),
           axis.title = element_blank() 
    )+
    scale_x_reverse()+
    coord_flip()+
    scale_color_gradientn(colors = c('black','lightgrey','tomato'))
  
  pp
  
}



pplot <- lapply(idx, pf )

pplot[[9]]

psave <- ggpubr::ggarrange(
  ncol=4,nrow=2,common.legend = T,
  pplot[[1]],pplot[[2]],pplot[[3]],pplot[[4]],pplot[[5]],
  pplot[[6]],pplot[[7]],pplot[[8]]
)
psave 

ggsave(psave,width = 10,height=6,
       file='D:/Aproj/aMR/Pt_lung1107/Spatial_res/RCTD_deconvo_p16.pdf')


p <- SpatialFeaturePlot(
  st,
  pt.size.factor = 1.5,
  features = "RACGAP1",
  alpha = 1,
  image.alpha = 0,
  ncol = 1 
)+scale_fill_gradientn(colors = c('lightcyan','yellow','red')) 

pdat <- p[[1]]$data

pra <- ggplot(pdat, aes(x = imagerow, y = imagecol, fill = RACGAP1)) +
  geom_point(size = 2,alpha = 1, shape = 21) +
  scale_fill_gradientn(colors = c('purple','white','yellow')) +
  theme_void() +
  coord_flip() +
  scale_x_reverse()+
  labs(title = "RACGAP1 Spatial Expression") +
  guides(color = guide_colorbar(title=NULL, barwidth=0.5, barheight=10)) +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    plot.title = element_text(color = "white", hjust=0.5),
    legend.position = "right"
  )

ggsave(pra,width = 5,height=4,
       file='D:/Aproj/aMR/Pt_lung1107/Spatial_res/sp_ra_Exp1.pdf')


p <- SpatialFeaturePlot(
  st,
  pt.size.factor = 1.5,
  features = c("RACGAP1"),
  alpha = 1,
  image.alpha = 0,
  
) #+scale_fill_gradientn(colors = c('black','grey','yellow')) 

p

ggsave(p,width = 4,height=4,
       file='D:/Aproj/aMR/Pt_lung1107/Spatial_res/sp_ra_Exp1.pdf')


theme(
  plot.background = element_rect(fill = "black", color = NA),   
  panel.background = element_rect(fill = "black", color = NA),  
  legend.background = element_rect(fill = "black", color = NA), 
  legend.text = element_text(color = "white"),                  
  axis.text = element_blank(),
  axis.ticks = element_blank()
)                 




#### RACGAP1 co   loca ####
pth1 <- "D:/Spatial/sp_and_luad_ctl/P16_T1/"
pth2 <- "D:/Spatial/sp_and_luad_ctl/P16_T1/spatial"

st <- Load10X_Spatial(data.dir =pth1, 
                      filename = paste0("P16_T1-filtered_feature_bc_matrix.",'h5'), 
                      slice =pth2)

SpatialFeaturePlot(
  st,pt.size.factor = 2,
  features = cycle,alpha = 1 ) +
  scale_fill_gradientn(colours = nature_colors1 )



easyspot <- function(gene){
  p <- SpatialFeaturePlot(
    st,pt.size.factor = 2,
    features = gene,alpha = 1 )
  pdat <- p[[1]]$data
  nature_colors1 <- c("black", 
                      "#edf8b1", "#ffffd9","tomato")
  ggplot(pdat,aes(imagerow,imagecol,fill=!!sym(gene)))+
    geom_point(shape=21,size=1)+scale_x_reverse()+coord_flip()+
    scale_fill_gradientn(colours = nature_colors1 )+
    theme_void()+
    theme(  panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = NA)  )
  
}

#### cycle 
cycle <- c('RACGAP1','CDK1','CDC20','CCND1','CCNB2','CDCA3','EPCAM','KRT8')


pp <- lapply(cycle, easyspot)

pp

psave <- ggpubr::ggarrange(
  align = c("hv"),
  pp[[1]],pp[[2]],pp[[3]],pp[[4]],
  pp[[5]],pp[[6]],pp[[7]],pp[[8]],ncol=4,nrow=2)

psave
ggsave(psave,width = 12,height = 5,
       file='D:/Aproj/aMR/Pt_lung1107/Spatial_res/RAcoloc_cellcycle_p24.pdf')

#### sac 
sac <- c('RACGAP1','PRC1','KIF23','KIF4A',
         'ANLN','CIT','ECT2','AURKB','PLK1','GMNN')

pp <- lapply(sac, easyspot)

psave <- ggpubr::ggarrange(
  align = c("hv"),
  pp[[1]],pp[[2]],pp[[3]],pp[[4]],
  pp[[5]],pp[[6]],pp[[7]],pp[[8]],pp[[9]],pp[[10]],ncol=5,nrow=2)

psave
ggsave(psave,width = 14,height = 5,
       file='D:/Aproj/aMR/Pt_lung1107/Spatial_res/RAcoloc_sac_p24.pdf')

#### CA 
CA_genes <- c(
  "AURKA",
  "CCNE2",
  "CEP63",
  "CEP152",
  "E2F1",
  "E2F2",
  "LMO4",
  "MDM2",
  "MYCN",
  "NDRG1",
  "NEK2",
  "PIN1",
  "PLK4",
  "SASS6",
  "STIL",
  "TUBG1"
)

pp <- lapply(CA_genes, easyspot)

psave <- ggpubr::ggarrange(
  align = c("hv"),
  pp[[1]],pp[[2]],pp[[3]],pp[[4]],
  pp[[5]],pp[[6]],pp[[7]],pp[[8]],pp[[9]],pp[[10]],
  pp[[11]],pp[[12]],pp[[13]],pp[[14]],pp[[15]],pp[[16]],
  
  ncol=4,nrow=4)

psave
ggsave(psave,width = 10,height = 10,
       file='D:/Aproj/aMR/Pt_lung1107/Spatial_res/RAcoloc_ca_p24.pdf')





#### ST diff genes ####

if(T){
  library(tidyverse)
  library(Seurat)
  library(SeuratObject)
  library(DropletUtils)
  library(org.Hs.eg.db) 
  library(clusterProfiler)
  library(stringr)
}

pth1 <- "D:/Spatial/sp_and_luad_ctl/P16_T1/"
pth2 <- "D:/Spatial/sp_and_luad_ctl/P16_T1/spatial"

st <- Load10X_Spatial(data.dir =pth1, 
                      filename = paste0("P16_T1-filtered_feature_bc_matrix.",'h5'),  #h5文件声明
                      slice =pth2)

if(T){
  st <- SCTransform(st,assay = 'Spatial',verbose = F)
  st <- RunPCA(st, assay = 'SCT',verbose = F)
  st <- FindNeighbors(st,reduction = 'pca',dims=1:30)
  st <- FindClusters(st,verbose = F)
  st <- RunUMAP(st,reduction = 'pca',dim=1:30)
}

num_idents <- length(levels(st))
col <- Getcolor::getrdmcol(9)
col <- RColorBrewer::brewer.pal(num_idents, "Set1")

colors <- c("grey", "lightcyan", "purple", "cyan", "#FF7F00", 
            "lightyellow", "#A65628", "#F781BF", "#999999")


p <- SpatialDimPlot(st, pt.size.factor = 1.2)+
  scale_fill_manual(values=colors)
p

ggsave(p,width = 5,height = 4.5,
       file='D:/Aproj/aMR/Pt_lung1107/Spatial_res/p16_cluster.pdf')


target_cluster <- "2"

cells_to_highlight <- WhichCells(st, idents = target_cluster)

p <- SpatialDimPlot(st, 
                    cells.highlight = cells_to_highlight,
                    cols.highlight = c("purple", "lightcyan"), 
                    pt.size.factor = 1.2)
p


ggsave(p,width = 5,height = 4.5,
       file='D:/Aproj/aMR/Pt_lung1107/Spatial_res/p16_cluster2.pdf')




gene <- c("RACGAP1", "PLK1", "ECT2", "CEACAM5")

target_gene <-  gene[3]


expression_data <- FetchData(st, vars = target_gene)


cells_to_highlight <- rownames(expression_data)[rowSums(expression_data) > 0]

p <- SpatialDimPlot(st, 
                    cells.highlight = cells_to_highlight,
                    cols.highlight = c("purple", "lightcyan"), 
                    pt.size.factor = 1.6)

p <- p + labs(title = "Cells expressing proliferation markers")

p



easyspot <- function(gene){
  p <- SpatialFeaturePlot(
    st,pt.size.factor = 2,
    features = gene,alpha = 1 )
  pdat <- p[[1]]$data
  nature_colors1 <- c("black", 
                      "#edf8b1", "#ffffd9","tomato")
  ggplot(pdat,aes(imagerow,imagecol,fill=!!sym(gene)))+
    geom_point(shape=21,size=1)+scale_x_reverse()+coord_flip()+
    scale_fill_gradientn(colours = nature_colors1 )+
    theme_void()+
    theme(  panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = NA)  )
  
}

gene <- c("RACGAP1", "PLK1", "ECT2", "CEACAM5")

pp <- lapply(gene, easyspot)

pp

psave <- ggpubr::ggarrange(
  align = c("hv"),
  pp[[1]],pp[[2]],pp[[3]],pp[[4]],
  ncol=4,nrow=1)

psave
ggsave(psave,width = 12,height = 2.5,
       file='D:/Aproj/aMR/Pt_lung1107/Spatial_res/RAcoloc_gene.pdf')







marker <- FindAllMarkers(st,
                         logfc.threshold = 0.20,
                         min.pct = 0.05
)


#### 
mar2 <- marker %>% filter(cluster=='2')

mar2$gene <- rownames(mar2) 

mar2$gene <- sub("\\..*", "", mar2$gene)

mar2 <- mar2[!duplicated(mar2$gene),]

rownames(mar2) <- mar2$gene

if(T){
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(tidyverse)
  library(GseaVis)
  library(msigdbr)
}

deff <- mar2 %>% dplyr::select(c('avg_log2FC','p_val')) %>% 
  rename(logFC='avg_log2FC',P.Value='p_val')

deff <- deff %>% arrange(desc(logFC))


if(T){
  reslimma = data.frame(gene =rownames(deff ),logFC=deff$logFC)
  ranklist = function(x){
    x= reslimma
    entrez <- bitr(x$gene,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
    
    entrez <- entrez[!duplicated(entrez$SYMBOL),]
    list = x[x$gene%in%entrez$SYMBOL,]
    if(identical(list$gene,entrez$SYMBOL)){
      print("construct gene list sucessfully")
    }
    list1 = cbind(list,entrez)
    if(identical(list1$gene,list1$SYMBOL)){
      print("construct merge list sucessfully")
    }
    
    list1 = list1[order(list$logFC,decreasing = T),]
    generank = as.vector(as.numeric(list1$logFC))
    names(generank) = list1$ENTREZID
    return(generank)
  }
  genelist = ranklist(reslimma)
  ranklist2 = function(x){
    x= reslimma
    entrez <- bitr(x$gene,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
    
    entrez <- entrez[!duplicated(entrez$SYMBOL),]
    list = x[x$gene%in%entrez$SYMBOL,]
    if(identical(list$gene,entrez$SYMBOL)){
      print("construct gene list sucessfully")
    }
    list1 = cbind(list,entrez)
    if(identical(list1$gene,list1$SYMBOL)){
      print("construct merge list sucessfully")
    }
    
    list1 = list1[order(list$logFC,decreasing = T),]
    generank = as.vector(as.numeric(list1$logFC))
    names(generank) = list1$SYMBOL
    return(generank)
  }
  list3 =  ranklist2(reslimma)
}

m_df = msigdbr(species = 'Homo sapiens' , category = "H") %>%
  dplyr::select(gs_name,gene_symbol) 

colnames(m_df)<-c("term","gene")


h_ges_tumor_region <- GSEA(list3, 
                           TERM2GENE = m_df,
                           minGSSize = 5,
                           maxGSSize = 500,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           verbose = FALSE,
                           eps = 0)

id <- h_ges_tumor_region@result %>% rownames()
print(id)

p <- gseaNb(object = h_ges_tumor_region,
            subPlot = 3,
            base_size = 14,
            addPval=T,
            geneSetID = id[7])

p     

ggsave(p,width=3,height=2.5,file=
         "D:/Aproj/aMR/Pt_lung1107/Spatial_res/p16_gsea_myc.pdf")

p <- gseaNb(object = h_ges_tumor_region,
            subPlot = 3,
            base_size = 14,
            addPval=T,
            geneSetID = id[21])

p     

ggsave(p,width=3,height=2.5,file=
         "D:/Aproj/aMR/Pt_lung1107/Spatial_res/p16_gsea_p53.pdf")



easyspot <- function(gene){
  p <- SpatialFeaturePlot(
    st,pt.size.factor = 2,
    features = gene,alpha = 1 )
  pdat <- p[[1]]$data
  nature_colors1 <- c("black", 
                      "#edf8b1", "#ffffd9","tomato")
  ggplot(pdat,aes(imagerow,imagecol,fill=!!sym(gene)))+
    geom_point(shape=21,size=1)+scale_x_reverse()+coord_flip()+
    scale_fill_gradientn(colours = nature_colors1 )+
    theme_void()+
    theme(  panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = NA)  )
  
}

cycle <- c('RACGAP1','CDK1','CDC20','CCND1','CCNB2','CDCA3','EPCAM','KRT8')

pp <- lapply(cycle, easyspot)

pp


easyspot <- function(gene){
  p <- SpatialFeaturePlot(
    st,pt.size.factor = 2,
    features = gene,alpha = 1 )
  pdat <- p[[1]]$data
  nature_colors1 <- c("black", 
                      "#edf8b1", "#ffffd9","tomato")
  ggplot(pdat,aes(imagerow,imagecol,fill=!!sym(gene)))+
    geom_point(shape=21,size=1)+scale_x_reverse()+coord_flip()+
    scale_fill_gradientn(colours = nature_colors1 )+
    theme_void()+
    theme(  panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = NA)  )
  
}


SpatialFeaturePlot(
  st,pt.size.factor = 2,
  features = 'CXCL12',alpha = 1 ) 

SpatialFeaturePlot(
  st,pt.size.factor = 2,
  features = 'CXCR4',alpha = 1 ) 



easyspot <- function(gene){
  p <- SpatialFeaturePlot(
    st,pt.size.factor = 2,
    features = gene,alpha = 1 )
  pdat <- p[[1]]$data
  nature_colors1 <- c("black", 
                      "#edf8b1", "#ffffd9","tomato")
  ggplot(pdat,aes(imagerow,imagecol,fill=!!sym(gene)))+
    geom_point(shape=21,size=1)+scale_x_reverse()+coord_flip()+
    scale_fill_gradientn(colours = nature_colors1 )+
    theme_void()+
    theme(  panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = NA)  )
  
}


easyspot(c('CXCL12'))
c('NPR1','VEGFA')

easyspot(c('NPR1','VEGFA'))

gs <- c('NPR1','VEGFA')
pp <- lapply(gs, easyspot)
psave <- ggpubr::ggarrange(
  align = c("hv"),
  pp[[1]],pp[[2]],ncol=2)

psave
ggsave(psave,width=5,height=2.5,file=
         "D:/Aproj/aMR/Pt_lung1107/Spatial_res/p16_coloc.pdf")



gs <- c('SPATA2','CD274') #PD1 PDL1
pp <- lapply(gs, easyspot)
psave <- ggpubr::ggarrange(
  align = c("hv"),
  pp[[1]],pp[[2]],ncol=2)

psave
ggsave(psave,width=5,height=2.5,file=
         "D:/Aproj/aMR/Pt_lung1107/Spatial_res/p16_colocpd1.pdf")










