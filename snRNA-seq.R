library(Seurat)
library(tidyverse)
library(mascarade)
library(tidydr)

load('00.rda/sn_s16_20260113.rda')

anno_color_d = c("Cor-Glu" = '#3a949e',"Subcor-Glu" = '#2a78c0',"GABA-N" = '#6f9b68',
    "AST" = '#f48fb1',"OLG" = '#b483c1',"OPC" = '#8f5eb5',"MG" = '#ca5655',
    "EC" = '#e67e22',"FB" = '#f2c94c',
    "T cell" = '#cc908d')

Region_color = c('LP1' = '#f6b2b2', 'LP2' = '#e78e8e', 'LP3' = '#d76f6f', 'LP4' = '#c94f4f',
                 'LP5' = '#e3a6b8', 'LP6' = '#d78aa0', 'LP7' = '#c56e88', 'LP8' = '#b1546f',
                 'CP1' = '#6fa6c2', 'CP2' = '#4b8db8', 'CP3' = '#2f6e9e', 'CP4' = '#1d4f77',
                 'CP5' = '#a8d8e9', 'CP6' = '#6db5c8', 'CP7' = '#3ca2b9', 'CP8' = '#1e7a92')


anno_umap_d_label = DimPlot(sn_s16, reduction = "umap", group.by = "anno",cols = anno_color_d,raster =F,pt.size = 0.1,label =T)+coord_equal()
ggsave('anno_umap_d_label.png',anno_umap_d_label,width = 8,height = 6,dpi = 300)

anno_umap_region = DimPlot(sn_s16, reduction = "umap", group.by = "Region",cols = Region_color,raster =F,pt.size = 0.05,label =F,alpha = 0.3)+coord_equal()
ggsave('anno_umap_region.png',anno_umap_region,width = 6,height = 6,dpi = 300)

cell_qc = VlnPlot(sn_s16,features = c('nCount_RNA', 'nFeature_RNA', 'percent_mt'), pt.size = 0,
                  group.by = 'anno',fill.by = 'ident',cols = anno_color_d, stack =T, flip =T)+
                  NoLegend()
ggsave('anno_qc.png',width=8,height = 4)
ggsave('anno_qc.pdf',width=8,height = 4)

mask_table = generateMask(dims=anno_umap_d_label@data[,1:2],cluster=anno_umap_d_label@data$anno,minDensity=1.5,smoothSigma=0.05)

mask_table_f = filter(mask_table,part %in% c('AST#1','Cor-Glu#3','EC#1','FB#1','GABA-N#2','MG#1','OLG#1','OPC#1','Subcor-Glu#1','T cell#1'))

anno_plot_data = anno_umap_d_label@data
anno_plot_label <- anno_plot_data %>%
    group_by(anno) %>%
    summarise(umap_1 = mean(umap_1),umap_2 = mean(umap_2))

anno_plot_label[8,'umap_1'] = -7
anno_plot_label[8,'umap_2'] = 3
anno_plot_label[9,'umap_1'] = -8
anno_plot_label[9,'umap_2'] = 1

anno_umap_fullname = ggplot()+
  geom_point(data=anno_umap_d_label@data,aes(x=umap_1,y =umap_2,color = anno),size = 0.2,shape = 16, stroke = 0,show.legend =F) +
  geom_path(data=mask_table_f,aes(x=umap_1,y =umap_2,group=cluster),linewidth=0.3,linetype=3,show.legend=F)+
  geom_text(data = anno_plot_label,aes(x= umap_1,y = umap_2,label = label),vjust = 0.5)+
  scale_color_manual(values = anno_color_d) +coord_equal()+
  theme_dr()+theme(axis.title = element_text(size = 8, hjust = 0.05),panel.grid = element_blank(),
                 legend.text = element_text(size = 12, face = "plain",color = "black"), 
                 plot.title = element_blank())
ggsave('anno_umap_fullname.png',anno_umap_fullname,width = 6,height = 6,dpi = 300)
ggsave('anno_umap_fullname.pdf',anno_umap_fullname,width = 6,height = 6)

marker_gene = c('CAMK2A','SLC17A7','SATB2','SLC17A6','TCF7L2','GAD1','GAD2','ALDH1L1','GFAP','MOG','MBP','COL9A1','TNR','ITGAM','C3','DOCK2','FLT1','RGS5','COL1A2','DCN','CD3D','CD3E','CD3G')
anno_marker_Dot = DotPlot(sn_s16,features = marker_gene,group.by = "anno",cols = 'RdGy')+theme(axis.text.x = element_text(angle = 90))
ggsave('anno_marker_Dot.png',anno_marker_Dot,width=8,height =4)
ggsave('anno_marker_Dot.pdf',anno_marker_Dot,width=8,height =4)

Region_proportions <- sn_s16@meta.data %>%
  group_by(Region, anno) %>%
  summarise(count = n()) %>% 
  mutate(cell_per = count / sum(count) * 100) %>%
  ungroup()
write.csv(Region_proportions,'Region_proportions.csv')

Annotation_proportions <- sn_s16@meta.data %>%
  group_by(anno, Region) %>%
  summarise(count = n()) %>% 
  mutate(cell_per = count / sum(count) * 100) %>%
  ungroup()
write.csv(Annotation_proportions,'Annotation_proportions.csv')