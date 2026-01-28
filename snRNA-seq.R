library(Seurat)
library(tidyverse)
library(mascarade)
library(tidydr)

load('00.rda/sn_s16_20260113.rda')

# Fig2c
anno_color_d = c("Cor-Glu" = '#3a949e',"Subcor-Glu" = '#2a78c0',"GABA-N" = '#6f9b68',
    "AST" = '#f48fb1',"OLG" = '#b483c1',"OPC" = '#8f5eb5',"MG" = '#ca5655',
    "EC" = '#e67e22',"FB" = '#f2c94c',
    "T cell" = '#cc908d')

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


# Fig2d
marker_gene = c('CAMK2A','SLC17A7','SATB2','SLC17A6','TCF7L2','GAD1','GAD2','ALDH1L1','GFAP','MOG','MBP','COL9A1','TNR','ITGAM','C3','DOCK2','FLT1','RGS5','COL1A2','DCN','CD3D','CD3E','CD3G')
anno_marker_Dot = DotPlot(sn_s16,features = marker_gene,group.by = "anno",cols = 'RdGy')+theme(axis.text.x = element_text(angle = 90))
ggsave('anno_marker_Dot.png',anno_marker_Dot,width=8,height =4)
ggsave('anno_marker_Dot.pdf',anno_marker_Dot,width=8,height =4)

# Fig2e

Anno_pro <- sn_s16@meta.data %>%
  group_by(anno, Region) %>%
  summarise(count = n()) %>% 
  mutate(cell_per = count / sum(count) * 100) %>%
  ungroup()

Anno_pro$anno = factor(Anno_pro$anno,levels = levels(sn_s16$anno))
Anno_pro$part = substr(Anno_pro$Region, start = 2, stop =3)
Anno_pro$side = substr(Anno_pro$Region, start = 1, stop =1)

cell_Region_plot = ggplot(Anno_pro,aes(x = anno ,y = cell_per,fill = Region ))+
  geom_col(width = 0.6)+
    geom_hline(yintercept = 50,linetype = 2,color = 'white')+
  scale_fill_manual(values = Region_color)+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text = element_text(color = 'black'),
        legend.position = 'bottom')+
  guides(fill=guide_legend(ncol=8))
ggsave('cell_Region_plot.png',cell_Region_plot,width = 6.5,height = 4)
ggsave('cell_Region_plot.pdf',cell_Region_plot,width = 6.5,height = 4)

anno_per_plot = ggplot(Anno_pro,aes(x = part,y = cell_per,fill = side))+
  geom_col(width = 0.6, position = position_dodge(width = 0.8))+
  scale_fill_manual(values = c(C = '#6fa6c2',L = '#e78e8e'))+
  scale_y_continuous(expand = c(0,0))+
  facet_wrap(~anno,scales = 'free',nrow =2)+
  theme_classic()
ggsave('anno_per_plot.png',anno_per_plot,width = 12,height = 4)
ggsave('anno_per_plot.pdf',anno_per_plot,width = 12,height = 4)

# Fig2i

Part_color = c('P1' = '#f2dd9f', 'P2' = '#f1b6b0', 'P3' = '#f4c2d3', 'P4' = '#cebdd6',
                 'P5' = '#a0a1c1', 'P6' = '#8bc0d7', 'P7' = '#84b3b2', 'P8' = '#b3ccb3')

Region_list = list.files(pattern = '_LC_DEG.csv')

plot_data = data.frame(matrix(nrow = 0,ncol = 7))
colnames(plot_data) = c('X','p_val','avg_log2FC','pct.1','pct.2','p_val_adj','celltype')
for(i in Region_list){
  tmp_data = read.csv(i)
  tmp_data$celltype = gsub('_LC_DEG.csv','',i)
  plot_data = rbind(plot_data,tmp_data)
}
plot_data = filter(plot_data, p_val_adj < 0.05, abs(avg_log2FC) > 0)
up_y = ceiling(max(plot_data$avg_log2FC))
down_y = floor(min(plot_data$avg_log2FC))

DE_up = plot_data %>%
  group_by(celltype) %>%
  summarize(count = sum(avg_log2FC > 0 & p_val_adj < 0.05 ))
DE_down = plot_data %>%
  group_by(celltype) %>%
  summarize(count = sum(avg_log2FC < 0 & p_val_adj < 0.05 ))

tmp_plot = ggplot(plot_data)+
  geom_jitter(aes(x=celltype,y=avg_log2FC,color=celltype),width=0.2,show.legend =FALSE)+
  geom_text(data = DE_up,aes(x=celltype,y=up_y,label=count),color="darkred")+
  geom_text(data = DE_down,aes(x=celltype,y=down_y,label=count),color="darkblue")+
  scale_color_manual()+
  labs(title = "DEG between lesioned and contralateral")+
  xlab("")+
  ylab(expression(log[2]("Fold Change")))+
  theme_classic(values = Part_color)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
        axis.text.x = element_text(color = "black",size = 10,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(color = "black",size=8),
        axis.title.y = element_text(color = "black",size = 12))
ggsave('Region_DEGs_jitter.pdf',tmp_plot,width = 6,height= 4)
ggsave('Region_DEGs_jitter.png',tmp_plot,width = 6,height= 4,dpi = 300)

# Fig2j
for(i in REgion_list){
  file_name = paste0('03.KEGG_res_0.5/',i,'_Lesioned_all_DEGs_KEGG_res.csv')
  if (file.exists(file_name)){
    tmp_data =  read.csv(file_name)
    if(nrow(tmp_data) > 0){
      tmp_data$anno = i
      KEGG_data_Regio_all = rbind(KEGG_data_Regio_all,tmp_data)
    }
  }else{
    message(i,' have no KEGG iterms eneich')
  }
}

KEGG_data_all_f10 = KEGG_data_Regio_all %>% 
  filter(category_map == 'Cellular Processes',pvalue < 0.05) %>%
  group_by(anno,category_map ) %>%
  slice_min(pvalue, n = 10, with_ties = FALSE)


KEGG_data_all_f10 = KEGG_data_Regio_all %>% 
  filter(category_map == 'Cellular Processes',pvalue < 0.05) %>%
  group_by(anno) %>%
  slice_min(pvalue, n = 10, with_ties = FALSE)

KEGG_data_all_f10$Description = factor(KEGG_data_all_f10$Description,levels = rev(unique(KEGG_data_all_f10$Description)))
KEGG_data_all_f10_ID = unique(KEGG_data_all_f10$ID)

KEGG_data_all_f10_ppd = KEGG_data_all_f10[,c('anno','subcategory_map','ID','Description')] %>%
  mutate(value = 1)
sankey_data = rbind(data.frame(x = 'Region',node = KEGG_data_all_f10_ppd$anno,
                               next_x = 'ID',next_node = KEGG_data_all_f10_ppd$ID),
                    data.frame(x = 'ID',node = KEGG_data_all_f10_ppd$ID,
                               next_x = NA,next_node = NA)) %>%
  as.tibble()
sankey_data$x = factor(sankey_data$x,levels = c('Region','ID'))
sankey_data$next_x = factor(sankey_data$next_x,levels = c('Region','ID'))
sankey_data$node = factor(sankey_data$node,levels = c(paste0('P',8:1),rev(KEGG_data_all_f10_ID)))
sankey_data$next_node = factor(sankey_data$next_node,levels = rev(KEGG_data_all_f10_ID))

all_KEGG_10_sankey = ggplot(sankey_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.8,width = 0.4) +
  geom_sankey_text(size = 3, color = "black",width = 1) +
  theme_sankey(base_size = 0) +
  scale_fill_manual(values = Region_color,na.value = "grey70")+
  labs(x = NULL) +
  theme(legend.position = "none")
ggsave('all_KEGG_10_CP_sankey.pdf',all_KEGG_10_sankey,width = 5,height = 5)
