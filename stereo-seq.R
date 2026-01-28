library(Seurat)
library(pheatmap)
library(viridis)
library(tidyverse)

# lesioned concentric zones
Y01405H1L4@meta.data <- Y01405H1L4@meta.data %>%
  mutate(ref_x = ifelse(x > 43000, 89000, -3000),
    ref_y = 52000,
    distance = sqrt((x - ref_x)^2 + (y - ref_y)^2),
    CZ_id = case_when(
      distance <= 20000 ~ "0",
      distance > 20000 & distance <= 25000 ~ "1",
      distance > 25000 & distance <= 30000 ~ "2",
      distance > 30000 & distance <= 35000 ~ "3",
      distance > 35000 & distance <= 40000 ~ "4",
      distance > 40000 & distance <= 45000 ~ "5",
      distance > 45000 ~ "6",
      TRUE ~ NA_character_
    ),
    side = ifelse(x > 43000, 'Lesioned', 'Contralateral'),
    CZ_part = ifelse(x > 43000, paste0('LCZ', CZ_id), paste0('CCZ', CZ_id))
  )


# cell density
cell_density_data = FetchData(object=Y01405H1L4, vars=c("spatialdz_1","spatialdz_2","orig.ident","anno"))
library(ggpointdensity)
p <- ggplot(data = dat, mapping = aes(x = spatialdz_1, y = spatialdz_2)) +
  geom_pointdensity(method='neighbors',size = 0.01) +
  scale_color_viridis_c(option="inferno") + 
  theme_classic(base_size = 15) 
ggsave('pointdensity.png',p,width = 13,height = 12,dpi = 300)

Neuron_hex_50 = ggplot(filter(cell_density_data, anno %in% c('Cor-Glu','Subcor-Glu','GABA-N')), aes(spatialdz_1, spatialdz_2)) +
  stat_bin_hex(bins = 50)+
  scale_fill_gradientn(
    colours = viridis(7, option = "G", direction = -1),
    name = "Density",
    trans = "log10",
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  ) +
  theme_void() +
  coord_equal()
ggsave('Neuron_hex_50.png',Neuron_hex_50,width = 6.5,height = 6,dpi = 300)


# bin heatmap

Y01405H1L4@meta.data = Y01405H1L4@meta.data %>%
  group_by(side) %>%
  mutate(rank = ceiling((distance - min(distance)) / (max(distance) - min(distance)) * 200)) %>%
  ungroup() %>%
  data.frame() %>%
  column_to_rownames('X_index')
Y01405H1L4$dis_group = paste(Y01405H1L4$side,Y01405H1L4$rank,sep = '_')
lesioned_bin_order = c(paste('Lesioned',1:200,sep = '_'))
filter_bin = names(table(Y01405H1L4$dis_group))[table(Y01405H1L4$dis_group) > 1000]
filter_bin = lesioned_bin_order[lesioned_bin_order %in% filter_bin]

filter_bin = names(table(Y01405H1L4$dis_group))[table(Y01405H1L4$dis_group) > 1000]
lesioned_filter_data = subset(Y01405H1L4,dis_group %in% filter_bin)

CZ_marker_gene = CZ_marker$X %>% unique()

CZ_marker_gene_Dot = DotPlot(lesioned_filter_data,features = CZ_marker_gene, group.by = "dis_group")
CZ_marker_gene_data = pivot_wider(CZ_marker_gene_Dot@data[,c(3:5)],names_from = features.plot, values_from = avg.exp.scaled)
CZ_marker_gene_data = column_to_rownames(CZ_marker_gene_data,'id')

bin_anno = Y01405H1L4@meta.data[,c('dis_group','circle_part')]  %>%
  filter(circle_part != "Left")%>%unique()
bin_anno = bin_anno %>%
  filter(dis_group %in% filter_bin)
bin_anno = bin_anno[!duplicated(bin_anno$dis_group),]
rownames(bin_anno) = bin_anno$dis_group
bin_anno = select(bin_anno,-'dis_group')

CZ_marker_gene_data = CZ_marker_gene_data[filter_bin,]

png(filename = 'lesioned_LCZ_marker_f_Bin200_order.png',width = 5,height = 6,units = 'in',res = 300)
pheatmap(heatmap_data,cluster_rows = F,cluster_cols = F,border_color = NA,show_colnames=F,
color = viridis(200, option = "G",direction = 1),annotation_col = bin_anno,show_rownames = T,fontsize_row=3)
dev.off()

pdf('lesioned_LCZ_marker_f_Bin200_order.pdf',width = 5,height = 6)
pheatmap(heatmap_data,cluster_rows = F,cluster_cols = F,border_color = NA,show_colnames=F,gaps_col = 162,
color = viridis(200, option = "G",direction = 1),annotation_col = bin_anno)
dev.off()

