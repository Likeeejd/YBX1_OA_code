#### YBX1 code ######

library(Seurat)
library(tidyverse)
library(magrittr)
library(ggsci)
library(RColorBrewer)
library(ggpubr)
library(data.table)
library(patchwork)
library(readxl)
set.seed(12345)


### Figure 1B
visCluster(object = cm,
           ncol = 1,
           #add.mline = FALSE,
           #ms.col = c("green","orange","red"),
           add_new_line = F,
           add.line = F,
           textbox.pos = c(0.5,0.8),
           textbox.size = 8,
           textbar.pos = c(0.5,0.5),
           mline.size = 1.5,
           plot.type = "line")

####  Figure 1C
ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description),fill = Cluster)) +
  scale_fill_manual(values =c('#a75eb6')) +  
  geom_bar(stat = "identity", width = 0.5, alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, 12)) + 
  labs(x = "-log10(pvalue)", y = "", title = " ") +
  geom_text(size=5, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  
  #### 基因名 是否显示 
  # geom_text(size=4, aes(x = 0.05, label = geneID), hjust = 0, vjust = 2.5, 
  #           #color=rep(c('#6bb9d2', '#d55640'),
  #            color=rep(c('black'),
  #                     length.out = nrow(dt))) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme +
  NoLegend()

### Figure 1D
ggplot(BRCA_Match_DEG, aes(x = `log2FoldChange`, y=log10FDR)) +
  geom_point(aes(fill = DEG),alpha= 0.85 , size=5,shape = 21 )+
  scale_fill_manual(values = c("steelblue", "#E5E5E5", "#DE4A4F"))+
  xlim(c(-4, 4)) + 
  geom_vline(xintercept=c(-0.5,0.5),lty=2,col="black",lwd=0.5) + 
  geom_hline(yintercept = -log10(0.05), lty=2,col="black",lwd=0.5) +  
  
  labs(x="log2FoldChange", y="-log10pval") + 
  ### 调整主题风格
  theme(plot.title = element_text(hjust = 0.5,size = 18,
  ),
  legend.position="bottom",
  legend.title = element_blank(),  # legend 的title 
  axis.text.x = element_text(size=18, hjust= 0.5, color="black"),    
  axis.text.y = element_text(size=18, color="black"),
  axis.title.y = element_text(size=18,colour = 'black',hjust = 0.5), 
  axis.title.x = element_text(size=18,colour = 'black',hjust = 0.5), 
  panel.grid = ggplot2::element_blank(), 
  legend.text = element_text(size = 18, family = "Arial"))


### Figure 1E 
Heatmap(selected_expression_data, 
        name = "Expression", 
        col = structure(my_palette, breaks = my_breaks),
        column_title = NULL,      # column 的title 设置  设置成NULL 头上方的1 2就会消失 
        
        # row_title = "Mitochondrial fission",
        
        width = unit(4, "cm"),
        height = unit(4,"cm"),
        
        column_order = colnames(selected_expression_data),
        column_gap = unit(0, "mm"),     #分隔间距  
        #border_gp = gpar(lwd = 2),    #边框粗细
        cluster_rows = T,
        cluster_columns = F,
        show_row_dend = F,
        show_column_dend = F,
        show_row_names = T,
        row_names_side = c("left"),  # 调整行名位置 ， 基因名位于那一侧
        row_names_rot = 0,    # 调整行名角度
        column_title_side = c("top"),
        show_column_names = F,      # 不显示列名 图下方
        column_names_centered = F,
        top_annotation = HeatmapAnnotation(df = data.frame(Group = col_anno),
                                           col = list(Group = col_anno_colors),
                                           gp = gpar(col = "white",lwd=1),
                                           show_annotation_name = FALSE,  # 不展示图例名  top
                                           show_legend = F,  # 不展示该部分 图例 
                                           simple_anno_size = unit(0.3, "cm")),
        # column_km = 2,
        
        show_heatmap_legend = T ,
        
        
        heatmap_legend_param = list(   
          at = c(-1, 0, 1),
          grid_width = unit(3, "mm"),
          # labels = c("low", "zero", "high"),
          title = " ",
          legend_height = unit(1, "cm"),       ###### 此处可以修改 legend的 长度
          legend_direction = "vertical"  # 设置图例为垂直方向，帮助控制放置在上方
        ))

### FigureS1 
pca.plot<-scatterplot3d(pca_mat[1:3], #数据
                        angle = 40, #角度
                        axis=T,
                        label.tick.marks=TRUE,
                        xlab="",
                        ylab="",
                        zlab="", # XYZ 轴名称
                        col.grid="black", 
                        color = c("black","black")[as.factor(pca_mat$group)],#创建颜色 设置圆圈的边框颜色
                        pch = c(21,21)[as.factor(pca_mat$group)], #形状 空心圆圈，
                        bg=c("#91bdd1","#ba93cf")[as.factor(pca_mat$group)],  # 设置圆圈的填充色 
                        box = F, 
                        grid=F, 
                        highlight.3d = F,
                        #  ellipsoid = TRUE,
                        lty.axis = 3,   #坐标轴线的类型
                        font.axis = 1,   # 1 正常 ，2 加粗
                        lty.grid = 6,
                        main = "Proteome",
                        cex.symbols = 1.8
)


### Figure S5
DimPlot(srsc,reduction = "tsne")+ 
  scale_fill_manual(values = allcolour) +
  scale_color_manual(values = allcolour) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(size = 10),  
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=10), 
        legend.key.size=unit(0.52,'cm') ) + 
  ggtitle(" ") +
  # theme_dr(xlength = 0.35, ylength = 0.35,
  #          arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  #) +
  theme_dr(xlength = 0.35, ylength = 0.35,  # 轴长度
           arrow = grid::arrow(length = unit(0.1, "inches"), 
                               type = "closed") 
  ) +
  theme(panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=3.5)))


jjDotPlot(object = srsc,
          ytree = F,
          id = "Celltype",  # seurat_clusters 
          x.text.angle = 90, # x 轴角度 
          x.text.vjust = 0.5,
          x.text.hjust = 1, # 对齐 
          gene = c("YBX1","CCL2","PTGS2","NOS2"))




