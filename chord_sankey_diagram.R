install.packages("readxl")
install.packages("tidyverse")
install.packages("devtools")
devtools::install_github("r-lib/conflicted")

library(readxl)
library(tidyverse)

data <- read_excel("C:/Users/CFL/Desktop/1.xlsx", sheet = "Sheet1")

# sankey diagram ####
install.packages("ggalluvial")
library(ggalluvial)

sankey_data <- data %>%
  group_by(Mobile_ARG, Class, Genera) %>%
  summarise(Count = n()) %>%
  ungroup()

group_colors <- c(
  "aminoglycoside" = "#E64B35",
  "beta-lactam" = "#4DBBD5",
  "elfamycins" = "#32CD32", 
  "glycopeptide" = "#3C5488", 
  "lincosamide" = "#FF69B4",  
  "macrolide" = "#8491B4",
  "multidrug" = "#91D1C2",
  "phenicol" = "#FFD700", 
  "streptothricin" = "#7E6148",
  "tetracycline" = "#9370DB",
  
  "Integron" = "#A52A2A",
  "IS/Transposon" = "#1E90FF",
  "Gene Island" = "#808000",
  "Plasmid" = "#228B22",
  
  "Bacteroidia" = "#1F77B4",
  "Gammaproteobacteria" = "#D62728",
  "Clostridia" = "#FF7F0E",
  "Bacilli" = "#BCBD22", 
  "Fibrobacteria" = "#17BECF",
  
  "Cellvibrio" = "#2CA02C",
  "Fibrobacter" = "#D62728", 
  "Dorea" = "#9467BD", 
  "Bacteroides" = "#8C564B", 
  "Streptococcus" = "#E377C2", 
  "Enterococcus" = "#98DF8A", 
  "Oscillibacter" = "#FFBB78", 
  "Citrobacter" = "#17BECF", 
  "Enterobacter" = "#FF9896",
  "Escherichia" = "#393B79", 
  "Shigella" = "#C7C7C7", 
  "Bacillus" = "#F7B6D2",
  "Lactococcus" = "#9EDAE5",
  
  "B. thetaiotaomicron NLAE-zl-G288" = "#1F77B4",  # 蓝色
  "S. gallolyticus ATCC 700065" = "#FF7F0E",       # 橙色
  "B. ovatus NLAE-zl-C11" = "#2CA02C",             # 绿色
  "B. ovatus NLAE-zl-C34" = "#D62728",             # 红色
  "B. xylanisolvens NLAE-zl-C182" = "#9467BD",     # 紫色
  "B. xylanisolvens NLAE-zl-G421" = "#8C564B",     # 棕色
  "D. longicatena AGR2136" = "#E377C2",            # 粉色
  "Fibrobacter sp. UWOV1" = "#7F7F7F",             # 灰色
  "Fibrobacter sp. UWB12" = "#BCBD22",             # 橄榄绿
  "Cellvibrio sp. BR" = "#17BECF",                 # 青色
  "B. thetaiotaomicron NLAE-zl-C523" = "#AEC7E8",  # 浅蓝
  "B. xylanisolvens NLAE-zl-C29" = "#FFBB78",      # 浅橙
  "Citrobacter sp. NLAE-zl-C269" = "#98DF8A",      # 浅绿
  "E. coli PA-3" = "#FF9896",                      # 浅红
  "E. gallinarum SKF1" = "#C5B0D5",               # 浅紫
  "Enterobacter sp. KPR-6" = "#C49C94",            # 浅棕
  "L. lactis DPC6856" = "#F7B6D2",                # 浅粉
  "S. bovis ATCC 33317" = "#C7C7C7",              # 浅灰
  "S. sonnei NLAE-zl-G496" = "#DBDB8D",           # 浅橄榄
  "O. ruminantium GH1" = "#9EDAE5",               # 浅青
  "B. nealsonii AAU1" = "#393B79",                 # 深蓝
  
  "other_category" = "gray"  # 备用颜色
)


p <- ggplot(sankey_data,
       aes(axis1 = Mobile_ARG, axis2 = Class, axis3 = Genera, y = Count)) +
  # 流动层：按Mobile_ARG填充颜色
  geom_alluvium(aes(fill = Mobile_ARG), width = 1/24, alpha = 0.65) +
  # 节点层：按group_colors填充颜色
  geom_stratum(
    width = 1/24,
    aes(fill = after_stat(stratum)),  # 使用节点名称映射颜色
    color = NA
  ) +
  scale_fill_manual(
    values = group_colors,
    guide = "none"
  ) +
  scale_x_discrete(limits = c("Mobile_ARG", "Class", "Genera"), expansion(add = c(1, 1))) + # 安排列的顺序并左右各扩展1个单位空间
  coord_cartesian(clip = "off") + # 允许绘图溢出面板
  theme_minimal() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
  )
    
print(p)

ggsave(
  "plot.png",       # 文件名
  plot = p,         # 图形对象
  width = 8,        # 宽度（英寸）
  height = 6,       # 高度（英寸）
  dpi = 300,        # 分辨率（默认 300，可调高至 600+）
  units = "in"      # 单位（英寸）
)

# chord diagram ####
install.packages("circlize")
library(circlize)

chord_data <- data %>%
  group_by(Mobile_ARG, MGE) %>%
  summarise(Count = n()) %>%
  ungroup()

png("chord_plot.png", 
    width = 8,         # 宽度（英寸）
    height = 7,        # 高度（英寸）
    units = "in",      # 单位（英寸）
    res = 300)         # 分辨率（DPI）
circos.par(
  canvas.xlim = c(-1.2, 1.2),  # 默认c(-1, 1)
  canvas.ylim = c(-1.2, 1.2),   # 扩展边界
  start.degree = 90
)

# 自定义颜色（可选）
group_colors <- c(
  "aminoglycoside" = "#E64B35",
  "beta-lactam" = "#4DBBD5",
  "elfamycins" = "#32CD32", 
  "glycopeptide" = "#3C5488", 
  "lincosamide" = "#FF69B4",  
  "macrolide" = "#8491B4",
  "multidrug" = "#91D1C2",
  "phenicol" = "#FFD700", 
  "streptothricin" = "#7E6148",
  "tetracycline" = "#9370DB",
  "Integron" = "#A52A2A",
  "IS/Transposon" = "#1E90FF",
  "Gene Island" = "#808000",
  "Plasmid" = "#228B22")


# 设置全局旋转参数
circos.par(start.degree = 0)  # 将整个图形旋转90度

# 绘制和弦图（标签自动跟随旋转）
chordDiagram(
  data[, c("Mobile_ARG", "MGE")],
  grid.col = group_colors,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  link.arr.type = "big.arrow",
  transparency = 0.25,
  annotationTrack = "grid",  # 显示名称和网格
  annotationTrackHeight = c(0.05, 0.02), # 调整标签和网格的高度
  preAllocateTracks = list(track.height = 0.01) # 为自定义标签预留空间
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xcenter = get.cell.meta.data("xcenter")
    ycenter = get.cell.meta.data("ycenter")
    # 方式1：使用基础文本（绝对垂直）
    circos.text(
      xcenter, 
      ycenter, 
      sector.index,
      facing = "clockwise",  # 强制向下垂直
      niceFacing = TRUE,   # 禁用自动调整
      adj = c(0, 0.5),
      cex = 0.8,
      col = "black"
    )
  },
  bg.border = NA
)
circos.clear()
dev.off()

