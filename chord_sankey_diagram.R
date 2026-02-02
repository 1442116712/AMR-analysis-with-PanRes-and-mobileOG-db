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
  
  "B. thetaiotaomicron NLAE-zl-G288" = "#1F77B4",  
  "S. gallolyticus ATCC 700065" = "#FF7F0E",       
  "B. ovatus NLAE-zl-C11" = "#2CA02C",             
  "B. ovatus NLAE-zl-C34" = "#D62728",             
  "B. xylanisolvens NLAE-zl-C182" = "#9467BD",     
  "B. xylanisolvens NLAE-zl-G421" = "#8C564B",     
  "D. longicatena AGR2136" = "#E377C2",            
  "Fibrobacter sp. UWOV1" = "#7F7F7F",             
  "Fibrobacter sp. UWB12" = "#BCBD22",             
  "Cellvibrio sp. BR" = "#17BECF",                 
  "B. thetaiotaomicron NLAE-zl-C523" = "#AEC7E8",  
  "B. xylanisolvens NLAE-zl-C29" = "#FFBB78",      
  "Citrobacter sp. NLAE-zl-C269" = "#98DF8A",      
  "E. coli PA-3" = "#FF9896",                     
  "E. gallinarum SKF1" = "#C5B0D5",               
  "Enterobacter sp. KPR-6" = "#C49C94",          
  "L. lactis DPC6856" = "#F7B6D2",               
  "S. bovis ATCC 33317" = "#C7C7C7",              
  "S. sonnei NLAE-zl-G496" = "#DBDB8D",          
  "O. ruminantium GH1" = "#9EDAE5",              
  "B. nealsonii AAU1" = "#393B79",                 
  
  "other_category" = "gray"  
)


p <- ggplot(sankey_data,
       aes(axis1 = Mobile_ARG, axis2 = Class, axis3 = Genera, y = Count)) +

  geom_alluvium(aes(fill = Mobile_ARG), width = 1/24, alpha = 0.65) +
  
  geom_stratum(
    width = 1/24,
    aes(fill = after_stat(stratum)),  
    color = NA
  ) +
  scale_fill_manual(
    values = group_colors,
    guide = "none"
  ) +
  scale_x_discrete(limits = c("Mobile_ARG", "Class", "Genera"), expansion(add = c(1, 1))) + 
  coord_cartesian(clip = "off") +
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
  "plot.png",       
  plot = p,         
  width = 8,        
  height = 6,    
  dpi = 300,      
  units = "in"      
)

# chord diagram ####
install.packages("circlize")
library(circlize)

chord_data <- data %>%
  group_by(Mobile_ARG, MGE) %>%
  summarise(Count = n()) %>%
  ungroup()

png("chord_plot.png", 
    width = 8,         
    height = 7,        
    units = "in",      
    res = 300)         
circos.par(
  canvas.xlim = c(-1.2, 1.2),  
  canvas.ylim = c(-1.2, 1.2), 
  start.degree = 90
)

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


# Set global rotation parameters
circos.par(start.degree = 0) 

chordDiagram(
  data[, c("Mobile_ARG", "MGE")],
  grid.col = group_colors,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  link.arr.type = "big.arrow",
  transparency = 0.25,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.02), 
  preAllocateTracks = list(track.height = 0.01)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xcenter = get.cell.meta.data("xcenter")
    ycenter = get.cell.meta.data("ycenter")
    circos.text(
      xcenter, 
      ycenter, 
      sector.index,
      facing = "clockwise",  
      niceFacing = TRUE,  
      adj = c(0, 0.5),
      cex = 0.8,
      col = "black"
    )
  },
  bg.border = NA
)
circos.clear()
dev.off()


