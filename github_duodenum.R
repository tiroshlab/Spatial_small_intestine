# Load libraries ----------------------------------------------------------
library(ggraph)
library(scales)
library(ggdark)
library(ComplexHeatmap)
library(tidyverse)
library(grid)
library(WeightedTreemaps)

# Load custom functions and levels ----------------------------------------
Savegg <- function(filename, plot) {
  ggsave(filename = filename, width = 20, height = 12, units = "cm", scale = 2, plot = plot)
}

Transform_cell_sum <- function(m) {
  count_sum <- rowSums(m)
  m_sum <- m/count_sum
  m_sum
}

`%ni%` <- Negate(`%in%`)

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

Count_cell_type <- function(cell_type = cell_type) {
  cells %>% count({{cell_type}}) %>% print(n=Inf)
}

hotmap <- readRDS('/home/labs/tirosh/alissag/st_glioma/scripts_workspaces/hotmap.rds')

Make_qpan <- function (export = "pg_cluster") {
  qpan <- cells %>% 
    mutate(image_name = paste0(sample, "_v2")) %>% 
    select(c(centroid_x, 
             centroid_y, 
             image_name, 
             !!sym(export)
    )
    )
  
  tib <- tibble(1:length(qpan %>% pull(!!sym(export)))) # create tibble to insert loop results
  
  for (i in unique(qpan %>% pull(!!sym(export)))) {
    loop <- qpan %>% mutate(i = if_else(!!sym(export) == i, true = 1, false = 0))
    tib <- cbind(tib, loop[[5]])
  }
  
  tib <- tib[,-1] # remove unnecessary column
  colnames(tib) <- unique(qpan %>% pull(!!sym(export))) %>% as.character()
  qpan <- cbind(qpan, tib)
  
  write_csv(qpan[,-4], file = paste0("qpan_all.tsv"))
  return(qpan)
}

Run_DEPs_hm <- function(matrix = rmat, cluster = "pg_cluster") {
  
  DEPs_hm <- tibble(dim(matrix)[1]:dim(matrix)[1]) # Create df in which the DEG for each cluster will be inserted
  
  for (i in sort(unique(cells %>% pull( !!sym(cluster) )))) {    # i iterates all clusters
    top_DEP_value <- (rowMeans(matrix[,which(i == cells %>% pull( !!sym(cluster) )) ]) - rowMeans(matrix[,which(i != cells %>% pull( !!sym(cluster) ))]))
    DEPs_hm <- cbind(DEPs_hm, top_DEP_value)   # select top gene names
  }
  
  DEPs_hm <- DEPs_hm[,-1] # Delete redundant first column
  colnames(DEPs_hm) <- rep( x = paste(sort(unique(cells %>% pull( !!sym(cluster) )))), each = 1)
  return(DEPs_hm)
  
}

Run_DEPs <- function(matrix = rmat, cell_table = cells, cluster = "pg_cluster") {
  DEPs <- tibble(1:dim(matrix)[1]) # Create df in which the DEG for each cluster will be inserted
  
  for (i in sort(unique(cell_table %>% pull( !!sym(cluster) )))) {    # i iterates all clusters
    top_DEP_ix <- sort( (rowMeans(matrix[,which(i == cell_table %>% pull( !!sym(cluster) ) ) ]) - rowMeans(matrix[,which(i != cell_table %>% pull( !!sym(cluster) ))]) ), decreasing = T, index.return = T )$ix[1:dim(matrix)[1]]   # sort gene ix
    top_DEP_value <- sort( (rowMeans(matrix[,which(i == cell_table %>% pull( !!sym(cluster) )) ]) - rowMeans(matrix[,which(i != cell_table %>% pull( !!sym(cluster) ))]) ), decreasing = T)
    DEPs <- cbind(DEPs, rownames(matrix)[top_DEP_ix[1:dim(matrix)[1]]], top_DEP_value)   # select top gene names
  }
  
  DEPs <- DEPs[,-1] # Delete redundant first column
  colnames(DEPs) <- rep( x = paste(sort(unique(cell_table %>% pull( !!sym(cluster) )))), each = 2)
  
  return(print.data.frame(DEPs))
}

Run_DEPs <- function(sample) {
  amat_sample <- amat[, cells %>% 
                        filter(sample %in% c({{sample}})) %>% 
                        filter(cell_type %ni% c("excluded")) %>%
                        pull(cell_name)
  ]
  
  cells_sample <- cells %>% 
    filter(sample %in% c({{sample}})) %>%  
    filter(cell_type %ni% c("excluded"))
  
  DEPs_hm_sample <-
    tibble(dim(amat_sample)[1]:dim(amat_sample)[1]) # Create df in which the DEG for each cluster will be inserted
  
  for (i in sort(unique(cells_sample$cell_type))) {
    # i iterates all clusters
    top_DEP_value <-
      (rowMeans(amat_sample[, which(i == cells_sample$cell_type)]) - rowMeans(amat_sample[, which(i != cells_sample$cell_type)]))
    DEPs_hm_sample <-
      cbind(DEPs_hm_sample, top_DEP_value)   # select top gene names
  }
  
  DEPs_hm_sample <-
    DEPs_hm_sample[, -1] # Delete redundant first column
  colnames(DEPs_hm_sample) <-
    paste0(sort(unique(cells_sample$cell_type)), "_", {
      {
        sample
      }
    })
  return(DEPs_hm_sample)
}

Min_dist_to_nth_cell_type_distribution <- function(cell_type = "epi", 
                                                   ROI = "P13_Serrated_villi_05", 
                                                   n = 1,
                                                   cell_type_column = "cell_type_2") {
  # Create distance matrix
  cells_clean <- cells %>% filter(cell_type %ni% c("excluded"))
  
  dist <- cells_clean %>% 
    filter(ROI == {{ROI}}) %>% 
    select(c(centroid_x,centroid_y)) %>% 
    dist() %>% 
    as.matrix()
  
  # Find the index of the closest spot of cell_type "interest
  ix_interest <- cells_clean %>% 
    filter(ROI == {{ROI}}) %>% 
    pull({{cell_type_column}}) == {{cell_type}}
  
  
  # Calculate min distance for each spot to the closest spot of cell type {{cell_type}}
  min_dist_to_interest <- apply(dist, 1, function(x) sort(x[ix_interest])[n])
  
  # Calculate average distance per cell type to closest interest cell
  cells_clean <- cells_clean %>% 
    filter(ROI == {{ROI}}) %>% 
    select(cell_name, {{cell_type_column}}, ROI) %>% 
    cbind(min_dist_to_interest) %>% 
    group_by(!!sym(cell_type_column)) #%>% 
  # summarise(mean_dist = mean(min_dist_to_interest), 
  #           median_dist = median(min_dist_to_interest)) %>% 
  # arrange(desc(mean_dist))
  
  rm(dist)
  gc()
  
  return(cells_clean)
}


levels_marker <- c(
  #immune
  "CD45",
  "CD14",
  "CD163",
  "CD206",
  "CD1c",
  "MHCII",
  "CD11c",
  "CD4",
  "DCLAMP",
  "CD8",
  "CD3",
  "CD279",
  "CD69",
  "CD19",
  "CD138",
  
  #Fibro
  "aSMA",
  "CD90",
  "VCAN",
  "MT1MMP",
  "VIM",
  #"DCX",
  "FN1",
  
  #vasc
  "CD31",
  "PDPN",
  
  #Epi
  "GLUT1",
  "NDRG1",
  "CA9",
  "Ki67",
  "POSTN",
  
  #Telocyte
  "PDGFRA",
  "TNR",
  "TNC",
  "HPLN1",
  
  #MALT
  "CD24",
  "IL1B",
  
  #rest
  "CD44",
  "ACAN",
  "GAP43",
  "area",
  "circularity",
  "solidity"

  
)

cols <- c(
  "mac" = "#d0003f",
  "Bcell" = "grey",
  "Pcell" = "#4c617c",
  "cDC2" = "#cf54c4",
  "cDC1" = "#7166d1",
  "FDC" = "#45005a",
  "Tcell_CD4" = "#465637",
  "Tcell_CD8" = "#77da4f",
  "fib" = "#ff9a84",
  "fib_BM" = "#ff5420",
  "endo_arterial" = "#69d29e",
  "endo_lymph" = "#7cc0c6",
  "endo_venous" = "#cfd351",
  "epi_GLUT1" = "#c9c296",
  "epi_NDRG1" = "#ca8c3f",
  "epi_CA9" = "#884834",
  "telocyte" = "#c1457a",
  "telocyte_tip" = "#ce94a9"
)

cols_fig <- c(
  "Macrophage" = "#d0003f",
  "B cell" = "grey",
  "Plasma cell" = "#4c617c",
  "cDC2" = "#cf54c4",
  "cDC1" = "#7166d1",
  "FDC" = "#45005a",
  "T cell CD4" = "#465637",
  "T cell CD8" = "#77da4f",
  "Fibro." = "#ff9a84",
  "Fibro. BM" = "#ff5420",
  "Endo. Art." = "#69d29e",
  "Endo. Lymp." = "#7cc0c6",
  "Endo. Ven." = "#cfd351",
  "Epi. GLUT1" = "#c9c296",
  "Epi. NDRG1" = "#ca8c3f",
  "Epi. CA9" = "#884834",
  "Telocyte" = "#c1457a",
  "Telocyte tip" = "#ce94a9"
)

# Import and tidy data ----------------------------------------------------
cells <- data.table::fread("cells.csv") %>% as_tibble()
amat <- data.table::fread("norm_mat.csv") %>% as.data.frame()
marker <- read_csv2("marker.csv", col_names = F)
marker <- marker$X1

rownames(amat) <- marker
# Overview ----------------------------------------------------------------
abundancies <- cells %>%
  #filter(cell_type %ni% c("excluded")) %>%
  group_by(cell_type) %>%
  summarise(percentage = round( n() / nrow(cells) * 100, digits = 1) ) %>% 
  mutate(cell_class = case_when(
    cell_type %in% c("Bcell", "cDC1", "cDC2", "FDC", "mac", "Pcell", "Tcell_CD4", "Tcell_CD8") ~ "Immune",
    cell_type %in% c("epi_CA9", "epi_NDRG1", "epi_GLUT1") ~ "Epithelial",
    cell_type %in% c("endo_arterial", "endo_lymph", "endo_venous", "fib_BM", "fib", "telocyte", "telocyte_tip") ~ "Mesenchyme",
    T ~ cell_type)
  ) %>% 
  mutate(cell_type_fig = case_when(
    cell_type == "Bcell" ~ "B cell",
    cell_type == "endo_arterial" ~ "Endo. Art.",
    cell_type == "endo_lymph" ~ "Endo. Lymp.",
    cell_type == "endo_venous" ~ "Endo. Ven.",
    cell_type == "epi_CA9" ~ "Epi. CA9",
    cell_type == "epi_NDRG1" ~ "Epi. NDRG1",
    cell_type == "epi_GLUT1" ~ "Epi. GLUT1",
    cell_type == "fib" ~ "Fibro.",
    cell_type == "fib_BM" ~ "Fibro. BM",
    cell_type == "mac" ~ "Macrophage",
    cell_type == "Pcell" ~ "Plasma cell",
    cell_type == "Tcell_CD4" ~ "T cell CD4",
    cell_type == "Tcell_CD8" ~ "T cell CD8",
    cell_type == "telocyte" ~ "Telocyte",
    cell_type == "telocyte_tip" ~ "Telocyte tip",
    T ~ cell_type)
  ) %>% 
  left_join(., tibble(cell_type = names(cols), col_code = cols), by = "cell_type") %>% 
  mutate(voronoi = paste0(cell_type_fig, "\n", percentage, "%"))

abundancies <- abundancies[with(abundancies, order(cell_class, voronoi)), ] #order for Voronoi




ggplot(abundancies, aes(area=percentage, 
                fill=cell_type_fig, 
                label=paste0(cell_type_fig,"\n", paste0(percentage), "%"), 
                subgroup=cell_class)) +
  treemapify::geom_treemap() +
  treemapify::geom_treemap_text(place = "centre") +
  treemapify::geom_treemap_subgroup_border(colour = "white", size = 10) +
  treemapify::geom_treemap_subgroup_text(place = "centre", grow = TRUE,
                             alpha = 0.25, colour = "black",
                             fontface = "italic") +
  scale_fill_manual(values = cols_fig) +
  labs(fill = "cell type")

# Voroni abundancy plot
# generate treemap; set seed to obtain same pattern every time
tm <- voronoiTreemap(
  data = abundancies,
  levels = c("cell_class","voronoi"),
  cell_size = "percentage",
  shape = "rounded_rect",
  seed = 100
)

drawTreemap(tm, 
            label_size = 2.5, 
            label_color = "white", 
            color_palette = abundancies$col_code, 
            color_level = 2, border_size = 10)  

tm2 <- sunburstTreemap(
  data = abundancies,
  levels = c("cell_class","cell_type_fig"),
  cell_size = "percentage"
)

drawTreemap(tm2, 
            label_size = 2.5, 
            color_palette = cols_fig, 
            label_level = c(2), 
            label_color = "grey" )

# Heatmap -----------------------------------------------------------------
tmp <- cells
tmp <- map(tmp$sample %>% unique, ~Run_DEPs(.))
tmp <- do.call(cbind, tmp)
tmp <- tmp %>% select(order(colnames(tmp)))

Heatmap(t(tmp[levels_marker,]),
        col = circlize::colorRamp2(breaks = seq(-2, 2, length.out = 11), hotmap),
        #col = viridis(option = "V", n = 10),
        row_names_side = "left",
        show_row_names = F,
        row_dend_side = NULL,
        column_order = levels_marker,
        row_gap = unit(1, "mm"),
        border = T,
        name = "z-score",
        cluster_row_slices = F,
        cluster_column_slices = F,
        row_split = c(
          rep("Bcell", 3),
          rep("cDC1", 3),
          rep("cDC2", 3),
          rep("endo_arterial", 3),
          rep("endo_lymph", 3),
          rep("endo_venous", 3),
          rep("epi_CA9", 3),
          rep("epi_GLUT1", 3),
          rep("epi_NDRG1", 3),
          rep("FDC", 3),
          rep("fib", 3),
          rep("fib_BM", 3),
          rep("mac", 3),
          rep("Pcell", 3),
          rep("Tcell_CD4", 3),
          rep("Tcell_CD8", 3),
          rep("telocyte", 3),
          rep("telocyte_tip", 3))
        %>% 
          factor(., levels = c("mac",
                               "cDC2",
                               "cDC1",
                               "FDC",
                               "Tcell_CD4",
                               "Tcell_CD8",
                               "Bcell",
                               "Pcell",
                               "fib",
                               "fib_BM",
                               "endo_arterial",
                               "endo_lymph",
                               "endo_venous",
                               "epi_GLUT1",
                               "epi_NDRG1",
                               "epi_CA9",
                               "telocyte",
                               "telocyte_tip")),
        left_annotation = rowAnnotation(cell_type = c(
          rep("Bcell", 3),
          rep("cDC1", 3),
          rep("cDC2", 3),
          rep("endo_arterial", 3),
          rep("endo_lymph", 3),
          rep("endo_venous", 3),
          rep("epi_CA9", 3),
          rep("epi_GLUT1", 3),
          rep("epi_NDRG1", 3),
          rep("FDC", 3),
          rep("fib", 3),
          rep("fib_BM", 3),
          rep("mac", 3),
          rep("Pcell", 3),
          rep("Tcell_CD4", 3),
          rep("Tcell_CD8", 3),
          rep("telocyte", 3),
          rep("telocyte_tip", 3)),
          col = list(cell_type = cols),
          show_legend = F),
        column_split = c(rep("immune", 15),
                         rep("fibro", 6),
                         rep("vasc", 2),
                         rep("epi", 5),
                         rep("telo", 4),
                         rep("MALT", 2),
                         rep("other", 6))  %>%
          factor(., levels = c("immune", "fibro", "vasc", "epi", "telo", "MALT", "other")),
        row_title_rot = 0,
        column_title_gp = gpar(fontsize = 13, fontface = "bold", col = "black"),
        row_title_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
        column_names_gp = gpar(col = "black", fontsize = 13, fontfamily = "Helvetica"),
        show_row_dend = F,
        heatmap_legend_param = list(title_gp = gpar(col = "black"), labels_gp = gpar(col = "black"), border = "black"))




# Spatial map -------------------------------------------------------------

### Specific cell types
ggplot(cells %>% filter(cell_type %in% c("endo_arterial", "endo_venous")), aes(x=centroid_x, y=centroid_y, color=cell_type)) +
  geom_point(size=0.1) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha = 1))) +
  labs(color = "cell_type") +
  facet_wrap(~sample, scales = "free") +
  dark_theme_gray() +
  scale_y_reverse()

### Specific cell types
ggplot(cells %>% filter(cell_type %in% c("telocyte", "telocyte_tip")), aes(x=centroid_x, y=centroid_y, color=cell_type)) +
  geom_point(size=0.1) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha = 1))) +
  labs(color = "cell_type") +
  facet_wrap(~sample, scales = "free") +
  dark_theme_gray() +
  scale_y_reverse()

### Villi example spatial maps
ggplot(cells %>% filter(ROI == "P9_Straight_villi_03") %>% filter(cell_type %ni% c("excluded")),
       aes(x=centroid_x, y=centroid_y, color=cell_type_fig)) +
  geom_point(size=1) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha =1))) +
  scale_color_manual(values = cols_fig) +
  labs(color = "cell_type") +
  theme_void() +
  scale_y_reverse()


### cDC example
ggplot(cells %>% filter(ROI == "P13_Straight_villi_02") %>% filter(cell_type_2 %in% c("epi", "cDC1", "cDC2", "fib_BM")),
       aes(x=centroid_x, y=centroid_y, color=cell_type_fig)) +
  geom_point(size=1) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha =1))) +
  scale_color_manual(values = cols_fig) +
  labs(color = "cell_type") +
  theme_void() +
  scale_y_reverse()

### Branched villi example
ggplot(cells %>% 
         filter(ROI %in% c("P8_Serrated_villi_05", "P8_Straight_villi_13")) %>% 
         filter(cell_type %ni% c("excluded")),
       aes(x=centroid_x, y=centroid_y, color=cell_type_fig)) +
  geom_point(size=1) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha =1))) +
  scale_color_manual(values = cols_fig) +
  labs(color = "cell_type") +
  theme_void() +
  scale_y_reverse()

### Branched villi example
ggplot(cells %>% filter(ROI == "P13_ILF_02") %>% filter(cell_type %ni% c("excluded")),
       aes(x=centroid_x, y=centroid_y, color=cell_type_fig)) +
  geom_point(size=0.5) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha =1))) +
  scale_color_manual(values = cols_fig) +
  labs(color = "cell_type") +
  theme_void() +
  scale_y_reverse()

# Min dist analyses -------------------------------------------------------

## To Epi for all villi types
tictoc::tic()
min_dist_to_epi <- lapply(X = cells %>% pull("ROI") %>% unique(),
                                   FUN = function(x) Min_dist_to_nth_cell_type_distribution(cell_type = "epi",
                                                                                            ROI = x,
                                                                                            cell_type_column = "cell_type_2"))
tictoc::toc()


min_dist_to_epi_df <- do.call(rbind, min_dist_to_epi)

## Exclude ILF villi and cells that are below the BM (e.g. min dist to epi)
min_dist_to_epi_df <- min_dist_to_epi_df %>% 
  filter(min_dist_to_interest < 100) %>% 
  filter(ROI %ni% c(cells %>% filter(villi_type %in% c("ILF")) %>% pull(ROI) %>% unique()))


## Add cell_class
min_dist_to_epi_df <- min_dist_to_epi_df %>% 
  mutate(cell_class = case_when(
    cell_type_2 %in% c("Bcell", "cDC1", "cDC2", "FDC", "mac", "Pcell", "Tcell_CD4", "Tcell_CD8") ~ "Immune",
    cell_type_2 %in% c("epi") ~ "Epithelial",
    cell_type_2 %in% c("endo_arterial", "endo_lymph", "endo_venous", "fib_BM", "fib", "telocyte", "telocyte_tip") ~ "Mesenchyme",
    T ~ cell_type_2)
  )

## Add cell_type_fig
min_dist_to_epi_df <- min_dist_to_epi_df %>% 
  mutate(cell_type_2_fig = case_when(
    cell_type_2 == "Bcell" ~ "B cell",
    cell_type_2 == "endo_arterial" ~ "Endo. Art.",
    cell_type_2 == "endo_lymph" ~ "Endo. Lymp.",
    cell_type_2 == "endo_venous" ~ "Endo. Ven.",
    cell_type_2 == "epi" ~ "Epi.",
    cell_type_2 == "fib" ~ "Fibro.",
    cell_type_2 == "fib_BM" ~ "Fibro. BM",
    cell_type_2 == "mac" ~ "Macrophage",
    cell_type_2 == "Pcell" ~ "Plasma cell",
    cell_type_2 == "Tcell_CD4" ~ "T cell CD4",
    cell_type_2 == "Tcell_CD8" ~ "T cell CD8",
    cell_type_2 == "telocyte" ~ "Telocyte",
    cell_type_2 == "telocyte_tip" ~ "Telocyte tip",
    T ~ cell_type_2)
  )

## Plot violin with cell class facets for all cells
ggplot(min_dist_to_epi_df %>%  filter(cell_type_2_fig %ni% c("excluded", "B cell", "FDC", "Epi.")), 
       aes(y=min_dist_to_interest, 
           x=fct_reorder(cell_type_2_fig, min_dist_to_interest, .desc = F),
           fill=cell_type_2_fig)) + 
  geom_violin(width=0.8, scale = "width", alpha = 0.8) +
  geom_boxplot(width=0.1, color="darkgray", outlier.shape = NA) +
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 20,
    size = 3,
    color = "white",
    position = position_dodge(width = 0.2)
  ) +
  scale_fill_manual(values = cols_fig) +
  facet_grid(~cell_class, scales = "free", space = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 12),
        strip.text = element_text(size=12)) +
  labs(x="Cell type", y="Distance to closest Enterocyte [um]")


## Plot violin closest to enterocyte T cells 
ggplot(min_dist_to_epi_df %>%  filter(cell_type_2 %in% c("Tcell_CD4", "Tcell_CD8")), 
       aes(x=min_dist_to_interest, 
           y=fct_reorder(cell_type_2_fig, min_dist_to_interest, .desc = F),
           fill=cell_type_2_fig)) + 
  geom_violin(width=0.8, scale = "width", alpha = 0.8) +
  geom_boxplot(width=0.1, color="darkgray", outlier.shape = NA) +
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 20,
    size = 3,
    color = "white",
    position = position_dodge(width = 0.2)
  ) +
  scale_fill_manual(values = cols_fig) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 12),
        strip.text = element_text(size=12)) +
  labs(y="Cell type", x="Distance to closest Enterocyte [um]") +
  #xlim(0,60) +
  ggpubr::stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("T cell CD4", "T cell CD8")),
                     label = "p.format")

## Plot violin closest to enterocyte cDCs
ggplot(min_dist_to_epi_df %>%  filter(cell_type_2 %in% c("cDC1", "cDC2")), 
       aes(x=min_dist_to_interest, 
           y=fct_reorder(cell_type_2_fig, min_dist_to_interest, .desc = F),
           fill=cell_type_2_fig)) + 
  geom_violin(width=0.8, scale = "width", alpha = 0.8) +
  geom_boxplot(width=0.1, color="darkgray", outlier.shape = NA) +
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 20,
    size = 3,
    color = "white",
    position = position_dodge(width = 0.2)
  ) +
  scale_fill_manual(values = cols_fig) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 12),
        strip.text = element_text(size=12)) +
  labs(y="Cell type", x="Distance to closest Enterocyte [um]") +
  #xlim(0,60) +
  ggpubr::stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("cDC1", "cDC2")),
                     label = "p.format")

## Crypt villus axis
tictoc::tic()
min_dist_to_BM <- lapply(X = cells %>% pull("ROI") %>% unique(),
                          FUN = function(x) Min_dist_to_nth_cell_type_distribution(cell_type = "fib_BM",
                                                                                   ROI = x,
                                                                                   cell_type_column = "cell_type"))
tictoc::toc()


min_dist_to_BM_df <- do.call(rbind, min_dist_to_BM)

## Exclude ILF villi and cells that are below the BM (e.g. min dist to epi)
min_dist_to_BM_df <- min_dist_to_BM_df %>% 
  filter(cell_name %in% min_dist_to_epi_df$cell_name)

## Add cell_class
min_dist_to_BM_df <- min_dist_to_BM_df %>% 
  mutate(cell_class = case_when(
    cell_type %in% c("Bcell", "cDC1", "cDC2", "FDC", "mac", "Pcell", "Tcell_CD4", "Tcell_CD8") ~ "Immune",
    cell_type %in% c("epi_CA9", "epi_GLUT1", "epi_NDRG1") ~ "Epithelial",
    cell_type %in% c("endo_arterial", "endo_lymph", "endo_venous", "fib_BM", "fib", "telocyte", "telocyte_tip") ~ "Mesenchyme",
    T ~ cell_type)
  )

## Add cell_type_fig
min_dist_to_BM_df <- min_dist_to_BM_df %>% 
  mutate(cell_type_fig = case_when(
    cell_type == "Bcell" ~ "B cell",
    cell_type == "endo_arterial" ~ "Endo. Art.",
    cell_type == "endo_lymph" ~ "Endo. Lymp.",
    cell_type == "endo_venous" ~ "Endo. Ven.",
    cell_type == "epi_CA9" ~ "Epi. CA9",
    cell_type == "epi_NDRG1" ~ "Epi. NDRG1",
    cell_type == "epi_GLUT1" ~ "Epi. GLUT1",
    cell_type == "fib" ~ "Fibro.",
    cell_type == "fib_BM" ~ "Fibro. BM",
    cell_type == "mac" ~ "Macrophage",
    cell_type == "Pcell" ~ "Plasma cell",
    cell_type == "Tcell_CD4" ~ "T cell CD4",
    cell_type == "Tcell_CD8" ~ "T cell CD8",
    cell_type == "telocyte" ~ "Telocyte",
    cell_type == "telocyte_tip" ~ "Telocyte tip",
    T ~ cell_type)
  )

## Plot violin with cell class facets
ggplot(min_dist_to_BM_df %>%  filter(cell_type %ni% c("excluded", "fib_BM", "FDC", "Bcell")), 
       aes(y=min_dist_to_interest, 
           x=fct_reorder(cell_type_fig, min_dist_to_interest, .desc = F),
           fill=cell_type_fig)) + 
  geom_violin(width=0.8, scale = "width", alpha = 0.8) +
  geom_boxplot(width=0.1, color="darkgray", outlier.shape = NA) +
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 20,
    size = 3,
    color = "white",
    position = position_dodge(width = 0.2)
  ) +
  scale_fill_manual(values = cols_fig) +
  facet_grid(~cell_class, scales = "free", space = "free") +
  theme_classic() +
  labs(x="") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 12),
        strip.text = element_text(size=12)) +
  labs(x="Cell type", y="distance from basal membrane [um]")

## Plot violin closest to enterocyte cDCs
ggplot(min_dist_to_BM_df %>%  filter(cell_type %in% c("cDC1", "cDC2")), 
       aes(y=min_dist_to_interest, 
           x=fct_reorder(cell_type, min_dist_to_interest, .desc = F),
           fill=cell_type)) + 
  geom_violin(width=0.8, scale = "width", alpha = 0.8) +
  geom_boxplot(width=0.1, color="darkgray", outlier.shape = NA) +
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 20,
    size = 3,
    color = "white",
    position = position_dodge(width = 0.2)
  ) +
  scale_fill_manual(values = cols_fig) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 12),
        strip.text = element_text(size=12)) +
  labs(x="Cell type", y="Distance to closest Enterocyte [um]") +
  #xlim(0,60) +
  ggpubr::stat_compare_means(method = "wilcox.test",
                             comparisons = list(c("cDC1", "cDC2")),
                             label = "p.format")

# Export ------------------------------------------------------------------
## QuPath
Make_qpan(export = "cell_type")
