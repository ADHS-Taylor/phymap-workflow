library(dplyr)
library(tidyr)
library(treeio)
library(ape)
library(phytools)

timetree_file <- "results_virus_x_non_pruned/treetime/timetree.nexus"
metadata_file <- "results_virus_x_non_pruned/phymapr_metadata.tsv"
col_specimen  <- "accessionVersion"
col_date      <- "sampleCollectionDate"
col_location  <- "location"

tt_tree <- treeio::read.beast(timetree_file) %>% treeio::drop.tip("Reference")
meta <- read.delim(metadata_file, sep = "\t")

meta <- meta %>%
  dplyr::rename(
    SPECIMEN_NUMBER = dplyr::all_of(col_specimen),
    location = dplyr::all_of(col_location)
  )

if (col_date %in% colnames(meta) && col_date != "date") {
  meta <- meta %>% dplyr::rename(date = dplyr::all_of(col_date))
}

td <- tibble::as_tibble(tt_tree) %>%
  dplyr::mutate(date = as.numeric(as.character(date)))

print(paste("Number of rows in td:", nrow(td)))
print("First few rows of td:")
print(head(td))

all_dates <- td$date
print(paste("Number of non-NA dates in all_dates:", sum(!is.na(all_dates))))

tree_numeric <- treeio::as.treedata(td)
phy_jittered <- tree_numeric@phylo
phy_jittered$edge.length[phy_jittered$edge.length == 0] <- 1e-06 

tree_labels <- phy_jittered$tip.label
print(paste("Number of tips in tree:", length(tree_labels)))

meta_cleaned <- meta %>%
  dplyr::filter(SPECIMEN_NUMBER %in% tree_labels) %>%
  dplyr::distinct(SPECIMEN_NUMBER, .keep_all = TRUE)
print(paste("Number of rows in meta_cleaned:", nrow(meta_cleaned)))

meta_ordered <- meta_cleaned[match(tree_labels, meta_cleaned$SPECIMEN_NUMBER), ]
print("meta_ordered check:")
print(head(meta_ordered))

tip_lat_vec <- setNames(as.numeric(meta_ordered$lat), meta_ordered$SPECIMEN_NUMBER)
tip_lon_vec <- setNames(as.numeric(meta_ordered$long), meta_ordered$SPECIMEN_NUMBER)

print(paste("tip_lat_vec non-NA count:", sum(!is.na(tip_lat_vec))))
print(paste("tip_lon_vec non-NA count:", sum(!is.na(tip_lon_vec))))

node_lat_reconstructed <- phytools::fastAnc(phy_jittered, tip_lat_vec)
node_lon_reconstructed <- phytools::fastAnc(phy_jittered, tip_lon_vec)

all_lat <- c(as.numeric(tip_lat_vec), as.numeric(node_lat_reconstructed))
all_lon <- c(as.numeric(tip_lon_vec), as.numeric(node_lon_reconstructed))

print(paste("Total lat count:", length(all_lat)))
print(paste("Total lon count:", length(all_lon)))

edges <- phy_jittered$edge
print(paste("Number of edges:", nrow(edges)))

mcc_tab_no_na <- data.frame(
  startYear = all_dates[edges[,1]],
  endYear   = all_dates[edges[,2]],
  startLat  = all_lat[edges[,1]],
  startLon  = all_lon[edges[,1]],
  endLat    = all_lat[edges[,2]],
  endLon    = all_lon[edges[,2]]
)
print("mcc_tab_no_na dimensions before drop_na:")
print(dim(mcc_tab_no_na))
print(head(mcc_tab_no_na))

mcc_tab <- mcc_tab_no_na %>% tidyr::drop_na()
print("mcc_tab dimensions after drop_na:")
print(dim(mcc_tab))
print(head(mcc_tab))
