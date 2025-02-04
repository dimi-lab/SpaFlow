library(data.table)
library(dplyr)

all_files <- list.files(".", pattern = "*.csv")

all_rois <- sub("som_metaclusters_", "", all_files)
all_rois <- unique(sub("_.clusters.csv", "", all_rois))

for (roi in all_rois) {
  all_files_roi <- grep(roi, all_files, value = T)
  
  dat <- fread(all_files_roi[1])
  name_parts <- strsplit(all_files_roi[1], "_")[[1]]
  som_number <- sub("clusters.csv", "", name_parts[length(name_parts)])
  new_column_name <- paste0("som_clusters_", som_number)
  col_ind <- which(colnames(dat) == "som_cluster")
  colnames(dat)[col_ind] <- new_column_name
  
  if (length(all_files_roi) > 1) {
    for (file in all_files_roi[2:length(all_files_roi)]) {
      name_parts <- strsplit(file, "_")[[1]]
      som_number <- sub("clusters.csv", "", name_parts[length(name_parts)])
      new_column_name <- paste0("som_clusters_", som_number)
      
      dat_new <- fread(file)
      col_ind <- which(colnames(dat_new) == "som_cluster")
      colnames(dat_new)[col_ind] <- new_column_name
      
      dat <- cbind(dat, dat_new[,..col_ind])
    }
  }
  
  dat <- dat %>% select(-cluster)
  
  fwrite(dat, paste0("som_metaclusters_", roi, ".csv"))
}

