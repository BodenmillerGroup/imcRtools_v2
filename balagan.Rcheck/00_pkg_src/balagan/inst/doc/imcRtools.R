## ---- echo=FALSE, results="hide"----------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE,
                        fig.retina = 0.75, crop = NULL)
library(BiocStyle)

## ----library, echo=FALSE------------------------------------------------------
library(imcRtools)

## ----quickstart-load-data-----------------------------------------------------
#library(data.table)
# load the single cell data
#cells <- as.data.frame(fread(file = "data/hubmap_processed_v2/cell.csv",stringsAsFactors = FALSE))

# load the image level metadata
#image_mat <- as.data.frame(read.csv(file = "data/hubmap_processed_v2/Image.csv",stringsAsFactors = FALSE))

# load the panel information
#panel_mat <- read.csv(file = "data/hubmap_processed_v2/panel.csv", sep= ",",  stringsAsFactors = FALSE )

# get an example file that contains the channel order
#tags <- read.csv( "data/protein/20191021_ZTMA256.1_slide3_TH_s0_p9_r1_a1_ac_full.csv", header = FALSE)

# load acqusition meta data
#acquisition_meta <- read.csv(file = "data/hubmap_processed_v2/Experiment.csv", stringsAsFactors = FALSE)

# load the clinical data. this is a clinical datatable that already contains the ImageNumbers. It has been prepared in the clinical metadata preparation.Rmd script (prepared for RNA and protein dataset separately)
#clinical_mat <- read.csv("data/protein/clinical_data_protein.csv",stringsAsFactors = FALSE)

## ----quickstart-plotPixels----------------------------------------------------
#sce <- generateSCE(cells = cells,
#                   cell_meta = c("ImageNumber", "ObjectNumber"),
#                   image_meta = image_mat,
#                   panel_meta = panel_mat,
#                   metal_tag = "Metal.Tag",
#                   channel_number = "channel", 
#                   target_names = "clean_target",
#                   imageNumber = "ImageNumber",
#                   objectNumber = "ObjectNumber",
#                   measurement_column = "Intensity_MeanIntensityComp_FullStackFiltered",
#                   scaling_column = "Scaling_FullStack")

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

