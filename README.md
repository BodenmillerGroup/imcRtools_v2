# Balagan (בלגן)

Balagan is an R-package dedicated to the study of large-scale highly-multiplexed imaging datasets. It contains several tools allowing to :
- Normalize and process the original data, i.e transforming and clustering the data.
- Vizualise the data with various plotting functions.
- Study and quantify cell-cell interactions using various interaction scores derived from the spatial point pattern theory and tensor decomposition approaches.
- Infer the best parameters for an optimal spatial sampling strategy.
- Predict sparse multi-range interaction models (Implementation of PenGE package)

Balagan is based on the **SingleCellExperiment** object structure and is therefore compatible with a variety of other single-cell analysis tools.


# Installation


# Data loading

Data can be loaded throug different functions. The simplest manner is to load the output of the CellProfiler **(CP)** pipeline as it will be automaticall processed and annotated. This can be done by simply typing the following command :

```r
List_data = Load_IMC_data(Path_to_cell_file = "Desktop/analysis/cpout/cell.csv",
                         Path_to_panel_file ="Desktop/analysis/cpout/panel.csv")
sce = Create_SCE(List_data,dimension = "2D",Bit_mode = 16,N_core = 6)
```

The first function loads the CP output (cell.csv) together with the annotation panel file (panel.csv file) while the second functions aggregates the data into a SingleCellExperiment (sce) object. Please carefully selects the number of cores that will be used by the parallelized functions of Balagan (N_core function) !
For a full description of how those functions can be tuned to load other types of files please see their associated documentation files.

# Data normalization and clustering

Highly multiplexed intensity data can be transformed using various functions. In Balagan two types transforms are available :
+ The inverse hyperbolic sinus transform. This transformation is the standard variance stabilizing transform (VST) for many distributions, including Negative Binomial distribution. The cofactor is either provided by the user or is automatically computed by calculting the relation between the mean and the variance across channels.
+ The residual count transform. This transformation is only possible if the data provided represents the total metal/fluorescence intensity (and not the mean) and if the size (in pixel) of each cell is available in the sce object. Briefly for each channel a log-Poisson regression between the total intensity and the size of the cell is computed. The residuals are then extracted and transformed and are then considered as the new expression value. 
