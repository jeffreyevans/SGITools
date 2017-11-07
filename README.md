# GeNetIt
GeNetIt R package for spatial graph-theoretic gravity modelling with implementation of spatial graph-theoretic genetic gravity models.
The model framework is applicable for other types of spatial flow questions. Includes functions for constructing spatial graphs, sampling and summarizing associated raster variables and building unconstrained and singly constrained gravity models.

Available functions in GeNetIt are:

​

         area.graph.statistics - Statistics for edges (lines) based on a defined scale (area).

         build.node.data - Build node data

         dmatrix.df - Distance matrix to data.frame

         dps - dps genetic distance matrix for Columbia spotted frog (Rana luteiventris)

         graph.statistics - Point sample and statistics for edges (lines)

         gravity - Gravity model

         knn.graph - Saturated or K Nearest Neighbor Graph

         plot.gravity - Generic plot function for a gravity model object

         predict.gravity - Predict gravity model

         print.gravity - Print gravity model

         ralu.model - Columbia spotted frog (Rana luteiventris) data for specifying gravity model. Note, the data.frame is already log transformed.

         ralu.site - Subset of site-level spatial point data for Columbia spotted frog (Rana luteiventris)

         rasters - Subset of raster data for Columbia spotted frog (Rana luteiventris)

         summary.gravity - generic summary function for gravity model objects

**Bugs**: Users are encouraged to report bugs here. Go to [issues](https://github.com/jeffreyevans/GeNetIt/issues) in the menu above, and press new issue to start a new bug report, documentation correction or feature request. You can direct questions to <jeffrey_evans@tnc.org>.

**To install `GeNetIt` in R use install.packages() to download curent stable release from CRAN** 

**or, for the development version, run the following (requires the devtools package):**
`devtools::install_github("jeffreyevans/GeNetIt")`
