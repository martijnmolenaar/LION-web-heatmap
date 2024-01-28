# LION/web: LION heatmap module

App available online via [lipidontology.com](http://heatmap.lipidontology.com)

<img src="https://raw.githubusercontent.com/martijnmolenaar/LION-web-heatmap/main/HeatmapApp/www/LIONicon%20heatmap.png" alt="LION logo">


Please [cite](https://academic.oup.com/gigascience/article/8/6/giz061/5505544):
> **LION/web: a web-based ontology enrichment tool for lipidomic data analysis.**
> Martijn R Molenaar,  Aike Jeucken,  Tsjerk A Wassenaar,  Chris H A van de Lest, Jos F Brouwers,  J Bernd Helms. 
> *GigaScience, Volume 8, Issue 6, June 2019, giz061, https://doi.org/10.1093/gigascience/giz061*

> **Lipidomic profiling of rat hepatic stellate cells during activation reveals a two-stage process accompanied by increased levels of lysosomal lipids.**
> Martijn R Molenaar, Maya W Haaker, A Bas Vaandrager, Martin Houweling, J Bernd Helms.
> *J Biol Chem. 2023 Feb 18;299(4):103042. doi: 10.1016/j.jbc.2023.103042.*

## installation of R-packages required for LION/web
```R
if("devtools" %in% rownames(installed.packages()) == FALSE) {install.packages("devtools",  repos = c(CRAN = "http://cran.rstudio.com"))}
library(devtools)


if("shiny" %in% rownames(installed.packages()) == FALSE) {install.packages("shiny",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinyBS" %in% rownames(installed.packages()) == FALSE) {install.packages("shinyBS",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinythemes" %in% rownames(installed.packages()) == FALSE) {install.packages("shinythemes",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinyTree" %in% rownames(installed.packages()) == FALSE) {install.packages("shinyTree",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinyWidgets" %in% rownames(installed.packages()) == FALSE) {install.packages("shinyWidgets",  repos = c(CRAN = "http://cran.rstudio.com"))}

if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("ggthemes" %in% rownames(installed.packages()) == FALSE) {install.packages("ggthemes",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("ggplotify" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplotify",  repos = c(CRAN = "http://cran.rstudio.com"))}

if("httr" %in% rownames(installed.packages()) == FALSE) {install.packages("httr",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("pheatmap" %in% rownames(installed.packages()) == FALSE) {install.packages("pheatmap",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("colortools" %in% rownames(installed.packages()) == FALSE) {install.packages("colortools",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("colourpicker" %in% rownames(installed.packages()) == FALSE) {install.packages("colourpicker",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("cowplot" %in% rownames(installed.packages()) == FALSE) {install.packages("cowplot",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("formattable" %in% rownames(installed.packages()) == FALSE) {install.packages("formattable",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("sortable" %in% rownames(installed.packages()) == FALSE) {install.packages("sortable",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("visNetwork" %in% rownames(installed.packages()) == FALSE) {install.packages("visNetwork",  repos = c(CRAN = "http://cran.rstudio.com"))}

if("topOnto" %in% rownames(installed.packages()) == FALSE) {install_github("martijnmolenaar/topOnto")}

##install topOnto.LION.db package:
if("topOnto.LION.db" %in% rownames(installed.packages()) == FALSE) {install_github("martijnmolenaar/topOnto.LION2.db/topOnto.LION.db")}

source("https://bioconductor.org/biocLite.R")
biocLite()


if("RSQLite" %in% rownames(installed.packages()) == FALSE) {install.packages("RSQLite",  repos = c(CRAN = "http://cran.rstudio.com"))}

if("ggrepel" %in% rownames(installed.packages()) == FALSE) {install.packages("ggrepel",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinycssloaders" %in% rownames(installed.packages()) == FALSE) {install.packages("shinycssloaders",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("jsonlite" %in% rownames(installed.packages()) == FALSE) {install.packages("jsonlite",  repos = c(CRAN = "http://cran.rstudio.com"))}

```

## Load all necessary libraries to check installation

```R
library(shiny)
library(shinyTree)
library(shinyWidgets)
library(shinyBS)
library(shinythemes)

library(data.table)
library(ggplot2)
library(ggthemes)
library(ggplotify)

library(httr)
library(pheatmap)
library(RColorBrewer)
library(colortools)
library(colourpicker)
library(cowplot)
library(formattable)
library(sortable)
library(visNetwork)


## loading lipid ontology data
library(RSQLite)
library(topOnto)
library('topOnto.LION.db')
topOnto::initONT('LION')
```

## run LION/web app
```R
runApp("HeatmapApp")
```
