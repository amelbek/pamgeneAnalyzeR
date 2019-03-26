## ---- include=TRUE-------------------------------------------------------

library(pamgeneAnalyzeR)
library(tiff) # load tiff library to read tiff images 

example_files_path <- system.file("extdata", package = "pamgeneAnalyzeR")

example_files <- list.files(example_files_path)

example_files

## ----plot target , include=TRUE, message=FALSE, warning=FALSE, fig.align='center', fig.width = 7, fig.height = 7----
target <- readTIFF(source=paste(example_files_path,'/', example_files[1], sep = ''))
plot_centers(image = target, centers_coords = NULL)

## ----plot centers , include=TRUE, fig.align='center', fig.width = 7, fig.height = 7----
cc <- find_centers(target) ## find centers coordinates

plot_centers(image = target, centers_coords = cc) #plot 


## ----plot circles , include=TRUE, message=FALSE, warning=FALSE, fig.align='center', fig.width = 7, fig.height = 7----
plot_circles(centers.coords = cc, img = target, radius = 5)

## ---- include=TRUE, fig.align='center', message=FALSE, warning=FALSE-----

STK_path = '../../STK_Array_Layout.txt' # use your own PamChip layout path 
STK.layout.file <- read.table(file = STK_path, header = TRUE, sep = '\t', quote = "")

reg <- regiter_one_image(source = readTIFF(source=paste(example_files_path,'/', example_files[2], sep = '')), target = target, centers.coords = cc, parallel.sig.extract = TRUE, pamgene.layout.file = STK.layout.file)

tail(reg$signal)

## ---- include=TRUE, fig.align='center', message=FALSE, warning=FALSE, results='asis', fig.width = 7, fig.height = 7----
check_registration(source = readTIFF(source=paste(example_files_path,'/', example_files[2], sep = '')), target = target)

## ---- include=TRUE, fig.align='center', message=FALSE, warning=FALSE, eval=FALSE----
#  #Takes time
#  registerAllImages(inDirectory = example_files_path, outDirectory = "outDir", target = target, centers.coords = cc, parallel.registration = TRUE, pamgene.layout.file = STK.layout.file)
#  

## ---- include=TRUE, fig.align='center', fig.width = 5, fig.height = 4----
plot_image_from_signal(reg$signal)

## ---- include=TRUE, fig.align='center', eval=FALSE-----------------------
#  
#  rawData <- mergeExperiments(path = "outDir/" )
#  normData <- substractBackground(rawData = rawData)
#  
#  ZpFactor <- ZpFactor_forAll(rawData = normData)

