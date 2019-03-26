# pamgeneAnalyzeR

---
title: "PamgeneAnalyzeR"
author: "Amel Bekkar"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{PamgeneAnalyzeR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

#PamgeneAnalyzeR example pipeline

##Load example data

We load here 12 PamGene PamChip raw images that are provided as examples of raw data.
```

library(pamgeneAnalyzeR)
library(tiff) # load tiff library to read tiff images 

example_files_path <- system.file("extdata", package = "pamgeneAnalyzeR")

example_files <- list.files(example_files_path)

example_files
```

##Choose a reference image 

We will use here first image as a reference. This image will be used to align the other images on it and extract the coordinates of peptides and controls spot centers. 
The reference image should have at least clear bright left and right control spots, which are used to identify the optimal parameter for the radius of the spot during quantification.
```{r plot target , include=TRUE, message=FALSE, warning=FALSE, fig.align='center', fig.width = 7, fig.height = 7}
target <- readTIFF(source=paste(example_files_path,'/', example_files[1], sep = ''))
plot_centers(image = target, centers_coords = NULL)
```

##Find spots centers coordinates 
`find_centers` functions is used to detect the centers coordinates of the peptides spots, the left and right controls spots and 8 points of background. 
`plot_centers` function is then used to plot the original image with the centers detected as a control. 
```{r plot centers , include=TRUE, fig.align='center', fig.width = 7, fig.height = 7}
cc <- find_centers(target) ## find centers coordinates

plot_centers(image = target, centers_coords = cc) #plot 

```

##Define a radius
`plot_circles` function can be used to visualize the area of signal that will be extracted from each spot. The default radius used from the centers detected is 5.
A QC image is generated. The centers of the spots are used thereafter with the chosen radius to draw a circle on each spot, inside which the pixels brightness values are captured. 
```{r plot circles , include=TRUE, message=FALSE, warning=FALSE, fig.align='center', fig.width = 7, fig.height = 7}
plot_circles(centers.coords = cc, img = target, radius = 5)
```


##Image registration - Align all the other images on the reference one (target)

`regiter_one_image` function is used to align one image on the reference (target) image. The function will give as output the aligned image as well as the extracted signal from this image. If the layout of the chip (STK or PTK) is given it will be also added to the extracted signal. 


```

STK_path = '../STK_Array_Layout.txt' # use your own PamChip layout path 
STK.layout.file <- read.table(file = STK_path, header = TRUE, sep = '\t', quote = "")

reg <- regiter_one_image(source = readTIFF(source=paste(example_files_path,'/', example_files[2], sep = '')), target = target, centers.coords = cc, parallel.sig.extract = TRUE, pamgene.layout.file = STK.layout.file)

tail(reg$signal)
```

`check_registration` can be used to check the alignment by displaying the edges of the image on the target.

```
check_registration(source = readTIFF(source=paste(example_files_path,'/', example_files[2], sep = '')), target = target)
```


`registerAllImages` aligned all images in a directory on the reference, extract the signal fom them in the form of text files written in the chosen ouput directory. The process can be parralized to reduce time.
Automatically aligns all the other images onto the chosen reference image and the signal of each spot of each image is extracted. Summary statistics of pixels brightness (mean, median, standard deviation and sum) are collected for each peptide spot, the left and right control point as well as 8 control background points.

```
#Takes time
registerAllImages(inDirectory = example_files_path, outDirectory = "outDir", target = target, centers.coords = cc, parallel.registration = TRUE, pamgene.layout.file = STK.layout.file)

```


## Quality control image 
As an assessment of the quality of the process, a synthetic image can be reconstructed from the extracted signal and compared to the original one to ensure that the pixel values where extracted correctly and no drift or shift is observable

```
plot_image_from_signal(reg$signal)
```

##Background normalization and Z’-factor
`mergeExperiments` puts together all the extracted data in a single dataframe.

`substractBackground` subtracts from each spot the mean of background control spots for normalization.

`ZpFactor_forAll` calculates the Z’-factor for the assay QC. Z’-factor values between 0.5 and 1 are considered excellent, values between 0 and 0.5 may be acceptable, and values less than 0 indicate the assay is unlikely to be usable, since positive and negative control readout overlap heavily in the readout, indicating that the chip for a given experiment has had some failure.

```

rawData <- mergeExperiments(path = "outDir/" )
normData <- substractBackground(rawData = rawData)

ZpFactor <- ZpFactor_forAll(rawData = normData)
```


