![TLSLeAF](TLSLeAF.png)

# TLSLeAF
Algorithms described in "PAPER TITLE" (Stovall et al. Submitted) for generating leaf angle distributions (LADs) from single-scan terrestrial laser scanning (TLS) data.

## INITIAL SETUP

Add your gridded TLS data in the "input" folder.

```{r,echo=FALSE}
files<-list.files(pattern = "ptx", 
                  recursive = TRUE, 
                  full.names = TRUE)
i=1

input_file = files[i]
output_file = gsub(".ptx",".asc", input_file)
```

## EDIT THE SETUP FILE FOR YOUR SYSTEM

Make sure the setup file is correct and save. The most important step is to have Cloud Compare installed and add the executable path into the setup file.

```{r,echo=FALSE}
file.edit('R/000_setup.R')
```

Add your operating system type and specify the path to CloudCompare. Further testing is necessary to ensure this works on a Windows machine.

Two setup option are available: 

```{r,echo=FALSE}
SCATTER_LIM<-85 #threshold of scattering angle to remove
correct.topography = TRUE #topographically normalize the voxels?
```
`SCATTER_LIM` is a filter based on the angle of reflection of the laser pulse, relative to the normal of the leaf surface. For instance, a 0 degree scattering angle would be found for surfaces perpendicular to the laser trajectory. A 90 degree angle would be a surface parallel to the laser trajectory. The algorithm filters scattering angles greater than 85 degrees by default.

`correct.topography` creates a surface model from the TLS data and normalizes returns and angles according to the ground surface.

## RUN THE PIPELINE

Load packages, functions, and input parameters.
```{r,echo=FALSE}
source('R/000_setup.R')
```
Calculate normals for gridded TLS point cloud
```{r,echo=FALSE}
source('../R/00_calculateNormals.R')
```
Calculate scattering angle and leaf angle
```{r,echo=FALSE}
source('../R/01_calculateAngle.R')
```
Classify wood and leaf from random forest classifier
```{r,echo=FALSE}
source('../R/02_classify_wood_leaves.R')
```
Correct topography
```{r,echo=FALSE}
source('../R/03_normalize_topography.R')
```
Leaf angle voxelation and density normalization
```{r,echo=FALSE}
source('../R/04_voxelize.R')
```
Simulate LAD from voxel statistics
```{r,echo=FALSE}
source('../R/05_simulate_LAD.R')
```

Fit beta function and get beta parameters from LAD
```{r,echo=FALSE}
source('../R/06_fit_beta.R')
```

Success (hopefully)! You will find all of the processed files in the `output` directory.
The `figures` folder should also include a LAD figure with beta parameters from TLSLeAF.
