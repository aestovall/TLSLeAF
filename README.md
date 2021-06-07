![TLSLeAF](TLSLeAF.png)

# TLSLeAF - Terrestrial Laser Scanning Leaf Angle Function
Algorithms described in "Automatic leaf angle estimates from single-scan terrestrial laser scanning with TLSLeAF" (Stovall et al. 2021. New Phytologist) for generating leaf angle distributions (LADs) from single-scan terrestrial laser scanning (TLS) data. 

TLSLeAF is based in the R programming language and takes advantage of CloudCompare commandline tools. 

## EDIT THE SETUP FILE FOR YOUR SYSTEM

Make sure the setup file is correct and save. The most important step is to have Cloud Compare installed and add the executable path into the setup file.

```{r,echo=FALSE}
file.edit('R/000_setup.R')
```

Add your operating system type and specify the path to CloudCompare. Further testing is necessary to ensure this works on a Windows machine.

Several setup option are available: 

```{r,echo=FALSE}
scatterLim = 85          #Threshold of scattering angle to remove'
correct.topography = TRUE #topographically normalize the voxels?
SS = 0.02                 #Spatial subsampling resolution
scales = c(0.1,0.5,0.75)  #3 scales of normal computation. Set for RF model
voxRes = 0.1             #voxel resolution for LAD normalization
minVoxDensity = 5         #minimum number of measurments per voxel
superDF = TRUE            #Merges output into a single TLSLeAF.class
clean = TRUE              #removes temporary files created during processing
```
`scatterLim` is a filter based on the angle of reflection of the laser pulse, relative to the normal of the leaf surface. For instance, a 0 degree scattering angle would be found for surfaces perpendicular to the laser trajectory. A 90 degree angle would be a surface parallel to the laser trajectory. The algorithm filters scattering angles greater than 85 degrees by default.

`correct.topography` creates a surface model from the TLS data and normalizes returns and angles according to the ground surface.

`SS` is the resolution of the spatial subsample. TLSLeAF defaults to 0.02 m or 2 cm point spacing since this resolution retains substantial detail, but minimizes computational needs.

`scales` are the three octree sizes at which normals are computed for the leaf wood classification. The values here MUST align with those of the classifier imported into the environment - in the case of TLSLeAF, `rf_model`.

`voxRes` is the resolution of the voxelization for normalizing density of leaf angle measurments. Laboratory test show 0.1 or 10 cm to be an appropriate resolution for this step, but the parameter may be set to a more coarse resolution (e.g. 1 m) if desired.

`minVoxDensity` is the minimum number of measurments per voxel that are required to retain a voxel. Default is 5 measurments, but should be adjusted depending on the voxel resolution (`voxRes`).

`superDF` is a logical option (`TRUE`/`FALSE`) to have the intermediate and final output of TLSLeAF be saved in the environment as a `TLSLeAF.class`. The `TLSLeAF.class` has the following structure (class):

```{r,echo=FALSE}
parameters ("data.frame")
dat ("data.frame")
voxels ("data.frame")
LAD ("data.frame")
Beta_parameters ("data.frame")
```
All slots can be accessed easily (e.g. `dat@LAD`) from a single object.

`clean` is a logical option. When set to `TRUE` the temporary files in the `input` folder will be removed.

## SET UP CLOUDCOMPARE

TLSLeAF relies on CloudCompare commandline for key portions of the processing. To setup CloudCompare for TLSLeAF you must first designate your operating system:

```{r,echo=FALSE}
OS<-"mac" 
OS<-"windows"
```
Depending on the operating system selected, note if CloudCompare is already in your `PATH`.

```{r,echo=FALSE}
PATH = FALSE
```

If not, please add the location of the CloudCompare executable/application. For example:

```{r,echo=FALSE}
#MAC
"/Applications/CloudCompare.app/Contents/MacOS/CloudCompare"

#WINDOWS
shQuote('C:\\Program Files\\CloudCompare\\CloudCompare.exe')
```

Now, TLSLeAF should be properly configured and ready to run!

## RUN THE PIPELINE
Load packages, functions, and input parameters.
```{r,echo=FALSE}
source('R/000_setup.R')
```

Load in your TLS file in PTX format (other gridded formats may work, but are untested).
```{r,echo=FALSE}
files<-list.files(pattern = "ptx", 
                  recursive = TRUE, 
                  full.names = TRUE)
input_file = files[1]

```

Find center coordinates. Modify the `center` variable according to file format.
```{r,echo=FALSE}

center<-fread(input_file, nrows=4, header=FALSE)[1,]
colnames(center)<-c("x","y","z")
```

FINALLY, Run TLSLeAF. The output from the following default parameters will be of the `TLSLeAF.class`, allowing access to all intermediate products used to create the leaf angle distribution.
```{r,echo=FALSE}
df<-TLSLeAF(input_file, 
            overwrite=TRUE,
            center, 
            scatterLim=85,
            SS=0.02, 
            scales=c(0.1,0.5,0.75),
            rf_model,
            voxRes=0.1,
            minVoxDensity=5,
            superDF=TRUE)
```

Success (hopefully)! You will find all of the processed files in the `output` directory.
The `figures` folder should also include a LAD figure with beta parameters from TLSLeAF:

![LAD](LAD.png)



## How fast is it?

TLSLeAF is an efficient means of calculating single-scan leaf angles and LADs. One can expect processing times in line with the figure below on a consumer-grade machine.

![Processing_speed](Processing_speed.png)
