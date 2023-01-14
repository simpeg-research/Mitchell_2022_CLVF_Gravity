# Mitchell_2022_CLVF_Gravity

This repository contains all of the files required to reproduce the SimPEG gravity inversions of the CLVF field dataset.


Gravity Dataset
-----------------------------------------------------------------------------------------------------------
File names: 

- **groundGrav_Combined_zEllipsoid_Full.pkl**
- **groundGrav_Combined_zEllipsoid_Full.csv**

For this study, a regional CLVF gravity dataset was compiled that includes 2,929 gravity stations. Measurements
used in this study include the dataset collected by Langenheim et al. (2006), datasets from
Chapman and Bishop (1974), Youngs et al. (1985), and Smith (1992) compiled by Langenheim et al. (2006), and a
2018 dataset collected in The Geysers (Peacock et al., 2020). All of the raw measurements were processed using
the standard gravity reduction methods outlined in Blakely (1995) to produce the complete Bouguer anomaly (CBA)
and isostatic residual gravity anomaly. The CBA includes the terrain correction, computed using BOUGUER
(Godson and Plouff, 1988), which accounts for the gravitational attraction of terrain above sea level. The
residual isostatic gravity anomaly, which accounts for regional variations in the upper mantle and lower crust
that compensate topographic loads, were calculated in the manner of Simpson et al. (1986). Additional details
regarding the standard USGS gravity reduction techniques used here can be found in Langenheim et al. (2006) and
Langenheim et al. (2007).

The field dataset is provided in two different formats (.csv and .pkl). Both datafiles were created from a Pandas
Dataframe with 12 data columns: 

 - **Station_ID**: Name/ID assigned to the gravity station
 - **lonWGS84**: Longitude of the gravity station (CRS: WGS84)
 - **latWGS84**: Latitude of the gravity station (CRS: WGS84)
 - **xWGS84_UTM10N**: Easting of the gravity station [m] (CRS: WGS84, Zone 10N)
 - **yWGS84_UTM10N**: Northing of the gravity station [m] (CRS: WGS84, Zone 10N)
 - **zWGS84**: Height of the gravity station above the WGS84 ellipsoid [m]
 - **OG**: Original (raw) gravity meansurement [mGal]
 - **FAA**: Free-Air anomaly [mGal]
 - **SBA**: Simple Bouguer Anomaly [mGal]
 - **TTC**: Total Terrian Correction [mGal]
 - **CBA**: Complete Bouguer Anomaly [mGal]
 - **ISO**: Residual isostatic gravity anomaly [mGal]

Gravity corrections were calculated using station locations based on the NAD27 CRS and the NVD29 vertical
datum. Absoluted gravity measurements were calculated using the 1967 formula based on the IGSN71 datum. For
consistency with the SRTM DEM gravity station locations were transformed into lat,lon and UTM coordinates
based on the WGS84 CRS and the WGS84 ellipsoid.


Mesh File
-----------------------------------------------------------------------------------------------------------
File name: **GroundGravCombined_InvMesh_BaseCell_200_200_50_50kmPad.msh**

OcTree mesh file saved in UBC-GIF format.

```
4096 4096 16384                         # Number of base cells in x, y, and z directions
107357.3290 3887596.6230 410325.0000    # Top SW corner of the mesh [m]
200.000 200.000 50.000                  # Base cell (smallest cell) size in x, y, and z directions [m]
3575380                                 # Total number of cells in mesh
1 1 1 2048                              # i,j,k (indices of cell location) and cell size relative to base cell
2049 1 1 2048
1 2049 1 2048
2049 2049 1 2048
```


Active Cells File
-----------------------------------------------------------------------------------------------------------
File name: **actCellTopoInd_200_200_50.npy**

(# Cells,) Numpy boolean array which is True if the cell is active (below the topographic surface) and False
if the cell is inactive (air cell above topography).



Cell Weights File
-----------------------------------------------------------------------------------------------------------
File name: **cellWeights_Depth.npy**

(# Active Cells,) Numpy array containing a cell weights for each active cell in the octree mesh. These cell
weights are based off a depth weighting which counteracts the natural ( 1/z^2) decay of the gravity kernel
function (Li and Oldenburg, 1998).



Starting model (m0) File
-----------------------------------------------------------------------------------------------------------
File name: **m0_CLV600m_Wells_Full.npy**

(# Cells,) Numpy array containing starting density contrast values for each cell in the octree mesh.

For the 5-unit PGI, L1-norm, and L2-norm inversions a zero starting model was used. This starting model
was used for the geologically constrained 5-unit PGI and contains an approximately 600 m-thick layer of
Clear Lake Volcanics at the surface and density constrast information from 4 borehole logs.



Cell Density Contrast Bound Files
-----------------------------------------------------------------------------------------------------------
File name: **lowerBoundFull_wells.npy**

(# Cells,) Numpy array containing a lower density contrast bound for each cell in the octree mesh.


File name: **upperBoundFull_wells.npy**

(# Cells,) Numpy array containing a upper density contrast bound for each cell in the octree mesh.

For the 5-unit PGI, L1-norm, and L2-norm inversions a these bounds are defined within the inversion scripts.
These bound files were only used for the geologically constrained 5-unit PGI and reflect that addition of an
approximately 600 m-thick layer of Clear Lake Volcanics at the surface and density constrast information
from 4 borehole logs to the starting model.


Inversion Scripts
-----------------------------------------------------------------------------------------------------------
File name: **Inv_GroundGrav_Combined_bounds_L2_depthWeight_dObsNeg.py**

Python script to setup and run the field dataset L2-norm (smooth) inversion.


File name: **Inv_GroundGrav_Combined_generalBounds_L1_DW_gradTotal.py**

Python script to setup and run the field dataset L1-norm (sparse/blocky) inversion.


File name: **Inv_PGI_grav_5units_depthWieght_TikhonovReg_dObsNeg.py**

Python script to setup and run the field dataset 5-unit Petrophysically and Geologically Guided Inversion
(PGI) See Astic and Oldenburg, (2019) for details.


File name: **Inv_PGI_grav_5units_DW_dObsNeg_m0CLV_wells.py**

Python script to setup and run the field dataset 5-unit PGI with additional geologic and borehole constraints.



Recovered Model Files
-----------------------------------------------------------------------------------------------------------
File name: **rhoInv_groundGrav_Combined_ISO_bounds_L2_depthWeight_dObsNeg.npy**

(# Cells,) Numpy array containing the recovered density contrast model from the L2-norm inversion.


File name: **rhoInv_groundGrav_Combined_ISO_generalBounds_L1_DW_gradTotal.npy**

(# Cells,) Numpy array containing the recovered density contrast model from the L1-norm inversion.


File name: **rhoInv_groundGravCombined_generalBounds_PGI5_depthWeight_dObsNeg.npy**

(# Cells,) Numpy array containing the recovered density contrast model from the 5-unit PGI.


File name: **rhoInv_groundGravCombined_generalBounds_PGI5_depthWeight_dObsNeg.npy**

(# Cells,) Numpy array containing the recovered density contrast model from the 5-unit PGI with additional
geologic and borehole constraints.



Predictive Data Files
-----------------------------------------------------------------------------------------------------------
File name: **dPred_groundGrav_Combined_ISO_bounds_L2_depthWeight_dObsNeg.npy**

(# Data,) Numpy array containing the predicted data from the L2-norm inversion's recovered density contrast model.


File name: **dPred_groundGrav_Combined_ISO_generalBounds_L1_DW_gradTotal.npy**

(# Data,) Numpy array containing the predicted data from the L1-norm inversion's recovered density contrast model.


File name: **dPred_mInv_groundGravCombined_generalBounds_PGI5_depthWeight_dObsNeg.npy**

(# Data,) Numpy array containing the predicted data from the 5-unit PGI's recovered density contrast model.


File name: **dPred_mInv_groundGravCombined_PGI5_DW_dObsNeg_m0CLV600m_wells.npy**

(# Data,) Numpy array containing the predicted data from the geologically constrained 5-unit PGI's recovered
density contrast model.

Data residuals can be calculated by subtracting these prediceted data from the observed data.



References:
-----------------------------------------------------------------------------------------------------------
Astic, T., Oldenburg, D.W., 2019. A framework for petrophysically and geologically guided
geophysical inversion using a dynamic Gaussian mixture model prior. Geophysical Journal
International 219, 1989–2012. doi:10.1093/gji/ggz389.

Blakely, R.J., 1995. Potential Theory in Gravity and Magnetic Applications. Cambridge
University Press, Cambridge. doi:10.1017/CBO9780511549816.

Chapman, R.H., Bishop, C.C., 1974. Bouguer gravity map of California, Santa Rosa sheet.
Technical Report.

Langenheim, V.E., Jachens, R.C., Morin, R.L., McCabe, C.A., 2007. Preliminary Gravity
and Magnetic Data of the Lake Pillsbury Region, Northern Coast Ranges, California.
Technical Report. U.S. Geological Survey.

Godson, R.H., Plouff, D., 1988. BOUGUER Version 1.0 : a microcomputer gravity-terrain-
correction program. Technical Report. U.S. Geological Survey. doi:10.3133/ofr88644B.

Langenheim, V.E., Roberts, C.W., McCabe, C.A., McPhee, D.K., Tilden, J.E., Jachens,
R.C., 2006. Preliminary Isostatic Gravity Map of the Sonoma Volcanic Field and Vicinity,
Sonoma and Napa Counties, California. Technical Report. U.S. Geological Survey. doi:10.
3133/ofr20061056.

Li, Y., Oldenburg, D.W., 1998. 3-D inversion of gravity data. Geophysics 63, 109–119.
doi:10.1190/1.1444302.

Peacock, J.R., Earney, T.E., Mangan, M.T., Schermerhorn, W.D., Glen, J.M., Walters, M.,
Hartline, C., 2020. Geophysical characterization of the Northwest Geysers geothermal
field, California. Journal of Volcanology and Geothermal Research 399, 106882. doi:10.
1016/j.jvolgeores.2020.106882.

Smith, N., 1992. Gravity interpretation of San Pablo Bay and vicinity, in: Wright, T.L. (Ed.),
Field trip guide to Late Cenozoic geology in the North Bay region. Northern California
Geologic Society, pp. 71–80.

Youngs, L.R., Chapman, R.H., Chase, G.W., 1985. Complete Bouguer gravity and aeromag-
netic maps with geology and thermal wells and springs of the Santa Rosa-Sonoma area, Sonoma
and Napa counties, California. Technical Report. California Division of Mines and Geology.
