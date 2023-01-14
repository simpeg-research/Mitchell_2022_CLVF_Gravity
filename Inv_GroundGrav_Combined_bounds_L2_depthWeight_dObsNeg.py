import sys

## Define local path of SimPEG and Discretize packages on your system
sys.path.insert(0,'/gp/mamitchell/bin/git/simpeg/discretize')
import discretize
print(discretize.__version__)
print(discretize.__path__)

sys.path.insert(1,'/gp/mamitchell/bin/git/simpeg/Mira_simpeg/simpeg')
import SimPEG
print(SimPEG.__version__)
print(SimPEG.__path__)

print(sys.path)


## Import required packages
from discretize import TreeMesh, utils
from discretize.utils import meshutils
from discretize.utils import active_from_xyz

from SimPEG import data, data_misfit, directives, maps, inverse_problem, optimization, inversion, regularization
from SimPEG.potential_fields import gravity

from scipy.interpolate import griddata
import pandas as pd

import pyproj
import time

import numpy as np


## Load octree mesh
mesh = TreeMesh.read_UBC('./GroundGravCombined_InvMesh_BaseCell_200_200_50_50kmPad.msh')


## Load active cell vector
actInd = np.load('./actCellTopoInd_200_200_50.npy')
print(actInd.shape)


## Create reference model (mRef)
mRef = np.zeros([mesh.nC,])
print('mRef shape=', mRef.shape)


## Create starting model (m0) by taking the active cells of mRef
m0 = mRef[actInd]
print('m0 shape=', m0.shape)
nC_act = m0.shape[0]
print(nC_act)


## Set lower and upper bounds for recovered density contrasts
lowerBound  = np.zeros_like(m0)
upperBound  = np.zeros_like(m0)

# Get cell centers from mesh
meshCC = mesh.cell_centers
actCC = meshCC[actInd]

# Find cell centers below z=-6km
CCbelow_Neg6000Ind = np.zeros_like(m0, dtype=bool)
CCbelow_Neg6000Ind[np.where(actCC[:,2] < -6000)[0]] = True

# Bounds for cells above z=-6km
lowerBound[~CCbelow_Neg6000Ind] = -0.25
upperBound[~CCbelow_Neg6000Ind] = 0.2

print('Above z = -6000m')
lowerBoundCheck = np.all(m0[~CCbelow_Neg6000Ind] > lowerBound[~CCbelow_Neg6000Ind])
print('lowerBoundCheck =', lowerBoundCheck)
upperBoundCheck = np.all(m0[~CCbelow_Neg6000Ind] < upperBound[~CCbelow_Neg6000Ind])
print('upperBoundCheck =', upperBoundCheck)

# Bounds for cells below z=-6km
lowerBound[CCbelow_Neg6000Ind] = -0.5
upperBound[CCbelow_Neg6000Ind] = 0.1

print('Below z = -6000m')
lowerBoundCheck = np.all(m0[CCbelow_Neg6000Ind] > lowerBound[CCbelow_Neg6000Ind])
print('lowerBoundCheck =', lowerBoundCheck)
upperBoundCheck = np.all(m0[CCbelow_Neg6000Ind] < upperBound[CCbelow_Neg6000Ind])
print('upperBoundCheck =', upperBoundCheck)


lowerBoundCheck = np.all(m0 > lowerBound)
print('Full lowerBoundCheck =', lowerBoundCheck)

upperBoundCheck = np.all(m0 < upperBound)
print('Full upperBoundCheck =', upperBoundCheck)


## Define Mappings
# Create active map to go from reduce set to full
actMap = maps.InjectActiveCells(mesh, actInd, np.nan)
idenMap = maps.IdentityMap(nP=nC_act)


## Load gravity data
gravData = pd.read_pickle('./groundGrav_Combined_zEllipsoid_Full.pkl')


## Create a Gravity Survey
# Create array using xyz locations of gravity stations
rxLoc = np.c_[gravData.xWGS84_UTM10N, gravData.yWGS84_UTM10N, gravData.zWGS84]

# Create a list of receivers from station point locations
rxList = gravity.receivers.Point(rxLoc)

# Compute the source field
srcField = gravity.sources.SourceField(receiver_list=[rxList])

# Create survey object containing rx and src information
survey = gravity.survey.Survey(srcField)


## Define a uniform uncertainty of 0.25 mGal for all gravity data
std = 0.25  # mGal
wd = np.ones_like(gravData.ISO) * std


## Create data object and assign data to the survey
dObs = -gravData.ISO # Flip sign of gravity data to account for SimPEG's z-positive coordinate system
dataObj_ISO = data.Data(survey, dobs=np.array(dObs), standard_deviation=wd)


## Setup simulation
simulation = gravity.simulation.Simulation3DIntegral(
    mesh=mesh, survey=survey, rhoMap=idenMap, actInd=actInd, store_sensitivities="ram")


## Create a regularization function
reg = regularization.Sparse(mesh, indActive=actInd, mapping=idenMap, alpha_s=1.0, alpha_x=1.0, alpha_y=1.0, alpha_z=1.0)

# Set the reference model
reg.mref = m0

# Depth weighting
# exponent=2.0
# wr_Depth = SimPEG.utils.depth_weighting(mesh, reference_locs=topoXYZ, indActive=actInd, exponent=exponent, threshold=0.1)
wr_Depth = np.load('./cellWeights_Depth.npy')
reg.cell_weights = wr_Depth


## Create optimization to determine which numerical methods will be used to solve the inversion
opt = optimization.ProjectedGNCG(maxIter=50, lower=lowerBound, upper=upperBound, maxIterLS=20, maxIterCG=100, tolCG=1e-4)


## Define data misfit function (phi_d)
dmis = data_misfit.L2DataMisfit(simulation=simulation, data=dataObj_ISO)


## Setup directive for the inverse problem
# Specify how the initial beta is found
betaest = directives.BetaEstimate_ByEig(beta0_ratio=2)

# Iteratively Reweighted Least Squares: try to fit as much as possible of the signal
IRLS = directives.Update_IRLS(f_min_change=1e-3, minGNiter=1, beta_tol=1e-1, max_irls_iterations=5)

# Updat Jacobi pre-conditioner
update_Jacobi = directives.UpdatePreconditioner()

# Save inversion output from each iteration to a log file
save = directives.SaveOutputEveryIteration()

# Save the model every iteration
save_model = directives.SaveModelEveryIteration()

## Create inverse problem
invProb = inverse_problem.BaseInvProblem(dmis, reg, opt)
inv = inversion.BaseInversion(
    invProb,
    # directives: set directives
    directiveList=[
        betaest,
        IRLS,
        update_Jacobi,
        save,
        save_model
    ]
)


## Run the inversion
time0 = time.time()
mInv = inv.run(m0)
time1 = time.time()
print('Inv time:' + str(time1-time0))


## Map active cells back to full domain
rhoInv = actMap*mInv


## Save the results
np.save('mInv_groundGrav_Combined_ISO_bounds_L2_depthWeight_dObsNeg.npy', mInv)
np.save('rhoInv_groundGrav_Combined_ISO_bounds_L2_depthWeight_dObsNeg.npy', rhoInv)


## Forward model predicted data from final model
time0 = time.time()
dPred = simulation.dpred(mInv)
time1 = time.time()
print('FwdMod time:' + str(time1-time0))

np.save('dPred_groundGrav_Combined_ISO_bounds_L2_depthWeight_dObsNeg.npy', dPred)
