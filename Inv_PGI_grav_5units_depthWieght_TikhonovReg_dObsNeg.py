import sys

## Define local path of SimPEG and Discretize packages on your system
sys.path.insert(0,'/user/bin/git/simpeg/discretize')
import discretize
print(discretize.__version__)
print(discretize.__path__)

sys.path.insert(1,'/user/bin/git/simpeg/simpeg')
import SimPEG
print(SimPEG.__version__)
print(SimPEG.__path__)

print(sys.path)


## Import required packages
from discretize import TreeMesh, utils
from discretize.utils import meshutils
from discretize.utils import active_from_xyz

from SimPEG.potential_fields import gravity
from SimPEG import (
    maps,
    utils,
    simulation,
    inverse_problem,
    inversion,
    optimization,
    regularization,
    data,
    data_misfit,
    directives,
)
from SimPEG.utils import io_utils
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Reproducible science
np.random.seed(518936)

import pandas as pd
import time


## Load octree mesh
mesh = TreeMesh.read_UBC('./GroundGravCombined_InvMesh_BaseCell_200_200_50_50kmPad.msh')


## Load active cell vector
actInd = np.load('./actCellTopoInd_200_200_50.npy')
print(actInd.shape)


## Create reference model (mRef)
mRef = 1e-8 * np.zeros([mesh.nC,])
print('mRef shape=', mRef.shape)


## Create starting model (m0) by taking the active cells of mRef
m0 = mRef[actInd]
print('m0 shape=', m0.shape)
nC_act = m0.shape[0]
print(nC_act)


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
dObs = -gravData.ISO   # Flip sign of gravity data to account for SimPEG's z-positive coordinate system
dataObj_ISO = data.Data(survey, dobs=np.array(dObs), standard_deviation=wd)


## Setup simulation
simulation = gravity.simulation.Simulation3DIntegral(
    mesh=mesh, survey=survey, rhoMap=idenMap, actInd=actInd, store_sensitivities="ram")


## Define data misfit function (phi_d)
dmis = data_misfit.L2DataMisfit(simulation=simulation, data=dataObj_ISO)


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


## Define reference Gaussian Mixture Model (GMM) for PGI
gmmref = utils.WeightedGaussianMixture(
    n_components=5,  # units: Franciscan, Great Valley/Geysers Plutonic Complex, dense volcanics/greenstone, low density volcanics, melt
    mesh=mesh,  # inversion mesh
    actv=actInd,  # actv cells
    covariance_type="diag",  # diagonal covariances
)

# Initialization of gmm with fit
# fake random samples, size of the mesh, number of physical properties: 2 (density and mag.susc)
gmmref.fit(np.random.randn(nC_act).reshape(-1,1))


## Set GMM parameters manually
# Set mean density contrasts for each unit
rhoFranciscan = 0.0
rhoGreatValley = -0.1
rhoVolcanics = 0.1
rhoLowDensityVolcanics = -0.2
rhoMelt = -0.3

gmmref.means_ = np.c_[
    [rhoFranciscan],   # Franciscan Assemblage
    [rhoGreatValley],  # Great Valley/Geysers Plutonic Complex
    [rhoVolcanics],    # dense volcanics/greenstone
    [rhoLowDensityVolcanics],    # low density volcanics
    [rhoMelt]          # partial melt zone
].T


## Set covariances for each unit
# Some experimentation required see article by Thibaut Astic for guideance

# Franciscan
minFranciscan = rhoFranciscan - 0.04
maxFranciscan = rhoFranciscan + 0.02
stdFranciscan = (maxFranciscan - minFranciscan)/500
print(stdFranciscan)

# CL volcanics
minVolcanics = rhoVolcanics - 0.1
maxVolcanics = rhoVolcanics + 0.1
stdVolcanics = (maxVolcanics - minVolcanics)/100
print(stdVolcanics)

# CL low density volcanics
minLowDensityVolcanics = rhoLowDensityVolcanics - 0.1
maxLowDensityVolcanics = rhoLowDensityVolcanics + 0.1
stdLowDensityVolcanics = (maxLowDensityVolcanics - minLowDensityVolcanics)/100
print(stdLowDensityVolcanics)

# Great Valley
minGreatValley = rhoGreatValley - 0.05
maxGreatValley = rhoGreatValley + 0.05
stdGreatValley = (maxGreatValley - minGreatValley)/500
print(stdGreatValley)

# Melt
minMelt = rhoMelt - 0.2
maxMelt = rhoMelt + 0.1
stdMelt = (maxMelt - minMelt)/50
print(stdMelt)

gmmref.covariances_ = np.array(
    [[stdFranciscan],  # Franciscan density contrast
    [stdGreatValley],  # Great Valley
    [stdVolcanics],    # CL volcanics
    [stdLowDensityVolcanics],   # low density volcanics
    [stdMelt]]         # Melt
)

check = gmmref.covariances_
print(check.shape)

# After setting covariances manually: compute precision matrices and cholesky
gmmref.compute_clusters_precisions()


## Set global proportions; low-impact as long as not 0 or 1 (total=1)
weightsFranciscan = 0.8 * np.ones_like(m0)

weightsGreatValley = 0.15 * np.ones_like(m0)
weightsGreatValley[CCbelow_Neg6000Ind] = 1e-256   # Great Valley not allowed below z =-6 km

weightsVolcanics = 0.025 * np.ones_like(m0)

weightsLowDensityVolcanics = 0.025 * np.ones_like(m0)
weightsLowDensityVolcanics[CCbelow_Neg6000Ind] = 1e-256   # CL Volcanics not allowed below z =-6 km

weightsMelt = 1e-256 * np.ones_like(m0)   # Melt not allowed above z =-6 km
weightsMelt[CCbelow_Neg6000Ind] = 0.175

gmmref.weights_ = np.vstack([weightsFranciscan, weightsGreatValley, weightsVolcanics, weightsLowDensityVolcanics, weightsMelt]).T

check = np.sum(gmmref.weights_,axis=1)
print(check)
print(check.shape)
print('min(check)=', np.min(check))
print('max(check)=', np.max(check))
print('mean(check)=', np.mean(check))


## Create PGI regularization
# Depth weighting
# wr_Depth = SimPEG.utils.depth_weighting(mesh, reference_locs=topoXYZ, indActive=actInd, exponent=2.0, threshold=0.1)
wr_Depth = np.load('./cellWeights_Depth.npy')

# Dummy wire setup for single property inversion
wires = maps.Wires(("m", m0.shape[0]))

# Create joint PGI regularization with smoothness
reg = utils.make_PGI_regularization(
    gmmref=gmmref,                  # Set reference Gaussian mixture model
    mesh=mesh,                      # Inversion mesh
    wiresmap=wires,                 # Wires of physical properties
    maplist=[idenMap],              # List of used mappings
    mref=m0,                        # Set reference/starting model
    indActive=actInd,               # Active cell vector
    alpha_s=1.0,                    # Model smallness weight
    alpha_x=1.0,                    # Model smoothness weight on first x-derivative
    alpha_y=1.0,                    # Model smoothness weight on first y-derivative
    alpha_z=1.0,                    # Model smoothness weight on first z-derivative
    alpha_xx=0.0,                   # Model smoothness weight on second x-derivative
    alpha_yy=0.0,                   # Model smoothness weight on second y-derivative
    alpha_zz=0.0,                   # Model smoothness weight on second z-derivative
    cell_weights_list=[wr_Depth],   # Cell weights (Depth weighting)
)


## Create optimization to determine which numerical methods will be used to solve the inversion
opt = optimization.ProjectedGNCG(
    maxIter=50,         # Max number of inversion iterations
    lower=lowerBound,   # Lower bound for density contrast
    upper=upperBound,   # Upper bound for density contrast
    maxIterLS=20,       # Max number of iterations for each line search
    maxIterCG=100,      # Max number of CG solves per inexact Gauss-Newton iteration
    tolCG=1e-4,         # Tolerance of CG solve
)


## Setup directives for the inverse problem
# Ratio to use for each phys prop. smoothness in each direction;
# roughly the ratio of the order of magnitude of each phys. prop.
alpha0_ratio = np.r_[
    np.zeros(len(reg.objfcts[0].objfcts)),
    1e-2 * np.ones(len(reg.objfcts[1].objfcts))]
print("alpha0_ratio =", alpha0_ratio)

Alphas = directives.AlphasSmoothEstimate_ByEig(alpha0_ratio=alpha0_ratio, verbose=True)

# Initialize beta and beta/alpha_s schedule
beta = directives.BetaEstimate_ByEig(beta0_ratio=1e-2)
betaIt = directives.PGI_BetaAlphaSchedule(
    verbose=True, coolingFactor=2.0, tolerance=0.2, progress=0.2,
)

# Geophysical and petrophysical target misfits
targets = directives.MultiTargetMisfits(verbose=True,)

# Add learned mref into smooth regularization once stable
MrefInSmooth = directives.PGI_AddMrefInSmooth(wait_till_stable=True, verbose=True,)

# Update the parameters in smallness (L2-approx of PGI)
update_smallness = directives.PGI_UpdateParameters(
    update_gmm=False  # keep GMM model fixed
)

# Updat Jacobi pre-conditioner
update_Jacobi = directives.UpdatePreconditioner()

# Save the model every iteration
save_model = directives.SaveModelEveryIteration()


## Create inverse problem
invProb = inverse_problem.BaseInvProblem(dmis, reg, opt)
inv = inversion.BaseInversion(
    invProb,
    # directives: set directives
    directiveList=[
        Alphas,
        beta,
        update_smallness,
        targets,
        betaIt,
        MrefInSmooth,
        update_Jacobi,
        save_model
    ]
)


## Run the inversion
pgi_model = inv.run(m0)


## Map active cells back to full domain
density_model = actMap * pgi_model
quasi_geology_model = actMap * reg.objfcts[0].membership(reg.objfcts[0].mref)
final_petrophysical_misfit = reg.objfcts[0].objfcts[0](pgi_model, externalW=False)


## Save the results
np.save('mInv_groundGravCombined_generalBounds_PGI5_depthWeight_dObsNeg.npy', pgi_model)
np.save('rhoInv_groundGravCombined_generalBounds_PGI5_depthWeight_dObsNeg.npy', density_model)
np.save('quasiGeologyMod_groundGravCombined_generalBounds_PGI5_depthWeight_dObsNeg.npy', quasi_geology_model)
np.save('final_petrophysical_misfit_PGI5_depthWeight_dObsNeg.npy', final_petrophysical_misfit)


## Forward model predicted data from final model
time0 = time.time()
dPred = simulation.dpred(pgi_model)
time1 = time.time()
print('FwdMod time:' + str(time1-time0))

np.save('dPred_mInv_groundGravCombined_generalBounds_PGI5_depthWeight_dObsNeg.npy', dPred)
