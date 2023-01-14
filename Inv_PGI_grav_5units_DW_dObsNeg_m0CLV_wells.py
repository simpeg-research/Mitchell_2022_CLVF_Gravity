import sys

## Define local path of SimPEG and Discretize packages on your system
sys.path.insert(0,'/gp/mamitchell/bin/git/simpeg/discretize')
import discretize
print(discretize.__version__)
print(discretize.__path__)

sys.path.insert(1,'/gp/mamitchell/bin/git/simpeg/simpeg')
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
mRef = 1e-8 * np.ones([mesh.nC,])
print('mRef shape=', mRef.shape)

## Load starting model with Clear Lake Volcanics unit and well info
m0_CLV600m_Wells_Full = np.load('./m0_CLV600m_Wells_Full.npy')
print(m0_CLV600m_Wells_Full.shape)
print('min(m0_CLV600m_Wells_Full)=', np.nanmin(m0_CLV600m_Wells_Full))
print('max(m0_CLV600m_Wells_Full)=', np.nanmax(m0_CLV600m_Wells_Full))

# Create starting model (m0) by taking the active cells of m0_CLV600m_Wells_Full
m0 = m0_CLV600m_Wells_Full[actInd]
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


## Load lower and upper bound files which account for welllog info
lowerBoundFull_wells = np.load('./lowerBoundFull_wells.npy')
print(lowerBoundFull_wells.shape)
print('min(lowerBoundFull_wells)=', np.nanmin(lowerBoundFull_wells))
print('max(lowerBoundFull_wells)=', np.nanmax(lowerBoundFull_wells))

upperBoundFull_wells = np.load('./upperBoundFull_wells.npy')
print(upperBoundFull_wells.shape)
print('min(upperBoundFull_wells)=', np.nanmin(upperBoundFull_wells))
print('max(upperBoundFull_wells)=', np.nanmax(upperBoundFull_wells))

lowerBound = lowerBoundFull_wells[actInd]
print('lowerBound size:', lowerBound.shape)
upperBound = upperBoundFull_wells[actInd]
print('upperBound size:', upperBound.shape)

# # Define Bounds
# lowerBoundFull  = np.zeros_like(m0_CLV600m_Wells_Full)
# upperBoundFull  = np.zeros_like(m0_CLV600m_Wells_Full)
# meshCC = mesh.cell_centers
# CCbelow_Neg6000Ind = np.zeros_like(m0_CLV600m_Wells_Full, dtype=bool)
# CCbelow_Neg6000Ind[np.where(meshCC[:,2] < -6000)[0]] = True

# # Above z = -6000m
# lowerBoundFull[~CCbelow_Neg6000Ind] = -0.25
# upperBoundFull[~CCbelow_Neg6000Ind] = 0.7

# # print('Above z = -6000m')
# # lowerBoundFullCheck = np.all(m0_CLV600m_Wells_Full[~CCbelow_Neg6000Ind] > lowerBoundFull[~CCbelow_Neg6000Ind])
# # print('lowerBoundFullCheck =', lowerBoundFullCheck)
# # upperBoundFullCheck = np.all(m0_CLV600m_Wells_Full[~CCbelow_Neg6000Ind] < upperBoundFull[~CCbelow_Neg6000Ind])
# # print('upperBoundFullCheck =', upperBoundFullCheck)

# # Below z = -6000m
# lowerBoundFull[CCbelow_Neg6000Ind] = -0.5
# upperBoundFull[CCbelow_Neg6000Ind] = 0.2

# # print('Below z = -6000m')
# # lowerBoundFullCheck = np.all(m0_CLV600m_Wells_Full[CCbelow_Neg6000Ind] > lowerBoundFull[CCbelow_Neg6000Ind])
# # print('lowerBoundFullCheck =', lowerBoundFullCheck)
# # upperBoundFullCheck = np.all(m0_CLV600m_Wells_Full[CCbelow_Neg6000Ind] < upperBoundFull[CCbelow_Neg6000Ind])
# # print('upperBoundFullCheck =', upperBoundFullCheck)

# # Set bounds of borehole cells to reflect std in measurements
# # The following minimum std values are used to replace unreasonalby low std values
# # minStd_deltaRho = 0.01
# # minStd_deltaRho_single = 0.02
# Wilson1_meshInd = np.load('./Wilson1_meshInd.npy')
# Neasham1_meshInd = np.load('./Neasham1_meshInd.npy')
# MagmaWatson1_meshInd = np.load('./MagmaWatson1_meshInd.npy')
# Kettenhofen1_meshInd = np.load('./Kettenhofen1_meshInd.npy')

# Wilson1_avgDeltaRho = np.load('./Wilson1_avgDeltaRho.npy')
# Neasham1_avgDeltaRho = np.load('./Neasham1_avgDeltaRho.npy')
# MagmaWatson1_avgDeltaRho = np.load('./MagmaWatson1_avgDeltaRho.npy')
# Kettenhofen1_avgDeltaRho = np.load('./Kettenhofen1_avgDeltaRho.npy')

# Wilson1_stdDeltaRho = np.load('./Wilson1_stdDeltaRho.npy')
# Neasham1_stdDeltaRho = np.load('./Neasham1_stdDeltaRho.npy')
# MagmaWatson1_stdDeltaRho = np.load('./MagmaWatson1_stdDeltaRho.npy')
# Kettenhofen1_stdDeltaRho = np.load('./Kettenhofen1_stdDeltaRho.npy')

# # Asssign bounds for Wilson1
# for ii in range(0, len(Wilson1_meshInd)):
#     lowerBoundFull[Wilson1_meshInd[ii]] = Wilson1_avgDeltaRho[ii] - Wilson1_stdDeltaRho[ii]
#     upperBoundFull[Wilson1_meshInd[ii]] = Wilson1_avgDeltaRho[ii] + Wilson1_stdDeltaRho[ii]

# # Asssign bounds for Neasham1
# for ii in range(0, len(Neasham1_meshInd)):
#     lowerBoundFull[Neasham1_meshInd[ii]] = Neasham1_avgDeltaRho[ii] - Neasham1_stdDeltaRho[ii]
#     upperBoundFull[Neasham1_meshInd[ii]] = Neasham1_avgDeltaRho[ii] + Neasham1_stdDeltaRho[ii]

# # Asssign bounds for MagmaWatson1
# for ii in range(0, len(MagmaWatson1_meshInd)):
#     lowerBoundFull[MagmaWatson1_meshInd[ii]] = MagmaWatson1_avgDeltaRho[ii] - MagmaWatson1_stdDeltaRho[ii]
#     upperBoundFull[MagmaWatson1_meshInd[ii]] = MagmaWatson1_avgDeltaRho[ii] + MagmaWatson1_stdDeltaRho[ii]

# # Asssign bounds for Kettenhofen1
# for ii in range(0, len(Kettenhofen1_meshInd)):
#     lowerBoundFull[Kettenhofen1_meshInd[ii]] = Kettenhofen1_avgDeltaRho[ii] - Kettenhofen1_stdDeltaRho[ii]
#     upperBoundFull[Kettenhofen1_meshInd[ii]] = Kettenhofen1_avgDeltaRho[ii] + Kettenhofen1_stdDeltaRho[ii]

# lowerBoundFull[~actInd] = np.nan
# upperBoundFull[~actInd] = np.nan

# lowerBound = lowerBoundFull[actInd]
# print('lowerBound size:', lowerBound.shape)
# upperBound = upperBoundFull[actInd]
# print('upperBound size:', upperBound.shape)

lowerBoundCheck = np.all(m0_CLV600m_Wells_Full[actInd] > lowerBound)
print('lowerBoundCheck =', lowerBoundCheck)

upperBoundCheck = np.all(m0_CLV600m_Wells_Full[actInd] < upperBound)
print('upperBoundCheck =', upperBoundCheck)


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
    [rhoFranciscan],   # Franciscan density contrast
    [rhoGreatValley],  # Great Valley/Geysers Plutonic Complex
    [rhoVolcanics],    # dense volcanics/greenstone
    [rhoLowDensityVolcanics],    # low density volcanics
    [rhoMelt]          # partial melt zone
].T


## Set covariances for each unit
# Some experimentation required see article by Thibaut Astic for guideance

#Franciscan
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
meshCC = mesh.cell_centers
actCC = meshCC[actInd]
CCbelow_Neg6000Ind = np.zeros_like(m0, dtype=bool)
CCbelow_Neg6000Ind[np.where(actCC[:,2] < -6000)[0]] = True

weightsFranciscan = 0.8 * np.ones_like(m0)

weightsGreatValley = 0.15 * np.ones_like(m0)
weightsGreatValley[CCbelow_Neg6000Ind] = 1e-256 # Great Valley not allowed below z =-6 km

weightsVolcanics = 0.025 * np.ones_like(m0)

weightsLowDensityVolcanics = 0.025 * np.ones_like(m0)
weightsLowDensityVolcanics[CCbelow_Neg6000Ind] = 1e-256 # CL Volcanics not allowed below z =-6 km

weightsMelt = 1e-256 * np.ones_like(m0) # Melt not allowed above z =-6 km
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
# exponent=2.0
# wr_Depth = SimPEG.utils.depth_weighting(mesh, reference_locs=topoXYZ, indActive=actInd, exponent=exponent, threshold=0.1)
wr_Depth = np.load('./cellWeights_Depth.npy')
# reg.cell_weights = wr_Depth

# Dummy wire setup for single property inversion
wires = maps.Wires(("m", m0.shape[0]))

# create joint PGI regularization with smoothness
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
np.save('mInv_groundGravCombined_PGI5_DW_dObsNeg_m0CLV600m_wells.npy', pgi_model)
np.save('rhoInv_groundGravCombined_PGI5_DW_dObsNeg_m0CLV600m_wells.npy', density_model)
np.save('quasiGeologyMod_groundGravCombined_PGI5_DW_dObsNeg_m0CLV600m_wells.npy', quasi_geology_model)


## Forward model predicted data from final model
time0 = time.time()
dPred = simulation.dpred(pgi_model)
time1 = time.time()
print('FwdMod time:' + str(time1-time0))

np.save('dPred_mInv_groundGravCombined_PGI5_DW_dObsNeg_m0CLV600m_wells.npy', dPred)


## Save final petrophysical misfit
np.save('final_petrophysical_misfit_PGI5_DW_dObsNeg_m0CLV600m_wells.npy', final_petrophysical_misfit)
