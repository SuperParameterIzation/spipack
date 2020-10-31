import pymuqModeling as mm
import pymuqSamplingAlgorithms as ms
import pymuqApproximationWrappers as ma
import pymuqUtilities as mu
import pymuqOptimization as mo

import matplotlib.pyplot as plt
import numpy as np

# Make sure the python path includes the DiffusionEquation model from the "Modeling" examples
import sys
sys.path.append('../../../../Modeling/CustomModPiece/python')

# Import a DiffusionEquation ModPiece from the CustomModPiece example
from CustomModPiece import DiffusionEquation

# The standard deviation of the additive Gaussian noise
noiseStd = 0.05

# Number of cells to use in Finite element discretization
numCells = 200

obsIndStart = 1
obsIndEnd = numCells
obsSkip = 20

obsMod = mm.SliceOperator(numCells+1,obsIndStart,obsIndEnd,obsSkip)

# The forward model
mod = DiffusionEquation(numCells)

# Create an exponential operator to transform log-conductivity into conductivity
rechargeVal = 10.0
recharge = mm.ConstantVector( rechargeVal*np.ones((numCells,)) )

# Combine the two models into a graph
graph = mm.WorkGraph()
graph.AddNode(mm.IdentityOperator(numCells), "Log-Conductivity")
graph.AddNode(mm.ExpOperator(numCells), "Conductivity")
graph.AddNode(recharge, "Recharge")
graph.AddNode(mod,"Forward Model")
graph.AddNode(obsMod, "Observables")
graph.AddEdge("Log-Conductivity", 0, "Conductivity", 0)
graph.AddEdge("Conductivity",0,"Forward Model",0)
graph.AddEdge("Recharge",0,"Forward Model",1)
graph.AddEdge("Forward Model",0,"Observables",0)


# Set up the Gaussian process prior on the log-conductivity
gpKern = ma.SquaredExpKernel(1, 1.0, 0.3) + ma.MaternKernel(1, 0.1, 0.2, 1.0/2.0)
gpMean = ma.ZeroMean(1,1)
prior = ma.GaussianProcess(gpMean,gpKern).Discretize(mod.xs[0:-1].reshape(1,-1))
graph.AddNode(prior.AsDensity(), "Prior")
graph.AddEdge("Log-Conductivity",0,"Prior",0)

# Generate a "true" log conductivity and generate synthetic observations
trueLogK = prior.Sample()

trueHead = mod.Evaluate([np.exp(trueLogK), rechargeVal*np.ones((numCells,))])[0]
obsHead = obsMod.Apply(trueHead)[:,0] + noiseStd*mu.RandomGenerator.GetNormal(obsMod.rows())[:,0]

# Set up the likelihood and posterior
likely = mm.Gaussian(obsHead, noiseStd*noiseStd*np.ones((obsMod.rows(),)))
graph.AddNode(likely.AsDensity(),"Likelihood")
graph.AddNode(mm.DensityProduct(2), "Posterior")

graph.AddEdge("Observables",0,"Likelihood",0)
graph.AddEdge("Prior",0,"Posterior",0)
graph.AddEdge("Likelihood",0,"Posterior",1)

graph.Visualize('ModelGraph.png')

#### MCMC Algorithm definition

# Basic options for MCMC and DILI
# NOTE: It might be cleaner in general to define these options in an external
#       YAML or XML file and read that file to construction the opts dictionary.

opts = dict()
opts['NumSamples'] = 20000 # Number of MCMC steps to take
opts['BurnIn'] = 0 # Number of steps to throw away as burn in
opts['PrintLevel'] = 3 # in {0,1,2,3} Verbosity of the output
opts['HessianType'] = 'GaussNewton'  # Type of Hessian to use in DILI.  Either "Exact" or "GaussNewton".  Note that the "Exact" Hessian is not always positive definite.
opts['Adapt Interval'] = 500  # How often the LIS should be adapted.  If negative, the LIS computed at the intial point will be used and held constant
opts['Adapt Start'] = 100 # The LIS will start being adapted after this many steps
opts['Adapt End'] = 10000 # The LIS will stop being adapted after this many steps
opts['Initial Weight'] = 100 # When the average Hessian is constructed, this represents the weight or "number of samples" given to the initial Hessian.

# Options for the LOBPCG generalized eigensolver used by DILI
eigOpts = dict()
eigOpts['NumEigs'] = 200 # Maximum number of generalized eigenvalues to compute (e.g., maximum LIS dimension)
eigOpts['RelativeTolerance'] = 1e-1 # Fraction of the largest eigenvalue used as stopping criteria on how many eigenvalues to compute
eigOpts['AbsoluteTolerance'] = -100 # Minimum allowed eigenvalue
#eigOpts['BlockSize'] = 30  # Number of eigenvalues computed simultaneously in LOBPCG.  Can be larger or smaller than NumEigs, but LOBPCG is generally more efficient if this is larger than NumEigs
#eigOpts['MaxIts'] = 30 # Maximum number of iterations taken by LOBPCG within each block
eigOpts['Verbosity'] = 3 # Controls how much information the solver prints to terminal

opts['Eigensolver Block'] = 'EigSolver' # Key in the opts dictionary where LOBPCG options are specified
opts['EigSolver'] = eigOpts

# Options for the transition kernel employed in the LIS
lisOpts = dict()
lisOpts['Method'] = 'MHKernel'
lisOpts['Proposal'] = 'PropOpts' # Key in the lisOpts dictionary where proposal options are specified
lisOpts['PropOpts.Method'] = 'MALAProposal'
lisOpts['PropOpts.StepSize'] = 0.1

opts['LIS Block'] = "LIS" # Dictionary key where the LIS options are specified
opts['LIS'] = lisOpts

# Options for the transition kernel employed in the CS
csOpts = dict()
csOpts['Method'] = 'MHKernel'
csOpts['Proposal'] = 'PropOpts' # Key in the csOpts dictionary where proposal options are specified
csOpts['PropOpts.Method'] = 'CrankNicolsonProposal'
csOpts['PropOpts.Beta'] = 0.3 # Crank-Nicolson beta
csOpts['PropOpts.PriorNode'] = 'Prior' # Name of the prior node in the graph constructed above

opts['CS Block'] =  "CS" # Dictionary key where the CS options are specified
opts['CS'] = csOpts


### Set up the sampling problem
postDens = graph.CreateModPiece("Posterior")
likely = graph.CreateModPiece("Likelihood")
prior = graph.CreateModPiece("Prior")
problem = ms.SamplingProblem(postDens)

# Construct the DILI Kernel based on the options specified above
kern = ms.DILIKernel(opts, problem)

# Construct the MCMC sampler using this transition kernel
sampler = ms.SingleChainMCMC(opts, [kern])

# Run the MCMC sampler
x0 = [trueLogK]
samps = sampler.Run(x0)

# Extract the posteior samples as a matrix and compute some posterior statistics
sampMat = samps.AsMatrix()

postMean = np.mean(sampMat,axis=1)
q05 = np.percentile(sampMat,5,axis=1)
q95 = np.percentile(sampMat,95,axis=1)

# Plot the results
fig, axs = plt.subplots(ncols=3,nrows=2, figsize=(12,8))
axs[0,0].plot(mod.xs, trueHead, label='True Head')
axs[0,0].plot(mod.xs[obsIndStart:obsIndEnd:obsSkip],obsHead,'.k',label='Observed Head')
axs[0,0].legend()
axs[0,0].set_title('Data')
axs[0,0].set_xlabel('Position, $x$')
axs[0,0].set_ylabel('Hydraulic head $h(x)$')

axs[0,1].fill_between(mod.xs[0:-1], q05, q95, alpha=0.5,label='5%-95% CI')
axs[0,1].plot(mod.xs[0:-1],postMean, label='Posterior Mean')
axs[0,1].plot(mod.xs[0:-1],trueLogK,label='Truth')
axs[0,1].legend()
axs[0,1].set_title('Posterior on log(K)')
axs[0,1].set_xlabel('Position, $x$')
axs[0,1].set_ylabel('Log-Conductivity $log(K(x))$')

axs[0,2].plot(sampMat[0,:], label='$log K_0$')
axs[0,2].plot(sampMat[100,:],label='$log K_{100}$')
axs[0,2].plot(sampMat[190,:],label='$log K_{190}$')
axs[0,2].set_title('MCMC Trace')
axs[0,2].set_xlabel('MCMC Iteration (after burnin)')
axs[0,2].set_ylabel('Log-Conductivity')

axs[1,0].scatter(sampMat[1,:],sampMat[0,:],edgeColor='None',alpha=0.02)
axs[1,0].set_xlabel('$log K_1$')
axs[1,0].set_ylabel('$log K_0$')
axs[1,1].scatter(sampMat[10,:],sampMat[0,:],edgeColor='None',alpha=0.02)
axs[1,1].set_xlabel('$log K_{10}$')
axs[1,2].scatter(sampMat[20,:],sampMat[0,:],edgeColor='None',alpha=0.02)
axs[1,2].set_xlabel('$log K_{20}$')


lisVecs = kern.LISVecs()
lisVals = kern.LISVals()

fig3 = plt.figure(figsize=(12,5))
gs = fig3.add_gridspec(1, 3)

ax1 = fig3.add_subplot(gs[0, 0])
ax2 = fig3.add_subplot(gs[0, 1:])

ax1.plot(lisVals,'.',markersize=10)
ax1.plot([0,len(lisVals)], [eigOpts['RelativeTolerance']*lisVals[0], eigOpts['RelativeTolerance']*lisVals[0]],'--k')

ax1.set_ylabel('$\lambda_i$')
ax1.set_xlabel('Index $i$')
ax1.set_title('Generalized Eigenvalues')

for i in range(lisVecs.shape[1]):
    ax2.plot(mod.xs[0:-1], lisVecs[:,i])
ax2.set_title('Generalized Eigenvectors')
ax2.set_xlabel('Position $x$')
ax2.set_ylabel('$v_i$')


plt.show()
