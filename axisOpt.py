#!/usr/bin/env python3

#### Axis Stellarator Optimization Code #####
##### combination of a suite of codes #######
###  Rogerio Jorge, March 8, 2021, v0.0  ####
from axisOptMain import initialize, optaxis, optaxiscoil, getResidue, pyoculusPoincare, runqsc, runSPEC, runVMEC, runBOOZXFORM, runREGCOIL, coils2simsgeo, CoilsBonAxis, runFIELDLINES, filesClean

### Magnetic axis to use as initial condition
NFP    =  2
p2     =  0.0
etabar =  0.8782361878175473
B2c    =  0.911877842043238
iota   =  0.3650662495270209
rc     = [ 1.0,-0.15539792101233177,0.018848828622397246,-0.002376647535978239,0.00032799264423656984,-4.6838441584813964e-05,4.056333840130882e-06 ]
zs     = [ 0.0,0.1670115095908119,-0.018820974430055693,0.002362225786275556,-0.0003434355033843257,5.131269616185478e-05,-4.611592081745671e-06 ]
## Extra inputs
qvfilename  = 'OptAxis'  # Overall filenames
ncores      = 8 # Number of cores to use
rr          = 0.18 # Radial location to send to VMEC
nIterations = 50 # Number of iterations in the optimization routine
coilSeparation = 0.22 # REGCOIL parameter (distance between plasma and coils)
targetValue = 0.006 # REGCOIL parameter (target value for normal Bfield on plasma surface)
nCoilsPerNFP = 5 # REGCOIL parameter (number of coils per NFP to cutcoilsfromregcoil)
nPCoils     = 22 # Number of Fourier coefficients to Fourier transform coils
accuracy    = 50000 # Resolution of integrals in Fourier transform
findResidue = 0 # Compute the residue with SPEC+pyoculus
pResidue    = 9 # Residue parameter iota=p/q
qResudue    = 25 # Residue parameter iota=p/q
s_guessResidue = 0.5 # Initial guess for the radial location of the residue
## Runs and codes
optAxisCoil = 0 # Optimize the axis and coils together
####################
optAxis     = 0 # Optimize the axis
nQSC        = 0 # Run quasisymmetry code
nSPEC       = 0 # Run SPEC
nVMEC       = 0 # Run VMEC
nBOOZXFORM  = 0 # Run BOOZ_XFORM
nREGCOIL    = 0 # Run REGCOIL
ncoilSims   = 0 # Convert coils to SIMSGEO
plotCoilsNB = 0 # Plot Fourier transformed coils and on-axis magnetic field
nFIELDLINES = 0 # Run FIELDLINES
ngetResidue = 0 # Compute residue from pyoculus and REGCOIL's coils
nPoincare   = 0 # Compute poincare plot from pyoculus and REGCOIL's coils
cleanFiles  = 1 # Clean all results files

## Actually run the code
stel, axis = initialize(rc,zs,NFP,etabar,B2c,p2)
if optAxis==1:     stel, axis = optaxis(stel,iota,nIterations)
if optAxisCoil==1: stel, axis = optaxiscoil(qvfilename,stel,iota,nIterations)
if nQSC==1:        runqsc(qvfilename,stel,rr)
if nSPEC==1:       runSPEC(qvfilename,findResidue, pResidue, qResudue, s_guessResidue)
if nVMEC==1:       runVMEC(qvfilename)
if nBOOZXFORM==1:  runBOOZXFORM(qvfilename)
if nREGCOIL==1:    runREGCOIL(qvfilename,coilSeparation,targetValue,nCoilsPerNFP)
if ncoilSims==1:   coils2simsgeo(qvfilename,NFP,nPCoils,accuracy)
if plotCoilsNB==1: CoilsBonAxis(axis,qvfilename,NFP)
if nFIELDLINES==1: runFIELDLINES(qvfilename)
if ngetResidue==1: getResidue(stel,axis,qvfilename)
if nPoincare==1:   pyoculusPoincare(stel,qvfilename)
if cleanFiles==1:  filesClean(qvfilename)