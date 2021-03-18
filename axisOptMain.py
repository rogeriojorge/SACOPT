#!/usr/bin/env python3

## Initialize near-axis stellarator
def initialize(rcInitial,zsInitial,NFP,etabarInitial,B2cInitial,p2,nphi=101,nquadrature=100):
	from numpy import concatenate
	from simsopt.geo.curverzfourier import CurveRZFourier
	from simsopt.core.optimizable import optimizable
	from qsc import Qsc
	nfourier = len(rcInitial)-1
	curveFourier = concatenate((rcInitial,zsInitial[1:]))
	axis = CurveRZFourier(nquadrature,nfourier,NFP,True)
	axis.set_dofs(curveFourier)
	stel = optimizable(Qsc(rc=rcInitial, zs=zsInitial, etabar=etabarInitial, nfp=NFP, nphi=nphi, B2c=B2cInitial, order='r2', p2=p2))
	return stel, axis
## Optimize axis
def optaxis(stel,iotaTarget,nIterations=30,nquadrature=100):
	from simsopt.core.least_squares_problem import LeastSquaresProblem
	from simsopt.solve.serial_solve import least_squares_serial_solve
	from glob import glob
	from os import remove
	from numpy import concatenate
	from simsopt.geo.curverzfourier import CurveRZFourier
	stel.all_fixed()
	stel.set_fixed('rc(1)', False) 
	stel.set_fixed('zs(1)', False)
	stel.set_fixed('rc(2)', False) 
	stel.set_fixed('zs(2)', False)
	stel.set_fixed('rc(3)', False) 
	stel.set_fixed('zs(3)', False)
	stel.set_fixed('rc(4)', False)
	stel.set_fixed('zs(4)', False)
	stel.set_fixed('rc(5)', False)
	stel.set_fixed('zs(5)', False)
	stel.set_fixed('rc(6)', False)
	stel.set_fixed('zs(6)', False)
	stel.set_fixed('etabar', False)
	stel.set_fixed('B2c', False)
	# stel.set_fixed('p2',False)
	## Form a nonlinear-least-squares optimization problem
	term = [(stel, 'iota', iotaTarget, 1e4),
			#(stel, 'p2', 0.0, 1e-1),
			(stel, 'max_elongation', 0.0, 1e-1),
			(stel, 'B20_anomaly', 0.0, 1e2),
			(stel, 'X20', 0.0, 1e-3),
			(stel, 'X2c', 0.0, 1e-3),
			(stel, 'X2s', 0.0, 1e-3),
			(stel, 'Y20', 0.0, 1e-3),
			(stel, 'Y2c', 0.0, 1e-3),
			(stel, 'Y2s', 0.0, 1e-3),
			(stel, 'Z20', 0.0, 1e-3),
			(stel, 'Z2c', 0.0, 1e-3),
			(stel, 'Z2s', 0.0, 1e-3),
			#(stel, 'DMerc_times_r2', 1e+0,1e-0), ### CHECK POSITIVE
			#(stel, 'grad_grad_B_inverse_scale_length', 0.0,1e-1)
			]
	prob = LeastSquaresProblem(term)
	## Print initial conditions
	print('Before optimization:')
	print('NFP    = ',stel.nfp)
	print('p2     = ',stel.p2)
	print('etabar = ',stel.etabar)
	print('B2c    = ',stel.B2c)
	print('iota   = ',stel.iota)
	print('rc     = [',','.join([str(elem) for elem in stel.rc]),']')
	print('zs     = [',','.join([str(elem) for elem in stel.zs]),']')
	print('DMerc  = ',stel.DMerc_times_r2)
	print('Max elongation  = ',stel.max_elongation)
	print('gradgradB inverse length: ', stel.grad_grad_B_inverse_scale_length)
	print('objective function: ', prob.objective())
	## Solve the minimization problem:
	try:
		least_squares_serial_solve(prob, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=nIterations, method='lm')
	except KeyboardInterrupt:
		print("Terminated optimization")
	for f in glob("residuals_2021*.dat"):
		remove(f)
	for f in glob("simsopt_2021*.dat"):
		remove(f)
	## Print final conditions
	print('After optimization:')
	print('NFP    = ',stel.nfp)
	print('p2     = ',stel.p2)
	print('etabar = ',stel.etabar)
	print('B2c    = ',stel.B2c)
	print('iota   = ',stel.iota)
	print('rc     = [',','.join([str(elem) for elem in stel.rc]),']')
	print('zs     = [',','.join([str(elem) for elem in stel.zs]),']')
	print('DMerc  = ',stel.DMerc_times_r2)
	print('Max elongation  = ',stel.max_elongation)
	print('gradgradB inverse length: ', stel.grad_grad_B_inverse_scale_length)
	print('objective function: ', prob.objective())
	nN=stel.iota-stel.iotaN
	if nN==0:
		print('Quasi-axisymmetric solution')
	else:
		print('Quasi-helically symmetric solution with N =',nN)
	axis = CurveRZFourier(nquadrature,len(stel.rc)-1,stel.nfp,True)
	axis.set_dofs(concatenate((stel.rc,stel.zs[1:])))
	return stel, axis
## Optimize the axis and coils together
from simsopt.core.optimizable import Optimizable
from simsopt.geo.biotsavart import BiotSavart
from numpy import concatenate,squeeze,asarray,sqrt,dot,mean,std,array
from simsopt.geo.objectives import CurveLength
from axisOptFuncs import SimsgeoBiotSavart
from pyoculus.solvers import FixedPoint
# optimization class
class objaxiscoil(Optimizable):
	def __init__(self,coils,current,axis,stel):
		self.axis      = axis
		self.nAxisFour = len(axis.get_dofs())
		self.stel      = stel
		self.nfp       = stel.nfp
		self.coils     = coils
		self.ncoils    = len(coils)
		self.nCoilFour = len(coils[0].get_dofs())
		self.current   = current
		self.currents  = [current for i in range(len(coils))]
		self.ndofs     = self.ncoils*self.nCoilFour+self.nAxisFour+4
		self.max_elongation = self.stel.max_elongation
		self._set_names()
	def get_dofs(self):
		allCoils=[coil.get_dofs() for coil in self.coils]
		allCoilsFlatten=asarray(allCoils).flatten()
		allAxis=self.axis.get_dofs()
		return concatenate((allAxis,allCoilsFlatten,[self.stel.etabar,self.stel.p2,self.stel.B2c],[self.current]))
	def set_dofs(self, dofs):
		self.axis.set_dofs(dofs[0:self.nAxisFour])
		self.stel.rc = dofs[0:int((self.nAxisFour+1)/2)]
		self.stel.zs = concatenate(([0],dofs[int((self.nAxisFour+1)/2):int(self.nAxisFour)]))
		[coil.set_dofs(dofs[self.nAxisFour+count*self.nCoilFour:self.nAxisFour+(count+1)*self.nCoilFour]) for count,coil in enumerate(self.coils)]
		self.stel.etabar = dofs[-4]
		self.stel.p2 = dofs[-3]
		self.stel.B2c = dofs[-2]
		self.stel.calculate()
		self.current=dofs[-1]
		self.currents=[dofs[-1] for i in range(self.ncoils)]
	def _set_names(self):
		names = []
		names += ['rc({})'.format(j) for j in range(self.stel.nfourier)]
		names += ['zs({})'.format(j) for j in range(1,self.stel.nfourier)]
		names += ['coil({})'.format(i)+'({})'.format(j) for i in range(self.ncoils) for j in range(self.nCoilFour)]
		names += ['etabar', 'p2', 'B2c', 'current']
		self.names = names
	def devBonAxis(self):
		bs = BiotSavart(self.coils,self.currents)
		bs.set_points(self.axis.gamma())
		Bfield=bs.B()
		devBstrength=[sqrt(dot(Bfieldi,Bfieldi))-self.stel.B0 for Bfieldi in Bfield]
		return array(devBstrength)
	def stdBonAxis(self):
		bs = BiotSavart(self.coils,self.currents)
		bs.set_points(self.axis.gamma())
		Bfield=bs.B()
		devBstrength=[sqrt(dot(Bfieldi,Bfieldi))-self.stel.B0 for Bfieldi in Bfield]
		return std(devBstrength)
	def BonAxis(self):
		bs = BiotSavart(self.coils,self.currents)
		bs.set_points(self.axis.gamma())
		Bfield=bs.B()
		Bstrength=[sqrt(dot(Bfieldi,Bfieldi)) for Bfieldi in Bfield]
		return Bstrength
	def coilOptLength(self):
		coil_lengths = sum([CurveLength(coil).J() for coil in self.coils])
		return coil_lengths
	def iota(self):
		return self.stel.iota
	def residue(self,guess=1.6,pp=3,qq=8,sbegin=1.58,send=1.62):
		bs = BiotSavart(self.coils,self.currents)
		sbsp = SimsgeoBiotSavart(bs, R0=sum(self.stel.rc), Z0=0, Nfp=self.nfp)
		fp = FixedPoint(sbsp, {"Z":0.0})
		output=fp.compute(guess=guess, pp=pp, qq=qq, sbegin=sbegin, send=send)
		try:
			residue=output.GreenesResidue
		except:
			residue=0
		return residue
# optimization function
def optaxiscoil(qvfilename,stel,iotaTarget,nIterations=30,nquadrature=150):
	from simsopt.core.least_squares_problem import LeastSquaresProblem
	from simsopt.solve.serial_solve import least_squares_serial_solve
	from simsopt.geo.curverzfourier import CurveRZFourier
	from axisOptFuncs import getFourierCurve, plot_stellarator, export_coils
	from numpy import sqrt, dot
	from os import remove
	from glob import glob
	import matplotlib.pyplot as plt
	outputFile  = qvfilename+"_coil_coeffs.dat"
	current=3.11049660e+05
	coils, _ = getFourierCurve(outputFile,current)
	axis = CurveRZFourier(nquadrature,len(stel.rc)-1,stel.nfp,True)
	axis.set_dofs(concatenate((stel.rc,stel.zs[1:])))
	obj=objaxiscoil(coils,current,axis,stel)
	term = [(obj.iota, iotaTarget, 1e4),
			(obj.devBonAxis, 0, 1e0),
			(obj, 'max_elongation', 0, 1e1),
			(obj.stdBonAxis, 0, 1e5),
			(obj.residue, 0, 1e5)
			]
	prob = LeastSquaresProblem(term)
	obj.all_fixed()
	# obj.set_fixed('rc(1)', False) 
	# obj.set_fixed('zs(1)', False) 
	# obj.set_fixed('rc(2)', False) 
	# obj.set_fixed('zs(2)', False) 
	obj.set_fixed('rc(3)', False) 
	obj.set_fixed('zs(3)', False) 
	obj.set_fixed('rc(4)', False) 
	obj.set_fixed('zs(4)', False) 
	obj.set_fixed('rc(5)', False) 
	obj.set_fixed('zs(5)', False) 
	#obj.set_fixed('etabar', False)
	#obj.set_fixed('current', False)
	#obj.fixed[obj.nAxisFour:obj.nAxisFour+(obj.ncoils+1)*obj.nCoilFour]=False
	################################################
	rcInit=obj.stel.rc
	zsInit=obj.stel.zs
	plt.figure()
	plt.plot(obj.BonAxis())
	plt.xlabel('phi')
	plt.ylabel('B')
	################  OPTIMIZE  ####################
	least_squares_serial_solve(prob,max_nfev=nIterations)#, method='lm')
	################################################
	for f in glob("residuals_2021*.dat"):
		remove(f)
	for f in glob("simsopt_2021*.dat"):
		remove(f)
	################################################
	rcFinal=obj.stel.rc
	zsFinal=obj.stel.zs
	rcDelta=[100*(rcf-rcInit[count])/rcInit[count] for count, rcf in enumerate(rcFinal)]
	zsDelta=[100*(zsf-zsInit[count])/zsInit[count] for count, zsf in enumerate(zsFinal)]
	print('Initial rc              = [',','.join([str(elem) for elem in rcInit]),']')
	print('Final rc                = [',','.join([str(elem) for elem in rcFinal]),']')
	print('Initial zs              = [',','.join([str(elem) for elem in zsInit]),']')
	print('Final zs                = [',','.join([str(elem) for elem in zsFinal]),']')
	print('Percentage change in rc = [',','.join([str(elem) for elem in rcDelta]),']')
	print('Percentage change in zs = [',','.join([str(elem) for elem in zsDelta]),']')
	print('Target iota = ',iotaTarget,', final iota = ',obj.stel.iota)
	################################################
	print("Plot on-axis B from these coils")
	plt.plot(obj.BonAxis())
	plt.savefig("optimizedBonaxis_"+qvfilename+'.pdf', bbox_inches = 'tight', pad_inches = 0)
	################################################
	print("Plot final coils and axis")
	plot_stellarator("coils_axis_Optimized_", qvfilename, obj.coils, obj.nfp, obj.axis)
	filename = "coils."+qvfilename+"Optimized"
	export_coils(obj.coils,filename,obj.currents,obj.nfp)
	return obj.stel, obj.axis
## Compute residue with Pyoculus
def getResidue(stel,axis,qvfilename,guess=1.6,pp=3,qq=8,sbegin=1.58,send=1.62):
	from axisOptFuncs import importCoils, getFourierCurve
	from pyoculus.solvers import FixedPoint
	_, current  = importCoils("coils."+qvfilename)
	outputFile  = qvfilename+"_coil_coeffs.dat"
	coils, currents = getFourierCurve(outputFile,current)
	bs = BiotSavart(coils,currents)
	sbsp = SimsgeoBiotSavart(bs, R0=sum(stel.rc), Z0=0, Nfp=stel.nfp)
	fp = FixedPoint(sbsp, {"Z":0.0})
	#sbegin=1.01*sum(stel.rc)
	output=fp.compute(guess=guess, pp=pp, qq=qq, sbegin=sbegin, send=send)
	residue=output.GreenesResidue
	print(residue)
## Poincare plot with Pyoculus
def pyoculusPoincare(stel,qvfilename,Rbegin=1.55,Rend=1.62,nPpts=50,nPtrj=5):
	from axisOptFuncs import importCoils, getFourierCurve
	from pyoculus.solvers import PoincarePlot
	import matplotlib.pyplot as plt
	_, current  = importCoils("coils."+qvfilename)
	outputFile  = qvfilename+"_coil_coeffs.dat"
	coils, currents = getFourierCurve(outputFile,current)
	bs = BiotSavart(coils,currents)
	sbsp = SimsgeoBiotSavart(bs, R0=sum(stel.rc), Z0=0, Nfp=stel.nfp)
	params = dict()
	Rbegin=1.01*sum(stel.rc)
	params["Rbegin"] = Rbegin
	params["Rend"] = Rend
	params["nPpts"] = nPpts
	params["nPtrj"] = nPtrj
	p = PoincarePlot(sbsp, params)
	poincare_output = p.compute()
	p.plot(s=1)
	plt.savefig("poincarePyoculus_"+qvfilename+'.pdf', bbox_inches = 'tight', pad_inches = 0)
	iota = p.compute_iota()
	p.plot_iota()
	plt.savefig("iotaPyoculus_"+qvfilename+'.pdf', bbox_inches = 'tight', pad_inches = 0)
	plt.show()
## Run quasisymmetry code
def runqsc(qvfilename,stel,rr):
	from axisOptFuncs import output2qsc
	from subprocess import run
	print("Output to quasisymmetry code")
	filename="quasisymmetry_in."+qvfilename
	output2qsc(stel,rr,filename)
	print("Run quasisymmetry code")
	run(["./quasisymmetry", filename])
	print("Plot quasisymmetry result")
	import quasisymmetryPlotSingle
	quasisymmetryPlotSingle.main("quasisymmetry_out."+qvfilename+".nc")
## Run SPEC
def runSPEC(qvfilename, findResidue=0, p=2, q=5, s_guess=0.5, RunSPEC=1, nmodes=25):
	from axisOptFuncs import output2spec
	from simsopt.mhd.spec import Spec
	from simsopt.mhd.spec import Residue
	print("Output to SPEC")
	qscfile="quasisymmetry_out."+qvfilename+".nc"
	output2spec(qvfilename,qscfile,nmodes,findResidue)
	if RunSPEC==1:
		print("Run SPEC")
		spec = Spec(qvfilename+".sp")
		spec.run()
	if findResidue==1:
		residue1 = Residue(spec, p, q, s_guess=s_guess)
		r1 = residue1.J()
		print("Residue = ",r1)
	print("Plot SPEC result")
	import specPlot
	specPlot.main('spec00000.sp.h5',qvfilename)
## Run VMEC
def runVMEC(qvfilename):
	from axisOptFuncs import replace
	import subprocess
	from subprocess import run
	import vmecPlot2
	print("Run VMEC")
	replace("input."+qvfilename,'NTOR = 0075','NTOR = 0012')
	#from simsopt.mhd.vmec import Vmec
	#v = Vmec("input."+qvfilename)
	bashCommand = "./xvmec2000 input."+qvfilename
	run(bashCommand.split())
	print("Plot VMEC result")
	vmecPlot2.main("wout_"+qvfilename+".nc")
## Run booz_xform
def runBOOZXFORM(qvfilename):
	from subprocess import run
	import boozPlot
	print("Run BOOZ_XFORM")
	#from simsopt.mhd.vmec import Vmec
	#v = Vmec("input."+qvfilename)
	#from simsopt.mhd.boozer import Boozer
	#b = Boozer(v, mpol=32, ntor=16)
	with open("in_booz."+qvfilename, "w") as txt_file:
		txt_file.write('120 35\n')
		txt_file.write("'wout_"+qvfilename+"'\n")
		txt_file.write('5 10 20 30 40 50 60 70 80 90 99\n')
	bashCommand = "./xbooz_xform in_booz."+qvfilename
	run(bashCommand.split())
	print("Plot BOOZ_XFORM result")
	boozPlot.main("boozmn_wout_"+qvfilename+".nc",qvfilename)
## Output to REGCOIL
def runREGCOIL(qvfilename,coilSeparation,targetValue,nCoilsPerNFP):
	from axisOptFuncs import output2regcoil
	from subprocess import run
	print("Output to REGCOIL")
	nescinfilename = qvfilename+"_nescin.out"
	output2regcoil("regcoil_in."+qvfilename,"wout_"+qvfilename+".nc",nescinfilename,coilSeparation,targetValue)
	print("Run REGCOIL")
	bashCommand = "./regcoil regcoil_in."+qvfilename
	run(bashCommand.split())
	print("Cut coils from regcoil")
	bashCommand = "./cutCoilsFromRegcoil regcoil_out."+qvfilename+".nc "+nescinfilename+" "+str(nCoilsPerNFP)+" 0 -1"
	run(bashCommand.split())
## Convert resulting coils to simsgeo curves
def coils2simsgeo(qvfilename,NFP,nPCoils,accuracy):
	print("Convert resulting coils to SIMSGEO curves")
	from axisOptFuncs import importCoils, cartesianCoils2fourier, getFourierCurve, export_coils
	coilsCartesian, current  = importCoils("coils."+qvfilename)
	outputFile      = qvfilename+"_coil_coeffs.dat"
	cartesianCoils2fourier(coilsCartesian,outputFile,nPCoils,accuracy)
	coils, currents = getFourierCurve(outputFile,current)
	filename = "coils."+qvfilename+"2"
	export_coils(coils,filename,currents,NFP)
	return coils, currents
## Look at field on-axis from these coils
def CoilsBonAxis(axis,qvfilename,NFP):
	from axisOptFuncs import importCoils, getFourierCurve, plot_stellarator
	from simsopt.geo.biotsavart import BiotSavart
	from numpy import sqrt, dot
	import matplotlib.pyplot as plt
	_, current  = importCoils("coils."+qvfilename)
	outputFile  = qvfilename+"_coil_coeffs.dat"
	coils, currents = getFourierCurve(outputFile,current)
	print("Look at resulting coils")
	plot_stellarator("coils_FOURIER_", qvfilename, coils, NFP, axis, "quasisymmetry_out."+qvfilename+".nc")
	print("Plot on-axis B from these coils")
	bs = BiotSavart(coils,currents)
	bs.set_points(axis.gamma())
	Bfield=bs.B()
	Bstrength = [sqrt(dot(Bfieldi,Bfieldi)) for Bfieldi in Bfield]
	plt.figure()
	plt.plot(Bstrength)
	plt.xlabel('phi')
	plt.ylabel('B')
	plt.savefig("Bonaxis"+qvfilename+'.pdf', bbox_inches = 'tight', pad_inches = 0)
	#plt.show()
## Run FIELDLINES
def runFIELDLINES(qvfilename):
	from axisOptFuncs import output2fieldlines
	from subprocess import run
	print("Output to FIELDLINES")
	output2fieldlines("input.wout_"+qvfilename)
	print("Run FIELDLINES")
	bashCommand = "./xfieldlines -vmec wout_OptAxis -coils coils.OptAxis -vac"
	run(bashCommand.split())
	print("Plot FIELDLINES")
	bashCommand = "./fieldlinesPlot.py fieldlines_wout_OptAxis.h5 save"
	run(bashCommand.split())
	bashCommand = "./fieldlinesPlot.py fieldlines_wout_OptAxis.h5 "
	run(bashCommand.split())
## Remove spurious files
def filesClean(qvfilename):
	from axisOptFuncs import silentremove
	print("Cleaning files")
	silentremove("in_booz."+qvfilename)
	silentremove("boozmn_wout_"+qvfilename+".nc")
	silentremove("input."+qvfilename)
	silentremove("jxbout_"+qvfilename+".nc")
	silentremove("mercier."+qvfilename)
	silentremove("parvmecinfo.txt")
	silentremove("quasisymmetry_in."+qvfilename)
	silentremove("quasisymmetry_out."+qvfilename+".nc.pdf")
	silentremove("quasisymmetry_out."+qvfilename+".nc")
	silentremove("threed1."+qvfilename)
	silentremove("timings.txt")
	silentremove("wout_"+qvfilename+".nc")
	silentremove("coils."+qvfilename)
	silentremove("coils."+qvfilename+"2")
	silentremove("regcoil_out."+qvfilename+".nc")
	silentremove("regcoil_in."+qvfilename)
	silentremove(qvfilename+"_nescin.out")
	silentremove("coils_REGCOIL_"+qvfilename+".pdf")
	silentremove("Bonaxis"+qvfilename+".pdf")
	silentremove("coils_FOURIER_"+qvfilename+".pdf")
	silentremove(qvfilename+"_coil_coeffs.dat")
	silentremove("wout_"+qvfilename+".ncVMEC3Dplot.pdf")
	silentremove("wout_"+qvfilename+".ncVMECsurfaces.pdf")
	silentremove("wout_"+qvfilename+".ncVMECparams.pdf")
	silentremove("fieldlines_wout_"+qvfilename+".h5")
	silentremove("fieldlines_wout_"+qvfilename+".h5.pdf")
	silentremove("coils_axis_Optimized_"+qvfilename+".pdf")
	silentremove("coils."+qvfilename+"Optimized")
	silentremove("optimizedBonaxis_"+qvfilename+".pdf")
	silentremove("iotaPyoculus_"+qvfilename+".pdf")
	silentremove("poincarePyoculus_"+qvfilename+".pdf")
	silentremove("quasisymmetry_out."+qvfilename+".nc_surfs.pdf")
	silentremove("quasisymmetry_out."+qvfilename+".nc_params.pdf")
	silentremove("SPECplot"+qvfilename+".pdf")
	silentremove("spec00000.sp")
	silentremove("spec00000.sp.end")
	silentremove("spec00000.sp.h5")
	silentremove(qvfilename+".sp")
	silentremove(qvfilename+".sp.end")
	silentremove(qvfilename+".sp.h5")
	silentremove("input.wout_"+qvfilename)
	silentremove(qvfilename+"_coil_coeffs2.dat")
	silentremove("BoozxformPlot1"+qvfilename+".pdf")
	silentremove("BoozxformPlot2"+qvfilename+".pdf")
	silentremove("BoozxformPlot3"+qvfilename+".pdf")