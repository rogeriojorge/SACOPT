#!/usr/bin/env python3
# import numpy as np
# from tempfile import mkstemp
# from os import path, fdopen, remove
# from shutil import move, copymode, copyfile
# import pandas as pd
# import simsopt as so
# from scipy import interpolate
# import matplotlib.pyplot as plt
# import errno
# from joblib import Parallel, delayed
# import multiprocessing
# from scipy.io import netcdf
# from numpy.matlib import repmat
# from matplotlib import cm

def output2qsc(stel,rr,filename):
	from numpy import savetxt
	vmecExample='"input.VMECexample"'
	with open(filename, "w") as txt_file:
		txt_file.write('&quasisymmetry\n')
		txt_file.write(' general_option="single"\n')
		txt_file.write(' finite_r_option="nonlinear"\n')
		txt_file.write(' order_r_option="r3_flux_constraint"\n')
		txt_file.write(' N_phi=151\n')
		txt_file.write(' nfp = '+str(stel.nfp)+'\n')
		txt_file.write(' eta_bar = '+str(stel.etabar)+'\n')
		txt_file.write(' R0c = ')
		savetxt(txt_file, stel.rc, newline=',')
		txt_file.write('\n')
		txt_file.write(' Z0s = ')
		savetxt(txt_file, stel.zs, newline=',')
		txt_file.write('\n')
		txt_file.write(' B2c = '+str(stel.B2c)+'\n')
		txt_file.write(' B2s = '+str(stel.B2s)+'\n')
		txt_file.write(' sigma_initial = '+str(stel.sigma0)+'\n')
		txt_file.write(' I2_over_B0 = '+str(stel.I2)+'\n')
		txt_file.write(' p2 = '+str(stel.p2)+'\n')
		txt_file.write(' r = '+str(rr)+'\n')
		txt_file.write(' vmec_template_filename = '+vmecExample+'\n')
		txt_file.write(' r_singularity_high_order = T'+'\n')
		txt_file.write(' r_singularity_Newton_iterations = 20'+'\n')
		txt_file.write(' r_singularity_line_search = 10'+'\n')
		txt_file.write('/\n')

## Function to replace text in files
def replace(file_path, pattern, subst):
    #Create temp file
    from tempfile import mkstemp
    from os import fdopen, remove
    from shutil import copymode, move
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    #Copy the file permissions from the old file to the new file
    copymode(file_path, abs_path)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)
    return 0

def output2spec(qvfilename,qscfile,nmodes,findResidue):
	from scipy.io import netcdf
	from shutil import copyfile
	f = netcdf.netcdf_file(qscfile,mode='r',mmap=False)
	RBC = f.variables['RBC'][()]
	ZBS = f.variables['ZBS'][()]
	rc = f.variables['R0c'][()]
	zs = f.variables['Z0s'][()]
	nfp = f.variables['nfp'][()]
	with open("input."+qvfilename, 'r') as read_obj:
		for line in read_obj:
			if "PHIEDGE" in line:
				PHIEDGE=line.split()[2]
	rbc=[0 for i in range(len(RBC))]
	zbs=[0 for i in range(len(RBC))]
	nmodesTot=len(RBC[1])
	for count in range(len(RBC)):
		if count==0:
			val = next((index for index,value in enumerate(RBC[0]) if value != 0), None)
			rbc[count]=RBC[count][val:val+2*nmodes]
			val = next((index for index,value in enumerate(ZBS[0]) if value != 0), None)
			zbs[count]=ZBS[count][val-1:val-1+2*nmodes]
		else:
			rbc[count]=RBC[count][int((nmodesTot-1)/2-nmodes):int((nmodesTot-1)/2+nmodes)]
			zbs[count]=ZBS[count][int((nmodesTot-1)/2-nmodes):int((nmodesTot-1)/2+nmodes)]
	text="! Axis shape\n"
	text=text+"rac(0:"+str(len(rc)-1)+")="+", ".join([str(elem) for elem in rc])+"\n"
	text=text+"zas(0:"+str(len(zs)-1)+")="+", ".join([str(elem) for elem in zs])+"\n"
	copyfile("inputExample.sp", qvfilename+".sp")
	replace(qvfilename+".sp","! Axis shape",text)
	text="!----- Boundary Parameters -----\n"
	for countn in range(len(rbc)-1):
		if countn==0:
			for countm in range(nmodes):
				text=text+"rbc("+str(countm)+","+str(countn)+")= "+str(rbc[countn][countm])+", zbs("+str(countm)+","+str(countn)+")= "+str(zbs[countn][countm])+",\n"
		elif countn==len(rbc)-2:
			for countm in range(2*nmodes):
				if countm==2*nmodes-1:
					text=text+"rbc("+str(countm-nmodes)+","+str(countn)+")= "+str(rbc[countn][countm])+", zbs("+str(countm-nmodes)+","+str(countn)+")= "+str(zbs[countn][countm])+"\n"
				else:
					text=text+"rbc("+str(countm-nmodes)+","+str(countn)+")= "+str(rbc[countn][countm])+", zbs("+str(countm-nmodes)+","+str(countn)+")= "+str(zbs[countn][countm])+",\n"
		else:
			for countm in range(2*nmodes):
				text=text+"rbc("+str(countm-nmodes)+","+str(countn)+")= "+str(rbc[countn][countm])+", zbs("+str(countm-nmodes)+","+str(countn)+")= "+str(zbs[countn][countm])+",\n"
	replace(qvfilename+".sp","!----- Boundary Parameters -----",text)
	replace(qvfilename+".sp"," Nfp         =         4"," Nfp         =         "+str(nfp))
	replace(qvfilename+".sp","phiedge     =   2.000000000000000E+00","phiedge     =   "+str(PHIEDGE))
	#replace(qvfilename+".sp","pressure    =   0.000000000000000E+00","pressure    =   "+p2)
	if findResidue==1:
		replace(qvfilename+".sp","nPpts       =        500","nPpts       =        0")
		replace(qvfilename+".sp","nPtrj       =        30","nPtrj       =        0")

def output2regcoil(regcoilFile,vmecFile,nescinfilename,coilSeparation,targetValue):
	with open(regcoilFile, "w") as txt_file:
		txt_file.write('&regcoil_nml\n')
		txt_file.write(' general_option=5\n')
		txt_file.write(' nlambda=100\n')
		txt_file.write(' lambda_min = 1e-19\n')
		txt_file.write(' lambda_max = 1e-12\n')
		txt_file.write(' target_option = "max_Bnormal"\n')
		txt_file.write(' target_value = '+str(targetValue)+'\n')
		txt_file.write(' ntheta_plasma=64\n')
		txt_file.write(' ntheta_coil  =64\n')
		txt_file.write(' nzeta_plasma =64\n')
		txt_file.write(' nzeta_coil   =64\n')
		txt_file.write(' mpol_potential = 8\n')
		txt_file.write(' ntor_potential = 8\n')
		txt_file.write(' geometry_option_plasma = 2\n')
		txt_file.write(" wout_filename = '"+vmecFile+"'\n")
		txt_file.write(' geometry_option_coil=2\n')
		txt_file.write(' separation = '+str(coilSeparation)+'\n')
		txt_file.write(" nescin_filename = '"+nescinfilename+"'\n")
		txt_file.write("/\n")

def output2fieldlines(fieldlinesFile):
	with open(fieldlinesFile, "w") as txt_file:
		txt_file.write('#!fortran\n&INDATA\n')
		txt_file.write(' EXTCUR = 10000.00\n')
		txt_file.write('/\n')
		txt_file.write('&FIELDLINES_INPUT\n')
		txt_file.write(' NR = 101\n')
		txt_file.write(' NPHI = 41\n')
		txt_file.write(' NZ = 101\n')
		txt_file.write(' RMIN = 0.4\n')
		txt_file.write(' RMAX = 1.4\n')
		txt_file.write(' ZMIN = -0.4\n')
		txt_file.write(' ZMAX = 0.4\n')
		txt_file.write(' PHIMIN = 0.0\n')
		txt_file.write(' PHIMAX = 2.094\n')
		txt_file.write(' MU = 0.0\n')
		txt_file.write(" R_START = 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.1 1.11 1.12 1.13 1.14 1.15\n")
		txt_file.write(' Z_START = 64*0.00\n')
		txt_file.write(' PHI_START = 64*0.00\n')
		txt_file.write(" PHI_END =  64*2600.\n")
		txt_file.write(" NPOINC = 8\n")
		txt_file.write(" INT_TYPE = 'LSODE'\n")
		txt_file.write(" FOLLOW_TOL = 1.0E-9\n")
		txt_file.write("/\n")
		txt_file.write("&END\n")

## Fourier coefficients calculation
from numpy import pi
# cos coefficient calculation.
def a(f, n, accuracy = 50, L=pi):
	from numpy import linspace,sum,cos
	a, b = 0, 2*L
	dx = (b - a) / accuracy
	integration = 0
	x=linspace(a, b, accuracy)
	integration=sum(f(x) * cos((n * pi * x) / L))
	# for x in np.linspace(a, b, accuracy):
	# 	integration += f(x) * np.cos((n * np.pi * x) / L)
	integration *= dx
	if n==0: integration=integration/2
	return (1 / L) * integration
# sin coefficient calculation.
def b(f, n, accuracy = 50, L=pi):
	from numpy import linspace,sum,sin
	a, b = 0, 2*L
	dx = (b - a) / accuracy
	integration = 0
	x=linspace(a, b, accuracy)
	integration=sum(f(x) * sin((n * pi * x) / L))
	# for x in np.linspace(a, b, accuracy):
	# 	integration += f(x) * np.sin((n * np.pi * x) / L)
	integration *= dx
	return (1 / L) * integration
# Fourier series.   
def Sf(f, x, n = 15, accuracy=50, L=pi):
	from numpy import zeros,cos,sin,arange,size
	a0 = a(f, 0, accuracy)
	sum = zeros(size(x))
	for i in arange(1, n + 1):
		sum += ((a(f, i, accuracy) * cos((i * pi * x) / L)) + (b(f, i, accuracy) * sin((i * pi * x) / L)))
	return (a0) + sum  

def importCoils(file):
	from pandas import read_csv
	allCoils = read_csv(file)
	allCoilsValues = allCoils.values
	coilN=0
	coilPosN=0
	coilPos=[[],[]]
	for nVals in range(len(allCoilsValues)-2):
		listVals=allCoilsValues[nVals+2][0]
		vals=listVals.split()
		try:
			floatVals = [float(nVals) for nVals in vals][0:3]
			coilPos[coilN].append(floatVals)
			coilPosN=coilPosN+1
		except:
			try:
				floatVals = [float(nVals) for nVals in vals[0:3]][0:3]
				coilPos[coilN].append(floatVals)
				coilN=coilN+1
				coilPos.append([])
			except:
				break
	current=allCoilsValues[6][0].split()
	current=float(current[3])
	coilPos=coilPos[:-2]
	return coilPos, current

def getCoilsFourier(coil0,nPCoils,accuracy):
	from numpy import linspace,concatenate,zeros,sqrt,array
	from scipy import interpolate
	xArr=[i[0] for i in coil0]
	yArr=[i[1] for i in coil0]
	zArr=[i[2] for i in coil0]
	# xthetaArr=linspace(0,2*pi,len(xArr))
	# ythetaArr=linspace(0,2*pi,len(yArr))
	# zthetaArr=linspace(0,2*pi,len(zArr))
	# xf  = interpolate.interp1d(xthetaArr,xArr, kind='cubic')
	# yf  = interpolate.interp1d(ythetaArr,yArr, kind='cubic')
	# zf  = interpolate.interp1d(zthetaArr,zArr, kind='cubic')
	### Uniform arclength along the coil reparametrization
	### independent variable ivariable = sum_i[sqrt(dx_i^2+dy_i^2+dz_i^2)]
	### Renormalize ivariable - 0 to 2pi
	### xf  = interpolate.interp1d(ivariable,xArr, kind='cubic')
	L = [0 for i in range(len(xArr))]
	for itheta in range(1,len(xArr)):
		dx = xArr[itheta]-xArr[itheta-1]
		dy = yArr[itheta]-yArr[itheta-1]
		dz = zArr[itheta]-zArr[itheta-1]
		dL = sqrt(dx*dx+dy*dy+dz*dz)
		L[itheta]=L[itheta-1]+dL
	L=(1+1e-12)*array(L)*2*pi/L[-1]
	xf  = interpolate.interp1d(L,xArr, kind='cubic')
	yf  = interpolate.interp1d(L,yArr, kind='cubic')
	zf  = interpolate.interp1d(L,zArr, kind='cubic')
	coilsFourierXS=[b(xf,j,accuracy) for j in range(nPCoils)]
	coilsFourierXC=[a(xf,j,accuracy) for j in range(nPCoils)]
	coilsFourierYS=[b(yf,j,accuracy) for j in range(nPCoils)]
	coilsFourierYC=[a(yf,j,accuracy) for j in range(nPCoils)]
	coilsFourierZS=[b(zf,j,accuracy) for j in range(nPCoils)]
	coilsFourierZC=[a(zf,j,accuracy) for j in range(nPCoils)]
	return concatenate([coilsFourierXS,coilsFourierXC,coilsFourierYS,coilsFourierYC,coilsFourierZS,coilsFourierZC])

def cartesianCoils2fourier(coilPos,outputFile,nPCoils=20,accuracy=500):
	from multiprocessing import cpu_count
	from joblib import Parallel, delayed
	from numpy import asarray,transpose,savetxt
	num_cores = cpu_count()
	coilsFourier = Parallel(n_jobs=num_cores)(delayed(getCoilsFourier)(coil0,nPCoils,accuracy) for coil0 in coilPos)
	#coilsFourier = [getCoilsFourier(coil0,nPCoils,accuracy) for coil0 in coilPos]
	coilsFourier = asarray(coilsFourier)
	coilsFourier = coilsFourier.reshape(6*len(coilPos),nPCoils)
	coilsFourier=transpose(coilsFourier)
	#coilsFourier = [np.ndarray.flatten(np.asarray(coilsF)) for coilsF in coilsFourier]
	savetxt(outputFile,coilsFourier, delimiter=',')
	return coilsFourier

def getFourierCurve(outputFile,current,ppp=10):
	from numpy import loadtxt, concatenate
	from simsopt.geo.curvexyzfourier import CurveXYZFourier
	from simsopt.core.optimizable import optimizable
	coil_data = loadtxt(outputFile, delimiter=',')
	Nt_coils=len(coil_data)-1
	num_coils = int(len(coil_data[0])/6)
	coils = [optimizable(CurveXYZFourier(Nt_coils*ppp, Nt_coils)) for i in range(num_coils)]
	for ic in range(num_coils):
		dofs = coils[ic].dofs
		dofs[0][0] = coil_data[0, 6*ic + 1]
		dofs[1][0] = coil_data[0, 6*ic + 3]
		dofs[2][0] = coil_data[0, 6*ic + 5]
		for io in range(0, Nt_coils):
			dofs[0][2*io+1] = coil_data[io+1, 6*ic + 0]
			dofs[0][2*io+2] = coil_data[io+1, 6*ic + 1]
			dofs[1][2*io+1] = coil_data[io+1, 6*ic + 2]
			dofs[1][2*io+2] = coil_data[io+1, 6*ic + 3]
			dofs[2][2*io+1] = coil_data[io+1, 6*ic + 4]
			dofs[2][2*io+2] = coil_data[io+1, 6*ic + 5]
		coils[ic].set_dofs(concatenate(dofs))
	currents = [current for i in range(num_coils)]
	return (coils, currents)

def export_coils(coils,filename,currents,NFP):
	from numpy import c_,savetxt,ones,append
	with open(filename, "w") as txt_file:
		txt_file.write("periods "+str(NFP)+"\n")
		txt_file.write("begin filament\n")
		txt_file.write("mirror NIL\n")
	for count,coil in enumerate(coils):
		with open(filename, "ab") as txt_file:
			coilE=coil.gamma()
			coilE=c_[coilE, currents[count]*ones(len(coilE))]
			coilLast=coilE[-1]
			coilE=coilE[:-1, :]
			savetxt(txt_file, coilE, fmt='%.8e')
			coilLast[3]=0.0
			coilLast=append(coilLast,"1")
			coilLast=append(coilLast,"Modular")
			savetxt(txt_file, coilLast, fmt='%.10s', newline=" ")
			txt_file.write(b"\n")
	with open(filename, "ab") as txt_file:
		txt_file.write(b"end\n")

def plot_stellarator(outName, qvfilename, coils, nfp, axis, qscfile=None):
	from numpy import zeros,asarray,cos,sin,vstack,linspace,meshgrid,interp,cos,sin
	from numpy.matlib import repmat
	import matplotlib.pyplot as plt
	from scipy.io import netcdf
	from matplotlib import cm
	## Get Coils
	gamma = coils[0].gamma()
	N = gamma.shape[0]
	l = len(coils)
	data = zeros((l*(N+1), 3))
	for i in range(l):
		data[(i*(N+1)):((i+1)*(N+1)-1), :] = coils[i].gamma()
		data[((i+1)*(N+1)-1), :] = coils[i].gamma()[0, :]
	## Get Axis
	if axis is not None:
		N = axis.gamma().shape[0]
		ma_ = zeros((nfp*N+1, 3))
		ma0 = axis.gamma().copy()
		theta = 2*pi/nfp
		rotmat = asarray([
			[cos(theta), -sin(theta), 0],
			[sin(theta), cos(theta), 0],
			[0, 0, 1]]).T

		for i in range(nfp):
			ma_[(i*N):(((i+1)*N)), :] = ma0
			ma0 = ma0 @ rotmat
		ma_[-1, :] = axis.gamma()[0, :]
		data = vstack((data, ma_))
	## Plot coils and axis
	fig=plt.figure(figsize=(7,7))
	fig.patch.set_facecolor('white')
	ax = fig.gca(projection='3d')
	maxR=max(data[:,0])
	ax.scatter(data[:,0], data[:,1], data[:,2], '.-',s=0.3)
	ax.auto_scale_xyz([-maxR,maxR],[-maxR,maxR],[-maxR,maxR])
	## Plot surface
	if qscfile is not None:
		N_theta = 40
		N_phi = 150
		f = netcdf.netcdf_file(qscfile,mode='r',mmap=False)
		B0 = f.variables['B0'][()]
		r = f.variables['r'][()]
		eta_bar = f.variables['eta_bar'][()]
		mpol = f.variables['mpol'][()]
		ntor = f.variables['ntor'][()]
		RBC = f.variables['RBC'][()]
		RBS = f.variables['RBS'][()]
		ZBC = f.variables['ZBC'][()]
		ZBS = f.variables['ZBS'][()]
		theta1D = linspace(0,2*pi,N_theta)
		phi1D = linspace(0,2*pi,N_phi)
		phi2D,theta2D = meshgrid(phi1D,theta1D)
		R = zeros((N_theta,N_phi))
		z = zeros((N_theta,N_phi))
		for m in range(mpol+1):
			for jn in range(ntor*2+1):
				n = jn-ntor
				angle = m * theta2D - nfp * n * phi2D
				sinangle = sin(angle)
				cosangle = cos(angle)
				R += RBC[m,jn] * cosangle + RBS[m,jn] * sinangle
				z += ZBC[m,jn] * cosangle + ZBS[m,jn] * sinangle
		x = R * cos(phi2D)
		y = R * sin(phi2D)
		B = B0 * (1 + r * eta_bar * cos(theta2D))
		def toString(ncVar):
			temp = [c.decode('UTF-8') for c in ncVar]
			return (''.join(temp)).strip()
		order_r_option = toString(f.variables["order_r_option"][()])
		order_r_squared = (order_r_option != 'r1' and order_r_option != 'r1_compute_B2')
		if order_r_squared:
			B20 = f.variables['B20'][()]
			B2s = f.variables['B2s'][()]
			B2c = f.variables['B2c'][()]
			phi = f.variables['phi'][()]
			B20_interpolated = interp(phi1D,phi,B20,period=2*pi/nfp)
			B20_2D = repmat(B20_interpolated,N_theta,1)
			B += r * r * (B2s * sin(2*theta2D) + B2c * cos(2*theta2D) + B20_2D)
		# Rescale to lie in [0,1]:
		B_rescaled = (B - B.min()) / (B.max() - B.min())
		#ax.set_aspect('equal')
		ax.plot_surface(x, y, z, facecolors = cm.viridis(B_rescaled), rstride=1, cstride=1, antialiased=False)
	# Hide grid lines
	ax.grid(False)
	# Hide axes ticks
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_zticks([])
	plt.axis('off')
	print("Saving coils PDF")
	plt.savefig(outName+qvfilename+'.pdf', bbox_inches = 'tight', pad_inches = 0)
	#plt.show()

def silentremove(filename):
	from os import path, fdopen, remove, mkdir
	from shutil import move
	from errno import ENOENT
	if not path.isdir('Results'):
		mkdir('Results')
	try:
		#remove(filename)
		move(filename,"Results/"+filename)
	except OSError as e: # this would be "except OSError, e:" before Python 2.6
		if e.errno != ENOENT: # errno.ENOENT = no such file or directory
			raise # re-raise exception if a different error occurred

from pyoculus.problems import CartesianBfield
class SimsgeoBiotSavart(CartesianBfield):
    def __init__(self, bs, R0, Z0, Nfp=1):
        """! Set up the problem to compute the magnetic field for simsopt.geo.BiotSavart
        @param bs an instance of simsopt.geo.BiotSavart, use to compute the magnetic field
        @param R0 the magnetic axis R coordinate at phi=0 plane
        @param Z0 the magnetic axis Z coordinate at phi=0 plane
        """
        from simsopt.geo.biotsavart import BiotSavart

        super().__init__(R0, Z0, Nfp)

        if not isinstance(bs, BiotSavart):
            raise TypeError("bs should be an instance of simsopt.geo.BiotSavart")

        self._bs = bs

    def B(self, xyz, args=None):
        """! The magnetic field, being used by parent class CartesianBfield
        @param xyz array with coordinates \f$(x, y z)\f$
        @returns \f$(B_x, B_y, B_z)\f$
        """
        point = [xyz]
        self._bs.set_points(point)
        Bfield=self._bs.B()
        return Bfield[0]

    def dBdX(self, xyz, args=None):
        """! The derivative of the magnetic field, being used by parent class CartesianBfield
        @param xyz array with coordinates \f$(x, y z)\f$
        @returns \f$\partial (B_x, B_y, B_z)/\partial (x,y,z)\f$
        """
        point = [xyz]
        self._bs.set_points(point)
        dB=self._bs.dB_by_dX()
        return dB[0]