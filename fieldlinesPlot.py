#!/opt/local/bin/python3

#myMarkerSize = 0.005
myMarkerSize = 1.0

print("Usage: "+__file__+" <fieldlines_*.h5> <Add any arguments to save a PDF>")

import sys, os
if len(sys.argv) < 2:
    print("Error! Wrong number of arguments.")
    exit(1)

makePDF = False
if len(sys.argv)>2:
    makePDF=True

if makePDF:
    import matplotlib
    matplotlib.use('PDF')

filename = sys.argv[1]

import h5py
import numpy as np
f = h5py.File(filename,'r')

raxis = f["raxis"][()]
zaxis = f["zaxis"][()]

R_lines = f["R_lines"][()]
Z_lines = f["Z_lines"][()]
print("R_lines.shape:",R_lines.shape)
M = R_lines.shape[0]
N = R_lines.shape[1]
print("Number of field lines found:",N)

import matplotlib.pyplot as plt

if makePDF:
	print("Saving PDF")
	fig = plt.figure(figsize=(14,7))
	fig.patch.set_facecolor('white')
	for j in range(N):
		data_R1 = R_lines[0:M:4,j]
		data_Z1 = Z_lines[0:M:4,j]
		mask = data_R1 != data_R1[-1]
		data_R1 = data_R1[mask]
		data_Z1 = data_Z1[mask]
		plt.plot(data_R1,data_Z1,'o',ms=myMarkerSize,markeredgecolor=None,mew=0)
		data_R2 = R_lines[2:M:4,j]
		data_Z2 = Z_lines[2:M:4,j]
		mask = data_R2 != data_R2[-1]
		data_R2 = data_R2[mask]
		data_Z2 = data_Z2[mask]
		plt.plot(data_R2,data_Z2,'o',ms=myMarkerSize,markeredgecolor=None,mew=0)
	plt.xlabel('R')
	plt.ylabel('Z')
	data_R2=data_R1
	data_Z2=data_Z1
	if np.min([np.min(data_R1),np.min(data_R2)]) != 0:
		plt.xlim((0.85*np.min([np.min(data_R1),np.min(data_R2)]),1.05*np.max([np.max(data_R1),np.max(data_R2)])))
		plt.ylim((-1.35*np.max([np.max(data_Z1),np.max(data_Z2)]),1.35*np.max([np.max(data_Z1),np.max(data_Z2)])))
	plt.gca().set_aspect('equal', adjustable='box')
	plt.savefig(filename+'.pdf', bbox_inches = 'tight', pad_inches = 0)
else:
	fig = plt.figure(figsize=(14,7))
	fig.patch.set_facecolor('white')

	numRows = 2
	numCols = 2

	plt.subplot(numRows,numCols,1)
	for j in range(N):
		data_R = R_lines[0:M:4,j]
		data_Z = Z_lines[0:M:4,j]
		mask = data_R != data_R[-1]
		data_R = data_R[mask]
		data_Z = data_Z[mask]
		plt.plot(data_R,data_Z,'o',ms=myMarkerSize,markeredgecolor=None,mew=0)
	plt.xlabel('R')
	plt.ylabel('Z')
	plt.axis('equal')
	plt.xlim((np.min(raxis),np.max(raxis)))
	plt.ylim((np.min(zaxis),np.max(zaxis)))

	plt.subplot(numRows,numCols,2)
	for j in range(N):
		data_R = R_lines[1:M:4,j]
		data_Z = Z_lines[1:M:4,j]
		mask = data_R != data_R[-1]
		data_R = data_R[mask]
		data_Z = data_Z[mask]
		print("Size of data_R:",len(data_R))
		plt.plot(data_R,data_Z,'o',ms=myMarkerSize,markeredgecolor=None,mew=0)
	plt.xlabel('R')
	plt.ylabel('Z')
	plt.axis('equal')
	plt.xlim((np.min(raxis),np.max(raxis)))
	plt.ylim((np.min(zaxis),np.max(zaxis)))

	plt.subplot(numRows,numCols,3)
	for j in range(N):
		data_R = R_lines[2:M:4,j]
		data_Z = Z_lines[2:M:4,j]
		mask = data_R != data_R[-1]
		data_R = data_R[mask]
		data_Z = data_Z[mask]
		plt.plot(data_R,data_Z,'o',ms=myMarkerSize,markeredgecolor=None,mew=0)
	plt.xlabel('R')
	plt.ylabel('Z')
	plt.axis('equal')
	plt.xlim((np.min(raxis),np.max(raxis)))
	plt.ylim((np.min(zaxis),np.max(zaxis)))

	plt.subplot(numRows,numCols,4)
	for j in range(N):
		data_R = R_lines[3:M:4,j]
		data_Z = Z_lines[3:M:4,j]
		mask = data_R != data_R[-1]
		data_R = data_R[mask]
		data_Z = data_Z[mask]
		plt.plot(data_R,data_Z,'o',ms=myMarkerSize,markeredgecolor=None,mew=0)
	plt.xlabel('R')
	plt.ylabel('Z')
	plt.axis('equal')
	plt.xlim((np.min(raxis),np.max(raxis)))
	plt.ylim((np.min(zaxis),np.max(zaxis)))

	plt.tight_layout()

	plt.figtext(0.5,0.99,os.path.abspath(filename),ha='center',va='top',fontsize=6)

	plt.show()
