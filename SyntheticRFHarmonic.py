# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Author: Tolulope Olugboji
# Date:  June 30, 2014
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# Objective:   Class for managing workflow involved with generate synthetic receiver function response.
#              object methods: loadSACData, buildModel, genSynthSACData, genRFResponse
#
#              class requires following kernel codes :
#						1. anirecSynth.f   ---- For generating synthetic SAC data giving a particular velocity model (iso/aniso - tropic)
#                       2. RecFuncHarmonic.cc  --- For computing synthetic response giving a particular SAC manifest ( real or synthetic )
#                       3. jp2tab - utility code to translate velocity model from jpark code format to table for numpy array load
#
#

import subprocess
import numpy as np
import time, sys, os, shutil

import matplotlib.pyplot as plt
from netCDF4 import Dataset

import obspyWaveFormCrawler as crwl
import runRFmodulesJaps as runRF
from matplotlib.ticker import AutoMinorLocator


from matplotlib import rcParams
rcParams['axes.labelsize'] = 14
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14
rcParams['legend.fontsize'] = 14

reload(runRF)

class SyntheticRFHarmonic():
	""" Python class to manage workflow for generating synthetic RF response given vel model """
	CODEDIR = "/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/Code/bin/"

	#KERNELSYNTH = "/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/SACFiles/anirec_synth_new"
	KERNELSYNTH =  CODEDIR + "anirec_synth_Tweak"
	KERENELRECFUNC = CODEDIR + "RecFuncHarmonic"
	JP2TAB = CODEDIR + "jp2tab"
	
	
	PROJDIR = "/Users/tmo22/Documents/JapaneseProject/"
	BUFDIR = PROJDIR + "ProcData/buffer/"
	SYNTHPARAMFILE =  PROJDIR + "ProcData/buffer/paramFile.txt"
	
	SYNTHRESULTDIR = PROJDIR + "RFDataStacks/synthStacks/"
	GRDRESULTDIR = PROJDIR + "RFDataStacks/gridStacks/"
	DATARESULTDIR = PROJDIR + "RFDataStacks/OBS/"
	
	SYNTHIMAGEDIR = PROJDIR + "Images/RFImages/synthetics/Solutions/synth_"
	INVIMAGEDIR = PROJDIR + "Images/RFImages/synthetics/GridInversions/"
	
	MODELLEDSTACK_MIG = "MigrateHarmonic.Modelled."
	MODELLEDSTACK_NOMIG = ".Modelled."
	
	M_MEAN = MODELLEDSTACK_MIG + "mean.xyz"
	M_PDTHETA = MODELLEDSTACK_MIG + "devMax.xyz"
	M_NDTHETA = MODELLEDSTACK_MIG + "devMin.xyz"
	
	M_MEAN_N = MODELLEDSTACK_NOMIG + "mean.xyz"
	M_PDTHETA_N = MODELLEDSTACK_NOMIG + "devMax.xyz"
	M_NDTHETA_N = MODELLEDSTACK_NOMIG + "devMin.xyz"
	
	
	PLTDRIVERS = PROJDIR + "RFScripts/"
	
	# No velocity models are stored internal to object, only references to the local files are stored
	fvelModel = ""				# velocity model linked to this particular object.
	fVelModelWithLAB = ""		# velocity model grafted with LAB reduction
	
	fSynthManifest = ""
	fDataManifest = ""
	
	manifestHead = ""
	
	
	# Internal arrays - Base Model and LABtweak model
	synthHarmonicRFArray = None
	synthHarmonicRFArrayLAB = None
	
	dataHarmonicRFArrayM = None
	dataHarmonicRFArrayP = None
	dataHarmonicRFArrayN = None

	velModelArray = None
	velModelArrayLAB = None
	
	nEvents = 0
		
	def __init__(self, msg):
		self.msg = msg
		print "Driver to automate workflow for generating RF synthetics ", self.msg
	
	def makeKernelSynth(self):
		codeDir = "/Users/tmo22/Documents/OlugbojiPark2013/Seismology/JeffRecFunc/"
		
		sourceFile = codeDir + "JJParkSeismoSrc/anirec_synth_Tweak.f"
		incFileA = codeDir + "plotit/plotlib.a"
		incFileB = codeDir + "plotit/eislib.a"
		incFileC = codeDir + "plotit/jlib.a"
		
		
		# gfortran -o  $CodeDir/bin/anirec_synth $CodeDir/JJParkSeismoSrc/anirec_synth.f $CodeDir/plotit/plotlib.a $CodeDir/plotit/eislib.a $CodeDir/plotit/jlib.a
		
		makeArgument = ['gfortran', '-o', self.KERNELSYNTH, sourceFile, incFileA, incFileB, incFileC]
		status = subprocess.call(makeArgument)

		if ( status == 0):
			print "*** Succesful Make! ",  status, "No Errors"
			if( status == 0):
				print " Successful updated anistoropic modeling code",  status
		else:
			print "Make Fail.. Revisit code and recompile !!!"



	# Code to control work flow for synthetic computation and visualization of computations
	def runWorkFlowSingle(self, prepBuf=True, getDataRF=True, buildSeismograms=True, buildRecFunc=True, useLAB=True):
		
		self.prepBuffer(copyManifest=True, useLAB=False)
		self.buildSynthSeismograms()
		self.buildRecFuncResponse(isMig=False, updateData=getDataRF, updateSynth=True, useLAB=False)
		
		## prepBuffer redefines model before building synthetic sesimogram: LAB model
		## Also since synthetic seismograms have changed, I need to redefine output RF data using new fileName
		self.prepBuffer(copyManifest=False, useLAB=True)
		self.buildSynthSeismograms()
		self.buildRecFuncResponse(isMig=False, updateData=getDataRF, updateSynth=True, useLAB=True)
	
	def runWorkFlowStack(self, nPts = 5, useH=50, useVsDrop=0.1, useDhGrad=0, lTime=2, rTime=9):
		hList = np.linspace(useH-10, useH+10, nPts)
		
		if useDhGrad > 0:
			vsDrop = useVsDrop
			xList = np.linspace(0,useDhGrad, nPts)
			# Use text buffer to store grid and then during save, transform to grd file
			tempXYZOut = open(self.BUFDIR+self.manifestHead+"_HDhGradGrd.txt", 'w')
			saveDH = True
		if useDhGrad == 0:
			xList = np.linspace(0.02, 0.15, nPts)
			dhGrad = 0
			# Use text buffer to store grid and then during save, transform to grd file
			tempXYZOut = open(self.BUFDIR+self.manifestHead+"_HVsDropGrd.txt", 'w')
			saveDH = False
		
		
		hLen = len(hList)
		xLen = len(xList)
		
		hGrid = np.zeros(shape=(hLen, xLen))
		xGrid = np.zeros(shape=(hLen, xLen))
		HXGrid = np.zeros(shape=(hLen, xLen))
		
		## Load Buffer Once with LABModel
		## Also since synthetic seismograms have changed, I need to redefine output RF data using new fileName
		self.prepBuffer(copyManifest=True, useLAB=True)
		print "runing grid search for LAB depth ..."
		
		
		iCnt = 0
		TotCnt = hLen*xLen
		
		if True:
			for iDep in range(hLen):
				for iX in range(xLen):
					LABdepth  = hList[iDep]
					hGrid[iDep, iX] = LABdepth
					
					if useDhGrad > 0:
						dhGrad = xList[iX]
						xGrid[iDep, iX] = dhGrad
					if useDhGrad == 0:
						vsDrop = xList[iX]
						xGrid[iDep, iX] = vsDrop * 100
						
					
					self.pertModel(depthLAB=LABdepth, vsDrop=vsDrop, dhGrad=dhGrad)
					self.prepBuffer(copyManifest=False, useLAB=True)
					
					self.buildSynthSeismograms()
					self.buildRecFuncResponse(isMig=False, updateData=False, updateSynth=True, useLAB=True)
					
					# Plot Data ...
					RFdata  = self.DATARESULTDIR + self.manifestHead
					RFsynth = self.SYNTHRESULTDIR + self.manifestHead
					RFsynthLAB = self.SYNTHRESULTDIR + self.manifestHead + "_LAB"
					# LoaD 3 receiver function arrays - 1. Data 2. Synthetic-Crust 3.Synthetic-CRUSTLAB
					self.loadRFArrays(RFdata, RFsynth, RFsynthLAB)
					
					misfitVal = self.calcMisfit(minTime=lTime, maxTime=rTime)
					HXGrid[iDep, iX] = misfitVal
					
					# Write grd search results to ascii file
					tempXYZOut.write(str(hGrid[iDep,iX])+ "\t"+ str(xGrid[iDep,iX]) + "\t"+ str(HXGrid[iDep,iX])+ "\n")
					
					print "************************************* Grid node computation done: ", iCnt, "of", TotCnt
					print LABdepth, vsDrop, misfitVal
					iCnt = iCnt + 1
		
		tempXYZOut.close()
		self.plotNsave(hGrid,xGrid,HXGrid, flagXDh=saveDH)
		self.flushBuffer()
		
		return hGrid, xGrid, HXGrid
	
	def loadManifestFile(self, manifestFile, manifestHead):
		self.fDataManifest = manifestFile
		self.manifestHead = manifestHead

		print "Manifest saved:", self.fDataManifest
		print "Manifest Head: ", self.manifestHead
	
		
	def buildModel(self):
		pass
	
	def setModel(self, modelDir, modelFile, depthLAB = 0, vsDrop = 0, dhGrad=0):
		self.fvelModel = modelDir + modelFile
		print "Model file updated ", self.fvelModel
		
		## Translate model into internal array and plot locally .
		# Remove  plot routine into visualize routine ... CONTINUE HERE
		# jptab converts velocity model to format:
		# z(i),vp(i),vp2(i),vp4(i),vs(i),vs2(i),rho(i)
		# z - depth (km)
		# vp - p wave velocity
		# vp2 -  pk-to-pk cos(2th) relative P pert
		# vp4 -  pk-to-pk cos(4th) relative P pert
		# vs - s wave velocity
		# vs2 -  pk-to-pk cos(2th) relative s pert
		## rho - density (kg/m3)
	
		# Copy velModel and JP2TAB into BUFFER...
		p = subprocess.Popen(['cp', self.fvelModel, self.BUFDIR], stdout=subprocess.PIPE)
		output, err = p.communicate()
		
		p = subprocess.Popen(['cp', self.JP2TAB, self.BUFDIR], stdout=subprocess.PIPE)
		output, err = p.communicate()
		
		os.chdir(self.BUFDIR)
		p = subprocess.Popen(['jp2tab', modelFile], stdout=subprocess.PIPE)
		output, err = p.communicate()
		
		self.velModelArray =  np.loadtxt(modelFile+'.tab')
		lastIndex = len(self.velModelArray[:,0]) - 1
		
		# Sharp velocity drop.
		if depthLAB > 0 and vsDrop > 0 and dhGrad == 0:
			# copy velocity arrray and scale depth and velocity drop ...
			self.velModelArrayLAB = np.copy(self.velModelArray)
			
			# First column, last layer, (2 last rows)
			self.velModelArrayLAB[lastIndex -1, 0] = depthLAB
			self.velModelArrayLAB[lastIndex -2, 0] = depthLAB
			
			# Vp last indexs
			Vp = self.velModelArrayLAB[lastIndex -2,1]
			Vs = self.velModelArrayLAB[lastIndex -2, 4]
			
			# Scale halfspace below
			self.velModelArrayLAB[lastIndex, 1] = Vp * (1)
			self.velModelArrayLAB[lastIndex, 4] = Vs * (1-vsDrop)
		
			self.velModelArrayLAB[lastIndex-1, 1] = Vp * (1)
			self.velModelArrayLAB[lastIndex-1, 4] = Vs * (1-vsDrop)
		
		# Define velocity gradient here. simple code enhancement
		if (depthLAB > 0 and vsDrop > 0 and dhGrad > 0):
			self.SYNTHIMAGEDIR = self.SYNTHIMAGEDIR + "_grad_"
			nPts = 20
			
			# copy velocity arrray and scale depth and velocity drop ...
			self.velModelArrayLAB = np.copy(self.velModelArray)
			
			# h, Vp, and Vs last indexs
			#hAbove = self.velModelArrayLAB[lastIndex -2,0]
			# The dv is centered on the depthLAB
			hAbove = depthLAB - (0.5)*dhGrad
			VpAbove = self.velModelArrayLAB[lastIndex -2,1]
			VsAbove = self.velModelArrayLAB[lastIndex -2, 4]
			
			dh = hAbove + np.linspace(0,dhGrad,nPts)
			Vpi = VpAbove * np.ones(nPts)
			Vsi = VsAbove - ((VsAbove*vsDrop)/dhGrad) * (dh - hAbove)
			
			# First column, last layer, (2 last rows)
			#self.velModelArrayLAB[lastIndex -1, 0] = depthLAB
			#self.velModelArrayLAB[lastIndex -2, 0] = depthLAB
			
			newCol = self.velModelArrayLAB[lastIndex -2, :]
			halfspace = self.velModelArrayLAB[lastIndex, :]
			
			#print "model b4 delete"
			#print self.velModelArrayLAB
			
			## Delete last  last three rows - layer above and halfspace below
			self.velModelArrayLAB =  np.delete(self.velModelArrayLAB, lastIndex, 0)
			self.velModelArrayLAB =  np.delete(self.velModelArrayLAB, lastIndex-1, 0)
			self.velModelArrayLAB = np.delete(self.velModelArrayLAB, lastIndex-2, 0)
			
			#print "model b4 update"
			#print self.velModelArrayLAB
			## Then update with a gradational model.
			for iZ in range(nPts):
				#update depth, vp, and vs
				#print dh, Vpi, Vsi
				newCol[0] = dh[iZ]
				newCol[1] = Vpi[iZ]
				newCol[4] = Vsi[iZ]
				
				#print "Row", iZ
				#print "newCol", newCol
				
				#self.velModelArrayLAB.append(newCol)
				#self.velModelArrayLAB.append(newCol)
				self.velModelArrayLAB =  np.vstack((self.velModelArrayLAB, newCol))
				self.velModelArrayLAB = np.vstack((self.velModelArrayLAB, newCol))
			
			
			print "Gradational updated"
			# Scale halfspace below
			halfspace[1] = Vpi[nPts -1]
			halfspace[4] = Vsi[nPts -1]
			#self.velModelArrayLAB.append(halfspace)
			self.velModelArrayLAB = np.vstack((self.velModelArrayLAB, halfspace))
			
			#print "model after update"
			#print self.velModelArrayLAB

			#self.velModelArrayLAB[lastIndex, 1] = Vp * (1)
			#self.velModelArrayLAB[lastIndex, 4] = Vs * (1-vsDrop)
			
			#self.velModelArrayLAB[lastIndex-1, 1] = Vp * (1)
			#self.velModelArrayLAB[lastIndex-1, 4] = Vs * (1-vsDrop)
		# Leave model unchanged ..
		if depthLAB == 0 and vsDrop == 0 and dhGrad == 0:
			print "Only crustal model provided. No LAB depth or velocity drop specified"
		
		#write tabular data to JP file format
		lastIndex = len(self.velModelArrayLAB[:,0]) - 1
		print "New Layer Count:", lastIndex
		self.writeTAB2JP(self.velModelArrayLAB, lastIndex)
	
	# Single model perturbation used for gridStackingRoutine
	def pertModel(self, depthLAB=0, vsDrop=0, dhGrad=0):
		lastIndex = len(self.velModelArray[:,0]) - 1
		if depthLAB > 0 and vsDrop > 0 and dhGrad == 0:
			# copy velocity arrray and scale depth and velocity drop ...
			self.velModelArrayLAB = np.copy(self.velModelArray)
			
			# First column, last layer, (2 last rows)
			self.velModelArrayLAB[lastIndex -1, 0] = depthLAB
			self.velModelArrayLAB[lastIndex -2, 0] = depthLAB
			
			# Vp last indexs
			Vp = self.velModelArrayLAB[lastIndex -2,1]
			Vs = self.velModelArrayLAB[lastIndex -2, 4]
			
			# Scale halfspace below
			#self.velModelArrayLAB[lastIndex, 1] = Vp * (1-vsDrop)
			self.velModelArrayLAB[lastIndex, 1] = Vp * (1)
			self.velModelArrayLAB[lastIndex, 4] = Vs * (1-vsDrop)
			
			#self.velModelArrayLAB[lastIndex-1, 1] = Vp * (1-vsDrop)
			self.velModelArrayLAB[lastIndex-1, 1] = Vp * (1)
			self.velModelArrayLAB[lastIndex-1, 4] = Vs * (1-vsDrop)
		# Define velocity gradient here. simple code enhancement
		if (depthLAB > 0 and vsDrop > 0 and dhGrad > 0):
			#self.SYNTHIMAGEDIR = self.SYNTHIMAGEDIR + "_grad_"
			nPts = 15
			
			# copy velocity arrray and scale depth and velocity drop ...
			self.velModelArrayLAB = np.copy(self.velModelArray)
			
			# h, Vp, and Vs last indexs
			#hAbove = self.velModelArrayLAB[lastIndex -2,0]
			# The dv is centered on the depthLAB
			hAbove = depthLAB - (0.5)*dhGrad
			VpAbove = self.velModelArrayLAB[lastIndex -2,1]
			VsAbove = self.velModelArrayLAB[lastIndex -2, 4]
			
			dh = hAbove + np.linspace(0,dhGrad,nPts)
			Vpi = VpAbove * np.ones(nPts)
			Vsi = VsAbove - ((VsAbove*vsDrop)/dhGrad) * (dh - hAbove)
			
			# First column, last layer, (2 last rows)
			#self.velModelArrayLAB[lastIndex -1, 0] = depthLAB
			#self.velModelArrayLAB[lastIndex -2, 0] = depthLAB
			
			newCol = self.velModelArrayLAB[lastIndex -2, :]
			halfspace = self.velModelArrayLAB[lastIndex, :]
			
			#print "model b4 delete"
			#print self.velModelArrayLAB
			
			## Delete last  last three rows - layer above and halfspace below
			self.velModelArrayLAB =  np.delete(self.velModelArrayLAB, lastIndex, 0)
			self.velModelArrayLAB =  np.delete(self.velModelArrayLAB, lastIndex-1, 0)
			self.velModelArrayLAB = np.delete(self.velModelArrayLAB, lastIndex-2, 0)
			
			#print "model b4 update"
			#print self.velModelArrayLAB
			## Then update with a gradational model.
			for iZ in range(nPts):
				#update depth, vp, and vs
				#print dh, Vpi, Vsi
				newCol[0] = dh[iZ]
				newCol[1] = Vpi[iZ]
				newCol[4] = Vsi[iZ]
				
				#print "Row", iZ
				#print "newCol", newCol
				
				#self.velModelArrayLAB.append(newCol)
				#self.velModelArrayLAB.append(newCol)
				self.velModelArrayLAB =  np.vstack((self.velModelArrayLAB, newCol))
				self.velModelArrayLAB = np.vstack((self.velModelArrayLAB, newCol))
			
			
			#print "Gradational updated"
			# Scale halfspace below
			halfspace[1] = Vpi[nPts -1]
			halfspace[4] = Vsi[nPts -1]
			#self.velModelArrayLAB.append(halfspace)
			self.velModelArrayLAB = np.vstack((self.velModelArrayLAB, halfspace))

		if depthLAB == 0 and vsDrop == 0 and dhGrad == 0:
			#print "Only crustal model provided. No LAB depth or velocity drop specified"
			pass

		lastIndex = len(self.velModelArrayLAB[:,0]) - 1
		self.writeTAB2JP(self.velModelArrayLAB, lastIndex)
	
	
	# Show model parameters
	def getModel(self):
		print "Original Crust + Mantle Model >>>>>>>>>>>>>>>>>>>>>>"
		print "Header: z(i),vp(i),vp2(i),vp4(i),vs(i),vs2(i),rho(i) "
		print self.velModelArray
		print "scaled Model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
		print self.velModelArrayLAB
	
	
	def loadRFArrays(self, fileData, fileSynth, fileSynthLAB=None):
		for iNme in range(3):
			if (iNme == 0):
				fileName = fileSynth + self.M_MEAN_N
				#print fileName
				#print "loading synthetic arrays internally"
				synthHarmonicRFArray = np.loadtxt(fileName, comments='>')
				
				const = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 40]
				cosT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 30]
				sinT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 20]
				cos2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 10]
				sin2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 0]
				# Time: column1, const: column2, cosT, sinT ...
				self.synthHarmonicRFArray =  np.column_stack(( const[:,0], const[:,2], cosT[:,2] , sinT[:,2], cos2T[:,2], sin2T[:,2] ))
				
				# Load LAB Synthetic Array... Compare with base crust mantle model above.
				if fileSynthLAB is not None:
					fileName = fileSynthLAB + self.M_MEAN_N
					synthHarmonicRFArray = np.loadtxt(fileName, comments='>')
					const = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 40]
					cosT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 30]
					sinT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 20]
					cos2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 10]
					sin2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 0]
					# Time: column1, const: column2, cosT, sinT ...
					self.synthHarmonicRFArrayLAB =  np.column_stack(( const[:,0], const[:,2], cosT[:,2] , sinT[:,2], cos2T[:,2], sin2T[:,2] ))
				
				# Load Synthetic if second fileName is note None
				print "should load synthetic array. just mean. since no errors in synthetics."
				fileName = fileData + self.M_MEAN_N
				#print fileName
				#print "loading synthetic arrays internally"
				dataHarmonicRFArray = np.loadtxt(fileName, comments='>')
				
				const = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 40]
				cosT = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 30]
				sinT = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 20]
				cos2T = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 10]
				sin2T = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 0]
				
				# Time: column1, const: column2, cosT, sinT ...
				self.dataHarmonicRFArrayM =  np.column_stack(( const[:,0], const[:,2], cosT[:,2] , sinT[:,2], cos2T[:,2], sin2T[:,2] ))
			
			if (iNme == 1):
				fileName = fileData + self.M_PDTHETA_N
				#print fileName
				#print "loading synthetic arrays internally"
				dataHarmonicRFArray = np.loadtxt(fileName, comments='>')
				
				const = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 40]
				cosT = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 30]
				sinT = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 20]
				cos2T = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 10]
				sin2T = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 0]
				
				# Time: column1, const: column2, cosT, sinT ...
				self.dataHarmonicRFArrayP =  np.column_stack(( const[:,0], const[:,2], cosT[:,2] , sinT[:,2], cos2T[:,2], sin2T[:,2] ))
			if (iNme == 2):
				fileName = fileData + self.M_NDTHETA_N
				#print fileName
				#print "loading synthetic arrays internally"
				dataHarmonicRFArray = np.loadtxt(fileName, comments='>')
				
				const = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 40]
				cosT = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 30]
				sinT = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 20]
				cos2T = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 10]
				sin2T = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 0]
				
				# Time: column1, const: column2, cosT, sinT ...
				self.dataHarmonicRFArrayN =  np.column_stack(( const[:,0], const[:,2], cosT[:,2] , sinT[:,2], cos2T[:,2], sin2T[:,2] ))

	
	# Internal routine for displaying RF...
	def displayRFOld(self, PassAx, fileHead, fileHead2 = None, scaleVal = None):
		for iNme in range(3):
			if (iNme == 0):
				fileName = fileHead + self.M_MEAN_N
				#print fileName
				#print "loading synthetic arrays internally"
				synthHarmonicRFArray = np.loadtxt(fileName, comments='>')
				
				const = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 40]
				cosT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 30]
				sinT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 20]
				cos2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 10]
				sin2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 0]
				# Time: column1, const: column2, cosT, sinT ...
				self.synthHarmonicRFArrayM =  np.column_stack(( const[:,0], const[:,2], cosT[:,2] , sinT[:,2], cos2T[:,2], sin2T[:,2] ))
			if (iNme == 1):
				fileName = fileHead + self.M_PDTHETA_N
				#print fileName
				#print "loading synthetic arrays internally"
				synthHarmonicRFArray = np.loadtxt(fileName, comments='>')
				
				const = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 40]
				cosT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 30]
				sinT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 20]
				cos2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 10]
				sin2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 0]
				
				# Time: column1, const: column2, cosT, sinT ...
				self.synthHarmonicRFArrayP =  np.column_stack(( const[:,0], const[:,2], cosT[:,2] , sinT[:,2], cos2T[:,2], sin2T[:,2] ))
			if (iNme == 2):
				fileName = fileHead + self.M_NDTHETA_N
				#print fileName
				#print "loading synthetic arrays internally"
				synthHarmonicRFArray = np.loadtxt(fileName, comments='>')
				
				const = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 40]
				cosT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 30]
				sinT = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 20]
				cos2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 10]
				sin2T = synthHarmonicRFArray[synthHarmonicRFArray[:,1] == 0]
				
				# Time: column1, const: column2, cosT, sinT ...
				self.synthHarmonicRFArrayN =  np.column_stack(( const[:,0], const[:,2], cosT[:,2] , sinT[:,2], cos2T[:,2], sin2T[:,2] ))
		# Load Synthetic if second fileName is note None
		if fileHead2 is not None:
			print "should load synthetic array. just mean. since no errors in synthetics."
			fileName = fileHead2 + self.M_MEAN_N
			#print fileName
			#print "loading synthetic arrays internally"
			dataHarmonicRFArray = np.loadtxt(fileName, comments='>')
			
			const = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 40]
			cosT = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 30]
			sinT = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 20]
			cos2T = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 10]
			sin2T = dataHarmonicRFArray[dataHarmonicRFArray[:,1] == 0]
			
			# Time: column1, const: column2, cosT, sinT ...
			self.dataHarmonicRFArray =  np.column_stack(( const[:,0], const[:,2], cosT[:,2] , sinT[:,2], cos2T[:,2], sin2T[:,2] ))
		#print "Dimension of synthetic array: ", self.synthHarmonicRFArrayM.shape
		# Trying to mimick gmt plot display ...
		#fig = plt.figure(figsize = (6,6))
		t = self.synthHarmonicRFArrayM
		t1 = self.synthHarmonicRFArrayN
		t2 = self.synthHarmonicRFArrayP
		
		if scaleVal is None:
			maxConst = np.max(t[:,1])
		else:
			maxConst = scaleVal
					
		offset = 0
		
		for i in [5,4,3,2, 1]:
			
			if (i == 5):
				x = t[:,0]
				y = t[:,5]/maxConst
				y1 = t1[:,5]/maxConst
				y2 = t2[:,5]/maxConst

				PassAx.fill_between(x, 0, y2, y2 >  0, color='#99FF99')
				PassAx.fill_between(x, 0, y1, y1 <  0, color='#99FF99')

				PassAx.fill_between(x, 0, y1, y1 >  0, color='blue')
				PassAx.fill_between(x, 0, y2, y2 <  0, color='red')

				PassAx.plot(x, y, 'k', linewidth =2)
				
				if fileHead2 is not None:
					ts = self.dataHarmonicRFArray
					x = ts[:,0]
					y = ts[:,5]/maxConst
					PassAx.plot(x, y, 'k--', linewidth=2)
					PassAx.fill_between(x, offset, y, y >  offset, color='white', alpha = 0.4)
					PassAx.fill_between(x, offset, y, y <  offset, color='white', alpha = 0.4)
			else:
				offset += 1.5 * maxConst
				
				x = t[:,0]
				y = t[:,i]/maxConst + offset

				y1 = t1[:,i]/maxConst + offset
				y2 = t2[:,i]/maxConst + offset

				PassAx.fill_between(x, offset, y2, y2 >  offset, color='#99FF99')
				PassAx.fill_between(x, offset, y1, y1 <  offset, color='#99FF99')

				PassAx.fill_between(x, offset, y1, y1 >  offset, color='blue')
				PassAx.fill_between(x, offset, y2, y2 <  offset, color='red')

				PassAx.plot(x, y, 'k', linewidth =2)
	
				if fileHead2 is not None:
					ts = self.dataHarmonicRFArray
					x = ts[:,0]
					y = ts[:,i]/maxConst + offset
					PassAx.plot(x, y, 'k--', linewidth=2)
					PassAx.fill_between(x, offset, y, y >  offset, color='white', alpha = 0.4)
					PassAx.fill_between(x, offset, y, y <  offset, color='white', alpha = 0.4)
				
		#plt.plot(t[:,0], t[:,2])
		offset += 1.5
		PassAx.set_xlim([-2,10])
		PassAx.set_ylim([-1,offset])
		PassAx.set_yticks([])
		return maxConst
		#plt.show()
	
	
	# Internal routine for displaying RF...
	def displayRF(self, PassAx, showData = True, showSynth = False, scaleVal = None, displayConstOnly = False):
		if showData:
			t = self.dataHarmonicRFArrayM
			t1 = self.dataHarmonicRFArrayN
			t2 = self.dataHarmonicRFArrayP
		if showSynth:
			ts = self.synthHarmonicRFArrayLAB
			ts1 = self.synthHarmonicRFArrayLAB
			ts2 = self.synthHarmonicRFArrayLAB
		if showData and showSynth:
			ts = self.synthHarmonicRFArray
			ts1 = self.synthHarmonicRFArray
			ts2 = self.synthHarmonicRFArray
			tsLAB = self.synthHarmonicRFArrayLAB
		
		if scaleVal is None:
			maxConst = np.max(t[:,1])
		else:
			maxConst = scaleVal
		


		if displayConstOnly == True:
			maxI = 1
			offset = 0
			offsetTick = []
		else:
			maxI = 5
			offset = 0
			offsetTick = [0]

		for i in np.arange(maxI, 0, -1):
			# Plot last harmonic stack here.
			if (i == 5):
				if showData:
					x = t[:,0]
					y = t[:,5]/maxConst
					y1 = t1[:,5]/maxConst
					y2 = t2[:,5]/maxConst
				if showSynth:
					x = ts[:,0]
					y = ts[:,5]/maxConst
					y1 = ts1[:,5]/maxConst
					y2 = ts2[:,5]/maxConst
				if showData and showSynth:
					x = t[:,0]
					y = t[:,5]/maxConst
					y1 = t1[:,5]/maxConst
					y2 = t2[:,5]/maxConst
				
				PassAx.fill_between(x, 0, y2, y2 >  0, color='#99FF99')
				PassAx.fill_between(x, 0, y1, y1 <  0, color='#99FF99')
				
				PassAx.fill_between(x, 0, y1, y1 >  0, color='blue')
				PassAx.fill_between(x, 0, y2, y2 <  0, color='red')
				
				PassAx.plot(x, y, 'k', linewidth =2)
				
				if showData and showSynth:
					x = ts[:,0]
					y = ts[:,5]/maxConst
					PassAx.plot(x, y, 'k--', linewidth=2)
					PassAx.fill_between(x, offset, y, y >  offset, color='white', alpha = 0.4)
					PassAx.fill_between(x, offset, y, y <  offset, color='white', alpha = 0.4)
			else:
				if displayConstOnly is True:
					offset = 0
				else:
					offset += 1.5 * maxConst

				offsetTick.append(offset)
				
				if showData:
					x = t[:,0]
					y = t[:,i]/maxConst + offset
					y1 = t1[:,i]/maxConst + offset
					y2 = t2[:,i]/maxConst + offset
				if showSynth:
					x = ts[:,0]
					y = ts[:,i]/maxConst + offset
					y1 = ts1[:,i]/maxConst + offset
					y2 = ts2[:,i]/maxConst + offset
				if showData and showSynth:
					x = t[:,0]
					y = t[:,i]/maxConst + offset
					y1 = t1[:,i]/maxConst + offset
					y2 = t2[:,i]/maxConst + offset
				
				# Plot Data (Or Synthetic) If Both. THen Data First
				PassAx.fill_between(x, offset, y2, y2 >  offset, color='#99FF99')
				PassAx.fill_between(x, offset, y1, y1 <  offset, color='#99FF99')
				
				PassAx.fill_between(x, offset, y1, y1 >  offset, color='blue')
				PassAx.fill_between(x, offset, y2, y2 <  offset, color='red')
				PassAx.plot(x, y, 'k', linewidth =2)

				# Then Synthetic Next, If Both Plots are requested, then redo ..
				if showData and showSynth:
					# Plot crust-mantle synthetics
					x = ts[:,0]
					y = ts[:,i]/maxConst + offset
					PassAx.plot(x, y, 'k-.', linewidth=2)
					PassAx.fill_between(x, offset, y, y >  offset, color='white', alpha = 0.4)
					PassAx.fill_between(x, offset, y, y <  offset, color='white', alpha = 0.4)
						
					# Replot Data Again ....
					offset += 1.5 * maxConst
					offsetTick.append(offset)
					
					x = t[:,0]
					y = t[:,i]/maxConst + offset
					y1 = t1[:,i]/maxConst + offset
					y2 = t2[:,i]/maxConst + offset
					
					PassAx.fill_between(x, offset, y2, y2 >  offset, color='#99FF99')
					PassAx.fill_between(x, offset, y1, y1 <  offset, color='#99FF99')

					PassAx.fill_between(x, offset, y1, y1 >  offset, color='blue')
					PassAx.fill_between(x, offset, y2, y2 <  offset, color='red')
					PassAx.plot(x, y, 'k', linewidth =2)

					
					# Plot LAB synthetics
					x = tsLAB[:,0]
					y = tsLAB[:,i]/maxConst + offset
					PassAx.plot(x, y, 'k--', linewidth=2)
					PassAx.fill_between(x, offset, y, y >  offset, color='white', alpha = 0.4)
					PassAx.fill_between(x, offset, y, y <  offset, color='white', alpha = 0.4)
					
					minorLocator   = AutoMinorLocator()
					PassAx.xaxis.set_minor_locator(minorLocator)
		
		#plt.plot(t[:,0], t[:,2])
		offset += 1.5
		PassAx.set_xlim([-2,10])
		PassAx.set_ylim([-1,offset])
		PassAx.set_yticks(offsetTick)

		if displayConstOnly is True:
			PassAx.set_yticklabels([])
			PassAx.text(6, offsetTick[0]-0.5, 'No LAB', fontsize = 15)
			PassAx.text(6, offsetTick[1]-0.5, 'LAB', fontsize = 15)
		else:
			PassAx.set_yticklabels(['$\sin2\pi$', '$\cos2\pi$', '$\sin\pi$', '$\cos\pi$', '$con$'], rotation=40,fontsize=14)

		return maxConst
	#plt.show()
	
	# Show model, dataRF and synthRF
	def showAll(self, format  = 'png'):
		fig = plt.figure(figsize = (11,8))
		
		
		
		axRatio  = plt.subplot2grid((2,3), (0,0), rowspan = 2)
		dataSynthAx = plt.subplot2grid((2,3), (0,1), colspan = 2)
		dataAx = plt.subplot2grid((2,3), (1,1) )
		synthAx = plt.subplot2grid((2,3), (1,2) )
		
		# Plot Velocity ratio
		print "vp/vs ratio", self.velModelArray[:,1] / self.velModelArray[:,4]
		axVel = axRatio.plot(self.velModelArray[:,1] / self.velModelArray[:,4], self.velModelArray[:,0], 'grey', linewidth=1)
		axRatio.set_xlim([1.5, 2.0])
		axRatio.set_xlabel("vp/vs")
		
		velAx = axRatio.twiny()
		
		## Plot VELOCITY MODEL CRUST ...>>>=
		velAx.plot(self.velModelArray[:,1], self.velModelArray[:,0], 'k-.', linewidth=2)
		velAx.plot(self.velModelArray[:,4], self.velModelArray[:,0], 'k-.', linewidth=2 )
		
		## Plot VELOCITY MODEL CRUST + LAB ...>>>=
		velAx.plot(self.velModelArrayLAB[:,1], self.velModelArrayLAB[:,0], 'k--', linewidth=2)
		velAx.plot(self.velModelArrayLAB[:,4], self.velModelArrayLAB[:,0], 'k--', linewidth=2 )
		
		velAx.set_xlim([2.0, 9.0])
		velAx.set_xlabel("vs and vp")
		velAx.set_ylabel("Depth (km)")
		
		

		
		velAx.invert_yaxis()
		
		# Plot Data ...
		RFdata  = self.DATARESULTDIR + self.manifestHead
		RFsynth = self.SYNTHRESULTDIR + self.manifestHead
		RFsynthLAB = self.SYNTHRESULTDIR + self.manifestHead + "_LAB"
		
		# LoaD 3 receiver function arrays - 1. Data 2. Synthetic-Crust 3.Synthetic-CRUSTLAB
		self.loadRFArrays(RFdata, RFsynth, RFsynthLAB)
		
		
		#scaleDat = self.displayRFOld(dataAx, RFdata)
		scaleDat = self.displayRF(dataAx, showData=True, showSynth=False, scaleVal=None, displayConstOnly= False)
		dataAx.set_title("Observed")
	
		
		
		# Plot Synthetics
		#scaleSyn = self.displayRFOld(synthAx, RFsynth, None, scaleDat)
		scaleSyn = self.displayRF(synthAx, showData=False, showSynth=True, scaleVal=scaleDat, displayConstOnly= False)
		synthAx.set_title("Predicted")
		
		# Plot Both Data and Synthetic
		#self.displayRFOld(dataSynthAx, RFdata, RFsynth)
		scaleSyn = self.displayRF(dataSynthAx, showData=True, showSynth=True, scaleVal=None, displayConstOnly = True)
		dataSynthAx.set_title(self.manifestHead + ": Observed + Predicted")
		
		plt.tick_params(which='both', width=2)
		plt.tick_params(which='major', length=7)
		plt.tick_params(which='minor', length=4, color='r')
		
		plt.savefig(self.SYNTHIMAGEDIR+'All_'+self.manifestHead+'.'+format, format=format)
		plt.show()
		
	
	def showModel(self):
		print "quick plot of velocity model, with other diagnostics: "
		fig = plt.figure(figsize = (6,6))
		
		plt.subplot(1,2,1)
		# Plot P velocity with depth
		axVel = plt.plot(self.velModelArray[:,1], self.velModelArray[:,0], self.velModelArray[:,4], self.velModelArray[:,0] )
		plt.gca().set_xlim([3, 8])
		plt.gca().set_xlabel("vs and vp")
		plt.gca().set_ylabel("Depth (km)")
		axRatio = plt.gca().twiny()
		
		print "vp/vs ratio moho:", self.velModelArray[:,1] / self.velModelArray[:,4]
		axRatio.plot(self.velModelArray[:,1] / self.velModelArray[:,4], self.velModelArray[:,0], 'r')
		axRatio.set_xlim([1, 8])
		axRatio.set_xlabel("moho vp/vs")

		plt.gca().invert_yaxis()
		
		plt.subplot(1,2,2)
		# Plot P velocity with depth
		axVel = plt.plot(self.velModelArrayLAB[:,1], self.velModelArrayLAB[:,0], self.velModelArrayLAB[:,4], self.velModelArrayLAB[:,0] )
		plt.gca().set_xlim([3, 8])
		plt.gca().set_xlabel("vs and vp")
		#plt.gca().set_ylabel("Depth (km)")
		axRatio = plt.gca().twiny()
		
		print "vp/vs ratio mohoLAB", self.velModelArrayLAB[:,1] / self.velModelArrayLAB[:,4]
		axRatio.plot(self.velModelArrayLAB[:,1] / self.velModelArrayLAB[:,4], self.velModelArrayLAB[:,0], 'r')
		axRatio.set_xlim([1, 8]), 	axRatio.set_yticklabels([])
		axRatio.set_xlabel("mohoLAB vp/vs")
		plt.gca().invert_yaxis()


		plt.show()
		
	def writeTAB2JP(self, modelArray2Save, lastIndex):
		fileSve = open(self.BUFDIR+'LABmodel.txt', 'w')
		#print "Generated in synthethicRFHarmonic.py moho model with LAB"
		fileSve.write("Generated in synthethicRFHarmonic.py moho model with LAB \n")
		listRows = np.arange(1, lastIndex+1, 2)
		noLayers =  len(listRows) - 1
		#print noLayers
		fileSve.write(str(noLayers) + "\n")
		for rowWrite in listRows:
			#print "0  0"
			fileSve.write("0  0\n")
			#print modelArray2Save[rowWrite,:] * 1000.0
			rowPrint = modelArray2Save[rowWrite,:] * 1000.0
			fileSve.write(' '.join(map(str, rowPrint)) + "\n")
		
		fileSve.close()
		pass
	
	def flushBuffer(self):
		#p = subprocess.Popen(['rm', self.BUFDIR+"*.*"], stdout=subprocess.PIPE)
		#output, err = p.communicate()
		try:
			shutil.rmtree(self.BUFDIR)
		except OSError, e:
			#print  output
			print("Error: %s - %s." % (e.filename, e.strerror))
		os.mkdir(self.BUFDIR)
	
	
	def prepBuffer(self, copyManifest = True, useLAB=False):
		if self.fvelModel == "" or self.fDataManifest == "" :
			print "No model file or data file set, can't run anirecsynth. Need Model & Data File"
		else:
			#print "running ", self.fvelModel, " on ", self.fDataManifest
			paramTable = []
			
			paramTable.append('1')  # Pulse Type - 1 sided pulse
			paramTable.append( '2.0')  # Incident Wave Period
			
			# Use velocity model provided or internally generated LAB model ..
			if useLAB == True:
				paramTable.append('LABmodel.txt')
			else:
				splitName = self.fvelModel.split("/")
				paramTable.append(splitName[len(splitName) -1])  # <<<< MODEL FILE
			
			
			paramTable.append('SAClist.txt')  # List of SAC Files
			#save parameter file ...
			np.savetxt(self.SYNTHPARAMFILE, paramTable, fmt='%s')
			
			fileNames = np.loadtxt(self.fDataManifest, dtype = 'string')
			self.nEvents = len(fileNames)
			
			print "Total Events", self.nEvents
			
			# Copy files into buffer: 3 component data seismograms, anirec code, and fortranCompatibleFileNames
			if copyManifest == True:
				forCompNames = []

				for fileName in fileNames:
					splitName = fileName.split("/")
					#print "file Directory array: ", splitName, len(splitName), splitName[len(splitName) -1]
					
					forCompNames.append(splitName[len(splitName) -1]+'Z')
					
					p = subprocess.Popen(['cp', fileName+'R', self.BUFDIR], stdout=subprocess.PIPE)
					output, err = p.communicate()
					
					p = subprocess.Popen(['cp', fileName+'T', self.BUFDIR], stdout=subprocess.PIPE)
					output, err = p.communicate()
					
					p = subprocess.Popen(['cp', fileName+'Z', self.BUFDIR], stdout=subprocess.PIPE)
					output, err = p.communicate()
					
					#print fileName, self.BUFDIR, output, err
					# Copy kernelSynth, and truncManifest (for fortranCompatibility) into BUFFER...
					
				p = subprocess.Popen(['cp', self.KERNELSYNTH, self.BUFDIR+'anirec'], stdout=subprocess.PIPE)
				output, err = p.communicate()
					
				# save fortran compatible manifest list in buffer directory
				forCompNames.append("stop")
				np.savetxt(self.BUFDIR+'SAClist.txt', forCompNames, fmt='%s')


	def buildSynthSeismograms(self):
		# run anirec_synth on velocity model for sac data
		
		os.chdir(self.BUFDIR)
		CommandArgumentList = ['gtime', 'anirec']

		f = open(self.SYNTHPARAMFILE)
		p = subprocess.Popen(CommandArgumentList, stdin=f, stdout=subprocess.PIPE)
		output, err = p.communicate()

		# append directory name to synthetic manifest list
		newManifestList = []
		fileNames = np.loadtxt(self.BUFDIR+'synthSAClist.txt', dtype = 'string')

		for fileName in fileNames:
				newManifestList.append(self.BUFDIR + fileName)
					
		np.savetxt(self.BUFDIR+'synthSAClist.txt', newManifestList, fmt='%s')
		#print output  # only output if you need to



	# Use runRFmodules to generate RF response, and then to print file
	def buildRecFuncResponse(self, isMig = True, updateData = True, updateSynth = True, useLAB=False):
		print "creating synthetic rf response using basic harmonic stacks, no migration ..."
		freqCut = 1
		migDepth = 3
		
		headCode = 1    # 1. for QC tagged files. 0. otherwise
		lqtRot = 0.5    # use to scale lqt rotation
			
		dataURL = self.DATARESULTDIR + self.manifestHead
		if useLAB == True:
			synthURL = self.SYNTHRESULTDIR + self.manifestHead + "_LAB"
		else:
			synthURL = self.SYNTHRESULTDIR + self.manifestHead
		
		
		# Generate both synthetics and RF data - runRFcodes Twice.
		for runtwice in range(2):
			if runtwice == 0 and updateData == True:
				print "Recomputing RF Response for Data"
				manifestURL = self.fDataManifest
				stackURL = dataURL
				migrateURL = self.fvelModel
				#1.0 Run and Visualize Single Harmonic Stack ******************#####################**********************
				if isMig:
					runRF.runHarmonicRFstackWithMigrate(manifestURL, stackURL, migrateURL, migDepth, freqCut, headCode,
														lqtRot);
				else:
					# 2.0 Run stack with no migration ******************#####################**********************
					runRF.runHarmonicRFstack(manifestURL, stackURL, freqCut, lqtRot)
			if runtwice == 1 and updateSynth == True:
				print "Recomputing RF Response for Synthetics"
				manifestURL = self.BUFDIR+'synthSAClist.txt'
				stackURL = synthURL
				print "file save", stackURL
				migrateURL = self.fvelModel
				#1.0 Run and Visualize Single Harmonic Stack ******************#####################**********************
				if isMig:
					runRF.runHarmonicRFstackWithMigrate(manifestURL, stackURL, migrateURL, migDepth, freqCut, headCode,
													lqtRot);
				else:
					# 2.0 Run stack with no migration ******************#####################**********************
					runRF.runHarmonicRFstack(manifestURL, stackURL, freqCut, lqtRot)
		manifestHead = self.manifestHead + '_synth'
		pltDir = self.PLTDRIVERS
		saveURL = self.SYNTHIMAGEDIR + self.manifestHead
		#runRF.plotHarmonicRFstackWithMigrate(manifestHead, pltDir, stackURL, saveURL, migDepth, freqCut);


	# Function that calculates misfit between synthetic and data (first version works with just the constant stacks)
	# maxTime is limit of data to fit.
	def calcMisfit(self, minTime, maxTime, showMisfit = False):
		# Plot Data ...
		#RFdata  = self.DATARESULTDIR + self.manifestHead
		#RFsynth = self.SYNTHRESULTDIR + self.manifestHead
		
		#self.loadRFArrays(RFdata, RFsynth)
		
		time = self.dataHarmonicRFArrayM[:,0]
		index = (time < maxTime) & (time > minTime)
		
		dataConst = self.dataHarmonicRFArrayM[:,1]
		synthConst = self.synthHarmonicRFArrayLAB[:,1]
		
		dataConst = dataConst[index]
		synthConst = synthConst[index]
		
		plusDat = self.dataHarmonicRFArrayP[:,1]
		negDat = self.dataHarmonicRFArrayN[:,1]
		
		plusDat = plusDat[index]
		negDat = negDat[index]
		
		if (showMisfit):
			plt.subplot(3,1,1)
			plt.plot(time[index], dataConst[index], time[index], synthConst[index])
			plt.plot(time[index], plusDat[index], 'r--')
			plt.plot(time[index], negDat[index], 'm--')
			
	
		# Plot misfit emulators - slice within time window ...
		
		error = (dataConst - synthConst) ** 2
		var = (plusDat - negDat) ** 2
		
		sd = np.sqrt(var)
		weight = ((dataConst - sd) ** 2/ var)
		weightError =  np.sqrt( error/var)
		
		
		#print "weightErorr", weightError
		#print "data", np.mean(dataConst)
		#print "synthetic", np.mean(synthConst)
		#print "misfit values [error], [variance]", np.mean(error), np.mean(var)
		
		if (showMisfit):
			plt.subplot(3,1,2)
			plt.plot(time[index], error[index], 'r--')
			plt.plot(time[index], var[index], 'b--')
			
			plt.subplot(3,1,3)
			plt.plot(time[index], weightError[index], 'k')
		
		rootMeanSquare = np.mean(weightError)
		print "Misfit value, root weighted square residual ", rootMeanSquare
		return rootMeanSquare

	def plotNsave(self, hGrid,xGrid,HXGrid, flagXDh = False, format='ps'):
		
		if True:
			if flagXDh == False:
				xyzFile = self.BUFDIR+self.manifestHead+"_HVsDropGrd.txt"
				#"RFDataStacks/gridStacks/HVsDrop/"
				grdFile = self.GRDRESULTDIR + "HVsDrop/" + self.manifestHead+'.nc'
				grdSave = self.INVIMAGEDIR +  "HVsDrop/grd_" + self.manifestHead
				pltYlabel = r"$\frac{\Delta Vs}{Vs}$ (%)"
			if flagXDh == True:
				xyzFile = self.BUFDIR+self.manifestHead+"_HDhGradGrd.txt"
				grdFile = self.GRDRESULTDIR + "HDhGrad/" + self.manifestHead+'.nc'
				grdSave = self.INVIMAGEDIR +  "HDhGrad/grd_" + self.manifestHead
				pltYlabel = "$\Delta$ H (km)"
				
			
			
			# x and y output of grid
			H_out = hGrid[:,1]
			X_out = xGrid[1,:]
			
			hstep = str(H_out[1] - H_out[0])
			xstep = str(X_out[1] - X_out[0])
			
			hMin = str(H_out[0])
			hMax = str(H_out[len(H_out)-1])
			
			xMin = str(X_out[0])
			xMax = str(X_out[len(X_out)-1])

			print hstep, xstep, hMin, hMax, xMin, xMax
				
			xyz2grdPrompt = ['xyz2grd', xyzFile, '-G'+grdFile, '-I'+hstep+'/'+xstep,'-R'+hMin+'/'+hMax+'/'+xMin+'/'+xMax, '-N0']
			p = subprocess.Popen(xyz2grdPrompt, stdout=subprocess.PIPE)
			output, err = p.communicate()
			print "*** Convert to grid ", "errors?: \n", output
					
			# the output array to write will be nH x nY
			#nH = len(hGrid[:,1]); nY = len(xGrid[1,:])

			# open a new netCDF file for writing.
			
			#ncfile = Dataset(self.GRDRESULTDIR+self.manifestHead+'.nc',mode='w',clobber=True)

			# output data.
			#grid_out = HXGrid.T # 2d array

			# create the X and Y dimensions.
			#ncfile.createDimension('LABDepth',nH)
			#ncfile.createDimension('VsDrop',nY)


			# Define the coordinate variables. They will hold the coordinate
			# information, that is, the LABDepth and VsDrop.
			#H = ncfile.createVariable('LABDepth',dtype('float32').char,('LABDepth',))
			#Y = ncfile.createVariable('VsDrop',dtype('float32').char,('VsDrop',))

			# Assign units attributes to coordinate var data. This attaches a
			# text attribute to each of the coordinate variables, containing the
			# units.
			#H.units = 'km'
			#Y.units = '%'
			# write data to coordinate vars.
			#H[:] = H_out
			#Y[:] = Y_out

			# create the Z variables (2D output grid)
			#RMSOut = ncfile.createVariable('RMS error Output',dtype('float32').char,('LABDepth','VsDrop'))

			# set the units attribute.
			#RMSOut.units =  'RFamp'

			# write data to variables.
			#RMSOut[:] = grid_out

			# close the file.
			#ncfile.close()
			print '*** SUCCESS writing' + self.GRDRESULTDIR+self.manifestHead	+'.nc'
		
		## Start of by writing file to netCdf format ..
		fig = plt.figure(figsize = (4,3))
		
		CS = plt.contourf(hGrid, xGrid, HXGrid,cmap=plt.cm.RdBu)
		CS2 = plt.contour(hGrid, xGrid, HXGrid, linewidth=0.5, colors='w')
		cbar = plt.colorbar(CS, format="%1.2f")
		
		plt.xlabel("H (km)", fontsize=16)
		plt.ylabel(pltYlabel, fontsize=16)
		cbar.ax.set_ylabel('misfit ', fontsize=16)
		
		plt.savefig(grdSave+'.'+format, format=format)
		plt.show()






