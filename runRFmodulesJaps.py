# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Author: Tolulope Olugboji
# Date:  January 30, 2013
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# Objective:   Debug, Test, and Run Multi-HK stacks in the frequency domain
#              Used in discussing OceanIsland Anistoropy, and modeling crustal structure..
#
#       
#
#

import subprocess
import time, sys, os
from scipy.io import netcdf 
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

import matplotlib as mpl
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['xtick.color'] = 'k'

def recompileCode(type):
    makeDir = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/Code/RFMigHarmonic'
    binDir = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/Code/bin'
    makeFile = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/Code/RFMigHarmonic/MakefileOnMAC'
    
    if (type == 'RFVelStck'):
        print "MultiHK grd search stack in progress.."
        makeArgument = ['make', 'RFVelStck', '-f' + makeFile]
        moveArgument = ['mv', './RFVelStck', binDir]
    else:
        print "RecFuncHarmonic Compilation and Link In Progress.."
        makeArgument = ['make', 'RecFunc', '-f' + makeFile]
        binDir = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/Code/bin/RecFuncHarmonic'
        moveArgument = ['mv', './RecFunc', binDir]
    
    # Recompile code ...
    os.chdir(makeDir)
    status = subprocess.call(makeArgument)
    
    # If successful move to bin Directory ...
    
    if ( status == 0):
        print "*** Succesful Make! ",  status, "No Errors"
    
        status = subprocess.call(moveArgument)
        
        if( status == 0): 
            print " Successful update of ", type, 'Code!', status
    else:
        print "Make Fail.. Revisit code and recompile !!!"
        
#0.0 Run Synthetics and Harmonics For A particular anisotropic model
def runSyntheticsHarmonic(synthURL, outURL, manifestHead, freqCut):
    
    
    #Flag used for debug purposes. close and it doesn't work
    if (1):
        print "Preparing input parameters ..."
        fileParams = open(synthURL+'runfile', 'wb')
        
        fileParams.write('1 \n')  # Pulse Type - 1 sided pulse
        fileParams.write(manifestHead + '_mod.txt \n')  # ANISOTROIC MODEL 
        fileParams.write(manifestHead + '.txt \n')  # List of SAC Files
        fileParams.write( '2.0 \n')  # Incident Wave Period
        
        manifestFile = manifestHead + '.txt'
        #wcOut = !wc -l  $manifestFile
    
        
        #for elem in wcOut:
        #   recNo = elem.split()
        
        #response  = int(recNo[0])
        
        p = subprocess.Popen(['wc', '-l', manifestFile], stdout=subprocess.PIPE)
        output, err = p.communicate()
        line =  output.split()
        response = int(line[0])
        
        for i in range(response):
            fileParams.write('1 \n')  # Pulse Type - Incident Wave
            fileParams.write('2 \n')  # Incident Wave Period
        
        fileParams.close()
        
    
    CommandArgumentList = ['anirec_synth_new']
    f = open("runfile")
    
    p = subprocess.Popen(CommandArgumentList, stdin=f, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    
    #output = !anirec_synth < runfile > /dev/null
    print "*** Running anirecsynth on ", manifestHead, "errors?: \n",  output

# 1.0 Epicentral Crude Stack
def runRFstackByEpic(manifestURL, outURL, azimWinMin, azimWinMax, azimWinStep, epicWinMin, epicWinMax, epicWinStep, freqCut):
    

    epicOption = str(azimWinMin) + '/' + str(azimWinMax) + '/' + str(epicWinMin) + '/' + str(epicWinMax) + '/' + str(epicWinStep)
    
    CommandArgumentList = ['RecFuncHarmonic', '-L' + manifestURL,
                           '-F' + str(freqCut), '-T60', '-O' + outURL,
                          '-B2', '-S0/50', '-E' + epicOption, '-H1' , '-V1', '-R0.7']
    
    #azimOption = str(azimWinMin) + '/' + str(azimWinMax) + '/'  + str(azimWinStep)
    
    #CommandArgumentList = ['RecFuncHarmonic', '-L' + manifestURL,
    #                       '-F' + str(freqCut), '-T60', '-O' + outURL,
    #                       '-B1', '-S0/50', '-A' + azimOption, '-H0', '-V1', '-R1.0']

    
    p = subprocess.Popen(CommandArgumentList, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running RecFuncHarmonic on ", outURL, "errors?: \n", err, output
    

# 2.0 Azimuth Crude Stack
def runRFstackByAzim(manifestURL, outURL, azimWinMin, azimWinMax, azimWinStep, freqCut, headCode):

    azimOption = str(azimWinMin) + '/' + str(azimWinMax) + '/'  + str(azimWinStep)
    
    CommandArgumentList = ['RecFuncHarmonic', '-L' + manifestURL,
                           '-F' + str(freqCut), '-T60', '-O' + outURL,
                           '-B1', '-S0/50', '-A' + azimOption, '-H'+str(headCode), '-V1', '-R0.55']

    
    p = subprocess.Popen(CommandArgumentList, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running RecFuncHarmonic on ", outURL, "errors?: \n", err, output

# 2.1 Azimuth Crude Stack - with horizontal rotation [Jap OBS Data]
def runRFstackByAzimWithRot(manifestURL, outURL, azimWinMin, azimWinMax, azimWinStep, freqCut, headCode,rotAngle):

    azimOption = str(azimWinMin) + '/' + str(azimWinMax) + '/'  + str(azimWinStep)
    
    CommandArgumentList = ['RecFuncHarmonic', '-L' + manifestURL,
                           '-F' + str(freqCut), '-T60', '-O' + outURL,
                           '-B1', '-S0/50', '-A' + azimOption, '-H'+str(headCode), '-V1', '-R0.55',
                           '-Rh/'+str(rotAngle)]

    
    p = subprocess.Popen(CommandArgumentList, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running RecFuncHarmonic on ", outURL, "errors?: \n", err, output
    
    
# 3.0 Harmonic Crude Stack
def runHarmonicRFstack(manifestURL, outURL, freqCut, lqtRot):
    
    CommandArgumentList = ['RecFuncHarmonic', '-L' + manifestURL,
                           '-F' + str(freqCut), '-T60', '-O' + outURL,
                           '-S1/100', '-H1', '-V1', '-R'+str(lqtRot)]

    
    p = subprocess.Popen(CommandArgumentList, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    #print "*** Running RecFuncHarmonic on ", outURL, "errors?: \n", err, output
    
# 3.1 Harmonic Crude Stack - With horizontal rotation [Jap OBS Data]
def runHarmonicRFstackWithRot(manifestURL, outURL, freqCut, rotAngle):
    
    CommandArgumentList = ['RecFuncHarmonic', '-L' + manifestURL,
                           '-F' + str(freqCut), '-T60', '-O' + outURL,
                           '-S1/100', '-H1', '-V1', '-R0.55', '-Rh/'+str(rotAngle)]

    
    p = subprocess.Popen(CommandArgumentList, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running RecFuncHarmonic on ", outURL, "errors?: \n", err, output
    

# 4.0 Azimuth Stack With Migration ...
def runRFstackByAzimWithMigrate(manifestURL, outURL, migrateURL, migrateDepth,
                                azimWinMin, azimWinMax, azimWinStep, freqCut, hCode, lqtRot):

    azimOption = str(azimWinMin) + '/' + str(azimWinMax) + '/'  + str(azimWinStep)
    
    CommandArgumentList = ['RecFuncHarmonic', '-L' + manifestURL,
                           '-F' + str(freqCut), '-T60', '-O' + outURL,
                           '-B1', '-S0/50', '-A' + azimOption, '-H'+str(hCode), '-V1', '-R'+str(lqtRot),
                           '-M0/'+str(migrateDepth), '-I'+migrateURL]

    
    p = subprocess.Popen(CommandArgumentList, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running RecFuncHarmonic on ", outURL, "errors?: \n", err, output
    
# 5.0 Harmonic Stack With Migration ...
def runHarmonicRFstackWithMigrate(manifestURL, outURL, migrateURL, migrateDepth, freqCut, hCode, lqtRot):
    
    CommandArgumentList = ['RecFuncHarmonic', '-L' + manifestURL,
                           '-F' + str(freqCut), '-T60', '-O' + outURL,
                           '-B1', '-S1/10', '-H'+str(hCode), '-V1', '-R'+str(lqtRot),
                           '-M0/'+str(migrateDepth), '-I'+migrateURL]

    
    p = subprocess.Popen(CommandArgumentList, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    #print "*** Running RecFuncHarmonic on ", outURL, "errors?: \n", err, output
    

# 6.0 Epicentral Stack With Migration
def runRFstackByEpicWithMigrate(manifestURL, outURL, migrateURL, migrateDepth, 
                     azimWinMin, azimWinMax, azimWinStep, epicWinMin, epicWinMax, epicWinStep, freqCut):
    

    epicOption = str(azimWinMin) + '/' + str(azimWinMax) + '/' + str(epicWinMin) + '/' + str(epicWinMax) + '/' + str(epicWinStep)
    
    CommandArgumentList = ['RecFuncHarmonic', '-L' + manifestURL,
                           '-F' + str(freqCut), '-T60', '-O' + outURL,
                          '-B2', '-S0/50', '-E' + epicOption, '-H1' , '-V1', '-R0.55',
                          '-M0/'+str(migrateDepth), '-I'+migrateURL]


    
    p = subprocess.Popen(CommandArgumentList, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running RecFuncHarmonic on ", outURL, "errors?: \n", err, output
    
# 7.0 run multiHK inversion    
def runRFmultiHKinversion(manifestURL, outURL, srchURL,freqCut ):
    
    CommandArgumentList = ['RFVelStck', '-L' + manifestURL,
                           '-F' + str(freqCut), '-T60', '-O' + outURL,
                           '-S' + srchURL, '-H1', '-V1', '-R1.0', '-I5']
    
    p = subprocess.Popen(CommandArgumentList, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running multiHK time domain inversion ", outURL, "errors?: \n", err, output

# 8.0 helper module used to run synthetic tests through 3 layer models
def parseSynthetic(noLayers, cSyn, cRFSyn, cHKTiming, rHKTime, rHKFreq):
    
    codeDir = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/Code/multiHKbenchmark/code/'
    
    bmarkScript = codeDir + 'benchmarkSyntheticMultiHK.sh'
    pltSynRFs = codeDir + 'plotSynRFs.sh'
    pltModelScript = codeDir + 'plotJPvelModel.sh'
    
    # Then run local 3 layer grid view inhouse ... depending on synthetic model to run ...

    
    # Choose Synthetic Model
    # 1. Generate Synthetic RF - Jeff's synthetic Code anirecsynth.f
    # 2. Run multiHK stack with Inverted Model over Synthetic Model
    # 3. Plot Snythetic Models - Synthetic RFs and Phase Timing
    # 4. Display Synthetic and Inverted Velocity models

    # Step 1 & 2
    scriptArgument = [bmarkScript,str(noLayers) , str(cSyn), str(cRFSyn), str(cHKTiming), str(rHKTime), 
                      str(rHKFreq) ]
    
    status = subprocess.call(scriptArgument)
    print "Synthetics and Inversion Complete: ", status
    
    # Step 3
    scriptArgument = [pltSynRFs,str(noLayers) ,  str(cHKTiming) ]
    
    status = subprocess.call(scriptArgument)
    print "View Synthetic RFs and Timing Here: ", status
    
    # Step 4
    scriptArgument = [pltModelScript,str(noLayers) ,  str(cHKTiming) ]
    
    status = subprocess.call(scriptArgument)
    print "View Optimal Inverted H-K Model Overlain with Synthetic Here: ", status

# 9.0 returns amplitude of RF stacks [harmonic] at zero time - used for multifreq. RF migration    
def returnZeroTimeConstRF(stackURL):
    
    Amplitudes = np.zeros(shape=(3,2))
    stackConst = ['MigrateHarmonic.Modelled.mean.xyz', 'MigrateHarmonic.UnModelled.mean.xyz']
    stackDip = ['MigrateHarmonic.Modelled.mean.rt', 'MigrateHarmonic.UnModelled.mean.rt']
    stackHor = ['MigrateHarmonic.Modelled.mean.rt', 'MigrateHarmonic.UnModelled.mean.rt']
    
    stackNames = np.array([stackConst, stackDip, stackHor])
    
    for istack in range(3):
        for itype in range(2):
            modeledStack = stackURL + stackNames[istack, itype]
            
            if istack == 0:
                columnTag = 40
            if istack == 1:
                columnTag = 0
            if istack == 2:
                columnTag = 2
    
            # read file for modeled stack and pick out just zero time and constant harmonic
            locFile = open(modeledStack, 'r');     
            locTable = locFile.readlines()

            #locZeroTimeAmp = []
            #print "length", len(locTable)
    
            lineFirst = locTable[0].rstrip().split("\t")
            lineSecond = locTable[1].rstrip().split("\t")
    
            #print "First line", lineSecond[0]
            deltaT = float(lineSecond[0]) - float(lineFirst[0])
            #print "Delta T", deltaT
    
            for line in locTable[0:len(locTable)]:
                linestrip = line.rstrip()
                columns = linestrip.split("\t")
        
                if len(columns) > 1 and int(columns[1]) == columnTag and abs(float(columns[0]) - 0.0) <= (0.5*deltaT):
                    #locZeroTimeAmp.append( [ float(columns[0]), float(columns[2]) ] )
                    Amplitudes[istack, itype] = float(columns[2])
                    #print "Zero Amplitude line: ", columns
                    break   
                   
    
    #return locZeroTimeAmp[0]
    return Amplitudes

# 10.0 saves 2-D multifrequency migrated RF stacks [xyz] file 
def saveHFstack(hGrid, fGrid, HFGrid, manifestHead, outDir, imgDir):
    #outDir = "/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/DepthFreqHFData/"
    
    stackConst = ['_Constant_HF.txt', '_Constant_res_HF.txt']
    stackDip = ['_Dip_HF.txt', '_Dip_res_HF.txt']
    stackHor = ['_Hor_HF.txt', '_Hor_res_HF.txt']
    
    stackNames = np.array([stackConst, stackDip, stackHor])
    

    
    for istack in range(3):
        for itype in range(2):  
            xyzFile = outDir + manifestHead + stackNames[istack,itype]
            xyzOut = open(xyzFile, 'w')
             
            for i in range(hGrid.shape[0]):
                for k in range(hGrid.shape[1]):
                    xyzOut.write(str(hGrid[i,k])+ "\t"+ str(fGrid[i,k]) + "00\t"+ str(HFGrid[istack,itype,i,k])+ "\n")
            
            xyzOut.close()
            
    print   "done saving file"     
    #loadPlotMultiFreqGrid(hGrid, fGrid, HFGrid, manifestHead, outDir, imgDir)
    #loadPlotMultiFreqGrid(xyBounds, manifestHead, outDir, imgDir)
    

# 10.0b plot xyz files in 10.0 above with color map
def loadPlotMultiFreqGrid(xyBounds, manifestHead, outDir, imagesDir):
    
    stackConst = ['_Constant_HF', '_Constant_res_HF']
    stackDip = ['_Dip_HF', '_Dip_res_HF']
    stackHor = ['_Hor_HF', '_Hor_res_HF']
    
    stackNames = np.array([stackConst, stackDip, stackHor])
    
    pngFile = imagesDir + manifestHead + "_Constant_HF.png"
    
    #hMin = str(hGrid[0,0])
    #hMax = str(hGrid[len(hGrid[:,0])-1, 0])
    hMin = str(xyBounds[0])
    hMax = str(xyBounds[1])
    
    #fMin = str(fGrid[0,0])
    #fMax = str(fGrid[0, len(hGrid[0,:])-1])
    
    fMin = str(xyBounds[2])
    fMax = str(xyBounds[3])
    
    #hstep = str(hGrid[1,1] - hGrid[0,1])
    #fstep = str(fGrid[0,1] - fGrid[0,0])
    
    hstep = str(xyBounds[4])
    fstep = str(xyBounds[5])
    
    print hstep, fstep, hMin, hMax, fMin, fMax
    
    #figs = plt.figure(figsize=(11,11))
    iCnt = 1
        
    
    im = []
    #gs1 = gridspec.GridSpec(3, 2)
    figs, ax = plt.subplots(3,2, sharex=True, sharey=True, figsize=(11,11))
    for istack in range(3):
        for itype in range(2):
            
            xyzFile = outDir + manifestHead + stackNames[istack, itype] + '.txt'
            grdFile = outDir + manifestHead + stackNames[istack, itype] + '.grd'
            
            xyz2grdPrompt = ['xyz2grd', xyzFile, '-G'+grdFile, '-I'+hstep+'/'+fstep,'-R'+hMin+'/'+hMax+'/'+fMin+'/'+fMax, '-N0']
            p = subprocess.Popen(xyz2grdPrompt, stdout=subprocess.PIPE) 
            output, err = p.communicate()
            print "*** Convert to grid ", "errors?: \n", output
            
            
            f = netcdf.netcdf_file(grdFile, version=2)
            x1 = f.variables['x'][:]
            y1 = f.variables['y'][:]
            z1 = f.variables['z'][:]
            X1,Y1 = meshgrid(x1,y1)
            
            valmax = abs(np.amax(z1))
            valmin = abs(np.amin(z1))
            # Load amplitudes for constant stack ...  and OTHERS?
            if (istack == 0 and itype == 0):
                
                # Clarify this logic ... Escpecially as it relates to just the last HFgrid? Not the isotropic?          
                valmax = abs(np.amax(z1))
                valmin = abs(np.amin(z1))
                
                valmax = np.amax([valmax,valmin])
                
                if (valmax == valmin):
                    valmax = abs(np.amax(z1))
                    valmin = -1*abs(np.amin(z1))
                    print "values", valmax, valmin
                    vStep  = (valmax - valmin) /5
                    ampLevs = np.arange(valmin, valmax, vStep)
                else:
                    valmin = -1 * valmax
                    vStep  = (valmax - valmin) /5
                    #print "values", valmax, valmin, z1
                    ampLevsA = np.arange(valmin, 0, vStep)
                    ampLevsB = np.arange(0, valmax, vStep/2)
                    
                    ampLevs = np.concatenate([ampLevsA, ampLevsB])
                
                vStep  = (valmax - valmin) / 4
                ampContour = ampLevs
                #print "amp Levels", ampLevs
                #print "contour levels", ampContour
            else:
                ampContour = [0.12, 0.18]
                levels = MaxNLocator(nbins=3).tick_values(z1.min(), z1.max())
                ampContour =levels

            
            #print x1, y1, z1
            # Define amplitude levels
            
            cmap = cm.coolwarm_r 
            cmap = cm.RdBu
            
            print "RF Levels", valmin, valmax
            print "amp Levels", ampLevs
            print "contour levels", ampContour
            
            
            ### NEW ADDITION ... Test new colormap seismic?            
            #ax = figs.add_subplot(3,2,iCnt)
            iCnt = iCnt + 1
            
            # toggle ampContour and ampLevs
                        
            im.append(ax[istack,itype].contourf(X1,Y1,z1,cmap=cm.RdBu,origin='lower',
            norm = colors.BoundaryNorm(ampLevs,cmap.N,clip=False)) )
            
            #ax[istack,itype].contour(X1,Y1,z1,linewidth=0.01,colors='grey', 
            #norm = colors.BoundaryNorm(ampContour,cmap.N,clip=False))
            
            ax[istack,itype].contour(X1,Y1,z1, ampContour, linewidth=0.01,colors='grey')


            
            ax[istack,itype].grid(which='major', color='white', linestyle='-')
            ax[istack,itype].set_xlim([int(hMin), int(hMax)])
            
            #if itype == 1:
            #    ax.yaxis.set_ticklabels([])
            #if istack < 2:
            #    ax.xaxis.set_ticklabels([])
    ax[2,0].set_xlabel('Depth (km)')
    ax[2,0].xaxis.set_label_coords(1.0, -0.15)
    
    ax[2,0].set_ylabel('Frequency (Hz)')
    ax[2,0].yaxis.set_label_coords(-0.15, 1.5)
    
    figs.subplots_adjust(right=0.84) 
    figs.subplots_adjust(wspace=0.05, hspace=0.05)
    #figs.tight_layout()   
    cax = figs.add_axes([0.85, 0.1, 0.03, 0.8])

    cb = figs.colorbar(im[0], cax=cax) 
    cb.ax.set_ylabel('RF avg. amp.', labelpad=-1)
    
    figs.suptitle(manifestHead, y= 0.94)  
    
    #ax.set_xlabel('Depth')
    
    #figs.colorbar(im[0])
     
    
   
    
    #imagesDir = "/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/RFImages/DepthFreqHFStacks/"
    plt.savefig(imagesDir+manifestHead+'.png', format='png', dpi=600)
    plt.savefig(imagesDir+manifestHead+'.ps', format='ps', dpi=600)
    
    plt.show()
    
# 11.0 plot single frequency Harmonic Stack With Migration ...
def plotHarmonicRFstackWithMigrate(manifestHead, pltScriptDir, stackURL, imageURL, migDepth, freqCut):
    
    # Harmonic Stack --- Display with shell script
    plotTitle = manifestHead+ ' Mig. to:' + str(migDepth) + 'km, @' + str(freqCut)  + 'Hz'


    CommandArgumentAzim = [pltScriptDir+'PlotRFHarmonicMig.sh', stackURL, imageURL, plotTitle]

    #Harmonic Stack --- Display with shell script
    p = subprocess.Popen(CommandArgumentAzim, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running Display Script on ", stackURL, "errors?: ", err, "Output: ", output

# 11.00 TRANSECT PLOT (multiple-station-on-line) single frequency Harmonic Stack With Migration - 
def plotTransectHarmonicWithMigrate(TransURL, TransName, TransDesc,pltScriptDir, stackURL, imageURL, migDepth, freqCut):
    
    # Construct Transect Plotting Script for Harmonic Stack 
    plotTitle = "Transect " +TransName+ ": " +TransDesc+ ' Mig. to:' + str(migDepth) + 'km, @' + str(freqCut)  + 'Hz'

    
    # This script doesn't yet exist but infrastructure already given ...
    CommandArgumentAzim = [pltScriptDir+'PlotTransectRFHarmonicMig.sh', TransURL+'.txt', imageURL, plotTitle]

    #Harmonic Stack --- Display with shell script
    p = subprocess.Popen(CommandArgumentAzim, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running Display Script 4 Transect ", stackURL, "errors?: ", err, "Output: ", output   

# 11.01 JAPAN - 5 TRANSECT PLOT (multiple-station-on-line) single frequency Harmonic Stack With Migration -
def plotAllTransectsHarmonicWithMigrate(TransURL, TransName, TransDesc,pltScriptDir, stackURL, imageURL, migDepth, freqCut, stckCde):
    
    # This script doesn't yet exist but infrastructure already given ...
    CommandArgumentAzim = [pltScriptDir+'PlotSummaryTransectRFHarmonicMig.sh', TransURL, imageURL, stckCde]

    #Harmonic Stack --- Display with shell script
    p = subprocess.Popen(CommandArgumentAzim, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    
    print "*** Running Display Script Summary Transects ", stackURL, "errors?: ", err, "Output: ", output 
    

# 11.0b plot single frequency Harmonic Stack - NO Migration ...    
def plotHarmonicRFstack(manifestHead, pltScriptDir, stackURL, imageURL, freqCut):
    
    # Harmonic Stack --- Display with shell script
    plotTitle = manifestHead+ ' No mig, @' + str(freqCut)  + 'Hz'


    CommandArgumentAzim = [pltScriptDir+'PlotRFHarmonic.sh', stackURL, imageURL, plotTitle]

    #Harmonic Stack --- Display with shell script
    p = subprocess.Popen(CommandArgumentAzim, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running Display Script on ", stackURL, "errors?: ", err, "Output: ", output
    
    
# 11.0 plot single frequency Azimuth Stack With Migration ...
def plotAzimRFstack(manifestHead, pltScriptDir, stackURL, imageURL, headCode):
    
    # Harmonic Stack --- Display with shell script
    plotTitle = manifestHead + '[H:' +  str(headCode) + ']'


    CommandArgumentAzim = [pltScriptDir+'PlotRFazim.sh', stackURL, imageURL, plotTitle]

    #Harmonic Stack --- Display with shell script
    p = subprocess.Popen(CommandArgumentAzim, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running Display Script on ", stackURL, "errors?: ", err, "Output: ", output
    
# 11.1 plot single frequency Azimuth Stack With Migration ...
def plotAzimRFstackWithMigrate(manifestHead, pltScriptDir, stackURL, imageURL, migDepth, freqCut):
    
    # Harmonic Stack --- Display with shell script
    plotTitle = manifestHead + ' Mig. to:' + str(migDepth) + 'km, @' + str(freqCut)  + 'Hz'


    CommandArgumentAzim = [pltScriptDir+'PlotRFazimMig.sh', stackURL, imageURL, plotTitle]

    #Harmonic Stack --- Display with shell script
    p = subprocess.Popen(CommandArgumentAzim, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running Display Script on ", stackURL, "errors?: ", err, "Output: ", output

# 12.1 plot single frequency Epic Stack With Migration ...
def plotEpicRFstackWithMigrate(manifestHead, pltScriptDir, stackURL, imageURL, migDepth, azimWinMin, azimWinMax, azimWinStep,freqCut):
    
    # Harmonic Stack --- Display with shell script
    plotTitle = manifestHead + ',' + str(migDepth) + 'km, @' + \
     '[' + str(azimWinMin) + ','+ str(azimWinMax)+ ','+ str(azimWinStep)+']'+ str(freqCut)  + ' Hz'


    CommandArgumentAzim = [pltScriptDir+'PlotRFepicMig.sh', stackURL, imageURL, plotTitle]

    #Harmonic Stack --- Display with shell script
    p = subprocess.Popen(CommandArgumentAzim, stdout=subprocess.PIPE) 
    output, err = p.communicate()
    print "*** Running Display Script on ", stackURL, "errors?: ", err, "Output: ", output

    