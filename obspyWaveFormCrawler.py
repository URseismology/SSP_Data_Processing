# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Author: Tolulope Olugboji
# Date:  April 11, 2014
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# Objective:   Helper modules for crawling waveforms [SAC etc.]
#              Used to read header, and waveform data etc.
#
#       
#
#

from obspy import UTCDateTime
from obspy import read
import obspy.signal
import time, sys, os
from IPython.display import display, clear_output
from obspy.taup.taup import getTravelTimes
from obspy.core.util.geodetics import gps2DistAzimuth

import numpy as np
import math
import matplotlib.pyplot as plt
import datetime
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
from numpy import sin, cos, exp, pi, arange, mean, array
from matplotlib.transforms import Affine2D

import matplotlib
matplotlib.rc('xtick', labelsize=16)

# Because we have hacked the transforms, you need a special method to
# set the voltage gain; this is a naive implementation of how you
# might want to do this in real life (eg make the scale changes
# exponential rather than linear) but it gives you the idea
def set_ygain(direction):
    set_ygain.scale += direction
    if set_ygain.scale <=0:
        set_ygain.scale -= direction
        return

    ax.set_ylim(0, set_ygain.scale)
    plt.draw()

# 1.0 load SAC headers into python dictionary...
def readHeaders(trace):
    "Read relevant header information ..."
    return {'gcarc': trace.stats.sac['gcarc'], 'evdp':trace.stats.sac['evdp'], 
            'stla': trace.stats['sac']['stla'], 'stlo': trace.stats['sac']['stlo'], 
            'evla': trace.stats['sac']['evla'], 'evlo': trace.stats['sac']['evlo'],
            'npts': trace.stats.npts,'baz':trace.stats.sac['baz'],
            'mag':trace.stats.sac['mag'], 'stTime':trace.stats.sac['b'], 'endTime': trace.stats.sac['e'],
            'delta': trace.stats['delta'], 'predTime': trace.stats.sac['t0'], 
            'predPhase': trace.stats.sac['kt0'], 'qcTime': trace.stats.sac['t1'],'eqTime':trace.stats['starttime'] }

# 2.0 return 3 component obsPy stream for fileName.[r,t,z]            
def read3CompStrm(fileName):
    fileHead = fileName.split("\n")
    stream = read(fileHead[0]+'R')
    stream += read(fileHead[0]+'T')
    stream += read(fileHead[0]+'Z')
    return stream

# 3.0 return n-by-3 two dimensional array of waveform data for each fileName in Manifest
def loadManifest(manifest, getTime=None):

    openManifest = open(manifest, 'r'); 
    manifestTable = openManifest.readlines()
    
    totNoEquakes = len(manifestTable)
    print "Total Equakes: ", totNoEquakes
    firstHead =  readHeaders(read3CompStrm(manifestTable[1])[0])
  
    dSze = (totNoEquakes, firstHead['npts'])
    vertWaveFm  = np.zeros(dSze)
    radWaveFm = np.zeros(dSze)
    transWaveFm = np.zeros(dSze)
    
    # use this to create time vector if requested ... New Hack
    Npts = firstHead['npts']
    stTime = firstHead['stTime']
    endTime = firstHead['endTime']
    
    timeSeq = np.linspace(stTime, endTime, Npts)
    
    ithWvFrm = 0
    
    # Timing info ...
    predTime = []
    qcTime = []
    predPhase = []
    eqInfo = []
    
    for fileName in manifestTable:
        
        #print fileName
        stream = read3CompStrm(fileName)
        #print stream[0].data
        #eqHeaders = readHeaders(stream[0])
        #print eqHeaders
        #print "len Data: ", len(stream[0].data)
        #print "npts: ", dSze
        
        radWaveFm[ithWvFrm,0:len(stream[0].data)] = stream[0].data
        transWaveFm[ithWvFrm,0:len(stream[1].data)] = stream[1].data
        vertWaveFm[ithWvFrm,0:len(stream[2].data)] = stream[2].data
        
        nextHdr = readHeaders(stream[0])
        predTime.append(nextHdr['predTime'])
        qcTime.append(nextHdr['qcTime'])
        predPhase.append(nextHdr['predPhase'])
        
        day = UTCDateTime(nextHdr["eqTime"]).strftime("%d. %B %Y %I:%M%p")
        mag = nextHdr['mag']
        eqInfo.append(day + ', Mag. '+str(mag))
        
        ithWvFrm = ithWvFrm + 1
    
    if getTime is None:
        return [radWaveFm, transWaveFm, vertWaveFm]
    else:
        timeInfo = [predTime, qcTime, predPhase, eqInfo]
        return [radWaveFm, transWaveFm, vertWaveFm], timeInfo, timeSeq

def loadManifestHeaders(manifest):
    evBaz = []
    evMag = []
    evDist = []
    
    openManifest = open(manifest, 'r'); 
    manifestTable = openManifest.readlines()
    
    totNoEquakes = len(manifestTable)
    
    for fileName in manifestTable:
        fileHead = fileName.split("\n")
        stream = read(fileHead[0]+'Z')
        
        eqHeaders = readHeaders(stream[0])
        
        evBaz.append(eqHeaders['baz'])
        evMag.append(eqHeaders['mag'])
        evDist.append(eqHeaders['gcarc'])
    
    #print "Total Equakes: ", totNoEquakes
    
    return [evBaz, evMag, evDist, totNoEquakes]

def getStationLoc(stationMetaDir, stationLocList):
    #stationMetaDir = "/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/StationMetaData/"
    #projects = ['Pacific.txt', 'Borehole.txt', 'OBS2005.txt', 'OBS2006.txt', 'OBS2007.txt']
    #projects = ['Pacific.txt']

    allCords = []

    for project in stationLocList:
        fName = stationMetaDir + project
        eventFile = open(fName, 'r'); 
        eventTable = eventFile.readlines()
        eventFile.close()

        stCords = []
            
        for line in eventTable:
            cols = line.split()    
            stCords.append( [float(cols[2]), float(cols[3]), float(cols[5])] )
    
        print "Project: ", project ,"No Stations: ", len(stCords)
        allCords.append(stCords)
    
    print "No of StationGroups: ", len(allCords)
    
    return allCords

# Use eventList file to cross reference 
def loadEventInfo(manifest, eventList):
    print "Cross-Referencing Manifests With Event List ..."
    
    evntLoc = []
    evntDep = []
    evntMag = []
    lYrs = []
    
    cntNoMatch = 0
    cntMatch = 0
    cntManifest = 0
    
    flag2005 = False
    flag2006 = False
    flag2007 = False
    flag2008 = False
    
    
    openManifest = open(manifest, 'r'); 
    manifestTable = openManifest.readlines()
    
    totManifest = len(manifestTable)
    
    for fileName in manifestTable[0:totManifest]:
    #for fileName in manifestTable[0:1]:
    
        
        fileNmeSplit = fileName.split(".")
        fileYearMonthDay =  fileNmeSplit[1]
        fileHrMinSec = fileNmeSplit[2]
        
        fileDate = fileYearMonthDay + fileHrMinSec[0:5]
        #print "Year Code on Manifest: ", fileYearMonthDay[0:2], fileYearMonthDay
        
        openEvents = open(eventList, 'r'); 
        evntTable = openEvents.readlines()
        allEvnts = len(evntTable)
        
        flagNoMatch = True
        cntManifest = cntManifest + 1
        
        for eventLine in evntTable[1:allEvnts]:
        
            evntYr = eventLine.split("\t")[0]
            evntMnth = eventLine.split("\t")[1]
            evntDay = eventLine.split("\t")[2]
            evntHr = eventLine.split("\t")[3]
            evntMin = eventLine.split("\t")[4]
            evntSec = eventLine.split("\t")[5]
            
            evntSec = str(int(float(evntSec)))
            evntSecRnd =  str(round(float(evntSec), -1))
            
            evLon = eventLine.split("\t")[7]
            evLat = eventLine.split("\t")[8]
            evMag = eventLine.split("\t")[10] 
            evDep = eventLine.split("\t")[9]
              
            
            #print "sec ", evntSec
            
            
            if (len(evntMnth) == 1):
                evntMnth = '0'+evntMnth
            if(len(evntDay) == 1):
                evntDay  = '0'+evntDay
            if(len(evntHr) == 1):
                evntHr = '0'+evntHr
            if(len(evntMin) == 1):
                evntMin = '0'+evntMin
            if(len(evntSec) == 1):
                evntSec = '0'+evntSec
            
            if (float(evMag) == 0.0):
                evMag = '5.0'
            
            
            evntDate = evntYr[2:4] + evntMnth + evntDay+ evntHr + evntMin + evntSec[0:1]
            evntDate2 = evntYr[2:4] + evntMnth + evntDay+ evntHr + evntMin + evntSecRnd[0:1]
            
            #print len(eventLine.split("\t")), evMag
            
            if (fileDate == evntDate or fileDate == evntDate2):
                #print  "Match", fileDate, evntDate, evLon, evLat, evMag, evDep
                evntLoc.append( [float(evLon), float(evLat)] )
                evntDep.append(float(evDep))
                evntMag.append(float(evMag))
                cntMatch = cntMatch + 1
                flagNoMatch = False
                break
            #else:
            #    print "No Match", fileDate, evntDate, len(fileDate), len(evntDate)
            
        if (flagNoMatch):
            print "No Match", cntManifest, "of", totManifest,  fileName, fileDate, evntDate, len(fileDate), len(evntDate)
            cntNoMatch = cntNoMatch + 1
        
        if(flag2005 == False and fileYearMonthDay[0:2] == '05'):
            flag2005 = True
        if(flag2006 == False and fileYearMonthDay[0:2] == '06'):
            flag2006 = True
        if(flag2007 == False and fileYearMonthDay[0:2] == '07'):
            flag2007 = True
        if(flag2008 == False and fileYearMonthDay[0:2] == '08'):
            flag2008 = True
                                
    print "All Match", cntMatch, "No Match", cntNoMatch, "tot Manifest", totManifest
    if(flag2005 == True ):
       lYrs.append('2005')
    if(flag2006 == True ):
       lYrs.append('2006')
    if(flag2007 == True ):
       lYrs.append('2007')
    if(flag2008 == True ):
       lYrs.append('2008')
       
        #print fileName
        #print "File Stub", fileDate
        #print "Event Stup", evntDate
    #print "Equakes : ", totNoEquakes
    
    return evntLoc, evntDep, evntMag, lYrs
    
def getStationName(stationMetaDir, stationLocList):
    allNames = []

    for project in stationLocList:
        fName = stationMetaDir + project
        eventFile = open(fName, 'r'); 
        eventTable = eventFile.readlines()
        eventFile.close()

        stNames = []   
        for line in eventTable:
            cols = line.split()    
            stNames.append( cols[1] )
    
        print "Project: ", project ,"No Stations: ", len(stNames)
        allNames.append(stNames)
    
    print "No of StationGroups: ", len(allNames)
    
    return allNames
    

def getDistance(statA, statB):
    latA = statA[0]; lonA = statA[1]
    latB = statB[0]; lonB = statB[1]
    
    #print latA, latB, lonA, lonB
    #print statA, statB
    
    dist, az, baz = gps2DistAzimuth( latA , lonA, latB, lonB )
    distKm = dist/1000
    
    #print dist, az, baz
    return distKm
     
def sortByLocation(stationGp, stationNm):
    
    locationTol = 50; # tolerance is 50 km
    stationPass = stationGp[0]
    namePass = stationNm[0]
    
    uniquePass = []
    uniqueNme = []
    
    stationCnt = [1] * len(stationPass) 
    #print "First station group", stationPass
    
    # check other station groups against 1st
    for indxGp in range(1,len(stationGp)):
        # Scan stations in current group
        if (len(uniquePass) > 0):
            stationPass = stationPass+uniquePass
            #print indxGp, stationPass, len(stationPass)
            
        for indxSt in range(len(stationGp[indxGp])):
            stLonLatA =  stationGp[indxGp][indxSt] 
            nme = stationNm[indxGp][indxSt]
            
            updateStation = True
            # Check against 1st station group
            for indxSt1 in range(len(stationPass)):
                
                stLonLatB = stationPass[indxSt1]
                dist = getDistance(stLonLatA, stLonLatB)
                #print "Gp", indxGp, "stat: ", indxSt, stLonLatA, "dist2A", dist, stLonLatB
                
                
                    
                if (dist <= locationTol):
                    updateStation = False
                    stationCnt[indxSt1] =  stationCnt[indxSt1] + 1
                    #print nme
                    #print "Gp", indxGp, "stat: ", indxSt, stLonLatA, "dist2A", dist
                    
                    
            if(updateStation):
                stationCnt.append(1)
                uniquePass.append(stLonLatA)
                uniqueNme.append(nme)
                #print "Unique station", stLonLatA
        
    namePass = namePass+uniqueNme   
    
    # Display Year, Station, Statistics     
    for iLen in range(len(stationCnt)):
        print namePass[iLen], ",",stationCnt[iLen]
    
    #print stationCnt , "number of stations", len(stationCnt)      
    #print "station name: ", namePass, len(namePass)     

def plot3CmpTrace(AllTraces, trIndx, titleInfo='', timeInfo=None,timeSeq=None,   fName=None, format='ps'):
    # Code to plot 3 component Trace for publication
    # AllTraces - 3-D array holding rad, trans, vert waveforms for all manifest..
    # indx - trace index for display ...
    # titleInfo: title for the plot - station, magnitude, depth bla bla bla
    # timeTag: timeTag for timing indication 
    # qcTimeTag: timeTag used by QC
    # tagPhase: character description for earthquake phase e.g PkIkP
    print "work in progress"
    predTime = timeInfo[0][trIndx]
    qcTime = timeInfo[1][trIndx]
    predPhase = timeInfo[2][trIndx]
    eqInfo = timeInfo[3][trIndx]
    
    figs = plt.figure(figsize=(11,6))
    
    colors = ['blue', '#1E90FF', 'red']
    for iCmp in range(3):
        # Pick rad, trans, or vert..
        sig = AllTraces[iCmp].tolist()
        
        # Normal particular trace for display use trIndx to identify trace... then bias by iCmp
        pltSig = np.array(sig[trIndx])/ np.max(np.abs(np.array(sig[trIndx]))) + (1.5)*iCmp
        
        if timeSeq is None:
            plt.plot(pltSig)
        else:
            plt.plot(timeSeq, pltSig, color= colors[iCmp])
            plt.axvline(x=predTime, color='grey')
            plt.axvline(x=qcTime, color='black', linewidth=2)
    
    plt.yticks([])
    plt.title(titleInfo + ', ' + eqInfo, fontsize=20)
    plt.xlabel('Time (sec)', fontweight = 'bold', fontsize=16)   
    
    if fName is None:
        print "file not saved. No file URL provided"
    else:
        print "saving file ..", fName+'.'+format
        plt.savefig(fName+'.'+format, format=format)
    plt.show()