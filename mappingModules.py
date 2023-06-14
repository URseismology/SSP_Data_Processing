from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import scipy.ndimage
import time, sys, os
from netCDF4 import Dataset

import matplotlib.mlab as ml

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker

from obspy.core.util.geodetics import gps2DistAzimuth, crossTrackDistance

import matplotlib as mpl
mpl.rcParams['font.size'] = 15
mpl.rcParams['axes.labelsize'] = 'large'

# Code to cross-reference with ocean age digital data - 
# Inherited from continental geotherm work from KaratoEtAl_2014
def findAge(lon, lat, inAgeFile, tol=50):
    
    AgeFile = open(inAgeFile, 'r'); 
    AgeTable = AgeFile.readlines()
    
    for ageVal in AgeTable[1:len(AgeTable)]:
        columns = ageVal.split()
        #print "Columns ", columns[0], columns[1], columns[2]
        
        tempLon = float( columns[0] )
        tempLat = float( columns[1] )
        age = float ( columns[2] )
    
        
        
        #print "Pos. File" , tempLon, tempLat, "Pos. Station:", lon, lat
        dist, az, baz = gps2DistAzimuth( float( lat) , float( lon ), tempLat, tempLon )
        #print "Pos. File" , tempLon, tempLat, "Pos. Station:", lon, lat, "Distance: ", dist
        distKm = dist/1000
        
        try:
                
            if (distKm < tol  ):
                print "Found!!"
                print "Dist: ", distKm, "Other Data", columns
                
                print "lon, lat ", lon, lat,  "age: ", age
                return age
                break
                
        except ValueError:
            print "Not a float"
            break
    
    if (tempMin < 1000):
        print  "Not Found: ", lon, lat
    
    return 0

def returnHotspotLocations(fName):
    # read hotspots.gmt
    hotspotFile = open(fName, 'r');     
    hotspotTable = hotspotFile.readlines()
    
    hotSpotCords = []
    
    for line in hotspotTable[1:len(hotspotTable)]:
        #clear_output()
        columns = line.split("\t")

        # append - lat, lon   
        #print len(columns), columns
        if ( len(columns) > 2 ):
            hotSpotCords.append( [  float(columns[0]), float(columns[1]) ] )
        
    return hotSpotCords


def normVector(vector):
    mag =  math.sqrt(vector[0]**2 + vector[1]**2)
    nmVector = vector / mag
    return nmVector


def returnAngle(vectorA, vectorB):
    normA = normVector(vectorA)
    normB = normVector(vectorB)
    
    dotProd = math.fabs(np.dot(normA, normB))
    angle = np.arccos(dotProd) * (180 / np.pi)
    return angle
    
def buildVectorGrid():
    # geographical locations ...staticOrientation[0] = lon ; staticOrientation[1] = lat 
    staticOrientation = [{'JOHN':-169.5292, 
                          'RAR':-159.7733,
                          'KIP':-158.0112,
                          'XMAS':-157.4457,
                          'POHA':-155.5326,
                          'PTCN':-130.0953,
                          'RPN':-109.3344,
                          'PAYG':-90.2861,
                          'WAKE':166.652,
                          'KWAJ':167.613,
                          'TARA':172.9229, }, 
                     
                         {'JOHN':16.7329, 
                          'RAR':-21.2125,
                          'KIP':21.42,
                          'XMAS':2.0448,
                          'POHA':19.7573,
                          'PTCN':-25.0713,
                          'RPN':-27.1267,
                          'PAYG':-0.6742,
                          'WAKE':19.2834,
                          'KWAJ':8.8019,
                          'TARA':1.3549, }
                     
                             ]
    dipOrientation = {'JOHN':[ [0.120407,-27.6038],	[0.0464743,-12.1245],	0.0249996,	3.25,	1,	1 ],
                          'RAR':[ [0.00001,0.0001],	[0.0320044,148.124],	1.45,	2.0,	1,	1 ],
                          'KIP':[ [0.0637498,-22.6283],	[0.050362,-179.716],	1.25,	2.225,	1,	1 ],
                          'XMAS':[ [0.124448,98.8942],	[0.0687762,-35.6039],	0.525,	4.0,	1,	1 ],
                          'POHA':[ [0.0632611,-137.334],	[0.0614051,72.8078],	1.475,	2.525,	1,	1 ],
                          'PTCN':[ [0.0734992,7.92846],	[0.0582371,159.558],	0.0249996,	2.025,	1,	1 ],
                          'RPN':[ [0.0806555,-44.1168],	[0.00001,0.0001],	0.65,	4.0,	1,	1 ],
                          'PAYG':[ [0.101381,-18.1291],	[0.0241671,39.207],	0.7,	2.95,	1,	1 ],
                          'WAKE':[ [0.00001,0.0001],	[0.0404969,165.284],	1.95,	2.45,	1,	1 ],
                          'KWAJ':[ [0.100257,-96.2637],	[0.0804664,133.365],	0.0249996,	2.025,	1,	1 ],
                          'TARA':[ [0.14015,-107.424],	[0.136934,85.1776],	0.775,	2.475,	1,	1 ],
                          }

    # Plate Vectors [0] Grips and Gordon [E, North] - Normalize for display
    Nuvel1A = {'JOHN':[ [-96,55] ],
                          'RAR':[ [-103,53] ],
                          'KIP':[[-89,52] ],
                          'XMAS':[ [-103,52] ],
                          'POHA':[[-90,52]  ],
                          'PTCN':[ [-112,36]  ],
                          'RPN':[ [31,-9]  ],
                          'PAYG':[ [21,0] ],
                          'WAKE':[[-102,55] ],
                          'KWAJ':[ [-104,55] ],
                          'TARA':[ [-104,56] ],
                          }
    
    # Shear wave splitting results from various authors [dt, theta] - vectorize for display
    # Results culled from shear wave splitting database IRIS
    # JOHN - (2003) Walker, K.T.; Bokelmann, G.H.R.; Klemperer, S.L.
    # RAR - (2007) Fontaine, F.R.; Barruol, G.; Tommasi, A.; Bokelmann, G.H.R.
    # KIP - 1.(1992) Vinnik, L.P., Makeyeva, L.I., Milev, A., Usenko, Y.
    #       2.(1999) Barruol, G., Hoffman, R.
    #       3.(2001) Walker, K. T., Bokelmann, G. H. R. and Klemperer, S.L.
    #       4.(2012) Collins, J.A., Wolfe, C.J. and Laske, G.
    # XMAS - None
    # POHA - 1.(2001) Walker, K. T., Bokelmann, G. H. R. and Klemperer, S.L.
    #        2.(2012) Collins, J.A., Wolfe, C.J. and Laske, G.
    # PTCN - 1.(1998) Russo, R.M.; Okal, E.A.
    #     - 2.(2007) Fontaine, F.R.; Barruol, G.; Tommasi, A.; Bokelmann, G.H.R.
    # RPN  - (2007) Fontaine, F.R.; Barruol, G.; Tommasi, A.; Bokelmann, G.H.R.
    # PAYG  - (2005) Fontaine, F. R., Hooft, E. E. E., Burkett, P. G., Toomey, D. R., Solomon, S. C., and Silver, P. G
    # WAKE - None
    # KWAJ - None
    # TARA - None
    SplitReslts = {'JOHN':[ [1.1,103.0] ],
                          'RAR':[ [1.71,297.0] ],
                          'KIP':[[1.5,45.0],[1.08,45.0],[1.5,60.0],[0.92,80.0] ],
                          'XMAS':[ [0.0,0.0] ],
                          'POHA':[[0.9,308.0], [0.55,315.0] ],
                          'PTCN':[ [1.1,322.0], [1.09,276.0]  ],
                          'RPN':[ [1.37,33.0]  ],
                          'PAYG':[ [0.0,0.0] ],
                          'WAKE':[[0.0,0.0] ],
                          'KWAJ':[ [0.0,0.0] ],
                          'TARA':[ [0.0,0.0] ],
                          }
    SplitReslts = {'JOHN':[ [1.1,103.0] ],
                          'RAR':[ [1.71,297.0] ],
                          'KIP':[[1.5,45.0],[1.08,45.0],[1.5,60.0],[0.92,80.0] ],
                          'XMAS':[ [0.0,0.0] ],
                          'POHA':[[0.9,308.0], [0.55,315.0] ],
                          'PTCN':[ [1.1,322.0], [1.09,276.0]  ],
                          'RPN':[ [1.37,33.0]  ],
                          'PAYG':[ [0.0,0.0] ],
                          'WAKE':[[0.0,0.0] ],
                          'KWAJ':[ [0.0,0.0] ],
                          'TARA':[ [0.0,0.0] ],
                          }
                          
    srcCnt = []
    for key in SplitReslts:
        srcCnt.append(len(SplitReslts[key]))

    #print srcCnt
    #print max(srcCnt)
    spltLen  = max(srcCnt)

    # Station crust age[0] and sea-floor age[1]
    CrustSeaAge = {'JOHN':[ 71, 118.7 ],
                          'RAR':[ 1.4, 96.1],
                          'KIP':[5, 89.0 ],
                          'XMAS':[ 36, 116.3 ],
                          'POHA':[ 0.4, 91.4 ],
                          'PTCN':[ 0.6, 24.9 ],
                          'RPN':[ 2, 5.9  ],
                          'PAYG':[ 1,11.1 ],
                          'WAKE':[ 90, 165.9],
                          'KWAJ':[ 88, 156.1 ],
                          'TARA':[ 88, 141.1 ],
                          }
                          
    ### Use Orientation Data to Build Grid ...
    #sorted coordinates ...
    lon = [-169.5292, -159.7733, -158.0112, -157.4457, -155.5326, -130.0953, 
           -109.3344, -90.2861, 166.652, 167.613, 172.9229]

    lat = [-27.1267, -25.0713, -21.2125, -0.6742, 1.3549, 2.0448, 
           8.8019, 16.7329, 19.2834, 19.7573, 21.42]

    lenlon = len(lon)
    lenlat = len(lat)


    # Initialize plate vectors
    ugridPlate = np.zeros(shape=(lenlon, lenlat))
    vgridPlate =  np.zeros(shape=(lenlon, lenlat))

    # Initialize top crust grid to zeros
    ugrid = np.zeros(shape=(lenlon, lenlat))
    vgrid = np.zeros(shape=(lenlon, lenlat))

    # Initialize low crust grid to zeros
    ugrid2 = np.zeros(shape=(lenlon, lenlat))
    vgrid2 = np.zeros(shape=(lenlon, lenlat))
    
    # Define grid for split results.
    nUgridSplit = np.empty(shape=(spltLen, lenlon, lenlat)); nUgridSplit[:] = np.NAN
    nVgridSplit = np.empty(shape=(spltLen, lenlon, lenlat)); nVgridSplit[:] = np.NAN

    # parse stations 
    cAge = []
    sAge = []
    prAngTop = []
    prAngBot = []
    stationName = []
    
    stations = ['JOHN', 'RAR', 'KIP', 'XMAS', 'POHA', 'PTCN', 'RPN', 'PAYG', 'WAKE', 'KWAJ', 'TARA']
    
    for station in stations:
        stalon = staticOrientation[0][station]
        stalat = staticOrientation[1][station]
    
        vector1 = dipOrientation[station][0]
        vector2 = dipOrientation[station][1]
        time1 =  dipOrientation[station][2]
        time2 =  dipOrientation[station][0]
        iSpltLen = len(SplitReslts[station])
    
        for ilon in range(lenlon):
            if stalon == lon[ilon]:
                for ilat in range(lenlat):
                    if stalat == lat[ilat]:
                    
                    
                        # Plate..
                        North = Nuvel1A[station][0][1] 
                        East = Nuvel1A[station][0][0]
                    
                        plateVec = normVector(np.array([East,North]))                       #Plate VectorHere!!!
                    
                        #normNorth = North / math.sqrt(North**2 + East**2)
                        #normEast = East / math.sqrt(North**2 + East**2)
                    
                        vgridPlate[ilat, ilon]= 10* plateVec[1]
                        ugridPlate[ilat,ilon]= 10 * plateVec[0]
                    
                        thetaRad = math.radians(vector1[1])
                        Rscaled = 100.0 * vector1[0]
                    
                        # top crust
                        ugrid[ilat, ilon] = Rscaled*np.cos(thetaRad)
                        vgrid[ilat, ilon] = Rscaled*np.sin(thetaRad)
                    
                        aniTopCrust = normVector(np.array([Rscaled*np.cos(thetaRad),Rscaled*np.sin(thetaRad)]))
                        predAngleTop = returnAngle(plateVec, aniTopCrust)
                    
                        # bot crust
                        thetaRad = math.radians(vector2[1])
                        Rscaled = 100.0 * vector2[0]
                    
                        ugrid2[ilat, ilon] = Rscaled*np.cos(thetaRad)
                        vgrid2[ilat, ilon] = Rscaled*np.sin(thetaRad)
                        
                        # split results..
                        for ilen in range(iSpltLen):
                            vec = SplitReslts[station][ilen]
                            dt = vec[0]
                            
                            # convert to radians
                            vec[1] = 90 - vec[1]
                            if vec[1] < - 180:
                                vec[1] = 360+vec[1]
                                
                            theta = math.radians(vec[1])
                            
                            nUgridSplit[ilen,ilat, ilon] = 10 * dt*np.cos(theta)
                            nVgridSplit[ilen,ilat, ilon] = 10 * dt*np.sin(theta)
                            
                    
                        aniBotCrust = normVector(np.array([Rscaled*np.cos(thetaRad),Rscaled*np.sin(thetaRad)]))
                        predAngleBot = returnAngle(plateVec, aniBotCrust)
                    
                        #print "Station:", station, stalon, ilon, ilat
                        print "Station:", station, "Pred. Angle, Top Crust: [", predAngleTop, predAngleBot, "]", \
                        "Crust Age: ",CrustSeaAge[station][0]
                    
                        cAge.append(CrustSeaAge[station][0])
                        sAge.append(CrustSeaAge[station][1])
                        prAngTop.append(predAngleTop)
                        prAngBot.append(predAngleBot)
                        stationName.append(station)
                    else:
                        # Plate..
                        vgridPlate[ilat, ilon]=np.nan
                        ugridPlate[ilat,ilon]=np.nan
                    
                        # top crust
                        vgrid[ilat, ilon]=np.nan
                        ugrid[ilat,ilon]=np.nan
                    
                        # bot crust
                        vgrid2[ilat, ilon]=np.nan
                        ugrid2[ilat,ilon]=np.nan
                        
                        # split results..
                        for ilen in range(iSpltLen):
                            nUgridSplit[ilen,ilat, ilon] = np.nan
                            nVgridSplit[ilen,ilat, ilon] = np.nan
       
    # transofrm to array ..
    arlon = np.array(lon)
    arlat = np.array(lat)
    
    resultGrid = [[ugridPlate,vgridPlate], [ugrid, vgrid], [ugrid2, vgrid2]]
    splitGrid = [nUgridSplit, nVgridSplit]
    
    resultLoc = [arlon, arlat]
    
    return resultGrid, splitGrid, resultLoc, CrustSeaAge   


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
            stCords.append( [float(cols[2]), float(cols[3]),  float(cols[5])] )
    
        print "Project: ", project ,"No Stations: ", len(stCords)
        allCords.append(stCords)
    
    print "No of StationGroups: ", len(allCords)
    
    return allCords
    
def getStationLocPlusName(stationMetaDir, stationLocList, ageFileLkUp=None):
    # Here station info is retrieved with name.
    # Write code that uses location to retrieve age?
    #
    allCords = []

    for project in stationLocList:
        fName = stationMetaDir + project
        eventFile = open(fName, 'r'); 
        eventTable = eventFile.readlines()
        eventFile.close()

        stCordsNme = []
            
        for line in eventTable:
            cols = line.split()    
            #print cols
            
            if ageFileLkUp is None :
                stCordsNme.append( [float(cols[2]), float(cols[3]), cols[1]] )
            else:
                print "No Age Look up requested ..."
    
        print "Project: ", project ,"No Stations: ", len(stCordsNme)
        allCords.append(stCordsNme)
    
    print "No of StationGroups: ", len(allCords)
    
    return allCords

def getStationDistOnTransect(locAnotPair, AATransect, tolerance):
    statOnTran = []     #Stations on transect
    
    begTranLat = AATransect[1]; begTranLon = AATransect[0]
    endTranLat = AATransect[3]; endTranLon = AATransect[2]
    
    for indxGp in range(len(locAnotPair)):
        for indxSt in range(len(locAnotPair[indxGp])):
            stLatLonName =  locAnotPair[indxGp][indxSt] 
            
            dxt, dat = crossTrackDistance(begTranLat, begTranLon, endTranLat, endTranLon,
             stLatLonName[0], stLatLonName[1])
             
            if (np.abs(dxt) <= tolerance):
                #print stLatLonName[2], dxt, dat
                statOnTran.append([stLatLonName[2], dxt, dat])
            #print "case:", stLatLonName
            
    
    return statOnTran
    
def plotRegionalMap(mapBounds, stationGp, imagesFileDir,  locAnotPair=None, transectLine=None, EquakeLoc=None):
    
    if locAnotPair is None:
        print "Annotation not set"
    else:
        print "Will Print Annotation on Map"
        
    
    # Region = [lowerLeftlonA, latA, upperRightlonB, latB, centLon, centLat]
    #Hawaii = [-179, -5, -138, 40, 180, 20  ]        #Hawaii
    #Hawaii = [135, -10, 180, 30, 160, 20]          #NwPacific
    Hawaii = [-179, -30, -100, 00, -120, -10]        #superSwell
    #Hawaii = [-111, -11, -75, 21, -91, 20] #Galapagos
    
    lonLeft = mapBounds[0]; lonRight = mapBounds[2]; lonStep=20
    latLeft = mapBounds[1]; latRight = mapBounds[3]; latStep=10
    fzHawaiiDir = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/PlateBoundaryShapeInfo/FractureZones/SHP/GSFML_SF_FZ_RM'
    fzHawaiiFile = 'GSFML_SF_FZ_RM'
    
    #ridges
    plateBoundaries = "/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/PlateBoundaryShapeInfo/PLATES_PlateBoundary_ArcGIS/ridge"
    plateFile = 'ridge'
    
    # Plot topography and sea-floor ...
    
    # read in etopo5 topography/bathymetry.
    #url = 'http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/etopo5.nc'
    url = '/Users/tmo22/Documents/OlugbojiKarato2011/AgeGridsMuller/Data/etopo5.nc'
    etopodata = Dataset(url)
    
    topoin = etopodata.variables['ROSE'][:]
    lons = etopodata.variables['ETOPO05_X'][:]
    lats = etopodata.variables['ETOPO05_Y'][:]

    # read in sea floor age
    url = '/Users/tmo22/Documents/OlugbojiKarato2011/AgeGridsMuller/Data/age.3.6.nc'
    agedata = Dataset(url)
    
    agein = agedata.variables['z'][:] / 100
    agelons = agedata.variables['x'][:]
    agelats = agedata.variables['y'][:]


    # shift data so lons go from -180 to 180 instead of 20 to 380.
    topoin,lons = shiftgrid(180.,topoin,lons,start=False)
    agein, agelons = shiftgrid(180.,agein,agelons,start=False)

    plt.figure(figsize=(11,11))
    # lon_0, lat_0 are the center point of the projection.
    # resolution = 'l' means use low resolution coastlines.
    # create polar stereographic Basemap instance.
    m = Basemap(projection='lcc',\
            llcrnrlat=latLeft,urcrnrlat=latRight,\
            llcrnrlon=lonLeft,urcrnrlon=lonRight, lat_0=mapBounds[5], lon_0=mapBounds[4], lat_1=latLeft)
    
    # transform to nx x ny regularly spaced 5km(2km) native projection grid
    nx = int((m.xmax-m.xmin)/5000.)+1; ny = int((m.ymax-m.ymin)/5000.)+1

    topodat = m.transform_scalar(topoin,lons,lats,nx,ny, masked=True)
    agedat, agex, agey =  m.transform_scalar(agein,agelons,agelats,nx,ny, masked=True, returnxy= True)
        
    m.drawcoastlines()
    m.drawparallels(np.arange(latLeft, latRight, latStep), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(lonLeft, lonRight, lonStep), labels=[0, 0, 0, 1], linewidth=0)
    
    
    #1.  Topography shade and contour...
    clevSwell = np.arange(-4100, -4000, 5)
    clevLand = np.arange(-100, 100,5)
    clevs = np.append(clevSwell, clevLand)
    
    im = m.imshow(topodat, cmap=cm.get_cmap('gray'), zorder=0)

    #2. Contour
    x, y = m(*np.meshgrid(lons,lats))
    m.contour(x, y, topoin,clevs, interpolation='guassian', linewidths=0.1,colors='white', zorder=2)
    
    discLevs = [0, 9.7, 20.1, 33.1, 40.1, 47.9, 55.9, 67.7, 83.5, 120.4, 126.7, 131.9, 139.6, 147.7,
            154.3, 180]

    #discLevs =  [120.4, 126.7, 131.9, 139.6, 147.7,154.3, 180]

    cmap = plt.cm.jet_r

    #3. Age colormap
    im = m.imshow(agedat,cmap=cm.get_cmap('jet_r'), alpha = 0.4, zorder=1, 
                  norm = colors.BoundaryNorm(discLevs,cmap.N,clip=False))
                  
    #4. Plot Fracture Zones... & Ridges ...
    shp_info_Hawaii = m.readshapefile(fzHawaiiDir,fzHawaiiFile,drawbounds=False)
    for shapedict, shape in zip(m.GSFML_SF_FZ_RM_info, m.GSFML_SF_FZ_RM):
        xx,yy = zip(*shape)
        m.plot(xx,yy,linewidth=2,color='k', zorder=1)
        
    shp_info = m.readshapefile(plateBoundaries,plateFile, drawbounds=False)
    for shapedict, shape in zip(m.ridge_info, m.ridge):
        desc = shapedict['geogdesc']
        #print desc[0:7]
        if desc.find('PACIFIC') != -1 or desc.find('EXPLORER') != -1 or desc.find('NAZCA') != -1 :
            xx,yy = zip(*shape)
            m.plot(xx,yy,linewidth=1.5,color='white')
        
    # Place hotspot locations on map
    #hotspotFileDir = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/CrustalAgeData/'
    #hotspotLocs = returnHotspotLocations(hotspotFileDir+'hotspots.gmt')
    #hotspotLonLat = np.array( hotspotLocs )
    #xhot, yhot = m(hotspotLonLat[:,0], hotspotLonLat[:,1])
    #m.scatter(xhot, yhot, marker='o', s=500, c= 'red', alpha=1, zorder= 4)
    
    #4. Plot station locations - single group, then test all groups.
    markType = ['*', 'o', '^']
    markColor = ['red', 'white', 'black']
    markEdge = ['black', 'black', 'white']
    
    for indxGp in range(len(stationGp)):
        for indxSt in range(len(stationGp[indxGp])):
            stLonLat = np.array( stationGp[indxGp][indxSt] )
            xsta, ysta = m(stLonLat[:,1], stLonLat[:,0])
            
            #if indxGp == 2:
            #print "stations: ", stLonLat
            yrColIndx = stLonLat[:,2] -2
            #print "No. of yrs operational:", yrColIndx
            #else:
            #   yrColIndx = indxGp
            for ithSta in range(len(stationGp[indxGp][indxSt])):
                m.scatter(xsta[ithSta], ysta[ithSta], marker=markType[indxGp], s=250, c= markColor[int(yrColIndx[ithSta])], alpha=1, 
                edgecolor = markEdge[int(yrColIndx[ithSta])], zorder= 4)
    
    #5.1 Annotate plot with station name using python text command
    for indxGp in range(len(locAnotPair)):
        for indxSt in range(len(locAnotPair[indxGp])):
            stLonLat = np.array( locAnotPair[indxGp][indxSt] )
            #print stLonLat
            xsta, ysta = m( float(stLonLat[1]), float(stLonLat[0])+0.5)
            name = stLonLat[2]
            plt.annotate(name, xy=(xsta, ysta), fontsize=14, fontweight='bold', zorder= 4, alpha=1)
            
    #5.2 Plot transects on station, visualize lateral coverage ..
    if transectLine is None:
        print "No transects defined"
    else:
        print "plotting transect as great circle"
        for indxTrans in range(len(transectLine)):
            #print transectLine[indxTrans]
            lonA = transectLine[indxTrans][0]; latA = transectLine[indxTrans][1] 
            lonB = transectLine[indxTrans][2]; latB =transectLine[indxTrans][3] 
        
            m.drawgreatcircle(lonA, latA, lonB, latB, color='white')
            
    #5.3. Plot Earthquakelocations - single group, then test all groups.
    if EquakeLoc is None:
        print "No Equake Data Set"
    else:
        print "plot equake data "
        
        for indxGp in range(len(EquakeLoc)):
                stLonLat = np.array( EquakeLoc[indxGp] )
                xsta, ysta = m(stLonLat[:,1], stLonLat[:,0])
                m.scatter(xsta, ysta, marker='o', s=40, c= 'black', 
                edgecolor = 'black', alpha=1, zorder= 3)
    
    m.drawmapboundary(fill_color='aqua')    
    
    # save figures ..
    plt.savefig(imagesFileDir, format='png', dpi=600)
    plt.show()

def plotRegionalMapWithTectonics(mapBounds, stationGp, tectonicDatDir, imagesFileDir,  locAnotPair=None, transectLine=None, EquakeLoc=None, EqMag=None, EqDep=None, form='ps'):
    # function plots map with focus on tectonics, earthquake distribution, and tectonic groupings... 
        
    
    lonLeft = mapBounds[0]; lonRight = mapBounds[2]; lonStep=20
    latLeft = mapBounds[1]; latRight = mapBounds[3]; latStep=10
    fzHawaiiDir = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/PlateBoundaryShapeInfo/FractureZones/SHP/GSFML_SF_FZ_RM'
    fzHawaiiFile = 'GSFML_SF_FZ_RM'
    
    #trench
    plateBoundaries = "/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/PlateBoundaryShapeInfo/PLATES_PlateBoundary_ArcGIS/trench"
    plateFile = 'trench'
    
    # Plot topography and sea-floor ...
    
    # read in etopo5 topography/bathymetry.
    #url = 'http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/etopo5.nc'
    url = '/Users/tmo22/Documents/OlugbojiKarato2011/AgeGridsMuller/Data/etopo5.nc'
    etopodata = Dataset(url)
    
    topoin = etopodata.variables['ROSE'][:]
    lons = etopodata.variables['ETOPO05_X'][:]
    lats = etopodata.variables['ETOPO05_Y'][:]

    # read in sea floor age
    url = '/Users/tmo22/Documents/OlugbojiKarato2011/AgeGridsMuller/Data/age.3.6.nc'
    agedata = Dataset(url)
    
    agein = agedata.variables['z'][:] / 100
    agelons = agedata.variables['x'][:]
    agelats = agedata.variables['y'][:]


    # shift data so lons go from -180 to 180 instead of 20 to 380.
    topoin,lons = shiftgrid(180.,topoin,lons,start=False)
    agein, agelons = shiftgrid(180.,agein,agelons,start=False)

    plt.figure(figsize=(11,11))
    # lon_0, lat_0 are the center point of the projection.
    # resolution = 'l' means use low resolution coastlines.
    # create polar stereographic Basemap instance.
    m = Basemap(projection='lcc',\
            llcrnrlat=latLeft,urcrnrlat=latRight,\
            llcrnrlon=lonLeft,urcrnrlon=lonRight, lat_0=mapBounds[5], lon_0=mapBounds[4], lat_1=latLeft)
    
    # transform to nx x ny regularly spaced 5km(2km) native projection grid
    nx = int((m.xmax-m.xmin)/5000.)+1; ny = int((m.ymax-m.ymin)/5000.)+1

    topodat = m.transform_scalar(topoin,lons,lats,nx,ny, masked=True)
    agedat, agex, agey =  m.transform_scalar(agein,agelons,agelats,nx,ny, masked=True, returnxy= True)
        
    m.drawcoastlines()
    m.drawparallels(np.arange(latLeft, latRight, latStep), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(lonLeft, lonRight, lonStep), labels=[0, 0, 0, 1], linewidth=0)
    
    
    #1.  Topography shade and contour...
    clevSwell = np.arange(-4100, -4000, 5)
    clevLand = np.arange(-100, 100,5)
    clevs = np.append(clevSwell, clevLand)
    
    cmap = plt.cm.gray_r
    discLevs =  [-100, 0, 5, 100, 200, 500, 1000, 2000, 3000, 40000]
    im = m.imshow(topodat, cmap=cm.get_cmap('gray_r'), zorder=0, norm = colors.BoundaryNorm(discLevs,cmap.N,clip=False))

    #2. Contour
    x, y = m(*np.meshgrid(lons,lats))
    m.contour(x, y, topoin,clevs, interpolation='guassian', linewidths=0.1,colors='black', zorder=2)
    
    discLevs = [0, 9.7, 20.1, 33.1, 40.1, 47.9, 55.9, 67.7, 83.5, 120.4, 126.7, 131.9, 139.6, 147.7,
            154.3, 180]

    

    cmap = plt.cm.jet_r

    #3. Age colormap
    #im = m.imshow(agedat,cmap=cm.get_cmap('jet_r'), alpha = 0.4, zorder=1, 
    #              norm = colors.BoundaryNorm(discLevs,cmap.N,clip=False))
    
    #im = m.contour(agex, agey, agedat,cmap=cm.get_cmap('jet_r'), zorder=1, 
    #              norm = colors.BoundaryNorm(discLevs,cmap.N,clip=False))
    
                  
    #4. Plot Fracture Zones... & Ridges ...
    shp_info_Hawaii = m.readshapefile(fzHawaiiDir,fzHawaiiFile,drawbounds=False)
    for shapedict, shape in zip(m.GSFML_SF_FZ_RM_info, m.GSFML_SF_FZ_RM):
        xx,yy = zip(*shape)
        m.plot(xx,yy,linewidth=0.5,color='k', zorder=1)
        
    shp_info = m.readshapefile(plateBoundaries,plateFile, drawbounds=False)
    for shapedict, shape in zip(m.trench_info, m.trench):
        desc = shapedict['geogdesc']
        #print desc
        #if desc.find('PACIFIC') != -1 or desc.find('EXPLORER') != -1 or desc.find('NAZCA') != -1 :
        xx,yy = zip(*shape)
        m.plot(xx,yy,linewidth=2,color='black')
        
    # 4. Continued.. use local data downloaded from http://www.mri-jma.go.jp/Dep/st/member/fhirose/en/en.PAC.html
    #    load contours and attempt to plot
    x,y,z = np.loadtxt(tectonicDatDir+'PAC/plate_combine.dat', unpack=True)

    #size of 1 m grid
    nx = (int(lonRight - lonLeft + 1))#CHANGE HERE
    ny = (int(latRight - latLeft + 1))#CHANGE HERE

    # Generate a regular grid to interpolate the data.
    xi = np.linspace(lonLeft, lonRight, nx)
    yi = np.linspace(latLeft, latRight, ny)
    
    xi, yi = m(*np.meshgrid(xi, yi))
    xx, yy = m(x,y)
    
    zi = ml.griddata(xx,yy,z,xi,yi,interp='nn')
    #print xi, yi, zi
    
    levels = [10,50,80,100,200,  250, 300,400, 450, 500, 550]
    #manLocs = [m(130,38), m(132,38), m(135,38), m(140,38), m(141,38), m(142,38) ]
    cs = m.contour(xi,yi,zi,levels, linewidths=2, colors='brown', linestyles='dashed', zorder=4)
    plt.clabel(cs, levels[1::3], inline=1, fmt='%1d', fontsize=15, fontweight='bold', colors='blue' )#, manual= manLocs)
    
    
    # Place hotspot locations on map
    #hotspotFileDir = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/CrustalAgeData/'
    #hotspotLocs = returnHotspotLocations(hotspotFileDir+'hotspots.gmt')
    #hotspotLonLat = np.array( hotspotLocs )
    #xhot, yhot = m(hotspotLonLat[:,0], hotspotLonLat[:,1])
    #m.scatter(xhot, yhot, marker='o', s=500, c= 'red', alpha=1, zorder= 4)
    
    #4. Plot station locations - single group, then test all groups.
    markType = ['*', 'o', '^']
    markColor = ['black', 'black', 'black']
    markEdge = ['white', 'white', 'white']
    
    for indxGp in range(len(stationGp)):
        for indxSt in range(len(stationGp[indxGp])):
            stLonLat = np.array( stationGp[indxGp][indxSt] )
            xsta, ysta = m(stLonLat[:,1], stLonLat[:,0])
            
            #if indxGp == 2:
            #print "stations: ", stLonLat
            yrColIndx = stLonLat[:,2] -2
            #print "No. of yrs operational:", yrColIndx
            #else:
            #   yrColIndx = indxGp
            for ithSta in range(len(stationGp[indxGp][indxSt])):
                m.scatter(xsta[ithSta], ysta[ithSta], marker=markType[indxGp], s=250, c= markColor[int(yrColIndx[ithSta])], alpha=1, 
                edgecolor = markEdge[int(yrColIndx[ithSta])], zorder= 4)
    
    #5.1 Annotate plot with station name using python text command
    if locAnotPair is None:
        print "Annotation not set"
    else:
        print "Will Print Annotation on Map"
        for indxGp in range(len(locAnotPair)):
            for indxSt in range(len(locAnotPair[indxGp])):
                stLonLat = np.array( locAnotPair[indxGp][indxSt] )
                #print stLonLat
                xsta, ysta = m( float(stLonLat[1]), float(stLonLat[0])+0.5)
                name = stLonLat[2]
                plt.annotate(name, xy=(xsta, ysta), fontsize=14, fontweight='bold', zorder= 4, alpha=1)
            
    #5.2 Plot transects on station, visualize lateral coverage ..
    if transectLine is None:
        print "No transects defined"
    else:
        print "plotting transect as great circle"
        for indxTrans in range(len(transectLine)):
            #print transectLine[indxTrans]
            lonA = transectLine[indxTrans][0]; latA = transectLine[indxTrans][1] 
            lonB = transectLine[indxTrans][2]; latB =transectLine[indxTrans][3] 
        
            m.drawgreatcircle(lonA, latA, lonB, latB, color='white')
            
    #5.3. Plot Earthquakelocations - single group, then test all groups.
    if EquakeLoc is None:
        print "No Equake Data Set"
    else:
        print "plot equake data "
        
        for indxGp in range(len(EquakeLoc)):
                stLonLat = np.array( EquakeLoc[indxGp] )
                xsta, ysta = m(stLonLat[:,1], stLonLat[:,0])
                eqInfo = np.array(EqMag[indxGp])
                eqInfo = np.array(EqDep[indxGp])
                #print "Depth ", eqInfo
                c = m.scatter(xsta, ysta, marker='o', s=100, c= eqInfo, 
                edgecolor = 'black', cmap=plt.cm.hot_r, alpha=0.5, zorder= 3)
        
        
        cbar = plt.colorbar(c, pad=0.10, orientation= 'horizontal', alpha=1)
        cbar.ax.set_xlabel('Earthquake Depth')
    
    #m.drawmapboundary(fill_color='aqua')    
    
    # save figures ..
    plt.savefig(imagesFileDir+'.'+form, format=form)
    plt.show()
    
def plotGlobalMap(coords, stationGp, imagesFileDir, EquakeLoc=None):
    lonIn = coords[0]; latIn = coords[1]
    
    # Plot topography and sea-floor ...
    
    # read in etopo5 topography/bathymetry.
    url = '/Users/tmo22/Documents/OlugbojiKarato2011/AgeGridsMuller/Data/etopo5.nc'
    etopodata = Dataset(url)
    
    topoin = etopodata.variables['ROSE'][:]
    lons = etopodata.variables['ETOPO05_X'][:]
    lats = etopodata.variables['ETOPO05_Y'][:]

    # read in sea floor age
    url = '/Users/tmo22/Documents/OlugbojiKarato2011/AgeGridsMuller/Data/age.3.6.nc'
    agedata = Dataset(url)
    
    agein = agedata.variables['z'][:] / 100
    agelons = agedata.variables['x'][:]
    agelats = agedata.variables['y'][:]


    # shift data so lons go from -180 to 180 instead of 20 to 380.
    topoin,lons = shiftgrid(180.,topoin,lons,start=False)
    agein, agelons = shiftgrid(180.,agein,agelons,start=False)

    plt.figure(figsize=(11,11))
    # lon_0, lat_0 are the center point of the projection.
    # resolution = 'l' means use low resolution coastlines.
    # create polar stereographic Basemap instance.
    m = Basemap(projection='ortho',lon_0=lonIn,lat_0=latIn,resolution='l')
        
    # transform to nx x ny regularly spaced 5km(2km) native projection grid
    nx = int((m.xmax-m.xmin)/5000.)+1; ny = int((m.ymax-m.ymin)/5000.)+1

    topodat = m.transform_scalar(topoin,lons,lats,nx,ny, masked=True)
    agedat, agex, agey =  m.transform_scalar(agein,agelons,agelats,nx,ny, masked=True, returnxy= True)
        
    m.drawcoastlines()
    
    #1.  Topography shade and contour...
    clevSwell = np.arange(-4100, -4000, 5)
    clevLand = np.arange(-100, 100,5)
    clevs = np.append(clevSwell, clevLand)
    
    im = m.imshow(topodat, cmap=cm.get_cmap('gray'), zorder=0)

    #2. Contour
    x, y = m(*np.meshgrid(lons,lats))
    m.contour(x, y, topoin,clevs, interpolation='guassian', linewidths=0.05,colors='white', zorder=2)
    
    discLevs = [0, 9.7, 20.1, 33.1, 40.1, 47.9, 55.9, 67.7, 83.5, 120.4, 126.7, 131.9, 139.6, 147.7,
            154.3, 180]

    #discLevs =  [120.4, 126.7, 131.9, 139.6, 147.7,154.3, 180]

    cmap = plt.cm.jet_r

    #3. Age colormap
    im = m.imshow(agedat,cmap=cm.get_cmap('jet_r'), alpha = 0.4, zorder=1, 
                  norm = colors.BoundaryNorm(discLevs,cmap.N,clip=False))
                  
    #4. Plot station locations - single group, then test all groups.
    #markType = ['*', 'o', '^']
    #markColor = ['red', 'black', 'white']
    #markEdge = ['black', 'white', 'black']
    
    #for indxGp in range(len(stationGp)):
    #    for indxSt in range(len(stationGp[indxGp])):
    #        stLonLat = np.array( stationGp[indxGp][indxSt] )
    #        xsta, ysta = m(stLonLat[:,1], stLonLat[:,0])
    #        m.scatter(xsta, ysta, marker=markType[indxGp], s=70, c= markColor[indxSt], 
    #        edgecolor = markEdge[indxSt], alpha=1, zorder= 3)
    
    #4. Plot station locations - single group, then test all groups.
    markType = ['*', 'o', '^']
    markColor = ['red', 'white', 'black']
    markEdge = ['black', 'black', 'white']
    
    for indxGp in range(len(stationGp)):
        for indxSt in range(len(stationGp[indxGp])):
            stLonLat = np.array( stationGp[indxGp][indxSt] )
            xsta, ysta = m(stLonLat[:,1], stLonLat[:,0])
            
            #if indxGp == 2:
            #print "stations: ", stLonLat
            yrColIndx = stLonLat[:,2] -2
            #print "No. of yrs operational:", yrColIndx
            #else:
            #   yrColIndx = indxGp
            for ithSta in range(len(stationGp[indxGp][indxSt])):
                m.scatter(xsta[ithSta], ysta[ithSta], marker=markType[indxGp], s=70, c= markColor[int(yrColIndx[ithSta])], alpha=1, 
                edgecolor = markEdge[int(yrColIndx[ithSta])], zorder= 3)        
        
    cb = m.colorbar(im,"right", size="5%", pad='2%')
    cb.ax.set_ylabel('Million Years B.P.')
    
    #4. Plot Earthquakelocations - single group, then test all groups.
    if EquakeLoc is None:
        print "No Equake Data Set"
    else:
        print "plot equake data "
        
        for indxGp in range(len(EquakeLoc)):
                stLonLat = np.array( EquakeLoc[indxGp] )
                xsta, ysta = m(stLonLat[:,1], stLonLat[:,0])
                m.scatter(xsta, ysta, marker='o', s=10, c= 'black', 
                edgecolor = 'black', alpha=1, zorder= 3)
          
        
                  
    #4. Plot Fracture Zones... & Ridges ...
    fzHawaiiDir = '/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/PlateBoundaryShapeInfo/FractureZones/SHP/GSFML_SF_FZ_RM'
    fzHawaiiFile = 'GSFML_SF_FZ_RM'
    
    #ridges
    plateBoundaries = "/Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012/SeismicData/SACbyStation/AllStationProcessed/PlateBoundaryShapeInfo/PLATES_PlateBoundary_ArcGIS/ridge"
    plateFile = 'ridge'

    #shp_info_Hawaii = m.readshapefile(fzHawaiiDir,fzHawaiiFile,drawbounds=False)
    #for shapedict, shape in zip(m.GSFML_SF_FZ_RM_info, m.GSFML_SF_FZ_RM):
    #    xx,yy = zip(*shape)
    #    m.plot(xx,yy,linewidth=2,color='k', zorder=1)
        
    #shp_info = m.readshapefile(plateBoundaries,plateFile, drawbounds=False)
    #for shapedict, shape in zip(m.ridge_info, m.ridge):
    #    desc = shapedict['geogdesc']
        #print desc[0:7]
    #    if desc.find('PACIFIC') != -1 or desc.find('EXPLORER') != -1 or desc.find('NAZCA') != -1 :
    #        xx,yy = zip(*shape)
    #        m.plot(xx,yy,linewidth=1.5,color='white')
    
    plt.savefig(imagesFileDir, format='png', dpi=600)
    plt.show()


###### Event Statistics Plot ...
#1.0 Polar Event -Single Station - Plot ...
def showPolarEventPlot(evntBaz, evntDist, areaMag, evMag, pltTitle, fName, format='ps'):
    #areaMag.append(3**eqHeaders['mag'])
    ax = plt.subplot(1,1,1, projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    c = plt.scatter((np.pi/180)*np.array(evntBaz), np.array(evntDist),  
                            s=70, c= np.array(evMag), cmap=plt.cm.hot)
    c.set_alpha(0.5)
    

    cbar = plt.colorbar(c, pad=0.10)
    
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    
    cbar.ax.set_ylabel('Earthquake Magnitude')
    
    plt.title(pltTitle, y= 1.08)
    plt.savefig(fName+'.'+format, format=format)
    plt.show()

#2.0 Histogram Event - Single Station - Plot
def showHistogramScatterPlot(evntBaz, evntDist, area, evMag,  pltTitle, fName, format='ps', nbins=4):
    # the random data
    x = np.array(evntBaz)
    y = np.array(evntDist)
    
    
    fig, axScatter = plt.subplots(figsize=(8.5,8.5))
    
    # the scatter plot:
    mapHolder = axScatter.scatter(x, y, s = 70, c= np.array(evMag), cmap=plt.cm.hot, alpha=0.5)
    #axScatter.set_aspect(1.)
    plt.xlabel("Back Azimuth")
    plt.ylabel("Great Circle Distance")
    
    
    # create new axes on the right and on the top of the current axes
    # The first argument of the new_vertical(new_horizontal) method is
    # the height (width) of the axes to be created in inches.
    divider = make_axes_locatable(axScatter)
    axHistx = divider.append_axes("top", 1.2, pad=0.18, sharex=axScatter)
    plt.title(pltTitle)
    axHisty = divider.append_axes("right", 1.2, pad=0.25, sharey=axScatter)
    
    # make some labels invisible
    plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
             visible=False)
    
    # now determine nice limits by hand:
    binwidth = 10
    xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
    lim = ( int(xymax/binwidth) + 1) * binwidth
    
    bins = np.arange(-lim, lim + binwidth, binwidth)
    valX = axHistx.hist(x, bins=bins)
    valY = axHisty.hist(y, bins=bins, orientation='horizontal')
    
    # the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
    # thus there is no need to manually adjust the xlim and ylim of these
    # axis.
    
    #axHistx.axis["bottom"].major_ticklabels.set_visible(False)
    for tl in axHistx.get_xticklabels():
        tl.set_visible(False)
    
    histX_Yticks = np.linspace(min(valX[0]),max(valX[0]),nbins)
    histX_Yticks = histX_Yticks.astype(int)
    axHistx.set_yticks(histX_Yticks)
    
    #axHisty.axis["left"].major_ticklabels.set_visible(False)
    for tl in axHisty.get_yticklabels():
        tl.set_visible(False)
    
    histY_Yticks = np.linspace(min(valY[0]),max(valY[0]),nbins-1)
    histY_Yticks = histY_Yticks.astype(int)
    axHisty.set_xticks(histY_Yticks)
    
    #plt.locator_params(nbins=4)
    
    plt.savefig(fName+'.'+format, format=format)

    plt.draw()
    plt.show()

#3.0 Bar chart for station event summary - This is for all Station Bar chart...
def showStationEventSummary(staNames, totGdEvnts, totAllEvnts, noYrsOper, statType, fName, format):
    # prints station event summmary as a bar chart color coded by noYrsOperation or OBS vs Borehole
    #
    # staNames : Names of stations used for annotation - maybe add age?
    # totGdEvnts: total good events used for first dimension bar chart and annotation
    # totBdEvnts: total bad events for second dimension bar chart and annotation
    # noYrsOper: could be used for color code 
    # statType: could be color coded or line coded
    print "not yet coded... work in progress"
    
    plt.figure(figsize=(12,7))

    X = np.arange(len(staNames))
    val = np.max(totAllEvnts)
    valA = 4
    totAllEvnts = (-1./valA) * np.array(totAllEvnts)

    val = np.max(totGdEvnts)
    valG = 1
    totGdEvnts =  (1./valG) * np.array(totGdEvnts)
    #print val, totAllEvnts
    
    barGd = plt.bar(X, totGdEvnts, facecolor='#9999ff', edgecolor='white')
    barAll = plt.bar(X, totAllEvnts, facecolor='#ff9999', edgecolor='white')
    
    for x in X:
        
        if (statType[x] == 1):
            rectO = barGd.patches[x]
            rectO.set_edgecolor('black')
            
            rectO = barAll.patches[x]
            rectO.set_edgecolor('black')
            #rectO.set_linestyle('dashed')
    
    for xPos, yPosn, yPos, staName in zip(X, totAllEvnts, totGdEvnts, staNames):
        plt.text(xPos+0.5, yPos+1, '%d' %yPos, fontsize=11, ha='center')
        plt.text(xPos, -20,  staName, fontsize=12, fontweight='bold')
        yPosnN = -1 * valA * yPosn
        plt.text(xPos+0.5, yPosn-20, '%d' %yPosnN, fontsize=11, ha='center')
    
    plt.xlim([-1,len(staNames)+1])
    plt.ylim([-450, 300])
    plt.yticks([])
    plt.xticks([])
    #plt.yticks([300, 200, 100, 0, -400, -800, 1200, -1600])
    
    plt.savefig(fName+'.'+format, format=format)
    plt.show()

#3.0 Summary Histogram Plots for All Stations ... All 19 stations. See Manuscript
def showEqDistAzHist(staNames, evBAZbySta, evDistbySta, fName, format):
    
    fig = plt.figure(figsize=(12,12))
    
    # Iterate through stations ..
    groupBounds = [[0,5], [5,12], [12,19]]
    cntFig = 1
    for groupBnd in groupBounds:
        #print groupBnd
        rt = groupBnd[0]
        lt = groupBnd[1]
        
        for staNme, evBaz, evDist in zip(staNames[rt:lt], evBAZbySta[rt:lt], evDistbySta[rt:lt] ):
        
            fig.add_subplot(3, 2,cntFig)
            
            #print "All", staNme, len(evBaz), len(evDist)
            y, binEdges = np.histogram(evBaz)
            bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
            plt.plot(bincenters, y, '-',linewidth=2,label=staNme)
            plt.scatter(bincenters, y)
            plt.xlim([1, 360])
            plt.ylim([0, 51])
            
            nextFig = cntFig + 1
            fig.add_subplot(3, 2 ,nextFig)
            y, binEdges = np.histogram(evDist)
            bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
            plt.plot(bincenters, y, '-', linewidth=2,label=staNme)
            plt.scatter(bincenters, y)
            plt.xlim([1, 165])
            plt.ylim([0, 73])
            
            #print nextFig
            
        cntFig = cntFig + 2
        #print cntFig
        
        #plt.legend()
        plt.legend(bbox_to_anchor=(0.93, 1), loc=2, borderaxespad=0.)
    
    plt.savefig(fName+'.'+format, format=format)
    plt.show()
    
    