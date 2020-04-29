import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from netCDF4 import Dataset
import pylab
import numpy.ma as ma
import pandas as pd
import math as math
import glob
import xarray
import os
import sys
from optparse import OptionParser
from wofs_sounding_variable_extractor_Lidar_Relative import *

#################### Command Options ####################
parser = OptionParser()
parser.add_option("-o", dest="indir", type="string", default= None, help="Input file of observation to be compared with WoFS")

(options, args) = parser.parse_args()

if (options.indir == None):
   parser.print_help()
   print()
   sys.exit(1)
else:
   indir = options.indir


############## read in observation files ################
files = os.listdir(indir)
files.sort()
print('Here are all the files in the folder you designated: ' + str(files))

for f, infile in enumerate(files):
	infile = os.path.join(indir, infile)
	Nfile = infile
	Lidar = Dataset(infile)
	date = []
	date.append(infile[43:51])
	print('Obs date: ' + str(date[0]))
	
	#Converts Lidar time output to HHMM format
	tLid = Lidar.variables['hour']
	tUTC = np.ones(len(tLid), dtype = np.int)
	minUTC = np.ones(len(tLid), dtype = np.int)
	for i in range(len(tLid)):
		minUTC[i] = int((tLid[i]-int(tLid[i]))*60)
		if minUTC[i] < 10:
			tUTC[i] = str(int(tLid[i])) + '0' + str(minUTC[i])
		else:
			tUTC[i] = str(int(tLid[i])) + str(minUTC[i])
# 	time = [2235,2240]
#  	time = [2245,2250,2255,2300,2305,2310]
 	time = [2315,2320]
#  	,2315,2320]
	print(time)
	print("Obs Start time: " + str(tUTC[0]))
	print("Times for which Lidar data is available: " + str(tUTC))

	
	#Need example grid in order to determine closest surface grid point to observation
	for x in range(len(time)):
		print(time[x])
		
		WoFSfile = sounding_extractor(str(time[x]),date[0],Nfile)

		print('WoFS Sounding file used in comparison: ' + WoFSfile)
	
		#files date/time
		WoFSrun = WoFSfile[21:25]
		WoFSvalid = WoFSfile[26:30]
		print( WoFSrun + 'UTC' + ' forecast, valid at ' + WoFSvalid +  'UTC')
		WoFS = Dataset('WoFS_Soundings/' + str(date[0][-4:]) + '/' + WoFSfile)

		#WoFS variables have dimensions [member, z-level, x, y]
		xdex, ydex = [1,1] #pulling middle grid point from the 9 provided in the sounding file
		Wlat = WoFS.variables['xlat']
		Wlon = WoFS.variables['xlon']
		uWoFS = WoFS.variables['u'][:,:,xdex,ydex]
		vWoFS = WoFS.variables['v'][:,:,xdex,ydex] 
		hWoFS = WoFS.variables['z_agl'][:,:,xdex,ydex]
		member_list = [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 2, 3, 4, 5, 6, 7, 8, 9]
		#ysu = [1,2,7,8,13,14]
		#myj = [3,4,9,10,15,16]            #pbl schemes based off of member
		#mynn = [5,6,11,12,17,18]
		ysu = [0,4,5,10,15,16]
		myj = [1,6,7,11,12,17]             #pbl indicies
		mynn = [2,3,8,9,13,14]
	
		#Lidar variables
		#Lidar wind variable have dimensions [t,z]
		uLid = Lidar.variables['u']
		vLid = Lidar.variables['v']
		hLid = Lidar.variables['height']
	
			
		# #Mean WoFS U,V,h
		uMeanWoFS = np.ones(len(uWoFS[1,:]))*np.nan
		for i in range(len(uWoFS[1,:])):
			uMeanWoFS[i] = np.average(uWoFS[:,i])
		vMeanWoFS = np.ones(len(vWoFS[1,:]))*np.nan
		for i in range(len(vWoFS[1,:])):
			vMeanWoFS[i] = np.average(vWoFS[:,i])
		hMeanWoFS = np.ones(len(hWoFS[1,:]))*np.nan
		for i in range(len(hWoFS[1,:])):
			hMeanWoFS[i] = np.average(hWoFS[:,i])
		
		#Averaging Lidar data based off of time
		Index = np.where((tUTC <= (int(time[x]) + 3)) & (tUTC >= (int(time[x]) - 3)))[0]
		print("Scan times used in averaging: " + str(tUTC[Index]))
		print(Index)

		#max height of obs, loop is for grabbing the lowest height of non-nan data for averaging
		hindex = []
		for j in Index:
			
			indh = np.nanargmax(uLid[j])
			print(indh)
			hindex.append(indh)
			maxhLid = hLid[indh]
			print("Maximum height of non-nan lidar data for time " + str(time[x]) + ' :' + str(maxhLid))
		maxhLid = np.min(hLid[hindex])
		indh = np.min(hindex)
		print("Maximum height of non-nan lidar data: " + str(maxhLid))
		print("Max height index for lidar: " + str(indh))

		#Index for WoFS data cutoff
		hWoFSco = []
		for j in range(50):
			if maxhLid > hWoFS[1,j]:
					hWoFSco.append(hWoFS[1,j])
		hdex = hWoFSco.index(max(hWoFSco)) + 1
		print("Max height index for WoFS: " + str(hdex))
	

		#Removing WoFS data above maximum observation height
		hMeanWoFSco = np.ones(len(hMeanWoFS))*np.nan
		for j in range(len(hMeanWoFS)):
			if maxhLid > hMeanWoFS[j]:
				hMeanWoFSco[j] = hMeanWoFS[j]
		WoFShdex = np.argmax(hMeanWoFSco, axis = None) + 1
		hMeanWoFSco = hMeanWoFSco[~np.isnan(hMeanWoFSco)]
		hMeanWoFS = hMeanWoFS[:hdex]
		print("Max Height for WoFS data: " +str(hWoFS[:,(WoFShdex-1)]))


		#Finding levels in lidar closest to WoFS grid
		distlvl = np.ones((len(hMeanWoFS),len(hLid[:])))*np.nan
		for i in range(len(hMeanWoFS)):
			for j in range(len(hLid[:])):
				distlvl[i,j] = abs((hMeanWoFS[i]-hLid[j]))

		decihgt = np.ones(len(distlvl), dtype = np.int)*np.nan
		decidex = np.ones(len(distlvl), dtype = np.int)*np.nan
		for i in range(len(distlvl)):
			decidex[i] = np.argmin(distlvl[i,:], axis = None) - 2
			decihgt[i] = hLid[decidex[i]]
	
		
		#Time Averaging Lidar Data, 6 min windows
		
		uLidAvg0 = np.ones(len(uLid[1,:]))*np.nan
		for i in range(len(uLid[1,:])):
			uLidAvg0[i] = np.average(uLid[Index,i])
		
		
		vLidAvg0 = np.ones(len(vLid[1,:]))*np.nan
		for i in range(len(vLid[1,:])):
			vLidAvg0[i] = np.average(vLid[Index,i])
	

	
		#Decimating wind onto decihgt (averaging nearby winds)
		indbin = 2

		deciU0 = np.ones(len(decidex))*np.nan
		deciU0[0] = np.average(uLid[:3])                                 #averaging lowest levels to get surface wind
		for i in range(len(decidex)-1):                                                #decimating using bins of 6 indicies(about 30 meters)
			deciU0[i+1] = np.average(uLidAvg0[int((decidex[i+1]-indbin)):int((decidex[i+1]+indbin))])

	
		deciV0 = np.ones(len(decidex))*np.nan
		deciV0[0] = uLidAvg0[0] 
		deciV0[1] = uLidAvg0[1]                               #averaging lowest levels to get surface wind
		for i in range(len(decidex)-1):                                                #decimating using bins of 6 indicies(about 30 meters)
			deciV0[i+1] = np.average(vLidAvg0[int((decidex[i+1]-indbin)):int((decidex[i+1]+indbin))])
		

	###################   STATISTICS   ########################
		#calculating difference between observation and each member, so diff should have dimensions of vertical level vs member
		#Differences
		#Zonal Component
		diffU0 = np.ones((18,len(deciU0)))*np.nan
		for i in range(18):
			for j in range(len(deciU0)):
				diffU0[i,j] = uWoFS[i,j]-deciU0[j]


		#Meridional Component
		diffV0 = np.ones((18,len(deciV0)))*np.nan
		for i in range(18):
			for j in range(len(deciV0)):
				diffV0[i,j] = vWoFS[i,j]-deciV0[j]
			


	############# Writing .nc file #########################3

		ne = 18
		nz = len(deciU0)
		outdir = "/home/jordan.laser/WoFS_stats/0517/Relative/Lidar"
		outname = "Lidar_WoFSMemberStats_" + date[0]+ "_2230_" + str(time[x]) + ".nc"         #output file
		output_path = os.path.join(outdir,outname)
		try:
			fout = netCDF4.Dataset(output_path, "w")
		except:
		   print("Could not create %s!\n" % output_path)


		### Create file and dimensions: ###
		fout.createDimension('NE', ne)
		fout.createDimension('NZ', nz)

		### 3-D Variables ###
		diffu_var = fout.createVariable('DU0', 'f4', ('NE','NZ',))
		diffu_var.long_name = "difference in U-component of wind"
		diffu_var.units = "m s**-1"
	

		diffv_var = fout.createVariable('DV0', 'f4', ('NE','NZ',))
		diffv_var.long_name = "difference in V-component of wind"
		diffv_var.units = "m s**-1"
	
	
		diffq_var = fout.createVariable('h', 'f4', ('NE','NZ',))
		diffq_var.long_name = "Height of Veritcal Levels in WoFS"
		diffq_var.units = "m"


		###Write Variables###
		fout.variables['DU0'][:] = diffU0

	
		fout.variables['DV0'][:] = diffV0

		print(WoFShdex)
		fout.variables['h'][:] = hWoFS[:,:(WoFShdex-1)]

	


		fout.close()
		del fout