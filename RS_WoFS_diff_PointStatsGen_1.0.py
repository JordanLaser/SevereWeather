import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from netCDF4 import Dataset
from netCDF4 import MFDataset
import pylab
import numpy.ma as ma
import pandas as pd
import math as math
import glob
import xarray
import os
import sys
from optparse import OptionParser
from wofs_sounding_variable_extractor_RS_Point import *

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
	df = pd.read_csv(infile, skiprows = 3, sep = ',', header = None)
	date = []
	time = []
	date.append(infile[-19:-11])
	time.append(infile[-10:-6])
	
	print('Obs date: ' + str(date[0]))
	print('Obs time: ' + str(time[0]) + ' UTC')
	
	#Need example grid in order to determine closest surface grid point to observation
	Nfile = infile
	WoFSfile = sounding_extractor(time[0],date[0],Nfile)
	print('WoFS Sounding file used in comparison: ' + WoFSfile)
	
	if date[0] == '20190518':
		date[0] = '20190517'

	#files date/times
	WoFSrun = WoFSfile[21:25]
	WoFSvalid = WoFSfile[26:30]
	print( WoFSrun + ' UTC' + ' forecast, valid at ' + WoFSvalid +  ' UTC')
	WoFS = Dataset('WoFS_Soundings/' + str(date[0][-4:]) + '/' + WoFSfile)

	#Pull wind components
	#WoFS variables have dimensions [member, z-level, x, y]
	Sxdex, Sydex = [1,1] #pulling middle grid point from the 9 provided in the sounding file
	Wlat = WoFS.variables['xlat']
	Wlon = WoFS.variables['xlon']
	uWoFS = WoFS.variables['u'][:,:,Sxdex,Sydex]
	vWoFS = WoFS.variables['v'][:,:,Sxdex,Sydex]
	tWoFS = WoFS.variables['t'][:,:,Sxdex,Sydex] 
	pWoFS = WoFS.variables['p'][:,:,Sxdex,Sydex] 
	hWoFS = WoFS.variables['z_agl'][:,:,Sxdex,Sydex]
	qWoFS = WoFS.variables['q'][:,:,Sxdex,Sydex]
	
	thetaWoFS = (tWoFS+273.15)*(1000/pWoFS)**(287/1004)
	thetaWoFS = thetaWoFS - 273.15
	
	
	member_list = [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 2, 3, 4, 5, 6, 7, 8, 9]
	#ysu = [1,2,7,8,13,14]
	#myj = [3,4,9,10,15,16]            #pbl schemes based off of member
	#mynn = [5,6,11,12,17,18]
	ysu = [0,4,5,10,15,16]
	myj = [1,6,7,11,12,17]             #pbl indicies
	mynn = [2,3,8,9,13,14]

	#radiosonde variables
	#set columns of df to variable names
	temp = df[2]
	rh = df[3]
	dew = df[4]
	pressure = df[5]
	winddir = df[6]
	windspd = df[7]
	elv = df[8]
	alt = elv-elv[0]
	Slat = df[9]
	Slon = df[10]
	windU = -windspd*np.cos((winddir*3.1415/180)-(3.1415/2))
	windV = windspd*np.sin((winddir*3.1415/180)-(3.1415/2))
	w = df[13]
	
	#Calculate Potential Temperature (temp = kelvin)
	thetaObs = temp *(1000/pressure)**(287/1004)
	thetaObs = thetaObs - 273.15

	# #Switch units back to Celcius
	# temp = temp - 273.15

	#Calculate water vapor mixing ratio from obs
	e_0 = 6.1173
	t_0 = 273.16
	Rv = 461.5
	Lv = 2501000
	es = e_0*np.exp((Lv/Rv)*((1/t_0)-(1/temp)))
	e = (rh/100)*es
	q = 0.622*(e/(pressure-e))


	# #Mean WoFS U,V,h,t,p
	uMeanWoFS = np.ones(len(uWoFS[1,:]))*np.nan
	for i in range(len(uWoFS[1,:])):
		uMeanWoFS[i] = np.average(uWoFS[:,i])
	vMeanWoFS = np.ones(len(vWoFS[1,:]))*np.nan
	for i in range(len(vWoFS[1,:])):
		vMeanWoFS[i] = np.average(vWoFS[:,i])
	hMeanWoFS = np.ones(len(hWoFS[1,:]))*np.nan
	for i in range(len(hWoFS[1,:])):
		hMeanWoFS[i] = np.average(hWoFS[:,i])
	# tMeanWoFS = np.ones(len(tWoFS[1,:]))*np.nan

		
	########################################################
	#          QUALITY CONTROL
	########################################################
		
	#set cut off index to either the first below 2000 m or the maximum height of radiosonde if
	#sonde did not make it to 2000 m

	defaultH = 2000
	# hcutoffind = np.argmax(alt)
# 	for i in range(len(alt)):
# 		if alt[i] > defaultH:
# 			hcutoffind = i -1
# 			break
# 
# 	#removing values above 2000m or indicies with rh > 95%
# 	# cloud equals RH > 95
# 	cloudlvls = []
# 	for i in range(len(rh)):
# 
# 		if i < hcutoffind:
# 
# 			if rh[i] > 95:
# 				cloudlvls.append(alt[i])
# 
# 				if all(rh[i:hcutoffind] > 95):
# 					hcutoffind = i
# 					hcutoff = alt[i]
# 
# 			else:
# 				hcutoff = defaultH
# 
# 
# 	if len(cloudlvls) == 0:
# 		print("No Clouds")
# 	else:
# 		cloudbase = round(np.min(cloudlvls),2)
# 		cloudtop = round(np.max(cloudlvls),2)
# 		print("*WARNING* Clouds possible from: " + str(cloudbase) + 'm to '+ str(cloudtop) + 'm')

	hcutoff = defaultH

	print("Maximum height of QCed observation data: " + str(hcutoff) + 'm')
	
	hMeanWoFSco = np.ones(len(hMeanWoFS))*np.nan
	for j in range(len(hMeanWoFS)):
		if hcutoff > hMeanWoFS[j]:
			hMeanWoFSco[j] = hMeanWoFS[j]
		else:
			break
	WoFShdex = np.argmax(hMeanWoFSco, axis = None) + 1
	hMeanWoFSco = hMeanWoFSco[~np.isnan(hMeanWoFSco)]
	print('Truncation index for WoFS: ' + str(WoFShdex))
	print('Highest WoFS grid pt: ' + str(round(np.max(hMeanWoFSco),2)) + 'm')

	Sco = np.ones(len(alt))*np.nan
	for j in range(len(alt)):
		if hcutoff > alt[j]:
				Sco[j] = alt[j]
	Shdex = np.argmax(Sco, axis = None) + 1
	
	

	###################################################################
	#             DECIMATING OBSERVATION ONTO WOFS GRID
	###################################################################

	#Finding levels in radiosonde closest to WoFS grid
	distlvl = np.ones((len(hMeanWoFSco),len(alt)))*np.nan
	for i in range(len(hMeanWoFSco)):
		for j in range(len(alt)):
			distlvl[i,j] = abs((hMeanWoFSco[i]-alt[j]))

	decihgt = np.ones(len(distlvl))*np.nan
	decidex = np.ones(len(distlvl))*np.nan
	for i in range(len(distlvl)):
		decidex[i] = np.argmin(distlvl[i,:], axis = None)
		decihgt[i] = alt[decidex[i]]
		

	#Decimating wind onto decihgt (averaging nearby winds)
	indbin = 3
	deciU = np.ones(len(decidex))*np.nan
	deciU[0] = np.average(windU[:int(decidex[0])])                                 #averaging lowest levels to get surface wind
	for i in range(len(decidex)-1):                                                #decimating using bins of 6 indicies(about 30 meters)
		deciU[i+1] = np.average(windU[int((decidex[i+1]-indbin)):int((decidex[i+1]+indbin))])

	deciV = np.ones(len(decidex))*np.nan
	deciV[0] = np.average(windV[:int(decidex[0])])                                 #averaging lowest levels to get surface wind
	for i in range(len(decidex)-1):                                                #decimating using bins of 6 indicies(about 30 meters)
		deciV[i+1] = np.average(windV[int((decidex[i+1]-indbin)):int((decidex[i+1]+indbin))])

	decit = np.ones(len(decidex))*np.nan
	decit[0] = np.average(thetaObs[:int(decidex[0])])                                 #averaging lowest levels to get surface wind
	for i in range(len(decidex)-1):                                                #decimating using bins of 6 indicies(about 30 meters)
		decit[i+1] = np.average(thetaObs[int((decidex[i+1]-indbin)):int((decidex[i+1]+indbin))])

	deciq = np.ones(len(decidex))*np.nan
	deciq[0] = np.average(q[:int(decidex[0])])                                 #averaging lowest levels to get surface wind
	for i in range(len(decidex)-1):                                                #decimating using bins of 6 indicies(about 30 meters)
		deciq[i+1] = np.average(q[int((decidex[i+1]-indbin)):int((decidex[i+1]+indbin))])


	###################   STATISTICS   ########################
	#calculating difference between observation and each member, so diff should have dimensions of vertical level vs member
	#Differences
	#Zonal Component
	diffU = np.ones((18,len(deciU)))*np.nan
	for i in range(18):
		for j in range(len(deciU)):
			diffU[i,j] = uWoFS[i,j]-deciU[j]

	#Meridional Component
	diffV = np.ones((18,len(deciV)))*np.nan
	for i in range(18):
		for j in range(len(deciV)):
			diffV[i,j] = vWoFS[i,j]-deciV[j]


	#Thermodynamic differences
	#Linear Bias in t/q
	difftheta = np.ones((18,len(decit)))*np.nan
	for i in range(18):
		for j in range(len(decit)):
			difftheta[i,j] = (thetaWoFS[i,j]-decit[j])

	diffq = np.ones((18,len(deciq)))*np.nan
	for i in range(18):
		for j in range(len(deciq)):
			diffq[i,j] = qWoFS[i,j]-deciq[j]


	###########################################################################
	###########################################################################
	############################## Writing output file ########################
	ne = 18
	nz = len(deciU)
	outdir = "/home/jordan.laser/WoFS_stats/0517/Relative"
	outname = "WoFSMemberStats_" + str(date[0]) + "_" + str(time[0]) + ".nc"         #output file
	output_path = os.path.join(outdir,outname)
	try:
		fout = netCDF4.Dataset(output_path, "w")
	except:
	   print("Could not create %s!\n" % output_path)


	### Create file and dimensions: ###
	no = 1
	fout.createDimension('NE', ne)
	fout.createDimension('NZ', nz)
	fout.createDimension('NO', no)

	### 3-D Variables ###
	diffu_var = fout.createVariable('DU', 'f4', ('NE','NZ',))
	diffu_var.long_name = "difference in U-component of Wind"
	diffu_var.units = "m s**-1"

	diffv_var = fout.createVariable('DV', 'f4', ('NE','NZ',))
	diffv_var.long_name = "Difference in V-component of Wind"
	diffv_var.units = "m s**-1"

	difft_var = fout.createVariable('DT', 'f4', ('NE','NZ',))
	difft_var.long_name = "Difference in Potential Temperature"
	difft_var.units = "\u00B0C"

	diffq_var = fout.createVariable('Dq', 'f4', ('NE','NZ',))
	diffq_var.long_name = "Difference in Mixing Ratio"
	diffq_var.units = "Kg/Kg"

	diffq_var = fout.createVariable('h', 'f4', ('NE','NZ',))
	diffq_var.long_name = "Height of Veritcal Levels in WoFS"
	diffq_var.units = "m"

	deciU_var = fout.createVariable('deciU', 'f4', ('NO','NZ',))
	deciU_var.long_name = "Decimated Observed U-component of Wind"
	deciU_var.units = "m/s"

	deciV_var = fout.createVariable('deciV', 'f4', ('NO','NZ',))
	deciV_var.long_name = "Decimated Observed V-component of Wind"
	deciV_var.units = "m/s"

	decit_var = fout.createVariable('decit', 'f4', ('NO','NZ',))
	decit_var.long_name = "Decimated Observed Potential Temperature"
	decit_var.units = "K"

	deciq_var = fout.createVariable('deciq', 'f4', ('NO','NZ',))
	deciq_var.long_name = "Decimated Observed Mixing Ratio"
	deciq_var.units = "Kg/Kg"

	thetaWoFS_var = fout.createVariable('tWoFS', 'f4', ('NE','NZ',))
	thetaWoFS_var.long_name = "WoFS Potential Temperature"
	thetaWoFS_var.units = "K"


	###Write Variables###
	fout.variables['DU'][:] = diffU
	fout.variables['DV'][:] = diffV
	fout.variables['DT'][:] = difftheta
	fout.variables['Dq'][:] = diffq
	fout.variables['h'][:] = hWoFS[:,:(WoFShdex-1)]
	fout.variables['deciU'][:] = deciU
	fout.variables['deciV'][:] = deciV
	fout.variables['decit'][:] = decit
	fout.variables['deciq'][:] = deciq
	fout.variables['tWoFS'][:] = thetaWoFS[:,:(WoFShdex-1)]

	fout.close()
	del fout
	
	#Hodograph
	a = 0.8
	lw = 0.7
	vertical_levels = hMeanWoFS[:11]
	xticks = [-30,-20,-10,10,20,30]
	yticks = [-30,-20,-10,10,20,30]
	fig = pylab.figure(figsize = (9,4.5))  
	ax1 = fig.add_subplot(1, 1, 1)
	ax1.set_xlim(-30,30)
	ax1.set_ylim(0,30)
	circle1 = plt.Circle((0,0), radius = 10, linestyle = '--', fill = False)
	circle2 = plt.Circle((0,0), radius = 20, linestyle = '--', fill = False)
	circle3 = plt.Circle((0,0), radius = 30, linestyle = '--', fill = False)
	ax1.add_artist(circle1)
	ax1.add_artist(circle2)
	ax1.add_artist(circle3)
	ax1.axhline(y=0, color = 'k')
	ax1.axvline(x=0, color = 'k')
	plt.axis('off')
	ax1.text(10,0, '10')
	ax1.text(20,0, '20')
	ax1.text(30,0, '30')
	ax1.text(-10,0, '-10')
	ax1.text(-20,0, '-20')
	ax1.text(-30,0, '-30')
	ax1.text(0,9, '10')
	ax1.text(0,19, '20')
	ax1.text(0,29, '30')
	ax1.text(0,-10, '-10')
	ax1.text(0,-20, '-20')
	ax1.text(0,-30, '-30')
	for i in ysu:
		ax1.plot(uWoFS[i,:WoFShdex],vWoFS[i,:WoFShdex], linewidth = lw, alpha = a, color = 'r')
	for i in myj:
		ax1.plot(uWoFS[i,:WoFShdex],vWoFS[i,:WoFShdex], linewidth = lw, alpha = a, color = 'y')
	for i in mynn:
		ax1.plot(uWoFS[i,:WoFShdex],vWoFS[i,:WoFShdex], linewidth = lw, alpha = a, color = 'c')
	
		
	p1, = ax1.plot(uMeanWoFS[:WoFShdex],vMeanWoFS[:WoFShdex], color = 'k', linewidth = 2.0, label = 'WoFS Mean')

	p1, = ax1.plot(deciU[:WoFShdex],deciV[:WoFShdex], '--', color = 'k', linewidth = 2.0, label = "Radiosonde: " + str(time[0]))
	ax1.legend(loc = "lower left")

	plt.savefig('plots/Point/Wind_Plots/RS_WoFS_hodo_' + str(date[0]) + '_' + WoFSrun + '_' + str(time[0]) + '.png', dpi = 300)
	
