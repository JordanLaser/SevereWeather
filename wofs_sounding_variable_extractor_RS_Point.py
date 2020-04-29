###################################################################################################

#from mpl_toolkits.basemap import Basemap
#import matplotlib
import math
from scipy import *
from scipy import spatial
import pylab as P
import numpy as np
import sys
import os
import netCDF4
#import argparse
from optparse import OptionParser
from news_e_post_cbook import *
import glob
import pandas as pd
from netCDF4 import Dataset
####################################### File Variables: ######################################################
#optparse is deprecated, should switch to argparse
def sounding_extractor(time,date,Nfile):
	
	#Reading observation date and time in order to select closest forecast time for comparison.
	hours = ['0000','0030','0100','0130','0200','0230','0300','1900','1930','2000','2030','2100','2130','2200','2230','2300','2330']
	delt = np.ones(17)*np.nan
	
	if date == '20190518':
		date = '20190517'
	print(int(time))
	if 500 < int(time) < 1900:
			hour = '1900'
			t = 0
	if time == '0000':
		hour = '0000'
		t = 0
	else:
			
		for i in range(17):
			if int(time) > int(hours[i]):
				delt[i] = int(time)-int(hours[i])

		ind = np.nanargmin(delt[:])
		hour = hours[ind]
		if int(hour) == 30:
			hour = '0030'
		if int(hour) < 1000:
			hour = str(hour)
	
	# define t so closest forecast to obs time is selected
		interval = [00,05,10,15,20,25,30,35,40,45,50,55]
		delt = np.ones(12)*np.nan
		if hour[-2:] == '00':
			for i in range(12):
				delt[i] = abs(int(time[2:4])-interval[i])
			closeind = np.argmin(delt[:])
			t = interval[closeind]/5
	
		if hour[-2:] == '30':
			for i in range(12):
				delt[i] = abs(int(time[2:4])-interval[i])
			closeind = np.argmin(delt[:])
			t =  (interval[closeind]/5) - 6
	print(hour)
	indir = '/scratch/wof/realtime/FCST/' + str(date) + '/' + str(hour)
	outdir = '/home/jordan.laser/WoFS_Soundings/0517'
    
	###################################

	#read in sounding file and defining variables
	#radiosondefile = 'Far_Field_MW41_output_20190517_230230.csv'
	
	if Nfile[-3:] == 'csv':
		radiosondefile = str(Nfile)
		df = pd.read_csv(radiosondefile, skiprows = 3, sep = ',', header = None)
		lat = df[9][0]
		lon = df[10][0]
	if Nfile[-2:] == 'nc':
		lidarfile = Dataset(Nfile)
		lat = lidarfile.variables['lat'][0]
		lon = lidarfile.variables['lon'][0]

	#read in example WoFS file used for extracting grid point closest to observation
	directory = os.listdir(indir + '/ENS_MEM_10/')
	for d, dir in enumerate(directory):
		if (dir[0:6] == 'wrfwof'):
			firstfile = directory[d]
			print('File used for lat/lon: ' + str(firstfile))
			break		
			
	WoFS = Dataset(indir + '/ENS_MEM_10/' + firstfile)

	#Defining WoFS lat&lon for comparison 
	Wlat = WoFS.variables['XLAT'][0,:,:]
	Wlon = WoFS.variables['XLONG'][0,:,:]
	print(Wlat[0,0])
	print(lat)

	#Find WoFS grid point closest to sounding surface level
	Sdist = np.ones(Wlat.shape)*np.nan
	for i in range(Sdist.shape[0]):
		for j in range(Sdist.shape[1]):
			Sdist[i,j] = ((Wlat[i,j]-lat)**2+(Wlon[i,j]-lon)**2)**.5

	#Find index/position of nearest grid point	
	Smindex = np.unravel_index(np.argmin(Sdist, axis = None), Sdist.shape)
	Sxdex = Smindex[0]
	Sydex = Smindex[1]

	print('lat/lon indicies: ' + str(Sxdex) + ',' + str(Sydex))

	############################ Find wrfwof files to process: #################################

	### Find member dirs ### 

	ne = 18
	member_dirs = []

	member_dirs_temp = os.listdir(indir)
	#member_test = glob.glob(os.path.join(indir,'ENS*'))
	for d, dir in enumerate(member_dirs_temp):
	   if (dir[0:3] == 'ENS'):
		  member_dirs.append(dir)

	member_dirs.sort() #sorts as members [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 2, 3, 4, 5, 6, 7, 8, 9]

	print('member dirs',member_dirs,sorted(member_dirs))
	#print('test glob',member_test,sorted(member_test))
	files = []

	for n in range(0, len(member_dirs)):
	   temp_dir = os.path.join(indir, member_dirs[n])

	   member_files = []
	   temp_files = os.listdir(temp_dir)

	   for f, file in enumerate(temp_files):
		  if (file[0:9] == 'wrfwof_d0'):                               #assumes filename format of: "wrfwof_d02_yyyy-mm-dd_hh:mm:ss
			 member_files.append(file)
	   member_files.sort()

	   files.append(os.path.join(temp_dir, member_files[t]))

	files.sort()  #should have sorted directory paths to each ensemble file to be processed


	##################### Process wrfwof Files ########################
	edge = 7

	for f, infile in enumerate(files):

	   try:                                                 #open wrfwof file
		  fin = netCDF4.Dataset(infile, "r")
		  print("Opening %s \n" % infile)
	   except:
		  print("%s does not exist! \n" %infile)
		  sys.exit(1)

	   if (f == 0):

	##################### Get dimensions and attributes using first wrfwof file ########################

		  ### Get init and valid times ###
		  start_date = fin.START_DATE                               #Read initilization time string
		  init_year = start_date[0:4]
		  init_mon = start_date[5:7]
		  init_day = start_date[8:10]
		  init_hr = start_date[11:13]
		  init_min = start_date[14:16]

		  init_date = init_year + init_mon + init_day               #YYYYMMDD string for output file
		  init_time = init_hr + init_min                            #HHMM string for initialization time

		  valid_hr = infile[-8:-6]                          #Parse valid time from wrfwof filename
		  valid_min = infile[-5:-3]

		  valid_time = valid_hr + valid_min                 #HHMM string for valid time

		  ### Set output path ###
		  timestep = str(t)
		  if (len(timestep) == 1):
			 timestep = '0' + timestep
		  outname = "wofs_SND_" + timestep + "_" + init_date + "_" + init_time + "_" + valid_time + "_" + time + ".nc"         #output file
		  output_path = os.path.join(outdir,outname)
		  print('outname is: ' + outname)
		   

		  ### Get grid/projection info ### 

		  dx = fin.DX                                             #east-west grid spacing (m)
		  dy = fin.DY                                             #north-south grid spacing (m)
		  cen_lat = fin.CEN_LAT                                   #center of domain latitude (dec deg)
		  cen_lon = fin.CEN_LON                                   #center of domain longitude (dec deg)
		  stand_lon = fin.STAND_LON                               #center lon of Lambert conformal projection
		  true_lat_1 = fin.TRUELAT1                               #true lat value 1 for Lambert conformal conversion (dec deg)
		  true_lat_2 = fin.TRUELAT2                               #true lat value 2 for Lambert conformal conversion (dec deg)

		  xlat = fin.variables["XLAT"][0,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]                     #latitude (dec deg; Lambert conformal)
		  xlon = fin.variables["XLONG"][0,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]                    #longitude (dec deg; Lambert conformal)
		  hgt = fin.variables["HGT"][0,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]                       #terrain height above MSL (m)

		  nz = fin.dimensions['bottom_top'].size
		  ny, nx = xlat.shape

		  ### Calculate initial and valid time in seconds ### 

		  init_time_seconds = int(init_hr) * 3600. + int(init_min) * 60.
		  valid_time_seconds = int(valid_hr) * 3600. + int(valid_min) * 60.

		  if (valid_time_seconds < 43000.):         #Convert values past 0000 UTC, assumes forecast not run past ~12 UTC 
			 valid_time_seconds = valid_time_seconds + 86400.

	################################ Initialize output variables ########################################
		  #3D variables
		  u = np.zeros((ne,nz,ny,nx))               #U component of wind (m/s)
		  v = np.zeros((ne,nz,ny,nx))               #V component of wind (m/s)
		  tmp = np.zeros((ne,nz,ny,nx))             #temperature (K) -->writes to netcdf as degC
		  p = np.zeros((ne,nz,ny,nx))               #pressure (hPa)
		  z_agl = np.zeros((ne,nz,ny,nx))           #height AGL (m)
		  q = np.zeros((ne,nz,ny,nx))               #water vaport mixing ratio (kg/kg)
		  omeg = np.zeros((ne,nz,ny,nx))            #pressure vertical velocity (Pa/s)

	##################### Process wrfwof file: ########################

	############################### Read wrfwof variables: ##################################

	   uc = fin.variables["U"][0,:,:,:]		                #expects var dimensions of (nt, nz, ny, nx) with nt = 1
	   vc = fin.variables["V"][0,:,:,:]
	   wc = fin.variables["W"][0,:,:,:]

	   ph = fin.variables["PH"][0,:,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]
	   phb = fin.variables["PHB"][0,:,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]
	   pr = fin.variables["P"][0,:,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]
	   pb = fin.variables["PB"][0,:,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]

	   u[f,:,:,:] = ((uc[:,:,:-1]+uc[:,:,1:])/2.)[:,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]                          #convert staggered grids to centered
	   v[f,:,:,:] = ((vc[:,:-1,:]+vc[:,1:,:])/2.)[:,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]
	   w = ((wc[:-1,:,:]+wc[1:,:,:])/2.)[:,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]

	   qv = fin.variables["QVAPOR"][0,:,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]

	   q[f,:,:,:] = np.where(qv < 0., 0.0001, qv)  #force qv to be positive definite

	   th = fin.variables["T"][0,:,Sxdex-1:Sxdex+1,Sydex-1:Sydex+1]
	   th = th + 300.    	#add base state temp (300 K) to potential temp

	### Close wrfwof file ###

	   fin.close()
	   del fin

	########################## Calculate derived output variables (using news_e_post_cbook.py): ##############################

	######### Calculate vertical grid values #########
	   z, dz = calc_height(ph, phb)                            #height and layer thickness (m)
	   z_agl[f,:,:,:] = z - hgt
	   p[f,:,:,:] = (pr + pb) / 100.                           #pressure (convert Pa to hPa)
	   tmp[f,:,:,:] = calc_t(th, p[f,:,:,:])

	   #compute omega
	   R = 287.058
	   rho = (p[f,:,:,:] * 100.)/(R * (tmp[f,:,:,:]))
	   omeg[f,:,:,:] = -w * rho * 9.80665

	##################################################################
	##################### Write Summary File: ########################
	##################################################################
	### Create file and dimensions: ###
	try:
	   fout = netCDF4.Dataset(output_path, "w")
	except:
	   print("Could not create %s!\n" % output_path)

	fout.createDimension('NE', ne)
	fout.createDimension('NZ', nz)
	fout.createDimension('NX', nx)
	fout.createDimension('NY', ny)

	### Set Attributes: ###

	setattr(fout,'DX',dx)
	setattr(fout,'DY',dy)
	setattr(fout,'CEN_LAT',cen_lat)
	setattr(fout,'CEN_LON',cen_lon)
	setattr(fout,'STAND_LON',stand_lon)
	setattr(fout,'TRUE_LAT1',true_lat_1)
	setattr(fout,'TRUE_LAT2',true_lat_2)
	setattr(fout,'PROJECTION','Lambert Conformal')
	setattr(fout,'START_DATE',start_date)
	setattr(fout,'INIT_TIME_SECONDS',init_time_seconds)
	setattr(fout,'VALID_TIME_SECONDS',valid_time_seconds)
	setattr(fout,'FORECAST_TIME_STEP',t)

	### Create variables ###

	xlat1 = fout.createVariable('xlat', 'f4', ('NY','NX',))
	xlat1.long_name = "Latitude"
	xlat1.units = "degrees North"

	xlon1 = fout.createVariable('xlon', 'f4', ('NY','NX',))
	xlon1.long_name = "Longitude"
	xlon1.units = "degrees West"

	hgt1 = fout.createVariable('hgt', 'f4', ('NY','NX',))
	hgt1.long_name = "Height AGL"
	hgt1.units = "m"

	### 3D variables ###

	u_var = fout.createVariable('u', 'f4', ('NE','NZ','NY','NX',))
	u_var.long_name = "U-component of wind"
	u_var.units = "m s**-1"

	v_var = fout.createVariable('v', 'f4', ('NE','NZ','NY','NX',))
	v_var.long_name = "V-component of wind"
	v_var.units = "m s**-1"

	p_var = fout.createVariable('p', 'f4', ('NE','NZ','NY','NX',))
	p_var.long_name = "Pressure"
	p_var.units = "hPa"

	t_var = fout.createVariable('t', 'f4', ('NE','NZ','NY','NX',))
	t_var.long_name = "Temperature"
	t_var.units = "C"

	q_var = fout.createVariable('q', 'f4', ('NE','NZ','NY','NX',))
	q_var.long_name = "Water vapor mixing ratio"
	q_var.units = "Kg/Kg"

	zagl_var = fout.createVariable('z_agl', 'f4', ('NE','NZ','NY','NX',))
	zagl_var.long_name = "Height above ground level"
	zagl_var.units = "m"

	omeg_var = fout.createVariable('omega', 'f4', ('NE','NZ','NY','NX',))
	omeg_var.long_name = "Pressure vertical velocity"
	omeg_var.units = "Pa s**-1"

	### Write variables ###

	fout.variables['xlat'][:] = xlat
	fout.variables['xlon'][:] = xlon
	fout.variables['hgt'][:] = hgt

	fout.variables['u'][:] = u 
	fout.variables['v'][:] = v
	fout.variables['p'][:] = p
	fout.variables['t'][:] = tmp - 273.15 #K to C
	fout.variables['q'][:] = q
	fout.variables['z_agl'][:] = z_agl
	fout.variables['omega'][:] = omeg

	### Close output file ### 
	fout.close()
	del fout
	return outname
	